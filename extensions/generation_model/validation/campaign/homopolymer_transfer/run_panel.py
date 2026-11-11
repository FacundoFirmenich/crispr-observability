"""Score the frozen public human germline panel using indexed CRAM reads."""

import collections
import concurrent.futures
import gzip
import hashlib
import json
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path


ROOT = Path('/content/hg006_hg007_homopolymer')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def save(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')
    tmp.replace(path)


def distances(query, targets, edlib):
    return [edlib.align(query, t, mode='HW', task='distance')['editDistance'] for t in targets]


def anchored_run(read, reference, position, site):
    base=site['base']
    start=stop=position+1
    while start>0 and reference[start-1]==base:
        start-=1
    while stop<len(reference) and reference[stop]==base:
        stop+=1
    assert stop-start==site['ref_copies']
    left=anchor(read,start-1)
    right=anchor(read,stop)
    if left is None or right is None or left>=right:
        return 'unmapped_flank',None
    seq=read.query_sequence
    if seq[left]!=reference[start-1] or seq[right]!=reference[stop]:
        return 'flank_mismatch',None
    motif=seq[left+1:right]
    length=len(motif)
    if set(motif)-{base}:
        return 'mixed_segment',length
    if length==site['ref_copies']:
        return '0',length
    if length==site['alt_copies']:
        return '1',length
    return 'other_length',length


def anchor(read, position):
    q, r = 0, read.reference_start
    for op, length in read.cigartuples or ():
        if op in (0,7,8):
            if r<=position<r+length:
                return q+position-r
            q+=length
            r+=length
        elif op in (1,4):
            q+=length
        elif op in (2,3):
            if r<=position<r+length:
                return None
            r+=length
    return None


def call(scores):
    n=len(scores)
    if n<10 or 5*sum(min(s)<=59 for s in scores)<4*n:
        return 'U'
    candidates=[str(a) for a in (0,1) if 2*sum(s[a]<=s[1-a] for s in scores)>n]
    return candidates[0] if len(candidates)==1 else 'U'


def locus(cram, site, reference, edlib):
    position=site['pos']
    end=position+len(site['ref'])
    assert reference[position:end]==site['ref']
    targets=[reference[position-500:position]+site[key]+reference[end:end+500] for key in ('ref','alt')]
    eligible={}
    exclusions=collections.Counter()
    for read in cram.fetch('chr20',position,position+1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality<20:
            exclusions['flags_or_mapq']+=1
            continue
        q=anchor(read,position)
        seq=read.query_sequence
        if q is None or seq is None or not 150<=q<=len(seq)-150:
            exclusions['anchor_or_context']+=1
            continue
        if set(seq[q-150:q+150])-set('ACGT'):
            exclusions['non_acgt']+=1
            continue
        if read.query_name in eligible:
            exclusions['duplicate_id']+=1
            continue
        eligible[read.query_name]=(read,q)
    reads=[]
    for name in sorted(eligible,key=lambda s:digest(s.encode()))[:20]:
        read,q=eligible[name]
        query=read.query_sequence[q-75:q+75]
        score=distances(query,targets,edlib)
        assert len(query)==150 and min(score)>=0
        flank,copies=anchored_run(read,reference,position,site)
        reads.append({'id_sha256':digest(name.encode()),'scores':score,'flank':flank,'copies':copies,'reverse':read.is_reverse})
    primary=call([r['scores'] for r in reads])
    guarded=primary if primary!='U' and 2*sum(r['flank']==primary for r in reads)>len(reads) else 'U'
    return site|{'eligible':len(eligible),'exclusions':dict(exclusions),'primary':primary,'exact_majority':guarded,'reads':reads}


def prepare_index(source, path):
    with urllib.request.urlopen(source['crai'],timeout=40) as response:
        original=response.read(5_000_001)
    assert digest(original)==source['crai_sha256']
    delta=source['index_offset_delta']
    rows=[]
    for line in gzip.decompress(original).decode().splitlines():
        fields=line.split('\t')
        fields[3]=str(int(fields[3])+delta)
        rows.append('\t'.join(fields))
    data=gzip.compress(('\n'.join(rows)+'\n').encode(),mtime=0)
    path.write_bytes(data)
    return digest(data)


def summarize(rows):
    primary,guard,withdrawn=collections.Counter(),collections.Counter(),collections.Counter()
    for row in rows:
        def status(value):
            return 'unresolved' if value=='U' else 'concordant' if value==row['expected'] else 'discordant'
        before,after=status(row['primary']),status(row['exact_majority'])
        primary[before]+=1
        guard[after]+=1
        if before!='unresolved' and after=='unresolved':
            withdrawn[before]+=1
    return {'primary':dict(primary),'exact_majority':dict(guard),'withdrawn':dict(withdrawn),
            'aciertos_retirados':withdrawn['concordant']/primary['concordant'] if primary['concordant'] else None,
            'errores_retirados':withdrawn['discordant']/primary['discordant'] if primary['discordant'] else None}


def run_donor(donor, strata, source, reference, fasta, bindings, pysam, edlib):
    directory=ROOT/donor
    directory.mkdir(exist_ok=True)
    index=directory/'source.crai'
    index_hash=prepare_index(source,index)
    records={}
    ordered=sorted((site['pos'],key,site) for key,sites in strata.items() for site in sites)
    with pysam.AlignmentFile(source['cram'],'rc',index_filename=str(index),reference_filename=str(fasta)) as cram:
        sq=next(x for x in cram.header.to_dict()['SQ'] if x['SN']=='chr20')
        assert sq['LN']==len(reference) and sq['M5']==source['chr20']['M5']
        for i,(_,key,site) in enumerate(ordered,1):
            path=directory/f"{site['pos']}.json"
            if path.exists():
                prior=json.loads(path.read_text())
                assert prior['binding']==bindings
                row=prior['row']
            else:
                row=locus(cram,site,reference,edlib)
                save(path,{'binding':bindings,'row':row})
            records.setdefault(key,[]).append(row)
            if i==1 or i%8==0:
                print(json.dumps({'donor':donor,'completed':i,'eligible_last':row['eligible']}),flush=True)
    assert {key:len(rows) for key,rows in records.items()}=={key:len(rows) for key,rows in strata.items() if rows}
    results={'source':source['cram'],'derived_index_sha256':index_hash,'strata':records,
             'summary':{key:summarize(rows) for key,rows in records.items()},
             'length_summary':{key:summarize([r for s,rows in records.items() if s.endswith(key) for r in rows]) for key in ('ge15','lt15')},
             'all':summarize([r for rows in records.values() for r in rows])}
    save(directory/'RESULT.json',results)
    return results


def main(input_zip='/content/hg006_hg007_inputs.zip'):
    try:
        import pysam
        import edlib
    except ModuleNotFoundError:
        subprocess.run([sys.executable,'-m','pip','install','--quiet','pysam==0.23.3','edlib==1.3.9.post1'],check=True)
        import pysam
        import edlib
    ROOT.mkdir(exist_ok=True)
    with zipfile.ZipFile(input_zip) as archive:
        inventory=archive.read('inventory.json')
        panel=archive.read('panel.json')
        reference=gzip.decompress(archive.read('chr20.txt.gz')).decode()
        runner_hash=digest(archive.read('run_panel.py'))
        alias_hash=digest(archive.read('aliases.json')) if 'aliases.json' in archive.namelist() else None
    sources=json.loads(inventory)
    selection=json.loads(panel)
    assert selection['source_inventory_sha256']==digest(inventory)
    assert hashlib.md5(reference.encode()).hexdigest()=='b18e6c531b0bd70e949a7fc20859cb01'
    bindings={'inventory':digest(inventory),'panel':digest(panel),'runner':runner_hash,'reference_sha256':digest(reference.encode())}
    if 'alias_readback_sha256' in selection:
        assert selection['alias_readback_sha256']==alias_hash
        bindings['aliases']=alias_hash
    fasta=ROOT/'chr20.fa'
    if not fasta.exists():
        with fasta.open('w') as f:
            f.write('>chr20\n')
            for i in range(0,len(reference),60):
                f.write(reference[i:i+60]+'\n')
    if not fasta.with_suffix('.fa.fai').exists():
        pysam.faidx(str(fasta))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        futures={donor:executor.submit(run_donor,donor,entry['strata'],sources[donor],reference,fasta,bindings,pysam,edlib)
                 for donor,entry in selection['donors'].items()}
        results={donor:future.result() for donor,future in futures.items()}
    save(ROOT/'RESULT.json',{'binding':bindings,'donors':results})
    print(json.dumps({'status':'COMPLETE','summary':{d:{'all':r['all'],'length':r['length_summary']} for d,r in results.items()}}),flush=True)


if __name__=='__main__':
    main()
