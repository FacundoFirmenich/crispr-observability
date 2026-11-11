"""Freeze a human germline panel before accessing its sequencing outcomes."""

import argparse
import bisect
import collections
import gzip
import hashlib
import json
import sys
import zipfile
from pathlib import Path


def eligible(truth, hg002, bed, reference, old, repeat_interval, simple, hom_alt):
    bed_starts = [lo for lo,_ in bed]
    positions = [int(row[1])-1 for row in truth]
    assert positions==sorted(positions)
    ends=[]
    max_end=0
    for pos,row in zip(positions,truth):
        max_end=max(max_end,pos+len(row[3]))
        ends.append(max_end)
    pools=collections.defaultdict(list)
    for label,records in (('1',truth),('0',hg002)):
        seen=set()
        for row in records:
            pos=int(row[1])-1
            ref,alt=row[3:5]
            if pos in old or pos in seen or not simple(ref,alt) or max(len(ref),len(alt))!=2 or ref[0]!=alt[0] or not hom_alt(row):
                continue
            lo,hi=pos-500,pos+len(ref)+500
            idx=bisect.bisect_right(bed_starts,lo)-1
            if idx<0 or hi>bed[idx][1] or reference[pos:pos+len(ref)]!=ref or set(reference[lo:hi])-set('ACGT'):
                continue
            left,right=bisect.bisect_right(ends,lo),bisect.bisect_left(positions,hi)
            if any(positions[j]+len(truth[j][3])>lo and (label=='0' or truth[j] is not row) for j in range(left,right)):
                continue
            base,_,_,copies_ref,copies_alt=repeat_interval(reference,pos,ref,alt)
            stratum=label+('_ge15' if max(copies_ref,copies_alt)>=15 else '_lt15')
            pools[stratum].append({'pos':pos,'ref':ref,'alt':alt,'expected':label,'base':base,
                                   'ref_copies':copies_ref,'alt_copies':copies_alt,'block':pos//1_000_000})
            seen.add(pos)
    return pools


def rank(sample,stratum,row):
    key=f"HG00X-homopolymer-v1|{sample}|{stratum}|{row['pos']}|{row['ref']}|{row['alt']}"
    return hashlib.sha256(key.encode()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workspace',type=Path,required=True)
    parser.add_argument('--inventory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--aliases',type=Path)
    parser.add_argument('--exclude-panel',type=Path)
    args=parser.parse_args()
    sys.path.insert(0,str(args.workspace.resolve()))
    import run_hg002_expanded_position_gate as base
    from audit_hg002_single_base_runs import repeat_interval

    inventory_bytes=args.inventory.read_bytes()
    inventory=json.loads(inventory_bytes)
    aliases=json.loads(args.aliases.read_bytes()) if args.aliases else None
    if aliases:
        assert aliases['source_inventory_sha256']==base.sha(inventory_bytes)
    old,_=base.old_records()
    for name in ('HG002_CHR20_POSITION_GATE.json','HG002_POSITION_GATE.json'):
        data=json.loads((args.workspace/name).read_text(encoding='utf-8'))
        old.update(row['pos1']-1 for row in data.get('loci',()))
        for platform in data.get('platforms',{}).values():
            old.update(row['position_1_based']-1 for row in platform['loci'])
    with zipfile.ZipFile(base.INPUT/'MP36_COMPLETE_AUDITED.zip') as archive:
        reference=gzip.decompress(archive.read('inputs/chr20_sequence.txt.gz')).decode()
        if aliases:
            catalog=[line.split('\t') for line in archive.read('inputs/truth_chr20.vcf').decode().splitlines()
                     if line.startswith('chr20\t')]
    assert hashlib.md5(reference.encode()).hexdigest()=='b18e6c531b0bd70e949a7fc20859cb01'
    if not aliases:
        catalog,_=base.vcf_chr20()
    if args.exclude_panel:
        previous=json.loads(args.exclude_panel.read_bytes())
        old.update(row['pos'] for entry in previous['donors'].values() for rows in entry['strata'].values() for row in rows)
    result={'scope':'Genotype-only, before sequencing outcomes. Coordinates are zero-based. One locus per 1 Mb block per stratum.',
            'source_inventory_sha256':base.sha(inventory_bytes),'donors':{}}
    used=set(old)
    if aliases:
        result['alias_readback_sha256']=base.sha(args.aliases.read_bytes())
        result['excluded_panel_sha256']=base.sha(args.exclude_panel.read_bytes())
        result['reference_candidate_catalog']='HG005 GIAB v4.2.1 chr20, archived in MP36_COMPLETE_AUDITED.zip'
        result['quota']='Up to 16 per stratum; keep all eligible blocks if fewer, without substitution.'
    for donor in (('HG003','HG004') if aliases else ('HG006','HG007')):
        source=inventory[donor]
        if aliases:
            proof=aliases['donors'][donor]
            assert proof['status']=='PUBLIC_RUN_LINK_VERIFIED' and proof['cram']==source['cram']
        else:
            assert source['status']=='HEADER_INDEX_AND_ACQUISITION_VERIFIED'
        base.CACHE=args.workspace/f'{donor}_HOMOPOLYMER_TRUTH'
        base.TRUTH=source['truth']['vcf']['url'].removesuffix('.vcf.gz')
        base.CACHE.mkdir(exist_ok=True)
        if not aliases:
            assert all((base.CACHE/n).exists() for n in ('GIAB.vcf.gz.tbi','GIAB_chr20.bgzf','GIAB_confidence.bed'))
        truth,_=base.vcf_chr20()
        base.cached('GIAB_confidence.bed',source['truth']['bed']['url'],25_000_000)
        bed=sorted((int(f[1]),int(f[2])) for line in (base.CACHE/'GIAB_confidence.bed').read_text().splitlines()
                   if (f:=line.split()) and f[0]=='chr20')
        pools=eligible(truth,catalog,bed,reference,old,repeat_interval,base.is_simple,base.hom_alt)
        print(donor,{key:len({r['block'] for r in rows if r['pos'] not in used}) for key,rows in pools.items()},flush=True)
        selected={}
        for stratum in ('0_ge15','1_ge15','0_lt15','1_lt15'):
            blocks=set()
            panel=[]
            for row in sorted(pools[stratum],key=lambda r:rank(donor,stratum,r)):
                if row['pos'] in used or row['block'] in blocks:
                    continue
                panel.append(row)
                blocks.add(row['block'])
                used.add(row['pos'])
                if len(panel)==16:
                    break
            if not aliases:
                assert len(panel)==16,(donor,stratum,len(panel))
            selected[stratum]=panel
        result['donors'][donor]={'eligible_before_cross_donor_exclusion':{k:len(v) for k,v in pools.items()},'strata':selected}
        print(donor,'selected',{key:len(rows) for key,rows in selected.items()},flush=True)
    encoded=json.dumps(result,indent=2)+'\n'
    if args.output.exists():
        assert args.output.read_text(encoding='utf-8')==encoded,'Frozen selection changed'
    else:
        args.output.write_text(encoded,encoding='utf-8')


if __name__=='__main__':
    main()
