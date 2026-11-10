"""Exploratory library audit; not a genotype or P1 accuracy validation.
Selection frozen before reading: Clone10, four DNA/cDNA mate files.
All raw records retained in denominators. Unique exact overlap >=30 bases,
ACGT-only; no mismatch tolerance, no tuning, no reference-derived labels.
Raw downloads stay in RAM; published MD5 and size checked before parsing.
"""
import urllib.request, json, gzip, hashlib, collections, pathlib, datetime

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / 'COL17A1_CLONE10_AUDIT_20261105.json'
def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()
def rc(s):
    return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def parse(b):
    lines=gzip.decompress(b).decode('ascii').splitlines()
    assert len(lines)%4==0
    out=[]
    for i in range(0,len(lines),4):
        h,s,p,q=lines[i:i+4]
        assert h.startswith('@') and p.startswith('+') and len(s)==len(q)
        out.append((h.split()[0].removesuffix('/1').removesuffix('/2'),s))
    assert len({h for h,s in out})==len(out)
    return out
def main():
    meta=json.loads(fetch('https://zenodo.org/api/records/6242279'))
    selected=sorted([f for f in meta['files'] if 'Clone10_' in f['key']],key=lambda f:f['key'])
    assert len(selected)==4
    result={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'scope':'Exploratory exact-overlap audit, one named clone; no P1 validation or independent molecular truth',
      'source':'https://zenodo.org/records/6242279','protocol':__doc__,'files':[],'libraries':{}}
    reads={}
    for f in selected:
        b=fetch(f['links']['self'])
        assert len(b)==f['size'] and 'md5:'+hashlib.md5(b).hexdigest()==f['checksum']
        reads[f['key']]=parse(b)
        result['files'].append(dict(f,sha256=hashlib.sha256(b).hexdigest(),records=len(reads[f['key']])))
        print('VERIFIED',f['key'],len(reads[f['key']]),flush=True)
    consensuses={}
    for typ in ['DNA','cDNA']:
        a=reads[f'C17Ex52_{typ} immoJEB_PPN_Clone10_fwd.fastq.gz']
        b=reads[f'C17Ex52_{typ} immoJEB_PPN_Clone10_rev.fastq.gz']
        assert len(a)==len(b)
        statuses=collections.Counter(); seqs=collections.Counter(); lengths=collections.Counter()
        for (h,s),(h2,t) in zip(a,b):
            assert h==h2
            t=rc(t)
            if set(s+t)-set('ACGT'):
                statuses['non_ACGT']+=1; continue
            ks=[k for k in range(30,min(len(s),len(t))+1) if s[-k:]==t[:k]]
            if len(ks)==1:
                statuses['unique_exact_overlap']+=1
                seqs[s+t[ks[0]:]]+=1; lengths[ks[0]]+=1
            else: statuses['no_exact_overlap' if not ks else 'multiple_exact_overlaps']+=1
        assert sum(statuses.values())==len(a)
        consensuses[typ]=seqs
        result['libraries'][typ]={'raw_pairs':len(a),'statuses':dict(statuses),'overlap_lengths':dict(lengths),'merged_sequences':seqs.most_common()}
    common=set(consensuses['DNA'])&set(consensuses['cDNA'])
    result['cross_library']={'shared_exact_full_merged_sequences':len(common),
       'DNA_pairs_in_shared_sequences':sum(consensuses['DNA'][s] for s in common),
       'cDNA_pairs_in_shared_sequences':sum(consensuses['cDNA'][s] for s in common),
       'caution':'Full amplicon equality may be inappropriate across DNA/cDNA primers; absence is not biological discordance.'}
    # Exploratory extension chosen AFTER inspecting the initial top sequences.
    # This is explicitly post hoc and must not be described as held-out validation.
    left,right='GATGCCTCCCACAGTC','ATGTCAGGCGGGGCAG'
    internal={'status':'EXPLORATORY_POSTHOC_ANCHORS_NOT_P1_VALIDATION',
      'left_anchor':left,'right_anchor':right,'libraries':{}}
    for typ,seqs in consensuses.items():
        c=collections.Counter(); rejected=0
        for s,n in seqs.items():
            if s.count(left)==1 and s.count(right)==1 and s.index(right)>=s.index(left)+len(left):
                c[s[s.index(left)+len(left):s.index(right)]]+=n
            else: rejected+=n
        internal['libraries'][typ]={'raw_pairs':result['libraries'][typ]['raw_pairs'],
          'admitted_internal':sum(c.values()),'merged_anchor_rejected':rejected,'internal_sequences':c.most_common()}
    result['exploratory_internal_comparison']=internal
    OUT.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:{a:b for a,b in v.items() if a!='merged_sequences'} for k,v in result['libraries'].items()}))
    print(json.dumps(result['cross_library']))
    for typ,c in consensuses.items(): print(typ,'TOP3',c.most_common(3))
if __name__=='__main__': main()
