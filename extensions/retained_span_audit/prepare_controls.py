"""Recover empirical control references from the fixed public Zenodo version."""
import pathlib,json,gzip,hashlib,collections
from col17a1_clone10_audit import fetch
ROOT=pathlib.Path(__file__).resolve().parent
meta=json.loads(fetch('https://zenodo.org/api/records/6242279'))
names={f'C17Ex52_DNA_immo{x}_Ctl_fwd.fastq.gz' for x in ['HKC','JEB']}
selected=[f for f in meta['files'] if f['key'] in names]
assert len(selected)==2
for f in selected:
    b=fetch(f['links']['self'])
    assert len(b)==f['size'] and 'md5:'+hashlib.md5(b).hexdigest()==f['checksum']
    lines=gzip.decompress(b).decode('ascii').splitlines();assert len(lines)%4==0
    for i in range(0,len(lines),4):
        assert lines[i].startswith('@') and lines[i+2].startswith('+') and len(lines[i+1])==len(lines[i+3])
    seqs=lines[1::4];counts=collections.Counter(seqs)
    result={'source':f,'sha256':hashlib.sha256(b).hexdigest(),'raw_reads':len(seqs),'top':counts.most_common(5)}
    (ROOT/(f['key']+'.top.json')).write_text(json.dumps(result),encoding='utf-8')
    print(f['key'],len(seqs))
