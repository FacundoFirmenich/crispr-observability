"""Compare HiFi with Nanopore on the same frozen HG003/HG004 loci."""

import concurrent.futures
import gzip
import hashlib
import json
import urllib.request
import zipfile
from pathlib import Path

import edlib
import pysam


ROOT=Path('/content/hg003_hg004_hifi')
BASE='https://storage.googleapis.com/deepvariant/pacbio-case-study-testdata/'


def run(donor, strata, reference, core):
    source=BASE+f'{donor}.pfda_challenge.grch38.phased.chr20.bam'
    with urllib.request.urlopen(urllib.request.Request(source,headers={'Range':'bytes=0-63'}),timeout=30) as response:
        assert response.status==206 and response.headers['Content-Range'].startswith('bytes 0-63/')
        prefix=response.read(65)
        assert len(prefix)==64
        source_headers=dict(response.headers)
    with urllib.request.urlopen(source+'.bai',timeout=40) as response:
        index=response.read(5_000_001)
    assert len(index)<=5_000_000 and index[:4]==b'BAI\1'
    directory=ROOT/donor
    directory.mkdir(exist_ok=True)
    index_path=directory/'source.bai'
    index_path.write_bytes(index)
    rows={}
    with pysam.AlignmentFile(source,'rb',index_filename=str(index_path)) as reads:
        header=reads.header.to_dict()
        sq=next(item for item in header['SQ'] if item['SN']=='chr20')
        assert sq['LN']==len(reference)
        if 'M5' in sq:
            assert sq['M5']==hashlib.md5(reference.encode()).hexdigest()
        aliases={'HG003':{'HG003','NA24149'},'HG004':{'HG004','NA24143'}}
        samples={rg['SM'].upper() for rg in header['RG']}
        assert samples and samples<=aliases[donor],samples
        ordered=sorted((s['pos'],key,s) for key,loci in strata.items() for s in loci)
        for i,(_,key,site) in enumerate(ordered,1):
            result=core['locus'](reads,site,reference,edlib)
            rows.setdefault(key,[]).append(result)
            if i==1 or i%8==0 or i==len(ordered):
                print(json.dumps({'platform':'HiFi','donor':donor,'completed':i}),flush=True)
    assert {k:len(v) for k,v in rows.items()}=={k:len(v) for k,v in strata.items()}
    return {'source':source,'source_headers':source_headers,'prefix_sha256':core['digest'](prefix),
            'index_sha256':core['digest'](index),'chr20':sq,'samples':sorted(samples),'strata':rows,
            'summary':{key:core['summarize'](loci) for key,loci in rows.items()}}


def main():
    ROOT.mkdir(exist_ok=True)
    with zipfile.ZipFile('/content/hg003_hg004_inputs.zip') as z:
        reference=gzip.decompress(z.read('chr20.txt.gz')).decode()
        panel=z.read('panel.json')
        runner=z.read('run_panel.py')
    assert hashlib.md5(reference.encode()).hexdigest()=='b18e6c531b0bd70e949a7fc20859cb01'
    core={'__name__':'frozen_runner'}
    exec(compile(runner,'run_panel.py','exec'),core)
    selection=json.loads(panel)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        futures={d:ex.submit(run,d,s['strata'],reference,core) for d,s in selection['donors'].items()}
        donors={d:f.result() for d,f in futures.items()}
    core['save'](ROOT/'RESULT.json',{'panel_sha256':core['digest'](panel),'runner_sha256':core['digest'](runner),'donors':donors})
    print(json.dumps({'status':'COMPLETE','summary':{d:r['summary'] for d,r in donors.items()}}),flush=True)


if __name__=='__main__':
    main()
