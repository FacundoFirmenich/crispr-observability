"""Verify public human germline read sources and genotype reference files."""

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import struct
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


S3 = 'https://ont-open-data.s3.amazonaws.com/'
NIST = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/'
DONORS = {
    'HG003': ('AshkenazimTrio/HG003_NA24149_father', 'PAY87794', '_noinconsistent.bed'),
    'HG004': ('AshkenazimTrio/HG004_NA24143_mother', 'PAY87778', '_noinconsistent.bed'),
    'HG006': ('ChineseTrio/HG006_NA24694_father', 'PAY77227', '.bed'),
    'HG007': ('ChineseTrio/HG007_NA24695_mother', 'PAY12990', '.bed'),
}


def fetch(url, limit, span=None, head=False):
    headers = {} if span is None else {'Range':f'bytes={span[0]}-{span[1]}'}
    request = urllib.request.Request(url, headers=headers, method='HEAD' if head else 'GET')
    with urllib.request.urlopen(request,timeout=20) as r:
        if span is not None:
            assert r.status==206 and r.headers['Content-Range'].startswith(f'bytes {span[0]}-{span[1]}/')
        data = r.read(limit+1)
        assert len(data)<=limit and (span is None or len(data)==span[1]-span[0]+1)
        return data, dict(r.headers)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def itf8(data, pos):
    first = data[pos]
    pos += 1
    n = 0 if first<128 else 1 if first<192 else 2 if first<224 else 3 if first<240 else 4
    value = first & (127>>n)
    for i in range(n):
        bits = 4 if n==4 and i==3 else 8
        value = (value<<bits) | (data[pos]&((1<<bits)-1))
        pos += 1
    return value,pos


def cram_dictionary(data):
    assert data[:6]==b'CRAM\x03\x00'
    pos=30
    for _ in range(4):
        _,pos=itf8(data,pos)
    assert data[pos:pos+2]==b'\x00\x00'
    _,pos=itf8(data,pos+2)
    landmarks,pos=itf8(data,pos)
    for _ in range(landmarks):
        _,pos=itf8(data,pos)
    pos+=4
    first_container = pos+struct.unpack_from('<i',data,26)[0]
    method,kind=data[pos:pos+2]
    assert kind==0
    _,pos=itf8(data,pos+2)
    packed,pos=itf8(data,pos)
    unpacked,pos=itf8(data,pos)
    block=data[pos:pos+packed]
    sam=gzip.decompress(block) if method==1 else block
    assert method in (0,1) and len(sam)==unpacked
    n=struct.unpack_from('<i',sam)[0]
    lines=sam[4:4+n].decode().splitlines()
    sq=[dict(field.split(':',1) for field in line.split('\t')[1:]) for line in lines if line.startswith('@SQ\t')]
    ref_id=next(i for i,row in enumerate(sq) if row['SN']=='chr20')
    return first_container,ref_id,sq[ref_id]


def inspect(donor, config):
    folder, flowcell, bed_suffix = config
    cram = S3+f'giab_2025.01/analysis/wf-human-variation/sup/{donor}/{flowcell}/output/SAMPLE.haplotagged.cram'
    prefix, h = fetch(cram,65536,(0,65535))
    total=int(h['Content-Range'].split('/')[-1])
    first,ref_id,sq=cram_dictionary(prefix)
    assert sq['LN']=='64444167' and sq['M5']=='b18e6c531b0bd70e949a7fc20859cb01'
    index, _ = fetch(cram+'.crai',5_000_000)
    rows=[tuple(map(int,line.split())) for line in gzip.decompress(index).splitlines()]
    offset=first-min(row[3] for row in rows)
    selected=sorted({row[3] for row in rows if row[0]==ref_id})
    assert selected
    checks=[]
    physical=sorted({row[3] for row in rows})
    for published in (selected[0],selected[len(selected)//2],selected[-1]):
        b,_=fetch(cram,64,(published+offset,published+offset+63))
        i=physical.index(published)
        length=struct.unpack_from('<i',b)[0]
        gap=physical[i+1]-published if i+1<len(physical) else None
        assert length>0 and (gap is None or 0<gap-length<256)
        checks.append({'original_offset':published,'physical_offset':published+offset,'sha256':sha(b)})
    flow_prefix=f'giab_2025.01/flowcells/{donor}/{flowcell}/'
    listing,_=fetch(S3+'?'+urllib.parse.urlencode({'list-type':'2','prefix':flow_prefix,'delimiter':'/'}),300000)
    root=ET.fromstring(listing)
    summaries=[key.text for key in root.findall('.//{*}Key') if '/final_summary_' in key.text]
    assert len(summaries)==1
    summary_url=S3+summaries[0]
    content,_=fetch(summary_url,30000)
    acquisition=dict(line.split('=',1) for line in content.decode().splitlines() if '=' in line)
    assert acquisition['flow_cell_id']==flowcell
    matches_sample=acquisition['sample_id'].upper()==donor
    base=NIST+folder+f'/NISTv4.2.1/GRCh38/{donor}_GRCh38_1_22_v4.2.1_benchmark'
    truth={}
    for key,suffix in (('vcf','.vcf.gz'),('tabix','.vcf.gz.tbi'),('bed',bed_suffix)):
        url=base+suffix
        _,headers=fetch(url,0,head=True)
        truth[key]={'url':url,'bytes':int(headers['Content-Length']) if 'Content-Length' in headers else None}
    status='HEADER_INDEX_AND_ACQUISITION_VERIFIED' if matches_sample else 'HOLD_SAMPLE_ALIAS'
    return {'status':status,'cram':cram,'cram_bytes':total,'cram_etag':h.get('ETag'),
            'header_prefix_sha256':sha(prefix),'chr20':sq,'crai':cram+'.crai','crai_sha256':sha(index),
            'index_offset_delta':offset,'chr20_container_checks':checks,'flowcell':flowcell,
            'acquisition_url':summary_url,'acquisition_sha256':sha(content),'acquisition':acquisition,'truth':truth}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    result={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures={donor:executor.submit(inspect,donor,config) for donor,config in DONORS.items()}
        for donor,future in futures.items():
            try:
                result[donor]=future.result()
            except Exception as error:
                result[donor]={'status':'PROBE_FAILED','error':repr(error)}
            print(donor,result[donor]['status'],result[donor].get('error',''),flush=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return int(any(v['status']=='PROBE_FAILED' for v in result.values()))


if __name__=='__main__':
    raise SystemExit(main())
