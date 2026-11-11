"""Link public acquisition aliases to the sample in the CRAM read group."""

import argparse
import csv
import datetime
import gzip
import io
import json
import struct
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

from discover import S3, fetch, itf8, sha


def read_group(header):
    pos=30
    for _ in range(4):
        _,pos=itf8(header,pos)
    _,pos=itf8(header,pos+2)
    n,pos=itf8(header,pos)
    for _ in range(n):
        _,pos=itf8(header,pos)
    pos+=4
    method=header[pos]
    _,pos=itf8(header,pos+2)
    packed,pos=itf8(header,pos)
    _,pos=itf8(header,pos)
    data=header[pos:pos+packed]
    sam=gzip.decompress(data) if method==1 else data
    lines=sam[4:4+struct.unpack_from('<i',sam)[0]].decode().splitlines()
    groups=[line for line in lines if line.startswith('@RG\t')]
    assert len(groups)==1
    return dict(x.split(':',1) for x in groups[0].split('\t')[1:])


def verify(sample, source):
    data,_=fetch(source['cram'],65536,(0,65535))
    assert sha(data)==source['header_prefix_sha256']
    group=read_group(data)
    acquisition=source['acquisition']
    run=group['ID'].split('_')[0]
    flow=group['PU']
    assert group['SM'].upper()==sample and flow==acquisition['flow_cell_id'] and run==acquisition['acquisition_run_id']
    parse=datetime.datetime.fromisoformat
    assert abs((parse(group['DT'])-parse(acquisition['started'])).total_seconds())<0.01
    prefix=f'giab_2025.01/flowcells/{sample}/{flow}/'
    listing,_=fetch(S3+'?'+urllib.parse.urlencode({'list-type':'2','prefix':prefix,'delimiter':'/'}),300000)
    keys=[k.text for k in ET.fromstring(listing).findall('.//{*}Key') if '/sample_sheet_' in k.text]
    assert len(keys)==1
    sheet,_=fetch(S3+keys[0],10000)
    rows=list(csv.DictReader(io.StringIO(sheet.decode('utf-8-sig'))))
    alias=acquisition['sample_id']
    assert any(alias in row.values() for row in rows),rows
    return {'status':'PUBLIC_RUN_LINK_VERIFIED','sample':group['SM'],'alias':alias,'run_id':run,'flowcell':flow,
            'read_group_time':group['DT'],'sample_sheet_url':S3+keys[0],'sample_sheet_sha256':sha(sheet),
            'cram':source['cram'],'summary':source['acquisition_url']}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('inventory',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    data=args.inventory.read_bytes()
    inventory=json.loads(data)
    donors={sample:verify(sample,inventory[sample]) for sample in ('HG003','HG004')}
    args.output.write_text(json.dumps({'source_inventory_sha256':sha(data),'donors':donors},indent=2)+'\n',encoding='utf-8')
    for sample,row in donors.items():
        print(sample,row['alias'],row['flowcell'],row['status'])


if __name__=='__main__':
    main()
