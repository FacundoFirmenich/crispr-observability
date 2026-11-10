"""Read back the cluster matches and audit publication-side evidence carriers."""
from collections import Counter
import csv
import gzip
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile

from audit_fragments import HERE, SOURCES, digest, fastq, raw_path


def checksums():
    sources = {}
    for receipt in SOURCES.glob('*RECEIPT.json'):
        for name, meta in json.loads(receipt.read_bytes()).items():
            if 'sha256' in meta:
                assert (SOURCES/name).stat().st_size == meta['bytes']
                assert digest(SOURCES/name) == meta['sha256']
                sources[name] = meta['url']
    return sources


def check_control():
    base = HERE.parent/'CRISPR_P1_MP08_VALIDACION/fastq'
    count = 0
    for a,b in zip(*(fastq(base/f'SRR32414092_{e}.fastq.gz') for e in (1,2)),strict=True):
        assert a[0] == b[0] and a[1] and b[1] and a[1][0] == b[1][0]
        assert (a[1][1],b[1][1]) == (b'1',b'2')
        count += 1
    return count


def check_repeats(audit):
    expected = Counter()
    for row in audit['cluster_comparisons']:
        assert row['type'] == 'identical_same_end'
        a,b = row['run_a'],row['run_b']
        for x,y in zip(fastq(raw_path(a)),fastq(raw_path(b)),strict=True):
            assert x[1] and x[1] == y[1] and x[2:] == y[2:]
            expected[(a,b,x[0],y[0],x[1][0].decode())] += 1
        assert sum(n for key,n in expected.items() if key[:2] == (a,b)) == row['n']
    ledger = HERE/'CLUSTER_LINKS.csv.gz'
    assert digest(ledger) == audit['cluster_ledger_sha256']
    with gzip.open(ledger,'rt',newline='') as f:
        for r in csv.DictReader(f):
            assert r['end_a'] == r['end_b'] == '1' and r['comparison'] == 'identical_same_end'
            key=tuple(r[k] for k in ('run_a','run_b','read_a','read_b','cluster'))
            assert expected.pop(key) == 1
    assert not expected


def workbook_inventory():
    workbooks = {}
    with ZipFile(SOURCES/'Sousa_supplementary.zip') as supplement:
        assert supplement.testzip() is None
        for name in sorted(n for n in supplement.namelist() if n.endswith('.xlsx')):
            with ZipFile(BytesIO(supplement.read(name))) as book:
                assert book.testzip() is None
                entries = ET.fromstring(book.read('xl/workbook.xml')).findall('.//{*}sheet')
                relations = {x.attrib['Id']:x.attrib['Target'] for x in
                             ET.fromstring(book.read('xl/_rels/workbook.xml.rels'))}
                sheets = {}
                for entry in entries:
                    dest = relations[entry.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
                    root = ET.fromstring(book.read(dest.lstrip('/') if dest.startswith('/') else 'xl/'+dest))
                    cells = root.findall('.//{*}c')
                    sheets[entry.attrib['name']] = {'rows':len(root.findall('.//{*}row')),'cells':len(cells)}
                shared = (ET.fromstring(book.read('xl/sharedStrings.xml'))
                          if 'xl/sharedStrings.xml' in book.namelist() else None)
                text = ' '.join(shared.itertext()) if shared is not None else ''
                for worksheet in (n for n in book.namelist() if n.startswith('xl/worksheets/sheet') and n.endswith('.xml')):
                    text += ' '.join(ET.fromstring(book.read(worksheet)).itertext())
                assert re.search(r'SRR\d+\.\d+',text) is None
                assert re.search(r'\w+:\d+:[A-Za-z0-9-]+:\d+:\d+:\d+:\d+',text) is None
                workbooks[name]={'sheets':sheets,'archive_record_or_cluster_identifiers_detected':0}
    return workbooks


def main():
    audit=json.loads((HERE/'FRAGMENT_AUDIT.json').read_bytes())
    check_repeats(audit)
    control=check_control()
    catalog=json.loads((SOURCES/'CRISPR-A_PLOS_figshare.json').read_bytes())
    assert len(catalog)==29 and {r['resource_doi'] for r in catalog}=={'10.1371/journal.pcbi.1011137'}
    figshare=json.loads((SOURCES/'CRISPR-A_S5_figshare.json').read_bytes())
    files=figshare['files']
    assert len(files)==1 and files[0]['name']=='pcbi.1011137.s005.jpg'
    assert digest(SOURCES/'CRISPR-A_S5.jpg','md5')==files[0]['computed_md5']
    with ZipFile(SOURCES/'CRISPR-A_Bitbucket.zip') as repo:
        assert repo.testzip() is None
        paths=[name for name in repo.namelist() if not name.endswith('/')]
        candidates=[name for name in paths if re.search(r'10285861|validat|\.xlsx?$|\.ods$',name,re.I)]
    assert not candidates
    output={'verified_source_receipt_files':checksums(),'source_workbooks':workbook_inventory(),
            'CRISPR_A_S5':{'figshare':figshare['url_public_html'],'counts':{'WT':328,'TMPL':49,'INDELS':6},
                           'records_linkable_from_table':0,'pinned_pipeline_files':len(paths),
                           'publisher_deposited_items':len(catalog)},
            'paired_parser_positive_control':{'run':'SRR32414092','verified_pairs':control,
                                             'R1_sha256':digest(HERE.parent/'CRISPR_P1_MP08_VALIDACION/fastq/SRR32414092_1.fastq.gz'),
                                             'R2_sha256':digest(HERE.parent/'CRISPR_P1_MP08_VALIDACION/fastq/SRR32414092_2.fastq.gz')},
            'same_end_rows_read_back':sum(r['n'] for r in audit['cluster_comparisons'])}
    (HERE/'READBACK_VERIFICATION.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(control,output['same_end_rows_read_back'],output['CRISPR_A_S5'])


if __name__=='__main__':
    main()
