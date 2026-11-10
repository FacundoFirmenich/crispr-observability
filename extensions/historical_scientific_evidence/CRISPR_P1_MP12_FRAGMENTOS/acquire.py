"""Retrieve original study metadata and accession-level evidence privately."""
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parent / 'sources'
ENA = 'https://www.ebi.ac.uk/ena/'
FIELDS = ('run_accession,experiment_accession,sample_accession,secondary_sample_accession,'
          'library_layout,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,'
          'submitted_ftp,experiment_title,run_alias')
HEK3 = ('SRR10285860', 'SRR10285861', 'SRR10285862')
CFTR = tuple(f'SRR296232{x}' for x in (65, 66, 67, 68, 70, 71, 72, 73, 74))
CFTR_EXTRA = ('SRR29623656','SRR29622959','SRR29622955','SRR29622845','SRR29622832','SRR29622827')


def report(project):
    return ENA+'portal/api/filereport?'+urlencode(dict(accession=project, result='read_run',
                                                        fields=FIELDS, format='tsv'))


def retrieve(item):
    name, url = item
    try:
        with urlopen(Request(url, headers={'User-Agent': 'CRISPR-P1-fragment-provenance/1.0'}), timeout=25) as f:
            body = f.read()
            final_url = f.geturl()
        if not body:
            raise ValueError('Empty response')
        (ROOT/name).write_bytes(body)
        return name, dict(url=url, retrieved_url=final_url, bytes=len(body),
                          sha256=hashlib.sha256(body).hexdigest())
    except (OSError, ValueError) as e:
        return name, dict(url=url, error=f'{type(e).__name__}: {e}')


def main():
    ROOT.mkdir(exist_ok=True)
    sources = {p+'.tsv': report(p) for p in ('PRJNA565979','PRJNA1055086')}
    sources.update({r+'.xml': ENA+'browser/api/xml/'+r for r in HEK3+CFTR+CFTR_EXTRA})
    sources.update({
        'Anzalone.xml': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6907074/fullTextXML',
        'Sousa.xml': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11754097/fullTextXML',
        'CRISPR-A.xml': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10256225/fullTextXML',
        'CRISPR-A_S5.jpg': 'https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1011137.s005&type=supplementary',
        'CFTR_Supplementary_Data.xlsx': 'https://static-content.springer-cdn.net/esm/art%3A10.1038%2Fs41551-024-01233-3/MediaObjects/41551_2024_1233_MOESM3_ESM.xlsx',
    })
    with ThreadPoolExecutor(max_workers=8) as pool:
        receipts = dict(pool.map(retrieve, sources.items()))
        experiments = {}
        samples = set()
        for project, wanted in [('PRJNA565979', HEK3), ('PRJNA1055086', CFTR+CFTR_EXTRA)]:
            name = project+'.tsv'
            if 'error' in receipts[name]:
                continue
            rows = list(csv.DictReader(io.StringIO((ROOT/name).read_text(encoding='utf-8')), delimiter='\t'))
            assert len({r['run_accession'] for r in rows}) == len(rows)
            selected = [r for r in rows if r['run_accession'] in wanted]
            assert len(selected) == len(wanted), (project, len(selected))
            experiments.update({row['experiment_accession']: ENA+'browser/api/xml/'+row['experiment_accession']
                                for row in selected})
            samples.update(row['sample_accession'] for row in selected)
        receipts.update(pool.map(retrieve, ((acc+'.xml', url) for acc,url in experiments.items())))
        receipts.update(pool.map(retrieve, ((acc+'.xml', ENA+'browser/api/xml/'+acc) for acc in sorted(samples))))
    (ROOT/'SOURCE_RECEIPT.json').write_text(json.dumps(receipts,indent=2)+'\n',encoding='utf-8')
    print(len(receipts), {n:r['error'] for n,r in receipts.items() if 'error' in r})


if __name__=='__main__':
    main()
