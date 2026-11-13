"""Recover the published CFTR sequence/function tables without inventing pairing.

Downloads two small public workbooks into memory; emits selected source cells,
full-object hashes and descriptive, explicitly unpaired function comparisons.
"""
from pathlib import Path
import argparse, hashlib, io, itertools, json, math, posixpath, urllib.request, zipfile
import xml.etree.ElementTree as ET

NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
BASE='https://static-content.springer.com/esm/art%3A10.1038%2Fs41551-024-01233-3/MediaObjects/'

def read_book(number):
    filename=f'41551_2024_1233_MOESM{number}_ESM.xlsx';url=BASE+filename
    with urllib.request.urlopen(url,timeout=30) as response:
        data=response.read(2_000_001)
    if len(data)>2_000_000:raise ValueError('Workbook exceeded bounded in-memory acquisition')
    z=zipfile.ZipFile(io.BytesIO(data))
    strings=[''.join(t.text or '' for t in e.findall('.//s:t',NS)) for e in ET.fromstring(z.read('xl/sharedStrings.xml'))]
    rels={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    sheets={}
    for s in ET.fromstring(z.read('xl/workbook.xml')).find('s:sheets',NS):
        target=rels[s.attrib[REL]];path=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
        cells={}
        for c in ET.fromstring(z.read(path)).findall('.//s:sheetData/s:row/s:c',NS):
            v=c.find('s:v',NS)
            if v is not None:cells[c.attrib['r']]=strings[int(v.text)] if c.attrib.get('t')=='s' else v.text
            elif c.attrib.get('t')=='inlineStr':cells[c.attrib['r']]=''.join(t.text or '' for t in c.findall('.//s:t',NS))
        sheets[s.attrib['name']]=cells
    return sheets,{'url':url,'filename':filename,'sha256':hashlib.sha256(data).hexdigest(),'MB':len(data)/1e6,'worksheets_parsed':list(sheets)}

def pearson(x,y):
    a=[v-sum(x)/len(x) for v in x];b=[v-sum(y)/len(y) for v in y]
    return sum(v*w for v,w in zip(a,b))/math.sqrt(sum(v*v for v in a)*sum(v*v for v in b))

def build():
    source,p4=read_book(4);meta,p3=read_book(3);cells=source['Figure 4'];st=meta['ST2 gRNA sequences + SRA files']
    assert cells['B4']=='Fig. 4a' and cells['B139']=='Fig. 4c'
    assert cells['C141']=='F&I' and cells['C142']=='CFTR-172'
    assert cells['J139']=='PE treated' and cells['D139']=='Mock'
    sequence=[float(cells[c+'8']) for c in 'DEF']
    endpoints={}
    for name,row,sign in [('F_and_I',141,1),('CFTR_inhibitor_response',142,-1)]:
        t=[sign*float(cells[c+str(row)]) for c in 'JKL'];c=[sign*float(cells[c+str(row)]) for c in 'DEF']
        permutations=[{'functional_replicate_order':list(order),'pearson':pearson(sequence,[t[i] for i in order])} for order in itertools.permutations(range(3))]
        endpoints[name]={'unit':'microampere per square centimetre','sign_convention':'published response' if sign==1 else 'negative of published inhibitor response',
            'treated':t,'mock':c,'unpaired_cross_difference_range':[min(t)-max(c),max(t)-min(c)],
            'all_observed_treated_exceed_all_observed_mock':min(t)>max(c),
            'conditional_column_order_pearson':pearson(sequence,t),
            'pairing_sensitivity_all_six':permutations,
            'pairing_sensitivity_range':[min(x['pearson'] for x in permutations),max(x['pearson'] for x in permutations)]}
    selected_meta={k:v for k,v in st.items() if k.lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ') in {'191','192','193'}}
    assert selected_meta['A191']==selected_meta['A192']==selected_meta['A193']=='F4a'
    selected_cells={k:v for k,v in cells.items() if k.lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ') in {'4','5','6','7','8','139','140','141','142'}}
    return {'source_doi':'10.1038/s41551-024-01233-3','source_attribution':'Sousa et al., original Figure 4 and Supplementary Table 2; CC BY 4.0 article. Derived analysis by the CRISPR Observability project.',
        'provenance':[p4,p3],'source_cells':selected_cells,'figure4a_metadata_cells':selected_meta,
        'published_CTT_insertion_percent':sequence,'endpoints':endpoints,
        'pairing_status':'UNAUTHENTICATED_ACROSS_ASSAYS','same_pool_status':'NOT_ESTABLISHED',
        'scope':'Three donor-level source observations per condition. Pairing sensitivity enumerates possible label assignments, not new biological observations; it is not a confidence interval or a validated sequence-function correlation.',
        'reason':'Both panels label Replicate 1-3, but the functional table lacks donor/culture identifiers authenticating correspondence to the sequencing file names. Column order alone is not a join key. The paper describes functional recording after three weeks of differentiation.',
        'required_link_fields':['donor_id','culture_or_pool_id','condition','sampling_time','sequencing_library_id','functional_assay_id','source_manifest_supporting_link'],
        'decision':'Admit condition-level descriptive functional context. Do not calibrate read-level certificates, pair records across assays, or infer sequence-function dose response from unauthenticated row order.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    if a.output.exists():raise ValueError('Preserve previous result; select a new output file')
    result=build();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'pairing_status':result['pairing_status'],'endpoints':{k:{f:v[f] for f in ['unpaired_cross_difference_range','pairing_sensitivity_range','conditional_column_order_pearson']} for k,v in result['endpoints'].items()}}))

if __name__=='__main__':main()
