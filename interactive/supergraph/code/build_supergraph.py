"""Build the standalone CRISPR adaptation of the published KwanTube spatial map.

SPDX-License-Identifier: GPL-3.0-only
The interface derives from Facundo Firmenich's KwanTube v3.5.1.1 map.
Scientific source data retain their own terms; see SOURCE_NOTICE.md.
"""
from pathlib import Path
import argparse, hashlib, json, re
ROOT = Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--model',type=Path,default=ROOT/'supergraph_model.json')
    p.add_argument('--template',type=Path,default=ROOT/'crispr_supergraph.template.html')
    p.add_argument('--output',type=Path,default=ROOT/'crispr_supergraph.html')
    a=p.parse_args()
    source=a.model.read_bytes();model=json.loads(source)
    raw=json.dumps(model,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    embedded=hashlib.sha256(raw.encode()).hexdigest()
    template=a.template.read_text()
    assert template.count('__SUPERGRAPH_DATA__')==1
    html=template.replace('__SUPERGRAPH_DATA__',raw).replace('__EMBEDDED_SHA256__',embedded).replace('__INPUT_SHA256__',hashlib.sha256(source).hexdigest())
    assert not re.search(r'<(?:script|link)[^>]+(?:src|href)=[\"\']https?://',html,re.I)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(html)
    receipt={'html':a.output.name,'html_sha256':hashlib.sha256(html.encode()).hexdigest(),
       'input_json_sha256':hashlib.sha256(source).hexdigest(),'embedded_json_sha256':embedded,
       'nodes':len(model['nodes']),'edges':len(model['edges']),'AND_gates':len(model['hyperedges']),
       'reference_architecture':'Published KwanTube v3.5.1.1 immersive spatial SVG map; adapted, not an identical scientific dataset.',
       'license':'GPL-3.0-only for derived interface/build code; data terms remain separate.',
       'dependencies':'None: data, CSS and JavaScript embedded; no CDN.',
       'scope':'Declared evidence navigation; no new biological experiment or molecular accuracy certification.'}
    (a.output.parent/'BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
