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
    p.add_argument('--model',type=Path,default=ROOT/'review_graph_model.json')
    p.add_argument('--tour',type=Path,default=ROOT/'guided_tour.json')
    p.add_argument('--template',type=Path,default=ROOT/'crispr_supergraph.template.html')
    p.add_argument('--output',type=Path,default=ROOT/'crispr_supergraph.html')
    a=p.parse_args()
    source=a.model.read_bytes();model=json.loads(source)
    raw=json.dumps(model,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    embedded=hashlib.sha256(raw.encode()).hexdigest()
    tour_source=a.tour.read_bytes()
    tour=json.loads(tour_source)
    graph_nodes={n['id'] for n in model['nodes']}
    graph_edges={e['id']:(e['source'],e['target']) for e in model['edges']}
    assert tour and len({block['title'] for block in tour})==len(tour)
    for block in tour:
        assert set(block['nodes']) <= graph_nodes and block['focus'] in block['nodes']
        assert block['edges'] and set(block['nodes'])=={node for edge in block['edges'] for node in graph_edges[edge]}
        unresolved=set(block['nodes'])
        while unresolved:
            ready=unresolved-{graph_edges[edge][1] for edge in block['edges']
                              if graph_edges[edge][0] in unresolved}
            assert ready, f"Cycle in tour block: {block['title']}"
            unresolved-=ready
        assert all((ROOT.parents[1]/path).is_file() for path in block.get('sources', []))
        assert all(block[field].strip() for field in ('title','question','finding','scope'))
    tour_raw=json.dumps(tour,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    template=a.template.read_text(encoding='utf-8')
    assert template.count('__SUPERGRAPH_DATA__')==1 and template.count('__GUIDED_TOUR_DATA__')==1
    html=template.replace('__SUPERGRAPH_DATA__',raw).replace('__EMBEDDED_SHA256__',embedded).replace('__INPUT_SHA256__',hashlib.sha256(source).hexdigest())
    html=html.replace('__GUIDED_TOUR_DATA__',tour_raw)
    assert not re.search(r'<(?:script|link)[^>]+(?:src|href)=[\"\']https?://',html,re.I)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(html, encoding='utf-8', newline='\n')
    receipt={'html':a.output.name,'html_sha256':hashlib.sha256(html.encode()).hexdigest(),
       'input_json_sha256':hashlib.sha256(source).hexdigest(),'embedded_json_sha256':embedded,
       'guided_tour_sha256':hashlib.sha256(tour_source).hexdigest(),'guided_tour_blocks':len(tour),
       'nodes':len(model['nodes']),'edges':len(model['edges']),'AND_gates':len(model['hyperedges']),
       'reference_architecture':'Published KwanTube v3.5.1.1 immersive spatial SVG map; adapted, not an identical scientific dataset.',
       'license':'GPL-3.0-only for derived interface/build code; data terms remain separate.',
       'dependencies':'None: data, CSS and JavaScript embedded; no CDN.',
       'scope':'Declared evidence navigation; no new biological experiment or molecular accuracy certification.'}
    (a.output.parent/'BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
