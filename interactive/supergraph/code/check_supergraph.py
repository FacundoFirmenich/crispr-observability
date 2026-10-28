"""Check standalone graph identities and exact rebuilding (not browser rendering).
SPDX-License-Identifier: GPL-3.0-only
"""
from pathlib import Path
import hashlib,json,re,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parents[1]

def sha(data):return hashlib.sha256(data).hexdigest()
def main():
    model_bytes=(ROOT/'supergraph_model.json').read_bytes();model=json.loads(model_bytes)
    html_bytes=(ROOT/'crispr_supergraph.html').read_bytes();html=html_bytes.decode()
    receipt=json.loads((ROOT/'BUILD_RECEIPT.json').read_text())
    match=re.search(r'<script type="application/json" id="supergraph-data">(.*?)</script>',html,re.S)
    assert match,'Missing embedded graph'
    raw=match.group(1)
    assert json.loads(raw)==model
    assert sha(model_bytes)==receipt['input_json_sha256']
    assert sha(raw.encode())==receipt['embedded_json_sha256']
    assert sha(html_bytes)==receipt['html_sha256']
    assert not re.search(r'<(?:script|link)[^>]+(?:src|href)=["\']https?://',html,re.I)
    ids={n['id'] for n in model['nodes']}
    assert len(ids)==len(model['nodes'])
    for e in model['edges']:assert e['source'] in ids and e['target'] in ids
    for gate in model['hyperedges']:assert gate['conclusion'] in ids and all(x in ids for x in gate['premises'])
    with tempfile.TemporaryDirectory(prefix='crispr-supergraph-') as f:
        out=Path(f)/'crispr_supergraph.html'
        subprocess.run([sys.executable,'-B',str(ROOT/'code/build_supergraph.py'),'--output',str(out)],check=True,stdout=subprocess.DEVNULL)
        assert out.read_bytes()==html_bytes
    print(json.dumps({'status':'PASS','model_object_equal':True,'byte_identical_rebuild':True,'nodes':len(ids),'edges':len(model['edges']),'AND_gates':len(model['hyperedges']),'html_sha256':sha(html_bytes),'scope':'Artifact integrity and rebuilding; browser interaction is checked separately.'},indent=2))
if __name__=='__main__':main()
