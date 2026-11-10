"""Validate release identity and claim/evidence consistency without a network."""
from pathlib import Path
import argparse, collections, hashlib, json

ROOT=Path(__file__).resolve().parent
def validate(root=ROOT):
    root=Path(root).resolve();manifest=json.loads((root/'MANIFEST.json').read_text())
    paths=set()
    for row in manifest['files']:
        p=(root/row['path']).resolve()
        assert p.is_relative_to(root) and row['path'] not in paths,row['path']
        assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
        paths.add(row['path'])
    graph=json.loads((root/'epistemic_graph.json').read_text(encoding='utf8'))
    byid={n['id']:n for n in graph['nodes']};assert len(byid)==len(graph['nodes'])
    allowed_relations={'requires','checks','limits','locates','repaired_before'}
    successors=collections.defaultdict(list);indegree=collections.Counter({k:0 for k in byid})
    for e in graph['edges']:
        assert e['source'] in byid and e['target'] in byid and e['relation'] in allowed_relations
        successors[e['source']].append(e['target']);indegree[e['target']]+=1
    queue=[k for k,n in indegree.items() if n==0];seen=[]
    while queue:
        k=queue.pop();seen.append(k)
        for t in successors[k]:
            indegree[t]-=1
            if indegree[t]==0:queue.append(t)
    assert len(seen)==len(byid),'Dependency cycle'
    for n in byid.values():
        assert n['scope'] and n['status']
        for a in n['artifacts']: assert a in paths,(n['id'],a)
        if n['kind']=='claim':
            assert n['falsifier'] and n['assumptions'] and n['artifacts'],n['id']
            for a in n['assumptions']:
                assert byid[a]['kind']=='premise'
                assert dict(source=a,target=n['id'],relation='requires') in graph['edges']
    html=(root/'epistemic_graph.html').read_text(encoding='utf8')
    encoded=html.split('<script id="graph-data" type="application/json">',1)[1].split('</script>',1)[0]
    assert json.loads(encoded)==graph,'HTML graph differs from machine-readable source'
    for name in graph['receipt_sources']: assert json.loads((root/name).read_text())['status']=='PASS'
    disagreements=json.loads((root/'CODING_DISAGREEMENT_RECORDS.json').read_text())
    assert len(disagreements)==83
    assert all(len(set(r['codings'].values()))>1 for r in disagreements)
    assert sum(len(set(r['native_labels']))==1 for r in disagreements)==18
    return dict(status='PASS',files=len(paths),nodes=len(byid),edges=len(graph['edges']),
                claims=sum(n['kind']=='claim' for n in byid.values()),
                status_counts=dict(collections.Counter(n['status'] for n in byid.values())),
                browser_interaction='NOT_EXECUTED',fresh_environment_install='NOT_EXECUTED')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);a=ap.parse_args();r=validate()
    if a.output:
        assert not a.output.exists();a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
