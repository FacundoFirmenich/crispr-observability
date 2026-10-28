"""Verify all final run hashes and record a closed recoverable file inventory."""
import gzip, json, os, zipfile
import run_grid as g

assert json.loads((g.BASE/'VERIFICATION.json').read_text())['stage']=='all'
records=list((g.BASE/'executions').rglob('*.json'));assert len(records)==144
checks=[]
for p in records:
    r=json.loads(p.read_text());assert r['state']=='COMPLETED'
    folder=g.BASE/'runs'/r['reference']/r['setting']/('CRISPResso_on_'+r['run'])
    for name,meta in r.get('files',{}).items():assert g.sha(folder/name)==meta['sha256'],str(folder/name)
    for z in folder.glob('*.zip'):
        with zipfile.ZipFile(z) as f:assert f.testzip() is None
    for z in folder.glob('*.gz'):gzip.open(z,'rb').read()
    checks.append({'id':r['id'],'state':'PASS'})
g.dump(g.BASE/'FINAL_OUTPUT_INTEGRITY.json',{'state':'PASS','utc':g.now(),'jobs':144,'checks':checks})
inventory={}
for p in sorted(g.BASE.rglob('*')):
    if not p.is_file() or p.name=='CLOSED_FILE_MANIFEST.json' or '__pycache__' in p.parts:continue
    assert not p.name.endswith('.tmp'), 'Incomplete temporary output: '+str(p)
    with p.open('rb') as f:os.fsync(f.fileno())
    inventory[str(p.relative_to(g.BASE))]={'bytes':p.stat().st_size,'sha256':g.sha(p)}
g.dump(g.BASE/'CLOSED_FILE_MANIFEST.json',{'state':'CLOSED_VERIFIED_CHECKPOINT','utc':g.now(),'engine_bundled':False,'files':inventory,'total_bytes':sum(x['bytes'] for x in inventory.values())})
for name,meta in inventory.items():assert g.sha(g.BASE/name)==meta['sha256']
print(json.dumps({'state':'CLOSED_VERIFIED_CHECKPOINT','files':len(inventory),'total_bytes':sum(x['bytes'] for x in inventory.values())}),flush=True)
