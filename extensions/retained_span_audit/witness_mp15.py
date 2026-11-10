import pathlib,json
from recompute_panel import solve
import optimal_lattice as lattice
R=pathlib.Path(__file__).resolve().parent;O=R/'ADMISSION_MP15B'
p=json.loads((O/'P1_POSITIVE_PROTOCOL.json').read_text());v=json.loads((O/'P1_POSITIVE_RESULT.json').read_text());out=[]
for x in v['results']:
 if x['status']!='COMPUTED':continue
 q=x['DNA_modal'][0]['sequence']
 for ref,expected in zip(p['references'],x['reference_deficits']):
  case={'reference_sequence':ref['sequence'],'observed_sequence':q,'lo':ref['lo'],'hi':ref['hi'],'target_kind':'HEK3_DNA','programmed_local':q[10:-10],'unedited_local':next(r['sequence'][10:-10] for r in p['references'] if r['name']=='JEB')}
  opt,w=solve(case,lattice,54);assert opt==expected['optimum']
  for k,z in w['programmed_local'].items():assert opt-z['score']==expected['event_deficits'][k]
  out.append({'clone':x['clone'],'event':x['event'],'reference':ref['name'],'optimum':opt,'witnesses':w})
(O/'P1_POSITIVE_WITNESSES.json').write_text(json.dumps({'reference_cases':len(out),'rescored_and_reprojected':True,'cases':out},indent=2)+'\n');print('WITNESS_CASES_PASS',len(out))
