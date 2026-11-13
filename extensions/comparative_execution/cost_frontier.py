"""Exact cost jurisdictions for donor-held-out routing; exploratory extension."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,argparse
from local_policy import ROOT,select,evaluate
MODES=['baseline','majority','pooled_empirical','local_empirical','local_beta11']
def cells():
 d=json.loads((ROOT/'extensions/generation_model/validation/results/HOMOPOLYMER_REPLICATION_COMPARISON.json').read_bytes())
 return [{'donor':donor,'tract':tract,**x} for donor,z in d['donors'].items() for tract,x in z['length'].items()]
def training(allcells,test,mode):
 return [c for c in allcells if c['donor']!=test['donor'] and (not mode.startswith('local') or c['tract']==test['tract'])]
def counts(allcells,r,mode):
 total=dict(C=0,D=0,U=0);actions=[]
 for test in allcells:
  use=mode=='majority' if mode in ['baseline','majority'] else select(training(allcells,test,mode),r,mode=='local_beta11')[0]
  for k,v in test['final' if use else 'initial'].items():total[k]+=v
  actions.append({'donor':test['donor'],'tract':test['tract'],'majority':use})
 assert sum(total.values())==241
 return {**total,'actions':actions}
def calculate():
 cs=cells();thresholds={F(0)}
 for mode in MODES[2:]:
  for test in cs:
   tr=training(cs,test,mode);d=sum(x['initial']['D']-x['final']['D'] for x in tr);n=sum(x['initial']['C']+x['initial']['D']-x['final']['C']-x['final']['D'] for x in tr)
   if n and (d or mode=='local_beta11'):thresholds.add(F(n+2,d+1) if mode=='local_beta11' else F(n,d))
 thresholds=sorted(thresholds);roots=set(thresholds)
 # Inside each routing cell loss is affine. Add every policy-loss crossing.
 for a,b in zip(thresholds,thresholds[1:]+[None]):
  r=(a+b)/2 if b else a+1;out={m:counts(cs,r,m) for m in MODES}
  for i,m in enumerate(MODES):
   for n in MODES[i+1:]:
    dd=out[m]['D']-out[n]['D'];du=out[m]['U']-out[n]['U']
    if dd:
     cross=F(-du,dd)
     if cross>a and (b is None or cross<b):roots.add(cross)
 roots=sorted(roots);intervals=[];points=[]
 def record(r):
  out={m:counts(cs,r,m) for m in MODES};loss={m:r*x['D']+x['U'] for m,x in out.items()};best=min(loss.values())
  return {'counts':out,'loss_at_evaluation':{m:str(v) for m,v in loss.items()},'best':[m for m,v in loss.items() if v==best], 'local_empirical_beats_both_fixed':loss['local_empirical']<min(loss['baseline'],loss['majority']), 'local_empirical_beats_pooled':loss['local_empirical']<loss['pooled_empirical']}
 for a,b in zip(roots,roots[1:]+[None]):
  r=(a+b)/2 if b else a+1;intervals.append({'lower_open':str(a),'upper_open':str(b) if b else None,'evaluation_ratio':str(r),**record(r)})
 for r in roots[1:]:points.append({'ratio':str(r),**record(r)})
 # Agreement with the separately implemented frozen-grid evaluation.
 for row in evaluate()['rows']:
  got=counts(cs,F(row['cost_ratio']),row['policy']);assert all(got[k]==row[k] for k in ['C','D','U'])
 # Probe additional exact rational costs on each cell and check no ranking change.
 for x in intervals:
  a=F(x['lower_open']);b=F(x['upper_open']) if x['upper_open'] else a+4
  for r in [(3*a+b)/4,(a+3*b)/4]:assert record(r)['best']==x['best']
 return {'status':'PASS','design':'Post-grid exploratory exact cost partition; no new biological observations','N':241,'routing_thresholds':[str(r) for r in thresholds[1:]],'intervals':intervals,'boundary_points':points,'claim_scope':'A complete finite-panel decision-cost map, not prospective dominance or calibrated probability.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();assert not a.output.exists();d=calculate();a.output.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps([{'lower':x['lower_open'],'upper':x['upper_open'],'best':x['best']} for x in d['intervals'] if x['local_empirical_beats_both_fixed']],indent=2));print('cells',len(d['intervals']),'boundaries',len(d['boundary_points']))
