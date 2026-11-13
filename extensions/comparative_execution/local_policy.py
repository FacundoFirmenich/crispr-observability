"""Donor-held-out decision routing on authentic frozen germline cells."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[2]
COSTS=(2,4,6,10,20)
def select(training,cost,prior):
 d=sum(x['initial']['D']-x['final']['D'] for x in training)
 c=sum(x['initial']['C']-x['final']['C'] for x in training)
 n=d+c
 if not n:return False,{'withdrawn':0,'discordant':0,'posterior_mean':None}
 p=F(d+(1 if prior else 0),n+(2 if prior else 0))
 return cost*p>1,{'withdrawn':n,'discordant':d,'estimated_discordant_fraction':str(p),'expected_unit_loss_change':str(1-cost*p)}
def evaluate():
 p=ROOT/'extensions/generation_model/validation/results/HOMOPOLYMER_REPLICATION_COMPARISON.json';raw=p.read_bytes();data=json.loads(raw);cells=[]
 for donor,info in data['donors'].items():
  for tract,x in info['length'].items():
   assert sum(x['initial'].values())==sum(x['final'].values())
   assert x['final']['D']<=x['initial']['D'] and x['final']['C']<=x['initial']['C']
   cells.append({'donor':donor,'tract':tract,**x})
 assert len(cells)==8 and sum(sum(c['initial'].values()) for c in cells)==241
 rows=[];details=[]
 for cost in COSTS:
  for mode in ['baseline','majority','pooled_empirical','local_empirical','local_beta11']:
   total={k:0 for k in ['C','D','U']}
   for test in cells:
    training=[c for c in cells if c['donor']!=test['donor'] and (not mode.startswith('local') or c['tract']==test['tract'])]
    assert test['donor'] not in {c['donor'] for c in training}
    use,stats=select(training,F(cost),mode=='local_beta11') if mode not in ['baseline','majority'] else (mode=='majority',{})
    chosen=test['final' if use else 'initial'];loss=cost*chosen['D']+chosen['U']
    details.append({'cost_ratio':cost,'policy':mode,'held_out_donor':test['donor'],'tract':test['tract'],'training_donors':sorted({x['donor'] for x in training}),'training':stats,'action':'majority' if use else 'retain','counts':chosen,'loss':loss})
    for k in total:total[k]+=chosen[k]
   assert sum(total.values())==241
   rows.append({'cost_ratio':cost,'policy':mode,**total,'loss':cost*total['D']+total['U'],'coverage':str(F(total['C']+total['D'],241)),'retained_discordance':str(F(total['D'],total['C']+total['D'])) if total['C']+total['D'] else None})
 return {'status':'COMPLETE','source_sha256':hashlib.sha256(raw).hexdigest(),'plan_sha256':hashlib.sha256((Path(__file__).parent/'ANALYSIS_PLAN.md').read_bytes()).hexdigest(),'timing':'Exploratory donor-held-out reanalysis of previously observed panels; policies fixed before this evaluation','N':241,'donors':sorted(data['donors']),'cost_grid':COSTS,'rows':rows,'local_decisions':details,'scope':'External germline benchmark concordance; no original-fragment editing claim. The folds reuse training donors and may share loci across donors.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists();d=evaluate();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d['rows'],indent=2))
