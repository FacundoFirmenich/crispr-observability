"""Prevent task invariance, detected editing and cross-model coverage being conflated."""
from pathlib import Path
from collections import Counter,defaultdict
import csv,json
HERE=Path(__file__).resolve().parent
def build():
 out={}
 for name,targets in [('cftr_lattice',['native','exact']),('hek3_lattice',['programmed_local','unedited_local'])]:
  cases=json.loads((HERE/'expected'/name/'OPTIMAL_ALIGNMENT_ENVELOPES.json').read_text())['cases'];out[name]={}
  for target in targets:
   cases_added=[c for c in cases if c[target+'_invariant'] and not c['all_optimal_strings_identical']]
   out[name][target]={'additional_invariant_case_labels':{'event':sum(c[target+'_values']==[True] for c in cases_added),'non_event':sum(c[target+'_values']==[False] for c in cases_added)},'scope':'Read-reference cases in the selected actual panel; not accuracy or additional editing discoveries.'}
 cases=json.loads((HERE/'expected/cftr_lattice/OPTIMAL_ALIGNMENT_ENVELOPES.json').read_text())['cases'];domains=defaultdict(set)
 for c in cases:domains[c['read_id']].update(c['valid_coding_strings'])
 focus=list(csv.DictReader((HERE/'inputs/FOCUS_89.csv').open()));mechanism=list(csv.DictReader((HERE/'expected/mechanism/FOCUS_89_MECHANISM.csv').open()));assert {x['read_id'] for x in focus}=={x['read_id'] for x in mechanism}
 compatible=sum(x['PE_default'] in domains[x['read_id']] for x in mechanism)
 out['selected_PE_vs_zero_incentive_two_reference_union']={'focus_records':len(mechanism),'PE_coding_in_union':compatible,'PE_coding_outside_union':len(mechanism)-compatible,'scope':'Compatibility with the union of the two declared fixed-reference zero-incentive optimal domains, not the selected PE reference alone. Outside-domain codings demonstrate that this optimal-path certificate cannot silently cover dedicated PE mode or other objectives/reference grammars.'}
 rows=list(csv.DictReader((HERE/'expected/mechanism/FAMILY_CERTIFICATE_COUNTS.csv').open()));out['F9_additional_task_determinations']={}
 for event in ['native','exact_programmed','native_nonexact']:
  def total(cert,key):return sum(int(r[key]) for r in rows if r['family']=='F9' and r['certificate']==cert and r['event']==event)
  out['F9_additional_task_determinations'][event]={'event_delta':total('target_relative','E')-total('string_witness','E'),'non_event_delta':total('target_relative','NE')-total('string_witness','NE')}
 assert compatible==44 and len(mechanism)==89
 assert all(v['event_delta']==0 for v in out['F9_additional_task_determinations'].values())
 out['scope']='Actual conditional negative determinations are useful but do not constitute additional positive edit calls or biological accuracy gains. Fixed-model completeness is distinct from cross-model coverage.'
 return out
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.write_text(json.dumps(build(),indent=2)+'\n');print(a.output.read_text())
