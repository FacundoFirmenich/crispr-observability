"""Post hoc constructive extension from the observed binary-label obstruction.
No new observations. Three nonempty actual target signatures establish the
minimum number of classes needed to preserve both events on the observed set.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R/'src'))
from observed_contracts.core import compile_partition, property_value

def main():
    data=json.loads((R/'results/CODING_CATALOGS.json').read_text());exact=data['exact_programmed_coding']
    counts=json.loads((R/'inputs/M8_WITNESS_COUNTS.json').read_text())
    idx={(r['run_accession'],'native' if r['event'].startswith('native') else 'exact'):r for r in counts}
    rows=[]
    for run,cat in data['catalogs']['witness'].items():
        p=compile_partition(cat,['native_event','exact_programmed_event'],exact)
        assert all(not g['target_signature'][1] or g['target_signature'][0] for g in p['classes'])
        assert all(not property_value(r['coding'],'exact_programmed_event',exact) or property_value(r['coding'],'native_event',exact) for r in cat)
        row={'run':run,'partition':p}
        if run!='ALL':
            n,e=idx[run,'native'],idx[run,'exact'];assert n['N']==e['N'] and n['U']==e['U']
            b=n['E']-e['E'];assert b>=0
            assert sum(g['records'] for g in p['classes'] if g['target_signature']==[True,False])==b
            row['raw_reporting_counts']={'N':n['N'],'exact_programmed':e['E'],'native_other_observed_coding':b,'non_native':n['NE'],'unknown':n['U']}
            row['native_minus_exact_within_run']={'naive_separate_marginal_hull_exact':[str(F(b-n['U'],n['N'])),str(F(b+n['U'],n['N']))],'joint_logically_compatible_hull_exact':[str(F(b,n['N'])),str(F(b+n['U'],n['N']))],'premise':'Exact programmed coding implies native local peptide; both predicates refer to the same witnessed coding or shared unresolved record. Not a functional benefit.'}
        rows.append(row)
    pooled=next(r for r in rows if r['run']=='ALL')
    assert pooled['partition']['minimum_observed_classes']==3
    (R/'results/JOINT_TARGET_PARTITIONS.json').write_text(json.dumps({'timing':'Post hoc constructive extension after the finite-catalog audit on already known data','minimum_classes_not_new_theorem':True,'rows':rows},ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(pooled))

if __name__=='__main__':main()
