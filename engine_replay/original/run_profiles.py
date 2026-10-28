"""Pre-output stock-profile extension and prioritized execution of frozen OAT grid."""
import argparse, json, subprocess
from concurrent.futures import ThreadPoolExecutor
import run_grid as g

ORIGINAL_COMMAND=g.command
OMIT=['--default_min_aln_score','--needleman_wunsch_gap_open','--needleman_wunsch_gap_extend','--needleman_wunsch_gap_incentive','--quantification_window_coordinates','--exclude_bp_from_left','--exclude_bp_from_right']
def command(job,spec):
    if job['setting']!='stock_default_no_guide': return ORIGINAL_COMMAND(job,spec)
    cmd=ORIGINAL_COMMAND({**job,'setting':'baseline'},spec)
    cmd[cmd.index('-o')+1]=str(g.BASE/'runs'/job['reference']/job['setting'])
    for opt in OMIT:
        i=cmd.index(opt); del cmd[i:i+2]
    return cmd
g.command=command

def prepare():
    assert not (g.BASE/'STOCK_PROFILE_FREEZE.json').exists()
    args=json.loads((g.SOURCE/'CRISPResso2/args.json').read_text())
    def seek(obj,key):
        if isinstance(obj,dict):
            if key in obj:return obj[key]
            for v in obj.values():
                found=seek(v,key)
                if found is not None:return found
        return None
    defaults={k:seek(args,k) for k in [x[2:] for x in OMIT]}
    spec={'id':'STOCK_DEFAULT_NO_GUIDE','created_utc':g.now(),'state':'FIXED_BEFORE_PROFILE_EXECUTION','parent_grid_sha256':g.sha(g.BASE/'GRID.json'),
       'jobs':18,'scope':'same nine authenticated Q30 inputs and same two frozen references; single end, no guide, no trimming',
       'omitted_explicit_options':OMIT,'upstream_default_definitions':defaults,
       'retained_options':'explicit zero internal quality filters; one process; detailed alleles/annotated FASTQ; suppress reports/plots; output/reference/name/input',
       'interpretation':'Stock alignment/window profile without guide; several settings change jointly, separate from OAT. With no guide, effective quantification includes reference positions except upstream terminal exclusions; effective positions verified in every run info. MODIFIED is not the CFTR native peptide event.',
       'limits':'Retrospective sensitivity on the same biological genealogy. Custom baseline is not standard defaults.'}
    g.dump(g.BASE/'STOCK_PROFILE.json',spec)
    g.dump(g.BASE/'STOCK_PROFILE_FREEZE.json',{'created_utc':g.now(),'files':{'STOCK_PROFILE.json':g.sha(g.BASE/'STOCK_PROFILE.json'),'run_profiles.py':g.sha(__file__)},'state':'FROZEN'})
    print('STOCK_PROFILE_FROZEN',flush=True)

def run(which,workers):
    for f in ['GRID_FREEZE.json','STOCK_PROFILE_FREEZE.json']:
        for p,h in json.loads((g.BASE/f).read_text())['files'].items():assert g.sha(g.BASE/p)==h
    assert subprocess.check_output(['git','-C',str(g.SOURCE),'rev-parse','HEAD'],text=True).strip()==g.COMMIT
    subprocess.run(['git','-C',str(g.SOURCE),'diff','--exit-code'],check=True)
    spec=json.loads((g.BASE/'GRID.json').read_text()); alljobs=json.loads((g.BASE/'JOBS.json').read_text())
    if which=='baseline':jobs=[j for j in alljobs if j['setting']=='baseline']
    elif which=='oat':jobs=[j for j in alljobs if j['setting']!='baseline']
    else:jobs=[{**j,'setting':'stock_default_no_guide','id':j['id'].replace('/baseline/','/stock_default_no_guide/')} for j in alljobs if j['setting']=='baseline']
    with ThreadPoolExecutor(max_workers=workers) as pool:outcomes=list(pool.map(lambda j:g.execute(j,spec),jobs))
    g.dump(g.BASE/(which.upper()+'_STATUS.json'),{'state':'COMPLETED' if all(x['state']=='COMPLETED' for x in outcomes) else 'FAILED','finished_utc':g.now(),'jobs':len(outcomes),'workers':workers})
    assert all(x['state']=='COMPLETED' for x in outcomes)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','baseline','stock','oat']);p.add_argument('--workers',type=int,default=3);a=p.parse_args()
    if a.mode=='prepare':prepare()
    else:run(a.mode,a.workers)
