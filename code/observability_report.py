"""Dependency-free user entry point for finite-record observation reports.

Counts are supplied observations, not inferred molecular truth. The existing
frozen event-contract implementation remains the provenance authority for CFTR.
"""
from pathlib import Path
from fractions import Fraction
import argparse, hashlib, html, json

ROOT = Path(__file__).resolve().parents[1]
VERSION = '1.0'

def require(condition, message):
    if not condition:
        raise ValueError(message)

def fraction(value):
    return {'exact':str(value), 'percent':float(100*value)}

def analyse(data):
    require(data.get('schema_version') == VERSION, 'Unsupported schema_version')
    contract = data['contract']
    for key in ['target','observation_unit','reference','coordinates','quality_policy','candidate_domain','accepted_label_assumption','decision_timing']:
        require(isinstance(contract.get(key),str) and bool(contract[key].strip()), 'Specify contract.'+key)
    require(data.get('comparability') == 'same_declared_contract', 'A contrast requires the same declared contract')
    groups=[]
    for name in ['treated','control']:
        arm=data[name]
        for key in ['sample_id','universe_id','source_locator']:
            require(isinstance(arm.get(key),str) and bool(arm[key].strip()), 'Specify '+name+'.'+key)
        n,e,ne,u=(arm[k] for k in ['N','E','NE','U'])
        require(all(type(v) is int and v>=0 for v in [n,e,ne,u]) and n>0 and e+ne+u==n, 'Counts must partition a positive fixed raw universe')
        lo,hi=Fraction(e,n),Fraction(e+u,n)
        groups.append({'arm':name, **arm,'bounds':{'lower':fraction(lo),'upper':fraction(hi)},'unresolved_fraction':fraction(Fraction(u,n)), 'resolved_only_event_fraction':fraction(Fraction(e,e+ne)) if e+ne else None})
    require(data['treated']['universe_id']!=data['control']['universe_id'], 'The same universe is not a treatment-control contrast')
    t,c=groups
    lo=Fraction(t['E'],t['N'])-Fraction(c['E']+c['U'],c['N'])
    hi=Fraction(t['E']+t['U'],t['N'])-Fraction(c['E'],c['N'])
    state='positive_for_all_binary_completions' if lo>0 else 'negative_for_all_binary_completions' if hi<0 else 'sign_not_resolved_by_these_records'
    return {'schema_version':VERSION,'contract':contract,'arms':groups,'contrast':{'lower':fraction(lo),'upper':fraction(hi),'decision':state},'interpretation':'Conditional finite-record contrast under the declared accepted-label assumption. Unknown records remain in the original universe. These bounds are not confidence intervals.','validation_scope':'Structural and arithmetic validation; input declarations and classification truth are not independently authenticated.'}

def demo_input():
    import importlib.util
    p=ROOT/'extensions/reproduction/contract/code/event_contract.py'
    spec=importlib.util.spec_from_file_location('frozen_contract',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    defs,cards=m.registry()
    t,c=cards['SRR29623272:witness:native'],cards['SRR29623273:witness:native']
    expected=m.contrast(t,c,defs,cards)
    event=defs[t['event_contract_sha256']]
    result={'schema_version':VERSION,'comparability':'same_declared_contract','contract':{
        'target':event['predicate'],'observation_unit':'sequenced raw record','reference':';'.join(event['references_sha256']),
        'coordinates':event['correspondence_level'],'quality_policy':event['quality'],'candidate_domain':'Frozen common observed-coding witness; complete-domain certificates are separate computations',
        'accepted_label_assumption':'Accepted E and NE labels are correct; unresolved labels may take either binary value',
        'decision_timing':'Retrospective; outcomes were known'},'provenance':{'frozen_contract_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'event_contract_sha256':t['event_contract_sha256']}}
    for key,card in [('treated',t),('control',c)]:
        result[key]={'sample_id':card['run_accession'],'universe_id':card['universe']['sequence_quality_stream_sha256'],'source_locator':card['source']['artifact'],'N':card['universe']['raw_N'],**card['counts']}
    actual=analyse(result)
    require(actual['contrast']['lower']['exact']==expected['L_exact'] and actual['contrast']['upper']['exact']==expected['H_exact'],'Independent calculation disagrees with frozen contract')
    return result

def utility():
    p=ROOT/'extensions/generation_model/validation/results/HOMOPOLYMER_REPLICATION_COMPARISON.json'
    data=json.loads(p.read_bytes());rows=[]
    for donor,v in data['donors'].items():
        for tract,x in v['length'].items():
            a,b=x['initial'],x['final'];n=sum(a.values())
            require(n==sum(b.values()),'Changing locus universe')
            require(a['C']-b['C']==x['withdrawn']['C'] and a['D']-b['D']==x['withdrawn']['D'],'Withdrawal accounting mismatch')
            rows.append({'donor':donor,'tract':tract,'N':n,'initial':a,'final':b,
                'discordances_removed':fraction(Fraction(a['D']-b['D'],a['D'])) if a['D'] else None,
                'concordances_sacrificed':fraction(Fraction(a['C']-b['C'],a['C'])) if a['C'] else None,
                'break_even_discordance_to_abstention_cost':str(Fraction(b['U']-a['U'],a['D']-b['D'])) if a['D']>b['D'] else None,
                'coverage_before':fraction(Fraction(a['C']+a['D'],n)), 'coverage_after':fraction(Fraction(b['C']+b['D'],n)),
                'discordance_among_retained_before':fraction(Fraction(a['D'],a['C']+a['D'])) if a['C']+a['D'] else None,
                'discordance_among_retained_after':fraction(Fraction(b['D'],b['C']+b['D'])) if b['C']+b['D'] else None})
    require(sum(x['N'] for x in rows)==data['total_distinct_loci'],'Overlapping or missing strata')
    totals={stage:{k:sum(r[stage][k] for r in rows) for k in ['C','D','U']} for stage in ['initial','final']}
    a,b=totals['initial'],totals['final']
    return {'source':p.relative_to(ROOT).as_posix(),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'scope':'Descriptive concordance against external germline genotype benchmarks in selected donor/tract jurisdictions; not editing accuracy or a threshold sweep.','local_rows':rows,'totals':totals,'decision_loss':{'definition':'L = c_D * discordant_calls + c_U * unresolved_calls; concordant_calls have zero loss','extra_unresolved':b['U']-a['U'],'discordances_removed':a['D']-b['D'],'break_even_c_D_over_c_U':str(Fraction(b['U']-a['U'],a['D']-b['D'])),'scope':'Decision sensitivity under user-specified costs; costs are not estimated or optimized from these outcomes.'},'interpretation':'Report both error withdrawal and information loss. A rule retaining no calls has undefined retained discordance, not zero error.'}

def render(result,local=None):
    h=html.escape;d=result['contrast']
    rows=''.join('<tr>'+''.join('<td>'+h(str(v))+'</td>' for v in [a['arm'],a['sample_id'],a['N'],a['E'],a['NE'],a['U'],f"{a['bounds']['lower']['percent']:.5f}–{a['bounds']['upper']['percent']:.5f}"])+'</tr>' for a in result['arms'])
    detail=''
    if local:
        detail='<h2>What does abstention buy—and cost?</h2><p>'+h(local['scope'])+'</p><table><tr><th>Donor</th><th>Tract</th><th>Initial C/D/U</th><th>Retained C/D/U</th><th>Discordances removed</th><th>Concordances sacrificed</th></tr>'
        for r in local['local_rows']:
            values=[r['donor'],r['tract'],'/'.join(str(r['initial'][k]) for k in ['C','D','U']),'/'.join(str(r['final'][k]) for k in ['C','D','U']),r['discordances_removed']['exact'] if r['discordances_removed'] else 'not estimable',r['concordances_sacrificed']['exact'] if r['concordances_sacrificed'] else 'not estimable']
            detail+='<tr>'+''.join('<td>'+h(str(v))+'</td>' for v in values)+'</tr>'
        detail+='</table><p>'+h(local['interpretation'])+'</p>'
    return '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>CRISPR Observability — decision report</title><style>body{font:17px/1.6 system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#283957;background:#fbfaff}h1,h2{color:#49388b}section{padding:22px;background:#eef2fc;border-radius:14px;margin:20px 0}table{border-collapse:collapse;width:100%}td,th{padding:10px;border-bottom:1px solid #cbd3e6;text-align:left}pre{white-space:pre-wrap}summary{cursor:pointer;color:#7055a4}</style><h1>What can these reads establish?</h1><section><h2>1 · Name the question</h2><p>'+h(result['contract']['target'])+'</p><p>'+h(result['contract']['quality_policy'])+'</p></section><section><h2>2 · Keep every acquired record</h2><table><tr><th>Arm</th><th>Sample</th><th>N</th><th>E</th><th>NE</th><th>U</th><th>Bounds (%)</th></tr>'+rows+'</table></section><section><h2>3 · Read the decision</h2><p><strong>'+h(d['decision'].replace('_',' '))+'</strong></p><p>Treated − control: '+f"{d['lower']['percent']:+.5f} to {d['upper']['percent']:+.5f}"+' percentage points.</p><p>'+h(result['interpretation'])+'</p></section><section><h2>4 · Choose the next measurement</h2><p>If the sign is unresolved, inspect which records are U and which flanks are absent. If it is resolved, test sensitivity to accepted-label error and connect functional measurements only at an authenticated shared unit. Use the 3D laboratory to test geometry, and the Supergraph to trace premises.</p></section>'+detail+'<details><summary>Inspect the complete machine-readable report</summary><pre>'+h(json.dumps(result,indent=2))+'</pre></details></html>'

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['demo','analyse']);p.add_argument('input',nargs='?',type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    data=demo_input() if a.action=='demo' else json.loads(a.input.read_bytes())
    report=analyse(data);report['input_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()
    a.output.mkdir(parents=True,exist_ok=True)
    for name in ['input.json','report.json','report.html','local_utility.json']:
        require(not (a.output/name).exists(),'Choose a new output directory; prior results are preserved')
    local=utility() if a.action=='demo' else None
    (a.output/'input.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    (a.output/'report.html').write_text(render(report,local),encoding='utf-8')
    if local:(a.output/'local_utility.json').write_text(json.dumps(local,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','output':str(a.output),'decision':report['contrast']['decision']}))

if __name__=='__main__':main()
