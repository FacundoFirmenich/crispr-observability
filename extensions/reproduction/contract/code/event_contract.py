"""Executable observed-event contract on the actual CFTR evidence.

Validates source identity, event/sequence contract, denominator and claim scope.
Metadata mutation checks are software checks; they add no biological records.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,copy,hashlib,json
B=Path(__file__).resolve().parents[1];S=B/'inputs/M7/inputs'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def read(s):return json.loads((S/s).read_text())
def put(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
class ContractError(ValueError):pass
def need(ok,message):
    if not ok:raise ContractError(message)
def registry():
    runs={r['run_accession']:r for r in read('M2/DESIGN.json')['runs']}
    refhashes=[sha(S/'M2/REFERENCE.json'),sha(S/'M2/REFERENCE_PROGRAMMED.json')]
    defs={};cards=[]
    for source,operator,defaultfamily in [('M3/WITNESS_COUNTS.json','witness',None),('M2/COUNTS.json','mutant_projection','native'),('M2/REFERENCE_COUNTS.json','programmed_projection','native')]:
        for row in read(source):
            if operator=='mutant_projection' and (row['observer']!='crispresso_projected' or not row['event'].startswith('native')):continue
            if operator=='programmed_projection' and (row['observer']!='programmed_reference' or not row['event'].startswith('native')):continue
            family=defaultfamily or ('native' if row['event'].startswith('native') else 'CTT_full_SE2')
            event={'family':family,'predicate':('in frame observed coding translates to GTIKENIIFGVSYD' if family=='native' else 'observed coding equals full programmed core plus AT'),'sequence_domain':operator,'correspondence_level':'identical observed coding string' if operator=='witness' else 'observed bases projected on declared reference interval','references_sha256':refhashes if operator=='witness' else [refhashes[0 if operator=='mutant_projection' else 1]],'quality':'arithmetic_mean_original_whole_read_Phred33_ge_30','raw_transform':'none','warning_rule':'unresolved required motif or flank conditions remain U','biological_truth_calibrated':False}
            eid=digest(event);defs[eid]=event;r=runs[row['run_accession']]
            card={'schema_version':'1.0','id':f"{r['run_accession']}:{operator}:{family}",'run_accession':r['run_accession'],'preparation_nominal':r['source_replica'],'condition':r['condition'],'event_contract_sha256':eid,'universe':{'raw_N':row['N'],'sequence_quality_stream_sha256':r['sequence_quality_stream_sha256'],'unit':'sequenced raw record'},'counts':{k:row[k] for k in ['E','NE','U']},'interval_hull_exact':[row['L_exact'],row['H_exact']],'source':{'bioproject':'PRJNA1055086','doi':'10.1038/s41551-024-01233-3','artifact':source,'sha256':sha(S/source)},'information_timing':'retrospective outcomes known','permitted_claim':'finite record binary completion under correct accepted labels','not_authenticated':['classification truth','independent cohort','unique physical copy','full allele','cell identity','per haplotype function','population efficacy']}
            cards.append(card)
    return defs,{r['id']:r for r in cards}
def validate(card,defs,expected):
    need(card['id'] in expected,'unknown original card identity');ref=expected[card['id']]
    need(set(card)==set(ref),'unexpected or missing contract fields')
    need(card['schema_version']=='1.0','schema version')
    need(card['event_contract_sha256'] in defs,'unknown event contract')
    need(digest(defs[card['event_contract_sha256']])==card['event_contract_sha256'],'event digest mismatch')
    need(card['source']['sha256']==sha(S/card['source']['artifact']),'source artifact digest mismatch')
    need(card['source']==ref['source'],'provenance declaration differs from the original evidence')
    need(card['universe']==ref['universe'],'raw universe or denominator mismatch')
    need(card['event_contract_sha256']==ref['event_contract_sha256'],'event or correspondence contract changed')
    need(card['counts']==ref['counts'],'counts differ from frozen evidence')
    need(card['run_accession']==ref['run_accession'] and card['condition']==ref['condition'] and card['preparation_nominal']==ref['preparation_nominal'],'source assignment mismatch')
    need(card['information_timing']=='retrospective outcomes known','no retrospective to prospective relabeling')
    need(card['permitted_claim']==ref['permitted_claim'] and card['not_authenticated']==ref['not_authenticated'],'unsupported scientific claim scope')
    n=card['universe']['raw_N'];e,ne,u=(card['counts'][k] for k in ['E','NE','U'])
    need(all(type(x)==int and x>=0 for x in [e,ne,u,n]) and n>0 and e+ne+u==n,'invalid count partition')
    need([F(e,n),F(e+u,n)]==[F(x) for x in card['interval_hull_exact']],'invalid sharp interval hull')
    return True
def contrast(t,c,defs,expected,claim='finite_record_contrast'):
    validate(t,defs,expected);validate(c,defs,expected)
    need(claim=='finite_record_contrast','requested claim requires evidence not authenticated here')
    need(t['event_contract_sha256']==c['event_contract_sha256'],'incompatible event or observed segment contracts')
    need(t['preparation_nominal']==c['preparation_nominal'],'nominal preparation mismatch')
    need(t['condition']=='PE treated' and c['condition'] in ['mock','untreated'],'unsupported source contrast')
    need(t['run_accession']!=c['run_accession'],'same run is not a treatment contrast')
    tl,th=map(F,t['interval_hull_exact']);cl,ch=map(F,c['interval_hull_exact']);lo,hi=tl-ch,th-cl
    return {'treated':t['id'],'control':c['id'],'event_contract_sha256':t['event_contract_sha256'],'L_exact':str(lo),'H_exact':str(hi),'L_pp':float(100*lo),'H_pp':float(100*hi),'strict_positive':lo>0,'scope':'Conditional finite raw record contrast; no claim of classification accuracy or biological efficacy.'}
def build():
    defs,cards=registry()
    for card in cards.values():validate(card,defs,cards)
    rows=[]
    for t in cards.values():
        if t['condition']!='PE treated':continue
        for c in cards.values():
            if c['event_contract_sha256']==t['event_contract_sha256'] and c['condition'] in ['mock','untreated'] and c['preparation_nominal']==t['preparation_nominal']:rows.append(contrast(t,c,defs,cards))
    negative=[]
    actual=cards['SRR29623272:witness:native'];control=cards['SRR29623273:witness:native']
    for label,field,value in [('wrong raw denominator','universe',dict(actual['universe'],raw_N=actual['universe']['raw_N']+1)),('retrospective presented as prospective','information_timing','prospective'),('claim elevated to efficacy','permitted_claim','biological efficacy')]:
        mutated=copy.deepcopy(actual);mutated[field]=value
        try:validate(mutated,defs,cards)
        except ContractError as exc:negative.append({'check':label,'status':'rejected','reason':str(exc)})
        else:raise AssertionError('invalid metadata accepted')
    for label,c,claim in [('cross observer segment join',cards['SRR29623273:mutant_projection:native'],'finite_record_contrast'),('unsupported accuracy claim',control,'biological_accuracy')]:
        try:contrast(actual,c,defs,cards,claim)
        except ContractError as exc:negative.append({'check':label,'status':'rejected','reason':str(exc)})
        else:raise AssertionError('unsupported comparison accepted')
    spec={'version':'1.0','status':'proposed reporting and verification contract, exercised on this case only','structural_fields':list(actual),'semantics':['fixed raw universe and original record identity','declared event and observed sequence correspondence level','quality and transform policy','complete E NE U partition and exact interval hull','source provenance and retrospective information timing','separate unverified claims','within contract comparison only'],'not_a_claim':['community standard adoption','independent ground truth','experimental validation','universal accuracy','clinical safety']}
    put(B/'delivery/OBSERVED_EVENT_CONTRACT_1_0.json',{'specification':spec,'event_definitions':defs,'actual_cards':list(cards.values()),'actual_contrasts':rows})
    put(B/'results/CONTRACT_VERIFICATION.json',{'status':'PASS','real_cards':len(cards),'within_contract_contrasts':len(rows),'semantic_rejections':negative,'no_biological_data_created':True})
    print(json.dumps({'status':'PASS','real_cards':len(cards),'contrasts':len(rows),'negative_checks':len(negative)}))
def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['build','contrast']);p.add_argument('--treated');p.add_argument('--control');a=p.parse_args()
    if a.action=='build':build()
    else:
        defs,cards=registry();need(a.treated in cards and a.control in cards,'actual card not found');print(json.dumps(contrast(cards[a.treated],cards[a.control],defs,cards),indent=2))
if __name__=='__main__':main()
