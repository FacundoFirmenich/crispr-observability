"""Independent boundary cases for the finite-record calculator and source model."""
from fractions import Fraction
from pathlib import Path
import importlib.util, json, math, hashlib, re

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('build',ROOT/'code/build_model.py')
build=importlib.util.module_from_spec(spec);spec.loader.exec_module(build)
rejections=[]
for name,values in [('zero universe',(0,0,0,0)),('negative category',(10,2,9,-1)),('nonpartition',(10,2,4,3)),('noninteger',(10,2.5,4.5,3)),('bool-as-int',(10,True,6,3))]:
    try: build.partition(*values)
    except ValueError: rejections.append(name)
    else: raise AssertionError(name+' was accepted')
assert build.partition(10,2,5,3)==(Fraction(1,5),Fraction(1,2))
attainable={Fraction(2+k,10) for k in range(4)}
assert attainable=={Fraction(1,5),Fraction(3,10),Fraction(2,5),Fraction(1,2)}
assert Fraction(1,4) not in attainable # the interval hull does not imply continuous attainability
assert build.partition(10,2,8,0)==(Fraction(1,5),Fraction(1,5))
t={'N':10,'E':2,'NE':5,'U':3};c={'N':8,'E':1,'NE':5,'U':2}
delta=build.contrast(t,c)
assert Fraction(delta['L_exact'])==Fraction(-7,40)
assert Fraction(delta['H_exact'])==Fraction(3,8)
model=json.loads((ROOT/'evidence_model.json').read_text())
ids={n['id'] for n in model['graph']['nodes']}
assert len(ids)==len(model['graph']['nodes'])
for e in model['graph']['hyperedges']:
    assert e['type']=='AND' and len(set(e['premises']))>=2
    assert all(p in ids for p in e['premises']) and e['conclusion'] in ids
for section in ['observations','sensitivityObservations']:
    for r in model[section]:
        assert all(type(r[k]) is int and r[k]>=0 for k in ['N','E','NE','U'])
        assert r['N']==r['E']+r['NE']+r['U'] and r['N']>0
        assert Fraction(r['L_exact'])==Fraction(r['E'],r['N'])
        assert Fraction(r['H_exact'])==Fraction(r['E']+r['U'],r['N'])
        assert all(math.isfinite(r[k]) for k in ['L','H','L_pp','H_pp'])
assert len(model['sensitivityContrasts'])==192
assert {r['case']:r['inputUnit'] for r in model['externalCases']}['NHEJ']=='paired input unit'
assert {r['case']:r['unresolvedTemplateWeight'] for r in model['externalCases']}['P23H/WT']==572
text=(ROOT/'evidence_model.json').read_text()
assert not re.search(r'/(?:blob|commit|raw)/[0-9a-f]{7,40}(?:/|\b)',text)
assert '\\begin{document}' not in text
report={'status':'PASS','invalid_calculator_inputs_rejected':rejections,'finite_attainability_checked':True,'zero_unknowns_checked':True,'unequal_denominator_two_arm_extrema_checked':True,'AND_hyperedge_referential_integrity':True,'raw_and_assignment_units_kept_distinct':True,'old_manuscript_commit_links_absent':True,'model_sha256':hashlib.sha256((ROOT/'evidence_model.json').read_bytes()).hexdigest()}
(ROOT/'CONTRACT_BOUNDARY_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
