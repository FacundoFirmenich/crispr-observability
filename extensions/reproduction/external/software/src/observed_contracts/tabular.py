"""Finite weighted-row contracts; no biological truth is inferred from a row.

Unlike the CFTR-specific properties in core.py, features are supplied explicitly.
None is unresolved for a target, never a new biological class. A certificate on
resolved rows cannot authorize the complete table when any target is unresolved.
"""
from collections import defaultdict
import json

SCOPE='finite_external_assignment_catalog'

def key(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False)

def validate(rows,names):
    if not rows:raise ValueError('A nonempty observed catalog is required')
    ids=set()
    for row in rows:
        if type(row.get('weight')) is not int or row['weight']<=0:
            raise ValueError('Positive integer assignment weight required')
        rid=key(row['locator'])
        if rid in ids:raise ValueError('Duplicate row locator')
        ids.add(rid)
        if any(n not in row['features'] for n in names):
            raise ValueError('Missing feature definition; use None for unresolved')

def audit(rows,representation,target):
    validate(rows,[representation,target])
    fibers=defaultdict(dict);unknown=0;target_support=set()
    for row in rows:
        value=row['features'][target]
        if value is None:
            unknown+=row['weight'];continue
        rk=key(row['features'][representation]);tk=key(value);target_support.add(tk)
        cell=fibers[rk].setdefault(tk,{'weight':0,'locator':row['locator']})
        cell['weight']+=row['weight']
    bad=[(r,ys) for r,ys in sorted(fibers.items()) if len(ys)>1]
    witnesses=[]
    for r,ys in bad[:8]:
        ks=sorted(ys)[:2]
        witnesses.append({'representation_value':json.loads(r),
          'target_values':[json.loads(k) for k in ks],
          'actual_rows':[ys[k]['locator'] for k in ks]})
    return {'representation':representation,'target':target,'scope':SCOPE,
      'assignment_weight':sum(r['weight'] for r in rows),
      'target_unresolved_weight':unknown,'observed_resolved_fibers':len(fibers),
      'distinct_resolved_target_values':len(target_support),
      'target_constant_on_resolved_domain':len(target_support)==1,
      'heterogeneous_fibers':len(bad),
      'assignment_weight_in_heterogeneous_fibers':sum(sum(v['weight'] for v in ys.values()) for _,ys in bad),
      'factors_on_resolved_subcatalog':not bad,
      'complete_catalog_authorized':not bad and unknown==0,
      'counterexamples':witnesses}

def require_claim(certificate,scope):
    if scope!=SCOPE:
        raise ValueError('Assignment-row contract cannot authenticate raw denominators, molecules, biological accuracy or future sufficiency')
    if not certificate['complete_catalog_authorized']:
        raise ValueError('Target is heterogeneous or unresolved on the supplied complete catalog')
    return True

def compile_joint(rows,targets):
    if not targets or len(set(targets))!=len(targets):raise ValueError('Distinct nonempty targets required')
    validate(rows,targets);cells={};unknown=0
    for row in rows:
        values=[row['features'][t] for t in targets]
        if any(v is None for v in values):unknown+=row['weight'];continue
        k=key(values);cell=cells.setdefault(k,{'signature':values,'weight':0,'rows':0,'actual_row':row['locator']})
        cell['weight']+=row['weight'];cell['rows']+=1
    return {'scope':SCOPE,'targets':targets,'minimum_classes_on_joint_resolved_domain':len(cells),
      'joint_unresolved_weight':unknown,'joint_resolved_weight':sum(c['weight'] for c in cells.values()),
      'classes':[cells[k] for k in sorted(cells)],
      'authority':'Coarsest observed partition for the specified targets and resolved domain; no out-of-support or truth claim.'}

def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('catalog');p.add_argument('--representation',required=True);p.add_argument('--target',required=True)
    p.add_argument('--require-scope')
    a=p.parse_args();rows=json.load(open(a.catalog))['rows'];c=audit(rows,a.representation,a.target)
    if a.require_scope:
        try:require_claim(c,a.require_scope);c['claim_authorization']='ACCEPTED_WITHIN_ASSIGNMENT_CATALOG'
        except ValueError as e:
            c.update(claim_authorization='REJECTED',reason=str(e));print(json.dumps(c,indent=2));return 2
    print(json.dumps(c,indent=2));return 0

if __name__=='__main__':raise SystemExit(main())
