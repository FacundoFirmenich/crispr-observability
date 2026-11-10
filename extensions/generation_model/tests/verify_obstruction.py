"""Independent model consequence checks; not a replacement for the oracle."""
import itertools,json,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
from crispr_extended import minimum_costs
n=0
for h in map(''.join,itertools.product('AC',repeat=4)):
    for q in map(''.join,itertools.product('AC',repeat=4)):
        costs=minimum_costs(h,1,3,q,'AC',4)
        for lab in ('E','NE'):
            c=costs[lab]
            if c is not None and c<=2:assert costs['U'] is not None and costs['U']<=c+2;n+=1
result={'status':'PASS','finite_class_inequalities_checked':n,'claim':'c_U <= c_determinate + 2 under the declared generation grammar'}
(R/'results/obstruction_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
