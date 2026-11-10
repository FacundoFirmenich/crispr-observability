"""Independent substring enumeration plus full edit distance; no scientific data."""
import itertools,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from crispr_extended.calibration import fragment_distance,block_quantile,prediction_set
from crispr_extended.core import minimum_costs

def edit(a,b):
    d=list(range(len(b)+1))
    for i,x in enumerate(a):
        v=[i+1]
        for j,y in enumerate(b):v.append(min(v[-1]+1,d[j+1]+1,d[j]+(x!=y)))
        d=v
    return d[-1]

def main():
    count=0
    for n in range(1,5):
        for h in map(''.join,itertools.product('AC',repeat=n)):
            for m in range(1,5):
                for r in map(''.join,itertools.product('AC',repeat=m)):
                    expected=min(edit(h[i:j],r) for i in range(n) for j in range(i+1,n+1))
                    for cap in (0,1,2,4):assert fragment_distance(h,r,cap)==(expected if expected<=cap else None),(h,r,cap)
                    if n==4:
                        values=minimum_costs(h,1,3,r,'AC',4).values()
                        assert min(x for x in values if x is not None)==expected
                    count+=1
    for read in ('ANCA','NNNN','ACNA'):
        assert fragment_distance('AACA',read,4)==min(edit('AACA'[i:j],read) for i in range(4) for j in range(i+1,5))
    q=block_quantile(list(range(256)),1,20);assert q['rank']==245 and q['q']==244
    assert block_quantile([0]*256,1,10)['rank']==232
    assert block_quantile([0]*8,1,20)['infinity']
    assert prediction_set({'0':None,'1':0},block_quantile([None]*256,1,20))==['0','1']
    assert prediction_set({'0':3,'1':None},{'infinity':False,'q':2})==[]
    print(json.dumps({'status':'PASS','independent_substring_cases':count+3,'caps':[0,1,2,4],'full_core_comparison_cases':480,'quantile_and_infinity_checks':'PASS'}))

if __name__=='__main__':main()
