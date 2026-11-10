"""Exhaustive finite-domain verification, not a biological experiment."""
import itertools,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from crispr_extended.core import minimum_costs,certify
from crispr_extended.oracle import enumerate_histories
def main():
    start=time.monotonic();count=histories=0
    for h in map(''.join,itertools.product('AC',repeat=4)):
        for m in (1,2,3,4):
            for q in map(''.join,itertools.product('AC',repeat=m)):
                for target in ('A','C','AC'):
                    actual=minimum_costs(h,1,3,q,target,2);expected,n=enumerate_histories(h,1,3,q,target,2)
                    assert actual==expected,(h,q,target,actual,expected)
                    count+=1;histories+=n
    for q in ('ANCA','NACA','ACNA','ACAN'):
        a=minimum_costs('AACA',1,3,q,'AC',3);b,n=enumerate_histories('AACA',1,3,q,'AC',3);assert a==b;(count:=count+1);histories+=n
    r=certify('ACGT','CG',[{'id':'h','sequence':'ACGT','lo':1,'hi':3,'allele':'ALT'}],(0,1,2))
    assert r['decisions'][0]['observed_tract_certificate']=='E'
    sets=[set(d['possible_observed_labels']) for d in r['decisions']];assert all(a<=b for a,b in zip(sets,sets[1:]))
    result={'status':'PASS','instances':count,'histories_enumerated':histories,'seconds':time.monotonic()-start,'scope':'Complete bounded input domain, not independent biological evidence'}
    Path('results').mkdir(exist_ok=True);Path('results/oracle_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
