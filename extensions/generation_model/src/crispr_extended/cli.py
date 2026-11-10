"""JSON/JSONL contract interface. Outputs preserve record identity."""
import argparse,json,sys
from .core import certify
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',help='JSONL requests, or - for stdin');p.add_argument('-o','--output',default='-');a=p.parse_args()
    inp=sys.stdin if a.input=='-' else open(a.input,encoding='utf8');out=sys.stdout if a.output=='-' else open(a.output,'w',encoding='utf8')
    try:
        for line in inp:
            if not line.strip():continue
            r=json.loads(line);result=certify(r['read'],r['target'],r['haplotypes'],r.get('budgets',(0,1,2,3,5,8,12)))
            out.write(json.dumps({'id':r['id'],**result})+'\n');out.flush()
    finally:
        if inp is not sys.stdin:inp.close()
        if out is not sys.stdout:out.close()
if __name__=='__main__':main()
