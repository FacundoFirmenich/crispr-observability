"""Independent arithmetic/coverage checks over actual paired-policy outputs."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,json,hashlib

def must(x,msg):
    if not x:raise ValueError(msg)
def main(directory):
    d=Path(directory);load=lambda n:json.loads((d/n).read_text())
    rows=load('COUNTS.json');idx={(r['run_accession'],r['trim'],r['policy'],r['observer'],r['event']):r for r in rows}
    must(len(rows)==960 and len(idx)==len(rows),'unexpected count coverage')
    rules=0
    for key,r in idx.items():
        n,e,ne,u=(r[x] for x in ['N','E','NE','U']);must(n>0 and e+ne+u==n,'partition')
        must(F(r['L_exact'])==F(e,n) and F(r['H_exact'])==1-F(ne,n),'exact endpoints')
        base=idx[key[:2]+('structural_only',)+key[3:]]
        both=idx[key[:2]+('common_and_own_meanQ30',)+key[3:]]
        must(r['N']==base['N'] and r['E']<=base['E'] and r['NE']<=base['NE'],'quality cannot add labels')
        must(both['E']<=r['E'] and both['NE']<=r['NE'],'intersection cannot add labels')
        # Removing reliable labels widens intervals: lower decreases and upper increases.
        must(F(r['L_exact'])<=F(base['L_exact']) and F(r['H_exact'])>=F(base['H_exact']),'interval widening')
        rules+=1
    contrasts=load('CONTRASTS.json')
    for r in contrasts:
        t=idx[r['treated'],r['trim'],r['policy'],r['observer'],r['event']]
        c=idx[r['control'],r['trim'],r['policy'],r['observer'],r['event']]
        low=F(t['E'],t['N'])+F(c['NE'],c['N'])-1
        high=1-F(t['NE'],t['N'])-F(c['E'],c['N'])
        must((low,high)==(F(r['L_exact']),F(r['H_exact'])),'contrast via complementary counts')
        must(r['sign']==('positive' if low>0 else 'negative' if high<0 else 'includes_zero'),'sign')
    pairs=load('PAIRS.json');loss=load('LOSS_INTERSECTIONS.json');totals=Counter()
    for r in loss:
        key=tuple(r[x] for x in ['run_accession','trim','policy','observer_A','observer_B','event','B_label'])
        totals[key]+=r['count']
    for r in pairs:
        key=tuple(r[x] for x in ['run_accession','trim','policy','observer_A','observer_B','event'])
        must(sum(r['transitions'].values())==r['N'],'pair partition')
        for label in ['E','N']:
            must(totals[key+(label,)]==r['transitions']['U>'+label],'loss diagnostics partition rescued labels')
    decomp=load('DECOMPOSITIONS.json')
    for r in decomp:
        g=r['gap_exact']
        must(F(g['own_meanQ30'])==F(g['structural_only'])+F(r['own_quality_net_exact']),'quality structural identity')
    v=load('VALIDATION.json');must(v['frozen_regression_checks']==480 and v['raw_reads']==429812 and v['unique_streams']==15,'baseline scope')
    result={'status':'PASS','actual_count_cases':rules,'actual_contrasts':len(contrasts),'paired_tables':len(pairs),
            'loss_partition_checks':2*len(pairs),'quality_decompositions':len(decomp),'frozen_regression_checks':v['frozen_regression_checks'],
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'independent arithmetic and actual-output invariants; not second FASTQ classification or biological truth validation'}
    (d/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory');main(p.parse_args().directory)
