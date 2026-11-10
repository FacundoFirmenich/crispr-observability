"""Independent bounded exhaustive software tests, never biological observations.

Enumerate legal alignment paths directly, without DP or production projection.
Compare exact scores, all-optimal image sets, dual boundary transport and
selected uncensored target maxima. Output is a new receipt, not a source edit.
"""
from pathlib import Path
import itertools, json, hashlib, random, time, sys
import optimal_lattice as production

R=Path(__file__).resolve().parent
OUT=R/'REVIEW_MP21'
COMP=str.maketrans('ACGTN-','TGCAN-')
def rc(x): return x.translate(COMP)[::-1]

def enumerate_paths(ref,read):
    n,m=len(ref),len(read)
    def walk(i,j,last,score,columns):
        if (i,j)==(n,m):
            yield score,columns
            return
        if i<n and j<m:
            a,b=ref[i],read[j]
            s=(5 if a==b else -4) if 'N' not in (a,b) else (-1 if a==b else -2)
            yield from walk(i+1,j+1,'M',score+s,columns+[(i,j)])
        if j<m and last in ('M','I'):
            cost=-2 if last=='I' or i==0 or i==n or j+1==m else -20
            yield from walk(i,j+1,'I',score+cost,columns+[(None,j)])
        if i<n and last in ('M','D'):
            cost=-2 if last=='D' or j==0 or j==m or i+1==n else -20
            yield from walk(i+1,j,'D',score+cost,columns+[(i,None)])
    return list(walk(0,0,'M',0,[]))

def oracle_projection(read,cols,lo,hi,right=False):
    positions={i:k for k,(i,j) in enumerate(cols) if i is not None}
    for flank in (lo-1,hi):
        j=cols[positions[flank]][1]
        if j is None or read[j] not in 'ACGT': return None
    left=positions[lo] if right else positions[lo-1]+1
    stop=positions[hi] if right else positions[hi-1]+1
    chars=[read[j] for i,j in cols[left:stop] if j is not None]
    return ''.join(chars) if chars and all(c in 'ACGT' for c in chars) else None

def image(result):
    return set(result['valid_coding_strings']) | ({None} if result['some_optimal_projection_unresolved'] else set())

def main():
    OUT.mkdir(exist_ok=True)
    started=time.monotonic(); counters=dict(sequence_pairs=0,paths_enumerated=0,window_contracts=0,dual_transport_checks=0,target_maxima_checks=0,modulo_checks=0)
    seqs=lambda alphabet,length: map(''.join,itertools.product(alphabet,repeat=length))
    pairs=[(r,q) for n in range(3,6) for r in seqs('AC',n) for m in range(6) for q in seqs('AC',m)]
    pairs += [(r,q) for r in seqs('AC',3) for m in range(1,5) for q in seqs('ACN',m) if 'N' in q]
    rng=random.Random(20261110)
    pairs += [(''.join(rng.choices('ACGT',k=rng.randrange(3,8))),''.join(rng.choices('ACGTN',k=rng.randrange(0,8)))) for _ in range(80)]
    from recompute_panel import solve
    for ci,(ref,read) in enumerate(pairs):
        paths=enumerate_paths(ref,read);best=max(s for s,c in paths);optimal=[c for s,c in paths if s==best]
        counters['sequence_pairs']+=1;counters['paths_enumerated']+=len(paths)
        if 'N' not in read:
            assert all((s+2*(len(ref)+len(read)))%9==0 for s,c in paths)
            counters['modulo_checks']+=len(paths)
        for lo in range(1,len(ref)-1):
            for hi in range(lo+1,len(ref)):
                for right in (False,True):
                    expected={oracle_projection(read,c,lo,hi,right) for c in optimal}
                    got=production.audit(ref,read,lo,hi,enumeration_limit=0,boundary_convention='right_insertions' if right else 'left_insertions')
                    assert got['optimal_integer_objective']==best and got['optimal_path_count']==len(optimal),(ref,read)
                    assert image(got)==expected,(ref,read,lo,hi,right,image(got),expected)
                    rev=production.audit(rc(ref),rc(read),len(ref)-hi,len(ref)-lo,enumeration_limit=0,boundary_convention='left_insertions' if right else 'right_insertions')
                    assert rev['optimal_integer_objective']==best and image(rev)=={rc(x) if x is not None else None for x in expected},(ref,read,lo,hi,right,'reverse')
                    counters['window_contracts']+=1;counters['dual_transport_checks']+=1
                if ci%79==0:
                    target=ref[lo:hi];alternative='A' if target!='A' else 'C'
                    case=dict(target_kind='HEK3_DNA',reference_sequence=ref,observed_sequence=read,lo=lo,hi=hi,programmed_local=target,unedited_local=alternative)
                    # The production maximum routine assumes nonempty read input;
                    # empty reads are still covered above as unresolved DAG cases.
                    if not read: continue
                    opt,w=solve(case,production,-1);assert opt==best
                    for name,event in [('programmed_local',target),('unedited_local',alternative)]:
                        expected={}
                        for score,cols in paths:
                            s=oracle_projection(read,cols,lo,hi)
                            lab=2 if s is None else int(s==event)
                            expected[lab]=max(expected.get(lab,-10**12),score)
                        actual={k:v['score'] for k,v in w[name].items()}
                        assert actual==expected,(ref,read,lo,hi,name,actual,expected)
                        counters['target_maxima_checks']+=1
        if ci%500==0:print('oracle progress',ci,len(pairs),flush=True)
    result=dict(status='PASS',scope='Bounded exhaustive computational tests; artificial test strings are not biological data or additional study observations.',alphabet='Exhaustive AC reference lengths 3..5 and read lengths 0..5; N-containing ACN reads to length 4 for length-3 AC references; 80 seeded ACGT/ACGTN cases.',**counters,seconds=round(time.monotonic()-started,3),source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),R/'optimal_lattice.py',R/'recompute_panel.py']})
    (OUT/'INDEPENDENT_ORACLE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
