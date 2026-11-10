"""Post hoc boundary sensitivity; does not replace the frozen baseline result.
Frozen grid after observing baseline: symmetric padding 0,3,6,9; deltas
0,9,18,27,36,45,54. Same raw pairs, admission, references and engine.
All windows reported; no winner selected as independent validation.
"""
import pathlib,json,gzip,hashlib,collections
from recompute_panel import encode,classify_intervals,maximize,NEG,solve
import optimal_lattice as lattice
ROOT=pathlib.Path(__file__).resolve().parent
def digest(b):return hashlib.sha256(b).hexdigest()
def main():
    freeze=json.loads((ROOT/'P1_FREEZE.json').read_text())
    old=json.loads((ROOT/'P1_SEGMENTS.json').read_text())
    refs=freeze['references'];targets=freeze['targets'];deltas=freeze['deltas']
    pads=[0,3,6,9]
    protocol={'scope':__doc__,'padding':pads,'parent_freeze_sha256':digest((ROOT/'P1_FREEZE.json').read_bytes()),
              'target_semantics':'Exact local string extended by exact observed admitted flanks; reference-projected event boundaries differ.'}
    (ROOT/'BOUNDARY_PROTOCOL.json').write_text(json.dumps(protocol,indent=2)+'\n')
    byhash={};witnesses=[]
    for idx,(q,prior) in enumerate(old.items()):
        cases={}
        for pad in pads:
            padded={name:(freeze['outer_left'][-pad:] if pad else '')+t+freeze['outer_right'][:pad] for name,t in targets.items()}
            lab=classify_intervals(q,padded);deficits=[]
            for ref in refs:
                opt,mx,_,_,suffix=maximize(encode(ref['sequence']),encode(q),ref['lo']-pad,ref['hi']+pad,lab,max(deltas))
                assert opt==suffix[0,0,0] and all(max(x)==opt for x in mx)
                deficits.append([[None if v==NEG else int(opt-v) for v in row] for row in mx])
            calls=[]
            for ti in range(2):
                tc=[]
                for delta in deltas:
                    sets=[{k for k,v in zip(('NE','E','U'),d[ti]) if v is not None and v<=delta} for d in deficits]
                    tc.append(next(iter(sets[0])) if len(sets[0])==1 and sets[0]==sets[1] else 'U')
                for a,b in zip(tc,tc[1:]):assert a==b or b=='U'
                calls.append(tc)
            if pad==0:assert calls==prior['calls'] and deficits==prior['deficits']
            cases[pad]={'calls':calls,'deficits':deficits}
        byhash[digest(q.encode())]=cases
        if idx%100==0:print('SEGMENTS',idx,len(old),flush=True)
    for pad in pads:
        for name,target in targets.items():
            q=freeze['outer_left']+target+freeze['outer_right'];assert q in old
            for ref in refs:
                case={'reference_sequence':ref['sequence'],'observed_sequence':q,'lo':ref['lo']-pad,'hi':ref['hi']+pad,
                    'target_kind':'HEK3_DNA','programmed_local':(freeze['outer_left'][-pad:] if pad else '')+target+freeze['outer_right'][:pad],
                    'unedited_local':ref['sequence'][ref['lo']-pad:ref['hi']+pad]}
                opt,ws=solve(case,lattice,cutoff=54)
                witnesses.append({'padding':pad,'target':name,'reference':ref['name'],'optimum':opt,'witnesses':ws})
    counts=collections.Counter();totals=collections.Counter()
    with gzip.open(ROOT/'P1_PAIRS.jsonl.gz','rt') as f:
        for line in f:
            row=json.loads(line);typ=row['library'];totals[typ]+=1
            for pad in pads:
                cs=[byhash[h][pad]['calls'] if h else [['U']*len(deltas) for _ in targets] for h in row['segment_hashes']]
                for ti,name in enumerate(targets):
                    for di,d in enumerate(deltas):
                        a,b=cs[0][ti][di],cs[1][ti][di];call=a if a==b else 'U'
                        counts[typ,pad,name,d,call]+=1
                        if a!='U' and b!='U' and a!=b:counts[typ,pad,name,d,'mate_conflict']+=1
    for typ,n in totals.items():
        for pad in pads:
            for name in targets:
                for d in deltas:assert sum(counts[typ,pad,name,d,c] for c in ['E','NE','U'])==n
    baseline=json.loads((ROOT/'P1_RESULT.json').read_text())
    for r in baseline['counts']:assert counts[r['library'],0,r['target'],r['delta'],r['call']]==r['pairs']
    result={'protocol_sha256':digest((ROOT/'BOUNDARY_PROTOCOL.json').read_bytes()),'totals':dict(totals),
         'baseline_reproduced':True,'witness_checks':len(witnesses),
         'counts':[{'library':t,'padding':p,'target':n,'delta':d,'call':c,'pairs':v} for (t,p,n,d,c),v in sorted(counts.items())]}
    for name,obj in [('BOUNDARY_RESULT.json',result),('BOUNDARY_SEGMENTS.json',byhash),('BOUNDARY_WITNESSES.json',witnesses)]:
        (ROOT/name).write_text(json.dumps(obj,indent=2)+'\n')
    for p in pads:
        for name in targets:
            print('RESULT',p,name,[(d,counts['DNA',p,name,d,'E'],counts['cDNA',p,name,d,'E']) for d in deltas],flush=True)
if __name__=='__main__':main()
