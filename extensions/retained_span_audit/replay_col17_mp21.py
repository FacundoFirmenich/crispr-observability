"""Read-only replay of retained COL17A1 spans and complete pair-ledger counts.

--mode example: authentic Clone3 sequence, exact deletion placements, margins.
--mode counts: every pair and target/deficit count, including strict admission.
--mode full: counts plus independently reconstructed attaining score witnesses.
Input hashes are checked before execution. Only --output is written.
"""
import argparse,collections,gzip,hashlib,json,pathlib,time
from recompute_panel import solve
import optimal_lattice as lattice

R=pathlib.Path(__file__).resolve().parent
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--mode',choices=['example','counts','full'],default='example');ap.add_argument('--output',type=pathlib.Path,required=True);args=ap.parse_args()
    assert not args.output.exists(),'Output must be a new file'
    manifest=json.loads((R/'MANIFEST.json').read_text())
    for x in manifest['files']:
        p=(R/x['path']).resolve();assert p.is_relative_to(R.resolve()) and sha(p)==x['sha256'],x['path']
    began=time.monotonic();P=R/'NONMODAL_MP16';A=R/'ADMISSION_MP15B';protocol=json.loads((P/'PROTOCOL.json').read_text())
    result={'mode':args.mode,'input_files_verified':len(manifest['files']),'biological_accuracy':'NOT_ESTIMATED'}
    if args.mode=='example':
        j=json.loads((P/'CLONE3_RECONCILIATION.json').read_text());ref=j['observed_reference_JEB'];q=j['observed_clone3_sequence']
        starts=[i for i in range(len(ref)-5) if ref[:i]+ref[i+6:]==q];assert starts==[40,41]
        targets=protocol['targets']['3'];cases=[]
        for ref in protocol['references']:
            case=dict(target_kind='HEK3_DNA',reference_sequence=ref['sequence'],observed_sequence=q,lo=ref['lo'],hi=ref['hi'],programmed_local=targets[0]['sequence'][10:-10],unedited_local=targets[-1]['sequence'][10:-10])
            opt,w=solve(case,lattice,54);deficits={str(k):opt-v['score'] for k,v in w['programmed_local'].items()}
            assert deficits['1']==0 and min(v for k,v in deficits.items() if k!='1')==45
            cases.append({'reference':ref['name'],'target_class_deficits':deficits,'witnesses':w['programmed_local']})
        result.update(exact_deletion_starts=starts,net_length=-6,source='Zenodo 6242279; Clone3 observed modal span; original source Table S2 +6 retained as inconsistent with Figure 6A.',cases=cases)
    else:
        physical=collections.Counter();clones=[];allpairs=0;totalw=0;totalcases=0
        for clone,targets in sorted(protocol['targets'].items(),key=lambda x:int(x[0])):
            seqs=json.loads((A/f'clone{clone}_segments.json').read_text());scores=json.loads((P/f'clone{clone}_scores.json').read_text());expected=json.loads((P/f'clone{clone}_counts.json').read_text())
            assert sha(P/f'clone{clone}_scores.json')==expected['scores_sha256']
            if args.mode=='full':
                for h,q in seqs.items():
                    for ref,expected_score in zip(protocol['references'],scores[h]['reference_scores'],strict=True):
                        case=dict(target_kind='HEK3_DNA',reference_sequence=ref['sequence'],observed_sequence=q,lo=ref['lo'],hi=ref['hi'],programmed_local=targets[0]['sequence'][10:-10],unedited_local=targets[-1]['sequence'][10:-10])
                        opt,w=solve(case,lattice,54);assert opt==expected_score['optimum'];totalcases+=1
                        for ti,name in enumerate(('programmed_local','unedited_local')):
                            actual=[opt-w[name][k]['score'] if k in w[name] else None for k in range(3)]
                            assert actual==expected_score['deficits'][ti],(clone,h,ref['name'],name)
                            totalw+=len(w[name])
            counts=collections.Counter();raw=collections.Counter()
            with gzip.open(A/f'clone{clone}_pairs.jsonl.gz','rt') as f:
                for line in f:
                    row=json.loads(line);lib=row['library'];raw[lib]+=1;hs=set(filter(None,row['segments']+[row['merged_segment']]))
                    group='no_span' if not hs else 'modal_only' if all(seqs[h] in [t['sequence'] for t in targets] for h in hs) else 'nonmodal'
                    strict=row['event_relative_JEB'] is not None
                    for d in protocol['deltas']:
                        calls=[]
                        for ti,t in enumerate(targets):
                            span_calls=[]
                            for h in hs:
                                labels=[]
                                for sc in scores[h]['reference_scores']:
                                    admitted=[lab for lab,z in zip(('NE','E','U'),sc['deficits'][ti]) if z is not None and z<=d]
                                    labels.append(admitted)
                                span_calls.append(labels[0][0] if len(labels[0])==1 and labels[0]==labels[1] else 'U')
                            call=span_calls[0] if span_calls and span_calls[0]!='U' and all(v==span_calls[0] for v in span_calls) else 'U'
                            counts[lib,t['label'],d,group,call]+=1;calls.append(call)
                        pc='any_E' if 'E' in calls else 'all_NE' if calls and all(v=='NE' for v in calls) else 'unresolved'
                        for policy in ('any_available_span','MP15B_gate'):
                            physical[lib,policy,d,group,pc if policy=='any_available_span' or strict else 'unresolved']+=1
            assert dict(raw)==expected['raw']
            expected_counts={(x['library'],x['target'],x['delta'],x['group'],x['call']):x['pairs'] for x in expected['counts']}
            # Frozen tables explicitly include zero cells; Counter omits them.
            assert dict(counts)=={k:v for k,v in expected_counts.items() if v},(clone,'counts')
            allpairs+=sum(raw.values());clones.append(dict(clone=clone,spans=len(seqs),pairs=sum(raw.values()),count_rows=len(counts)))
            print('verified',clone,sum(raw.values()),len(seqs),flush=True)
        exp=json.loads((P/'ADMISSION_SENSITIVITY.json').read_text())['counts']
        expected_physical={(x['library'],x['policy'],x['delta'],x['group'],x['call']):x['pairs'] for x in exp}
        assert dict(physical)==expected_physical and allpairs==439589
        strict36={}
        for lib in ('DNA','cDNA'):
            strict36[lib]={call:sum(v for k,v in physical.items() if k[0]==lib and k[1]=='MP15B_gate' and k[2]==36 and k[4]==call) for call in ('any_E','all_NE','unresolved')}
        result.update(clones=clones,raw_pairs=allpairs,reference_cases=totalcases,attaining_witnesses=totalw,strict_deficit36=strict36)
    result.update(status='PASS',seconds=round(time.monotonic()-began,3))
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','clones')}))
if __name__=='__main__':main()
