"""Exact projection sets over every optimal path of a declared affine alignment DAG.

The objective reproduces the zero-gap-incentive recurrence and terminal-gap
convention of CRISPResso2 2.3.4. It does not reproduce its single-traceback tie
choice or confuse its reported identity percentage with this integer objective.
No biological truth or out-of-model optimality is claimed.
"""
M,I,J=0,1,2
NEG=-10**12
def substitution(a,b):return 5 if a==b and a!='N' else -2 if 'N' in (a,b) and a!=b else -1 if a==b=='N' else -4
def project(aref,aseq,lo,hi,boundary_convention='left_insertions'):
 columns=[k for k,c in enumerate(aref) if c!='-']
 if any(aseq[columns[k]] not in 'ACGT' for k in [lo-1,hi]):return None
 s=(aseq[columns[lo-1]+1:columns[hi-1]+1] if boundary_convention=='left_insertions' else aseq[columns[lo]:columns[hi]]).replace('-','')
 return s if s and not set(s)-set('ACGT') else None
def rescore(aref,aseq,gap_open=-20,gap_extend=-2):
 assert len(aref)==len(aseq);n=len(aref.replace('-',''));m=len(aseq.replace('-',''));i=j=0;prev=M;score=0
 for a,b in zip(aref,aseq):
  assert a!='-' or b!='-'
  if a!='-' and b!='-':i+=1;j+=1;score+=substitution(a,b);prev=M
  elif a=='-':
   j+=1;assert prev in [M,I] or i==0;score+=gap_extend if i==0 or prev==I or i==n or j==m else gap_open;prev=I
  else:
   i+=1;assert prev in [M,J] or j==0;score+=gap_extend if j==0 or prev==J or i==n or j==m else gap_open;prev=J
 return score
def audit(reference,observed,lo,hi,gap_open=-20,gap_extend=-2,enumeration_limit=1000,reverse_domain=None,boundary_convention='left_insertions'):
 assert set(reference)<=set('ACGT') and set(observed)<=set('ACGTN') and 0<lo<hi<len(reference)
 n,m=len(reference),len(observed);w=m+1;size=(n+1)*w;score=[[NEG]*size for _ in range(3)];pred=[[0]*size for _ in range(3)];score[M][0]=0
 for j in range(1,m+1):score[I][j]=gap_extend*j;pred[I][j]=1<< (M if j==1 else I)
 for i in range(1,n+1):score[J][i*w]=gap_extend*i;pred[J][i*w]=1<<(M if i==1 else J)
 for i in range(1,n+1):
  for j in range(1,m+1):
   k=i*w+j;diag=k-w-1;left=k-1;up=k-w;opening=gap_extend if i==n or j==m else gap_open
   a,b=score[M][left]+opening,score[I][left]+gap_extend;best=max(a,b);score[I][k]=best;pred[I][k]=(1<<M if a==best else 0)|(1<<I if b==best else 0)
   a,b=score[M][up]+opening,score[J][up]+gap_extend;best=max(a,b);score[J][k]=best;pred[J][k]=(1<<M if a==best else 0)|(1<<J if b==best else 0)
   v=[score[s][diag] for s in range(3)];best=max(v);score[M][k]=best+substitution(reference[i-1],observed[j-1]);pred[M][k]=sum(1<<s for s in range(3) if v[s]==best)
 end=n*w+m;best=max(s[end] for s in score)
 def parents(s,k):
  if k==0:return []
  delta=w+1 if s==M else 1 if s==I else w
  return [(t,k-delta) for t in range(3) if pred[s][k]&(1<<t)]
 final=[(s,end) for s in range(3) if score[s][end]==best];good=set(final);stack=list(final)
 while stack:
  s,k=stack.pop()
  for p in parents(s,k):
   if p not in good:good.add(p);stack.append(p)
 ordered=sorted(good,key=lambda x:x[1]);ways={(M,0):1}
 for node in ordered:
  if node!=(M,0):ways[node]=sum(ways.get(p,0) for p in parents(*node))
 paths=sum(ways[node] for node in final);starts={};intervals=set();invalid=False
 left_valid={}
 for s,k in ordered:
  i,j=divmod(k,w)
  if i in [lo,hi+1] and (s==J or s==M and (j==0 or observed[j-1] not in 'ACGT')):invalid=True
  if boundary_convention=='right_insertions':
   if i==lo:
    left_valid[s,k]=(s==M and j>0 and observed[j-1] in 'ACGT') or (s==I and any(left_valid.get(p,False) for p in parents(s,k)))
   if not lo<i<=hi+1:continue
   if i==lo+1 and s!=I:bits=(1<<(j-1 if s==M else j)) if any(left_valid.get(p,False) for p in parents(s,k)) else 0
   else:
    bits=0
    for p in parents(s,k):bits|=starts.get(p,0)
   starts[s,k]=bits
   if i==hi+1 and s==M and j and observed[j-1] in 'ACGT':
    while bits:
     b=bits&-bits;start=b.bit_length()-1;bits-=b;coding=observed[start:j-1]
     if not coding or set(coding)-set('ACGT'):invalid=True
     else:intervals.add((start,j-1))
   continue
  if i<lo or i>hi:continue
  if i==lo and s==M:bits=1<<j if j and observed[j-1] in 'ACGT' else 0
  else:bits=0
  if not (i==lo and s==M):
   for p in parents(s,k):bits|=starts.get(p,0)
  starts[s,k]=bits
  if i!=hi or s==I or not bits:continue
  state=s;col=j;right_ok=False
  while col<m:
   nxt=(hi+1)*w+col+1
   if (M,nxt) in good and pred[M][nxt]&(1<<state) and observed[col] in 'ACGT':right_ok=True;break
   ins=hi*w+col+1
   if (I,ins) not in good or not pred[I][ins]&(1<<state):break
   state=I;col+=1
  if right_ok:
   while bits:
    b=bits&-bits;start=b.bit_length()-1;bits-=b;coding=observed[start:j]
    if not coding or set(coding)-set('ACGT'):invalid=True
    else:intervals.add((start,j))
 codings={observed[a:b] for a,b in intervals};enumerated=False;reverse_codes=set()
 if paths<=enumeration_limit:
  brute=set();count=0
  def walk(s,k,rr,ss):
   nonlocal count
   if k==0:
    count+=1;brute.add(project(rr[::-1],ss[::-1],lo,hi,boundary_convention))
    if reverse_domain is not None:
     comp=str.maketrans('ACGTN','TGCAN');reverse_codes.add(project(rr.translate(comp),ss.translate(comp),*reverse_domain))
    return
   i,j=divmod(k,w);a=reference[i-1] if s!=I else '-';b=observed[j-1] if s!=J else '-'
   for t,p in parents(s,k):walk(t,p,rr+a,ss+b)
  for s,k in final:walk(s,k,'','')
  assert count==paths and brute==codings|({None} if invalid else set()),'Interval propagation differs from complete actual-path enumeration'
  enumerated=True
 result={'optimal_integer_objective':best,'optimal_path_count':paths,'optimal_DAG_nodes':len(good),'valid_projection_intervals':[list(x) for x in sorted(intervals)],'valid_coding_strings':sorted(codings),'some_optimal_projection_unresolved':invalid,'all_optimal_strings_identical':not invalid and len(codings)==1,'all_paths_independently_enumerated':enumerated,'enumeration_limit':enumeration_limit,'scope':'All optimal paths for the declared affine objective, zero gap incentive and terminal convention; not all scoring models, alignments or biological truth.'}
 if reverse_domain is not None:
  assert enumerated,'Reverse-domain check needs full enumeration for this actual case'
  result['canonical_reverse_coding_strings']=sorted(x for x in reverse_codes if x is not None);result['canonical_reverse_some_unresolved']=None in reverse_codes
 return result
