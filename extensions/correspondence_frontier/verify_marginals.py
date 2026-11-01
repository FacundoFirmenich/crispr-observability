"""Verify sharp aggregate agreement bounds by integer min-cost transport.

Matrices are mathematical completions of published marginals, not patient/read
truth labels. No particular completion is asserted to be the observed one.
"""
from pathlib import Path
import json
def transport(rows,cols,sign):
 n=len(rows);source=2*n;sink=source+1;g=[[] for _ in range(sink+1)]
 def edge(u,v,c,w):
  g[u].append([v,c,w,len(g[v])]);g[v].append([u,0,-w,len(g[u])-1]);return len(g[u])-1
 for i,v in enumerate(rows):edge(source,i,v,0)
 for j,v in enumerate(cols):edge(n+j,sink,v,0)
 ids=[[edge(i,n+j,sum(rows),sign*int(i==j)) for j in range(n)] for i in range(n)]
 flow=0
 while flow<sum(rows):
  d=[10**9]*len(g);d[source]=0;prev={}
  for _ in range(len(g)-1):
   changed=False
   for u,es in enumerate(g):
    for k,(v,c,w,rev) in enumerate(es):
     if c and d[u]+w<d[v]:d[v]=d[u]+w;prev[v]=(u,k);changed=True
   if not changed:break
  assert sink in prev
  v=sink;amount=sum(rows)-flow
  while v!=source:u,k=prev[v];amount=min(amount,g[u][k][1]);v=u
  v=sink
  while v!=source:
   u,k=prev[v];e=g[u][k];e[1]-=amount;g[v][e[3]][1]+=amount;v=u
  flow+=amount
 matrix=[[sum(rows)-g[i][ids[i][j]][1] for j in range(n)] for i in range(n)]
 assert [sum(row) for row in matrix]==rows
 assert [sum(matrix[i][j] for i in range(n)) for j in range(n)]==cols
 return sum(matrix[i][i] for i in range(n))
def verify(path):
 data=json.loads(Path(path).read_text());human=data['marginals']['human'];keys=list(human)
 for b in data['sharp_agreement_bounds_from_marginals_only']:
  rows=[human[k] for k in keys];cols=[data['marginals'][b['tool']][k] for k in keys]
  assert transport(rows,cols,1)==b['minimum_possible_agreements']
  assert transport(rows,cols,-1)==b['maximum_possible_agreements']
 print('Exact integer transport verifies both sharp endpoints for both published marginal tables.')
if __name__=='__main__':
 import sys
 verify(sys.argv[1])
