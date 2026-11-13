from pathlib import Path
from fractions import Fraction as F
import argparse,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path);output=parser.parse_args().output
assert not output.exists(), 'Choose a new output path'
d=json.loads((HERE/'results/COST_FRONTIER.json').read_bytes());p=json.loads((HERE/'results/LOCAL_POLICY.json').read_bytes());colors=['#426db5','#b879a7','#9d88d1','#593e99','#df996d'];modes=['baseline','majority','pooled_empirical','local_empirical','local_beta11'];labels=['Baseline','Always majority','Pooled routing','Local routing','Local Beta(1,1)']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'#273b62','axes.labelcolor':'#273b62','xtick.color':'#273b62','ytick.color':'#273b62','axes.edgecolor':'#aab7d5','pdf.fonttype':42})
fig,(ax,bx)=plt.subplots(1,2,figsize=(11.7,5.0),gridspec_kw={'width_ratios':[1.35,1]})
for j,mode in enumerate(modes):
 for i,x in enumerate(d['intervals']):
  a=max(2,float(F(x['lower_open'])));b=min(12,float(F(x['upper_open']))) if x['upper_open'] else 12
  if a>=b:continue
  c=x['counts'][mode];ax.plot([a,b],[a*c['D']+c['U'],b*c['D']+c['U']],color=colors[j],lw=2 if mode=='local_empirical' else 1.4,linestyle='-' if j in [0,3] else '--',label=labels[j] if a==2 else None)
ax.axvspan(17/3,20/3,color='#b0cff0',alpha=.3);ax.set(xlabel='Discordance cost / abstention cost',ylabel='Held-out loss',xlim=(2,12),title='A  Decision regimes, including adverse transfer');ax.legend(fontsize=8,loc='upper left');ax.grid(axis='y',color='#e3e8f3',lw=.6)
rows=sorted([x for x in p['local_decisions'] if x['policy']=='local_empirical' and x['cost_ratio']==6],key=lambda x:(x['held_out_donor'],x['tract']));v=[[x['counts'][k] for k in ['C','D','U']] for x in rows];cmap=LinearSegmentedColormap.from_list('blueviolet',['#f6f3fb','#b6c9ee','#5c58a4']);bx.imshow(v,cmap=cmap,vmin=0,vmax=32,aspect='auto')
for i,row in enumerate(v):
 for j,n in enumerate(row):bx.text(j,i,str(n),ha='center',va='center',color='white' if n>17 else '#273b62')
bx.set_xticks([0,1,2],['Concordant','Discordant','Unresolved']);bx.set_yticks(range(len(rows)),[x['held_out_donor']+(' ≥15 nt' if x['tract']=='ge15' else ' <15 nt') for x in rows]);bx.set_title('B  Local outcomes at cost ratio six');bx.tick_params(length=0);fig.tight_layout();fig.savefig(output);plt.close(fig)
