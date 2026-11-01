"""Exact scientific plots from the MP05 finite-record outputs; no simulated data."""
from pathlib import Path
from collections import Counter
import csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[2];D=R/'extensions/reproduction/derived';O=R/'figures/frontier';O.mkdir(exist_ok=True)
def rows(n):return list(csv.DictReader((D/n).open()))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','font.family':'DejaVu Sans'})
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(O/(name+'.'+ext),bbox_inches='tight',dpi=190)
 plt.close(fig)
names={1:'CFB-26-16',2:'CFB-26-02',3:'CFB-06-04'};col={'PE_selected_reference':'#146888','old_three_operator_witness':'#744585','PE_corroborated_witness':'#b96426'}
fig=plt.figure(figsize=(8.5,7.8));grid=fig.add_gridspec(2,2,height_ratios=[1.6,1]);ax=[fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,:])]
t=[r for r in rows('PE_TOOL_VS_EVENT_COMPARISON.csv') if r['condition']=='PE treated']
for r in t:
 p=int(r['preparation']);y=3-p
 ax[0].plot(float(r['Prime_edited_assignment_percent_of_all_aligned']),y+.14,'D',color='#595959',ms=5)
 lo=float(r['PE_projected_native_L_percent_raw_N']);hi=float(r['PE_projected_native_H_percent_raw_N']);ax[0].plot([lo,hi],[y-.14]*2,lw=4,color=col['PE_selected_reference']);ax[0].plot([lo,hi],[y-.14]*2,'|',color=col['PE_selected_reference'],ms=9)
ax[0].set(yticks=[2,1,0],yticklabels=list(names.values()),xlabel='Percent (denominators differ)',title='A  Assignment and local event');ax[0].set_xlim(10,42)
ax[0].legend(handles=[Line2D([],[],marker='D',ls='',color='#595959',label='Prime-edited / all aligned'),Line2D([],[],lw=4,color=col['PE_selected_reference'],label='Native-local hull / raw N')],loc='upper left',bbox_to_anchor=(-.03,-.20),fontsize=8,frameon=False)
c=[r for r in rows('PE_EVENT_CONTRASTS.csv') if r['control']=='mock' and r['event']=='native_local']
for r in c:
 p=int(r['preparation']);op=r['operator'];offset={'old_three_operator_witness':.19,'PE_selected_reference':0,'PE_corroborated_witness':-.19}[op];y=3-p+offset;lo=float(r['L_pp']);hi=float(r['H_pp']);ax[1].plot([lo,hi],[y,y],lw=2.5,color=col[op]);ax[1].plot([lo,hi],[y,y],'|',color=col[op],ms=7)
ax[1].axvline(0,color='#777777',lw=.8,ls=':');ax[1].set(yticks=[2,1,0],yticklabels=list(names.values()),xlabel='PE minus mock (pp)',title='B  Shared predicate, raw N');ax[1].set_xlim(-1,38)
ax[1].legend(handles=[Line2D([],[],lw=3,color=col[k],label=v) for k,v in [('old_three_operator_witness','Original T+M+P witness'),('PE_selected_reference','Dedicated PE projection'),('PE_corroborated_witness','Witness + PE correspondence')]],loc='upper left',bbox_to_anchor=(-.03,-.20),fontsize=8,frameon=False)
d=rows('PE_STRING_DISAGREEMENTS.csv');classes=Counter()
for r in d:
 k='Same two labels' if r['old_native']==r['PE_native'] and r['old_exact']==r['PE_exact'] else 'Native label only' if r['old_exact']==r['PE_exact'] else 'Native + exact labels';classes[k]+=1
labels=['Same two labels','Native label only','Native + exact labels'];vals=[classes[k] for k in labels];bars=ax[2].barh(range(3),vals,color=['#6f7b87','#c78635','#ad4f48'],height=.55)
for bar,v in zip(bars,vals):ax[2].text(v+.6,bar.get_y()+bar.get_height()/2,str(v),va='center')
ax[2].set(yticks=range(3),yticklabels=labels,xlabel='Records with unequal coding strings',title='C  89 strings disagree; 36 labels change',xlim=(0,62));ax[2].invert_yaxis()
fig.subplots_adjust(left=.20,right=.97,top=.95,bottom=.10,wspace=.75,hspace=1.05);save(fig,'fig5_PE_comparison')
# Resolution-cost plot; every mark represents a declared string certificate.
cs=[r for r in rows('CORRESPONDENCE_ABLATION_CONTRASTS.csv') if r['control']=='mock' and r['event']=='native_local'];subs=['target','mutant_projection','programmed_projection','target+mutant_projection','target+programmed_projection','mutant_projection+programmed_projection','target+mutant_projection+programmed_projection'];short=['T','M','P','T+M','T+P','M+P','T+M+P']
fig,ax=plt.subplots(figsize=(8.2,4.6))
for prep,color in [(1,'#146888'),(2,'#744585'),(3,'#b96426')]:
 ys=[float(next(r for r in cs if int(r['preparation'])==prep and r['certificate']==s)['width_pp']) for s in subs];ax.plot([i+(prep-2)*.13 for i in range(7)],ys,'o',color=color,label=names[prep],ms=5)
ax.set(xticks=range(7),xticklabels=short,xlabel='Required identical-string operators',ylabel='Native PE minus mock hull width (pp)',title='Correspondence requirements have an explicit resolution cost')
ax.legend(frameon=False);fig.tight_layout();save(fig,'figS3_correspondence_cost')
used=['PE_TOOL_VS_EVENT_COMPARISON.csv','PE_EVENT_CONTRASTS.csv','PE_STRING_DISAGREEMENTS.csv','CORRESPONDENCE_ABLATION_CONTRASTS.csv']
(O/'SOURCE_IDENTITIES.json').write_text(json.dumps({'edition':'MP05','figure_type':'Deterministic scientific plots of actual finite-record results; not AI-generated artwork or simulated biological data.','sources':{n:hashlib.sha256((D/n).read_bytes()).hexdigest() for n in used},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print(json.dumps({'status':'GENERATED','figures':['fig5_PE_comparison','figS3_correspondence_cost']}))
