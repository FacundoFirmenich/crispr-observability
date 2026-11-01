"""Deterministic article figures from real, frozen MP06 scientific outputs."""
from pathlib import Path
import csv,json,os
os.environ.setdefault("SOURCE_DATE_EPOCH","1793491200")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Rectangle
R=Path(__file__).resolve().parents[2];E=R/'extensions/correspondence_frontier/expected';O=R/'figures/correspondence_frontier'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'crispr-p1-mp06','pdf.fonttype':42,'font.family':'DejaVu Sans'})
BLUE='#185f8d';GREEN='#24785c';RED='#b34c46';GREY='#a7b2bd'
def load(p):return json.loads(p.read_text())
def label(ax,letter,title):ax.set_title(title,loc='left',fontsize=10,pad=10);ax.text(-.14,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=12)
def export(fig,name):
 O.mkdir(parents=True,exist_ok=True)
 for ext in ['svg','pdf','png']:fig.savefig(O/(name+'.'+ext),dpi=180,bbox_inches='tight',metadata={'Creator':'CRISPR P1; deterministic scientific plotting'})
 plt.close(fig)
def polygon(ax,p,n,targets):
 points=[(100*x/n,100*y/n) for x,y in p];xs,ys=zip(*points);x0,x1=min(xs),max(xs);y0,y1=min(ys),max(ys)
 ax.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fill=False,ls='--',edgecolor=GREY,lw=1.5,label='Marginal rectangle'))
 ax.add_patch(Polygon(points,closed=True,alpha=.25,facecolor=BLUE,edgecolor=BLUE,lw=2,label='Joint convex hull'))
 ax.plot([*xs,xs[0]],[*ys,ys[0]],color=BLUE,lw=1.5)
 padx=(x1-x0)*.12;pady=(y1-y0)*.12;ax.set_xlim(x0-padx,x1+padx);ax.set_ylim(y0-pady,y1+pady);ax.set_xlabel(targets[0]+' (% raw N)');ax.set_ylabel(targets[1]+' (% raw N)');ax.legend(fontsize=7,loc='upper left')
def main():
 m=load(E/'mechanism/MECHANISM_AND_CERTIFICATE_VERIFICATION.json');l=load(E/'cftr_lattice/LATTICE_VERIFICATION.json')
 fig,ax=plt.subplots(2,2,figsize=(7.6,6.9),layout='constrained')
 a=ax[0,0];vals=[l[k] for k in ['whole_string_invariant_cases','native_invariant_cases','exact_invariant_cases']];a.bar(['Exact\nstring','Native\npeptide','Programmed\nDNA'],vals,color=[GREY,GREEN,BLUE]);a.set_ylim(0,420);a.set_ylabel('Invariant cases / 398');label(a,'A','All optimal paths: CFTR')
 for i,v in enumerate(vals):a.text(i,v+6,str(v),ha='center')
 a.text(.02,.95,'199 actual reads × 2 references\nPurpose-enriched panel',transform=a.transAxes,va='top',fontsize=8)
 a=ax[0,1];a.bar(['Orientation\nchange','PE gap\nincentive'],[44,45],color=[BLUE,RED]);a.set_ylim(0,55);a.set_ylabel('Records among 89 discrepancies');label(a,'B','Fixed partial factorial')
 for i,v in enumerate([44,45]):a.text(i,v+1,str(v),ha='center')
 a.text(.04,.95,'PE incentive 0 = competing RC\nfor all 89 records',transform=a.transAxes,va='top',fontsize=8)
 a=ax[1,0];u=m['unknown_counts']['F9'];events=['native','exact_programmed','native_nonexact']
 for i,k in enumerate(events):
  v=u['target_relative'][k];s=u['string_witness'][k];a.plot([v,s],[i,i],color=GREY,lw=2);a.scatter([s],[i],color=GREY,label='Identical-string requirement' if i==0 else None);a.scatter([v],[i],color=GREEN,label='Target-relative certificate' if i==0 else None);a.text(v-5,i-.12,str(v),ha='right',fontsize=8);a.text(s+5,i,str(s),va='center',fontsize=8)
 a.set_yticks([0,1,2],['Native','Programmed','Native only']);a.set_xlim(6170,6455);a.set_ylim(2.8,-.6);a.set_xlabel('Unknown records / 177,006');label(a,'C','Finite family: retain each target');a.legend(fontsize=7,loc='lower left')
 p=next(p for p in load(E/'mechanism/JOINT_IDENTIFICATION_POLYGONS.json')['polygons'] if p['family']=='F9' and p['run']=='SRR29623272');a=ax[1,1];polygon(a,p['integer_count_hull_vertices'],p['N'],['Native local','Exact programmed']);label(a,'D','Joint inclusion: CFTR preparation 1')
 export(fig,'fig6_complete_correspondence')
 fig,ax=plt.subplots(2,2,figsize=(7.6,6.9),layout='constrained')
 rows=list(csv.DictReader((E/'external/HEK3_EVENT_BOUNDS.csv').open()));a=ax[0,0]
 for i,rid in enumerate(['SRR10285862','SRR10285861','SRR10285860']):
  for j,(op,cert,color) in enumerate([('M_fw','single_projection',GREY),('P_rc','single_projection',RED),('F5','target_relative',BLUE)]):
   r=next(r for r in rows if r['run']==rid and r['operator']==op and r['certificate']==cert and r['event']=='programmed_local');lo,hi=float(r['L_percent']),float(r['H_percent']);y=i+(j-1)*.16;a.plot([lo,hi],[y,y],color=color,lw=3,label=['M forward','P reverse','F5 target-relative'][j] if i==0 else None)
 a.set_yticks([0,1,2],['rep1 · 86,980','rep2 · 73,360','rep3 · 58,590']);a.set_ylim(3.0,-.4);a.set_xlabel('Programmed 10-nt tract (% raw N)');a.legend(fontsize=7,loc='lower right');label(a,'A','HEK3: fixed-rule external application')
 a=ax[0,1];h=load(E/'hek3_lattice/LATTICE_VERIFICATION.json');vals=[h[k] for k in ['whole_string_invariant_cases','programmed_local_invariant_cases','unedited_local_invariant_cases']];a.bar(['Exact\nstring','Programmed\ntract','Unedited\ntract'],vals,color=[GREY,BLUE,GREEN]);a.set_ylim(0,72);a.set_ylabel('Invariant cases / 66');label(a,'B','All optimal paths: HEK3')
 for i,v in enumerate(vals):a.text(i,v+1,str(v),ha='center')
 a.text(.02,.96,'33 actual reads × 2 references\nPurpose-enriched panel',transform=a.transAxes,va='top',fontsize=8)
 p=load(E/'external/HEK3_JOINT_POLYGONS.json')['polygons'][0];a=ax[1,0];polygon(a,p['vertices_integer_counts'],p['N'],['Programmed tract','Unedited tract']);a.legend(fontsize=7,loc='upper right');label(a,'C','Joint exclusion: HEK3 rep3')
 a=ax[1,1];manual=load(E/'external/PUBLISHED_HUMAN_MARGINALS_IDENTIFICATION.json')['sharp_agreement_bounds_from_marginals_only']
 for i,r in enumerate(manual):
  lo,hi=r['minimum_agreement_percent'],r['maximum_agreement_percent'];a.plot([lo,hi],[i,i],lw=6,color=[GREEN,BLUE][i]);a.text(lo,i+.16,f'{lo:.2f}–{hi:.2f}%',fontsize=8)
 a.set_yticks([0,1],['CRISPR-A','CRISPResso2']);a.set_ylim(-.35,1.6);a.set_xlim(68,102);a.set_xlabel('Possible agreement from marginals (%)');a.text(.02,.95,'383 selected human-reviewed records\nNo linked per-read labels available here\nIdentification bounds; not accuracy',transform=a.transAxes,va='top',fontsize=8);label(a,'D','Published benchmark: what is identified?')
 export(fig,'fig7_external_transport')
 print('Figures 6 and 7 exported in SVG, PDF and PNG; no manuscript PDF.')
if __name__=='__main__':main()
