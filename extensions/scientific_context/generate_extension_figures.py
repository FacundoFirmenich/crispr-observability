"""Scientific figures from exact reference and raw quality audits; no new data."""
from pathlib import Path
import json,csv,math,argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--results',type=Path,default=ROOT/'results');ap.add_argument('--output',type=Path,default=ROOT/'figures');ap.add_argument('--molecular',type=Path);args=ap.parse_args();DATA=args.results;MOL=args.molecular or DATA;OUT=args.output;OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','pdf.fonttype':42})
audit=json.loads((MOL/'MOLECULAR_REFERENCE_AUDIT.json').read_text())
with (MOL/'CANONICAL_AND_PROGRAMMED_CODONS.csv').open()as f:codons=list(csv.DictReader(f))
fig=plt.figure(figsize=(12,7.3));gs=fig.add_gridspec(2,1,height_ratios=[1.1,1],hspace=.65)
ax=fig.add_subplot(gs[0]);ax.set_xlim(-3,14);ax.set_ylim(-.8,3);ax.axis('off')
ax.text(-3,2.8,'A  Same local peptide; five nucleotide differences',fontsize=12,weight='bold')
labels=[('Residue',2),('Peptide',1.45),('Canonical DNA',.85),('Programmed DNA',.25)]
for text,y in labels:ax.text(-3,y,text,va='center',fontsize=9)
for i,row in enumerate(codons):
 color='#006d77'if row['synonymous_difference']=='True'else'#273746'
 if row['canonical_residue_1']=='508':
  ax.add_patch(Rectangle((i-.43,-.05),.86,2.35,facecolor='#f8edf0',edgecolor='#9a3651',linewidth=1))
 for text,y in [(row['canonical_residue_1'],2),(row['amino_acid'],1.45),(row['canonical_codon'],.85),(row['programmed_codon'],.25)]:
  ax.text(i,y,text,ha='center',va='center',color=color,fontfamily='DejaVu Sans Mono',fontsize=10)
ax.text(-3,-.55,'CDS c.1498-1539 | c.1521_1523del spans codons 507/508; SE2 changes are synonymous at 510-512.',fontsize=9)
ax=fig.add_subplot(gs[1]);records=audit['reference_set_local_peptide_counterexample']
for i,r in enumerate(records):
 n=r['sequence_length'];s=r['native_peptide_start_1'][0];ax.barh(i,n,height=.5,color='#dbe5e8',edgecolor='#778899')
 ax.barh(i,14,left=s-1,height=.5,color='#006d77')
 ax.text(n+18,i,f"{n} aa; local peptide {s}-{s+13}",va='center',fontsize=9)
ax.set_yticks(range(len(records)),[r['uniprot_accession']for r in records]);ax.invert_yaxis();ax.set_xlim(0,1900);ax.set_xlabel('Position in each published reference sequence (1-based)');ax.set_title('B  Equal local peptide does not identify the complete reference sequence',loc='left',fontsize=12,weight='bold',pad=15);ax.spines[['top','right']].set_visible(False)
fig.suptitle('CRISPR P1 · supplementary molecular reference audit',x=.08,y=.99,ha='left',fontsize=14,weight='bold')
fig.text(.08,.02,'Reference context only: these isoforms are not assigned to the study samples, and no expression or channel function is inferred.',fontsize=9)
fig.subplots_adjust(left=.15,right=.96,top=.9,bottom=.12)
for ext in ['svg','pdf','png']:fig.savefig(OUT/f'Figure_S1_Molecular_Reference.{ext}',dpi=180,metadata={'Date':None}if ext=='svg'else({'CreationDate':None}if ext=='pdf'else None))
plt.close(fig)

qc=json.loads((DATA/'raw_qc/RAW_QC_SUMMARY.json').read_text())
with (DATA/'raw_qc/RAW_QC_BY_POSITION.csv').open()as f:positions=list(csv.DictReader(f))
rows=qc['sample_rows'];model=json.loads((ROOT/'sources/SAMPLE_METADATA.json').read_text());meta={r['run']:r for r in model['records']}
fig,(a,b)=plt.subplots(1,2,figsize=(12.8,6),gridspec_kw={'width_ratios':[1.1,1]},layout='constrained')
colors={'CFB-26-16':'#006d77','CFB-26-02':'#b36620','CFB-06-04':'#7059a0'};styles={'untreated':':','mock':'--','PE treated':'-'}
labels=[]
for i,r in enumerate(rows):
 labels.append(r['run']+' / '+meta[r['run']]['condition']);h=r['read_length_histogram']
 for length,count in h.items():a.scatter(int(length),i,s=5+330*count/r['records'],color=colors[meta[r['run']]['preparation']],alpha=.65,edgecolors='none')
 p=[x for x in positions if x['run']==r['run']]
 b.plot([int(x['position_1'])for x in p],[float(x['mean_recorded_Phred'])for x in p],color=colors[meta[r['run']]['preparation']],linestyle=styles[meta[r['run']]['condition']],linewidth=1.3,label=r['run'])
a.set_yticks(range(len(rows)),labels);a.invert_yaxis();a.set_xlim(20,315);a.set_xticks([35,100,150,200,224,250,300]);a.set_xlabel('Original recorded read length (nt)');a.set_title('A  Length distribution in each raw stream\nPoint area scales with within-stream frequency',loc='left',fontsize=11);a.axvline(224,color='#555',linestyle='--',linewidth=.8)
b.set_xlabel('Position in the original read (1-based)');b.set_ylabel('Mean recorded base Phred score');b.set_title('B  Position-specific quality summaries\nColour: preparation; line: experimental condition',loc='left',fontsize=11);b.spines[['top','right']].set_visible(False)
from matplotlib.lines import Line2D
handles=[Line2D([0],[0],color=c,label=p)for p,c in colors.items()]+[Line2D([0],[0],color='#444',linestyle=s,label=c)for c,s in styles.items()]
b.legend(handles=handles,fontsize=8,loc='lower left',frameon=False,ncol=2)
fig.suptitle('CRISPR P1 · all nine original FASTQ streams, no trimming',fontsize=14,weight='bold')
fig.supxlabel('177,006 records retained; original whole-record meanQ30 gate: 175,010 pass / 1,996 fail. Scores do not calibrate molecular truth.',fontsize=9)
for ext in ['svg','pdf','png']:fig.savefig(OUT/f'Figure_S2_Raw_Quality_Context.{ext}',dpi=180,metadata={'Date':None}if ext=='svg'else({'CreationDate':None}if ext=='pdf'else None))
plt.close(fig)
print(json.dumps({'status':'RENDERED','figures':2,'formats':['svg','pdf','png'],'matplotlib':matplotlib.__version__}))
