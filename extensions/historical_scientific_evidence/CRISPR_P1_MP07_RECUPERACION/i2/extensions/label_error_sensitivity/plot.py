"""Plot exact conditional sensitivity; no simulated measurements."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import contrast, SOURCE


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    with SOURCE.open() as f:
        rows = [r for r in csv.DictReader(f) if (r['family'],r['certificate'],r['event']) == ('F9','target_relative','native')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.hashsalt':'crispr-i2-label-errors'})
    fig, ax = plt.subplots(figsize=(8.3,4.7),layout='constrained')
    colors = ('#176B91','#B85B13','#388351')
    names = ('CFB-26-16','CFB-26-02','CFB-06-04')
    npoints = {}
    for prep, name, color in zip(('1','2','3'),names,colors):
        arms = {r['condition']:{**r,**{k:int(r[k]) for k in ('E','NE','U')}} for r in rows if r['preparation']==prep}
        t,c = arms['PE treated'], arms['mock']
        cuts={F(0),F(1,4)}
        for cap in (t['E'],c['NE']):
            cuts.update(F(k,cap) for k in range(cap//4+1))
        rates=sorted(cuts); lows=[contrast(t,c,r)[0] for r in rates]
        ax.step([100*float(r) for r in rates],[100*float(v) for v in lows],where='post',color=color,label=name,lw=1.6)
        npoints[prep]=len(rates)
    ax.axhline(0,color='#444444',lw=.9)
    ax.set(xlim=(0,25),xlabel='Hypothetical cap within each accepted class (%)',
           ylabel='PE minus mock lower endpoint (percentage points)',
           title='Accepted-label sensitivity of the native-local certificate')
    ax.text(.02,.035,'F9 target-relative domain; raw denominator retained\nCaps are assumptions, not measured error rates',transform=ax.transAxes,fontsize=9,
            bbox={'facecolor':'white','alpha':.9,'edgecolor':'none'})
    ax.legend(loc='upper right',frameon=False)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.18)
    for fmt in ('svg','png','pdf'):
        metadata={'Date':None} if fmt=='svg' else {'CreationDate':None,'ModDate':None} if fmt=='pdf' else {}
        fig.savefig(out/f'figS4_label_error_bounds.{fmt}',dpi=180,metadata=metadata)
    plt.close(fig)
    (out/'FIGURE_SCOPE.json').write_text(json.dumps({'kind':'exact_conditional_function_of_actual_counts','model':'same hypothetical cap r within accepted E and NE classes in each arm','target':'native','family':'F9','certificate':'target_relative','evaluated_breakpoints_by_preparation':npoints,'new_measurements':False,'rates_calibrated':False},indent=2)+'\n')


if __name__=='__main__':
    main()
