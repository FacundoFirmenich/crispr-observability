"""Recalculate a paired reporting comparison from authenticated historical outputs.
This replays output accounting, not a fresh execution of CRISPResso2.
"""
from pathlib import Path
from fractions import Fraction
import argparse,csv,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def compare():
 d=ROOT/'extensions/reproduction/derived'
 names=['PE_TOOL_VS_EVENT_COMPARISON.csv','PE_NATIVE_TOOL_METRICS.csv','PE_EVENT_BOUNDS.csv','PE_ANALYSIS_VERIFICATION.json']
 read=lambda n:list(csv.DictReader((d/n).open(encoding='utf-8')))
 native=read(names[1]);bounds=read(names[2]);rows=[]
 for r in read(names[0]):
  n=int(r['raw_N']);aligned=int(r['all_amplicons_aligned']);q=int(r['Q30_input']);rid=r['run']
  nr=next(x for x in native if x['run']==rid and x['Amplicon']=='Prime-edited')
  assert int(nr['Reads_in_input'])==q and int(nr['Reads_aligned_all_amplicons'])==aligned
  prime=int(nr['Reads_aligned']);assert prime==int(r['Prime_edited_assigned'])
  row={'run':rid,'condition':r['condition'],'raw_N':n,'Q30':q,'aligned':aligned,'outside_aligned_denominator':n-aligned,'prime_assigned':prime,'prime_assignment_over_aligned':str(Fraction(prime,aligned)),'prime_assignment_over_raw':str(Fraction(prime,n))}
  for event in ['native_local','exact_programmed']:
   for op in ['PE_selected_reference','old_three_operator_witness']:
    x=next(x for x in bounds if x['run']==rid and x['operator']==op and x['event']==event)
    e,ne,u=[int(x[k]) for k in ['E','NE','U']];assert e+ne+u==n
    lo,hi=Fraction(e,n),Fraction(e+u,n);assert lo==Fraction(x['L_exact']) and hi==Fraction(x['H_exact'])
    row[event+':'+op]={'E':e,'NE':ne,'U':u,'lower':str(lo),'upper':str(hi)}
  rows.append(row)
 assert len(rows)==9 and sum(x['raw_N'] for x in rows)==177006
 return {'status':'PASS','engine':'CRISPResso2 2.3.4','scope':'Historical same-record output-accounting replay; nine libraries; no fresh caller execution or independent editing truth','source_sha256':{n:hashlib.sha256((d/n).read_bytes()).hexdigest() for n in names},'rows':rows,'interpretation':'A native Prime-edited assignment and an exact-DNA or native-local event have different definitions. Denominator-normalized assignment remains an assignment, not event accuracy. No ampliCan/CRISPR-A execution is represented by these rows.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(compare(),indent=2)+'\n');print('PASS: 9 libraries, 177006 original records, exact rational partitions')
