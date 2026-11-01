"""Post hoc mechanism localization; not a preregistered new biological outcome."""
from pathlib import Path
from collections import Counter
import json
ROOT=Path(__file__).resolve().parents[1]
def main():
    examples=[];unmodified=[]
    for case,chosen in [('nhej',[58]),('allele_specific',[10,18])]:
        cat=json.loads((ROOT/'results'/f'{case}_CATALOG.json').read_text())['rows']
        info=json.loads((ROOT/'inputs'/f'{case}_CRISPResso2_info.json').read_text())
        outside_weight=0;unmod_nonexact=0;boundary_or_inside=0
        for row in cat:
            f=row['features'];refname=f['source_reference']
            if f['source_status']!='UNMODIFIED' or f['exact_declared_template'] is not False:continue
            mask=info['results']['refs'][refname]['include_idxs']['value'];allowed=set(mask)
            s,r=f['full_alignment'][:2];p=-1;diff=[]
            for x,y in zip(s,r):
                if y!='-':p+=1
                if x!=y:diff.append({'reference_position_zero_based':p,'read_symbol':x,'reference_symbol':y,'inside_declared_mask_by_left_coordinate':p in allowed})
            unmod_nonexact+=row['weight']
            outside=all(not d['inside_declared_mask_by_left_coordinate'] for d in diff)
            if outside:outside_weight+=row['weight']
            else:boundary_or_inside+=row['weight']
            if row['locator']['line_including_header'] in chosen:
                # Chosen deletion witnesses have no insertion-boundary convention.
                assert all(d['read_symbol']=='-' and d['reference_symbol']!='-' for d in diff)
                assert outside
                examples.append({'locator':row['locator'],'assignment_weight':row['weight'],
                  'source_status':f['source_status'],'exact_declared_template':False,
                  'declared_included_positions_zero_based':mask,'differences':diff,
                  'explanation':'Deletion outside the source quantification mask; source UNMODIFIED does not state full-template equality.'})
        unmodified.append({'case':case,'unmodified_nonexact_unambiguous_weight':unmod_nonexact,
          'all_differences_outside_mask_under_left_coordinate_weight':outside_weight,
          'inside_or_boundary_candidate_weight':boundary_or_inside,
          'scope':'Post hoc descriptive coordinate audit of supplied alignments; insertion-boundary conventions are not adjudicated as source software errors.'})
    cat=json.loads((ROOT/'results/base_editor_CATALOG.json').read_text())['rows'];overhangs=Counter()
    for row in cat:
        ref=row['features']['full_alignment'][1]
        leading=len(ref)-len(ref.lstrip('-'));trailing=len(ref)-len(ref.rstrip('-'))
        overhangs[(leading,trailing)]+=row['weight']
    out={'timing':'POST_HOC_MECHANISM_EXPLANATION_AFTER_FIXED_EXTERNAL_ATLAS',
      'deletion_witnesses':examples,'unmodified_diagnostics':unmodified,
      'base_editor_terminal_reference_gap_weight':[{'leading':a,'trailing':b,'weight':n} for (a,b),n in sorted(overhangs.items())],
      'limits':'Alignment gaps include terminal overhangs and are not interchangeable with source quantified biological indels. The software does not authenticate the biochemical origin of any sequence.'}
    (ROOT/'results/MECHANISM_DIAGNOSTICS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'witnesses':len(examples),'diagnostics':unmodified,'base_editor_overhang_classes':len(overhangs)},indent=2))
if __name__=='__main__':main()
