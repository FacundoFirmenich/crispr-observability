"""MIT. Generate every declared comparison from verified MP28 result objects."""
import csv,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parents[1];d=json.loads((R/'results/MP28_RESULT.json').read_text());v=json.loads((R/'results/MP28_VERIFICATION.json').read_text());assert d['status']=='COMPLETE' and v['status']=='PASS_FULL_RECORD_REPLAY_AND_RECOUNT'
names={'baseline_matched':'Matched baseline','noise_shift':'Noise shift, unchanged calibration','catalogue_shift':'Missing background, unchanged calibration','catalogue_recalibrated':'Missing background, recalibrated','catalogue_expanded_fixed_threshold':'Expanded catalogue, fixed threshold','catalogue_expanded_recalibrated':'Expanded catalogue, recalibrated'}
rows=[]
for r in d['comparisons']:
    rows.append({'comparison':r['comparison'],'nominal_percent':100*(1-r['alpha'][0]/r['alpha'][1]),'threshold':r['threshold']['q'],'covered_blocks':r['covered_blocks'],'blocks':r['blocks'],'block_coverage_percent':100*r['covered_blocks']/r['blocks'],'included_records':r['included'],'records':r['records'],'record_coverage_percent':100*r['included']/r['records'],**{k:r.get(k,0) for k in ('correct_singleton','wrong_singleton','ambiguous','empty')},'scope':r['scope']})
with (R/'results/MP28_COMPARISONS.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
table='| Comparison | Nominal | q | Covered blocks /256 | Correct singleton /3072 | Wrong singleton | Both labels | Empty |\n|---|---:|---:|---:|---:|---:|---:|---:|\n'
for r in rows:table+=f"| {names[r['comparison']]} | {r['nominal_percent']:.0f}% | {r['threshold']} | {r['covered_blocks']} ({r['block_coverage_percent']:.2f}%) | {r['correct_singleton']} | {r['wrong_singleton']} | {r['ambiguous']} | {r['empty']} |\n"
text=f'''# MP28: locus-block calibration, catalogue coverage and noise shift

The frozen Kaggle campaign completed in **{d['seconds']:.2f} seconds** of campaign
time, excluding queue and image setup. It used 256 new calibration loci and 256
new test loci from the fixed 2024 chr20 draft benchmark domain. There were 6,144
calibration records and 9,216 test records. All were controlled simulations;
there are no additional biological observations. Six paired comparisons reuse
the relevant 3,072-record test regime and must not be summed as independent data.

The score is the minimum unit-edit distance from a read to any positive-length
source crop for an allele. A block score is the maximum true-allele score across
12 generated records at one locus. Split-conformal thresholds were frozen before
test scoring. The conventional marginal guarantee requires exchangeable locus
blocks and a fixed matched scoring/generation regime; site disjointness alone
does not prove these premises. The noise-shift and fixed-threshold diagnostics
have no nominal coverage guarantee. See [Angelopoulos and Bates](https://arxiv.org/abs/2107.07511).

{table}

## Scientific interpretation

The matched baseline at nominal 95% covered all 12 records in 255/256 loci,
but only 312/3072 records gave a correct singleton. Of 1,024 centered reads,
261 were correct singletons; the two partial-event geometries contributed only
48 and 3. High coverage is therefore achieved largely through ambiguity and
must not be presented as high classification accuracy or operational resolution.

With the same threshold q=5, noise shift left only 3/256 whole blocks covered,
and unrepresented haplotype backgrounds left 1/256. Their 3 and 11 wrong singleton
calls are retained alongside empty sets. A larger number of singleton calls
under shift is not evidence of improved reliability.

For the same nuisance test records, recalibrating the canonical catalogue gave
q=9 and 254/256 covered blocks, with 393 correct singletons. Expanding the
catalogue and recalibrating gave q=4 and 245/256 covered blocks, with 504 correct
singletons. These are different empirical coverage/resolution operating points,
not a dominance claim. Catalogue expansion at the fixed q=5 gave 254/256 blocks,
315 correct singletons and no wrong singleton in this sample. Candidate-set
inclusion under catalogue expansion was checked for every record.

The added backgrounds were deliberately constructed and then supplied to the
expanded catalogue. This is an oracle-catalogue stress experiment, not discovery
of unknown variants from reads. It localizes the effect of catalogue omission,
but leaves practical haplotype discovery and real-data calibration unresolved.
All six comparisons and both alpha values remain in the report; no threshold
was selected after test outcomes. Zero observed wrong singletons in matched
panels does not imply zero future or biological error.

## Verification and reintegration

All 15,360 generated strings were replayed from the independent source/edit
ledger; every primary count and all calibration quantiles were independently
recounted. An independent dynamic-programming implementation rescored 36 test
records. Prior exhaustive small-domain score verification covered 903 cases,
including 480 comparisons with the full generation-history engine. This is
not an independent rescoring of the entire campaign. There were no detected
true-history inclusion failures within the applicable catalogue and error cap.

The complete 2.01 MB result archive was uploaded and fully read back from the
versioned GitHub release. The user authorized reintegration with the principal
line. MP27's failed literal-flank nonvacuity gate remains unchanged: calibrated
source-allele sets do not establish the ancestry of individual observed bases.
The new contribution is a measured separation of model coverage, calibration,
set resolution and transport failure. Independent biological validation remains
open. A next confirmatory campaign should use genuinely new loci and calibrate
under the intended noise regime, while reporting resolution at comparable
coverage and keeping haplotype-discovery performance separate.
'''
(R/'docs/RESULTS_MP28.md').write_text(text,encoding='utf8',newline='\n')
g=json.loads((R/'results/EPISTEMIC_GRAPH.json').read_text());g['version']=2
for key,path in [('mp28_protocol','protocols/MP28.json'),('mp28_result','results/MP28_RESULT.json'),('mp28_verification','results/MP28_VERIFICATION.json'),('mp28_calibration','results/MP28_CALIBRATION_FREEZE.json')]:
    g['nodes'].append({'id':key,'kind':'artifact','path':path,'sha256':hashlib.sha256((R/path).read_text(encoding='utf8').encode()).hexdigest(),'hash_representation':'UTF-8 text with LF line endings'})
g['nodes'] += [{'id':'C6','kind':'claim','text':'Locus-block calibration attains useful coverage under matched controlled regimes','status':'EMPIRICALLY_SUPPORTED_WITH_LOW_SINGLETON_RESOLUTION'}, {'id':'C7','kind':'claim','text':'Unchanged calibration transports across the declared noise and catalogue shifts','status':'ADVERSE_TEST_FAILED'}, {'id':'C8','kind':'claim','text':'Oracle catalogue expansion and matched recalibration repair coverage on the nuisance panel','status':'SUPPORTED_CONTROLLED_SIMULATION_NOT_HAPLOTYPE_DISCOVERY'}]
g['edges'] += [{'from':a,'to':b,'relation':rel} for a,b,rel in [('mp28_protocol','mp28_result','freezes_design'),('mp28_calibration','mp28_result','fixes_thresholds_before_test'),('mp28_verification','mp28_result','replays_and_recounts'),('mp28_result','C6','supports_conditionally'),('mp28_result','C7','refutes_in_declared_stress_panel'),('mp28_result','C8','supports_conditionally'),('C4','C6','limits_transport'),('C8','C5','does_not_establish')]]
ids={x['id'] for x in g['nodes']};assert len(ids)==len(g['nodes']) and all(e['from'] in ids and e['to'] in ids for e in g['edges'])
(R/'results/EPISTEMIC_GRAPH_MP28.json').write_text(json.dumps(g,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps({'comparisons':len(rows),'graph_nodes':len(g['nodes']),'graph_edges':len(g['edges']),'report':'docs/RESULTS_MP28.md'}))
