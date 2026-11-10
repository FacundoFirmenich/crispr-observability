# MP28: locus-block calibration, catalogue coverage and noise shift

The frozen Kaggle campaign completed in **164.66 seconds** of campaign
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

| Comparison | Nominal | q | Covered blocks /256 | Correct singleton /3072 | Wrong singleton | Both labels | Empty |
|---|---:|---:|---:|---:|---:|---:|---:|
| Matched baseline | 90% | 4 | 246 (96.09%) | 512 | 0 | 2550 | 10 |
| Matched baseline | 95% | 5 | 255 (99.61%) | 312 | 0 | 2759 | 1 |
| Noise shift, unchanged calibration | 90% | 4 | 1 (0.39%) | 614 | 3 | 1091 | 1364 |
| Noise shift, unchanged calibration | 95% | 5 | 3 (1.17%) | 687 | 3 | 1515 | 867 |
| Missing background, unchanged calibration | 90% | 4 | 0 (0.00%) | 201 | 8 | 1185 | 1678 |
| Missing background, unchanged calibration | 95% | 5 | 1 (0.39%) | 682 | 11 | 1485 | 894 |
| Missing background, recalibrated | 90% | 8 | 233 (91.02%) | 622 | 0 | 2427 | 23 |
| Missing background, recalibrated | 95% | 9 | 254 (99.22%) | 393 | 0 | 2677 | 2 |
| Expanded catalogue, fixed threshold | 90% | 4 | 245 (95.70%) | 504 | 0 | 2557 | 11 |
| Expanded catalogue, fixed threshold | 95% | 5 | 254 (99.22%) | 315 | 0 | 2755 | 2 |
| Expanded catalogue, recalibrated | 90% | 4 | 245 (95.70%) | 504 | 0 | 2557 | 11 |
| Expanded catalogue, recalibrated | 95% | 4 | 245 (95.70%) | 504 | 0 | 2557 | 11 |


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
