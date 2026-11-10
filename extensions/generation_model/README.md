# Generation-model and calibrated source-allele extension

This is the author-authorized reintegration of MP27 and MP28 into CRISPR P1.
The prior releases and scientific observations remain unchanged. Original code
in this directory is MIT; the principal repository's other license scopes remain.
Manuscripts and the formal supplement are not distributed.

- [MP27 findings and adverse gate](docs/RESULTS_MP27.md).
- [MP28 calibration and all shift comparisons](docs/RESULTS_MP28.md).
- [Scientific reintegration boundaries](docs/REINTEGRATION_MP27_MP28.md).
- [Evidence graph](results/EPISTEMIC_GRAPH_MP28.json).
- [Full result archives](https://github.com/FacundoFirmenich/crispr-observability/releases/tag/v0.8.0).

## Reproduce finite checks

From this directory, with Python 3.11 or newer:

```sh
python -m pip install .
python tests/verify_oracle.py
python tests/verify_regressions.py
python tests/verify_obstruction.py
python tests/verify_calibration.py
python analysis/report_mp28.py
```

The frozen MP28 Kaggle source and input bindings are in `campaign/mp28/kaggle`.
Its outputs contain the calibration/test selection, every generated sequence,
literal edit ledger, source-allele scores, per-locus resumable checkpoints,
calibration freeze and manifest. MP27 outputs are retained separately.

The MP28 guarantee is conditional on exchangeable locus blocks and a matched
regime; fixed chr20 site disjointness does not establish biological calibration.
At nominal 95%, baseline whole-block coverage was 255/256 with 312/3072 correct
singletons. Noise shift and omitted haplotype backgrounds reduced coverage to
3/256 and 1/256. Matched recalibration and an oracle expanded catalogue repair
coverage in this controlled panel, with substantial ambiguity remaining.

Literal observed-base ancestry remains unidentifiable under the permissive
history model at budget 2 on clean centered reads. Neither source-allele sets nor
high set coverage negate that MP27 result. No new wet-lab data or universal
accuracy claim is supplied by these campaigns.
