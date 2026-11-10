# MP27: completed computational campaign

The complete frozen campaign processed **61,056 records** at **818 loci** in
541.89 seconds of campaign time (queue and image setup excluded):

| Panel | Loci | Records | Jurisdiction |
|---|---:|---:|---|
| Development | 209 | 33,440 | Known MP24 simulation, regression only |
| Real HG002 | 489 | 16,816 | Previously inspected read-locus records, conditional genotype comparison |
| Held-out | 120 | 10,800 | New loci and a separate controlled simulation generator |

Both padding endpoints (5 and 15) and all seven frozen error budgets were
completed. The independent terminal monitor re-counted 854,784 decision rows,
replayed the held-out sources/edits, and verified full remote custody of the
7.62 MB output archive. It does not independently reimplement the full DP for
all records. The separate enumerator exhaustively checked 24,817 histories on
1,444 small instances before the campaign; 30 known adverse records passed
regression, and 206 finite-class obstruction inequalities were checked.

## Gates: retain the failure

- **In-domain computational soundness: PASS.** No known source label or allele
  was excluded when the literal generating history was within the declared
  error budget. This is conditional model validation, not a biological guarantee.
- **Observed-tract nonvacuity at budget 2, clean centered held-out reads: FAIL.**
  All 720 such records are U. The frozen failure is not replaced by a different
  endpoint or a favorable budget.
- **Independent biological accuracy: NOT ESTABLISHED.**

The failure agrees with the [literal-flank obstruction](ORIGIN_IDENTIFIABILITY_OBSTRUCTION.md):
deletion plus identical reinsertion can preserve the read while removing the
original source-flank base. Broadening the model repaired the known missing-flank
overcertifications but cannot make original-base ancestry identifiable from a
sequence that admits indistinguishable histories.

## A distinct endpoint retains information

Within the same 720 clean centered held-out records, the **source-allele** endpoint
at budget 2 resolves 594 correctly and leaves 126 unresolved. This is 82.5% of
that controlled subset, not biological accuracy. At budget 1, the observed-tract
endpoint resolves 717 correctly and leaves three U; at budget 0 it resolves all
720. No operational threshold is selected from these results.

For all 10,800 held-out records at budget 2, 6,703 have literal error counts within
the budget. Among that in-domain subset, observed-tract certificates resolve 966
correctly; source-allele certificates resolve 4,230 correctly; all other decisions
for each endpoint are U and none is discordant. Across the full panel, including
4,097 out-of-budget histories, there are 11 incompatible observed-tract
certifications and eight wrong source-allele certifications. Underestimating the
error budget can therefore produce incorrect determinate calls, not merely
detectable model mismatch. The literal error count is known only in simulation.

For real HG002 at budget 2, source-allele calls give 7,826 concordances, 29
discordances and 7,098 abstentions among 14,953 conditionally labeled records;
another 1,863 records lack that external label. This endpoint differs from the
earlier exact-observed-tract comparison: their counts must not be treated as an
accuracy improvement over the historical algorithm. Only canonical reference and
alternate templates are supplied; other real haplotype variation and unknown
noise remain model limitations. Read-locus/mate observations are not independent
molecules. See [all budgets](MP27_ALL_BUDGETS.md).

## Scientific disposition

The complete v1 supplies a verified finite-model certificate and an explicit
identifiability counterexample. The strong original-flank claim fails its utility
gate under the broad noise model. The conditional source-allele endpoint remains
nontrivial and is a candidate for further model-coverage/calibration work. No
probability, universal accuracy, state-of-the-art superiority or successful
biological validation is inferred. Original P1 manuscripts and repository were
not modified. Future versions must retain these results and separate development
from new confirmation data.
