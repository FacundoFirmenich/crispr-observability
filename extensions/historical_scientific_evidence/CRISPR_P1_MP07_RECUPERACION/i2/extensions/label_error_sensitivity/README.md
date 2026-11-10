# Accepted-label sensitivity, MP06 integrity successor I2

This extension addresses a conditional premise of the finite-record bounds:
accepted event/non-event labels need not all be correct. It does not calibrate
errors or provide new biological observations. It preserves all frozen MP06
tables, computational models and existing graph baselines.

For one binary target and one arm, E/NE/U are frozen counts, N=E+NE+U.
If at most a accepted events are false positives and at most b accepted
non-events are false negatives, the sharp interval hull is
[(E-a)/N, (E+U+b)/N]. For treated and mock arms, subtract the corresponding
opposite endpoints. Unknowns remain unrestricted. No statistical independence
of reads or estimated error probability is assumed.

The common class cap r permits floor(r E) and floor(r NE) errors in each arm.
It is a hypothetical deterministic allowance within each accepted class.
It is not an observed misclassification rate, a confidence level or a bound on
errors before sequencing. Targets are analyzed separately; these marginal
intervals cannot be combined into a new simultaneous joint-target region.
Shared-family mechanisms remain uncalibrated.

Run outside the repository output tree:

    python -B extensions/label_error_sensitivity/analyze.py --output /tmp/crispr-i2-sensitivity

The script authenticates all 162 actual source count rows, computes 54 paired
contrasts over nine explicitly declared caps, and finds exact first nonpositive
breakpoints. It verifies sharpness against exhaustive mathematical binary
assignments for all 791 admissible small count/budget configurations with
1 <= N <= 7. These are logical checks, not synthetic biological observations.
The zero-budget endpoints equal the frozen MP06 endpoints exactly.

`expected/` contains the actual output of this execution. Reproduction can
compare its three files byte-for-byte with the command output. Figure source:

    python extensions/label_error_sensitivity/plot.py --output /tmp/crispr-i2-figures

Plotting requires Matplotlib. Curves evaluate every exact error-budget breakpoint
through r=1/4; they are mathematical sensitivity functions of actual count data,
not additional observed rates or an interpolation of biological results.

The complete article and formal supplement remain private. This public code
and numerical output are subordinate outputs of CRISPR P1. Broader finite-state,
acoustic and quantum projects supply methodological context, not additional
editing or functional evidence. No general mathematical priority is claimed.
