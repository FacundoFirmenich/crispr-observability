# Literal-flank identifiability obstruction in generation model v1

This deduction was recorded while MP27 was running, before inspecting held-out
outcomes. It does not alter the frozen model, budgets or nonvacuity gate.

**Proposition.** For any admissible history with a determinate observed-tract
label and edit cost c, the same template and the same observed read admit a
history with missing original-flank observation and edit cost at most c+2.

**Construction.** A determinate history emits an ACGT base by consuming the
source coordinate lo-1. Replace that match/substitution by insertion of exactly
the emitted observed base at boundary lo-1 followed by deletion of source base
lo-1. The rest of the history is unchanged. Insertions at that boundary and I-D
transitions are admitted by model v1. The read is identical. The source fragment
is still positive length and consumes the same source coordinates. However, the
base at lo-1 has now been deleted rather than observed, so the flank endpoint is
U. The original step cost 0 or 1; its replacement costs 2. Thus the added cost is
at most 2. This is a model-specific identifiability limit, not a software defect
or a claim of mathematical priority.

Consequently c_U <= min(c_E,c_NE)+2 whenever the right-hand finite class exists.
A budget reaching that bound cannot yield a determinate literal-flank
certificate. In particular, clean fully observed histories have c=0, so the
budget-2 clean-center observed-tract nonvacuity gate is predicted to fail under
the deliberately broad edit-history model. The frozen criterion is retained;
it must not be retrospectively replaced by source-allele utility.

Source allele remains a different endpoint. The construction does not change
the allele, so it does not imply that all allele candidate sets are ambiguous.
Their resolution is evaluated separately, conditional on the candidate catalogue
and edit bound. A probabilistic treatment could make the substitute history
unlikely, but would supply model-dependent probability, not recover a lost
deterministic guarantee. Additional measurements can help only with justified
linkage and a declared joint error model; they are not assumed independent.
