# Bounded generation model v1

## Contract

Inputs: a nonempty observed ACGTN string, a nonempty ACGT target string, a finite
catalogue of ACGT source templates, and an absolute integer edit budget B. Each
template has an interval [lo,hi) with existing immediate source flanks lo-1 and
hi. Templates can have different intervals and lengths for different alleles.

A history selects a positive-length source fragment [a,b), then consumes it by
match/substitution (one source and observed base), deletion (one source base),
or insertion (one observed base). Insertions may occur at fragment boundaries.
Substitution, insertion and deletion each cost one. Both I-D and D-I transitions
are allowed. Source cropping is free; it is an explicit choice, not a gap score.
N is an observed unknown symbol with substitution cost one relative to ACGT.

An inserted base at source boundary i belongs immediately before source base i.
The observed tract comprises emitted bases whose source coordinate is in
[lo,hi). The flank test requires an emitted non-insertion ACGT base from each
of lo-1 and hi. A missing flank, an N inside the tract, or an empty tract gives
U. Otherwise exact equality with the target gives E and inequality gives NE.
Outside-target errors need not invalidate the tract. This observed-base event
differs from source allele identity, which is reported separately.

## Algorithm

The min-plus dynamic program carries source index i, read index j and a small
event automaton: unavailable; exact target prefix; mismatch; completed E; or
completed NE. A source-flank emission enters/completes the prefix automaton;
deletion of a flank makes the observation unavailable. Tract emissions compare
against the target; gaps and crops are permitted by the generation contract.
Starting histories are introduced with j initial insertions followed by the
first source-consuming operation. Only active histories can terminate, so an
empty source fragment is excluded. Termination is allowed at every source
position after consuming the full read. Endpoint insertions are included.

Each graph edge strictly advances i+j, so the history graph is acyclic. Minimum
cost per state suffices: two histories at the same state have identical possible
future transitions and event outcomes, and the lower-cost one dominates for all
budgets. Induction in graph order therefore gives the minimum error cost of every
possible event class. The automaton records precisely the information used by
the endpoint. Taking labels whose minimum cost is at most B gives the complete
label set of the declared history domain. Across templates the sets are united.

For a candidate set C, a determinate certificate is emitted iff C is nonempty
and all its labels equal E or all equal NE. If the real generating history is in
C, the certificate equals its endpoint by inclusion and unanimity. This is the
existing target-invariance principle applied to a different domain; mathematical
priority is not claimed. Empty C, possible U and event disagreement are distinct
reasons for no determinate certificate. No inference that C contains reality is
made from the fact that C is nonempty.

Space: O(read_length * target_length), with two source rows. Time:
O(source_length * read_length * target_length) per template. Costs above the
largest requested budget are pruned. All requested budgets reuse one DP.

## Verification jurisdiction

The independent oracle enumerates source starts and all literal M/S/I/D paths,
reconstructs origins, and computes labels without using the production states.
Small instances are legitimate exhaustive software verification, not a surrogate
for the full scientific panels. Exhaustive tests do not constitute a formal
machine-checked proof, nor validation of the model's biological assumptions.

## Quantities that must remain distinct

1. All-model observed-tract label set.
2. Source-allele candidate set within the supplied catalogue.
3. Genotype of the biological sample.
4. Population editing frequency and independent molecular accuracy.

Only the first two are computed here. Arbitrarily large error budgets or broad
catalogues can make these sets uninformative. Utility must be assessed jointly
with abstention; no posterior is implied by a Phred score or by a budget.
