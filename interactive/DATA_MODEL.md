# Numerical and epistemic data for the offline companion

`evidence_model.json` is the machine-readable model for the interactive evidence
explorer. It is built only from 19 frozen CSV tables, two numerical JSON artifacts
and the three-row coordinate table. No manuscript or supplement prose is an input.
The source files are retained under `sources/`; the model records their SHA-256
fingerprints and clean public repository paths.

Rebuild and check with the Python standard library:

```sh
python3 code/build_model.py
python3 code/check_contracts.py
```

The model contains 144 observation rows, 96 source contrast rows, 288 alignment
event/job rows and 192 derived alignment contrasts. These are explicitly different
catalogs, not independent replicates. All proportions retain their raw record
denominators. External catalog assignment weights, paired input units and raw
single-end records retain their source units.

The event definitions `native_peptide` and `exact_programmed` remain separate.
Software `MODIFIED/aligned` is an additional statistic with its own denominator;
it is not relabeled as either observed local DNA event. The common-string witness
continues to require the original target extraction and both baseline projections;
it is not retrospectively enlarged to the fourteen grid profiles.

`graph.nodes` and `graph.hyperedges` implement an inspectable evidence dependency
map. Each typed AND hyperedge has multiple required premises. Premise selection
can change which support path is displayed; it cannot alter the archived data or
turn a hypothetical toggle into a scientific refutation. Calibration and future
transport remain separate scope boundaries.

For the calculator, accepted labels must be correct for the declared event.
Unknown labels may vary independently. The hull is `[E/N,(E+U)/N]`, but attainable
finite fractions are `(E+k)/N` for **integer** `k=0..U`. A fraction slider should
therefore map to an integer `k` and display the selected count; a continuous relaxed
interpolation is not necessarily an attainable completion. Edited inputs and
completion choices are hypothetical and never new observations.

`MODEL_VERIFICATION.json` records exact source arithmetic checks.
`CONTRACT_BOUNDARY_CHECKS.json` records invalid-input rejection, zero-unknown,
unequal-denominator and finite-attainability boundary cases, and graph integrity.
Successful checks establish internal numerical consistency within these declared
contracts; they do not establish molecular accuracy or independent human replication.
