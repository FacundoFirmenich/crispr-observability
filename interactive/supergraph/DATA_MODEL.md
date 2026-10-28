# Source-bound data for the reference scientific map

This graph model is a separate conversion of the existing public numerical
companion. The previous evidence explorer and its numerical model are untouched.
No manuscript or formal-supplement prose is read during conversion.

The six layers contain 21 existing core nodes, eight explicit AND gates, 22 frozen
source-artifact nodes, and six declared correction-test nodes. The artifact nodes
can be hidden initially for legibility; their paths and exact fingerprints remain
available for inspection. Correction tests are declared conditions for revision,
not newly executed experiments or new biological evidence.

The graph has 57 nodes and 59 typed edges. Required-premise edges terminate at an
AND gate; its outgoing edge is conditional support. Artifact edges document
provenance. Revision-condition edges state what would require revising a local
conclusion. Neither layout nor these navigation relations measures causal effect,
evidential strength or molecular accuracy.

The context selector can display two declared events and three preparations.
Its 18 common-string-witness partitions and 12 compatible two-arm contrasts are
copied as data objects from the frozen input model. The event, raw denominator,
reference contract and exact rational endpoints remain attached. No value from a
different software-defined event is relabeled as this event.

Build and verify from a repository checkout, adjusting script and output paths to
their installed location:

```sh
python3 -B build_reference_graph_model.py \
  --model /path/to/repository/interactive/evidence_model.json \
  --repository /path/to/repository \
  --output reference_graph_model.json
python3 -B verify_reference_graph_model.py \
  --model /path/to/repository/interactive/evidence_model.json \
  --repository /path/to/repository \
  --graph reference_graph_model.json \
  --output MODEL_CONVERSION_CHECK.json
```

The builder checks the actual bytes of all 22 canonical source artifacts before
writing. The checker additionally requires a byte-identical temporary rebuild,
unchanged core source fields and hyperedge premises, identical witness rows,
valid exact interval endpoints, and clearly unexecuted future-test status.

These are computational source and conversion checks. Accepted-label correctness,
orthogonal molecular truth, prospective transport and independent human
replication remain separate scientific questions.
