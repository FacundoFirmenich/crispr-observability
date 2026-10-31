# CRISPR observability: computational companion

Version **0.4.1 / MP04** (31 October 2026). Public code, numerical data and editable figures for
a retrospective genome-editing readout audit. The complete manuscript and
supplement are held privately pending the author's chosen preprint deposition.

## Reproduce the computational checks

With Python 3.11 or later, from this directory:

```sh
python3 code/run_review_checks.py
```

The standard-library runner verifies exact derived-table arithmetic and
reconstructs observers on nine authentic public CFTR FASTQ streams using
eighteen frozen alignment tables. It checks 177,006 original records, 170,725
common-string witnesses and 6,281 unknowns against preserved result objects.
It does not calibrate molecular classification accuracy or establish independent
biological replication. Source experiments remain attributable to Sousa et al.

## Public assets

- [PaperGraph](interactive/papergraph/PaperGraph_CRISPR.html): a third, standalone
  reading and scientific-argument tool. Twelve stations distinguish definitions,
  source facts, formal conditions, adverse cases and consequences. Student,
  Technical and Audit depths preserve the same quantities and claim boundaries.
  The public model contains scientific summaries and numerical-source bindings;
  the complete article and supplement, their snapshots and quoted prose are withheld.
- `interactive/supergraph/crispr_supergraph.html`: separate dark spatial graph,
  adapted from the actual published *Beyond Redfield* KwanTube companion.
  Its side filters, navigable SVG canvas and inspector retain that interface
  while presenting CRISPR-specific claims, evidence, falsifiers and boundaries.
- `interactive/epistemic_graph.html`: existing self-contained white evidence explorer,
  preserved byte for byte.
  Download and open locally; GitHub's source viewer does not execute HTML.
  Four views expose conjunctive support paths, event partitions, conditional
  hulls, sensitivity profiles, source fingerprints and CSV/JSON/SVG exports.
- `figures/`: editable SVG/PDF/PNG exports, plotting code and exact input tables.
- `tables/`: nineteen authentic numerical CSVs, retained byte for byte.
- `raw_replay/`: nine genuine FASTQ inputs, eighteen frozen alignment tables,
  original comparators and portable replay code.
- `engine_replay/`: portable source and configurations for the historical
  144-job alignment sensitivity experiment, with clear preparation/execution scopes.
- `metadata/`: source identities, provenance, runtime records and release status.

`docs/REPRODUCTION_SCOPE.md` and `docs/SCIENTIFIC_FRONTIER.md` state what each
operation verifies. Conditional finite-read intervals, source software percentages,
sequence correspondence and biological truth are distinct quantities.

## Figures

Install `figures/requirements-figures.txt` and run:

```sh
python3 figures/code/generate_figures.py
```

The plotting source uses an explicit standalone coordinate contract and frozen
numerical tables. It does not require manuscript prose.

PaperGraph has separate portrait, landscape and complete argument-atlas figures
in [figures/papergraph/](figures/papergraph/). Editable SVGs, rendered PNG/PDF
figures, source and rebuild/check instructions are documented in
[interactive/papergraph/README.md](interactive/papergraph/README.md). Figure PDFs
are illustration assets; the article and formal supplement are not distributed.

## Interactive companions

The PaperGraph is a separate third companion. Its dotted reading route is a
pedagogical sequence, while typed scientific links expose support, qualification,
refuted inference and unestablished claims. Seven AND groups preserve required
premises. The model contains 92 nodes and 86 links: 75 scientific relations and
11 separately labelled reading-order links. Navigation does not create evidence
or establish causality. All three tools work offline without a CDN.

Argument graphs have direct prior art, including
[Micropublications (Clark, Ciccarese and Goble, 2014)](https://doi.org/10.1186/2041-1480-5-28).
This project contributes a source-bound CRISPR reading interface, without claiming
a new scientific standard or measured educational benefit.

The dark supergraph is a separate tool. Its corresponding editable source,
model builder, provenance and GPLv3 notice are in
[interactive/supergraph/README.md](interactive/supergraph/README.md). It does not
replace or restyle the white explorer. The original two tools retain their own source models and licenses.

The preserved white explorer contains 21 nodes and eight typed AND hyperedges. Its source model
keeps accepted events, non-events and unknown records separate under each declared
event, observer and reference. Completion controls choose integer hypothetical
unknown labels without changing the data. See [interactive/README.md](interactive/README.md)
for offline use and standard-library rebuild/check commands.

The exact source inventory is `metadata/REPOSITORY_MANIFEST.json`. Run
`python3 -B code/check_public_distribution.py` to check its fingerprints and text
withholding boundary. Generated runtime outputs and browser screenshots are
excluded from the public source inventory.

## Attribution and license

Current personal and institutional acknowledgments are in [ACKNOWLEDGMENTS.md](ACKNOWLEDGMENTS.md).
The verified author identity is in `metadata/AUTHOR_PROFILE_VERIFIED.json`.
MIT applies only to the original-code allowlist in `metadata/LICENSE_SCOPE.json`.
The separate KwanTube-derived supergraph code is distributed under GPL-3.0-only,
with the complete license and modification notice in `interactive/supergraph/`.
The independently implemented PaperGraph interface code is also in the exact
MIT allowlist; generated HTML includes scientific content whose source terms
remain separate. Scientific observations, data, figures and source publications retain their
distinct attribution and applicable terms; neither code license relicenses them.
See `raw_replay/DATA_NOTICE.md`.

No preprint DOI is assigned by this repository. The private manuscript and
supplement will be linked only after the author elects to deposit them.

## Added reference and recording-context checks (31 October 2026)

This repository accompanies the CRISPR P1 article; all diagrams and interactive
views are outputs of that paper. The supplementary extension in
[extensions/scientific_context/](extensions/scientific_context/README.md) adds
versioned molecular-reference checks, an independent raw quality/length audit
and accessible reference-sequence carriers. It preserves every original count,
input and earlier interactive artifact. Two new figure sets are in
`figures/scientific_context/`. Run
`python3 -B extensions/scientific_context/check_extension.py` offline.
These additions neither identify sample isoforms nor establish function,
molecular accuracy or educational efficacy. Full formal texts remain private.
