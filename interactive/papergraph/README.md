# CRISPR PaperGraph: public scientific diagram

The PaperGraph connects the study question, observation contract, scientific
evidence, conditional derivations, adverse cases and remaining frontier. It is
a third, separate companion: the white evidence explorer and the KwanTube-derived
dark spatial tool are preserved. This interface is original code, independently
implemented, and does not copy the KwanTube interface.

Download `PaperGraph_CRISPR.html` and open it locally. This single HTML file
works offline, with no account, network request, CDN, server or installation.
GitHub's source viewer does not execute HTML.

## Read the scientific argument

The **Reading guide** follows twelve numbered stations. Its dotted links mean
"read next"; they are not proof or causal arrows. **Scientific argument** shows
the selected station's typed dependency closure, including every required input
to its named AND conjunctions. Relations and nodes can be inspected by keyboard
or touch. Refuted inferences are identified explicitly in text, not by colour
alone. Student, Technical and Audit depths explain the same objects and values
at increasing detail; navigation does not change scientific observations.

The public model retains 12 stations, 49 facts, 92 nodes, 86 declared links
(75 scientific relations and 11 separate reading-order links), and 7 named AND
groups. The finite-record hull is conditional on the event, original recorded
universe, accepted-label correctness and unrestricted combinatorial completion
of unknowns. It is not a confidence interval or a measured molecular error rate.
Read records are not independent cells, donors or molecules. EMX1 remains an
interpretive control; neither the graph nor its checks create biological data.

## Formal text remains withheld

This public release distributes graphical summaries, their authored model,
interface source and derived figures. It excludes the complete article and
formal supplement, full-text source snapshots, literal private prose anchors,
the full private subsection-routing map and private integration instructions.
Two source records retain only the byte/hash identities of the withheld
pre-MP04 baseline. Those identities are not downloadable files, current article
hashes, or a claim that the public checker independently verified private text.

Eleven numerical/model/geometry source records point to exact files already
in this repository. `references.json` contains bibliographic metadata only.
`PUBLIC_PROJECTION.json` records the public projection and the unchanged
scientific-content fingerprints, excluding provenance locators. It does not
publish the private input. The generated HTML model exports and argument-atlas
SVG metadata contain the same sanitized public model.

## Editable figures

The three figure sets are in [`../../figures/papergraph/`](../../figures/papergraph/):

- `PaperGraph_Article.svg`: portrait reading figure, with source identifiers,
  evidence badges and an explicit conditional-hull premise/conclusion inset.
- `PaperGraph_Overview.svg`: landscape reading guide.
- `PaperGraph_Argument.svg`: complete zoomable atlas containing all actual
  graph nodes and declared links. Role columns are layout, not causal ranking.

PNG previews and vector PDFs are figure exports; no article or supplement PDF
is included. At 235 mm printed height, the portrait figure's smallest essential
32 px text is approximately 8.0 pt. The full argument atlas is intended for
zoomed inspection, not reduction to a single A4 page. SVG and HTML exports are
editable, and retain declared scientific roles and limitations.

## Rebuild and check

From the repository root, use Python 3.11 or later and install the pinned
rendering dependency:

```sh
python3 -m pip install -r interactive/papergraph/requirements-figures.txt
python3 -B interactive/papergraph/code/build_papergraph.py
python3 -B interactive/papergraph/code/check_papergraph.py
```

The builder uses the DejaVu Sans regular and bold font files from
`/usr/share/fonts/truetype/dejavu/`, supplied by the Linux `fonts-dejavu-core`
package. It validates public source files in place and never copies or reads
the withheld full text. To also regenerate PNG/PDF figures, install Inkscape and
run:

```sh
python3 -B interactive/papergraph/code/build_papergraph.py --render
```

`papergraph.json` is the canonical public authored model. The template and
`code/argument_runtime.js` are editable interface sources. The checker verifies
public source fingerprints, withheld-source separation, graph/reference
integrity, the exact embedded model, atlas node/edge identities and a temporary
byte-identical HTML/model/SVG rebuild. PDF metadata may vary between renders;
cross-environment PDF-byte identity is not claimed. The numerical raw-read replay
remains the separate root command `python3 code/run_review_checks.py`.

## Antecedents and evidence limits

Claim/evidence/argument graphs already have direct precedents, including Clark,
Ciccarese and Goble's *Micropublications* (2014),
[doi:10.1186/2041-1480-5-28](https://doi.org/10.1186/2041-1480-5-28).
The contribution here is this source-bound CRISPR implementation and its linked
reading interface. It is a proposed communication artifact, not a new accepted
standard or a demonstrated improvement in reader comprehension. No PROV, RDF,
RO-Crate or complete WCAG conformance is asserted.

Local browser checks and repository CI test rendered behaviour and source
integrity. They do not establish independent biological replication, comparative
molecular accuracy, haplotype-specific function, prospective transport, or
pedagogical efficacy. Scientific claims, source publications and observations
retain their own attribution and terms. MIT is scoped to the enumerated
original interface/build/check code in `metadata/LICENSE_SCOPE.json`; it does not
relicense scientific content or the separate GPL-3.0-only supergraph.
