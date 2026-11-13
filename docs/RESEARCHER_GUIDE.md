# Researcher and community guide

Start with the [four-step decision guide](START_HERE.md).

[Castellano](RESEARCHER_GUIDE.es.md) · [Leading case](LEADING_CASE.md)

## Choose the scientific task

| Task | Entry point | Output |
| --- | --- | --- |
| Understand or challenge a claim | Epistemic Supergraph and PaperGraph | Traceable premises, sources and scope; exportable map |
| Inspect exact historical counts | Numerical epistemic explorer | Original finite-record counts and conditional intervals |
| Manipulate a sequence/structure or inspect an alignment | Read Observability Lab | Interactive 3D scene, coordinate report, result/replay JSON, PNG and alignment TSV where applicable |
| Run several named inputs or connect an IDE | Multi-input engine CLI/API/MCP | Versioned manifest and structured source-bound report |
| Recompute published calculations | Native launcher, Python runner, Streamlit or notebook | Scoped calculation and execution receipt |

These interfaces answer different questions. A graphical scene is not an
alignment oracle; a source checksum authenticates bytes, not a biological label.
Start with the leading case, then bring a case whose coordinate and evidence
contracts you can state explicitly.

## Open the tools

The repository landing page links the standalone tools. Download the HTML if
offline work is preferable. Use a modern WebGL 2 browser for the laboratory.
Orbit by dragging the scene; zoom with the wheel; use the camera and view controls
to restore orientation. The Supergraph supports individual cards, columns,
whole-map movement, highlighted relations and a conceptual play/pause tour.
PaperGraph provides Student, Technical and Audit explanations.

The browser laboratory presents five modes: anchor containment, competing spans,
ideal DNA, imported structures and read alignments. Change one parameter at a
time and inspect the explanation and numerical result, then save the replay.
Source-bound views and deliberate what-if changes remain distinguishable.

## Bring your inputs

| Input | Contract |
| --- | --- |
| FASTA | One sequence, 2–360 ACGTN bases; short local inspection |
| PDB / textual mmCIF / XYZ | Normalized atom coordinates; at most 30,000 atoms and 6 MB text; declared structural provenance |
| SAM | Record/CIGAR inspection with explicit reference and locus |
| BAM plus BAI/CSI | Local indexed access; BAI tested with public fixtures; CSI fixture coverage remains limited |
| Portable JSON | Versioned sequence, structure, alignment or replay schema; incompatible schemas rejected |
| Project manifest | Up to 32 named inputs and 64 ordered operations; SHA-256 checks available |

Coordinates are zero-based, half-open. Structural distances use angstroms.
The alignment window is at most 50,000 reference bases and 2,500 records.
The browser holds one active case; the engine supports multiple named inputs.
CRAM and whole-genome streaming are outside the present implementation. Read
the format-specific [laboratory guide](../interactive/lab/README.md) and
[junction contract](../engine/GUIA_EVIDENCIA_SPLICING_v0.5.md) before conversion.

Import normalizes supported formats into portable JSON. Export preserves a
replayable computational case; it is not a lossless converter for every field
of every source format. Keep the original file and hash alongside the normalized
case. Absence from a bounded input window means not observed there.

## Programmatic work and MCP

Install the pinned engine dependencies with `npm ci --prefix engine`. Use
`node engine/project-cli.mjs check <manifest>` followed by `run <manifest>`.
The same kernel is available through JavaScript and the stdio MCP server:
`node engine/project-mcp.mjs <authorized-study-directory>`. The server exposes
`lab_project_run`; paths remain inside that directory. The separate browser MCP
bridge supports normalization, boundary and span calculations, alignment reports
and rendered exports. See [MCP instructions](../interactive/lab/MCP.md).

An agent should submit an explicit manifest, retain the returned parameters and
source identities, and distinguish a requested operation from a completed result.
Never infer sample or molecule identity from a matching filename alone.

## Reproduction, Streamlit and Docker

Native v0.10.0 packages include their Python runtime and frozen inputs:
Windows x86-64 `START_HERE.bat`, Linux x86-64 `START_HERE.sh`, macOS Intel/Apple
silicon `START_HERE.command`. These are platform-specific launchers, not one
binary for every operating system. Current source updates postdate that fixed
release; consult the release commit before attributing new tools to a package.

The source dashboard has English/Spanish controls, interactive evidence,
Supergraph, laboratory, scoped replay, execution ledger and community guides:

```sh
python -m pip install -r reproducibility/requirements-reviewer.txt
python -m streamlit run reproducibility/reviewer_app.py
```

Alternatively:

```sh
docker build -f Dockerfile.reviewer -t crispr-reviewer .
docker run --rm -p 127.0.0.1:8501:8501 -v crispr-review-results:/results crispr-reviewer
```

Open `http://localhost:8501`. The named volume retains outputs. The Docker
dashboard includes browser tools and Python reproduction; the Node engine has
its separately pinned installation. Quick/counts/full replay have different
scopes; read the receipt. Choose a persistent output location rather than an
ephemeral container or a crowded local disk. No launcher starts Colab or Kaggle.

## Hardware and realistic accessibility

Local development and editorial verification used a conventional 2016 laptop
with an Intel Core i5-6200U at 2.30 GHz, as reported by the author. Larger scientific
campaigns used hosted resources; native builds use platform CI. This establishes
the development context, not a benchmark that every full workflow runs on a
ten-year-old computer. The standalone laboratory is small, but rendering still
depends on browser/WebGL support. No dedicated GPU is required by the frozen
Python replay; runtime and memory depend on the selected scope.

## Contribute to the research

Three programme lines already use or adapt the infrastructure; see their
[names and integration states](RESEARCH_LINES.md). To contribute, report a
localized discrepancy with the commit, command/manifest, source hash, coordinate
system, expected result, observed result and evidence unit. Share only data you
are permitted to publish. Use an issue for reproducible defects, a discussion
for scientific questions or proposed use, and a pull request for a tested change.
The [community roadmap](COMMUNITY_ROADMAP.md) lists concrete participation routes.

Useful contributions include a new independently documented case, an adversarial
coordinate fixture, a verified format adapter, an accessible explanation or an
external reproduction receipt. Report unsuccessful cases too. Citations,+independent use and reproduction are distinct outcomes; none is manufactured
through badges, promotional claims or issue counts.

## Publish and cite a result

Retain original inputs, code commit, manifest, complete report, evidence unit,
unresolved outcomes and renderer/replay version. Cite the exact component and
version. DOI 10.5281/zenodo.24121127 identifies v0.8.1, not later laboratory
development. Original current software is GPL-3.0-only; the independent
generation-model extension retains MIT and dependencies retain their licences.
The private manuscript and supplement are not included in public source releases.
