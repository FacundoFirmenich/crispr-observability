# CRISPR Observability

**Start with a usable answer:** [four-step guide](docs/START_HERE.md) · [authentic input JSON](examples/observation-report/input.json) · [decision report](examples/observation-report/report.html). Python-only CLI and guided Streamlit entry; no graph expertise required.

[English](README.md) · [Castellano](README.es.md) · **[Open the research environment](https://facundofirmenich.github.io/crispr-observability/)**

**Trace the claim. Inspect the evidence. Change the parameters. Replay the result.**

CRISPR Observability connects a method for auditing genome-editing readouts to an
interactive evidence environment and a reusable computational laboratory.
A product can move the evidence needed to detect it beyond a finite read.
Different rules can agree on a label while identifying different observed bases.
The method makes those failures visible and certifies targets over explicitly
declared interpretation domains.

## Three connected research capabilities

| Start here | What you can do | Scientific purpose |
| --- | --- | --- |
| **[Epistemic Supergraph](interactive/supergraph/crispr_supergraph.html)** | Follow a conceptual tour; inspect claims, required premises and sources; move cards or columns; trace highlighted connections; export the map | Understand and challenge how a conclusion follows from its evidence |
| **[Read Observability Lab](interactive/lab/CFTR_OBSERVABILITY_LAB.html)** | Rotate sequences and structures in 3D; alter boundaries and spans; inspect SAM/BAM coordinates; import cases; export results, replay and detailed renders | Test the consequences of a specified computational intervention |
| **[Reproduce the analysis](reproducibility/README.md)** | Launch native Windows, macOS or Linux packages; select quick, counts or full replay; inspect the execution ledger in Streamlit or a notebook | Recompute the published calculations at an explicit evidence scope |

Download an HTML file using GitHub's **Download raw file** control and open it in a
modern browser. The Supergraph and laboratory are standalone applications; the
GitHub source preview does not execute them. The laboratory requires WebGL 2.
Imports and calculations stay in the browser until you explicitly export a file.

### Explore the argument at the resolution you need

The **Supergraph** connects the core CFTR findings to later COL17A1, germline and
external paired-library checks. Its tour advances by conceptual block, with
play, pause, previous and next controls and a scientific explanation at each step.
Typed links and conjunction gates preserve the distinction between a premise,
an observation, a calculation, a limitation and a future revision test.

The **[numerical epistemic explorer](interactive/epistemic_graph.html)** gives four
focused views of the original CFTR counts, bounds, sensitivity and contracts.
**[PaperGraph](interactive/papergraph/PaperGraph_CRISPR.html)** provides a separate
guided argument with Student, Technical and Audit reading depths. These are
complementary views with declared scopes, not interchangeable versions of one graph.

### Bring your own scientific case

Read Observability Lab accepts FASTA, PDB, textual mmCIF, XYZ, SAM, indexed BAM
and versioned JSON within its [documented format limits](interactive/lab/README.md).
Its experiments connect coordinates to observed sequence spans and a manipulable
3D scene. Export the result JSON to preserve parameters and camera state, or a
PNG for presentation. Structural inspection uses imported coordinates or explicitly
ideal DNA geometry; it does not predict molecular dynamics or editing efficiency.

The **[multi-input engine](engine/README.md)** adds study manifests, reusable
operations, a JavaScript API, a CLI and an MCP connection for IDE workspaces.
It can compare multiple named inputs while preserving source provenance and
project namespaces. The browser currently displays one active case.
See the [capability and evidence map](docs/COMPANION_CAPABILITIES.md).

## Start with a scientific question

**[Leading case: when one missing base changes what can be observed](docs/LEADING_CASE.md)**
connects real CFTR reads, a geometric intervention, a seven-operation manifest
and exact expected outcomes. **[Researcher and community guide](docs/RESEARCHER_GUIDE.md)**
explains how to import a case, connect MCP, export evidence and contribute.

The infrastructure is used and adapted across **CRISPR Observability**, **Biogen /
BioGen S1 (local splicing observability)** and **ChronoWitness (mathematical biology
and local temporal dynamics)**. [Their integration states differ](docs/RESEARCH_LINES.md):
native laboratory, applied junction engine and a separate temporal adapter,
respectively. Shared browser integration for the temporal adapter is pending.

Local development used a conventional **2016 Intel Core i5-6200U laptop**.
Larger campaigns used hosted compute. This is a hardware-provenance statement,
not a claim that every full campaign runs on an old laptop. Browser inspection,
native launchers, Docker and Streamlit provide distinct access paths.

## A short first reproduction: three commands

With Python 3.11+ and Git available, choose a persistent research directory, then:

```sh
git clone https://github.com/FacundoFirmenich/crispr-observability.git
cd crispr-observability
python -B reproducibility/reproduce.py --mode quick --output replay-quick
```

The output directory must be new. Inspect `replay-quick/REPRODUCTION.json` for the
executed scope and outcome. Quick mode reconstructs exact CFTR endpoints from
archived counts, resolves the actual Clone3 example and replays preserved
ONT/HiFi locus decisions. It is not a fresh sequencing or alignment campaign.
Elapsed time depends on your system and initial repository transfer.

For a graphical start without installing Python, use the packages in
**[native release v0.10.0](https://github.com/FacundoFirmenich/crispr-observability/releases/tag/v0.10.0)**:
`START_HERE.bat` on Windows, `START_HERE.command` on macOS or `START_HERE.sh`
on Linux. The [reviewer workspace](reproducibility/README.md) also documents
Docker, Streamlit, the interactive notebook and granular replay.

## What the evidence establishes

| Finding | Observed result | Scope |
| --- | --- | --- |
| Product-dependent anchor observability | A +3-base shift takes a required anchor to [213,225) in a 224-base read; one CFTR lower contrast endpoint changes from −10.09177 to +13.23196 percentage points under a common quality mask | Finite recorded universe and declared label assumptions; the earlier interval includes zero |
| Label agreement can conceal unequal strings | 83 unequal strings among 170,808 jointly resolved records, including 18 with agreeing native-local labels | A correspondence counterexample, not 83 adjudicated molecular errors |
| Target-relative certification | The target is resolved when every candidate in a complete declared domain gives the same determinate value | Model-relative certificate with explicit abstention |
| The cost of withdrawal | A fixed four-donor germline rule withdraws 25/26 discordant and 113/207 concordant calls across 241 loci | Locus-level benchmark agreement under the frozen selection and rule |

The [evidence map](reproducibility/EVIDENCE.md) connects each result to its unit,
input and replay operation. Local disagreements and failed gates remain visible.
The [scientific scope](docs/SCIENTIFIC_FRONTIER.md) and
[practical workflow](docs/PRACTICAL_WORKFLOW.md) explain how to reuse the method.

## Versions, citation and access

The **v0.8.1 source snapshot** is archived at
**[10.5281/zenodo.24121127](https://doi.org/10.5281/zenodo.24121127)**.
The **v0.10.0 native packages** have their own release and recorded build commit.
Current graph and laboratory development postdates those archives: cite its exact
Git commit alongside the component. A development-branch name is not a substitute
for a versioned citation. New source components are not retroactively included in
the existing DOI or native binaries.

The repository is public. The complete manuscript and formal supplement remain
private for author deposit. Public code, numerical evidence and interactive tools
are accessible here; manuscript privacy does not imply that their source code is private.

## Contribute a scientific check

Use [scientific feedback](https://github.com/FacundoFirmenich/crispr-observability/issues/new/choose)
for a localized discrepancy, failed reproduction or supported improvement.
Include the commit, input identity, command, observed result and affected claim.
See [CONTRIBUTING.md](CONTRIBUTING.md). Independent use and community adoption
will be reported when measured; an empty issue tracker is not a performance test.

## Documentation and provenance

- [Repository guide](docs/REPOSITORY_GUIDE.md): choose a workflow and evidence scope.
- [Research environment](docs/COMPANION_CAPABILITIES.md): graph, lab, API, MCP and replay.
- [Native and dashboard reproduction](reproducibility/README.md).
- [Evidence and reproduction map](reproducibility/EVIDENCE.md).
- [Source and historical release coverage](extensions/retained_span_audit/RELEASE_COVERAGE.json).
- [Integrated generation-model evidence](extensions/generation_model/README.md).
- [Visual conventions](docs/VISUAL_STYLE.md), [publication policy](PUBLICATION_POLICY.md)
  and [acknowledgments](ACKNOWLEDGMENTS.md).

## Licensing

Current original observability software, the portable runner and integrated
laboratory use **GPL-3.0-only**. The independent source-allele extension under
`extensions/generation_model/` retains **MIT**. Earlier MIT grants and third-party
notices are preserved. Source data and publications retain their own terms.
See [LICENSE_SCOPE.md](LICENSE_SCOPE.md); no code licence assigns rights over
the private manuscript or relabels a third-party dataset.
