# Two executable scopes and the remaining pipeline boundary

The default command is `python3 code/run_review_checks.py`. Both scopes execute in
temporary working directories; retained reports go under `results/review/`.

| Level | Operation | Status in this repository |
| --- | --- | --- |
| Derived-table arithmetic | Verify four CSV identities and the inherited checker; recompute balances, bounds, contrasts and K/K−1 certificates | Executable offline |
| Raw target observer | Reapply historical target motifs, two events, trims 0/8 and four quality/structural policies to nine original streams | Executable offline; 144 original-result rows |
| Common-string coding witness | Reconstruct motif extraction, local translation and both coordinate projections using frozen detailed alignments; reconcile raw identities and allele abundances | Executable offline; 18 counts and 12 contrasts |
| All historical observers and programme stages | Recreate every operator, intervention analysis and previous package | Not included |
| Fresh CFTR engine alignment | Run separately installed CRISPResso2 on the original quality-filtered reads under both declared references | Not executed by these checks; engine absent |
| External-example source processing/alignment | Recreate source preprocessing, NHEJ merging and alignment from tutorial FASTQ | Not included |
| Molecular or biological validation | Compare with suitably authenticated independent orthogonal evidence | Not supplied by a computational replay |

## What is genuinely reconstructed from raw

`raw_replay/run_raw_replay.py` verifies the identities of its 35 frozen inputs and
recomputes observations on all 177006 records. It checks unique within-run read
identifiers, ordered sequence/quality stream identities and conservation of both
sets of frozen allele-table abundances. It retains the declared target policies
separately. The common witness accepts only a coding string shared by motifs and
both reference projections; unresolved or discordant correspondence stays unknown.

The three generated scientific files are compared with preserved original-result
subsets under `raw_replay/expected/`. Their bytes are also compared independently
by `code/run_review_checks.py`, which retains the generated files and a receipt.
The inputs, reference definitions and comparators are not modified. Identity
checking reads comparator bytes for hashes; the observer computation does not
interpret their scientific contents as an input to its calculations.

The raw companion adapts historical contracts from M1
`separate_quality_geometry.py`, M2 `verify_m2.py` and M3 `build_witness.py`. Its input
manifest records those sources and the original result provenance. This is a
retrospective implementation on the same evidence, not another human analyst,
biological cohort or experiment. The source companion's documents and manifests
are preserved unchanged even where their wording describes delivery before MP02
integration. `metadata/RAW_REPLAY_INTEGRATION.json` records the later integration.

## What is still required for fresh alignment

The inspected historical entry point `advance_m2/code/run_external.py` filters
original FASTQ according to the frozen mean-Q30 policy and invokes a separately
installed `env/bin/CRISPResso`. Fresh reference-sensitivity execution also needs
the programmed-reference command/configuration. These invocations are distinct
from projecting coding coordinates through an existing alignment.

A full fresh-alignment distribution would require the exact external engine
source and commit acquired separately under its terms, its compatible compiled
runtime/dependencies, both command sets and parameter freezes, the source quality
policy, complete detailed engine outputs, and new execution/exit-status receipts.
The historical engine attribution is CRISPResso2 2.3.4, commit
`1fc4b629e5400a48d6569b0ab4b8d1766f33636f`. This repository redistributes no engine
code and does not claim to run it. A separate private engine execution or grid has
its own evidence and cannot be credited to this default command or CI.

## Other preserved replay levels

The inspected S1-02 `code/reproduce.py` regenerates six scientific outputs from a
complete frozen input tree and runs external/inherited contract tests. Its recorded
scope expressly excludes source alignment and NHEJ merging. That entire external
source collection and software module are not part of this minimal MP02 repository.
A successful raw CFTR witness reconstruction therefore does not establish complete
S1-02 tutorial processing replay, historical training replay or the full programme.

## Interpretation and rights

Arithmetic remains conditional on accepted labels, the event definition and the
finite observed universe. Raw reconstruction shows that the same records under
the preserved observation rules reproduce the original outputs. Neither step
estimates molecular error rates, establishes independent biological replication,
or proves full-allele identity or function. Runtime is reported for each actual
execution, not guaranteed across environments. Docker and remote CI remain
unexecuted until their own receipts exist.

The original experiments belong to the source authors. The raw data notice separates FASTQ, frozen tables, comparators and metadata
from the scoped MIT license for original code. Its MP03 update documents attributed
public SRA snapshot redistribution under the primary SRA/INSDC access policies.
Integration does not create third-party rights or relicense source data as software.
