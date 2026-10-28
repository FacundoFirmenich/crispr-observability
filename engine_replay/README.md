# Full 144-job CRISPResso2 sensitivity replay

This compact source package makes the executed MP02 sensitivity experiment
relaunchable from a clone of the repository. It preserves the original scripts
and scientific metadata without bundling CRISPResso2 or the 144 large alignment
output directories. `original/VERIFICATION.json` and `original/FINDINGS.json`
describe the completed historical executions. They are not receipts for a new
run or for independent reproduction.

The source study used nine authentic CFTR FASTQ streams (177,006 raw records),
two frozen references, a custom baseline and six one-at-a-time parameter changes
per reference (126 jobs), plus an upstream stock profile without a guide (18
jobs). The profiles are retrospective sensitivity checks on the same records.
The standard-library raw replay at the repository root is a separate, shorter
operation that uses frozen alignment tables.

## 1. Prepare and authenticate inputs without installing the engine

From the repository root, choose a **new directory outside the repository**:

```bash
python3 engine_replay/prepare.py --repo . --work /tmp/crispr-engine-replay
```

This checks the included source and dependency manifests, authenticates the nine
original sequence/quality streams, reconstructs the common original whole-read
Q30 inputs, and checks every filtered stream against its preserved original
hash and count. It stages 126 grid jobs and the scripts needed for the additional
18 stock-profile jobs. No alignment engine is required or executed. An existing
work directory is rejected rather than overwritten. The result receipt is
`/tmp/crispr-engine-replay/PREPARATION_VALIDATION.json`.

The `GRID_FREEZE.json` created in this new work directory freezes its **new
computational replay**. Historical results are already known; this operation
does not establish prospective validation or a previously unseen holdout.

## 2. Obtain and install the official engine separately

Use the official source <https://github.com/pinellolab/CRISPResso2>, version 2.3.4,
commit `1fc4b629e5400a48d6569b0ab4b8d1766f33636f`, following its installation and
license instructions. This repository neither redistributes nor relicenses
the engine. A separate Python environment and the upstream dependencies are
required; the native checks in step 1 do not install them.

Attach the separately installed source and executable:

```bash
python3 engine_replay/attach_engine.py \
  --work /tmp/crispr-engine-replay \
  --engine-source /absolute/path/to/CRISPResso2 \
  --engine-executable /absolute/path/to/engine-env/bin/CRISPResso
```

This verifies the source commit and absence of tracked modifications, executable
version, script command-line interfaces and actual engine `--help`. It prepares
the stock-profile freeze from the installed upstream `args.json` and requires
those default definitions to equal the recorded original definitions. It starts
no alignment jobs. Its receipt is `ENGINE_ATTACHMENT_VALIDATION.json` in the
chosen work directory. Engine attachment was **not executed** when preparing this
MP03 compact source package, because no engine is bundled or reinstalled here.

## 3. Execute and verify the 144 new jobs

After successful attachment, run from the new work directory:

```bash
cd /tmp/crispr-engine-replay
python3 publication_02/sensitivity/checked_run.py baseline --workers 3
python3 publication_02/sensitivity/summarize_grid.py baseline
python3 publication_02/sensitivity/staged_run.py stock --workers 3
python3 publication_02/sensitivity/staged_run.py oat --workers 3
python3 publication_02/sensitivity/summarize_grid.py all
python3 publication_02/sensitivity/build_findings.py
python3 publication_02/sensitivity/close_manifest.py
```

The first 18 jobs must reproduce the preserved original projected-event counts
and native software quantification tables exactly. The complete verifier then
checks all 144 outputs using original read identities, annotated FASTQ and
allele-table abundance conservation, two coordinate-projection implementations,
native quantification classes and raw-denominator interval arithmetic. Input
hashes are verified around execution blocks. The staged runner closes and checks
gzip/ZIP files, copies verified outputs and fsyncs them before accepting them.
Preparation gives each new work directory its own scratch staging location;
no scientific parameter is changed.

These commands start new alignments only when a person actually executes them.
No additional engine job was run while producing this MP03 port. A technical
failure must be preserved and diagnosed rather than counted as a completed
scientific execution. The original `checked_run.py`, `staged_run.py` and frozen
metadata are retained unchanged under `original/`; only the staged working copy
receives a unique temporary directory.

## Included and excluded objects

- `original/`: all top-level source scripts, metadata, historical summaries,
  contrasts and distinct-record label/stability artifacts from the closed MP02
  sensitivity archive. These files are copied byte for byte.
- `dependencies/m2/`: the exact original design, two references, baseline event
  counts, nine quality-count receipts and eighteen native quantification tables
  read by the original scripts. No additional historical workspace is needed.
- The nine raw FASTQ inputs remain in `../raw_replay/inputs/fastq/`.
- `ORIGINAL_SOURCE_MANIFEST.json` and `DEPENDENCY_MANIFEST.json` authenticate the
  compact package's selected scientific inputs and scripts.
- `original/CLOSED_FILE_MANIFEST.json` records the **larger historical archive**.
  Most of its 2,390 closed files are intentionally absent here. It is retained for
  provenance, not presented as this compact package's inventory.
- The original 144 complete alignment directories, failed attempts and full
  historical byte-level audit remain in the previously delivered MP02 sensitivity
  archive. This compact source package does not claim to include those bytes.

Replaying the same records adds computational traceability. It does not create a
new biological cohort, molecular ground truth, functional validation or a blinded
independent analyst. Observed label changes and aggregate interval stability
remain separate conclusions, conditional on the declared event and acceptance
contracts.
