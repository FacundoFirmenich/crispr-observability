# Read Observability Lab: multi-input research engine

[Castellano: guía general](GUIA_USUARIO_MOTOR_MULTIENTRADA_v0.4.md) ·
[Castellano: uniones de splicing](GUIA_EVIDENCIA_SPLICING_v0.5.md)

The v0.5 engine executes versioned study manifests containing named inputs and
ordered operations. It extends the v0.4 contract with exact same-record junction
evidence. CLI, JavaScript and MCP use the same scientific kernel.

## Start

With Node.js 22+ in a research workspace:

```sh
npm ci --prefix engine --no-fund --no-audit
node engine/project-cli.mjs check examples/cftr-leading-case/project.json
node engine/project-cli.mjs run examples/cftr-leading-case/project.json
```

The leading case uses documented CFTR coordinates and an authentic archived read.
See [the complete demonstration](../docs/LEADING_CASE.md).

## Inputs and operations

Up to 32 inputs and 64 operations are accepted per manifest. Named inputs may be
single-record FASTA, portable sequence/structure JSON, PDB, textual mmCIF, XYZ,
SAM or indexed BAM. Paths are relative to the manifest directory, remain inside
the authorized workspace and may include an expected SHA-256.

| Operation | Result |
| --- | --- |
| `inspect` | Format, sequence length, structural atom count or alignment summary |
| `anchor-containment` | Shifted interval, observed positions and missing boundary |
| `compare-sequence-spans` | Oriented observed strings and equality; overlap only within the same input coordinate space |
| `atom-distance` | Distance in angstroms between two normalized atoms |
| `alignment-locus` | CIGAR projection, reference coverage, selected record and evidence TSV |
| `junction-evidence` | Conjunction of exact CIGAR N intervals in one record of the specified mate; separate any/unique/primary endpoints |

Coordinates are zero-based and half-open. FASTA accepts 2–360 ACGTN bases. Text
input is limited to 6 MB; structural cases to 30,000 atoms. Loci are limited to
50,000 bases and 2,500 records, with indexed-reading bounds documented in the lab.
An `inspect` or `check` result does not establish that unused inputs were read.
Missing, incompatible, out-of-workspace or checksum-mismatched inputs fail explicitly.

For junction evidence, each required interval must occur in the same alignment
record of the chosen mate. Separate mates or alternative alignments cannot be
combined into a witness. `exactAnyReported`, `exactUniqueReported` and
`exactUniquePrimary` retain distinct semantics. NH determines the declared
uniqueness check; MAPQ 255 remains unavailable rather than silently rejecting
a junction witness. Base quality, NH and MAPQ are separate fields. A local
not-observed result concerns the supplied locus, not an entire library.

## API and MCP

```js
import {readManifest, runProject} from './engine/project-engine.mjs';
const report = await runProject(await readManifest(root, 'project.json'), root);
```

```sh
node engine/project-mcp.mjs /absolute/path/to/studies
```

The stdio server exposes `lab_project_run` with `manifest_path` and optional
`check_only`. Separate study directories retain their own inputs and namespace.
Reports contain normalized manifest identity, effective parameters, source
provenance and ordered results. Persist the report with the inputs and code commit.
Rendering and browser replay use the [laboratory MCP bridge](../interactive/lab/MCP.md).

## Validation and interpretation

Public-fixture checks compare API, CLI and MCP outputs, BAM/SAM decoding,
namespace isolation, rejected invalid inputs and exact junction semantics.
The junction fixture is a labelled software test, not biological evidence.
Biogen has separately exercised the junction operation on its retained biological
records; its frozen sample and endpoint scope remain that study's responsibility.
The shared engine is not a guarantee of identical scientific targets across studies.
