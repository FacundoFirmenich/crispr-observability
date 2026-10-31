# CRISPR P1: supplementary reference and recording-context audit

This is an output of the CRISPR paper, alongside its tables and figures. It
adds post hoc public-reference checks and a separately implemented raw FASTQ
audit. The complete article and formal supplement remain private pending the
author's chosen preprint deposit. The earlier three interactive companions
and all original numerical results remain frozen.

## Results and scope

- `derived/molecular/`: exact canonical/proposed coding context, the local
  peptide, codon table and a three-record reference-set identity counterexample.
- `derived/raw_qc/`: all nine original streams, 177,006 records; the historical
  whole-record meanQ30 gate again yields 175,010 passes and 1,996 failures.
  No trimming, new alignment, altered event or molecular accuracy estimate.
- `sources/`: versioned public retrieval snapshots and compact sample metadata.
- `accessible/`: Seq2Music 0.3.2 carriers and actual local decode reports for
  three public reference inputs. MIDI, WAV and MusicXML statuses concern exact
  normalized symbols, not DNA recovery from a peptide or biological acoustics.
- `SCIENTIFIC_PLUGIN_SCOPE.json`, `CAPABILITY_ADJUDICATION.json` and
  `LOOKUP_OUTCOMES.json`: active scientific capability inventory, task decisions,
  and successful, empty and failed source lookups. Overlapping plugins and
  sources are not independent evidence channels.
- `../../figures/scientific_context/`: two new supplementary figure sets,
  editable SVG plus PDF and PNG exports. These are figure files, not a public
  article or supplementary-manuscript PDF.

The same declared local peptide occurs in published reference records of unequal
complete lengths. This bounds the target that local identity can determine; it
does not assign those isoforms to the study's samples or establish their function.
AlphaFold copies reference sequences and is not an independent molecular truth
channel. Canonical and programmed DNA differ at five bases while preserving that
local peptide. The coding deletion c.1521_1523delCTT crosses codons 507 and 508.

## Rebuild numerical checks offline

From the repository root with Python 3.11 or later:

```sh
python3 -B extensions/scientific_context/check_extension.py
```

The standard-library checker reconstructs the added reference and raw-QC outputs
in a temporary directory and compares their exact bytes. It checks frozen carrier
files and the actual saved local round-trip reports; it does not rerun the
third-party carrier decoders in CI. Raw inputs are the original repository files.

The reference snapshots are the inputs, not current API state. Refreshing them
would create a new dated analysis, not silently change this one. Installed
source-specific skills performed the logged lookups. No new model generation,
training, docking, molecular dynamics or clinical experiment was performed.

For optional figure rendering, install `requirements-figures.txt` and run:

```sh
python3 extensions/scientific_context/generate_extension_figures.py \
  --results extensions/scientific_context/derived \
  --molecular extensions/scientific_context/derived/molecular \
  --output results/scientific_context_figures
```

## Attribution

Original CFTR experiments remain attributable to Sousa et al.,
doi:10.1038/s41551-024-01233-3. Reference data retain their provider attribution
and applicable terms: [UniProt P13569](https://www.uniprot.org/uniprotkb/P13569/entry),
[Ensembl CFTR transcript](https://www.ensembl.org/id/ENST00000003084),
[AlphaFold reference records](https://alphafold.ebi.ac.uk/api/prediction/P13569),
[RCSB 5UAK](https://www.rcsb.org/structure/5UAK), and
[ClinVar 7105](https://www.ncbi.nlm.nih.gov/clinvar/variation/7105/).
Targeted source hygiene checked PubMed and PMC; IndelQuant
(doi:10.3390/ijms27209054) is additional context, not a pipeline run here.

Seq2Music is by Jacob Vogan, distributed under MIT; its full notice is retained
in `SEQ2MUSIC_LICENSE`. Generated scientific content and public reference
sequences retain their separate attribution and terms. Only the repository's
explicit original-code allowlist carries the author's MIT grant.
