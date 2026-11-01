# CRISPR P1: complete finite-record reproduction extension (MP05)

This directory is a subordinate computational output of **CRISPR P1**. The article and formal supplement remain private pending the author's preprint deposit.

It closes source-distribution gaps in the former companion: all fifteen nonduplicate streams can now exercise the original edge/inner/target geometry and quality rules; the original ambiguity-veto code, external catalogue parser and 36-card/24-contrast checker are present. Six additional original FASTQ containers supplement the nine already under `raw_replay/inputs/fastq/`. Repeated deposits are not counted as new streams or independent biological units. Original code and decoded scientific comparator bytes are preserved exactly. Large comparator JSON files use lossless gzip containers with original decoded hashes retained in the source manifest. Third-party data retain their source terms.

From the repository root, use a **new directory outside the repository**:

```sh
python -B extensions/reproduction/replay.py --work /tmp/crispr-p1-full-replay
python -B extensions/reproduction/analyze_correspondence.py --output /tmp/crispr-p1-certificate-analysis
```

The first command authenticates its scientific source files, rebuilds geometry/quality/veto outputs over 429,812 records, reconstructs the three external assignment catalogues and builds the reporting contract. Scientific outputs are compared to their historical bytes. It also reconstructs the original CFTR coding catalogues, all 600 finite factorization cases, 12 whole-coding-family withdrawal certificates and 10 joint-target partitions, using the original source and a second grouping/translation verifier over real record identifiers. It does not rerun source NHEJ merging or external tutorial alignment. The second reconstructs seven correspondence certificates from the nine primary streams and frozen allele tables, and characterizes 128 native-local but nonexact-programmed witness records. It does not authenticate a full allele or a sequencing-error mechanism.

The output directory contains actual command logs, new scientific files and a verification receipt. `--only geometry`, `--only external` and `--only cftr` and `--only contract` allow separate scopes. A successful internal replay is distinct from an outside person's independent reproduction and from molecular truth.

The original `observed-contracts` package is installable without third-party runtime dependencies:

```sh
python -m pip install --no-deps extensions/reproduction/external/software
observed-contracts-table --help
```

Its source contains the finite-catalogue factorization, joint-target partition and claim-scope checks. Passing on observed support is not proof of future sufficiency, and a constant target supplies no discrimination.

## Dedicated prime-editing comparison

`pe_mode/CONSTRUCT.json` transcribes **epeg147 / LrAS0147**, the source nicking guide and the separately recorded dead guide from Sousa et al., Supplementary Table 2, Figure 4a, Excel row 193. The source workbook identity and DOI are retained. RTT followed by PBS is passed as the extension; the stabilizing motif is not mistaken for the RTT/PBS extension. The dead guide is not passed as the nicking guide.

Obtain CRISPResso2 **2.3.4**, upstream commit `1fc4b629e5400a48d6569b0ab4b8d1766f33636f`, separately under its terms. The exact source ZIP at `https://github.com/pinellolab/CRISPResso2/archive/1fc4b629e5400a48d6569b0ab4b8d1766f33636f.zip` has SHA256 `7e4f2ae9ce9bad27f00afd53abf195431f68e02fc06b39da30d4f7696fca6468`. The engine source and binaries are not redistributed here.

```sh
python -B extensions/reproduction/pe_mode/prepare.py --work /tmp/crispr-p1-pe --engine /absolute/path/to/CRISPResso --source-archive /absolute/path/to/CRISPResso2-source.zip
python -B extensions/reproduction/pe_mode/run.py --work /tmp/crispr-p1-pe
python -B extensions/reproduction/pe_mode/analyze.py --runs /tmp/crispr-p1-pe/runs --output /tmp/crispr-p1-pe-analysis
```

Preparation verifies the archive, executable version and raw stream identities, closes and reads back each filtered gzip, and freezes commands before execution. The supplied amplicon is the reverse complement of the frozen forward reference to satisfy the engine's PE orientation convention. RNA guides are unchanged. There is no edited-reference override and no bypass of sequence checks. The automatically inferred edited reference must equal the reverse complement of the frozen programmed reference during analysis.

The dated MP05 execution completed nine dedicated PE jobs plus one baseline installation-control job; the latter reproduced all 6,668 frozen allele rows and abundances in SRR29623272. A failed orientation attempt and a failed execution using a truncated derivative gzip were retained privately. The two affected derivative files were reconstructed from unchanged raw input and restored the already frozen gzip hashes. Failed attempts are not counted among completed jobs.

`derived/PE_TOOL_VS_EVENT_COMPARISON.csv` separates native Prime-edited assignment percentages over all aligned reads from event bounds over original raw records. `PE_NATIVE_TOOL_METRICS.csv` retains window-relative modification percentages with their own denominators. They are different estimands. Native assignment is not native-peptide identity.

The dedicated PE projection differed from the former three-operator common string in 89 records (19/27/43 by preparation): 53 retained both binary labels, while 36 changed native-event labels and 32 of those changed exact-programmed labels. A further four records had unresolved selected-reference correspondence. Requiring agreement with this added projection returns all 93 records to U, retaining 170,632 witnessed records and 6,374 unknowns. All three native PE-minus-mock lower endpoints remain unchanged and positive; upper endpoints widen. This is conditional robustness of the stated finite-record conclusion, not comparative accuracy or universal correspondence.

The new analyses are retrospective extensions on the same biological genealogy. They are not a held-out cohort, an independent truth channel or prospective validation. Frozen historical results are retained under their original definitions rather than silently redefined by the new certificate.

## Source and distribution boundaries

`SOURCE_MANIFEST.json` authenticates original scientific sources and comparator files. `derived/` contains brief machine-readable outputs, not full manuscript or supplement text. Original external HTML reports, FASTQ files and allele tables come from the public CRISPResso2 examples and retain upstream terms. The ZIP allowlist contains only exact allele-table archives; private historical dependency bundles are excluded. The runtime engine and dependencies remain separate.
