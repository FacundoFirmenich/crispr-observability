# Public archive inputs, attribution and separation of rights

Original biological experiments: Sousa et al., Nature Biomedical Engineering, DOI `10.1038/s41551-024-01233-3`, BioProject `PRJNA1055086`. The public SRA accessions are SRR29623265, SRR29623266, SRR29623267, SRR29623268, SRR29623270, SRR29623271, SRR29623272, SRR29623273 and SRR29623274. These are genuine source reads, not simulations.

Bundled FASTQ files are the attributed public-archive snapshot used for the preserved analysis. NCBI SRA processing documentation states that public sequence data may be retrieved and redistributed. INSDC policy provides access without restrictions on use or redistribution; INSDC does not accept controlled-access submissions. These policies support redistribution of this public SRA snapshot. They do not turn source data into original software or grant MIT terms over the source authors' work. NCBI's general policy distinguishes its access policy from third-party rights.

Primary policies checked 2026-10-28:

- https://www.ncbi.nlm.nih.gov/sra/docs/sequence-data-processing/
- https://www.insdc.org/policy/
- https://www.insdc.org/about-insdc/
- https://www.ncbi.nlm.nih.gov/home/about/policies/

The 18 detailed allele tables are frozen programme outputs under CRISPResso2 2.3.4, commit `1fc4b629e5400a48d6569b0ab4b8d1766f33636f`, using observed-mutant and programmed references. They are attributed derived alignment data. This package does not contain or license the CRISPResso2 engine. Pertinent frozen configurations remain under `metadata/`.

Replay adapts the programme's M1 `separate_quality_geometry.py`, M2 `verify_m2.py` and M3 `build_witness.py` contracts. Manifests preserve source implementation and original-result hashes. The analysis does not reassign experimental improvements to this programme.

The original-code MIT allowlist excludes FASTQ, alignment tables, metadata and scientific results. This notice records provenance and public-archive redistribution scope; it does not replace third-party terms. The manuscript remains a separate scientific text object.

Re-extracting or recompressing SRA records may change bytes or order. This replay verifies bundled bytes and logical sequence/quality streams; it does not silently replace discrepancies. Accessions support independent retrieval without promising byte identity from another extraction method.
