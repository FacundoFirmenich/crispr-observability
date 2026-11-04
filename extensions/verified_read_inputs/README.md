# Verified public KEAP1 sequencing inputs

This extension pins six public FASTQ files for KEAP1 from the CRISPRnano study. The [original publication](https://doi.org/10.1093/nar/gkac440) and [ENA study PRJNA828591](https://www.ebi.ac.uk/ena/browser/view/PRJNA828591) provide provenance for the sequencing. The manifest links each accession to the platform, original submitted alias and exact ENA export. It records the ENA length and MD5 plus a SHA-256 computed from the verified download.

From the repository root, use Python 3.10 or later.

    python extensions/verified_read_inputs/verify_fastq.py data/keap1

The script writes each download to a temporary file, checks its length, hashes, FASTQ structure, accession prefixes, unique identifiers and record count, and then moves it into place. Append --offline to check files already present without a network request. All six inputs are available from their original archive and retain the source terms.

This check establishes the identity and readability of the source files. Biological interpretation depends on the published experiment and its source linkage.
