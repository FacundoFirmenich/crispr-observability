# CRISPR P1: complete record-ledger audit, MP06 integrity addendum I1

This subordinate output retains every actual primary CFTR and HEK3 record in a
derived CSV: source accession, identifier, complete raw sequence, original quality,
quality decision and local coding from every declared operator. Empty coding is
unresolved. Quality failures and unresolved alignments remain in the denominator.

Run from the repository root into a new directory outside the checkout:

```sh
python -B extensions/evidence_integrity/export_record_ledgers.py --output /tmp/crispr-p1-complete-ledgers
```

The program authenticates the existing inputs, reconciles native allele-family
multiplicities during export, closes the gzip files, reads the saved CSVs back and
joins every saved sequence/quality/identifier to the original FASTQ in order.
From the saved coding columns, it independently rebuilds all 162 CFTR family-count
rows and all 42 HEK3 count rows and compares them to MP06. Both compressed tables
are actual data products, not lists of fingerprints. `LEDGER_VERIFICATION.json`
records their sizes, identities and readback results. The completed internal
execution is recorded in `VERIFIED_EXECUTION.json`.

There are 177,006 CFTR records (175,010 pass Q30) and 218,930 HEK3 records
(216,376 pass Q30). CFTR uses nine declared operators and HEK3 five. These cohorts
are distinct; the fifteen-stream 429,812-record geometry analysis includes the
nine primary CFTR streams and must not be added to them as independent data.

The ledgers are newly generated derived tables. They do not impersonate the
incomplete native annotated FASTQ outputs preserved in the original execution
archive. Those failures remain failures. The analytic channel uses complete
native allele tables with exact raw-family reconciliation; no absent annotation
is imputed. An empty coding cannot be turned into a negative observation.

This addendum changes no frozen biological count, event contract, confidence
claim or graph. Internal reconstruction is not independent human reproduction,
a new alignment-engine execution, comparative molecular accuracy or functional
validation. Complete article and supplementary text remain private.

Original CFTR experiments: Sousa et al., DOI 10.1038/s41551-024-01233-3.
Original HEK3 experiments: Anzalone et al., DOI 10.1038/s41586-019-1711-4.
Original observations and generated scientific data retain their source terms;
the explicitly allowlisted original audit code uses the scoped code license.
