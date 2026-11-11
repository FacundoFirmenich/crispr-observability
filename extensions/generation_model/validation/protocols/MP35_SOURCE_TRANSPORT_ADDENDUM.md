# Source transport correction before MP35 scientific outcomes

The published PAW87816 CRAI does not point at valid CRAM container headers in the retrieved object. Parsing the CRAM v3 file-header container gives the first data-container offset9308, whereas the published CRAI lists9494. The derived offset correction is minus186 file-position units. Three independent positions (first two containers and the infrastructure pilot chr22 container) satisfy the container data-length plus header-length relation after this correction. The original CRAI is preserved; a separate rebased CRAI records the derived correction. Source-author confirmation has not been obtained and no universal integrity claim is made.

This is an access-layer repair, not an allele, coordinate, sequence, label, cohort-selection or threshold change. Every used container must decode successfully with CRAM integrity checks and the exact reference. A failed decode is a source failure, never a biological negative. Raw range diagnostics, original/rebased hashes and decoding receipt are retained. Specification: https://samtools.github.io/hts-specs/CRAMv3.pdf .

MP35 scientific analysis remains forbidden until the pilot decode succeeds. The frozen selection and evaluation rules remain unchanged. The 186-unit offset is recorded as a machine coordinate, not a storage-size report; storage is reported in MB/GB.
