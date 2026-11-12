# Local alignment auditing

Open **Read alignments → Open public HG00096 example**. The record has FLAG 145, MAPQ 0 and CIGAR `183S67M`. Stored SEQ index 0 lies in soft clipping. Toggle original orientation and index 0 maps to reference index 10065. Changing minimum MAPQ from 0 to 20 changes the observed reference positions from 67 to 0 while holding the source read and CIGAR fixed.

Load a SAM slice or a BAM plus its matching BAI/CSI. BAM reads a local indexed locus. The inspection table distinguishes stored query, original query, reference index and zero-width insertion anchor. The 3D scene is a schematic of up to 100 query positions. The depth plot shows local maxima in display bins, with exact depth runs in the TSV and replay JSON.

## Evidence contract

All coordinates are zero-based and half-open. SAM POS is converted once on import. For FLAG 0x10, stored SEQ follows reference traversal. Original orientation reverses the query index and complements the base. Hard clips contribute to the unclipped-original coordinate offset and supply no bases.

M, = and X contribute depth. D and N consume reference without contributing a base. An insertion has an anchor boundary and no reference base. A soft clip belongs to stored SEQ without a mapped reference base. M alone leaves match or mismatch status unspecified. Coverage counts alignment records. MAPQ 255 is unavailable and excluded separately. Secondary and supplementary records are excluded by default. Duplicate and QC-failure flags are preserved without silently excluding their records.

An indexed query accepts up to 50,000 reference bases, 2,500 mapped records and 24 MiB of indexed reads. Coverage stays within the queried locus, with a visible error if a limit is exceeded. Small BAMs carry a full SHA-256. Larger inputs retain slice and index digests without a full-file hash. The tested binary input is BAM/BAI. CSI is accepted by the reader and awaits an independent fixture. CRAM, index-to-file authenticity and reference-backed M mismatch inference are outside the verified scope.

## Output and validation

The browser exports the case, exact CIGAR/coverage TSV, view PNG and replay. The MCP tool `lab_alignment_locus` produces TSV, replay and a SHA-256 manifest from a workspace SAM/JSON or BAM with its index. `lab_render` accepts that replay and produces the scene PNG and the same exact TSV. The TSV's query intervals remain in stored orientation and the replay records the display orientation.

The Python oracle independently decodes public BAM fixtures using gzip and struct, then accumulates exact depth per position. The JavaScript regression suite compares its results to the indexed BAM and SAM paths. The three windows cover a reverse soft clip, a long CIGAR with 5,227 operations, and an overlapping 12-record window. These are transport and coordinate checks for existing alignments.

- [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
- [Pinned public HG00096 and CHM1 fixtures](https://github.com/GMOD/bam-js/tree/ca0836f866c65e61c70cc89199cf93032c27de1f/test/data)
- [GMOD BAM reader](https://github.com/GMOD/bam-js)
