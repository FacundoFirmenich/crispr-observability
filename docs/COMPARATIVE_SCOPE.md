# Comparative scope and documented capabilities

This is a source-grounded capability comparison, not a performance ranking. Unreported features are not proven absent. CRISPResso2 reports input, preprocessed and aligned counts; it can retain/handle ambiguous assignments. Its native output is not characterized as hiding all uncertainty. ampliCan control normalization removes background events without necessarily removing the reads carrying those events. Different event definitions cannot be ranked by subtracting percentages.

| Tool | Documented input and operation | Denominator/uncertainty context | Evidence in this work |
|---|---|---|---|
| CRISPResso2 | Amplicon alignment, editing windows, prime-editing references | Input/preprocessing/alignment statistics; configurable ambiguous assignments | Nine fresh 2.3.4 jobs, authenticated same-record projections and preserved historical replay |
| ampliCan | Amplicon alignment, filtering, control normalization | Quality/primer filtering and event-level control correction | All nine CFTR libraries under three declared computational anchor configurations; exact admission transitions |
| CRISPR-A | Analytics, normalization and simulation pipeline | Pipeline classifications and controls; published human adjudications | Existing HEK3 analysis plus nine new CFTR variant stages; explicit native-interface audit and repair intervention |
| CRISPRLungo | Long-read amplicon editing, UMI consensus and outcome classification | Consensus/filtering and structural-event analysis | Contextual comparator; not a short-read CFTR benchmark |
| ALPINE (Maggi-Chen) | Long-read amplicon editing outcome classification | Alignment-based outcome categories | Contextual comparator; no execution claimed |
| CRISPR Observability | Declared event, observed-string correspondence, conditional finite-record bounds | E/NE/U retained on acquired N; certificates relative to the declared candidate domain | Authentic CFTR FASTQ entry, fresh native comparisons, local donor-held-out routing and an exact decision-cost map |

Complete all-optimal-domain target certification and the exact E/NE/U reporting contract are our explicitly implemented layers. The cited comparator sources have not been audited exhaustively for every possible equivalent feature; no categorical universal absence claim is made. The new FASTQ entry is an exact-anchor observer, not that complete alignment certificate.

Primary sources (checked 13 November 2026):
- CRISPResso2 documentation: https://crispresso2.pinellolab.org/help
- ampliCan: https://pmc.ncbi.nlm.nih.gov/articles/PMC6499316/
- CRISPR-A: https://doi.org/10.1371/journal.pcbi.1011137
- CRISPRLungo: https://doi.org/10.1038/s41551-026-01776-7 and https://github.com/pinellolab/CRISPRLungo
- ALPINE: https://github.com/Maggi-Chen/ALPINE and https://doi.org/10.64898/2026.03.27.714831

The original comparison's CRISPResso2 version is pinned to 2.3.4; current documentation provides context and does not retroactively change the executable version. Reproduce the historical accounting with `python code/compare_frozen_callers.py --output /persistent/new-comparison.json`. This recalculates archived outputs and is not a new engine execution.

The [current comparative extension](COMPARATIVE_UPDATE.md) preserves all native categories and every adverse decision-cost regime. [Native execution adjudication](NATIVE_EXECUTION_ADJUDICATION.md) distinguishes successful computation from a validated return contract. The historical replay command above retains its original scope.
