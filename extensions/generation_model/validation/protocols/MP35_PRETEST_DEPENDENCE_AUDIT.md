# Pre-test dependence audit

Before executing any test batch, require a report of query-name overlap between calibration and test and across test loci. The original protocol forbids shared loci, but long reads can span multiple loci; it does not guarantee independent molecules or blocks. Do not silently change primary retention or recalibrate after seeing test outcomes. If shared read IDs occur, preserve the original result and mark the independence limitation; any disjoint-read sensitivity must be explicitly secondary and separately denominated. Even zero query-ID overlap does not establish independent libraries or independent biological replicates.

This adds an audit and an interpretation boundary, not an outcome-dependent selection rule. Nominal conformal coverage remains conditional on assumptions; observed test coverage is reported descriptively.
