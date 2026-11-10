# Reproduction ledger

This file is generated from the same frozen execution receipts linked by the graph. It is not the paper supplement and contains no private manuscript text.

| Check | Result | Scope |
|---|---|---|
| Independent path enumeration | PASS | 4328 sequence pairs; 176305 paths; 33928 window/ownership contracts; 410 selected target maxima |
| Retained COL17A1 replay | PASS | 439589 raw pairs; 20260 reference cases; 81112 attaining witnesses |
| Biological accuracy | NOT ESTIMATED | No new same-molecule truth for the original CFTR/HEK3 records |
| Fresh extraction of all 72 COL17A1 FASTQ files | NOT EXECUTED | Reproduction starts with retained spans and complete pair ledgers |
| Second OS / fresh dependency installation | NOT EXECUTED | Installed Windows/Python environment only |
| Interactive browser testing | NOT EXECUTED | Browser use unavailable; static graph, data and model checks are separate |

The zero-count serialization mismatch and its repair remain in `checks/REPLAY_FAILURE_AND_REPAIR.json`. No observed counts were changed. Do not translate the number of passing computational checks into a biological accuracy score.

## One entry point

Install `requirements.txt`, then run `python reproduce.py --mode quick --output new-quick-output`. This verifies every distributed hash, checks the graph and reproduces the authentic Clone3 example. Run `--mode full --output new-full-output` for all retained-pair counts, attaining witnesses and the independent bounded path oracle. Both commands require a new output directory and leave source data unchanged. Runtime is environment-dependent; the full components took about two and a half minutes on the authoring machine.

For all earlier CFTR/HEK3 analyses, use the pinned source repository in `DEPENDENCY_BASELINE.json` and its documented entry points. See `RELEASE_COVERAGE.json` for the distribution boundary of every evidence family.
