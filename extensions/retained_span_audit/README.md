# COL17A1 retained-span reproducibility extension

This scientific companion extends `crispr-observability` with the retained
COL17A1 spans, complete 439,589-pair identity ledgers, fixed targets, score
results and a read-only reproduction runner. The paper and supplementary text
remain private. This is computational reproduction, not new sequencing,
independent human replication, or a calibrated biological accuracy benchmark.

Tested environment: CPython 3.14, with the exact package versions in
`requirements.txt`. Install them in a separate environment. A clean extraction
was checked with the same installed interpreter; no second OS or fresh package
installation is claimed. No internet access is required for the three replays.

```
python -m pip install -r requirements.txt
python replay_col17_mp21.py --mode example --output example_receipt.json
python replay_col17_mp21.py --mode counts --output counts_receipt.json
python replay_col17_mp21.py --mode full --output full_receipt.json
```

Outputs must not already exist. The runner hashes every distributed input
before running and writes only the specified new output (Python may also write
normal module/compiler caches). The example uses the actual Clone3 modal span:
net -6 versus JEB, two equivalent starts 40/41, and first competing target-class
deficit 45 under both references. It does not infer a unique genomic coordinate.

`counts` independently aggregates every retained pair using all per-target
class-score deficits and verifies both all-available-span and strict
paired-fragment admission policies. `full` additionally reconstructs, rescores
and reprojects every attaining witness for 10,130 clone/span cases under two
references. Expected strict deficit-36 counts:

| Library | Any E | All NE | Unresolved | Raw pairs |
| --- | ---: | ---: | ---: | ---: |
| DNA | 172149 | 7185 | 35004 | 214338 |
| cDNA | 180758 | 4141 | 40352 | 225251 |

The full replay checks 20,260 reference cases and 81,112 attaining witnesses.
This reuses retained exact source spans and pair ledgers; it does not download
or re-extract all original FASTQ files. `ADMISSION_MP15B/PROTOCOL.json` lists all
72 public source URLs, file checksums, anchors and admission rules. Original
data: Bischof et al., DOI 10.1016/j.ymthe.2022.04.020, Zenodo record 6242279.
The original study retains attribution for its biological experiments.

`NONMODAL_MP16` includes the separately frozen recurrent nonmodal panel,
sequence-support results and explicit source-sign reconciliation. Repeated
DNA/cDNA support is not same-molecule truth. Historical version identifiers are
preserved as source locators. The earlier CFTR/HEK3 companion is pinned in
`DEPENDENCY_BASELINE.json`; this extension is not a blanket assertion that every
historical analysis in the manuscript has been freshly rerun or fully archived
here. A release coverage inventory states remaining distribution boundaries.

## Unified reproduction and evidence navigation

The preferred commands are now:

```
python reproduce.py --mode quick --output new-quick-output
python reproduce.py --mode full --output new-full-output
```

The full command combines every retained-pair count and witness replay with
the independent bounded path oracle. The quick command uses an authentic
Clone3 example and checks every distributed hash and graph dependency.
Open `epistemic_graph.html` locally for the original standalone viewer, or read
`epistemic_graph.json` and `CLAIMS_LEDGER.md`. The HTML contains no remote script,
tracking or third-party runtime. Artifact links work inside the extracted archive.
The graph is an evidence navigator, not an inference of biological accuracy.

`RELEASE_COVERAGE.json` identifies the separate historical scientific archive,
the pinned earlier repository, the current reruns and unexecuted checks.
`REPRODUCTION_LEDGER.md` is generated from the same receipts as the graph.
`ENVIRONMENT.json` records the tested interpreter and direct dependencies;
it is not represented as a fully resolved cross-platform lockfile.
