# P1 figures: editable source and vector exports

The four figures are re-rendered from frozen MP02 tables and the independent coordinate contract in `sources/anchor_geometry_contract.json`. The new SVG/PDF outputs consist of native vector paths and text; they contain no embedded raster images. The original MP02 PNGs and the corrected Figure 2 PDF are preserved byte-for-byte in `originals/`. This is a documented rendering revision, not a claim to have recovered the original plotting code for Figures 1, 3 and 4.

## Files

- `figures/fig1.pdf` through `fig4.pdf`: vector figures for LaTeX inclusion.
- `figures/fig1.svg` through `fig4.svg`: editable vector figures with live text.
- `figures/fig1.png` through `fig4.png`: 300 dpi previews.
- `code/generate_figures.py`: the complete plotting and exact-input extraction code.
- `inputs/`: small, extracted figure tables, including exact rational interval endpoints.
- `sources/`: unchanged MP02 source tables, witness replay output and the standalone Figure1 coordinate contract.
- `SOURCE_MANIFEST.json`: source paths, sizes and SHA-256 hashes.
- `verification/`: input/identity checks, vector structure checks and PDF page renders used for visual review.

## Rebuild

Tested with Python 3.12.14, matplotlib 3.10.8 and numpy 2.3.5. Python's standard library handles the exact rational computations.

```sh
python3 -m pip install -r requirements-figures.txt
python3 code/generate_figures.py
```

The code first verifies every frozen source hash. It then extracts the figure input tables and renders all four figures as PDF, SVG and PNG. Run it from any working directory. No network access or parent workspace is required after dependencies are installed.

## Scientific provenance and interpretation

| Figure | Source and operation | Scope |
| --- | --- | --- |
| 1 | Standalone coordinate contract: half-open anchor coordinates `[210,222)`, `[213,225)`, `[206,218)` and read end 224 nt. | A coordinate illustration. It is not a new observation or an assertion that every product has this geometry. |
| 2 | `Table_1_Native_Contrast_Intervals.csv`, all 21 exact two-arm intervals. | Conditional identification hulls in raw denominators, not sampling confidence intervals. All endpoints fit within the corrected common x-axis domain `[-12,45]`; the upper endpoints above 38 in the middle panel remain visible. |
| 3 | Six `common_raw_meanQ30` native peptide intervals for CFB-26-16, PE versus mock, from Table S3. The seventh uses the original-ambiguity-veto treated and mock E/U/N counts in Table S4 and the exact two-arm hull formula. | A finite truncation probe. The original veto prevents a shortened read from resolving uncertainty that was already present. It does not establish biological truth. The lower panel is a zoom of the three target rows, not extra data. |
| 4 | Table S9 yields 34, 58 and 45 mutant non-event to programmed event label changes. `WITNESS_BOUNDS.json` supplies nine raw-stream domain partitions, summed to 83 different-coding, 1,996 quality-failure, 4,202 missing-correspondence and 170,725 common-string-witness records. | The four categories partition 177,006 raw records. Label disagreement, sequence correspondence and unknown status remain distinct. The right panel has a log count axis. These are computational checks of the same data, not independent biological replication or comparative classification accuracy. |

The rendering revision preserves the numerical inputs and labels while making every figure available as editable vectors. It does not broaden any scientific claim in MP02.
