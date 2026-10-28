# CRISPR observability: computational companion

Version **0.3.2 / MP03D**. Public code, numerical data and editable figures for
a retrospective genome-editing readout audit. The complete manuscript and
supplement are held privately pending the author's chosen preprint deposition.

## Reproduce the computational checks

With Python 3.11 or later, from this directory:

```sh
python3 code/run_review_checks.py
```

The standard-library runner verifies exact derived-table arithmetic and
reconstructs observers on nine authentic public CFTR FASTQ streams using
eighteen frozen alignment tables. It checks 177,006 original records, 170,725
common-string witnesses and 6,281 unknowns against preserved result objects.
It does not calibrate molecular classification accuracy or establish independent
biological replication. Source experiments remain attributable to Sousa et al.

## Public assets

- `figures/`: editable SVG/PDF/PNG exports, plotting code and exact input tables.
- `tables/`: nineteen authentic numerical CSVs, retained byte for byte.
- `raw_replay/`: nine genuine FASTQ inputs, eighteen frozen alignment tables,
  original comparators and portable replay code.
- `engine_replay/`: portable source and configurations for the historical
  144-job alignment sensitivity experiment, with clear preparation/execution scopes.
- `metadata/`: source identities, provenance, runtime records and release status.

`docs/REPRODUCTION_SCOPE.md` and `docs/SCIENTIFIC_FRONTIER.md` state what each
operation verifies. Conditional finite-read intervals, source software percentages,
sequence correspondence and biological truth are distinct quantities.

## Figures

Install `figures/requirements-figures.txt` and run:

```sh
python3 figures/code/generate_figures.py
```

The plotting source uses an explicit standalone coordinate contract and frozen
numerical tables. It does not require manuscript prose.

## Attribution and license

Current personal and institutional acknowledgments are in [ACKNOWLEDGMENTS.md](ACKNOWLEDGMENTS.md).
The verified author identity is in `metadata/AUTHOR_PROFILE_VERIFIED.json`.
MIT applies only to the original-code allowlist in `metadata/LICENSE_SCOPE.json`;
source observations, scientific data, graphics and third-party software retain
their distinct attribution and applicable terms. See `raw_replay/DATA_NOTICE.md`.

No preprint DOI is assigned by this repository. The private manuscript and
supplement will be linked only after the author elects to deposit them.
