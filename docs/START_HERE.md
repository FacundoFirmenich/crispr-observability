# Start with a biological question, leave with a decision report

**Question:** do these acquired reads resolve a positive treated–mock contrast for my explicitly defined sequence event?

Open Streamlit and choose **Start here**. The authentic CFTR example is already loaded. No graph terminology, Node installation or 3D controls are needed for this first result. Read the four steps, download the HTML decision report and retain the JSON input/output. The laboratory and evidence maps are optional deeper inspections.

For a terminal, Python 3.10+ and the source repository are sufficient; no Python packages are required:

```sh
python code/observability_report.py demo --output /your/persistent/location/cftr-review
```

Open `report.html`. It shows the event, the original record counts, the exact contrast, the unresolved fraction and the next measurement question. `input.json` is an actual CFTR input contract; `report.json` contains exact fractions. The command refuses to replace an existing report.

## Four steps

1. **Define the event.** Exact programmed DNA, local in-silico peptide and a caller's generic modification label are different targets. State reference, coordinates, quality policy and candidate domain.
2. **Account for observations.** Every acquired record belongs to E, NE or U. Counts must sum to the original N; filtering must not silently shrink the denominator.
3. **Interpret the decision.** A lower contrast bound above zero means every binary completion of U remains positive under the accepted-label assumption. An interval crossing zero identifies a measurement/interpretation gap to investigate. Neither result is a confidence interval.
4. **Select the next action.** Inspect absent flanks or ambiguous strings, test accepted-label sensitivity, or acquire independently linked evidence. The 3D laboratory tests coordinate interventions; the Supergraph identifies the premises behind a conclusion.

## Your own counts

Edit the downloaded input contract and run:

```sh
python code/observability_report.py analyse your-study.json --output /your/persistent/location/your-study-review
```

This entry point accepts already classified counts. It does not classify arbitrary FASTQ files or certify the completeness of an alignment-path domain. Those operations use the documented caller, alignment engine and certificate workflows. The report validates arithmetic and declared metadata; the researcher supplies the evidence supporting those declarations. Keep each library/event jurisdiction separate before aggregating disjoint units.

## What does abstention cost?

The demo includes eight donor-by-tract jurisdictions from the existing four-donor germline comparison. It reports discordances removed alongside concordances sacrificed. Zero retained calls means retained-call discordance is undefined. Under the explicit decision loss `L = c_D D + c_U U`, the aggregate majority rule improves loss only when `c_D / c_U > 138/25 = 5.52` (for positive abstention cost). The appropriate costs are a user decision, not fitted from the benchmark. The report preserves local differences; the aggregate threshold does not authorize a universal filter.

## Functional evidence

Published CFTR ion-channel measurements provide biological context. To correlate them with sequence bounds, first authenticate donor, culture/pool, condition, time and assay/library identifiers. Repeated `Replicate 1` labels do not establish that link. See `FUNCTIONAL_BRIDGE.md` for the recovered values and the pairing sensitivity.

## From FASTQ to a report

For a declared exact DNA insert between unique exact anchors, run the standard-library entry point on your persistent input files:

```sh
python code/fastq_decision.py --treated treated.fastq.gz --control mock.fastq.gz --config examples/observation-report/cftr_fastq_config.json --output /persistent/new-report
```

The supplied configuration is specifically CFTR. For another locus, explicitly replace both anchors, the expected insert, reference description, maximum span and quality threshold before examining the outcome. The observer checks both orientations, retains low-quality/missing/ambiguous-anchor records as U, rejects malformed FASTQ and duplicate read identifiers, and refuses to overwrite reports. The parser accepts four-line FASTQ with uppercase ACGTN and Phred+33; other encodings need explicit conversion. It streams the records and keeps an identifier set in memory. It does not imply unique original molecules or certify all alignment paths.

The CFTR exact-DNA example gives 13.11530–22.69773 percentage points for the treated–mock contrast. Its target is different from the native-local three-operator witness used by the count-only demo. Do not mix their labels or compare them as accuracy estimates.

Open [the guided tutorial](../interactive/CFTR_TUTORIAL.html). Its 15–20 minute duration is a design target, not a measured usability result. See [comparative scope](COMPARATIVE_SCOPE.md) and [prospective validation](PROSPECTIVE_VALIDATION.md).

## Direct sequence input in Streamlit

In **Start here**, expand **Analyse your FASTQ**. Upload treated and control FASTQ or gzip files, inspect or edit the CFTR configuration, then select **Compute sequence evidence**. Download the HTML report and the combined JSON evidence/contract. Recompute after changing inputs. Files are processed in host memory without persistent upload files; a hosted app processes data on its host. Limits are 20 MB and 250,000 records per file. The CLI streams larger inputs. The exact-anchor observer is not all-alignment certification.
