# CFTR sequence–function bridge: source recovery and admission

We recovered the complete published source workbook and supplementary workbook in memory from [Sousa et al.](https://doi.org/10.1038/s41551-024-01233-3). `examples/observation-report/functional_bridge.json` preserves source URLs, full-object hashes, selected cell coordinates, values and the explicit linkage decision. `code/functional_bridge.py` reproduces the extraction with Python's standard library; it does not download sequencing data.

## Observed biological context

Figure 4c reports three observations per condition. For forskolin/IBMX (F&I), the PE values are 10.0214797, 19.7655267 and 12.8944473 microampere/cm²; mock values are −0.3421667, −1.0758233 and −1.0234183. Every recorded PE value exceeds every recorded mock value. Across all nine cross-condition differences, the range is 10.3636464–20.84135 microampere/cm². This is an observed unpaired difference range, not a confidence interval or an estimate of a population causal effect.

After sign reversal of the published CFTR-inhibitor response, every PE value also exceeds every mock value; the nine cross-condition differences span 4.9167197–11.82787711 microampere/cm². The two pharmacological endpoints provide functional context for the source experiment. They do not measure the correctness of individual sequenced records.

## The link that is still missing

Supplementary Table 2 rows 191–193 identify the nine Figure 4a sequencing filenames and donors CFB-26-16, CFB-26-02 and CFB-06-04. Figure 4c labels columns only Replicate 1–3. The article describes functional recordings after three weeks of differentiation. These source cells do not independently authenticate a donor/culture/library-to-functional-assay crosswalk. We therefore do not assume that matching column order establishes the same pool or time point.

If Figure 4a's published CTT-insertion percentages are paired with F&I values in displayed column order, Pearson r is 0.916136. Across all six possible assignments of the three functional labels, r ranges from −0.998609 to 0.916136. For the sign-reversed inhibitor endpoint the corresponding range is −0.795163 to 0.957725. These are permutation sensitivity diagnostics, not six new experiments or evidence for any particular pairing. CTT insertion is also a different endpoint from our exact programmed sequence and native-local peptide contracts.

## Admission decision

**Admitted:** authentic functional values and condition-level descriptive comparison. **Not admitted:** paired sequence–function correlation, pool-level calibration of observation bounds, or original-fragment editing accuracy. Required crosswalk fields are donor, culture/pool, condition, collection time, sequencing library, functional assay and the source supporting each link. A verified join would enable a small, explicitly descriptive pilot; with three donors it would still not establish a general dose–response model.
