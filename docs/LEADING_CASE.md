# Leading case: when one missing base changes what can be observed

[Castellano](LEADING_CASE.es.md) · [Researcher guide](RESEARCHER_GUIDE.md)

This worked example takes the reader from a physical read boundary to a
sequence-correspondence failure, then to an auditable computational intervention.
It uses two actual archived CFTR reads and the published analysis's anchor
coordinates. All input identities are retained in
[`PROVENANCE.json`](../examples/cftr-leading-case/PROVENANCE.json).

## 1. See the mechanism

Open [Read Observability Lab](../interactive/lab/CFTR_OBSERVABILITY_LAB.html) in a
WebGL 2 browser. Select **Read boundary**. With anchor [210,222), shift +3 and
read end 224, the shifted interval is [213,225): only 11 of 12 anchor positions
are observed. Change the shift to zero; containment returns. Restore +3 and
change the read end to 225; containment also returns.

The last change is a counterfactual geometry calculation. It does not invent a
measured 225th base or demonstrate its sequence or quality. In the full recorded
CFTR analysis, changing the observation rule under a common quality mask changes
one lower contrast endpoint from −10.09177 to +13.23196 percentage points. That
empirical contrast is a separate frozen result, not an efficiency predicted by
the 3D scene.

## 2. Inspect the actual bases

In **Competing spans**, choose `SRR29623270.2091`. Compare [114,156) with
[114,201): 42 versus 87 bases, with 45 bases exclusive to the longer span.
These strings differ. This read has different native-local labels under the two
rules; it is not the equal-label counterexample.

Choose `SRR29623271.868` for that second question. The [114,186) and [114,153)
projections contain 72 and 39 bases while their archived native-local labels
agree (both false). An alternative optimal projection [114,186) restores string
agreement with the first rule. This shows why inspecting one traceback is
different from certifying every candidate in a complete declared domain.

The **Atomic model** view allows rotation, selection and measurement of ideal
geometry associated with the declared sequence. Imported experimental coordinates
have their own provenance. Neither representation establishes that the imaged
conformation occurred in the sequenced material.

## 3. Execute and preserve the result

From the repository root with Node.js 22 or later:

```sh
npm ci --prefix engine --no-fund --no-audit
node engine/project-cli.mjs check examples/cftr-leading-case/project.json
node engine/project-cli.mjs run examples/cftr-leading-case/project.json > leading-case-result.json
node examples/cftr-leading-case/verify.mjs
```

Choose a persistent output directory appropriate to your workstation. Installation
requires network access; the demonstration itself runs offline. `check` validates
the manifest contract. `run` reads the named files, checks their SHA-256 values and
computes seven operations. The final command verifies the expected outcomes.

| Operation | Expected outcome | Interpretation |
| --- | --- | --- |
| Unshifted anchor | 12/12 positions; contained | Original geometric condition |
| +3 shift, read end 224 | 11/12; one outside | Product-dependent observation loss |
| +3 shift, read end 225 | 12/12; contained | Explicit geometric counterfactual |
| First authentic read | 42 versus 87; unequal | Different observed strings |
| Second authentic read | 72 versus 39; unequal | Label agreement does not identify a string |
| Alternative optimum | Identical 72-base strings | Candidate choice matters; domain certification is a separate operation |

## 4. Connect the result to its claim

Open the [Supergraph](../interactive/supergraph/crispr_supergraph.html). Follow the
conceptual tour from observation geometry through correspondence and conditional
bounds. Inspect the input and premise links for each claim. Use the numerical
epistemic explorer for exact recorded counts. A graph link expresses a declared
dependency; the numerical computation and source record supply its evidence.

## 5. Test a new idea responsibly

Copy the manifest into a new study directory, keep the original unchanged, and
alter one named operation. Record whether the change is a measurement, an analysis
choice or a counterfactual. Compare outputs within the same coordinate system;
preserve unresolved outcomes and the original denominator. Export JSON/replay
alongside any rendered PNG. A picture alone cannot reproduce an intervention.

The example demonstrates observation geometry, local correspondence, provenance
and replay. It does not independently measure original-molecule editing accuracy.
That question requires an independently linked observation at the relevant unit.
