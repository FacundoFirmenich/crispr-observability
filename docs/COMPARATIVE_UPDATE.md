# From observation contracts to local decisions

The 13 November 2026 comparative extension adds fresh native executions, record-level admission transitions and a complete finite-panel decision-cost map. It preserves the preceding releases and every adverse regime. The analyses are retrospective computational experiments on authenticated public data.

## What changes for a researcher?

**A stable percentage can conceal a changed population.** On the same 175,010 Q30 records from nine CFTR libraries, ampliCan 1.32.1 assigns 112,822 records with the declared reference-[0,20) computational anchor, 174,439 at [4,24), and 174,510 at [8,28). These are computational assignment anchors, not recovered experimental PCR-primer sequences. Relocation also changes the native primer-excluded event boundary. It therefore changes both admission and the observation domain.

**A small net change can conceal opposing transitions.** The second relocation newly admits 209 records and excludes 138, for a net gain of 71. Full observed sequence families and unchanged multiplicities identify these transitions within each library. This is not a new molecular identity claim. A researcher can inspect which evidence changed instead of interpreting the gain as a uniformly improved observation set.

**A local policy can help in one cost regime and harm in another.** We hold out each of four donors and train only on the other donors' two tract strata. At discordance cost six and abstention cost one, local routing retains 124 concordant and four discordant calls and leaves 113 unresolved. Loss is 137, versus 164 for baseline, 152 for always applying the fixed abstention action, and 165 for pooled routing. The same local policy is adverse at other costs. Exact rational enumeration finds its favorable region against both fixed policies and pooled routing: 17/3 < cost ratio <= 20/3, on this finite panel only.

These results connect observation geometry, evidence identity and selective decision-making. They do not establish a universal accuracy ranking or a population-optimal policy. Independent original-fragment editing accuracy remains outside this release's completion gates.

## Evidence and reproduction

| Layer | Evidence | Reproduce or inspect |
|---|---|---|
| Fresh CRISPResso2 | Nine native 2.3.4 jobs; all acquired-record projections authenticated | `extensions/comparative_execution/results/native_projection/` |
| ampliCan admission | Nine libraries, three anchors, 27 within-library transitions | `results/AMP_ANCHOR_SUMMARY.json`; `collect_native.py --output /persistent/new.json` |
| Donor-held-out decisions | 241 observations, eight donor/tract cells, all five frozen costs | `local_policy.py --output /persistent/new-policy.json` |
| Exact applicability | 24 open cost cells and 23 boundary points | `cost_frontier.py --output /persistent/new-frontier.json` |
| Integrity checks | Recompute source-linked decisions and verify transitions | `python extensions/comparative_execution/check_results.py` |
| CRISPR-A | Nine native variant stages plus nine typed-return reruns; linked counts reconcile and six control tables remain identical | `NATIVE_EXECUTION_ADJUDICATION.md` |

The command paths in the middle rows are relative to `extensions/comparative_execution/`. Output paths must be new. Native archives remain in the [comparative evidence release](https://github.com/FacundoFirmenich/crispr-observability/releases/tag/comparison-evidence-43493655364). The collector reads the ampliCan archives in memory. Native engine images are run on hosted Linux workers, not installed on the development laptop. Runtime observations are stage-specific and do not rank tool speed across different operations.

## Start with your own files

Open the Streamlit reviewer application and expand **Analyse your FASTQ**. Upload treated and control files, inspect the editable exact-anchor configuration, and generate the HTML/JSON report. The host processes uploads in memory; a remote host receives those files. The browser route rejects more than 20 MB or 250,000 records per file rather than truncating data. Use the streaming CLI for larger inputs. The exact-anchor observer and complete alignment-domain certification have separate contracts.

The Supergraph explains the new applicability and admission nodes. Numerical Graph and PaperGraph retain their historical evidence identities; their current companion documents explain the broader extension. The tool is an inspectable analysis environment, not a replacement for deciding what event an experiment should measure.
