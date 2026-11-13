"""Build a typed, source-bound graph without reading manuscript prose.

The input is the previously verified numerical companion. This conversion adds
navigation objects and declared checks, never new observations or experiments.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json


LAYERS = [
    "Premises", "Observation mechanisms", "Quantitative constraints",
    "Evidence artifacts", "Research tools", "External validation", "Decision tests", "Scope boundaries",
]
LENSES = [
    {"id": "geometry", "label": "Product and read geometry"},
    {"id": "correspondence", "label": "Sequence correspondence"},
    {"id": "identification", "label": "Unknowns and conditional bounds"},
    {"id": "sensitivity", "label": "Reference and alignment sensitivity"},
    {"id": "domain", "label": "Reference and target domains"},
    {"id": "reproducibility", "label": "Source identities and replay"},
    {"id": "scope", "label": "Scope and future tests"},
    {"id": "transport", "label": "External panels and local outcomes"},
]
NODE_LENSES = {
    "E_geometry": ["geometry"], "E_anchor": ["geometry", "identification"],
    "E_strings": ["correspondence"], "E_counts": ["identification"],
    "E_contrasts": ["identification"], "E_grid": ["sensitivity"],
    "E_labels": ["sensitivity", "correspondence"], "E_external": ["domain"],
    "A_quality": ["geometry", "correspondence", "sensitivity", "identification"],
    "A_labels": ["identification", "scope"],
    "A_completion": ["identification"], "A_domain": ["domain", "scope"],
    "A_validation_scope": ["scope", "reproducibility"],
    "C_anchor": ["geometry", "identification"],
    "C_strings": ["correspondence"], "C_hull": ["identification"],
    "C_witness": ["correspondence", "identification"],
    "C_grid": ["sensitivity", "identification"], "C_domain": ["domain"],
    "B_accuracy": ["scope", "identification"],
    "B_transport": ["scope", "sensitivity"],
}
STATUS = {
    "observed_configuration": "supported", "replayed": "supported",
    "recomputed": "supported", "executed_retrospective": "supported",
    "audited": "supported", "declared": "constraint",
    "conditional": "conditional", "mathematical": "conditional",
    "supported_local": "supported", "conditional_exact": "conditional",
    "supported_conditional": "conditional", "supported_finite_grid": "conditional",
    "supported_catalog": "supported", "open": "blocked",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tidy(text):
    # Typographic spacing only; the numerical and scientific contracts are frozen.
    return text.replace("at225", "at 225").replace("at224", "at 224").replace("with163", "with 163").replace("across14", "across 14")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # Installed layout: interactive/supergraph/code/<this script>.
    component = Path(__file__).resolve().parents[1]
    parser.add_argument("--model", type=Path, default=component.parent / "evidence_model.json")
    parser.add_argument("--repository", type=Path, default=component.parent.parent)
    parser.add_argument("--output", type=Path, default=component / "review_graph_model.json")
    args = parser.parse_args()
    model = json.loads(args.model.read_bytes())
    model_sha = digest(args.model)
    sources = {row["id"]: row for row in model["sources"]}
    for row in sources.values():
        relative = Path(row["path"])
        assert not relative.is_absolute() and ".." not in relative.parts
        assert digest(args.repository / relative) == row["sha256"], row["path"]

    nodes, edges = [], []
    source_lenses = {key: {"reproducibility"} for key in sources}
    for original in model["graph"]["nodes"]:
        ident = original["id"]
        if original["type"] == "assumption":
            layer = "Premises"
        elif original["type"] == "boundary":
            layer = "Scope boundaries"
        elif original["type"] == "evidence":
            layer = "Evidence artifacts"
        elif ident in {"C_anchor", "C_strings"}:
            layer = "Observation mechanisms"
        else:
            layer = "Quantitative constraints"
        refs = [sources[key]["path"] for key in original["evidenceIds"]]
        for key in original["evidenceIds"]:
            source_lenses[key].update(NODE_LENSES[ident])
        nodes.append({
            "id": ident, "label": original["label"], "type": original["type"],
            "layer": layer, "status": STATUS[original["status"]],
            "status_original": original["status"],
            "summary": tidy(original["mechanism"]),
            "details": original["scope"], "scope": original["scope"],
            "falsification": original["falsification"],
            "artifacts": refs, "evidence_ids": original["evidenceIds"],
            "lenses": NODE_LENSES[ident], "metrics": {},
            "optional_detail": False,
        })

    by_id = {node["id"]: node for node in nodes}
    by_id["E_geometry"]["metrics"] = {
        item["label"]: f"[{item['start']}, {item['end']}); read end {item['read_end']}"
        for item in model["geometry"]
    }
    summary = model["summary"]
    by_id["E_strings"]["metrics"] = {
        "Different observed coding strings": summary["codingDisagreement"],
        "Equal labels but different coding strings": summary["equalLabelsDifferentCoding"],
    }
    by_id["E_counts"]["metrics"] = {
        "Original raw records": summary["rawRecords"],
        "Common-string witnesses": summary["witnessRetained"],
        "Unknown records retained in denominator": summary["unknown"],
    }
    by_id["E_grid"]["metrics"] = {
        "Completed historical alignment jobs": summary["alignmentJobs"],
        "Reference-by-OAT profiles": 14,
        "Separate stock profiles": 2,
    }
    by_id["E_labels"]["metrics"] = {
        "Baseline reference-dependent native-label changes": summary["referenceNativeLabelChanges"],
    }
    by_id["C_hull"]["metrics"] = {
        "Partition": model["calculationContract"]["partition"],
        "Finite-record hull": model["calculationContract"]["hull"],
        "Two-arm hull": model["calculationContract"]["contrast"],
    }
    for record in model["sensitivityEnvelopes"]:
        by_id["C_grid"]["metrics"][f"{record['preparation']} / {record['event']} / minimum lower endpoint (pp)"] = record["L_min_pp"]
    for case in model["externalCases"]:
        by_id["E_external"]["metrics"][f"{case['case']} / unresolved template assignment weight"] = case["unresolvedTemplateWeight"]

    def edge(source, target, relation, label, status, lenses, description, **extra):
        edges.append({
            "id": f"{source}->{target}:{relation}", "source": source, "target": target,
            "relation": relation, "label": label, "status": status,
            "lenses": lenses, "summary": description, **extra,
        })

    for hyperedge in model["graph"]["hyperedges"]:
        ident = hyperedge["id"]
        conclusion = by_id[hyperedge["conclusion"]]
        lenses = conclusion["lenses"]
        nodes.append({
            "id": ident, "label": "AND: " + conclusion["label"],
            "type": "conjunction", "layer": "Quantitative constraints",
            "status": "conditional", "summary": hyperedge["label"],
            "details": "Every listed premise is required for this declared support path. Withholding a premise changes displayed support, not recorded evidence or scientific truth.",
            "scope": hyperedge["scope"], "falsification": "Inspect the required premises and the stated conclusion; a navigation edge is not a causal effect or an experimental validation.",
            "artifacts": [], "lenses": lenses,
            "metrics": {"Required premises": len(hyperedge["premises"]), "Conclusion": conclusion["label"]},
            "required_premises": hyperedge["premises"], "conclusion": hyperedge["conclusion"],
            "optional_detail": False,
        })
        for premise in hyperedge["premises"]:
            edge(premise, ident, "requires", "required premise", "conditional", lenses,
                 "This premise is jointly required with the other premises of the AND gate.", conjunction=ident)
        edge(ident, hyperedge["conclusion"], "conditional_support", "joint support", "conditional", lenses,
             "The displayed conclusion is supported only within the complete declared conjunction and its scope.", conjunction=ident)

    for source in sources.values():
        ident = "S_" + source["id"]
        lenses = sorted(source_lenses[source["id"]])
        metrics = {"SHA-256": source["sha256"], "Source kind": source["kind"]}
        if "rows" in source:
            metrics["Source rows"] = source["rows"]
        nodes.append({
            "id": ident, "label": source["id"].replace("_", " "),
            "type": "artifact", "layer": "Evidence artifacts", "status": "supported",
            "status_original": "source_fingerprint_verified",
            "summary": "Frozen source artifact; its bytes match the recorded SHA-256 fingerprint.",
            "details": source["path"],
            "scope": "Byte identity and source attribution. This is not molecular truth calibration or independent human replication.",
            "falsification": "A byte mismatch changes the source identity and requires a new version and recheck.",
            "artifacts": [source["path"]], "lenses": lenses, "metrics": metrics,
            "optional_detail": True, "source_identity": source,
        })
        for original in model["graph"]["nodes"]:
            if source["id"] in original["evidenceIds"]:
                edge(ident, original["id"], "documents", "source artifact", "supported", NODE_LENSES[original["id"]],
                     "This frozen artifact documents the displayed node; the edge does not establish molecular accuracy.")

    for original in model["graph"]["nodes"]:
        if original["type"] != "claim":
            continue
        ident = "F_" + original["id"]
        claim = by_id[original["id"]]
        lenses = sorted(set(claim["lenses"] + ["scope"]))
        nodes.append({
            "id": ident, "label": "Test: " + claim["label"],
            "type": "falsifier", "layer": "Decision tests", "status": "open",
            "status_original": "declared_correction_condition_not_new_experiment",
            "summary": original["falsification"],
            "details": "Declared check or future evidence that would require correcting or narrowing the linked conclusion. Listing it does not mean a new experiment has been executed.",
            "scope": claim["scope"], "falsification": original["falsification"],
            "artifacts": claim["artifacts"], "lenses": lenses,
            "metrics": {"Linked conclusion": claim["label"], "Execution claim": "No new experiment is asserted"},
            "optional_detail": False,
        })
        edge(ident, claim["id"], "would_require_revision", "revision condition", "open", lenses,
             "If this stated condition is established, the linked local conclusion must be revised; this is not a reported experimental outcome.")

    validation_root = "extensions/generation_model/validation/results/"
    def read_result(name):
        return json.loads((args.repository / (validation_root + name)).read_bytes())

    def current_node(ident, title, summary, scope, artifacts, metrics, status="supported"):
        nodes.append({
            "id": ident, "label": title, "type": "evidence", "layer": "External validation",
            "status": status, "status_original": "public_source_bound_review_checkpoint",
            "summary": summary, "details": scope, "scope": scope, "artifacts": artifacts,
            "lenses": ["transport", "reproducibility", "scope"], "metrics": metrics,
            "optional_detail": False,
        })

    for case in ("MP37", "MP38"):
        filename = case + "_RESULT.json"
        result = read_result(filename)
        current_node("V_" + case, case + " held-out transfer",
                     "Transfer of the scoring rules to a reserved germline panel with local outcomes.",
                     "The selected loci and source conditions define the transfer test. The source retains unresolved calls and stratifies by event size and region.",
                     [validation_root + filename, validation_root + case + "_LOCAL_OBSERVATIONS.csv"],
                     {"Planned loci": result["planned_loci"],
                      "Zero-threshold consensus at 300 bases": result["arms"]["consensus0_300"]["outcomes"],
                      "Median rule at 300 bases": result["arms"]["median0_300"]["outcomes"]})
        edge("V_" + case, "B_transport", "tests", "Reserved external transport", "supported",
             ["transport", "scope"], "A reserved germline panel tests transfer at the locus level under its own domain.")

    hifi = read_result("HOMOPOLYMER_HIFI_CROSSCHECK.json")["donors"]
    current_node("V_HIFI", "HiFi and Nanopore at matched loci",
                 "A genotype benchmark distinguishes a stable readout from an accurate call across acquisitions.",
                 "This comparison concerns germline genotype at matched loci. The donor and locus remain distinct evidence units.",
                 [validation_root + "HOMOPOLYMER_HIFI_CROSSCHECK.json"],
                 {donor: value["summary"] for donor, value in hifi.items()})
    edge("V_HIFI", "B_accuracy", "limits", "Different accuracy estimand", "conditional",
         ["transport", "scope"], "Germline agreement does not authenticate editing in the original fragment.")

    paired_path = "interactive/supergraph/paired_library_summary.json"
    paired = json.loads((args.repository / paired_path).read_bytes())
    current_node("V_PAIRED", "Paired external library audit",
                 "All deposited mate identifiers link. Joint read lengths differ by library.",
                 "The unit is the amplified insert. Event coverage and original-molecule editing accuracy remain unestimated.",
                 [paired_path], paired["runs"])
    edge("V_PAIRED", "B_accuracy", "limits", "Insert linkage boundary", "conditional",
         ["transport", "scope"], "Mate linkage locates the same library insert and leaves original-molecule identity open.")

    current_node("V_CLONE", "COL17A1 retained-span replay",
                 "Clone3, boundary sensitivity and complete retained-span counts link to the public replay.",
                 "Clone-level and paired-fragment evidence answer different questions from independent truth on HEK3 or CFTR.",
                 ["extensions/retained_span_audit/README.md", "reproducibility/EVIDENCE.md"], {})
    edge("V_CLONE", "B_accuracy", "limits", "Clone and fragment scope", "conditional",
         ["transport", "scope"], "The COL17A1 retained spans support local correspondence, under the named admission rules.")

    current_node("V_FUNCTION", "CFTR functional context and linkage gate",
                 "Recovered Figure 4 source cells show condition-level channel responses. Across six possible sequence/function pairings, F&I Pearson correlation ranges from -0.998609 to 0.916136.",
                 "Three donor-level observations per condition. Replicate column order does not authenticate a shared culture, library or original fragment. Paired correlation remains unadmitted.",
                 ["examples/observation-report/functional_bridge.json", "docs/FUNCTIONAL_BRIDGE.md", "code/functional_bridge.py"],
                 {"Functional pairing": "UNAUTHENTICATED_ACROSS_ASSAYS", "Same pool": "NOT_ESTABLISHED"}, status="conditional")
    nodes[-1]["layer"] = "Evidence artifacts"
    nodes[-1]["type"] = "artifact"
    edge("V_FUNCTION", "B_accuracy", "limits", "Functional response is not fragment truth", "conditional",
         ["transport", "scope"], "Condition-level functional context does not independently label individual sequenced fragments.")

    for ident, title, description, artifacts, target in [
        ("T_FASTQ", "FASTQ to an exact-DNA decision", "Streaming exact-anchor observer, direct in-memory Streamlit FASTQ upload, raw denominator retention and HTML/JSON reports. Regression agrees with nine archived CFTR libraries; this is not complete alignment-domain certification.",
         ["code/fastq_decision.py", "code/test_fastq_decision.py", "examples/observation-report/cftr_fastq_config.json"], "T_REPORT"),
        ("T_TUTORIAL", "Seven-block CFTR learning route", "Bilingual conceptual tutorial with explicit observation units and interpretation checkpoints. Duration is a design target; independent usability testing remains pending.",
         ["interactive/CFTR_TUTORIAL.html", "docs/PROSPECTIVE_VALIDATION.md", "docs/COMPARATIVE_SCOPE.md"], "T_REPORT"),
        ("T_REPORT", "Four-step observation report", "Define, account, decide and investigate. A dependency-free Python entry point and guided Streamlit view emit exact JSON and readable HTML; local withdrawal costs remain visible.",
         ["code/observability_report.py", "docs/START_HERE.md", "examples/observation-report/report.json", "examples/observation-report/local_utility.json"], "T_CASE"),
        ("T_LAB", "Read Observability Lab", "Manipulate read boundaries and sequence spans in a source-bound 3D laboratory; import cases and export replayable results.",
         ["interactive/lab/README.md", "interactive/lab/CFTR_OBSERVABILITY_LAB.html"], "E_geometry"),
        ("T_CASE", "Executable CFTR leading case", "Two authentic reads, seven operations, exact expected outcomes and explicit geometric counterfactuals.",
         ["examples/cftr-leading-case/project.json", "examples/cftr-leading-case/PROVENANCE.json", "docs/LEADING_CASE.md"], "E_strings"),
        ("T_ENGINE", "Multi-input CLI, API and MCP", "Named studies share one computational kernel. Six operation types include exact same-record junction evidence; researcher documentation is bilingual.",
         ["engine/README.md", "docs/RESEARCH_LINES.md", "docs/RESEARCHER_GUIDE.md"], "T_LAB"),
    ]:
        nodes.append({"id": ident, "label": title, "type": "artifact", "layer": "Research tools",
                      "status": "supported", "status_original": "implemented_research_interface",
                      "summary": description, "details": description,
                      "scope": "Implemented access and computation. Tool availability is not independent validation or community adoption.",
                      "falsification": "Replay the named operation and compare its report with the source contract.",
                      "artifacts": artifacts, "lenses": ["reproducibility", "geometry", "correspondence"],
                      "metrics": {}, "optional_detail": False})
        edge(ident, target, "documents", "Inspect and replay", "supported", ["reproducibility"],
             "This tool exposes or recomputes the linked object; it does not add an independent biological observation.")

    comparison_root = "extensions/comparative_execution/results/"
    policy = json.loads((args.repository / (comparison_root + "LOCAL_POLICY.json")).read_bytes())
    current_node("V_LOCAL_POLICY", "Donor-held-out local decision routing",
                 "At discordance cost six, local loss is 137 versus 164 for baseline, 152 for always-majority and 165 for pooled routing.",
                 "All five costs and eight donor/tract cells are retained. Local routing is worse than baseline at cost four and worse than always-majority at ten and twenty. This is retrospective genotype transport, not fragment-editing accuracy.",
                 [comparison_root + "LOCAL_POLICY.json", "extensions/comparative_execution/local_policy.py"],
                 {"N": policy["N"], "Cost grid": policy["cost_grid"], "Cost six C/D/U": [124,4,113]})
    edge("V_LOCAL_POLICY", "B_transport", "tests", "Hold out each donor", "supported", ["transport","scope"], "Training excludes the evaluated donor; shared loci and relatedness remain dependence considerations.")
    current_node("V_COST_MAP", "Exact applicability map of local routing",
                 "The empirical local policy beats both fixed policies and pooled routing precisely for 17/3 < cost ratio <= 20/3 on this panel.",
                 "Exploratory enumeration after the frozen cost grid; 24 open cells and 23 boundary points preserve ties and adverse regimes. This is not prospective dominance.",
                 [comparison_root + "COST_FRONTIER.json", "extensions/comparative_execution/cost_frontier.py"], {"Open cells":24,"Boundary points":23})
    edge("V_COST_MAP", "V_LOCAL_POLICY", "limits", "Bound the favorable jurisdiction", "supported", ["transport","sensitivity"], "The complete decision-cost map prevents selection of one favorable cost from becoming a universal claim.")
    current_node("V_NATIVE_ANCHOR", "Executed native assignment-anchor intervention",
                 "On 175010 Q30 records, ampliCan admission is 112822, 174439 and 174510 at anchor offsets zero, four and eight.",
                 "All nine libraries use the same native version and input. Relocation also changes the primer-excluded event boundary. The offset-four to offset-eight net gain of 71 comprises 209 admissions and 138 exclusions. More admission is not demonstrated accuracy; native HDR and generic edited categories remain distinct.",
                 [comparison_root + "AMP_ANCHOR_SUMMARY.json", "extensions/comparative_execution/ANCHOR_PROBE_ADDENDUM.md"], {"Q30 input":175010,"Assigned":[112822,174439,174510]})
    edge("V_NATIVE_ANCHOR", "E_anchor", "tests", "Intervene on an observation boundary", "supported", ["geometry","sensitivity"], "This new native-engine sensitivity complements the finite-read downstream-anchor mechanism; it is not the same experiment.")

    current_node("V_NATIVE_RETURN", "Native return-contract intervention",
                 "Three treated libraries return multiple counts before the event type. The first counts are 2, 4 and 1; linked-ID totals are 897, 1398 and 882.",
                 "The original objects and native outputs are preserved. Scalar-interface validity is distinct from biological truth. All nine typed-return reruns pass, and the six control tables are unchanged. The intervention changes aggregation, not aligned data or matching criteria.",
                 [comparison_root + "CRISPRA_INTERFACE_RESULT.json", "extensions/comparative_execution/INTERFACE_REPAIR_ADDENDUM.md"], {"Affected treated libraries":3,"Controls with scalar interface":6})
    edge("V_NATIVE_RETURN", "T_REPORT", "tests", "Check the event contract at an interface", "supported", ["reproducibility","scope"], "Execution success does not guarantee that count, type and identifiers retain their declared meanings.")

    witness = [row for row in model["observations"] if row["operator"] == "triple_observed_coding_witness"]
    witness_contrasts = [row for row in model["contrasts"] if row["operator"] == "triple_observed_coding_witness"]
    for row in witness:
        assert row["E"] + row["NE"] + row["U"] == row["N"]
        assert Fraction(row["L_exact"]) == Fraction(row["E"], row["N"])
        assert Fraction(row["H_exact"]) == Fraction(row["E"] + row["U"], row["N"])
    assert len({node["id"] for node in nodes}) == len(nodes)
    ids = {node["id"] for node in nodes}
    assert all(row["source"] in ids and row["target"] in ids for row in edges)
    assert len({row["id"] for row in edges}) == len(edges)

    graph = {
        "schema_version": "2.1-crispr-reference-map",
        "generated_date": model["generatedDate"],
        "public_evidence_revision": "2026-11-13",
        "title": "CRISPR observability: interactive scientific map",
        "subtitle": "Premises, mechanisms, conditional constraints, source artifacts, revision tests and scope",
        "summary": {
            "source_fingerprints_verified": len(sources),
            "validation": f"{len(sources)}/{len(sources)} source identities",
            "sha256_prefix": model_sha[:16], "evidence_model_sha256": model_sha,
            "raw_records": summary["rawRecords"], "witness_records": summary["witnessRetained"],
            "unknown_records": summary["unknown"], "historical_alignment_jobs": summary["alignmentJobs"],
            "core_nodes": len(model["graph"]["nodes"]), "AND_gates": len(model["graph"]["hyperedges"]),
            "source_artifact_nodes": len(sources), "declared_revision_test_nodes": 6,
            "total_nodes": len(nodes), "total_edges": len(edges),
        },
        "map_semantics": "Typed navigation of frozen evidence and declared conditional support. AND gates preserve joint premises; ordinary source and revision-condition edges are not causal effects. Node positions, counts and filters do not quantify evidential strength.",
        "source_scope": model["scope"],
        "calculation_contract": model["calculationContract"],
        "paper_refs": {
            "fails": "Local obstructions and recorded disagreement; no universal observer ranking.",
            "candidate": "Declared future checks; no new biological experiment is asserted.",
            "open": "Correction conditions and future evidence requirements.",
            "blocked": "Molecular truth calibration and prospective transport remain outside this released evidence.",
            "supported": "Frozen observations, numerical artifacts and verified source identities.",
            "all": "Current public numerical companion. The full manuscript and formal supplement remain private.",
        },
        "lenses": LENSES, "layers": LAYERS,
        "target_domains": [dict(row) for row in model["catalogs"]["events"]],
        "preparations": model["catalogs"]["preparations"],
        "context_metrics": {"witness_observations": witness, "witness_contrasts": witness_contrasts},
        "nodes": nodes, "edges": edges, "hyperedges": model["graph"]["hyperedges"],
        "sources": model["sources"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(graph, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "PASS", "model_sha256": model_sha,
                      "graph_sha256": digest(args.output), "nodes": len(nodes),
                      "edges": len(edges), "source_identities_verified": len(sources),
                      "witness_observation_partitions_verified": len(witness),
                      "scope": "Deterministic graph conversion and source identities; no new experimental result."}))


if __name__ == "__main__":
    main()
