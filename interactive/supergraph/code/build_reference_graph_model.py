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
    "Evidence artifacts", "Decision tests", "Scope boundaries",
]
LENSES = [
    {"id": "geometry", "label": "Product and read geometry"},
    {"id": "correspondence", "label": "Sequence correspondence"},
    {"id": "identification", "label": "Unknowns and conditional bounds"},
    {"id": "sensitivity", "label": "Reference and alignment sensitivity"},
    {"id": "domain", "label": "Reference and target domains"},
    {"id": "reproducibility", "label": "Source identities and replay"},
    {"id": "scope", "label": "Scope and future tests"},
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
    parser.add_argument("--output", type=Path, default=component / "supergraph_model.json")
    args = parser.parse_args()
    model = json.loads(args.model.read_text())
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
    assert len(nodes) == 57
    graph = {
        "schema_version": "2.0-crispr-reference-map",
        "generated_date": model["generatedDate"],
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
    args.output.write_text(json.dumps(graph, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": "PASS", "model_sha256": model_sha,
                      "graph_sha256": digest(args.output), "nodes": len(nodes),
                      "edges": len(edges), "source_identities_verified": len(sources),
                      "witness_observation_partitions_verified": len(witness),
                      "scope": "Deterministic graph conversion and source identities; no new experimental result."}))


if __name__ == "__main__":
    main()
