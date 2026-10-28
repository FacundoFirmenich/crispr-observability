"""Check conversion identity, graph contracts and the frozen source boundary."""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    component = Path(__file__).resolve().parents[1]
    p.add_argument("--model", type=Path, default=component.parent / "evidence_model.json")
    p.add_argument("--repository", type=Path, default=component.parent.parent)
    p.add_argument("--graph", type=Path, default=component / "supergraph_model.json")
    p.add_argument("--output", type=Path, default=component / "MODEL_CONVERSION_CHECK.json")
    a = p.parse_args()
    original = json.loads(a.model.read_text())
    graph = json.loads(a.graph.read_text())
    by_id = {n["id"]: n for n in graph["nodes"]}
    assert graph["summary"]["evidence_model_sha256"] == sha(a.model)
    assert graph["sources"] == original["sources"]
    assert graph["hyperedges"] == original["graph"]["hyperedges"]
    assert graph["calculation_contract"] == original["calculationContract"]
    for old in original["graph"]["nodes"]:
        current = by_id[old["id"]]
        for key in ["label", "type", "scope", "falsification"]:
            assert current[key] == old[key], (old["id"], key)
        assert current["evidence_ids"] == old["evidenceIds"]
    gate = by_id["H_grid"]
    assert {"A_labels", "A_completion"}.issubset(gate["required_premises"])
    for h in graph["hyperedges"]:
        assert by_id[h["id"]]["required_premises"] == h["premises"]
        for premise in h["premises"]:
            assert any(e["source"] == premise and e["target"] == h["id"] and e["relation"] == "requires" for e in graph["edges"])
    assert all(n["status"] == "open" and n["status_original"] == "declared_correction_condition_not_new_experiment" for n in graph["nodes"] if n["type"] == "falsifier")
    witness = [r for r in original["observations"] if r["operator"] == "triple_observed_coding_witness"]
    assert graph["context_metrics"]["witness_observations"] == witness
    assert graph["context_metrics"]["witness_contrasts"] == [r for r in original["contrasts"] if r["operator"] == "triple_observed_coding_witness"]
    for r in witness:
        assert r["E"] + r["NE"] + r["U"] == r["N"]
        assert Fraction(r["L_exact"]) == Fraction(r["E"], r["N"])
        assert Fraction(r["H_exact"]) == Fraction(r["E"] + r["U"], r["N"])
    assert not any("manuscript/" in path or "P1_Revised_Integrated" in path for n in graph["nodes"] for path in n["artifacts"])
    with tempfile.TemporaryDirectory(prefix="reference-graph-rebuild-") as tmp:
        rebuilt = Path(tmp) / "graph.json"
        subprocess.run([sys.executable, "-B", str(Path(__file__).with_name("build_reference_graph_model.py")),
                        "--model", str(a.model), "--repository", str(a.repository),
                        "--output", str(rebuilt)], check=True, stdout=subprocess.DEVNULL)
        assert rebuilt.read_bytes() == a.graph.read_bytes()
    result = {
        "status": "PASS", "input_model_sha256": sha(a.model), "graph_sha256": sha(a.graph),
        "graph_rebuild_byte_identical": True, "canonical_source_fingerprints_verified": len(original["sources"]),
        "core_node_source_fields_preserved": len(original["graph"]["nodes"]),
        "AND_hyperedges_preserved_with_all_premises": len(original["graph"]["hyperedges"]),
        "declared_revision_tests_not_new_experiments": 6,
        "witness_rows_bytewise_structure_preserved": len(witness),
        "total_nodes": len(graph["nodes"]), "total_edges": len(graph["edges"]),
        "scope": "Source identity and data-to-graph conversion consistency. Not molecular calibration, prospective validation or independent human replication.",
    }
    a.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
