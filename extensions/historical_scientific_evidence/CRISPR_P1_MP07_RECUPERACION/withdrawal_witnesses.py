"""Materialize exact witnesses for the first observed withdrawal in each CFTR class."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import sys

from recompute_panel import solve
from recompute_universes import allele_orientations, canonical


def flank_bases(witness, lo, hi):
    columns = [i for i, base in enumerate(witness["aligned_reference"]) if base != "-"]
    return [witness["aligned_read"][columns[i]] for i in (lo-1, hi)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, required=True)
    p.add_argument("--ledger", type=Path, required=True)
    p.add_argument("--diagnosis", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    frontier = a.repository / "extensions/correspondence_frontier"
    sys.path.insert(0, str(frontier.resolve()))
    import optimal_lattice as lattice

    findings = json.loads(a.diagnosis.read_bytes())
    selected = {r["example"]: r for r in findings if r["condition"] == "PE treated"
                and r["target"] == "native" and r["label_at_27"] == "E" and r["outcome_at_36"] != "retained"}
    with gzip.open(a.ledger, "rt", encoding="ascii", newline="") as stream:
        raw = {row["read_id"]: row["raw_sequence"] for row in csv.DictReader(stream)
               if row["read_id"] in selected}
    panel = json.loads((frontier / "inputs/cftr_REAL_ALIGNMENT_CASES.json").read_bytes())["cases"]
    specs = {case["reference"]: case for case in panel}
    maps = {}
    output = []
    for read_id, selection in selected.items():
        accession = selection["accession"]
        case_results = {}
        for name, case in specs.items():
            if (name, accession) not in maps:
                table = a.repository / "raw_replay/inputs/alleles" / name / f"{accession}.zip"
                maps[name, accession], _ = allele_orientations(table, case["reference_sequence"])
            oriented, = maps[name, accession][canonical(raw[read_id])]
            optimum, targets = solve({**case, "observed_sequence": oriented}, lattice, cutoff=-1)
            witnesses = targets["native"]
            for witness in witnesses.values():
                witness.update(deficit=optimum-witness["score"], flank_bases=flank_bases(
                    witness, case["lo"], case["hi"]))
            assert witnesses[1]["deficit"] == 0
            case_results[name] = {"optimum": optimum, "witnesses": witnesses}
        opposite = any(0 in r["witnesses"] and r["witnesses"][0]["deficit"] == 36 for r in case_results.values())
        unknown = any(2 in r["witnesses"] and r["witnesses"][2]["deficit"] == 36 for r in case_results.values())
        assert (opposite, unknown) == {"contradiction_only": (True, False),
            "unobservable_only": (False, True), "both": (True, True)}[selection["outcome_at_36"]]
        output.append({"selected": selection, "references": case_results})
        print(read_id, selection["outcome_at_36"], {
            name: [(k, w["deficit"], w["flank_bases"]) for k, w in r["witnesses"].items()]
            for name, r in case_results.items()
        })
    a.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
