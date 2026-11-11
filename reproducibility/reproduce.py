"""Offline reproduction of frozen public evidence. SPDX-License-Identifier: GPL-3.0-only."""

import argparse
import csv
import hashlib
import json
import platform
import sys
from fractions import Fraction
from pathlib import Path


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def load(path):
    return json.loads(path.read_bytes())


def verify_files(root, manifest):
    seen = set()
    for entry in manifest["files"]:
        relative = entry["path"]
        path = (root / relative).resolve()
        require(path.is_relative_to(root) and relative not in seen, relative)
        require(path.is_file() and path.stat().st_size == entry["size"], relative)
        require(sha256(path) == entry["sha256"], relative)
        seen.add(relative)
    return len(seen)


def verify_validation(root):
    sys.path.insert(0, str(root / "campaign" / "homopolymer_transfer"))
    try:
        from verify_replication import verify
    finally:
        sys.path.pop(0)
    manifest = load(root / "MANIFEST.json")
    checked = verify_files(root, manifest)
    results = root / "results"
    legacy_graph = (results / "MP36_EPISTEMIC_GRAPH.json").read_bytes().decode("cp1252")
    require(json.loads(legacy_graph) == load(root / "MP36_EPISTEMIC_GRAPH.utf8.json"),
            "UTF-8 graph view changed the archived graph")
    older_panel = results / "HOMOPOLYMER_FROZEN_PANEL.json"
    newer_panel = results / "HOMOPOLYMER_HG003_HG004_FROZEN_PANEL.json"
    earlier = load(results / "HOMOPOLYMER_HG006_HG007_RESULT.json")
    later = load(results / "HOMOPOLYMER_HG003_HG004_RESULT.json")
    hifi = load(results / "HOMOPOLYMER_HIFI_CROSSCHECK.json")
    runner = root / "campaign" / "homopolymer_transfer" / "run_panel.py"

    require(earlier["binding"]["panel"] == sha256(older_panel), "HG006/HG007 panel identity")
    require(later["binding"]["panel"] == hifi["panel_sha256"] == sha256(newer_panel),
            "HG003/HG004 shared panel identity")
    require(load(newer_panel)["excluded_panel_sha256"] == sha256(older_panel),
            "Panel exclusion identity")
    require(later["binding"]["runner"] == hifi["runner_sha256"] == sha256(runner),
            "Scoring source identity")
    for field, name in (("inventory", "HOMOPOLYMER_PUBLIC_SOURCE_INVENTORY.json"),
                        ("aliases", "HOMOPOLYMER_ALIAS_READBACK.json")):
        require(later["binding"][field] == sha256(results / name), f"Source identity: {field}")

    before, old_positions = verify(earlier["donors"], load(older_panel))
    after, new_positions = verify(later["donors"], load(newer_panel))
    paired, hifi_positions = verify(hifi["donors"], load(newer_panel))
    require(not old_positions & new_positions and new_positions == hifi_positions,
            "Panel loci and cross-acquisition pairing")
    comparison = {"total_distinct_loci": len(old_positions | new_positions), "donors": before | after}
    require(comparison == load(results / "HOMOPOLYMER_REPLICATION_COMPARISON.json"),
            "Reconstructed versus frozen comparison")
    initial = {state: sum(s["initial"][state] for d in comparison["donors"].values()
                          for s in d["strata"].values()) for state in "CDU"}
    withdrawn = {state: sum(s["withdrawn"][state] for d in comparison["donors"].values()
                            for s in d["strata"].values()) for state in "CD"}
    require(initial == {"C": 207, "D": 26, "U": 8}, "Four-donor starting decisions")
    require(withdrawn == {"C": 113, "D": 25}, "Cost of the exact-length majority gate")

    hifi_counts = {state: sum(s["initial"][state] for d in paired.values()
                              for s in d["strata"].values()) for state in "CDU"}
    require(hifi_counts == {"C": 113, "D": 0, "U": 0}, "HiFi finite panel")
    retained = [row for donor in later["donors"].values() for rows in donor["strata"].values()
                for row in rows if row["exact_majority"] not in ("U", row["expected"])]
    require([row["pos"] for row in retained] == [52989163], "Persistent ONT discordance")
    return {
        "status": "PASS",
        "scope": "Reconstructed decisions from archived scores and flank calls at fixed loci",
        "verified_files": checked,
        "ont": comparison,
        "ont_initial": initial,
        "ont_withdrawn": withdrawn,
        "hifi_by_donor_and_stratum": paired,
        "hifi_initial": hifi_counts,
        "persistent_ont_discordance": {"donor": "HG004", "contig": "chr20", "position_1_based": retained[0]["pos"] + 1},
        "raw_alignment_scoring": "NOT_EXECUTED",
    }



def verify_cftr_bounds(root):
    """Reconstruct exact contrasts from the archived per-arm counts."""
    tables = root / "tables"

    def rows(name):
        with (tables / name).open(encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream))

    def index_counts(name, dimensions):
        indexed = {}
        for row in rows(name):
            key = tuple(row[k] for k in dimensions)
            require(key not in indexed, f"Repeated CFTR row: {key}")
            n, e, ne, u = (int(row[k]) for k in ("N", "E", "NE", "U"))
            require(n > 0 and min(e, ne, u) >= 0 and e + ne + u == n, key)
            interval = Fraction(e, n), Fraction(e + u, n)
            require(interval == (Fraction(row["L_exact"]), Fraction(row["H_exact"])), key)
            indexed[key] = interval
        return indexed

    def contrast(counts, treated, control, context, row):
        tlo, thi = counts[(treated, *context)]
        clo, chi = counts[(control, *context)]
        interval = tlo - chi, thi - clo
        require(interval == (Fraction(row["L_exact"]), Fraction(row["H_exact"])),
                f"Exact CFTR contrast: {treated}, {control}, {context}")
        return interval

    dims = ("trim", "policy", "observer", "event")
    geometry = index_counts("Table_S2_M1_All_Counts.csv", ("run_accession", *dims))
    contrasts = rows("Table_S3_M1_All_Contrasts.csv")
    primary = {}
    for row in contrasts:
        interval = contrast(geometry, row["treated"], row["control"],
                            tuple(row[k] for k in dims), row)
        if (row["prep"], row["control_condition"], row["trim"], row["policy"], row["event"]) == (
                "CFB-26-16", "mock", "0", "common_and_own_meanQ30",
                "native_local_peptide_all_observed_sequences"):
            primary[row["observer"]] = interval
    require(primary["broad_edge"][0] < 0 < primary["broad_edge"][1]
            and primary["target"][0] > 0, "Primary CFTR positivity contrast")

    witness = index_counts("Table_S7_Coding_Witness_Counts.csv",
                           ("run_accession", "observer", "event"))
    witness_contrasts = rows("Table_S8_Coding_Witness_Contrasts.csv")
    for row in witness_contrasts:
        contrast(witness, row["treated"], row["comparator"], (row["observer"], row["event"]), row)

    projections = index_counts("Table_S5_Projection_Counts.csv",
                               ("run_accession", "observer", "event"))
    design = {(row["preparation"], row["condition"]): row["run_accession"]
              for row in rows("Table_S1_Primary_Run_Manifest.csv")}
    table_one = rows("Table_1_Native_Contrast_Intervals.csv")
    reported = {}
    for row in table_one:
        prep, operator = row["preparation"], row["operator"]
        treated, control = design[(prep, "PE treated")], design[(prep, "mock")]
        if operator == "witness":
            interval = contrast(witness, treated, control,
                                ("triple_observed_coding_witness",
                                 "native_local_peptide_on_witnessed_observed_coding"), row)
        else:
            interval = contrast(projections, treated, control,
                                (operator, "native_local_peptide_all_observed_sequences"), row)
        if prep == "CFB-26-16":
            reported[operator] = interval
    require(reported["broad_edge"][0] < 0 < reported["broad_edge"][1]
            and reported["target"][0] > 0, "Table 1 positivity contrast")

    def endpoints(intervals):
        return {observer: {"exact": [str(v) for v in values],
                           "percentage_points": [round(100 * float(v), 5) for v in values]}
                for observer, values in sorted(intervals.items())}

    return {
        "status": "PASS",
        "scope": "Exact rational reconstruction from archived finite-record counts",
        "geometry_count_rows": len(geometry),
        "geometry_contrasts": len(contrasts),
        "witness_count_rows": len(witness),
        "witness_contrasts": len(witness_contrasts),
        "projection_count_rows": len(projections),
        "table_1_contrasts": len(table_one),
        "table_1_CFB_26_16_PE_minus_mock": endpoints(reported),
        "M1_CFB_26_16_common_and_own_meanQ30": endpoints(primary),
    }


def reproduce_retained(root, output, mode):
    sys.path.insert(0, str(root))
    original_argv = sys.argv
    try:
        from verify_release import validate
        from replay_col17_mp21 import main as replay

        integrity = validate(root)
        result = output / "COL17_REPLAY.json"
        sys.argv = ["replay_col17_mp21.py", "--mode", mode, "--output", str(result)]
        replay()
        report = {"integrity": integrity, "cftr_bounds": verify_cftr_bounds(root), "replay": load(result)}
        if mode == "full":
            import verify_alignment_oracle_mp21 as oracle
            oracle.OUT = output
            oracle.main()
            report["oracle"] = load(output / "INDEPENDENT_ORACLE.json")
        return report
    finally:
        sys.argv = original_argv
        sys.path.pop(0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("quick", "counts", "full"), default="quick")
    parser.add_argument("--output", type=Path, default=Path("crispr-reproduction-output"),
                        help="New directory for the evidence receipt and computation outputs")
    args = parser.parse_args()
    require(__debug__, "Run with Python assertions enabled")
    frozen = getattr(sys, "frozen", False)
    root = (Path(sys.executable).resolve().parent / "evidence" if frozen
            else Path(__file__).resolve().parents[1])
    out = args.output.resolve()
    require(not out.is_relative_to(root / "extensions"), "Keep outputs outside the evidence directories")
    validation = verify_validation(root / "extensions" / "generation_model" / "validation")
    out.mkdir(parents=True, exist_ok=False)
    retained = reproduce_retained(root / "extensions" / "retained_span_audit", out,
                                  "example" if args.mode == "quick" else args.mode)
    receipt = {
        "status": "PASS",
        "mode": args.mode,
        "environment": {
            "os": platform.system(), "architecture": platform.machine(),
            "python": platform.python_version(), "frozen": bool(frozen)
        },
        "retained_spans": retained,
        "germline_validation": validation,
        "original_fragment_editing_accuracy": "NOT_ESTIMATED",
        "original_fastq_extraction": "NOT_EXECUTED",
    }
    (out / "REPRODUCTION.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"PASS {args.mode} | {validation['ont']['total_distinct_loci']} ONT loci | "
          f"{validation['hifi_initial']['C']} HiFi concordances | {out / 'REPRODUCTION.json'}")


if __name__ == "__main__":
    main()
