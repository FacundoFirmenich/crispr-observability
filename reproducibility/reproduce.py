"""Offline reproduction of frozen public evidence. SPDX-License-Identifier: GPL-3.0-only."""

import argparse
import hashlib
import json
import platform
import sys
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
        report = {"integrity": integrity, "replay": load(result)}
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
