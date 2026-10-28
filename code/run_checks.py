#!/usr/bin/env python3
"""Run supplied-byte and inherited arithmetic checks; no raw processing."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=root / "results")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((root / "metadata/DERIVED_INPUT_MANIFEST.json").read_text())
    records = []
    for row in manifest["files"] + [manifest["inherited_checker"]]:
        path = root / row["path"]
        data = path.read_bytes()
        actual = hashlib.sha256(data).hexdigest()
        records.append({"path": row["path"], "bytes": len(data), "sha256": actual,
                        "matches_frozen_identity": actual == row["sha256"] and len(data) == row["bytes"]})
    identities_ok = all(r["matches_frozen_identity"] for r in records)
    report = args.output_dir / "SELECTED_P1_ARITHMETIC_CHECK.json"
    command = [sys.executable, "-B", str(root / "code/check_selected_p1_tables.py"),
               "--tables-dir", str(root / "tables"), "--output", str(report)]
    if identities_ok:
        process = subprocess.run(command, check=False)
        arithmetic_ok = process.returncode == 0
    else:
        arithmetic_ok = False
    result = {
        "status": "PASS" if identities_ok and arithmetic_ok else "FAIL",
        "version": "0.2.0", "edition": "MP02", "python": platform.python_version(),
        "scope": "Supplied-byte identity and exact arithmetic of selected derived tables only.",
        "identity_checks": records, "arithmetic_checker_passed": arithmetic_ok,
        "arithmetic_report": report.name if identities_ok else None,
        "not_claimed": {"FASTQ_pipeline": False, "observation_replay": False,
                        "fresh_alignment": False, "independent_label_truth": False,
                        "new_biological_validation": False},
    }
    target = args.output_dir / "MP02_CHECK_SUMMARY.json"
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "scope": result["scope"],
                      "files_checked": len(records), "summary": str(target)}))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
