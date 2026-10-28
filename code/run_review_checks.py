#!/usr/bin/env python3
"""Run derived arithmetic, then raw-observer/frozen-alignment reconstruction.

Both scopes use temporary working output. Retained reports go to --output-dir.
No source data, expected result, source manifest or engine is modified or run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time

OUTPUTS = ("M1_TARGET_COUNTS.json", "RAW_UNIVERSE_IDENTITIES.json", "WITNESS_BOUNDS.json")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=root / "results/review")
    args = parser.parse_args()
    destination = args.output_dir.resolve()
    for protected in (root / "raw_replay", root / "tables", root / "metadata", root / "code"):
        if destination == protected or destination.is_relative_to(protected):
            parser.error("Output must remain outside frozen data, code and metadata directories")
    started = time.monotonic()
    retained = []
    with tempfile.TemporaryDirectory(prefix="mp02-review-checks-") as scratch:
        scratch = Path(scratch)
        arithmetic = scratch / "arithmetic"
        process = subprocess.run([sys.executable, "-B", str(root / "code/run_checks.py"),
                                  "--output-dir", str(arithmetic)], check=False)
        if process.returncode:
            return process.returncode
        raw = scratch / "raw"
        process = subprocess.run([sys.executable, "-B", str(root / "raw_replay/run_raw_replay.py"),
                                  "--output", str(raw)], check=False)
        if process.returncode:
            return process.returncode
        matches = []
        for name in OUTPUTS:
            generated, expected = raw / name, root / "raw_replay/expected" / name
            equal = generated.read_bytes() == expected.read_bytes()
            matches.append({"output": name, "byte_identical": equal,
                            "sha256": digest(generated), "expected_sha256": digest(expected)})
        if not all(row["byte_identical"] for row in matches):
            print("Generated raw-replay output bytes differ from preserved comparators", file=sys.stderr)
            return 1
        arithmetic_result = json.loads((arithmetic / "MP02_CHECK_SUMMARY.json").read_text())
        raw_result = json.loads((raw / "VERIFICATION.json").read_text())
        destination.mkdir(parents=True, exist_ok=True)
        for folder, files in [(arithmetic, ("MP02_CHECK_SUMMARY.json", "SELECTED_P1_ARITHMETIC_CHECK.json")),
                              (raw, OUTPUTS + ("VERIFICATION.json",))]:
            label = "arithmetic" if folder == arithmetic else "raw_replay"
            target = destination / label
            target.mkdir(parents=True, exist_ok=True)
            for name in files:
                shutil.copyfile(folder / name, target / name)
                retained.append(label + "/" + name)
    elapsed = time.monotonic() - started
    result = {
        "status": "PASS", "version": "0.2.0", "edition": "MP02",
        "python": platform.python_version(), "runtime_seconds": elapsed,
        "under_five_minutes_on_this_execution": elapsed < 300,
        "working_output_was_temporary_and_removed": True,
        "scope_1": {"status": arithmetic_result["status"],
                    "scope": "Supplied-table/code identities and exact derived-table arithmetic"},
        "scope_2": {"status": raw_result["status"],
                    "scope": "Observers reconstructed from authentic raw reads with frozen alignment projections",
                    "raw_reads": raw_result["raw_reads"],
                    "witnessed_records": raw_result["witnessed_records"],
                    "unresolved_records": raw_result["unresolved_records"],
                    "M1_target_rows": raw_result["M1_target_rows"],
                    "aligned_by_reference": raw_result["aligned_by_reference"],
                    "original_comparator_outputs": matches},
        "retained_reports": retained,
        "not_claimed": {"fresh_external_alignment": False, "full_historical_pipeline": False,
                        "independent_human_or_biological_validation": False,
                        "molecular_accuracy_calibration": False,
                        "universal_five_minute_runtime": False},
    }
    path = destination / "REVIEW_CHECK_SUMMARY.json"
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "runtime_seconds": elapsed,
                      "completed_scopes": ["exact derived-table arithmetic", "raw observers with frozen alignments"],
                      "raw_records": raw_result["raw_reads"], "unknown_records": raw_result["unresolved_records"],
                      "raw_outputs_byte_identical": len(matches),
                      "external_alignment_engine_executed": False, "receipt": str(path)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
