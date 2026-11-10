"""Rejoin complete per-reference ledgers and calculate exact finite-record contrasts."""

import argparse
from collections import Counter
import csv
from fractions import Fraction
import gzip
import json
from pathlib import Path

from recompute_panel import file_hash
from recompute_universes import DELTAS, label_at


def support(entry, target, delta):
    if entry["status"] != "evaluated":
        return "U"
    return label_at(entry["deficits"][target], delta)


def paired_rows(directory, left, right):
    files = [directory / f'{r["run"]}_{r["reference"]}.jsonl.gz' for r in (left, right)]
    for file, receipt in zip(files, (left, right), strict=True):
        assert file_hash(file) == receipt["records_sha256"]
    counts = Counter()
    marginal = [Counter(), Counter()]
    total = 0
    with gzip.open(files[0], "rt", encoding="ascii") as a, gzip.open(files[1], "rt", encoding="ascii") as b:
        for la, lb in zip(a, b, strict=True):
            x, y = json.loads(la), json.loads(lb)
            assert x["id"] == y["id"]
            total += 1
            for index, target in enumerate(left["targets"]):
                for delta in DELTAS:
                    lx, ly = support(x, index, delta), support(y, index, delta)
                    counts[target, delta, lx if lx == ly else "U"] += 1
                    marginal[0][target, delta, lx] += 1
                    marginal[1][target, delta, ly] += 1
    assert total == left["raw_N"] == right["raw_N"]
    for i, receipt in enumerate((left, right)):
        for r in receipt["rows"]:
            for label in ("E", "NE", "U"):
                assert marginal[i][r["target"], r["delta"], label] == r[label]
    return [{"cohort": left["cohort"], "run": left["run"], "reference": "both", "target": target, "delta": delta,
             "N": total, **{label: counts[target, delta, label] for label in ("E", "NE", "U")}}
            for target in left["targets"] for delta in DELTAS]


def bounds(row):
    assert row["E"] + row["NE"] + row["U"] == row["N"]
    return Fraction(row["E"], row["N"]), Fraction(row["E"] + row["U"], row["N"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--universes", type=Path, required=True)
    args = parser.parse_args()
    full = json.loads((args.universes / "full_receipt.json").read_bytes())
    assert full["status"] == "COMPLETE_DELTA_0_TO_54"
    grouped = {}
    for receipt in full["runs_references"]:
        grouped.setdefault(receipt["run"], []).append(receipt)
    rows = []
    for pair in grouped.values():
        left, right = pair
        assert left["targets"] == right["targets"]
        for r in pair:
            rows.extend({"cohort": r["cohort"], "run": r["run"], "reference": r["reference"],
                         "N": r["raw_N"], **e} for e in r["rows"])
        rows.extend(paired_rows(args.universes, left, right))
    for row in rows:
        low, high = bounds(row)
        row.update(lower=str(low), upper=str(high))
    indexed = {(r["run"], r["reference"], r["target"], r["delta"]): r for r in rows}
    contrasts = []
    runs = json.loads((args.repository / "raw_replay/metadata/RUNS.json").read_bytes())["runs"]
    for replica in (1, 2, 3):
        by_condition = {r["condition"]: r["run_accession"] for r in runs if r["source_replica"] == replica}
        for reference in ("mutant", "programmed", "both"):
            for target in ("native", "exact"):
                for delta in DELTAS:
                    treatment = indexed[by_condition["PE treated"], reference, target, delta]
                    mock = indexed[by_condition["mock"], reference, target, delta]
                    tlo, thi = bounds(treatment)
                    mlo, mhi = bounds(mock)
                    lo, hi = tlo-mhi, thi-mlo
                    contrasts.append({"preparation": replica, "reference": reference, "target": target,
                                      "delta": delta, "lower": str(lo), "upper": str(hi),
                                      "lower_pp": float(100*lo), "upper_pp": float(100*hi)})
    for name, data in (("FULL_EVENT_BOUNDS.csv", rows), ("CFTR_CONTRASTS.csv", contrasts)):
        with (args.universes / name).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    for row in contrasts:
        if row["reference"] == "both" and row["delta"] in (0, 9, 54):
            print(row)


if __name__ == "__main__":
    main()
