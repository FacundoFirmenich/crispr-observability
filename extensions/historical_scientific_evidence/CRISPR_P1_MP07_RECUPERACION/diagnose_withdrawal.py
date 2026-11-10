"""Separate label contradiction from loss of projection at the CFTR transition."""

import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path

from recompute_universes import label_at


def classify(entry, target, delta):
    if entry["status"] != "evaluated":
        return "U"
    return label_at(entry["deficits"][target], delta)


def compare(first, second, target):
    a, b = (classify(record, target, 27) for record in (first, second))
    if a != b or a == "U":
        return
    original = int(a == "E")
    contradictory = False
    unobservable = False
    for record in (first, second):
        opposite, unknown = (record["deficits"][target][i] for i in (1-original, 2))
        contradictory |= opposite is not None and opposite <= 36
        unobservable |= unknown is not None and unknown <= 36
    cause = {(False, False): "retained", (True, False): "contradiction_only",
             (False, True): "unobservable_only", (True, True): "both"}[contradictory, unobservable]
    assert (cause == "retained") == all(classify(record, target, 36) == a for record in (first, second))
    return a, cause


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, required=True)
    p.add_argument("--universes", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    runs = json.loads((a.repository / "raw_replay/metadata/RUNS.json").read_bytes())["runs"]
    counts = []
    for run in runs:
        if run["condition"] not in ("PE treated", "mock"):
            continue
        accession = run["run_accession"]
        statuses = Counter()
        example = defaultdict(dict)
        with gzip.open(a.universes / f"{accession}_mutant.jsonl.gz", "rt", encoding="ascii") as x, gzip.open(
                a.universes / f"{accession}_programmed.jsonl.gz", "rt", encoding="ascii") as y:
            n = 0
            for i, j in zip(x, y, strict=True):
                r, s = json.loads(i), json.loads(j)
                assert r["id"] == s["id"]
                n += 1
                for idx, name in enumerate(("native", "exact")):
                    result = compare(r, s, idx)
                    if result is not None:
                        label, cause = result
                        statuses[name, label, cause] += 1
                        example[name, label].setdefault(cause, r["id"])
            assert n == run["raw_reads"]
        counts.extend({"accession": accession, "preparation": run["source_replica"],
                       "condition": run["condition"], "target": t, "label_at_27": label,
                       "outcome_at_36": cause, "n": count, "example": example[t, label][cause]}
                      for (t, label, cause), count in sorted(statuses.items()))
    a.output.write_text(json.dumps(counts, indent=2) + "\n", encoding="utf-8")
    for row in counts:
        if row["target"] == "native" and row["condition"] == "PE treated" and row["label_at_27"] == "E":
            print(row)


if __name__ == "__main__":
    main()
