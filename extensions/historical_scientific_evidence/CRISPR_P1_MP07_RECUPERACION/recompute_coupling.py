"""Rebuild the CFTR operator-coupling comparison from the complete I1 ledger."""

import argparse
from collections import Counter
import csv
from fractions import Fraction as F
import gzip
import json
from pathlib import Path

from recompute_panel import file_hash, is_event
from recompute_universes import LEDGER_HASHES


OPERATORS = ("T", "M_fw", "P_fw", "M_rc", "P_rc", "MP_fw", "MP_rc", "PE_default", "PE_zero")
EVENTS = ("native", "exact_programmed", "native_nonexact")


def events(coding, native, programmed):
    if not coding:
        return ("U",) * len(EVENTS)
    n = is_event("native", native, coding)
    e = coding == programmed
    assert not e or n
    return tuple("E" if flag else "NE" for flag in (n, e, n and not e))


def build(ledger, native, programmed):
    assert file_hash(ledger) == LEDGER_HASHES["cftr"]
    counts = Counter()
    totals = Counter()
    with gzip.open(ledger, "rt", encoding="ascii", newline="") as f:
        for row in csv.DictReader(f):
            run = row["run"]
            totals[run] += 1
            labels = [events(row[operator], native, programmed) for operator in OPERATORS]
            for target, column in zip(EVENTS, zip(*labels), strict=True):
                for operator, label in zip(OPERATORS, column, strict=True):
                    counts[run, operator, target, label] += 1
                counts[run, "F9", target, column[0] if len(set(column)) == 1 else "U"] += 1
    assert sum(totals.values()) == 177006
    return counts, totals


def validate(path, counts, totals, is_family):
    checked = 0
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if is_family and (row["family"], row["certificate"]) != ("F9", "target_relative"):
                continue
            run, event = row["run"], row["event"]
            op = "F9" if is_family else row["mode"]
            assert totals[run] == int(row["N"])
            assert all(counts[run, op, event, k] == int(row[k]) for k in ("E", "NE", "U"))
            checked += 1
    return checked


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, required=True)
    p.add_argument("--ledger", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    raw = a.repository / "raw_replay/metadata"
    spec = json.loads((raw / "REFERENCE.json").read_bytes())
    counts, totals = build(a.ledger, spec["native_peptide"], spec["exact_programmed_core"] + "AT")
    frozen = a.repository / "extensions/correspondence_frontier/expected/mechanism"
    assert validate(frozen / "PROFILE_EVENT_COUNTS.csv", counts, totals, False) == 243
    assert validate(frozen / "FAMILY_CERTIFICATE_COUNTS.csv", counts, totals, True) == 27
    runs = json.loads((raw / "RUNS.json").read_bytes())["runs"]

    def bound(run, op, event):
        e, u = (counts[run, op, event, k] for k in ("E", "U"))
        return F(e, totals[run]), F(e+u, totals[run])

    rows = []
    for replica in (1, 2, 3):
        ids = {r["condition"]: r["run_accession"] for r in runs if r["source_replica"] == replica}
        treatment = ids["PE treated"]
        for control in ("mock", "untreated"):
            for event in EVENTS:
                t = {op: bound(treatment, op, event) for op in (*OPERATORS, "F9")}
                c = {op: bound(ids[control], op, event) for op in (*OPERATORS, "F9")}
                shared = (min(t[op][0]-c[op][1] for op in OPERATORS),
                          max(t[op][1]-c[op][0] for op in OPERATORS))
                separate = (min(t[op][0] for op in OPERATORS)-max(c[op][1] for op in OPERATORS),
                            max(t[op][1] for op in OPERATORS)-min(c[op][0] for op in OPERATORS))
                recordwise = t["F9"][0]-c["F9"][1], t["F9"][1]-c["F9"][0]
                assert shared == separate == (t["T"][0]-c["T"][1], t["T"][1]-c["T"][0])
                assert recordwise[0] <= shared[0] <= shared[1] <= recordwise[1]
                rows.append({"preparation": replica, "control": control, "event": event,
                             "shared": list(map(str, shared)), "separate": list(map(str, separate)),
                             "recordwise": list(map(str, recordwise)),
                             "recordwise_width_increment_pp": float(100 * (
                                 recordwise[1]-recordwise[0]-shared[1]+shared[0])),
                             "recordwise_positive": recordwise[0] > 0})
    receipt = {"records": sum(totals.values()), "verified_operator_rows": 243,
               "verified_F9_rows": 27, "contrasts": rows}
    a.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    for row in rows:
        if row["event"] == "native" and row["control"] == "mock":
            print(row)


if __name__ == "__main__":
    main()
