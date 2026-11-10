"""Reconcile all raw identifiers and report nested, cost-stratified admissions."""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
result = json.loads((ROOT/"ANCHOR_RESULT.json").read_bytes())
rows_path = ROOT/"ANCHOR_RECORDS.jsonl.gz"
assert hashlib.sha256(rows_path.read_bytes()).hexdigest() == result["output_sha256"]
with gzip.open(rows_path, "rt", encoding="ascii") as stream:
    rows = list(map(json.loads, stream))
assert len(rows) == len({r["id"] for r in rows}) == 2095
deltas = (0, 9, 18, 27, 36, 45, 54)
counter = Counter()
for row in rows:
    observed = [c for c in row["calls"] if c != "U"]
    assert row["calls"] == observed + ["U"]*(len(deltas)-len(observed))
    assert len(set(observed)) <= 1
    status = row["admission"]
    group = "original" if status == "exact" else (
        "ordered_exact_pair" if status == "added" and row["cost"] == 0 else
        "substitution_rescue" if status == "added" else status
    )
    for delta, call in zip(deltas, row["calls"], strict=True):
        outcome = "U" if call == "U" else "concordant" if call == row["expected"] else "discordant"
        counter[row["clone"], group, delta, outcome] += 1
report = [
    {"clone": clone, "group": group, "deficit": delta, "outcome": outcome, "n": n}
    for (clone,group,delta,outcome),n in sorted(counter.items())
]
for clone, n in (("WT",504), ("HETS",797), ("MT",794)):
    for delta in deltas:
        assert sum(r["n"] for r in report if r["clone"] == clone and r["deficit"] == delta) == n
(ROOT/"STRATIFIED_RESULT.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
