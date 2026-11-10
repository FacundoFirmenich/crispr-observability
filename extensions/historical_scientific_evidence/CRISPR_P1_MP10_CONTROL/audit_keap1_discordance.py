"""Recover and rescore all disagreements with the published clone classes."""
from collections import Counter
import gzip
import hashlib
import json
import sys

from keap1_control import ROOT, fastq, rc, segment, write_json
from recompute_panel import solve

sys.path.insert(0, str(ROOT.parent / "CRISPR_P1_MP07_RECUPERACION/i2/extensions/correspondence_frontier"))
import optimal_lattice as lattice

rule = json.loads((ROOT / "KEAP1_FREEZE.json").read_bytes())
evaluation = json.loads((ROOT / "KEAP1_ONT_EVALUATION.json").read_bytes())
ledger = ROOT / "KEAP1_ONT_READS.jsonl.gz"
assert hashlib.sha256(ledger.read_bytes()).hexdigest() == evaluation["read_rows_sha256"]
raw = {spec["ont"]: {identity: (seq, quality) for identity, seq, quality in fastq(spec["ont"])[0]}
       for spec in rule["clones"].values()}
wt, mt = rule["reference_wt"], rule["reference_mt"]
lo, hi = rule["local_lo"], rule["local_hi"]
witnesses = []
counts = Counter()
with gzip.open(ledger, "rt", encoding="ascii") as stream:
    for line in stream:
        row = json.loads(line)
        counts[row["clone"], row["expected"], row["baseline"], row["calls"][0]] += 1
        if row["calls"][0] in ("U", row["expected"]):
            continue
        seq, quality = raw[row["run"]][row["id"]]
        if row["strand"] == "-":
            seq, quality = rc(seq), quality[::-1]
        a, b = segment(seq, rule["outer_left"], rule["outer_right"])
        query = seq[a:b]
        reports = []
        for reference, end in ((wt, hi), (mt, hi-17)):
            case = {
                "reference_sequence": reference, "observed_sequence": query,
                "lo": lo, "hi": end, "target_kind": "HEK3_DNA",
                "programmed_local": rule["target"], "unedited_local": wt[lo:hi],
            }
            optimum, paths = solve(case, lattice, cutoff=-1)
            reports.append({"reference": "WT" if reference == wt else "MT",
                            "optimum": optimum, "local_paths": paths["programmed_local"]})
        row.update({"outer_observed": query, "outer_quality": quality[a:b],
                    "exact_full_c2": query == mt, "paths": reports})
        witnesses.append(row)
write_json("KEAP1_DISCORDANCE_WITNESSES.json", {
    "all_optimal_discordant_records": witnesses,
    "paired_baseline": [{"clone": c, "expected": t, "baseline": b, "optimum": o, "n": n}
                        for (c,t,b,o),n in sorted(counts.items())],
})
print("Discordances", Counter(row["clone"] for row in witnesses))
for row in witnesses:
    if row["clone"] == "HETS" or row["calls"][3] != "U":
        print(row["id"], row["clone"], "full exact C2", row["exact_full_c2"],
              "anchored length",len(row["outer_observed"]), "calls",row["calls"])
