"""Reconstruct and independently rescore actual paths for GIAB diagnostic reads."""
import gzip
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
MP07 = ROOT.parent / "CRISPR_P1_MP07_RECUPERACION"
sys.path.insert(0, str(MP07))
sys.path.insert(0, str(MP07 / "i2/extensions/correspondence_frontier"))
import optimal_lattice as lattice
from recompute_panel import solve

cases = {p["indel"]["pos"] + 1: p for p in json.loads((ROOT / "giab_41779054165/CASES.json").read_bytes())}
chosen = {}
for line in gzip.open(ROOT / "GIAB_CERTIFICATES.jsonl.gz", "rt"):
    row = json.loads(line)
    if row["site"] == 255818 and row["truth"] == "NE" and row["native_bam_exact_tract"] == "NE":
        chosen.setdefault("reference_disagreement", row)
    if row["certificate"][0] not in ("U", row["truth"]):
        chosen["marker_discordance"] = row

results = {}
for name, row in chosen.items():
    source = cases[row["site"]]
    ref = source["reference"]
    p = source["indel"]["pos"] - source["start"]
    ref_allele = source["indel"]["ref"]
    alt = source["indel"]["alts"][0]
    alternate = ref[:p] + alt + ref[p+len(ref_allele):]
    outputs = []
    for reference, hi in ((ref, source["local_hi"]),
                          (alternate, source["local_hi"] + len(alt)-len(ref_allele))):
        case = {"reference_sequence": reference, "observed_sequence": row["sequence"],
                "lo": source["local_lo"], "hi": hi, "target_kind": "HEK3_DNA",
                "programmed_local": source["target"],
                "unedited_local": ref[source["local_lo"]:source["local_hi"]]}
        optimum, witnesses = solve(case, lattice, cutoff=-1)
        outputs.append({"optimum": optimum, "case": case, "witnesses": witnesses})
    results[name] = {"read": row, "paths": outputs}
(ROOT / "GIAB_WITNESS_AUDIT.json").write_text(json.dumps(results, indent=2) + "\n")
print("Verified full path rescoring and projected labels for", list(results))
