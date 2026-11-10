"""Localize external certificate withdrawals and audit real attaining paths."""

from collections import Counter
from fractions import Fraction
import gzip
import json
import sys

from freeze import ROOT, DELTAS, sha
from run_external import PREVIOUS, agree, label_at
from recompute_panel import solve

LATTICE = PREVIOUS / "i2/extensions/correspondence_frontier"
sys.path.insert(0, str(LATTICE))
import optimal_lattice


def withdrawal(refs):
    label = agree(*[label_at(r["deficits"][0], 0) for r in refs])
    if label != "E":
        return None
    entries = [(r["deficits"][0][k], k) for r in refs for k in (0,2)
               if r["deficits"][0][k] is not None]
    if not entries:
        return None
    d = min(x[0] for x in entries)
    kinds = {k for x,k in entries if x == d}
    return d, "both" if len(kinds) == 2 else "unobservable" if 2 in kinds else "contradiction"


def flank_bases(ref, query, lo, hi):
    columns = [i for i,b in enumerate(ref) if b != "-"]
    return [query[columns[k]] for k in (lo-1,hi)]


def materialize(family, spec):
    value = []
    for k, hi_key, reference in (("wt","hi_wt","wt"),("edited","hi_edited","edited")):
        case = {**spec, "target_kind":"HEK3_DNA", "reference_sequence":spec[reference],
                "observed_sequence":family["sequence"],"hi":spec[hi_key]}
        optimum, targets = solve(case, optimal_lattice, cutoff=-1)
        for labels in targets.values():
            for w in labels.values():
                w["deficit"] = optimum-w["score"]
                w["flank_bases"] = flank_bases(w["aligned_reference"],w["aligned_read"],spec["lo"],spec[hi_key])
        value.append({"reference":k,"score":optimum,"targets":targets})
    return value


def main():
    freeze = json.loads((ROOT/"EXTERNAL_FREEZE.json").read_bytes())
    folder = ROOT/"external_results"
    data = json.loads((folder/"EXTERNAL_RESULTS.json").read_bytes())
    counts, examples, totals, sign = [], {}, [], []
    for r in data["samples"]:
        path = folder/f'{r["run"]}_families.jsonl.gz'
        assert sha(path) == r["output_sha256"][path.name]
        n, result = 0, Counter()
        with gzip.open(path,"rt",encoding="ascii") as f:
            for row in map(json.loads,f):
                n += row["N"]
                cause = withdrawal(row["references"])
                if cause:
                    result[cause] += row["N"]
                    if r["arm"] == "PE_unmodified" and cause[0] == 45:
                        key = r["locus"]
                        rank = (row["N"],row["key"])
                        if key not in examples or rank > examples[key][0]:
                            examples[key] = (rank,r["run"],row)
        assert n == r["exclusions"]["evaluated"]
        base = next(x for x in r["counts"] if (x["target"],x["delta"]) == ("programmed_local",0))
        last = next(x for x in r["counts"] if (x["target"],x["delta"]) == ("programmed_local",54))
        assert sum(result.values()) == base["E"]-last["E"]
        counts.extend({"run":r["run"],"locus":r["locus"],"arm":r["arm"],"delta":d,
                       "cause":cause,"N":n} for (d,cause),n in sorted(result.items()))
        upstream = r["N"]-n
        totals.append({"run":r["run"],"locus":r["locus"],"arm":r["arm"],
                       "raw_N":r["N"],"upstream_U":upstream,"alignment_U_at_optimum":base["U"]-upstream})
    by_rep = {(r["replicate"],r["arm"]):r for r in data["samples"] if r["locus"] == "PCSK9"}
    for i in (1,2,3):
        t,c = (by_rep[i,a] for a in ("PE_unmodified","untreated"))
        rows = []
        for sample in (t,c):
            cnt = Counter()
            for x in sample["same_score_single_trace_vs_certificate"]:
                if x["target"]=="programmed_local" and x["delta"]==0:
                    cnt[x["single_trace"]] += x["N"]
            assert sum(cnt.values()) == sample["N"]
            rows.append(cnt)
        mt,mc = rows
        point_complete_case = Fraction(mt["E"],mt["E"]+mt["NE"])-Fraction(mc["E"],mc["E"]+mc["NE"])
        point_raw_lower = Fraction(mt["E"],t["N"])-Fraction(mc["E"]+mc["U"],c["N"])
        sign.append({"replicate":i,"single_trace_complete_case_pp":float(100*point_complete_case),
                     "single_trace_raw_lower_pp":float(100*point_raw_lower)})
    paths = []
    for locus,(_,run,family) in sorted(examples.items()):
        maximum = materialize(family,freeze["targets"][locus])
        stored = family["references"]
        for true, limited in zip(maximum,stored,strict=True):
            assert true["score"] == limited["score"]
            for j,target in enumerate(true["targets"].values()):
                for k in (0,1,2):
                    expected = target.get(k,{}).get("deficit")
                    assert limited["deficits"][j][k] == (expected if expected in DELTAS else None)
        paths.append({"locus":locus,"run":run,"family_sha256":family["key"],"N":family["N"],
                      "sequence":family["sequence"],"paths":maximum})
    assert paths
    (ROOT/"EXTERNAL_DIAGNOSIS.json").write_text(json.dumps(
        {"upstream_vs_alignment":totals,"first_positive_withdrawal":counts,
         "pcsk9_single_trace":sign,"witnesses":paths},indent=2)+"\n",encoding="utf-8")
    print(sign)
    print([(p["locus"],p["run"],p["N"]) for p in paths])


if __name__ == "__main__":
    main()
