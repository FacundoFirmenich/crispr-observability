"""Apply a single-substitution outer-anchor sensitivity test to the MP10 reads."""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
MP10 = PROJECT / "CRISPR_P1_MP10_CONTROL"
sys.path.insert(0, str(MP10))
from keap1_control import DELTAS, fastq, rc, segment
from recompute_panel import NEG, classify_intervals, encode, maximize


def neighborhood(anchor):
    variants = {anchor: 0}
    for i, original in enumerate(anchor):
        for base in "ACGT":
            if base != original:
                variants[anchor[:i] + base + anchor[i+1:]] = 1
    return variants


def hits(sequence, choices, width):
    return [(i, choices[sequence[i:i+width]]) for i in range(len(sequence)-width+1)
            if sequence[i:i+width] in choices]


def admission(sequence, left, right):
    orientations = (("+", sequence), ("-", rc(sequence)))
    exact = [(strand, query, span) for strand, query in orientations
             if (span := segment(query, left, right)) is not None]
    if len(exact) == 1:
        strand, query, (a, b) = exact[0]
        return "exact", strand, query[a:b], (a, b), 0
    lm, rm = neighborhood(left), neighborhood(right)
    candidates = [(x+y, strand, a, b+len(right), query)
                  for strand, query in orientations
                  for a, x in hits(query, lm, len(left))
                  for b, y in hits(query, rm, len(right)) if a+len(left) <= b]
    if not candidates:
        return "no_ordered_pair", None, None, None, None
    best = min(row[0] for row in candidates)
    best_pairs = [row for row in candidates if row[0] == best]
    if len(best_pairs) != 1:
        return "tied", None, None, None, best
    _, strand, a, b, query = best_pairs[0]
    return "added", strand, query[a:b], (a, b), best


def labels(query, references, lo, hi, target):
    if query is None:
        return ["U"] * len(DELTAS), []
    encoded = encode(query)
    classes = classify_intervals(query, {"programmed_local": target, "unedited_local": references[0][lo:hi]})
    deficits = []
    for ref, end in zip(references, (hi, hi-17), strict=True):
        optimum, maxima, _, _, suffix = maximize(encode(ref), encoded, lo, end, classes, max(DELTAS))
        assert optimum == suffix[0, 0, 0] and all(row.max() == optimum for row in maxima)
        deficits.append([None if s == NEG else int(optimum-s) for s in maxima[0]])
    calls = []
    for delta in DELTAS:
        admitted = [{k for k, score in zip(("NE","E","U"), scores, strict=True)
                     if score is not None and score <= delta} for scores in deficits]
        calls.append(next(iter(admitted[0])) if admitted[0] == admitted[1] and len(admitted[0]) == 1 else "U")
    return calls, deficits


def main():
    rule_bytes = (MP10 / "KEAP1_FREEZE.json").read_bytes()
    rule = json.loads(rule_bytes)
    base_path = MP10 / "KEAP1_ONT_READS.jsonl.gz"
    prior = json.loads((MP10 / "KEAP1_ONT_EVALUATION.json").read_bytes())
    assert hashlib.sha256(rule_bytes).hexdigest() == prior["freeze_sha256"]
    assert hashlib.sha256(base_path.read_bytes()).hexdigest() == prior["read_rows_sha256"]
    with gzip.open(base_path, "rt", encoding="ascii") as stream:
        original = {row["id"]: row for row in map(json.loads, stream)}
    assert len(original) == 2095
    stats, disposition = Counter(), Counter()
    seen = set()
    refs = (rule["reference_wt"], rule["reference_mt"])
    output = ROOT / "ANCHOR_RECORDS.jsonl.gz"
    with gzip.open(output, "wt", encoding="ascii") as stream:
        for name, spec in rule["clones"].items():
            rows, sha = fastq(spec["ont"])
            assert sha == prior["ont_sha256"][spec["ont"]]
            for identity, sequence, _ in rows:
                seen.add(identity)
                old = original[identity]
                assert old["clone"] == name and old["expected"] == spec["expected"]
                status, strand, query, span, cost = admission(sequence, rule["outer_left"], rule["outer_right"])
                calls, deficits = labels(query, refs, rule["local_lo"], rule["local_hi"], rule["target"])
                if status == "exact":
                    assert strand == old["strand"] and calls == old["calls"] and deficits == old["deficits"]
                else:
                    assert old["strand"] is None and set(old["calls"]) == {"U"}
                disposition[name, status] += 1
                for delta, call in zip(DELTAS, calls, strict=True):
                    result = "U" if call == "U" else "agree" if call == spec["expected"] else "disagree"
                    stats[name, status, delta, result] += 1
                stream.write(json.dumps({"id": identity, "clone": name, "expected": spec["expected"],
                                         "admission": status, "strand": strand, "span": span, "cost": cost,
                                         "calls": calls, "deficits": deficits}) + "\n")
    assert seen == set(original)
    result = {
        "protocol_sha256": hashlib.sha256((ROOT/"ANCHOR_PROTOCOL.md").read_bytes()).hexdigest(),
        "source_freeze_sha256": prior["freeze_sha256"],
        "source_ledger_sha256": prior["read_rows_sha256"],
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "admission": [{"clone": c, "status": a, "n": n} for (c,a),n in sorted(disposition.items())],
        "counts": [{"clone": c, "status": a, "deficit": d, "result": r, "n": n}
                   for (c,a,d,r),n in sorted(stats.items())],
    }
    (ROOT / "ANCHOR_RESULT.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(result["admission"])


if __name__ == "__main__":
    main()
