"""Measure KEAP1 Nanopore concordance with separately sequenced Illumina clones."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "CRISPR_P1_MP07_RECUPERACION"))
from recompute_panel import NEG, classify_intervals, encode, maximize

RUNS = {
    "WT": ("SRR18884574", "SRR18840484", "F6", "NE"),
    "HETS": ("SRR18884573", "SRR18840476", "B2", "NE"),
    "MT": ("SRR18884572", "SRR18840477", "C2", "E"),
}
DELTAS = (0, 9, 18, 27, 36, 45, 54)
GUIDE = "CCGTGGGAGCCGCCGACGG"
TRANS = str.maketrans("ACGTN", "TGCAN")


def rc(sequence):
    return sequence.translate(TRANS)[::-1]


def runs_metadata():
    with (ROOT / "PRJNA828591.tsv").open(encoding="utf-8") as file:
        return {x["run_accession"]: x for x in csv.DictReader(file, delimiter="\t")}


def fastq(run):
    path = ROOT / (run + ".fastq.gz")
    metadata = runs_metadata()[run]
    data = path.read_bytes()
    assert len(data) == int(metadata["fastq_bytes"])
    assert hashlib.md5(data).hexdigest() == metadata["fastq_md5"]
    rows = []
    with gzip.open(path, "rt", encoding="ascii") as file:
        for head, seq, separator, quality in zip(file, file, file, file, strict=True):
            head, seq, separator, quality = (x.rstrip("\r\n") for x in (head, seq, separator, quality))
            assert head.startswith("@") and separator.startswith("+")
            assert len(seq) == len(quality) and set(seq) <= set("ACGTN")
            rows.append((head[1:].split()[0], seq, quality))
    assert len(rows) == int(metadata["read_count"])
    assert len({row[0] for row in rows}) == len(rows)
    return rows, hashlib.sha256(data).hexdigest()


def segment(sequence, left, right):
    if sequence.count(left) == 1 and sequence.count(right) == 1:
        a, b = sequence.index(left), sequence.index(right)
        if a + len(left) <= b:
            return a, b + len(right)
    return None


def write_json(name, value):
    (ROOT / name).write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def freeze():
    wt_rows, _ = fastq(RUNS["WT"][0])
    seed, _ = Counter(seq for _, seq, _ in wt_rows if rc(GUIDE) in seq).most_common(1)[0]
    start = seed.index(rc(GUIDE))
    left, right = seed[start-40:start-25], seed[start+45:start+60]
    truth = {}
    for name, (illumina, ont, clone, expected) in RUNS.items():
        rows, digest = fastq(illumina)
        counts, q30 = Counter(), Counter()
        for _, seq, qual in rows:
            span = segment(seq, left, right)
            if span is not None:
                a, b = span
                counts[seq[a:b]] += 1
                if min(map(ord, qual[a:b])) >= 63:
                    q30[seq[a:b]] += 1
        leading = counts.most_common(2 if name == "HETS" else 1)
        assert all(q30[s] > 0 for s, _ in leading)
        truth[name] = {
            "illumina": illumina, "ont": ont, "clone": clone, "expected": expected,
            "illumina_sha256": digest, "raw": len(rows), "anchored": sum(counts.values()),
            "q30_all_bases": sum(q30.values()),
            "alleles": [{"sequence": s, "all": n, "q30_all_bases": q30[s]} for s, n in leading],
            "all_observed": dict(counts), "all_q30": dict(q30),
        }
        for run in (illumina, ont):
            alias = ET.parse(ROOT / (run + "_runbundle.xml")).getroot().find("RUN").get("alias")
            assert alias.lower().endswith(name.lower() + ".fastq")
            truth[name][("illumina" if run == illumina else "ont") + "_alias"] = alias
    wt = truth["WT"]["alleles"][0]["sequence"]
    mt = truth["MT"]["alleles"][0]["sequence"]
    assert len(wt) - len(mt) == 17
    assert {len(wt)-len(x["sequence"]) for x in truth["HETS"]["alleles"]} == {1, 10}
    prefix = next(i for i, (a, b) in enumerate(zip(wt, mt)) if a != b)
    assert wt[:prefix] + wt[prefix+17:] == mt
    lo, hi = prefix - 5, prefix + 17 + 5
    assert len(left) <= lo-7 and hi+7 <= len(wt)-len(right)
    contract = {
        "rule": "Published clones F6, B2, C2 and aliases, corroborated by distinct Illumina FASTQs. Exact C2 tract is E. F6 and both B2 modes are NE. Tests are clone-conditional, not same-molecule orthogonal measurements.",
        "admission": "Unique ordered exact 15-base outer anchors in exactly one orientation. Retain every raw ONT record, with U for no qualifying span. No ONT quality filter. Only Illumina opened for allele definition.",
        "truth_source": "Nguyen et al. 10.1093/nar/gkac440 Figure 2D and Supplementary information, PRJNA828591 run bundles",
        "outer_left": left, "outer_right": right, "reference_wt": wt, "reference_mt": mt,
        "local_lo": lo, "local_hi": hi, "target": mt[lo:hi-17],
        "deltas": DELTAS, "clones": truth,
    }
    write_json("KEAP1_FREEZE.json", contract)
    (ROOT / "KEAP1_CLONE_REFERENCES.fasta").write_text(
        "".join(f">{name}_{data['clone']}_Illumina_rank{rank}\n{allele['sequence']}\n"
                for name, data in truth.items() for rank, allele in enumerate(data["alleles"], 1)),
        encoding="ascii")
    print("Frozen", {name: (d["raw"], d["anchored"], [a["all"] for a in d["alleles"]])
                    for name, d in truth.items()})


def evaluate():
    frozen = (ROOT / "KEAP1_FREEZE.json").read_bytes()
    rule = json.loads(frozen)
    wt, mt = rule["reference_wt"], rule["reference_mt"]
    lo, hi, target = rule["local_lo"], rule["local_hi"], rule["target"]
    refs = ((encode(wt), hi), (encode(mt), hi-17))
    counts, admission, hashes = Counter(), Counter(), {}
    out = ROOT / "KEAP1_ONT_READS.jsonl.gz"
    with gzip.open(out, "wt", encoding="ascii") as output:
        for name, spec in rule["clones"].items():
            rows, hashes[spec["ont"]] = fastq(spec["ont"])
            for identity, seq, _ in rows:
                options = [(s, strand, segment(s, rule["outer_left"], rule["outer_right"]))
                           for s, strand in ((seq, "+"), (rc(seq), "-"))]
                hits = [(s, strand, span) for s, strand, span in options if span is not None]
                calls = ["U"] * len(DELTAS)
                baseline, strand, deficits = "U", None, []
                if len(hits) == 1:
                    s, strand, (a, b) = hits[0]
                    query = s[a:b]
                    encoded = encode(query)
                    classes = classify_intervals(query, {"programmed_local": target, "unedited_local": wt[lo:hi]})
                    for ref, end in refs:
                        optimum, maxima, _, _, suffix = maximize(ref, encoded, lo, end, classes, max(DELTAS))
                        assert optimum == suffix[0, 0, 0] and all(row.max() == optimum for row in maxima)
                        deficits.append([None if value == NEG else int(optimum-value) for value in maxima[0]])
                    for j, delta in enumerate(DELTAS):
                        labels = [{k for k, v in zip(("NE", "E", "U"), ds, strict=True) if v is not None and v <= delta}
                                  for ds in deficits]
                        if labels[0] == labels[1] and len(labels[0]) == 1:
                            calls[j] = next(iter(labels[0]))
                    inner_left, inner_right = wt[lo-7:lo], wt[hi:hi+7]
                    inner = segment(query, inner_left, inner_right)
                    if inner is not None:
                        observed = query[inner[0]+len(inner_left):inner[1]-len(inner_right)]
                        if observed and set(observed) <= set("ACGT"):
                            baseline = "E" if observed == target else "NE"
                admission[name, "admitted" if len(hits)==1 else "unresolved"] += 1
                counts[name, "exact_inner_anchors", baseline] += 1
                for delta, call in zip(DELTAS, calls, strict=True):
                    counts[name, f"delta_{delta}", call] += 1
                output.write(json.dumps({"run": spec["ont"], "id": identity, "clone": name,
                                         "expected": spec["expected"], "strand": strand,
                                         "baseline": baseline, "deficits": deficits, "calls": calls})+"\n")
    result = {"freeze_sha256": hashlib.sha256(frozen).hexdigest(), "ont_sha256": hashes,
              "read_rows_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
              "admission": [{"clone": c, "status": s, "n": n} for (c,s),n in sorted(admission.items())],
              "counts": [{"clone": c, "profile": p, "label": l, "n": n} for (c,p,l),n in sorted(counts.items())]}
    write_json("KEAP1_ONT_EVALUATION.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "evaluate"))
    args = parser.parse_args()
    (freeze if args.stage == "freeze" else evaluate)()
