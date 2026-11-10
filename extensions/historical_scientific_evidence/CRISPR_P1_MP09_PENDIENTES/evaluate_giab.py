"""Evaluate frozen MP07 exact-tract certificates against external phased markers."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "CRISPR_P1_MP07_RECUPERACION"))
from recompute_panel import NEG, classify_intervals, encode, maximize
from reconcile_crispra import observed_projection

DELTAS = (0, 9, 18, 27, 36, 45, 54)


def label(deficits, delta):
    admitted = [i for i, x in enumerate(deficits) if x is not None and x <= delta]
    return ("NE", "E", "U")[admitted[0]] if len(admitted) == 1 else "U"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    evidence = args.evidence
    checks = json.loads((evidence / "SHA256.json").read_bytes())
    for name, expected in checks.items():
        path = evidence / name
        if name in ("benchmark.vcf.gz", "benchmark.bed"):
            continue
        with path.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == expected
    cases = {p["indel"]["pos"] + 1: p for p in json.loads((evidence / "CASES.json").read_bytes())}
    counts = Counter()
    sites = Counter()
    identities = set()
    outpath = ROOT / "GIAB_CERTIFICATES.jsonl.gz"
    with gzip.open(evidence / "reads.jsonl.gz", "rt") as source, gzip.open(outpath, "wt") as output:
        for line in source:
            row = json.loads(line)
            case = cases[row["site"]]
            key = (row["site"], row["read_id"], row["flag"] & 192)
            assert key not in identities
            identities.add(key)
            seq = row["sequence"]
            assert set(seq) <= set("ACGTN")
            p = case["indel"]["pos"] - case["start"]
            nref = len(case["indel"]["ref"])
            alt = case["indel"]["alts"][0]
            reference = case["reference"]
            alternate = reference[:p] + alt + reference[p+nref:]
            lo, hi = case["local_lo"], case["local_hi"]
            targets = {"programmed_local": case["target"], "unedited_local": reference[lo:hi]}
            classifications = classify_intervals(seq, targets)
            deficits = []
            for ref, end in ((reference, hi), (alternate, hi + len(alt) - nref)):
                optimum, maxima, endpoints, prefix, suffix = maximize(
                    encode(ref), encode(seq), lo, end, classifications, 54)
                assert optimum == suffix[0, 0, 0] and all(maximum.max() == optimum for maximum in maxima)
                deficits.append([None if x == NEG else int(optimum-x) for x in maxima[0]])
            projected = observed_projection(row["position"] - case["start"] + 1,
                                            row["cigar"], seq, lo, hi)
            native = "U" if projected is None else "E" if projected == case["target"] else "NE"
            decisions = []
            for delta in DELTAS:
                a, b = (label(d, delta) for d in deficits)
                decision = a if a == b else "U"
                decisions.append(decision)
                counts[delta, row["truth"], native, decision] += 1
                sites[row["site"], delta, row["truth"], decision] += 1
            row.update({"deficits": deficits, "native_bam_exact_tract": native,
                        "deltas": DELTAS, "certificate": decisions})
            output.write(json.dumps(row) + "\n")
    result = {
        "N_reads": len(identities), "N_sites": len(cases),
        "scope": "External HG002 read-level phase-marker concordance, conditional on benchmark phase and correct marker identification. BAM baseline is its published alignment. CRISPR-A was not run on this control. Counts are paired reads and are not independent molecule trials.",
        "counts": [{"delta": d, "truth": t, "bam": b, "certificate": c, "n": n}
                   for (d, t, b, c), n in sorted(counts.items())],
        "sites": [{"site": s, "delta": d, "truth": t, "certificate": c, "n": n}
                  for (s, d, t, c), n in sorted(sites.items())],
    }
    with outpath.open("rb") as stream:
        result["records_sha256"] = hashlib.file_digest(stream, "sha256").hexdigest()
    (ROOT / "GIAB_EVALUATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"N": len(identities), "sites": len(cases)}), flush=True)


if __name__ == "__main__":
    main()
