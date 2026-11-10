"""Apply the frozen MP08 audit to the twelve publisher-declared ENA runs."""

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sys

from numba import njit

from freeze import ROOT, DELTAS, rc, sha

PREVIOUS = ROOT.parent / "CRISPR_P1_MP07_RECUPERACION"
sys.path.insert(0, str(PREVIOUS))
from recompute_panel import M, I, D, NEG, classify_intervals, encode, maximize, opening, substitution
from recompute_universes import label_at


def fastq(path):
    with gzip.open(path, "rt", encoding="ascii") as f:
        while header := f.readline():
            sequence, plus, quality = (f.readline().rstrip("\r\n") for _ in range(3))
            assert header.startswith("@") and plus.startswith("+") and sequence and len(sequence) == len(quality)
            name = header.split()[0][1:]
            if name.endswith(("/1", "/2")):
                name = name[:-2]
            assert set(sequence) <= set("ACGTN")
            yield name, sequence, quality


def mean_q30(quality):
    return sum(map(ord, quality)) >= 63 * len(quality)


def primer_clip(sequence, primer, opposing):
    start = sequence.find(primer, 0, len(primer) + 8)
    if start < 0:
        return None
    oriented = sequence[start:]
    end = oriented.find(opposing, 24)
    return oriented[:end+len(opposing)] if end >= 0 else oriented


def merge(seq_a, seq_b, head, tail):
    pairs = []
    for a, b in ((seq_a, seq_b), (seq_b, seq_a)):
        left, reverse = primer_clip(a, head, rc(tail)), primer_clip(b, tail, rc(head))
        if left is not None and reverse is not None:
            pairs.append((left, rc(reverse)))
    if len(pairs) != 1:
        return "primer_orientation", None
    a, b = pairs[0]
    starts = []
    at = a.find(b[:24])
    while at >= 0:
        overlap = len(a) - at
        if 24 <= overlap <= len(b) and a[at:] == b[:overlap]:
            starts.append(at)
        at = a.find(b[:24], at+1)
    if len(starts) != 1:
        return "overlap_absent" if not starts else "overlap_ambiguous", None
    sequence = a[:starts[0]] + b
    if "N" in sequence:
        return "non_acgt", None
    return "evaluated", sequence


@njit(nogil=True, cache=True)
def deterministic_projection(reference, read, lo, hi, prefix, optimum, labels):
    n, m = len(reference), len(read)
    i, j = n, m
    state = M if prefix[M, i, j] == optimum else I if prefix[I, i, j] == optimum else D
    start = end = -1
    left_flank = right_flank = False
    while i or j:
        value = prefix[state, i, j]
        if state != I:
            if i == lo:
                start = j
                left_flank = state == M
            elif i == hi:
                end = j
            elif i == hi+1:
                right_flank = state == M
        if state == M:
            previous = value - substitution(reference[i-1], read[j-1])
            i -= 1
            j -= 1
            state = M if prefix[M, i, j] == previous else I if prefix[I, i, j] == previous else D
        elif state == I:
            cost = opening(i, j, n, m)
            j -= 1
            state = M if prefix[M, i, j] + cost == value else I
        else:
            cost = opening(i, j, n, m)
            i -= 1
            state = M if prefix[M, i, j] + cost == value else D
    if not left_flank or not right_flank:
        return 2, 2
    assert 0 <= start <= end
    return labels[0, start, end], labels[1, start, end]


def analyze(sequence, references, lo, hi_values, targets):
    read = encode(sequence)
    labels = classify_intervals(sequence, targets)
    values = []
    for ref, hi in zip(references, hi_values, strict=True):
        optimum, maxima, endpoints, prefix, suffix = maximize(ref, read, lo, hi, labels, DELTAS[-1])
        assert optimum == suffix[M, 0, 0]
        deficits = [[None if x == NEG else int(optimum-x) for x in row] for row in maxima]
        point = tuple(map(int, deterministic_projection(ref, read, lo, hi, prefix, optimum, labels)))
        for target, kind in enumerate(point):
            assert deficits[target][kind] == 0
            assert all(x is None or x in DELTAS for x in deficits[target])
        values.append({"score": int(optimum), "deficits": deficits, "single_trace": point,
                       "endpoints": endpoints.tolist()})
    return values


def agree(a, b):
    return a if a == b else "U"


def metrics(results, weight, names, counts, paired):
    for index, name in enumerate(names):
        point = agree(*[("NE", "E", "U")[r["single_trace"][index]] for r in results])
        for delta in DELTAS:
            certified = agree(*[label_at(r["deficits"][index], delta) for r in results])
            counts[name, delta, certified] += weight
            paired[name, delta, point, certified] += weight


def evaluate(sample, spec, directory, pool):
    run = sample["run"]
    path = directory / f"{run}.json"
    if path.exists():
        value = json.loads(path.read_bytes())
        assert all(sha(directory / n) == h for n, h in value["output_sha256"].items())
        return value
    paths = [ROOT / "fastq" / x["url"].split("/")[-1] for x in sample["files"]]
    groups = defaultdict(list)
    reasons = Counter()
    seen = set()
    out_reads = directory / f"{run}_records.tsv.gz"
    out_families = directory / f"{run}_families.jsonl.gz"
    with gzip.open(out_reads.with_suffix(".tmp"), "wt", encoding="ascii", newline="\n") as f:
        f.write("read_id\tstatus\toriented_sequence_sha256\n")
        for a, b in zip(fastq(paths[0]), fastq(paths[1]), strict=True):
            name, sa, qa = a
            assert name == b[0] and name not in seen
            seen.add(name)
            status, sequence = ("quality_fail", None)
            if mean_q30(qa) and mean_q30(b[2]):
                status, sequence = merge(sa, b[1], spec["reference_head_primer"], spec["reference_tail_primer"])
            reasons[status] += 1
            key = ""
            if sequence is not None:
                if sequence not in groups:
                    groups[sequence] = [hashlib.sha256(sequence.encode()).hexdigest(), 0]
                key = groups[sequence][0]
                groups[sequence][1] += 1
            f.write(f"{name}\t{status}\t{key}\n")
    n = sum(reasons.values())
    assert n == sample["ENA_read_count"], (run, n, sample["ENA_read_count"])
    print("START", run, sample["locus"], n, dict(reasons), "families", len(groups), flush=True)
    refs = tuple(encode(spec[k]) for k in ("wt", "edited"))
    his = spec["hi_wt"], spec["hi_edited"]
    targets = {k: spec[k] for k in ("programmed_local", "unedited_local")}

    def score(seq):
        return analyze(seq, refs, spec["lo"], his, targets)

    counts, paired = Counter(), Counter()
    seqs = list(groups)
    with gzip.open(out_families.with_suffix(".tmp"), "wt", encoding="ascii", newline="\n") as f:
        for start in range(0, len(seqs), 256):
            batch = seqs[start:start+256]
            for seq, values in zip(batch, pool.map(score, batch), strict=True):
                key, weight = groups[seq]
                metrics(values, weight, targets, counts, paired)
                f.write(json.dumps({"key": key, "sequence": seq, "N": weight, "references": values},
                                   separators=(",", ":")) + "\n")
    recount, repaired = Counter(), Counter()
    with gzip.open(out_families.with_suffix(".tmp"), "rt", encoding="ascii") as f:
        for line in f:
            row = json.loads(line)
            assert [row["key"],row["N"]] == groups[row["sequence"]]
            metrics(row["references"], row["N"], targets, recount, repaired)
    assert recount == counts and repaired == paired
    with gzip.open(out_reads.with_suffix(".tmp"), "rt", encoding="ascii") as f:
        next(f)
        saved_reasons, saved_families = Counter(), Counter()
        for line in f:
            _, status, key = line.rstrip("\n").split("\t")
            saved_reasons[status] += 1
            if key:
                saved_families[key] += 1
        assert saved_reasons == reasons and saved_families == Counter({k:w for k,w in groups.values()})
    for name in targets:
        for delta in DELTAS:
            omitted = n - reasons["evaluated"]
            counts[name, delta, "U"] += omitted
            paired[name, delta, "U", "U"] += omitted
            assert sum(counts[name, delta, k] for k in ("E","NE","U")) == n
    out_reads.with_suffix(".tmp").replace(out_reads)
    out_families.with_suffix(".tmp").replace(out_families)
    result = {
        **{k: sample[k] for k in ("locus", "replicate", "arm", "run")}, "N": n, "exclusions": dict(reasons),
        "families": len(groups), "output_sha256": {p.name: sha(p) for p in (out_reads, out_families)},
        "counts": [{"target": t, "delta": d, **{k: counts[t,d,k] for k in ("E","NE","U")},
                    "lower": str(Fraction(counts[t,d,"E"], n)), "upper": str(Fraction(n-counts[t,d,"NE"], n))}
                   for t in targets for d in DELTAS],
        "same_score_single_trace_vs_certificate": [
            {"target": t, "delta": d, "single_trace": a, "certificate": b, "N": c}
            for (t,d,a,b),c in sorted(paired.items())],
    }
    path.with_suffix(".part").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    path.with_suffix(".part").replace(path)
    print("VERIFIED", run, flush=True)
    return result


def summarize(results):
    pcsk9 = {(r["replicate"],r["arm"]):r for r in results if r["locus"] == "PCSK9"}
    contrasts = []
    for replicate in (1,2,3):
        t = pcsk9[replicate,"PE_unmodified"]
        c = pcsk9[replicate,"untreated"]
        for a,b in zip(t["counts"],c["counts"],strict=True):
            assert (a["target"],a["delta"]) == (b["target"],b["delta"])
            if a["target"] != "programmed_local":
                continue
            low = Fraction(a["lower"]) - Fraction(b["upper"])
            high = Fraction(a["upper"]) - Fraction(b["lower"])
            contrasts.append({"replicate":replicate,"delta":a["delta"],"lower":str(low),"upper":str(high),
                              "lower_pp":float(100*low),"upper_pp":float(100*high)})
    return contrasts


if __name__ == "__main__":
    manifest = ROOT / "EXTERNAL_FREEZE.json"
    freeze = json.loads(manifest.read_bytes())
    acquisition = json.loads((ROOT / "ACQUISITION.json").read_bytes())
    assert acquisition["freeze_sha256"] == sha(manifest)
    for entry in acquisition["inputs"]:
        assert sha(ROOT / "fastq" / entry["file"]) == entry["sha256"]
    directory = ROOT / "external_results"
    directory.mkdir(exist_ok=True)
    with ThreadPoolExecutor(4) as pool:
        results = [evaluate(s, freeze["targets"][s["locus"]], directory, pool) for s in freeze["samples"]]
    assert len(results) == 12
    (directory / "EXTERNAL_RESULTS.json").write_text(json.dumps({
        "freeze_sha256": sha(manifest), "input_receipt_sha256": sha(ROOT / "ACQUISITION.json"),
        "N": sum(r["N"] for r in results), "samples": results,
        "PCSK9_PE_minus_untreated": summarize(results)}, indent=2) + "\n", encoding="utf-8")
