"""Recompute MP07 target maxima from the authentic I2 read/reference panel."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

from Bio.Seq import Seq
from numba import njit
import numpy as np


NEG = -1_000_000
M, I, D = 0, 1, 2
STATES = (M, I, D)
ALPHABET = "ACGTN"
TRANSLATION = {
    a + b + c: str(Seq(a + b + c).translate())
    for a in ALPHABET[:4] for b in ALPHABET[:4] for c in ALPHABET[:4]
}


@njit
def substitution(a, b):
    if a == 4 or b == 4:
        return -1 if a == b else -2
    return 5 if a == b else -4


@njit
def opening(i, j, n, m):
    return -2 if i == 0 or j == 0 or i == n or j == m else -20


@njit
def forward(reference, read, origin_i, origin_j, end_i, end_j):
    n, m = len(reference), len(read)
    scores = np.full((3, end_i - origin_i + 1, end_j - origin_j + 1), NEG, np.int32)
    scores[M, 0, 0] = 0
    for a in range(scores.shape[1]):
        i = origin_i + a
        for b in range(scores.shape[2]):
            j = origin_j + b
            cost = opening(i, j, n, m)
            if a and b:
                scores[M, a, b] = max(scores[M, a-1, b-1], scores[I, a-1, b-1],
                                       scores[D, a-1, b-1]) + substitution(reference[i-1], read[j-1])
            if b:
                scores[I, a, b] = max(scores[M, a, b-1] + cost, scores[I, a, b-1] - 2)
            if a:
                scores[D, a, b] = max(scores[M, a-1, b] + cost, scores[D, a-1, b] - 2)
    return scores


@njit
def backward(reference, read):
    n, m = len(reference), len(read)
    scores = np.full((3, n+1, m+1), NEG, np.int32)
    scores[:, n, m] = 0
    for i in range(n, -1, -1):
        for j in range(m, -1, -1):
            if i == n and j == m:
                continue
            diagonal = (scores[M, i+1, j+1] + substitution(reference[i], read[j])
                        if i < n and j < m else NEG)
            ins = scores[I, i, j+1] if j < m else NEG
            delete = scores[D, i+1, j] if i < n else NEG
            scores[M, i, j] = max(diagonal, ins + opening(i, j+1, n, m),
                                   delete + opening(i+1, j, n, m))
            scores[I, i, j] = max(diagonal, ins - 2)
            scores[D, i, j] = max(diagonal, delete - 2)
    return scores


@njit(nogil=True, cache=True)
def maximize(reference, read, lo, hi, labels, cutoff):
    n, m = len(reference), len(read)
    prefix = forward(reference, read, 0, 0, n, m)
    suffix = backward(reference, read)
    optimum = max(prefix[M, n, m], prefix[I, n, m], prefix[D, n, m])
    threshold = NEG // 2 if cutoff < 0 else optimum - cutoff
    maxima = np.full((2, 3), NEG, np.int32)
    endpoints = np.full((2, 3, 5), -1, np.int32)
    right = np.full((3, m+1), NEG, np.int32)
    right_end = np.full((3, m+1), -1, np.int32)
    for j in range(m-1, -1, -1):
        diagonal = (suffix[M, hi+1, j+1] + substitution(reference[hi], read[j])
                    if read[j] != 4 else NEG)
        for state in STATES:
            right[state, j] = diagonal
            right_end[state, j] = j + 1
        for state in (M, I):
            cost = opening(hi, j+1, n, m) if state == M else -2
            if right[I, j+1] + cost > diagonal:
                right[state, j] = right[I, j+1] + cost
                right_end[state, j] = right_end[I, j+1]
    for flank in (lo, hi+1):
        for j in range(m+1):
            for state in (M, D):
                if state == M and (j == 0 or read[j-1] != 4):
                    continue
                score = prefix[state, flank, j] + suffix[state, flank, j]
                if score >= threshold and score > maxima[0, 2]:
                    for target in range(2):
                        maxima[target, 2] = score
                        endpoints[target, 2, 1] = flank
                        endpoints[target, 2, 2] = j
                        endpoints[target, 2, 3] = state
    for start in range(1, m):
        if read[start-1] == 4 or prefix[M, lo, start] + suffix[M, lo, start] < threshold:
            continue
        local = forward(reference, read, lo, start, hi, m-1)
        for end in range(start, m):
            for state in (M, D):
                score = prefix[M, lo, start] + local[state, hi-lo, end-start] + right[state, end]
                if score < threshold:
                    continue
                for target in range(2):
                    label = labels[target, start, end]
                    if score > maxima[target, label]:
                        maxima[target, label] = score
                        endpoints[target, label, 0] = start
                        endpoints[target, label, 1] = hi
                        endpoints[target, label, 2] = end
                        endpoints[target, label, 3] = state
                        endpoints[target, label, 4] = right_end[state, end]
    return optimum, maxima, endpoints, prefix, suffix


def trace_prefix(reference, read, scores, origin, end):
    start_i, start_j = origin
    i, j, state = end
    path = []
    while (i, j) != origin:
        a, b = i-start_i, j-start_j
        current = int(scores[state, a, b])
        if state == M:
            path.append((reference[i-1], read[j-1]))
            cost = int(substitution(ALPHABET.index(reference[i-1]), ALPHABET.index(read[j-1])))
            i, j = i-1, j-1
            previous = STATES
        elif state == I:
            path.append(("-", read[j-1]))
            cost = int(opening(i, j, len(reference), len(read)))
            j -= 1
            previous = (M, I)
        else:
            path.append((reference[i-1], "-"))
            cost = int(opening(i, j, len(reference), len(read)))
            i -= 1
            previous = (M, D)
        state = next(p for p in previous if scores[p, i-start_i, j-start_j]
                     + (-2 if state != M and p == state else cost) == current)
    assert state == M
    return path[::-1]


def trace_suffix(reference, read, scores, i, j, state):
    n, m = len(reference), len(read)
    path = []
    while (i, j) != (n, m):
        current = int(scores[state, i, j])
        if i < n and j < m and current == scores[M, i+1, j+1] + substitution(
                ALPHABET.index(reference[i]), ALPHABET.index(read[j])):
            path.append((reference[i], read[j]))
            i, j, state = i+1, j+1, M
        elif j < m and state != D and current == scores[I, i, j+1] + (
                -2 if state == I else opening(i, j+1, n, m)):
            path.append(("-", read[j]))
            j, state = j+1, I
        else:
            assert i < n and state != I and current == scores[D, i+1, j] + (
                -2 if state == D else opening(i+1, j, n, m))
            path.append((reference[i], "-"))
            i, state = i+1, D
    return path


def encode(s):
    return np.array([ALPHABET.index(b) for b in s], dtype=np.int8)


def make_targets(case):
    if case.get("target_kind") == "HEK3_DNA":
        return {"programmed_local": case["programmed_local"], "unedited_local": case["unedited_local"]}
    return {"native": case["native_peptide"], "exact": case["exact_programmed_coding"]}


def is_event(name, target, coding):
    if name == "native":
        return len(coding) == 3 * len(target) and all(
            TRANSLATION[coding[i:i+3]] == target[i//3] for i in range(0, len(coding), 3))
    return coding == target


def classify_intervals(read, targets):
    m = len(read)
    labels = np.zeros((2, m+1, m+1), np.int8)
    for start in range(m+1):
        labels[:, start, start] = 2
        n_position = read.find("N", start)
        if n_position >= 0:
            labels[:, start, n_position+1:] = 2
        for index, (name, target) in enumerate(targets.items()):
            end = start + len(target) * (3 if name == "native" else 1)
            if end <= m and labels[index, start, end] != 2:
                labels[index, start, end] = is_event(name, target, read[start:end])
    return labels


def solve(case, lattice, cutoff):
    reference, read = case["reference_sequence"], case["observed_sequence"]
    lo, hi = case["lo"], case["hi"]
    targets = make_targets(case)
    r, q = encode(reference), encode(read)
    optimum, maxima, endpoints, prefix, suffix = maximize(r, q, lo, hi, classify_intervals(read, targets), cutoff)
    assert optimum == suffix[M, 0, 0]
    results = {}
    for index, (name, target) in enumerate(targets.items()):
        witnesses = {}
        for outcome, value in enumerate(maxima[index]):
            if value == NEG:
                continue
            start, i, j, state, right_end = map(int, endpoints[index, outcome])
            if start < 0:
                path = trace_prefix(reference, read, prefix, (0, 0), (i, j, state))
                path += trace_suffix(reference, read, suffix, i, j, state)
            else:
                path = trace_prefix(reference, read, prefix, (0, 0), (lo, start, M))
                local = forward(r, q, lo, start, i, j)
                path += trace_prefix(reference, read, local, (lo, start), (i, j, state))
                path += [("-", b) for b in read[j:right_end-1]]
                path.append((reference[i], read[right_end-1]))
                path += trace_suffix(reference, read, suffix, i+1, right_end, M)
            aligned_ref, aligned_read = ("".join(x) for x in zip(*path))
            assert aligned_ref.replace("-", "") == reference and aligned_read.replace("-", "") == read
            assert lattice.rescore(aligned_ref, aligned_read) == value
            coding = lattice.project(aligned_ref, aligned_read, lo, hi)
            assert (2 if coding is None else int(is_event(name, target, coding))) == outcome
            witnesses[outcome] = {"score": int(value), "projection": coding,
                                  "aligned_reference": aligned_ref, "aligned_read": aligned_read}
        assert max(maxima[index]) == optimum
        results[name] = witnesses
    return int(optimum), results


def file_hash(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.repository / "extensions" / "correspondence_frontier"
    sys.path.insert(0, str(root.resolve()))
    import optimal_lattice as lattice

    args.output.mkdir(exist_ok=False)
    summary = {}
    for cohort in ("cftr", "hek3"):
        cases_path = root / "inputs" / f"{cohort}_REAL_ALIGNMENT_CASES.json"
        envelope_path = root / "expected" / f"{cohort}_lattice" / "OPTIMAL_ALIGNMENT_ENVELOPES.json"
        cases = json.loads(cases_path.read_bytes())["cases"]
        envelopes = json.loads(envelope_path.read_bytes())["cases"]
        assert len(cases) == len(envelopes)
        counts = Counter()
        with (args.output / f"{cohort}_witnesses.jsonl").open("x", encoding="utf-8") as stream:
            for case, envelope in zip(cases, envelopes, strict=True):
                for key in ("read_id", "reference", "run"):
                    assert case[key] == envelope[key]
                optimum, targets = solve(case, lattice, cutoff=-1)
                assert optimum == envelope["optimal_integer_objective"]
                for name, witnesses in targets.items():
                    labels = {outcome for outcome, w in witnesses.items() if w["score"] == optimum}
                    expected = set(map(int, envelope[f"{name}_values"]))
                    if envelope["some_optimal_projection_unresolved"]:
                        expected.add(2)
                    assert labels == expected, (case["read_id"], name, labels, expected)
                    if len(labels) == 1 and 2 not in labels:
                        label = next(iter(labels))
                        for delta in (0, 9, 18, 27, 36, 45, 54):
                            if all(optimum - w["score"] > delta for k, w in witnesses.items() if k != label):
                                counts[name, label, delta] += 1
                    counts["attaining_witnesses"] += len(witnesses)
                record = {key: case[key] for key in ("run", "read_id", "reference", "selection")}
                record.update(optimum=optimum, targets=targets)
                stream.write(json.dumps(record) + "\n")
        summary[cohort] = {
            "cases": len(cases), "targets": len(cases)*2,
            "inputs": {p.name: file_hash(p) for p in (cases_path, envelope_path)},
            "attaining_witnesses": counts["attaining_witnesses"],
            "certified": [{"target": name, "label": label, "delta": delta, "n": count}
                          for (name, label, delta), count in sorted(
                              (k, v) for k, v in counts.items() if isinstance(k, tuple))],
        }
        print(cohort, len(cases), counts["attaining_witnesses"], flush=True)
    (args.output / "panel_receipt.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
