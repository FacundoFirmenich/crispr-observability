"""Verify the clipped engine on every exact panel case and enumerate selected real paths."""

import argparse
import json
from pathlib import Path
import sys

from recompute_panel import (ALPHABET, D, I, M, backward, encode, is_event, make_targets)
from recompute_universes import DELTAS, process_sequence


def enumerate_paths(case, lattice, optimum):
    ref, read = case["reference_sequence"], case["observed_sequence"]
    n, m = len(ref), len(read)
    suffix = backward(encode(ref), encode(read))
    best = {name: {} for name in make_targets(case)}
    count = 0

    def walk(i, j, state, score, aligned_ref, aligned_read):
        nonlocal count
        if score + suffix[state, i, j] < optimum - 9:
            return
        if (i, j) == (n, m):
            assert score == lattice.rescore(aligned_ref, aligned_read)
            count += 1
            coding = lattice.project(aligned_ref, aligned_read, case["lo"], case["hi"])
            for name, target in make_targets(case).items():
                label = 2 if coding is None else int(is_event(name, target, coding))
                best[name][label] = max(best[name].get(label, score), score)
            return
        if i < n and j < m:
            walk(i+1, j+1, M, score+lattice.substitution(ref[i], read[j]),
                 aligned_ref+ref[i], aligned_read+read[j])
        if j < m and state != D:
            cost = -2 if state == I or i in (0, n) or j+1 == m else -20
            walk(i, j+1, I, score+cost, aligned_ref+"-", aligned_read+read[j])
        if i < n and state != I:
            cost = -2 if state == D or j in (0, m) or i+1 == n else -20
            walk(i+1, j, D, score+cost, aligned_ref+ref[i], aligned_read+"-")

    assert not set(ref+read) - set(ALPHABET)
    walk(0, 0, M, 0, "", "")
    return count, best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frontier = args.repository / "extensions/correspondence_frontier"
    sys.path.insert(0, str(frontier.resolve()))
    import optimal_lattice as lattice

    verified = 0
    enumerations = []
    for cohort in ("cftr", "hek3"):
        cases = json.loads((frontier / "inputs" / f"{cohort}_REAL_ALIGNMENT_CASES.json").read_bytes())["cases"]
        with (args.panel / f"{cohort}_witnesses.jsonl").open(encoding="utf-8") as stream:
            for index, (case, line) in enumerate(zip(cases, stream, strict=True)):
                witness = json.loads(line)
                for key in ("run", "read_id", "reference"):
                    assert case[key] == witness[key]
                targets = make_targets(case)
                optimum, deficits, _ = process_sequence(case["observed_sequence"],
                    encode(case["reference_sequence"]), case["lo"], case["hi"], targets)
                assert optimum == witness["optimum"]
                for name, values in zip(targets, deficits, strict=True):
                    for outcome, value in enumerate(values):
                        w = witness["targets"][name].get(str(outcome))
                        expected = None if w is None or optimum-w["score"] > DELTAS[-1] else optimum-w["score"]
                        assert value == expected
                    verified += 1
                is_example = case["read_id"] == "SRR29623270.2091" and case["reference"] == "programmed"
                if index < 4 or is_example:
                    count, maxima = enumerate_paths(case, lattice, optimum)
                    for name, found in maxima.items():
                        expected = {int(k): w["score"] for k, w in witness["targets"][name].items()
                                    if optimum-w["score"] <= 9}
                        assert found == expected
                    enumerations.append({"cohort": cohort, "read_id": case["read_id"],
                                         "reference": case["reference"], "paths": count})
    assert verified == 928 and len(enumerations) == 9
    receipt = {"targets_verified": verified, "cutoff": DELTAS[-1],
               "complete_path_enumeration_cutoff": 9, "enumerations": enumerations,
               "enumerated_paths": sum(row["paths"] for row in enumerations)}
    args.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(receipt, flush=True)


if __name__ == "__main__":
    main()
