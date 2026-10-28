#!/usr/bin/env python3
"""Replay three P1 results from authentic raw records and frozen allele tables.

Standard library only. No CRISPResso executable is bundled, imported or called.
The target observer follows M1 separate_quality_geometry.py; the witnessed
coding contract follows M3 build_witness.py and M2 verify_m2.py. These rules
reconstruct an existing retrospective analysis, not a new biological cohort.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parent
COMP = str.maketrans("ACGTN", "TGCAN")
CODONS = dict(zip(("".join(c) for c in itertools.product("TCAG", repeat=3)),
                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
LEFT, RIGHT = "TGGATTATGCCT", "ATGAATATAGATACAGAAGC"
M1_EVENTS = ["CTT_plus_full_SE2", "native_local_peptide_all_observed_sequences"]
WITNESS_EVENTS = ["CTT_full_SE2_on_witnessed_observed_coding",
                  "native_local_peptide_on_witnessed_observed_coding"]
POLICIES = ["structural_only", "own_meanQ30", "common_raw_meanQ30", "common_and_own_meanQ30"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def verify_inputs():
    manifest = load(ROOT / "INPUT_MANIFEST.json")
    for item in manifest["files"]:
        path = ROOT / item["path"]
        require(path.is_file(), "Missing input: " + item["path"])
        require(path.stat().st_size == item["size"] and sha(path) == item["sha256"],
                "Input identity mismatch: " + item["path"])
    return manifest


def positions(sequence, motif):
    start = 0
    while True:
        start = sequence.find(motif, start)
        if start < 0:
            return
        yield start
        start += 1


def target_observer(oriented, exact):
    """M1 unique observed motif pair, no imputation; target interval only."""
    hits = []
    for sequence, quality, orientation in oriented:
        for a in positions(sequence, LEFT):
            for b in positions(sequence, RIGHT):
                if a + len(LEFT) <= b <= a + len(LEFT) + 350:
                    hits.append((sequence, quality, a, b, orientation))
    if len(hits) != 1:
        return {"labels": ("U", "U"), "qpass": False, "coding": None}
    sequence, quality, a, b, orientation = hits[0]
    core = sequence[a + len(LEFT):b]
    span = quality[a:b + len(RIGHT)]
    require(bool(span), "Empty target quality span")
    qpass = sum(ord(c) - 33 for c in span) >= 30 * len(span)
    if set(core) - set("ACGT"):
        return {"labels": ("U", "U"), "qpass": qpass, "coding": None}
    coding = core + RIGHT[:2]
    peptide = "".join(CODONS[coding[i:i + 3]] for i in range(0, len(coding), 3)) if len(coding) % 3 == 0 else None
    return {"labels": ("E" if core == exact else "N", "E" if peptide == "GTIKENIIFGVSYD" else "N"),
            "qpass": qpass, "coding": coding}


def gated_labels(info, policy, common):
    gate = (policy == "structural_only" or
            policy == "own_meanQ30" and info["qpass"] or
            policy == "common_raw_meanQ30" and common or
            policy == "common_and_own_meanQ30" and common and info["qpass"])
    return info["labels"] if gate else ("U", "U")


def project_coding(aligned_reference, aligned_read, reference):
    """M2 frozen alignment columns; immediate observed flanks must be ACGT.

    Include insertion at the left coding boundary, exclude insertion after the
    coding interval, and remove internal deletion columns. No bases are imputed.
    """
    require(len(aligned_reference) == len(aligned_read), "Unequal aligned lengths")
    require(aligned_reference.replace("-", "") == reference["sequence"], "Frozen alignment/reference mismatch")
    columns = [i for i, base in enumerate(aligned_reference) if base != "-"]
    lo, hi = reference["core_start_0"], reference["coding_end_exclusive_0"]
    if any(aligned_read[columns[i]] not in "ACGT" for i in (lo - 1, hi)):
        return None, ("U", "U")
    coding = aligned_read[columns[lo - 1] + 1:columns[hi - 1] + 1].replace("-", "")
    if not coding or set(coding) - set("ACGT"):
        return None, ("U", "U")
    peptide = "".join(CODONS[coding[i:i + 3]] for i in range(0, len(coding), 3)) if len(coding) % 3 == 0 else None
    return coding, ("E" if coding == reference["exact_programmed_core"] + "AT" else "N",
                    "E" if peptide == reference["native_peptide"] else "N")


def allele_map(path, reference):
    mapping, expected = {}, Counter()
    with zipfile.ZipFile(path) as archive:
        names = [x for x in archive.namelist() if not x.endswith("/")]
        require(len(names) == 1, "Expected exactly one detailed allele table")
        with archive.open(names[0]) as binary:
            for row in csv.DictReader(io.TextIOWrapper(binary, encoding="utf-8"), delimiter="\t"):
                key = row["Reference_Sequence"], row["Aligned_Sequence"]
                sequence = row["Aligned_Sequence"].replace("-", "")
                value = key, project_coding(*key, reference)
                require(sequence not in mapping or mapping[sequence] == value, "Ambiguous table sequence")
                mapping[sequence] = value
                abundance = int(row["#Reads"])
                require(abundance >= 0, "Negative allele abundance")
                expected[key] += abundance
    return mapping, expected


def fastq(path):
    with gzip.open(path, "rt", encoding="ascii", newline="") as f:
        while header := f.readline():
            sequence = f.readline().rstrip("\n")
            plus = f.readline()
            quality = f.readline().rstrip("\n")
            require(header.startswith("@") and plus.startswith("+") and bool(sequence), "Malformed FASTQ")
            require(len(sequence) == len(quality), "FASTQ sequence/quality length mismatch")
            require(not set(sequence) - set("ACGTN"), "Unexpected raw sequence alphabet")
            require(all(33 <= ord(c) <= 126 for c in quality), "Invalid Phred+33 quality")
            yield header.split()[0][1:], sequence, quality


def summarize(counter, n):
    require(sum(counter.values()) == n, "Label partition does not exhaust raw universe")
    e, ne, u = (counter[k] for k in "ENU")
    return {"N": n, "E": e, "NE": ne, "U": u,
            "L_exact": str(Fraction(e, n)), "H_exact": str(Fraction(e + u, n)),
            "L_pp": float(100 * Fraction(e, n)), "H_pp": float(100 * Fraction(e + u, n))}


def replay():
    design = load(ROOT / "metadata/RUNS.json")
    references = [load(ROOT / "metadata/REFERENCE.json"), load(ROOT / "metadata/REFERENCE_PROGRAMMED.json")]
    m1_rows, identities, witness_rows, domains = [], [], [], []
    aligned_totals = [0, 0]
    for run in design["runs"]:
        rid = run["run_accession"]
        mappings, expected = [], []
        for folder, reference in zip(("mutant", "programmed"), references):
            mapping, abundance = allele_map(ROOT / "inputs/alleles" / folder / (rid + ".zip"), reference)
            mappings.append(mapping)
            expected.append(abundance)
        used = [Counter(), Counter()]
        target_counts = {(trim, policy, ei): Counter() for trim in (0, 8) for policy in POLICIES for ei in range(2)}
        witness_counts = [Counter(), Counter()]
        domain, flags, ids = Counter(), Counter(), set()
        digest, n = hashlib.sha256(), 0
        raw_path = ROOT / "inputs/fastq" / (rid + ".fastq.gz")
        for read_id, sequence, quality in fastq(raw_path):
            require(read_id not in ids, "Duplicate read ID within run " + rid)
            ids.add(read_id)
            n += 1
            digest.update((sequence + "\n" + quality + "\n").encode("ascii"))
            common = sum(ord(c) - 33 for c in quality) >= 30 * len(quality)
            reverse, reverse_q = sequence.translate(COMP)[::-1], quality[::-1]
            info0 = None
            for trim in (0, 8):
                oriented = ([(sequence, quality, "forward"), (reverse, reverse_q, "reverse")] if trim == 0 else
                            [(sequence[:-trim], quality[:-trim], "forward"), (reverse[trim:], reverse_q[trim:], "reverse")])
                info = target_observer(oriented, references[0]["exact_programmed_core"])
                if trim == 0:
                    info0 = info
                for policy in POLICIES:
                    for ei, label in enumerate(gated_labels(info, policy, common)):
                        target_counts[trim, policy, ei][label] += 1
            labels, reason = ("U", "U"), "common_quality_fail"
            if common:
                target = info0["coding"] if info0["labels"][1] != "U" else None
                observed = []
                for ri, mapping in enumerate(mappings):
                    matches = []
                    for candidate in (sequence, reverse):
                        if candidate in mapping and mapping[candidate] not in matches:
                            matches.append(mapping[candidate])
                    require(len(matches) <= 1, "Ambiguous raw-to-frozen-table correspondence")
                    if matches:
                        key, value = matches[0]
                        used[ri][key] += 1
                        aligned_totals[ri] += 1
                        observed.append(value)
                    else:
                        observed.append((None, ("U", "U")))
                flags["target" + str(int(target is not None)) + "_mut" + str(int(observed[0][0] is not None)) +
                      "_prog" + str(int(observed[1][0] is not None))] += 1
                if target is not None and all(x[0] is not None for x in observed):
                    domain["triple_resolved"] += 1
                    if target == observed[0][0] == observed[1][0]:
                        labels = info0["labels"]
                        require(labels == observed[0][1] == observed[1][1], "Same-coding event label mismatch")
                        reason = "identical_observed_coding_witness"
                    else:
                        reason = "different_observed_coding"
                        domain["different_coding"] += 1
                        if info0["labels"][1] == observed[0][1][1] == observed[1][1][1]:
                            domain["equal_native_labels_different_coding"] += 1
                else:
                    reason = "missing_correspondence_witness"
            domain[reason] += 1 if reason != "different_observed_coding" else 0
            for ei, label in enumerate(labels):
                witness_counts[ei][label] += 1
        require(n == run["raw_reads"], "Raw-record count mismatch " + rid)
        require(digest.hexdigest() == run["sequence_quality_stream_sha256"], "Logical raw stream changed " + rid)
        require(used == expected, "Frozen allele abundances do not match authentic raw sequences " + rid)
        identities.append({"run_accession": rid, "N": n, "sequence_quality_stream_sha256": digest.hexdigest(),
                           "fastq_file_sha256": sha(raw_path)})
        for trim in (0, 8):
            for policy in sorted(POLICIES):
                for ei, event in enumerate(M1_EVENTS):
                    meta = {k: run[k] for k in ("run_accession", "figure", "condition", "source_replica", "source_filename")}
                    m1_rows.append({**meta, "trim": trim, "policy": policy, "observer": "target", "event": event,
                                    **summarize(target_counts[trim, policy, ei], n)})
        meta = {k: run[k] for k in ("run_accession", "condition", "source_replica")}
        for ei, event in enumerate(WITNESS_EVENTS):
            witness_rows.append({**meta, "event": event, "observer": "triple_observed_coding_witness",
                                 "universe_sha256": digest.hexdigest(), **summarize(witness_counts[ei], n)})
        domains.append({**meta, "N_raw": n, "counts": dict(domain), "resolved_intersections": dict(flags)})
        print(json.dumps({"run": rid, "raw_reads": n, "witness_U": witness_counts[0]["U"]}), flush=True)
    index = {(row["run_accession"], row["event"]): row for row in witness_rows}
    contrasts = []
    for treated in design["runs"]:
        if treated["condition"] != "PE treated":
            continue
        for control in ("mock", "untreated"):
            comparator = next(r for r in design["runs"] if r["source_replica"] == treated["source_replica"] and r["condition"] == control)
            for event in WITNESS_EVENTS:
                a, b = index[treated["run_accession"], event], index[comparator["run_accession"], event]
                lo = Fraction(a["L_exact"]) - Fraction(b["H_exact"])
                hi = Fraction(a["H_exact"]) - Fraction(b["L_exact"])
                contrasts.append({"preparation": treated["source_replica"], "control": control,
                                  "treated": treated["run_accession"], "comparator": comparator["run_accession"],
                                  "event": event, "observer": "triple_observed_coding_witness",
                                  "L_exact": str(lo), "H_exact": str(hi), "L_pp": float(100 * lo), "H_pp": float(100 * hi),
                                  "treated_counts": {k: a[k] for k in ("N", "E", "NE", "U")},
                                  "control_counts": {k: b[k] for k in ("N", "E", "NE", "U")}})
    require(sum(x["N"] for x in identities) == 177006, "Wrong total raw universe")
    require(sum(x["U"] for x in witness_rows if x["event"] == WITNESS_EVENTS[0]) == 6281, "Wrong unresolved witness total")
    return {"M1_TARGET_COUNTS.json": m1_rows, "RAW_UNIVERSE_IDENTITIES.json": identities,
            "WITNESS_BOUNDS.json": {"counts": witness_rows, "contrasts": contrasts, "domains": domains,
                                    "aligned_by_reference": aligned_totals}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/generated")
    args = parser.parse_args()
    out = args.output.resolve()
    for protected in (ROOT / "inputs", ROOT / "metadata", ROOT / "expected"):
        require(not out.is_relative_to(protected), "Output must not modify frozen inputs/expected files")
    started = time.monotonic()
    manifest = verify_inputs()
    outputs = replay()
    checks = {}
    # Expected files are read only after all three outputs have been computed.
    for name, value in outputs.items():
        require(value == load(ROOT / "expected" / name), "Original scientific output mismatch: " + name)
        checks[name] = True
    out.mkdir(parents=True, exist_ok=True)
    for name, value in outputs.items():
        dump(out / name, value)
    result = {"status": "PASS", "input_files_checked": len(manifest["files"]), "outputs_checked": checks,
              "raw_reads": 177006, "witnessed_records": 170725, "unresolved_records": 6281,
              "M1_target_rows": len(outputs["M1_TARGET_COUNTS.json"]), "witness_count_rows": 18,
              "witness_contrast_rows": 12, "aligned_by_reference": outputs["WITNESS_BOUNDS.json"]["aligned_by_reference"],
              "observers_recomputed_from_raw": True, "projection_from_frozen_alignments": True,
              "external_alignment_engine_rerun": False, "historical_training_rerun": False,
              "independent_human_verification": False, "new_biological_observations": False,
              "runtime_seconds": time.monotonic() - started,
              "scope": "Retrospective raw-record observer and frozen-table correspondence reconstruction; no molecular truth or function calibration."}
    dump(out / "VERIFICATION.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, StopIteration, zipfile.BadZipFile) as error:
        print("REPLAY FAILED: " + str(error), file=sys.stderr)
        sys.exit(1)
