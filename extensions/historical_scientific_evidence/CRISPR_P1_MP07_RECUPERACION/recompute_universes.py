"""Exact target certificate survival at delta 0 through 54 in CFTR and HEK3."""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import io
from itertools import islice
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from recompute_panel import NEG, classify_intervals, encode, file_hash, make_targets, maximize


DELTAS = (0, 9, 18, 27, 36, 45, 54)
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
LEDGER_HASHES = {
    "cftr": "48e897dcf4e848f7c1327caf400119387dedf917438845f01924c95903ac8c3c",
    "hek3": "e23d7f441938c020c368838257baca1c63ce861e4e023572bfac50b285b72f56",
}


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def canonical(sequence):
    return min(sequence, rc(sequence))


def allele_orientations(path, reference):
    options = defaultdict(set)
    weights = Counter()
    with ZipFile(path) as archive:
        members = [entry for entry in archive.infolist() if not entry.is_dir()]
        assert len(members) == 1
        with io.TextIOWrapper(archive.open(members[0]), encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                assert row["Reference_Sequence"].replace("-", "") == reference
                sequence = row["Aligned_Sequence"].replace("-", "")
                key = canonical(sequence)
                options[key].add(sequence)
                weights[key] += int(row["#Reads"])
    return options, weights


def load_inputs(repository, ledger_dir):
    frontier = repository / "extensions" / "correspondence_frontier" / "inputs"
    for cohort in ("cftr", "hek3"):
        cases_path = frontier / f"{cohort}_REAL_ALIGNMENT_CASES.json"
        cases = json.loads(cases_path.read_bytes())["cases"]
        specs = {}
        for case in cases:
            name = case["reference"]
            spec = (case["reference_sequence"], case["lo"], case["hi"], make_targets(case))
            assert specs.setdefault(name, spec) == spec
        ledger = ledger_dir / f"{cohort.upper()}_COMPLETE_RECORD_LEDGER.csv.gz"
        assert file_hash(ledger) == LEDGER_HASHES[cohort]
        fastq_dir = (repository / "raw_replay/inputs/fastq" if cohort == "cftr"
                     else frontier / "hek3/inputs")
        records = defaultdict(list)
        with gzip.open(ledger, "rt", encoding="ascii", newline="") as stream:
            for row in csv.DictReader(stream):
                records[row["run"]].append(row)
        for run, rows in records.items():
            fastq = fastq_dir / f"{run}.fastq.gz"
            inputs = {fastq.name: file_hash(fastq), ledger.name: LEDGER_HASHES[cohort]}
            seen_ids = set()
            with gzip.open(fastq, "rt", encoding="ascii") as source:
                for row in rows:
                    header, sequence, plus, quality = (source.readline().strip() for _ in range(4))
                    assert header.startswith("@") and plus.startswith("+")
                    assert (header[1:].split()[0], sequence, quality) == (
                        row["read_id"], row["raw_sequence"], row["raw_quality"])
                    assert row["read_id"] not in seen_ids and len(sequence) == len(quality)
                    assert bool(int(row["Q30"])) == (sum(ord(q)-33 for q in quality) >= 30*len(quality))
                    seen_ids.add(row["read_id"])
                assert not source.readline()
            tables = {}
            for name, (ref, _, _, _) in specs.items():
                if cohort == "cftr":
                    table = repository / "raw_replay/inputs/alleles" / name / f"{run}.zip"
                else:
                    mode = "M_fw" if name == "reference" else "P_fw"
                    table = frontier / "hek3/runs" / mode / f"CRISPResso_on_{run}" / "Alleles_frequency_table.zip"
                options, weights = allele_orientations(table, ref)
                observed = Counter(canonical(row["raw_sequence"]) for row in rows if row["Q30"] == "1")
                assert all(observed[key] == count for key, count in weights.items()), (run, name)
                inputs[table.relative_to(repository).as_posix()] = file_hash(table)
                tables[name] = options
            differences = 0
            for case in (c for c in cases if c["run"] == run):
                sequence = case["observed_sequence"]
                options = tables[case["reference"]].get(canonical(sequence), set())
                differences += options != {sequence}
            yield cohort, run, specs, rows, tables, inputs, differences


def process_sequence(sequence, reference, lo, hi, targets):
    optimum, maxima, endpoints, prefix, suffix = maximize(
        reference, encode(sequence), lo, hi, classify_intervals(sequence, targets), DELTAS[-1])
    assert optimum == suffix[0, 0, 0]
    assert np.all(maxima.max(axis=1) == optimum)
    deficits = [[None if value == NEG else int(optimum - value) for value in values]
                for values in maxima]
    assert all(d is None or (d in DELTAS) for values in deficits for d in values)
    return int(optimum), deficits, endpoints.tolist()


def chunks(values, size):
    iterator = iter(values)
    while batch := list(islice(iterator, size)):
        yield batch


def label_at(deficits, delta):
    admitted = [outcome for outcome, value in enumerate(deficits) if value is not None and value <= delta]
    return ("NE", "E", "U")[admitted[0]] if len(admitted) == 1 else "U"


def run_reference(output, cohort, run, name, spec, rows, orientations, pool, inputs, panel_differences):
    ref, lo, hi, targets = spec
    basename = f"{run}_{name}"
    final_rows = output / f"{basename}.jsonl.gz"
    final_receipt = output / f"{basename}.json"
    if final_rows.exists() and final_receipt.exists():
        receipt = json.loads(final_receipt.read_bytes())
        assert receipt["inputs"] == inputs and receipt["records_sha256"] == file_hash(final_rows)
        return receipt
    assert not final_rows.exists() and not final_receipt.exists()
    reference = encode(ref)
    grouped = defaultdict(list)
    statuses = Counter()
    table = {}
    for row in rows:
        read_id, raw = row["read_id"], row["raw_sequence"]
        options = orientations.get(canonical(raw), set())
        if row["Q30"] != "1":
            status = "quality_fail"
        elif not options:
            status = "orientation_absent"
        elif len(options) != 1:
            status = "orientation_ambiguous"
        else:
            sequence = next(iter(options))
            assert "N" not in sequence
            grouped[sequence].append(read_id)
            status = "evaluated"
        statuses[status] += 1
        table[read_id] = {"id": read_id, "status": status}

    def calculate(sequence):
        return process_sequence(sequence, reference, lo, hi, targets)

    counts = Counter()
    print(f"START {cohort} {run} {name} raw={len(rows)} families={len(grouped)}", flush=True)
    for batch in chunks(grouped, 256):
        for sequence, (optimum, deficits, endpoints) in zip(batch, pool.map(calculate, batch), strict=True):
            read_ids = grouped[sequence]
            for read_id in read_ids:
                raw_entry = table[read_id]
                raw_entry.update(score=optimum, deficits=deficits, endpoints=endpoints,
                                 orientation_sequence_sha256=__import__("hashlib").sha256(sequence.encode()).hexdigest())
            for target, target_deficits in zip(targets, deficits, strict=True):
                for delta in DELTAS:
                    counts[target, delta, label_at(target_deficits, delta)] += len(read_ids)
    for target in targets:
        for delta in DELTAS:
            counts[target, delta, "U"] += len(rows) - statuses["evaluated"]
            assert sum(counts[target, delta, k] for k in ("E", "NE", "U")) == len(rows)
    temp_rows = final_rows.with_suffix(".tmp")
    with temp_rows.open("xb") as raw_stream:
        with gzip.GzipFile(fileobj=raw_stream, mode="wb", filename="", mtime=0) as compressed:
            for entry in table.values():
                compressed.write((json.dumps(entry, separators=(",", ":")) + "\n").encode())
    records_hash = file_hash(temp_rows)
    receipt = {
        "cohort": cohort, "run": run, "reference": name, "raw_N": len(rows), "statuses": dict(statuses),
        "unique_oriented_sequences": len(grouped), "panel_orientation_differences_across_both_references": panel_differences,
        "targets": list(targets), "inputs": inputs, "records_sha256": records_hash,
        "rows": [{"target": target, "delta": delta,
                  **{k: counts[target, delta, k] for k in ("E", "NE", "U")}}
                 for target in targets for delta in DELTAS],
    }
    with gzip.open(temp_rows, "rt", encoding="ascii") as readback:
        observed = Counter()
        for entry, expected in zip(readback, table.values(), strict=True):
            got = json.loads(entry)
            assert got == expected
            for index, target in enumerate(targets):
                for delta in DELTAS:
                    label = label_at(got["deficits"][index], delta) if got["status"] == "evaluated" else "U"
                    observed[target, delta, label] += 1
        assert observed == counts
    temp_receipt = final_receipt.with_suffix(".tmp")
    temp_receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    temp_rows.replace(final_rows)
    temp_receipt.replace(final_receipt)
    print(f"VERIFIED {cohort} {run} {name} {dict(statuses)}", flush=True)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--ledgers", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=True)
    receipts = []
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        for cohort, run, specs, rows, tables, inputs, differences in load_inputs(args.repository, args.ledgers):
            for name, spec in specs.items():
                receipts.append(run_reference(args.output, cohort, run, name, spec, rows, tables[name],
                                              pool, inputs, differences))
    assert len(receipts) == 24 and sum(r["raw_N"] for r in receipts) == 2 * 395936
    final = args.output / "full_receipt.json"
    final.write_text(json.dumps({"status": "COMPLETE_DELTA_0_TO_54", "runs_references": receipts}, indent=2)
                     + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
