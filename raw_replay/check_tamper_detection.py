#!/usr/bin/env python3
"""Negative controls on temporary copies of authentic inputs; never modify ROOT."""
from pathlib import Path
import gzip
import hashlib
import io
import json
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.monotonic()
    source_manifest_hash = digest(ROOT / "INPUT_MANIFEST.json")
    source_manifest = json.loads((ROOT / "INPUT_MANIFEST.json").read_text())
    results = []
    with tempfile.TemporaryDirectory(prefix="p1-raw-negative-") as scratch:
        temp = Path(scratch) / "package"
        shutil.copytree(ROOT, temp, ignore=shutil.ignore_patterns("results", "*.stdout", "*.stderr", "__pycache__"))
        original_manifest = (temp / "INPUT_MANIFEST.json").read_bytes()

        def restore(relative):
            shutil.copyfile(ROOT / relative, temp / relative)
            (temp / "INPUT_MANIFEST.json").write_bytes(original_manifest)

        def refresh_transport_hash(relative):
            manifest = json.loads(original_manifest)
            item = next(x for x in manifest["files"] if x["path"] == relative)
            item.update(size=(temp / relative).stat().st_size, sha256=digest(temp / relative))
            (temp / "INPUT_MANIFEST.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")

        def expect_rejection(label, diagnostic):
            t = time.monotonic()
            run = subprocess.run([sys.executable, str(temp / "run_raw_replay.py"), "--output", str(Path(scratch) / "generated")],
                                 text=True, capture_output=True, timeout=290)
            if run.returncode == 0 or diagnostic not in run.stderr:
                raise AssertionError(label + " did not reject as required: " + run.stderr)
            results.append({"case": label, "status": "REJECTED_AS_EXPECTED", "exit_code": run.returncode,
                            "diagnostic": run.stderr.strip(), "runtime_seconds": time.monotonic() - t})

        raw = "inputs/fastq/SRR29623274.fastq.gz"
        lines = gzip.decompress((temp / raw).read_bytes()).splitlines(keepends=True)
        base = b"C" if lines[1][:1] != b"C" else b"A"
        lines[1] = base + lines[1][1:]
        (temp / raw).write_bytes(gzip.compress(b"".join(lines), mtime=0))
        expect_rejection("authentic_raw_base_changed", "Input identity mismatch")
        refresh_transport_hash(raw)
        expect_rejection("raw_base_changed_even_after_file_hash_update", "Logical raw stream changed")
        restore(raw)

        table = "inputs/alleles/mutant/SRR29623274.zip"
        with zipfile.ZipFile(temp / table) as archive:
            name = archive.namelist()[0]
            text = archive.read(name).decode("utf-8")
        records = text.splitlines(keepends=True)
        header = records[0].rstrip("\r\n").split("\t")
        row = records[1].rstrip("\r\n").split("\t")
        column = header.index("#Reads")
        row[column] = str(int(row[column]) + 1)
        records[1] = "\t".join(row) + "\n"
        with zipfile.ZipFile(temp / table, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr(name, "".join(records))
        refresh_transport_hash(table)
        expect_rejection("frozen_allele_abundance_changed_even_after_file_hash_update", "Frozen allele abundances do not match")
        restore(table)

        for relative, change in [
            ("expected/M1_TARGET_COUNTS.json", lambda x: x[0].update(E=x[0]["E"] + 1)),
            ("expected/WITNESS_BOUNDS.json", lambda x: x["contrasts"][0].update(L_exact="0")),
        ]:
            expected = json.loads((temp / relative).read_text())
            change(expected)
            (temp / relative).write_text(json.dumps(expected, ensure_ascii=False, indent=2) + "\n")
            refresh_transport_hash(relative)
            expect_rejection(relative + "_changed_after_file_hash_update", "Original scientific output mismatch")
            restore(relative)

    for item in source_manifest["files"]:
        if digest(ROOT / item["path"]) != item["sha256"]:
            raise AssertionError("Original source package changed: " + item["path"])
    if digest(ROOT / "INPUT_MANIFEST.json") != source_manifest_hash:
        raise AssertionError("Original manifest changed")
    result = {"status": "PASS", "cases": results, "original_inputs_unchanged": True,
              "runtime_seconds": time.monotonic() - started,
              "scope": "Temporary corruption controls on authentic raw/frozen/golden files; no new biological evidence and no external alignment-engine execution."}
    out = ROOT / "results/TAMPER_VERIFICATION.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
