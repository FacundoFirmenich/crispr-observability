"""Acquire exactly the ENA mate files fixed by the external design manifest."""

from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import urllib.request

from freeze import ROOT, sha


def acquire(record):
    url, expected, size = record["url"], record["md5"], record["bytes"]
    destination = ROOT / "fastq" / url.split("/")[-1]
    if not destination.exists():
        temp = destination.with_suffix(".part")
        with urllib.request.urlopen(url, timeout=50) as source, temp.open("wb") as output:
            while chunk := source.read(1024*1024):
                output.write(chunk)
        with temp.open("rb") as f:
            assert temp.stat().st_size == size and hashlib.file_digest(f, "md5").hexdigest() == expected
        temp.replace(destination)
    with destination.open("rb") as f:
        assert destination.stat().st_size == size and hashlib.file_digest(f, "md5").hexdigest() == expected
    return {"file": destination.name, **record, "sha256": sha(destination)}


if __name__ == "__main__":
    freeze = ROOT / "EXTERNAL_FREEZE.json"
    contract = json.loads(freeze.read_bytes())
    (ROOT / "fastq").mkdir(exist_ok=True)
    inputs = [x for s in contract["samples"] for x in s["files"]]
    assert len(inputs) == len({x["url"] for x in inputs}) == 24
    with ThreadPoolExecutor(4) as pool:
        verified = list(pool.map(acquire, inputs))
    (ROOT / "ACQUISITION.json").write_text(json.dumps(
        {"freeze_sha256": sha(freeze), "inputs": verified}, indent=2) + "\n")
    print(len(verified), sum(x["bytes"] for x in verified))
