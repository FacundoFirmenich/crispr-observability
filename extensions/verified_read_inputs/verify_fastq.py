"""Acquire and verify a pinned set of publicly archived FASTQ inputs."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import tempfile
from urllib.parse import urlsplit
from urllib.request import urlopen

MANIFEST = Path(__file__).with_name("keap1_ena.json")


def verify(path, entry):
    content = path.read_bytes()
    if len(content) != entry["bytes"] or hashlib.md5(content).hexdigest() != entry["md5"]:
        raise ValueError(f"ENA size or MD5 mismatch: {path.name}")
    digest = hashlib.sha256(content).hexdigest()
    if digest != entry["sha256"]:
        raise ValueError(f"SHA-256 mismatch: {path.name}")
    identifiers = set()
    with gzip.open(path, "rt", encoding="ascii") as stream:
        for header, seq, plus, qual in zip(stream, stream, stream, stream, strict=True):
            header, seq, plus, qual = (x.rstrip("\r\n") for x in (header,seq,plus,qual))
            if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Invalid FASTQ record in {path.name}")
            name = header[1:].split()[0]
            if not name.startswith(entry["accession"]+".") or name in identifiers:
                raise ValueError(f"Unexpected or duplicate identifier in {path.name}: {name}")
            identifiers.add(name)
    if len(identifiers) != entry["reads"]:
        raise ValueError(f"Unexpected FASTQ record count: {path.name}")
    return {"accession": entry["accession"], "reads": len(identifiers), "sha256": digest}


def acquire(directory, entry, offline):
    accession = entry["accession"]
    if not re.fullmatch(r"SRR\d+", accession):
        raise ValueError(f"Unexpected accession {accession}")
    path = directory / (accession+".fastq.gz")
    if path.is_file() or offline:
        return verify(path, entry)
    source = urlsplit(entry["url"])
    if source.scheme != "https" or source.netloc != "ftp.sra.ebi.ac.uk":
        raise ValueError(f"Unexpected FASTQ provider for {accession}")
    with tempfile.NamedTemporaryFile(dir=directory, suffix=".part", delete=False) as stage:
        temporary = Path(stage.name)
    try:
        with urlopen(entry["url"], timeout=45) as response:
            content = response.read(entry["bytes"]+1)
        temporary.write_bytes(content)
        result = verify(temporary, entry)
        temporary.replace(path)
        return result
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_bytes())
    verified = [acquire(args.directory, entry, args.offline) for entry in manifest["runs"]]
    print(json.dumps({"study": manifest["study"], "verified": verified}, indent=2))


if __name__ == "__main__":
    main()
