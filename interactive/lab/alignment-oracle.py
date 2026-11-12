"""Independent fixture oracle. Read small public BAMs sequentially with stdlib."""
import gzip
import hashlib
import json
import re
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = "https://github.com/GMOD/bam-js/tree/ca0836f866c65e61c70cc89199cf93032c27de1f/test/data"


def read_bam(path):
    b = gzip.decompress(path.read_bytes())
    assert b[:4] == b"BAM\x01"
    header_length, = struct.unpack_from("<i", b, 4)
    header = b[8:8 + header_length].decode().rstrip("\n")
    off = 8 + header_length
    n_ref, = struct.unpack_from("<i", b, off)
    off += 4
    refs = []
    for _ in range(n_ref):
        n, = struct.unpack_from("<i", b, off)
        refs.append(b[off + 4:off + 3 + n].decode())
        off += 8 + n
    rows = []
    while off < len(b):
        size, = struct.unpack_from("<i", b, off)
        entry = off + 4
        ref_id, pos, bin_mq_nl, flag_nc, seq_len = struct.unpack_from("<iiIIi", b, entry)
        off = entry + size
        flag, count = flag_nc >> 16, flag_nc & 65535
        if flag & 4:
            continue
        name_len = bin_mq_nl & 255
        name_end = entry + 32 + name_len
        name = b[entry + 32:name_end - 1].decode()
        packed = struct.unpack_from(f"<{count}I", b, name_end)
        ops = [(item >> 4, "MIDNSHP=X"[item & 15]) for item in packed]
        cigar = "".join(f"{length}{op}" for length, op in ops)
        seq_start = name_end + 4 * count
        half = (seq_len + 1) // 2
        seq_codes = b[seq_start:seq_start + half]
        seq = "".join("=ACMGRSVTWYHKDBN"[(seq_codes[i // 2] >> (4 if not i % 2 else 0)) & 15] for i in range(seq_len))
        quals = b[seq_start + half:seq_start + half + seq_len]
        qual = "".join(chr(x + 33) for x in quals) if quals and quals[0] != 255 else "*"
        rows.append({"name": name, "flag": flag, "mapq": (bin_mq_nl >> 8) & 255, "ref": refs[ref_id],
                     "start": pos, "cigar": cigar, "sequence": seq or "*", "quality": qual})
    return header, rows


def ref_end(row):
    return row["start"] + sum(int(n) for n, op in re.findall(r"(\d+)(\D)", row["cigar"]) if op in "MDN=X")


def oracle(filename, threshold, dense=False):
    header, rows = read_bam(ROOT / filename)
    first = next((row for row in rows if "S" in row["cigar"] and row["flag"] & 16), rows[0])
    if dense:
        first = max(rows, key=lambda r: sum(
            candidate["ref"] == r["ref"] and candidate["start"] <= r["start"] < ref_end(candidate)
            for candidate in rows))
    ref = first["ref"]
    lo = first["start"]
    hi = min(lo + 200, ref_end(first))
    selected = [r for r in rows if r["ref"] == ref and r["start"] < hi and ref_end(r) > lo]
    exact = [0] * (hi - lo)
    for row in selected:
        if row["mapq"] == 255 or row["mapq"] < threshold or row["flag"] & (256 | 2048):
            continue
        pos = row["start"]
        for length, op in re.findall(r"(\d+)(\D)", row["cigar"]):
            length = int(length)
            if op in "M=X":
                for site in range(max(pos, lo), min(pos + length, hi)):
                    exact[site - lo] += 1
            if op in "MDN=X":
                pos += length
    sam = header + "\n" + "".join("\t".join(map(str, (
        r["name"], r["flag"], r["ref"], r["start"] + 1, r["mapq"], r["cigar"],
        "*", 0, 0, r["sequence"], r["quality"]))) + "\n" for r in selected)
    sam_name = filename.replace(".bam", ".dense.sam" if dense else ".sam")
    (ROOT / sam_name).write_text(sam, encoding="utf-8", newline="\n")
    return {"bam": filename, "sam": sam_name, "bamSHA256": hashlib.sha256((ROOT / filename).read_bytes()).hexdigest(),
            "ref": ref, "lo": lo, "hi": hi, "threshold": threshold, "records": selected, "depth": exact}


if __name__ == "__main__":
    checks = [oracle("1000genomes_hg00096_chr1.bam", 20), oracle("CHM1_pacbio_clip2.bam", 10),
              oracle("CHM1_pacbio_clip2.bam", 10, dense=True)]
    (ROOT / "ALIGNMENT_FIXTURE_ORACLE.json").write_text(json.dumps({"source": BASE, "checks": checks}) + "\n")
    for case in checks:
        print(case["bam"], case["ref"], case["lo"], case["hi"], len(case["records"]), max(case["depth"]))
