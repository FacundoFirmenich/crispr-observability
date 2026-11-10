"""Recover actual HG002 reads with an external, phased marker outside the indel."""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

import pysam

OUT = Path("giab-evidence")
OUT.mkdir(exist_ok=True)
BASE = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/"
PREFIX = "GRCh38_HG2-T2TQ100-V1.1_smvar"
BAM = "https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"


def fetch(url, name, md5=None):
    path = OUT / name
    urllib.request.urlretrieve(url, path)
    if md5:
        with path.open("rb") as stream:
            assert hashlib.file_digest(stream, "md5").hexdigest() == md5
    return path


def reference(start, end):
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr20;start={start};end={end}"
    path = fetch(url, f"reference_chr20_{start}_{end}.json")
    seq = json.loads(path.read_bytes())["dna"].upper()
    assert len(seq) == end - start and set(seq) <= set("ACGT")
    return seq


def main():
    fetch(BASE + "NIST_HG002_DraftBenchmark_defrabbV0.019-20241113_README.md", "BENCHMARK_README.md")
    fetch(BASE + "checksum.md5", "original_checksum.md5")
    vcf = fetch(BASE + PREFIX + ".vcf.gz", "benchmark.vcf.gz", "a77065ecb99ece71939be763d7811209")
    bed = fetch(BASE + PREFIX + ".benchmark.bed", "benchmark.bed", "77a7f0b49d8bbbab16d346cf1a086ba4")
    intervals = [tuple(map(int, row.split()[1:3])) for row in bed.read_text().splitlines()
                 if row.startswith("chr20\t")]
    records = []
    with pysam.VariantFile(vcf) as variants:
        for v in variants:
            if v.chrom != "chr20":
                continue
            sample = next(iter(v.samples.values()))
            records.append({"pos": v.start, "ref": v.ref, "alts": v.alts,
                            "gt": sample["GT"], "phased": sample.phased,
                            "filter": list(v.filter), "vcf": str(v).strip()})
    pairs = []
    for i, v in enumerate(records):
        if not v["phased"] or v["gt"] not in ((0, 1), (1, 0)) or len(v["alts"] or ()) != 1:
            continue
        alt = v["alts"][0]
        if min(len(v["ref"]), len(alt)) != 1 or not 2 <= abs(len(v["ref"]) - len(alt)) <= 10:
            continue
        neighbors = [r for r in records[max(0, i-20):i+21]
                     if abs(r["pos"] - v["pos"]) <= 110 and r is not v]
        if len(neighbors) != 1:
            continue
        s = neighbors[0]
        if not (s["phased"] and s["gt"] in ((0, 1), (1, 0)) and
                len(s["ref"]) == 1 and len(s["alts"] or ()) == 1 and len(s["alts"][0]) == 1 and
                25 <= abs(s["pos"] - v["pos"]) <= 60):
            continue
        lo = min(v["pos"], s["pos"]) - 15
        hi = max(v["pos"] + len(v["ref"]), s["pos"] + 1) + 15
        if not any(a <= lo and hi <= b for a, b in intervals):
            continue
        pairs.append({"indel": v, "marker": s})
        if len(pairs) == 12:
            break
    assert pairs, "No candidate passes the outcome-independent selection rule"
    contract = {
        "scope": "External HG002 phased-marker control. Truth is conditional on the assembly phase, marker assignment and source identity. It does not adjudicate HEK3 or CFTR molecules.",
        "selection": "First 12 chr20 biallelic phased heterozygous indels of 2-10 bp with exactly one other variant within 110 bp, a phased heterozygous SNP 25-60 bp away. Joint core and marker plus 15 bp must be inside the benchmark BED. Selection precedes reading BAM observations.",
        "marker": "Unique exact 21-base marker match outside the indel, minimum base quality 30 across all 21 bases. Marker phase supplies the indel label. BAM target bases do not supply truth.",
        "reads": "All regional primary mapped nonduplicate QC-pass alignments with MAPQ >=30, keeping each mate as a distinct observation. Counts are reads, not independent molecules. No target sequence or predicted label filter.",
        "reference_source": "UCSC hg38 sequence API, exact response bytes retained",
        "benchmark": BASE, "bam": BAM, "pairs": pairs,
    }
    (OUT / "SELECTION_FREEZE.json").write_text(json.dumps(contract, indent=2) + "\n")
    index = fetch(BAM + ".bai", "input.bam.bai")
    counts = Counter()
    with pysam.AlignmentFile(BAM, "rb", index_filename=str(index)) as bam, \
         gzip.open(OUT / "reads.jsonl.gz", "wt") as output:
        (OUT / "bam_header.sam").write_text(str(bam.header))
        for pair in pairs:
            v, s = pair["indel"], pair["marker"]
            start, end = v["pos"] - 170, v["pos"] + len(v["ref"]) + 170
            ref = reference(start, end)
            p, q = v["pos"] - start, s["pos"] - start
            assert ref[p:p+len(v["ref"])] == v["ref"] and ref[q] == s["ref"]
            marker0 = ref[q-10:q+11]
            marker1 = marker0[:10] + s["alts"][0] + marker0[11:]
            local_lo, local_hi = p-5, p+len(v["ref"])+5
            target = ref[local_lo:p] + v["alts"][0] + ref[p+len(v["ref"]):local_hi]
            pair.update({"start": start, "end": end, "reference": ref,
                         "local_lo": local_lo, "local_hi": local_hi, "target": target,
                         "marker0": marker0, "marker1": marker1})
            seen = set()
            for read in bam.fetch("chr20", min(v["pos"], s["pos"])-11, max(v["pos"], s["pos"])+12):
                counts["fetched"] += 1
                if read.flag & (4 | 256 | 512 | 1024 | 2048) or read.mapping_quality < 30:
                    counts["mapping_excluded"] += 1
                    continue
                seq, quals = read.query_sequence, read.query_qualities
                hits = [(a, seq.find(marker)) for a, marker in enumerate((marker0, marker1)) if seq.count(marker)]
                if len(hits) != 1 or seq.count((marker0, marker1)[hits[0][0]]) != 1:
                    counts["marker_unresolved"] += 1
                    continue
                allele, offset = hits[0]
                if quals is None or min(quals[offset:offset+21]) < 30:
                    counts["marker_quality_excluded"] += 1
                    continue
                haplotype = s["gt"].index(allele)
                truth = "E" if v["gt"][haplotype] == 1 else "NE"
                identity = (read.query_name, read.is_read1, read.is_read2)
                assert identity not in seen
                seen.add(identity)
                counts["tagged_" + truth] += 1
                row = {"site": v["pos"] + 1, "read_id": read.query_name, "flag": read.flag,
                       "position": read.reference_start, "cigar": read.cigarstring,
                       "mapq": read.mapping_quality, "sequence": seq, "quality": read.qual,
                       "marker_allele": allele, "marker_offset": offset, "truth": truth}
                output.write(json.dumps(row) + "\n")
    (OUT / "CASES.json").write_text(json.dumps(pairs, indent=2) + "\n")
    (OUT / "COUNTS.json").write_text(json.dumps(counts, indent=2) + "\n")
    manifest = {}
    for path in OUT.iterdir():
        if path.is_file():
            with path.open("rb") as f:
                manifest[path.name] = hashlib.file_digest(f, "sha256").hexdigest()
    (OUT / "SHA256.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"sites": len(pairs), "counts": counts}), flush=True)


if __name__ == "__main__":
    main()
