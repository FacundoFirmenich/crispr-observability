"""Join an actual CRISPR-A native artifact to the frozen MP07 HEK3 ledgers."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
MP07 = ROOT.parent/"CRISPR_P1_MP07_RECUPERACION"
RUN = "SRR10285861"
DELTAS = (0,9,18,27,36,45,54)


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def verify_artifact(path):
    assert (path/"nextflow.exit").read_text().strip() == "0"
    for line in (path/"output.sha256").read_text().splitlines():
        identity,name = line.split("  ",1)
        assert name.startswith("../evidence/")
        p = (path/name.removeprefix("../evidence/")).resolve()
        assert p.is_relative_to(path.resolve())
        assert digest(p) == identity


def covers(pos, cigar, low, high):
    i = pos-1
    covered = 0
    ops = re.findall(r"(\d+)([MIDNSHP=X])",cigar)
    assert "".join(n+k for n,k in ops) == cigar
    for n,k in ops:
        n = int(n)
        if k in "M=X":
            covered += max(0,min(i+n,high)-max(i,low))
        if k in "MDN=X":
            i += n
    return covered == high-low


def observed_projection(pos, cigar, query, lo, hi):
    ref_pos, read_pos = pos-1, 0
    flank = {}
    segments = []
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])",cigar):
        n = int(length)
        if op in "M=X":
            for k in (lo-1,hi):
                if ref_pos<=k<ref_pos+n:
                    flank[k] = query[read_pos+k-ref_pos]
            start = max(lo,ref_pos)
            end = min(hi,ref_pos+n)
            if start<end:
                segments.append(query[read_pos+start-ref_pos:read_pos+end-ref_pos])
        elif op=="I" and lo<=ref_pos<hi:
            segments.append(query[read_pos:read_pos+n])
        if op in "MDN=X":
            ref_pos += n
        if op in "MIS=X":
            read_pos += n
    assert read_pos == len(query)
    coding = "".join(segments)
    if set(flank)!={lo-1,hi} or not coding or set(coding+"".join(flank.values()))-set("ACGT"):
        return None
    return coding


def label(row, delta):
    if row["status"] != "evaluated":
        return "U"
    flags = [i for i,d in enumerate(row["deficits"][0]) if d is not None and d <= delta]
    return ("NE","E","U")[flags[0]] if len(flags) == 1 else "U"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    a = parser.parse_args()
    verify_artifact(a.evidence)
    native = a.evidence/"native-results"
    rule = json.loads((MP07/"i2/extensions/correspondence_frontier/inputs/hek3/RULE_FREEZE.json").read_bytes())
    reference = a.evidence/f"{RUN}_reference-correctOrient.fasta"
    seq = "".join(line.strip() for line in reference.read_text().splitlines() if not line.startswith(">")).upper()
    assert seq == rule["reference"], "Native orientation needs a verified coordinate transform"
    with (native/f"{RUN}_template-reads.csv").open(encoding="ascii",newline="") as f:
        rows = csv.reader(f)
        assert len(next(rows)) == 2
        ids = [row[1] for row in rows]
    assert all(re.fullmatch(r"SRR10285861\.\d+",x) for x in ids)
    positive = set(ids)
    assert len(positive)==len(ids)
    alignments = defaultdict(list)
    raw_path = MP07/"i2/extensions/correspondence_frontier/inputs/hek3/inputs/SRR10285861.fastq.gz"
    assert digest(raw_path) == "e89a193fe84eda143cbb2c1df154e44157d467a578aeaf077b13f9703d0a816a"
    raw_sequences = {}
    with gzip.open(raw_path,"rt",encoding="ascii") as f:
        while header := f.readline().strip():
            sequence,plus,quality = (f.readline().strip() for _ in range(3))
            assert header.startswith("@") and plus.startswith("+") and len(sequence)==len(quality)
            name = header[1:].split()[0]
            assert name not in raw_sequences
            raw_sequences[name] = sequence
    sequence_identity = Counter()
    with gzip.open(a.evidence/"primary-alignments.tsv.gz","rt",encoding="ascii") as f:
        for line in f:
            name,flag,ref,pos,mapq,cigar,query = line.rstrip("\n").split("\t")
            assert int(flag)&2308 == 0
            raw_sequence = raw_sequences[name]
            oriented = raw_sequence.translate(str.maketrans("ACGTN","TGCAN"))[::-1] if int(flag)&16 else raw_sequence
            identity = "exact_raw_oriented" if query==oriented else "trimmed_substring_of_raw" if query in oriented else "other_sequence_change"
            sequence_identity[identity] += 1
            alignments[name].append((int(pos),cigar,query))
    assert positive <= set(alignments)
    index = MP07/"universes_MP07"
    checks = {(r["run"],r["reference"]):r for r in json.loads((index/"full_receipt.json").read_bytes())["runs_references"]}
    paths = [index/f"{RUN}_{k}.jsonl.gz" for k in ("reference","programmed_reference")]
    for path,k in zip(paths,("reference","programmed_reference"),strict=True):
        assert digest(path)==checks[RUN,k]["records_sha256"]
    counts, admission = Counter(), Counter()
    identifiers = set()
    out_rows = a.output.with_suffix(".csv.gz")
    with gzip.open(paths[0],"rt",encoding="ascii") as x,gzip.open(paths[1],"rt",encoding="ascii") as y, \
         gzip.open(out_rows,"wt",encoding="ascii",newline="") as output:
        w = csv.writer(output)
        w.writerow(("id","crispra_native_template_call","unique_local_span","crispra_span_conditioned","crispra_exact_tract",
                    *(f"mp07_delta_{d}" for d in DELTAS)))
        for u,v in zip(x,y,strict=True):
            l,r = json.loads(u),json.loads(v)
            name = l["id"]
            assert name==r["id"] and name not in identifiers
            identifiers.add(name)
            aln = alignments.get(name,[])
            span = len(aln)==1 and covers(*aln[0][:2],rule["local_lo"],rule["local_hi"])
            observed = observed_projection(*aln[0],rule["local_lo"],rule["local_hi"]) if len(aln)==1 else None
            projected = "U" if observed is None else "E" if observed==rule["programmed_local"] else "NE"
            p = name in positive
            native_label = "E" if p else "NE" if aln else "U"
            conditioned = native_label if span else "U"
            admission[native_label,span] += 1
            certs = []
            for delta in DELTAS:
                first,second = label(l,delta),label(r,delta)
                certificate = first if first==second else "U"
                certs.append(certificate)
                counts["native_template_span_conditioned",delta,conditioned,certificate] += 1
                counts["same_tract_from_native_alignment",delta,projected,certificate] += 1
            w.writerow((name,native_label,int(span),conditioned,projected,*certs))
    assert len(identifiers)==73360 and set(alignments)<=identifiers
    result = {
        "N":len(identifiers),"native_template_IDs":len(positive),"native_primary_alignment_IDs":len(alignments),
        "native_sequence_identity":dict(sequence_identity),
        "paired_endpoint":"Two profiles: native CRISPR-A template event conditioned on a unique primary alignment spanning the target; and the exact ten-base tract projected from that native alignment using MP07 flank ownership. MP07 scores all admitted paths under its own affine objective.",
        "admission":[{"native_call":k,"full_local_span":s,"n":n} for (k,s),n in sorted(admission.items())],
        "counts":[{"profile":p,"delta":d,"crispra":c,"mp07":m,"n":n} for (p,d,c,m),n in sorted(counts.items())],
        "row_sha256":digest(out_rows),
    }
    a.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(len(positive),len(alignments),result["admission"])


if __name__=="__main__":
    main()
