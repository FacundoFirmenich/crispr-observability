"""Measure technical concordance of exact local strings in the two original mates."""

from collections import Counter
import gzip
import json

from freeze import ROOT, rc, sha
from run_external import fastq, mean_q30, primer_clip, agree, label_at


def literal_call(sequence, motifs):
    hits = [kind for kind, motif in motifs.items()
            if sequence.find(motif) != -1 and sequence.find(motif) == sequence.rfind(motif)]
    return hits[0] if len(hits) == 1 else "U"


def technical(a, b, spec, motifs):
    if not mean_q30(a[2]) or not mean_q30(b[2]):
        return "U", "U"
    head, tail = spec["reference_head_primer"], spec["reference_tail_primer"]
    candidates = []
    for first, second in ((a, b), (b, a)):
        left = primer_clip(first[1], head, rc(tail))
        right = primer_clip(second[1], tail, rc(head))
        if left is not None and right is not None:
            candidates.append((left, rc(right)))
    if len(candidates) != 1:
        return "U", "U"
    return tuple(literal_call(s, motifs) for s in candidates[0])


def motifs_for(spec):
    left = spec["wt"][spec["lo"]-10:spec["lo"]]
    right = spec["wt"][spec["hi_wt"]:spec["hi_wt"]+10]
    assert len(left) == len(right) == 10
    motifs = {key:left+spec[target]+right for key,target in (("edited","programmed_local"),("wt","unedited_local"))}
    assert spec["wt"].count(motifs["wt"]) == spec["edited"].count(motifs["edited"]) == 1
    assert motifs["edited"] not in spec["wt"] and motifs["wt"] not in spec["edited"]
    return motifs


def evaluate(sample, spec, folder):
    run = sample["run"]
    summary = json.loads((folder/f"{run}.json").read_bytes())
    read_path = folder/f"{run}_records.tsv.gz"
    family_path = folder/f"{run}_families.jsonl.gz"
    assert all(sha(p) == summary["output_sha256"][p.name] for p in (read_path,family_path))
    files = [ROOT/"fastq"/x["url"].split("/")[-1] for x in sample["files"]]
    motifs = motifs_for(spec)
    with gzip.open(family_path,"rt",encoding="ascii") as f:
        families = {row["key"]:row["references"] for row in map(json.loads,f)}
    outcomes = Counter()
    with gzip.open(read_path,"rt",encoding="ascii") as f:
        next(f)
        for a,b,line in zip(*(fastq(p) for p in files),f,strict=True):
            name,status,key = line.rstrip("\n").split("\t")
            assert name == a[0] == b[0]
            left,right = technical(a,b,spec,motifs)
            evidence = left if left==right else "discordant" if "U" not in (left,right) else "one_mate_only"
            outcomes["evidence",evidence] += 1
            if evidence in ("edited","wt"):
                for delta in (0,27,36,45,54):
                    label = agree(*[label_at(r["deficits"][0],delta) for r in families[key]]) if key else "U"
                    outcomes["paired",evidence,delta,label] += 1
                outcomes["statuses",evidence,status] += 1
    assert sum(n for k,n in outcomes.items() if k[0]=="evidence") == summary["N"]
    return {
        "run":run,"locus":sample["locus"],"arm":sample["arm"],"N":summary["N"],
        "motifs":motifs,
        "mate_evidence":[{"label":key[1],"n":n} for key,n in outcomes.items() if key[0]=="evidence"],
        "comparison":[{"mate_label":key[1],"delta":key[2],"certificate":key[3],"n":n}
                      for key,n in sorted(outcomes.items()) if key[0]=="paired"],
        "merge_status":[{"mate_label":key[1],"status":key[2],"n":n}
                      for key,n in sorted(outcomes.items()) if key[0]=="statuses"],
    }


if __name__ == "__main__":
    primary = ROOT/"EXTERNAL_FREEZE.json"
    contract = json.loads(primary.read_bytes())
    folder = ROOT/"external_results"
    specification = ROOT/"MATE_SUPPORT_CONTRACT.json"
    # Record this supplementary analysis before scanning the paired reads.
    with specification.open("x",encoding="utf-8") as f:
        json.dump({"scope":"Supplementary technical corroboration specified after observing merge coverage, before evaluating mate agreement. Same fragment, no orthogonal molecular labels.",
                   "rule":"Use two primer-authenticated raw mates, both mean Q30. Unique exact motif in each mate, built from event DNA plus ten unchanged bases on each side. Classify as corroborated only when both mates agree on the same unique motif.",
                   "comparison":"Join by original read ID and retain all primary outcomes; compare programmed-event certificates at 0,27,36,45,54.",
                   "primary_freeze_sha256":sha(primary),"code_sha256":sha(ROOT/"paired_mate_support.py")},
                  f,indent=2)
    result = [evaluate(sample,contract["targets"][sample["locus"]],folder) for sample in contract["samples"]]
    (ROOT/"MATE_SUPPORT.json").write_text(json.dumps({"contract_sha256":sha(specification),"runs":result},
                                                    indent=2)+"\n",encoding="utf-8")
