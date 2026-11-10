"""Adjudicate the observable local motif in the three original untreated libraries."""

from collections import Counter, defaultdict
import gzip
import json

from freeze import ROOT, sha
from paired_mate_support import motifs_for, technical
from run_external import agree, fastq, label_at


def main():
    contract = json.loads((ROOT/"EXTERNAL_FREEZE.json").read_bytes())
    spec = contract["targets"]["PCSK9"]
    motifs = motifs_for(spec)
    folder = ROOT/"external_results"
    results = []
    for sample in contract["samples"]:
        if (sample["locus"],sample["arm"]) != ("PCSK9","untreated"):
            continue
        run = sample["run"]
        summary = json.loads((folder/f"{run}.json").read_bytes())
        reads = folder/f"{run}_records.tsv.gz"
        families = folder/f"{run}_families.jsonl.gz"
        assert all(sha(path) == summary["output_sha256"][path.name] for path in (reads,families))
        files = [ROOT/"fastq"/p["url"].split("/")[-1] for p in sample["files"]]
        with gzip.open(families,"rt",encoding="ascii") as f:
            positives = {row["key"]:row for row in map(json.loads,f)
                         if agree(*[label_at(ref["deficits"][0],0) for ref in row["references"]])=="E"}
        counts = Counter()
        examples = defaultdict(list)
        with gzip.open(reads,"rt",encoding="ascii") as f:
            next(f)
            for a,b,line in zip(*(fastq(p) for p in files),f,strict=True):
                name,status,key = line.rstrip("\n").split("\t")
                assert name == a[0] == b[0]
                if key not in positives:
                    continue
                assert status=="evaluated"
                left,right = technical(a,b,spec,motifs)
                counts[left,right] += 1
                if not examples[left,right]:
                    examples[left,right].append({"id":name,"sequence_a":a[1],"sequence_b":b[1],
                                                "merged_sequence":positives[key]["sequence"]})
        assert sum(counts.values()) == sum(row["N"] for row in positives.values())
        assert sum(counts.values()) == next(row["E"] for row in summary["counts"]
                                          if (row["target"],row["delta"])==("programmed_local",0))
        results.append({"run":run,"N":summary["N"],"mp07_E_at_optimum":sum(counts.values()),
                        "literal_mate_calls":[{"head_mate":a,"tail_mate":b,"N":n,"examples":examples[a,b]}
                                              for (a,b),n in sorted(counts.items())]})
    assert len(results)==3
    (ROOT/"PCSK9_CONTROL_OBSERVATIONS.json").write_text(json.dumps({
        "scope":"Follow-up after observing the untreated-control positives. The raw mates belong to the same sequenced fragment.",
        "motifs":motifs,"runs":results},indent=2)+"\n",encoding="utf-8")
    print([(r["run"],r["mp07_E_at_optimum"],
            [(c["head_mate"],c["tail_mate"],c["N"]) for c in r["literal_mate_calls"]]) for r in results])


if __name__=="__main__":
    main()
