"""Join preserved CRISPResso2 projections with exact MP07 labels by original read ID."""

from collections import Counter
import csv
from fractions import Fraction
import gzip
import json
from pathlib import Path
import sys

from freeze import ROOT, DELTAS, sha

PREVIOUS = ROOT.parent / "CRISPR_P1_MP07_RECUPERACION"
sys.path.insert(0, str(PREVIOUS))
from recompute_panel import is_event
from recompute_universes import LEDGER_HASHES, label_at


def projection(column, coding, target, expected):
    return "U" if not coding[column] else "E" if is_event(target, expected, coding[column]) else "NE"


def certificate(record, index, delta):
    return label_at(record["deficits"][index], delta) if record["status"] == "evaluated" else "U"


def main():
    universes = PREVIOUS / "universes_MP07"
    contract = json.loads((universes / "full_receipt.json").read_bytes())
    assert len(contract["runs_references"]) == 24
    indexed = {(r["run"], r["reference"]): r for r in contract["runs_references"]}
    loci = PREVIOUS / "i2/extensions/correspondence_frontier/inputs"
    raw, comparisons = [], []
    seen = Counter()
    for cohort in ("cftr", "hek3"):
        caselist = json.loads((loci / f"{cohort}_REAL_ALIGNMENT_CASES.json").read_bytes())["cases"]
        if cohort == "cftr":
            refs = "mutant", "programmed"
            targets = {"native": caselist[0]["native_peptide"], "exact": caselist[0]["exact_programmed_coding"]}
        else:
            refs = "reference", "programmed_reference"
            targets = {k: caselist[0][k] for k in ("programmed_local","unedited_local")}
        ledger = PREVIOUS / "complete_record_ledgers" / f"{cohort.upper()}_COMPLETE_RECORD_LEDGER.csv.gz"
        assert sha(ledger) == LEDGER_HASHES[cohort]
        counts = Counter()
        with gzip.open(ledger, "rt", encoding="ascii", newline="") as stream:
            from itertools import groupby
            for accession, rows in groupby(csv.DictReader(stream), key=lambda x: x["run"]):
                left, right = [universes / f"{accession}_{ref}.jsonl.gz" for ref in refs]
                for path, ref in zip((left, right), refs, strict=True):
                    assert sha(path) == indexed[accession, ref]["records_sha256"]
                with gzip.open(left, "rt", encoding="ascii") as a, gzip.open(right, "rt", encoding="ascii") as b:
                    for row, first, second in zip(rows,a,b,strict=True):
                        x, y = json.loads(first), json.loads(second)
                        assert row["read_id"] == x["id"] == y["id"]
                        seen[cohort,accession] += 1
                        for idx,(name,expected) in enumerate(targets.items()):
                            m = projection("M_fw",row,name,expected)
                            p = projection("P_fw",row,name,expected)
                            both = m if m == p else "U"
                            old = {"M_fw":m, "P_fw":p, "M_and_P_consensus":both}
                            if cohort == "cftr":
                                old["PE_default"] = projection("PE_default",row,name,expected)
                            for delta in DELTAS:
                                cx, cy = certificate(x,idx,delta), certificate(y,idx,delta)
                                cb = cx if cx == cy else "U"
                                for operator,label in old.items():
                                    comparator = cx if operator == "M_fw" else cy if operator == "P_fw" else cb
                                    counts[accession,name,operator,delta,label,comparator] += 1
        for (run,name,operator,delta,old,new),n in sorted(counts.items()):
            raw.append({"cohort":cohort,"run":run,"target":name,"operator":operator,"delta":delta,
                        "baseline":old,"mp07":new,"n":n})
        compare = Counter()
        for row in raw:
            if row["cohort"] != cohort:
                continue
            k = tuple(row[x] for x in ("run","target","operator","delta"))
            compare[k,"paired_N"] += row["n"]
            a,b = row["baseline"],row["mp07"]
            if (a,b) == ("E","U"):
                compare[k,"baseline_event_audit_unknown"] += row["n"]
            if (a,b) == ("U","E"):
                compare[k,"baseline_unknown_audit_event"] += row["n"]
            if {a,b} == {"E","NE"}:
                compare[k,"opposite_determinate"] += row["n"]
            if a == "E":
                compare[k,"baseline_E"] += row["n"]
            if b == "E":
                compare[k,"audit_E"] += row["n"]
        for key, label in sorted(compare):
            if label != "paired_N":
                continue
            run,name,operator,delta = key
            assert compare[key,label] == seen[cohort,run]
            comparisons.append({"cohort":cohort,"run":run,"target":name,"operator":operator,"delta":delta,
                                **{k:compare[key,k] for k in ("paired_N","baseline_E","audit_E",
                                    "baseline_event_audit_unknown","baseline_unknown_audit_event",
                                    "opposite_determinate")}})
    assert sum(seen.values()) == 395936 and len(seen) == 12
    for name,rows in (("PAIRED_C2_CONFUSION.csv",raw),("PAIRED_C2_SUMMARY.csv",comparisons)):
        with (ROOT / name).open("w",encoding="utf-8",newline="") as f:
            w = csv.DictWriter(f,fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    (ROOT / "PAIRED_C2_RECEIPT.json").write_text(json.dumps({
        "scope":"Actual retained CRISPResso2 local projections, joined on authentic read IDs. Exact events and quality are harmonized. No new CRISPResso2 run; the dedicated PE_default is distinct from the single-reference modes.",
        "evidence":"Comparisons by locus, reference profile, target and deficit. U is abstention; opposite labels are methodological disagreements and are not adjudicated against molecular truth.",
        "total_raw_pairs":sum(seen.values()),"matched_runs":len(seen),
        "mp07_receipt_sha256":sha(universes/"full_receipt.json"),
        "outputs":{n:sha(ROOT/n) for n in ("PAIRED_C2_CONFUSION.csv","PAIRED_C2_SUMMARY.csv")}},
        indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    main()
