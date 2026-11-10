"""Locate and independently rescore the CFTR determinate discordances at delta zero."""

from collections import defaultdict
import csv
import gzip
from itertools import groupby
import json
import sys

from freeze import ROOT, sha
from run_external import PREVIOUS, agree, label_at
from recompute_panel import is_event, solve
from recompute_universes import allele_orientations, canonical

frontier = PREVIOUS/"i2/extensions/correspondence_frontier"
sys.path.insert(0,str(frontier))
import optimal_lattice


def main():
    cases = json.loads((frontier/"inputs/cftr_REAL_ALIGNMENT_CASES.json").read_bytes())["cases"]
    references = {c["reference"]:c for c in cases}
    peptide = cases[0]["native_peptide"]
    universe = PREVIOUS/"universes_MP07"
    source = PREVIOUS/"complete_record_ledgers/CFTR_COMPLETE_RECORD_LEDGER.csv.gz"
    output = []
    with gzip.open(source,"rt",encoding="ascii",newline="") as f:
        for run, rows in groupby(csv.DictReader(f),key=lambda row:row["run"]):
            files = [universe/f"{run}_{ref}.jsonl.gz" for ref in references]
            with gzip.open(files[0],"rt",encoding="ascii") as a,gzip.open(files[1],"rt",encoding="ascii") as b:
                differences = []
                for row,la,lb in zip(rows,a,b,strict=True):
                    x,y = json.loads(la),json.loads(lb)
                    assert row["read_id"]==x["id"]==y["id"]
                    if not row["PE_default"] or not all(z["status"]=="evaluated" for z in (x,y)):
                        continue
                    prior = "E" if is_event("native",peptide,row["PE_default"]) else "NE"
                    now = agree(*[label_at(z["deficits"][0],0) for z in (x,y)])
                    if {prior,now} == {"E","NE"}:
                        differences.append((row,prior,now))
            if not differences:
                continue
            tables = {}
            for ref,c in references.items():
                source = PREVIOUS/"i2/raw_replay/inputs/alleles"/ref/f"{run}.zip"
                tables[ref],_ = allele_orientations(source,c["reference_sequence"])
            for row,prior,now in differences:
                entry = {"run":run,"read_id":row["read_id"],"PE_default_label":prior,
                         "PE_default_projection":row["PE_default"],"MP07_delta0":now,"references":{}}
                for ref,template in references.items():
                    oriented, = tables[ref][canonical(row["raw_sequence"])]
                    optimum, targets = solve({**template,"observed_sequence":oriented},optimal_lattice,cutoff=-1)
                    witness = targets["native"]
                    assert {k for k,w in witness.items() if w["score"]==optimum} == {int(now=="E")}
                    entry["references"][ref] = {"optimum":optimum,"witnesses":witness}
                output.append(entry)
    assert len(output) == 9
    (ROOT/"C2_DISCORDANT_WITNESSES.json").write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    print([(e["read_id"],e["PE_default_label"],e["MP07_delta0"]) for e in output])


if __name__ == "__main__":
    main()
