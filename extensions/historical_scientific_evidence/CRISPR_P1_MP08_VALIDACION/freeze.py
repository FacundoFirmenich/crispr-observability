"""Freeze selection and sequence identity from publisher/ENA metadata before reading FASTQ."""

import csv
import hashlib
import json
from pathlib import Path
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
LOCI = ("PCSK9", "CTLA4", "HBEGF")
COMP = str.maketrans("ACGTN", "TGCAN")
DELTAS = [0, 9, 18, 27, 36, 45, 54]


def rc(s):
    return s.translate(COMP)[::-1]


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def workbook(index):
    file = INPUT / f"41467_2025_58653_MOESM{index}_ESM.xlsx"
    with ZipFile(file) as archive:
        assert archive.testzip() is None
        strings = ["".join(e.itertext()) for e in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        for row in sheet.findall(".//{*}row"):
            cells = {}
            for cell in row.findall("{*}c"):
                value = cell.findtext("{*}v")
                if value is not None:
                    name = re.match("[A-Z]+", cell.get("r"))[0]
                    cells[name] = strings[int(value)] if cell.get("t") == "s" else value
            yield int(row.get("r")), cells


def construct_target(wt, edited, expected):
    assert len(edited) - len(wt) == expected > 0
    prefix = next(i for i, (a, b) in enumerate(zip(wt, edited)) if a != b)
    suffix = next(i for i, (a, b) in enumerate(zip(wt[::-1], edited[::-1])) if a != b)
    leftmost, rightmost = sorted((prefix, len(wt) - suffix))
    lo, hi = leftmost - 8, rightmost + 8
    assert wt[:lo] == edited[:lo] and wt[hi:] == edited[hi+expected:]
    assert all(0 < lo < end < len(ref) for end, ref in ((hi, wt), (hi+expected, edited)))
    return lo, hi, hi + expected


def freeze():
    primers = {r["A"]: (n, r["B"], r["C"]) for n, r in workbook(4) if n > 2 and r.get("A") in LOCI}
    guides = {r["A"].strip(): (n, r["D"].replace("U", "T"), r["F"].replace("U", "T"), r["C"])
              for n, r in workbook(3) if r.get("B") == "pegRNA"}
    targets = {}
    for locus in LOCI:
        pn, forward, reverse = primers[locus]
        gn, guide, extension, edit = guides[locus]
        matches = []
        for row, data in workbook(5):
            if row <= 2 or not data.get("C"):
                continue
            wt = data["C"]
            primer_match = any(wt.startswith(a) and wt.endswith(rc(b)) for a, b in
                               ((forward, reverse), (reverse, forward)))
            if primer_match and data.get("B") == guide and data.get("D") == extension:
                matches.append((row, data))
        (row, data), = matches
        wt, edited = data["C"], data["L"]
        expected = int(re.match(r"(\d+)-bp insertion", edit)[1])
        lo, hi_wt, hi_edited = construct_target(wt, edited, expected)
        head, tail = (forward, reverse) if wt.startswith(forward) else (reverse, forward)
        targets[locus] = {
            "primer_row": pn, "guide_row": gn, "parameter_row": row,
            "parameter_printed_label": data["A"], "spacer": guide, "extension": extension,
            "reference_head_primer": head, "reference_tail_primer": tail, "wt": wt, "edited": edited,
            "lo": lo, "hi_wt": hi_wt, "hi_edited": hi_edited,
            "programmed_local": edited[lo:hi_edited], "unedited_local": wt[lo:hi_wt],
            "published_crispresso2_parameters": data,
        }
    with (INPUT / "PRJNA1226148_ENA.tsv").open(encoding="utf-8", newline="") as source:
        by_alias = {}
        for row in csv.DictReader(source, delimiter="\t"):
            by_alias.setdefault(row["experiment_alias"], []).append(row)
    samples = []
    locus = group = ""
    for n, row in workbook(7):
        if n < 4:
            continue
        if row.get("A"):
            locus = row["A"].split()[0]
        if row.get("D"):
            group = row["D"]
        if locus not in LOCI or not re.fullmatch(r"K562 x[123]", group):
            continue
        if row.get("B") == "PEn" and row.get("C") == "pegRNA":
            arm = "PE_unmodified"
        elif row.get("B") == row.get("C") == "-":
            arm = "untreated"
        else:
            continue
        data, = by_alias[row["F"]]
        assert data["library_layout"] == "PAIRED" and locus in data["run_alias"] and "K562" in data["run_alias"]
        files = [{"url": "https://" + u, "md5": h, "bytes": int(b)} for u,h,b in
                 zip(*(data[k].split(";") for k in ("fastq_ftp","fastq_md5","fastq_bytes")), strict=True)]
        assert len(files) == 2
        samples.append({"locus": locus, "replicate": int(group[-1]), "arm": arm,
                        "metadata_row": n, "library_alias": row["F"], "run": data["run_accession"],
                        "ENA_read_count": int(data["read_count"]), "files": files})
    expected = {(loc, rep, arm) for loc in LOCI for rep in (1,2,3)
                for arm in (("untreated", "PE_unmodified") if loc == "PCSK9" else ("PE_unmodified",))}
    assert len(samples) == len(expected) == 12
    assert {(s["locus"],s["replicate"],s["arm"]) for s in samples} == expected
    return {
        "study": "Antoniou et al., 10.1038/s41467-025-58653-1; ENA PRJNA1226148",
        "selection": "All K562 replicates from first three human-locus blocks of Supplementary Data 5, PEn with unmodified pegRNA and paired untreated controls when deposited. No editing frequencies used.",
        "targets": targets, "samples": samples, "score_deltas": DELTAS,
        "observation": "One original read pair per unit. Both raw mate means >=30. Require complete opposing locus-specific 5-prime primer in each mate with offset <=8, trim its upstream prefix, trim readthrough following the first complete downstream opposite primer if it occurs after position 24. Unique exact overlap >=24 bases; both mates must agree throughout their overlap. Reject any N in merged sequence. No fill-in. Every rejection is U in raw N.",
        "analysis": "Both wild-type and programmed amplicons in primer-authenticated orientation, frozen MP07 affine score and flank ownership. Core includes insertion equivalence interval with eight unchanged nucleotides each side. Exact local DNA, both programmed and unedited events. Agreement on determinate label between references required. Report E/NE/U by accession; untreated contrasts only when matched, using [E/N,(E+U)/N]. Report 0,9,18,27,36,45,54; no pooling across loci.",
        "decision": "A change from positive to nonidentified lower bound on the local programmed event is a change in the reportable sign. Count direct label contradictions separately from projection failure. Joint reference-independent molecular accuracy requires linked independent truth.",
        "freeze_sha256": {p.name: sha(p) for p in sorted(INPUT.iterdir()) if p.is_file()},
    }


if __name__ == "__main__":
    result = freeze()
    out = ROOT / "EXTERNAL_FREEZE.json"
    with out.open("x", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(sha(out))
    for locus, spec in result["targets"].items():
        print(locus, spec["parameter_printed_label"], spec["parameter_row"],
              spec["programmed_local"], spec["unedited_local"])
    print([(s["run"], s["locus"], s["arm"]) for s in result["samples"]])
