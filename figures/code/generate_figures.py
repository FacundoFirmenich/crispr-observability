"""Re-render the four P1 figures from frozen MP02 tables and stated geometry.

No image tracing, inferred measurements or raster embedding. Original MP02
rasters remain in originals/. SVG and PDF outputs contain native vectors.
"""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import csv
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "sources"
INPUTS = ROOT / "inputs"
FIGURES = ROOT / "figures"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(name):
    with (SOURCES / name).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def save_csv(name, rows):
    with (INPUTS / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def verify_sources():
    manifest = json.loads((ROOT / "SOURCE_MANIFEST.json").read_text())
    for entry in manifest["files"]:
        p = ROOT / entry["local_path"]
        assert p.stat().st_size == entry["size_bytes"], str(p)
        assert sha(p) == entry["sha256"], str(p)


def prepare_inputs():
    # Figure 1 uses its explicit public data contract, independently of private text.
    contract = json.loads((SOURCES / "anchor_geometry_contract.json").read_text())
    geometry = contract["rows"]
    assert [(r["start"], r["end"], r["read_end"]) for r in geometry] == [(210, 222, 224), (213, 225, 224), (206, 218, 224)]
    save_csv("fig1_coordinate_geometry.csv", geometry)

    intervals = read_csv("Table_1_Native_Contrast_Intervals.csv")
    assert len(intervals) == 21
    for row in intervals:
        assert abs(float(row["L_pp"]) - float(Fraction(row["L_exact"]) * 100)) < 1e-12
        assert abs(float(row["H_pp"]) - float(Fraction(row["H_exact"]) * 100)) < 1e-12
        assert -12 < float(row["L_pp"]) <= float(row["H_pp"]) < 45
    save_csv("fig2_native_contrast_intervals.csv", intervals)

    contrasts = read_csv("Table_S3_M1_All_Contrasts.csv")
    selected = {}
    for row in contrasts:
        if (row["prep"] == "CFB-26-16" and row["control_condition"] == "mock"
                and row["event"] == "native_local_peptide_all_observed_sequences"
                and row["policy"] == "common_raw_meanQ30"
                and row["observer"] in ("broad_edge", "broad_current", "target")):
            key = (row["observer"], int(row["trim"]))
            assert key not in selected
            selected[key] = row
    assert len(selected) == 6
    labels = {"broad_edge": "Edge anchor", "broad_current": "Inner anchor", "target": "Target motifs"}
    truncation = []
    for observer in ("broad_edge", "broad_current", "target"):
        for trim in (0, 8):
            r = selected[(observer, trim)]
            truncation.append({"label": f"{labels[observer]} trim {trim}",
                               "observer": observer, "trim": trim, "original_veto": False,
                               "L_exact": r["L_exact"], "H_exact": r["H_exact"],
                               "L_pp": r["L_pp"], "H_pp": r["H_pp"]})
    veto = [r for r in read_csv("Table_S4_Original_Ambiguity_Veto.csv")
            if r["policy"] == "common_raw_meanQ30"
            and r["event"] == "native_local_peptide_all_observed_sequences"
            and r["observer"] == "target" and r["trim"] == "8"
            and r["run_accession"] in ("SRR29623272", "SRR29623273")]
    assert len(veto) == 2
    t = next(r for r in veto if r["run_accession"] == "SRR29623272")
    c = next(r for r in veto if r["run_accession"] == "SRR29623273")
    # Exact two-arm identification hull; no sampling interval is computed.
    low = Fraction(int(t["E"]), int(t["N"])) - Fraction(int(c["E"]) + int(c["U"]), int(c["N"]))
    high = Fraction(int(t["E"]) + int(t["U"]), int(t["N"])) - Fraction(int(c["E"]), int(c["N"]))
    truncation.append({"label": "Target trim 8 with original veto", "observer": "target", "trim": 8,
                       "original_veto": True, "L_exact": str(low), "H_exact": str(high),
                       "L_pp": float(low * 100), "H_pp": float(high * 100)})
    save_csv("fig3_truncation_intervals.csv", truncation)

    changes = read_csv("Table_S9_Reference_Changes.csv")
    assert len(changes) == 137
    assert all(r["mutant_projection_native"] == "N" and r["programmed_projection_native"] == "E" for r in changes)
    counts = Counter(r["run"] for r in changes)
    prep_runs = [("CFB-26-16", "SRR29623272"), ("CFB-26-02", "SRR29623271"), ("CFB-06-04", "SRR29623270")]
    reference = [{"preparation": prep, "run_accession": run, "NE_to_E": counts[run]} for prep, run in prep_runs]
    assert [r["NE_to_E"] for r in reference] == [34, 58, 45]
    save_csv("fig4_reference_changes.csv", reference)

    domains = json.loads((SOURCES / "WITNESS_BOUNDS.json").read_text())["domains"]
    assert len(domains) == 9
    total = sum(r["N_raw"] for r in domains)
    catkeys = [("Different coding", "different_coding"), ("Quality failure", "common_quality_fail"),
               ("Missing correspondence", "missing_correspondence_witness"),
               ("Coding witness", "identical_observed_coding_witness")]
    categories = [{"category": name, "count": sum(r["counts"].get(key, 0) for r in domains)} for name, key in catkeys]
    assert total == sum(r["count"] for r in categories) == 177006
    assert [r["count"] for r in categories] == [83, 1996, 4202, 170725]
    save_csv("fig4_witness_partition.csv", categories)
    return geometry, intervals, truncation, reference, categories


def save(fig, stem):
    for suffix in ("pdf", "svg", "png"):
        kwargs = {"dpi": 300} if suffix == "png" else {}
        # Fixed metadata avoids wall-clock timestamps in vector deliverables.
        metadata = {"Creator": "P1 MP03 figure regeneration", "Date": "2026-10-28"} if suffix == "svg" else None
        if suffix == "pdf":
            metadata = {"Creator": "P1 MP03 figure regeneration", "CreationDate": None, "ModDate": None}
        fig.savefig(FIGURES / f"{stem}.{suffix}", bbox_inches="tight", metadata=metadata, **kwargs)
    plt.close(fig)


def figure1(rows):
    fig, ax = plt.subplots(figsize=(7.3, 3.0))
    for y, row in zip((2, 1, 0), rows):
        ax.hlines(y, 202, row["read_end"], color="#999999", lw=2)
        ax.barh(y, row["end"] - row["start"], left=row["start"], height=.30,
                color="#AA3377" if row["end"] > row["read_end"] else "#4477AA", alpha=.85)
        ax.text(row["start"], y + .22, f'[{row["start"]},{row["end"]})', fontsize=10)
    ax.axvline(224, color="black", ls="--", lw=1.5)
    ax.text(224, 2.53, "224 nt read end", ha="right", va="top", fontsize=10)
    ax.set(xlim=(201, 228), ylim=(-.6, 3), xlabel="Zero based read coordinate", yticks=(2, 1, 0),
           yticklabels=[r["label"] for r in rows], xticks=(205, 210, 215, 220, 225))
    for s in ("left", "top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    save(fig, "fig1")


def figure2(rows):
    names = {"broad_edge": "Edge anchor", "broad_current": "Inner anchor", "target": "Target motifs",
             "crispresso_projected": "Mutant projection", "programmed_reference": "Programmed projection",
             "reference_agreement": "Reference label agreement", "witness": "Coding witness"}
    colors = ["#7A3E00", "#4477AA", "#228833", "#AA3377", "#66AABB", "#777777", "#222222"]
    fig, axes = plt.subplots(1, 3, figsize=(7.3, 4.3), sharey=True)
    for ax, prep in zip(axes, ("CFB-26-16", "CFB-26-02", "CFB-06-04")):
        selected = [r for r in rows if r["preparation"] == prep]
        assert len(selected) == 7
        for y, row, color in zip(range(6, -1, -1), selected, colors):
            lo, hi = float(row["L_pp"]), float(row["H_pp"])
            ax.plot((lo, hi), (y, y), color=color, lw=2.8)
            ax.plot((lo, hi), (y, y), "|", color=color, ms=10)
        ax.axvline(0, color="#555555", ls="--", lw=.8)
        ax.set(title=prep, xlim=(-12, 45), xticks=(-10, 0, 10, 20, 30, 40), xlabel="PE minus mock (pp)")
        ax.grid(axis="x", alpha=.18)
    axes[0].set_yticks(range(7), [names[r["operator"]] for r in reversed(rows[:7])])
    save(fig, "fig2")


def figure3(rows):
    fig, axes = plt.subplots(2, 1, figsize=(7.3, 5.8), gridspec_kw={"height_ratios": (1.8, 1)})
    for ax, selected in ((axes[0], rows), (axes[1], rows[4:])):
        for y, row in enumerate(selected):
            lo, hi = float(row["L_pp"]), float(row["H_pp"])
            color = "#228833" if row["observer"] == "target" else "#4477AA"
            ax.plot((lo, hi), (y, y), color=color, lw=2.6)
            ax.plot((lo, hi), (y, y), "|", color="#222222", ms=10)
        ax.set_yticks(range(len(selected)), [r["label"] for r in selected])
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=.18)
    axes[0].set(xlim=(-110, 110), xlabel="PE minus mock native interval (pp)")
    axes[0].axvline(0, color="#777777", ls="--", lw=1.3)
    axes[1].set(xlim=(12.8, 23.1), xlabel="Target interval detail (pp)")
    fig.subplots_adjust(hspace=.5)
    save(fig, "fig3")


def figure4(reference, categories):
    fig, axes = plt.subplots(1, 2, figsize=(7.7, 3.6), gridspec_kw={"width_ratios": (1, 1.3)})
    counts = [r["NE_to_E"] for r in reference]
    axes[0].bar(range(3), counts, color="#AA3377")
    for i, n in enumerate(counts):
        axes[0].text(i, n + 1, str(n), ha="center", fontsize=10)
    axes[0].set(title="Reference sensitivity", ylabel="NE to E read changes", ylim=(0, 66),
                xticks=range(3), xticklabels=[r["preparation"] for r in reference])
    plt.setp(axes[0].get_xticklabels(), rotation=25, ha="right")
    counts = [r["count"] for r in categories]
    axes[1].barh(range(4), counts, color=("#AA3377", "#4477AA", "#777777", "#228833"))
    axes[1].set_xscale("log")
    for i, n in enumerate(counts):
        axes[1].text(n * 1.10, i, str(n), va="center", fontsize=9)
    axes[1].set(title="Raw universe 177006", yticks=range(4),
                yticklabels=[r["category"] for r in categories], xlim=(40, 700000), xlabel="Records (log scale)")
    axes[1].invert_yaxis()
    fig.subplots_adjust(wspace=.95)
    save(fig, "fig4")


def main():
    for folder in (INPUTS, FIGURES, ROOT / "verification"):
        folder.mkdir(exist_ok=True, parents=True)
    verify_sources()
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none",
                         "svg.hashsalt": "crispr-observability-mp03d"})
    geometry, intervals, truncation, reference, categories = prepare_inputs()
    figure1(geometry)
    figure2(intervals)
    figure3(truncation)
    figure4(reference, categories)
    result = {"status": "PASS", "source_version": "MP02", "rendering_version": "MP03",
              "scope": "New vector rendering of frozen source values; originals preserved byte-for-byte.",
              "source_hashes_verified": True, "figure2_xlim": [-12, 45], "figure2_intervals": 21,
              "figure3_exact_fraction_hulls": True, "figure4_partition_total": 177006,
              "matplotlib_version": matplotlib.__version__,
              "files": [{"path": str(p.relative_to(ROOT)), "sha256": sha(p), "size_bytes": p.stat().st_size}
                        for folder in (INPUTS, FIGURES) for p in sorted(folder.iterdir()) if p.is_file()]}
    (ROOT / "verification/FIGURE_REGENERATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "files"}, indent=2))


if __name__ == "__main__":
    main()
