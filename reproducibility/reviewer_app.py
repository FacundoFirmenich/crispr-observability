"""Local reviewer dashboard for the frozen public evidence."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

import altair as alt
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT / "interactive" / "supergraph"
OUTPUT = Path(os.environ.get("CRISPR_REVIEW_OUTPUT", Path.home() / "CRISPR-Observability"))
SCOPE = {
    "All available replay sections": "all",
    "CFTR and COL17A1": "editing",
    "Germline genotype decisions": "germline",
}
COLORS = ["#527ed7", "#9664cd", "#e38fae"]
st.set_page_config(page_title="CRISPR Observability", page_icon="🧬", layout="wide")
st.title("CRISPR Observability")
st.caption("Explore the public evidence, reproduce a section, and inspect its execution record.")
model = json.loads((GRAPH / "review_graph_model.json").read_bytes())
evidence, replay, ledger = st.tabs(["Evidence and figures", "Reproduce", "Execution ledger"])

with evidence:
    st.subheader("Conditional finite-record bounds")
    rows = pd.DataFrame(model["context_metrics"]["witness_observations"])
    preparation = st.selectbox("Preparation", list(dict.fromkeys(rows["preparation"])))
    subset = rows[rows["preparation"] == preparation]
    event = st.selectbox("Event", list(dict.fromkeys(subset["event"])))
    selected = subset[subset["event"] == event]
    bars = alt.Chart(selected).mark_rule(strokeWidth=5).encode(
        x=alt.X("L_pp:Q", title="Percent of the raw records"),
        x2="H_pp:Q",
        y=alt.Y("condition:N", title=None),
        color=alt.Color("condition:N", legend=None, scale=alt.Scale(range=COLORS)),
        tooltip=["run:N", "condition:N", "N:Q", "E:Q", "NE:Q", "U:Q", "L_pp:Q", "H_pp:Q"],
    ).interactive()
    st.altair_chart(bars, width="stretch")
    st.dataframe(selected[["run", "condition", "N", "E", "NE", "U", "L_pp", "H_pp"]], hide_index=True)
    st.caption("The endpoints are bounds conditional on accepted labels. Unresolved records remain in the raw denominator.")
    st.subheader("Evidence relationships")
    components.html((GRAPH / "crispr_supergraph.html").read_text(encoding="utf-8"), height=740, scrolling=True)

with replay:
    scope = st.selectbox("Choose the evidence to reproduce", list(SCOPE))
    section = SCOPE[scope]
    modes = ["quick"] if section == "germline" else ["quick", "counts", "full"]
    mode = st.radio("Depth of replay", modes, horizontal=True)
    st.caption("The offline runner replays the distributed inputs and writes its scope to a new execution receipt.")
    if st.button("Run selected replay", type="primary"):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        destination = OUTPUT / f"{stamp}-{section}-{mode}"
        command = [sys.executable, "-B", str(ROOT / "reproducibility" / "reproduce.py"),
                   "--section", section, "--mode", mode, "--output", str(destination)]
        with st.spinner("Replaying frozen evidence…", show_time=True):
            result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                                    encoding="utf-8", errors="replace", timeout=1800)
        st.code(result.stdout + result.stderr)
        if result.returncode:
            st.error(f"The runner exited with status {result.returncode}. See the execution ledger for details.")
        else:
            st.success(f"Replay completed. Evidence and receipt saved in {destination}.")

with ledger:
    st.subheader("Local run history")
    receipts = sorted(OUTPUT.glob("*/REPRODUCTION.json"), reverse=True)
    if not receipts:
        st.info("Run a replay to create the first receipt.")
    else:
        receipt = st.selectbox("Execution", receipts, format_func=lambda p: p.parent.name)
        data = receipt.read_bytes()
        st.json(json.loads(data), expanded=1)
        st.download_button("Download execution receipt", data, "REPRODUCTION.json", "application/json")
