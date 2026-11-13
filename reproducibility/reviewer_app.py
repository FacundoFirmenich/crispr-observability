"""Local reviewer dashboard for the frozen public evidence."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from observability_report import analyse, demo_input, render, utility

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
language = st.sidebar.selectbox("Language / Idioma", ["English", "Castellano"])
def T(en, es):
    return es if language == "Castellano" else en
st.title("CRISPR Observability")
st.caption(T("Explore the public evidence, reproduce a section, and inspect its execution record.","Explore la evidencia pública, reproduzca una sección y examine su registro de ejecución."))
model = json.loads((GRAPH / "review_graph_model.json").read_bytes())
start, evidence, lab, replay, ledger, guides = st.tabs([
    T("Start here", "Empezar aquí"),
    T("Evidence and figures", "Evidencia y figuras"), T("3D research laboratory", "Laboratorio 3D"),
    T("Reproduce", "Reproducir"), T("Execution ledger", "Ledger de ejecución"),
    T("Researcher guides", "Guías de investigación")])

with start:
    st.subheader(T("One question. Four steps. A reproducible answer.", "Una pregunta. Cuatro pasos. Una respuesta reproducible."))
    st.markdown(T("**Question:** do the acquired CFTR records resolve a positive treated–mock contrast for the declared native-local sequence event?", "**Pregunta:** ¿los registros CFTR adquiridos resuelven un contraste positivo entre tratado y mock para el evento local declarado?"))
    st.markdown(T("**1. Define** the event and quality policy. **2. Account** for E, NE and U. **3. Interpret** the exact bounds. **4. Choose** the next check or measurement.", "**1. Defina** evento y calidad. **2. Conserve** E, NE y U. **3. Interprete** las cotas exactas. **4. Elija** la comprobación o medición siguiente."))
    data=demo_input()
    upload=st.file_uploader(T("Optional: load your observation contract (JSON)","Opcional: cargue su contrato de observación (JSON)"),type=['json'])
    try:
        if upload:data=json.loads(upload.getvalue())
        result=analyse(data)
        st.info(T("The output validates declared counts and arithmetic. It does not independently authenticate your labels.", "El resultado valida los recuentos declarados y la aritmética. No autentica independientemente sus etiquetas."))
        components.html(render(result,utility() if upload is None else None),height=900,scrolling=True)
        st.download_button(T("Download decision report (HTML)","Descargar informe de decisión (HTML)"),render(result,utility() if upload is None else None),'observability-report.html','text/html')
        st.download_button(T("Download exact results (JSON)","Descargar resultados exactos (JSON)"),json.dumps(result,indent=2),'observability-report.json','application/json')
        st.download_button(T("Download editable input contract","Descargar contrato de entrada editable"),json.dumps(data,indent=2),'observation-input.json','application/json')
    except (ValueError,KeyError,TypeError) as exc:
        st.error(T("Contract rejected: ","Contrato rechazado: ")+str(exc))

with evidence:
    st.subheader(T("Conditional finite-record bounds","Cotas condicionadas sobre registros finitos"))
    rows = pd.DataFrame(model["context_metrics"]["witness_observations"])
    preparation = st.selectbox(T("Preparation","Preparación"), list(dict.fromkeys(rows["preparation"])))
    subset = rows[rows["preparation"] == preparation]
    event = st.selectbox(T("Event","Evento"), list(dict.fromkeys(subset["event"])))
    selected = subset[subset["event"] == event]
    bars = alt.Chart(selected).mark_rule(strokeWidth=5).encode(
        x=alt.X("L_pp:Q", title=T("Percent of the raw records","Porcentaje de registros originales")),
        x2="H_pp:Q",
        y=alt.Y("condition:N", title=None),
        color=alt.Color("condition:N", legend=None, scale=alt.Scale(range=COLORS)),
        tooltip=["run:N", "condition:N", "N:Q", "E:Q", "NE:Q", "U:Q", "L_pp:Q", "H_pp:Q"],
    ).interactive()
    st.altair_chart(bars, width="stretch")
    st.dataframe(selected[["run", "condition", "N", "E", "NE", "U", "L_pp", "H_pp"]], hide_index=True)
    st.caption(T("The endpoints are bounds conditional on accepted labels. Unresolved records remain in the raw denominator.","Los extremos son cotas condicionadas a las etiquetas aceptadas. Los registros irresueltos permanecen en el denominador original."))
    st.subheader(T("Evidence relationships","Relaciones entre evidencias"))
    components.html((GRAPH / "crispr_supergraph.html").read_text(encoding="utf-8"), height=740, scrolling=True)

with replay:
    scope = st.selectbox(T("Choose the evidence to reproduce","Elija la evidencia que desea reproducir"), list(SCOPE),
        format_func=lambda x: T(x, {"all":"Todas las secciones", "editing":"CFTR y COL17A1", "germline":"Decisiones germinales"}[SCOPE[x]]))
    section = SCOPE[scope]
    modes = ["quick"] if section == "germline" else ["quick", "counts", "full"]
    mode = st.radio(T("Depth of replay","Profundidad de reproducción"), modes, horizontal=True)
    st.caption(T("The offline runner replays the distributed inputs and writes its scope to a new execution receipt.","El motor reproduce las entradas distribuidas y registra su alcance en un recibo nuevo."))
    if st.button(T("Run selected replay","Ejecutar reproducción seleccionada"), type="primary"):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        destination = OUTPUT / f"{stamp}-{section}-{mode}"
        command = [sys.executable, "-B", str(ROOT / "reproducibility" / "reproduce.py"),
                   "--section", section, "--mode", mode, "--output", str(destination)]
        with st.spinner(T("Replaying frozen evidence…","Reproduciendo evidencia congelada…"), show_time=True):
            result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                                    encoding="utf-8", errors="replace", timeout=1800)
        st.code(result.stdout + result.stderr)
        if result.returncode:
            st.error(T(f"Runner exit status: {result.returncode}. See the ledger.", f"Código de salida: {result.returncode}. Consulte el ledger."))
        else:
            st.info(T(f"Replay completed. Evidence and receipt: {destination}", f"Reproducción terminada. Evidencia y recibo: {destination}"))

with ledger:
    st.subheader(T("Local run history","Historial de ejecuciones"))
    receipts = sorted(OUTPUT.glob("*/REPRODUCTION.json"), reverse=True)
    if not receipts:
        st.info(T("Run a replay to create the first receipt.","Ejecute una reproducción para crear el primer recibo."))
    else:
        receipt = st.selectbox(T("Execution","Ejecución"), receipts, format_func=lambda p: p.parent.name)
        data = receipt.read_bytes()
        st.json(json.loads(data), expanded=1)
        st.download_button(T("Download execution receipt","Descargar recibo de ejecución"), data, "REPRODUCTION.json", "application/json")


with lab:
    st.subheader(T("Read Observability Lab", "Read Observability Lab"))
    st.caption(T("Rotate, inspect, change a parameter and export its computational consequences. Start with the leading case in Researcher guides.",
                 "Gire, inspeccione, cambie un parámetro y exporte sus consecuencias computacionales. Empiece por el caso guía en Guías de investigación."))
    html = (ROOT / "interactive/lab/CFTR_OBSERVABILITY_LAB.html").read_text(encoding="utf-8")
    components.html(html, height=1080, scrolling=True)
    st.download_button(T("Download standalone laboratory", "Descargar laboratorio autónomo"), html,
                       "Read_Observability_Lab.html", "text/html")
with guides:
    suffix = ".es.md" if language == "Castellano" else ".md"
    choice = st.selectbox(T("Guide", "Guía"), ["RESEARCHER_GUIDE", "LEADING_CASE", "RESEARCH_LINES", "COMMUNITY_ROADMAP"],
        format_func=lambda x: {"RESEARCHER_GUIDE": T("Researcher and community guide", "Guía de investigadores y comunidad"),
            "LEADING_CASE": T("Executable leading case", "Caso guía ejecutable"),
            "RESEARCH_LINES": T("Three research lines", "Tres líneas de investigación"),
            "COMMUNITY_ROADMAP": T("Contribute and develop", "Contribuir y desarrollar")}[x])
    text = (ROOT / "docs" / (choice + suffix)).read_text(encoding="utf-8")
    st.markdown(text)
    st.download_button(T("Download editable guide", "Descargar guía editable"), text, choice + suffix, "text/markdown")
