# Published reference and adaptation

This component adapts the spatial SVG interface of **Interactive Epistemic Map of
Claims, Evidence, and Falsification**, accompanying *Beyond Redfield: Thermodynamic
Bounds and Non-Perturbative Quantum Dynamics in Tubulin Networks*, by Facundo
Firmenich, Pau Firmenich and León Firmenich (bioRxiv DOI
10.64898/2026.05.10.724047).

The exact reference was recovered from the published KwanTube v3.5.1.1 archive:

- Zenodo record: https://zenodo.org/records/20102651
- Archive: `KwanTube-v.3.5.1.1.zip`
- Archived interface: `FacundoFirmenich-KwanTube-620bdeb/outputs_data/interactive/epistemic_graph.html`
- Interface SHA-256: `fd59a3f2409bf6b752eee81a203ffc312880f8d664ce95553df4687b40107c42`
- Repository: https://github.com/FacundoFirmenich/KwanTube

The archived interface uses the title above, a dark full-viewport spatial canvas,
collapsible filter and inspector panels, SVG nodes and typed edges, lenses,
search, node/edge inspection, pan, zoom and fit controls. This CRISPR component
preserves that architecture and its visual language. It replaces the tubulin
content with a verified CRISPR evidence model, adds explicit AND-gate navigation,
preparation/event witness readouts, keyboard and touch navigation, exports and
embedded-data verification. It does not reproduce the original scientific data.

The original KwanTube project is distributed under the GNU General Public License
version 3. The complete upstream license text is retained in `LICENSE`.
The derived interface, template and build/check code in this component are
provided under **GPL-3.0-only**. They are separate from the original-code MIT
allowlist elsewhere in this repository. The model JSON, scientific tables,
observations, figure data and third-party assets retain their separate attribution
and applicable terms; this notice does not relicense them.

Changes made for this CRISPR adaptation are dated 28 October 2026. Corresponding
editable source is supplied in `crispr_supergraph.template.html` and `code/`.
No manuscript or formal supplement text is embedded. Source links refer to
public numerical artifacts only. The older four-view CRISPR explorer remains a
separate, unchanged output.
