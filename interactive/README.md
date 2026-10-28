# Offline evidence explorer

Download `epistemic_graph.html` and open it in a browser. The HTML contains its
own data, scripts and styles; it does not contact a CDN or require a server.
GitHub's source-code viewer displays the source rather than running the tool.

The four views provide an evidence dependency map, finite-record readout hulls,
the historical alignment-parameter sensitivity catalog, and calculation/source
contracts. Search, keyboard node selection, premise toggles, mobile text
navigation, view links and CSV/JSON/SVG exports support inspection.

An AND hyperedge declares which premises are jointly required for a displayed
support path. Toggling a premise changes that display; it neither removes source
observations nor establishes or refutes a scientific claim. The unknown-completion
slider selects an integer number of unknown records hypothetically assigned to
the event. It cannot add biological evidence.

## Rebuild and verify

Run from the repository root, using Python 3.11 or later:

```sh
python3 -B interactive/code/build_model.py
python3 -B interactive/code/check_contracts.py
python3 -B interactive/code/build_explorer.py
python3 -B interactive/code/check_explorer.py
```

`sources/` is a minimal frozen copy of the 19 canonical numerical tables and three
additional numerical artifacts. Verification requires every copy's bytes to match
the canonical public path recorded in the model. It also rebuilds the model and
HTML in a temporary directory and compares their bytes with this release. See
`DATA_MODEL.md` for catalog units, event definitions and finite attainability.

`BUILD_RECEIPT.json` distinguishes the SHA-256 of the input JSON file from that of
the exact JSON serialization embedded in HTML. Both represent the same parsed
object. The browser verifies the embedded serialization when its Web Crypto API
is available; the Python checker verifies all three identities without a browser.

The original interface and build code are independently implemented. KwanTube's
separation of public numerical artifacts from formal text and its inspectable
evidence map informed the architecture; no KwanTube source code is copied here.
The scoped MIT grant covers only the enumerated original code. Source data and
scientific content retain their attribution and applicable terms.

The formal article and supplement are absent. Conditional identification hulls
assume correct accepted labels; they are not confidence intervals or comparative
molecular-accuracy estimates. Reproduction here remains computational and
retrospective.
