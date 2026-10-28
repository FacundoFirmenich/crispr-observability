# CRISPR scientific supergraph

Download `crispr_supergraph.html` and open it in a browser. The file contains its
styles, script and frozen graph model and works offline. GitHub's code viewer
displays source rather than running the interface.

This is the spatial-map adaptation of the published *Beyond Redfield* companion:
a dark full-viewport graph, collapsible filter and inspector panels, scientific
lenses, status views, node and edge inspection, pan/zoom/fit, and SVG/JSON exports.
It is a second output; the original four-view numerical explorer is retained.
`Focus graph` collapses both panels. `Filters` and `Inspector` restore them.

On a touchscreen, drag to pan and pinch to zoom. Select a node to open its
inspector. The filter panel also offers a native node selector. With keyboard
focus on the canvas, arrow keys pan, `+` and `-` zoom, and `0` fits the graph;
Enter or Space activates a focused node or link. Escape closes a mobile panel.

The map has 57 nodes and 59 typed links: 21 core nodes, eight AND gates, 22 source
artifacts and six declared revision tests. The source-file nodes are initially
hidden for legibility; their links remain available in the inspector. Revision
tests are conditions requiring evidence, not newly performed experiments.
Node positions and visible counts do not measure evidential strength.

Preparation and event select preserved common-string witness partitions and
conditional hulls shown in the inspector. The dependency map is shared; these
selectors do not discover preparation-specific causal mechanisms. AND gates
require all listed premises, including accepted-label correctness where needed.
Their display is not an automatic proof of molecular accuracy or function.

Rebuild from this component directory with Python 3.11 or later:

```sh
python3 -B code/build_reference_graph_model.py
python3 -B code/verify_reference_graph_model.py
python3 -B code/build_supergraph.py
python3 -B code/check_supergraph.py
```

`supergraph_model.json` is built separately from the canonical numerical model
and verified source files. The receipt distinguishes the source JSON file's
SHA-256 from the exact JSON serialization embedded in HTML. The parsed objects
are equal. Web Crypto verifies the embedded serialization when available; the
Python checker verifies both representations and byte-identical rebuilding.

The scientific observations and previous alignment calculations remain frozen.
This tool supplies navigation and computational audit, not new biological
replication, independent molecular truth or prospective transport validation.
See `SOURCE_NOTICE.md` and `LICENSE` for the exact published interface reference
and the GPL scope of this derived component.
