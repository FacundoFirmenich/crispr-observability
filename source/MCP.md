# Read Observability Lab · MCP

The local stdio bridge runs four tools and renders the same standalone HTML used by the web application. It takes an explicit workspace directory. File arguments are relative to that directory and resolved before use. PNG, replay and manifest are written to a unique subdirectory. The renderer loads the local HTML and blocks HTTP/HTTPS requests.

## Start from an IDE

1. Extract the full package into a directory such as `C:/research/read-observability-lab`.
2. Run `npm ci` in that directory. Use Node.js 22 or newer and an installed Google Chrome for the render tool. The calculation tools do not launch Chrome.
3. Configure an MCP **stdio** server with command `node` and arguments `[absolute path to mcp-server.mjs, absolute path to your data workspace]`.

In Codex, the equivalent CLI registration is:

```sh
codex mcp add read-observability -- node /absolute/path/mcp-server.mjs /absolute/path/data-workspace
```

For an IDE with a JSON MCP configuration:

```json
{
  "mcpServers": {
    "read-observability": {
      "command": "node",
      "args": ["C:/research/read-observability-lab/mcp-server.mjs", "C:/research/data"]
    }
  }
}
```

| Tool | Input and output |
| --- | --- |
| `lab_boundary` | Anchor start/end, displacement and read end. Returns exact interval containment and observed-base count. |
| `lab_compare_spans` | Two half-open intervals on an archived read or supplied sequence. Returns strings, overlap and exclusive lengths. |
| `lab_normalize_case` | Workspace-relative `case_path`, or `case_name` with inline `case_text`. Returns the SHA-256, format and counts. |
| `lab_render` | A built-in preset, inline replay or supplied case, and scale 1, 2 or 4. Returns the output folder, dimensions, result and hashes. |

To reproduce the archived one-base witness, call `lab_boundary` with `{"start":210,"end":222,"shift":3,"read_end":224}`. The returned count is 11. Set `read_end` to 225 for the counterfactual.

To create a 4× detail from the archived CFTR model, call `lab_render` with:

```json
{
  "scale": 4,
  "replay": {
    "schema": "read-observability-result/1",
    "kind": "ideal_duplex_inspection",
    "inputs": {
      "start": 129,
      "end": 141,
      "strand": "both",
      "representation": "bond"
    }
  }
}
```

For a custom angle or zoom, first save the replay JSON from the browser, then pass that JSON object as the `replay` argument. The tool writes `scene.png`, `replay.json` and `manifest.json`. These are independent of the browser's screen resolution up to the render cap. The manifest links each file to its SHA-256. Imported scientific content remains in the workspace used by the MCP host. This implementation uses a local stdio transport and does not expose a public MCP server.
