# Integrated source identity

The browser source is the locally preserved v0.3 alignment-capable laboratory,
integrated with its source modules and validation scripts. The multi-input engine
under `../../engine/` is imported from commit
`d6a12e4dae968250d6791fdee281ea15c1cfd25b` and includes the v0.5 junction operation.
The engine and browser retain distinct schemas and evidence scopes.

Current integration CI overlays these source modules onto the frozen public
fixtures before rebuilding and testing. A historical test of code inside a
fixture ZIP is not used to certify a different source tree.

Original laboratory code uses GPL-3.0-only under the repository licence.
Bundled Three.js and other dependencies retain their own licence notices.
