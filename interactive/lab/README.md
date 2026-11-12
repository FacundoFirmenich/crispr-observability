# Read Observability Lab

Open `CFTR_OBSERVABILITY_LAB.html` directly in a modern WebGL 2 browser. The single file embeds the application, its graphics library and its CFTR data. It can run offline or from any static HTTPS host. Imported files are processed in browser memory. The application makes no network requests for imports or calculations.

**Open walkthrough** starts five guided activities. Load an example, vary its controls and check the live checkpoint. **What these controls change** explains each experiment’s parameters. In the coordinate rail, moving the position slider highlights the corresponding object in 3D.

## Experiments

| Experiment | Change | Inspect |
| --- | --- | --- |
| Read boundary | Anchor start and end, displacement, exclusive read end | Anchor coverage, first missing index, minimum required length, a threshold sweep and a pinned prior result |
| Competing spans | Archived read, extraction rule, alternate optimal span, interval coordinates, orientation | Exact overlap, each extracted string and the status of the archived labels |
| Atomic model | Crop, strand visibility and atom representation | Original model coordinates, visible atoms and covalent links |
| Imported structure | File, chain, two selected atoms and distance threshold | The measured distance, threshold decision and portable structure |
| Read alignments | SAM or indexed BAM locus, record, orientation and MAPQ threshold | CIGAR mapping, clipping, reference depth and source provenance |

In the read-boundary experiment the baseline is the unshifted interval. The frozen treated-minus-mock contrasts are shown separately from the coordinate counterfactual. In the archived comparison the first record produces 42 and 87 nt under two rules. On the second record, the default 72- and 39-nt spans differ although the labels agree. A second optimal placement makes the compared spans coincide. The structural CFTR model uses the observed 42-base string and AmberTools ideal B-DNA coordinates. The two scaffold experiments use their own coordinate geometry. The structure is inspectable in three dimensions, and the model yields no new sequencing or editing observations.

## Open and convert a case

- **DNA FASTA** (`.fa`, `.fasta`, `.fna`) accepts one record, 2–360 bases in A/C/G/T/N. Its sequence enters the same span-comparison engine. Custom intervals are zero-based and half-open. N complements to N.
- **PDB** accepts the first model and blank or A alternate-location records. The lab imports ATOM/HETATM Cartesian positions and explicit CONECT relationships. PDB files often omit polymer connectivity. The portable export preserves imported positions and explicit bonds, with no inferred topology.
- **PDBx/mmCIF** (`.cif`, `.mmcif`) reads the first model in the textual `atom_site` loop, in ångströms. It uses author IDs when present and blank/unknown or A alternate locations. Import currently covers coordinates and atom identity; it does not apply assembly operators or recover connectivity from other categories.
- **XYZ** accepts one frame and Cartesian positions in ångströms. Files declare the atom count on line one, a comment on line two, then an element and three coordinates per atom. XYZ supplies no bonds here.
- **Portable case JSON** accepts the versioned contract in `case.schema.json`. The browser converter emits that contract and accepts it back. Import limits are 6 MB and 30,000 atoms. The converter preserves the supported fields and is not a lossless archiver of all PDB or XYZ metadata.

Click a point to inspect its position. In imported structures, choose two atoms by clicking or by visible index, then vary the distance threshold. Use **Download portable case** to retain the normalized input and **Download current experiment** to retain parameters and results. An archived result includes source hashes. Imported data, conversion and calculations stay on the device until the user explicitly downloads or shares a file.

## Replay and export

Open an exported experiment JSON with **Open a case or replay** to reconstruct the controls, case and camera. The app recalculates the outputs from those inputs. Imported structures embed their normalized coordinates in the replay. Keep the file private if its case is private.

**Export PNG** samples the current WebGL view at 2× or 4×, with a 4,096-pixel long-side cap. Orbit and zoom before exporting, or use the atom start/end controls to isolate a segment. The JSON carries calculation provenance separately. Saving the same zoom in replay allows the view to be generated again.

See [MCP.md](MCP.md) for local IDE integration, workspace-restricted case loading and render exports with hashes.

See [ALIGNMENT_AUDIT.md](ALIGNMENT_AUDIT.md) for indexed BAM, SAM, precise coordinates, coverage, export schemas and the independent oracle.

## Build and source

Install Node dependencies with `npm ci`, then run `node build.mjs` in this directory. The frozen inputs in `data.json` make this standalone. `build_data.py` records how the first CFTR case was derived from the companion project artifacts. Run the nine scientific checks from the browser. The dependency versions are pinned in `package-lock.json`.

The experiment kernels in `science.js` are independent of the browser UI. `import-case.js` defines the case conversion and uses the same normalized record contract across sources. `molecular-view.js` owns the 3D renderer and navigation. Each layer can be reused independently by a future host.

- [Prime editing mechanism, Anzalone et al. (2019)](https://doi.org/10.1038/s41586-019-1711-4)
- [AmberTools B-DNA constructor implementation](https://github.com/Amber-MD/AmberClassic/blob/main/src/nabc/fd_helix.c)
- [Three.js documentation](https://threejs.org/docs/)
