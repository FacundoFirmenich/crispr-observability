# Offline reproduction

Extract the package for your operating system and architecture. Open a terminal
in the extracted folder and run:

| Platform | Command |
| --- | --- |
| Windows x86-64 (PowerShell) | `.\crispr-reproduce.exe --output .\replay-quick` |
| Linux x86-64 | `./crispr-reproduce --output ./replay-quick` |
| macOS Apple silicon or Intel | `./crispr-reproduce --output ./replay-quick` |

The package includes its Python runtime and frozen public inputs. Keep its
`evidence` and runtime directories alongside the executable. The output
directory must be new. Validate the archive against `SHA256SUMS.txt` from
the same release. Each successful replay writes `REPRODUCTION.json` with
the executed mode, environment, count and evidence scope.

- **quick** checks frozen files and graph links, reconstructs the exact CFTR
  contrast endpoints from archived counts, reruns the authentic Clone3
  example, and reconstructs 241 ONT locus decisions and 113 matched HiFi
  decisions from stored scores and flank calls.
- **counts** adds the 439,589 retained COL17A1 pair ledger counts in place of
  the single example.
- **full** replays the retained COL17A1 counts, 20,260 reference cases and
  81,112 attaining witnesses and runs the bounded independent path oracle.
  The germline reconstruction has the same scope in all modes.

Original FASTQ extraction, fresh ONT/HiFi read scoring and independent
original-fragment editing accuracy have their own evidence units and are
outside these offline commands. MP33–MP38 protocols and results in the
bundle are checked for file identity. Their full remote acquisitions and
executions are identified separately.

The earlier v0.8.1 source archive has its own [Zenodo DOI](https://doi.org/10.5281/zenodo.24121127). The native release is identified by its version and the source commit in BUILD.json.

## Run or build from source

From the repository root, with CPython 3.14 and a fresh environment.
On Intel macOS use CPython 3.13 and the pinned requirements-macos-intel.txt instead.

```sh
python -m pip install -r reproducibility/requirements-native.txt
python reproducibility/reproduce.py --mode full --output replay-full
python reproducibility/build_native.py --source-commit "$(git rev-parse HEAD)" --output native-build --test-mode full
```

Build on the target OS and architecture. The workflow uses native hosts,
tests the executable from the archive and records exact output digests.
Linux builds target glibc on Ubuntu 22.04. macOS uses the GitHub Actions
macOS 15 runners. The executables have no project-issued Windows code
signature or Apple notarization. Source execution is also supported.
