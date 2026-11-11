"""Build and test a native offline package on the current system. SPDX-License-Identifier: GPL-3.0-only."""

import argparse
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

from reproduce import sha256


ROOT = Path(__file__).resolve().parents[1]
PROJECT = Path(__file__).resolve().parent
PACKAGED_EXTENSIONS = ("retained_span_audit", "generation_model/validation")
RUNTIME_DISTRIBUTIONS = ("numpy", "numba", "llvmlite", "biopython", "pyinstaller")


def copy_notices(destination):
    destination.mkdir()
    shutil.copy2(ROOT / "LICENSE", destination / "GPL-3.0-only.txt")
    shutil.copy2(ROOT / "extensions" / "generation_model" / "LICENSE", destination / "MIT-extension.txt")
    for distribution_name in RUNTIME_DISTRIBUTIONS:
        dist = importlib.metadata.distribution(distribution_name)
        target = destination / distribution_name
        target.mkdir()
        for record in dist.files:
            if ".." not in record.parts and any(part.upper().startswith(("LICENSE", "COPYING", "NOTICE")) for part in record.parts):
                source = dist.locate_file(record)
                if source.is_file():
                    copy = target.joinpath(*record.parts)
                    copy.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, copy)
    for name in ("LICENSE.txt", "LICENSE"):
        license_path = Path(sys.base_prefix) / name
        if license_path.is_file():
            shutil.copy2(license_path, destination / "CPython-LICENSE.txt")
            break
    (destination / "RUNTIME.md").write_text(
        "# Runtime notices\n\n"
        "This package includes CPython (Python Software Foundation License), NumPy, "
        "Numba, llvmlite and Biopython. The PyInstaller bootloader exception applies "
        "to the wrapper. License notices installed by the dependencies are included here.\n\n"
        "CPython source and license: https://github.com/python/cpython/tree/v" +
        platform.python_version() + "\n"
        "Scientific data retain the source attribution and terms named in the evidence.\n",
        encoding="utf-8",
    )


def archive_package(bundle, output, platform_name):
    if platform_name.startswith("windows-"):
        archive = output / (bundle.name + ".zip")
        with zipfile.ZipFile(archive, "x", zipfile.ZIP_DEFLATED) as package:
            for path in sorted(bundle.rglob("*")):
                if path.is_file():
                    package.write(path, path.relative_to(bundle.parent))
        with zipfile.ZipFile(archive) as package:
            if package.testzip() is not None:
                raise ValueError("Native ZIP failed archive verification")
    else:
        archive = output / (bundle.name + ".tar.gz")
        with tarfile.open(archive, "x:gz") as package:
            package.add(bundle, arcname=bundle.name)
    return archive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--test-mode", choices=("quick", "counts", "full"), default="quick")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    operating_system = platform.system().lower()
    machine = {"amd64": "x86_64", "aarch64": "arm64"}.get(platform.machine().lower(), platform.machine().lower())
    platform_name = f"{operating_system}-{machine}"
    binary = "crispr-reproduce.exe" if operating_system == "windows" else "crispr-reproduce"
    with tempfile.TemporaryDirectory(prefix="crispr-native-") as temporary:
        temp = Path(temporary)
        subprocess.run(
            [sys.executable, "-m", "PyInstaller", "--noconfirm", "--onedir", "--noupx",
             "--name", "crispr-reproduce", "--distpath", str(temp / "dist"),
             "--workpath", str(temp / "build"), "--specpath", str(temp),
             "--hidden-import", "Bio.Seq", "--hidden-import", "numba",
             "--hidden-import", "numpy", str(PROJECT / "reproduce.py")],
            check=True, cwd=ROOT,
        )
        bundle = temp / f"crispr-reproduce-{platform_name}"
        (temp / "dist" / "crispr-reproduce").rename(bundle)
        for extension in PACKAGED_EXTENSIONS:
            shutil.copytree(ROOT / "extensions" / extension, bundle / "evidence" / "extensions" / extension,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for name in ("README.md", "LEEME_ES.md", "EVIDENCE.md"):
            shutil.copy2(PROJECT / name, bundle / name)
        source = bundle / "source"
        source.mkdir()
        for name in ("reproduce.py", "build_native.py", "requirements-native.txt", "requirements-macos-intel.txt"):
            shutil.copy2(PROJECT / name, source / name)
        copy_notices(bundle / "licenses")
        metadata = {
            "source_commit": args.source_commit,
            "platform": platform_name,
            "platform_version": platform.platform(),
            "python": platform.python_version(),
            "build_distributions": {d: importlib.metadata.version(d) for d in RUNTIME_DISTRIBUTIONS},
            "source_repository": "https://github.com/FacundoFirmenich/crispr-observability",
        }
        (bundle / "BUILD.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        archive = archive_package(bundle, output, platform_name)
        # Test the downloaded form in a fresh path containing spaces.
        test_root = temp / "extracted package"
        test_root.mkdir()
        shutil.unpack_archive(archive, test_root, filter="data" if operating_system != "windows" else None)
        executable = test_root / bundle.name / binary
        report_path = output / "replay"
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONHOME", None)
        env["PYTHONNOUSERSITE"] = "1"
        subprocess.run([str(executable), "--mode", args.test_mode, "--output", str(report_path)],
                       cwd=test_root, env=env, check=True)
        report = json.loads((report_path / "REPRODUCTION.json").read_bytes())
        if report["status"] != "PASS" or report["environment"]["frozen"] is not True:
            raise ValueError("Packaged executable did not complete the selected replay")
    (output / "SHA256SUMS.txt").write_text(f"{sha256(archive)}  {archive.name}\n", encoding="utf-8")
    print(json.dumps({"archive": archive.name, "platform": platform_name, "test_mode": args.test_mode}))


if __name__ == "__main__":
    main()
