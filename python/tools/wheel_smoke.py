"""One search through an installed pylosat wheel (.github/workflows/python-wheels.yml).

Run it with the Python of a fresh virtual environment that holds the wheel alone, with `-I` and
from outside the checkout, so that `import pylosat` finds the wheel's package and not
python/pylosat:

    <venv>/bin/python -I <checkout>/python/tools/wheel_smoke.py --repo <checkout> --wheel <wheel>

The search is the row `blastp pairwise_default_serial` of the command line's frozen hashes
(docs/evidence/losat_web_e1a/baseline/hashes.tsv; ci.yml checks the command line against them on
every pull request). Its words go to `Search.from_argv` as they are, but for the inputs and the
format, which are arguments of `run`, and `-out`. The report must have the row's frozen hash and
the diagnostics that of the row's standard error. Prints one Markdown table row for the job
summary: the Python, the wheel's file name, tags and size, and the hash.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import sys
import zipfile
from pathlib import Path

HASHES = Path("docs/evidence/losat_web_e1a/baseline/hashes.tsv")
CASE = ("blastp", "pairwise_default_serial")
FIXTURES = "LOSAT/tests/fasta/"
ARGUMENTS_OF_RUN = ("-query", "-subject", "-outfmt", "-out")


def frozen_row(repo: Path) -> dict[str, str]:
    with open(repo / HASHES, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row["program"], row["case_id"]) == CASE:
                if row["exit"] != "0" or row["matches_expected"] != "true":
                    raise SystemExit(f"{HASHES}: {CASE} is not a passing frozen case")
                return row
    raise SystemExit(f"{HASHES}: no row {CASE}")


def wheel_metadata(wheel: Path) -> tuple[str, list[str]]:
    """The wheel's version and tags, from its .dist-info."""
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()

        def read(suffix: str) -> list[str]:
            (name,) = [name for name in names if name.endswith(f".dist-info/{suffix}")]
            return archive.read(name).decode().splitlines()

        tags = [line.split(":", 1)[1].strip() for line in read("WHEEL") if line.startswith("Tag:")]
        (version,) = [line.split(":", 1)[1].strip() for line in read("METADATA") if line.startswith("Version:")]
    return version, tags


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repo", type=Path, required=True, help="the checkout: fixtures and frozen hashes")
    parser.add_argument("--wheel", type=Path, required=True, help="the wheel installed in this environment")
    args = parser.parse_args()

    row = frozen_row(args.repo)
    program, *words = json.loads(row["command"])
    named: dict[str, str] = {}
    options: list[str] = []
    rest = iter(words)
    for word in rest:
        if word in ARGUMENTS_OF_RUN:
            named[word] = next(rest)
        else:
            options.append(word)

    def fixture(path: str) -> Path:
        return args.repo / FIXTURES / path.rsplit(FIXTURES, 1)[1]

    version, tags = wheel_metadata(args.wheel)
    if sys.prefix == sys.base_prefix:
        raise SystemExit("not in a virtual environment")

    import pylosat

    module = Path(pylosat.__file__).resolve()
    if not module.is_relative_to(Path(sys.prefix).resolve()):
        raise SystemExit(f"pylosat was imported from {module}, not from this environment ({sys.prefix})")
    installed = importlib.metadata.version("pylosat")
    if installed != version:
        raise SystemExit(f"pylosat {installed} is installed, the wheel is {version}")

    search = pylosat.Search.from_argv(program, options)
    outcome = search.run(
        pylosat.Sequences.from_fasta(fixture(named["-query"])),
        pylosat.Sequences.from_fasta(fixture(named["-subject"])),
        outfmt=named["-outfmt"],
    )
    report = hashlib.sha256(outcome.reports[0].bytes).hexdigest()
    diagnostics = hashlib.sha256(outcome.diagnostics.encode()).hexdigest()
    failures = []
    if report != row["expected_sha256"]:
        failures.append(f"report sha256 {report}, frozen {row['expected_sha256']}")
    if diagnostics != row["stderr_sha256"]:
        failures.append(f"diagnostics {outcome.diagnostics!r}, frozen stderr sha256 {row['stderr_sha256']}")
    size = args.wheel.stat().st_size
    verdict = "**differs**" if failures else "= frozen"
    print(
        f"| {platform.python_implementation()} {platform.python_version()} | `{args.wheel.name}` | "
        f"{', '.join(f'`{tag}`' for tag in tags)} | {size:,} B ({size / 2**20:.1f} MiB) | "
        f"`{report[:12]}` {verdict} |"
    )
    if failures:
        raise SystemExit(f"{HASHES.as_posix()}, {' '.join(CASE)}: " + "; ".join(failures))


if __name__ == "__main__":
    main()
