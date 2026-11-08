#!/usr/bin/env python3
"""Merge the lcov files of the nightly coverage shards and print the job summary.

No NCBI counterpart: CI tooling for `.github/workflows/nightly.yml` (jobs `coverage` and
`coverage report`). Each shard runs part of the test set (`cargo llvm-cov nextest --partition
count:k/N`) and exports lcov (`cargo llvm-cov report --lcov`) with the outcome of its test step.
This script adds the shards' counts, writes the merged lcov file and a table per file, and prints
the totals and a table by module of `LOSAT/src` as Markdown.

- A line is a `DA` line of the file, counted once, covered when its count is above 0 in a shard;
  the merged count is the sum. llvm-cov's own table counts the lines of a closure for the closure
  and again for the function around it (the whole test set of 2026-11-07: 248,001 lines here,
  256,105 in `cargo llvm-cov report`).
- A function is the group of lcov records (instantiations of one generic function, the copies in
  each test binary) that start on one line of one file, covered when any of them ran in any shard.
  llvm-cov groups them by line and column; lcov gives no column, so functions that start on one
  line count once (2026-11-07: 19,939 here, 20,053 in `cargo llvm-cov report`).

Usage: coverage_summary.py --src LOSAT/src --shards 4 --out DIR [--commit SHA] SHARD_DIR...
Each SHARD_DIR (`coverage-shard-k`) may hold `lcov.info` and `outcome` (the test step's outcome,
one word); the summary names each of `coverage-shard-1` to `-N` that is absent or has no lcov.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class FileCoverage:
    """Counts of one source file: line -> count, function name -> [first line, count]."""

    lines: dict[int, int] = field(default_factory=dict)
    functions: dict[str, list[int]] = field(default_factory=dict)

    def line_counts(self) -> tuple[int, int]:
        """(covered, total) lines."""
        return sum(1 for c in self.lines.values() if c > 0), len(self.lines)

    def function_counts(self) -> tuple[int, int]:
        """(covered, total) functions, grouped by first line."""
        groups: dict[int, bool] = {}
        for line, count in self.functions.values():
            groups[line] = groups.get(line, False) or count > 0
        return sum(groups.values()), len(groups)


def read_lcov(text: str, files: dict[str, FileCoverage]) -> None:
    """Add the records of one lcov file to `files` (keyed by the SF path)."""
    current: FileCoverage | None = None
    named: set[str] = set()  # the functions this record's FN lines name
    for raw in text.splitlines():
        tag, _, value = raw.partition(":")
        if raw == "end_of_record":
            current = None
        elif tag == "SF":
            current = files.setdefault(value, FileCoverage())
            named = set()
        elif current is None:
            continue
        elif tag == "DA":
            line, count = value.split(",")[:2]
            current.lines[int(line)] = current.lines.get(int(line), 0) + int(count)
        elif tag == "FN":
            line, name = value.split(",", 1)
            current.functions.setdefault(name, [int(line), 0])
            named.add(name)
        elif tag == "FNDA":
            # llvm-cov writes a function's FN line before its FNDA line; a count without one has
            # no line to group it by.
            count, name = value.split(",", 1)
            if name in named:
                current.functions[name][1] += int(count)


def write_lcov(files: dict[str, FileCoverage]) -> str:
    """The merged counts as lcov, LF/LH and FNF/FNH counted as above."""
    out = ["TN:"]
    for path in sorted(files):
        cov = files[path]
        out.append(f"SF:{path}")
        for name, (line, _) in sorted(cov.functions.items(), key=lambda kv: (kv[1][0], kv[0])):
            out.append(f"FN:{line},{name}")
        for name, (_, count) in sorted(cov.functions.items(), key=lambda kv: (kv[1][0], kv[0])):
            out.append(f"FNDA:{count},{name}")
        fh, fnf = cov.function_counts()
        out += [f"FNF:{fnf}", f"FNH:{fh}"]
        for line in sorted(cov.lines):
            out.append(f"DA:{line},{cov.lines[line]}")
        lh, lf = cov.line_counts()
        out += [f"LF:{lf}", f"LH:{lh}", "end_of_record"]
    return "\n".join(out) + "\n"


def module_of(path: str, src: Path) -> str:
    """`algorithm/<program>` under `algorithm/`, else the first directory under `src`."""
    try:
        parts = Path(path).resolve().relative_to(src).parts
    except ValueError:
        return "(outside src)"
    if len(parts) == 1:
        return "src/*.rs"
    if parts[0] == "algorithm" and len(parts) > 2:
        return f"algorithm/{parts[1]}"
    return parts[0]


def shard_order(item: tuple[str, str]) -> tuple[int, str]:
    """`coverage-shard-2` before `coverage-shard-10`."""
    name = item[0]
    number = name.rsplit("-", 1)[-1]
    return (int(number) if number.isdigit() else 0, name)


def pct(covered: int, total: int) -> str:
    return f"{100 * covered / total:.1f}" if total else "-"


def summary(files: dict[str, FileCoverage], src: Path, commit: str, shards: list[str]) -> str:
    """The Markdown of the job summary."""
    out = ["### Coverage of LOSAT (report, not a gate)", ""]
    head = f"Commit {commit[:9]}; " if commit else ""
    out += [f"{head}shards: {', '.join(shards)}.", "",
            "Test files are left out; inline `#[cfg(test)]` modules of `src` files are counted.", ""]
    if any(not s.endswith(": success") for s in shards):
        out += ["A shard did not pass its tests or gave no report (see its log): the numbers "
                "count the tests that ran.", ""]
    lines = [c.line_counts() for c in files.values()]
    funcs = [c.function_counts() for c in files.values()]
    out += ["| | covered | total | % |", "|---|---:|---:|---:|"]
    for kind, counts in (("lines", lines), ("functions", funcs)):
        covered, total = sum(c for c, _ in counts), sum(t for _, t in counts)
        out.append(f"| {kind} | {covered} | {total} | {pct(covered, total)} |")
    groups: dict[str, list[int]] = {}
    for path, cov in files.items():
        g = groups.setdefault(module_of(path, src), [0, 0, 0, 0, 0])
        (lh, lf), (fh, fnf) = cov.line_counts(), cov.function_counts()
        g[0] += 1
        g[1] += lh
        g[2] += lf
        g[3] += fh
        g[4] += fnf
    out += ["", "By module of `LOSAT/src`, most uncovered lines first (per file: artifact "
            "`coverage`, `files.tsv`):", "",
            "| module | files | lines % | uncovered lines | functions % |",
            "|---|---:|---:|---:|---:|"]
    for key, g in sorted(groups.items(), key=lambda kv: (kv[1][1] - kv[1][2], kv[0])):
        out.append(f"| `{key}` | {g[0]} | {pct(g[1], g[2])} | {g[2] - g[1]} | {pct(g[3], g[4])} |")
    return "\n".join(out) + "\n"


def files_table(files: dict[str, FileCoverage], src: Path) -> str:
    """One row per file: path under `src` (or as given), lines and functions covered / total."""
    out = ["file\tlines_covered\tlines\tfunctions_covered\tfunctions"]
    rows = []
    for path, cov in files.items():
        try:
            name = Path(path).resolve().relative_to(src).as_posix()
        except ValueError:
            name = path
        (lh, lf), (fh, fnf) = cov.line_counts(), cov.function_counts()
        rows.append(f"{name}\t{lh}\t{lf}\t{fh}\t{fnf}")
    return "\n".join(out + sorted(rows)) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--src", type=Path, required=True, help="LOSAT/src")
    parser.add_argument("--shards", type=int, required=True, help="number of shards run")
    parser.add_argument("--out", type=Path, required=True, help="directory for lcov.info, files.tsv")
    parser.add_argument("--commit", default="")
    parser.add_argument("shard_dirs", nargs="*", type=Path)
    args = parser.parse_args(argv)
    src = args.src.resolve()
    files: dict[str, FileCoverage] = {}
    found: dict[str, str] = {}
    for shard_dir in args.shard_dirs:
        if not shard_dir.is_dir():  # an unmatched glob reaches here as the pattern itself
            continue
        lcov = shard_dir / "lcov.info"
        outcome = shard_dir / "outcome"
        state = outcome.read_text().strip() if outcome.is_file() else "unknown"
        if lcov.is_file():
            read_lcov(lcov.read_text(), files)
        else:
            state = "no report"
        found[shard_dir.name] = state
    for k in range(1, args.shards + 1):
        found.setdefault(f"coverage-shard-{k}", "no report")
    shards = [f"{name}: {state}" for name, state in sorted(found.items(), key=shard_order)]
    print(summary(files, src, args.commit, shards), end="")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "lcov.info").write_text(write_lcov(files))
    (args.out / "files.tsv").write_text(files_table(files, src))
    return 0 if files else 1


if __name__ == "__main__":
    sys.exit(main())
