#!/usr/bin/env python3
"""Which side of each conditional performance switch the CI fixtures reached (LOSAT_X_REACH).

A conditional switch of the strict set (LOSAT/src/utils/xswitch.rs STRICT_SET) takes its own path
only where a rule holds (input sizes, thread count, options, memory); elsewhere the reference path
runs although the switch is on. A conditional or platform switch that graduated (GRADUATED in the
same file) has no variable any more but keeps its rule. Runs with LOSAT_X_REACH=<directory> write
one file per LOSAT process (LOSAT/src/utils/xreach.rs). This script folds the files of a directory
and compares them with switch_reach_expected.tsv, one row per switch of the set and per graduated
switch:

  switch    the variable, as in STRICT_SET or GRADUATED
  programs  the LOSAT subcommands whose runs reach it (comma-separated)
  class     always   acts on every run of its programs; not recorded
            platform chooses a kernel by CPU or target only; covered by the platform-checks matrix
            conditional  recorded at its decision points (x_reach!)
  new, ref  for a conditional switch: required (some fixture must reach the side) or gap (no fixture
            reaches it yet; the note says why)
  note      free text

It fails when a switch of the set or a graduated one has no row, a row names neither, a graduated
switch's row says always, or a required side of a conditional switch was not reached. A conditional row is checked only when every program it
names was run (--programs); a gap that was reached is reported, so the row can say required.

Usage:
  switch_reach_check.py --dir DIR --programs all|LIST [--expected TSV]
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parents[1]
XSWITCH = ROOT / "LOSAT/src/utils/xswitch.rs"
EXPECTED = TESTS / "switch_reach_expected.tsv"
PROGRAMS = ("blastn", "blastp", "blastx", "tblastn", "tblastx")
CLASSES = ("always", "platform", "conditional")
SIDES = ("required", "gap")


def names_of(list_name: str) -> list[str]:
    text = XSWITCH.read_text()
    body = re.search(list_name + r": &\[&str\] = &\[(.*?)\];", text, re.S)
    if body is None:
        raise SystemExit(f"{list_name} not found in {XSWITCH}")
    return re.findall(r'"(LOSAT_X_[A-Z0-9_]+)"', body.group(1))


def strict_set() -> list[str]:
    return names_of("STRICT_SET")


def graduated() -> list[str]:
    return names_of("GRADUATED")


def opt_in() -> list[str]:
    """The opt-in switches (OPT_IN in xswitch.rs): off unless their variable is set, no row until they join the set."""
    return names_of("OPT_IN")


def read_expected(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))


def read_reach(directory: Path) -> dict[str, dict[str, set[str]]]:
    """switch -> side ("new"/"ref") -> programs that reached it."""
    reached: dict[str, dict[str, set[str]]] = {}
    for path in sorted(directory.glob("reach-*.tsv")):
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sides = reached.setdefault(row["switch"], {"new": set(), "ref": set()})
                for side in ("new", "ref"):
                    if row[side] == "1":
                        sides[side].add(row["program"])
    return reached


def check(names: list[str], rows: list[dict[str, str]], reached, programs: set[str],
          graduated_names: tuple[str, ...] = ()) -> tuple[list[str], list[str]]:
    failures: list[str] = []
    notes: list[str] = []
    names = [*names, *graduated_names]
    by_name = {}
    for row in rows:
        name = row["switch"]
        if name in by_name:
            failures.append(f"{name}: two rows")
        by_name[name] = row
        if name not in names:
            failures.append(f"{name}: row for a switch in neither STRICT_SET nor GRADUATED (remove the row)")
        if name in graduated_names and row["class"] == "always":
            failures.append(f"{name}: graduated with class always (an always-on switch leaves the file when it graduates)")
        if row["class"] not in CLASSES:
            failures.append(f"{name}: class {row['class']!r} is not one of {CLASSES}")
        listed = {p for p in row["programs"].split(",") if p}
        if not listed or not listed <= set(PROGRAMS):
            failures.append(f"{name}: programs {row['programs']!r} must name LOSAT subcommands {PROGRAMS}")
        if row["class"] == "conditional" and not {row["new"], row["ref"]} <= set(SIDES):
            failures.append(f"{name}: new/ref must be one of {SIDES}")
    for name in names:
        if name not in by_name:
            failures.append(f"{name}: no row in {EXPECTED.name} (add one; a conditional switch also needs x_reach! "
                            "at its decision points)")
    for name, row in sorted(by_name.items()):
        sides = reached.get(name, {"new": set(), "ref": set()})
        if row["class"] != "conditional":
            if sides["new"] or sides["ref"]:
                notes.append(f"{name}: recorded although its class is {row['class']}")
            continue
        listed = set(row["programs"].split(","))
        if not listed <= programs:
            continue
        for side in ("new", "ref"):
            if row[side] == "required" and not sides[side]:
                what = "its own path" if side == "new" else "the reference path (rule not met)"
                failures.append(f"{name}: no fixture reached {what}")
            elif row[side] == "gap" and sides[side]:
                notes.append(f"{name}: the {side} side was reached ({','.join(sorted(sides[side]))}); "
                             "the row can say required")
    unknown = sorted(set(reached) - set(names) - set(opt_in()))
    for name in unknown:
        failures.append(f"{name}: recorded but in neither STRICT_SET nor GRADUATED")
    return failures, notes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", type=Path, required=True)
    parser.add_argument("--programs", default="all")
    parser.add_argument("--expected", type=Path, default=EXPECTED)
    args = parser.parse_args()
    if args.programs == "all":
        programs = set(PROGRAMS)
    else:
        programs = {p for p in args.programs.replace(" ", ",").split(",") if p}
    if not args.dir.is_dir():
        print(f"switch_reach_check: {args.dir} is not a directory (no LOSAT run wrote a reach file)")
        return 1
    failures, notes = check(strict_set(), read_expected(args.expected), read_reach(args.dir), programs,
                            tuple(graduated()))
    for note in notes:
        print(f"note: {note}")
    for failure in failures:
        print(f"FAIL: {failure}")
    print(f"switch_reach_check: {len(failures)} failure(s), {len(notes)} note(s); programs {','.join(sorted(programs))}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
