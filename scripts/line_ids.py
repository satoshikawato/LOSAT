#!/usr/bin/env python3
"""Line and stage IDs that the repository cites, and their rows in docs/decisions/LINE_IDS.md.

Code comments, tests, changelog fragments and documents cite the maintainer's work lines and their
stages by short IDs ("line DBID, N1b", "DI-8f", "OF8a", "plan v4 line R", "SF S3"). The work records
behind them are outside the repository; `docs/decisions/LINE_IDS.md` gives each ID one row.

Usage:
  line_ids.py list [--all]   every cited ID with its number of citing lines and first location
                             (CHANGELOG.md, changelog.d/ and docs/; --all: the whole tree)
  line_ids.py check [--strict]
                             IDs cited in CHANGELOG.md, changelog.d/ and docs/ that have no row in
                             LINE_IDS.md, as GitHub warnings; exit status 0 (1 with --strict)

Not counted: families with a registry of their own in the tree (DW-n and TD-n in
docs/losat_web_gui_plan.md, REQ-nn in docs/web/requirements_trace.tsv, PD-* in docs/product_decisions/,
LIB-n in docs/losat_library_surface_plan.md).
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = ROOT / "docs" / "decisions"  # the registries themselves; not citations
REGISTRY = DECISIONS / "LINE_IDS.md"
EXTS = (".md", ".tsv", ".txt", ".json", ".py", ".sh", ".yml", ".yaml", ".rs", ".ts", ".vue", ".toml",
        ".js", ".mjs", ".html", ".css")
SKIP_DIRS = {"node_modules", ".git", "target", "dist"}
SKIP_FILES = {"package-lock.json", "Cargo.lock"}

# Lines of plan v4, cited as "line R", "plan v4 line W" or "v4-MEM".
V4 = {"A", "B", "C", "D", "E", "F", "G", "R", "S", "W", "MF", "MEM", "NL", "BPF", "P0",
      "A10", "A12", "A14", "A15", "B10", "D12", "D14", "R6"}
LINES = {"SF", "NP", "DI", "DBID", "BDC", "OF", "CLI", "PS", "GZ", "WF", "WA", "WR", "WO", "WG", "WFL",
         "SX", "TS", "ST", "EB", "BW", "LIB", "RL", "NCMEM", "SW", "MK"}


def _line_mention(m: re.Match) -> str | None:
    """'line X', 'Line X, S1', '線 X', 'plan v4 line W', 'structural line S3'."""
    pre, code, stage = m.group("pre") or "", m.group("code"), m.group("stage")
    if "structural" in pre and re.fullmatch(r"S[0-9]+", code):
        return "SF " + code
    if re.fullmatch(r"DI-?[0-9]+[a-z]?", code):
        return "DI-" + code.replace("DI-", "").replace("DI", "")
    if re.fullmatch(r"OF[0-9]+[a-z]?|RL[0-9]+|BDC[0-9]+[a-z]?|NCMEM[0-9]+|WF-[0-9]+", code):
        return code
    if code not in V4 and code not in LINES:
        return None  # prose: "line N of the file", "line NCBI prints"
    base = "v4-" + code if code in V4 else code
    if not stage or stage.startswith(("s", "item")) or stage == "V8":  # V8: the JavaScript engine
        return base
    if code in ("OF", "PS", "BDC", "RL") and stage.startswith(code):
        return stage
    return base + " " + stage


PATTERNS: list[tuple[re.Pattern, object]] = [
    (r"(?P<pre>(?:plan v4 |v4 |structural )?)(?:\b[Ll]ines?|線) (?P<code>[A-Z][A-Za-z0-9]*(?:-[0-9]+[a-z]?)?)"
     r"(?:,? (?:item )?(?P<stage>[A-Z]{1,4}[0-9]+[a-z]?|s[0-9]+|session [0-9]+))?", _line_mention),
    (r"\bDI-?([0-9]+[a-z]?)\b", lambda m: "DI-" + m.group(1)),
    (r"\bDBID (?:item )?([A-Z]{1,2}[0-9]+[a-z]?)\b", lambda m: "DBID " + m.group(1)),
    (r"\bBDC([0-9]+[a-z]?)\b", lambda m: "BDC" + m.group(1)),
    (r"\bOF([0-9]+[a-z]?)\b", lambda m: "OF" + m.group(1)),
    (r"\bSF S([0-9]+)\b", lambda m: "SF S" + m.group(1)),
    (r"\bNP (M[0-9]+|WW[0-9]+|SC[0-9]+[a-z]?|SM[0-9]+)\b", lambda m: "NP " + m.group(1)),
    (r"(?:\bNP,? |\bline NP, |parity line )(O[0-9]+[a-z](?:-[0-9]+)?)\b", lambda m: "NP " + m.group(1).split("-")[0]),
    (r"\b(O[0-9]+[a-z])-[0-9]+\b", lambda m: "NP " + m.group(1)),
    (r"\bSM-([0-9]+)\b", lambda m: "NP SM-" + m.group(1)),
    (r"(?:series SC, |\bNP SC|\(NP )(SC[0-9]+[a-z]?)\b", lambda m: "NP " + m.group(1)),
    (r"\bSC([0-9]+[a-z]?) (?:of line NP|\(line NP)", lambda m: "NP SC" + m.group(1)),
    (r"\bNCMEM([0-9]+)\b", lambda m: "NCMEM" + m.group(1)),
    (r"\bWF-([0-9]+)\b", lambda m: "WF-" + m.group(1)),
    (r"\bWA2?-([A-Z]{2})-([0-9]+)\b", lambda m: m.group(0)),
    (r"\bWA[- ](W[0-9]+)\b", lambda m: "WA " + m.group(1)),
    (r"\bSX([a-d][0-9]?)\b", lambda m: "SX" + m.group(1)),
    (r"\bv4-([A-Z][A-Za-z0-9]*)\b", lambda m: "v4-" + m.group(1) if m.group(1) in V4 else None),
    (r"\bCLI (C[0-9]+[a-z]?)\b", lambda m: "CLI " + m.group(1)),
    (r"(?:line PS,? |PS item |\(PS)(PS[0-9]+)\b", lambda m: m.group(1)),
    (r"\bRL([0-9]+)\b", lambda m: "RL" + m.group(1)),
    (r"\bmemory line \(v4-MEM\)|\bline MEM\b", lambda m: "v4-MEM"),
]
PATTERNS = [(re.compile(rx), fn) for rx, fn in PATTERNS]


def files(whole_tree: bool) -> list[Path]:
    if not whole_tree:
        bases = [ROOT / "changelog.d", ROOT / "docs"]
        out = [ROOT / "CHANGELOG.md"]
    else:
        bases, out = [ROOT], []
    for base in bases:
        for d, ds, fs in os.walk(base):
            rel = Path(d).relative_to(ROOT).as_posix()
            if SKIP_DIRS & set(rel.split("/")) or rel.startswith("LOSAT/tests/fixtures"):
                ds[:] = []
                continue
            out += [Path(d) / f for f in fs if f.endswith(EXTS) and f not in SKIP_FILES]
    return sorted(p for p in out if p.is_file() and DECISIONS not in p.parents)


def cited(paths: list[Path]) -> dict[str, list[str]]:
    """ID -> ['path:line', ...] (one entry per citing line)."""
    found: dict[str, list[str]] = {}
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="ignore")
        for number, line in enumerate(text.splitlines(), 1):
            ids = set()
            for rx, fn in PATTERNS:
                for m in rx.finditer(line):
                    if (cid := fn(m)):
                        ids.add(cid)
            for cid in ids:
                found.setdefault(cid, []).append(f"{path.relative_to(ROOT).as_posix()}:{number}")
    return found


def registered() -> set[str]:
    """IDs in the first column of the tables of LINE_IDS.md (`| \\`ID\\` | ...`)."""
    if not REGISTRY.exists():
        return set()
    return set(re.findall(r"^\| `([^`]+)` \|", REGISTRY.read_text(encoding="utf-8"), re.M))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)
    sub.add_parser("list").add_argument("--all", action="store_true")
    sub.add_parser("check").add_argument("--strict", action="store_true")
    args = parser.parse_args()
    if args.mode == "list":
        found = cited(files(args.all))
        for cid in sorted(found):
            print(f"{cid}\t{len(found[cid])}\t{found[cid][0]}")
        print(f"{len(found)} IDs", file=sys.stderr)
        return 0
    found = cited(files(False))
    missing = sorted(set(found) - registered())
    for cid in missing:
        path, line = found[cid][0].rsplit(":", 1)
        print(f"::warning file={path},line={line}::{cid} has no row in docs/decisions/LINE_IDS.md")
    print(f"{len(found)} line IDs cited in CHANGELOG.md, changelog.d/ and docs/; {len(missing)} without a row")
    return 1 if missing and args.strict else 0


if __name__ == "__main__":
    sys.exit(main())
