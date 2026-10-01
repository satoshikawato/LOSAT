#!/usr/bin/env python3
"""E2g stage 1: map LOSAT's NCBI reference annotations on the BLASTN path to NCBI functions.

1. Files: the Rust modules reachable from LOSAT/src/algorithm/blastn/ through `use`
   paths and `crate::`/`super::`/`self::` paths (and paths that start with a child
   module), not entering the directories of the other programs
   (algorithm/{blastp,blastx,tblastn,tblastx}).
2. Annotations: every line that contains "NCBI reference". The referenced NCBI files
   and line ranges are the `<file>.<ext>:<lines>` tokens of that line, or, if it has
   none, of the next three comment lines (the snippet usually follows).
3. Functions: each referenced NCBI file (resolved under the pinned NCBI checkout;
   a bare file name is looked up by name, preferring algo/blast) is scanned for
   function bodies, and the first line of each range is mapped to the enclosing
   function (or to the enclosing struct, macro or file scope).

The first check reproduces the stage-1 counts pinned in EXPECTED (see README.md).

Usage: inventory_refs.py [--ncbi DIR] [--out DIR] [--no-check]
Writes stage1_refs.tsv, stage1_functions.tsv and stage1_summary.json to --out.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SRC = REPO / "LOSAT/src"
OTHER_PROGRAMS = tuple(f"algorithm/{p}/" for p in ("blastp", "blastx", "tblastn", "tblastx"))
DEFAULT_NCBI = Path("/mnt/c/Users/genom/GitHub/ncbi-blast")
# Pinned on 2026-10-02 at the S07++ engine (LOSAT/ unchanged since 7693a9c73).
EXPECTED = {"files": 54, "annotations": 1829, "ncbi_functions": 390}
REF = re.compile(r"((?:[A-Za-z0-9_.+\-]+/)*[A-Za-z0-9_+\-]+\.(?:cpp|hpp|inl|prt|c|h))(?![A-Za-z0-9_])"
                 r"(?::(\d+)(?:\s*[-–]\s*(\d+))?)?")


# ---------------------------------------------------------------- reachability
def module_file(parts: list[str]) -> Path | None:
    for k in range(len(parts), 0, -1):
        base = SRC.joinpath(*parts[:k])
        for candidate in (base.with_suffix(".rs"), base / "mod.rs"):
            if candidate.is_file():
                return candidate
    return None


def module_parts(path: Path) -> list[str]:
    parts = list(path.relative_to(SRC).with_suffix("").parts)
    if parts[-1] == "mod":
        parts = parts[:-1]
    return [] if parts in (["lib"], ["main"]) else parts


def expand_use(body: str) -> list[list[str]]:
    body = re.sub(r"\s+", "", body)
    out: list[list[str]] = []

    def rec(prefix: list[str], text: str) -> None:
        depth, current, items = 0, "", []
        for ch in text:
            depth += ch == "{"
            depth -= ch == "}"
            if ch == "," and depth == 0:
                items.append(current)
                current = ""
            else:
                current += ch
        if current:
            items.append(current)
        for item in items:
            if "{" in item:
                i = item.index("{")
                rec(prefix + [x for x in item[:i].rstrip(":").split("::") if x], item[i + 1:-1])
            else:
                item = item.split("as")[0] if re.search(r"[a-z]as[A-Z_a-z]", item) is None else item
                out.append(prefix + [x for x in item.split("::") if x and x not in ("*", "self")])

    rec([], body)
    return out


def references(path: Path) -> set[Path]:
    code = re.sub(r"//[^\n]*", "", path.read_text(errors="replace"))
    me = module_parts(path)
    is_mod = path.name in ("mod.rs", "lib.rs", "main.rs")
    paths = []
    for match in re.finditer(r"\buse\s+([^;]+);", code):
        paths += expand_use(match.group(1))
    for match in re.finditer(r"\b((?:crate|super|self)(?:::[A-Za-z_][A-Za-z0-9_]*)+)", code):
        paths.append(match.group(1).split("::"))
    found = set()
    for p in paths:
        if not p:
            continue
        if p[0] == "crate":
            parts = p[1:]
        elif p[0] in ("super", "self"):
            base = me if is_mod else me[:-1]
            rest = list(p)
            while rest and rest[0] == "super":
                base, rest = base[:-1], rest[1:]
            if rest and rest[0] == "self":
                rest = rest[1:]
            parts = base + rest
        else:
            child = SRC.joinpath(*(me + [p[0]]))
            if child.with_suffix(".rs").is_file() or (child / "mod.rs").is_file():
                target = module_file(me + p)
                if target:
                    found.add(target)
            continue
        target = module_file(parts)
        if target:
            found.add(target)
    return found


def blastn_files() -> list[Path]:
    start = sorted((SRC / "algorithm/blastn").rglob("*.rs"))
    seen, todo = set(start), list(start)
    while todo:
        for target in references(todo.pop()):
            if str(target.relative_to(SRC)).startswith(OTHER_PROGRAMS) or target in seen:
                continue
            seen.add(target)
            todo.append(target)
    return sorted(seen)


# ---------------------------------------------------------------- annotations
def rust_items(lines: list[str]) -> list[str]:
    """For each line, the next `fn` within 60 lines, else the previous `fn`."""
    fn = re.compile(r"\bfn\s+([A-Za-z_][A-Za-z0-9_]*)")
    names = [m.group(1) if (m := fn.search(line)) and not line.lstrip().startswith("//") else None for line in lines]
    result, previous = [], "<module>"
    for i in range(len(lines)):
        if names[i]:
            previous = names[i]
        upcoming = next((names[j] for j in range(i, min(i + 60, len(lines))) if names[j]), None)
        result.append(upcoming or previous)
    return result


def annotations(path: Path) -> list[dict]:
    lines = path.read_text(errors="replace").splitlines()
    items = rust_items(lines)
    out = []
    for i, line in enumerate(lines):
        if "NCBI reference" not in line:
            continue
        tokens = REF.findall(line.split("NCBI reference", 1)[1])
        j = i + 1
        while not tokens and j < min(i + 4, len(lines)) and lines[j].lstrip().startswith("//"):
            tokens = REF.findall(lines[j])
            j += 1
        out.append({"rust_file": str(path.relative_to(REPO)), "rust_line": i + 1, "rust_item": items[i],
                    "tokens": [(name, int(a) if a else None, int(b) if b else (int(a) if a else None))
                               for name, a, b in tokens]})
    return out


# ---------------------------------------------------------------- NCBI functions
class NcbiTree:
    def __init__(self, root: Path):
        self.root = root
        self.cxx = root / "c++"
        self.by_name: dict[str, list[Path]] = collections.defaultdict(list)
        for sub in ("src", "include"):
            for path in (self.cxx / sub).rglob("*"):
                if path.suffix in (".c", ".cpp", ".h", ".hpp", ".inl", ".prt") and path.is_file():
                    self.by_name[path.name].append(path)
        self.cache: dict[Path, list[tuple[int, int, str]]] = {}

    def resolve(self, token: str) -> Path | None:
        if "/" in token:
            tail = token.split("c++/", 1)[1] if "c++/" in token else token
            candidate = self.cxx / tail
            if candidate.is_file():
                return candidate
        matches = self.by_name.get(Path(token).name, [])
        if not matches:
            return None
        rank = lambda p: (0 if "algo/blast/core" in str(p) else 1 if "algo/blast" in str(p) else 2, str(p))
        if "/" in token:
            suffix = [p for p in matches if str(p).endswith(token.split("c++/", 1)[-1])]
            matches = suffix or matches
        return sorted(matches, key=rank)[0]

    def blocks(self, path: Path) -> list[tuple[int, int, str]]:
        if path not in self.cache:
            self.cache[path] = scan_blocks(path.read_text(errors="replace"))
        return self.cache[path]

    def enclosing(self, path: Path, first: int, last: int) -> tuple[str, int, int]:
        """The innermost block around the first line of the range that is in a block
        (a range often starts with the comment above a function)."""
        blocks = self.blocks(path)
        for line in range(first, max(first, last) + 1):
            best = None
            for start, end, name in blocks:
                if start <= line <= end and (best is None or start >= best[0]):
                    best = (start, end, name)
            if best:
                return best[2], best[0], best[1]
        return "<file scope>", 0, 0


def strip_code(text: str) -> str:
    """Blanks comments, literals and preprocessor lines, keeping line breaks."""
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i)); i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append(re.sub(r"[^\n]", " ", text[i:j])); i = j
        elif c in "\"'":
            j = i + 1
            while j < n and text[j] != c and text[j] != "\n":
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            out.append(re.sub(r"[^\n]", " ", text[i:j])); i = j
        elif c == "#" and (i == 0 or text[text.rfind("\n", 0, i) + 1:i].strip() == ""):
            j = i
            while True:
                k = text.find("\n", j)
                k = n if k < 0 else k
                if k > 0 and text[k - 1] == "\\" and k < n:
                    j = k + 1
                    continue
                break
            directive = text[i:k]
            define = re.match(r"#\s*define\s+([A-Za-z_][A-Za-z0-9_]*)", directive)
            out.append(f"\x01{define.group(1)}\x01" if define else "")
            out.append(re.sub(r"[^\n]", " ", directive)[len(f"\x01{define.group(1)}\x01") if define else 0:])
            i = k
        else:
            out.append(c); i += 1
    return "".join(out)


def scan_blocks(text: str) -> list[tuple[int, int, str]]:
    """(start line, end line, name) of function bodies, aggregates and #define lines."""
    code = strip_code(text)
    blocks = []
    line_of = []
    line = 1
    for ch in code:
        line_of.append(line)
        line += ch == "\n"
    for match in re.finditer(r"\x01([A-Za-z_][A-Za-z0-9_]*)\x01", code):
        start = line_of[match.start()]
        end = start
        raw_lines = text.splitlines()
        while end - 1 < len(raw_lines) and raw_lines[end - 1].rstrip().endswith("\\"):
            end += 1
        blocks.append((start, end, f"#define {match.group(1)}"))
    stack: list[tuple[str, str, int]] = []  # (kind, name, start line)
    header_start = 0
    i = 0
    while i < len(code):
        ch = code[i]
        if ch in ";":
            header_start = i + 1
        elif ch == "{":
            header = code[header_start:i]
            kind, name = classify(header, stack)
            name_line = line_of[header_start + len(header) - len(header.lstrip())] if header.strip() else line_of[i]
            stack.append((kind, name, name_line))
            header_start = i + 1
        elif ch == "}":
            if stack:
                kind, name, start = stack.pop()
                if kind in ("function", "aggregate") and name:
                    if not any(k == "function" for k, _, _ in stack):
                        blocks.append((start, line_of[i], name if kind == "function" else f"<{name}>"))
            header_start = i + 1
        i += 1
    return blocks


def classify(header: str, stack) -> tuple[str, str]:
    h = " ".join(header.replace("\x01", " ").split())
    if any(k == "function" for k, _, _ in stack):
        return "inner", ""
    if re.search(r"\bnamespace\b", h) or re.match(r'^extern\s*$', h) or h.endswith('extern'):
        return "transparent", ""
    if re.search(r'extern\s*$', h) or h == "":
        return "transparent", ""
    aggregate = re.search(r"\b(struct|class|union|enum)\s+([A-Za-z_][A-Za-z0-9_]*)?[^()]*$", h)
    if aggregate and "(" not in h.split(aggregate.group(0))[0][-1:]:
        name = aggregate.group(2) or "anonymous " + aggregate.group(1)
        if "=" in h and "(" not in h:
            return "inner", ""
        return "aggregate", name
    if "=" in h.split("(")[0] if "(" in h else "=" in h:
        return "inner", ""
    if "(" not in h:
        return "inner", ""
    depth, first = 0, None
    for k, c in enumerate(h):
        if c == "(":
            if depth == 0 and first is None:
                first = k
            depth += 1
        elif c == ")":
            depth -= 1
    before = h[:first].rstrip()
    m = re.search(r"((?:[A-Za-z_~][A-Za-z0-9_]*\s*(?:<[^<>]*>)?\s*::\s*)*(?:operator\s*\S+|~?[A-Za-z_][A-Za-z0-9_]*))$", before)
    if not m:
        return "inner", ""
    name = re.sub(r"\s+", "", m.group(1))
    if name in ("if", "for", "while", "switch", "catch", "return", "sizeof", "BEGIN_SCOPE", "BEGIN_NCBI_SCOPE"):
        return "inner", ""
    classes = [n for k, n, _ in stack if k == "aggregate" and n]
    if classes and "::" not in name:
        name = f"{classes[-1]}::{name}"
    return "function", name


# ---------------------------------------------------------------- main
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ncbi", type=Path, default=DEFAULT_NCBI)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--no-check", action="store_true")
    args = parser.parse_args()

    files = blastn_files()
    notes = [a for path in files for a in annotations(path)]
    counts = {"files": len(files), "annotations": len(notes)}
    print(f"files {counts['files']} (expected {EXPECTED['files']}), annotations {counts['annotations']} "
          f"(expected {EXPECTED['annotations']})", flush=True)
    if not args.no_check and any(counts[k] != EXPECTED[k] for k in counts):
        print("FIRST CHECK FAILED: the counts differ from the pinned stage-1 counts", file=sys.stderr)
        return 1

    tree = NcbiTree(args.ncbi)
    rows = []
    unresolved = collections.Counter()
    functions: dict[tuple[str, str], dict] = {}
    for note in notes:
        if not note["tokens"]:
            rows.append({**{k: note[k] for k in ("rust_file", "rust_line", "rust_item")},
                         "ncbi_token": "", "ncbi_file": "", "ncbi_lines": "", "ncbi_function": "<no file:line>"})
            continue
        for token, first, last in note["tokens"]:
            path = tree.resolve(token)
            if path is None:
                unresolved[token] += 1
                function, start, end, relative = "<unresolved file>", 0, 0, token
            else:
                relative = str(path.relative_to(args.ncbi))
                if first is None:
                    function, start, end = "<no line>", 0, 0
                else:
                    function, start, end = tree.enclosing(path, first, last)
            rows.append({"rust_file": note["rust_file"], "rust_line": note["rust_line"], "rust_item": note["rust_item"],
                         "ncbi_token": token, "ncbi_file": relative,
                         "ncbi_lines": "" if first is None else (f"{first}-{last}" if last != first else str(first)),
                         "ncbi_function": function})
            if path is not None and first is not None:
                key = (relative, function)
                entry = functions.setdefault(key, {"ncbi_file": relative, "ncbi_function": function,
                                                   "ncbi_start": start, "ncbi_end": end, "refs": 0, "rust": []})
                entry["refs"] += 1
                location = f"{note['rust_file'].removeprefix('LOSAT/src/')}:{note['rust_line']}({note['rust_item']})"
                if location not in entry["rust"]:
                    entry["rust"].append(location)

    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / "stage1_refs.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    real = [f for f in functions.values() if not f["ncbi_function"].startswith(("<", "#define"))]
    with open(args.out / "stage1_functions.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["ncbi_file", "ncbi_function", "ncbi_start", "ncbi_end", "refs", "losat_locations"])
        for f in sorted(functions.values(), key=lambda f: (f["ncbi_file"], f["ncbi_start"], f["ncbi_function"])):
            writer.writerow([f["ncbi_file"], f["ncbi_function"], f["ncbi_start"], f["ncbi_end"], f["refs"], "; ".join(f["rust"])])
    summary = {
        **counts,
        "reference_tokens": sum(1 for r in rows if r["ncbi_token"]),
        "annotations_without_file_line": sum(1 for r in rows if r["ncbi_function"] == "<no file:line>"),
        "ncbi_files": len({r["ncbi_file"] for r in rows if r["ncbi_file"] and not r["ncbi_function"].startswith("<unresolved")}),
        "ncbi_functions": len(real),
        "ncbi_non_function_targets": len(functions) - len(real),
        "unresolved_tokens": dict(unresolved.most_common()),
        "rust_files": [str(p.relative_to(REPO)) for p in files],
    }
    (args.out / "stage1_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "rust_files"}, indent=2))
    if not args.no_check and summary["ncbi_functions"] != EXPECTED["ncbi_functions"]:
        print(f"FIRST CHECK FAILED: {summary['ncbi_functions']} NCBI functions, expected {EXPECTED['ncbi_functions']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
