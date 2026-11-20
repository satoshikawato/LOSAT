#!/usr/bin/env python3
"""Keep the source files that only unit tests read out of the published `losat` crate.

No NCBI counterpart: this checks LOSAT's own packaging, not ported behavior.

crates.io refuses a `.crate` file over 10 MiB (the "Crate size" step of ci.yml). A file that only
`#[cfg(test)]` code reaches (`#[cfg(test)] mod x_tests;`, a `tests.rs`, a table that a test
`include!`s) is never read by `cargo build`, by the build with which `cargo publish` verifies the
package, or by docs.rs, so packaging it only adds size: on 2026-11-20 the 118 such files were 0.58
MiB of the 7.84 MiB crate. `exclude` in LOSAT/Cargo.toml leaves them out by their names
(`*_tests.rs`, `tests.rs`, `tests_*.rs`, a `*_tests/` directory) and names the few others.

The check follows, from the crate roots (`src/lib.rs`, the `[[bin]]` paths and `src/main.rs`,
`build.rs`, `examples/`), the `mod x;` declarations (with `#[path]`) and the `include!`,
`include_str!` and `include_bytes!` of a string path (also
`concat!(env!("CARGO_MANIFEST_DIR"), "...")`). A file that one of them reaches outside an item under
`#[cfg(test)]` or `#[cfg(all(test, ...))]` is read by the build; a file reached only from such items
is a test file. It fails, naming the file and the `exclude` pattern to change,
- when the build reads a file that is not packaged (an `exclude` pattern matches it, or it lies
  outside the crate): the published crate would not build;
- when a `.rs` file that only tests reach is packaged: name it as above, or add its path to
  `exclude`;
- when a packaged `.rs` file is reached by nothing that the check reads: delete a dead file, or
  teach the check the form that declares it.

The packaged files are computed as Cargo lists them, without cargo: `git ls-files` of LOSAT/
(tracked files, and untracked ones that .gitignore does not ignore), less a nested package and
the root `target`, less what `exclude` matches with the gitignore rules of the `ignore` crate that
Cargo uses (`matched_path_or_any_parents`: the file, then each parent directory; the last
matching pattern decides). `cargo package --list --offline` is not used: it needs the registry
index, which a runner with an empty cache lacks (2026-11-20: "no matching package named `bzip2`"
offline). `--list` prints the computed list, to compare with `cargo package --list`.

Usage: python3 LOSAT/tests/check_crate_test_files.py --root . [--list]
"""

from __future__ import annotations

import argparse
import posixpath
import re
import subprocess
import sys
import tomllib
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

CRATE = "LOSAT"
FIX_TEST = ("name it <name>_tests.rs, tests.rs or tests_<name>.rs, or put it in a <name>_tests/ "
            "directory, or add its path to `exclude` in LOSAT/Cargo.toml")


# ---------------------------------------------------------------- Cargo's `exclude`
@dataclass(frozen=True)
class Pattern:
    text: str
    negated: bool
    dir_only: bool
    regex: re.Pattern[str]


def glob_regex(glob: str) -> str:
    """The regex of a gitignore glob in which `*` and `?` do not match `/` (the `ignore` crate)."""
    parts = glob.split("/")
    out = []
    for i, part in enumerate(parts):
        last = i == len(parts) - 1
        if part == "**":
            if last:
                out.append(".*")
            else:
                out.append("(?:.*/)?")
            continue
        j = 0
        piece = ""
        while j < len(part):
            c = part[j]
            if c == "*":
                while j + 1 < len(part) and part[j + 1] == "*":
                    j += 1
                piece += "[^/]*"
            elif c == "?":
                piece += "[^/]"
            elif c == "[":
                start = j + 1 + (part[j + 1:j + 2] in ("!", "^"))
                close = part.find("]", start + (part[start:start + 1] == "]"))
                if close < 0:
                    piece += re.escape(c)
                else:
                    body = part[j + 1:close]
                    if body[:1] in ("!", "^"):
                        body = "^" + body[1:]
                    piece += "[" + body.replace("\\", "\\\\") + "]"
                    j = close
            elif c == "\\" and j + 1 < len(part):
                j += 1
                piece += re.escape(part[j])
            else:
                piece += re.escape(c)
            j += 1
        out.append(piece + ("" if last else "/"))
    return "".join(out)


def compile_pattern(line: str) -> Pattern:
    """One `exclude` entry, read as a .gitignore line is (`!`, a leading or inner `/`, a trailing `/`)."""
    text = line
    negated = line.startswith("!")
    if negated:
        line = line[1:]
    elif line.startswith("\\!") or line.startswith("\\#"):
        line = line[1:]
    dir_only = line.endswith("/")
    line = line.rstrip("/")
    anchored = line.startswith("/")
    line = line.lstrip("/")
    if not anchored and "/" not in line and not line.startswith("**/"):
        line = "**/" + line
    if line.endswith("/**"):
        line += "/*"
    return Pattern(text, negated, dir_only, re.compile(glob_regex(line)))


class Exclude:
    """The `exclude` list of a Cargo.toml, as Cargo applies it to a path under the package."""

    def __init__(self, lines: Iterable[str]) -> None:
        self.patterns = [compile_pattern(line) for line in lines]

    def _match(self, path: str, is_dir: bool) -> Pattern | None:
        for pattern in reversed(self.patterns):
            if pattern.dir_only and not is_dir:
                continue
            if pattern.regex.fullmatch(path):
                return pattern
        return None

    @lru_cache(maxsize=None)
    def _dir(self, path: str) -> Pattern | None:
        hit = self._match(path, True)
        if hit is not None:
            return hit
        parent = posixpath.dirname(path)
        return self._dir(parent) if parent else None

    def excluding(self, path: str) -> Pattern | None:
        """The pattern that leaves `path` (a file) out of the package, or None when it is packaged."""
        hit = self._match(path, False)
        if hit is None:
            parent = posixpath.dirname(path)
            hit = self._dir(parent) if parent else None
        return None if hit is None or hit.negated else hit


def git_files(root: Path) -> list[str]:
    """The files under LOSAT/ that Git tracks or would add (not ignored), relative to LOSAT/."""
    def ls_files(*options: str) -> set[str]:
        out = subprocess.run(["git", "-C", str(root), "ls-files", "-z", *options, "--", CRATE],
                             capture_output=True, check=True).stdout.decode("utf-8")
        return {path[len(CRATE) + 1:] for path in out.split("\0") if path.startswith(CRATE + "/")}

    return sorted(ls_files("--cached", "--others", "--exclude-standard") - ls_files("--deleted"))


def package_files(files: Sequence[str], exclude: Exclude) -> tuple[list[str], dict[str, Pattern]]:
    """Cargo's package list: no nested package, no root `target`, nothing `exclude` matches."""
    nested = [posixpath.dirname(f) for f in files if f.endswith("/Cargo.toml")]
    packaged = []
    excluded = {}
    for f in files:
        if f.startswith("target/") or any(f.startswith(d + "/") for d in nested):
            continue
        pattern = exclude.excluding(f)
        if pattern is None:
            packaged.append(f)
        else:
            excluded[f] = pattern
    return packaged, excluded


# ---------------------------------------------------------------- Rust sources
TOKEN_RE = re.compile(r'//|/\*|(?<![A-Za-z0-9_])(?:br|cr|r)(#*)"|"|\'')
STRING_BODY_RE = re.compile(r'(?:[^"\\]|\\.)*"', re.S)
CHAR_RE = re.compile(r"'(?:[^\\'\n]|\\(?:x[0-9a-fA-F]{2}|u\{[0-9a-fA-F_]{1,6}\}|[nrt\\0'\"]))'")
CFG_RE = re.compile(r"#(!?)\[\s*cfg\s*\(")
MOD_RE = re.compile(r"(?<![A-Za-z0-9_])mod\s+(?:r#)?([A-Za-z_][A-Za-z0-9_]*)\s*([;{])")
PATH_ATTR_RE = re.compile(r'#\[\s*path\s*=\s*"([^"]*)"\s*\]')
INCLUDE_RE = re.compile(r"(?<![A-Za-z0-9_$])(include|include_str|include_bytes)!\s*\(\s*")
STRING_ARG_RE = re.compile(r'"([^"]*)"\s*,?\s*\)')
MANIFEST_ARG_RE = re.compile(r'concat!\s*\(\s*env!\s*\(\s*"([^"]*)"\s*\)\s*,\s*"([^"]*)"\s*,?\s*\)\s*,?\s*\)')


def blank(text: str) -> str:
    return re.sub(r"[^\n]", " ", text) if "\n" in text else " " * len(text)


def mask(text: str) -> str:
    """`text` with comments and the insides of string literals blanked (newlines kept)."""
    out = []
    pos = 0
    n = len(text)
    while True:
        m = TOKEN_RE.search(text, pos)
        if not m:
            break
        s = m.start()
        token = m.group(0)
        if token == "//":
            e = text.find("\n", s)
            e = n if e < 0 else e
            out.append(text[pos:s] + blank(text[s:e]))
        elif token == "/*":
            depth, i = 1, s + 2
            while depth and i < n:
                inner, end = text.find("/*", i), text.find("*/", i)
                if end < 0:
                    i = n
                    break
                if 0 <= inner < end:
                    depth, i = depth + 1, inner + 2
                else:
                    depth, i = depth - 1, end + 2
            e = i
            out.append(text[pos:s] + blank(text[s:e]))
        elif token == "'":
            c = CHAR_RE.match(text, s)
            e = c.end() if c else s + 1
            out.append(text[pos:e])
        else:
            if token == '"':
                body = STRING_BODY_RE.match(text, s + 1)
                start, end = s + 1, (body.end() - 1 if body else n)
                e = body.end() if body else n
            else:  # a raw string, `r#"..."#`
                start = m.end()
                end = text.find('"' + m.group(1), start)
                end, e = (n, n) if end < 0 else (end, end + 1 + len(m.group(1)))
            out.append(text[pos:start] + blank(text[start:end]) + text[end:e])
        pos = e
    out.append(text[pos:])
    return "".join(out)


def closing(masked: str, i: int, open_ch: str, close_ch: str) -> int:
    """The index of the bracket that closes the one at `i`."""
    depth = 0
    for j in range(i, len(masked)):
        if masked[j] == open_ch:
            depth += 1
        elif masked[j] == close_ch:
            depth -= 1
            if depth == 0:
                return j
    return len(masked) - 1


def item_end(masked: str, i: int) -> int:
    """The end of the item at `i`: its `;` or the `}` of its block, or the end of the enclosing one."""
    depth = 0
    braced = False
    for j in range(i, len(masked)):
        c = masked[j]
        if c in "([{":
            braced = braced or (c == "{" and depth == 0)
            depth += 1
        elif c in ")]}":
            depth -= 1
            if depth < 0:
                return j
            if depth == 0 and c == "}" and braced:
                return j + 1
        elif c == ";" and depth == 0:
            return j + 1
    return len(masked)


def after_attributes(masked: str, i: int) -> int:
    while True:
        while i < len(masked) and masked[i].isspace():
            i += 1
        if not masked.startswith("#[", i):
            return i
        i = closing(masked, i + 1, "[", "]") + 1


def split_args(text: str) -> list[str]:
    args, depth, current = [], 0, ""
    for c in text:
        depth += (c == "(") - (c == ")")
        if c == "," and depth == 0:
            args.append(current)
            current = ""
        else:
            current += c
    return args + [current]


def test_only(predicate: str) -> bool:
    """Whether a `cfg(...)` predicate holds only in a test build: `test`, or `all(...)` with such an argument."""
    predicate = re.sub(r"\s+", "", predicate)
    if predicate == "test":
        return True
    if predicate.startswith("all(") and predicate.endswith(")"):
        return any(test_only(arg) for arg in split_args(predicate[4:-1]))
    return False


@dataclass(frozen=True)
class Edge:
    kind: str  # "mod" or "include"
    name: str  # module name, or the include's path
    path_attr: str | None  # the module's `#[path]`
    inline: tuple[str, ...]  # the inline modules around a `mod` declaration
    test: bool
    line: int
    manifest: bool = False  # an include path relative to CARGO_MANIFEST_DIR


def edges(text: str) -> tuple[list[Edge], int]:
    """The `mod x;` declarations and string-path includes of a source file, and its unread includes."""
    masked = mask(text)
    inline = []  # (open, close, name) of every `mod name { ... }`
    declared = []
    for m in MOD_RE.finditer(masked):
        if m.group(2) == "{":
            brace = m.end() - 1
            inline.append((brace, closing(masked, brace, "{", "}"), m.group(1)))
        else:
            declared.append(m)
    gated = []
    for m in CFG_RE.finditer(masked):
        paren = m.end() - 1
        if not test_only(masked[paren + 1:closing(masked, paren, "(", ")")]):
            continue
        if m.group(1):  # `#![cfg(test)]` gates the enclosing inline module, or the file
            around = [(o, c) for o, c, _ in inline if o < m.start() < c]
            gated.append(max(around) if around else (0, len(masked)))
        else:
            end = closing(masked, m.start() + 1, "[", "]")
            gated.append((m.start(), item_end(masked, after_attributes(masked, end + 1))))

    def is_test(pos: int) -> bool:
        return any(s <= pos < e for s, e in gated)

    def line(pos: int) -> int:
        return masked.count("\n", 0, pos) + 1

    out = []
    for m in declared:
        s = m.start()
        window = max(masked.rfind(c, 0, s) for c in ";{}") + 1
        attrs = list(PATH_ATTR_RE.finditer(masked, window, s))
        path_attr = None
        if attrs:  # the masked string is blank: read it from the source
            path_attr = text[attrs[-1].start(1):attrs[-1].end(1)]
        around = tuple(name for o, c, name in sorted(inline) if o < s < c)
        out.append(Edge("mod", m.group(1), path_attr, around, is_test(s), line(s)))
    unread = 0
    for m in INCLUDE_RE.finditer(masked):
        arg = STRING_ARG_RE.match(masked, m.end())
        if arg:
            target = text[arg.start(1):arg.end(1)]
            out.append(Edge("include", target, None, (), is_test(m.start()), line(m.start())))
            continue
        arg = MANIFEST_ARG_RE.match(masked, m.end())
        if arg and text[arg.start(1):arg.end(1)] == "CARGO_MANIFEST_DIR":
            target = text[arg.start(2):arg.end(2)].lstrip("/")
            out.append(Edge("include", target, None, (), is_test(m.start()), line(m.start()), True))
            continue
        unread += 1
    return out, unread


# ---------------------------------------------------------------- the module tree
def crate_roots(manifest: dict, files: set[str]) -> list[str]:
    """The files Cargo compiles as crate roots for the library, binaries, build script and examples."""
    package = manifest.get("package", {})
    roots = [manifest.get("lib", {}).get("path", "src/lib.rs")]
    roots += [b["path"] for b in manifest.get("bin", []) if "path" in b]
    roots += [e["path"] for e in manifest.get("example", []) if "path" in e]
    build = package.get("build", "build.rs")
    if isinstance(build, str):
        roots.append(build)
    if package.get("autobins", True):
        roots.append("src/main.rs")
        roots += [f for f in files if re.fullmatch(r"src/bin/(?:[^/]+\.rs|[^/]+/main\.rs)", f)]
    if package.get("autoexamples", True):
        roots += [f for f in files if re.fullmatch(r"examples/(?:[^/]+\.rs|[^/]+/main\.rs)", f)]
    return sorted({posixpath.normpath(r) for r in roots if posixpath.normpath(r) in files})


def module_candidates(path: str, mod_rs: bool, edge: Edge) -> list[str]:
    """Where rustc looks for the file of `mod name;` declared in `path` (the Rust reference, "Modules")."""
    directory = posixpath.dirname(path)
    owner = directory if mod_rs else posixpath.join(directory, posixpath.basename(path)[:-3])
    if edge.path_attr is not None:
        base = posixpath.join(owner, *edge.inline) if edge.inline else directory
        return [posixpath.normpath(posixpath.join(base, edge.path_attr))]
    base = posixpath.join(owner, *edge.inline)
    return [posixpath.join(base, edge.name + ".rs"), posixpath.join(base, edge.name, "mod.rs")]


@dataclass
class Reach:
    build: bool  # reached outside `#[cfg(test)]`
    via: str  # the first declaration that reached it, `file:line`


def walk(crate: Path, files: set[str], roots: Sequence[str]) -> tuple[dict[str, Reach], list[str], int]:
    """Every file reached from the roots, whether the build reads it, the unresolved declarations,
    and the number of includes whose path is not a string the check reads."""
    reach: dict[str, Reach] = {}
    mod_rs: dict[str, bool] = {}
    parsed: dict[str, tuple[list[Edge], int]] = {}
    unresolved = []
    stack = []
    for root in roots:
        reach[root] = Reach(True, "crate root")
        mod_rs[root] = True
        stack.append(root)
    while stack:
        path = stack.pop()
        if not path.endswith(".rs") or path not in files:
            continue
        if path not in parsed:
            parsed[path] = edges((crate / path).read_text(encoding="utf-8", errors="replace"))
        for edge in parsed[path][0]:
            build = reach[path].build and not edge.test
            via = f"{path}:{edge.line}"
            if edge.kind == "mod":
                candidates = module_candidates(path, mod_rs[path], edge)
                target = next((c for c in candidates if c in files), None)
                if target is None:
                    unresolved.append(f"{via}: mod {edge.name} (looked for {', '.join(candidates)})")
                    continue
                child_mod_rs = edge.path_attr is not None or posixpath.basename(target) == "mod.rs"
            else:
                base = "" if edge.manifest else posixpath.dirname(path)
                target = posixpath.normpath(posixpath.join(base, edge.name))
                child_mod_rs = False
            seen = reach.get(target)
            if seen is not None and (seen.build or not build):
                continue
            reach[target] = Reach(build, via)
            if edge.kind == "mod":  # an included file is read as text or as part of its includer
                mod_rs[target] = child_mod_rs
                stack.append(target)
    return reach, unresolved, sum(unread for _, unread in parsed.values())


# ---------------------------------------------------------------- the check
def check(root: Path) -> tuple[list[str], str, list[str]]:
    """The findings, a one-line summary, and the packaged files."""
    crate = root / CRATE
    manifest = tomllib.loads((crate / "Cargo.toml").read_text(encoding="utf-8"))
    package = manifest.get("package", {})
    if "include" in package:
        raise ValueError("LOSAT/Cargo.toml has `include`: this check reads only `exclude`")
    exclude = Exclude(package.get("exclude", []))
    files = git_files(root)
    packaged, excluded = package_files(files, exclude)
    in_package = set(packaged)
    roots = crate_roots(manifest, set(files))
    if not roots:
        raise ValueError("LOSAT/Cargo.toml: no crate root found")
    reach, unresolved, unread = walk(crate, set(files), roots)
    findings = [f"LOSAT/{u}: the check finds no such file; fix the declaration or teach the check its form"
                for u in sorted(set(unresolved))]
    for path, r in sorted(reach.items()):
        if not r.build or path in in_package:
            continue
        if path in excluded:
            findings.append(f"LOSAT/{path}\tthe build reads it ({r.via}) but `exclude` pattern "
                            f"{excluded[path].text!r} leaves it out: the published crate would not build; "
                            f"rename the file or narrow the pattern in LOSAT/Cargo.toml")
        else:
            findings.append(f"LOSAT/{path}\tthe build reads it ({r.via}) but it is not in the package "
                            f"(outside the crate, or not in Git)")
    for path in packaged:
        if not path.endswith(".rs"):
            continue
        r = reach.get(path)
        if r is None:
            findings.append(f"LOSAT/{path}\tpackaged, but no `mod` or `include!` that the check reads "
                            f"reaches it: delete a dead file, or teach the check the form that declares it")
        elif not r.build:
            findings.append(f"LOSAT/{path}\tonly tests read it ({r.via}) but it is packaged: {FIX_TEST}")
    test_files = sorted(p for p, r in reach.items() if not r.build and p.endswith(".rs") and p in excluded)
    by_pattern = Counter(excluded[p].text for p in test_files)
    summary = (f"crate test-file check: {len(packaged)} files packaged; the build reads "
               f"{sum(1 for r in reach.values() if r.build)}, all packaged; the {len(test_files)} .rs files "
               f"that only tests read are left out ("
               + ", ".join(f"{pattern} {count}" for pattern, count in by_pattern.most_common())
               + f"); includes of a computed path, not followed: {unread}")
    return findings, summary, packaged


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path("."), help="the repository root (holds LOSAT/)")
    parser.add_argument("--list", action="store_true", help="print the packaged files and exit")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    try:
        findings, summary, packaged = check(root)
    except (OSError, ValueError, subprocess.CalledProcessError, tomllib.TOMLDecodeError) as error:
        print(f"crate test-file check error: {error}", file=sys.stderr)
        return 2
    if args.list:
        print("\n".join(packaged))
        return 0
    for finding in findings:
        print(finding)
    if findings:
        print(f"crate test-file check: {len(findings)} findings", file=sys.stderr)
        return 1
    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
