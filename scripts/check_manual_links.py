#!/usr/bin/env python3
"""Checks LOSAT user documentation sources before publishing.

Manual sources (SRC_DIR, built by mdBook):
- Every page that SUMMARY.md lists exists, and every page in the tree is listed.
- Every relative link, image, and <source srcset> points at a file in the tree, and a
  `#fragment` names a heading of that page (mdBook heading ids).
- Text outside code says none of "verification", "verify", "verified", "certified",
  "certification": the Owner's rule for text users read (2026-11-04).

Repository Markdown files (--repo-readme PATH):
- Relative links, <img>, <source srcset>, and <a> targets stay inside --repo-root and
  exist (with valid `#fragment` heading anchors on `.md` targets).
- Text outside code contains none of the forbidden words above.

Package READMEs (--package-readme PATH, for crates.io / PyPI / npm):
- Same forbidden-word check outside code.
- Only https:// URLs or same-page `#fragment` heading anchors are allowed; relative file
  paths fail because package registries do not serve files outside the README.

Usage: check_manual_links.py SRC_DIR [--repo-root DIR] [--repo-readme PATH]... [--package-readme PATH]...
Exit status 1 with one line per problem.
"""

from __future__ import annotations

import argparse
import re
import sys
import unicodedata
from pathlib import Path

LINK = re.compile(r"!?\[(?:[^\]\[]|\[[^\]]*\])*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
HTML_SRC = re.compile(r"""<img\s[^>]*src=["']([^"']+)["']""", re.IGNORECASE)
HTML_SRCSET = re.compile(r"""<source\s[^>]*srcset=["']([^"']+)["']""", re.IGNORECASE)
HTML_HREF = re.compile(r"""<a\s[^>]*href=["']([^"']+)["']""", re.IGNORECASE)
HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
FENCE = re.compile(r"^\s*(```|~~~)")
INLINE_CODE = re.compile(r"`+[^`]*`+")
ATTR_URL = re.compile(r"""(?:src|srcset|href)=["'][^"']*["']""", re.IGNORECASE)
FORBIDDEN = re.compile(r"\b(verif(?:y|ied|ies|ication|ying)|certif(?:ied|ication|y))\b", re.IGNORECASE)


def heading_id(text: str) -> str:
    """mdBook's `normalize_id`: keep letters, digits, `_` and `-`, spaces become `-`, lower case."""
    text = re.sub(r"<[^>]+>", "", text)
    text = INLINE_CODE.sub(lambda m: m.group(0).strip("`"), text)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    out = []
    for ch in text.strip():
        if ch.isalnum() or ch in "_-":
            out.append(ch.lower())
        elif ch.isspace():
            out.append("-")
    return unicodedata.normalize("NFC", "".join(out))


def prose_lines(text: str):
    """(line number, line) of the lines outside fenced code blocks."""
    in_fence = False
    for number, line in enumerate(text.splitlines(), 1):
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if not in_fence:
            yield number, line


def anchors(path: Path) -> set[str]:
    ids: dict[str, int] = {}
    found = set()
    for _, line in prose_lines(path.read_text(encoding="utf-8")):
        match = HEADING.match(line)
        if not match:
            continue
        base = heading_id(match.group(2))
        count = ids.get(base, 0)
        ids[base] = count + 1
        found.add(base if count == 0 else f"{base}-{count}")
    return found


def _extract_srcset_urls(prose: str) -> list[str]:
    urls: list[str] = []
    for raw in HTML_SRCSET.findall(prose):
        for candidate in raw.split(","):
            part = candidate.strip().split()
            if part:
                urls.append(part[0])
    return urls


def _extract_targets(prose: str) -> list[str]:
    return (
        LINK.findall(prose)
        + HTML_SRC.findall(prose)
        + _extract_srcset_urls(prose)
        + HTML_HREF.findall(prose)
    )


def _check_forbidden(label: str, number: int, prose: str, problems: list[str]) -> None:
    cleaned = ATTR_URL.sub("", re.sub(r"\]\([^)]*\)", "]", prose))
    for word in FORBIDDEN.findall(cleaned):
        problems.append(f"{label}:{number}: says '{word}' (use 'check', 'Supported' or 'Not yet supported')")


def _resolve_anchor_check(
    label: str,
    number: int,
    target: str,
    dest: Path,
    anchor_cache: dict[Path, set[str]],
    problems: list[str],
) -> None:
    _, _, fragment = target.partition("#")
    if fragment and dest.suffix == ".md":
        if dest not in anchor_cache:
            anchor_cache[dest] = anchors(dest)
        if fragment not in anchor_cache[dest]:
            problems.append(f"{label}:{number}: {target} names no heading of {dest.name}")


def check(src: Path) -> list[str]:
    problems: list[str] = []
    pages = {p.relative_to(src).as_posix() for p in src.rglob("*.md")}
    summary = src / "SUMMARY.md"
    if not summary.is_file():
        return [f"{summary}: missing"]
    listed = set()
    for number, line in prose_lines(summary.read_text(encoding="utf-8")):
        for target in LINK.findall(line):
            listed.add(target)
            if target not in pages:
                problems.append(f"SUMMARY.md:{number}: lists {target}, which does not exist")
    for page in sorted(pages - listed - {"SUMMARY.md"}):
        problems.append(f"{page}: not listed in SUMMARY.md")

    anchor_cache: dict[Path, set[str]] = {}
    src_resolved = src.resolve()
    for page in sorted(pages - {"SUMMARY.md"}):
        path = src / page
        for number, line in prose_lines(path.read_text(encoding="utf-8")):
            prose = INLINE_CODE.sub("", line)
            _check_forbidden(page, number, prose, problems)
            for target in _extract_targets(prose):
                if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.IGNORECASE):
                    continue
                file_part, _, _ = target.partition("#")
                dest = (path.parent / file_part).resolve() if file_part else path.resolve()
                try:
                    dest.relative_to(src_resolved)
                except ValueError:
                    problems.append(f"{page}:{number}: {target} leaves the manual (use a full https:// link)")
                    continue
                if not dest.exists():
                    problems.append(f"{page}:{number}: {target} does not exist")
                    continue
                _resolve_anchor_check(page, number, target, dest, anchor_cache, problems)
    return problems


def check_repo_markdown(path: Path, repo_root: Path | None = None) -> list[str]:
    problems: list[str] = []
    if not path.is_file():
        return [f"{path}: missing"]
    root_resolved = (repo_root or path.parent).resolve()
    path_resolved = path.resolve()
    try:
        label = path_resolved.relative_to(root_resolved).as_posix()
    except ValueError:
        label = str(path)

    anchor_cache: dict[Path, set[str]] = {}
    for number, line in prose_lines(path_resolved.read_text(encoding="utf-8")):
        prose = INLINE_CODE.sub("", line)
        _check_forbidden(label, number, prose, problems)
        for target in _extract_targets(prose):
            if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.IGNORECASE):
                continue
            file_part, _, _ = target.partition("#")
            dest = (path_resolved.parent / file_part).resolve() if file_part else path_resolved
            try:
                dest.relative_to(root_resolved)
            except ValueError:
                problems.append(f"{label}:{number}: {target} leaves the repository root (use a full https:// link)")
                continue
            if not dest.exists():
                problems.append(f"{label}:{number}: {target} does not exist")
                continue
            _resolve_anchor_check(label, number, target, dest, anchor_cache, problems)
    return problems


def check_package_readme(path: Path, display_name: str | None = None) -> list[str]:
    problems: list[str] = []
    if not path.is_file():
        return [f"{path}: missing"]
    label = display_name or str(path)
    anchor_cache: dict[Path, set[str]] = {}
    path_resolved = path.resolve()
    for number, line in prose_lines(path_resolved.read_text(encoding="utf-8")):
        prose = INLINE_CODE.sub("", line)
        _check_forbidden(label, number, prose, problems)
        for target in _extract_targets(prose):
            if target.startswith("#"):
                _resolve_anchor_check(label, number, target, path_resolved, anchor_cache, problems)
                continue
            if not re.match(r"^https?://", target, re.IGNORECASE):
                problems.append(
                    f"{label}:{number}: relative or non-HTTP target '{target}' in package README (must use full https:// URL)"
                )
    return problems


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Checks LOSAT user manual, repository README/docs, and package READMEs.",
    )
    parser.add_argument("src_dir", nargs="?", type=Path, help="Manual source directory (docs/manual/src)")
    parser.add_argument("--repo-root", type=Path, default=None, help="Repository root for --repo-readme checks")
    parser.add_argument("--repo-readme", action="append", default=[], type=Path, help="Repo-scoped Markdown file to check")
    parser.add_argument("--package-readme", action="append", default=[], type=Path, help="Package README (no relative links)")
    return parser


def main(argv: list[str]) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv[1:])
    if args.src_dir is None and not args.repo_readme and not args.package_readme:
        parser.print_usage(sys.stderr)
        return 2

    problems: list[str] = []
    if args.src_dir is not None:
        problems.extend(check(args.src_dir))
    repo_root = args.repo_root or Path.cwd()
    for readme in args.repo_readme:
        problems.extend(check_repo_markdown(readme, repo_root=repo_root))
    for pkg in args.package_readme:
        problems.extend(check_package_readme(pkg))

    for problem in problems:
        print(problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
