#!/usr/bin/env python3
"""Check that every NCBI quotation of AUTHORITY.md is on the cited lines.

A table row whose second cell is `path:lines` (relative to the NCBI `c++/` directory;
lines `a-b` or `a`, several spans separated by commas) and whose third cell holds one or
more back-quoted snippets (separated by ` / `) is checked: each snippet, with white space
collapsed, must occur in one of the cited lines (CR removed, white space collapsed).
Prints each failing row and exits 1 when any fails.

Usage: check_authority.py [AUTHORITY.md]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

NCBI = Path("/mnt/c/Users/genom/GitHub/ncbi-blast/c++")
HERE = Path(__file__).resolve().parent


def norm(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "AUTHORITY.md"
    checked = failed = 0
    cache: dict[str, list[str]] = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        cells = [cell.strip() for cell in line.split(" | ")]
        if len(cells) < 3:
            continue
        ref = re.fullmatch(r"`([^`:]+\.(?:c|cpp|h|hpp)):([\d,\- ]+)`", cells[1])
        if not ref:
            continue
        source = ref.group(1)
        if source not in cache:
            cache[source] = (NCBI / source).read_text(errors="replace").replace("\r", "").splitlines()
        text = cache[source]
        window = []
        for span in ref.group(2).replace(" ", "").split(","):
            lo, _, hi = span.partition("-")
            window += [norm(x) for x in text[int(lo) - 1:int(hi or lo)]]
        snippets = re.findall(r"`([^`]+)`", cells[2])
        checked += 1
        missing = [s for s in snippets if not any(norm(s) in w for w in window)]
        if not snippets or missing:
            failed += 1
            print(f"{path.name}:{number}: {source}:{ref.group(2)}: not found: {missing or 'no snippet'}")
    print(f"{checked} quotations checked, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
