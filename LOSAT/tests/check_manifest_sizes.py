#!/usr/bin/env python3
"""Keep each file of `LOSAT/tests/db_search_manifest/` at 512 KiB or less.

No NCBI counterpart: this checks the size of LOSAT's own fixture lists, not ported behavior.

The `-db` fixture rows are one file per program (db_search_manifest/<program>.tsv), but the row
families moved to their own file (<program>.<family>.tsv: the rows `db.<family>.*` of that program,
`db_search_fixtures.family_files`). The rule "no committed file over about 1 MB without the Owner's
consent" would stop every line that adds rows once a file nears it, so a file that passes 512 KiB
fails here: move a row family of it to its own file before it grows further (the loader reads
every such file; no other change is needed).

Usage: python3 LOSAT/tests/check_manifest_sizes.py --root .
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

LIMIT = 512 * 1024
MANIFEST_DIR = Path("LOSAT") / "tests" / "db_search_manifest"


def oversized(directory: Path, limit: int = LIMIT) -> list[tuple[str, int]]:
    """The files of `directory` over `limit` bytes, with their sizes."""
    return [(path.name, path.stat().st_size) for path in sorted(directory.glob("*.tsv"))
            if path.stat().st_size > limit]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args(argv)
    directory = args.root.resolve() / MANIFEST_DIR
    if not directory.is_dir():
        print(f"manifest size check error: no directory {MANIFEST_DIR.as_posix()}", file=sys.stderr)
        return 2
    big = oversized(directory)
    for name, size in big:
        print(f"{MANIFEST_DIR.as_posix()}/{name}\t{size} bytes\tover {LIMIT}: move a row family to its own "
              f"file, <program>.<family>.tsv (db_search_fixtures.py reads it)")
    if big:
        return 1
    print(f"manifest size check: every file of {MANIFEST_DIR.as_posix()} is {LIMIT} bytes or less")
    return 0


if __name__ == "__main__":
    sys.exit(main())
