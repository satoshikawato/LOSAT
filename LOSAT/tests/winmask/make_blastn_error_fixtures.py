#!/usr/bin/env python3
"""Freeze tests/fixtures/winmask_blastn/: blastn -window_masker_db with counts files that cannot be read.

NCBI BLAST+ is the comparison oracle only. The exception reports of these runs name the source tree of
the NCBI build, so the frozen bytes come from the official Linux x64 release of BLAST+ 2.17.0
(`ncbi-blast-2.17.0+-x64-linux.tar.gz`), the tree LOSAT names (LOSAT/src/utils/winmask/ncbi_report.rs,
line NP M/DECISIONS.md 14); other builds differ in that path only.

- make_counts(): writes the counts files of the cases to LOSAT/target/winmask_blastn/ from
  tests/winmask/counts/mg1655_u9.obinary (tests/blastn_window_masker.rs makes the same bytes).
- freeze --ncbi-bin DIR: runs `DIR/blastn <args> -num_threads 1` from LOSAT/ for each row of
  cases.tsv with the row's environment and writes <case>.out, <case>.err and <case>.rc.

Usage: make_blastn_error_fixtures.py freeze --ncbi-bin /path/to/ncbi-blast-2.17.0+/bin
"""
from __future__ import annotations

import argparse
import os
import struct
import subprocess
import sys
from pathlib import Path

ENGINE = Path(__file__).resolve().parents[2]
FIXTURES = ENGINE / "tests/fixtures/winmask_blastn"
COUNTS = ENGINE / "tests/winmask/counts/mg1655_u9.obinary"
WORK = ENGINE / "target/winmask_blastn"
# NCBI reads these and each one changes the batches or the report (tests/blastn_regression_fixtures.py).
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE", "ADAPTIVE_CBS",
              "OVERLAP_CHUNK_SIZE", "PRE_FETCH_SEQS_LIMIT")


def words(*values: int) -> bytes:
    return b"".join(struct.pack("<I", value) for value in values)


def make_counts() -> None:
    """The counts files of the cases. The optimized binary layout: the metadata block (3, its
    length, its lines), then the format word, unit size, M, k, roff, bc, the four thresholds, the
    flag word of the cache bit array, the array (4^unit / 32 words), the hash table (2^k words) and
    the values table (M half-words)."""
    WORK.mkdir(parents=True, exist_ok=True)
    data = COUNTS.read_bytes()
    md = struct.unpack_from("<I", data, 4)[0]
    head, body = data[:8 + md], data[8 + md:]
    _fmt, unit, m, k, _roff, bc = struct.unpack_from("<6I", body, 0)
    ht_at = 40 + 4 + 4 * ((1 << (2 * unit)) // 32)
    bad = bytearray(body)
    for index in range(1 << k):
        # three collisions from index M: past the values table (`bad index`)
        struct.pack_into("<I", bad, ht_at + 4 * index, (m << bc) | 3)
    pct_md = b"##windowmasker-statistics-format-version:obinary 2.0.0\0##pct:\0"
    files = {
        "unit17.obinary": head + body[:4] + words(17) + body[8:],
        "short_cba.obinary": head + body[:48],
        "bad_index.obinary": head + bytes(bad),
        "pct6.obinary": words(3, len(pct_md)) + pct_md + body,
        "binary_unit0.binary": words(0, 0) + words(5, 1) + words(1, 2, 3, 4),
        "ascii.counts": b"##windowmasker-statistics-format-version:ascii 1.1.0\n9\nAAAAAAAAA 5\n",
        "empty.obinary": b"",
    }
    for name, content in files.items():
        (WORK / name).write_bytes(content)
    unreadable = WORK / "unreadable.obinary"
    if unreadable.exists():
        unreadable.chmod(0o600)
    unreadable.write_bytes(data)
    unreadable.chmod(0)


def read_cases() -> list[tuple[str, str, list[str]]]:
    rows = []
    for line in (FIXTURES / "cases.tsv").read_text().splitlines():
        if line.startswith("#") or not line:
            continue
        case, env, args = line.split("\t")
        rows.append((case, env, ["" if token == "''" else token.replace("%20", " ") for token in args.split(" ")]))
    return rows


def command_freeze(args) -> int:
    if any(key in os.environ for key in REPORT_ENV) or (Path.home() / ".ncbirc").exists():
        raise SystemExit(f"unset {REPORT_ENV} and remove ~/.ncbirc before freezing")
    make_counts()
    blastn = str(Path(args.ncbi_bin).resolve() / "blastn")
    for case, env_text, argv in read_cases():
        env = {key: value for key, value in os.environ.items() if not key.startswith(("LOSAT_", "RAYON_"))}
        for item in env_text.split():
            key, _, value = item.partition("=")
            env[key] = value
        result = subprocess.run([blastn, *argv, "-num_threads", "1"], cwd=ENGINE, capture_output=True, env=env)
        (FIXTURES / f"{case}.out").write_bytes(result.stdout)
        (FIXTURES / f"{case}.err").write_bytes(result.stderr)
        (FIXTURES / f"{case}.rc").write_text(f"{result.returncode}\n")
        print(f"{case}\trc {result.returncode}\t{len(result.stdout)} stdout bytes\t{len(result.stderr)} stderr bytes")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--ncbi-bin", required=True)
    args = parser.parse_args()
    return command_freeze(args)


if __name__ == "__main__":
    sys.exit(main())
