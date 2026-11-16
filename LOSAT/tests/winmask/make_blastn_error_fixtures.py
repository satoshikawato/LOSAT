#!/usr/bin/env python3
"""Freeze tests/fixtures/winmask_blastn/: blastn -window_masker_db with counts files that cannot be read,
and -window_masker_taxid with the directories of counts files made at run time.

NCBI BLAST+ is the comparison oracle only. The exception reports of these runs name the source tree of
the NCBI build, so the frozen bytes come from the official Linux x64 release of BLAST+ 2.17.0
(`ncbi-blast-2.17.0+-x64-linux.tar.gz`), the tree LOSAT names (LOSAT/src/utils/winmask/ncbi_report.rs,
line NP M/DECISIONS.md 14); other builds differ in that path only.

- make_counts(): writes the counts files of the cases to LOSAT/target/winmask_blastn/ from
  tests/winmask/counts/mg1655_u9.obinary, and the -window_masker_taxid directory taxid/ (562 and -1:
  that file, 7: the unit size 17, 8: unreadable, 9606: a directory in its place, 10: ascii counts, 11-59:
  the ascii and optimized ascii counts of tests/winmask/text_counts.tsv)
  with the registry files home_rc/.ncbirc, cfg_inidont/{blastn.ini,.ncbirc} and the malformed
  .ncbirc files of cfg_rcbad, cfg_rcbadapp, cfg_rcpre and cfg_rcbom, the files with an entry set again
  of cfg_rcdup and cfg_rcdupbad, and short.fa, subjects with fewer letters than subjects
  (tests/blastn_window_masker.rs makes the same bytes, with ncbi.ini for .ncbirc on Windows).
- freeze --ncbi-bin DIR [--cases NAME ...]: runs `DIR/blastn <args> -num_threads 1` from LOSAT/ for
  each row of cases.tsv (or the rows named) with the row's environment and writes <case>.out,
  <case>.err and <case>.rc.

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
TEXT_COUNTS = ENGINE / "tests/winmask/text_counts.tsv"
WORK = ENGINE / "target/winmask_blastn"
# NCBI reads these and each one changes the batches or the report (tests/blastn_regression_fixtures.py).
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE", "ADAPTIVE_CBS",
              "OVERLAP_CHUNK_SIZE", "PRE_FETCH_SEQS_LIMIT")


def words(*values: int) -> bytes:
    return b"".join(struct.pack("<I", value) for value in values)


def unescape(text: str) -> bytes:
    """The escapes of text_counts.tsv: \\n, \\r, \\t and \\\\."""
    out = bytearray()
    i = 0
    while i < len(text):
        if text[i] == "\\":
            out += {"n": b"\n", "r": b"\r", "t": b"\t", "\\": b"\\"}[text[i + 1]]
            i += 2
        else:
            out += text[i].encode()
            i += 1
    return bytes(out)


def text_counts() -> dict[str, bytes]:
    """The files of tests/winmask/text_counts.tsv, by their path under WORK."""
    files = {}
    for row in TEXT_COUNTS.read_text().splitlines():
        if not row or row.startswith("#"):
            continue
        taxid, base, *edits = row.split("\t")
        lines = [b""] if base == "-" else (COUNTS.parent / base).read_bytes().split(b"\n")
        cut, crlf = None, False
        for edit in edits:
            if edit == "crlf":
                crlf = True
            elif edit.startswith("cut="):
                cut = int(edit[4:])
            else:
                number, _, text = edit.partition("=")
                lines[int(number) - 1] = unescape(text)
        if cut is not None:
            lines = lines[:cut] + [b""]
        files[f"taxid/{taxid}/wmasker.obinary"] = (b"\r\n" if crlf else b"\n").join(lines)
    return files


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
    # the same with the entry of key 0 kept: a query without A, C, G or T only reads that one
    bad_after0 = bytearray(bad)
    bad_after0[ht_at:ht_at + 4] = body[ht_at:ht_at + 4]
    pct_md = b"##windowmasker-statistics-format-version:obinary 2.0.0\0##pct:\0"
    files = {
        "unit17.obinary": head + body[:4] + words(17) + body[8:],
        "short_cba.obinary": head + body[:48],
        "bad_index.obinary": head + bytes(bad),
        "bad_index_after0.obinary": head + bytes(bad_after0),
        "pct6.obinary": words(3, len(pct_md)) + pct_md + body,
        "binary_unit0.binary": words(0, 0) + words(5, 1) + words(1, 2, 3, 4),
        "ascii.counts": b"##windowmasker-statistics-format-version:ascii 1.1.0\n9\nAAAAAAAAA 5\n",
        "empty.obinary": b"",
    }
    for name, content in files.items():
        (WORK / name).write_bytes(content)
    unreadable = WORK / "unreadable.obinary"
    if not unreadable.exists():
        unreadable.write_bytes(data)
        unreadable.chmod(0)
    # -window_masker_taxid: <dir>/<taxid>/wmasker.obinary, the directory named by the registry files
    # relative to LOSAT/ (NCBI joins the strings and opens the path from the current directory).
    rc = b"[WINDOW_MASKER]\nWINDOW_MASKER_PATH = target/winmask_blastn/taxid\n"
    taxid_files = {
        "taxid/562/wmasker.obinary": data,
        "taxid/-1/wmasker.obinary": data,
        "taxid/7/wmasker.obinary": files["unit17.obinary"],
        "taxid/10/wmasker.obinary": files["ascii.counts"],
        "home_rc/.ncbirc": rc,
        "cfg_inidont/blastn.ini": b"[NCBI]\nDONT_USE_NCBIRC = 1\n",
        "cfg_inidont/.ncbirc": rc,
        # malformed .ncbirc files (line NP M3c): the entry above the line NCBI reports, the line
        # before the entry, a file of one byte whose stream fails
        "cfg_rcbad/blastn.ini": b"[NCBI]\nDONT_USE_NCBIRC = 1\n",
        "cfg_rcbad/.ncbirc": rc + b"[my blast]\nx = 1\n",
        "cfg_rcbadapp/.ncbirc": rc + b"[my blast]\nx = 1\n",
        "cfg_rcpre/blastn.ini": b"[NCBI]\nDONT_USE_NCBIRC = 1\n",
        "cfg_rcpre/.ncbirc": b"x y = 1\n" + rc,
        "cfg_rcbom/blastn.ini": b"[NCBI]\nDONT_USE_NCBIRC = 1\n",
        "cfg_rcbom/.ncbirc": b"\xef",
        # an entry set again (line NP M3c review): in a well-formed file, and above a malformed line
        "cfg_rcdup/blastn.ini": b"[NCBI]\nDONT_USE_NCBIRC = 1\n",
        "cfg_rcdup/.ncbirc": rc + rc.split(b"\n", 1)[1],
        "cfg_rcdupbad/blastn.ini": b"[NCBI]\nDONT_USE_NCBIRC = 1\n",
        "cfg_rcdupbad/.ncbirc": rc + rc.split(b"\n", 1)[1] + b"[my blast]\nx = 1\n",
        # subjects with fewer letters than subjects (the average subject length is 0)
        "short.fa": b">s1\nA\n>s2\n\n",
        # the ascii and optimized ascii counts (line NP M3d)
        **text_counts(),
    }
    for name, content in taxid_files.items():
        (WORK / name).parent.mkdir(parents=True, exist_ok=True)
        (WORK / name).write_bytes(content)
    (WORK / "taxid/9606/wmasker.obinary").mkdir(parents=True, exist_ok=True)
    unreadable = WORK / "taxid/8/wmasker.obinary"
    if not unreadable.exists():
        unreadable.parent.mkdir(parents=True, exist_ok=True)
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
    if any(key in os.environ for key in (*REPORT_ENV, "NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH")) \
            or (Path.home() / ".ncbirc").exists():
        raise SystemExit(f"unset {REPORT_ENV} and remove ~/.ncbirc before freezing")
    make_counts()
    blastn = str(Path(args.ncbi_bin).resolve() / "blastn")
    rows = read_cases()
    if args.cases:
        unknown = set(args.cases) - {case for case, _, _ in rows}
        if unknown:
            raise SystemExit(f"no such case in cases.tsv: {sorted(unknown)}")
        rows = [row for row in rows if row[0] in args.cases]
    for case, env_text, argv in rows:
        env = {key: value for key, value in os.environ.items()
               if not key.startswith(("LOSAT_", "RAYON_")) and key != "WINDOW_MASKER_PATH"}
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
    freeze.add_argument("--cases", nargs="+", help="freeze only these rows of cases.tsv")
    args = parser.parse_args()
    return command_freeze(args)


if __name__ == "__main__":
    sys.exit(main())
