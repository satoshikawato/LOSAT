#!/usr/bin/env python3
"""
Detect "overlong tail" regressions in TBLASTX tabular outputs.

This is a safety check for changes that affect stop handling / seeding / extension.
It is intentionally simple and conservative: it flags only *very* long HSPs that
are unlikely to be legitimate in typical bacterial self-comparisons.

Input:
  - BLAST+ outfmt 6/7 or LOSAT tblastx tabular output with the 12 standard columns
    (comment lines beginning with '#' are ignored).
Output:
  - Prints max length, some summary stats, and top-N longest HSPs.
  - Exits non-zero when thresholds are exceeded.
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd


COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def load_outfmt(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", names=COLUMNS, header=None)
    # Ensure numeric types.
    for c in ["pident", "length", "evalue", "bitscore"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["pident", "length"])
    df["length"] = df["length"].astype(int)
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="TBLASTX tabular output file (outfmt 6/7 style)")
    ap.add_argument("--max-length-threshold", type=int, default=10000, help="Fail if max(length) exceeds this (default: 10000)")
    ap.add_argument("--high-ident-threshold", type=float, default=95.0, help="High-identity threshold (default: 95.0)")
    ap.add_argument(
        "--high-ident-length-threshold",
        type=int,
        default=5000,
        help="Fail if any HSP with pident >= high-ident-threshold has length >= this (default: 5000)",
    )
    ap.add_argument("--top", type=int, default=20, help="Print top-N longest HSPs (default: 20)")
    args = ap.parse_args()

    if not os.path.exists(args.path):
        raise SystemExit(f"File not found: {args.path}")

    df = load_outfmt(args.path)
    if df.empty:
        print("No data rows parsed (empty after filtering).")
        return

    max_len = int(df["length"].max())
    n = len(df)
    hi = df[df["pident"] >= float(args.high_ident_threshold)]
    hi_max = int(hi["length"].max()) if not hi.empty else 0

    print("=== Overlong tail check ===")
    print(f"file: {args.path}")
    print(f"rows: {n}")
    print(f"max_length: {max_len}")
    print(f"high_identity: pident >= {args.high_ident_threshold}")
    print(f"high_identity_rows: {len(hi)}")
    print(f"high_identity_max_length: {hi_max}")
    print("")

    topn = df.sort_values(by=["length", "pident", "bitscore"], ascending=[False, False, False]).head(max(1, int(args.top)))
    print(f"=== Top {min(len(topn), int(args.top))} by length ===")
    for r in topn.itertuples(index=False):
        print(
            f"len={int(r.length)} pid={float(r.pident):.3f} bits={float(r.bitscore):.1f} e={float(r.evalue):g} "
            f"q={r.qseqid}:{int(r.qstart)}-{int(r.qend)} s={r.sseqid}:{int(r.sstart)}-{int(r.send)}"
        )

    failed = False
    reasons = []
    if max_len > int(args.max_length_threshold):
        failed = True
        reasons.append(f"max_length {max_len} > max_length_threshold {args.max_length_threshold}")

    if not hi.empty and hi_max >= int(args.high_ident_length_threshold):
        failed = True
        reasons.append(
            f"high_identity_max_length {hi_max} >= high_ident_length_threshold {args.high_ident_length_threshold}"
        )

    if failed:
        print("")
        print("=== FAIL ===")
        for x in reasons:
            print(f"- {x}")
        sys.exit(2)

    print("")
    print("=== OK ===")


if __name__ == "__main__":
    main()


