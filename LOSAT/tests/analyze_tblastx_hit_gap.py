#!/usr/bin/env python3
"""
Quantify BLAST+ vs LOSAT TBLASTX output differences ("fewer hits" debugging).

Input files:
  - BLAST+ tblastx outfmt 6/7 (tab-separated; comment lines start with '#')
  - LOSAT tblastx output with the same 12 columns (tab-separated; may have no '#')

Outputs:
  - Prints overall hit counts and distribution deltas
  - Prints top deficit bins (where BLAST+ has more hits than LOSAT)
  - Optionally writes a CSV summary of 2D bin deficits

This script is designed to be the first step in selecting a representative
"missing HSP" for deeper pipeline diagnostics.
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from typing import List

import numpy as np
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


def _load_outfmt(path: str, tool_name: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", names=COLUMNS, header=None)

    # Convert numeric columns; drop rows that don't parse cleanly.
    numeric_cols = [
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
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(
        subset=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ]
    )

    # Types we care about later.
    df["length"] = df["length"].astype(int)
    df["qstart"] = df["qstart"].astype(int)
    df["qend"] = df["qend"].astype(int)
    df["sstart"] = df["sstart"].astype(int)
    df["send"] = df["send"].astype(int)

    df["Tool"] = tool_name
    return df


def _edges_identity(step: float) -> np.ndarray:
    # Include 100 as the last edge.
    n = int(math.ceil(100.0 / step))
    edges = np.linspace(0.0, step * n, n + 1)
    edges[-1] = 100.0
    return edges


def _edges_len_log2(max_len: int) -> np.ndarray:
    # Log2-ish bins: [1,2,4,...] with 0 as underflow bin start.
    edges: List[int] = [0, 1, 2]
    x = 4
    while x < max_len + 1:
        edges.append(x)
        x *= 2
    if edges[-1] < max_len + 1:
        edges.append(max_len + 1)
    return np.array(edges, dtype=int)


def _edges_bitscore_log2(max_bits: float) -> np.ndarray:
    # Bitscore can be float in LOSAT output; bin in powers of 2.
    edges: List[float] = [0.0, 1.0, 2.0]
    x = 4.0
    while x < max_bits + 1.0:
        edges.append(x)
        x *= 2.0
    if edges[-1] < max_bits + 1.0:
        edges.append(max_bits + 1.0)
    return np.array(edges, dtype=float)


@dataclass(frozen=True)
class DeficitBin:
    pident_bin: str
    len_bin: str
    blast_count: int
    losat_count: int
    deficit: int
    blast_acc_len: int
    losat_acc_len: int
    deficit_acc_len: int


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--blast",
        default="./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out",
        help="BLAST+ tblastx output (outfmt 6/7; tab-separated)",
    )
    ap.add_argument(
        "--losat",
        default="./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out",
        help="LOSAT tblastx output (tab-separated, same 12 columns)",
    )
    ap.add_argument(
        "--identity-step",
        type=float,
        default=5.0,
        help="Identity bin step (percentage points, default: 5.0)",
    )
    ap.add_argument(
        "--min-bin-blast-count",
        type=int,
        default=50,
        help=(
            "Only consider 2D bins with at least this many BLAST+ hits when selecting "
            "the top deficit bin (default: 50)"
        ),
    )
    ap.add_argument(
        "--top",
        type=int,
        default=20,
        help="How many deficit bins to print (default: 20)",
    )
    ap.add_argument(
        "--out-csv",
        default="",
        help="Optional: write 2D deficit table as CSV to this path",
    )
    args = ap.parse_args()

    if not os.path.exists(args.blast):
        raise SystemExit(f"BLAST+ file not found: {args.blast}")
    if not os.path.exists(args.losat):
        raise SystemExit(f"LOSAT file not found: {args.losat}")

    df_b = _load_outfmt(args.blast, "BLAST+")
    df_l = _load_outfmt(args.losat, "LOSAT")

    print("=== Input summary ===")
    print(f"BLAST+: {args.blast}")
    print(f"LOSAT : {args.losat}")
    print("")

    print("=== Overall hit counts ===")
    n_b = len(df_b)
    n_l = len(df_l)
    print(f"BLAST+ hits: {n_b}")
    print(f"LOSAT  hits: {n_l}")
    print(
        f"Delta (BLAST+ - LOSAT): {n_b - n_l} ({(n_b - n_l) / max(1, n_b) * 100:.2f}% of BLAST+)"
    )
    print("")

    # Common bin edges derived from union of both datasets.
    max_len = int(max(df_b["length"].max(), df_l["length"].max()))
    max_bits = float(max(df_b["bitscore"].max(), df_l["bitscore"].max()))
    ident_edges = _edges_identity(args.identity_step)
    len_edges = _edges_len_log2(max_len)
    bits_edges = _edges_bitscore_log2(max_bits)

    def add_bins(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["pident_bin"] = pd.cut(
            out["pident"], bins=ident_edges, include_lowest=True, right=False
        )
        out["len_bin"] = pd.cut(
            out["length"], bins=len_edges, include_lowest=True, right=False
        )
        out["bits_bin"] = pd.cut(
            out["bitscore"], bins=bits_edges, include_lowest=True, right=False
        )
        return out

    df_b2 = add_bins(df_b)
    df_l2 = add_bins(df_l)

    print("=== 1D distribution deltas (counts) ===")
    # Identity
    b_ident = df_b2["pident_bin"].value_counts().sort_index()
    l_ident = df_l2["pident_bin"].value_counts().sort_index()
    ident = pd.DataFrame({"BLAST+": b_ident, "LOSAT": l_ident}).fillna(0).astype(int)
    ident["Delta"] = ident["BLAST+"] - ident["LOSAT"]
    print("\n-- pident bins --")
    print(ident.to_string())

    # Length
    b_len = df_b2["len_bin"].value_counts().sort_index()
    l_len = df_l2["len_bin"].value_counts().sort_index()
    lb = pd.DataFrame({"BLAST+": b_len, "LOSAT": l_len}).fillna(0).astype(int)
    lb["Delta"] = lb["BLAST+"] - lb["LOSAT"]
    print("\n-- length bins (log2-ish) --")
    print(lb.to_string())

    # Bitscore
    b_bits = df_b2["bits_bin"].value_counts().sort_index()
    l_bits = df_l2["bits_bin"].value_counts().sort_index()
    bb = pd.DataFrame({"BLAST+": b_bits, "LOSAT": l_bits}).fillna(0).astype(int)
    bb["Delta"] = bb["BLAST+"] - bb["LOSAT"]
    print("\n-- bitscore bins (log2-ish) --")
    print(bb.to_string())
    print("")

    # 2D deficit table: pident_bin x len_bin (counts and accumulated length)
    def pivot_counts(df: pd.DataFrame) -> pd.DataFrame:
        return pd.pivot_table(
            df,
            index="pident_bin",
            columns="len_bin",
            values="length",
            aggfunc="size",
            fill_value=0,
            observed=False,
        )

    def pivot_acc_len(df: pd.DataFrame) -> pd.DataFrame:
        return pd.pivot_table(
            df,
            index="pident_bin",
            columns="len_bin",
            values="length",
            aggfunc="sum",
            fill_value=0,
            observed=False,
        )

    b2d = pivot_counts(df_b2)
    l2d = pivot_counts(df_l2)
    d2d = b2d.subtract(l2d, fill_value=0).astype(int)

    b2d_len = pivot_acc_len(df_b2)
    l2d_len = pivot_acc_len(df_l2)
    d2d_len = b2d_len.subtract(l2d_len, fill_value=0).astype(int)

    # Flatten and pick top deficit bins (BLAST+ > LOSAT)
    rows: List[DeficitBin] = []
    for pbin in d2d.index:
        for lbin in d2d.columns:
            blast_count = int(b2d.at[pbin, lbin]) if (pbin in b2d.index and lbin in b2d.columns) else 0
            losat_count = int(l2d.at[pbin, lbin]) if (pbin in l2d.index and lbin in l2d.columns) else 0
            deficit = int(d2d.at[pbin, lbin])
            if deficit <= 0:
                continue
            blast_acc_len = int(b2d_len.at[pbin, lbin]) if (pbin in b2d_len.index and lbin in b2d_len.columns) else 0
            losat_acc_len = int(l2d_len.at[pbin, lbin]) if (pbin in l2d_len.index and lbin in l2d_len.columns) else 0
            deficit_acc_len = int(d2d_len.at[pbin, lbin])
            rows.append(
                DeficitBin(
                    pident_bin=str(pbin),
                    len_bin=str(lbin),
                    blast_count=blast_count,
                    losat_count=losat_count,
                    deficit=deficit,
                    blast_acc_len=blast_acc_len,
                    losat_acc_len=losat_acc_len,
                    deficit_acc_len=deficit_acc_len,
                )
            )

    rows_sorted = sorted(rows, key=lambda r: (r.deficit, r.deficit_acc_len), reverse=True)
    print("=== Top deficit bins (2D: pident x length) ===")
    if not rows_sorted:
        print("No deficit bins found (BLAST+ > LOSAT).")
    else:
        top = rows_sorted[: max(1, int(args.top))]
        for i, r in enumerate(top, 1):
            print(
                f"{i:2d}. pident={r.pident_bin}, len={r.len_bin} :: "
                f"BLAST+={r.blast_count}, LOSAT={r.losat_count}, deficit={r.deficit} "
                f"(acc_len_deficit={r.deficit_acc_len})"
            )

        # Pick "best" bin with enough BLAST support to be stable.
        stable = [r for r in rows_sorted if r.blast_count >= args.min_bin_blast_count]
        chosen = stable[0] if stable else rows_sorted[0]
        print("")
        print("=== Suggested focus bin (largest stable deficit) ===")
        print(
            f"pident={chosen.pident_bin}, len={chosen.len_bin} :: "
            f"BLAST+={chosen.blast_count}, LOSAT={chosen.losat_count}, deficit={chosen.deficit}"
        )

    if args.out_csv:
        out_path = args.out_csv
        out_dir = os.path.dirname(out_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        # Write flattened table for easy inspection.
        df_out = pd.DataFrame([r.__dict__ for r in rows_sorted])
        df_out.to_csv(out_path, index=False)
        print("")
        print(f"Wrote deficit CSV: {out_path}")


if __name__ == "__main__":
    main()


