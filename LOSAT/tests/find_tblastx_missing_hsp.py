#!/usr/bin/env python3
"""
Find representative BLAST+ TBLASTX HSPs that are missing in LOSAT output.

This script is step (2) in the "fewer hits" workflow:
  1) Quantify deficit bins (see analyze_tblastx_hit_gap.py)
  2) Pick one missing HSP (this script)
  3) Use the chosen HSP to trace where it is dropped inside LOSAT

We work with outfmt 6/7 tabular files (comment lines starting with '#').
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

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


def read_fasta_length(path: str) -> int:
    # Minimal FASTA reader: count non-header bases.
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                continue
            n += len(line.strip())
    if n <= 0:
        raise ValueError(f"Could not determine sequence length from FASTA: {path}")
    return n


def load_outfmt(path: str, tool_name: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", names=COLUMNS, header=None)

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
    df = df.dropna(subset=["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

    df["length"] = df["length"].astype(int)
    df["qstart"] = df["qstart"].astype(int)
    df["qend"] = df["qend"].astype(int)
    df["sstart"] = df["sstart"].astype(int)
    df["send"] = df["send"].astype(int)

    df["Tool"] = tool_name
    return df


def add_intervals_and_frames(df: pd.DataFrame, dna_len: int) -> pd.DataFrame:
    out = df.copy()
    out["q_lo"] = out[["qstart", "qend"]].min(axis=1)
    out["q_hi"] = out[["qstart", "qend"]].max(axis=1)
    out["s_lo"] = out[["sstart", "send"]].min(axis=1)
    out["s_hi"] = out[["sstart", "send"]].max(axis=1)
    out["q_strand"] = np.where(out["qstart"] <= out["qend"], 1, -1).astype(int)
    out["s_strand"] = np.where(out["sstart"] <= out["send"], 1, -1).astype(int)

    # Frame estimation (LOSAT convert_coords semantics).
    # Plus: start_bp = aa_start*3 + shift + 1  => shift = (q_start - 1) % 3
    # Minus: start_bp = dna_len - (aa_start*3 + shift) => shift = (dna_len - q_start) % 3
    q_shift = np.where(out["q_strand"] == 1, (out["qstart"] - 1) % 3, (dna_len - out["qstart"]) % 3).astype(int)
    s_shift = np.where(out["s_strand"] == 1, (out["sstart"] - 1) % 3, (dna_len - out["sstart"]) % 3).astype(int)
    out["q_frame"] = np.where(out["q_strand"] == 1, q_shift + 1, -(q_shift + 1)).astype(int)
    out["s_frame"] = np.where(out["s_strand"] == 1, s_shift + 1, -(s_shift + 1)).astype(int)

    out["q_bin"] = (out["q_lo"] // 5000).astype(int)  # default bin (can override later)
    return out


def overlap_ratio(a1: int, b1: int, a2: int, b2: int) -> float:
    # inclusive intervals
    lo = max(a1, a2)
    hi = min(b1, b2)
    if hi < lo:
        return 0.0
    ov = hi - lo + 1
    l1 = b1 - a1 + 1
    l2 = b2 - a2 + 1
    return ov / float(min(l1, l2))


@dataclass(frozen=True)
class CandidateMatch:
    idx: int
    q_ov: float
    s_ov: float
    score: float  # combined overlap score
    dist: int


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast", default="./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out")
    ap.add_argument("--losat", default="./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out")
    ap.add_argument("--fasta", default="./fasta/NZ_CP006932.fasta", help="FASTA used to infer genome length (needed for frame estimation)")

    ap.add_argument("--pident-min", type=float, default=40.0)
    ap.add_argument("--pident-max", type=float, default=45.0)
    ap.add_argument("--len-min", type=int, default=16)
    ap.add_argument("--len-max", type=int, default=32)

    ap.add_argument("--q-overlap", type=float, default=0.80, help="Query interval overlap ratio threshold (default: 0.80)")
    ap.add_argument("--s-overlap", type=float, default=0.80, help="Subject interval overlap ratio threshold (default: 0.80)")
    ap.add_argument("--bin-size", type=int, default=5000, help="Index bin size in bp for candidate search (default: 5000)")
    ap.add_argument("--nearby-slack", type=int, default=300, help="bp slack when selecting nearby candidates (default: 300)")
    ap.add_argument("--top", type=int, default=20, help="Print top-N missing BLAST+ HSPs (default: 20)")
    ap.add_argument("--show-nearest", type=int, default=10, help="For chosen HSP, show this many nearest LOSAT hits (default: 10)")
    ap.add_argument("--choose", type=int, default=1, help="Choose which missing HSP to show details for (1-based; default: 1)")
    ap.add_argument("--out-csv", default="", help="Optional: write missing BLAST+ HSP list to CSV")
    args = ap.parse_args()

    for p in [args.blast, args.losat, args.fasta]:
        if not os.path.exists(p):
            raise SystemExit(f"File not found: {p}")

    dna_len = read_fasta_length(args.fasta)
    df_b = add_intervals_and_frames(load_outfmt(args.blast, "BLAST+"), dna_len)
    df_l = add_intervals_and_frames(load_outfmt(args.losat, "LOSAT"), dna_len)

    # Override binning with user-selected size.
    df_l["q_bin"] = (df_l["q_lo"] // args.bin_size).astype(int)

    # Focus bin filter (BLAST+ side)
    b_focus = df_b[
        (df_b["pident"] >= args.pident_min)
        & (df_b["pident"] < args.pident_max)
        & (df_b["length"] >= args.len_min)
        & (df_b["length"] < args.len_max)
    ].copy()

    print("=== Focus bin ===")
    print(f"pident in [{args.pident_min}, {args.pident_max})")
    print(f"length in [{args.len_min}, {args.len_max}) (aa)")
    print(f"BLAST+ focus hits: {len(b_focus)}")
    print("")

    # Build losat index: (q_strand, s_strand, q_bin) -> list of row indices
    index: Dict[Tuple[int, int, int], List[int]] = {}
    for i, row in df_l.iterrows():
        key = (int(row["q_strand"]), int(row["s_strand"]), int(row["q_bin"]))
        index.setdefault(key, []).append(int(i))

    def losat_candidates_for(bl_row: pd.Series) -> List[int]:
        q_bin_lo = int((int(bl_row["q_lo"]) - args.nearby_slack) // args.bin_size)
        q_bin_hi = int((int(bl_row["q_hi"]) + args.nearby_slack) // args.bin_size)
        out_idx: List[int] = []
        for qb in range(q_bin_lo, q_bin_hi + 1):
            key = (int(bl_row["q_strand"]), int(bl_row["s_strand"]), qb)
            out_idx.extend(index.get(key, []))
        return out_idx

    # Determine which BLAST HSPs are "missing" (no close LOSAT match)
    missing_rows: List[Tuple[int, float, int]] = []  # (blast_row_idx, best_score, best_dist)
    for i, row in b_focus.iterrows():
        cand_idx = losat_candidates_for(row)
        if not cand_idx:
            missing_rows.append((int(i), 0.0, 10**18))
            continue

        best_score = 0.0
        best_dist = 10**18
        q_lo = int(row["q_lo"])
        q_hi = int(row["q_hi"])
        s_lo = int(row["s_lo"])
        s_hi = int(row["s_hi"])
        for j in cand_idx:
            lr = df_l.loc[j]
            q_ov = overlap_ratio(q_lo, q_hi, int(lr["q_lo"]), int(lr["q_hi"]))
            s_ov = overlap_ratio(s_lo, s_hi, int(lr["s_lo"]), int(lr["s_hi"]))
            score = min(q_ov, s_ov)
            dist = abs(q_lo - int(lr["q_lo"])) + abs(q_hi - int(lr["q_hi"])) + abs(s_lo - int(lr["s_lo"])) + abs(s_hi - int(lr["s_hi"]))
            if score > best_score or (score == best_score and dist < best_dist):
                best_score = float(score)
                best_dist = int(dist)

        if best_score < min(args.q_overlap, args.s_overlap):
            missing_rows.append((int(i), float(best_score), int(best_dist)))

    print("=== Missing summary (within focus bin) ===")
    print(f"Missing (no LOSAT match with overlap >= {min(args.q_overlap, args.s_overlap):.2f}): {len(missing_rows)} / {len(b_focus)}")

    if not missing_rows:
        print("No missing HSPs found under the current thresholds. Try lowering --q-overlap/--s-overlap or increasing --nearby-slack.")
        return

    # Sort missing by BLAST bitscore (desc), then by best_score (asc), then by distance (asc)
    miss_df = pd.DataFrame(missing_rows, columns=["blast_idx", "best_overlap", "best_dist"])
    miss_df = miss_df.merge(df_b[["bitscore", "evalue", "pident", "length", "qstart", "qend", "sstart", "send", "q_frame", "s_frame"]], left_on="blast_idx", right_index=True)
    miss_df = miss_df.sort_values(by=["bitscore", "best_overlap", "best_dist"], ascending=[False, True, True])

    if args.out_csv:
        out_dir = os.path.dirname(args.out_csv)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        miss_df.to_csv(args.out_csv, index=False)
        print(f"Wrote missing list CSV: {args.out_csv}")

    print("")
    print(f"=== Top {min(args.top, len(miss_df))} missing BLAST+ HSPs (sorted by bitscore) ===")
    head = miss_df.head(max(1, int(args.top)))
    for rank, r in enumerate(head.itertuples(index=False), 1):
        print(
            f"{rank:2d}. bits={r.bitscore:.1f} e={r.evalue:g} pid={r.pident:.2f} len={r.length} "
            f"q={r.qstart}-{r.qend} (q_frame={r.q_frame}) "
            f"s={r.sstart}-{r.send} (s_frame={r.s_frame}) "
            f"best_overlap={r.best_overlap:.2f} best_dist={r.best_dist}"
        )

    choose_idx = int(args.choose) - 1
    if choose_idx < 0 or choose_idx >= len(miss_df):
        raise SystemExit(f"--choose out of range: {args.choose} (missing count={len(miss_df)})")

    chosen = miss_df.iloc[choose_idx]
    blast_idx = int(chosen["blast_idx"])
    bl = df_b.loc[blast_idx]

    print("")
    print("=== Chosen missing HSP (BLAST+) ===")
    print(
        f"bits={bl['bitscore']:.1f} e={bl['evalue']:g} pid={bl['pident']:.2f} len={bl['length']} "
        f"q={bl['qstart']}-{bl['qend']} (q_frame={bl['q_frame']}) "
        f"s={bl['sstart']}-{bl['send']} (s_frame={bl['s_frame']})"
    )
    print(f"q_interval=[{int(bl['q_lo'])}, {int(bl['q_hi'])}] s_interval=[{int(bl['s_lo'])}, {int(bl['s_hi'])}]")
    print(f"Best LOSAT overlap score (min(q_ov,s_ov)) observed: {float(chosen['best_overlap']):.2f}")

    # Nearest LOSAT hits (same strands, nearby bins)
    cand_idx = losat_candidates_for(bl)
    matches: List[CandidateMatch] = []
    for j in cand_idx:
        lr = df_l.loc[j]
        q_ov = overlap_ratio(int(bl["q_lo"]), int(bl["q_hi"]), int(lr["q_lo"]), int(lr["q_hi"]))
        s_ov = overlap_ratio(int(bl["s_lo"]), int(bl["s_hi"]), int(lr["s_lo"]), int(lr["s_hi"]))
        score = float(min(q_ov, s_ov))
        dist = (
            abs(int(bl["q_lo"]) - int(lr["q_lo"]))
            + abs(int(bl["q_hi"]) - int(lr["q_hi"]))
            + abs(int(bl["s_lo"]) - int(lr["s_lo"]))
            + abs(int(bl["s_hi"]) - int(lr["s_hi"]))
        )
        matches.append(CandidateMatch(idx=int(j), q_ov=float(q_ov), s_ov=float(s_ov), score=score, dist=int(dist)))
    matches_sorted = sorted(matches, key=lambda m: (-m.score, m.dist))

    print("")
    print(f"=== Nearest LOSAT candidates (same strands; showing top {min(args.show_nearest, len(matches_sorted))}) ===")
    for k, m in enumerate(matches_sorted[: max(1, int(args.show_nearest))], 1):
        lr = df_l.loc[m.idx]
        print(
            f"{k:2d}. score={m.score:.2f} q_ov={m.q_ov:.2f} s_ov={m.s_ov:.2f} dist={m.dist} "
            f"LOSAT bits={lr['bitscore']:.1f} e={lr['evalue']:g} pid={lr['pident']:.2f} len={lr['length']} "
            f"q={lr['qstart']}-{lr['qend']} (q_frame={lr['q_frame']}) "
            f"s={lr['sstart']}-{lr['send']} (s_frame={lr['s_frame']})"
        )

    # Simple classification hint
    best = matches_sorted[0] if matches_sorted else None
    if best is None or best.score < 0.20:
        cls = "complete_missing"
    elif best.score < min(args.q_overlap, args.s_overlap):
        cls = "nearby_but_not_matching"
    else:
        cls = "has_matching_losat_hit"

    print("")
    print(f"=== Classification hint ===")
    print(f"{cls} (threshold={min(args.q_overlap, args.s_overlap):.2f})")


if __name__ == "__main__":
    main()


