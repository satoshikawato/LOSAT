#!/usr/bin/env python3
"""
Find representative LOSAT TBLASTX HSPs that are extra (present in LOSAT output, missing in BLAST+).

This is the inverse of `find_tblastx_missing_hsp.py`.

We work with outfmt 6/7 tabular files (comment lines starting with '#').
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

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
    q_shift = np.where(
        out["q_strand"] == 1,
        (out["qstart"] - 1) % 3,
        (dna_len - out["qstart"]) % 3,
    ).astype(int)
    s_shift = np.where(
        out["s_strand"] == 1,
        (out["sstart"] - 1) % 3,
        (dna_len - out["sstart"]) % 3,
    ).astype(int)
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
    ap.add_argument("--blast", default="./blast_out/AP027280.AP027280.tblastx.n1.out")
    ap.add_argument("--losat", default="./losat_out/AP027280.AP027280.tlosatx.n1.out")
    ap.add_argument(
        "--fasta",
        default="./fasta/AP027280.fasta",
        help="FASTA used to infer genome length (needed for frame estimation)",
    )

    # Optional focus bin (defaults cover all).
    ap.add_argument("--pident-min", type=float, default=0.0)
    ap.add_argument("--pident-max", type=float, default=101.0)
    ap.add_argument("--len-min", type=int, default=0)
    ap.add_argument("--len-max", type=int, default=10**9)

    ap.add_argument("--q-overlap", type=float, default=0.80, help="Query interval overlap ratio threshold (default: 0.80)")
    ap.add_argument("--s-overlap", type=float, default=0.80, help="Subject interval overlap ratio threshold (default: 0.80)")
    ap.add_argument("--bin-size", type=int, default=5000, help="Index bin size in bp for candidate search (default: 5000)")
    ap.add_argument("--nearby-slack", type=int, default=300, help="bp slack when selecting nearby candidates (default: 300)")
    ap.add_argument("--top", type=int, default=20, help="Print top-N extra LOSAT HSPs (default: 20)")
    ap.add_argument("--show-nearest", type=int, default=10, help="For chosen HSP, show this many nearest BLAST+ hits (default: 10)")
    ap.add_argument("--choose", type=int, default=1, help="Choose which extra HSP to show details for (1-based; default: 1)")
    ap.add_argument("--out-csv", default="", help="Optional: write extra HSP list to CSV")
    args = ap.parse_args()

    for p in [args.blast, args.losat, args.fasta]:
        if not os.path.exists(p):
            raise SystemExit(f"File not found: {p}")

    dna_len = read_fasta_length(args.fasta)
    df_b = add_intervals_and_frames(load_outfmt(args.blast, "BLAST+"), dna_len)
    df_l = add_intervals_and_frames(load_outfmt(args.losat, "LOSAT"), dna_len)

    # Override binning with user-selected size (used for candidate indexing).
    df_b["q_bin"] = (df_b["q_lo"] // args.bin_size).astype(int)
    df_b["s_bin"] = (df_b["s_lo"] // args.bin_size).astype(int)

    # Optional focus filter (LOSAT side)
    l_focus = df_l[
        (df_l["pident"] >= args.pident_min)
        & (df_l["pident"] < args.pident_max)
        & (df_l["length"] >= args.len_min)
        & (df_l["length"] < args.len_max)
    ].copy()

    print("=== Focus filter (LOSAT side) ===")
    print(f"pident in [{args.pident_min}, {args.pident_max})")
    print(f"length in [{args.len_min}, {args.len_max}) (aa)")
    print(f"LOSAT focus hits: {len(l_focus)} / {len(df_l)}")
    print("")

    # ---------------------------------------------------------------------
    # Fast path: exact coordinate match
    # ---------------------------------------------------------------------
    # Most LOSAT hits are identical (by coordinates) to BLAST+; only a small tail is "extra".
    # Build a hash-set of exact BLAST+ coordinate tuples so we can skip the vast majority fast.
    blast_exact = set(
        zip(
            df_b["q_strand"].astype(int),
            df_b["s_strand"].astype(int),
            df_b["q_lo"].astype(int),
            df_b["q_hi"].astype(int),
            df_b["s_lo"].astype(int),
            df_b["s_hi"].astype(int),
        )
    )

    # ---------------------------------------------------------------------
    # Candidate index for overlap matching (only used for the few non-exact cases)
    # ---------------------------------------------------------------------
    # Build BLAST index: (q_strand, s_strand, q_bin, s_bin) -> list of row indices
    index: Dict[Tuple[int, int, int, int], List[int]] = {}
    for row in df_b.itertuples(index=True):
        key = (int(row.q_strand), int(row.s_strand), int(row.q_bin), int(row.s_bin))
        index.setdefault(key, []).append(int(row.Index))

    def blast_candidates_for(q_strand: int, s_strand: int, q_lo: int, q_hi: int, s_lo: int, s_hi: int) -> List[int]:
        q_bin_lo = int((q_lo - args.nearby_slack) // args.bin_size)
        q_bin_hi = int((q_hi + args.nearby_slack) // args.bin_size)
        s_bin_lo = int((s_lo - args.nearby_slack) // args.bin_size)
        s_bin_hi = int((s_hi + args.nearby_slack) // args.bin_size)
        out_idx: List[int] = []
        for qb in range(q_bin_lo, q_bin_hi + 1):
            for sb in range(s_bin_lo, s_bin_hi + 1):
                key = (q_strand, s_strand, qb, sb)
                out_idx.extend(index.get(key, []))
        return out_idx

    # Determine which LOSAT HSPs are "extra" (no close BLAST+ match)
    extra_rows: List[Tuple[int, float, int]] = []  # (losat_idx, best_overlap, best_dist)
    thresh = float(min(args.q_overlap, args.s_overlap))
    for row in l_focus.itertuples(index=True):
        q_strand = int(row.q_strand)
        s_strand = int(row.s_strand)
        q_lo = int(row.q_lo)
        q_hi = int(row.q_hi)
        s_lo = int(row.s_lo)
        s_hi = int(row.s_hi)

        # Fast exact-match skip
        if (q_strand, s_strand, q_lo, q_hi, s_lo, s_hi) in blast_exact:
            continue

        cand_idx = blast_candidates_for(q_strand, s_strand, q_lo, q_hi, s_lo, s_hi)
        if not cand_idx:
            extra_rows.append((int(row.Index), 0.0, 10**18))
            continue

        best_score = 0.0
        best_dist = 10**18
        for j in cand_idx:
            br = df_b.loc[j]
            q_ov = overlap_ratio(q_lo, q_hi, int(br["q_lo"]), int(br["q_hi"]))
            s_ov = overlap_ratio(s_lo, s_hi, int(br["s_lo"]), int(br["s_hi"]))
            score = float(min(q_ov, s_ov))
            dist = (
                abs(q_lo - int(br["q_lo"]))
                + abs(q_hi - int(br["q_hi"]))
                + abs(s_lo - int(br["s_lo"]))
                + abs(s_hi - int(br["s_hi"]))
            )
            if score > best_score or (score == best_score and dist < best_dist):
                best_score = score
                best_dist = int(dist)
                if best_score >= thresh:
                    break

        if best_score < thresh:
            extra_rows.append((int(row.Index), float(best_score), int(best_dist)))

    print("=== Extra summary (within focus filter) ===")
    print(f"Extra (no BLAST+ match with overlap >= {min(args.q_overlap, args.s_overlap):.2f}): {len(extra_rows)} / {len(l_focus)}")

    if not extra_rows:
        print("No extra HSPs found under the current thresholds. Try lowering --q-overlap/--s-overlap or increasing --nearby-slack.")
        return

    extra_df = pd.DataFrame(extra_rows, columns=["losat_idx", "best_overlap", "best_dist"])
    extra_df = extra_df.merge(
        df_l[["bitscore", "evalue", "pident", "length", "qstart", "qend", "sstart", "send", "q_frame", "s_frame"]],
        left_on="losat_idx",
        right_index=True,
    )
    # Sort by LOSAT bitscore (desc), then by best_overlap (asc), then by best_dist (asc)
    extra_df = extra_df.sort_values(by=["bitscore", "best_overlap", "best_dist"], ascending=[False, True, True])

    if args.out_csv:
        out_dir = os.path.dirname(args.out_csv)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        extra_df.to_csv(args.out_csv, index=False)
        print(f"Wrote extra list CSV: {args.out_csv}")

    print("")
    print(f"=== Top {min(args.top, len(extra_df))} extra LOSAT HSPs (sorted by bitscore) ===")
    head = extra_df.head(max(1, int(args.top)))
    for rank, r in enumerate(head.itertuples(index=False), 1):
        print(
            f"{rank:2d}. bits={r.bitscore:.1f} e={r.evalue:g} pid={r.pident:.2f} len={r.length} "
            f"q={r.qstart}-{r.qend} (q_frame={r.q_frame}) "
            f"s={r.sstart}-{r.send} (s_frame={r.s_frame}) "
            f"best_overlap={r.best_overlap:.2f} best_dist={r.best_dist}"
        )

    choose_idx = int(args.choose) - 1
    if choose_idx < 0 or choose_idx >= len(extra_df):
        raise SystemExit(f"--choose out of range: {args.choose} (extra count={len(extra_df)})")

    chosen = extra_df.iloc[choose_idx]
    losat_idx = int(chosen["losat_idx"])
    lr = df_l.loc[losat_idx]

    print("")
    print("=== Chosen extra HSP (LOSAT) ===")
    print(
        f"bits={lr['bitscore']:.1f} e={lr['evalue']:g} pid={lr['pident']:.2f} len={lr['length']} "
        f"q={lr['qstart']}-{lr['qend']} (q_frame={lr['q_frame']}) "
        f"s={lr['sstart']}-{lr['send']} (s_frame={lr['s_frame']})"
    )
    print(f"q_interval=[{int(lr['q_lo'])}, {int(lr['q_hi'])}] s_interval=[{int(lr['s_lo'])}, {int(lr['s_hi'])}]")
    print(f"Best BLAST+ overlap score (min(q_ov,s_ov)) observed: {float(chosen['best_overlap']):.2f}")

    # Nearest BLAST+ hits (same strands, nearby bins)
    cand_idx = blast_candidates_for(
        int(lr["q_strand"]),
        int(lr["s_strand"]),
        int(lr["q_lo"]),
        int(lr["q_hi"]),
        int(lr["s_lo"]),
        int(lr["s_hi"]),
    )
    matches: List[CandidateMatch] = []
    for j in cand_idx:
        br = df_b.loc[j]
        q_ov = overlap_ratio(int(lr["q_lo"]), int(lr["q_hi"]), int(br["q_lo"]), int(br["q_hi"]))
        s_ov = overlap_ratio(int(lr["s_lo"]), int(lr["s_hi"]), int(br["s_lo"]), int(br["s_hi"]))
        score = float(min(q_ov, s_ov))
        dist = (
            abs(int(lr["q_lo"]) - int(br["q_lo"]))
            + abs(int(lr["q_hi"]) - int(br["q_hi"]))
            + abs(int(lr["s_lo"]) - int(br["s_lo"]))
            + abs(int(lr["s_hi"]) - int(br["s_hi"]))
        )
        matches.append(CandidateMatch(idx=int(j), q_ov=float(q_ov), s_ov=float(s_ov), score=score, dist=int(dist)))
    matches_sorted = sorted(matches, key=lambda m: (-m.score, m.dist))

    print("")
    print(f"=== Nearest BLAST+ candidates (same strands; showing top {min(args.show_nearest, len(matches_sorted))}) ===")
    for k, m in enumerate(matches_sorted[: max(1, int(args.show_nearest))], 1):
        br = df_b.loc[m.idx]
        print(
            f"{k:2d}. score={m.score:.2f} q_ov={m.q_ov:.2f} s_ov={m.s_ov:.2f} dist={m.dist} "
            f"BLAST bits={br['bitscore']:.1f} e={br['evalue']:g} pid={br['pident']:.2f} len={br['length']} "
            f"q={br['qstart']}-{br['qend']} (q_frame={br['q_frame']}) "
            f"s={br['sstart']}-{br['send']} (s_frame={br['s_frame']})"
        )

    # Simple classification hint
    best = matches_sorted[0] if matches_sorted else None
    if best is None or best.score < 0.20:
        cls = "complete_extra"
    elif best.score < min(args.q_overlap, args.s_overlap):
        cls = "nearby_but_not_matching"
    else:
        cls = "has_matching_blast_hit"

    print("")
    print("=== Classification hint ===")
    print(f"{cls} (threshold={min(args.q_overlap, args.s_overlap):.2f})")


if __name__ == "__main__":
    main()


