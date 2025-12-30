#!/usr/bin/env python3
"""
Analyze whether TBLASTX hits overlap annotated sequence features (GFF3).

This script is intentionally dependency-free (stdlib only).

It parses:
  - NCBI BLAST+ TBLASTX output (outfmt 6 or 7; lines starting with '#' are ignored)
  - GFF3 feature annotations (NCBI-style GFF3 is supported)

Then it reports:
  - Overlap rate (any overlap, and >=80% coverage) per feature type
  - Fraction of hits that are mostly non-CDS
  - Top-N longest hits with the list of overlapped CDS names

Notes:
  - BLAST coordinates are 1-based and can be reversed for minus-strand HSPs.
    We interpret intervals as [min(start,end), max(start,end)].
  - Even if the genome is circular, BLAST treats the provided sequence as linear
    unless you explicitly duplicate/rotate it, so we do NOT treat start>end as wrap.

Example:
  python3 tests/analyze_tblastx_feature_overlap.py \
    --gff3 /mnt/d/study/2025-05-04_gbdraw_test/NZ_CP006932.gff3 \
    --blast tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out \
    --top 10
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


FEATURE_TYPES_DEFAULT = ["CDS", "gene", "tRNA", "rRNA", "tmRNA"]


def parse_attrs(attr_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for item in attr_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def parse_genome_len_from_gff3(gff3_path: str) -> Optional[int]:
    with open(gff3_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##sequence-region"):
                # ##sequence-region <seqid> <start> <end>
                parts = line.strip().split()
                if len(parts) >= 4:
                    return int(parts[3])
    return None


def build_prefix_sum(mask: bytearray, genome_len: int) -> List[int]:
    # 1-based prefix sum (ps[0] = 0)
    ps = [0] * (genome_len + 2)
    s = 0
    for i in range(1, genome_len + 1):
        s += mask[i]
        ps[i] = s
    return ps


def interval_cov(ps: List[int], a: int, b: int) -> int:
    # a,b are 1-based inclusive, a<=b
    return ps[b] - ps[a - 1]


@dataclass(frozen=True)
class CdsFeature:
    start: int
    end: int
    strand: str
    locus_tag: str
    gene: str
    product: str

    def display_name(self) -> str:
        return self.gene or self.locus_tag or "(unknown)"


@dataclass(frozen=True)
class BlastHit:
    alen_aa: int
    q_start: int
    q_end: int
    s_start: int
    s_end: int

    @property
    def q_interval(self) -> Tuple[int, int]:
        return (min(self.q_start, self.q_end), max(self.q_start, self.q_end))

    @property
    def s_interval(self) -> Tuple[int, int]:
        return (min(self.s_start, self.s_end), max(self.s_start, self.s_end))


def parse_gff3_masks_and_cds(
    gff3_path: str, genome_len: int, feature_types: List[str]
) -> Tuple[Dict[str, bytearray], List[CdsFeature]]:
    masks: Dict[str, bytearray] = {t: bytearray(genome_len + 2) for t in feature_types}
    cds_list: List[CdsFeature] = []

    with open(gff3_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            _seqid, _source, ftype, start, end, _score, strand, _phase, attrs = cols
            if ftype not in masks:
                continue
            start_i = int(start)
            end_i = int(end)
            if start_i < 1 or end_i > genome_len or start_i > end_i:
                continue
            masks[ftype][start_i : end_i + 1] = b"\x01" * (end_i - start_i + 1)

            if ftype == "CDS":
                a = parse_attrs(attrs)
                cds_list.append(
                    CdsFeature(
                        start=start_i,
                        end=end_i,
                        strand=strand,
                        locus_tag=a.get("locus_tag", ""),
                        gene=a.get("gene", ""),
                        product=a.get("product", ""),
                    )
                )

    return masks, cds_list


def parse_blast_hits(blast_path: str) -> List[BlastHit]:
    hits: List[BlastHit] = []
    with open(blast_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            # outfmt 6/7 with fields:
            # qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
            if len(parts) < 12:
                continue
            try:
                alen = int(parts[3])
                q_start = int(parts[6])
                q_end = int(parts[7])
                s_start = int(parts[8])
                s_end = int(parts[9])
            except ValueError:
                continue
            hits.append(
                BlastHit(
                    alen_aa=alen,
                    q_start=q_start,
                    q_end=q_end,
                    s_start=s_start,
                    s_end=s_end,
                )
            )
    return hits


def build_cds_bin_index(cds: List[CdsFeature], bin_size: int) -> Dict[int, List[int]]:
    bins: Dict[int, List[int]] = defaultdict(list)
    for idx, feat in enumerate(cds):
        b0 = feat.start // bin_size
        b1 = feat.end // bin_size
        for b in range(b0, b1 + 1):
            bins[b].append(idx)
    return bins


def overlapping_cds_indices(
    cds: List[CdsFeature], cds_bins: Dict[int, List[int]], bin_size: int, a: int, b: int
) -> List[int]:
    # a<=b, 1-based inclusive
    out = set()
    for bin_id in range(a // bin_size, b // bin_size + 1):
        for idx in cds_bins.get(bin_id, []):
            feat = cds[idx]
            if feat.start <= b and feat.end >= a:
                out.add(idx)
    return sorted(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff3", required=True, help="GFF3 annotation file")
    ap.add_argument("--blast", required=True, help="BLAST+ tblastx output (outfmt 6 or 7)")
    ap.add_argument(
        "--feature-types",
        default=",".join(FEATURE_TYPES_DEFAULT),
        help="Comma-separated feature types to analyze (default: CDS,gene,tRNA,rRNA,tmRNA)",
    )
    ap.add_argument(
        "--top",
        type=int,
        default=10,
        help="Show top-N longest hits with overlapped CDS names (default: 10)",
    )
    ap.add_argument(
        "--bin-size",
        type=int,
        default=500,
        help="Bin size for CDS overlap enumeration (default: 500 bp)",
    )
    ap.add_argument(
        "--cds-coverage-threshold",
        type=float,
        default=0.8,
        help="Coverage threshold (0-1) for reporting '>=X%% covered' (default: 0.8)",
    )
    ap.add_argument(
        "--mostly-non-cds-threshold",
        type=float,
        default=0.2,
        help="Threshold (0-1) for counting 'mostly non-CDS' hits (default: 0.2)",
    )
    args = ap.parse_args()

    feature_types = [x.strip() for x in args.feature_types.split(",") if x.strip()]

    genome_len = parse_genome_len_from_gff3(args.gff3)
    if genome_len is None:
        raise SystemExit("Could not parse genome length from GFF3 header (##sequence-region).")

    masks, cds_list = parse_gff3_masks_and_cds(args.gff3, genome_len, feature_types)
    prefixes = {t: build_prefix_sum(m, genome_len) for t, m in masks.items()}

    hits = parse_blast_hits(args.blast)
    if not hits:
        raise SystemExit("No BLAST hits parsed (check that the file is outfmt 6/7 and tab-separated).")

    def summarize(ftype: str, thresh: float) -> Tuple[float, float, float, float]:
        any_q = 0
        any_s = 0
        high_q = 0
        high_s = 0
        for h in hits:
            qa, qb = h.q_interval
            sa, sb = h.s_interval
            qlen = qb - qa + 1
            slen = sb - sa + 1
            qcov = interval_cov(prefixes[ftype], qa, qb)
            scov = interval_cov(prefixes[ftype], sa, sb)
            if qcov > 0:
                any_q += 1
            if scov > 0:
                any_s += 1
            if qlen > 0 and (qcov / qlen) >= thresh:
                high_q += 1
            if slen > 0 and (scov / slen) >= thresh:
                high_s += 1
        n = float(len(hits))
        return any_q / n, any_s / n, high_q / n, high_s / n

    print(f"Genome length: {genome_len} bp")
    print(f"BLAST hits: {len(hits)}")
    print("---")
    for t in feature_types:
        any_q, any_s, high_q, high_s = summarize(t, args.cds_coverage_threshold)
        print(
            f"{t}: query any {any_q*100:.1f}% (>= {args.cds_coverage_threshold*100:.0f}% {high_q*100:.1f}%), "
            f"subject any {any_s*100:.1f}% (>= {args.cds_coverage_threshold*100:.0f}% {high_s*100:.1f}%)"
        )

    if "CDS" in prefixes:
        non_cds = 0
        for h in hits:
            qa, qb = h.q_interval
            qlen = qb - qa + 1
            qcov = interval_cov(prefixes["CDS"], qa, qb)
            if qlen > 0 and (qcov / qlen) < args.mostly_non_cds_threshold:
                non_cds += 1
        print("---")
        print(
            f"Hits with < {args.mostly_non_cds_threshold*100:.0f}% CDS coverage in query: "
            f"{non_cds} / {len(hits)} ({non_cds/len(hits)*100:.1f}%)"
        )

    # Top-N hits: enumerate overlapped CDS names for interpretability
    if cds_list:
        cds_bins = build_cds_bin_index(cds_list, args.bin_size)
        hits_sorted = sorted(hits, key=lambda h: h.alen_aa, reverse=True)
        print("---")
        print(f"Top {min(args.top, len(hits_sorted))} hits by alignment length:")
        for i, h in enumerate(hits_sorted[: args.top], 1):
            qa, qb = h.q_interval
            span_bp = qb - qa + 1
            cds_idx = overlapping_cds_indices(cds_list, cds_bins, args.bin_size, qa, qb)
            # Show up to 6 names to keep output readable
            shown = 6
            parts = []
            for idx in cds_idx[:shown]:
                feat = cds_list[idx]
                parts.append(f"{feat.display_name()}({feat.start}-{feat.end}{feat.strand})")
            more = "" if len(cds_idx) <= shown else f" ... +{len(cds_idx)-shown}"
            cds_str = ", ".join(parts) + more
            print(
                f"{i:2d}. alen={h.alen_aa} aa, q={h.q_start}-{h.q_end} span={span_bp} bp, CDS_count={len(cds_idx)}"
                + (f" :: {cds_str}" if cds_str.strip() else "")
            )


if __name__ == "__main__":
    main()



