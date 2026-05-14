#!/usr/bin/env python3
"""Classify NCBI-vs-LOSAT default outfmt 6 TBLASTX differences.

NCBI reference comments:
- c++/src/objtools/align_format/format_flags.cpp:39-40 defines default
  tabular fields as qaccver, saccver, pident, length, mismatch, gapopen,
  qstart, qend, sstart, send, evalue, bitscore.
- c++/src/algo/blast/core/blast_engine.c:1515-1541 calculates TBLASTX
  E-values and bit scores before output formatting.
- c++/src/algo/blast/core/blast_hits.c:1811-1905 assigns HSP E-values, and
  blast_hits.c:1907-1928 assigns bit scores. This diagnostic reports both the
  "remove only E-value" key and the coordinate/identity key that also removes
  bit-score text, so bit-score drift is not misclassified as HSP-set drift.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_FIELDS = (
    "qaccver",
    "saccver",
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
)
KEY_FIELDS = DEFAULT_FIELDS[:10] + (DEFAULT_FIELDS[11],)
COORD_KEY_FIELDS = DEFAULT_FIELDS[:10]


@dataclass(frozen=True)
class HspRow:
    line: str
    fields: tuple[str, ...]

    @property
    def key_without_evalue(self) -> tuple[str, ...]:
        return self.fields[:10] + (self.fields[11],)

    @property
    def coord_key_without_scores(self) -> tuple[str, ...]:
        return self.fields[:10]

    @property
    def evalue_text(self) -> str:
        return self.fields[10]

    @property
    def bitscore_text(self) -> str:
        return self.fields[11]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Classify default outfmt 6 differences by full row, HSP key, and E-value text."
    )
    parser.add_argument("--ncbi", required=True, type=Path, help="NCBI BLAST+ outfmt 6 output")
    parser.add_argument("--native", required=True, type=Path, help="LOSAT native outfmt 6 output")
    parser.add_argument("--out-dir", required=True, type=Path, help="Scratch output directory")
    return parser.parse_args()


def split_outfmt6_line(line: str, path: Path, line_no: int) -> tuple[str, ...]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) == 1:
        fields = line.split()
    if len(fields) != len(DEFAULT_FIELDS):
        raise ValueError(
            f"{path}:{line_no}: expected {len(DEFAULT_FIELDS)} default outfmt 6 fields, got {len(fields)}"
        )
    return tuple(fields)


def read_rows(path: Path) -> list[HspRow]:
    rows: list[HspRow] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = split_outfmt6_line(line, path, line_no)
            rows.append(HspRow(line=line.rstrip("\n"), fields=fields))
    return rows


def counter_items(counter: Counter[tuple[str, ...]]) -> Iterable[tuple[str, ...]]:
    for key in sorted(counter):
        for _ in range(counter[key]):
            yield key


def write_key_rows(path: Path, rows: Iterable[tuple[str, ...]]) -> int:
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(KEY_FIELDS) + "\n")
        for key in rows:
            handle.write("\t".join(key) + "\n")
            count += 1
    return count


def write_coord_key_rows(path: Path, rows: Iterable[tuple[str, ...]]) -> int:
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(COORD_KEY_FIELDS) + "\n")
        for key in rows:
            handle.write("\t".join(key) + "\n")
            count += 1
    return count


def parse_float(text: str) -> float | None:
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def bit_bin(key: tuple[str, ...]) -> int | None:
    value = parse_float(key[-1])
    if value is None:
        return None
    return math.floor(value)


def write_score_bins(path: Path, rows: Iterable[tuple[str, ...]]) -> tuple[int | None, int | None]:
    counts: Counter[int] = Counter()
    for key in rows:
        bin_value = bit_bin(key)
        if bin_value is not None:
            counts[bin_value] += 1

    with path.open("w", encoding="utf-8") as handle:
        handle.write("bit_score_floor\tcount\n")
        for bin_value in sorted(counts):
            handle.write(f"{bin_value}\t{counts[bin_value]}\n")

    if not counts:
        return None, None
    return min(counts), max(counts)


def write_row_score_bins(path: Path, rows: Iterable[HspRow]) -> tuple[int | None, int | None]:
    counts: Counter[int] = Counter()
    for row in rows:
        value = parse_float(row.bitscore_text)
        if value is not None:
            counts[math.floor(value)] += 1

    with path.open("w", encoding="utf-8") as handle:
        handle.write("bit_score_floor\tcount\n")
        for bin_value in sorted(counts):
            handle.write(f"{bin_value}\t{counts[bin_value]}\n")

    if not counts:
        return None, None
    return min(counts), max(counts)


def select_rows_for_coord_counter(
    rows: list[HspRow], key_counter: Counter[tuple[str, ...]]
) -> list[HspRow]:
    remaining = key_counter.copy()
    selected: list[HspRow] = []
    for row in rows:
        key = row.coord_key_without_scores
        if remaining[key] > 0:
            selected.append(row)
            remaining[key] -= 1
    return selected


def key_to_output(key: tuple[str, ...]) -> str:
    return "\t".join(key)


def write_shared_evalue_mismatches(
    path: Path,
    ncbi_by_key: dict[tuple[str, ...], Counter[str]],
    native_by_key: dict[tuple[str, ...], Counter[str]],
) -> int:
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                (
                    *KEY_FIELDS,
                    "ncbi_evalue",
                    "native_evalue",
                    "native_over_ncbi_ratio",
                )
            )
            + "\n"
        )
        for key in sorted(set(ncbi_by_key) & set(native_by_key)):
            ncbi_unmatched = ncbi_by_key[key] - native_by_key[key]
            native_unmatched = native_by_key[key] - ncbi_by_key[key]
            ncbi_values = sorted(ncbi_unmatched.elements())
            native_values = sorted(native_unmatched.elements())
            for idx in range(max(len(ncbi_values), len(native_values))):
                ncbi_evalue = ncbi_values[idx] if idx < len(ncbi_values) else ""
                native_evalue = native_values[idx] if idx < len(native_values) else ""
                ratio = ""
                ncbi_float = parse_float(ncbi_evalue) if ncbi_evalue else None
                native_float = parse_float(native_evalue) if native_evalue else None
                if (
                    ncbi_float is not None
                    and native_float is not None
                    and ncbi_float != 0.0
                    and native_float != 0.0
                ):
                    ratio = f"{native_float / ncbi_float:.17g}"
                handle.write(
                    "\t".join((key_to_output(key), ncbi_evalue, native_evalue, ratio)) + "\n"
                )
                count += 1
    return count


def write_shared_bitscore_mismatches(
    path: Path,
    ncbi_by_key: dict[tuple[str, ...], Counter[str]],
    native_by_key: dict[tuple[str, ...], Counter[str]],
) -> int:
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "\t".join((*COORD_KEY_FIELDS, "ncbi_bitscore", "native_bitscore"))
            + "\n"
        )
        for key in sorted(set(ncbi_by_key) & set(native_by_key)):
            ncbi_unmatched = ncbi_by_key[key] - native_by_key[key]
            native_unmatched = native_by_key[key] - ncbi_by_key[key]
            ncbi_values = sorted(ncbi_unmatched.elements())
            native_values = sorted(native_unmatched.elements())
            for idx in range(max(len(ncbi_values), len(native_values))):
                ncbi_bitscore = ncbi_values[idx] if idx < len(ncbi_values) else ""
                native_bitscore = native_values[idx] if idx < len(native_values) else ""
                handle.write(
                    "\t".join((key_to_output(key), ncbi_bitscore, native_bitscore)) + "\n"
                )
                count += 1
    return count


def write_sorted_keys(path: Path, rows: list[HspRow]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(row.key_without_evalue for row in rows):
            handle.write(key_to_output(key) + "\n")


def write_sorted_coord_keys(path: Path, rows: list[HspRow]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(row.coord_key_without_scores for row in rows):
            handle.write(key_to_output(key) + "\n")


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    ncbi_rows = read_rows(args.ncbi)
    native_rows = read_rows(args.native)

    ncbi_full = Counter(row.line for row in ncbi_rows)
    native_full = Counter(row.line for row in native_rows)
    ncbi_keys = Counter(row.key_without_evalue for row in ncbi_rows)
    native_keys = Counter(row.key_without_evalue for row in native_rows)
    ncbi_coord_keys = Counter(row.coord_key_without_scores for row in ncbi_rows)
    native_coord_keys = Counter(row.coord_key_without_scores for row in native_rows)

    shared_full = sum((ncbi_full & native_full).values())
    shared_keys = sum((ncbi_keys & native_keys).values())
    shared_coord_keys = sum((ncbi_coord_keys & native_coord_keys).values())

    ncbi_only_keys = ncbi_keys - native_keys
    native_only_keys = native_keys - ncbi_keys
    ncbi_only_coord_keys = ncbi_coord_keys - native_coord_keys
    native_only_coord_keys = native_coord_keys - ncbi_coord_keys

    ncbi_only_key_list = list(counter_items(ncbi_only_keys))
    native_only_key_list = list(counter_items(native_only_keys))
    ncbi_only_coord_key_list = list(counter_items(ncbi_only_coord_keys))
    native_only_coord_key_list = list(counter_items(native_only_coord_keys))
    ncbi_only_coord_rows = select_rows_for_coord_counter(ncbi_rows, ncbi_only_coord_keys)
    native_only_coord_rows = select_rows_for_coord_counter(native_rows, native_only_coord_keys)

    write_sorted_keys(args.out_dir / "ncbi.hsp_keys.no_evalue.sorted.tsv", ncbi_rows)
    write_sorted_keys(args.out_dir / "native.hsp_keys.no_evalue.sorted.tsv", native_rows)
    write_sorted_coord_keys(
        args.out_dir / "ncbi.alignment_keys.no_evalue_or_bitscore.sorted.tsv", ncbi_rows
    )
    write_sorted_coord_keys(
        args.out_dir / "native.alignment_keys.no_evalue_or_bitscore.sorted.tsv", native_rows
    )

    ncbi_only_count = write_key_rows(
        args.out_dir / "ncbi_only.hsp_keys.no_evalue.tsv", ncbi_only_key_list
    )
    native_only_count = write_key_rows(
        args.out_dir / "native_only.hsp_keys.no_evalue.tsv", native_only_key_list
    )
    ncbi_only_coord_count = write_coord_key_rows(
        args.out_dir / "ncbi_only.alignment_keys.no_evalue_or_bitscore.tsv",
        ncbi_only_coord_key_list,
    )
    native_only_coord_count = write_coord_key_rows(
        args.out_dir / "native_only.alignment_keys.no_evalue_or_bitscore.tsv",
        native_only_coord_key_list,
    )

    ncbi_min_bin, ncbi_max_bin = write_score_bins(
        args.out_dir / "ncbi_only.score_bins.tsv", ncbi_only_key_list
    )
    native_min_bin, native_max_bin = write_score_bins(
        args.out_dir / "native_only.score_bins.tsv", native_only_key_list
    )
    ncbi_coord_min_bin, ncbi_coord_max_bin = write_row_score_bins(
        args.out_dir / "ncbi_only.alignment_key_score_bins.tsv", ncbi_only_coord_rows
    )
    native_coord_min_bin, native_coord_max_bin = write_row_score_bins(
        args.out_dir / "native_only.alignment_key_score_bins.tsv", native_only_coord_rows
    )

    ncbi_by_key: dict[tuple[str, ...], Counter[str]] = defaultdict(Counter)
    native_by_key: dict[tuple[str, ...], Counter[str]] = defaultdict(Counter)
    ncbi_bits_by_coord_key: dict[tuple[str, ...], Counter[str]] = defaultdict(Counter)
    native_bits_by_coord_key: dict[tuple[str, ...], Counter[str]] = defaultdict(Counter)
    for row in ncbi_rows:
        ncbi_by_key[row.key_without_evalue][row.evalue_text] += 1
        ncbi_bits_by_coord_key[row.coord_key_without_scores][row.bitscore_text] += 1
    for row in native_rows:
        native_by_key[row.key_without_evalue][row.evalue_text] += 1
        native_bits_by_coord_key[row.coord_key_without_scores][row.bitscore_text] += 1

    evalue_mismatch_count = write_shared_evalue_mismatches(
        args.out_dir / "shared_hsp_evalue_mismatches.tsv",
        ncbi_by_key,
        native_by_key,
    )
    bitscore_mismatch_count = write_shared_bitscore_mismatches(
        args.out_dir / "shared_hsp_bitscore_mismatches.tsv",
        ncbi_bits_by_coord_key,
        native_bits_by_coord_key,
    )

    all_bins = [
        bin_value
        for bin_value in (ncbi_min_bin, ncbi_max_bin, native_min_bin, native_max_bin)
        if bin_value is not None
    ]
    lowest_diff_bin = min(all_bins) if all_bins else None
    highest_diff_bin = max(all_bins) if all_bins else None
    all_coord_bins = [
        bin_value
        for bin_value in (
            ncbi_coord_min_bin,
            ncbi_coord_max_bin,
            native_coord_min_bin,
            native_coord_max_bin,
        )
        if bin_value is not None
    ]
    lowest_coord_diff_bin = min(all_coord_bins) if all_coord_bins else None
    highest_coord_diff_bin = max(all_coord_bins) if all_coord_bins else None

    summary = {
        "ncbi_line_count": len(ncbi_rows),
        "native_line_count": len(native_rows),
        "exact_full_line_matches_after_sorting": shared_full,
        "matches_ignoring_only_evalue": shared_keys,
        "matches_ignoring_evalue_and_bitscore": shared_coord_keys,
        "shared_hsp_evalue_mismatch_rows": evalue_mismatch_count,
        "shared_hsp_bitscore_mismatch_rows": bitscore_mismatch_count,
        "ncbi_only_hsp_keys_ignoring_evalue": ncbi_only_count,
        "native_only_hsp_keys_ignoring_evalue": native_only_count,
        "ncbi_only_alignment_keys_ignoring_evalue_and_bitscore": ncbi_only_coord_count,
        "native_only_alignment_keys_ignoring_evalue_and_bitscore": native_only_coord_count,
        "hsp_key_diff_lowest_bit_score_floor": lowest_diff_bin,
        "hsp_key_diff_highest_bit_score_floor": highest_diff_bin,
        "first_descending_score_bin_with_hsp_key_diff": highest_diff_bin,
        "alignment_key_diff_lowest_bit_score_floor": lowest_coord_diff_bin,
        "alignment_key_diff_highest_bit_score_floor": highest_coord_diff_bin,
        "first_descending_score_bin_with_alignment_key_diff": highest_coord_diff_bin,
    }

    with (args.out_dir / "classification.summary.tsv").open("w", encoding="utf-8") as handle:
        for key, value in summary.items():
            handle.write(f"{key}\t{'' if value is None else value}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
