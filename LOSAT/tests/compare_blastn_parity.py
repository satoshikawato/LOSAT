#!/usr/bin/env python3
"""Compare LOSAT BLASTN tabular output against an NCBI BLASTN oracle.

This script intentionally treats NCBI output as the oracle and reports where
LOSAT diverges by final coordinate key and by same-coordinate values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import shlex
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-40
# ```c
# const char* kDfltArgTabularOutputFmt =
#     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
#     "evalue bitscore";
# ```
# NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
# ```c
# ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
#     if (iter != m_FieldsToShow.begin())
#         m_Ostream << m_FieldDelimiter;
#     x_PrintField(*iter);
# }
# m_Ostream << "\n";
# ```
DEFAULT_FIELDS = (
    "query",
    "subject",
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


@dataclass(frozen=True)
class HspRow:
    query: str
    subject: str
    pident: str
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: str
    bitscore: str
    ordinal: int

    @property
    def coord_key(self) -> tuple[str, str, int, int, int, int]:
        return (self.query, self.subject, self.qstart, self.qend, self.sstart, self.send)

    @property
    def value_key(self) -> tuple[str, int, int, int, str, str]:
        return (self.pident, self.length, self.mismatch, self.gapopen, self.evalue, self.bitscore)

    # NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
    # ```c
    # ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    #     if (iter != m_FieldsToShow.begin())
    #         m_Ostream << m_FieldDelimiter;
    #     x_PrintField(*iter);
    # }
    # m_Ostream << "\n";
    # ```
    # All emitted default fields participate in structured equality.  Keeping
    # this key separate from ordinal lets a Counter retain duplicate rows while
    # the ordered list comparison below retains their emission order.
    @property
    def structured_key(self) -> tuple[str, str, str, int, int, int, int, int, int, int, str, str]:
        return (
            self.query,
            self.subject,
            self.pident,
            self.length,
            self.mismatch,
            self.gapopen,
            self.qstart,
            self.qend,
            self.sstart,
            self.send,
            self.evalue,
            self.bitscore,
        )

    @property
    def bitscore_float(self) -> float:
        try:
            return float(self.bitscore)
        except ValueError:
            return float("-inf")

    @property
    def pident_float(self) -> float:
        try:
            return float(self.pident)
        except ValueError:
            return float("nan")

    @property
    def diagnostic_impact(self) -> float:
        bitscore = self.bitscore_float
        if math.isinf(bitscore) or math.isnan(bitscore):
            return 0.0
        return self.length * bitscore


def parse_outfmt(path: Path) -> list[HspRow]:
    rows: list[HspRow] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < len(DEFAULT_FIELDS):
                raise ValueError(f"{path}: expected at least 12 tabular fields, got {len(parts)}: {line!r}")
            rows.append(
                HspRow(
                    query=parts[0],
                    subject=parts[1],
                    pident=parts[2],
                    length=int(parts[3]),
                    mismatch=int(parts[4]),
                    gapopen=int(parts[5]),
                    qstart=int(parts[6]),
                    qend=int(parts[7]),
                    sstart=int(parts[8]),
                    send=int(parts[9]),
                    evalue=parts[10],
                    bitscore=parts[11],
                    ordinal=len(rows),
                )
            )
    return rows


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
# ```c
# m_Ostream << "\n";
# // The tabular writer emits the selected fields and delimiters directly to
# // the output stream; byte comparison is therefore the final output gate.
# ```
def compare_raw_bytes(ncbi_path: Path, losat_path: Path) -> bool:
    ncbi_bytes = ncbi_path.read_bytes()
    losat_bytes = losat_path.read_bytes()
    if ncbi_bytes == losat_bytes:
        print(f"  raw_bytes: exact ({len(ncbi_bytes)} bytes)")
        print(f"  raw_sha256_ncbi: {hashlib.sha256(ncbi_bytes).hexdigest()}")
        print(f"  raw_sha256_losat: {hashlib.sha256(losat_bytes).hexdigest()}")
        return True

    common_len = min(len(ncbi_bytes), len(losat_bytes))
    first_diff = next(
        (index for index in range(common_len) if ncbi_bytes[index] != losat_bytes[index]),
        common_len,
    )
    print(
        "  raw_bytes: DIFFER"
        f" (ncbi={len(ncbi_bytes)} bytes, losat={len(losat_bytes)} bytes,"
        f" first_diff={first_diff})"
    )
    print(f"  raw_sha256_ncbi: {hashlib.sha256(ncbi_bytes).hexdigest()}")
    print(f"  raw_sha256_losat: {hashlib.sha256(losat_bytes).hexdigest()}")
    return False


def rows_by_coord(rows: Iterable[HspRow]) -> dict[tuple[str, str, int, int, int, int], HspRow]:
    by_key: dict[tuple[str, str, int, int, int, int], HspRow] = {}
    duplicates: Counter[tuple[str, str, int, int, int, int]] = Counter()
    for row in rows:
        if row.coord_key in by_key:
            duplicates[row.coord_key] += 1
            continue
        by_key[row.coord_key] = row
    if duplicates:
        dup_count = sum(duplicates.values())
        print(
            f"warning: coordinate-key diagnostics collapse {dup_count} duplicate rows; "
            "structured equality still compares them",
            file=sys.stderr,
        )
    return by_key


# NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:811-832
# ```c
# if (results.HasAlignments()) {
#     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_HitlistSize);
#     ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
#         const CSeq_align& s = **itr;
#         tabinfo.SetFields(s, *m_Scope, &m_ScoringMatrix);
#         tabinfo.Print();
#     }
# }
# ```
# NCBI emits every retained alignment in CSeq_align_set iteration order.  These
# counters expose multiplicity without replacing the existing coordinate-key
# diagnostics, which intentionally retain only the first row for each key.
def duplicate_coord_rows(rows: Iterable[HspRow]) -> dict[tuple[str, str, int, int, int, int], int]:
    counts = Counter(row.coord_key for row in rows)
    return {key: count for key, count in counts.items() if count > 1}


# NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:828-832
# ```c
# ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
#     const CSeq_align& s = **itr;
#     tabinfo.SetFields(s, *m_Scope, &m_ScoringMatrix);
#     tabinfo.Print();
# }
# ```
def first_order_difference(ncbi_rows: list[HspRow], losat_rows: list[HspRow]) -> int | None:
    for index, (ncbi_row, losat_row) in enumerate(zip(ncbi_rows, losat_rows)):
        if ncbi_row.structured_key != losat_row.structured_key:
            return index
    if len(ncbi_rows) != len(losat_rows):
        return min(len(ncbi_rows), len(losat_rows))
    return None


def score_bins(rows: Iterable[HspRow]) -> Counter[str]:
    bins: Counter[str] = Counter()
    for row in rows:
        bits = row.bitscore_float
        if math.isinf(bits):
            bins["nan"] += 1
        else:
            start = int(bits // 10) * 10
            bins[f"{start:03d}-{start + 9:03d}"] += 1
    return bins


def length_bins(rows: Iterable[HspRow]) -> Counter[str]:
    bins: Counter[str] = Counter()
    for row in rows:
        start = int(row.length // 25) * 25
        bins[f"{start:04d}-{start + 24:04d}"] += 1
    return bins


def identity_bins(rows: Iterable[HspRow]) -> Counter[str]:
    bins: Counter[str] = Counter()
    for row in rows:
        pident = row.pident_float
        if math.isnan(pident):
            bins["nan"] += 1
            continue
        start = int(pident // 5) * 5
        end = min(start + 4, 100)
        bins[f"{start:03d}-{end:03d}"] += 1
    return bins


def sorted_counter(counter: Counter[str]) -> dict[str, int]:
    return {key: counter[key] for key in sorted(counter)}


def first_n(rows: list[HspRow], keys: set[tuple[str, str, int, int, int, int]], limit: int) -> list[HspRow]:
    return [row for row in rows if row.coord_key in keys][:limit]


def rows_for_keys(rows: Iterable[HspRow], keys: set[tuple[str, str, int, int, int, int]]) -> list[HspRow]:
    return [row for row in rows if row.coord_key in keys]


def ranked_by_impact(
    rows: Iterable[HspRow],
    keys: set[tuple[str, str, int, int, int, int]],
    limit: int,
) -> list[HspRow]:
    return sorted(
        rows_for_keys(rows, keys),
        key=lambda row: (row.diagnostic_impact, row.bitscore_float, row.length, -row.ordinal),
        reverse=True,
    )[:limit]


def total_alignment_length(rows: Iterable[HspRow]) -> int:
    return sum(row.length for row in rows)


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
# ```c
# ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
#     if (iter != m_FieldsToShow.begin())
#         m_Ostream << m_FieldDelimiter;
#     x_PrintField(*iter);
# }
# m_Ostream << "\n";
# ```
def compare_case(case_id: str, ncbi_path: Path, losat_path: Path, limit: int) -> dict[str, object]:
    ncbi_rows = parse_outfmt(ncbi_path)
    losat_rows = parse_outfmt(losat_path)
    ncbi_structured_rows = [row.structured_key for row in ncbi_rows]
    losat_structured_rows = [row.structured_key for row in losat_rows]
    row_count_equal = len(ncbi_rows) == len(losat_rows)
    row_multiset_equal = Counter(ncbi_structured_rows) == Counter(losat_structured_rows)
    row_order_equal = ncbi_structured_rows == losat_structured_rows
    structured_equal = row_count_equal and row_multiset_equal and row_order_equal
    order_difference = first_order_difference(ncbi_rows, losat_rows)
    ncbi_duplicates = duplicate_coord_rows(ncbi_rows)
    losat_duplicates = duplicate_coord_rows(losat_rows)
    ncbi_by_key = rows_by_coord(ncbi_rows)
    losat_by_key = rows_by_coord(losat_rows)
    ncbi_keys = set(ncbi_by_key)
    losat_keys = set(losat_by_key)
    common = ncbi_keys & losat_keys
    ncbi_only = ncbi_keys - losat_keys
    losat_only = losat_keys - ncbi_keys
    ncbi_only_rows = rows_for_keys(ncbi_rows, ncbi_only)
    losat_only_rows = rows_for_keys(losat_rows, losat_only)

    bitscore_diffs = []
    evalue_diffs = []
    pident_diffs = []
    other_value_diffs = []
    for key in sorted(common):
        ncbi = ncbi_by_key[key]
        losat = losat_by_key[key]
        if ncbi.bitscore != losat.bitscore:
            bitscore_diffs.append((key, ncbi, losat))
        if ncbi.evalue != losat.evalue:
            evalue_diffs.append((key, ncbi, losat))
        if ncbi.pident != losat.pident:
            pident_diffs.append((key, ncbi, losat))
        if (
            ncbi.length,
            ncbi.mismatch,
            ncbi.gapopen,
        ) != (
            losat.length,
            losat.mismatch,
            losat.gapopen,
        ):
            other_value_diffs.append((key, ncbi, losat))

    print(f"case: {case_id}")
    print(f"  ncbi_hits: {len(ncbi_rows)}")
    print(f"  losat_hits: {len(losat_rows)}")
    print(f"  structured_row_count_equal: {row_count_equal}")
    print(f"  structured_duplicate_rows_ncbi: {sum(count - 1 for count in ncbi_duplicates.values())}")
    print(f"  structured_duplicate_rows_losat: {sum(count - 1 for count in losat_duplicates.values())}")
    print(f"  structured_row_multiset_equal: {row_multiset_equal}")
    print(f"  structured_row_order_equal: {row_order_equal}")
    print(f"  structured_equal: {structured_equal}")
    print(f"  common_coord_keys: {len(common)}")
    print(f"  ncbi_only: {len(ncbi_only)}")
    print(f"  losat_only: {len(losat_only)}")
    print(f"  ncbi_only_accumulated_length: {total_alignment_length(ncbi_only_rows)}")
    print(f"  losat_only_accumulated_length: {total_alignment_length(losat_only_rows)}")
    print(f"  same_coordinate_bitscore_diffs: {len(bitscore_diffs)}")
    print(f"  same_coordinate_evalue_diffs: {len(evalue_diffs)}")
    print(f"  same_coordinate_pident_diffs: {len(pident_diffs)}")
    print(f"  same_coordinate_length_mismatch_gapopen_diffs: {len(other_value_diffs)}")

    if order_difference is not None:
        print(f"  first_structured_order_diff: {order_difference}")
        if order_difference < len(ncbi_rows):
            print(format_row("    ncbi=", ncbi_rows[order_difference]))
        else:
            print("    ncbi=<missing>")
        if order_difference < len(losat_rows):
            print(format_row("    losat=", losat_rows[order_difference]))
        else:
            print("    losat=<missing>")

    if ncbi_only:
        print("  first_ncbi_only:")
        for row in first_n(ncbi_rows, ncbi_only, limit):
            print(format_row("    ", row))
    if losat_only:
        print("  first_losat_only:")
        for row in first_n(losat_rows, losat_only, limit):
            print(format_row("    ", row))
    if ncbi_only:
        print("  top_ncbi_only_by_length_bitscore:")
        for row in ranked_by_impact(ncbi_rows, ncbi_only, limit):
            print(format_row("    ", row))
    if losat_only:
        print("  top_losat_only_by_length_bitscore:")
        for row in ranked_by_impact(losat_rows, losat_only, limit):
            print(format_row("    ", row))
    if bitscore_diffs:
        print("  first_bitscore_diffs:")
        for key, ncbi, losat in bitscore_diffs[:limit]:
            print(f"    key={key}")
            print(f"      ncbi={ncbi.bitscore}")
            print(f"      losat={losat.bitscore}")
    if evalue_diffs:
        print("  first_evalue_diffs:")
        for key, ncbi, losat in evalue_diffs[:limit]:
            print(f"    key={key}")
            print(f"      ncbi={ncbi.evalue}")
            print(f"      losat={losat.evalue}")
    if pident_diffs:
        print("  first_pident_diffs:")
        for key, ncbi, losat in pident_diffs[:limit]:
            print(f"    key={key}")
            print(f"      ncbi={ncbi.pident}")
            print(f"      losat={losat.pident}")
    if other_value_diffs:
        print("  first_length_mismatch_gapopen_diffs:")
        for key, ncbi, losat in other_value_diffs[:limit]:
            print(f"    key={key}")
            print(f"      ncbi={(ncbi.length, ncbi.mismatch, ncbi.gapopen)}")
            print(f"      losat={(losat.length, losat.mismatch, losat.gapopen)}")

    print("  score_bins_ncbi_only:", sorted_counter(score_bins(ncbi_only_rows)))
    print("  score_bins_losat_only:", sorted_counter(score_bins(losat_only_rows)))
    print("  identity_bins_ncbi_only:", sorted_counter(identity_bins(ncbi_only_rows)))
    print("  identity_bins_losat_only:", sorted_counter(identity_bins(losat_only_rows)))
    print("  length_bins_ncbi_only:", sorted_counter(length_bins(ncbi_only_rows)))
    print("  length_bins_losat_only:", sorted_counter(length_bins(losat_only_rows)))

    return {
        "case_id": case_id,
        "ncbi_rows": ncbi_rows,
        "losat_rows": losat_rows,
        "row_count_equal": row_count_equal,
        "row_multiset_equal": row_multiset_equal,
        "row_order_equal": row_order_equal,
        "structured_equal": structured_equal,
        "ncbi_duplicates": ncbi_duplicates,
        "losat_duplicates": losat_duplicates,
        "ncbi_by_key": ncbi_by_key,
        "losat_by_key": losat_by_key,
        "ncbi_only": ncbi_only,
        "losat_only": losat_only,
        "bitscore_diffs": bitscore_diffs,
        "evalue_diffs": evalue_diffs,
        "pident_diffs": pident_diffs,
        "other_value_diffs": other_value_diffs,
    }


def format_row(prefix: str, row: HspRow) -> str:
    return (
        f"{prefix}{row.query}\t{row.subject}\tq={row.qstart}-{row.qend}\t"
        f"s={row.sstart}-{row.send}\tpident={row.pident}\tlen={row.length}\t"
        f"mismatch={row.mismatch}\tgapopen={row.gapopen}\te={row.evalue}\tbit={row.bitscore}"
    )


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#") and line.strip()]
    return list(csv.DictReader(lines, delimiter="\t"))


def manifest_bool(row: dict[str, str], field: str) -> bool:
    value = row[field].strip().lower()
    if value in {"1", "true", "yes", "on"}:
        return True
    if value in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"{row['case_id']}: invalid boolean {field}={row[field]!r}")


def require_manifest_condition(row: dict[str, str], condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(f"{row['case_id']}: {message}")


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:365-368
# ```c
# arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
#                 "Apply filtering locations as soft masks",
#                 CArgDescriptions::eBoolean,
#                 kDfltArgLookupTableMaskingOnlyNucl);
# ```
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:475-480,505-511
# ```c
# arg_desc.AddOptionalKey(kArgWindowSize, "int_value",
#                         "Multiple hits window size, use 0 to specify 1-hit algorithm",
#                         CArgDescriptions::eInteger);
# arg_desc.AddDefaultKey(kArgOffDiagonalRange, "int_value",
#                         "Number of off-diagonals to search for the 2nd hit, use 0 to turn off",
#                         CArgDescriptions::eInteger, NStr::IntToString(kDfltOffDiagonalRange));
# ```
# LOSAT does not expose every one of these NCBI switches yet.  A paired oracle
# run is valid only for the manifest values matching the supported Rust path.
def validate_manifest_row(row: dict[str, str]) -> None:
    require_manifest_condition(
        row,
        manifest_bool(row, "dust"),
        "paired/refresh execution currently supports only dust=true",
    )
    require_manifest_condition(
        row,
        manifest_bool(row, "soft_masking"),
        "paired/refresh execution currently supports only soft_masking=true",
    )
    require_manifest_condition(
        row,
        row["strand"] == "both",
        "paired/refresh execution currently supports only strand=both",
    )
    require_manifest_condition(
        row,
        int(row["window_size"]) == 0,
        "paired/refresh execution currently supports only manifest window_size=0",
    )
    require_manifest_condition(
        row,
        int(row["scan_range"]) == 0,
        "paired/refresh execution currently supports only manifest scan_range=0",
    )
    require_manifest_condition(
        row,
        int(row["num_threads"]) == 1,
        "native serial parity execution requires num_threads=1",
    )
    require_manifest_condition(row, row["outfmt"] in {"6", "7"}, "only outfmt 6/7 are parseable")


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-107,157-168,208
# ```c
# const string kArgQuery("query");
# const string kArgOutput("out");
# const string kArgSubject("subject");
# const string kTask("task");
# const string kArgNumThreads("num_threads");
# const string kArgEvalue("evalue");
# const string kArgMaxTargetSequences("max_target_seqs");
# const string kArgGapOpen("gapopen");
# const string kArgGapExtend("gapextend");
# const string kArgMismatch("penalty");
# const string kArgMatch("reward");
# const string kArgWordSize("word_size");
# const string kArgDustFiltering("dust");
# const string kArgLookupTableMaskingOnly("soft_masking");
# const string kArgMaxHSPsPerSubject("max_hsps");
# ```
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:317-319,2960-2968
# ```c
# if (args.Exist(kArgMaxHSPsPerSubject) && args[kArgMaxHSPsPerSubject]) {
#     opt.SetMaxHspsPerSubject(args[kArgMaxHSPsPerSubject].AsInteger());
# }
# if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
#     hitlist_size = args[kArgMaxTargetSequences].AsInteger();
# }
# m_NumDescriptions = hitlist_size;
# m_NumAlignments = hitlist_size;
# ```
def build_losat_command(
    row: dict[str, str],
    losat_bin: Path,
    output_path: Path | str | None = None,
) -> list[str]:
    validate_manifest_row(row)
    output = str(output_path) if output_path is not None else row["losat_out"]

    command = [
        str(losat_bin),
        "blastn",
        "-q",
        row["query"],
        "-s",
        row["subject"],
        "--task",
        row["task"],
        "--reward",
        row["reward"],
        f"--penalty={row['penalty']}",
        "--gap-open",
        row["gap_open"],
        "--gap-extend",
        row["gap_extend"],
        "--word-size",
        row["word_size"],
        "-n",
        row["num_threads"],
        "--evalue",
        row["evalue"],
        "--max-target-seqs",
        row["max_target_seqs"],
        "--max-hsps-per-subject",
        row["max_hsps_per_subject"],
        "--outfmt",
        row["outfmt"],
        "-o",
        output,
    ]
    # LOSAT exposes dust as a clap flag; its default is the NCBI nucleotide
    # default (enabled), so only the enabled manifest value is representable.
    if manifest_bool(row, "dust"):
        command.insert(command.index("--outfmt"), "--dust")
    return command


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-107,157-168,208
# ```c
# const string kArgQuery("query");
# const string kArgOutput("out");
# const string kArgSubject("subject");
# const string kTask("task");
# const string kArgNumThreads("num_threads");
# const string kArgEvalue("evalue");
# const string kArgWordSize("word_size");
# const string kArgDustFiltering("dust");
# const string kArgLookupTableMaskingOnly("soft_masking");
# const string kArgMaxHSPsPerSubject("max_hsps");
# ```
def build_ncbi_command(row: dict[str, str], ncbi_bin: str, output_path: Path | str) -> list[str]:
    validate_manifest_row(row)
    command = [
        ncbi_bin,
        "-query",
        row["query"],
        "-subject",
        row["subject"],
        "-task",
        row["task"],
        "-reward",
        row["reward"],
        "-penalty",
        row["penalty"],
        "-gapopen",
        row["gap_open"],
        "-gapextend",
        row["gap_extend"],
        "-word_size",
        row["word_size"],
        "-dust",
        "yes" if manifest_bool(row, "dust") else "no",
        "-soft_masking",
        "true" if manifest_bool(row, "soft_masking") else "false",
        "-strand",
        row["strand"],
        "-evalue",
        row["evalue"],
        "-max_target_seqs",
        row["max_target_seqs"],
        "-num_threads",
        row["num_threads"],
        "-outfmt",
        row["outfmt"],
        "-out",
        str(output_path),
    ]
    # NCBI's CLI declares max_hsps as an optional integer constrained to >= 1;
    # manifest zero represents the unset/unlimited default and must be omitted.
    if int(row["max_hsps_per_subject"]) > 0:
        command.extend(["-max_hsps", row["max_hsps_per_subject"]])
    return command


def refresh_losat_output(root: Path, row: dict[str, str], losat_bin: Path, print_command: bool) -> None:
    losat_out = root / row["losat_out"]
    losat_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = build_losat_command(row, losat_bin)
    if print_command:
        print("refresh_losat:", shlex.join(cmd))
    subprocess.run(cmd, cwd=root, check=True)


# NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-782,828-832
# ```c
# if (m_FormatType == CFormattingArgs::eTabular ||
#     m_FormatType == CFormattingArgs::eTabularWithComments) {
#     CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
#     ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
#         tabinfo.SetFields(**itr, *m_Scope, &m_ScoringMatrix);
#         tabinfo.Print();
#     }
# }
# ```
# Fresh paired artifacts isolate both writers from stale checked-in output.  The
# NCBI executable is an explicitly supplied/local comparison oracle only; this
# helper does not locate, download, install, or link it into LOSAT.
def prepare_paired_output_dir(root: Path, requested: Path | None) -> Path:
    if requested is None:
        return Path(tempfile.mkdtemp(prefix="losat-blastn-parity-"))

    output_dir = requested if requested.is_absolute() else root / requested
    if output_dir.exists():
        if not output_dir.is_dir():
            raise ValueError(f"paired output path is not a directory: {output_dir}")
        if any(output_dir.iterdir()):
            raise ValueError(f"paired output directory must be empty: {output_dir}")
    else:
        output_dir.mkdir(parents=True)
    return output_dir


# NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-782,828-832
# ```c
# if (m_FormatType == CFormattingArgs::eTabular ||
#     m_FormatType == CFormattingArgs::eTabularWithComments) {
#     CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
#     ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
#         tabinfo.SetFields(**itr, *m_Scope, &m_ScoringMatrix);
#         tabinfo.Print();
#     }
# }
# ```
def run_fresh_paired(
    root: Path,
    rows: list[dict[str, str]],
    ncbi_bin: str,
    losat_bin: Path,
    requested_output_dir: Path | None,
) -> tuple[list[tuple[str, Path, Path]], Path]:
    output_dir = prepare_paired_output_dir(root, requested_output_dir)
    print(f"paired_cwd: {root}")
    print(f"paired_output_dir: {output_dir}")

    cases: list[tuple[str, Path, Path]] = []
    for row in rows:
        case_id = row["case_id"]
        ncbi_out = output_dir / f"{case_id}.ncbi.out"
        losat_out = output_dir / f"{case_id}.losat.out"
        ncbi_command = build_ncbi_command(row, ncbi_bin, ncbi_out)
        losat_command = build_losat_command(row, losat_bin, losat_out)
        print("paired_ncbi:", shlex.join(ncbi_command))
        subprocess.run(ncbi_command, cwd=root, check=True)
        print("paired_losat:", shlex.join(losat_command))
        subprocess.run(losat_command, cwd=root, check=True)
        cases.append((case_id, ncbi_out, losat_out))
    return cases, output_dir


def select_targets(result: dict[str, object]) -> list[tuple[str, HspRow]]:
    ncbi_rows: list[HspRow] = result["ncbi_rows"]  # type: ignore[assignment]
    losat_rows: list[HspRow] = result["losat_rows"]  # type: ignore[assignment]
    ncbi_only: set[tuple[str, str, int, int, int, int]] = result["ncbi_only"]  # type: ignore[assignment]
    losat_only: set[tuple[str, str, int, int, int, int]] = result["losat_only"]  # type: ignore[assignment]
    targets: list[tuple[str, HspRow]] = []

    ncbi_ranked = ranked_by_impact(ncbi_rows, ncbi_only, 1)
    if ncbi_ranked:
        targets.append(("ncbi_only_top_length_bitscore", ncbi_ranked[0]))

    losat_ranked = ranked_by_impact(losat_rows, losat_only, 1)
    if losat_ranked:
        targets.append(("losat_only_top_length_bitscore", losat_ranked[0]))

    losat_perfect = [row for row in losat_rows if row.coord_key in losat_only and row.pident_float == 100.0]
    if losat_perfect:
        targets.append(("losat_only_low_score_perfect", min(losat_perfect, key=lambda row: row.bitscore_float)))

    return targets


def write_targets(path: Path, results: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "case_id",
                "target_kind",
                "query_id",
                "subject_id",
                "qstart",
                "qend",
                "sstart",
                "send",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "evalue",
                "bitscore",
                "trace_env",
            ]
        )
        for result in results:
            case_id = str(result["case_id"])
            for kind, row in select_targets(result):
                writer.writerow(
                    [
                        case_id,
                        kind,
                        row.query,
                        row.subject,
                        row.qstart,
                        row.qend,
                        row.sstart,
                        row.send,
                        row.pident,
                        row.length,
                        row.mismatch,
                        row.gapopen,
                        row.evalue,
                        row.bitscore,
                        f"LOSAT_TRACE_BLASTN_HSP={row.qstart},{row.qend},{row.sstart},{row.send}",
                    ]
                )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("tests/blastn_parity_manifest.tsv"))
    parser.add_argument("--case-id", action="append", help="Run only the named manifest case. May be repeated.")
    parser.add_argument("--ncbi", type=Path, help="Direct NCBI outfmt 6/7 path.")
    parser.add_argument("--losat", type=Path, help="Direct LOSAT outfmt 6/7 path.")
    parser.add_argument("--limit", type=int, default=5)
    parser.add_argument("--write-targets", type=Path)
    parser.add_argument(
        "--refresh-losat",
        action="store_true",
        help="Regenerate selected LOSAT output paths from the manifest before comparison.",
    )
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-51
    # ```c
    # const string kArgQuery("query");
    # const string kArgOutput("out");
    # const string kArgSubject("subject");
    # ```
    # Paired mode invokes an already-available blastn only as a test oracle and
    # gives both processes the same manifest paths and cwd.
    parser.add_argument(
        "--fresh-paired",
        action="store_true",
        help=(
            "Run the selected manifest cases through the local NCBI blastn oracle and LOSAT "
            "into a fresh output directory before comparison. No executable is installed or downloaded."
        ),
    )
    parser.add_argument(
        "--paired-output-dir",
        type=Path,
        help="New or empty output directory for --fresh-paired (default: a new persistent temp directory).",
    )
    parser.add_argument(
        "--ncbi-bin",
        default="blastn",
        help="Existing local NCBI blastn oracle executable used only by --fresh-paired.",
    )
    parser.add_argument(
        "--losat-bin",
        type=Path,
        default=Path("target/release/LOSAT"),
        help="LOSAT executable for --refresh-losat/--fresh-paired, relative to cwd unless absolute.",
    )
    parser.add_argument("--print-refresh-commands", action="store_true")
    parser.add_argument("--fail-on-diff", action="store_true")
    parser.add_argument(
        "--fail-on-byte-diff",
        action="store_true",
        help="Compare complete NCBI/LOSAT files byte-for-byte and fail on any difference.",
    )
    args = parser.parse_args()

    root = Path.cwd()
    cases: list[tuple[str, Path, Path]] = []
    refresh_rows: list[dict[str, str]] = []
    selected_rows: list[dict[str, str]] = []
    if args.ncbi or args.losat:
        if not (args.ncbi and args.losat):
            parser.error("--ncbi and --losat must be provided together")
        if args.refresh_losat:
            parser.error("--refresh-losat requires manifest-selected cases, not direct --ncbi/--losat paths")
        if args.fresh_paired:
            parser.error("--fresh-paired requires manifest-selected cases, not direct --ncbi/--losat paths")
        cases.append(("direct", args.ncbi, args.losat))
    else:
        wanted = set(args.case_id or [])
        for row in read_manifest(args.manifest):
            case_id = row["case_id"]
            if wanted and case_id not in wanted:
                continue
            selected_rows.append(row)

        if not selected_rows:
            parser.error("no cases selected")

        if args.fresh_paired:
            if args.refresh_losat:
                parser.error("--fresh-paired and --refresh-losat are mutually exclusive")
            losat_bin = args.losat_bin if args.losat_bin.is_absolute() else root / args.losat_bin
            try:
                cases, _ = run_fresh_paired(
                    root,
                    selected_rows,
                    args.ncbi_bin,
                    losat_bin,
                    args.paired_output_dir,
                )
            except ValueError as error:
                parser.error(str(error))
        else:
            fresh_only: list[str] = []
            for row in selected_rows:
                ncbi_output = row["ncbi_out"].strip()
                losat_output = row["losat_out"].strip()
                if ncbi_output in {"", "-"} or losat_output in {"", "-"}:
                    require_manifest_condition(
                        row,
                        ncbi_output in {"", "-"} and losat_output in {"", "-"},
                        "ncbi_out and losat_out must both be artifact paths or both be '-'",
                    )
                    fresh_only.append(row["case_id"])
                    continue
                refresh_rows.append(row)
                cases.append((row["case_id"], root / ncbi_output, root / losat_output))

            if fresh_only:
                names = ", ".join(fresh_only)
                if wanted:
                    parser.error(f"selected cases require --fresh-paired: {names}")
                print(
                    f"note: skipped fresh-paired-only manifest cases: {names}; use --fresh-paired to execute them",
                    file=sys.stderr,
                )

    if not cases:
        parser.error("no cases selected")

    if args.paired_output_dir and not args.fresh_paired:
        parser.error("--paired-output-dir requires --fresh-paired")

    if args.refresh_losat:
        losat_bin = args.losat_bin if args.losat_bin.is_absolute() else root / args.losat_bin
        for row in refresh_rows:
            refresh_losat_output(root, row, losat_bin, args.print_refresh_commands)

    results = []
    raw_failed = False
    structured_failed = False
    for case_id, ncbi_path, losat_path in cases:
        raw_equal = compare_raw_bytes(ncbi_path, losat_path) if args.fail_on_byte_diff or args.fresh_paired else True
        result = compare_case(case_id, ncbi_path, losat_path, args.limit)
        results.append(result)
        raw_failed |= not raw_equal
        structured_failed |= not bool(result["structured_equal"])

    if args.write_targets:
        write_targets(args.write_targets, results)
        print(f"wrote trace targets: {args.write_targets}")

    # NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
    # ```c
    # ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    #     x_PrintField(*iter);
    # }
    # m_Ostream << "\n";
    # ```
    # Raw and structured gates are independent: requesting either one makes its
    # own observed divergence determine the process status.
    return 1 if (args.fail_on_byte_diff and raw_failed) or (args.fail_on_diff and structured_failed) else 0


if __name__ == "__main__":
    raise SystemExit(main())
