#!/usr/bin/env python3
"""Compare LOSAT BLASTN tabular output against an NCBI BLASTN oracle.

This script intentionally treats NCBI output as the oracle and reports where
LOSAT diverges by final coordinate key and by same-coordinate values.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
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
        print(f"warning: ignored {dup_count} duplicate coordinate-key rows", file=sys.stderr)
    return by_key


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


def first_n(rows: list[HspRow], keys: set[tuple[str, str, int, int, int, int]], limit: int) -> list[HspRow]:
    return [row for row in rows if row.coord_key in keys][:limit]


def compare_case(case_id: str, ncbi_path: Path, losat_path: Path, limit: int) -> dict[str, object]:
    ncbi_rows = parse_outfmt(ncbi_path)
    losat_rows = parse_outfmt(losat_path)
    ncbi_by_key = rows_by_coord(ncbi_rows)
    losat_by_key = rows_by_coord(losat_rows)
    ncbi_keys = set(ncbi_by_key)
    losat_keys = set(losat_by_key)
    common = ncbi_keys & losat_keys
    ncbi_only = ncbi_keys - losat_keys
    losat_only = losat_keys - ncbi_keys

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
    print(f"  common_coord_keys: {len(common)}")
    print(f"  ncbi_only: {len(ncbi_only)}")
    print(f"  losat_only: {len(losat_only)}")
    print(f"  same_coordinate_bitscore_diffs: {len(bitscore_diffs)}")
    print(f"  same_coordinate_evalue_diffs: {len(evalue_diffs)}")
    print(f"  same_coordinate_pident_diffs: {len(pident_diffs)}")
    print(f"  same_coordinate_length_mismatch_gapopen_diffs: {len(other_value_diffs)}")

    if ncbi_only:
        print("  first_ncbi_only:")
        for row in first_n(ncbi_rows, ncbi_only, limit):
            print(format_row("    ", row))
    if losat_only:
        print("  first_losat_only:")
        for row in first_n(losat_rows, losat_only, limit):
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

    print("  score_bins_ncbi_only:", dict(score_bins(ncbi_by_key[key] for key in ncbi_only)))
    print("  score_bins_losat_only:", dict(score_bins(losat_by_key[key] for key in losat_only)))
    print("  length_bins_ncbi_only:", dict(length_bins(ncbi_by_key[key] for key in ncbi_only)))
    print("  length_bins_losat_only:", dict(length_bins(losat_by_key[key] for key in losat_only)))

    return {
        "case_id": case_id,
        "ncbi_rows": ncbi_rows,
        "losat_rows": losat_rows,
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


def build_losat_command(row: dict[str, str], losat_bin: Path) -> list[str]:
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-97
    # ```c
    # const string kArgQuery("query");
    # const string kArgSubject("subject");
    # const string kTask("task");
    # const string kArgNumThreads("num_threads");
    # const string kArgEvalue("evalue");
    # const string kArgMaxTargetSequences("max_target_seqs");
    # const string kArgGapOpen("gapopen");
    # const string kArgGapExtend("gapextend");
    # const string kArgMismatch("penalty");
    # const string kArgMatch("reward");
    # ```
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:203-207
    # ```c
    # arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value",
    #                    "Set maximum number of HSPs per subject sequence to save for each query",
    #                    CArgDescriptions::eInteger);
    # ```
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2960-2968
    # ```c
    # if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    #     hitlist_size = args[kArgMaxTargetSequences].AsInteger();
    # }
    # m_NumDescriptions = hitlist_size;
    # m_NumAlignments = hitlist_size;
    # ```
    require_manifest_condition(row, manifest_bool(row, "dust"), "LOSAT CLI refresh currently supports only dust=true")
    require_manifest_condition(
        row,
        manifest_bool(row, "soft_masking"),
        "LOSAT CLI refresh currently supports only soft_masking=true",
    )
    require_manifest_condition(row, row["strand"] == "both", "LOSAT CLI refresh currently supports only strand=both")
    require_manifest_condition(
        row,
        int(row["window_size"]) == 0,
        "LOSAT CLI refresh currently supports only manifest window_size=0",
    )
    require_manifest_condition(
        row,
        int(row["scan_range"]) == 0,
        "LOSAT CLI refresh currently supports only manifest scan_range=0",
    )

    return [
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
        row["losat_out"],
    ]


def refresh_losat_output(root: Path, row: dict[str, str], losat_bin: Path, print_command: bool) -> None:
    losat_out = root / row["losat_out"]
    losat_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = build_losat_command(row, losat_bin)
    if print_command:
        print("refresh_losat:", " ".join(cmd))
    subprocess.run(cmd, cwd=root, check=True)


def select_targets(result: dict[str, object]) -> list[tuple[str, HspRow]]:
    ncbi_rows: list[HspRow] = result["ncbi_rows"]  # type: ignore[assignment]
    losat_rows: list[HspRow] = result["losat_rows"]  # type: ignore[assignment]
    ncbi_only: set[tuple[str, str, int, int, int, int]] = result["ncbi_only"]  # type: ignore[assignment]
    losat_only: set[tuple[str, str, int, int, int, int]] = result["losat_only"]  # type: ignore[assignment]
    targets: list[tuple[str, HspRow]] = []

    ncbi_first = first_n(ncbi_rows, ncbi_only, 50)
    if ncbi_first:
        targets.append(("ncbi_only_highish_first50", max(ncbi_first, key=lambda row: row.bitscore_float)))

    losat_first = first_n(losat_rows, losat_only, 50)
    if losat_first:
        targets.append(("losat_only_highish_first50", max(losat_first, key=lambda row: row.bitscore_float)))

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
    parser.add_argument(
        "--losat-bin",
        type=Path,
        default=Path("target/release/LOSAT"),
        help="LOSAT executable to use with --refresh-losat, relative to cwd unless absolute.",
    )
    parser.add_argument("--print-refresh-commands", action="store_true")
    parser.add_argument("--fail-on-diff", action="store_true")
    args = parser.parse_args()

    root = Path.cwd()
    cases: list[tuple[str, Path, Path]] = []
    refresh_rows: list[dict[str, str]] = []
    if args.ncbi or args.losat:
        if not (args.ncbi and args.losat):
            parser.error("--ncbi and --losat must be provided together")
        if args.refresh_losat:
            parser.error("--refresh-losat requires manifest-selected cases, not direct --ncbi/--losat paths")
        cases.append(("direct", args.ncbi, args.losat))
    else:
        wanted = set(args.case_id or [])
        for row in read_manifest(args.manifest):
            case_id = row["case_id"]
            if wanted and case_id not in wanted:
                continue
            refresh_rows.append(row)
            cases.append((case_id, root / row["ncbi_out"], root / row["losat_out"]))

    if not cases:
        parser.error("no cases selected")

    if args.refresh_losat:
        losat_bin = args.losat_bin if args.losat_bin.is_absolute() else root / args.losat_bin
        for row in refresh_rows:
            refresh_losat_output(root, row, losat_bin, args.print_refresh_commands)

    results = []
    failed = False
    for case_id, ncbi_path, losat_path in cases:
        result = compare_case(case_id, ncbi_path, losat_path, args.limit)
        results.append(result)
        failed |= bool(
            result["ncbi_only"]
            or result["losat_only"]
            or result["bitscore_diffs"]
            or result["evalue_diffs"]
            or result["pident_diffs"]
            or result["other_value_diffs"]
        )

    if args.write_targets:
        write_targets(args.write_targets, results)
        print(f"wrote trace targets: {args.write_targets}")

    return 1 if args.fail_on_diff and failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
