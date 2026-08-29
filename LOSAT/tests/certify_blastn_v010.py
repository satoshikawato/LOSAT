#!/usr/bin/env python3
"""Certify fresh BLASTN v0.1.0 paired outputs against the release contract."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


# NCBI reference: ncbi-blast-2.17.0+-src/c++/src/objtools/align_format/format_flags.cpp:38-40
# ```c
# const char* kDfltArgTabularOutputFmt =
#     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
#     "evalue bitscore";
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
FIELD_INDEX = {field: index for index, field in enumerate(DEFAULT_FIELDS)}
COORD_FIELDS = ("query", "subject", "qstart", "qend", "sstart", "send")
PERMITTED_SOURCE_UNDETERMINED_FIELDS = frozenset(
    {"pident", "length", "mismatch", "gapopen"}
)
EXPECTED_MANIFEST_CASE_COUNT = 14
EXPECTED_SOURCE_EXCEPTION_CASES = frozenset({"Sakai.MG1655.megablast"})
EXCEPTION_COLUMNS = (
    "case_id",
    "allowed_fields",
    "expected_differing_coord_keys",
    "expected_pident_diffs",
    "expected_length_diffs",
    "expected_mismatch_diffs",
    "expected_gapopen_diffs",
)


class CertificationError(RuntimeError):
    """The paired output suite is outside the frozen v0.1.0 contract."""


@dataclass(frozen=True)
class SourceException:
    case_id: str
    allowed_fields: frozenset[str]
    expected_differing_coord_keys: int
    expected_field_diffs: dict[str, int]


@dataclass(frozen=True)
class CaseResult:
    case_id: str
    classification: str
    differing_coord_keys: int = 0
    field_diffs: dict[str, int] | None = None


def _data_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    if not lines:
        raise CertificationError(f"no tabular header found: {path}")
    return list(csv.DictReader(lines, delimiter="\t"))


def read_manifest_case_ids(path: Path) -> list[str]:
    rows = _data_rows(path)
    case_ids: list[str] = []
    for row in rows:
        case_id = (row.get("case_id") or "").strip()
        if not case_id:
            raise CertificationError(f"manifest row has no case_id: {path}")
        case_ids.append(case_id)
    duplicates = [case_id for case_id, count in Counter(case_ids).items() if count > 1]
    if duplicates:
        raise CertificationError(f"duplicate manifest case_id values: {sorted(duplicates)}")
    return case_ids


def _nonnegative_int(row: dict[str, str], field: str, path: Path) -> int:
    try:
        value = int(row[field])
    except (KeyError, ValueError) as error:
        raise CertificationError(f"invalid {field} in {path}: {row.get(field)!r}") from error
    if value < 0:
        raise CertificationError(f"negative {field} in {path}: {value}")
    return value


def read_source_exceptions(path: Path) -> dict[str, SourceException]:
    rows = _data_rows(path)
    exceptions: dict[str, SourceException] = {}
    for row in rows:
        missing = [field for field in EXCEPTION_COLUMNS if field not in row]
        if missing:
            raise CertificationError(f"missing exception columns in {path}: {missing}")
        case_id = row["case_id"].strip()
        if not case_id or case_id in exceptions:
            raise CertificationError(f"invalid or duplicate exception case_id in {path}: {case_id!r}")
        allowed_fields = frozenset(
            field.strip() for field in row["allowed_fields"].split(",") if field.strip()
        )
        if allowed_fields != PERMITTED_SOURCE_UNDETERMINED_FIELDS:
            raise CertificationError(
                f"{case_id}: allowed_fields must be exactly "
                f"{sorted(PERMITTED_SOURCE_UNDETERMINED_FIELDS)}, got {sorted(allowed_fields)}"
            )
        exceptions[case_id] = SourceException(
            case_id=case_id,
            allowed_fields=allowed_fields,
            expected_differing_coord_keys=_nonnegative_int(
                row, "expected_differing_coord_keys", path
            ),
            expected_field_diffs={
                field: _nonnegative_int(row, f"expected_{field}_diffs", path)
                for field in sorted(PERMITTED_SOURCE_UNDETERMINED_FIELDS)
            },
        )
    return exceptions


def paired_paths(output_dir: Path, case_id: str) -> tuple[Path, Path]:
    return (
        output_dir / f"{case_id}.ncbi.out",
        output_dir / f"{case_id}.losat.out",
    )


def _line_body_and_ending(line: bytes) -> tuple[bytes, bytes]:
    if line.endswith(b"\r\n"):
        return line[:-2], b"\r\n"
    if line.endswith(b"\n"):
        return line[:-1], b"\n"
    if line.endswith(b"\r"):
        return line[:-1], b"\r"
    return line, b""


def _parse_tabular_fields(line: bytes, path: Path, line_number: int) -> tuple[list[bytes], bytes]:
    body, ending = _line_body_and_ending(line)
    fields = body.split(b"\t")
    if len(fields) != len(DEFAULT_FIELDS):
        raise CertificationError(
            f"{path}:{line_number}: expected exactly {len(DEFAULT_FIELDS)} fields, got {len(fields)}"
        )
    return fields, ending


def _coord_key(fields: list[bytes]) -> tuple[bytes, ...]:
    return tuple(fields[FIELD_INDEX[field]] for field in COORD_FIELDS)


# NCBI reference: ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/blast_hits.c:2268-2387
# ```c
# if (h1->subject.end < h2->subject.end)
#    return 1;
# if (h1->subject.end > h2->subject.end)
#    return -1;
# return 0;
# ```
# The comparator returns equality without examining gap_info. Certification
# therefore permits only the four reviewed edit-script-derived output fields;
# coordinates, scores, order, headers, delimiters, and line endings stay exact.
def certify_source_exception(
    exception: SourceException, ncbi_path: Path, losat_path: Path
) -> CaseResult:
    ncbi_lines = ncbi_path.read_bytes().splitlines(keepends=True)
    losat_lines = losat_path.read_bytes().splitlines(keepends=True)
    if len(ncbi_lines) != len(losat_lines):
        raise CertificationError(
            f"{exception.case_id}: physical line count differs "
            f"(NCBI={len(ncbi_lines)}, LOSAT={len(losat_lines)})"
        )

    ncbi_coord_keys: Counter[tuple[bytes, ...]] = Counter()
    losat_coord_keys: Counter[tuple[bytes, ...]] = Counter()
    differing_coord_keys: set[tuple[bytes, ...]] = set()
    field_diffs: Counter[str] = Counter()
    ncbi_row_count = 0
    losat_row_count = 0

    for line_number, (ncbi_line, losat_line) in enumerate(
        zip(ncbi_lines, losat_lines), start=1
    ):
        ncbi_body, _ = _line_body_and_ending(ncbi_line)
        losat_body, _ = _line_body_and_ending(losat_line)
        ncbi_is_data = bool(ncbi_body) and not ncbi_body.startswith(b"#")
        losat_is_data = bool(losat_body) and not losat_body.startswith(b"#")
        if ncbi_is_data != losat_is_data:
            raise CertificationError(
                f"{exception.case_id}: output structure differs at line {line_number}"
            )
        if not ncbi_is_data:
            if ncbi_line != losat_line:
                raise CertificationError(
                    f"{exception.case_id}: header/formatter difference at line {line_number}"
                )
            continue

        ncbi_fields, ncbi_ending = _parse_tabular_fields(
            ncbi_line, ncbi_path, line_number
        )
        losat_fields, losat_ending = _parse_tabular_fields(
            losat_line, losat_path, line_number
        )
        ncbi_row_count += 1
        losat_row_count += 1
        if ncbi_ending != losat_ending:
            raise CertificationError(
                f"{exception.case_id}: newline difference at line {line_number}"
            )

        ncbi_key = _coord_key(ncbi_fields)
        losat_key = _coord_key(losat_fields)
        ncbi_coord_keys[ncbi_key] += 1
        losat_coord_keys[losat_key] += 1
        line_differs = False
        for field, index in FIELD_INDEX.items():
            if ncbi_fields[index] == losat_fields[index]:
                continue
            if field not in exception.allowed_fields:
                raise CertificationError(
                    f"{exception.case_id}: unexpected {field} difference at line {line_number}"
                )
            field_diffs[field] += 1
            line_differs = True
        if line_differs:
            if ncbi_key != losat_key:
                raise CertificationError(
                    f"{exception.case_id}: differing row changed its coordinate/HSP key"
                )
            differing_coord_keys.add(ncbi_key)

    if ncbi_row_count != losat_row_count:
        raise CertificationError(
            f"{exception.case_id}: HSP row count differs "
            f"(NCBI={ncbi_row_count}, LOSAT={losat_row_count})"
        )
    if ncbi_coord_keys != losat_coord_keys:
        ncbi_only = ncbi_coord_keys - losat_coord_keys
        losat_only = losat_coord_keys - ncbi_coord_keys
        raise CertificationError(
            f"{exception.case_id}: coordinate/HSP key multiset differs "
            f"(NCBI-only={sum(ncbi_only.values())}, LOSAT-only={sum(losat_only.values())})"
        )

    observed_field_diffs = {
        field: field_diffs[field]
        for field in sorted(PERMITTED_SOURCE_UNDETERMINED_FIELDS)
    }
    if len(differing_coord_keys) != exception.expected_differing_coord_keys:
        raise CertificationError(
            f"{exception.case_id}: differing coordinate footprint changed "
            f"(expected={exception.expected_differing_coord_keys}, "
            f"observed={len(differing_coord_keys)})"
        )
    if observed_field_diffs != exception.expected_field_diffs:
        raise CertificationError(
            f"{exception.case_id}: field-difference footprint changed "
            f"(expected={exception.expected_field_diffs}, observed={observed_field_diffs})"
        )
    return CaseResult(
        case_id=exception.case_id,
        classification="SOURCE_UNDETERMINED_ACCEPTED",
        differing_coord_keys=len(differing_coord_keys),
        field_diffs=observed_field_diffs,
    )


def certify_suite(
    manifest_path: Path,
    output_dir: Path,
    exceptions_path: Path,
    *,
    expected_manifest_case_count: int = EXPECTED_MANIFEST_CASE_COUNT,
    expected_exception_cases: frozenset[str] = EXPECTED_SOURCE_EXCEPTION_CASES,
) -> list[CaseResult]:
    case_ids = read_manifest_case_ids(manifest_path)
    if len(case_ids) != expected_manifest_case_count:
        raise CertificationError(
            f"manifest case count changed (expected={expected_manifest_case_count}, "
            f"observed={len(case_ids)})"
        )
    exceptions = read_source_exceptions(exceptions_path)
    if frozenset(exceptions) != expected_exception_cases:
        raise CertificationError(
            f"source-exception cases changed (expected={sorted(expected_exception_cases)}, "
            f"observed={sorted(exceptions)})"
        )
    unknown_exceptions = set(exceptions) - set(case_ids)
    if unknown_exceptions:
        raise CertificationError(
            f"source exceptions are absent from the manifest: {sorted(unknown_exceptions)}"
        )

    results: list[CaseResult] = []
    for case_id in case_ids:
        ncbi_path, losat_path = paired_paths(output_dir, case_id)
        missing = [str(path) for path in (ncbi_path, losat_path) if not path.is_file()]
        if missing:
            raise CertificationError(f"{case_id}: missing paired output: {missing}")
        if case_id in exceptions:
            results.append(
                certify_source_exception(exceptions[case_id], ncbi_path, losat_path)
            )
            continue
        if ncbi_path.read_bytes() != losat_path.read_bytes():
            raise CertificationError(
                f"{case_id}: source-defined output is not byte-exact"
            )
        results.append(CaseResult(case_id=case_id, classification="EXACT_TEXT"))
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("tests/blastn_parity_manifest.tsv"),
    )
    parser.add_argument(
        "--paired-output-dir",
        type=Path,
        required=True,
        help="Directory produced by compare_blastn_parity.py --fresh-paired.",
    )
    parser.add_argument(
        "--exceptions",
        type=Path,
        default=Path("tests/blastn_v010_source_exceptions.tsv"),
    )
    args = parser.parse_args()

    try:
        results = certify_suite(args.manifest, args.paired_output_dir, args.exceptions)
    except (CertificationError, OSError, UnicodeError) as error:
        print(f"BLASTN v0.1.0 certification: FAIL: {error}", file=sys.stderr)
        return 1

    for result in results:
        print(f"{result.case_id}\t{result.classification}")
        if result.field_diffs is not None:
            print(f"  differing_coordinate_keys={result.differing_coord_keys}")
            for field in sorted(result.field_diffs):
                print(f"  {field}_differences={result.field_diffs[field]}")
    exact_count = sum(result.classification == "EXACT_TEXT" for result in results)
    exception_count = len(results) - exact_count
    print(
        "BLASTN v0.1.0 certification: PASS "
        f"({len(results)} cases; {exact_count} exact; {exception_count} source-underdetermined)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
