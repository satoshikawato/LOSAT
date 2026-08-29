#!/usr/bin/env python3
"""Run the narrow gbdraw BLASTP parity matrix against the fixed oracle."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Iterable


DEFAULT_ORACLE = Path(
    "/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastp"
)
EXPECTED_ORACLE_SHA256 = (
    "5ce267c04e4988c265357bfbedc64e809545b6fcfae7ff6775266fabbee8ba0e"
)
STANDARD_FIELDS = (
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
)
HSP_KEY_INDEXES = (0, 1, 6, 7, 8, 9)
EXACT_CLASSIFICATIONS = {"EXACT_TEXT"}


@dataclass(frozen=True)
class Case:
    case_id: str
    profile_ids: str
    use_case: str
    query: Path
    subject: Path
    max_hsps_per_subject: int | None
    max_target_seqs: int | None
    num_threads: int
    outfmt: str
    coverage_note: str


@dataclass(frozen=True)
class CommandResult:
    returncode: int
    stdout: str
    stderr: str


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def optional_positive_int(value: str, *, field: str, case_id: str) -> int | None:
    if not value:
        return None
    parsed = int(value)
    if parsed <= 0:
        raise ValueError(f"{case_id}: {field} must be blank or positive")
    return parsed


def load_manifest(path: Path, repo_root: Path) -> list[Case]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    cases: list[Case] = []
    seen: set[str] = set()
    for row in rows:
        case_id = row["case_id"].strip()
        if not case_id or case_id in seen:
            raise ValueError(f"manifest case_id must be non-empty and unique: {case_id!r}")
        seen.add(case_id)
        outfmt = row["outfmt"].strip()
        if outfmt != "6":
            raise ValueError(f"{case_id}: gbdraw certification only permits outfmt 6")
        num_threads = int(row["num_threads"])
        if num_threads <= 0:
            raise ValueError(f"{case_id}: num_threads must be positive")
        query = repo_root / row["query"]
        subject = repo_root / row["subject"]
        if not query.is_file() or not subject.is_file():
            raise ValueError(f"{case_id}: committed query/subject fixture is missing")
        cases.append(
            Case(
                case_id=case_id,
                profile_ids=row["profile_ids"].strip(),
                use_case=row["use_case"].strip(),
                query=query,
                subject=subject,
                max_hsps_per_subject=optional_positive_int(
                    row["max_hsps_per_subject"].strip(),
                    field="max_hsps_per_subject",
                    case_id=case_id,
                ),
                max_target_seqs=optional_positive_int(
                    row["max_target_seqs"].strip(),
                    field="max_target_seqs",
                    case_id=case_id,
                ),
                num_threads=num_threads,
                outfmt=outfmt,
                coverage_note=row["coverage_note"].strip(),
            )
        )
    return cases


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgQuery("query");
# const string kArgSubject("subject");
# const string kArgNumThreads("num_threads");
# const string kArgMaxTargetSequences("max_target_seqs");
# ```
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:190-210
# ```c
# const string kArgMaxHSPsPerSubject("max_hsps");
# ```
def build_commands(case: Case, oracle: Path, losat_bin: Path, outputs: Path) -> tuple[list[str], list[str]]:
    ncbi_output = outputs / f"{case.case_id}.ncbi.tsv"
    losat_output = outputs / f"{case.case_id}.losat.tsv"
    ncbi = [
        str(oracle),
        "-query",
        str(case.query),
        "-subject",
        str(case.subject),
        "-outfmt",
        case.outfmt,
        "-num_threads",
        str(case.num_threads),
    ]
    losat = [
        str(losat_bin),
        "blastp",
        "--query",
        str(case.query),
        "--subject",
        str(case.subject),
        "--outfmt",
        case.outfmt,
        "--num-threads",
        str(case.num_threads),
    ]
    if case.max_hsps_per_subject is not None:
        ncbi.extend(["-max_hsps", str(case.max_hsps_per_subject)])
        losat.extend(["--max-hsps-per-subject", str(case.max_hsps_per_subject)])
    if case.max_target_seqs is not None:
        ncbi.extend(["-max_target_seqs", str(case.max_target_seqs)])
        losat.extend(["--max-target-seqs", str(case.max_target_seqs)])
    ncbi.extend(["-out", str(ncbi_output)])
    losat.extend(["--out", str(losat_output)])
    return ncbi, losat


def run_command(command: list[str], log_path: Path) -> CommandResult:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    log_path.write_text(completed.stderr, encoding="utf-8")
    return CommandResult(completed.returncode, completed.stdout, completed.stderr)


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:39-45
# ```c++
# const char* kDfltArgTabularOutputFmt =
#     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
#     "evalue bitscore";
# ```
def parse_outfmt6(path: Path) -> list[tuple[str, ...]]:
    rows: list[tuple[str, ...]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        columns = tuple(line.split("\t"))
        if len(columns) != len(STANDARD_FIELDS):
            raise ValueError(
                f"{path}:{line_number}: expected {len(STANDARD_FIELDS)} fields, got {len(columns)}"
            )
        rows.append(columns)
    return rows


def hsp_key(row: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(row[index] for index in HSP_KEY_INDEXES)


def classify_outputs(ncbi_path: Path, losat_path: Path) -> tuple[str, dict[str, object]]:
    if not ncbi_path.is_file() or not losat_path.is_file():
        return "EXECUTION_ERROR", {"detail": "one or both output files are absent"}
    ncbi_bytes = ncbi_path.read_bytes()
    losat_bytes = losat_path.read_bytes()
    ncbi_rows = parse_outfmt6(ncbi_path)
    losat_rows = parse_outfmt6(losat_path)
    metrics: dict[str, object] = {
        "ncbi_rows": len(ncbi_rows),
        "losat_rows": len(losat_rows),
        "ncbi_sha256": hashlib.sha256(ncbi_bytes).hexdigest(),
        "losat_sha256": hashlib.sha256(losat_bytes).hexdigest(),
        "ncbi_hsp_keys": len(Counter(map(hsp_key, ncbi_rows))),
        "losat_hsp_keys": len(Counter(map(hsp_key, losat_rows))),
    }
    if ncbi_bytes == losat_bytes:
        return "EXACT_TEXT", metrics
    if ncbi_rows == losat_rows:
        return "EXACT_DATA", metrics
    if Counter(ncbi_rows) == Counter(losat_rows):
        return "ORDER_ONLY", metrics
    ncbi_keys = Counter(map(hsp_key, ncbi_rows))
    losat_keys = Counter(map(hsp_key, losat_rows))
    if ncbi_keys == losat_keys:
        return "VALUE_DIFF", metrics
    ncbi_key_set = set(ncbi_keys)
    losat_key_set = set(losat_keys)
    metrics["missing_key_count"] = len(ncbi_key_set - losat_key_set)
    metrics["extra_key_count"] = len(losat_key_set - ncbi_key_set)
    if losat_key_set < ncbi_key_set:
        return "MISSING", metrics
    return "HSP_SET_DIFF", metrics


def fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split(maxsplit=1)[0])
    return ids


def write_results(path: Path, results: Iterable[dict[str, object]]) -> None:
    rows = list(results)
    fieldnames = [
        "case_id",
        "profile_ids",
        "classification",
        "ncbi_exit",
        "losat_exit",
        "ncbi_rows",
        "losat_rows",
        "ncbi_hsp_keys",
        "losat_hsp_keys",
        "ncbi_no_hit_queries",
        "losat_no_hit_queries",
        "repeatability_runs",
        "repeatability",
        "repeatability_sha256",
        "missing_key_count",
        "extra_key_count",
        "ncbi_sha256",
        "losat_sha256",
        "detail",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str]) -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=script_path.with_name("blastp_v010_parity_manifest.tsv"),
    )
    parser.add_argument("--repo-root", type=Path, default=repo_root)
    parser.add_argument("--oracle", type=Path, default=DEFAULT_ORACLE)
    parser.add_argument(
        "--losat-bin",
        type=Path,
        default=repo_root / "LOSAT" / "target" / "release" / "LOSAT",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/tmp/losat-blastp-v010-audit/fresh"),
    )
    parser.add_argument(
        "--allow-differences",
        action="store_true",
        help="return success after classification even when a case is non-exact",
    )
    parser.add_argument(
        "--repeatability-runs",
        type=int,
        default=3,
        help="number of native LOSATP runs retained for byte-repeatability evidence",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    repo_root = args.repo_root.resolve()
    oracle = args.oracle.resolve()
    losat_bin = args.losat_bin.resolve()
    if not oracle.is_file() or not oracle.stat().st_mode & 0o111:
        raise SystemExit(f"fixed oracle is missing or not executable: {oracle}")
    oracle_sha256 = sha256_path(oracle)
    if oracle_sha256 != EXPECTED_ORACLE_SHA256:
        raise SystemExit(
            f"fixed oracle SHA-256 mismatch: expected {EXPECTED_ORACLE_SHA256}, got {oracle_sha256}"
        )
    version = subprocess.run(
        [str(oracle), "-version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    if "blastp: 2.17.0+" not in version or "Package: blast 2.17.0" not in version:
        raise SystemExit(f"fixed oracle version mismatch: {version}")
    if not losat_bin.is_file() or not losat_bin.stat().st_mode & 0o111:
        raise SystemExit(f"LOSATP release binary is missing or not executable: {losat_bin}")
    if args.repeatability_runs <= 0:
        raise SystemExit("--repeatability-runs must be positive")

    cases = load_manifest(args.manifest.resolve(), repo_root)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "oracle_identity.json").write_text(
        json.dumps(
            {
                "path": str(oracle),
                "sha256": oracle_sha256,
                "version": version.splitlines(),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    results: list[dict[str, object]] = []
    for case in cases:
        ncbi_command, losat_command = build_commands(case, oracle, losat_bin, output_dir)
        ncbi_path = output_dir / f"{case.case_id}.ncbi.tsv"
        losat_path = output_dir / f"{case.case_id}.losat.tsv"
        ncbi_path.unlink(missing_ok=True)
        losat_path.unlink(missing_ok=True)
        ncbi_result = run_command(ncbi_command, output_dir / f"{case.case_id}.ncbi.stderr")
        losat_result = run_command(losat_command, output_dir / f"{case.case_id}.losat.stderr")
        repeat_results = [losat_result]
        repeat_paths = [losat_path]
        for run_index in range(2, args.repeatability_runs + 1):
            repeat_path = output_dir / f"{case.case_id}.losat.run{run_index}.tsv"
            repeat_path.unlink(missing_ok=True)
            repeat_command = list(losat_command)
            repeat_command[-1] = str(repeat_path)
            repeat_results.append(
                run_command(
                    repeat_command,
                    output_dir / f"{case.case_id}.losat.run{run_index}.stderr",
                )
            )
            repeat_paths.append(repeat_path)
        if ncbi_result.returncode or losat_result.returncode:
            classification = "EXECUTION_ERROR"
            metrics: dict[str, object] = {
                "detail": "oracle or LOSATP returned a non-zero exit status"
            }
        else:
            try:
                classification, metrics = classify_outputs(ncbi_path, losat_path)
            except (OSError, UnicodeError, ValueError) as exc:
                classification = "EXECUTION_ERROR"
                metrics = {"detail": str(exc)}
        query_ids = set(fasta_ids(case.query))
        for owner, path in (("ncbi", ncbi_path), ("losat", losat_path)):
            try:
                hit_query_ids = {row[0] for row in parse_outfmt6(path)}
                metrics[f"{owner}_no_hit_queries"] = len(query_ids - hit_query_ids)
            except (OSError, UnicodeError, ValueError):
                metrics[f"{owner}_no_hit_queries"] = ""
        repeatability_hashes = [
            sha256_path(path) for path in repeat_paths if path.is_file()
        ]
        if any(result.returncode for result in repeat_results) or len(repeatability_hashes) != len(
            repeat_paths
        ):
            repeatability = "EXECUTION_ERROR"
        elif len(repeatability_hashes) < 2:
            repeatability = "NOT_CHECKED"
        elif len(set(repeatability_hashes)) == 1:
            repeatability = "REPEATABLE"
        else:
            repeatability = "NONDETERMINISTIC"
        row: dict[str, object] = {
            "case_id": case.case_id,
            "profile_ids": case.profile_ids,
            "classification": classification,
            "ncbi_exit": ncbi_result.returncode,
            "losat_exit": losat_result.returncode,
            "repeatability_runs": args.repeatability_runs,
            "repeatability": repeatability,
            "repeatability_sha256": ",".join(repeatability_hashes),
            **metrics,
        }
        results.append(row)
        print(
            f"{case.case_id}\t{case.profile_ids}\t{classification}\t"
            f"rows={metrics.get('ncbi_rows', '?')}/{metrics.get('losat_rows', '?')}\t"
            f"repeatability={repeatability}"
        )

    write_results(output_dir / "classifications.tsv", results)
    totals = Counter(str(row["classification"]) for row in results)
    (output_dir / "classification_totals.json").write_text(
        json.dumps(dict(sorted(totals.items())), indent=2) + "\n",
        encoding="utf-8",
    )
    non_exact = [
        row
        for row in results
        if row["classification"] not in EXACT_CLASSIFICATIONS
        or (args.repeatability_runs >= 2 and row["repeatability"] != "REPEATABLE")
    ]
    if non_exact and not args.allow_differences:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
