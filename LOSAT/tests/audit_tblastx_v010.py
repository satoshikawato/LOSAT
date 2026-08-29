#!/usr/bin/env python3
"""Audit the narrow gbdraw TBLASTX v0.1.0 profile against the fixed oracle."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Iterable


DEFAULT_ORACLE = Path(
    "/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/tblastx"
)
EXPECTED_ORACLE_SHA256 = (
    "583e5d60bbd444ac455d20e0956c5aa0aeef675da8daee8204d8f9376ddb8804"
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
PARITY_CONTRACT = "parity"
DEVIATION_CONTRACT = "approved_db_gencode_deviation"
EXACT_CLASSIFICATIONS = {"EXACT_TEXT"}
APPROVED_DEVIATION_CLASSIFICATIONS = {"HSP_SET_DIFF"}


@dataclass(frozen=True)
class Case:
    case_id: str
    profile_ids: str
    contract: str
    query: Path
    subject: Path
    query_gencode: int
    db_gencode: int
    gencode_args: str
    num_threads: int
    outfmt: str
    repeatability_required: bool
    coverage_note: str


@dataclass(frozen=True)
class CommandResult:
    returncode: int
    stderr: str


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_yes_no(value: str, *, field: str, case_id: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "yes":
        return True
    if normalized == "no":
        return False
    raise ValueError(f"{case_id}: {field} must be yes or no")


def load_manifest(path: Path, repo_root: Path) -> list[Case]:
    resolved_root = repo_root.resolve()
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    cases: list[Case] = []
    seen: set[str] = set()
    for row in rows:
        case_id = row["case_id"].strip()
        if not case_id or case_id in seen:
            raise ValueError(f"manifest case_id must be non-empty and unique: {case_id!r}")
        seen.add(case_id)
        contract = row["contract"].strip()
        if contract not in {PARITY_CONTRACT, DEVIATION_CONTRACT}:
            raise ValueError(f"{case_id}: unsupported contract {contract!r}")
        query = (resolved_root / row["query"]).resolve()
        subject = (resolved_root / row["subject"]).resolve()
        if not query.is_relative_to(resolved_root) or not subject.is_relative_to(resolved_root):
            raise ValueError(f"{case_id}: fixtures must remain inside the repository")
        if not query.is_file() or not subject.is_file():
            raise ValueError(f"{case_id}: committed query/subject fixture is missing")
        query_gencode = int(row["query_gencode"])
        db_gencode = int(row["db_gencode"])
        gencode_args = row["gencode_args"].strip()
        if gencode_args not in {"implicit", "explicit"}:
            raise ValueError(f"{case_id}: gencode_args must be implicit or explicit")
        if gencode_args == "implicit" and (query_gencode != 1 or db_gencode != 1):
            raise ValueError(f"{case_id}: implicit genetic-code arguments require code 1/1")
        if contract == PARITY_CONTRACT and db_gencode != 1:
            raise ValueError(f"{case_id}: parity cases must retain the default subject code")
        if contract == DEVIATION_CONTRACT and (gencode_args != "explicit" or db_gencode == 1):
            raise ValueError(
                f"{case_id}: the approved deviation requires explicit non-default db_gencode"
            )
        num_threads = int(row["num_threads"])
        if num_threads != 1:
            raise ValueError(f"{case_id}: gbdraw TBLASTX certification requires one thread per job")
        outfmt = row["outfmt"].strip()
        if outfmt != "6":
            raise ValueError(f"{case_id}: gbdraw TBLASTX certification only permits outfmt 6")
        repeatability_required = parse_yes_no(
            row["repeatability_required"],
            field="repeatability_required",
            case_id=case_id,
        )
        if contract == PARITY_CONTRACT and not repeatability_required:
            raise ValueError(f"{case_id}: every parity case requires repeatability")
        cases.append(
            Case(
                case_id=case_id,
                profile_ids=row["profile_ids"].strip(),
                contract=contract,
                query=query,
                subject=subject,
                query_gencode=query_gencode,
                db_gencode=db_gencode,
                gencode_args=gencode_args,
                num_threads=num_threads,
                outfmt=outfmt,
                repeatability_required=repeatability_required,
                coverage_note=row["coverage_note"].strip(),
            )
        )
    if not cases:
        raise ValueError("manifest must contain at least one case")
    if not any(case.contract == DEVIATION_CONTRACT and case.repeatability_required for case in cases):
        raise ValueError("at least one approved-deviation case must require repeatability")
    return cases


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
# ```c++
# const string kArgQuery("query");
# const string kArgSubject("subject");
# const string kArgQueryGeneticCode("query_gencode");
# const string kArgDbGeneticCode("db_gencode");
# const string kArgNumThreads("num_threads");
# ```
def build_commands(
    case: Case,
    oracle: Path,
    losat_bin: Path,
    ncbi_output: Path,
    losat_output: Path,
) -> tuple[list[str], list[str]]:
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
        "tblastx",
        "--query",
        str(case.query),
        "--subject",
        str(case.subject),
        "--outfmt",
        case.outfmt,
        "--num-threads",
        str(case.num_threads),
    ]
    if case.gencode_args == "explicit":
        ncbi.extend(
            [
                "-query_gencode",
                str(case.query_gencode),
                "-db_gencode",
                str(case.db_gencode),
            ]
        )
        losat.extend(
            [
                "--query-gencode",
                str(case.query_gencode),
                "--db-gencode",
                str(case.db_gencode),
            ]
        )
    ncbi.extend(["-out", str(ncbi_output)])
    losat.extend(["--out", str(losat_output)])
    return ncbi, losat


def run_command(command: list[str], log_path: Path, *, native: bool = False) -> CommandResult:
    environment = os.environ.copy()
    if native:
        environment["RAYON_NUM_THREADS"] = "1"
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    log_path.write_text(completed.stderr, encoding="utf-8")
    return CommandResult(completed.returncode, completed.stderr)


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


def data_key_without_scores(row: tuple[str, ...]) -> tuple[str, ...]:
    return row[:10]


def classify_outputs(ncbi_path: Path, losat_path: Path) -> tuple[str, dict[str, object]]:
    if not ncbi_path.is_file() or not losat_path.is_file():
        return "EXECUTION_ERROR", {"detail": "one or both output files are absent"}
    ncbi_bytes = ncbi_path.read_bytes()
    losat_bytes = losat_path.read_bytes()
    ncbi_rows = parse_outfmt6(ncbi_path)
    losat_rows = parse_outfmt6(losat_path)
    ncbi_hsp_keys = Counter(map(hsp_key, ncbi_rows))
    losat_hsp_keys = Counter(map(hsp_key, losat_rows))
    metrics: dict[str, object] = {
        "ncbi_rows": len(ncbi_rows),
        "losat_rows": len(losat_rows),
        "ncbi_sha256": hashlib.sha256(ncbi_bytes).hexdigest(),
        "losat_sha256": hashlib.sha256(losat_bytes).hexdigest(),
        "missing_key_count": len(set(ncbi_hsp_keys) - set(losat_hsp_keys)),
        "extra_key_count": len(set(losat_hsp_keys) - set(ncbi_hsp_keys)),
    }
    if ncbi_bytes == losat_bytes:
        return "EXACT_TEXT", metrics
    if ncbi_rows == losat_rows:
        return "EXACT_DATA", metrics
    if Counter(ncbi_rows) == Counter(losat_rows):
        return "ORDER_ONLY", metrics
    if Counter(map(data_key_without_scores, ncbi_rows)) == Counter(
        map(data_key_without_scores, losat_rows)
    ) or ncbi_hsp_keys == losat_hsp_keys:
        return "VALUE_DIFF", metrics
    if set(losat_hsp_keys) < set(ncbi_hsp_keys):
        return "MISSING", metrics
    return "HSP_SET_DIFF", metrics


def fasta_ids(path: Path) -> set[str]:
    identifiers: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.add(line[1:].strip().split(maxsplit=1)[0])
    return identifiers


def output_ids_are_valid(case: Case, output: Path) -> bool:
    query_ids = fasta_ids(case.query)
    subject_ids = fasta_ids(case.subject)
    return all(row[0] in query_ids and row[1] in subject_ids for row in parse_outfmt6(output))


def contract_accepts(case: Case, classification: str) -> bool:
    if case.contract == PARITY_CONTRACT:
        return classification in EXACT_CLASSIFICATIONS
    return classification in APPROVED_DEVIATION_CLASSIFICATIONS


def write_results(path: Path, results: Iterable[dict[str, object]]) -> None:
    rows = list(results)
    fieldnames = [
        "case_id",
        "profile_ids",
        "contract",
        "classification",
        "contract_status",
        "ncbi_exit",
        "losat_exit",
        "ncbi_rows",
        "losat_rows",
        "missing_key_count",
        "extra_key_count",
        "output_ids_valid",
        "repeatability_runs",
        "repeatability",
        "repeatability_sha256",
        "ncbi_sha256",
        "losat_sha256",
        "detail",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate_output_dir(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(Path("/tmp")):
        raise ValueError("certification output must remain under /tmp")
    return resolved


def run_default_code_equivalence_probe(
    baseline: Case,
    oracle: Path,
    losat_bin: Path,
    output_dir: Path,
    baseline_ncbi: Path,
    baseline_losat: Path,
) -> dict[str, object]:
    explicit = Case(
        case_id=f"{baseline.case_id}_explicit_code1",
        profile_ids=baseline.profile_ids,
        contract=PARITY_CONTRACT,
        query=baseline.query,
        subject=baseline.subject,
        query_gencode=1,
        db_gencode=1,
        gencode_args="explicit",
        num_threads=1,
        outfmt="6",
        repeatability_required=False,
        coverage_note="explicit default-code companion",
    )
    ncbi_output = output_dir / "default_code_explicit.ncbi.tsv"
    losat_output = output_dir / "default_code_explicit.losat.tsv"
    ncbi_output.unlink(missing_ok=True)
    losat_output.unlink(missing_ok=True)
    ncbi_command, losat_command = build_commands(
        explicit, oracle, losat_bin, ncbi_output, losat_output
    )
    ncbi_result = run_command(ncbi_command, output_dir / "default_code_explicit.ncbi.stderr")
    losat_result = run_command(
        losat_command,
        output_dir / "default_code_explicit.losat.stderr",
        native=True,
    )
    passed = (
        ncbi_result.returncode == 0
        and losat_result.returncode == 0
        and baseline_ncbi.read_bytes() == ncbi_output.read_bytes()
        and baseline_losat.read_bytes() == losat_output.read_bytes()
        and ncbi_output.read_bytes() == losat_output.read_bytes()
    )
    result = {
        "baseline_case_id": baseline.case_id,
        "ncbi_exit": ncbi_result.returncode,
        "losat_exit": losat_result.returncode,
        "implicit_ncbi_sha256": sha256_path(baseline_ncbi),
        "explicit_ncbi_sha256": sha256_path(ncbi_output) if ncbi_output.is_file() else "",
        "implicit_losat_sha256": sha256_path(baseline_losat),
        "explicit_losat_sha256": sha256_path(losat_output) if losat_output.is_file() else "",
        "status": "PASS" if passed else "FAIL",
    }
    (output_dir / "default_code_equivalence.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


def parse_args(argv: list[str]) -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=script_path.with_name("tblastx_v010_parity_manifest.tsv"),
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
        default=Path("/tmp/losat-tblastx-v010-certification/final-native"),
    )
    parser.add_argument(
        "--repeatability-runs",
        type=int,
        default=3,
        help="native runs retained for every parity case and the marked deviation regression",
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
    if "tblastx: 2.17.0+" not in version or "Package: blast 2.17.0" not in version:
        raise SystemExit(f"fixed oracle version mismatch: {version}")
    if not losat_bin.is_file() or not losat_bin.stat().st_mode & 0o111:
        raise SystemExit(f"LOSAT release binary is missing or not executable: {losat_bin}")
    if args.repeatability_runs < 3:
        raise SystemExit("--repeatability-runs must be at least 3")

    cases = load_manifest(args.manifest.resolve(), repo_root)
    try:
        output_dir = validate_output_dir(args.output_dir)
    except ValueError as error:
        raise SystemExit(str(error)) from error
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
    paths_by_case: dict[str, tuple[Path, Path]] = {}
    for index, case in enumerate(cases, start=1):
        case_dir = output_dir / case.case_id
        case_dir.mkdir(exist_ok=True)
        ncbi_path = case_dir / "ncbi.tsv"
        losat_path = case_dir / "losat.tsv"
        ncbi_path.unlink(missing_ok=True)
        losat_path.unlink(missing_ok=True)
        ncbi_command, losat_command = build_commands(
            case, oracle, losat_bin, ncbi_path, losat_path
        )
        (case_dir / "commands.json").write_text(
            json.dumps({"ncbi": ncbi_command, "losat": losat_command}, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"[{index}/{len(cases)}] {case.case_id}: oracle", flush=True)
        ncbi_result = run_command(ncbi_command, case_dir / "ncbi.stderr")
        print(f"[{index}/{len(cases)}] {case.case_id}: native run 1", flush=True)
        losat_result = run_command(losat_command, case_dir / "losat.stderr", native=True)
        if ncbi_result.returncode or losat_result.returncode:
            classification = "EXECUTION_ERROR"
            metrics: dict[str, object] = {"detail": "oracle or LOSAT returned non-zero"}
        else:
            try:
                classification, metrics = classify_outputs(ncbi_path, losat_path)
            except (OSError, UnicodeError, ValueError) as error:
                classification = "EXECUTION_ERROR"
                metrics = {"detail": str(error)}

        output_ids_valid = False
        if losat_path.is_file():
            try:
                output_ids_valid = output_ids_are_valid(case, losat_path)
            except (OSError, UnicodeError, ValueError):
                output_ids_valid = False

        repeat_paths = [losat_path]
        repeat_results = [losat_result]
        if case.repeatability_required:
            for run_index in range(2, args.repeatability_runs + 1):
                repeat_path = case_dir / f"losat.run{run_index}.tsv"
                repeat_path.unlink(missing_ok=True)
                _, repeat_command = build_commands(
                    case,
                    oracle,
                    losat_bin,
                    case_dir / "unused.ncbi.tsv",
                    repeat_path,
                )
                print(
                    f"[{index}/{len(cases)}] {case.case_id}: native run {run_index}",
                    flush=True,
                )
                repeat_results.append(
                    run_command(
                        repeat_command,
                        case_dir / f"losat.run{run_index}.stderr",
                        native=True,
                    )
                )
                repeat_paths.append(repeat_path)
        repeat_hashes = [sha256_path(path) for path in repeat_paths if path.is_file()]
        if not case.repeatability_required:
            repeatability = "NOT_REQUIRED"
        elif any(result.returncode for result in repeat_results) or len(repeat_hashes) != len(
            repeat_paths
        ):
            repeatability = "EXECUTION_ERROR"
        elif len(set(repeat_hashes)) == 1:
            repeatability = "REPEATABLE"
        else:
            repeatability = "NONDETERMINISTIC"

        contract_status = (
            "PASS"
            if contract_accepts(case, classification)
            and output_ids_valid
            and (not case.repeatability_required or repeatability == "REPEATABLE")
            else "FAIL"
        )
        row: dict[str, object] = {
            "case_id": case.case_id,
            "profile_ids": case.profile_ids,
            "contract": case.contract,
            "classification": classification,
            "contract_status": contract_status,
            "ncbi_exit": ncbi_result.returncode,
            "losat_exit": losat_result.returncode,
            "output_ids_valid": output_ids_valid,
            "repeatability_runs": len(repeat_paths),
            "repeatability": repeatability,
            "repeatability_sha256": ",".join(repeat_hashes),
            **metrics,
        }
        results.append(row)
        paths_by_case[case.case_id] = (ncbi_path, losat_path)
        print(
            f"[{index}/{len(cases)}] {case.case_id}: {classification} "
            f"rows={metrics.get('ncbi_rows', '?')}/{metrics.get('losat_rows', '?')} "
            f"contract={contract_status} repeatability={repeatability}",
            flush=True,
        )

    baseline = next((case for case in cases if case.gencode_args == "implicit"), None)
    default_probe: dict[str, object] = {"status": "FAIL", "detail": "implicit baseline absent"}
    if baseline is not None:
        baseline_ncbi, baseline_losat = paths_by_case[baseline.case_id]
        if baseline_ncbi.is_file() and baseline_losat.is_file():
            default_probe = run_default_code_equivalence_probe(
                baseline,
                oracle,
                losat_bin,
                output_dir,
                baseline_ncbi,
                baseline_losat,
            )

    write_results(output_dir / "classifications.tsv", results)
    totals = Counter(str(row["classification"]) for row in results)
    (output_dir / "classification_totals.json").write_text(
        json.dumps(dict(sorted(totals.items())), indent=2) + "\n",
        encoding="utf-8",
    )
    failed = [row for row in results if row["contract_status"] != "PASS"]
    if failed or default_probe.get("status") != "PASS":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
