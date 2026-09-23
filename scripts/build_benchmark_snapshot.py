#!/usr/bin/env python3
"""Build a deterministic benchmark snapshot from completed comparison runs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import platform
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
TESTS = REPO / "LOSAT" / "tests"
ALIGNMENT_HEADER = (
    "program",
    "case_id",
    "implementation",
    "classification",
    "primary_for_distribution",
    "source_row",
    "qseqid",
    "sseqid",
    "pident",
    "length",
)
TIMING_HEADER = (
    "provenance_id",
    "program",
    "mode",
    "case_id",
    "implementation",
    "thread_count",
    "seconds",
    "value_kind",
    "producing_sha",
    "ncbi_version",
    "recorded_date",
    "source_path",
    "source_sha256",
    "raw_elapsed",
    "warmup_count",
    "sample_index",
    "wall_seconds",
    "output_sha256",
    "benchmark_losat_sha",
    "environment_id",
    "tool",
    "contract",
    "benchmark_timestamp",
    "material_environment_id",
    "collection_segment",
    "boot_id",
    "effective_thread_label",
    "effective_thread_evidence",
)


@dataclass(frozen=True)
class TimingCase:
    case_id: str
    program: str
    native_stem: str
    ncbi_stem: str
    contract: str = "EXACT_TEXT"


TIMING_CASES = (
    TimingCase(
        "PesePMNV.MjPMNV.task_blastn",
        "blastn",
        "PesePMNV.MjPMNV.losatn.blastn",
        "PesePMNV.MjPMNV.blastn",
    ),
    TimingCase(
        "Sakai.MG1655.megablast",
        "blastn",
        "Sakai.MG1655.losatn.megablast",
        "Sakai.MG1655.blastn.megablast",
    ),
    TimingCase(
        "WSSV.PajaWSV.blastp",
        "blastp",
        "WSSV.PajaWSV.losatp",
        "WSSV.PajaWSV.BLASTP",
        "TARGET_SEMANTICS_DIFFER",
    ),
    TimingCase(
        "p03_mela_pemojnva",
        "tblastx",
        "MelaMJNV.PemoMJNVA.tlosatx",
        "MelaMJNV.PemoMJNVA.tblastx",
    ),
    TimingCase(
        "d04_ap027131_ap027133_code4",
        "tblastx",
        "AP027131.AP027133.tlosatx",
        "AP027131.AP027133.tblastx",
    ),
    TimingCase(
        "p11_avclpv_psclpv",
        "tblastx",
        "AvCLPV.PsCLPV.tlosatx",
        "AvCLPV.PsCLPV.tblastx",
    ),
)

MODES = (
    ("ncbi_n1", "NCBI BLAST+", "1", "blast_out", "ncbi", ".n1"),
    ("ncbi_n8", "NCBI BLAST+", "8", "blast_out", "ncbi", ".n8"),
    ("losat_native_n1", "LOSAT native", "1", "losat_out", "native", ""),
    ("losat_native_n8", "LOSAT native", "8", "losat_out", "native", ".n8"),
    (
        "losat_wasm_serial",
        "LOSAT serial Wasm",
        "1",
        "losat_out",
        "native",
        ".wasm",
    ),
    (
        "losat_wasm_threads_requested_n8",
        "LOSAT threaded Wasm requested n8",
        "8",
        "losat_out",
        "native",
        ".wasm.n8",
    ),
)


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def command_output(argv: list[str], cwd: Path = REPO) -> str:
    return subprocess.check_output(argv, cwd=cwd, text=True).strip()


def load_verified_output(path: Path) -> dict[str, object]:
    record_path = path.with_suffix(".run.json")
    log_path = path.with_suffix(".log")
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if record.get("status") != "PASS" or record.get("exit_status") != 0:
        raise ValueError(f"incomplete output record: {record_path}")
    if record.get("output_sha256") != sha256(path):
        raise ValueError(f"output checksum mismatch: {path}")
    if record.get("log_sha256") != sha256(log_path):
        raise ValueError(f"log checksum mismatch: {log_path}")
    return record


def parse_wall_seconds(path: Path) -> float:
    matches = re.findall(
        r"^real\s+(?:(\d+)m)?(\d+(?:\.\d+)?)s?\s*$",
        path.read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    if not matches:
        raise ValueError(f"missing wall time: {path}")
    minutes, seconds = matches[-1]
    return float(minutes or 0) * 60 + float(seconds)


# NCBI reference: ncbi-blast/c++/src/app/blastdb/makeblastdb.cpp:236-247
# dbtype and parse_seqids are explicit database-build inputs. Preserve every
# warmup/timed build as separate provenance while excluding it from search time.
def load_database_build(
    timing_root: Path,
    directory: Path,
    case: TimingCase,
    phase: str,
    sample_index: int,
) -> dict[str, object]:
    paths = sorted((directory / "blast_out" / "db").glob("*/*.makeblastdb.json"))
    if len(paths) != 1:
        raise ValueError(f"expected one database build record: {directory}")
    path = paths[0]
    log_path = path.with_suffix(".log")
    record = json.loads(path.read_text(encoding="utf-8"))
    manifest = json.loads((directory / "run.json").read_text(encoding="utf-8"))
    input_path = Path(str(record.get("input", "")))
    executable = Path(str(record.get("executable", "")))
    if manifest.get("status") != "COMPLETE":
        raise ValueError(f"incomplete database-build run: {directory}")
    if record.get("run_id") != manifest.get("run_id"):
        raise ValueError(f"database-build run mismatch: {path}")
    if record.get("exit_status") != 0 or not record.get("ordered_argv"):
        raise ValueError(f"failed database build: {path}")
    if not input_path.is_file() or record.get("input_sha256") != sha256(input_path):
        raise ValueError(f"database input checksum mismatch: {path}")
    if not executable.is_file() or record.get("executable_sha256") != sha256(executable):
        raise ValueError(f"database executable checksum mismatch: {path}")
    if not log_path.is_file() or record.get("log_sha256") != sha256(log_path):
        raise ValueError(f"database log checksum mismatch: {path}")
    if not str(record.get("version", "")).startswith("makeblastdb: 2.17.0+"):
        raise ValueError(f"unexpected database builder version: {path}")
    if float(record.get("wall_seconds", 0)) <= 0:
        raise ValueError(f"invalid database build time: {path}")
    return {
        "case_id": case.case_id,
        "phase": phase,
        "sample_index": sample_index,
        "source_path": str(path.relative_to(timing_root)),
        "record_sha256": sha256(path),
        **record,
    }


def read_cases() -> list[dict[str, str]]:
    with (TESTS / "comparison_cases.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-40
# qaccver saccver pident length mismatch gapopen qstart qend sstart send
# evalue bitscore
# Preserve recorded outfmt 6 pident/length values without biological filtering.
def build_alignment_dataset(
    distribution_run: Path, snapshot: Path
) -> tuple[int, list[dict[str, object]]]:
    target = snapshot / "alignment_results.tsv.gz"
    temporary = snapshot / ".alignment_results.tsv.gz.tmp"
    case_metadata: list[dict[str, object]] = []
    total_rows = 0
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(ALIGNMENT_HEADER)
                for case in read_cases():
                    task = case["task"]
                    program = "blastn" if task == "megablast" else task
                    case_id = f"{task}:{case['name']}"
                    native_stem = case["losat_stem"] + (".n1" if task == "tblastx" else "")
                    ncbi_stem = case["ncbi_stem"] + (
                        ".n1" if task == "tblastx" else ".subject.n1"
                    )
                    sources = (
                        ("NCBI BLAST+", distribution_run / "blast_out" / f"{ncbi_stem}.out"),
                        ("LOSAT", distribution_run / "losat_out" / f"{native_stem}.out"),
                    )
                    source_bytes = [path.read_bytes() for _, path in sources]
                    if source_bytes[0] != source_bytes[1]:
                        raise ValueError(f"distribution output mismatch: {case_id}")
                    source_metadata: dict[str, object] = {}
                    for implementation, path in sources:
                        record = load_verified_output(path)
                        row_count = 0
                        with path.open(encoding="utf-8", newline="") as handle:
                            for source_row, row in enumerate(
                                csv.reader(handle, delimiter="\t"), start=1
                            ):
                                if len(row) != 12:
                                    raise ValueError(f"unexpected outfmt 6 row: {path}:{source_row}")
                                writer.writerow(
                                    (
                                        program,
                                        case_id,
                                        implementation,
                                        "EXACT_TEXT",
                                        "true",
                                        source_row,
                                        row[0],
                                        row[1],
                                        row[2],
                                        row[3],
                                    )
                                )
                                row_count += 1
                                total_rows += 1
                        source_metadata[implementation] = {
                            "alignment_rows": row_count,
                            "argv": record["ordered_argv"],
                            "data_sha256": record["output_sha256"],
                            "record_sha256": sha256(path.with_suffix(".run.json")),
                            "target": (
                                "database"
                                if implementation == "NCBI BLAST+" and task == "tblastx"
                                else "local_subject"
                            ),
                        }
                    case_metadata.append(
                        {
                            "case_id": case_id,
                            "program": program,
                            "task": task,
                            "classification": "EXACT_TEXT",
                            "primary_for_distribution": True,
                            "sources": source_metadata,
                        }
                    )
    temporary.replace(target)
    return total_rows, case_metadata


def timing_output_path(directory: Path, case: TimingCase, mode: tuple[str, ...]) -> Path:
    _, _, _, subdirectory, stem_kind, suffix = mode
    stem = case.ncbi_stem if stem_kind == "ncbi" else case.native_stem
    if stem_kind == "native" and mode[0] == "losat_native_n1" and case.program == "tblastx":
        suffix = ".n1"
    return directory / subdirectory / f"{stem}{suffix}.out"


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3225-3236
# NCBI forces local-subject searches to its minimum thread count. Timing inputs
# therefore select only the separately recorded database-search n1/n8 outputs.
def build_timing_dataset(
    timing_root: Path, snapshot: Path, provenance_id: str, head: str
) -> tuple[int, dict[str, object], list[dict[str, object]]]:
    target = snapshot / "execution_times.tsv"
    boot_id_path = Path("/proc/sys/kernel/random/boot_id")
    boot_id = boot_id_path.read_text().strip() if boot_id_path.is_file() else "unavailable"
    environment_id = f"{platform.system().lower()}-{platform.machine()}-{datetime.now().date()}"
    sample_commands: dict[str, dict[str, object]] = {}
    database_builds: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []
    for case in TIMING_CASES:
        database_builds.append(
            load_database_build(
                timing_root,
                timing_root / "warmup-0" / case.case_id,
                case,
                "warmup",
                0,
            )
        )
        hashes_by_mode: dict[str, set[str]] = {mode[0]: set() for mode in MODES}
        for sample_index in range(1, 4):
            directory = timing_root / f"sample-{sample_index}" / case.case_id
            manifest = json.loads((directory / "run.json").read_text(encoding="utf-8"))
            if manifest.get("status") != "COMPLETE":
                raise ValueError(f"incomplete timing run: {directory}")
            database_builds.append(
                load_database_build(
                    timing_root, directory, case, "timed", sample_index
                )
            )
            for mode in MODES:
                mode_name, implementation, threads, *_ = mode
                output = timing_output_path(directory, case, mode)
                record = load_verified_output(output)
                wall_seconds = parse_wall_seconds(output.with_suffix(".log"))
                output_hash = str(record["output_sha256"])
                hashes_by_mode[mode_name].add(output_hash)
                case_commands = sample_commands.setdefault(case.case_id, {})
                mode_commands = case_commands.setdefault(mode_name, [])
                if not isinstance(mode_commands, list):
                    raise TypeError(f"invalid command collection: {case.case_id} {mode_name}")
                mode_commands.append(
                    {
                        "sample_index": sample_index,
                        "argv": record["ordered_argv"],
                        "output_sha256": output_hash,
                        "record_sha256": sha256(output.with_suffix(".run.json")),
                    }
                )
                effective = (
                    "serial command-Wasm"
                    if mode_name == "losat_wasm_serial"
                    else f"requested/configured {threads}"
                )
                evidence = (
                    "NCBI prebuilt -db command record"
                    if mode_name.startswith("ncbi_")
                    else "recorded command and artifact identity"
                )
                rows.append(
                    dict(
                        zip(
                            TIMING_HEADER,
                            (
                                provenance_id,
                                case.program,
                                mode_name,
                                case.case_id,
                                implementation,
                                threads,
                                wall_seconds,
                                "wall_clock_sample",
                                head,
                                "2.17.0+",
                                str(datetime.now().date()),
                                str(output.relative_to(timing_root)),
                                sha256(output.with_suffix(".run.json")),
                                wall_seconds,
                                1,
                                sample_index,
                                wall_seconds,
                                output_hash,
                                head,
                                environment_id,
                                implementation,
                                case.contract,
                                record.get("ended_utc", ""),
                                environment_id,
                                1,
                                boot_id,
                                effective,
                                evidence,
                            ),
                        )
                    )
                )
        for mode_name, hashes in hashes_by_mode.items():
            if len(hashes) != 1:
                raise ValueError(f"unstable timing output: {case.case_id} {mode_name}")
        if case.contract == "EXACT_TEXT":
            all_hashes = set().union(*hashes_by_mode.values())
            if len(all_hashes) != 1:
                raise ValueError(f"cross-mode timing mismatch: {case.case_id}")

    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=TIMING_HEADER, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows), sample_commands, database_builds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distribution-run", type=Path, required=True)
    parser.add_argument("--timing-root", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    snapshot = args.snapshot.resolve()
    snapshot.mkdir(parents=True, exist_ok=True)
    head = command_output(["git", "rev-parse", "HEAD"])
    short_head = head[:8]
    provenance_id = f"current_{short_head}_{datetime.now().date()}"
    alignment_rows, alignment_cases = build_alignment_dataset(
        args.distribution_run.resolve(), snapshot
    )
    timing_rows, timing_commands, database_builds = build_timing_dataset(
        args.timing_root.resolve(), snapshot, provenance_id, head
    )
    alignment_path = snapshot / "alignment_results.tsv.gz"
    timing_path = snapshot / "execution_times.tsv"
    metadata = {
        "schema_version": 1,
        "snapshot_id": provenance_id,
        "benchmark_completed_at": datetime.now().astimezone().isoformat(),
        "repository_base": head,
        "datasets": {
            "alignment_results": {
                "file": alignment_path.name,
                "sha256": sha256(alignment_path),
                "status": "AVAILABLE_EXACT",
                "row_count": alignment_rows,
                "case_count": len(alignment_cases),
                "plot_weighting": "alignment length",
                "target_protocol": {
                    "tblastx_ncbi": "prebuilt -db with matching -db_gencode",
                    "blastn_blastp_ncbi": "local -subject",
                    "losat": "local -subject",
                },
                "cases": alignment_cases,
            },
            "execution_times": {
                "file": timing_path.name,
                "sha256": sha256(timing_path),
                "status": "AVAILABLE_CURRENT",
                "row_count": timing_rows,
                "current_plot_provenance_id": provenance_id,
                "current_plot_statistic": "median",
                "provenance_groups": {
                    provenance_id: {
                        "used_for_current_plot": True,
                        "status": "AVAILABLE_CURRENT",
                        "benchmark_losat_sha": head,
                        "ncbi_version": "2.17.0+",
                        "protocol": {
                            "warmup_count_per_tool_case": 1,
                            "timed_repetitions_per_tool_case": 3,
                            "output_sink": "regular file",
                            "wall_clock": "Bash time around one process invocation",
                            "database_preparation_included": False,
                            "ncbi_timing_target": "prebuilt -db for every program",
                        },
                        "environment": {
                            "platform": platform.platform(),
                            "machine": platform.machine(),
                            "logical_cpu_count": os.cpu_count(),
                            "node": command_output(["node", "--version"]),
                            "node_versions": json.loads(
                                command_output(["node", "-p", "JSON.stringify(process.versions)"])
                            ),
                            "rustc": command_output(["rustc", "--version"]),
                            "cargo": command_output(["cargo", "--version"]),
                        },
                        "case_commands": timing_commands,
                        "database_builds": database_builds,
                        "target_semantics_note": (
                            "BLASTP database timing intentionally follows the requested NCBI -db target; "
                            "its multi-record result set differs from LOSAT local-subject output and is "
                            "not used for the hit-distribution parity plot."
                        ),
                    }
                },
            },
        },
    }
    (snapshot / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"alignment_rows": alignment_rows, "timing_rows": timing_rows, "snapshot_id": provenance_id}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
