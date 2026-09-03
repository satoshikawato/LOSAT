#!/usr/bin/env python3
"""Certify the integrated LOSAT v0.1.0 Linux native/serial-Wasm runtime."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tomllib
from types import ModuleType
from typing import Iterable, Sequence


EXPECTED_NATIVE_COUNTS = {"blastn": 14, "blastp": 9, "tblastx": 20}
EXPECTED_WASM_COUNTS = {"blastn": 14, "blastp": 7, "tblastx": 20}
EXPECTED_BLASTN_CLASSIFICATIONS = {
    "EXACT_TEXT": 13,
    "SOURCE_UNDETERMINED_ACCEPTED": 1,
}
EXPECTED_TBLASTX_CLASSIFICATIONS = {"EXACT_TEXT": 14, "HSP_SET_DIFF": 6}
BLASTP_THREAD_PAIRS = {
    "pairwise_default_thread4": "pairwise_default_serial",
    "all_hsp_forward_thread4": "all_hsp_forward_serial",
}
REPEATABILITY_REPRESENTATIVES = (
    ("blastn", "PesePMNV.MjPMNV.task_blastn", "ordinary_exact"),
    ("blastn", "Sakai.MG1655.megablast", "equal_hsp"),
    ("blastp", "pairwise_default_serial", "representative"),
    ("tblastx", "p11_avclpv_psclpv", "linking_heavy_exact"),
    ("tblastx", "p14_ap027131_ap027133_query4", "query_gencode_exact"),
    ("tblastx", "d06_ap027131_ap027133_db4", "db_gencode_deviation"),
)


@dataclass(frozen=True)
class OracleSpec:
    program: str
    path: Path
    sha256: str
    version_needles: tuple[str, ...]


class CertificationFailure(RuntimeError):
    """The integrated evidence is outside the frozen certification contract."""


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_authority(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise CertificationFailure(f"cannot load certification authority: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def validate_output_dir(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(Path("/tmp")):
        raise CertificationFailure("integrated evidence must remain under /tmp")
    if resolved.exists() and any(resolved.iterdir()):
        raise CertificationFailure(f"evidence directory must be new or empty: {resolved}")
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def run_capture(command: Sequence[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command), cwd=cwd, capture_output=True, text=True, check=False
    )


def run_logged(
    command: Sequence[str],
    cwd: Path,
    evidence_prefix: Path,
    *,
    environment: dict[str, str] | None = None,
) -> None:
    evidence_prefix.parent.mkdir(parents=True, exist_ok=True)
    print(f"[certification] {evidence_prefix}", flush=True)
    Path(f"{evidence_prefix}.command.json").write_text(
        json.dumps({"argv": list(command), "cwd": str(cwd)}, indent=2) + "\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    Path(f"{evidence_prefix}.stdout").write_text(
        completed.stdout, encoding="utf-8"
    )
    Path(f"{evidence_prefix}.stderr").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode != 0:
        raise CertificationFailure(
            f"command failed ({completed.returncode}): {shlex.join(command)}"
        )


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    materialized = list(rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(materialized)


def replace_output(command: Sequence[str], option: str, output: Path) -> list[str]:
    updated = list(command)
    try:
        option_index = updated.index(option)
    except ValueError as error:
        raise CertificationFailure(f"output option absent from command: {option}") from error
    updated[option_index + 1] = str(output)
    return updated


def wasm_command(node: str, runner: Path, wasm_bin: Path, losat_command: Sequence[str]) -> list[str]:
    if not losat_command:
        raise CertificationFailure("cannot convert an empty LOSAT command to Wasm")
    return [node, str(runner), str(wasm_bin), *losat_command[1:]]


def assert_count(actual: int, expected: int, label: str) -> None:
    if actual != expected:
        raise CertificationFailure(f"{label}: expected {expected}, observed {actual}")


def assert_counter(actual: Counter[str], expected: dict[str, int], label: str) -> None:
    if dict(actual) != expected:
        raise CertificationFailure(
            f"{label}: expected {expected}, observed {dict(actual)}"
        )


def validate_git_identity(repo_root: Path, expected_sha: str) -> str:
    head = run_capture(["git", "rev-parse", "HEAD"], repo_root)
    if head.returncode != 0:
        raise CertificationFailure(head.stderr.strip() or "cannot resolve git HEAD")
    observed = head.stdout.strip()
    if observed != expected_sha:
        raise CertificationFailure(
            f"certification SHA mismatch: expected {expected_sha}, observed {observed}"
        )
    output_paths = [
        "LOSAT/src",
        "LOSAT/Cargo.toml",
        "LOSAT/Cargo.lock",
        "LOSAT/.cargo/config.toml",
    ]
    production_diff = subprocess.run(
        [
            "git",
            "diff",
            "--ignore-cr-at-eol",
            "--quiet",
            expected_sha,
            "--",
            *output_paths,
        ],
        cwd=repo_root,
        check=False,
    )
    if production_diff.returncode != 0:
        raise CertificationFailure(
            "output-affecting production/build files differ from the certification SHA"
        )
    return observed


def record_toolchain(repo_root: Path, output_dir: Path) -> dict[str, object]:
    rustc = run_capture(["rustc", "-vV"], repo_root)
    cargo = run_capture(["cargo", "-V"], repo_root)
    node = run_capture(["node", "--version"], repo_root)
    for name, completed in (("rustc", rustc), ("cargo", cargo), ("node", node)):
        if completed.returncode != 0:
            raise CertificationFailure(f"cannot record {name} identity")
    host = next(
        (line.removeprefix("host: ") for line in rustc.stdout.splitlines() if line.startswith("host: ")),
        "",
    )
    if not host:
        raise CertificationFailure("rustc -vV did not report a host triple")
    cargo_config_path = repo_root / "LOSAT" / ".cargo" / "config.toml"
    cargo_config = tomllib.loads(cargo_config_path.read_text(encoding="utf-8"))
    wasm_rustflags = (
        cargo_config.get("target", {})
        .get("wasm32-wasip1", {})
        .get("rustflags", [])
    )
    relevant_env = {
        key: os.environ.get(key, "")
        for key in (
            "RUSTFLAGS",
            "CARGO_ENCODED_RUSTFLAGS",
            "CARGO_BUILD_RUSTFLAGS",
            "RUSTC_LINKER",
            "CC",
            "CFLAGS",
            "CPPFLAGS",
            "LDFLAGS",
        )
    }
    record: dict[str, object] = {
        "rustc_vV": rustc.stdout.splitlines(),
        "cargo_V": cargo.stdout.strip(),
        "node_version": node.stdout.strip(),
        "host_triple": host,
        "native": {
            "target_triple": host,
            "cargo_features": "default",
            "profile": "release",
            "locked": True,
            "build_command": "cargo build --release --locked",
        },
        "serial_wasm": {
            "target_triple": "wasm32-wasip1",
            "cargo_features": "no-default-features",
            "profile": "release",
            "locked": True,
            "build_command": (
                "cargo build --release --target wasm32-wasip1 "
                "--no-default-features --locked"
            ),
            "target_rustflags": wasm_rustflags,
        },
        "environment_flags": relevant_env,
        "cargo_config": str(cargo_config_path.relative_to(repo_root)),
        "cargo_config_sha256": sha256_path(cargo_config_path),
        "output_affecting_linker_build_flags": [],
    }
    (output_dir / "cert_toolchain.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    return record


def validate_oracles(
    repo_root: Path,
    output_dir: Path,
    specs: Sequence[OracleSpec],
    source_archive: Path,
    source_archive_sha256: str,
) -> dict[str, object]:
    identities: dict[str, object] = {}
    for oracle in specs:
        if not oracle.path.is_file() or not oracle.path.stat().st_mode & 0o111:
            raise CertificationFailure(f"official oracle is missing: {oracle.path}")
        observed_hash = sha256_path(oracle.path)
        if observed_hash != oracle.sha256:
            raise CertificationFailure(
                f"{oracle.program} oracle SHA-256 mismatch: {observed_hash}"
            )
        version = run_capture([str(oracle.path), "-version"], repo_root)
        version_text = (version.stdout + version.stderr).strip()
        if version.returncode != 0 or any(
            needle not in version_text for needle in oracle.version_needles
        ):
            raise CertificationFailure(
                f"{oracle.program} oracle version mismatch: {version_text}"
            )
        identities[oracle.program] = {
            "path": str(oracle.path),
            "sha256": observed_hash,
            "version": version_text.splitlines(),
        }
    if not source_archive.is_file():
        raise CertificationFailure(f"official source archive is missing: {source_archive}")
    observed_source_hash = sha256_path(source_archive)
    if observed_source_hash != source_archive_sha256:
        raise CertificationFailure(
            f"official source archive SHA-256 mismatch: {observed_source_hash}"
        )
    identities["official_source_archive"] = {
        "path": str(source_archive),
        "sha256": observed_source_hash,
    }
    (output_dir / "oracle_identities.json").write_text(
        json.dumps(identities, indent=2) + "\n", encoding="utf-8"
    )
    return identities


def run_boundary(repo_root: Path, output_dir: Path) -> None:
    checker = repo_root / "LOSAT" / "tests" / "check_pure_rust_runtime_boundary.py"
    completed = run_capture([sys.executable, str(checker)], repo_root)
    (output_dir / "zero_delegation.stdout").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output_dir / "zero_delegation.stderr").write_text(
        completed.stderr, encoding="utf-8"
    )
    required_lines = {
        "temporary_production_findings\t0",
        "cargo_links_findings\t2",
        "result\tPASS",
    }
    observed_lines = set(completed.stdout.splitlines())
    summary = next(
        (line for line in completed.stdout.splitlines() if line.startswith("summary\t")),
        "",
    )
    required_summary_parts = {
        "temporary_allowed=0",
        "temporary_unexpected=0",
        "temporary_stale=0",
        "temporary_invalid=0",
        "dependencies_reviewed=2",
        "dependencies_unexpected=0",
        "dependencies_stale=0",
        "dependencies_invalid=0",
    }
    if (
        completed.returncode != 0
        or not required_lines.issubset(observed_lines)
        or not required_summary_parts.issubset(set(summary.split("\t")))
    ):
        raise CertificationFailure("zero-delegation boundary evidence is not exact")


def certify(
    args: argparse.Namespace,
    output_dir: Path,
    blastn_compare: ModuleType,
    blastn_cert: ModuleType,
    blastp_audit: ModuleType,
    tblastx_audit: ModuleType,
) -> dict[str, object]:
    repo_root = args.repo_root.resolve()
    tests_dir = repo_root / "LOSAT" / "tests"
    native_bin = args.native_bin.resolve()
    wasm_bin = args.wasm_bin.resolve()
    runner = args.wasm_runner.resolve()
    for artifact in (native_bin, wasm_bin, runner):
        if not artifact.is_file():
            raise CertificationFailure(f"required certification artifact is missing: {artifact}")

    node_environment = os.environ.copy()
    node_environment["NODE_NO_WARNINGS"] = "1"
    native_single_environment = os.environ.copy()
    native_single_environment["RAYON_NUM_THREADS"] = "1"

    native_contracts: list[dict[str, object]] = []
    wasm_equalities: list[dict[str, object]] = []
    native_outputs: dict[tuple[str, str], Path] = {}
    wasm_outputs: dict[tuple[str, str], Path] = {}
    native_commands: dict[tuple[str, str], list[str]] = {}
    wasm_commands: dict[tuple[str, str], list[str]] = {}
    native_environments: dict[tuple[str, str], dict[str, str] | None] = {}

    # NCBI reference:
    # ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-832
    # Existing program authorities construct and classify the exact raw output;
    # this orchestration layer never rewrites or sorts it.
    blastn_manifest = tests_dir / "blastn_parity_manifest.tsv"
    blastn_exceptions = tests_dir / "blastn_v010_source_exceptions.tsv"
    blastn_rows = blastn_compare.read_manifest(blastn_manifest)
    assert_count(len(blastn_rows), EXPECTED_NATIVE_COUNTS["blastn"], "BLASTN manifest")
    blastn_dir = output_dir / "native" / "blastn"
    blastn_dir.mkdir(parents=True)
    for row in blastn_rows:
        case_id = row["case_id"]
        command_row = dict(row)
        command_row["query"] = str((repo_root / "LOSAT" / row["query"]).resolve())
        command_row["subject"] = str(
            (repo_root / "LOSAT" / row["subject"]).resolve()
        )
        ncbi_output = blastn_dir / f"{case_id}.ncbi.out"
        losat_output = blastn_dir / f"{case_id}.losat.out"
        ncbi_command = blastn_compare.build_ncbi_command(
            command_row, str(args.blastn_oracle), ncbi_output
        )
        losat_command = blastn_compare.build_losat_command(
            command_row, native_bin, losat_output
        )
        run_logged(ncbi_command, repo_root, blastn_dir / f"{case_id}.ncbi")
        run_logged(losat_command, repo_root, blastn_dir / f"{case_id}.losat")
        native_outputs[("blastn", case_id)] = losat_output
        native_commands[("blastn", case_id)] = losat_command
        native_environments[("blastn", case_id)] = None
    blastn_results = blastn_cert.certify_suite(
        blastn_manifest, blastn_dir, blastn_exceptions
    )
    assert_counter(
        Counter(result.classification for result in blastn_results),
        EXPECTED_BLASTN_CLASSIFICATIONS,
        "BLASTN native classifications",
    )
    for result in blastn_results:
        output = native_outputs[("blastn", result.case_id)]
        native_contracts.append(
            {
                "program": "blastn",
                "case_id": result.case_id,
                "contract": result.classification,
                "classification": result.classification,
                "losat_sha256": sha256_path(output),
            }
        )

    blastp_manifest = tests_dir / "blastp_v010_parity_manifest.tsv"
    blastp_cases = blastp_audit.load_manifest(blastp_manifest, repo_root)
    assert_count(len(blastp_cases), EXPECTED_NATIVE_COUNTS["blastp"], "BLASTP manifest")
    blastp_dir = output_dir / "native" / "blastp"
    blastp_dir.mkdir(parents=True)
    for case in blastp_cases:
        ncbi_command, losat_command = blastp_audit.build_commands(
            case, args.blastp_oracle, native_bin, blastp_dir
        )
        run_logged(ncbi_command, repo_root, blastp_dir / f"{case.case_id}.ncbi")
        run_logged(losat_command, repo_root, blastp_dir / f"{case.case_id}.losat")
        ncbi_output = blastp_dir / f"{case.case_id}.ncbi.tsv"
        losat_output = blastp_dir / f"{case.case_id}.losat.tsv"
        classification, _ = blastp_audit.classify_outputs(ncbi_output, losat_output)
        if classification not in blastp_audit.EXACT_CLASSIFICATIONS:
            raise CertificationFailure(
                f"BLASTP {case.case_id}: authority rejected {classification}"
            )
        native_outputs[("blastp", case.case_id)] = losat_output
        native_commands[("blastp", case.case_id)] = losat_command
        native_environments[("blastp", case.case_id)] = None
        native_contracts.append(
            {
                "program": "blastp",
                "case_id": case.case_id,
                "contract": "source_defined_exact",
                "classification": classification,
                "losat_sha256": sha256_path(losat_output),
            }
        )

    tblastx_manifest = tests_dir / "tblastx_v010_parity_manifest.tsv"
    tblastx_cases = tblastx_audit.load_manifest(tblastx_manifest, repo_root)
    assert_count(len(tblastx_cases), EXPECTED_NATIVE_COUNTS["tblastx"], "TBLASTX manifest")
    tblastx_dir = output_dir / "native" / "tblastx"
    tblastx_dir.mkdir(parents=True)
    tblastx_paths: dict[str, tuple[Path, Path]] = {}
    tblastx_classifications: Counter[str] = Counter()
    for case in tblastx_cases:
        case_dir = tblastx_dir / case.case_id
        case_dir.mkdir()
        ncbi_output = case_dir / "ncbi.tsv"
        losat_output = case_dir / "losat.tsv"
        ncbi_command, losat_command = tblastx_audit.build_commands(
            case, args.tblastx_oracle, native_bin, ncbi_output, losat_output
        )
        run_logged(ncbi_command, repo_root, case_dir / "ncbi")
        run_logged(
            losat_command,
            repo_root,
            case_dir / "losat",
            environment=native_single_environment,
        )
        classification, _ = tblastx_audit.classify_outputs(
            ncbi_output, losat_output
        )
        if not tblastx_audit.contract_accepts(case, classification):
            raise CertificationFailure(
                f"TBLASTX {case.case_id}: contract rejected {classification}"
            )
        if not tblastx_audit.output_ids_are_valid(case, losat_output):
            raise CertificationFailure(f"TBLASTX {case.case_id}: output IDs are invalid")
        tblastx_classifications[classification] += 1
        tblastx_paths[case.case_id] = (ncbi_output, losat_output)
        native_outputs[("tblastx", case.case_id)] = losat_output
        native_commands[("tblastx", case.case_id)] = losat_command
        native_environments[("tblastx", case.case_id)] = native_single_environment
        native_contracts.append(
            {
                "program": "tblastx",
                "case_id": case.case_id,
                "contract": case.contract,
                "classification": classification,
                "losat_sha256": sha256_path(losat_output),
            }
        )
    assert_counter(
        tblastx_classifications,
        EXPECTED_TBLASTX_CLASSIFICATIONS,
        "TBLASTX native classifications",
    )
    implicit_case = next(
        (case for case in tblastx_cases if case.gencode_args == "implicit"), None
    )
    if implicit_case is None:
        raise CertificationFailure("TBLASTX implicit default-code case is absent")
    implicit_ncbi, implicit_losat = tblastx_paths[implicit_case.case_id]
    default_probe = tblastx_audit.run_default_code_equivalence_probe(
        implicit_case,
        args.tblastx_oracle,
        native_bin,
        tblastx_dir,
        implicit_ncbi,
        implicit_losat,
    )
    if default_probe.get("status") != "PASS":
        raise CertificationFailure("TBLASTX implicit/explicit default-code probe failed")

    assert_count(len(native_contracts), 43, "native contract total")
    assert_counter(
        Counter(str(row["program"]) for row in native_contracts),
        EXPECTED_NATIVE_COUNTS,
        "native program counts",
    )

    serial_blastp_cases = [case for case in blastp_cases if case.num_threads == 1]
    assert_count(
        len(serial_blastp_cases), EXPECTED_WASM_COUNTS["blastp"], "serial BLASTP rows"
    )
    wasm_case_groups = {
        "blastn": blastn_rows,
        "blastp": serial_blastp_cases,
        "tblastx": tblastx_cases,
    }
    for program, cases in wasm_case_groups.items():
        program_dir = output_dir / "serial_wasm" / program
        program_dir.mkdir(parents=True)
        for case in cases:
            case_id = case["case_id"] if program == "blastn" else case.case_id
            native_output = native_outputs[(program, case_id)]
            if program == "blastn":
                command_case = dict(case)
                command_case["query"] = str(
                    (repo_root / "LOSAT" / case["query"]).resolve()
                )
                command_case["subject"] = str(
                    (repo_root / "LOSAT" / case["subject"]).resolve()
                )
                wasm_output = program_dir / f"{case_id}.losat.out"
                losat_template = blastn_compare.build_losat_command(
                    command_case, wasm_bin, wasm_output
                )
            elif program == "blastp":
                _, losat_template = blastp_audit.build_commands(
                    case, args.blastp_oracle, wasm_bin, program_dir
                )
                wasm_output = program_dir / f"{case_id}.losat.tsv"
            else:
                case_dir = program_dir / case_id
                case_dir.mkdir()
                wasm_output = case_dir / "losat.tsv"
                _, losat_template = tblastx_audit.build_commands(
                    case,
                    args.tblastx_oracle,
                    wasm_bin,
                    case_dir / "unused.ncbi.tsv",
                    wasm_output,
                )
            command = wasm_command(args.node, runner, wasm_bin, losat_template)
            log_prefix = (
                program_dir / case_id
                if program != "tblastx"
                else program_dir / case_id / "losat"
            )
            run_logged(
                command,
                repo_root,
                log_prefix,
                environment=node_environment,
            )
            native_bytes = native_output.read_bytes()
            wasm_bytes = wasm_output.read_bytes()
            status = "EXACT_BYTES" if native_bytes == wasm_bytes else "DIFFER"
            wasm_equalities.append(
                {
                    "program": program,
                    "case_id": case_id,
                    "native_sha256": hashlib.sha256(native_bytes).hexdigest(),
                    "wasm_sha256": hashlib.sha256(wasm_bytes).hexdigest(),
                    "status": status,
                }
            )
            if status != "EXACT_BYTES":
                raise CertificationFailure(f"{program} {case_id}: native/Wasm divergence")
            wasm_outputs[(program, case_id)] = wasm_output
            wasm_commands[(program, case_id)] = command

    assert_count(len(wasm_equalities), 41, "native/Wasm equality total")
    assert_counter(
        Counter(str(row["program"]) for row in wasm_equalities),
        EXPECTED_WASM_COUNTS,
        "native/Wasm program counts",
    )

    thread_equalities: list[dict[str, object]] = []
    for threaded_case, serial_case in BLASTP_THREAD_PAIRS.items():
        threaded = native_outputs[("blastp", threaded_case)].read_bytes()
        serial = native_outputs[("blastp", serial_case)].read_bytes()
        status = "EXACT_BYTES" if threaded == serial else "DIFFER"
        thread_equalities.append(
            {
                "threaded_case": threaded_case,
                "serial_case": serial_case,
                "threaded_sha256": hashlib.sha256(threaded).hexdigest(),
                "serial_sha256": hashlib.sha256(serial).hexdigest(),
                "status": status,
            }
        )
        if status != "EXACT_BYTES":
            raise CertificationFailure(
                f"BLASTP thread equality failed: {threaded_case} vs {serial_case}"
            )

    repeatability: list[dict[str, object]] = []
    for program, case_id, semantic_class in REPEATABILITY_REPRESENTATIVES:
        key = (program, case_id)
        if key not in native_outputs or key not in wasm_outputs:
            raise CertificationFailure(f"repeatability representative is absent: {key}")
        for target, first_output, template, option, environment in (
            (
                "native",
                native_outputs[key],
                native_commands[key],
                "-o" if program == "blastn" else "--out",
                native_environments[key],
            ),
            (
                "serial_wasm",
                wasm_outputs[key],
                wasm_commands[key],
                "-o" if program == "blastn" else "--out",
                node_environment,
            ),
        ):
            hashes = [sha256_path(first_output)]
            repeat_dir = output_dir / "repeatability" / program / case_id / target
            repeat_dir.mkdir(parents=True)
            for run_index in (2, 3):
                repeat_output = repeat_dir / f"run{run_index}.out"
                command = replace_output(template, option, repeat_output)
                run_logged(
                    command,
                    repo_root,
                    repeat_dir / f"run{run_index}",
                    environment=environment,
                )
                hashes.append(sha256_path(repeat_output))
            status = "REPEATABLE" if len(set(hashes)) == 1 else "NONDETERMINISTIC"
            repeatability.append(
                {
                    "program": program,
                    "case_id": case_id,
                    "semantic_class": semantic_class,
                    "target": target,
                    "runs": 3,
                    "sha256_run1": hashes[0],
                    "sha256_run2": hashes[1],
                    "sha256_run3": hashes[2],
                    "status": status,
                }
            )
            if status != "REPEATABLE":
                raise CertificationFailure(
                    f"{program} {case_id} is nondeterministic on {target}"
                )
    assert_count(
        len({(row["program"], row["case_id"]) for row in repeatability}),
        6,
        "repeatability semantic classes",
    )
    assert_count(len(repeatability), 12, "repeatability target rows")

    write_tsv(
        output_dir / "native_contracts.tsv",
        native_contracts,
        ["program", "case_id", "contract", "classification", "losat_sha256"],
    )
    write_tsv(
        output_dir / "native_wasm_equalities.tsv",
        wasm_equalities,
        ["program", "case_id", "native_sha256", "wasm_sha256", "status"],
    )
    write_tsv(
        output_dir / "blastp_native_thread_equalities.tsv",
        thread_equalities,
        [
            "threaded_case",
            "serial_case",
            "threaded_sha256",
            "serial_sha256",
            "status",
        ],
    )
    write_tsv(
        output_dir / "repeatability.tsv",
        repeatability,
        [
            "program",
            "case_id",
            "semantic_class",
            "target",
            "runs",
            "sha256_run1",
            "sha256_run2",
            "sha256_run3",
            "status",
        ],
    )
    return {
        "native_contracts": {
            "total": len(native_contracts),
            "by_program": dict(Counter(row["program"] for row in native_contracts)),
            "blastn_classifications": dict(
                Counter(
                    row["classification"]
                    for row in native_contracts
                    if row["program"] == "blastn"
                )
            ),
            "blastp_classifications": dict(
                Counter(
                    row["classification"]
                    for row in native_contracts
                    if row["program"] == "blastp"
                )
            ),
            "tblastx_classifications": dict(
                Counter(
                    row["classification"]
                    for row in native_contracts
                    if row["program"] == "tblastx"
                )
            ),
        },
        "native_wasm_equalities": {
            "total": len(wasm_equalities),
            "by_program": dict(Counter(row["program"] for row in wasm_equalities)),
            "all_exact": True,
        },
        "blastp_native_thread_equalities": {
            "total": len(thread_equalities),
            "all_exact": True,
        },
        "repeatability": {
            "semantic_classes": 6,
            "targets": 2,
            "runs_per_target": 3,
            "all_repeatable": True,
        },
        "tblastx_default_code_equivalence": default_probe,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=repo_root)
    parser.add_argument("--expected-sha", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--native-bin",
        type=Path,
        default=repo_root / "LOSAT" / "target" / "release" / "LOSAT",
    )
    parser.add_argument(
        "--wasm-bin",
        type=Path,
        default=repo_root
        / "LOSAT"
        / "target"
        / "wasm32-wasip1"
        / "release"
        / "LOSAT.wasm",
    )
    parser.add_argument(
        "--wasm-runner",
        type=Path,
        default=script_path.with_name("run_losat_wasi.js"),
    )
    parser.add_argument("--node", default="node")
    parser.add_argument(
        "--blastn-oracle",
        type=Path,
        default=Path(
            "/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastn"
        ),
    )
    parser.add_argument(
        "--blastp-oracle",
        type=Path,
        default=Path(
            "/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastp"
        ),
    )
    parser.add_argument(
        "--tblastx-oracle",
        type=Path,
        default=Path(
            "/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/tblastx"
        ),
    )
    parser.add_argument(
        "--source-archive",
        type=Path,
        default=Path(
            "/home/kawato/tools/ncbi-blast-source/ncbi-blast-2.17.0+-src.tar.gz"
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        output_dir = validate_output_dir(args.output_dir)
        repo_root = args.repo_root.resolve()
        certified_sha = validate_git_identity(repo_root, args.expected_sha)
        toolchain = record_toolchain(repo_root, output_dir)
        oracle_specs = (
            OracleSpec(
                "blastn",
                args.blastn_oracle.resolve(),
                "33b64bc67d3149cee2459b2f7766b363323df632cf12c099546de00aea9698b5",
                ("blastn: 2.17.0+", "Package: blast 2.17.0"),
            ),
            OracleSpec(
                "blastp",
                args.blastp_oracle.resolve(),
                "5ce267c04e4988c265357bfbedc64e809545b6fcfae7ff6775266fabbee8ba0e",
                ("blastp: 2.17.0+", "Package: blast 2.17.0"),
            ),
            OracleSpec(
                "tblastx",
                args.tblastx_oracle.resolve(),
                "583e5d60bbd444ac455d20e0956c5aa0aeef675da8daee8204d8f9376ddb8804",
                ("tblastx: 2.17.0+", "Package: blast 2.17.0"),
            ),
        )
        oracles = validate_oracles(
            repo_root,
            output_dir,
            oracle_specs,
            args.source_archive.resolve(),
            "502057a88e9990e34e62758be21ea474cc0ad68d6a63a2e37b2372af1e5ea147",
        )
        run_boundary(repo_root, output_dir)
        tests_dir = repo_root / "LOSAT" / "tests"
        blastn_compare = load_authority(
            "integrated_blastn_compare", tests_dir / "compare_blastn_parity.py"
        )
        blastn_cert = load_authority(
            "integrated_blastn_cert", tests_dir / "certify_blastn_v010.py"
        )
        blastp_audit = load_authority(
            "integrated_blastp_audit", tests_dir / "audit_blastp_v010.py"
        )
        tblastx_audit = load_authority(
            "integrated_tblastx_audit", tests_dir / "audit_tblastx_v010.py"
        )
        result = certify(
            args,
            output_dir,
            blastn_compare,
            blastn_cert,
            blastp_audit,
            tblastx_audit,
        )
        summary = {
            "result": "INTEGRATED_RUNTIME_CERTIFIED",
            "certified_sha": certified_sha,
            "cert_toolchain": toolchain,
            "oracle_identities": oracles,
            "zero_delegation": {
                "temporary_project_owned_production_findings": 0,
                "reviewed_dependency_metadata_signals": 2,
                "unexpected_stale_invalid_findings": 0,
            },
            **result,
        }
        (output_dir / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
    except (CertificationFailure, OSError, ValueError) as error:
        try:
            if "output_dir" in locals():
                (output_dir / "failure.json").write_text(
                    json.dumps({"result": "FAIL", "error": str(error)}, indent=2)
                    + "\n",
                    encoding="utf-8",
                )
        except OSError:
            pass
        print(f"integrated runtime certification: FAIL: {error}", file=sys.stderr)
        return 1
    print("INTEGRATED_RUNTIME_CERTIFIED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
