#!/usr/bin/env python3
"""Certify LOSAT v0.1.0 native output on Windows and macOS."""

from __future__ import annotations

import argparse
from collections import Counter
import contextlib
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
from types import ModuleType
from typing import Iterable, Sequence


EXPECTED_COUNTS = {"blastn": 14, "blastp": 9, "tblastx": 20}
EXPECTED_EXECUTIONS = {"matrix": 43, "oracle": 6, "repeatability": 12}
RUNTIME_CERT_SHA = "5845d22ed9842449628a647f29b8c6762511ca59"
PR5_EVIDENCE_MANIFEST_SHA256 = (
    "b9fc98a376d2849274c86b4e4769d2ee38b76025adbb4d63d3da4a5e3e7cdb5c"
)
CERT_TOOLCHAIN = "1.92.0"
EVIDENCE_SCHEMA = 1
SAFE_CASE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
SHA_PATTERN = re.compile(r"^[0-9a-f]{40}$")


@dataclass(frozen=True)
class PlatformSpec:
    platform_id: str
    os_name: str
    machine: str
    target_triple: str
    binary_arch: str
    runner_label: str
    archive_name: str
    archive_md5: str

    @property
    def archive_url(self) -> str:
        quoted = self.archive_name.replace("+", "%2B")
        return (
            "https://ftp.ncbi.nlm.nih.gov/blast/executables/"
            f"blast%2B/2.17.0/{quoted}"
        )

    @property
    def checksum_url(self) -> str:
        return f"{self.archive_url}.md5"


PLATFORMS = {
    "windows-x64": PlatformSpec(
        platform_id="windows-x64",
        os_name="Windows",
        machine="x86_64",
        target_triple="x86_64-pc-windows-msvc",
        binary_arch="x86_64",
        runner_label="windows-2025",
        archive_name="ncbi-blast-2.17.0+-x64-win64.tar.gz",
        archive_md5="dcd973097407a2910061ff4fb51b09fb",
    ),
    "macos-arm64": PlatformSpec(
        platform_id="macos-arm64",
        os_name="Darwin",
        machine="arm64",
        target_triple="aarch64-apple-darwin",
        binary_arch="arm64",
        runner_label="macos-15",
        archive_name="ncbi-blast-2.17.0+-aarch64-macosx.tar.gz",
        archive_md5="8dc685eb284713db76de41a4dabf96ad",
    ),
    "macos-x64": PlatformSpec(
        platform_id="macos-x64",
        os_name="Darwin",
        machine="x86_64",
        target_triple="x86_64-apple-darwin",
        binary_arch="x86_64",
        runner_label="macos-15-intel",
        archive_name="ncbi-blast-2.17.0+-x64-macosx.tar.gz",
        archive_md5="cde7979c0edca21da0567ab172b68b0e",
    ),
}

REPRESENTATIVES = (
    ("blastn", "PesePMNV.MjPMNV.task_blastn", "ordinary_exact"),
    ("blastn", "Sakai.MG1655.megablast", "source_undetermined"),
    ("blastp", "pairwise_default_serial", "representative"),
    ("tblastx", "p03_mela_pemojnva", "ordinary_exact"),
    ("tblastx", "p14_ap027131_ap027133_query4", "query_gencode_exact"),
    ("tblastx", "d06_ap027131_ap027133_db4", "db_gencode_deviation"),
)


@dataclass(frozen=True)
class SearchStep:
    execution_index: int
    step_id: str
    kind: str
    program: str
    case_id: str
    semantic_class: str
    run_index: int
    command: tuple[str, ...]
    environment: tuple[tuple[str, str], ...]
    output_rel: str
    final_output_rel: str
    expected_losat_sha256: str


@dataclass
class Catalog:
    canonical_rows: list[dict[str, str]]
    canonical: dict[tuple[str, str], dict[str, str]]
    blastn_rows: dict[str, dict[str, str]]
    blastp_cases: dict[str, object]
    tblastx_cases: dict[str, object]
    blastn_compare: ModuleType
    blastn_cert: ModuleType
    blastp_audit: ModuleType
    tblastx_audit: ModuleType


class CertificationFailure(RuntimeError):
    """Platform evidence is outside the frozen PR 6 contract."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: object) -> None:
    atomic_write_bytes(
        path,
        (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_path(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
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


def run_capture(command: Sequence[str], cwd: Path) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        list(command), cwd=cwd, capture_output=True, text=False, check=False
    )


def normalize_machine(value: str) -> str:
    normalized = value.strip().lower()
    return {
        "amd64": "x86_64",
        "x64": "x86_64",
        "aarch64": "arm64",
    }.get(normalized, normalized)


def binary_architecture(path: Path) -> str:
    data = path.read_bytes()[:4096]
    if data.startswith(b"MZ") and len(data) >= 64:
        pe_offset = int.from_bytes(data[60:64], "little")
        with path.open("rb") as handle:
            handle.seek(pe_offset)
            header = handle.read(6)
        if header[:4] != b"PE\0\0":
            raise CertificationFailure(f"invalid PE header: {path}")
        return {0x8664: "x86_64", 0xAA64: "arm64"}.get(
            int.from_bytes(header[4:6], "little"), "unknown"
        )
    if len(data) >= 8 and data[:4] in {b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe"}:
        cpu_type = int.from_bytes(data[4:8], "little")
        return {0x01000007: "x86_64", 0x0100000C: "arm64"}.get(
            cpu_type, "unknown"
        )
    if len(data) >= 8 and data[:4] in {b"\xfe\xed\xfa\xcf", b"\xfe\xed\xfa\xce"}:
        cpu_type = int.from_bytes(data[4:8], "big")
        return {0x01000007: "x86_64", 0x0100000C: "arm64"}.get(
            cpu_type, "unknown"
        )
    raise CertificationFailure(f"unsupported native executable format: {path}")


def verify_archive(
    platform_spec: PlatformSpec, archive: Path, checksum: Path
) -> dict[str, str]:
    if archive.name != platform_spec.archive_name:
        raise CertificationFailure(
            f"wrong NCBI archive for {platform_spec.platform_id}: {archive.name}"
        )
    if not archive.is_file() or not checksum.is_file():
        raise CertificationFailure("NCBI archive and published checksum are required")
    checksum_parts = checksum.read_text(encoding="ascii").strip().split()
    if checksum_parts != [platform_spec.archive_md5, platform_spec.archive_name]:
        raise CertificationFailure(
            f"published checksum record is not the pinned {platform_spec.platform_id} record"
        )
    observed = md5_path(archive)
    if observed != platform_spec.archive_md5:
        raise CertificationFailure(
            f"NCBI archive checksum mismatch: expected {platform_spec.archive_md5}, "
            f"observed {observed}"
        )
    return {
        "archive": str(archive.resolve()),
        "archive_name": platform_spec.archive_name,
        "archive_md5": observed,
        "archive_url": platform_spec.archive_url,
        "checksum": str(checksum.resolve()),
        "checksum_url": platform_spec.checksum_url,
    }


def replace_output(command: Sequence[str], output: Path) -> list[str]:
    updated = list(command)
    for option in ("--out", "-o", "-out"):
        if option in updated:
            option_index = updated.index(option)
            updated[option_index + 1] = str(output)
            return updated
    raise CertificationFailure(f"output option absent from command: {command}")


def _canonical_lines(path: Path) -> tuple[list[str], list[str]]:
    comments: list[str] = []
    data: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith("#"):
            comments.append(line.rstrip("\r\n"))
        elif line.strip():
            data.append(line)
    return comments, data


def load_catalog(repo_root: Path) -> Catalog:
    tests_dir = repo_root / "LOSAT" / "tests"
    blastn_compare = load_authority(
        "platform_blastn_compare", tests_dir / "compare_blastn_parity.py"
    )
    blastn_cert = load_authority(
        "platform_blastn_cert", tests_dir / "certify_blastn_v010.py"
    )
    blastp_audit = load_authority(
        "platform_blastp_audit", tests_dir / "audit_blastp_v010.py"
    )
    tblastx_audit = load_authority(
        "platform_tblastx_audit", tests_dir / "audit_tblastx_v010.py"
    )

    canonical_path = tests_dir / "platform_native_v010_canonical.tsv"
    comments, data_lines = _canonical_lines(canonical_path)
    required_comments = {
        "# Canonical LOSAT native output hashes from PR 5 integrated certification.",
        f"# RUNTIME_CERT_SHA={RUNTIME_CERT_SHA}",
        f"# PR5_EVIDENCE_MANIFEST_SHA256={PR5_EVIDENCE_MANIFEST_SHA256}",
    }
    if set(comments) != required_comments:
        raise CertificationFailure("canonical manifest provenance comments changed")
    canonical_rows = list(csv.DictReader(data_lines, delimiter="\t"))
    expected_fields = (
        "program",
        "case_id",
        "contract",
        "classification",
        "losat_sha256",
    )
    if not canonical_rows or tuple(canonical_rows[0]) != expected_fields:
        raise CertificationFailure("canonical manifest columns changed")
    if len(canonical_rows) != 43:
        raise CertificationFailure(
            f"canonical manifest must contain 43 rows, found {len(canonical_rows)}"
        )
    canonical: dict[tuple[str, str], dict[str, str]] = {}
    for row in canonical_rows:
        key = (row["program"], row["case_id"])
        if key in canonical or not SAFE_CASE_ID.fullmatch(row["case_id"]):
            raise CertificationFailure(f"invalid or duplicate canonical case: {key}")
        if not re.fullmatch(r"[0-9a-f]{64}", row["losat_sha256"]):
            raise CertificationFailure(f"invalid canonical SHA-256: {key}")
        canonical[key] = row
    if dict(Counter(row["program"] for row in canonical_rows)) != EXPECTED_COUNTS:
        raise CertificationFailure("canonical program counts changed")

    blastn_manifest = tests_dir / "blastn_parity_manifest.tsv"
    blastn_rows_list = blastn_compare.read_manifest(blastn_manifest)
    blastn_rows = {row["case_id"]: row for row in blastn_rows_list}
    blastp_cases_list = blastp_audit.load_manifest(
        tests_dir / "blastp_v010_parity_manifest.tsv", repo_root
    )
    blastp_cases = {case.case_id: case for case in blastp_cases_list}
    tblastx_cases_list = tblastx_audit.load_manifest(
        tests_dir / "tblastx_v010_parity_manifest.tsv", repo_root
    )
    tblastx_cases = {case.case_id: case for case in tblastx_cases_list}
    authority_ids = {
        "blastn": set(blastn_rows),
        "blastp": set(blastp_cases),
        "tblastx": set(tblastx_cases),
    }
    for program, expected_count in EXPECTED_COUNTS.items():
        canonical_ids = {
            row["case_id"] for row in canonical_rows if row["program"] == program
        }
        if len(authority_ids[program]) != expected_count or canonical_ids != authority_ids[program]:
            raise CertificationFailure(
                f"{program} canonical cases differ from the existing authority"
            )

    blastn_classes = Counter(
        row["classification"]
        for row in canonical_rows
        if row["program"] == "blastn"
    )
    if blastn_classes != Counter(
        {"EXACT_TEXT": 13, "SOURCE_UNDETERMINED_ACCEPTED": 1}
    ):
        raise CertificationFailure("BLASTN canonical classifications changed")
    sakai = canonical[("blastn", "Sakai.MG1655.megablast")]
    if sakai["contract"] != "SOURCE_UNDETERMINED_ACCEPTED":
        raise CertificationFailure("Sakai must not be represented as NCBI byte-exact")

    blastp_classes = Counter(
        row["classification"]
        for row in canonical_rows
        if row["program"] == "blastp"
    )
    if blastp_classes != Counter({"EXACT_TEXT": 9}):
        raise CertificationFailure("BLASTP exact authority changed")
    for case_id, case in tblastx_cases.items():
        row = canonical[("tblastx", case_id)]
        expected_classification = (
            "EXACT_TEXT"
            if case.contract == tblastx_audit.PARITY_CONTRACT
            else "HSP_SET_DIFF"
        )
        if row["contract"] != case.contract or row["classification"] != expected_classification:
            raise CertificationFailure(
                f"TBLASTX canonical contract differs from authority: {case_id}"
            )
    if Counter(
        row["classification"]
        for row in canonical_rows
        if row["program"] == "tblastx"
    ) != Counter({"EXACT_TEXT": 14, "HSP_SET_DIFF": 6}):
        raise CertificationFailure("TBLASTX canonical classifications changed")

    representative_keys = {(program, case_id) for program, case_id, _ in REPRESENTATIVES}
    if len(representative_keys) != 6 or not representative_keys <= set(canonical):
        raise CertificationFailure("the exact six PR 6 representatives are unavailable")
    if ("tblastx", "p11_avclpv_psclpv") in representative_keys:
        raise CertificationFailure("PR 5-only p11 repeatability selection leaked into PR 6")

    return Catalog(
        canonical_rows=canonical_rows,
        canonical=canonical,
        blastn_rows=blastn_rows,
        blastp_cases=blastp_cases,
        tblastx_cases=tblastx_cases,
        blastn_compare=blastn_compare,
        blastn_cert=blastn_cert,
        blastp_audit=blastp_audit,
        tblastx_audit=tblastx_audit,
    )


def _blastn_command_row(repo_root: Path, row: dict[str, str]) -> dict[str, str]:
    command_row = dict(row)
    command_row["query"] = str((repo_root / "LOSAT" / row["query"]).resolve())
    command_row["subject"] = str((repo_root / "LOSAT" / row["subject"]).resolve())
    return command_row


def _base_commands(
    repo_root: Path,
    catalog: Catalog,
    program: str,
    case_id: str,
    native_bin: Path,
    oracles: dict[str, Path],
    output_dir: Path,
) -> tuple[list[str], list[str]]:
    unused_ncbi = output_dir / "unused.ncbi.out"
    unused_losat = output_dir / "unused.losat.out"
    if program == "blastn":
        row = _blastn_command_row(repo_root, catalog.blastn_rows[case_id])
        return (
            catalog.blastn_compare.build_ncbi_command(
                row, str(oracles[program]), unused_ncbi
            ),
            catalog.blastn_compare.build_losat_command(row, native_bin, unused_losat),
        )
    if program == "blastp":
        return catalog.blastp_audit.build_commands(
            catalog.blastp_cases[case_id],
            oracles[program],
            native_bin,
            output_dir,
        )
    if program == "tblastx":
        return catalog.tblastx_audit.build_commands(
            catalog.tblastx_cases[case_id],
            oracles[program],
            native_bin,
            unused_ncbi,
            unused_losat,
        )
    raise CertificationFailure(f"unsupported program: {program}")


# NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-832
# ```c++
# ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
#     tabinfo.SetFields(**itr, *m_Scope, &m_ScoringMatrix);
#     tabinfo.Print();
# }
# ```
# Existing program authorities own command construction and output classification.
# This platform layer records and invokes those commands without sorting or rewriting
# LOSAT output; the PR 5 raw SHA-256 is the cross-platform gate.
def build_steps(
    repo_root: Path,
    output_dir: Path,
    catalog: Catalog,
    native_bin: Path,
    oracles: dict[str, Path],
) -> list[SearchStep]:
    steps: list[SearchStep] = []
    matrix_commands: dict[tuple[str, str], list[str]] = {}

    def append_step(
        *,
        step_id: str,
        kind: str,
        program: str,
        case_id: str,
        semantic_class: str,
        run_index: int,
        command: Sequence[str],
        output_rel: str,
    ) -> None:
        canonical_hash = catalog.canonical[(program, case_id)]["losat_sha256"]
        environment = (("RAYON_NUM_THREADS", "1"),) if program == "tblastx" else ()
        steps.append(
            SearchStep(
                execution_index=len(steps) + 1,
                step_id=step_id,
                kind=kind,
                program=program,
                case_id=case_id,
                semantic_class=semantic_class,
                run_index=run_index,
                command=tuple(command),
                environment=environment,
                output_rel=output_rel,
                final_output_rel=output_rel.removesuffix(".partial"),
                expected_losat_sha256=canonical_hash,
            )
        )

    for row in catalog.canonical_rows:
        program = row["program"]
        case_id = row["case_id"]
        rel = f"matrix/{program}/{case_id}/losat.out.partial"
        _, losat_command = _base_commands(
            repo_root,
            catalog,
            program,
            case_id,
            native_bin,
            oracles,
            output_dir / f"matrix/{program}/{case_id}",
        )
        losat_command = replace_output(losat_command, output_dir / rel)
        matrix_commands[(program, case_id)] = losat_command
        append_step(
            step_id=f"matrix:{program}:{case_id}",
            kind="matrix",
            program=program,
            case_id=case_id,
            semantic_class=row["contract"],
            run_index=1,
            command=losat_command,
            output_rel=rel,
        )

    for program, case_id, semantic_class in REPRESENTATIVES:
        rel = f"oracle/{program}/{case_id}/ncbi.out.partial"
        ncbi_command, _ = _base_commands(
            repo_root,
            catalog,
            program,
            case_id,
            native_bin,
            oracles,
            output_dir / f"oracle/{program}/{case_id}",
        )
        ncbi_command = replace_output(ncbi_command, output_dir / rel)
        append_step(
            step_id=f"oracle:{program}:{case_id}",
            kind="oracle",
            program=program,
            case_id=case_id,
            semantic_class=semantic_class,
            run_index=1,
            command=ncbi_command,
            output_rel=rel,
        )

    for program, case_id, semantic_class in REPRESENTATIVES:
        for run_index in (2, 3):
            rel = (
                f"repeatability/{program}/{case_id}/run{run_index}.out.partial"
            )
            command = replace_output(
                matrix_commands[(program, case_id)], output_dir / rel
            )
            append_step(
                step_id=f"repeatability:{program}:{case_id}:run{run_index}",
                kind="repeatability",
                program=program,
                case_id=case_id,
                semantic_class=semantic_class,
                run_index=run_index,
                command=command,
                output_rel=rel,
            )

    counts = Counter(step.kind for step in steps)
    if dict(counts) != EXPECTED_EXECUTIONS or len(steps) != 61:
        raise CertificationFailure(
            f"execution plan must be exactly 43+6+12=61, observed {dict(counts)}"
        )
    return steps


def canonical_newlines(data: bytes) -> bytes:
    return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def newline_profile(data: bytes) -> dict[str, int]:
    crlf = data.count(b"\r\n")
    without_crlf = data.replace(b"\r\n", b"")
    return {
        "crlf": crlf,
        "lf": data.count(b"\n") - crlf,
        "cr": without_crlf.count(b"\r"),
    }


def line_endings_only(left: Path, right: Path) -> bool:
    left_bytes = left.read_bytes()
    right_bytes = right.read_bytes()
    return left_bytes != right_bytes and canonical_newlines(left_bytes) == canonical_newlines(
        right_bytes
    )


def _normalized_copy(source: Path, destination: Path) -> Path:
    atomic_write_bytes(destination, canonical_newlines(source.read_bytes()))
    return destination


def classify_oracle(
    catalog: Catalog,
    step: SearchStep,
    ncbi_output: Path,
    losat_output: Path,
    diagnostic_dir: Path,
) -> dict[str, object]:
    diagnostic_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, object] = {
        "program": step.program,
        "case_id": step.case_id,
        "ncbi_sha256": sha256_path(ncbi_output),
        "losat_sha256": sha256_path(losat_output),
        "ncbi_newlines": newline_profile(ncbi_output.read_bytes()),
        "losat_newlines": newline_profile(losat_output.read_bytes()),
        "platform_diagnostic": "none",
    }
    if step.program == "blastn":
        with contextlib.redirect_stdout(io.StringIO()) as captured:
            structured = catalog.blastn_compare.compare_case(
                step.case_id, ncbi_output, losat_output, 5
            )
            raw_equal = catalog.blastn_compare.compare_raw_bytes(
                ncbi_output, losat_output
            )
        atomic_write_bytes(
            diagnostic_dir / "authority.txt", captured.getvalue().encode("utf-8")
        )
        if step.case_id != "Sakai.MG1655.megablast":
            if raw_equal:
                result["raw_authority_classification"] = "EXACT_TEXT"
                result["authority_classification"] = "EXACT_TEXT"
                return result
            if line_endings_only(ncbi_output, losat_output):
                result["raw_authority_classification"] = (
                    "STRUCTURED_EQUAL"
                    if structured["structured_equal"]
                    else "RAW_LINE_ENDING_DIFFERENCE"
                )
                result["authority_classification"] = "EXACT_DATA"
                result["platform_diagnostic"] = "LINE_ENDINGS_ONLY"
                return result
            raise CertificationFailure(
                f"BLASTN direct oracle rejected {step.case_id}"
            )

        exceptions = catalog.blastn_cert.read_source_exceptions(
            Path(catalog.blastn_cert.__file__).with_name(
                "blastn_v010_source_exceptions.tsv"
            )
        )
        exception = exceptions[step.case_id]
        try:
            authority_result = catalog.blastn_cert.certify_source_exception(
                exception, ncbi_output, losat_output
            )
        except catalog.blastn_cert.CertificationError as raw_error:
            normalized_ncbi = _normalized_copy(
                ncbi_output, diagnostic_dir / "ncbi.normalized.out"
            )
            normalized_losat = _normalized_copy(
                losat_output, diagnostic_dir / "losat.normalized.out"
            )
            if newline_profile(ncbi_output.read_bytes()) == newline_profile(
                losat_output.read_bytes()
            ):
                raise CertificationFailure(str(raw_error)) from raw_error
            try:
                authority_result = catalog.blastn_cert.certify_source_exception(
                    exception, normalized_ncbi, normalized_losat
                )
            except catalog.blastn_cert.CertificationError as normalized_error:
                raise CertificationFailure(str(normalized_error)) from normalized_error
            result["platform_diagnostic"] = "LINE_ENDINGS_ONLY"
        result["authority_classification"] = authority_result.classification
        result["differing_coord_keys"] = authority_result.differing_coord_keys
        result["field_diffs"] = authority_result.field_diffs
        return result

    if step.program == "blastp":
        classification, metrics = catalog.blastp_audit.classify_outputs(
            ncbi_output, losat_output
        )
        result["raw_authority_classification"] = classification
        result["authority_metrics"] = metrics
        if classification in catalog.blastp_audit.EXACT_CLASSIFICATIONS:
            result["authority_classification"] = classification
            return result
        if classification == "EXACT_DATA" and line_endings_only(
            ncbi_output, losat_output
        ):
            result["authority_classification"] = classification
            result["platform_diagnostic"] = "LINE_ENDINGS_ONLY"
            return result
        raise CertificationFailure(f"BLASTP direct oracle rejected {step.case_id}")

    case = catalog.tblastx_cases[step.case_id]
    classification, metrics = catalog.tblastx_audit.classify_outputs(
        ncbi_output, losat_output
    )
    result["raw_authority_classification"] = classification
    result["authority_metrics"] = metrics
    if catalog.tblastx_audit.contract_accepts(case, classification):
        result["authority_classification"] = classification
        return result
    if (
        case.contract == catalog.tblastx_audit.PARITY_CONTRACT
        and classification == "EXACT_DATA"
        and line_endings_only(ncbi_output, losat_output)
    ):
        result["authority_classification"] = classification
        result["platform_diagnostic"] = "LINE_ENDINGS_ONLY"
        return result
    raise CertificationFailure(
        f"TBLASTX direct oracle rejected {step.case_id}: {classification}"
    )


def validate_git_identity(repo_root: Path, expected_sha: str) -> str:
    if not SHA_PATTERN.fullmatch(expected_sha):
        raise CertificationFailure("expected SHA must be exactly 40 lowercase hex digits")
    head = run_capture(["git", "rev-parse", "HEAD"], repo_root)
    if head.returncode != 0:
        raise CertificationFailure("cannot resolve git HEAD")
    observed = head.stdout.decode().strip()
    if observed != expected_sha:
        raise CertificationFailure(
            f"certification SHA mismatch: expected {expected_sha}, observed {observed}"
        )
    dirty = run_capture(
        [
            "git",
            "diff",
            "--ignore-cr-at-eol",
            "--quiet",
            "HEAD",
            "--",
            "LOSAT/src",
            "LOSAT/Cargo.toml",
            "LOSAT/Cargo.lock",
            "LOSAT/.cargo/config.toml",
        ],
        repo_root,
    )
    if dirty.returncode != 0:
        raise CertificationFailure("tracked output-affecting files are dirty")
    return observed


def _tool_output(command: Sequence[str], repo_root: Path, label: str) -> str:
    completed = run_capture(command, repo_root)
    if completed.returncode != 0:
        raise CertificationFailure(f"cannot record {label} identity")
    return (completed.stdout + completed.stderr).decode("utf-8", errors="replace").strip()


def record_identity(
    repo_root: Path,
    output_dir: Path,
    expected_sha: str,
    platform_spec: PlatformSpec,
    native_bin: Path,
    archive_record: dict[str, str],
    oracles: dict[str, Path],
    canonical_manifest: Path,
) -> dict[str, object]:
    observed_os = platform.system()
    observed_machine = normalize_machine(platform.machine())
    if observed_os != platform_spec.os_name or observed_machine != platform_spec.machine:
        raise CertificationFailure(
            f"runner architecture mismatch: expected {platform_spec.os_name}/"
            f"{platform_spec.machine}, observed {observed_os}/{observed_machine}"
        )
    rustc = _tool_output(["rustc", "-vV"], repo_root, "rustc")
    cargo = _tool_output(["cargo", "-V"], repo_root, "cargo")
    host = next(
        (line.removeprefix("host: ") for line in rustc.splitlines() if line.startswith("host: ")),
        "",
    )
    if (
        f"rustc {CERT_TOOLCHAIN}" not in rustc
        or not cargo.startswith(f"cargo {CERT_TOOLCHAIN}")
        or host != platform_spec.target_triple
    ):
        raise CertificationFailure("Rust/Cargo toolchain differs from PR 6 authority")
    observed_binary_arch = binary_architecture(native_bin)
    if observed_binary_arch != platform_spec.binary_arch:
        raise CertificationFailure(
            f"LOSAT binary architecture mismatch: {observed_binary_arch}"
        )
    losat_version = _tool_output([str(native_bin), "--version"], repo_root, "LOSAT")
    oracle_identities: dict[str, object] = {}
    for program, oracle in oracles.items():
        if not oracle.is_file():
            raise CertificationFailure(f"official {program} oracle is missing: {oracle}")
        version = _tool_output([str(oracle), "-version"], repo_root, program)
        if f"{program}: 2.17.0+" not in version or "Package: blast 2.17.0" not in version:
            raise CertificationFailure(f"{program} oracle version mismatch: {version}")
        oracle_identities[program] = {
            "path": str(oracle.resolve()),
            "sha256": sha256_path(oracle),
            "version": version.splitlines(),
        }
    relevant_environment = {
        key: os.environ.get(key, "")
        for key in (
            "RUNNER_OS",
            "RUNNER_ARCH",
            "ImageOS",
            "ImageVersion",
            "PROCESSOR_ARCHITECTURE",
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
    identity: dict[str, object] = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "created_at": utc_now(),
        "source_sha": expected_sha,
        "platform_id": platform_spec.platform_id,
        "runner_label": platform_spec.runner_label,
        "os": {
            "system": observed_os,
            "machine": observed_machine,
            "platform": platform.platform(),
            "release": platform.release(),
            "version": platform.version(),
        },
        "toolchain": {
            "rustc_vV": rustc.splitlines(),
            "cargo_V": cargo,
            "python": sys.version,
            "target_triple": host,
            "profile": "release",
            "features": "default",
            "locked": True,
            "build_command": "cargo build --release --locked",
            "environment": relevant_environment,
        },
        "losat": {
            "path": str(native_bin.resolve()),
            "sha256": sha256_path(native_bin),
            "architecture": observed_binary_arch,
            "version": losat_version.splitlines(),
        },
        "oracle_archive": archive_record,
        "oracles": oracle_identities,
        "canonical_manifest": {
            "path": str(canonical_manifest.resolve()),
            "sha256": sha256_path(canonical_manifest),
            "runtime_cert_sha": RUNTIME_CERT_SHA,
            "pr5_evidence_manifest_sha256": PR5_EVIDENCE_MANIFEST_SHA256,
        },
        "execution_contract": {
            "matrix": 43,
            "oracle": 6,
            "repeatability": 12,
            "per_platform": 61,
        },
    }
    atomic_write_json(output_dir / "identity.json", identity)
    return identity


def _step_record_path(output_dir: Path, step: SearchStep) -> Path:
    safe = step.step_id.replace(":", "__")
    return output_dir / "completions" / f"{safe}.json"


def _write_state(
    output_dir: Path,
    platform_id: str,
    completed: dict[str, dict[str, object]],
    *,
    status: str,
    failure: str | None = None,
) -> None:
    by_kind = Counter(str(record["kind"]) for record in completed.values())
    state: dict[str, object] = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "updated_at": utc_now(),
        "platform_id": platform_id,
        "status": status,
        "completed_execution_count": len(completed),
        "completed_by_kind": dict(sorted(by_kind.items())),
        "completed_step_ids": sorted(completed),
        "expected_execution_count": 61,
    }
    if failure is not None:
        state["failure"] = failure
    atomic_write_json(output_dir / "state.json", state)


def _identity_resume_key(identity: dict[str, object]) -> dict[str, object]:
    return {
        "evidence_schema": identity["evidence_schema"],
        "source_sha": identity["source_sha"],
        "platform_id": identity["platform_id"],
        "runner_label": identity["runner_label"],
        "toolchain": identity["toolchain"],
        "losat_sha256": identity["losat"]["sha256"],
        "oracle_archive_md5": identity["oracle_archive"]["archive_md5"],
        "oracle_sha256": {
            program: value["sha256"] for program, value in identity["oracles"].items()
        },
        "canonical_manifest_sha256": identity["canonical_manifest"]["sha256"],
    }


def _safe_resume_path(root: Path, relative: str) -> Path:
    candidate = (root / relative).resolve()
    if not candidate.is_relative_to(root.resolve()):
        raise CertificationFailure(f"resume path escapes evidence root: {relative}")
    return candidate


def _validate_completion_shape(
    record: dict[str, object], step: SearchStep, identity: dict[str, object]
) -> None:
    expected = {
        "step_id": step.step_id,
        "execution_index": step.execution_index,
        "kind": step.kind,
        "program": step.program,
        "case_id": step.case_id,
        "run_index": step.run_index,
        "expected_losat_sha256": step.expected_losat_sha256,
        "source_sha": identity["source_sha"],
        "platform_id": identity["platform_id"],
        "losat_binary_sha256": identity["losat"]["sha256"],
    }
    for key, value in expected.items():
        if record.get(key) != value:
            raise CertificationFailure(
                f"resume completion mismatch for {step.step_id}: {key}"
            )
    if record.get("status") != "VERIFIED":
        raise CertificationFailure(f"resume step is not verified: {step.step_id}")


def import_resume(
    resume_dir: Path,
    output_dir: Path,
    identity: dict[str, object],
    steps: Sequence[SearchStep],
    catalog: Catalog,
) -> dict[str, dict[str, object]]:
    if not resume_dir.exists() or not any(resume_dir.iterdir()):
        return {}
    identity_path = resume_dir / "identity.json"
    if not identity_path.is_file():
        raise CertificationFailure("resume artifact has no identity.json")
    prior_identity = json.loads(identity_path.read_text(encoding="utf-8"))
    if _identity_resume_key(prior_identity) != _identity_resume_key(identity):
        raise CertificationFailure("resume artifact identity does not match this build")

    completed: dict[str, dict[str, object]] = {}
    for step in steps:
        prior_record_path = _step_record_path(resume_dir, step)
        if not prior_record_path.is_file():
            continue
        record = json.loads(prior_record_path.read_text(encoding="utf-8"))
        _validate_completion_shape(record, step, identity)
        prior_output = _safe_resume_path(resume_dir, str(record["output_rel"]))
        if not prior_output.is_file() or sha256_path(prior_output) != record["output_sha256"]:
            raise CertificationFailure(f"resume output hash mismatch: {step.step_id}")
        if step.kind in {"matrix", "repeatability"} and (
            record["output_sha256"] != step.expected_losat_sha256
        ):
            raise CertificationFailure(f"resume LOSAT output changed: {step.step_id}")
        current_output = _safe_resume_path(output_dir, str(record["output_rel"]))
        current_output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(prior_output, current_output)
        if step.kind == "oracle":
            losat_output = output_dir / (
                f"matrix/{step.program}/{step.case_id}/losat.out"
            )
            if not losat_output.is_file():
                raise CertificationFailure(
                    f"resume oracle lacks verified matrix run 1: {step.step_id}"
                )
            classification = classify_oracle(
                catalog,
                step,
                current_output,
                losat_output,
                output_dir / f"oracle/{step.program}/{step.case_id}/resume-diagnostic",
            )
            if classification != record.get("verification"):
                raise CertificationFailure(
                    f"resume oracle classification changed: {step.step_id}"
                )
        record["resumed_at"] = utc_now()
        record["resumed_from"] = str(resume_dir.resolve())
        atomic_write_json(_step_record_path(output_dir, step), record)
        completed[step.step_id] = record
        _write_state(
            output_dir,
            str(identity["platform_id"]),
            completed,
            status="RUNNING",
        )
    return completed


def _execute_step(
    repo_root: Path,
    output_dir: Path,
    identity: dict[str, object],
    catalog: Catalog,
    step: SearchStep,
) -> dict[str, object]:
    partial_output = output_dir / step.output_rel
    final_output = output_dir / step.final_output_rel
    partial_output.parent.mkdir(parents=True, exist_ok=True)
    if partial_output.exists():
        partial_output.unlink()
    environment = os.environ.copy()
    environment.update(dict(step.environment))
    started_at = utc_now()
    completed = subprocess.run(
        list(step.command),
        cwd=repo_root,
        env=environment,
        capture_output=True,
        text=False,
        check=False,
    )
    log_dir = partial_output.parent
    atomic_write_bytes(log_dir / "stdout.log", completed.stdout)
    atomic_write_bytes(log_dir / "stderr.log", completed.stderr)
    if completed.returncode != 0:
        raise CertificationFailure(
            f"search execution failed ({completed.returncode}): {step.step_id}"
        )
    if not partial_output.is_file():
        raise CertificationFailure(f"search produced no output: {step.step_id}")
    os.replace(partial_output, final_output)
    output_hash = sha256_path(final_output)
    verification: dict[str, object]
    if step.kind in {"matrix", "repeatability"}:
        if output_hash != step.expected_losat_sha256:
            raise CertificationFailure(
                f"CROSS_PLATFORM_LOSAT_DIVERGENCE: {step.program} {step.case_id} "
                f"expected {step.expected_losat_sha256}, observed {output_hash}"
            )
        verification = {
            "classification": "CANONICAL_PR5_RAW_BYTES",
            "canonical_sha256": step.expected_losat_sha256,
        }
    else:
        losat_output = output_dir / f"matrix/{step.program}/{step.case_id}/losat.out"
        if not losat_output.is_file():
            raise CertificationFailure(
                f"direct oracle lacks verified matrix run 1: {step.step_id}"
            )
        verification = classify_oracle(
            catalog,
            step,
            final_output,
            losat_output,
            final_output.parent / "diagnostic",
        )
    record: dict[str, object] = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "status": "VERIFIED",
        "step_id": step.step_id,
        "execution_index": step.execution_index,
        "kind": step.kind,
        "program": step.program,
        "case_id": step.case_id,
        "semantic_class": step.semantic_class,
        "run_index": step.run_index,
        "started_at": started_at,
        "completed_at": utc_now(),
        "source_sha": identity["source_sha"],
        "platform_id": identity["platform_id"],
        "losat_binary_sha256": identity["losat"]["sha256"],
        "expected_losat_sha256": step.expected_losat_sha256,
        "output_rel": step.final_output_rel,
        "output_sha256": output_hash,
        "returncode": completed.returncode,
        "verification": verification,
    }
    atomic_write_json(_step_record_path(output_dir, step), record)
    return record


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(materialized)
    os.replace(temporary, path)


def finalize_evidence(
    output_dir: Path,
    identity: dict[str, object],
    steps: Sequence[SearchStep],
    completed: dict[str, dict[str, object]],
    executed_this_attempt: int,
) -> None:
    if len(completed) != 61:
        raise CertificationFailure(
            f"cannot finalize incomplete platform evidence: {len(completed)}/61"
        )
    matrix_rows = []
    oracle_rows = []
    for step in steps:
        record = completed[step.step_id]
        if step.kind == "matrix":
            matrix_rows.append(
                {
                    "program": step.program,
                    "case_id": step.case_id,
                    "expected_sha256": step.expected_losat_sha256,
                    "observed_sha256": record["output_sha256"],
                    "status": "EXACT_BYTES",
                }
            )
        elif step.kind == "oracle":
            verification = record["verification"]
            oracle_rows.append(
                {
                    "program": step.program,
                    "case_id": step.case_id,
                    "semantic_class": step.semantic_class,
                    "classification": verification["authority_classification"],
                    "platform_diagnostic": verification["platform_diagnostic"],
                    "ncbi_sha256": record["output_sha256"],
                    "losat_sha256": step.expected_losat_sha256,
                    "status": "ACCEPTED",
                }
            )
    repeatability_rows = []
    for program, case_id, semantic_class in REPRESENTATIVES:
        run_hashes = [
            completed[f"matrix:{program}:{case_id}"]["output_sha256"],
            completed[f"repeatability:{program}:{case_id}:run2"]["output_sha256"],
            completed[f"repeatability:{program}:{case_id}:run3"]["output_sha256"],
        ]
        repeatability_rows.append(
            {
                "program": program,
                "case_id": case_id,
                "semantic_class": semantic_class,
                "sha256_run1": run_hashes[0],
                "sha256_run2": run_hashes[1],
                "sha256_run3": run_hashes[2],
                "status": "REPEATABLE" if len(set(run_hashes)) == 1 else "DIFFER",
            }
        )
    if any(row["status"] != "REPEATABLE" for row in repeatability_rows):
        raise CertificationFailure("repeatability aggregation failed")
    write_tsv(
        output_dir / "cross_platform_equalities.tsv",
        matrix_rows,
        ["program", "case_id", "expected_sha256", "observed_sha256", "status"],
    )
    write_tsv(
        output_dir / "direct_ncbi_checks.tsv",
        oracle_rows,
        [
            "program",
            "case_id",
            "semantic_class",
            "classification",
            "platform_diagnostic",
            "ncbi_sha256",
            "losat_sha256",
            "status",
        ],
    )
    write_tsv(
        output_dir / "repeatability.tsv",
        repeatability_rows,
        [
            "program",
            "case_id",
            "semantic_class",
            "sha256_run1",
            "sha256_run2",
            "sha256_run3",
            "status",
        ],
    )
    summary = {
        "decision": "PLATFORM_NATIVE_CERTIFIED",
        "platform_id": identity["platform_id"],
        "source_sha": identity["source_sha"],
        "completed_at": utc_now(),
        "search_executions": 61,
        "executed_this_attempt": executed_this_attempt,
        "resumed_verified_executions": 61 - executed_this_attempt,
        "cross_platform_losat": {"total": 43, "all_exact_bytes": True},
        "direct_ncbi": {"total": 6, "all_authority_accepted": True},
        "repeatability": {
            "representatives": 6,
            "runs_per_representative": 3,
            "additional_executions": 12,
            "all_repeatable": True,
        },
        "excluded": [
            "Linux matrix",
            "Wasm",
            "performance benchmark",
            "soak testing",
            "expanded NCBI matrix",
        ],
    }
    atomic_write_json(output_dir / "summary.json", summary)
    _write_state(
        output_dir,
        str(identity["platform_id"]),
        completed,
        status="CERTIFIED",
    )
    manifest_rows = []
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "FINAL_EVIDENCE_MANIFEST.sha256":
            manifest_rows.append(
                f"{sha256_path(path)}  {path.relative_to(output_dir).as_posix()}\n"
            )
    atomic_write_bytes(
        output_dir / "FINAL_EVIDENCE_MANIFEST.sha256",
        "".join(manifest_rows).encode("utf-8"),
    )


def certify(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    preexisting = {path.name for path in output_dir.iterdir()}
    if preexisting - {"workflow-started.txt"}:
        raise CertificationFailure(f"output directory must be new or empty: {output_dir}")
    platform_spec = PLATFORMS[args.platform_id]
    initial_state = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "created_at": utc_now(),
        "platform_id": platform_spec.platform_id,
        "status": "INITIALIZING",
        "expected_execution_count": 61,
    }
    atomic_write_json(output_dir / "state.json", initial_state)
    certified_sha = validate_git_identity(repo_root, args.expected_sha)
    archive_record = verify_archive(
        platform_spec, args.oracle_archive.resolve(), args.oracle_checksum.resolve()
    )
    oracles = {
        "blastn": args.blastn_oracle.resolve(),
        "blastp": args.blastp_oracle.resolve(),
        "tblastx": args.tblastx_oracle.resolve(),
    }
    canonical_manifest = (
        repo_root / "LOSAT" / "tests" / "platform_native_v010_canonical.tsv"
    )
    identity = record_identity(
        repo_root,
        output_dir,
        certified_sha,
        platform_spec,
        args.native_bin.resolve(),
        archive_record,
        oracles,
        canonical_manifest,
    )
    catalog = load_catalog(repo_root)
    steps = build_steps(
        repo_root,
        output_dir,
        catalog,
        args.native_bin.resolve(),
        oracles,
    )
    atomic_write_json(
        output_dir / "command_plan.json",
        {
            "source_sha": certified_sha,
            "platform_id": platform_spec.platform_id,
            "created_at": utc_now(),
            "execution_count": len(steps),
            "executions": [asdict(step) for step in steps],
        },
    )
    completed: dict[str, dict[str, object]] = {}
    if args.resume_dir is not None:
        completed = import_resume(
            args.resume_dir.resolve(), output_dir, identity, steps, catalog
        )
    _write_state(output_dir, platform_spec.platform_id, completed, status="RUNNING")
    executed_this_attempt = 0
    for step in steps:
        if step.step_id in completed:
            print(f"[platform-certification] resume {step.execution_index}/61 {step.step_id}")
            continue
        print(f"[platform-certification] run {step.execution_index}/61 {step.step_id}", flush=True)
        record = _execute_step(repo_root, output_dir, identity, catalog, step)
        completed[step.step_id] = record
        executed_this_attempt += 1
        _write_state(
            output_dir, platform_spec.platform_id, completed, status="RUNNING"
        )
    finalize_evidence(output_dir, identity, steps, completed, executed_this_attempt)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    archive_parser = subparsers.add_parser(
        "verify-archive", help="verify a pinned official NCBI archive"
    )
    archive_parser.add_argument("--platform-id", choices=sorted(PLATFORMS), required=True)
    archive_parser.add_argument("--archive", type=Path, required=True)
    archive_parser.add_argument("--checksum", type=Path, required=True)

    certify_parser = subparsers.add_parser(
        "certify", help="run or resume the exact 61-execution platform gate"
    )
    certify_parser.add_argument("--repo-root", type=Path, default=repo_root)
    certify_parser.add_argument("--expected-sha", required=True)
    certify_parser.add_argument("--platform-id", choices=sorted(PLATFORMS), required=True)
    certify_parser.add_argument("--output-dir", type=Path, required=True)
    certify_parser.add_argument("--resume-dir", type=Path)
    certify_parser.add_argument("--native-bin", type=Path, required=True)
    certify_parser.add_argument("--oracle-archive", type=Path, required=True)
    certify_parser.add_argument("--oracle-checksum", type=Path, required=True)
    certify_parser.add_argument("--blastn-oracle", type=Path, required=True)
    certify_parser.add_argument("--blastp-oracle", type=Path, required=True)
    certify_parser.add_argument("--tblastx-oracle", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = getattr(args, "output_dir", None)
    try:
        if args.mode == "verify-archive":
            record = verify_archive(
                PLATFORMS[args.platform_id], args.archive.resolve(), args.checksum.resolve()
            )
            print(json.dumps(record, sort_keys=True))
        else:
            certify(args)
            print("PLATFORM_NATIVE_CERTIFIED")
        return 0
    except (CertificationFailure, OSError, ValueError) as error:
        if output_dir is not None:
            resolved = output_dir.resolve()
            resolved.mkdir(parents=True, exist_ok=True)
            atomic_write_json(
                resolved / "failure.json",
                {"failed_at": utc_now(), "error": str(error), "status": "FAILED"},
            )
            state_path = resolved / "state.json"
            if state_path.is_file():
                state = json.loads(state_path.read_text(encoding="utf-8"))
                state["status"] = "FAILED"
                state["failure"] = str(error)
                state["updated_at"] = utc_now()
                atomic_write_json(state_path, state)
        print(f"PLATFORM_CERTIFICATION_FAILED: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
