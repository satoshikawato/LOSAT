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
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import subprocess
import sys
import tarfile
from types import ModuleType
from typing import Iterable, Sequence


EXPECTED_COUNTS = {"blastn": 14, "blastp": 9, "tblastx": 20}
EXPECTED_EXECUTIONS = {"matrix": 43, "oracle": 6, "repeatability": 12}
EXPECTED_STAGED_FIXTURE_COUNT = 34
RUNTIME_CERT_SHA = "5845d22ed9842449628a647f29b8c6762511ca59"
PR5_EVIDENCE_MANIFEST_SHA256 = (
    "b9fc98a376d2849274c86b4e4769d2ee38b76025adbb4d63d3da4a5e3e7cdb5c"
)
CERT_TOOLCHAIN = "1.92.0"
EVIDENCE_SCHEMA = 1
NATIVE_AUTHORITY_VERSION = "ncbi-platform-variance-v0.1.0"
SAFE_CASE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
SHA_PATTERN = re.compile(r"^[0-9a-f]{40}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
MD5_PATTERN = re.compile(r"^[0-9a-f]{32}$")
HISTORICAL_LEXICAL_ROOT = "/tmp/losat-pr5-runtime-cert-5845d22/LOSAT"
INPUT_OPTIONS = {
    "blastn": {
        "query": ("-q", "-query"),
        "subject": ("-s", "-subject"),
    },
    "blastp": {
        "query": ("--query", "-query"),
        "subject": ("--subject", "-subject"),
    },
    "tblastx": {
        "query": ("--query", "-query"),
        "subject": ("--subject", "-subject"),
    },
}


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
            f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B/2.17.0/{quoted}"
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


@dataclass(frozen=True)
class StagedFixture:
    source_repo_relative: str
    source_path: str
    lexical_target_path: str
    physical_target_path: str
    sha256: str


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


@dataclass(frozen=True)
class NativeAuthority:
    path: Path
    authority_version: str
    file_sha256: str
    document: dict[str, object]
    representatives: dict[tuple[str, str], dict[str, object]]
    platforms: dict[str, dict[str, object]]
    outputs: dict[tuple[str, str, str], dict[str, object]]
    diagnostics: dict[tuple[str, str, str], dict[str, object]]


class CertificationFailure(RuntimeError):
    """Platform evidence is outside the frozen PR 6 contract."""


def native_authority_miss(message: str) -> CertificationFailure:
    return CertificationFailure(f"NATIVE_NCBI_AUTHORITY_MISS: {message}")


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


def _strict_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise CertificationFailure(f"duplicate authority JSON key: {key}")
        result[key] = value
    return result


def _object(value: object, context: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise CertificationFailure(f"{context} must be an object")
    return value


def _array(value: object, context: str) -> list[object]:
    if not isinstance(value, list):
        raise CertificationFailure(f"{context} must be an array")
    return value


def _exact_fields(value: object, expected: set[str], context: str) -> dict[str, object]:
    record = _object(value, context)
    observed = set(record)
    if observed != expected:
        unknown = sorted(observed - expected)
        missing = sorted(expected - observed)
        raise CertificationFailure(
            f"{context} fields changed: unknown={unknown}, missing={missing}"
        )
    return record


def _string(value: object, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise CertificationFailure(f"{context} must be a non-empty string")
    return value


def _integer(value: object, context: str, *, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise CertificationFailure(f"{context} must be an integer >= {minimum}")
    return value


def _boolean(value: object, context: str) -> bool:
    if type(value) is not bool:
        raise CertificationFailure(f"{context} must be a boolean")
    return value


def _sha256(value: object, context: str) -> str:
    digest = _string(value, context)
    if not SHA256_PATTERN.fullmatch(digest):
        raise CertificationFailure(f"{context} must be 64 lowercase hex digits")
    return digest


def _sha40(value: object, context: str) -> str:
    digest = _string(value, context)
    if not SHA_PATTERN.fullmatch(digest):
        raise CertificationFailure(f"{context} must be 40 lowercase hex digits")
    return digest


def _no_wildcard(value: object, context: str) -> str:
    identity = _string(value, context)
    if any(character in identity for character in "*?[]"):
        raise CertificationFailure(f"{context} contains a wildcard")
    return identity


def _validate_authority_newline_profile(
    value: object, context: str
) -> dict[str, object]:
    profile = _exact_fields(
        value,
        {"lf_count", "crlf_count", "bare_lf_count", "bare_cr_count"},
        context,
    )
    counts = {key: _integer(profile[key], f"{context}.{key}") for key in profile}
    if counts["lf_count"] != counts["crlf_count"] + counts["bare_lf_count"]:
        raise CertificationFailure(f"{context} has an inconsistent LF count")
    return profile


def _validate_authority_input(value: object, context: str, role: str) -> None:
    record = _exact_fields(
        value,
        {
            "role",
            "option",
            "repository_relative",
            "lexical_path",
            "sha256",
            "byte_length",
            "newline_profile",
        },
        context,
    )
    if record["role"] != role:
        raise CertificationFailure(f"{context}.role must be {role}")
    if record["option"] not in {"-query", "-subject"}:
        raise CertificationFailure(f"{context}.option is not an NCBI input option")
    repository_relative = _no_wildcard(
        record["repository_relative"], f"{context}.repository_relative"
    )
    if not repository_relative.startswith("LOSAT/tests/fasta/"):
        raise CertificationFailure(f"{context} is outside the controlled fixtures")
    lexical_path = _no_wildcard(record["lexical_path"], f"{context}.lexical_path")
    expected_lexical = (
        f"{HISTORICAL_LEXICAL_ROOT}/{repository_relative.removeprefix('LOSAT/')}"
    )
    if lexical_path != expected_lexical:
        raise CertificationFailure(f"{context} lexical identity changed")
    _sha256(record["sha256"], f"{context}.sha256")
    _integer(record["byte_length"], f"{context}.byte_length", minimum=1)
    newline = _validate_authority_newline_profile(
        record["newline_profile"], f"{context}.newline_profile"
    )
    if newline["crlf_count"] != 0 or newline["bare_cr_count"] != 0:
        raise CertificationFailure(f"{context} must retain exact LF Git-blob bytes")


def _validate_native_authority_document(document: dict[str, object]) -> NativeAuthority:
    top = _exact_fields(
        document,
        {
            "authority_version",
            "schema_version",
            "normative_acceptance",
            "diagnostic_metadata",
            "characterization_provenance",
        },
        "authority",
    )
    if top["authority_version"] != NATIVE_AUTHORITY_VERSION:
        raise CertificationFailure(
            f"unknown native authority version: {top['authority_version']}"
        )
    if top["schema_version"] != 1:
        raise CertificationFailure("unknown native authority schema version")
    normative = _exact_fields(
        top["normative_acceptance"],
        {"ncbi_release", "representatives", "platforms"},
        "normative_acceptance",
    )
    if normative["ncbi_release"] != "2.17.0":
        raise CertificationFailure("native authority NCBI release changed")

    expected_representatives = {
        (program, case_id) for program, case_id, _ in REPRESENTATIVES
    }
    representatives: dict[tuple[str, str], dict[str, object]] = {}
    for index, value in enumerate(
        _array(normative["representatives"], "normative_acceptance.representatives")
    ):
        context = f"normative_acceptance.representatives[{index}]"
        record = _exact_fields(
            value,
            {
                "program",
                "case_id",
                "losat_contract",
                "losat_parity_class",
                "query",
                "subject",
                "ordered_argv_template",
                "authoritative_outfmt",
                "working_directory_contract",
                "lexical_fixture_root",
                "shell",
                "environment_overrides",
            },
            context,
        )
        program = _no_wildcard(record["program"], f"{context}.program")
        case_id = _no_wildcard(record["case_id"], f"{context}.case_id")
        key = (program, case_id)
        if key in representatives:
            raise CertificationFailure(f"duplicate native representative: {key}")
        if key not in expected_representatives:
            raise CertificationFailure(f"unexpected native representative: {key}")
        if record["losat_contract"] not in {
            "EXACT_TEXT",
            "SOURCE_UNDETERMINED_ACCEPTED",
            "source_defined_exact",
            "parity",
            "approved_db_gencode_deviation",
        }:
            raise CertificationFailure(f"{context}.losat_contract is invalid")
        if record["losat_parity_class"] not in {
            "EXACT_TEXT",
            "SOURCE_UNDETERMINED_ACCEPTED",
            "HSP_SET_DIFF",
        }:
            raise CertificationFailure(f"{context}.losat_parity_class is invalid")
        _validate_authority_input(record["query"], f"{context}.query", "query")
        _validate_authority_input(record["subject"], f"{context}.subject", "subject")
        argv = _array(
            record["ordered_argv_template"], f"{context}.ordered_argv_template"
        )
        if not argv or any(not isinstance(item, str) or not item for item in argv):
            raise CertificationFailure(f"{context}.ordered_argv_template is invalid")
        if argv[0] != "{oracle}" or argv.count("{oracle}") != 1:
            raise CertificationFailure(
                f"{context} must have one exact oracle placeholder"
            )
        if argv.count("{output}") != 1 or argv[-1] != "{output}":
            raise CertificationFailure(
                f"{context} must have one exact output placeholder"
            )
        for item in argv:
            if ("{" in item or "}" in item) and item not in {"{oracle}", "{output}"}:
                raise CertificationFailure(f"{context} contains an unknown placeholder")
            if item not in {"{oracle}", "{output}"}:
                _no_wildcard(item, f"{context}.ordered_argv_template")
        if (
            argv.count("-outfmt") != 1
            or argv[argv.index("-outfmt") + 1] != record["authoritative_outfmt"]
        ):
            raise CertificationFailure(f"{context} authoritative outfmt changed")
        if record["authoritative_outfmt"] not in {"6", "7"}:
            raise CertificationFailure(f"{context} authoritative outfmt is invalid")
        for role in ("query", "subject"):
            input_record = _object(record[role], f"{context}.{role}")
            option = str(input_record["option"])
            if (
                argv.count(option) != 1
                or argv[argv.index(option) + 1] != input_record["lexical_path"]
            ):
                raise CertificationFailure(f"{context} {role} argv identity changed")
        if record["working_directory_contract"] != "repository_root":
            raise CertificationFailure(f"{context} working-directory contract changed")
        if record["lexical_fixture_root"] != HISTORICAL_LEXICAL_ROOT:
            raise CertificationFailure(f"{context} lexical fixture root changed")
        if _boolean(record["shell"], f"{context}.shell") is not False:
            raise CertificationFailure(f"{context} must use shell=false")
        if record["environment_overrides"] != {}:
            raise CertificationFailure(f"{context} native environment must be empty")
        representatives[key] = record
    if set(representatives) != expected_representatives or len(representatives) != 6:
        raise CertificationFailure(
            "native authority must contain exactly six representatives"
        )

    expected_platforms = set(PLATFORMS)
    platforms: dict[str, dict[str, object]] = {}
    outputs: dict[tuple[str, str, str], dict[str, object]] = {}
    for index, value in enumerate(
        _array(normative["platforms"], "normative_acceptance.platforms")
    ):
        context = f"normative_acceptance.platforms[{index}]"
        record = _exact_fields(
            value,
            {
                "platform_id",
                "platform_system",
                "normalized_machine",
                "architecture",
                "archive",
                "executables",
                "expected_outputs",
            },
            context,
        )
        platform_id = _no_wildcard(record["platform_id"], f"{context}.platform_id")
        if platform_id in platforms:
            raise CertificationFailure(f"duplicate native platform: {platform_id}")
        if platform_id not in expected_platforms:
            raise CertificationFailure(f"unknown native platform: {platform_id}")
        spec = PLATFORMS[platform_id]
        if (
            record["platform_system"] != spec.os_name
            or record["normalized_machine"] != spec.machine
            or record["architecture"] != spec.binary_arch
        ):
            raise CertificationFailure(f"{context} platform sanity identity changed")
        archive = _exact_fields(
            record["archive"],
            {"release", "filename", "checksum", "retrieval_url"},
            f"{context}.archive",
        )
        if archive["release"] != "2.17.0":
            raise CertificationFailure(f"{context} archive release changed")
        if (
            _no_wildcard(archive["filename"], f"{context}.archive.filename")
            != spec.archive_name
        ):
            raise CertificationFailure(f"{context} archive filename changed")
        checksum = _exact_fields(
            archive["checksum"], {"algorithm", "value"}, f"{context}.archive.checksum"
        )
        if checksum["algorithm"] != "md5" or not MD5_PATTERN.fullmatch(
            str(checksum["value"])
        ):
            raise CertificationFailure(
                f"{context} archive checksum identity is invalid"
            )
        if checksum["value"] != spec.archive_md5:
            raise CertificationFailure(f"{context} archive checksum changed")
        retrieval_url = _no_wildcard(
            archive["retrieval_url"], f"{context}.archive.retrieval_url"
        )
        if retrieval_url != spec.archive_url:
            raise CertificationFailure(
                f"{context} archive retrieval provenance changed"
            )

        executables: dict[str, dict[str, object]] = {}
        for executable_index, executable_value in enumerate(
            _array(record["executables"], f"{context}.executables")
        ):
            executable_context = f"{context}.executables[{executable_index}]"
            executable = _exact_fields(
                executable_value,
                {"program", "sha256", "architecture"},
                executable_context,
            )
            program = _no_wildcard(
                executable["program"], f"{executable_context}.program"
            )
            if program in executables or program not in {"blastn", "blastp", "tblastx"}:
                raise CertificationFailure(
                    f"duplicate or invalid executable: {program}"
                )
            _sha256(executable["sha256"], f"{executable_context}.sha256")
            if executable["architecture"] != spec.binary_arch:
                raise CertificationFailure(f"{executable_context} architecture changed")
            executables[program] = executable
        if set(executables) != {"blastn", "blastp", "tblastx"}:
            raise CertificationFailure(f"{context} executable set changed")

        platform_output_keys = set()
        for output_index, output_value in enumerate(
            _array(record["expected_outputs"], f"{context}.expected_outputs")
        ):
            output_context = f"{context}.expected_outputs[{output_index}]"
            output = _exact_fields(
                output_value,
                {
                    "program",
                    "case_id",
                    "raw_sha256",
                    "byte_length",
                    "newline_profile",
                    "data_row_count",
                },
                output_context,
            )
            key = (
                platform_id,
                _no_wildcard(output["program"], f"{output_context}.program"),
                _no_wildcard(output["case_id"], f"{output_context}.case_id"),
            )
            representative_key = key[1:]
            if key in outputs or representative_key not in representatives:
                raise CertificationFailure(
                    f"duplicate or unexpected output tuple: {key}"
                )
            _sha256(output["raw_sha256"], f"{output_context}.raw_sha256")
            _integer(output["byte_length"], f"{output_context}.byte_length", minimum=1)
            _validate_authority_newline_profile(
                output["newline_profile"], f"{output_context}.newline_profile"
            )
            _integer(
                output["data_row_count"], f"{output_context}.data_row_count", minimum=1
            )
            outputs[key] = output
            platform_output_keys.add(representative_key)
        if platform_output_keys != expected_representatives:
            raise CertificationFailure(
                f"{context} must contain all six exact output tuples"
            )
        platforms[platform_id] = record
    if set(platforms) != expected_platforms or len(outputs) != 18:
        raise CertificationFailure(
            "native authority must contain exactly 3 platforms/18 outputs"
        )

    diagnostics: dict[tuple[str, str, str], dict[str, object]] = {}
    allowed_native_classes = {
        "EXACT_TO_RETAINED_LINUX_REFERENCE",
        "OUTPUT_NEWLINE_ONLY",
        "ROW_ORDER_VARIANCE",
        "HSP_STATISTIC_VARIANCE",
    }
    for index, value in enumerate(
        _array(top["diagnostic_metadata"], "diagnostic_metadata")
    ):
        context = f"diagnostic_metadata[{index}]"
        record = _exact_fields(
            value,
            {
                "platform_id",
                "program",
                "case_id",
                "native_ncbi_reference_class",
                "lf_diagnostic_sha256",
                "retained_linux_raw_sha256",
                "retained_linux_row_count",
                "row_key_multiset_equal",
                "row_key_order_equal",
                "complete_row_order_equal",
                "differing_row_count",
                "runtime_enforced",
            },
            context,
        )
        key = (
            _no_wildcard(record["platform_id"], f"{context}.platform_id"),
            _no_wildcard(record["program"], f"{context}.program"),
            _no_wildcard(record["case_id"], f"{context}.case_id"),
        )
        if key in diagnostics or key not in outputs:
            raise CertificationFailure(
                f"duplicate or unexpected diagnostic tuple: {key}"
            )
        if record["native_ncbi_reference_class"] not in allowed_native_classes:
            raise CertificationFailure(f"{context} native reference class is invalid")
        _sha256(record["lf_diagnostic_sha256"], f"{context}.lf_diagnostic_sha256")
        _sha256(
            record["retained_linux_raw_sha256"],
            f"{context}.retained_linux_raw_sha256",
        )
        _integer(
            record["retained_linux_row_count"],
            f"{context}.retained_linux_row_count",
            minimum=1,
        )
        _integer(record["differing_row_count"], f"{context}.differing_row_count")
        for field in (
            "row_key_multiset_equal",
            "row_key_order_equal",
            "complete_row_order_equal",
        ):
            _boolean(record[field], f"{context}.{field}")
        if (
            _boolean(record["runtime_enforced"], f"{context}.runtime_enforced")
            is not False
        ):
            raise CertificationFailure(
                f"{context} diagnostic metadata became normative"
            )
        diagnostics[key] = record
    if set(diagnostics) != set(outputs) or len(diagnostics) != 18:
        raise CertificationFailure(
            "native authority must contain exactly 18 diagnostics"
        )

    provenance = _exact_fields(
        top["characterization_provenance"],
        {
            "runtime_enforced",
            "characterization_run_id",
            "candidate_sha",
            "analysis_fix_sha",
            "corrected_replay_sha256",
            "pr6_v3_sha",
            "pr6_v3_run_id",
            "retained_linux_lineage",
            "execution_counts",
            "platform_runs",
            "rich_group_counts",
            "rich_groups",
        },
        "characterization_provenance",
    )
    if (
        _boolean(
            provenance["runtime_enforced"],
            "characterization_provenance.runtime_enforced",
        )
        is not False
    ):
        raise CertificationFailure("characterization provenance became normative")
    if (
        _integer(
            provenance["characterization_run_id"],
            "characterization_provenance.characterization_run_id",
            minimum=1,
        )
        != 33582928402
    ):
        raise CertificationFailure("characterization run provenance changed")
    if (
        _sha40(provenance["candidate_sha"], "characterization_provenance.candidate_sha")
        != "be67253156ce10b852166709e9519ab54709fd80"
    ):
        raise CertificationFailure("characterization candidate provenance changed")
    if (
        _sha40(
            provenance["analysis_fix_sha"],
            "characterization_provenance.analysis_fix_sha",
        )
        != "cfd213f57691607429c4e21edbf314899310d805"
    ):
        raise CertificationFailure("characterization analyzer provenance changed")
    if (
        _sha256(
            provenance["corrected_replay_sha256"],
            "characterization_provenance.corrected_replay_sha256",
        )
        != "37e95cf9252ede022cc0d3f3f6a43f4776f78b3e31febd2a76bf2d2d667be2db"
    ):
        raise CertificationFailure("characterization replay provenance changed")
    if (
        _sha40(provenance["pr6_v3_sha"], "characterization_provenance.pr6_v3_sha")
        != "0e29e2201e2d1b03124c0b9d6698a81bfed8cec0"
    ):
        raise CertificationFailure("PR 6 V3 provenance changed")
    if (
        _integer(
            provenance["pr6_v3_run_id"],
            "characterization_provenance.pr6_v3_run_id",
            minimum=1,
        )
        != 33503928773
    ):
        raise CertificationFailure("PR 6 V3 run provenance changed")
    linux = _exact_fields(
        provenance["retained_linux_lineage"],
        {"runtime_cert_sha", "evidence_manifest_sha256"},
        "characterization_provenance.retained_linux_lineage",
    )
    if (
        _sha40(linux["runtime_cert_sha"], "retained_linux_lineage.runtime_cert_sha")
        != RUNTIME_CERT_SHA
    ):
        raise CertificationFailure("retained Linux runtime provenance changed")
    if (
        _sha256(
            linux["evidence_manifest_sha256"],
            "retained_linux_lineage.evidence_manifest_sha256",
        )
        != PR5_EVIDENCE_MANIFEST_SHA256
    ):
        raise CertificationFailure("retained Linux evidence provenance changed")
    counts = _exact_fields(
        provenance["execution_counts"],
        {
            "authoritative_ncbi",
            "diagnostic_ncbi",
            "total_ncbi",
            "losat",
            "linux_ncbi_reruns",
        },
        "characterization_provenance.execution_counts",
    )
    if counts != {
        "authoritative_ncbi": 18,
        "diagnostic_ncbi": 18,
        "total_ncbi": 36,
        "losat": 0,
        "linux_ncbi_reruns": 0,
    }:
        raise CertificationFailure("characterization execution accounting changed")
    platform_run_ids = set()
    runner_fields = {
        "ImageOS",
        "ImageVersion",
        "RUNNER_ARCH",
        "RUNNER_OS",
        "machine",
        "platform",
        "release",
        "system",
        "version",
    }
    for index, value in enumerate(
        _array(provenance["platform_runs"], "characterization_provenance.platform_runs")
    ):
        context = f"characterization_provenance.platform_runs[{index}]"
        record = _exact_fields(
            value,
            {
                "platform_id",
                "runner_label",
                "runner",
                "fixture_manifest_sha256",
                "identity_sha256",
            },
            context,
        )
        platform_id = _no_wildcard(record["platform_id"], f"{context}.platform_id")
        if platform_id in platform_run_ids or platform_id not in expected_platforms:
            raise CertificationFailure(
                f"duplicate or unknown provenance platform: {platform_id}"
            )
        _no_wildcard(record["runner_label"], f"{context}.runner_label")
        runner = _exact_fields(record["runner"], runner_fields, f"{context}.runner")
        for key, item in runner.items():
            _string(item, f"{context}.runner.{key}")
        _sha256(record["fixture_manifest_sha256"], f"{context}.fixture_manifest_sha256")
        _sha256(record["identity_sha256"], f"{context}.identity_sha256")
        platform_run_ids.add(platform_id)
    if platform_run_ids != expected_platforms:
        raise CertificationFailure("characterization platform provenance changed")
    rich_counts = _exact_fields(
        provenance["rich_group_counts"],
        {"exhaustive", "authoritative_seeded", "rich_only"},
        "characterization_provenance.rich_group_counts",
    )
    if rich_counts != {"exhaustive": 11, "authoritative_seeded": 7, "rich_only": 4}:
        raise CertificationFailure("rich diagnostic group counts changed")
    rich_case_counts: Counter[tuple[str, str]] = Counter()
    seeded = 0
    rich_only = 0
    alignment_fields = {
        "program",
        "case_id",
        "query_id",
        "subject_id",
        "qstart",
        "qend",
        "sstart",
        "send",
        "score",
        "raw_score",
        "length",
        "evalue",
        "bitscore",
    }
    for index, value in enumerate(
        _array(provenance["rich_groups"], "characterization_provenance.rich_groups")
    ):
        context = f"characterization_provenance.rich_groups[{index}]"
        group = _exact_fields(
            value,
            {
                "alignment_identity",
                "authoritative_seeded",
                "rich_only",
                "platform_representation_multisets",
            },
            context,
        )
        identity = _exact_fields(
            group["alignment_identity"],
            alignment_fields,
            f"{context}.alignment_identity",
        )
        for key, item in identity.items():
            _no_wildcard(item, f"{context}.alignment_identity.{key}")
        case_key = (str(identity["program"]), str(identity["case_id"]))
        if case_key not in representatives:
            raise CertificationFailure(f"{context} has an unknown representative")
        authoritative_seeded = _boolean(
            group["authoritative_seeded"], f"{context}.authoritative_seeded"
        )
        is_rich_only = _boolean(group["rich_only"], f"{context}.rich_only")
        if authoritative_seeded == is_rich_only:
            raise CertificationFailure(f"{context} rich scope flags are invalid")
        seeded += authoritative_seeded
        rich_only += is_rich_only
        rich_case_counts[case_key] += 1
        multiset_platforms = set()
        for multiset_index, multiset_value in enumerate(
            _array(
                group["platform_representation_multisets"],
                f"{context}.platform_representation_multisets",
            )
        ):
            multiset_context = (
                f"{context}.platform_representation_multisets[{multiset_index}]"
            )
            multiset = _exact_fields(
                multiset_value,
                {
                    "platform_id",
                    "representation_multiset_sha256",
                    "representations",
                    "row_count",
                    "row_indices",
                },
                multiset_context,
            )
            platform_id = _no_wildcard(
                multiset["platform_id"], f"{multiset_context}.platform_id"
            )
            if (
                platform_id in multiset_platforms
                or platform_id not in expected_platforms
            ):
                raise CertificationFailure(f"{multiset_context} platform set changed")
            _sha256(
                multiset["representation_multiset_sha256"],
                f"{multiset_context}.representation_multiset_sha256",
            )
            row_count = _integer(
                multiset["row_count"], f"{multiset_context}.row_count", minimum=1
            )
            row_indices = _array(
                multiset["row_indices"], f"{multiset_context}.row_indices"
            )
            if len(row_indices) != row_count or any(
                type(item) is not int or item < 1 for item in row_indices
            ):
                raise CertificationFailure(f"{multiset_context} row identity changed")
            multiplicity = 0
            for representation_index, representation_value in enumerate(
                _array(
                    multiset["representations"], f"{multiset_context}.representations"
                )
            ):
                representation_context = (
                    f"{multiset_context}.representations[{representation_index}]"
                )
                representation = _exact_fields(
                    representation_value,
                    {"multiplicity", "btop_sha256", "qseq_sha256", "sseq_sha256"},
                    representation_context,
                )
                multiplicity += _integer(
                    representation["multiplicity"],
                    f"{representation_context}.multiplicity",
                    minimum=1,
                )
                for field in ("btop_sha256", "qseq_sha256", "sseq_sha256"):
                    _sha256(representation[field], f"{representation_context}.{field}")
            if multiplicity != row_count:
                raise CertificationFailure(f"{multiset_context} multiplicity changed")
            multiset_platforms.add(platform_id)
        if multiset_platforms != expected_platforms:
            raise CertificationFailure(f"{context} must bind all three platforms")
    if (
        len(
            _array(provenance["rich_groups"], "characterization_provenance.rich_groups")
        )
        != 11
        or seeded != 7
        or rich_only != 4
    ):
        raise CertificationFailure(
            "rich diagnostic evidence must remain exactly 11/7/4"
        )
    if rich_case_counts != Counter(
        {
            ("blastn", "PesePMNV.MjPMNV.task_blastn"): 3,
            ("blastn", "Sakai.MG1655.megablast"): 7,
            ("blastp", "pairwise_default_serial"): 1,
        }
    ):
        raise CertificationFailure("rich diagnostic representative counts changed")

    return NativeAuthority(
        path=Path(),
        authority_version=NATIVE_AUTHORITY_VERSION,
        file_sha256="",
        document=document,
        representatives=representatives,
        platforms=platforms,
        outputs=outputs,
        diagnostics=diagnostics,
    )


# NCBI references:
# - ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-832
#   dispatches the selected output format.
# - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
#   writes the selected fields and row terminator directly to the output stream.
def load_native_authority(path: Path) -> NativeAuthority:
    if not path.is_file():
        raise CertificationFailure(f"native authority file is missing: {path}")
    raw = path.read_bytes()
    try:
        document = json.loads(raw, object_pairs_hook=_strict_json_object)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise CertificationFailure(f"invalid native authority JSON: {error}") from error
    validated = _validate_native_authority_document(_object(document, "authority"))
    return NativeAuthority(
        path=path.resolve(),
        authority_version=validated.authority_version,
        file_sha256=hashlib.sha256(raw).hexdigest(),
        document=validated.document,
        representatives=validated.representatives,
        platforms=validated.platforms,
        outputs=validated.outputs,
        diagnostics=validated.diagnostics,
    )


def resolve_native_authority(
    authority: NativeAuthority, platform_id: str, program: str, case_id: str
) -> tuple[dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    if authority.authority_version != NATIVE_AUTHORITY_VERSION:
        raise native_authority_miss("unknown native authority version")
    platform_record = authority.platforms.get(platform_id)
    representative = authority.representatives.get((program, case_id))
    output = authority.outputs.get((platform_id, program, case_id))
    diagnostic = authority.diagnostics.get((platform_id, program, case_id))
    if any(
        item is None for item in (platform_record, representative, output, diagnostic)
    ):
        raise CertificationFailure(
            "NATIVE_NCBI_AUTHORITY_MISS: no exact authority tuple for "
            f"{authority.authority_version}/{platform_id}/{program}/{case_id}"
        )
    return platform_record, representative, output, diagnostic


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


def run_capture(
    command: Sequence[str], cwd: Path
) -> subprocess.CompletedProcess[bytes]:
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
        return {0x01000007: "x86_64", 0x0100000C: "arm64"}.get(cpu_type, "unknown")
    if len(data) >= 8 and data[:4] in {b"\xfe\xed\xfa\xcf", b"\xfe\xed\xfa\xce"}:
        cpu_type = int.from_bytes(data[4:8], "big")
        return {0x01000007: "x86_64", 0x0100000C: "arm64"}.get(cpu_type, "unknown")
    raise CertificationFailure(f"unsupported native executable format: {path}")


def _verify_archive_identity(
    platform_spec: PlatformSpec,
    archive: Path,
    checksum: Path,
    authority_platform: dict[str, object] | None = None,
) -> dict[str, str]:
    authority_archive = (
        _object(authority_platform["archive"], "authority platform archive")
        if authority_platform is not None
        else None
    )
    expected_name = (
        str(authority_archive["filename"])
        if authority_archive is not None
        else platform_spec.archive_name
    )
    expected_md5 = (
        str(_object(authority_archive["checksum"], "archive checksum")["value"])
        if authority_archive is not None
        else platform_spec.archive_md5
    )
    expected_url = (
        str(authority_archive["retrieval_url"])
        if authority_archive is not None
        else platform_spec.archive_url
    )
    if archive.name != expected_name:
        raise CertificationFailure(
            f"wrong NCBI archive for {platform_spec.platform_id}: {archive.name}"
        )
    if not archive.is_file() or not checksum.is_file():
        raise CertificationFailure("NCBI archive and published checksum are required")
    checksum_parts = checksum.read_text(encoding="ascii").strip().split()
    if checksum_parts != [expected_md5, expected_name]:
        raise CertificationFailure(
            f"published checksum record is not the pinned {platform_spec.platform_id} record"
        )
    observed = md5_path(archive)
    if observed != expected_md5:
        raise CertificationFailure(
            f"NCBI archive checksum mismatch: expected {expected_md5}, "
            f"observed {observed}"
        )
    return {
        "archive": str(archive.resolve()),
        "archive_name": expected_name,
        "archive_md5": observed,
        "archive_checksum_algorithm": "md5",
        "archive_url": expected_url,
        "checksum": str(checksum.resolve()),
        "checksum_url": f"{expected_url}.md5",
    }


def verify_archive(
    platform_spec: PlatformSpec,
    archive: Path,
    checksum: Path,
    authority_platform: dict[str, object] | None = None,
) -> dict[str, str]:
    try:
        return _verify_archive_identity(
            platform_spec, archive, checksum, authority_platform
        )
    except CertificationFailure as error:
        if "NATIVE_NCBI_AUTHORITY_MISS" in str(error):
            raise
        raise native_authority_miss(str(error)) from error
    except (KeyError, OSError, TypeError, UnicodeError, ValueError) as error:
        raise native_authority_miss(f"invalid archive identity: {error}") from error


def extract_verified_archive(
    platform_spec: PlatformSpec,
    archive: Path,
    checksum: Path,
    destination: Path,
    authority_platform: dict[str, object] | None = None,
) -> dict[str, object]:
    archive_record = verify_archive(
        platform_spec, archive, checksum, authority_platform
    )
    destination.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(archive, mode="r:gz") as bundle:
            members = bundle.getmembers()
            if sys.version_info >= (3, 12):
                bundle.extractall(path=destination, filter="fully_trusted")
            else:
                bundle.extractall(path=destination)
    except tarfile.TarError as error:
        raise native_authority_miss(
            f"verified NCBI archive extraction failed: {error}"
        ) from error
    return {
        "status": "EXTRACTED_VERIFIED_ARCHIVE",
        "archive": archive_record,
        "destination": str(destination.absolute()),
        "member_count": len(members),
        "implementation": "python.tarfile.open(mode=r:gz).extractall",
    }


def _controlled_relative(value: str | PurePosixPath) -> PurePosixPath:
    relative = PurePosixPath(value)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise CertificationFailure(f"invalid controlled fixture path: {value}")
    return relative


def historical_lexical_path(relative: str | PurePosixPath) -> str:
    controlled = _controlled_relative(relative)
    return f"{HISTORICAL_LEXICAL_ROOT}/{controlled.as_posix()}"


def controlled_physical_path(lexical_path: str) -> Path:
    prefix = f"{HISTORICAL_LEXICAL_ROOT}/"
    if lexical_path != HISTORICAL_LEXICAL_ROOT and not lexical_path.startswith(prefix):
        raise CertificationFailure(
            f"path is outside the controlled lexical root: {lexical_path}"
        )
    return Path(os.path.abspath(lexical_path))


def preflight_controlled_fixture_path(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "controlled-path-preflight.json"
    lexical_sentinel = historical_lexical_path(
        ".platform-native-certification-path-preflight"
    )
    physical_sentinel = controlled_physical_path(lexical_sentinel)
    payload = b"LOSAT_PLATFORM_CERTIFICATION_PATH_PREFLIGHT\n"
    record: dict[str, object] = {
        "status": "IN_PROGRESS",
        "lexical_root": HISTORICAL_LEXICAL_ROOT,
        "lexical_sentinel_path": lexical_sentinel,
        "physical_sentinel_path": str(physical_sentinel),
        "expected_sha256": hashlib.sha256(payload).hexdigest(),
        "platform": platform.system(),
        "subprocess_contract": {
            "api": "subprocess.run",
            "shell": False,
            "msys_or_git_bash_involved": False,
            "search_driver_uses_literal_argv": True,
        },
    }
    try:
        if platform.system() == "Windows" and not physical_sentinel.drive:
            raise CertificationFailure(
                "Windows controlled physical path is not drive-qualified"
            )
        physical_sentinel.parent.mkdir(parents=True, exist_ok=True)
        physical_sentinel.write_bytes(payload)
        record["physical_path_created"] = physical_sentinel.is_file()
        with open(lexical_sentinel, "rb") as handle:
            direct_bytes = handle.read()
        if direct_bytes != payload:
            raise CertificationFailure("controlled lexical sentinel read changed bytes")
        probe = (
            "import hashlib,json,sys;"
            "p=sys.argv[1];"
            "data=open(p,'rb').read();"
            "print(json.dumps({'argv':p,'sha256':hashlib.sha256(data).hexdigest()}))"
        )
        completed = subprocess.run(
            [sys.executable, "-c", probe, lexical_sentinel],
            capture_output=True,
            text=False,
            check=False,
            shell=False,
        )
        record["subprocess"] = {
            "argv": [sys.executable, "-c", "<probe>", lexical_sentinel],
            "returncode": completed.returncode,
            "stdout": completed.stdout.decode("utf-8", errors="replace"),
            "stderr": completed.stderr.decode("utf-8", errors="replace"),
        }
        if completed.returncode != 0:
            raise CertificationFailure(
                f"controlled path subprocess failed: {completed.returncode}"
            )
        child_record = json.loads(completed.stdout.decode("utf-8"))
        if child_record != {
            "argv": lexical_sentinel,
            "sha256": hashlib.sha256(payload).hexdigest(),
        }:
            raise CertificationFailure(
                "shell-free subprocess did not preserve the controlled lexical argv"
            )
        record["direct_read_sha256"] = hashlib.sha256(direct_bytes).hexdigest()
        record["lexical_read_verified"] = True
        record["child_observation"] = child_record
        record["status"] = "VERIFIED"
        atomic_write_json(evidence_path, record)
        return record
    except (CertificationFailure, OSError, ValueError) as error:
        record["status"] = "FAILED"
        record["error"] = str(error)
        atomic_write_json(evidence_path, record)
        if isinstance(error, CertificationFailure):
            raise
        raise CertificationFailure(
            f"controlled path preflight failed: {error}"
        ) from error
    finally:
        physical_sentinel.unlink(missing_ok=True)


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
        if (
            len(authority_ids[program]) != expected_count
            or canonical_ids != authority_ids[program]
        ):
            raise CertificationFailure(
                f"{program} canonical cases differ from the existing authority"
            )

    blastn_classes = Counter(
        row["classification"] for row in canonical_rows if row["program"] == "blastn"
    )
    if blastn_classes != Counter({"EXACT_TEXT": 13, "SOURCE_UNDETERMINED_ACCEPTED": 1}):
        raise CertificationFailure("BLASTN canonical classifications changed")
    sakai = canonical[("blastn", "Sakai.MG1655.megablast")]
    if sakai["contract"] != "SOURCE_UNDETERMINED_ACCEPTED":
        raise CertificationFailure("Sakai must not be represented as NCBI byte-exact")

    blastp_classes = Counter(
        row["classification"] for row in canonical_rows if row["program"] == "blastp"
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
        if (
            row["contract"] != case.contract
            or row["classification"] != expected_classification
        ):
            raise CertificationFailure(
                f"TBLASTX canonical contract differs from authority: {case_id}"
            )
    if Counter(
        row["classification"] for row in canonical_rows if row["program"] == "tblastx"
    ) != Counter({"EXACT_TEXT": 14, "HSP_SET_DIFF": 6}):
        raise CertificationFailure("TBLASTX canonical classifications changed")

    representative_keys = {
        (program, case_id) for program, case_id, _ in REPRESENTATIVES
    }
    if len(representative_keys) != 6 or not representative_keys <= set(canonical):
        raise CertificationFailure("the exact six PR 6 representatives are unavailable")
    if ("tblastx", "p11_avclpv_psclpv") in representative_keys:
        raise CertificationFailure(
            "PR 5-only p11 repeatability selection leaked into PR 6"
        )

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
    command_row["query"] = str(repo_root / "LOSAT" / row["query"])
    command_row["subject"] = str(repo_root / "LOSAT" / row["subject"])
    return command_row


def _input_arguments(
    command: Sequence[str], program: str
) -> list[tuple[str, str, str]]:
    arguments: list[tuple[str, str, str]] = []
    for role, aliases in INPUT_OPTIONS[program].items():
        present = [option for option in aliases if option in command]
        if len(present) != 1 or command.count(present[0]) != 1:
            raise CertificationFailure(
                f"{program} command must contain exactly one {role} input: {command}"
            )
        option = present[0]
        option_index = command.index(option)
        if option_index + 1 >= len(command):
            raise CertificationFailure(f"{program} command has no value for {option}")
        arguments.append((role, option, command[option_index + 1]))
    return arguments


# NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:794-808
# ```c++
# if (m_FormatType == CFormattingArgs::eTabularWithComments) {
#     if (m_IsDbScan)
#         dbname = string("User specified sequence set (Input: ")
#                  + m_SubjectTag + string(")");
#     tabinfo.PrintHeader(strProgVersion, *(bhandle.GetBioseqCore()),
#                         dbname, results.GetRID(), itr_num, aln_set,
#                         subject_bioseq);
# }
# ```
# The input argument is therefore certification provenance for commented tabular
# output. Only its filesystem layer is redirected to the immutable PR 5 spelling.
def _adapt_fixture_arguments(
    command: Sequence[str], program: str, repo_root: Path
) -> list[str]:
    updated = list(command)
    source_root = repo_root / "LOSAT"
    for _, option, source_value in _input_arguments(updated, program):
        source_path = Path(source_value)
        try:
            relative = source_path.relative_to(source_root)
        except ValueError as error:
            raise CertificationFailure(
                f"authoritative {program} fixture is outside LOSAT: {source_path}"
            ) from error
        updated[updated.index(option) + 1] = historical_lexical_path(
            PurePosixPath(*relative.parts)
        )
    return updated


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
        commands = (
            catalog.blastn_compare.build_ncbi_command(
                row, str(oracles[program]), unused_ncbi
            ),
            catalog.blastn_compare.build_losat_command(row, native_bin, unused_losat),
        )
    elif program == "blastp":
        commands = catalog.blastp_audit.build_commands(
            catalog.blastp_cases[case_id],
            oracles[program],
            native_bin,
            output_dir,
        )
    elif program == "tblastx":
        commands = catalog.tblastx_audit.build_commands(
            catalog.tblastx_cases[case_id],
            oracles[program],
            native_bin,
            unused_ncbi,
            unused_losat,
        )
    else:
        raise CertificationFailure(f"unsupported program: {program}")
    ncbi_command, losat_command = commands
    return (
        _adapt_fixture_arguments(ncbi_command, program, repo_root),
        _adapt_fixture_arguments(losat_command, program, repo_root),
    )


def _option_value(command: Sequence[str], option: str) -> str:
    if command.count(option) != 1:
        raise CertificationFailure(
            f"command must contain exactly one {option}: {command}"
        )
    option_index = command.index(option)
    if option_index + 1 >= len(command):
        raise CertificationFailure(f"command has no value for {option}")
    return command[option_index + 1]


def validate_planned_input_contract(steps: Sequence[SearchStep]) -> None:
    matrix_steps = [step for step in steps if step.kind == "matrix"]
    if len(matrix_steps) != 43:
        raise CertificationFailure(
            f"controlled path validation expected 43 matrix cases, found {len(matrix_steps)}"
        )
    prefix = f"{HISTORICAL_LEXICAL_ROOT}/"
    for step in steps:
        for _, option, value in _input_arguments(step.command, step.program):
            if not value.startswith(prefix):
                raise CertificationFailure(
                    f"planned {step.step_id} {option} is outside controlled root: {value}"
                )

    blastn_matrix = [step for step in matrix_steps if step.program == "blastn"]
    path_sensitive = [
        step for step in blastn_matrix if _option_value(step.command, "--outfmt") == "7"
    ]
    headerless = [
        step for step in blastn_matrix if _option_value(step.command, "--outfmt") == "6"
    ]
    if len(path_sensitive) != 13 or len(headerless) != 1:
        raise CertificationFailure(
            "BLASTN path contract must remain 13 outfmt 7 plus one outfmt 6 case"
        )
    for step in path_sensitive:
        subject = {
            role: value
            for role, _, value in _input_arguments(step.command, step.program)
        }["subject"]
        if not subject.startswith(prefix):
            raise CertificationFailure(
                f"BLASTN outfmt 7 subject is outside controlled root: {step.case_id}"
            )


def _fixture_relative_from_lexical(lexical_path: str) -> PurePosixPath:
    prefix = f"{HISTORICAL_LEXICAL_ROOT}/"
    if not lexical_path.startswith(prefix):
        raise CertificationFailure(
            f"fixture argument is outside controlled root: {lexical_path}"
        )
    return _controlled_relative(lexical_path.removeprefix(prefix))


def required_fixture_relatives(steps: Sequence[SearchStep]) -> list[PurePosixPath]:
    required: set[PurePosixPath] = set()
    for step in steps:
        if step.kind != "matrix":
            continue
        for _, _, value in _input_arguments(step.command, step.program):
            required.add(_fixture_relative_from_lexical(value))
    return sorted(required, key=PurePosixPath.as_posix)


def stage_required_fixtures(
    repo_root: Path, steps: Sequence[SearchStep], output_dir: Path
) -> list[StagedFixture]:
    evidence_path = output_dir / "fixture-staging.json"
    records: list[StagedFixture] = []
    status = "IN_PROGRESS"
    failure: str | None = None
    try:
        for relative in required_fixture_relatives(steps):
            source_repo_relative = PurePosixPath("LOSAT") / relative
            source = repo_root.joinpath(*source_repo_relative.parts)
            lexical_target = historical_lexical_path(relative)
            physical_target = controlled_physical_path(lexical_target)
            if not source.is_file():
                raise CertificationFailure(
                    f"authoritative fixture is missing: {source_repo_relative.as_posix()}"
                )
            source_hash = sha256_path(source)
            physical_target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, physical_target)
            target_hash = sha256_path(physical_target)
            if (
                target_hash != source_hash
                or physical_target.read_bytes() != source.read_bytes()
            ):
                raise CertificationFailure(
                    f"staged fixture bytes changed: {source_repo_relative.as_posix()}"
                )
            records.append(
                StagedFixture(
                    source_repo_relative=source_repo_relative.as_posix(),
                    source_path=str(source.absolute()),
                    lexical_target_path=lexical_target,
                    physical_target_path=str(physical_target),
                    sha256=source_hash,
                )
            )
        status = "VERIFIED"
        return records
    except (CertificationFailure, OSError) as error:
        failure = str(error)
        if isinstance(error, CertificationFailure):
            raise
        raise CertificationFailure(f"fixture staging failed: {error}") from error
    finally:
        evidence: dict[str, object] = {
            "status": status if failure is None else "FAILED",
            "lexical_root": HISTORICAL_LEXICAL_ROOT,
            "staged_input_count": len(records),
            "inputs": [asdict(record) for record in records],
        }
        if failure is not None:
            evidence["error"] = failure
        atomic_write_json(evidence_path, evidence)


def fixture_staging_identity(records: Sequence[StagedFixture]) -> dict[str, object]:
    return {
        "lexical_root": HISTORICAL_LEXICAL_ROOT,
        "staged_input_count": len(records),
        "inputs": [
            {
                "source_repo_relative": record.source_repo_relative,
                "lexical_target_path": record.lexical_target_path,
                "sha256": record.sha256,
            }
            for record in records
        ],
    }


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
        environment = (
            (("RAYON_NUM_THREADS", "1"),)
            if kind in {"matrix", "repeatability"} and program == "tblastx"
            else ()
        )
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
            rel = f"repeatability/{program}/{case_id}/run{run_index}.out.partial"
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
    validate_planned_input_contract(steps)
    return steps


def canonical_newlines(data: bytes) -> bytes:
    return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def newline_profile(data: bytes) -> dict[str, int]:
    crlf = data.count(b"\r\n")
    without_crlf = data.replace(b"\r\n", b"")
    return {
        "lf_count": data.count(b"\n"),
        "crlf_count": crlf,
        "bare_lf_count": data.count(b"\n") - crlf,
        "bare_cr_count": without_crlf.count(b"\r"),
    }


def output_data_row_count(data: bytes) -> int:
    return sum(
        bool(line.strip()) and not line.lstrip().startswith(b"#")
        for line in data.splitlines()
    )


def validate_native_platform_sanity(
    authority: NativeAuthority,
    platform_id: str,
    *,
    observed_system: str | None = None,
    observed_machine: str | None = None,
) -> dict[str, object]:
    platform_record = authority.platforms.get(platform_id)
    if platform_record is None:
        raise native_authority_miss(f"unknown native authority platform: {platform_id}")
    system = platform.system() if observed_system is None else observed_system
    raw_machine = platform.machine() if observed_machine is None else observed_machine
    machine = normalize_machine(raw_machine)
    if (
        system != platform_record["platform_system"]
        or machine != platform_record["normalized_machine"]
    ):
        raise native_authority_miss(
            "native authority platform mismatch: expected "
            f"{platform_record['platform_system']}/{platform_record['normalized_machine']}, "
            f"observed {system}/{raw_machine} ({machine})"
        )
    return platform_record


def validate_native_authority_catalog(
    authority: NativeAuthority, catalog: Catalog
) -> None:
    for key, representative in authority.representatives.items():
        canonical = catalog.canonical.get(key)
        if canonical is None:
            raise CertificationFailure(
                f"native representative is absent from Gate A: {key}"
            )
        if (
            representative["losat_contract"] != canonical["contract"]
            or representative["losat_parity_class"] != canonical["classification"]
        ):
            raise CertificationFailure(
                f"native authority changed the independent LOSAT parity axis: {key}"
            )
    sakai = authority.representatives[("blastn", "Sakai.MG1655.megablast")]
    if (
        sakai["losat_contract"] != "SOURCE_UNDETERMINED_ACCEPTED"
        or sakai["losat_parity_class"] != "SOURCE_UNDETERMINED_ACCEPTED"
    ):
        raise CertificationFailure("Sakai LOSAT parity authority changed")
    d06 = authority.representatives[("tblastx", "d06_ap027131_ap027133_db4")]
    if (
        d06["losat_contract"] != "approved_db_gencode_deviation"
        or d06["losat_parity_class"] != "HSP_SET_DIFF"
    ):
        raise CertificationFailure("d06 LOSAT parity authority changed")
    accepted_tblastx = sum(
        row["contract"] == "approved_db_gencode_deviation"
        for row in catalog.canonical_rows
        if row["program"] == "tblastx"
    )
    if accepted_tblastx != 6:
        raise CertificationFailure(
            f"TBLASTX accepted-deviation ceiling changed: {accepted_tblastx}"
        )


def _authority_executable(
    platform_record: dict[str, object], program: str
) -> dict[str, object]:
    matches = [
        _object(value, "native executable")
        for value in _array(platform_record["executables"], "native executables")
        if _object(value, "native executable").get("program") == program
    ]
    if len(matches) != 1:
        raise native_authority_miss(f"native executable authority miss: {program}")
    return matches[0]


# NCBI references:
# - ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-832
#   selects the requested tabular formatter before the search output is written.
# - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
#   preserves the requested field order and row terminator in the raw stream.
def _verify_native_invocation_identity(
    repo_root: Path,
    output_dir: Path,
    authority: NativeAuthority,
    identity: dict[str, object],
    step: SearchStep,
) -> dict[str, object]:
    platform_id = str(identity["platform_id"])
    platform_record, representative, _, _ = resolve_native_authority(
        authority, platform_id, step.program, step.case_id
    )
    if step.kind != "oracle":
        raise CertificationFailure(
            "Gate B invocation identity used outside an oracle step"
        )
    if (
        identity.get("authority_version") != authority.authority_version
        or identity.get("authority_file_sha256") != authority.file_sha256
    ):
        raise CertificationFailure("native invocation identity is not authority-bound")
    archive = _object(platform_record["archive"], "native archive authority")
    checksum = _object(archive["checksum"], "native archive checksum")
    recorded_archive = _object(identity["oracle_archive"], "recorded oracle archive")
    if (
        recorded_archive.get("archive_name") != archive["filename"]
        or recorded_archive.get("archive_checksum_algorithm") != checksum["algorithm"]
        or recorded_archive.get("archive_md5") != checksum["value"]
    ):
        raise CertificationFailure("native invocation archive identity changed")
    executable = _authority_executable(platform_record, step.program)
    recorded_oracle = _object(
        _object(identity["oracles"], "recorded oracles")[step.program],
        f"recorded {step.program} oracle",
    )
    if (
        recorded_oracle.get("sha256") != executable["sha256"]
        or recorded_oracle.get("architecture") != executable["architecture"]
        or step.command[0] != recorded_oracle.get("path")
    ):
        raise CertificationFailure("native invocation executable identity changed")
    if step.environment:
        raise CertificationFailure(
            "native official NCBI invocation received LOSAT environment overrides"
        )
    command = list(step.command)
    if command.count("-out") != 1:
        raise CertificationFailure("native invocation must have one output option")
    output_index = command.index("-out") + 1
    expected_output = (output_dir / step.output_rel).absolute()
    if Path(command[output_index]).absolute() != expected_output:
        raise CertificationFailure("native invocation output path identity changed")
    command[0] = "{oracle}"
    command[output_index] = "{output}"
    if command != representative["ordered_argv_template"]:
        raise CertificationFailure(
            f"native invocation argv identity changed: {step.program}/{step.case_id}"
        )
    if (
        representative["shell"] is not False
        or representative["environment_overrides"] != {}
        or representative["working_directory_contract"] != "repository_root"
        or representative["lexical_fixture_root"] != HISTORICAL_LEXICAL_ROOT
    ):
        raise CertificationFailure("native invocation process contract changed")
    input_identities = {}
    for role in ("query", "subject"):
        expected_input = _object(representative[role], f"native {role} authority")
        option = str(expected_input["option"])
        if command.count(option) != 1:
            raise CertificationFailure(f"native invocation {role} option changed")
        lexical = str(expected_input["lexical_path"])
        if command[command.index(option) + 1] != lexical:
            raise CertificationFailure(f"native invocation {role} lexical path changed")
        repository_path = repo_root / str(expected_input["repository_relative"])
        controlled_path = controlled_physical_path(lexical)
        for label, path in (
            ("repository", repository_path),
            ("controlled", controlled_path),
        ):
            if not path.is_file():
                raise CertificationFailure(f"native {role} {label} input is missing")
            data = path.read_bytes()
            observed = {
                "sha256": hashlib.sha256(data).hexdigest(),
                "byte_length": len(data),
                "newline_profile": newline_profile(data),
            }
            if observed != {
                "sha256": expected_input["sha256"],
                "byte_length": expected_input["byte_length"],
                "newline_profile": expected_input["newline_profile"],
            }:
                raise CertificationFailure(
                    f"native invocation {role} input identity changed"
                )
        if repository_path.read_bytes() != controlled_path.read_bytes():
            raise CertificationFailure(
                f"native invocation {role} staging changed bytes"
            )
        input_identities[role] = {
            "role": role,
            "repository_relative": expected_input["repository_relative"],
            "lexical_path": lexical,
            "sha256": expected_input["sha256"],
            "byte_length": expected_input["byte_length"],
            "newline_profile": expected_input["newline_profile"],
        }
    return {
        "authority_version": authority.authority_version,
        "authority_file_sha256": authority.file_sha256,
        "selector": {
            "authority_version": authority.authority_version,
            "platform_id": platform_id,
            "program": step.program,
            "case_id": step.case_id,
        },
        "ncbi_release": archive["release"],
        "archive_name": archive["filename"],
        "archive_checksum": checksum,
        "executable_sha256": executable["sha256"],
        "architecture": executable["architecture"],
        "inputs": input_identities,
        "ordered_argv_template_sha256": hashlib.sha256(
            json.dumps(command, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "authoritative_outfmt": representative["authoritative_outfmt"],
        "working_directory_contract": "repository_root",
        "lexical_fixture_root": HISTORICAL_LEXICAL_ROOT,
        "shell": False,
        "environment_overrides": {},
    }


def verify_native_invocation_identity(
    repo_root: Path,
    output_dir: Path,
    authority: NativeAuthority,
    identity: dict[str, object],
    step: SearchStep,
) -> dict[str, object]:
    try:
        return _verify_native_invocation_identity(
            repo_root, output_dir, authority, identity, step
        )
    except CertificationFailure as error:
        if "NATIVE_NCBI_AUTHORITY_MISS" in str(error):
            raise
        raise native_authority_miss(str(error)) from error
    except (IndexError, KeyError, OSError, TypeError, ValueError) as error:
        raise native_authority_miss(
            f"invalid native invocation identity: {error}"
        ) from error


# NCBI references:
# - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:59-91
#   expands the selected tabular fields in their requested order.
# - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
#   writes those exact rows and terminators; Gate B therefore hashes raw bytes.
def verify_native_raw_fingerprint(
    authority: NativeAuthority,
    platform_id: str,
    program: str,
    case_id: str,
    output: Path,
) -> dict[str, object]:
    _, representative, expected, diagnostic = resolve_native_authority(
        authority, platform_id, program, case_id
    )
    try:
        data = output.read_bytes()
    except OSError as error:
        raise native_authority_miss(
            f"cannot read native output fingerprint: {error}"
        ) from error
    observed = {
        "raw_sha256": hashlib.sha256(data).hexdigest(),
        "byte_length": len(data),
        "newline_profile": newline_profile(data),
        "data_row_count": output_data_row_count(data),
    }
    expected_fingerprint = {
        "raw_sha256": expected["raw_sha256"],
        "byte_length": expected["byte_length"],
        "newline_profile": expected["newline_profile"],
        "data_row_count": expected["data_row_count"],
    }
    if observed != expected_fingerprint:
        mismatches = [
            key
            for key in expected_fingerprint
            if observed[key] != expected_fingerprint[key]
        ]
        raise CertificationFailure(
            "NATIVE_NCBI_AUTHORITY_MISS: exact raw fingerprint mismatch for "
            f"{platform_id}/{program}/{case_id}: {mismatches}"
        )
    return {
        "classification": "NATIVE_NCBI_FINGERPRINT_EXACT",
        "authority_version": authority.authority_version,
        "authority_file_sha256": authority.file_sha256,
        "losat_contract": representative["losat_contract"],
        "losat_parity_class": representative["losat_parity_class"],
        "native_ncbi_reference_class": diagnostic["native_ncbi_reference_class"],
        **observed,
    }


def line_endings_only(left: Path, right: Path) -> bool:
    left_bytes = left.read_bytes()
    right_bytes = right.read_bytes()
    return left_bytes != right_bytes and canonical_newlines(
        left_bytes
    ) == canonical_newlines(right_bytes)


# NCBI references:
# - ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:828-832
#   emits retained alignments before this post-acceptance comparison exists.
# - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
#   writes exact field values and terminators already bound by Gate B's raw hash.
def report_native_vs_losat_diagnostic(
    catalog: Catalog,
    step: SearchStep,
    ncbi_output: Path,
    losat_output: Path,
    diagnostic_dir: Path,
) -> dict[str, object]:
    result: dict[str, object] = {
        "runtime_enforced": False,
        "program": step.program,
        "case_id": step.case_id,
        "platform_diagnostic": "none",
    }
    try:
        diagnostic_dir.mkdir(parents=True, exist_ok=True)
        result.update(
            {
                "ncbi_sha256": sha256_path(ncbi_output),
                "losat_sha256": sha256_path(losat_output),
                "ncbi_newlines": newline_profile(ncbi_output.read_bytes()),
                "losat_newlines": newline_profile(losat_output.read_bytes()),
            }
        )
        if step.program == "blastn":
            with contextlib.redirect_stdout(io.StringIO()) as captured:
                structured = catalog.blastn_compare.compare_case(
                    step.case_id, ncbi_output, losat_output, 5
                )
                raw_equal = catalog.blastn_compare.compare_raw_bytes(
                    ncbi_output, losat_output
                )
            atomic_write_bytes(
                diagnostic_dir / "comparison.txt",
                captured.getvalue().encode("utf-8"),
            )
            ncbi_keys = set(structured["ncbi_by_key"])
            losat_keys = set(structured["losat_by_key"])
            if raw_equal:
                classification = "EXACT_TEXT"
            elif structured["structured_equal"]:
                classification = "EXACT_DATA"
            elif structured["row_multiset_equal"]:
                classification = "ORDER_ONLY"
            elif ncbi_keys == losat_keys:
                classification = "VALUE_DIFF"
            else:
                classification = "HSP_SET_DIFF"
            if line_endings_only(ncbi_output, losat_output):
                result["platform_diagnostic"] = "LINE_ENDINGS_ONLY"
            result["classification"] = classification
            result["metrics"] = {
                "ncbi_rows": len(structured["ncbi_rows"]),
                "losat_rows": len(structured["losat_rows"]),
                "row_count_equal": structured["row_count_equal"],
                "row_multiset_equal": structured["row_multiset_equal"],
                "row_order_equal": structured["row_order_equal"],
                "structured_equal": structured["structured_equal"],
                "ncbi_only_coord_keys": len(structured["ncbi_only"]),
                "losat_only_coord_keys": len(structured["losat_only"]),
                "bitscore_differences": len(structured["bitscore_diffs"]),
                "evalue_differences": len(structured["evalue_diffs"]),
                "pident_differences": len(structured["pident_diffs"]),
                "other_value_differences": len(structured["other_value_diffs"]),
            }
            return result

        if step.program == "blastp":
            classification, metrics = catalog.blastp_audit.classify_outputs(
                ncbi_output, losat_output
            )
            result["classification"] = classification
            result["metrics"] = metrics
            if line_endings_only(ncbi_output, losat_output):
                result["platform_diagnostic"] = "LINE_ENDINGS_ONLY"
            return result

        classification, metrics = catalog.tblastx_audit.classify_outputs(
            ncbi_output, losat_output
        )
        result["classification"] = classification
        result["metrics"] = metrics
        if line_endings_only(ncbi_output, losat_output):
            result["platform_diagnostic"] = "LINE_ENDINGS_ONLY"
        return result
    except Exception as error:  # Diagnostic reporting cannot veto exact Gate B.
        result["classification"] = "DIAGNOSTIC_ERROR"
        result["diagnostic_error"] = f"{type(error).__name__}: {error}"
        return result


def validate_git_identity(repo_root: Path, expected_sha: str) -> str:
    if not SHA_PATTERN.fullmatch(expected_sha):
        raise CertificationFailure(
            "expected SHA must be exactly 40 lowercase hex digits"
        )
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
    return (
        (completed.stdout + completed.stderr).decode("utf-8", errors="replace").strip()
    )


def record_identity(
    repo_root: Path,
    output_dir: Path,
    expected_sha: str,
    platform_spec: PlatformSpec,
    native_bin: Path,
    archive_record: dict[str, str],
    oracles: dict[str, Path],
    canonical_manifest: Path,
    controlled_fixtures: dict[str, object],
    authority: NativeAuthority,
) -> dict[str, object]:
    observed_os = platform.system()
    observed_machine = normalize_machine(platform.machine())
    authority_platform = validate_native_platform_sanity(
        authority,
        platform_spec.platform_id,
        observed_system=observed_os,
        observed_machine=platform.machine(),
    )
    rustc = _tool_output(["rustc", "-vV"], repo_root, "rustc")
    cargo = _tool_output(["cargo", "-V"], repo_root, "cargo")
    host = next(
        (
            line.removeprefix("host: ")
            for line in rustc.splitlines()
            if line.startswith("host: ")
        ),
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
            raise native_authority_miss(
                f"official {program} oracle is missing: {oracle}"
            )
        try:
            version = _tool_output([str(oracle), "-version"], repo_root, program)
        except CertificationFailure as error:
            raise native_authority_miss(str(error)) from error
        if (
            f"{program}: 2.17.0+" not in version
            or "Package: blast 2.17.0" not in version
        ):
            raise native_authority_miss(f"{program} oracle version mismatch: {version}")
        expected_oracle = _authority_executable(authority_platform, program)
        try:
            observed_oracle_sha256 = sha256_path(oracle)
            observed_oracle_architecture = binary_architecture(oracle)
        except (CertificationFailure, OSError) as error:
            raise native_authority_miss(
                f"official {program} executable identity is unreadable: {error}"
            ) from error
        if (
            observed_oracle_sha256 != expected_oracle["sha256"]
            or observed_oracle_architecture != expected_oracle["architecture"]
        ):
            raise native_authority_miss(
                f"official {program} executable identity differs from native authority"
            )
        oracle_identities[program] = {
            "path": str(oracle.resolve()),
            "sha256": observed_oracle_sha256,
            "architecture": observed_oracle_architecture,
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
        "authority_version": authority.authority_version,
        "authority_file_sha256": authority.file_sha256,
        "authority_file": str(authority.path),
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
        "controlled_fixtures": controlled_fixtures,
        "execution_contract": {
            "matrix": 43,
            "oracle": 6,
            "repeatability": 12,
            "per_platform": 61,
            "gate_a_losat_canonical": 43,
            "gate_b_native_ncbi_reference": 6,
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
    if "controlled_fixtures" not in identity:
        raise CertificationFailure(
            "resume identity lacks the controlled fixture-path contract"
        )
    return {
        "evidence_schema": identity["evidence_schema"],
        "authority_version": identity["authority_version"],
        "authority_file_sha256": identity["authority_file_sha256"],
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
        "controlled_fixtures": identity["controlled_fixtures"],
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
        "authority_version": identity["authority_version"],
        "authority_file_sha256": identity["authority_file_sha256"],
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
    repo_root: Path,
    identity: dict[str, object],
    steps: Sequence[SearchStep],
    catalog: Catalog,
    authority: NativeAuthority,
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
        if (
            not prior_output.is_file()
            or sha256_path(prior_output) != record["output_sha256"]
        ):
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
            invocation = verify_native_invocation_identity(
                repo_root,
                output_dir,
                authority,
                identity,
                step,
            )
            fingerprint = verify_native_raw_fingerprint(
                authority,
                str(identity["platform_id"]),
                step.program,
                step.case_id,
                current_output,
            )
            diagnostic = report_native_vs_losat_diagnostic(
                catalog,
                step,
                current_output,
                losat_output,
                output_dir / f"oracle/{step.program}/{step.case_id}/resume-diagnostic",
            )
            verification = {
                "gate": "PLATFORM_NATIVE_NCBI_REFERENCE",
                "invocation_identity": invocation,
                "native_fingerprint": fingerprint,
                "native_vs_losat_diagnostic": diagnostic,
            }
            prior_verification = _object(
                record.get("verification"), "resume oracle verification"
            )
            if any(
                prior_verification.get(key) != verification[key]
                for key in ("gate", "invocation_identity", "native_fingerprint")
            ):
                raise CertificationFailure(
                    f"resume oracle normative verification changed: {step.step_id}"
                )
            record["verification"] = verification
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
    authority: NativeAuthority | None = None,
) -> dict[str, object]:
    partial_output = output_dir / step.output_rel
    final_output = output_dir / step.final_output_rel
    partial_output.parent.mkdir(parents=True, exist_ok=True)
    if partial_output.exists():
        partial_output.unlink()
    invocation_identity = None
    if step.kind == "oracle":
        if authority is None:
            raise CertificationFailure("Gate B requires the bounded native authority")
        invocation_identity = verify_native_invocation_identity(
            repo_root, output_dir, authority, identity, step
        )
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
        shell=False,
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
            "gate": (
                "LOSAT_CANONICAL" if step.kind == "matrix" else "LOSAT_REPEATABILITY"
            ),
            "classification": "CANONICAL_PR5_RAW_BYTES",
            "canonical_sha256": step.expected_losat_sha256,
        }
    else:
        assert authority is not None and invocation_identity is not None
        losat_output = output_dir / f"matrix/{step.program}/{step.case_id}/losat.out"
        if not losat_output.is_file():
            raise CertificationFailure(
                f"direct oracle lacks verified matrix run 1: {step.step_id}"
            )
        native_fingerprint = verify_native_raw_fingerprint(
            authority,
            str(identity["platform_id"]),
            step.program,
            step.case_id,
            final_output,
        )
        native_vs_losat_diagnostic = report_native_vs_losat_diagnostic(
            catalog,
            step,
            final_output,
            losat_output,
            final_output.parent / "diagnostic",
        )
        verification = {
            "gate": "PLATFORM_NATIVE_NCBI_REFERENCE",
            "invocation_identity": invocation_identity,
            "native_fingerprint": native_fingerprint,
            "native_vs_losat_diagnostic": native_vs_losat_diagnostic,
        }
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
        "authority_version": identity.get("authority_version"),
        "authority_file_sha256": identity.get("authority_file_sha256"),
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
            native_fingerprint = verification["native_fingerprint"]
            native_vs_losat = verification["native_vs_losat_diagnostic"]
            oracle_rows.append(
                {
                    "program": step.program,
                    "case_id": step.case_id,
                    "semantic_class": step.semantic_class,
                    "classification": native_fingerprint["losat_parity_class"],
                    "native_vs_losat_diagnostic": native_vs_losat["classification"],
                    "platform_diagnostic": native_vs_losat["platform_diagnostic"],
                    "native_fingerprint_classification": native_fingerprint[
                        "classification"
                    ],
                    "native_ncbi_reference_class": native_fingerprint[
                        "native_ncbi_reference_class"
                    ],
                    "ncbi_sha256": record["output_sha256"],
                    "losat_sha256": step.expected_losat_sha256,
                    "status": "EXACT_AUTHORITY_FINGERPRINT",
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
            "native_vs_losat_diagnostic",
            "platform_diagnostic",
            "native_fingerprint_classification",
            "native_ncbi_reference_class",
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
    sakai = next(
        row
        for row in oracle_rows
        if row["program"] == "blastn" and row["case_id"] == "Sakai.MG1655.megablast"
    )
    d06 = next(
        row
        for row in oracle_rows
        if row["program"] == "tblastx" and row["case_id"] == "d06_ap027131_ap027133_db4"
    )
    if sakai["classification"] != "SOURCE_UNDETERMINED_ACCEPTED":
        raise CertificationFailure("Sakai LOSAT parity ratchet changed")
    if d06["classification"] != "HSP_SET_DIFF":
        raise CertificationFailure("d06 LOSAT parity ratchet changed")
    summary = {
        "decision": "PLATFORM_NATIVE_CERTIFIED",
        "platform_id": identity["platform_id"],
        "source_sha": identity["source_sha"],
        "authority_version": identity["authority_version"],
        "authority_file_sha256": identity["authority_file_sha256"],
        "completed_at": utc_now(),
        "search_executions": 61,
        "executed_this_attempt": executed_this_attempt,
        "resumed_verified_executions": 61 - executed_this_attempt,
        "gate_a_losat_canonical": {
            "total": 43,
            "passed": 43,
            "all_exact_pr5_raw_bytes": True,
        },
        "gate_b_platform_native_ncbi_reference": {
            "total": 6,
            "passed": 6,
            "all_exact_authority_fingerprints": True,
        },
        "cross_platform_losat": {"total": 43, "all_exact_bytes": True},
        "direct_ncbi": {
            "total": 6,
            "all_authority_accepted": True,
            "all_exact_native_fingerprints": True,
        },
        "repeatability": {
            "representatives": 6,
            "runs_per_representative": 3,
            "additional_executions": 12,
            "all_repeatable": True,
        },
        "sakai_ratchet": {
            "expected_losat_parity_class": "SOURCE_UNDETERMINED_ACCEPTED",
            "observed_losat_parity_class": sakai["classification"],
            "passed": True,
        },
        "tblastx_deviation_ratchet": {
            "accepted_deviation_ceiling": 6,
            "d06_expected_losat_contract": "approved_db_gencode_deviation",
            "d06_expected_losat_parity_class": "HSP_SET_DIFF",
            "d06_observed_losat_parity_class": d06["classification"],
            "passed": True,
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


def _strict_json_path(path: Path, context: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_bytes(), object_pairs_hook=_strict_json_object)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise CertificationFailure(f"invalid {context} JSON: {error}") from error
    return _object(value, context)


def verify_final_evidence_manifest(root: Path) -> str:
    manifest = root / "FINAL_EVIDENCE_MANIFEST.sha256"
    if not manifest.is_file():
        raise CertificationFailure(f"final evidence manifest is missing: {root}")
    recorded: dict[str, str] = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if "  " not in line:
            raise CertificationFailure(f"malformed final evidence manifest: {root}")
        digest, relative = line.split("  ", 1)
        if not SHA256_PATTERN.fullmatch(digest) or relative in recorded:
            raise CertificationFailure(
                f"invalid final evidence manifest row: {relative}"
            )
        path = _safe_resume_path(root, relative)
        if not path.is_file() or sha256_path(path) != digest:
            raise CertificationFailure(f"final evidence integrity mismatch: {relative}")
        recorded[relative] = digest
    observed = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path != manifest
    }
    if set(recorded) != observed:
        raise CertificationFailure(f"final evidence manifest coverage changed: {root}")
    return sha256_path(manifest)


def expected_completion_topology(
    catalog: Catalog,
) -> dict[str, dict[str, object]]:
    topology: dict[str, dict[str, object]] = {}

    def add(
        step_id: str,
        kind: str,
        program: str,
        case_id: str,
        semantic_class: str,
        run_index: int,
        output_rel: str,
    ) -> None:
        topology[step_id] = {
            "step_id": step_id,
            "execution_index": len(topology) + 1,
            "kind": kind,
            "program": program,
            "case_id": case_id,
            "semantic_class": semantic_class,
            "run_index": run_index,
            "output_rel": output_rel,
            "expected_losat_sha256": catalog.canonical[(program, case_id)][
                "losat_sha256"
            ],
        }

    for row in catalog.canonical_rows:
        add(
            f"matrix:{row['program']}:{row['case_id']}",
            "matrix",
            row["program"],
            row["case_id"],
            row["contract"],
            1,
            f"matrix/{row['program']}/{row['case_id']}/losat.out",
        )
    for program, case_id, semantic_class in REPRESENTATIVES:
        add(
            f"oracle:{program}:{case_id}",
            "oracle",
            program,
            case_id,
            semantic_class,
            1,
            f"oracle/{program}/{case_id}/ncbi.out",
        )
    for program, case_id, semantic_class in REPRESENTATIVES:
        for run_index in (2, 3):
            add(
                f"repeatability:{program}:{case_id}:run{run_index}",
                "repeatability",
                program,
                case_id,
                semantic_class,
                run_index,
                f"repeatability/{program}/{case_id}/run{run_index}.out",
            )
    if len(topology) != 61:
        raise CertificationFailure("aggregate expected topology is not 61 steps")
    return topology


def expected_native_invocation_evidence(
    authority: NativeAuthority,
    platform_id: str,
    program: str,
    case_id: str,
) -> dict[str, object]:
    platform_record, representative, _, _ = resolve_native_authority(
        authority, platform_id, program, case_id
    )
    archive = _object(platform_record["archive"], "native archive authority")
    checksum = _object(archive["checksum"], "native archive checksum")
    executable = _authority_executable(platform_record, program)
    inputs = {}
    for role in ("query", "subject"):
        expected_input = _object(representative[role], f"native {role} authority")
        inputs[role] = {
            "role": role,
            "repository_relative": expected_input["repository_relative"],
            "lexical_path": expected_input["lexical_path"],
            "sha256": expected_input["sha256"],
            "byte_length": expected_input["byte_length"],
            "newline_profile": expected_input["newline_profile"],
        }
    return {
        "authority_version": authority.authority_version,
        "authority_file_sha256": authority.file_sha256,
        "selector": {
            "authority_version": authority.authority_version,
            "platform_id": platform_id,
            "program": program,
            "case_id": case_id,
        },
        "ncbi_release": archive["release"],
        "archive_name": archive["filename"],
        "archive_checksum": checksum,
        "executable_sha256": executable["sha256"],
        "architecture": executable["architecture"],
        "inputs": inputs,
        "ordered_argv_template_sha256": hashlib.sha256(
            json.dumps(
                representative["ordered_argv_template"], separators=(",", ":")
            ).encode("utf-8")
        ).hexdigest(),
        "authoritative_outfmt": representative["authoritative_outfmt"],
        "working_directory_contract": "repository_root",
        "lexical_fixture_root": HISTORICAL_LEXICAL_ROOT,
        "shell": False,
        "environment_overrides": {},
    }


def expected_native_fingerprint_evidence(
    authority: NativeAuthority,
    platform_id: str,
    program: str,
    case_id: str,
) -> dict[str, object]:
    _, representative, expected, diagnostic = resolve_native_authority(
        authority, platform_id, program, case_id
    )
    return {
        "classification": "NATIVE_NCBI_FINGERPRINT_EXACT",
        "authority_version": authority.authority_version,
        "authority_file_sha256": authority.file_sha256,
        "losat_contract": representative["losat_contract"],
        "losat_parity_class": representative["losat_parity_class"],
        "native_ncbi_reference_class": diagnostic["native_ncbi_reference_class"],
        "raw_sha256": expected["raw_sha256"],
        "byte_length": expected["byte_length"],
        "newline_profile": expected["newline_profile"],
        "data_row_count": expected["data_row_count"],
    }


def aggregate_platform_evidence(
    input_root: Path,
    output: Path,
    expected_sha: str,
    authority: NativeAuthority,
    catalog: Catalog,
) -> dict[str, object]:
    if not SHA_PATTERN.fullmatch(expected_sha):
        raise CertificationFailure(
            "aggregate candidate SHA must be 40 lowercase hex digits"
        )
    validate_native_authority_catalog(authority, catalog)
    expected_topology = expected_completion_topology(catalog)
    expected_step_ids = set(expected_topology)
    expected_platforms = set(PLATFORMS)
    expected_directories = {
        f"platform-native-{platform_id}" for platform_id in expected_platforms
    }
    observed_directories = {path.name for path in input_root.iterdir() if path.is_dir()}
    if observed_directories != expected_directories:
        raise CertificationFailure(
            "aggregate must consume exactly the three expected platform artifacts"
        )
    platform_summaries = []
    for platform_id in sorted(expected_platforms):
        root = input_root / f"platform-native-{platform_id}"
        final_manifest_sha256 = verify_final_evidence_manifest(root)
        required_evidence = {
            "identity.json",
            "command_plan.json",
            "state.json",
            "summary.json",
            "cross_platform_equalities.tsv",
            "direct_ncbi_checks.tsv",
            "repeatability.tsv",
            "fixture-staging.json",
            "controlled-path-preflight.json",
        }
        observed_evidence = {
            path.relative_to(root).as_posix()
            for path in root.rglob("*")
            if path.is_file()
        }
        if not required_evidence <= observed_evidence:
            raise CertificationFailure(
                f"aggregate required evidence is incomplete: {platform_id}"
            )
        identity = _strict_json_path(root / "identity.json", "platform identity")
        command_plan = _strict_json_path(root / "command_plan.json", "command plan")
        summary = _strict_json_path(root / "summary.json", "platform summary")
        state = _strict_json_path(root / "state.json", "platform state")
        expected_identity = {
            "source_sha": expected_sha,
            "platform_id": platform_id,
            "authority_version": authority.authority_version,
            "authority_file_sha256": authority.file_sha256,
        }
        if any(identity.get(key) != value for key, value in expected_identity.items()):
            raise CertificationFailure(
                f"aggregate platform identity mismatch: {platform_id}"
            )
        if any(summary.get(key) != value for key, value in expected_identity.items()):
            raise CertificationFailure(
                f"aggregate platform summary mismatch: {platform_id}"
            )
        identity_losat = _object(identity.get("losat"), "platform LOSAT identity")
        authority_platform = authority.platforms[platform_id]
        identity_os = _object(identity.get("os"), "platform OS identity")
        recorded_archive = _object(
            identity.get("oracle_archive"), "platform archive identity"
        )
        expected_archive = _object(
            authority_platform["archive"], "authority archive identity"
        )
        expected_checksum = _object(
            expected_archive["checksum"], "authority archive checksum"
        )
        recorded_oracles = _object(
            identity.get("oracles"), "platform oracle identities"
        )
        if (
            not SHA256_PATTERN.fullmatch(str(identity_losat.get("sha256", "")))
            or identity_losat.get("architecture") != PLATFORMS[platform_id].binary_arch
            or identity_os.get("system") != authority_platform["platform_system"]
            or normalize_machine(str(identity_os.get("machine", "")))
            != authority_platform["normalized_machine"]
            or recorded_archive.get("archive_name") != expected_archive["filename"]
            or recorded_archive.get("archive_checksum_algorithm")
            != expected_checksum["algorithm"]
            or recorded_archive.get("archive_md5") != expected_checksum["value"]
            or set(recorded_oracles) != {"blastn", "blastp", "tblastx"}
        ):
            raise CertificationFailure(
                f"aggregate native identity failed: {platform_id}"
            )
        for program in ("blastn", "blastp", "tblastx"):
            recorded_oracle = _object(
                recorded_oracles[program], f"platform {program} identity"
            )
            expected_oracle = _authority_executable(authority_platform, program)
            if (
                recorded_oracle.get("sha256") != expected_oracle["sha256"]
                or recorded_oracle.get("architecture")
                != expected_oracle["architecture"]
            ):
                raise CertificationFailure(
                    f"aggregate native executable identity failed: "
                    f"{platform_id}/{program}"
                )
        if identity.get("execution_contract") != {
            "matrix": 43,
            "oracle": 6,
            "repeatability": 12,
            "per_platform": 61,
            "gate_a_losat_canonical": 43,
            "gate_b_native_ncbi_reference": 6,
        }:
            raise CertificationFailure(
                f"aggregate execution identity failed: {platform_id}"
            )
        if any(
            command_plan.get(key) != value
            for key, value in {
                "source_sha": expected_sha,
                "platform_id": platform_id,
                "authority_version": authority.authority_version,
                "authority_file_sha256": authority.file_sha256,
                "execution_count": 61,
                "controlled_fixture_root": HISTORICAL_LEXICAL_ROOT,
            }.items()
        ):
            raise CertificationFailure(
                f"aggregate command plan identity failed: {platform_id}"
            )
        planned_executions = _array(
            command_plan.get("executions"), "command plan executions"
        )
        if len(planned_executions) != 61:
            raise CertificationFailure(
                f"aggregate command plan topology failed: {platform_id}"
            )
        for planned, expected_step in zip(
            planned_executions, expected_topology.values(), strict=True
        ):
            planned_step = _object(planned, "planned execution")
            expected_environment = (
                [["RAYON_NUM_THREADS", "1"]]
                if expected_step["kind"] in {"matrix", "repeatability"}
                and expected_step["program"] == "tblastx"
                else []
            )
            if any(
                planned_step.get(key) != value
                for key, value in expected_step.items()
                if key != "output_rel"
            ) or (
                planned_step.get("output_rel")
                != f"{expected_step['output_rel']}.partial"
                or planned_step.get("final_output_rel") != expected_step["output_rel"]
                or planned_step.get("environment") != expected_environment
            ):
                raise CertificationFailure(
                    f"aggregate command plan step failed: {platform_id}"
                )
        if (
            state.get("platform_id") != platform_id
            or state.get("status") != "CERTIFIED"
            or state.get("completed_execution_count") != 61
            or state.get("expected_execution_count") != 61
            or state.get("completed_by_kind")
            != {"matrix": 43, "oracle": 6, "repeatability": 12}
            or state.get("completed_step_ids") != sorted(expected_step_ids)
            or summary.get("decision") != "PLATFORM_NATIVE_CERTIFIED"
            or summary.get("search_executions") != 61
        ):
            raise CertificationFailure(
                f"aggregate platform completion failed: {platform_id}"
            )
        if summary.get("gate_a_losat_canonical") != {
            "total": 43,
            "passed": 43,
            "all_exact_pr5_raw_bytes": True,
        }:
            raise CertificationFailure(f"aggregate Gate A failed: {platform_id}")
        if summary.get("gate_b_platform_native_ncbi_reference") != {
            "total": 6,
            "passed": 6,
            "all_exact_authority_fingerprints": True,
        }:
            raise CertificationFailure(f"aggregate Gate B failed: {platform_id}")
        repeatability = _object(summary.get("repeatability"), "repeatability summary")
        if repeatability != {
            "representatives": 6,
            "runs_per_representative": 3,
            "additional_executions": 12,
            "all_repeatable": True,
        }:
            raise CertificationFailure(f"aggregate repeatability failed: {platform_id}")
        if summary.get("sakai_ratchet") != {
            "expected_losat_parity_class": "SOURCE_UNDETERMINED_ACCEPTED",
            "observed_losat_parity_class": "SOURCE_UNDETERMINED_ACCEPTED",
            "passed": True,
        }:
            raise CertificationFailure(f"aggregate Sakai ratchet failed: {platform_id}")
        if summary.get("tblastx_deviation_ratchet") != {
            "accepted_deviation_ceiling": 6,
            "d06_expected_losat_contract": "approved_db_gencode_deviation",
            "d06_expected_losat_parity_class": "HSP_SET_DIFF",
            "d06_observed_losat_parity_class": "HSP_SET_DIFF",
            "passed": True,
        }:
            raise CertificationFailure(
                f"aggregate TBLASTX deviation ratchet failed: {platform_id}"
            )
        completion_paths = sorted((root / "completions").glob("*.json"))
        expected_completion_names = {
            f"{step_id.replace(':', '__')}.json" for step_id in expected_step_ids
        }
        if (
            len(completion_paths) != 61
            or {path.name for path in completion_paths} != expected_completion_names
        ):
            raise CertificationFailure(
                f"aggregate completion count failed: {platform_id}"
            )
        kinds: Counter[str] = Counter()
        gate_a = 0
        gate_b = 0
        repeatability_count = 0
        observed_step_ids: set[str] = set()
        for completion_path in completion_paths:
            completion = _strict_json_path(completion_path, "completion")
            step_id = str(completion.get("step_id", ""))
            expected_step = expected_topology.get(step_id)
            if expected_step is None or step_id in observed_step_ids:
                raise CertificationFailure(
                    f"aggregate exact step topology failed: {completion_path.name}"
                )
            observed_step_ids.add(step_id)
            if (
                completion.get("status") != "VERIFIED"
                or completion.get("source_sha") != expected_sha
                or completion.get("platform_id") != platform_id
                or completion.get("authority_version") != authority.authority_version
                or completion.get("authority_file_sha256") != authority.file_sha256
                or completion.get("losat_binary_sha256") != identity_losat["sha256"]
                or any(
                    completion.get(key) != value for key, value in expected_step.items()
                )
            ):
                raise CertificationFailure(
                    f"aggregate completion identity failed: {completion_path.name}"
                )
            output_rel = str(completion.get("output_rel", ""))
            output_path = _safe_resume_path(root, output_rel)
            if not output_path.is_file() or sha256_path(output_path) != completion.get(
                "output_sha256"
            ):
                raise CertificationFailure(
                    f"aggregate completion output failed: {completion_path.name}"
                )
            kind = str(completion.get("kind"))
            kinds[kind] += 1
            verification = _object(
                completion.get("verification"), "completion verification"
            )
            if kind == "matrix":
                expected_verification = {
                    "gate": "LOSAT_CANONICAL",
                    "classification": "CANONICAL_PR5_RAW_BYTES",
                    "canonical_sha256": expected_step["expected_losat_sha256"],
                }
                if (
                    verification != expected_verification
                    or completion.get("output_sha256")
                    != expected_step["expected_losat_sha256"]
                ):
                    raise CertificationFailure(
                        f"aggregate Gate A completion failed: {completion_path.name}"
                    )
                gate_a += 1
            elif kind == "oracle":
                fingerprint = _object(
                    verification.get("native_fingerprint"), "native fingerprint"
                )
                _object(
                    verification.get("native_vs_losat_diagnostic"),
                    "native-versus-LOSAT diagnostic",
                )
                program = str(expected_step["program"])
                case_id = str(expected_step["case_id"])
                expected_invocation = expected_native_invocation_evidence(
                    authority, platform_id, program, case_id
                )
                expected_fingerprint = expected_native_fingerprint_evidence(
                    authority, platform_id, program, case_id
                )
                if (
                    set(verification)
                    != {
                        "gate",
                        "invocation_identity",
                        "native_fingerprint",
                        "native_vs_losat_diagnostic",
                    }
                    or verification.get("gate") != "PLATFORM_NATIVE_NCBI_REFERENCE"
                    or verification.get("invocation_identity") != expected_invocation
                    or fingerprint != expected_fingerprint
                ):
                    raise CertificationFailure(
                        f"aggregate Gate B completion failed: {completion_path.name}"
                    )
                if (
                    completion.get("output_sha256")
                    != expected_fingerprint["raw_sha256"]
                ):
                    raise CertificationFailure(
                        f"aggregate native raw fingerprint failed: {completion_path.name}"
                    )
                gate_b += 1
            elif kind == "repeatability":
                expected_verification = {
                    "gate": "LOSAT_REPEATABILITY",
                    "classification": "CANONICAL_PR5_RAW_BYTES",
                    "canonical_sha256": expected_step["expected_losat_sha256"],
                }
                if (
                    verification != expected_verification
                    or completion.get("output_sha256")
                    != expected_step["expected_losat_sha256"]
                ):
                    raise CertificationFailure(
                        f"aggregate repeatability completion failed: {completion_path.name}"
                    )
                repeatability_count += 1
            else:
                raise CertificationFailure(
                    f"aggregate found unknown completion kind: {kind}"
                )
        if kinds != Counter({"matrix": 43, "oracle": 6, "repeatability": 12}):
            raise CertificationFailure(
                f"aggregate execution topology failed: {platform_id}"
            )
        if (gate_a, gate_b, repeatability_count) != (43, 6, 12):
            raise CertificationFailure(
                f"aggregate subgate totals failed: {platform_id}"
            )
        if observed_step_ids != expected_step_ids:
            raise CertificationFailure(
                f"aggregate exact step set failed: {platform_id}"
            )
        platform_summaries.append(
            {
                "platform_id": platform_id,
                "source_sha": expected_sha,
                "authority_version": authority.authority_version,
                "authority_file_sha256": authority.file_sha256,
                "search_executions": 61,
                "gate_a_passed": 43,
                "gate_b_passed": 6,
                "repeatability_passed": 12,
                "final_evidence_manifest_sha256": final_manifest_sha256,
            }
        )
    aggregate = {
        "decision": "CROSS_PLATFORM_NATIVE_CERTIFIED",
        "created_at": utc_now(),
        "source_sha": expected_sha,
        "authority_version": authority.authority_version,
        "authority_file_sha256": authority.file_sha256,
        "expected_platform_ids": sorted(expected_platforms),
        "platform_count": 3,
        "search_executions_per_platform": 61,
        "total_search_executions": 183,
        "platforms": platform_summaries,
    }
    atomic_write_json(output, aggregate)
    return aggregate


def certify(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    authority = load_native_authority(args.authority.resolve())
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    preexisting = {path.name for path in output_dir.iterdir()}
    if preexisting - {"workflow-started.txt", "controlled-path-preflight.json"}:
        raise CertificationFailure(
            f"output directory must be new or empty: {output_dir}"
        )
    platform_spec = PLATFORMS[args.platform_id]
    authority_platform = validate_native_platform_sanity(
        authority, platform_spec.platform_id
    )
    initial_state = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "created_at": utc_now(),
        "platform_id": platform_spec.platform_id,
        "authority_version": authority.authority_version,
        "authority_file_sha256": authority.file_sha256,
        "status": "INITIALIZING",
        "expected_execution_count": 61,
    }
    atomic_write_json(output_dir / "state.json", initial_state)
    preflight_controlled_fixture_path(output_dir)
    certified_sha = validate_git_identity(repo_root, args.expected_sha)
    archive_record = verify_archive(
        platform_spec,
        args.oracle_archive.resolve(),
        args.oracle_checksum.resolve(),
        authority_platform,
    )
    oracles = {
        "blastn": args.blastn_oracle.resolve(),
        "blastp": args.blastp_oracle.resolve(),
        "tblastx": args.tblastx_oracle.resolve(),
    }
    canonical_manifest = (
        repo_root / "LOSAT" / "tests" / "platform_native_v010_canonical.tsv"
    )
    catalog = load_catalog(repo_root)
    validate_native_authority_catalog(authority, catalog)
    steps = build_steps(
        repo_root,
        output_dir,
        catalog,
        args.native_bin.resolve(),
        oracles,
    )
    staged_fixtures = stage_required_fixtures(repo_root, steps, output_dir)
    if len(staged_fixtures) != EXPECTED_STAGED_FIXTURE_COUNT:
        raise CertificationFailure(
            "controlled fixture set changed: "
            f"expected {EXPECTED_STAGED_FIXTURE_COUNT}, found {len(staged_fixtures)}"
        )
    controlled_fixtures = fixture_staging_identity(staged_fixtures)
    identity = record_identity(
        repo_root,
        output_dir,
        certified_sha,
        platform_spec,
        args.native_bin.resolve(),
        archive_record,
        oracles,
        canonical_manifest,
        controlled_fixtures,
        authority,
    )
    atomic_write_json(
        output_dir / "command_plan.json",
        {
            "source_sha": certified_sha,
            "platform_id": platform_spec.platform_id,
            "authority_version": authority.authority_version,
            "authority_file_sha256": authority.file_sha256,
            "created_at": utc_now(),
            "execution_count": len(steps),
            "controlled_fixture_root": HISTORICAL_LEXICAL_ROOT,
            "staged_fixture_count": len(staged_fixtures),
            "executions": [asdict(step) for step in steps],
        },
    )
    completed: dict[str, dict[str, object]] = {}
    if args.resume_dir is not None:
        completed = import_resume(
            args.resume_dir.resolve(),
            output_dir,
            repo_root,
            identity,
            steps,
            catalog,
            authority,
        )
    _write_state(output_dir, platform_spec.platform_id, completed, status="RUNNING")
    executed_this_attempt = 0
    for step in steps:
        if step.step_id in completed:
            print(
                f"[platform-certification] resume {step.execution_index}/61 {step.step_id}"
            )
            continue
        print(
            f"[platform-certification] run {step.execution_index}/61 {step.step_id}",
            flush=True,
        )
        record = _execute_step(
            repo_root, output_dir, identity, catalog, step, authority
        )
        completed[step.step_id] = record
        executed_this_attempt += 1
        _write_state(output_dir, platform_spec.platform_id, completed, status="RUNNING")
    finalize_evidence(output_dir, identity, steps, completed, executed_this_attempt)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]
    default_authority = script_path.with_name("ncbi_platform_variance_v010.json")
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    archive_parser = subparsers.add_parser(
        "verify-archive", help="verify a pinned official NCBI archive"
    )
    archive_parser.add_argument(
        "--platform-id", choices=sorted(PLATFORMS), required=True
    )
    archive_parser.add_argument("--archive", type=Path, required=True)
    archive_parser.add_argument("--checksum", type=Path, required=True)
    archive_parser.add_argument("--authority", type=Path, default=default_authority)

    extraction_parser = subparsers.add_parser(
        "extract-verified-archive",
        help="extract a pinned official NCBI archive with Python tarfile",
    )
    extraction_parser.add_argument(
        "--platform-id", choices=sorted(PLATFORMS), required=True
    )
    extraction_parser.add_argument("--archive", type=Path, required=True)
    extraction_parser.add_argument("--checksum", type=Path, required=True)
    extraction_parser.add_argument("--destination", type=Path, required=True)
    extraction_parser.add_argument("--authority", type=Path, default=default_authority)

    preflight_parser = subparsers.add_parser(
        "preflight-path", help="verify the controlled historical fixture path"
    )
    preflight_parser.add_argument("--output-dir", type=Path, required=True)

    certify_parser = subparsers.add_parser(
        "certify", help="run or resume the exact 61-execution platform gate"
    )
    certify_parser.add_argument("--repo-root", type=Path, default=repo_root)
    certify_parser.add_argument("--authority", type=Path, default=default_authority)
    certify_parser.add_argument("--expected-sha", required=True)
    certify_parser.add_argument(
        "--platform-id", choices=sorted(PLATFORMS), required=True
    )
    certify_parser.add_argument("--output-dir", type=Path, required=True)
    certify_parser.add_argument("--resume-dir", type=Path)
    certify_parser.add_argument("--native-bin", type=Path, required=True)
    certify_parser.add_argument("--oracle-archive", type=Path, required=True)
    certify_parser.add_argument("--oracle-checksum", type=Path, required=True)
    certify_parser.add_argument("--blastn-oracle", type=Path, required=True)
    certify_parser.add_argument("--blastp-oracle", type=Path, required=True)
    certify_parser.add_argument("--tblastx-oracle", type=Path, required=True)
    aggregate_parser = subparsers.add_parser(
        "aggregate", help="aggregate exactly three completed platform artifacts"
    )
    aggregate_parser.add_argument("--input-root", type=Path, required=True)
    aggregate_parser.add_argument("--output", type=Path, required=True)
    aggregate_parser.add_argument("--expected-sha", required=True)
    aggregate_parser.add_argument("--authority", type=Path, default=default_authority)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = getattr(args, "output_dir", None)
    try:
        if args.mode == "verify-archive":
            authority = load_native_authority(args.authority.resolve())
            record = verify_archive(
                PLATFORMS[args.platform_id],
                args.archive.resolve(),
                args.checksum.resolve(),
                authority.platforms[args.platform_id],
            )
            print(json.dumps(record, sort_keys=True))
        elif args.mode == "extract-verified-archive":
            authority = load_native_authority(args.authority.resolve())
            record = extract_verified_archive(
                PLATFORMS[args.platform_id],
                args.archive.resolve(),
                args.checksum.resolve(),
                args.destination.absolute(),
                authority.platforms[args.platform_id],
            )
            print(json.dumps(record, sort_keys=True))
        elif args.mode == "preflight-path":
            record = preflight_controlled_fixture_path(args.output_dir.resolve())
            print(json.dumps(record, sort_keys=True))
        elif args.mode == "certify":
            certify(args)
            print("PLATFORM_NATIVE_CERTIFIED")
        else:
            authority = load_native_authority(args.authority.resolve())
            catalog = load_catalog(Path(__file__).resolve().parents[2])
            aggregate_platform_evidence(
                args.input_root.resolve(),
                args.output.resolve(),
                args.expected_sha,
                authority,
                catalog,
            )
            print("CROSS_PLATFORM_NATIVE_CERTIFIED")
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
