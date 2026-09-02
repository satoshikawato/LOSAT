#!/usr/bin/env python3
"""Characterize official NCBI BLAST+ 2.17.0 output on three native platforms.

This driver is diagnostic-only.  On each platform it executes the six retained
direct-oracle representatives as six authoritative/diagnostic pairs, never
executes LOSAT, and never changes an authority classification.

NCBI references:
- ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-45
  defines the twelve default tabular fields used by outfmt 6 and 7.
- ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
  emits those fields in order and terminates each row without reordering it.
- ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-832
  iterates the retained alignment set in formatter order.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import subprocess
import sys
from types import ModuleType
from typing import Iterable, Sequence


PR6_CERT_SHA_V3 = "0e29e2201e2d1b03124c0b9d6698a81bfed8cec0"
PR5_EVIDENCE_MANIFEST_SHA256 = (
    "b9fc98a376d2849274c86b4e4769d2ee38b76025adbb4d63d3da4a5e3e7cdb5c"
)
PR5_EVIDENCE_ROOT = Path(
    "/mnt/c/Users/genom/LOSAT-certification-evidence/"
    "losat-pr5-integrated-5845d22-20260831-resumed-final"
)
EVIDENCE_SCHEMA = 2
EXPECTED_FIXTURE_COUNT = 10
EXPECTED_REPRESENTATIVE_COUNT = 6
EXPECTED_PLATFORM_AUTHORITATIVE_EXECUTIONS = 6
EXPECTED_PLATFORM_DIAGNOSTIC_EXECUTIONS = 6
EXPECTED_PLATFORM_EXECUTIONS = 12
EXPECTED_TOTAL_AUTHORITATIVE_EXECUTIONS = 18
EXPECTED_TOTAL_DIAGNOSTIC_EXECUTIONS = 18
EXPECTED_TOTAL_EXECUTIONS = 36
SHA40 = re.compile(r"^[0-9a-f]{40}$")
ALLOWED_DIAGNOSTIC_CHANGES = {
    ".github/workflows/ncbi-only-platform-characterization.yml",
    "LOSAT/tests/characterize_ncbi_platform_variance_v010.py",
    "LOSAT/tests/test_characterize_ncbi_platform_variance_v010.py",
}


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-45
# ```c++
# const char* kDfltArgTabularOutputFmt =
#     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
#     "evalue bitscore";
# ```
STANDARD_FIELDS = (
    "query_id",
    "subject_id",
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
KEY_FIELDS = ("query_id", "subject_id", "qstart", "qend", "sstart", "send")
# NCBI references:
# - ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-41
#   defines ``std`` as the twelve authoritative tabular fields above.
# - ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:96-101, 114-119,
#   144-146 defines qseq, sseq, score, and btop.
# - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:59-91 expands ``std``
#   in place and appends explicitly requested fields in requested order.
DIAGNOSTIC_OUTFMT = "6 std score btop qseq sseq"
DIAGNOSTIC_FIELDS = (*STANDARD_FIELDS, "score", "btop", "qseq", "sseq")


@dataclass(frozen=True)
class Representative:
    program: str
    case_id: str
    semantic_class: str
    linux_relative: str


REPRESENTATIVES = (
    Representative(
        "blastn",
        "PesePMNV.MjPMNV.task_blastn",
        "ordinary_exact",
        "native/blastn/PesePMNV.MjPMNV.task_blastn.ncbi.out",
    ),
    Representative(
        "blastn",
        "Sakai.MG1655.megablast",
        "source_undetermined",
        "native/blastn/Sakai.MG1655.megablast.ncbi.out",
    ),
    Representative(
        "blastp",
        "pairwise_default_serial",
        "representative",
        "native/blastp/pairwise_default_serial.ncbi.tsv",
    ),
    Representative(
        "tblastx",
        "p03_mela_pemojnva",
        "ordinary_exact",
        "native/tblastx/p03_mela_pemojnva/ncbi.tsv",
    ),
    Representative(
        "tblastx",
        "p14_ap027131_ap027133_query4",
        "query_gencode_exact",
        "native/tblastx/p14_ap027131_ap027133_query4/ncbi.tsv",
    ),
    Representative(
        "tblastx",
        "d06_ap027131_ap027133_db4",
        "db_gencode_deviation",
        "native/tblastx/d06_ap027131_ap027133_db4/ncbi.tsv",
    ),
)


class CharacterizationFailure(RuntimeError):
    """The diagnostic contract failed closed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: object) -> None:
    atomic_write_bytes(
        path, (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise CharacterizationFailure(f"cannot load authority module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_platform_authority(repo_root: Path) -> tuple[ModuleType, object]:
    path = repo_root / "LOSAT/tests/certify_platform_native_v010.py"
    authority = load_module("ncbi_characterization_platform_authority", path)
    catalog = authority.load_catalog(repo_root)
    observed = tuple(
        Representative(program, case_id, semantic_class, "")
        for program, case_id, semantic_class in authority.REPRESENTATIVES
    )
    expected = tuple(
        Representative(item.program, item.case_id, item.semantic_class, "")
        for item in REPRESENTATIVES
    )
    if observed != expected:
        raise CharacterizationFailure(
            "the retained six-representative authority changed"
        )
    return authority, catalog


def run_bytes(command: Sequence[str], cwd: Path) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        list(command),
        cwd=cwd,
        capture_output=True,
        text=False,
        check=False,
        shell=False,
    )


def git_output(repo_root: Path, arguments: Sequence[str]) -> bytes:
    completed = run_bytes(["git", *arguments], repo_root)
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise CharacterizationFailure(f"git {' '.join(arguments)} failed: {detail}")
    return completed.stdout


def verify_candidate(repo_root: Path, candidate_sha: str) -> None:
    if not SHA40.fullmatch(candidate_sha):
        raise CharacterizationFailure("candidate SHA must be 40 lowercase hex digits")
    head = git_output(repo_root, ["rev-parse", "HEAD"]).decode("ascii").strip()
    if head != candidate_sha:
        raise CharacterizationFailure(
            f"checked-out candidate mismatch: expected {candidate_sha}, observed {head}"
        )
    ancestor = run_bytes(
        ["git", "merge-base", "--is-ancestor", PR6_CERT_SHA_V3, candidate_sha],
        repo_root,
    )
    if ancestor.returncode != 0:
        raise CharacterizationFailure("candidate is not descended from exact PR 6 V3")
    changed = set(
        git_output(
            repo_root,
            ["diff", "--name-only", f"{PR6_CERT_SHA_V3}..{candidate_sha}"],
        )
        .decode("utf-8")
        .splitlines()
    )
    unexpected = changed - ALLOWED_DIAGNOSTIC_CHANGES
    if unexpected:
        raise CharacterizationFailure(
            "diagnostic candidate changes protected paths: "
            + ", ".join(sorted(unexpected))
        )


def newline_profile(data: bytes) -> dict[str, int]:
    crlf = data.count(b"\r\n")
    without_crlf = data.replace(b"\r\n", b"")
    return {
        "lf_count": data.count(b"\n"),
        "bare_lf_count": data.count(b"\n") - crlf,
        "crlf_count": crlf,
        "bare_cr_count": without_crlf.count(b"\r"),
    }


def lf_diagnostic_bytes(data: bytes) -> bytes:
    return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def require_lf_fixture(data: bytes, label: str) -> dict[str, int]:
    """Fail before search planning unless a staged fixture is LF-only."""

    profile = newline_profile(data)
    if (
        profile["lf_count"] == 0
        or profile["crlf_count"] != 0
        or profile["bare_cr_count"] != 0
        or profile["bare_lf_count"] != profile["lf_count"]
    ):
        raise CharacterizationFailure(f"fixture is not LF-only: {label}")
    return profile


def _replace_option(command: Sequence[str], option: str, value: str) -> list[str]:
    updated = list(command)
    if updated.count(option) != 1:
        raise CharacterizationFailure(f"command must contain one {option}: {command}")
    index = updated.index(option)
    if index + 1 >= len(updated):
        raise CharacterizationFailure(f"command lacks value for {option}: {command}")
    updated[index + 1] = value
    return updated


def _option_value(command: Sequence[str], option: str) -> str:
    if command.count(option) != 1:
        raise CharacterizationFailure(f"command must contain one {option}: {command}")
    index = command.index(option)
    if index + 1 >= len(command):
        raise CharacterizationFailure(f"command lacks value for {option}: {command}")
    return command[index + 1]


def _normalized_command_without_outfmt(command: Sequence[str]) -> list[str]:
    return _replace_option(command, "-outfmt", "<OUTFMT>")


def compare_pair_commands(
    authoritative: Sequence[str], diagnostic: Sequence[str]
) -> dict[str, object]:
    """Prove that a search pair differs only in the outfmt value.

    NCBI references:
    - ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:35-41 names
      ``outfmt`` and defines the default tabular field set.
    - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:59-91 expands the
      requested tabular field tokens without changing search inputs/options.
    """

    authoritative_outfmt = _option_value(authoritative, "-outfmt")
    diagnostic_outfmt = _option_value(diagnostic, "-outfmt")
    authoritative_normalized = _normalized_command_without_outfmt(authoritative)
    diagnostic_normalized = _normalized_command_without_outfmt(diagnostic)
    differences = [
        {
            "argument_index": index,
            "authoritative": authoritative_value,
            "diagnostic": diagnostic_value,
        }
        for index, (authoritative_value, diagnostic_value) in enumerate(
            zip(authoritative, diagnostic, strict=True)
        )
        if authoritative_value != diagnostic_value
    ] if len(authoritative) == len(diagnostic) else []
    identical_non_output_arguments = (
        len(authoritative) == len(diagnostic)
        and authoritative_normalized == diagnostic_normalized
        and len(differences) == 1
        and differences[0]["argument_index"]
        == authoritative.index("-outfmt") + 1
    )
    if not identical_non_output_arguments:
        raise CharacterizationFailure(
            "authoritative/diagnostic commands differ outside -outfmt"
        )
    if diagnostic_outfmt != DIAGNOSTIC_OUTFMT:
        raise CharacterizationFailure(
            f"diagnostic outfmt changed: {diagnostic_outfmt!r}"
        )
    if authoritative_outfmt == diagnostic_outfmt:
        raise CharacterizationFailure("diagnostic outfmt did not change")
    return {
        "identical_non_output_arguments": True,
        "excluded_output_format_option": "-outfmt",
        "authoritative_outfmt": authoritative_outfmt,
        "diagnostic_outfmt": diagnostic_outfmt,
        "argument_differences": differences,
        "normalized_authoritative_argv": authoritative_normalized,
        "normalized_diagnostic_argv": diagnostic_normalized,
    }


def _input_options(program: str) -> dict[str, tuple[str, ...]]:
    if program == "blastn":
        return {"query": ("-query", "-q"), "subject": ("-subject", "-s")}
    return {"query": ("-query", "--query"), "subject": ("-subject", "--subject")}


def input_arguments(command: Sequence[str], program: str) -> list[tuple[str, str, str]]:
    result: list[tuple[str, str, str]] = []
    for role, aliases in _input_options(program).items():
        present = [option for option in aliases if command.count(option) == 1]
        if len(present) != 1:
            raise CharacterizationFailure(
                f"{program} command lacks one {role}: {command}"
            )
        option = present[0]
        index = command.index(option)
        if index + 1 >= len(command):
            raise CharacterizationFailure(f"{program} command lacks {option} value")
        result.append((role, option, command[index + 1]))
    return result


def _authoritative_command(
    repo_root: Path,
    authority: ModuleType,
    catalog: object,
    representative: Representative,
    oracle: Path,
    output: Path,
) -> list[str]:
    program = representative.program
    case_id = representative.case_id
    if program == "blastn":
        row = authority._blastn_command_row(repo_root, catalog.blastn_rows[case_id])
        command = catalog.blastn_compare.build_ncbi_command(row, str(oracle), output)
    elif program == "blastp":
        command, _ = catalog.blastp_audit.build_commands(
            catalog.blastp_cases[case_id],
            oracle,
            Path("UNUSED_LOSAT_BINARY_MUST_NOT_EXECUTE"),
            output.parent,
        )
    elif program == "tblastx":
        command, _ = catalog.tblastx_audit.build_commands(
            catalog.tblastx_cases[case_id],
            oracle,
            Path("UNUSED_LOSAT_BINARY_MUST_NOT_EXECUTE"),
            output,
            output.parent / "UNUSED_LOSAT_OUTPUT_MUST_NOT_EXIST",
        )
    else:
        raise CharacterizationFailure(f"unsupported representative program: {program}")
    return authority.replace_output(command, output)


def _repo_relative(repo_root: Path, source: str) -> PurePosixPath:
    try:
        relative = Path(source).resolve().relative_to(repo_root.resolve())
    except ValueError as error:
        raise CharacterizationFailure(
            f"fixture escaped repository: {source}"
        ) from error
    result = PurePosixPath(*relative.parts)
    if not result.parts or result.parts[0] != "LOSAT":
        raise CharacterizationFailure(f"fixture is outside LOSAT: {result}")
    return result


def _historical_lexical_path(authority: ModuleType, relative: PurePosixPath) -> str:
    within_crate = PurePosixPath(*relative.parts[1:])
    return authority.historical_lexical_path(within_crate)


def build_command_plan(
    repo_root: Path,
    authority: ModuleType,
    catalog: object,
    oracles: dict[str, Path],
    evidence_root: Path,
) -> list[dict[str, object]]:
    plan: list[dict[str, object]] = []
    for pair_index, representative in enumerate(REPRESENTATIVES, start=1):
        case_dir = evidence_root / representative.program / representative.case_id
        # Both NCBI searches use one transient sink.  Each completed file is
        # moved byte-for-byte into its role-specific directory before the next
        # search.  Consequently the actual argv differs only at -outfmt.
        execution_output = case_dir / "ncbi.execution.raw.out"
        authoritative_command = _authoritative_command(
            repo_root,
            authority,
            catalog,
            representative,
            oracles[representative.program],
            execution_output,
        )
        inputs: dict[str, dict[str, str]] = {}
        for role, option, source in input_arguments(
            authoritative_command, representative.program
        ):
            relative = _repo_relative(repo_root, source)
            lexical = _historical_lexical_path(authority, relative)
            authoritative_command[authoritative_command.index(option) + 1] = lexical
            inputs[role] = {
                "option": option,
                "repository_relative": relative.as_posix(),
                "lexical_path": lexical,
            }
        diagnostic_command = _replace_option(
            authoritative_command, "-outfmt", DIAGNOSTIC_OUTFMT
        )
        pair_comparison = compare_pair_commands(
            authoritative_command, diagnostic_command
        )

        searches: dict[str, dict[str, object]] = {}
        for offset, (role, command) in enumerate(
            (
                ("authoritative", authoritative_command),
                ("diagnostic", diagnostic_command),
            )
        ):
            if (len(command) - 1) % 2:
                raise CharacterizationFailure(
                    f"unexpected BLAST argument topology: {command}"
                )
            option_pairs = [
                {"option": command[position], "value": command[position + 1]}
                for position in range(1, len(command), 2)
            ]
            if any(not item["option"].startswith("-") for item in option_pairs):
                raise CharacterizationFailure(
                    f"unexpected positional BLAST argument: {command}"
                )
            option_map = {
                str(item["option"]): str(item["value"]) for item in option_pairs
            }
            searches[role] = {
                "execution_index": (pair_index - 1) * 2 + offset + 1,
                "search_role": role,
                "argv": command,
                "option_pairs": option_pairs,
                "semantics": {
                    "task": option_map.get("-task", representative.program),
                    "outfmt": option_map["-outfmt"],
                    "threads": option_map["-num_threads"],
                    "query_gencode": option_map.get("-query_gencode", "1"),
                    "db_gencode": option_map.get("-db_gencode", "1"),
                },
            }
        plan.append(
            {
                "pair_index": pair_index,
                "program": representative.program,
                "case_id": representative.case_id,
                "semantic_class": representative.semantic_class,
                "inputs": inputs,
                "authoritative": searches["authoritative"],
                "diagnostic": searches["diagnostic"],
                "pair_command_comparison": pair_comparison,
            }
        )
    if len(plan) != EXPECTED_REPRESENTATIVE_COUNT:
        raise CharacterizationFailure(
            "execution plan is not exactly six representatives"
        )
    authoritative_count = sum("authoritative" in pair for pair in plan)
    diagnostic_count = sum("diagnostic" in pair for pair in plan)
    if (
        authoritative_count != EXPECTED_PLATFORM_AUTHORITATIVE_EXECUTIONS
        or diagnostic_count != EXPECTED_PLATFORM_DIAGNOSTIC_EXECUTIONS
        or authoritative_count + diagnostic_count != EXPECTED_PLATFORM_EXECUTIONS
    ):
        raise CharacterizationFailure("paired execution count ratchet failed")
    return plan


def _git_blob(
    repo_root: Path,
    relative: PurePosixPath,
    source_sha: str = PR6_CERT_SHA_V3,
) -> tuple[str, bytes]:
    spec = f"{source_sha}:{relative.as_posix()}"
    oid = git_output(repo_root, ["rev-parse", spec]).decode("ascii").strip()
    data = git_output(repo_root, ["show", spec])
    return oid, data


def stage_fixtures(
    repo_root: Path, candidate_sha: str, platform_id: str, output_dir: Path
) -> None:
    verify_candidate(repo_root, candidate_sha)
    authority, catalog = load_platform_authority(repo_root)
    if platform_id not in authority.PLATFORMS:
        raise CharacterizationFailure(f"unsupported platform: {platform_id}")
    placeholders = {
        program: Path(f"${{NCBI_BIN}}/{program}")
        for program in ("blastn", "blastp", "tblastx")
    }
    plan = build_command_plan(
        repo_root, authority, catalog, placeholders, Path("${EVIDENCE_ROOT}")
    )
    fixture_relatives = sorted(
        {
            PurePosixPath(str(details["repository_relative"]))
            for pair in plan
            for details in pair["inputs"].values()
        },
        key=PurePosixPath.as_posix,
    )
    if len(fixture_relatives) != EXPECTED_FIXTURE_COUNT:
        raise CharacterizationFailure(
            f"expected {EXPECTED_FIXTURE_COUNT} unique fixtures, found {len(fixture_relatives)}"
        )
    records = []
    for relative in fixture_relatives:
        oid, data = _git_blob(repo_root, relative)
        profile = require_lf_fixture(data, relative.as_posix())
        staged_relative = PurePosixPath("fixtures") / relative
        staged = output_dir.joinpath(*staged_relative.parts)
        atomic_write_bytes(staged, data)
        observed = staged.read_bytes()
        if observed != data:
            raise CharacterizationFailure(
                f"binary staging changed {relative.as_posix()}"
            )
        records.append(
            {
                "repository_relative": relative.as_posix(),
                "git_blob_oid": oid,
                "expected_git_blob_raw_sha256": sha256_bytes(data),
                "staged_relative": staged_relative.as_posix(),
                "staged_sha256": sha256_bytes(observed),
                "byte_length": len(data),
                "newline_profile": profile,
                "lf_only": True,
            }
        )
    manifest = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "status": "FIXTURES_STAGED_FROM_GIT_OBJECTS",
        "created_at": utc_now(),
        "platform_id": platform_id,
        "candidate_sha": candidate_sha,
        "fixture_source_sha": PR6_CERT_SHA_V3,
        "fixture_count": len(records),
        "inputs": records,
    }
    atomic_write_json(output_dir / "fixture_manifest.json", manifest)
    atomic_write_json(
        output_dir / "command_plan.json",
        {
            "evidence_schema": EVIDENCE_SCHEMA,
            "status": "PLANNED_NO_EXECUTION",
            "representative_count": len(plan),
            "authoritative_count": EXPECTED_PLATFORM_AUTHORITATIVE_EXECUTIONS,
            "diagnostic_count": EXPECTED_PLATFORM_DIAGNOSTIC_EXECUTIONS,
            "execution_count": EXPECTED_PLATFORM_EXECUTIONS,
            "pairs": plan,
        },
    )


def _invariant_input(record: dict[str, object]) -> dict[str, object]:
    return {
        key: record[key]
        for key in (
            "repository_relative",
            "git_blob_oid",
            "expected_git_blob_raw_sha256",
            "staged_sha256",
            "byte_length",
            "newline_profile",
            "lf_only",
        )
    }


def verify_fixture_invariance(manifests: Sequence[Path], output: Path) -> None:
    if len(manifests) != 3:
        raise CharacterizationFailure(
            "fixture invariance requires exactly three manifests"
        )
    loaded = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
    platform_ids = {manifest["platform_id"] for manifest in loaded}
    if platform_ids != {"windows-x64", "macos-arm64", "macos-x64"}:
        raise CharacterizationFailure(
            f"fixture platform topology changed: {platform_ids}"
        )
    candidate_shas = {manifest["candidate_sha"] for manifest in loaded}
    if len(candidate_shas) != 1 or not SHA40.fullmatch(next(iter(candidate_shas))):
        raise CharacterizationFailure(
            "fixture manifests do not share one candidate SHA"
        )
    if any(manifest["fixture_source_sha"] != PR6_CERT_SHA_V3 for manifest in loaded):
        raise CharacterizationFailure("fixture source is not exact PR 6 V3")
    baseline = [_invariant_input(item) for item in loaded[0]["inputs"]]
    for manifest in loaded[1:]:
        observed = [_invariant_input(item) for item in manifest["inputs"]]
        if observed != baseline:
            raise CharacterizationFailure(
                f"cross-platform fixture bytes differ for {manifest['platform_id']}"
            )
    record = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "status": "CROSS_PLATFORM_FIXTURES_BYTE_IDENTICAL",
        "verified_at": utc_now(),
        "platform_ids": sorted(platform_ids),
        "candidate_sha": next(iter(candidate_shas)),
        "fixture_source_sha": PR6_CERT_SHA_V3,
        "fixture_count": len(baseline),
        "inputs": baseline,
        "manifest_sha256": {
            manifest["platform_id"]: sha256_path(path)
            for manifest, path in zip(loaded, manifests, strict=True)
        },
    }
    atomic_write_json(output, record)


def _parse_manifest(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        digest, separator, relative = line.partition("  ")
        if not separator or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise CharacterizationFailure(f"invalid retained manifest line: {line!r}")
        records[relative.removeprefix("./")] = digest
    return records


def export_linux_reference(source_root: Path, output_dir: Path) -> None:
    manifest_path = source_root / "FINAL_EVIDENCE_MANIFEST.sha256"
    if not manifest_path.is_file():
        raise CharacterizationFailure("retained PR 5 Linux evidence is unavailable")
    observed_manifest = sha256_path(manifest_path)
    if observed_manifest != PR5_EVIDENCE_MANIFEST_SHA256:
        raise CharacterizationFailure(
            "retained PR 5 Linux evidence manifest differs from accepted evidence"
        )
    manifest = _parse_manifest(manifest_path)
    records = []
    for representative in REPRESENTATIVES:
        source = source_root / representative.linux_relative
        expected = manifest.get(representative.linux_relative)
        if expected is None or not source.is_file() or sha256_path(source) != expected:
            raise CharacterizationFailure(
                f"retained Linux output is unusable: {representative.linux_relative}"
            )
        destination_relative = PurePosixPath(
            representative.program, representative.case_id, "ncbi.raw.out"
        )
        destination = output_dir.joinpath(*destination_relative.parts)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        if sha256_path(destination) != expected:
            raise CharacterizationFailure("copy changed retained Linux evidence")
        records.append(
            {
                "program": representative.program,
                "case_id": representative.case_id,
                "semantic_class": representative.semantic_class,
                "source_relative": representative.linux_relative,
                "artifact_relative": destination_relative.as_posix(),
                "sha256": expected,
                "byte_length": destination.stat().st_size,
                "newline_profile": newline_profile(destination.read_bytes()),
            }
        )
    atomic_write_json(
        output_dir / "linux_reference_manifest.json",
        {
            "evidence_schema": EVIDENCE_SCHEMA,
            "status": "RETAINED_PR5_LINUX_EVIDENCE_VERIFIED",
            "created_at": utc_now(),
            "source_manifest_sha256": observed_manifest,
            "case_count": len(records),
            "cases": records,
        },
    )


def _validate_fixture_bundle(
    authority: ModuleType,
    candidate_sha: str,
    platform_id: str,
    bundle_root: Path,
    invariance_path: Path,
    evidence_dir: Path,
) -> dict[str, dict[str, object]]:
    manifest = json.loads((bundle_root / "fixture_manifest.json").read_text("utf-8"))
    invariance = json.loads(invariance_path.read_text("utf-8"))
    if (
        manifest["candidate_sha"] != candidate_sha
        or manifest["platform_id"] != platform_id
    ):
        raise CharacterizationFailure("fixture artifact identity mismatch")
    if invariance["status"] != "CROSS_PLATFORM_FIXTURES_BYTE_IDENTICAL":
        raise CharacterizationFailure("cross-platform fixture gate is not verified")
    if set(invariance["platform_ids"]) != {"windows-x64", "macos-arm64", "macos-x64"}:
        raise CharacterizationFailure("cross-platform fixture gate topology changed")
    common = {item["repository_relative"]: item for item in invariance["inputs"]}
    identities: dict[str, dict[str, object]] = {}
    staged_records = []
    for record in manifest["inputs"]:
        relative = record["repository_relative"]
        if _invariant_input(record) != common.get(relative):
            raise CharacterizationFailure(f"fixture gate mismatch: {relative}")
        artifact_source = bundle_root.joinpath(
            *PurePosixPath(record["staged_relative"]).parts
        )
        data = artifact_source.read_bytes()
        observed_profile = require_lf_fixture(data, relative)
        if (
            sha256_bytes(data) != record["staged_sha256"]
            or len(data) != record["byte_length"]
            or observed_profile != record["newline_profile"]
            or record["lf_only"] is not True
        ):
            raise CharacterizationFailure(f"fixture artifact changed: {relative}")
        repo_relative = PurePosixPath(relative)
        lexical = _historical_lexical_path(authority, repo_relative)
        physical = authority.controlled_physical_path(lexical)
        atomic_write_bytes(physical, data)
        staged_data = physical.read_bytes()
        if staged_data != data:
            raise CharacterizationFailure(f"command-path staging changed: {relative}")
        identity = {
            **record,
            "lexical_command_path": lexical,
            "physical_command_path": str(physical.absolute()),
            "command_staged_sha256": sha256_bytes(staged_data),
        }
        identities[relative] = identity
        staged_records.append(identity)
    final_manifest = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "status": "COMMAND_FIXTURES_BYTE_IDENTICAL_AND_STAGED",
        "platform_id": platform_id,
        "candidate_sha": candidate_sha,
        "fixture_source_sha": PR6_CERT_SHA_V3,
        "fixture_count": len(staged_records),
        "cross_platform_gate_sha256": sha256_path(invariance_path),
        "inputs": staged_records,
    }
    atomic_write_json(evidence_dir / "fixture_manifest.json", final_manifest)
    return identities


def _validate_linux_bundle(
    bundle_root: Path,
) -> tuple[dict[tuple[str, str], Path], dict]:
    manifest_path = bundle_root / "linux_reference_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest["status"] != "RETAINED_PR5_LINUX_EVIDENCE_VERIFIED"
        or manifest["source_manifest_sha256"] != PR5_EVIDENCE_MANIFEST_SHA256
        or manifest["case_count"] != EXPECTED_REPRESENTATIVE_COUNT
    ):
        raise CharacterizationFailure(
            "retained Linux evidence artifact is not authoritative"
        )
    paths: dict[tuple[str, str], Path] = {}
    for record in manifest["cases"]:
        path = bundle_root.joinpath(*PurePosixPath(record["artifact_relative"]).parts)
        if not path.is_file() or sha256_path(path) != record["sha256"]:
            raise CharacterizationFailure(
                f"retained Linux evidence changed: {record['case_id']}"
            )
        paths[(record["program"], record["case_id"])] = path
    expected = {(item.program, item.case_id) for item in REPRESENTATIVES}
    if set(paths) != expected:
        raise CharacterizationFailure("retained Linux evidence case inventory changed")
    return paths, manifest


def _record_identity(
    repo_root: Path,
    candidate_sha: str,
    platform_id: str,
    authority: ModuleType,
    archive: Path,
    checksum: Path,
    oracles: dict[str, Path],
    fixture_manifest: Path,
    linux_manifest: dict,
    output_dir: Path,
) -> dict:
    spec = authority.PLATFORMS[platform_id]
    observed_os = platform.system()
    observed_machine = authority.normalize_machine(platform.machine())
    if observed_os != spec.os_name or observed_machine != spec.machine:
        raise CharacterizationFailure(
            f"runner mismatch: expected {spec.os_name}/{spec.machine}, "
            f"observed {observed_os}/{observed_machine}"
        )
    try:
        archive_record = authority.verify_archive(spec, archive, checksum)
    except authority.CertificationFailure as error:
        raise CharacterizationFailure(str(error)) from error
    oracle_records = {}
    for program, executable in oracles.items():
        if not executable.is_file():
            raise CharacterizationFailure(f"official {program} executable is absent")
        if authority.binary_architecture(executable) != spec.binary_arch:
            raise CharacterizationFailure(f"official {program} architecture mismatch")
        oracle_records[program] = {
            "path": str(executable.resolve()),
            "sha256": sha256_path(executable),
            "architecture": spec.binary_arch,
            "version_evidence": {
                "package": "NCBI BLAST+ 2.17.0",
                "source": "verified pinned official archive",
                "archive_name": spec.archive_name,
                "archive_published_md5": spec.archive_md5,
            },
        }
    identity = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "status": "PREEXECUTION_IDENTITY_VERIFIED",
        "created_at": utc_now(),
        "candidate_sha": candidate_sha,
        "diagnostic_base_sha": PR6_CERT_SHA_V3,
        "platform_id": platform_id,
        "runner_label": spec.runner_label,
        "runner": {
            "system": observed_os,
            "machine": observed_machine,
            "platform": platform.platform(),
            "release": platform.release(),
            "version": platform.version(),
            "RUNNER_OS": os.environ.get("RUNNER_OS", ""),
            "RUNNER_ARCH": os.environ.get("RUNNER_ARCH", ""),
            "ImageOS": os.environ.get("ImageOS", ""),
            "ImageVersion": os.environ.get("ImageVersion", ""),
        },
        "oracle_archive": {
            **archive_record,
            "published_md5": spec.archive_md5,
            "observed_archive_md5": archive_record["archive_md5"],
        },
        "oracles": oracle_records,
        "fixture_manifest_sha256": sha256_path(fixture_manifest),
        "retained_linux_evidence": {
            "source_manifest_sha256": linux_manifest["source_manifest_sha256"],
            "case_count": linux_manifest["case_count"],
        },
        "execution_contract": {
            "losat_executions": 0,
            "authoritative_ncbi_searches_per_platform": (
                EXPECTED_PLATFORM_AUTHORITATIVE_EXECUTIONS
            ),
            "diagnostic_ncbi_searches_per_platform": (
                EXPECTED_PLATFORM_DIAGNOSTIC_EXECUTIONS
            ),
            "ncbi_executions_per_platform": EXPECTED_PLATFORM_EXECUTIONS,
            "platform_count": 3,
            "total_authoritative_ncbi_searches": (
                EXPECTED_TOTAL_AUTHORITATIVE_EXECUTIONS
            ),
            "total_diagnostic_ncbi_searches": EXPECTED_TOTAL_DIAGNOSTIC_EXECUTIONS,
            "total_ncbi_executions": EXPECTED_TOTAL_EXECUTIONS,
            "repeatability_executions": 0,
            "ncbi_version_probe_executions": 0,
        },
    }
    atomic_write_json(output_dir / "identity.json", identity)
    return identity


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
# ```c++
# ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
#     if (iter != m_FieldsToShow.begin()) m_Ostream << m_FieldDelimiter;
#     x_PrintField(*iter);
# }
# m_Ostream << "\n";
# ```
def parse_structured_rows(data: bytes) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    for line in data.splitlines():
        if not line or line.startswith(b"#"):
            continue
        values = line.decode("utf-8").split("\t")
        if len(values) != len(STANDARD_FIELDS):
            raise CharacterizationFailure(
                f"expected {len(STANDARD_FIELDS)} tabular fields, observed {len(values)}"
            )
        row: dict[str, str | int] = {"row_index": len(rows) + 1}
        row.update(dict(zip(STANDARD_FIELDS, values, strict=True)))
        rows.append(row)
    return rows


def row_key(row: dict[str, str | int]) -> tuple[str, ...]:
    return tuple(str(row[field]) for field in KEY_FIELDS)


def structured_comparison(
    platform_rows: Sequence[dict[str, str | int]],
    linux_rows: Sequence[dict[str, str | int]],
    platform_raw: bytes,
    linux_raw: bytes,
) -> dict[str, object]:
    differences = []
    for index in range(max(len(platform_rows), len(linux_rows))):
        platform_row = platform_rows[index] if index < len(platform_rows) else None
        linux_row = linux_rows[index] if index < len(linux_rows) else None
        field_differences = []
        for field in STANDARD_FIELDS:
            platform_value = None if platform_row is None else platform_row[field]
            linux_value = None if linux_row is None else linux_row[field]
            if platform_value != linux_value:
                field_differences.append(
                    {"field": field, "linux": linux_value, "platform": platform_value}
                )
        if field_differences:
            differences.append(
                {
                    "row_index": index + 1,
                    "linux_row": linux_row,
                    "platform_row": platform_row,
                    "field_differences": field_differences,
                }
            )
    platform_keys = [row_key(row) for row in platform_rows]
    linux_keys = [row_key(row) for row in linux_rows]
    return {
        "linux_raw_sha256": sha256_bytes(linux_raw),
        "platform_raw_sha256": sha256_bytes(platform_raw),
        "platform_lf_diagnostic_sha256": sha256_bytes(
            lf_diagnostic_bytes(platform_raw)
        ),
        "raw_bytes_equal": platform_raw == linux_raw,
        "lf_diagnostic_bytes_equal": lf_diagnostic_bytes(platform_raw)
        == lf_diagnostic_bytes(linux_raw),
        "linux_row_count": len(linux_rows),
        "platform_row_count": len(platform_rows),
        "row_order_equal": list(platform_rows) == list(linux_rows),
        "row_key_order_equal": platform_keys == linux_keys,
        "row_key_multiset_equal": Counter(platform_keys) == Counter(linux_keys),
        "differing_row_count": len(differences),
        "differing_rows": differences,
    }


def load_platform_authoritative_raw(path: Path) -> bytes:
    """Admit only the role-specific authoritative raw artifact.

    NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:
    1098-1108 writes the selected fields directly to the output stream.  This
    gate ensures diagnostic fields never enter the authoritative comparison.
    """

    if path.name != "ncbi.raw.out" or path.parent.name != "authoritative":
        raise CharacterizationFailure(
            f"authoritative comparison rejected non-authoritative path: {path}"
        )
    return path.read_bytes()


def parse_diagnostic_rows(data: bytes) -> list[dict[str, str | int]]:
    """Parse ``6 std score btop qseq sseq`` in exact NCBI field order.

    NCBI references:
    - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:59-91 expands std
      and appends score, btop, qseq, and sseq in requested order.
    - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:973-1027 derives
      qseq, sseq, and BTOP from the same alignment representation.
    """

    rows: list[dict[str, str | int]] = []
    for line in data.splitlines():
        if not line:
            continue
        if line.startswith(b"#"):
            raise CharacterizationFailure("diagnostic outfmt 6 emitted comments")
        values = line.decode("utf-8").split("\t")
        if len(values) != len(DIAGNOSTIC_FIELDS):
            raise CharacterizationFailure(
                f"expected {len(DIAGNOSTIC_FIELDS)} diagnostic fields, "
                f"observed {len(values)}"
            )
        row: dict[str, str | int] = {"row_index": len(rows) + 1}
        row.update(dict(zip(DIAGNOSTIC_FIELDS, values, strict=True)))
        rows.append(row)
    return rows


def standard_projection(
    rows: Iterable[dict[str, str | int]],
) -> list[dict[str, str | int]]:
    return [
        {"row_index": index, **{field: row[field] for field in STANDARD_FIELDS}}
        for index, row in enumerate(rows, start=1)
    ]


def alignment_representation_identity(
    program: str, case_id: str, row: dict[str, str | int]
) -> tuple[str, ...]:
    return (
        program,
        case_id,
        *(str(row[field]) for field in KEY_FIELDS),
        str(row["score"]),
    )


def rich_representation_finding(
    platform_representations: Sequence[dict[str, str | int]],
    expected_platforms: set[str],
) -> tuple[str, int, set[str]]:
    """Use BTOP/qseq/sseq evidence for the narrow diagnostic finding only.

    NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:
    973-1027 obtains qseq/sseq and constructs BTOP from the alignment.  A claim
    of alternate representation therefore requires distinct values in those
    fields for rows already grouped by equal raw score and endpoints.
    """

    distinct_representations = {
        (
            str(row["btop"]),
            str(row["qseq"]),
            str(row["sseq"]),
        )
        for row in platform_representations
    }
    matched_platforms = {
        str(row["platform_id"]) for row in platform_representations
    }
    if len(distinct_representations) > 1:
        finding = (
            "EVIDENCE_CONSISTENT_WITH_ALTERNATIVE_EQUAL_SCORING_"
            "ALIGNMENT_REPRESENTATIONS"
        )
    elif matched_platforms == expected_platforms:
        finding = (
            "SAME_SCORE_AND_ENDPOINTS_HAVE_ONE_OBSERVED_RICH_"
            "REPRESENTATION_ACROSS_NEW_PLATFORMS"
        )
    else:
        finding = "RICH_REPRESENTATION_COMPARISON_INCOMPLETE"
    return finding, len(distinct_representations), matched_platforms


def _case_authority_context(representative: Representative) -> dict[str, object]:
    if representative.case_id == "Sakai.MG1655.megablast":
        return {
            "existing_authority": "SOURCE_UNDETERMINED_ACCEPTED",
            "separate_from_platform_variance": True,
            "not_reinterpreted_as": ["exact", "intentional_LOSAT_divergence"],
        }
    if representative.case_id == "d06_ap027131_ap027133_db4":
        return {
            "existing_authority": "approved_db_gencode_deviation",
            "separate_from_platform_variance": True,
            "does_not_create_additional_tblastx_deviation": True,
        }
    return {"existing_authority_unchanged": True}


def _pair_identity_record(
    pair: dict[str, object], oracle_sha256: str
) -> dict[str, object]:
    authoritative = pair["authoritative"]
    diagnostic = pair["diagnostic"]
    comparison = compare_pair_commands(
        authoritative["argv"], diagnostic["argv"]
    )
    program = str(pair["program"])
    authoritative_inputs = {
        role: (option, value)
        for role, option, value in input_arguments(authoritative["argv"], program)
    }
    diagnostic_inputs = {
        role: (option, value)
        for role, option, value in input_arguments(diagnostic["argv"], program)
    }
    input_records = {}
    for role, details in pair["inputs"].items():
        authoritative_input = authoritative_inputs.get(role)
        diagnostic_input = diagnostic_inputs.get(role)
        if authoritative_input != diagnostic_input:
            raise CharacterizationFailure(
                f"paired {role} command paths differ for {pair['case_id']}"
            )
        if authoritative_input != (details["option"], details["lexical_path"]):
            raise CharacterizationFailure(
                f"paired {role} command escaped staged fixture for {pair['case_id']}"
            )
        input_records[role] = {
            "repository_relative": details["repository_relative"],
            "authoritative_command_path": authoritative_input[1],
            "diagnostic_command_path": diagnostic_input[1],
            "authoritative_sha256": details["sha256"],
            "diagnostic_sha256": details["sha256"],
            "identical_sha256": True,
            "byte_length": details["byte_length"],
            "newline_profile": details["newline_profile"],
            "lf_only": details["lf_only"],
        }
    if set(input_records) != {"query", "subject"}:
        raise CharacterizationFailure("pair input identity is incomplete")
    record = {
        "evidence_schema": EVIDENCE_SCHEMA,
        "status": "PREEXECUTION_PAIR_IDENTITY_VERIFIED",
        "pair_index": pair["pair_index"],
        "program": program,
        "case_id": pair["case_id"],
        **comparison,
        "same_executable_path": (
            authoritative["argv"][0] == diagnostic["argv"][0]
        ),
        "authoritative_executable_sha256": oracle_sha256,
        "diagnostic_executable_sha256": oracle_sha256,
        "identical_executable_sha256": True,
        "identical_input_sha256": all(
            record["identical_sha256"] for record in input_records.values()
        ),
        "inputs": input_records,
        "authoritative_output_is_canonical_for_this_characterization": True,
        "diagnostic_output_is_authoritative": False,
        "diagnostic_output_may_satisfy_authoritative_comparison": False,
    }
    if not (
        record["identical_non_output_arguments"]
        and record["same_executable_path"]
        and record["identical_executable_sha256"]
        and record["identical_input_sha256"]
        and all(item["lf_only"] for item in input_records.values())
    ):
        raise CharacterizationFailure(
            f"pair identity ratchet failed for {pair['case_id']}"
        )
    return record


def _command_record(
    repo_root: Path,
    pair: dict[str, object],
    role: str,
) -> dict[str, object]:
    search = pair[role]
    representative = next(
        item
        for item in REPRESENTATIVES
        if item.program == pair["program"] and item.case_id == pair["case_id"]
    )
    return {
        "evidence_schema": EVIDENCE_SCHEMA,
        "execution_index": search["execution_index"],
        "pair_index": pair["pair_index"],
        "search_role": role,
        "program": pair["program"],
        "case_id": pair["case_id"],
        "semantic_class": pair["semantic_class"],
        "argv": search["argv"],
        "cwd": str(repo_root),
        "shell": False,
        "environment_overrides": {},
        "option_pairs": search["option_pairs"],
        "semantics": search["semantics"],
        "inputs": pair["inputs"],
        "authority_context": _case_authority_context(representative),
        "authoritative_comparison_eligible": role == "authoritative",
        "canonical_baseline_eligible": role == "authoritative",
        "independent_biological_replicate": False,
    }


def _execute_paired_search(
    repo_root: Path,
    pair: dict[str, object],
    role: str,
) -> tuple[bytes, dict[str, object]]:
    """Execute one authorized search and retain its output without rewriting.

    NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:
    1098-1108 writes fields and the row terminator directly to the configured
    stream.  ``os.replace`` retains those exact bytes under the role directory.
    """

    pair_dir = Path(str(_option_value(pair[role]["argv"], "-out"))).parent
    role_dir = pair_dir / role
    role_dir.mkdir(parents=True, exist_ok=True)
    transient_raw = Path(str(_option_value(pair[role]["argv"], "-out")))
    retained_raw = role_dir / "ncbi.raw.out"
    if transient_raw.exists() or retained_raw.exists():
        raise CharacterizationFailure(
            f"refusing to overwrite {role} raw output: {pair['case_id']}"
        )
    atomic_write_json(
        role_dir / "command.json", _command_record(repo_root, pair, role)
    )
    started_at = utc_now()
    completed = run_bytes(pair[role]["argv"], repo_root)
    completed_at = utc_now()
    atomic_write_bytes(role_dir / "stdout.txt", completed.stdout)
    atomic_write_bytes(role_dir / "stderr.txt", completed.stderr)
    if completed.returncode != 0 or not transient_raw.is_file():
        raise CharacterizationFailure(
            f"official NCBI {role} search failed for {pair['case_id']}: "
            f"exit {completed.returncode}"
        )
    os.replace(transient_raw, retained_raw)
    raw = retained_raw.read_bytes()
    return raw, {
        "returncode": completed.returncode,
        "started_at": started_at,
        "completed_at": completed_at,
        "capture_method": "byte-preserving os.replace from shared pair sink",
        "retained_relative": f"{role}/ncbi.raw.out",
    }


def changed_hsp_diagnostic_records(
    comparison: dict[str, object],
    diagnostic_rows: Sequence[dict[str, str | int]],
) -> list[dict[str, object]]:
    records = []
    for difference in comparison["differing_rows"]:
        row_index = int(difference["row_index"])
        diagnostic_row = (
            diagnostic_rows[row_index - 1]
            if row_index <= len(diagnostic_rows)
            else None
        )
        records.append(
            {
                "row_index": row_index,
                "linux_retained_authoritative": difference["linux_row"],
                "platform_authoritative": difference["platform_row"],
                "platform_diagnostic": diagnostic_row,
                "authoritative_field_differences": difference[
                    "field_differences"
                ],
                "diagnostic_used_for_authoritative_pass_fail": False,
            }
        )
    return records


def _prepare_characterization(
    args: argparse.Namespace,
) -> tuple[Path, Path, list[dict[str, object]], dict, dict[tuple[str, str], Path]]:
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    verify_candidate(repo_root, args.candidate_sha)
    authority, catalog = load_platform_authority(repo_root)
    if args.platform_id not in authority.PLATFORMS:
        raise CharacterizationFailure(f"unsupported platform: {args.platform_id}")
    fixture_identities = _validate_fixture_bundle(
        authority,
        args.candidate_sha,
        args.platform_id,
        args.fixture_bundle.resolve(),
        args.invariance_manifest.resolve(),
        output_dir,
    )
    linux_paths, linux_manifest = _validate_linux_bundle(args.linux_bundle.resolve())
    oracles = {
        "blastn": args.blastn_oracle.resolve(),
        "blastp": args.blastp_oracle.resolve(),
        "tblastx": args.tblastx_oracle.resolve(),
    }
    identity = _record_identity(
        repo_root,
        args.candidate_sha,
        args.platform_id,
        authority,
        args.oracle_archive.resolve(),
        args.oracle_checksum.resolve(),
        oracles,
        output_dir / "fixture_manifest.json",
        linux_manifest,
        output_dir,
    )
    plan = build_command_plan(repo_root, authority, catalog, oracles, output_dir)
    for pair in plan:
        for details in pair["inputs"].values():
            relative = details["repository_relative"]
            details.update(
                {
                    "sha256": fixture_identities[relative]["command_staged_sha256"],
                    "byte_length": fixture_identities[relative]["byte_length"],
                    "newline_profile": fixture_identities[relative]["newline_profile"],
                    "lf_only": fixture_identities[relative]["lf_only"],
                }
            )
    atomic_write_json(
        output_dir / "command_plan.json",
        {
            "evidence_schema": EVIDENCE_SCHEMA,
            "status": "PREEXECUTION_PLAN_VERIFIED",
            "representative_count": len(plan),
            "authoritative_count": EXPECTED_PLATFORM_AUTHORITATIVE_EXECUTIONS,
            "diagnostic_count": EXPECTED_PLATFORM_DIAGNOSTIC_EXECUTIONS,
            "execution_count": EXPECTED_PLATFORM_EXECUTIONS,
            "losat_executions": 0,
            "pairs": plan,
        },
    )
    return repo_root, output_dir, plan, identity, linux_paths


def analyze_platforms(
    platform_roots: Sequence[Path], linux_bundle: Path, output: Path
) -> None:
    """Correlate rich diagnostic rows without changing authoritative results.

    NCBI references:
    - ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-41,
      96-101, 114-119, 144-146 defines the standard and added rich fields.
    - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:895-1027 derives
      coordinates, identity, aligned sequences, and BTOP from each alignment.
    """

    if len(platform_roots) != 3:
        raise CharacterizationFailure("rich analysis requires exactly three platforms")
    roots: dict[str, Path] = {}
    summaries: dict[str, dict[str, object]] = {}
    for root in platform_roots:
        identity = json.loads((root / "identity.json").read_text(encoding="utf-8"))
        summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
        platform_id = identity["platform_id"]
        if (
            identity["status"] != "CHARACTERIZATION_COMPLETE"
            or summary["platform_id"] != platform_id
            or summary["authoritative_search_executions"]
            != EXPECTED_PLATFORM_AUTHORITATIVE_EXECUTIONS
            or summary["diagnostic_search_executions"]
            != EXPECTED_PLATFORM_DIAGNOSTIC_EXECUTIONS
            or summary["ncbi_search_executions"] != EXPECTED_PLATFORM_EXECUTIONS
            or summary["losat_executions"] != 0
        ):
            raise CharacterizationFailure(
                f"platform evidence count/identity failed: {platform_id}"
            )
        if platform_id in roots:
            raise CharacterizationFailure(f"duplicate platform evidence: {platform_id}")
        roots[platform_id] = root
        summaries[platform_id] = summary
    expected_platforms = {"windows-x64", "macos-arm64", "macos-x64"}
    if set(roots) != expected_platforms:
        raise CharacterizationFailure(
            f"rich analysis platform topology changed: {set(roots)}"
        )
    linux_paths, linux_manifest = _validate_linux_bundle(linux_bundle)

    case_analyses = []
    total_changed_authoritative_hsps = 0
    total_unmatched_changed_authoritative_hsps = 0
    alternative_representation_groups = 0
    for representative in REPRESENTATIVES:
        linux_raw = linux_paths[
            (representative.program, representative.case_id)
        ].read_bytes()
        linux_rows = parse_structured_rows(linux_raw)
        diagnostics_by_platform: dict[
            str, list[dict[str, str | int]]
        ] = {}
        triggers_by_identity: dict[tuple[str, ...], list[dict[str, object]]] = {}
        unmatched = []
        platform_authoritative_results = []
        for platform_id in sorted(roots):
            case_dir = (
                roots[platform_id]
                / representative.program
                / representative.case_id
            )
            authoritative_path = case_dir / "authoritative" / "ncbi.raw.out"
            diagnostic_path = case_dir / "diagnostic" / "ncbi.raw.out"
            pair_identity = json.loads(
                (case_dir / "pair_identity.json").read_text(encoding="utf-8")
            )
            authoritative_command = json.loads(
                (case_dir / "authoritative" / "command.json").read_text(
                    encoding="utf-8"
                )
            )
            diagnostic_command = json.loads(
                (case_dir / "diagnostic" / "command.json").read_text(
                    encoding="utf-8"
                )
            )
            observed_pair_comparison = compare_pair_commands(
                authoritative_command["argv"], diagnostic_command["argv"]
            )
            if (
                pair_identity["status"] != "PAIR_EXECUTIONS_COMPLETE"
                or not pair_identity["identical_non_output_arguments"]
                or not pair_identity["identical_executable_sha256"]
                or not pair_identity["identical_input_sha256"]
                or pair_identity["diagnostic_outfmt"] != DIAGNOSTIC_OUTFMT
                or not observed_pair_comparison["identical_non_output_arguments"]
                or authoritative_command["authoritative_comparison_eligible"]
                is not True
                or diagnostic_command["authoritative_comparison_eligible"]
                is not False
            ):
                raise CharacterizationFailure(
                    f"pair identity artifact failed for "
                    f"{platform_id}/{representative.case_id}"
                )
            authoritative_raw = load_platform_authoritative_raw(authoritative_path)
            diagnostic_raw = diagnostic_path.read_bytes()
            if (
                sha256_bytes(authoritative_raw)
                != pair_identity["authoritative_raw_sha256"]
                or sha256_bytes(diagnostic_raw)
                != pair_identity["diagnostic_raw_sha256"]
            ):
                raise CharacterizationFailure(
                    f"paired raw output hash failed for "
                    f"{platform_id}/{representative.case_id}"
                )
            authoritative_rows = parse_structured_rows(authoritative_raw)
            diagnostic_rows = parse_diagnostic_rows(diagnostic_raw)
            if standard_projection(diagnostic_rows) != authoritative_rows:
                raise CharacterizationFailure(
                    f"diagnostic projection changed authoritative rows for "
                    f"{platform_id}/{representative.case_id}"
                )
            diagnostics_by_platform[platform_id] = diagnostic_rows
            comparison = structured_comparison(
                authoritative_rows, linux_rows, authoritative_raw, linux_raw
            )
            changed = changed_hsp_diagnostic_records(comparison, diagnostic_rows)
            total_changed_authoritative_hsps += len(changed)
            platform_authoritative_results.append(
                {
                    "platform_id": platform_id,
                    "authoritative_raw_sha256": sha256_bytes(authoritative_raw),
                    "diagnostic_raw_sha256": sha256_bytes(diagnostic_raw),
                    "authoritative_row_count": len(authoritative_rows),
                    "diagnostic_row_count": len(diagnostic_rows),
                    "changed_authoritative_hsp_count": len(changed),
                    "raw_authoritative_bytes_equal_to_linux": comparison[
                        "raw_bytes_equal"
                    ],
                }
            )
            for changed_record in changed:
                diagnostic_row = changed_record["platform_diagnostic"]
                authoritative_row = changed_record["platform_authoritative"]
                trigger = {
                    "platform_id": platform_id,
                    **changed_record,
                }
                if diagnostic_row is None or authoritative_row is None:
                    unmatched.append(trigger)
                    total_unmatched_changed_authoritative_hsps += 1
                    continue
                identity_key = alignment_representation_identity(
                    representative.program,
                    representative.case_id,
                    diagnostic_row,
                )
                triggers_by_identity.setdefault(identity_key, []).append(trigger)

        rich_groups = []
        for identity_key, triggers in sorted(triggers_by_identity.items()):
            platform_representations = []
            for platform_id, rows in sorted(diagnostics_by_platform.items()):
                for row in rows:
                    if (
                        alignment_representation_identity(
                            representative.program, representative.case_id, row
                        )
                        == identity_key
                    ):
                        platform_representations.append(
                            {
                                "platform_id": platform_id,
                                "row_index": row["row_index"],
                                **{field: row[field] for field in DIAGNOSTIC_FIELDS},
                            }
                        )
            finding, distinct_count, matched_platforms = (
                rich_representation_finding(
                    platform_representations, expected_platforms
                )
            )
            if finding.startswith("EVIDENCE_CONSISTENT_WITH_ALTERNATIVE"):
                alternative_representation_groups += 1
            rich_groups.append(
                {
                    "alignment_identity": {
                        "program": identity_key[0],
                        "case_id": identity_key[1],
                        "query_id": identity_key[2],
                        "subject_id": identity_key[3],
                        "qstart": identity_key[4],
                        "qend": identity_key[5],
                        "sstart": identity_key[6],
                        "send": identity_key[7],
                        "raw_score": identity_key[8],
                    },
                    "triggering_authoritative_differences": triggers,
                    "platform_representations": platform_representations,
                    "matched_new_platforms": sorted(matched_platforms),
                    "distinct_btop_qseq_sseq_count": distinct_count,
                    "finding": finding,
                    "diagnostic_used_for_authoritative_pass_fail": False,
                }
            )
        case_analyses.append(
            {
                "program": representative.program,
                "case_id": representative.case_id,
                "semantic_class": representative.semantic_class,
                "authority_context": _case_authority_context(representative),
                "linux_retained": {
                    "source": "retained PR 5 authoritative output",
                    "raw_sha256": sha256_bytes(linux_raw),
                    "row_count": len(linux_rows),
                    "rich_fields_available": False,
                    "unavailable_fields": ["score", "btop", "qseq", "sseq"],
                    "rerun_performed": False,
                },
                "platform_authoritative_results": platform_authoritative_results,
                "rich_diagnostic_groups": rich_groups,
                "unmatched_changed_authoritative_hsps": unmatched,
            }
        )

    authoritative_count = sum(
        int(summary["authoritative_search_executions"])
        for summary in summaries.values()
    )
    diagnostic_count = sum(
        int(summary["diagnostic_search_executions"])
        for summary in summaries.values()
    )
    total_count = sum(
        int(summary["ncbi_search_executions"]) for summary in summaries.values()
    )
    if (
        authoritative_count != EXPECTED_TOTAL_AUTHORITATIVE_EXECUTIONS
        or diagnostic_count != EXPECTED_TOTAL_DIAGNOSTIC_EXECUTIONS
        or total_count != EXPECTED_TOTAL_EXECUTIONS
    ):
        raise CharacterizationFailure("three-platform 18+18 count ratchet failed")
    atomic_write_json(
        output,
        {
            "evidence_schema": EVIDENCE_SCHEMA,
            "status": "RICH_DIAGNOSTIC_ANALYSIS_COMPLETE",
            "created_at": utc_now(),
            "platform_ids": sorted(roots),
            "retained_linux_manifest_sha256": linux_manifest[
                "source_manifest_sha256"
            ],
            "authoritative_search_count": authoritative_count,
            "diagnostic_search_count": diagnostic_count,
            "total_ncbi_search_count": total_count,
            "losat_execution_count": 0,
            "diagnostic_outfmt": DIAGNOSTIC_OUTFMT,
            "authoritative_pass_fail_source": "authoritative/ncbi.raw.out only",
            "diagnostic_output_is_canonical_baseline": False,
            "diagnostic_output_used_for_authoritative_pass_fail": False,
            "authority_modified": False,
            "total_changed_authoritative_hsps": total_changed_authoritative_hsps,
            "total_unmatched_changed_authoritative_hsps": (
                total_unmatched_changed_authoritative_hsps
            ),
            "alternative_equal_scoring_representation_groups": (
                alternative_representation_groups
            ),
            "cases": case_analyses,
        },
    )

def characterize(args: argparse.Namespace) -> None:
    repo_root, output_dir, plan, identity, linux_paths = _prepare_characterization(
        args
    )
    summaries = []
    for pair in plan:
        representative = next(
            item
            for item in REPRESENTATIVES
            if item.program == pair["program"] and item.case_id == pair["case_id"]
        )
        case_dir = output_dir / representative.program / representative.case_id
        case_dir.mkdir(parents=True, exist_ok=True)
        pair_identity = _pair_identity_record(
            pair, identity["oracles"][representative.program]["sha256"]
        )
        atomic_write_json(case_dir / "pair_identity.json", pair_identity)

        # The authoritative search is deliberately first in every pair.
        platform_raw, authoritative_execution = _execute_paired_search(
            repo_root, pair, "authoritative"
        )
        diagnostic_raw, diagnostic_execution = _execute_paired_search(
            repo_root, pair, "diagnostic"
        )
        authoritative_dir = case_dir / "authoritative"
        diagnostic_dir = case_dir / "diagnostic"
        retained_platform_raw = load_platform_authoritative_raw(
            authoritative_dir / "ncbi.raw.out"
        )
        if retained_platform_raw != platform_raw:
            raise CharacterizationFailure("authoritative raw bytes changed after capture")
        linux_raw = linux_paths[
            (representative.program, representative.case_id)
        ].read_bytes()
        platform_rows = parse_structured_rows(platform_raw)
        linux_rows = parse_structured_rows(linux_raw)
        diagnostic_rows = parse_diagnostic_rows(diagnostic_raw)
        if standard_projection(diagnostic_rows) != platform_rows:
            raise CharacterizationFailure(
                f"diagnostic std projection differs from authoritative rows: "
                f"{representative.case_id}"
            )
        raw_record = {
            "evidence_schema": EVIDENCE_SCHEMA,
            "raw_sha256": sha256_bytes(platform_raw),
            "byte_length": len(platform_raw),
            **newline_profile(platform_raw),
            "lf_diagnostic_sha256": sha256_bytes(lf_diagnostic_bytes(platform_raw)),
            "lf_diagnostic_materialized": False,
            **authoritative_execution,
            "source_role": "authoritative",
            "raw_bytes_rewritten": False,
        }
        atomic_write_json(authoritative_dir / "raw_profile.json", raw_record)
        atomic_write_json(
            authoritative_dir / "structured_rows.json",
            {
                "source": "authoritative/ncbi.raw.out",
                "source_sha256": raw_record["raw_sha256"],
                "fields": list(STANDARD_FIELDS),
                "row_count": len(platform_rows),
                "rows": platform_rows,
                "authoritative_comparison_eligible": True,
            },
        )
        comparison = structured_comparison(
            platform_rows, linux_rows, platform_raw, linux_raw
        )
        comparison["authority_context"] = _case_authority_context(representative)
        comparison["platform_source"] = "authoritative/ncbi.raw.out"
        comparison["diagnostic_output_used"] = False
        atomic_write_json(
            authoritative_dir / "linux_comparison.json", comparison
        )
        diagnostic_record = {
            "evidence_schema": EVIDENCE_SCHEMA,
            "raw_sha256": sha256_bytes(diagnostic_raw),
            "byte_length": len(diagnostic_raw),
            **newline_profile(diagnostic_raw),
            **diagnostic_execution,
            "source_role": "diagnostic",
            "outfmt": DIAGNOSTIC_OUTFMT,
            "authoritative_comparison_eligible": False,
            "canonical_baseline_eligible": False,
            "independent_biological_replicate": False,
        }
        atomic_write_json(diagnostic_dir / "raw_profile.json", diagnostic_record)
        atomic_write_json(
            diagnostic_dir / "alignment_rows.json",
            {
                "source": "diagnostic/ncbi.raw.out",
                "source_sha256": diagnostic_record["raw_sha256"],
                "outfmt": DIAGNOSTIC_OUTFMT,
                "fields": list(DIAGNOSTIC_FIELDS),
                "row_count": len(diagnostic_rows),
                "rows": diagnostic_rows,
                "standard_projection_matches_authoritative": True,
                "authoritative_comparison_eligible": False,
            },
        )
        changed_hsps = changed_hsp_diagnostic_records(comparison, diagnostic_rows)
        atomic_write_json(
            diagnostic_dir / "changed_authoritative_hsps.json",
            {
                "source_authoritative_comparison": (
                    "authoritative/linux_comparison.json"
                ),
                "diagnostic_source": "diagnostic/alignment_rows.json",
                "changed_authoritative_hsp_count": len(changed_hsps),
                "diagnostic_used_for_authoritative_pass_fail": False,
                "rows": changed_hsps,
            },
        )
        pair_identity["status"] = "PAIR_EXECUTIONS_COMPLETE"
        pair_identity["authoritative_raw_sha256"] = raw_record["raw_sha256"]
        pair_identity["diagnostic_raw_sha256"] = diagnostic_record["raw_sha256"]
        pair_identity["authoritative_executed_first"] = True
        atomic_write_json(case_dir / "pair_identity.json", pair_identity)
        summaries.append(
            {
                "program": representative.program,
                "case_id": representative.case_id,
                "semantic_class": representative.semantic_class,
                "authoritative_raw_sha256": raw_record["raw_sha256"],
                "authoritative_lf_diagnostic_sha256": raw_record[
                    "lf_diagnostic_sha256"
                ],
                "authoritative_row_count": len(platform_rows),
                "diagnostic_raw_sha256": diagnostic_record["raw_sha256"],
                "diagnostic_row_count": len(diagnostic_rows),
                "diagnostic_outfmt": DIAGNOSTIC_OUTFMT,
                "linux_raw_sha256": comparison["linux_raw_sha256"],
                "linux_row_count": comparison["linux_row_count"],
                "differing_row_count": comparison["differing_row_count"],
                "row_key_order_equal": comparison["row_key_order_equal"],
                "authority_context": _case_authority_context(representative),
            }
        )
    if len(summaries) != EXPECTED_REPRESENTATIVE_COUNT:
        raise CharacterizationFailure(
            "platform did not complete exactly six paired representatives"
        )
    identity["status"] = "CHARACTERIZATION_COMPLETE"
    atomic_write_json(output_dir / "identity.json", identity)
    atomic_write_json(
        output_dir / "summary.json",
        {
            "evidence_schema": EVIDENCE_SCHEMA,
            "decision": "DIAGNOSTIC_CHARACTERIZATION_COMPLETE",
            "completed_at": utc_now(),
            "candidate_sha": args.candidate_sha,
            "platform_id": args.platform_id,
            "representative_count": len(summaries),
            "authoritative_search_executions": len(summaries),
            "diagnostic_search_executions": len(summaries),
            "ncbi_search_executions": len(summaries) * 2,
            "diagnostic_outfmt": DIAGNOSTIC_OUTFMT,
            "losat_executions": 0,
            "repeatability_executions": 0,
            "authority_modified": False,
            "identity_sha256": sha256_path(output_dir / "identity.json"),
            "cases": summaries,
            "excluded": [
                "LOSAT execution",
                "43-case matrix",
                "repeatability",
                "performance",
                "Wasm",
                "Linux rerun",
                "authority decision",
            ],
        },
    )
    atomic_write_bytes(
        output_dir / "workflow_state.txt",
        (
            "status=CHARACTERIZATION_COMPLETE\n"
            f"candidate_sha={args.candidate_sha}\n"
            f"platform_id={args.platform_id}\n"
            "losat_executions=0\n"
            f"authoritative_ncbi_searches={len(summaries)}\n"
            f"diagnostic_ncbi_searches={len(summaries)}\n"
            f"total_ncbi_searches={len(summaries) * 2}\n"
        ).encode("utf-8"),
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    default_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)

    stage = subparsers.add_parser("stage-fixtures")
    stage.add_argument("--repo-root", type=Path, default=default_root)
    stage.add_argument("--candidate-sha", required=True)
    stage.add_argument("--platform-id", required=True)
    stage.add_argument("--output-dir", type=Path, required=True)

    invariant = subparsers.add_parser("verify-fixture-invariance")
    invariant.add_argument("--manifest", type=Path, action="append", required=True)
    invariant.add_argument("--output", type=Path, required=True)

    linux = subparsers.add_parser("export-linux-reference")
    linux.add_argument("--source-root", type=Path, default=PR5_EVIDENCE_ROOT)
    linux.add_argument("--output-dir", type=Path, required=True)

    run = subparsers.add_parser("characterize")
    run.add_argument("--repo-root", type=Path, default=default_root)
    run.add_argument("--candidate-sha", required=True)
    run.add_argument("--platform-id", required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--fixture-bundle", type=Path, required=True)
    run.add_argument("--invariance-manifest", type=Path, required=True)
    run.add_argument("--linux-bundle", type=Path, required=True)
    run.add_argument("--oracle-archive", type=Path, required=True)
    run.add_argument("--oracle-checksum", type=Path, required=True)
    run.add_argument("--blastn-oracle", type=Path, required=True)
    run.add_argument("--blastp-oracle", type=Path, required=True)
    run.add_argument("--tblastx-oracle", type=Path, required=True)

    analyze = subparsers.add_parser("analyze-platforms")
    analyze.add_argument(
        "--platform-evidence", type=Path, action="append", required=True
    )
    analyze.add_argument("--linux-bundle", type=Path, required=True)
    analyze.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = getattr(args, "output_dir", None)
    if output_dir is None and hasattr(args, "output"):
        output_dir = args.output.parent
    try:
        if args.mode == "stage-fixtures":
            stage_fixtures(
                args.repo_root.resolve(),
                args.candidate_sha,
                args.platform_id,
                args.output_dir.resolve(),
            )
        elif args.mode == "verify-fixture-invariance":
            verify_fixture_invariance(args.manifest, args.output.resolve())
        elif args.mode == "export-linux-reference":
            export_linux_reference(
                args.source_root.resolve(), args.output_dir.resolve()
            )
        elif args.mode == "analyze-platforms":
            analyze_platforms(
                [path.resolve() for path in args.platform_evidence],
                args.linux_bundle.resolve(),
                args.output.resolve(),
            )
        else:
            characterize(args)
        return 0
    except (
        CharacterizationFailure,
        OSError,
        ValueError,
        KeyError,
        json.JSONDecodeError,
    ) as error:
        if output_dir is not None:
            resolved = output_dir.resolve()
            resolved.mkdir(parents=True, exist_ok=True)
            atomic_write_json(
                resolved / "failure.json",
                {
                    "status": "FAILED_CLOSED",
                    "failed_at": utc_now(),
                    "error": str(error),
                },
            )
        print(f"NCBI_PLATFORM_CHARACTERIZATION_FAILED: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
