#!/usr/bin/env python3
"""Assemble and verify exact-SHA LOSAT v0.2.0 release-candidate artifacts.

The v0.2.0 release candidate reuses no LOSAT-side result of v0.1.0.  The contract
(`docs/release/v0.2.0_rc_contract.json`, schema 2) names the candidate commit C that every
certification run used; the release-candidate commit R built here may differ from C only
outside the runtime inputs.  Archive assembly, extraction and the source-crate checks are the
v0.1.0 helpers (`prepare_release_candidate_v010.py`, unchanged); the contract, lineage,
version and smoke rules below are the v0.2.0 ones.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence


def _load_v010():
    path = Path(__file__).with_name("prepare_release_candidate_v010.py")
    spec = importlib.util.spec_from_file_location("_losat_rc_v010_helpers", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_v010 = _load_v010()

# v0.1.0 helpers that do not name the release.
ReleaseFailure = _v010.ReleaseFailure
OUTPUT_ENV_KEYS = _v010.OUTPUT_ENV_KEYS
sha256_path = _v010.sha256_path
json_bytes = _v010.json_bytes
write_json = _v010.write_json
run_capture = _v010.run_capture
require_success = _v010.require_success
require_sha = _v010.require_sha
package_version = _v010.package_version
capture_toolchain = _v010.capture_toolchain
binary_architecture = _v010.binary_architecture
expected_architecture = _v010.expected_architecture
tar_bytes = _v010.tar_bytes
zip_bytes = _v010.zip_bytes
safe_member_name = _v010.safe_member_name
extract_archive = _v010.extract_archive
archive_entries = _v010.archive_entries
write_sidecars = _v010.write_sidecars
crate_members = _v010.crate_members
read_sidecar = _v010.read_sidecar

CONTRACT_PATH = Path("docs/release/v0.2.0_rc_contract.json")
RELEASE = "v0.2.0"
PACKAGE_VERSION = "0.2.0"
SCHEMA_VERSION = 2
SIDECAR_SCHEMA_VERSION = 2
EXPECTED_ARTIFACTS = (
    ("native", "x86_64-unknown-linux-gnu"),
    ("native", "aarch64-unknown-linux-gnu"),
    ("native", "x86_64-pc-windows-msvc"),
    ("native", "aarch64-apple-darwin"),
    ("native", "x86_64-apple-darwin"),
    ("wasm", "wasm32-wasip1"),
    ("wasm", "wasm32-wasip1-threads"),
    ("source", "source"),
)
SMOKE_PROGRAMS = frozenset({"blastn", "blastp", "blastx", "tblastn", "tblastx"})
SMOKE_PLACEHOLDERS = ("{query}", "{subject}", "{output}")
WASM_KIND_BY_TARGET = {
    "wasm32-wasip1": "serial-command",
    "wasm32-wasip1-threads": "threaded-command",
}
EXPECTED_RUN_KEYS = ("B", "C", "F", "G", "H", "I")
GATED_RUN_KEYS = ("B", "C", "F", "H")
LOCAL_EVIDENCE_KEYS = ("A", "D", "E")
HEX = "0123456789abcdef"


def is_hex(value: object, length: int) -> bool:
    return (
        isinstance(value, str)
        and len(value) == length
        and all(character in HEX for character in value)
    )


def load_contract(repo_root: Path) -> dict[str, object]:
    path = repo_root / CONTRACT_PATH
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseFailure(f"cannot read RC contract {path}: {error}") from error
    if not isinstance(contract, dict):
        raise ReleaseFailure("RC contract is not an object")
    if contract.get("schema_version") != SCHEMA_VERSION:
        raise ReleaseFailure("unsupported RC contract schema")
    if contract.get("release") != RELEASE or contract.get("package_version") != PACKAGE_VERSION:
        raise ReleaseFailure("RC contract release/version changed")
    artifacts = contract.get("artifacts")
    if not isinstance(artifacts, list) or len(artifacts) != len(EXPECTED_ARTIFACTS):
        raise ReleaseFailure("RC contract must declare exactly eight artifacts")
    if not all(isinstance(item, dict) for item in artifacts):
        raise ReleaseFailure("RC contract artifacts must be objects")
    filenames = [item.get("filename") for item in artifacts]
    if (
        len(filenames) != len(set(filenames))
        or any(not isinstance(name, str) or not name for name in filenames)
    ):
        raise ReleaseFailure("RC artifact filenames must be unique and non-empty")
    observed = sorted((str(item.get("kind")), str(item.get("target"))) for item in artifacts)
    if observed != sorted(EXPECTED_ARTIFACTS):
        raise ReleaseFailure("RC contract artifacts must be one per release target")
    hash_policy = contract.get("binary_hash_policy")
    if (
        not isinstance(hash_policy, dict)
        or hash_policy.get("cross_runner_binary_equality_required") is not False
    ):
        raise ReleaseFailure("RC binary hash policy must reject cross-runner equality claims")
    for item in artifacts:
        historical_hash = item.get("historical_certification_binary_sha256")
        if historical_hash is not None and not is_hex(historical_hash, 64):
            raise ReleaseFailure(
                "historical certification hash must be 64 lowercase hex digits or null"
            )
        if item["kind"] == "wasm":
            runner = item.get("runner_script")
            if not isinstance(runner, str) or not runner:
                raise ReleaseFailure("RC wasm artifacts require a runner_script")
        if item["kind"] != "source" and not item.get("binary_name"):
            raise ReleaseFailure("RC executable artifacts require a binary_name")
    validate_smoke_definitions(contract.get("smoke"))
    validate_lineage_shape(contract.get("certification_lineage"))
    publication = contract.get("publication")
    if (
        not isinstance(publication, dict)
        or not publication
        or any(value is not False for value in publication.values())
    ):
        raise ReleaseFailure("RC contract publication flags must all be false")
    toolchain = contract.get("cert_toolchain")
    if not isinstance(toolchain, dict) or any(
        not isinstance(toolchain.get(key), str) or not toolchain[key]
        for key in ("rust", "cargo", "node", "profile")
    ):
        raise ReleaseFailure("RC contract toolchain is incomplete")
    return contract


def validate_smoke_definitions(smoke: object) -> None:
    if not isinstance(smoke, list) or len(smoke) != len(SMOKE_PROGRAMS):
        raise ReleaseFailure("RC contract must declare exactly five smoke cases")
    programs = []
    for case in smoke:
        if not isinstance(case, dict):
            raise ReleaseFailure("RC smoke case is not an object")
        programs.append(case.get("program"))
        arguments = case.get("arguments")
        if not isinstance(arguments, list) or not all(isinstance(a, str) for a in arguments):
            raise ReleaseFailure("RC smoke arguments must be a list of strings")
        if arguments[:1] != [case.get("program")]:
            raise ReleaseFailure("RC smoke arguments must start with the program name")
        for placeholder in SMOKE_PLACEHOLDERS:
            if arguments.count(placeholder) != 1:
                raise ReleaseFailure(f"RC smoke arguments need exactly one {placeholder}")
        if (
            "-num_threads" not in arguments
            or arguments[arguments.index("-num_threads") + 1 : arguments.index("-num_threads") + 2]
            != ["1"]
        ):
            raise ReleaseFailure("RC smoke arguments must pass -num_threads 1")
        if not is_hex(case.get("expected_output_sha256"), 64):
            raise ReleaseFailure("RC smoke expected_output_sha256 must be 64 hex digits")
        for key in ("case_id", "query", "subject"):
            if not isinstance(case.get(key), str) or not case[key]:
                raise ReleaseFailure(f"RC smoke case needs {key}")
    if set(programs) != SMOKE_PROGRAMS or len(set(programs)) != len(programs):
        raise ReleaseFailure("RC smoke cases must cover each of the five programs once")


def validate_lineage_shape(lineage: object) -> None:
    if not isinstance(lineage, dict):
        raise ReleaseFailure("RC contract has no certification_lineage")
    paths = lineage.get("runtime_input_paths")
    excludes = lineage.get("runtime_input_excludes")
    if (
        not isinstance(paths, list)
        or not paths
        or not all(isinstance(item, str) and item for item in paths)
        or not isinstance(excludes, list)
        or not all(isinstance(item, str) and item for item in excludes)
    ):
        raise ReleaseFailure("RC lineage runtime input paths are malformed")
    ancestors = lineage.get("v010_ancestors")
    if not isinstance(ancestors, dict):
        raise ReleaseFailure("RC lineage has no v010_ancestors")
    runs = lineage.get("runs")
    if not isinstance(runs, dict) or sorted(runs) != list(EXPECTED_RUN_KEYS):
        raise ReleaseFailure("RC lineage runs must be B, C, F, G, H, I")
    group = runs["G"]
    group_runs = group.get("runs") if isinstance(group, dict) else None
    required = group.get("required_consecutive_green") if isinstance(group, dict) else None
    if (
        not isinstance(group_runs, list)
        or not isinstance(required, int)
        or isinstance(required, bool)
        or required != 3
        or len(group_runs) != required
    ):
        raise ReleaseFailure("RC lineage run G needs three consecutive green runs")
    local = lineage.get("local_evidence")
    if not isinstance(local, dict) or sorted(local) != list(LOCAL_EVIDENCE_KEYS):
        raise ReleaseFailure("RC lineage local_evidence must be A, D, E")
    ws = lineage.get("ws_evidence")
    if not isinstance(ws, dict) or not isinstance(ws.get("J"), dict) or not ws["J"]:
        raise ReleaseFailure("RC lineage ws_evidence needs J")


def artifact_spec(
    contract: dict[str, object], kind: str, target: str
) -> dict[str, object]:
    return _v010.artifact_spec(contract, kind, target)


def require_run(label: str, run: object) -> None:
    run_id = run.get("run_id") if isinstance(run, dict) else None
    if not isinstance(run_id, int) or isinstance(run_id, bool) or run_id <= 0:
        raise ReleaseFailure(f"contract not filled: run {label} has no positive run_id")
    if run.get("conclusion") != "success":
        raise ReleaseFailure(f"run {label} ({run_id}) did not conclude with success")


def git_is_ancestor(repo_root: Path, ancestor: str, descendant: str) -> bool:
    result = run_capture(
        ["git", "merge-base", "--is-ancestor", ancestor, descendant], repo_root
    )
    if result.returncode not in (0, 1):
        raise ReleaseFailure(
            f"git cannot compare {ancestor} with {descendant}: {result.stderr.strip()}"
        )
    return result.returncode == 0


def validate_candidate(
    repo_root: Path,
    candidate_sha: str,
    contract: dict[str, object],
    *,
    allow_unfilled_lineage: bool = False,
) -> bool:
    """Check the commit and the lineage; return whether the lineage was checked."""
    require_sha(candidate_sha, "candidate SHA")
    head = run_capture(["git", "rev-parse", "HEAD"], repo_root)
    require_success(head, "resolve candidate HEAD")
    if head.stdout.strip() != candidate_sha:
        raise ReleaseFailure(
            f"candidate SHA mismatch: expected {candidate_sha}, observed {head.stdout.strip()}"
        )
    if package_version(repo_root) != contract["package_version"]:
        raise ReleaseFailure("Cargo package version differs from the RC contract")

    if allow_unfilled_lineage:
        if os.environ.get("GITHUB_ACTIONS"):
            raise ReleaseFailure(
                "--allow-unfilled-lineage is a local dry-run flag and is refused under GITHUB_ACTIONS"
            )
        print(
            "WARNING: certification lineage and run identifiers are NOT checked "
            "(--allow-unfilled-lineage); the output is not a release candidate",
            file=sys.stderr,
        )
        return False

    lineage = contract["certification_lineage"]
    lineage_sha = lineage.get("candidate_sha")
    if not is_hex(lineage_sha, 40):
        raise ReleaseFailure(
            "contract not filled: certification_lineage.candidate_sha must be 40 lowercase hex digits"
        )
    if not git_is_ancestor(repo_root, str(lineage_sha), candidate_sha):
        raise ReleaseFailure(
            f"candidate does not descend from the certified candidate {lineage_sha}"
        )
    excludes = [f":(exclude){item}" for item in lineage["runtime_input_excludes"]]
    runtime_diff = run_capture(
        [
            "git",
            "diff",
            "--ignore-cr-at-eol",
            "--quiet",
            str(lineage_sha),
            candidate_sha,
            "--",
            *lineage["runtime_input_paths"],
            *excludes,
        ],
        repo_root,
    )
    if runtime_diff.returncode == 1:
        raise ReleaseFailure(
            "candidate changes a runtime/build/contract input after the certified candidate; "
            "certification reuse is invalid"
        )
    require_success(runtime_diff, "runtime input diff")

    for name, ancestor in sorted(lineage["v010_ancestors"].items()):
        if name == "note":
            continue
        require_sha(str(ancestor), f"v0.1.0 ancestor {name}")
        if not git_is_ancestor(repo_root, str(ancestor), candidate_sha):
            raise ReleaseFailure(
                f"candidate does not descend from v0.1.0 provenance ancestor {name} {ancestor}"
            )

    runs = lineage["runs"]
    for label in GATED_RUN_KEYS:
        require_run(label, runs[label])
    for index, run in enumerate(runs["G"]["runs"], 1):
        require_run(f"G{index}", run)
    for label, item in sorted(lineage["local_evidence"].items()):
        if not isinstance(item, dict) or not is_hex(item.get("evidence_sha256"), 64):
            raise ReleaseFailure(
                f"contract not filled: local evidence {label} has no 64-hex evidence_sha256"
            )
    for label, item in sorted(lineage["ws_evidence"]["J"].items()):
        if not isinstance(item, dict) or not item.get("record"):
            raise ReleaseFailure(f"contract not filled: ws evidence J {label} has no record")
    return True


def version_output(command_prefix: Sequence[str], repo_root: Path, contract: dict[str, object]) -> str:
    result = run_capture(
        [*command_prefix, "--version"], repo_root, clean_losat_environment=True
    )
    require_success(result, "LOSAT --version")
    if result.stderr:
        raise ReleaseFailure(f"LOSAT --version wrote stderr: {result.stderr.strip()}")
    expected = f"losat {contract['package_version']}"
    observed = result.stdout.strip()
    if observed != expected:
        raise ReleaseFailure(
            f"LOSAT version mismatch: expected {expected}, observed {observed}"
        )
    return observed


def wasm_command_prefix(
    node: str, runner: Path, binary: Path
) -> list[str]:
    return [node, "--no-warnings", str(runner), str(binary)]


def smoke_command(
    command_prefix: Sequence[str], repo_root: Path, case: dict[str, object], output: Path
) -> list[str]:
    # The WASI runners preopen "/" so the inputs are passed as absolute paths.
    replacements = {
        "{query}": str((repo_root / str(case["query"])).resolve()),
        "{subject}": str((repo_root / str(case["subject"])).resolve()),
        "{output}": str(output),
    }
    arguments = [replacements.get(str(item), str(item)) for item in case["arguments"]]
    return [*command_prefix, *arguments]


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-107
# The RC smoke invokes already-certified LOSAT CLI arguments for each program; this
# release helper does not implement, normalize, or classify BLAST behavior.
def run_smoke(
    command_prefix: Sequence[str], repo_root: Path, contract: dict[str, object]
) -> dict[str, object]:
    records: list[dict[str, str]] = []
    for case in contract["smoke"]:
        program = str(case["program"])
        for key in ("query", "subject"):
            if not (repo_root / str(case[key])).is_file():
                raise ReleaseFailure(f"RC smoke fixture is missing: {case[key]}")
        with tempfile.TemporaryDirectory(prefix="losat-rc-smoke-") as temp_dir:
            output = Path(temp_dir) / "smoke.out"
            result = run_capture(
                smoke_command(command_prefix, repo_root, case, output),
                repo_root,
                clean_losat_environment=True,
            )
            require_success(result, f"representative {program} smoke")
            if result.stdout or result.stderr:
                raise ReleaseFailure(
                    f"representative {program} smoke wrote unexpected console output"
                )
            if not output.is_file():
                raise ReleaseFailure(f"representative {program} smoke did not create output")
            observed = sha256_path(output)
        expected = str(case["expected_output_sha256"])
        if observed != expected:
            raise ReleaseFailure(
                f"representative {program} smoke hash mismatch: "
                f"expected {expected}, observed {observed}"
            )
        records.append(
            {
                "program": program,
                "case_id": str(case["case_id"]),
                "output_sha256": observed,
                "status": "PASS",
            }
        )
    return {"status": "PASS", "cases": records}


def check_smoke_record(record: object, contract: dict[str, object], label: str) -> None:
    cases = record.get("cases") if isinstance(record, dict) else None
    expected = [
        {
            "program": case["program"],
            "case_id": case["case_id"],
            "output_sha256": case["expected_output_sha256"],
            "status": "PASS",
        }
        for case in contract["smoke"]
    ]
    if not isinstance(record, dict) or record.get("status") != "PASS" or cases != expected:
        raise ReleaseFailure(f"smoke record of {label} differs from the contract")


# NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
# CATCH_ALL(status); return status;
# Host compatibility only: the threaded Wasm artifact is published with the shared-memory
# guard (`LOSAT/tests/wasi_shared_memory.js guard`, release v0.2.0 DECISIONS 52), so that
# browser hosts run guarded bytes too. An artifact whose guard does not verify is refused.
def threaded_guard_record(node: str, repo_root: Path, binary: Path) -> dict[str, object]:
    result = subprocess.run(
        [node, str(repo_root / "LOSAT/tests/wasi_shared_memory.js"), "verify", str(binary)],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        raise ReleaseFailure(
            f"the threaded Wasm artifact has no verified shared-memory guard: {result.stderr.strip()}"
        )
    return json.loads(result.stdout)


def assemble_binary(args: argparse.Namespace, kind: str) -> None:
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    allow_unfilled = bool(getattr(args, "allow_unfilled_lineage", False))
    contract = load_contract(repo_root)
    lineage_checked = validate_candidate(
        repo_root, args.candidate_sha, contract, allow_unfilled_lineage=allow_unfilled
    )
    spec = artifact_spec(contract, kind, args.target)
    binary = args.binary.resolve()
    if not binary.is_file():
        raise ReleaseFailure(f"built artifact is missing: {binary}")
    toolchain = capture_toolchain(
        repo_root,
        args.target,
        contract,
        require_native_host=kind == "native",
        node=args.node if kind == "wasm" else None,
    )
    wasm_identity: object = None
    guard_record: object = None
    runner: Path | None = None
    if kind == "native":
        observed_architecture = binary_architecture(binary)
        required_architecture = expected_architecture(args.target)
        if observed_architecture != required_architecture:
            raise ReleaseFailure(
                f"binary architecture mismatch: expected {required_architecture}, "
                f"observed {observed_architecture}"
            )
        command_prefix = [str(binary)]
    else:
        runner_option = getattr(args, "runner", None)
        runner = (
            runner_option.resolve()
            if runner_option is not None
            else (repo_root / str(spec["runner_script"])).resolve()
        )
        if not runner.is_file():
            raise ReleaseFailure(f"WASI runner script is missing: {runner}")
        # NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
        # CATCH_ALL(status); return status;
        # Reject a reactor or a wrong-threading artifact before running its entry point.
        inspection = subprocess.run(
            [
                args.node,
                str(repo_root / "LOSAT/tests/wasi_artifact.js"),
                str(binary),
                WASM_KIND_BY_TARGET[args.target],
            ],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
        if inspection.returncode:
            raise ReleaseFailure(
                f"invalid {WASM_KIND_BY_TARGET[args.target]} artifact: {inspection.stderr.strip()}"
            )
        wasm_identity = json.loads(inspection.stdout)
        if args.target == "wasm32-wasip1-threads":
            guard_record = threaded_guard_record(args.node, repo_root, binary)
        observed_architecture = "wasm32"
        command_prefix = wasm_command_prefix(args.node, runner, binary)

    observed_version = version_output(command_prefix, repo_root, contract)
    smoke = run_smoke(command_prefix, repo_root, contract)
    binary_sha256 = sha256_path(binary)
    historical_binary_sha256 = spec["historical_certification_binary_sha256"]

    filename = str(spec["filename"])
    archive = output_dir / filename
    root_name = filename.removesuffix(".tar.gz").removesuffix(".zip")
    metadata: dict[str, object] = {
        "schema_version": SIDECAR_SCHEMA_VERSION,
        "release": contract["release"],
        "candidate_sha": args.candidate_sha,
        "lineage_checked": lineage_checked,
        "artifact": {
            "kind": kind,
            "target": args.target,
            "filename": filename,
            "archive_format": spec["archive_format"],
            "features": spec["features"],
            "build_command": spec["build_command"],
        },
        "binary": {
            "filename": spec["binary_name"],
            "sha256": binary_sha256,
            "size": binary.stat().st_size,
            "architecture": observed_architecture,
            "version": observed_version,
        },
        "certification_binary_reference": {
            "sha256": historical_binary_sha256,
            "same_bytes_in_this_build": binary_sha256 == historical_binary_sha256,
            "role": "historical certification provenance; not a cross-runner equality gate",
        },
        "binary_hash_policy": contract["binary_hash_policy"],
        "toolchain": toolchain,
        "pre_archive_smoke": smoke,
        "certification_lineage": contract["certification_lineage"],
        "publication": contract["publication"],
    }
    if kind == "wasm":
        metadata["binary"]["wasi_identity"] = wasm_identity
        metadata["binary"]["runner_script"] = spec["runner_script"]
        if args.target == "wasm32-wasip1-threads":
            metadata["binary"]["shared_memory_guard"] = guard_record
    entries = archive_entries(
        repo_root, root_name, str(spec["binary_name"]), binary, metadata
    )
    if spec["archive_format"] == "zip":
        archive.write_bytes(zip_bytes(entries))
    elif spec["archive_format"] == "tar.gz":
        archive.write_bytes(tar_bytes(entries))
    else:
        raise ReleaseFailure(f"unsupported archive format: {spec['archive_format']}")

    with tempfile.TemporaryDirectory(prefix="losat-rc-extract-") as temp_dir:
        extracted_root = extract_archive(archive, Path(temp_dir), root_name)
        extracted_binary = extracted_root / str(spec["binary_name"])
        if kind == "native":
            extracted_binary.chmod(0o755)
            extracted_prefix = [str(extracted_binary)]
        else:
            assert runner is not None
            extracted_prefix = wasm_command_prefix(args.node, runner, extracted_binary)
        extracted_version = version_output(extracted_prefix, repo_root, contract)
        extracted_smoke = run_smoke(extracted_prefix, repo_root, contract)
        extracted_binary_sha256 = sha256_path(extracted_binary)
        if extracted_binary_sha256 != binary_sha256:
            raise ReleaseFailure("extracted binary differs from the assembled input")
    metadata["extracted_artifact"] = {
        "binary_sha256": extracted_binary_sha256,
        "version": extracted_version,
        "smoke": extracted_smoke,
        "status": "PASS",
    }
    metadata["artifact_sha256"] = sha256_path(archive)
    metadata["status"] = "PASS"
    write_sidecars(archive, output_dir, metadata)
    print(f"{kind.upper()}_ARTIFACT_READY {archive}")


def prepare_source(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    allow_unfilled = bool(getattr(args, "allow_unfilled_lineage", False))
    contract = load_contract(repo_root)
    lineage_checked = validate_candidate(
        repo_root, args.candidate_sha, contract, allow_unfilled_lineage=allow_unfilled
    )
    spec = artifact_spec(contract, "source", "source")
    toolchain = capture_toolchain(
        repo_root,
        "x86_64-unknown-linux-gnu",
        contract,
        require_native_host=True,
    )
    source_crate = args.crate.resolve()
    if source_crate.name != spec["filename"] or not source_crate.is_file():
        raise ReleaseFailure("source crate filename or path differs from RC contract")
    expected_root = f"LOSAT-{contract['package_version']}"
    members = crate_members(source_crate, expected_root)
    with tempfile.TemporaryDirectory(prefix="losat-rc-source-") as temp_dir:
        temp_root = Path(temp_dir)
        extracted = extract_archive(source_crate, temp_root / "extract", expected_root)
        install_root = temp_root / "install"
        build_root = temp_root / "target"
        result = run_capture(
            [
                "cargo",
                "install",
                "--path",
                str(extracted),
                "--locked",
                "--offline",
                "--root",
                str(install_root),
                "--target-dir",
                str(build_root),
            ],
            repo_root,
        )
        require_success(result, "clean source-crate installation")
        installed_binary = install_root / "bin/losat"
        if not installed_binary.is_file():
            raise ReleaseFailure("clean source-crate installation produced no LOSAT binary")
        installed_version = version_output([str(installed_binary)], repo_root, contract)
        smoke = run_smoke([str(installed_binary)], repo_root, contract)

    copied = output_dir / source_crate.name
    shutil.copyfile(source_crate, copied)
    metadata: dict[str, object] = {
        "schema_version": SIDECAR_SCHEMA_VERSION,
        "release": contract["release"],
        "candidate_sha": args.candidate_sha,
        "lineage_checked": lineage_checked,
        "artifact": {
            "kind": "source",
            "target": "source",
            "filename": copied.name,
            "archive_format": "crate",
            "features": spec["features"],
            "build_command": spec["build_command"],
            "sha256": sha256_path(copied),
            "size": copied.stat().st_size,
            "file_count": len(members),
        },
        "toolchain": toolchain,
        "clean_install": {
            "command": "cargo install --path <extracted-crate> --locked --offline",
            "version": installed_version,
            "smoke": smoke,
            "status": "PASS",
        },
        "certification_lineage": contract["certification_lineage"],
        "binary_hash_policy": contract["binary_hash_policy"],
        "publication": contract["publication"],
        "status": "PASS",
    }
    write_sidecars(copied, output_dir, metadata)
    print(f"SOURCE_ARTIFACT_READY {copied}")


def run_id_text(run: object) -> str:
    run_id = run.get("run_id") if isinstance(run, dict) else None
    return str(run_id) if run_id is not None else "not recorded"


def handoff_markdown(
    contract: dict[str, object], candidate_sha: str, workflow_run_id: str
) -> str:
    lineage = contract["certification_lineage"]
    runs = lineage["runs"]
    toolchain = contract["cert_toolchain"]
    lines = [
        "# LOSAT v0.2.0 release-candidate handoff",
        "",
        "`RC_HANDOFF_READY`",
        "",
        f"- Release-candidate SHA (R): `{candidate_sha}`",
        f"- Certified candidate SHA (C): `{lineage['candidate_sha']}`",
        f"- Release-readiness workflow run (I): `{workflow_run_id}`",
        "- Certification runs:",
    ]
    for label in GATED_RUN_KEYS:
        lines.append(f"  - {label} `{runs[label]['workflow']}`: `{run_id_text(runs[label])}`")
    group_ids = ", ".join(run_id_text(run) for run in runs["G"]["runs"])
    lines.append(f"  - G `{runs['G']['workflow']}`: `{group_ids}`")
    lines += [
        f"- Toolchain: Rust `{toolchain['rust']}`, Cargo `{toolchain['cargo']}`, "
        f"Node `{toolchain['node']}`, profile `{toolchain['profile']}`, `--locked`",
        "- Signing/notarization: not performed",
        "- No tag / release / upload performed (tag, GitHub Release, registry, Bioconda)",
        "",
        "## Platform claims",
        "",
        "| Platform | Claim | Evidence sets | Wording |",
        "| --- | --- | --- | --- |",
    ]
    for claim in contract["platform_claims"]:
        evidence = ", ".join(str(item) for item in claim["evidence_sets"])
        lines.append(
            f"| {claim['platform']} | {claim['claim']} | {evidence} | {claim['wording']} |"
        )
    lines += [
        "",
        "Artifact hashes are in `SHA256SUMS`; complete machine-readable provenance is in `RC-HANDOFF.json`.",
        "",
    ]
    return "\n".join(lines)


ARCHIVE_SUFFIXES = (".tar.gz", ".zip", ".crate")


def aggregate(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()
    contract = load_contract(repo_root)
    validate_candidate(repo_root, args.candidate_sha, contract)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ReleaseFailure("aggregate output directory must be empty")

    expected_names = sorted(str(spec["filename"]) for spec in contract["artifacts"])
    observed_names = sorted(
        path.name
        for path in input_root.rglob("*")
        if path.is_file() and path.name.endswith(ARCHIVE_SUFFIXES)
    )
    if observed_names != expected_names:
        raise ReleaseFailure(
            f"downloaded artifact set differs from the eight contract files: {observed_names}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, object]] = []
    checksum_lines: list[str] = []
    for spec in sorted(contract["artifacts"], key=lambda item: str(item["filename"])):
        filename = str(spec["filename"])
        matches = list(input_root.rglob(filename))
        if len(matches) != 1:
            raise ReleaseFailure(f"expected exactly one downloaded {filename}, found {len(matches)}")
        source = matches[0]
        metadata_matches = list(input_root.rglob(f"{filename}.metadata.json"))
        checksum_matches = list(input_root.rglob(f"{filename}.sha256"))
        if len(metadata_matches) != 1 or len(checksum_matches) != 1:
            raise ReleaseFailure(f"{filename} does not have exactly one metadata/checksum sidecar")
        metadata = read_sidecar(metadata_matches[0])
        observed_sha256 = sha256_path(source)
        expected_line = f"{observed_sha256}  {filename}\n"
        if checksum_matches[0].read_text(encoding="utf-8") != expected_line:
            raise ReleaseFailure(f"checksum sidecar mismatch for {filename}")
        if metadata.get("lineage_checked") is not True:
            raise ReleaseFailure(
                f"{filename} was assembled without the lineage check (lineage_checked is not true)"
            )
        if (
            metadata.get("candidate_sha") != args.candidate_sha
            or metadata.get("release") != contract["release"]
            or metadata.get("status") != "PASS"
            or metadata.get("publication") != contract["publication"]
            or metadata.get("binary_hash_policy") != contract["binary_hash_policy"]
        ):
            raise ReleaseFailure(f"artifact metadata identity/status mismatch for {filename}")
        artifact_metadata = metadata.get("artifact")
        if not isinstance(artifact_metadata, dict) or any(
            artifact_metadata.get(key) != spec[key]
            for key in ("kind", "target", "filename")
        ):
            raise ReleaseFailure(f"artifact metadata contract mismatch for {filename}")
        metadata_sha256 = (
            artifact_metadata.get("sha256")
            if spec["kind"] == "source"
            else metadata.get("artifact_sha256")
        )
        if metadata_sha256 != observed_sha256:
            raise ReleaseFailure(f"artifact metadata checksum mismatch for {filename}")
        if metadata.get("certification_lineage") != contract["certification_lineage"]:
            raise ReleaseFailure(f"artifact metadata lineage mismatch for {filename}")
        smoke_holder = (
            metadata.get("clean_install")
            if spec["kind"] == "source"
            else metadata.get("extracted_artifact")
        )
        check_smoke_record(
            smoke_holder.get("smoke") if isinstance(smoke_holder, dict) else None,
            contract,
            filename,
        )
        shutil.copyfile(source, output_dir / filename)
        shutil.copyfile(metadata_matches[0], output_dir / metadata_matches[0].name)
        records.append(
            {
                "filename": filename,
                "kind": spec["kind"],
                "target": spec["target"],
                "sha256": observed_sha256,
                "size": source.stat().st_size,
                "metadata": metadata_matches[0].name,
            }
        )
        checksum_lines.append(f"{observed_sha256}  {filename}\n")

    (output_dir / "SHA256SUMS").write_text(
        "".join(checksum_lines), encoding="utf-8", newline="\n"
    )
    handoff: dict[str, object] = {
        "schema_version": SIDECAR_SCHEMA_VERSION,
        "decision": "RC_HANDOFF_READY",
        "release": contract["release"],
        "candidate_sha": args.candidate_sha,
        "certified_candidate_sha": contract["certification_lineage"]["candidate_sha"],
        "workflow_run_id": str(args.workflow_run_id),
        "artifacts": records,
        "cert_toolchain": contract["cert_toolchain"],
        "platform_claims": contract["platform_claims"],
        "certification_lineage": contract["certification_lineage"],
        "binary_hash_policy": contract["binary_hash_policy"],
        "publication": contract["publication"],
    }
    write_json(output_dir / "RC-HANDOFF.json", handoff)
    (output_dir / "RC-HANDOFF.md").write_text(
        handoff_markdown(contract, args.candidate_sha, str(args.workflow_run_id)),
        encoding="utf-8",
        newline="\n",
    )
    print("RC_HANDOFF_READY")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    def common(subparser: argparse.ArgumentParser, *, dry_run_flag: bool = True) -> None:
        subparser.add_argument("--repo-root", type=Path, required=True)
        subparser.add_argument("--candidate-sha", required=True)
        subparser.add_argument("--output-dir", type=Path, required=True)
        if dry_run_flag:
            subparser.add_argument(
                "--allow-unfilled-lineage",
                action="store_true",
                help="local dry run: skip the lineage and run checks (refused under GITHUB_ACTIONS)",
            )

    native = subparsers.add_parser("native", help="assemble and smoke a native artifact")
    common(native)
    native.add_argument("--target", required=True)
    native.add_argument("--binary", type=Path, required=True)

    wasm = subparsers.add_parser("wasm", help="assemble and smoke a command-Wasm artifact")
    common(wasm)
    wasm.add_argument("--target", default="wasm32-wasip1")
    wasm.add_argument("--binary", type=Path, required=True)
    wasm.add_argument("--node", default="node")
    wasm.add_argument(
        "--runner",
        type=Path,
        default=None,
        help="override the contract's runner_script for the target",
    )

    source = subparsers.add_parser("source", help="verify and smoke a Cargo source crate")
    common(source)
    source.add_argument("--crate", type=Path, required=True)

    aggregate_parser = subparsers.add_parser(
        "aggregate", help="verify and aggregate the exact eight-artifact handoff"
    )
    common(aggregate_parser, dry_run_flag=False)
    aggregate_parser.add_argument("--input-root", type=Path, required=True)
    aggregate_parser.add_argument("--workflow-run-id", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "native":
            assemble_binary(args, "native")
        elif args.command == "wasm":
            if args.target not in WASM_KIND_BY_TARGET:
                raise ReleaseFailure(f"unsupported wasm target: {args.target}")
            assemble_binary(args, "wasm")
        elif args.command == "source":
            prepare_source(args)
        elif args.command == "aggregate":
            aggregate(args)
        else:
            raise ReleaseFailure(f"unknown command: {args.command}")
    except (OSError, ReleaseFailure, subprocess.SubprocessError) as error:
        print(f"RC_HANDOFF_FAILED: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
