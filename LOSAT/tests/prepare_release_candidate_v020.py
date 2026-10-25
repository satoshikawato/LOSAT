#!/usr/bin/env python3
"""Assemble and verify exact-SHA LOSAT v0.2.0 release-candidate artifacts.

Lineage model (docs/release/v0.2.0_certification_plan.md, section 1): every
LOSAT-side certification run uses the candidate commit C. The release-candidate
commit R (HEAD, the commit the RC workflow builds) may differ from C only
outside the runtime inputs listed in the contract. The contract records C, the
GitHub run identifiers and the local evidence hashes; this tool refuses an
unfilled contract, a candidate that is not a descendant of C, any difference
between C and R in a runtime input, and any run that is not a success.

The v0.1.0 assembler (prepare_release_candidate_v010.py) is frozen and is only
imported here for its generic helpers (hashing, deterministic archives, safe
extraction, toolchain capture, crate inspection).

``--allow-unfilled-lineage`` skips only the "filled" checks (candidate SHA,
run identifiers and conclusions, evidence hashes) so that local dry runs and
tests can use the committed, still unfilled contract. ``--allow-node-difference``
(wasm only) accepts a Node release other than the contract's. Both flags are
recorded in the artifact metadata, and ``aggregate`` rejects any artifact that
was assembled with either of them unless aggregate itself is told so. THE
RELEASE WORKFLOW MUST NEVER PASS EITHER FLAG.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence


def _load_v010():
    name = "prepare_release_candidate_v010"
    if name in sys.modules:
        return sys.modules[name]
    path = Path(__file__).with_name(f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


_v010 = _load_v010()

ReleaseFailure = _v010.ReleaseFailure
sha256_path = _v010.sha256_path
json_bytes = _v010.json_bytes
write_json = _v010.write_json
run_capture = _v010.run_capture
require_success = _v010.require_success
require_sha = _v010.require_sha
tar_bytes = _v010.tar_bytes
zip_bytes = _v010.zip_bytes
safe_member_name = _v010.safe_member_name
extract_archive = _v010.extract_archive
archive_entries = _v010.archive_entries
write_sidecars = _v010.write_sidecars
crate_members = _v010.crate_members
read_sidecar = _v010.read_sidecar
expected_architecture = _v010.expected_architecture
version_output = _v010.version_output
package_version = _v010.package_version
artifact_spec = _v010.artifact_spec
capture_toolchain = _v010.capture_toolchain
OUTPUT_ENV_KEYS = _v010.OUTPUT_ENV_KEYS

CONTRACT_PATH = Path("docs/release/v0.2.0_rc_contract.json")
RELEASE = "v0.2.0"
PACKAGE_VERSION = "0.2.0"
BINARY_STEM = "losat"
PROGRAMS = ("blastn", "blastp", "blastx", "tblastn", "tblastx")
PLACEHOLDERS = ("{query}", "{subject}", "{output}")
ARTIFACT_KINDS = {"native": 5, "wasm": 2, "source": 1}
WASM_ARTIFACT_KIND = {
    "wasm32-wasip1": "serial-command",
    "wasm32-wasip1-threads": "threaded-command",
}
# The threaded artifact also runs one smoke with this many threads.
THREAD_SMOKE_PROGRAM = "blastp"
THREAD_SMOKE_THREADS = "2"
SINGLE_RUN_KEYS = ("B", "C", "F", "H", "I")
REQUIRED_SUCCESS_KEYS = ("B", "C", "F", "H")
LOCAL_EVIDENCE_KEYS = ("A", "D", "E")
V010_ANCESTOR_KEYS = ("integrated_runtime", "cross_platform_native", "platform_native_v011")
ARCHIVE_SUFFIX = {"tar.gz": ".tar.gz", "zip": ".zip", "crate": ".crate"}
HEX = "0123456789abcdef"


def is_hex(value: object, length: int) -> bool:
    return isinstance(value, str) and len(value) == length and all(c in HEX for c in value)


def is_run_id(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


# --------------------------------------------------------------------------
# Binary format and architecture
# --------------------------------------------------------------------------

def binary_format(path: Path) -> str:
    """Return elf, pe, macho or wasm from the file's magic bytes."""
    head = path.read_bytes()[:4]
    if head == b"\0asm":
        return "wasm"
    if head == b"\x7fELF":
        return "elf"
    if head[:2] == b"MZ":
        return "pe"
    if head in {b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xfe\xed\xfa\xce"}:
        return "macho"
    raise ReleaseFailure(f"unsupported executable format: {path}")


def binary_architecture(path: Path) -> str:
    """Architecture of an ELF, PE, Mach-O or WebAssembly binary.

    The v0.1.0 function covers the three native formats; the WebAssembly
    magic (a version-1 module) is added here.
    """
    data = path.read_bytes()[:8]
    if data[:4] == b"\0asm":
        if data[4:8] != b"\x01\0\0\0":
            raise ReleaseFailure(f"unsupported WebAssembly module version: {path}")
        return "wasm32"
    return _v010.binary_architecture(path)


def expected_format(target: str) -> str:
    if target.endswith("-windows-msvc"):
        return "pe"
    if target.endswith("-apple-darwin"):
        return "macho"
    if target.endswith("-linux-gnu"):
        return "elf"
    raise ReleaseFailure(f"unsupported native release target: {target}")


# --------------------------------------------------------------------------
# Contract schema
# --------------------------------------------------------------------------

def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ReleaseFailure(f"RC contract: {message}")


def _check_run(run: object, label: str) -> None:
    _require(isinstance(run, dict), f"{label} must be an object")
    run_id = run.get("run_id")
    _require(run_id is None or is_run_id(run_id), f"{label}.run_id must be null or a positive integer")
    conclusion = run.get("conclusion")
    _require(conclusion is None or isinstance(conclusion, str), f"{label}.conclusion must be null or a string")


def _check_artifacts(contract: dict, repo_root: Path) -> None:
    artifacts = contract.get("artifacts")
    _require(isinstance(artifacts, list) and len(artifacts) == 8, "must declare exactly eight artifacts")
    _require(all(isinstance(item, dict) for item in artifacts), "artifacts must be objects")
    filenames = [item.get("filename") for item in artifacts]
    _require(
        all(isinstance(name, str) and name for name in filenames) and len(set(filenames)) == len(filenames),
        "artifact filenames must be unique and non-empty",
    )
    counts = {kind: 0 for kind in ARTIFACT_KINDS}
    seen: set[tuple[object, object]] = set()
    for item in artifacts:
        kind, target = item.get("kind"), item.get("target")
        _require(kind in ARTIFACT_KINDS, f"unknown artifact kind {kind!r}")
        counts[kind] += 1
        _require(isinstance(target, str) and target, "artifact target must be a non-empty string")
        _require((kind, target) not in seen, f"duplicate artifact {kind}/{target}")
        seen.add((kind, target))
        filename, binary_name = item["filename"], item.get("binary_name")
        suffix = ARCHIVE_SUFFIX.get(item.get("archive_format"))
        _require(suffix is not None, f"{filename}: unsupported archive_format")
        _require(filename.endswith(suffix), f"{filename}: filename does not match archive_format")
        _require(
            filename.startswith(f"{contract['package_name']}-{contract['package_version']}"),
            f"{filename}: filename must start with the package name and version",
        )
        for key in ("runner", "features", "build_command"):
            _require(isinstance(item.get(key), str) and item[key], f"{filename}: {key} must be a non-empty string")
        historical = item.get("historical_certification_binary_sha256")
        _require(
            historical is None or is_hex(historical, 64),
            f"{filename}: historical_certification_binary_sha256 must be 64 hex digits or null",
        )
        if kind == "source":
            _require(binary_name is None and item["archive_format"] == "crate", "source artifact has no binary")
            continue
        _require(
            isinstance(binary_name, str) and binary_name == binary_name.lower() and binary_name.startswith(BINARY_STEM),
            f"{filename}: binary_name must be lowercase and start with {BINARY_STEM!r}",
        )
        if kind == "native":
            expected_name = f"{BINARY_STEM}.exe" if target.endswith("-windows-msvc") else BINARY_STEM
            _require(binary_name == expected_name, f"{filename}: native binary_name must be {expected_name}")
            _require(item["archive_format"] == ("zip" if target.endswith("-windows-msvc") else "tar.gz"),
                     f"{filename}: archive_format does not suit the target")
            expected_format(target)
        else:
            _require(target in WASM_ARTIFACT_KIND, f"{filename}: unknown wasm target")
            _require(binary_name == f"{BINARY_STEM}.wasm" and item["archive_format"] == "tar.gz",
                     f"{filename}: wasm binary_name/archive_format")
            runner = item.get("runner_script")
            _require(isinstance(runner, str) and (repo_root / runner).is_file(),
                     f"{filename}: runner_script must name an existing file")
        _require(filename == f"{contract['package_name']}-{contract['package_version']}-{target}{suffix}",
                 f"{filename}: filename must be <package>-<version>-<target><suffix>")
    _require(counts == ARTIFACT_KINDS, f"artifact kinds must be {ARTIFACT_KINDS}, found {counts}")


def _check_smoke(contract: dict, repo_root: Path) -> None:
    smoke = contract.get("smoke")
    _require(isinstance(smoke, list), "smoke must be a list")
    _require(sorted(str(item.get("program")) for item in smoke if isinstance(item, dict)) == sorted(PROGRAMS)
             and len(smoke) == len(PROGRAMS),
             "smoke must have exactly one entry per program " + "/".join(PROGRAMS))
    for item in smoke:
        program = item["program"]
        _require(isinstance(item.get("case_id"), str) and item["case_id"], f"smoke {program}: case_id missing")
        for key in ("query", "subject"):
            value = item.get(key)
            _require(isinstance(value, str) and (repo_root / value).is_file(),
                     f"smoke {program}: {key} must name an existing file")
        _require(is_hex(item.get("expected_output_sha256"), 64),
                 f"smoke {program}: expected_output_sha256 must be 64 hex digits")
        arguments = item.get("arguments")
        _require(isinstance(arguments, list) and all(isinstance(a, str) for a in arguments),
                 f"smoke {program}: arguments must be a list of strings")
        _require(arguments and arguments[0] == program, f"smoke {program}: arguments must start with the program")
        for placeholder in PLACEHOLDERS:
            _require(arguments.count(placeholder) == 1, f"smoke {program}: {placeholder} must appear exactly once")
        _require(arguments[arguments.index("{query}") - 1] == "-query"
                 and arguments[arguments.index("{subject}") - 1] == "-subject"
                 and arguments[arguments.index("{output}") - 1] == "-out",
                 f"smoke {program}: placeholders must follow -query/-subject/-out")
        _require("-num_threads" in arguments
                 and arguments[arguments.index("-num_threads") + 1:][:1] == ["1"],
                 f"smoke {program}: arguments must contain -num_threads 1")


def _check_lineage(contract: dict) -> None:
    lineage = contract.get("certification_lineage")
    _require(isinstance(lineage, dict), "certification_lineage must be an object")
    candidate = lineage.get("candidate_sha")
    _require(candidate == "" or is_hex(candidate, 40), "certification_lineage.candidate_sha must be empty or 40 hex digits")
    for key in ("runtime_input_paths", "runtime_input_excludes"):
        value = lineage.get(key)
        _require(isinstance(value, list) and all(isinstance(p, str) and p for p in value),
                 f"certification_lineage.{key} must be a list of paths")
    _require(lineage["runtime_input_paths"], "certification_lineage.runtime_input_paths must not be empty")
    ancestors = lineage.get("v010_ancestors")
    _require(isinstance(ancestors, dict) and all(is_hex(ancestors.get(k), 40) for k in V010_ANCESTOR_KEYS),
             "certification_lineage.v010_ancestors must hold 40-hex " + ", ".join(V010_ANCESTOR_KEYS))
    runs = lineage.get("runs")
    _require(isinstance(runs, dict) and all(k in runs for k in (*SINGLE_RUN_KEYS, "G")),
             "certification_lineage.runs must have B, C, F, G, H, I")
    for key in SINGLE_RUN_KEYS:
        _check_run(runs[key], f"runs.{key}")
    group = runs["G"]
    _require(isinstance(group, dict), "runs.G must be an object")
    required = group.get("required_consecutive_green")
    _require(isinstance(required, int) and not isinstance(required, bool) and required >= 1,
             "runs.G.required_consecutive_green must be a positive integer")
    _require(isinstance(group.get("runs"), list) and len(group["runs"]) == required,
             "runs.G.runs must hold required_consecutive_green entries")
    for index, run in enumerate(group["runs"]):
        _check_run(run, f"runs.G.runs[{index}]")
    evidence = lineage.get("local_evidence")
    _require(isinstance(evidence, dict) and all(k in evidence for k in LOCAL_EVIDENCE_KEYS),
             "certification_lineage.local_evidence must have A, D, E")
    for key in LOCAL_EVIDENCE_KEYS:
        value = evidence[key].get("evidence_sha256") if isinstance(evidence[key], dict) else 0
        _require(value is None or is_hex(value, 64), f"local_evidence.{key}.evidence_sha256 must be null or 64 hex digits")


def check_contract_schema(contract: object, repo_root: Path) -> None:
    """Schema-only check; accepts the unfilled contract."""
    _require(isinstance(contract, dict), "must be a JSON object")
    _require(contract.get("schema_version") == 2, "unsupported schema_version (expected 2)")
    _require(contract.get("release") == RELEASE and contract.get("package_version") == PACKAGE_VERSION,
             f"release/package_version must be {RELEASE}/{PACKAGE_VERSION}")
    _require(contract.get("package_name") == "LOSAT", "package_name must be LOSAT")
    policy = contract.get("binary_hash_policy")
    _require(isinstance(policy, dict) and policy.get("cross_runner_binary_equality_required") is False,
             "binary_hash_policy.cross_runner_binary_equality_required must be false")
    toolchain = contract.get("cert_toolchain")
    _require(isinstance(toolchain, dict)
             and all(isinstance(toolchain.get(k), str) and toolchain[k] for k in ("rust", "cargo", "node", "profile"))
             and isinstance(toolchain.get("locked"), bool)
             and isinstance(toolchain.get("output_affecting_environment"), dict)
             and sorted(toolchain["output_affecting_environment"]) == sorted(OUTPUT_ENV_KEYS),
             "cert_toolchain is incomplete")
    publication = contract.get("publication")
    _require(isinstance(publication, dict) and publication and not any(publication.values()),
             "publication flags must all be false")
    _check_artifacts(contract, repo_root)
    _check_smoke(contract, repo_root)
    _check_lineage(contract)


def load_contract(repo_root: Path) -> dict:
    path = repo_root / CONTRACT_PATH
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseFailure(f"cannot read RC contract {path}: {error}") from error
    check_contract_schema(contract, repo_root)
    return contract


# --------------------------------------------------------------------------
# Candidate validation
# --------------------------------------------------------------------------

def _run_claims(lineage: dict) -> list[tuple[str, dict]]:
    runs = lineage["runs"]
    claims = [(key, runs[key]) for key in REQUIRED_SUCCESS_KEYS if key in runs]
    claims += [(f"G[{i}]", run) for i, run in enumerate(runs["G"]["runs"])]
    return claims


def check_lineage_filled(contract: dict) -> None:
    """Every certification run and evidence hash must be recorded (not run I)."""
    lineage = contract["certification_lineage"]
    if not is_hex(lineage["candidate_sha"], 40):
        raise ReleaseFailure("certification_lineage.candidate_sha is not filled")
    for label, run in _run_claims(lineage):
        if not is_run_id(run.get("run_id")):
            raise ReleaseFailure(f"certification run {label} has no run_id")
        if run.get("conclusion") != "success":
            raise ReleaseFailure(f"certification run {label} conclusion is {run.get('conclusion')!r}, not 'success'")
    for key in LOCAL_EVIDENCE_KEYS:
        if not is_hex(lineage["local_evidence"][key].get("evidence_sha256"), 64):
            raise ReleaseFailure(f"local evidence {key} has no evidence_sha256")


def _git_flag(repo_root: Path, command: Sequence[str], label: str) -> bool:
    """Run a git command whose exit status 0/1 means yes/no; anything else is an error."""
    result = run_capture(["git", *command], repo_root)
    if result.returncode not in (0, 1):
        require_success(result, label)
    return result.returncode == 0


def validate_candidate(
    repo_root: Path,
    candidate_sha: str,
    contract: dict,
    *,
    allow_unfilled_lineage: bool = False,
) -> None:
    """Check R (HEAD) against the contract. See the module docstring."""
    require_sha(candidate_sha, "candidate SHA")
    head = run_capture(["git", "rev-parse", "HEAD"], repo_root)
    require_success(head, "resolve candidate HEAD")
    if head.stdout.strip() != candidate_sha:
        raise ReleaseFailure(
            f"candidate SHA mismatch: expected {candidate_sha}, observed {head.stdout.strip()}"
        )
    lineage = contract["certification_lineage"]
    if not allow_unfilled_lineage:
        check_lineage_filled(contract)

    for name, ancestor in sorted(lineage["v010_ancestors"].items()):
        if name not in V010_ANCESTOR_KEYS:
            continue
        if not _git_flag(repo_root, ["merge-base", "--is-ancestor", ancestor, candidate_sha], f"ancestry of {name}"):
            raise ReleaseFailure(f"candidate does not descend from the v0.1.0 ancestor {name} ({ancestor})")

    certified = lineage["candidate_sha"]
    if certified:
        if not _git_flag(repo_root, ["merge-base", "--is-ancestor", certified, candidate_sha], "ancestry of C"):
            raise ReleaseFailure(f"candidate does not descend from the certified commit C ({certified})")
        unchanged = _git_flag(
            repo_root,
            [
                "diff",
                "--ignore-cr-at-eol",
                "--quiet",
                certified,
                candidate_sha,
                "--",
                *lineage["runtime_input_paths"],
                *(f":(exclude){path}" for path in lineage["runtime_input_excludes"]),
            ],
            "runtime input diff",
        )
        if not unchanged:
            raise ReleaseFailure(
                "candidate changes a runtime input after the certified commit C; "
                "certification does not cover it"
            )

    if package_version(repo_root) != contract["package_version"]:
        raise ReleaseFailure("Cargo package version differs from the RC contract")


# --------------------------------------------------------------------------
# Smoke
# --------------------------------------------------------------------------

def thread_variant(arguments: Sequence[str], threads: str) -> list[str]:
    result = list(arguments)
    result[result.index("-num_threads") + 1] = threads
    return result


def run_one_smoke(
    command_prefix: Sequence[str],
    repo_root: Path,
    entry: dict,
    *,
    arguments: Sequence[str] | None = None,
    threads: str = "1",
) -> dict:
    query = (repo_root / entry["query"]).resolve()
    subject = (repo_root / entry["subject"]).resolve()
    if not query.is_file() or not subject.is_file():
        raise ReleaseFailure(f"RC smoke fixture is missing for {entry['program']}")
    label = f"{entry['program']} smoke ({entry['case_id']}, {threads} thread(s))"
    with tempfile.TemporaryDirectory(prefix="losat-rc-smoke-") as temp_dir:
        output = Path(temp_dir) / "smoke.out"
        replacements = {"{query}": str(query), "{subject}": str(subject), "{output}": str(output)}
        template = entry["arguments"] if arguments is None else arguments
        resolved = [replacements.get(item, item) for item in template]
        result = run_capture([*command_prefix, *resolved], repo_root, clean_losat_environment=True)
        require_success(result, label)
        if result.stdout or result.stderr:
            raise ReleaseFailure(f"{label} wrote unexpected console output")
        if not output.is_file():
            raise ReleaseFailure(f"{label} did not create its output")
        observed = sha256_path(output)
    expected = entry["expected_output_sha256"]
    if observed != expected:
        raise ReleaseFailure(f"{label} hash mismatch: expected {expected}, observed {observed}")
    return {
        "program": entry["program"],
        "case_id": entry["case_id"],
        "num_threads": int(threads),
        "output_sha256": observed,
        "status": "PASS",
    }


def run_smokes(
    command_prefix: Sequence[str], repo_root: Path, contract: dict, *, threaded: bool = False
) -> list[dict]:
    """Every contract smoke as written; the threaded artifact adds a 2-thread blastp."""
    results = [run_one_smoke(command_prefix, repo_root, entry) for entry in contract["smoke"]]
    if threaded:
        entry = next(e for e in contract["smoke"] if e["program"] == THREAD_SMOKE_PROGRAM)
        results.append(
            run_one_smoke(
                command_prefix,
                repo_root,
                entry,
                arguments=thread_variant(entry["arguments"], THREAD_SMOKE_THREADS),
                threads=THREAD_SMOKE_THREADS,
            )
        )
    return results


# --------------------------------------------------------------------------
# Native and wasm artifacts
# --------------------------------------------------------------------------

def dry_run_flags(args: argparse.Namespace) -> dict[str, bool]:
    return {
        "allow_unfilled_lineage": bool(getattr(args, "allow_unfilled_lineage", False)),
        "allow_node_difference": bool(getattr(args, "allow_node_difference", False)),
    }


def wasm_prefix(node: str, runner: Path, binary: Path) -> list[str]:
    return [node, "--no-warnings", str(runner), str(binary)]


def capture_wasm_toolchain(
    repo_root: Path, target: str, contract: dict, node: str, *, allow_node_difference: bool
) -> dict:
    """capture_toolchain for Wasm; a dry run may accept another Node release."""
    if allow_node_difference:
        result = run_capture([node, "--version"], repo_root)
        require_success(result, "Node identity")
        contract = json.loads(json.dumps(contract))
        contract["cert_toolchain"]["node"] = result.stdout.strip().removeprefix("v")
    return capture_toolchain(repo_root, target, contract, require_native_host=False, node=node)


def assemble_binary(args: argparse.Namespace, kind: str) -> None:
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    contract = load_contract(repo_root)
    validate_candidate(
        repo_root, args.candidate_sha, contract,
        allow_unfilled_lineage=bool(getattr(args, "allow_unfilled_lineage", False)),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    spec = artifact_spec(contract, kind, args.target)
    binary = args.binary.resolve()
    if not binary.is_file():
        raise ReleaseFailure(f"built artifact is missing: {binary}")
    wasm_identity = None
    if kind == "native":
        toolchain = capture_toolchain(repo_root, args.target, contract, require_native_host=True)
        if binary_format(binary) != expected_format(args.target):
            raise ReleaseFailure(
                f"binary format mismatch: expected {expected_format(args.target)}, observed {binary_format(binary)}"
            )
        observed_architecture = binary_architecture(binary)
        required_architecture = expected_architecture(args.target)
        if observed_architecture != required_architecture:
            raise ReleaseFailure(
                f"binary architecture mismatch: expected {required_architecture}, "
                f"observed {observed_architecture}"
            )
        command_prefix = [str(binary)]
    else:
        toolchain = capture_wasm_toolchain(
            repo_root, args.target, contract, args.node,
            allow_node_difference=bool(getattr(args, "allow_node_difference", False)),
        )
        runner = (repo_root / str(spec["runner_script"])).resolve()
        if not args.runner.resolve().samefile(runner):
            raise ReleaseFailure(f"--runner differs from the contract's runner_script {spec['runner_script']}")
        observed_architecture = binary_architecture(binary)
        if observed_architecture != "wasm32":
            raise ReleaseFailure(f"not a WebAssembly module: {binary}")
        inspection = subprocess.run(
            [args.node, str(repo_root / "LOSAT/tests/wasi_artifact.js"), str(binary), WASM_ARTIFACT_KIND[args.target]],
            cwd=repo_root, capture_output=True, text=True, check=False,
        )
        if inspection.returncode:
            raise ReleaseFailure(
                f"invalid {WASM_ARTIFACT_KIND[args.target]} artifact: {inspection.stderr.strip()}"
            )
        wasm_identity = json.loads(inspection.stdout)
        command_prefix = wasm_prefix(args.node, runner, binary)

    expected_version = f"{BINARY_STEM} {contract['package_version']}"
    observed_version = version_output(command_prefix, repo_root)
    if observed_version != expected_version:
        raise ReleaseFailure(
            f"LOSAT version mismatch: expected {expected_version}, observed {observed_version}"
        )
    threaded = args.target == "wasm32-wasip1-threads"
    smokes = run_smokes(command_prefix, repo_root, contract, threaded=threaded)
    binary_sha256 = sha256_path(binary)
    historical = spec["historical_certification_binary_sha256"]

    filename = str(spec["filename"])
    archive = output_dir / filename
    root_name = filename.removesuffix(".tar.gz").removesuffix(".zip")
    metadata: dict = {
        "schema_version": 2,
        "release": contract["release"],
        "candidate_sha": args.candidate_sha,
        "integration_base_sha": contract["integration_base_sha"],
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
            "sha256": historical,
            "same_bytes_in_this_build": historical is not None and binary_sha256 == historical,
            "role": "historical certification provenance; not a cross-runner equality gate",
        },
        "binary_hash_policy": contract["binary_hash_policy"],
        "toolchain": toolchain,
        "pre_archive_smokes": smokes,
        "certification_lineage": contract["certification_lineage"],
        "lineage_validation": dry_run_flags(args),
        "publication": contract["publication"],
    }
    if wasm_identity is not None:
        metadata["binary"]["wasi_identity"] = wasm_identity
        metadata["wasm_runner"] = {"script": spec["runner_script"], "node": toolchain.get("node_V")}
    entries = archive_entries(repo_root, root_name, str(spec["binary_name"]), binary, metadata)
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
            extracted_prefix = wasm_prefix(args.node, runner, extracted_binary)
        extracted_version = version_output(extracted_prefix, repo_root)
        if extracted_version != expected_version:
            raise ReleaseFailure("extracted binary reports a different version")
        extracted_smokes = run_smokes(extracted_prefix, repo_root, contract, threaded=threaded)
        extracted_binary_sha256 = sha256_path(extracted_binary)
        if extracted_binary_sha256 != binary_sha256:
            raise ReleaseFailure("extracted binary differs from the assembled input")
    metadata["extracted_artifact"] = {
        "binary_sha256": extracted_binary_sha256,
        "version": extracted_version,
        "smokes": extracted_smokes,
        "status": "PASS",
    }
    metadata["artifact_sha256"] = sha256_path(archive)
    metadata["status"] = "PASS"
    write_sidecars(archive, output_dir, metadata)
    print(f"{kind.upper()}_ARTIFACT_READY {archive}")


# --------------------------------------------------------------------------
# Source crate
# --------------------------------------------------------------------------

def prepare_source(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    contract = load_contract(repo_root)
    validate_candidate(
        repo_root, args.candidate_sha, contract,
        allow_unfilled_lineage=bool(getattr(args, "allow_unfilled_lineage", False)),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    spec = artifact_spec(contract, "source", "source")
    toolchain = capture_toolchain(repo_root, "x86_64-unknown-linux-gnu", contract, require_native_host=True)
    source_crate = args.crate.resolve()
    if source_crate.name != spec["filename"] or not source_crate.is_file():
        raise ReleaseFailure("source crate filename or path differs from RC contract")
    expected_root = f"{contract['package_name']}-{contract['package_version']}"
    members = crate_members(source_crate, expected_root)
    with tempfile.TemporaryDirectory(prefix="losat-rc-source-") as temp_dir:
        temp_root = Path(temp_dir)
        extracted = extract_archive(source_crate, temp_root / "extract", expected_root)
        install_root = temp_root / "install"
        result = run_capture(
            [
                "cargo", "install", "--path", str(extracted), "--locked", "--offline",
                "--root", str(install_root), "--target-dir", str(temp_root / "target"),
            ],
            repo_root,
        )
        require_success(result, "clean source-crate installation")
        installed_binary = install_root / "bin" / BINARY_STEM
        if not installed_binary.is_file():
            raise ReleaseFailure("clean source-crate installation produced no LOSAT binary")
        installed_version = version_output([str(installed_binary)], repo_root)
        expected_version = f"{BINARY_STEM} {contract['package_version']}"
        if installed_version != expected_version:
            raise ReleaseFailure(
                f"LOSAT version mismatch: expected {expected_version}, observed {installed_version}"
            )
        smokes = run_smokes([str(installed_binary)], repo_root, contract)

    copied = output_dir / source_crate.name
    shutil.copyfile(source_crate, copied)
    metadata: dict = {
        "schema_version": 2,
        "release": contract["release"],
        "candidate_sha": args.candidate_sha,
        "integration_base_sha": contract["integration_base_sha"],
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
            "smokes": smokes,
            "status": "PASS",
        },
        "certification_lineage": contract["certification_lineage"],
        "lineage_validation": dry_run_flags(args),
        "binary_hash_policy": contract["binary_hash_policy"],
        "publication": contract["publication"],
        "status": "PASS",
    }
    write_sidecars(copied, output_dir, metadata)
    print(f"SOURCE_ARTIFACT_READY {copied}")


# --------------------------------------------------------------------------
# Aggregate
# --------------------------------------------------------------------------

def handoff_markdown(
    contract: dict, candidate_sha: str, workflow_run_id: str, records: list[dict]
) -> str:
    lineage = contract["certification_lineage"]
    runs = lineage["runs"]

    def run_line(key: str) -> str:
        run = runs[key]
        return f"| {key} | {run['workflow']} | {run.get('run_id')} | {run.get('conclusion')} |"

    lines = [
        f"# LOSAT {RELEASE} release-candidate handoff",
        "",
        "`RC_HANDOFF_READY`",
        "",
        f"- Candidate SHA (R): `{candidate_sha}`",
        f"- Workflow run: `{workflow_run_id}`",
        "- Signing/notarization: not performed",
        "- Tag/release/registry publication: not performed",
        "",
        "## Artifacts",
        "",
        "| File | Kind | Target | Size | SHA-256 |",
        "| --- | --- | --- | --- | --- |",
        *(
            f"| {r['filename']} | {r['kind']} | {r['target']} | {r['size']} | `{r['sha256']}` |"
            for r in records
        ),
        "",
        "Hashes are also in `SHA256SUMS`; machine-readable provenance is in `RC-HANDOFF.json`.",
        "",
        "## Certification lineage",
        "",
        f"- Certified commit (C): `{lineage['candidate_sha']}`",
        "- Runtime inputs of R equal those of C (checked by `git diff` over: "
        + ", ".join(f"`{p}`" for p in lineage["runtime_input_paths"])
        + "; excluded: "
        + ", ".join(f"`{p}`" for p in lineage["runtime_input_excludes"])
        + ")",
        "- v0.1.0 ancestors (provenance only): "
        + ", ".join(f"{k} `{v}`" for k, v in sorted(lineage["v010_ancestors"].items()) if k in V010_ANCESTOR_KEYS),
        "",
        "| Set | Workflow | Run id | Conclusion |",
        "| --- | --- | --- | --- |",
        *(run_line(key) for key in ("B", "C", "F", "H")),
        *(
            f"| G[{i + 1}] | {runs['G']['workflow']} | {run.get('run_id')} | {run.get('conclusion')} |"
            for i, run in enumerate(runs["G"]["runs"])
        ),
        run_line("I"),
        "",
        "| Local evidence | Tool | evidence_sha256 |",
        "| --- | --- | --- |",
        *(
            f"| {key} ({item['name']}) | {item['tool']} | `{item.get('evidence_sha256')}` |"
            for key, item in sorted(lineage["local_evidence"].items())
        ),
        "",
        "## Publication flags",
        "",
        *(f"- {name}: {str(value).lower()}" for name, value in sorted(contract["publication"].items())),
        "",
    ]
    return "\n".join(lines)


def aggregate(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()
    allow_unfilled = bool(getattr(args, "allow_unfilled_lineage", False))
    contract = load_contract(repo_root)
    validate_candidate(repo_root, args.candidate_sha, contract, allow_unfilled_lineage=allow_unfilled)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ReleaseFailure("aggregate output directory must be empty")

    specs = sorted(contract["artifacts"], key=lambda item: str(item["filename"]))
    allowed = {f"{s['filename']}{suffix}" for s in specs for suffix in ("", ".sha256", ".metadata.json")}
    stray = sorted(
        p.name for p in input_root.rglob("*")
        if p.is_file() and p.name.startswith(f"{contract['package_name']}-") and p.name not in allowed
    )
    if stray:
        raise ReleaseFailure(f"unexpected artifact files in the input: {stray}")

    records: list[dict] = []
    checksum_lines: list[str] = []
    plan: list[tuple[Path, Path, str]] = []
    for spec in specs:
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
        if checksum_matches[0].read_text(encoding="utf-8") != f"{observed_sha256}  {filename}\n":
            raise ReleaseFailure(f"checksum sidecar mismatch for {filename}")
        if (
            metadata.get("candidate_sha") != args.candidate_sha
            or metadata.get("status") != "PASS"
            or metadata.get("release") not in (None, contract["release"])
            or metadata.get("publication") != contract["publication"]
            or metadata.get("binary_hash_policy") != contract["binary_hash_policy"]
        ):
            raise ReleaseFailure(f"artifact metadata identity/status mismatch for {filename}")
        flags = metadata.get("lineage_validation")
        if not isinstance(flags, dict):
            raise ReleaseFailure(f"artifact metadata for {filename} has no lineage_validation record")
        if not allow_unfilled and any(flags.get(key) is not False for key in dry_run_flags(argparse.Namespace())):
            raise ReleaseFailure(f"{filename} was assembled with a dry-run flag; not releasable")
        artifact_metadata = metadata.get("artifact")
        if not isinstance(artifact_metadata, dict) or any(
            artifact_metadata.get(key) != spec[key] for key in ("kind", "target", "filename")
        ):
            raise ReleaseFailure(f"artifact metadata contract mismatch for {filename}")
        metadata_sha256 = (
            artifact_metadata.get("sha256") if spec["kind"] == "source" else metadata.get("artifact_sha256")
        )
        if metadata_sha256 != observed_sha256:
            raise ReleaseFailure(f"artifact metadata checksum mismatch for {filename}")
        if (
            metadata.get("integration_base_sha") != contract["integration_base_sha"]
            or metadata.get("certification_lineage") != contract["certification_lineage"]
        ):
            raise ReleaseFailure(f"artifact metadata lineage mismatch for {filename}")
        plan.append((source, metadata_matches[0], filename))
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

    output_dir.mkdir(parents=True, exist_ok=True)
    for source, metadata_path, filename in plan:
        shutil.copyfile(source, output_dir / filename)
        shutil.copyfile(metadata_path, output_dir / metadata_path.name)
    (output_dir / "SHA256SUMS").write_text("".join(checksum_lines), encoding="utf-8", newline="\n")
    handoff = {
        "schema_version": 2,
        "decision": "RC_HANDOFF_READY",
        "release": contract["release"],
        "candidate_sha": args.candidate_sha,
        "integration_base_sha": contract["integration_base_sha"],
        "workflow_run_id": str(args.workflow_run_id),
        "artifacts": records,
        "certification_lineage": contract["certification_lineage"],
        "binary_hash_policy": contract["binary_hash_policy"],
        "publication": contract["publication"],
    }
    write_json(output_dir / "RC-HANDOFF.json", handoff)
    (output_dir / "RC-HANDOFF.md").write_text(
        handoff_markdown(contract, args.candidate_sha, str(args.workflow_run_id), records),
        encoding="utf-8",
        newline="\n",
    )
    print("RC_HANDOFF_READY")


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    subparsers = parser.add_subparsers(dest="command", required=True)

    def lineage_flag(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument(
            "--allow-unfilled-lineage", action="store_true",
            help="local dry runs and tests only: skip the 'filled' checks; never used by the workflow",
        )

    def common(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--repo-root", type=Path, required=True)
        subparser.add_argument("--candidate-sha", required=True)
        lineage_flag(subparser)

    def with_output(subparser: argparse.ArgumentParser) -> None:
        common(subparser)
        subparser.add_argument("--output-dir", type=Path, required=True)

    check = subparsers.add_parser("check-contract", help="schema check of the contract (works unfilled)")
    check.add_argument("--repo-root", type=Path, required=True)

    validate = subparsers.add_parser("validate", help="schema plus candidate and lineage checks")
    common(validate)

    native = subparsers.add_parser("native", help="assemble and smoke a native artifact")
    with_output(native)
    native.add_argument("--target", required=True)
    native.add_argument("--binary", type=Path, required=True)

    wasm = subparsers.add_parser("wasm", help="assemble and smoke a command-Wasm artifact")
    with_output(wasm)
    wasm.add_argument("--target", choices=sorted(WASM_ARTIFACT_KIND), default="wasm32-wasip1")
    wasm.add_argument("--binary", type=Path, required=True)
    wasm.add_argument("--node", default="node")
    wasm.add_argument("--runner", type=Path, required=True)
    wasm.add_argument(
        "--allow-node-difference", action="store_true",
        help="local dry runs only: accept a Node release other than the contract's; never used by the workflow",
    )

    source = subparsers.add_parser("source", help="verify and smoke a Cargo source crate")
    with_output(source)
    source.add_argument("--crate", type=Path, required=True)

    aggregate_parser = subparsers.add_parser("aggregate", help="verify and aggregate the exact eight-artifact handoff")
    with_output(aggregate_parser)
    aggregate_parser.add_argument("--input-root", type=Path, required=True)
    aggregate_parser.add_argument("--workflow-run-id", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "check-contract":
            load_contract(args.repo_root.resolve())
            print("RC_CONTRACT_OK")
        elif args.command == "validate":
            repo_root = args.repo_root.resolve()
            contract = load_contract(repo_root)
            validate_candidate(
                repo_root, args.candidate_sha, contract, allow_unfilled_lineage=args.allow_unfilled_lineage
            )
            print("RC_CANDIDATE_OK")
        elif args.command in ("native", "wasm"):
            assemble_binary(args, args.command)
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
