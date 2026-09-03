#!/usr/bin/env python3
"""Assemble and verify exact-SHA LOSAT v0.1.0 release-candidate artifacts."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import tomllib
import zipfile
from pathlib import Path, PurePosixPath
from typing import Iterable, Sequence


CONTRACT_PATH = Path("docs/release/v0.1.0_rc_contract.json")
RUNTIME_INPUT_PATHS = (
    "LOSAT/src",
    "LOSAT/Cargo.toml",
    "LOSAT/Cargo.lock",
    "LOSAT/build.rs",
    "LOSAT/.cargo",
    "LOSAT/tests/blastn_parity_manifest.tsv",
    "LOSAT/tests/blastn_v010_source_exceptions.tsv",
    "LOSAT/tests/blastp_v010_parity_manifest.tsv",
    "LOSAT/tests/tblastx_v010_parity_manifest.tsv",
    "LOSAT/tests/ncbi_platform_variance_v010.json",
    "LOSAT/tests/platform_native_v010_canonical.tsv",
    "docs/product_decisions",
)
OUTPUT_ENV_KEYS = (
    "CARGO_BUILD_RUSTFLAGS",
    "CARGO_ENCODED_RUSTFLAGS",
    "CC",
    "CFLAGS",
    "CPPFLAGS",
    "LDFLAGS",
    "RUSTC_LINKER",
    "RUSTFLAGS",
)


class ReleaseFailure(RuntimeError):
    """A fail-closed release-candidate gate failure."""


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(json_bytes(value))


def run_capture(
    command: Sequence[str], cwd: Path, *, clean_losat_environment: bool = False
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    if clean_losat_environment:
        for name in tuple(environment):
            if name.startswith("LOSAT_"):
                del environment[name]
    return subprocess.run(
        list(command),
        cwd=cwd,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def require_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode != 0:
        details = (result.stderr or result.stdout).strip()
        raise ReleaseFailure(f"{label} failed ({result.returncode}): {details}")


def require_sha(value: str, label: str) -> str:
    if len(value) != 40 or any(character not in "0123456789abcdef" for character in value):
        raise ReleaseFailure(f"{label} must be exactly 40 lowercase hexadecimal digits")
    return value


def load_contract(repo_root: Path) -> dict[str, object]:
    path = repo_root / CONTRACT_PATH
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseFailure(f"cannot read RC contract {path}: {error}") from error
    if contract.get("schema_version") != 1:
        raise ReleaseFailure("unsupported RC contract schema")
    if contract.get("release") != "v0.1.0" or contract.get("package_version") != "0.1.0":
        raise ReleaseFailure("RC contract release/version changed")
    artifacts = contract.get("artifacts")
    if not isinstance(artifacts, list) or len(artifacts) != 6:
        raise ReleaseFailure("RC contract must declare exactly six artifacts")
    filenames = [item.get("filename") for item in artifacts if isinstance(item, dict)]
    if len(filenames) != len(set(filenames)) or any(not name for name in filenames):
        raise ReleaseFailure("RC artifact filenames must be unique and non-empty")
    return contract


def artifact_spec(
    contract: dict[str, object], kind: str, target: str
) -> dict[str, object]:
    matches = [
        item
        for item in contract["artifacts"]
        if isinstance(item, dict)
        and item.get("kind") == kind
        and item.get("target") == target
    ]
    if len(matches) != 1:
        raise ReleaseFailure(f"RC contract has no unique {kind}/{target} artifact")
    return matches[0]


def package_version(repo_root: Path) -> str:
    cargo_toml = repo_root / "LOSAT/Cargo.toml"
    with cargo_toml.open("rb") as handle:
        manifest = tomllib.load(handle)
    value = manifest.get("package", {}).get("version")
    if not isinstance(value, str):
        raise ReleaseFailure("LOSAT package version is missing")
    return value


def validate_candidate(
    repo_root: Path, candidate_sha: str, contract: dict[str, object]
) -> None:
    require_sha(candidate_sha, "candidate SHA")
    head = run_capture(["git", "rev-parse", "HEAD"], repo_root)
    require_success(head, "resolve candidate HEAD")
    if head.stdout.strip() != candidate_sha:
        raise ReleaseFailure(
            f"candidate SHA mismatch: expected {candidate_sha}, observed {head.stdout.strip()}"
        )

    lineage = contract["certification_lineage"]
    required_ancestors = (
        contract["integration_base_sha"],
        lineage["integrated_runtime"]["sha"],
        lineage["cross_platform_native"]["sha"],
    )
    for ancestor in required_ancestors:
        require_sha(str(ancestor), "certification ancestor SHA")
        result = run_capture(
            ["git", "merge-base", "--is-ancestor", str(ancestor), candidate_sha],
            repo_root,
        )
        if result.returncode != 0:
            raise ReleaseFailure(
                f"candidate does not descend from required certification ancestor {ancestor}"
            )

    certified_runtime_sha = str(lineage["cross_platform_native"]["sha"])
    runtime_diff = run_capture(
        [
            "git",
            "diff",
            "--ignore-cr-at-eol",
            "--quiet",
            certified_runtime_sha,
            candidate_sha,
            "--",
            *RUNTIME_INPUT_PATHS,
        ],
        repo_root,
    )
    if runtime_diff.returncode != 0:
        raise ReleaseFailure(
            "candidate changes a runtime/build/contract authority after PR 6; "
            "certification reuse is invalid"
        )
    if package_version(repo_root) != contract["package_version"]:
        raise ReleaseFailure("Cargo package version differs from the RC contract")


def capture_toolchain(
    repo_root: Path,
    target: str,
    contract: dict[str, object],
    *,
    require_native_host: bool,
    node: str | None = None,
) -> dict[str, object]:
    rustc = run_capture(["rustc", "-vV"], repo_root)
    cargo = run_capture(["cargo", "-V"], repo_root)
    require_success(rustc, "rustc identity")
    require_success(cargo, "Cargo identity")
    authority = contract["cert_toolchain"]
    rust_version = str(authority["rust"])
    cargo_version = str(authority["cargo"])
    if not rustc.stdout.startswith(f"rustc {rust_version} "):
        raise ReleaseFailure("rustc differs from CERT_TOOLCHAIN")
    if not cargo.stdout.startswith(f"cargo {cargo_version} "):
        raise ReleaseFailure("Cargo differs from CERT_TOOLCHAIN")
    host = next(
        (
            line.removeprefix("host: ")
            for line in rustc.stdout.splitlines()
            if line.startswith("host: ")
        ),
        "",
    )
    if require_native_host and host != target:
        raise ReleaseFailure(f"native target/host mismatch: target={target}, host={host}")
    expected_environment = authority["output_affecting_environment"]
    environment = {key: os.environ.get(key, "") for key in OUTPUT_ENV_KEYS}
    if environment != expected_environment:
        raise ReleaseFailure(
            f"output-affecting build environment differs from RC contract: {environment}"
        )
    metadata: dict[str, object] = {
        "rustc_vV": rustc.stdout.splitlines(),
        "cargo_V": cargo.stdout.strip(),
        "host": host,
        "target": target,
        "profile": authority["profile"],
        "locked": authority["locked"],
        "environment": environment,
        "python": sys.version,
        "runner": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "RUNNER_OS": os.environ.get("RUNNER_OS", ""),
            "RUNNER_ARCH": os.environ.get("RUNNER_ARCH", ""),
            "ImageOS": os.environ.get("ImageOS", ""),
            "ImageVersion": os.environ.get("ImageVersion", ""),
        },
    }
    if node is not None:
        node_result = run_capture([node, "--version"], repo_root)
        require_success(node_result, "Node identity")
        if node_result.stdout.strip() != f"v{authority['node']}":
            raise ReleaseFailure("Node differs from CERT_TOOLCHAIN")
        metadata["node_V"] = node_result.stdout.strip()
    return metadata


def binary_architecture(path: Path) -> str:
    data = path.read_bytes()[:4096]
    if data.startswith(b"\x7fELF") and len(data) >= 20:
        byteorder = "little" if data[5] == 1 else "big"
        machine = int.from_bytes(data[18:20], byteorder)
        return {62: "x86_64", 183: "arm64"}.get(machine, "unknown")
    if data.startswith(b"MZ") and len(data) >= 64:
        pe_offset = int.from_bytes(data[60:64], "little")
        with path.open("rb") as handle:
            handle.seek(pe_offset)
            header = handle.read(6)
        if header[:4] != b"PE\0\0":
            raise ReleaseFailure(f"invalid PE header: {path}")
        return {0x8664: "x86_64", 0xAA64: "arm64"}.get(
            int.from_bytes(header[4:6], "little"), "unknown"
        )
    if len(data) >= 8 and data[:4] in {b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe"}:
        cpu_type = int.from_bytes(data[4:8], "little")
        return {0x01000007: "x86_64", 0x0100000C: "arm64"}.get(cpu_type, "unknown")
    if len(data) >= 8 and data[:4] in {b"\xfe\xed\xfa\xcf", b"\xfe\xed\xfa\xce"}:
        cpu_type = int.from_bytes(data[4:8], "big")
        return {0x01000007: "x86_64", 0x0100000C: "arm64"}.get(cpu_type, "unknown")
    raise ReleaseFailure(f"unsupported native executable format: {path}")


def expected_architecture(target: str) -> str:
    if target.startswith("x86_64-"):
        return "x86_64"
    if target.startswith("aarch64-"):
        return "arm64"
    raise ReleaseFailure(f"unsupported native release target: {target}")


def version_output(command_prefix: Sequence[str], repo_root: Path) -> str:
    result = run_capture(
        [*command_prefix, "--version"], repo_root, clean_losat_environment=True
    )
    require_success(result, "LOSAT --version")
    if result.stderr:
        raise ReleaseFailure(f"LOSAT --version wrote stderr: {result.stderr.strip()}")
    return result.stdout.strip()


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-107
# The RC smoke invokes the already-certified LOSAT BLASTP CLI arguments; this
# release helper does not implement, normalize, or classify BLAST behavior.
def run_smoke(
    command_prefix: Sequence[str], repo_root: Path, contract: dict[str, object]
) -> dict[str, object]:
    smoke = contract["smoke"]
    query = repo_root / str(smoke["query"])
    subject = repo_root / str(smoke["subject"])
    if not query.is_file() or not subject.is_file():
        raise ReleaseFailure("RC smoke fixture is missing")
    with tempfile.TemporaryDirectory(prefix="losat-rc-smoke-") as temp_dir:
        output = Path(temp_dir) / "smoke.out"
        replacements = {
            "{query}": str(query),
            "{subject}": str(subject),
            "{output}": str(output),
        }
        arguments = [replacements.get(str(item), str(item)) for item in smoke["arguments"]]
        result = run_capture(
            [*command_prefix, *arguments],
            repo_root,
            clean_losat_environment=True,
        )
        require_success(result, "representative BLASTP smoke")
        if result.stdout or result.stderr:
            raise ReleaseFailure("representative BLASTP smoke wrote unexpected console output")
        if not output.is_file():
            raise ReleaseFailure("representative BLASTP smoke did not create output")
        observed = sha256_path(output)
    expected = str(smoke["expected_output_sha256"])
    if observed != expected:
        raise ReleaseFailure(
            f"representative BLASTP smoke hash mismatch: expected {expected}, observed {observed}"
        )
    return {
        "program": smoke["program"],
        "case_id": smoke["case_id"],
        "output_sha256": observed,
        "status": "PASS",
    }


def tar_bytes(entries: Iterable[tuple[str, bytes, int]]) -> bytes:
    with tempfile.TemporaryFile() as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            with tarfile.open(
                fileobj=compressed, mode="w", format=tarfile.USTAR_FORMAT
            ) as archive:
                for name, data, mode in sorted(entries):
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mode = mode
                    info.mtime = 0
                    info.uid = 0
                    info.gid = 0
                    info.uname = ""
                    info.gname = ""
                    archive.addfile(info, fileobj=_BytesReader(data))
        raw.seek(0)
        return raw.read()


class _BytesReader:
    def __init__(self, data: bytes):
        self.data = data
        self.offset = 0

    def read(self, size: int = -1) -> bytes:
        if size < 0:
            size = len(self.data) - self.offset
        chunk = self.data[self.offset : self.offset + size]
        self.offset += len(chunk)
        return chunk


def zip_bytes(entries: Iterable[tuple[str, bytes, int]]) -> bytes:
    with tempfile.TemporaryFile() as raw:
        with zipfile.ZipFile(
            raw, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
        ) as archive:
            for name, data, mode in sorted(entries):
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.create_system = 3
                info.external_attr = (stat.S_IFREG | mode) << 16
                archive.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
        raw.seek(0)
        return raw.read()


def safe_member_name(name: str, expected_root: str) -> None:
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ReleaseFailure(f"unsafe archive member: {name}")
    if path.parts[0] != expected_root:
        raise ReleaseFailure(f"archive member is outside expected root {expected_root}: {name}")


def extract_archive(archive: Path, destination: Path, expected_root: str) -> Path:
    if archive.suffix == ".zip":
        with zipfile.ZipFile(archive) as handle:
            for info in handle.infolist():
                safe_member_name(info.filename, expected_root)
                mode = info.external_attr >> 16
                if stat.S_ISLNK(mode):
                    raise ReleaseFailure(f"archive contains a symlink: {info.filename}")
            handle.extractall(destination)
    else:
        with tarfile.open(archive, "r:gz") as handle:
            for member in handle.getmembers():
                safe_member_name(member.name, expected_root)
                if not member.isfile():
                    raise ReleaseFailure(f"archive contains a non-file member: {member.name}")
            handle.extractall(destination, filter="data")
    root = destination / expected_root
    if not root.is_dir():
        raise ReleaseFailure("extracted artifact root is missing")
    return root


def archive_entries(
    repo_root: Path,
    root_name: str,
    binary_name: str,
    binary: Path,
    metadata: dict[str, object],
) -> list[tuple[str, bytes, int]]:
    return [
        (f"{root_name}/{binary_name}", binary.read_bytes(), 0o755),
        (f"{root_name}/LICENSE", (repo_root / "LICENSE").read_bytes(), 0o644),
        (f"{root_name}/README.md", (repo_root / "README.md").read_bytes(), 0o644),
        (f"{root_name}/RELEASE-METADATA.json", json_bytes(metadata), 0o644),
    ]


def write_sidecars(
    artifact: Path, output_dir: Path, metadata: dict[str, object]
) -> None:
    write_json(output_dir / f"{artifact.name}.metadata.json", metadata)
    checksum = sha256_path(artifact)
    (output_dir / f"{artifact.name}.sha256").write_text(
        f"{checksum}  {artifact.name}\n", encoding="utf-8", newline="\n"
    )


def assemble_binary(args: argparse.Namespace, kind: str) -> None:
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = load_contract(repo_root)
    validate_candidate(repo_root, args.candidate_sha, contract)
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
        if binary.read_bytes()[:4] != b"\0asm":
            raise ReleaseFailure("serial-Wasm artifact has no WebAssembly magic")
        observed_architecture = "wasm32"
        command_prefix = [args.node, "--no-warnings", str(args.runner.resolve()), str(binary)]

    expected_version = f"losat {contract['package_version']}"
    observed_version = version_output(command_prefix, repo_root)
    if observed_version != expected_version:
        raise ReleaseFailure(
            f"LOSAT version mismatch: expected {expected_version}, observed {observed_version}"
        )
    smoke = run_smoke(command_prefix, repo_root, contract)
    binary_sha256 = sha256_path(binary)
    certified_binary_sha256 = str(spec["certified_binary_sha256"])
    if binary_sha256 != certified_binary_sha256:
        raise ReleaseFailure(
            f"{args.target} binary does not reproduce its certified hash: "
            f"expected {certified_binary_sha256}, observed {binary_sha256}"
        )

    filename = str(spec["filename"])
    archive = output_dir / filename
    root_name = filename.removesuffix(".tar.gz").removesuffix(".zip")
    metadata: dict[str, object] = {
        "schema_version": 1,
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
            "certified_binary_sha256": certified_binary_sha256,
            "certified_hash_reproduced": True,
        },
        "toolchain": toolchain,
        "pre_archive_smoke": smoke,
        "certification_lineage": contract["certification_lineage"],
        "publication": contract["publication"],
    }
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
            extracted_prefix = [
                args.node,
                "--no-warnings",
                str(args.runner.resolve()),
                str(extracted_binary),
            ]
        extracted_version = version_output(extracted_prefix, repo_root)
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


def crate_members(crate: Path, expected_root: str) -> list[str]:
    members: list[str] = []
    manifest_bytes: bytes | None = None
    with tarfile.open(crate, "r:gz") as handle:
        for member in handle.getmembers():
            safe_member_name(member.name, expected_root)
            if not member.isfile():
                raise ReleaseFailure(f"source crate contains a non-file member: {member.name}")
            members.append(member.name)
            if member.name == f"{expected_root}/Cargo.toml":
                extracted = handle.extractfile(member)
                if extracted is None:
                    raise ReleaseFailure("source crate Cargo.toml cannot be read")
                manifest_bytes = extracted.read()
    required = {
        f"{expected_root}/Cargo.toml",
        f"{expected_root}/Cargo.lock",
        f"{expected_root}/README.md",
    }
    missing = required.difference(members)
    if missing:
        raise ReleaseFailure(f"source crate is missing required files: {sorted(missing)}")
    try:
        manifest = tomllib.loads((manifest_bytes or b"").decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise ReleaseFailure(f"source crate Cargo.toml is invalid: {error}") from error
    package = manifest.get("package", {})
    license_id = package.get("license")
    license_file = package.get("license-file") or package.get("license_file")
    if not isinstance(license_id, str) or not license_id.strip():
        if not isinstance(license_file, str) or not license_file.strip():
            raise ReleaseFailure("source crate has no license metadata")
        if f"{expected_root}/{license_file}" not in members:
            raise ReleaseFailure("source crate license-file is missing")
    forbidden_parts = {"target", "__pycache__", ".tmp", "losat_out", "blast_out", "ncbi_out"}
    for name in members:
        if forbidden_parts.intersection(PurePosixPath(name).parts):
            raise ReleaseFailure(f"source crate contains generated/scratch content: {name}")
        if name.endswith((".pyc", ".tmp")) or (
            name.endswith(".orig") and name != f"{expected_root}/Cargo.toml.orig"
        ):
            raise ReleaseFailure(f"source crate contains generated/scratch content: {name}")
    return sorted(members)


def prepare_source(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = load_contract(repo_root)
    validate_candidate(repo_root, args.candidate_sha, contract)
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
        installed_binary = install_root / "bin/LOSAT"
        if not installed_binary.is_file():
            raise ReleaseFailure("clean source-crate installation produced no LOSAT binary")
        installed_version = version_output([str(installed_binary)], repo_root)
        smoke = run_smoke([str(installed_binary)], repo_root, contract)

    copied = output_dir / source_crate.name
    shutil.copyfile(source_crate, copied)
    metadata: dict[str, object] = {
        "schema_version": 1,
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
            "smoke": smoke,
            "status": "PASS",
        },
        "certification_lineage": contract["certification_lineage"],
        "publication": contract["publication"],
        "status": "PASS",
    }
    write_sidecars(copied, output_dir, metadata)
    print(f"SOURCE_ARTIFACT_READY {copied}")


def read_sidecar(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseFailure(f"cannot read artifact metadata {path}: {error}") from error
    if not isinstance(value, dict):
        raise ReleaseFailure(f"artifact metadata is not an object: {path}")
    return value


def aggregate(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()
    contract = load_contract(repo_root)
    validate_candidate(repo_root, args.candidate_sha, contract)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ReleaseFailure("aggregate output directory must be empty")
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
        if (
            metadata.get("candidate_sha") != args.candidate_sha
            or metadata.get("status") != "PASS"
            or metadata.get("publication") != contract["publication"]
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
        if (
            metadata.get("integration_base_sha") != contract["integration_base_sha"]
            or metadata.get("certification_lineage")
            != contract["certification_lineage"]
        ):
            raise ReleaseFailure(f"artifact metadata lineage mismatch for {filename}")
        destination = output_dir / filename
        shutil.copyfile(source, destination)
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
    rerun = args.integrated_certification_rerun == "true"
    handoff: dict[str, object] = {
        "schema_version": 1,
        "decision": "RC_HANDOFF_READY",
        "release": contract["release"],
        "candidate_sha": args.candidate_sha,
        "integration_base_sha": contract["integration_base_sha"],
        "workflow_run_id": str(args.workflow_run_id),
        "artifacts": records,
        "certification_lineage": contract["certification_lineage"],
        "integrated_certification_rerun": rerun,
        "integrated_certification_rerun_reason": (
            "explicit workflow input requested a fresh integrated run"
            if rerun
            else contract["certification_lineage"]["rerun_reason"]
        ),
        "publication": contract["publication"],
    }
    write_json(output_dir / "RC-HANDOFF.json", handoff)
    markdown = [
        "# LOSAT v0.1.0 release-candidate handoff",
        "",
        "`RC_HANDOFF_READY`",
        "",
        f"- Candidate SHA: `{args.candidate_sha}`",
        f"- Workflow run: `{args.workflow_run_id}`",
        f"- Integrated certification rerun: `{'YES' if rerun else 'NO'}`",
        "- Signing/notarization: not performed",
        "- Tag/release/registry publication: not performed",
        "",
        "Artifact hashes are in `SHA256SUMS`; complete machine-readable provenance is in `RC-HANDOFF.json`.",
        "",
    ]
    (output_dir / "RC-HANDOFF.md").write_text(
        "\n".join(markdown), encoding="utf-8", newline="\n"
    )
    print("RC_HANDOFF_READY")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    def common(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--repo-root", type=Path, required=True)
        subparser.add_argument("--candidate-sha", required=True)
        subparser.add_argument("--output-dir", type=Path, required=True)

    native = subparsers.add_parser("native", help="assemble and smoke a native artifact")
    common(native)
    native.add_argument("--target", required=True)
    native.add_argument("--binary", type=Path, required=True)

    wasm = subparsers.add_parser("wasm", help="assemble and smoke serial command-Wasm")
    common(wasm)
    wasm.add_argument("--target", default="wasm32-wasip1")
    wasm.add_argument("--binary", type=Path, required=True)
    wasm.add_argument("--node", default="node")
    wasm.add_argument("--runner", type=Path, required=True)

    source = subparsers.add_parser("source", help="verify and smoke a Cargo source crate")
    common(source)
    source.add_argument("--crate", type=Path, required=True)

    aggregate_parser = subparsers.add_parser(
        "aggregate", help="verify and aggregate the exact six-artifact handoff"
    )
    common(aggregate_parser)
    aggregate_parser.add_argument("--input-root", type=Path, required=True)
    aggregate_parser.add_argument("--workflow-run-id", required=True)
    aggregate_parser.add_argument(
        "--integrated-certification-rerun", choices=("true", "false"), required=True
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "native":
            assemble_binary(args, "native")
        elif args.command == "wasm":
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
