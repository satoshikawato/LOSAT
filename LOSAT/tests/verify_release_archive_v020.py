#!/usr/bin/env python3
"""Verify one downloaded LOSAT v0.2.0 release-candidate archive on a fresh runner.

Checks the archive's sha256 against its `.sha256` sidecar, extracts it with the assembler's
path-safe extractor, runs `losat --version` (must be `losat 0.2.0`) and the contract's five
smoke searches from the extracted binary, and prints the smoke table.  The source crate has no
binary: its hash and archive layout are checked.  Any failure exits 1 with a clear message.
Contract, extraction, version and smoke rules are those of `prepare_release_candidate_v020.py`.
"""

from __future__ import annotations

import argparse
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

import prepare_release_candidate_v020 as rc  # noqa: E402


def read_checksum_sidecar(archive: Path) -> str:
    sidecar = archive.with_name(f"{archive.name}.sha256")
    try:
        text = sidecar.read_text(encoding="utf-8")
    except OSError as error:
        raise rc.ReleaseFailure(f"cannot read checksum sidecar {sidecar}: {error}") from error
    parts = text.split()
    if len(parts) != 2 or not rc.is_hex(parts[0], 64) or parts[1] != archive.name:
        raise rc.ReleaseFailure(
            f"checksum sidecar {sidecar.name} is not '<sha256>  {archive.name}'"
        )
    return parts[0]


def check_metadata(archive: Path, archive_sha256: str, binary_sha256: str | None) -> None:
    """Cross-check the metadata sidecar when it travels with the archive."""
    path = archive.with_name(f"{archive.name}.metadata.json")
    if not path.is_file():
        return
    metadata = rc.read_sidecar(path)
    artifact = metadata.get("artifact")
    recorded = metadata.get("artifact_sha256")
    if recorded is None and isinstance(artifact, dict):
        recorded = artifact.get("sha256")
    if recorded != archive_sha256:
        raise rc.ReleaseFailure(
            f"{path.name} records archive sha256 {recorded}, observed {archive_sha256}"
        )
    binary = metadata.get("binary")
    if binary_sha256 is not None and isinstance(binary, dict):
        if binary.get("sha256") != binary_sha256:
            raise rc.ReleaseFailure(
                f"{path.name} records binary sha256 {binary.get('sha256')}, "
                f"extracted {binary_sha256}"
            )


def verify(
    repo_root: Path, archive: Path, target: str, node: str, runner: Path | None
) -> dict[str, object]:
    contract = rc.load_contract(repo_root)
    specs = [item for item in contract["artifacts"] if item["target"] == target]
    if len(specs) != 1:
        raise rc.ReleaseFailure(f"target {target} is not one artifact of the RC contract")
    spec = specs[0]
    kind = str(spec["kind"])
    if archive.name != spec["filename"]:
        raise rc.ReleaseFailure(
            f"archive name {archive.name} differs from the contract's {spec['filename']}"
        )
    if not archive.is_file():
        raise rc.ReleaseFailure(f"archive is missing: {archive}")

    observed = rc.sha256_path(archive)
    expected = read_checksum_sidecar(archive)
    if observed != expected:
        raise rc.ReleaseFailure(
            f"archive sha256 differs from its sidecar: sidecar {expected}, observed {observed}"
        )

    filename = str(spec["filename"])
    root_name = filename.removesuffix(".tar.gz").removesuffix(".zip").removesuffix(".crate")
    result: dict[str, object] = {"target": target, "archive_sha256": observed}
    with tempfile.TemporaryDirectory(prefix="losat-rc-verify-") as temp_dir:
        try:
            extracted_root = rc.extract_archive(archive, Path(temp_dir), root_name)
        except (tarfile.TarError, zipfile.BadZipFile) as error:
            raise rc.ReleaseFailure(f"cannot extract {archive.name}: {error}") from error
        if kind == "source":
            check_metadata(archive, observed, None)
            result["status"] = "PASS"
            return result
        binary = extracted_root / str(spec["binary_name"])
        if not binary.is_file():
            raise rc.ReleaseFailure(f"extracted archive has no {spec['binary_name']}")
        if kind == "native":
            binary.chmod(0o755)
            prefix = [str(binary)]
        else:
            script = runner if runner is not None else repo_root / str(spec["runner_script"])
            if not script.is_file():
                raise rc.ReleaseFailure(f"WASI runner script is missing: {script}")
            prefix = rc.wasm_command_prefix(node, script.resolve(), binary)
        result["version"] = rc.version_output(prefix, repo_root, contract)
        result["smoke"] = rc.run_smoke(prefix, repo_root, contract)
        binary_sha256 = rc.sha256_path(binary)
        result["binary_sha256"] = binary_sha256
        check_metadata(archive, observed, binary_sha256)
    result["status"] = "PASS"
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--node", default="node")
    parser.add_argument("--runner", type=Path, default=None)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = verify(
            args.repo_root.resolve(),
            args.archive.resolve(),
            args.target,
            args.node,
            args.runner,
        )
    except (OSError, rc.ReleaseFailure, rc.subprocess.SubprocessError) as error:
        print(f"RC_ARCHIVE_VERIFY_FAILED: {error}", file=sys.stderr)
        return 1
    print(f"archive  {args.archive.name}")
    print(f"sha256   {result['archive_sha256']} (matches sidecar)")
    if "version" in result:
        print(f"version  {result['version']}")
        print(f"binary   {result['binary_sha256']}")
        print("program   case_id                    output_sha256                                                     status")
        for case in result["smoke"]["cases"]:
            print(
                f"{case['program']:<9} {case['case_id']:<26} "
                f"{case['output_sha256']} {case['status']}"
            )
    print(f"RC_ARCHIVE_VERIFIED {args.archive.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
