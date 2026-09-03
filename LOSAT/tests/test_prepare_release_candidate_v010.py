#!/usr/bin/env python3
"""Focused tests for the v0.1.0 exact-SHA release-candidate handoff."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import subprocess
import tarfile
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from unittest import mock


MODULE_PATH = Path(__file__).with_name("prepare_release_candidate_v010.py")
SPEC = importlib.util.spec_from_file_location("prepare_release_candidate_v010", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
release = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release)

REPO_ROOT = Path(__file__).resolve().parents[2]


class ReleaseCandidateTests(unittest.TestCase):
    def test_contract_has_one_exact_artifact_per_supported_target(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        observed = {
            (item["kind"], item["target"], item["filename"])
            for item in contract["artifacts"]
        }
        self.assertEqual(
            observed,
            {
                (
                    "native",
                    "x86_64-unknown-linux-gnu",
                    "LOSAT-0.1.0-x86_64-unknown-linux-gnu.tar.gz",
                ),
                (
                    "native",
                    "x86_64-pc-windows-msvc",
                    "LOSAT-0.1.0-x86_64-pc-windows-msvc.zip",
                ),
                (
                    "native",
                    "aarch64-apple-darwin",
                    "LOSAT-0.1.0-aarch64-apple-darwin.tar.gz",
                ),
                (
                    "native",
                    "x86_64-apple-darwin",
                    "LOSAT-0.1.0-x86_64-apple-darwin.tar.gz",
                ),
                (
                    "wasm",
                    "wasm32-wasip1",
                    "LOSAT-0.1.0-wasm32-wasip1.tar.gz",
                ),
                ("source", "source", "LOSAT-0.1.0.crate"),
            },
        )
        self.assertEqual(
            contract["certification_lineage"]["integrated_runtime"]["native_contracts"],
            43,
        )
        self.assertEqual(
            contract["certification_lineage"]["integrated_runtime"][
                "serial_wasm_equalities"
            ],
            41,
        )
        self.assertEqual(
            contract["certification_lineage"]["cross_platform_native"][
                "platform_executions"
            ],
            183,
        )
        self.assertFalse(any(contract["publication"].values()))
        self.assertFalse(
            contract["binary_hash_policy"]["cross_runner_binary_equality_required"]
        )
        for item in contract["artifacts"]:
            historical_hash = item["historical_certification_binary_sha256"]
            if item["kind"] == "source":
                self.assertIsNone(historical_hash)
            else:
                self.assertRegex(historical_hash, r"^[0-9a-f]{64}$")

    def test_tar_and_zip_assembly_are_deterministic(self) -> None:
        entries = [
            ("LOSAT-0.1.0-test/LOSAT", b"binary\0bytes", 0o755),
            ("LOSAT-0.1.0-test/README.md", b"readme\n", 0o644),
        ]
        first_tar = release.tar_bytes(entries)
        second_tar = release.tar_bytes(reversed(entries))
        self.assertEqual(hashlib.sha256(first_tar).digest(), hashlib.sha256(second_tar).digest())
        first_zip = release.zip_bytes(entries)
        second_zip = release.zip_bytes(reversed(entries))
        self.assertEqual(hashlib.sha256(first_zip).digest(), hashlib.sha256(second_zip).digest())

    def test_archive_extraction_rejects_path_escape(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            archive = root / "unsafe.tar.gz"
            with tarfile.open(archive, "w:gz") as handle:
                info = tarfile.TarInfo("LOSAT-0.1.0-test/../outside")
                data = b"unsafe"
                info.size = len(data)
                handle.addfile(info, io.BytesIO(data))
            with self.assertRaisesRegex(release.ReleaseFailure, "unsafe archive member"):
                release.extract_archive(archive, root / "extract", "LOSAT-0.1.0-test")

    def test_binary_architecture_recognizes_supported_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            elf = bytearray(64)
            elf[:6] = b"\x7fELF\x02\x01"
            elf[18:20] = (62).to_bytes(2, "little")
            elf_path = root / "elf"
            elf_path.write_bytes(elf)
            self.assertEqual(release.binary_architecture(elf_path), "x86_64")

            pe = bytearray(128)
            pe[:2] = b"MZ"
            pe[60:64] = (80).to_bytes(4, "little")
            pe[80:84] = b"PE\0\0"
            pe[84:86] = (0x8664).to_bytes(2, "little")
            pe_path = root / "pe.exe"
            pe_path.write_bytes(pe)
            self.assertEqual(release.binary_architecture(pe_path), "x86_64")

            macho = bytearray(32)
            macho[:4] = b"\xcf\xfa\xed\xfe"
            macho[4:8] = (0x0100000C).to_bytes(4, "little")
            macho_path = root / "macho"
            macho_path.write_bytes(macho)
            self.assertEqual(release.binary_architecture(macho_path), "arm64")

            macho_be = bytearray(32)
            macho_be[:4] = b"\xfe\xed\xfa\xcf"
            macho_be[4:8] = (0x01000007).to_bytes(4, "big")
            macho_be_path = root / "macho-big-endian"
            macho_be_path.write_bytes(macho_be)
            self.assertEqual(release.binary_architecture(macho_be_path), "x86_64")

    def test_smoke_uses_absolute_fixture_paths_for_wasi_preopen(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        contract["smoke"]["expected_output_sha256"] = hashlib.sha256(
            b"smoke-output\n"
        ).hexdigest()
        observed_command: list[str] = []

        def fake_run_capture(command, cwd, *, clean_losat_environment=False):
            del cwd, clean_losat_environment
            observed_command.extend(command)
            output = Path(command[command.index("--out") + 1])
            output.write_bytes(b"smoke-output\n")
            return subprocess.CompletedProcess(command, 0, "", "")

        with mock.patch.object(release, "run_capture", side_effect=fake_run_capture):
            release.run_smoke(["node", "runner.js", "LOSAT.wasm"], REPO_ROOT, contract)

        query = observed_command[observed_command.index("--query") + 1]
        subject = observed_command[observed_command.index("--subject") + 1]
        self.assertEqual(Path(query), (REPO_ROOT / contract["smoke"]["query"]).resolve())
        self.assertEqual(
            Path(subject), (REPO_ROOT / contract["smoke"]["subject"]).resolve()
        )

    def test_candidate_binary_hash_is_recorded_without_cross_runner_equality(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        target = "x86_64-unknown-linux-gnu"
        spec = release.artifact_spec(contract, "native", target)
        spec["historical_certification_binary_sha256"] = "f" * 64
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            binary = root / "LOSAT"
            image = bytearray(64)
            image[:6] = b"\x7fELF\x02\x01"
            image[18:20] = (62).to_bytes(2, "little")
            binary.write_bytes(image)
            output_dir = root / "output"
            args = argparse.Namespace(
                repo_root=REPO_ROOT,
                candidate_sha="a" * 40,
                target=target,
                binary=binary,
                output_dir=output_dir,
            )
            with (
                mock.patch.object(release, "load_contract", return_value=contract),
                mock.patch.object(release, "validate_candidate"),
                mock.patch.object(release, "capture_toolchain", return_value={}),
                mock.patch.object(release, "version_output", return_value="losat 0.1.0"),
                mock.patch.object(release, "run_smoke", return_value={"status": "PASS"}),
            ):
                release.assemble_binary(args, "native")

            metadata = json.loads(
                (output_dir / f"{spec['filename']}.metadata.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(metadata["candidate_sha"], args.candidate_sha)
            self.assertEqual(metadata["binary"]["sha256"], release.sha256_path(binary))
            self.assertFalse(
                metadata["certification_binary_reference"]["same_bytes_in_this_build"]
            )
            self.assertEqual(
                metadata["certification_binary_reference"]["sha256"], "f" * 64
            )

    def test_source_crate_rejects_generated_content(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            crate = Path(temp_dir) / "LOSAT-0.1.0.crate"
            with tarfile.open(crate, "w:gz") as handle:
                for name in (
                    "LOSAT-0.1.0/Cargo.toml",
                    "LOSAT-0.1.0/Cargo.lock",
                    "LOSAT-0.1.0/LICENSE",
                    "LOSAT-0.1.0/README.md",
                    "LOSAT-0.1.0/target/debug/LOSAT",
                ):
                    data = (
                        b'[package]\nname = "LOSAT"\nversion = "0.1.0"\nlicense = "MIT"\n'
                        if name.endswith("Cargo.toml")
                        else b"fixture"
                    )
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mode = 0o644
                    handle.addfile(info, io.BytesIO(data))
            with self.assertRaisesRegex(release.ReleaseFailure, "generated/scratch"):
                release.crate_members(crate, "LOSAT-0.1.0")

    def test_source_crate_requires_license_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            crate = Path(temp_dir) / "LOSAT-0.1.0.crate"
            with tarfile.open(crate, "w:gz") as handle:
                contents = {
                    "LOSAT-0.1.0/Cargo.toml": (
                        b'[package]\nname = "LOSAT"\nversion = "0.1.0"\n'
                    ),
                    "LOSAT-0.1.0/Cargo.lock": b"fixture",
                    "LOSAT-0.1.0/README.md": b"fixture",
                }
                for name, data in contents.items():
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mode = 0o644
                    handle.addfile(info, io.BytesIO(data))
            with self.assertRaisesRegex(release.ReleaseFailure, "license metadata"):
                release.crate_members(crate, "LOSAT-0.1.0")

    def test_aggregate_requires_and_hashes_the_exact_contract_set(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        candidate_sha = "a" * 40
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inputs = root / "inputs"
            outputs = root / "outputs"
            inputs.mkdir()
            for index, spec in enumerate(contract["artifacts"]):
                artifact = inputs / spec["filename"]
                artifact.write_bytes(f"artifact-{index}\n".encode())
                digest = release.sha256_path(artifact)
                (inputs / f"{artifact.name}.sha256").write_text(
                    f"{digest}  {artifact.name}\n", encoding="utf-8"
                )
                release.write_json(
                    inputs / f"{artifact.name}.metadata.json",
                    {
                        "candidate_sha": candidate_sha,
                        "integration_base_sha": contract["integration_base_sha"],
                        "artifact": {
                            "kind": spec["kind"],
                            "target": spec["target"],
                            "filename": spec["filename"],
                            **(
                                {"sha256": digest}
                                if spec["kind"] == "source"
                                else {}
                            ),
                        },
                        **(
                            {"artifact_sha256": digest}
                            if spec["kind"] != "source"
                            else {}
                        ),
                        "certification_lineage": contract["certification_lineage"],
                        "binary_hash_policy": contract["binary_hash_policy"],
                        "publication": contract["publication"],
                        "status": "PASS",
                    },
                )
            args = argparse.Namespace(
                repo_root=REPO_ROOT,
                candidate_sha=candidate_sha,
                input_root=inputs,
                output_dir=outputs,
                workflow_run_id="123456",
                integrated_certification_rerun="false",
            )
            with mock.patch.object(release, "validate_candidate"):
                release.aggregate(args)
            handoff = json.loads((outputs / "RC-HANDOFF.json").read_text(encoding="utf-8"))
            self.assertEqual(handoff["decision"], "RC_HANDOFF_READY")
            self.assertEqual(handoff["candidate_sha"], candidate_sha)
            self.assertEqual(len(handoff["artifacts"]), 6)
            self.assertFalse(handoff["integrated_certification_rerun"])
            self.assertEqual(len((outputs / "SHA256SUMS").read_text().splitlines()), 6)

    def test_release_workflow_is_manual_read_only_and_non_publishing(self) -> None:
        workflow = (REPO_ROOT / ".github/workflows/release-readiness.yml").read_text(
            encoding="utf-8"
        )
        self.assertIn("workflow_dispatch:", workflow)
        self.assertIn("contents: read", workflow)
        self.assertNotIn("actions/create-release", workflow)
        self.assertNotIn("softprops/action-gh-release", workflow)
        self.assertNotIn("cargo publish\n", workflow)
        self.assertNotIn("git tag", workflow)
        for target in (
            "x86_64-unknown-linux-gnu",
            "x86_64-pc-windows-msvc",
            "aarch64-apple-darwin",
            "x86_64-apple-darwin",
            "wasm32-wasip1",
        ):
            self.assertIn(target, workflow)
        self.assertIn("prepare_release_candidate_v010.py", workflow)


if __name__ == "__main__":
    unittest.main()
