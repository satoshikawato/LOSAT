#!/usr/bin/env python3
"""Focused tests for the v0.2.0 exact-SHA release-candidate assembler.

No network and no cargo: git is used only on temporary repositories.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import subprocess
import tarfile
import tempfile
import unittest
import zipfile
from copy import deepcopy
from pathlib import Path
from unittest import mock


MODULE_PATH = Path(__file__).with_name("prepare_release_candidate_v020.py")
SPEC = importlib.util.spec_from_file_location("prepare_release_candidate_v020", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
release = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release)

REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_FILE = REPO_ROOT / "docs/release/v0.2.0_rc_contract.json"

EXPECTED_ARTIFACTS = {
    ("native", "x86_64-unknown-linux-gnu", "LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz"),
    ("native", "aarch64-unknown-linux-gnu", "LOSAT-0.2.0-aarch64-unknown-linux-gnu.tar.gz"),
    ("native", "x86_64-pc-windows-msvc", "LOSAT-0.2.0-x86_64-pc-windows-msvc.zip"),
    ("native", "aarch64-apple-darwin", "LOSAT-0.2.0-aarch64-apple-darwin.tar.gz"),
    ("native", "x86_64-apple-darwin", "LOSAT-0.2.0-x86_64-apple-darwin.tar.gz"),
    ("wasm", "wasm32-wasip1", "LOSAT-0.2.0-wasm32-wasip1.tar.gz"),
    ("wasm", "wasm32-wasip1-threads", "LOSAT-0.2.0-wasm32-wasip1-threads.tar.gz"),
    ("source", "source", "LOSAT-0.2.0.crate"),
}


def committed_contract() -> dict:
    return json.loads(CONTRACT_FILE.read_text(encoding="utf-8"))


class ContractSchemaTests(unittest.TestCase):
    def test_committed_contract_passes_the_schema(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        self.assertEqual(contract["schema_version"], 2)

    def test_one_artifact_per_declared_target(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        observed = {(i["kind"], i["target"], i["filename"]) for i in contract["artifacts"]}
        self.assertEqual(observed, EXPECTED_ARTIFACTS)
        self.assertEqual(len(contract["artifacts"]), 8)

    def test_one_smoke_per_program_with_existing_fixtures(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        self.assertEqual(sorted(e["program"] for e in contract["smoke"]), sorted(release.PROGRAMS))
        for entry in contract["smoke"]:
            self.assertTrue((REPO_ROOT / entry["query"]).is_file(), entry["query"])
            self.assertTrue((REPO_ROOT / entry["subject"]).is_file(), entry["subject"])
            self.assertRegex(entry["expected_output_sha256"], r"^[0-9a-f]{64}$")
            for placeholder in release.PLACEHOLDERS:
                self.assertEqual(entry["arguments"].count(placeholder), 1)
            threads = entry["arguments"].index("-num_threads")
            self.assertEqual(entry["arguments"][threads + 1], "1")

    def test_wasm_entries_have_runner_scripts(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        runners = {
            i["target"]: i["runner_script"] for i in contract["artifacts"] if i["kind"] == "wasm"
        }
        self.assertEqual(
            runners,
            {
                "wasm32-wasip1": "LOSAT/tests/run_losat_wasi.js",
                "wasm32-wasip1-threads": "LOSAT/tests/run_losat_wasi_threads.js",
            },
        )
        for script in runners.values():
            self.assertTrue((REPO_ROOT / script).is_file())

    def test_unfilled_lineage_shape_and_publication_flags(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        self.assertFalse(any(contract["publication"].values()))
        self.assertFalse(contract["binary_hash_policy"]["cross_runner_binary_equality_required"])
        for item in contract["artifacts"]:
            historical = item["historical_certification_binary_sha256"]
            if item["target"] in ("aarch64-unknown-linux-gnu", "wasm32-wasip1-threads", "source"):
                self.assertIsNone(historical)

    def assertRejected(self, mutate, pattern: str) -> None:
        contract = committed_contract()
        mutate(contract)
        with self.assertRaisesRegex(release.ReleaseFailure, pattern):
            release.check_contract_schema(contract, REPO_ROOT)

    def test_schema_rejections(self) -> None:
        self.assertRejected(lambda c: c.update(schema_version=1), "schema_version")
        self.assertRejected(lambda c: c.update(release="v0.1.0"), "release")
        self.assertRejected(lambda c: c.update(package_version="0.1.0"), "package_version")
        self.assertRejected(lambda c: c["artifacts"].pop(), "eight artifacts")
        self.assertRejected(
            lambda c: c["artifacts"][1].update(filename=c["artifacts"][0]["filename"]), "unique"
        )
        self.assertRejected(lambda c: c["artifacts"][0].update(binary_name="LOSAT"), "binary_name")
        self.assertRejected(lambda c: c["artifacts"][5].update(runner_script="LOSAT/tests/none.js"), "runner_script")
        self.assertRejected(
            lambda c: c["artifacts"][0].update(historical_certification_binary_sha256="abc"), "historical"
        )
        self.assertRejected(
            lambda c: c["binary_hash_policy"].update(cross_runner_binary_equality_required=True), "equality"
        )
        self.assertRejected(lambda c: c["smoke"].pop(), "one entry per program")
        self.assertRejected(lambda c: c["smoke"][0].update(query="LOSAT/tests/fasta/none.fa"), "query")
        self.assertRejected(lambda c: c["smoke"][1].update(expected_output_sha256="x"), "expected_output_sha256")
        self.assertRejected(lambda c: c["smoke"][2]["arguments"].remove("{output}"), r"\{output\}")
        self.assertRejected(
            lambda c: c["smoke"][3]["arguments"].__setitem__(
                c["smoke"][3]["arguments"].index("-num_threads") + 1, "2"
            ),
            "num_threads",
        )
        self.assertRejected(lambda c: c["certification_lineage"].update(candidate_sha="abc"), "candidate_sha")
        self.assertRejected(
            lambda c: c["certification_lineage"]["runs"]["B"].update(run_id=0), "run_id"
        )
        self.assertRejected(
            lambda c: c["certification_lineage"]["runs"]["G"]["runs"].pop(), "required_consecutive_green"
        )
        self.assertRejected(
            lambda c: c["certification_lineage"]["local_evidence"]["A"].update(evidence_sha256="1"),
            "evidence_sha256",
        )
        self.assertRejected(lambda c: c["publication"].update(tag_created=True), "publication")


class BinaryAndArchiveTests(unittest.TestCase):
    def write(self, root: Path, name: str, data: bytes) -> Path:
        path = root / name
        path.write_bytes(data)
        return path

    def test_binary_architecture_recognises_every_format(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)

            def elf(machine: int) -> bytes:
                image = bytearray(64)
                image[:6] = b"\x7fELF\x02\x01"
                image[18:20] = machine.to_bytes(2, "little")
                return bytes(image)

            def pe(machine: int) -> bytes:
                image = bytearray(128)
                image[:2] = b"MZ"
                image[60:64] = (80).to_bytes(4, "little")
                image[80:84] = b"PE\0\0"
                image[84:86] = machine.to_bytes(2, "little")
                return bytes(image)

            def macho(cpu: int) -> bytes:
                image = bytearray(32)
                image[:4] = b"\xcf\xfa\xed\xfe"
                image[4:8] = cpu.to_bytes(4, "little")
                return bytes(image)

            cases = [
                ("elf-x64", elf(62), "x86_64", "elf"),
                ("elf-arm64", elf(183), "arm64", "elf"),
                ("pe-x64", pe(0x8664), "x86_64", "pe"),
                ("pe-arm64", pe(0xAA64), "arm64", "pe"),
                ("macho-x64", macho(0x01000007), "x86_64", "macho"),
                ("macho-arm64", macho(0x0100000C), "arm64", "macho"),
                ("module.wasm", b"\0asm\x01\0\0\0" + b"\0" * 8, "wasm32", "wasm"),
            ]
            for name, data, architecture, fmt in cases:
                path = self.write(root, name, data)
                self.assertEqual(release.binary_architecture(path), architecture, name)
                self.assertEqual(release.binary_format(path), fmt, name)
            with self.assertRaises(release.ReleaseFailure):
                release.binary_architecture(self.write(root, "text", b"#!/bin/sh\nexit 0\n" * 4))
            with self.assertRaisesRegex(release.ReleaseFailure, "module version"):
                release.binary_architecture(self.write(root, "old.wasm", b"\0asm\x02\0\0\0"))

    def test_expected_format_by_target(self) -> None:
        self.assertEqual(release.expected_format("x86_64-pc-windows-msvc"), "pe")
        self.assertEqual(release.expected_format("aarch64-apple-darwin"), "macho")
        self.assertEqual(release.expected_format("aarch64-unknown-linux-gnu"), "elf")
        with self.assertRaises(release.ReleaseFailure):
            release.expected_format("wasm32-wasip1")

    def test_tar_and_zip_assembly_are_deterministic(self) -> None:
        entries = [
            ("LOSAT-0.2.0-test/losat", b"binary\0bytes", 0o755),
            ("LOSAT-0.2.0-test/README.md", b"readme\n", 0o644),
        ]
        self.assertEqual(release.tar_bytes(entries), release.tar_bytes(reversed(entries)))
        self.assertEqual(release.zip_bytes(entries), release.zip_bytes(reversed(entries)))

    def test_archive_extraction_rejects_path_escape(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            tar_path = root / "unsafe.tar.gz"
            with tarfile.open(tar_path, "w:gz") as handle:
                info = tarfile.TarInfo("LOSAT-0.2.0-test/../outside")
                info.size = 6
                handle.addfile(info, io.BytesIO(b"unsafe"))
            with self.assertRaisesRegex(release.ReleaseFailure, "unsafe archive member"):
                release.extract_archive(tar_path, root / "x1", "LOSAT-0.2.0-test")
            zip_path = root / "unsafe.zip"
            with zipfile.ZipFile(zip_path, "w") as handle:
                handle.writestr("LOSAT-0.2.0-test/../outside", b"unsafe")
            with self.assertRaisesRegex(release.ReleaseFailure, "unsafe archive member"):
                release.extract_archive(zip_path, root / "x2", "LOSAT-0.2.0-test")

    def test_thread_variant_changes_only_the_thread_count(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        entry = next(e for e in contract["smoke"] if e["program"] == "blastp")
        changed = release.thread_variant(entry["arguments"], "2")
        self.assertEqual(changed[changed.index("-num_threads") + 1], "2")
        differing = [i for i, (a, b) in enumerate(zip(entry["arguments"], changed)) if a != b]
        self.assertEqual(len(differing), 1)


class SmokeTests(unittest.TestCase):
    def fake_contract(self, payload: bytes) -> dict:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        for entry in contract["smoke"]:
            entry["expected_output_sha256"] = hashlib.sha256(payload).hexdigest()
        return contract

    def test_every_smoke_runs_with_absolute_paths_and_empty_console(self) -> None:
        contract = self.fake_contract(b"smoke-output\n")
        commands: list[list[str]] = []

        def fake_run_capture(command, cwd, *, clean_losat_environment=False):
            del cwd
            self.assertTrue(clean_losat_environment)
            commands.append(list(command))
            Path(command[command.index("-out") + 1]).write_bytes(b"smoke-output\n")
            return subprocess.CompletedProcess(command, 0, "", "")

        with mock.patch.object(release, "run_capture", side_effect=fake_run_capture):
            results = release.run_smokes(["node", "runner.js", "losat.wasm"], REPO_ROOT, contract)
        self.assertEqual([r["program"] for r in results], [e["program"] for e in contract["smoke"]])
        self.assertEqual(len(commands), 5)
        for command in commands:
            for flag in ("-query", "-subject"):
                self.assertTrue(Path(command[command.index(flag) + 1]).is_absolute())
            self.assertEqual(command[command.index("-num_threads") + 1], "1")

    def test_threaded_artifact_adds_one_two_thread_blastp(self) -> None:
        contract = self.fake_contract(b"smoke-output\n")
        commands: list[list[str]] = []

        def fake_run_capture(command, cwd, *, clean_losat_environment=False):
            del cwd, clean_losat_environment
            commands.append(list(command))
            Path(command[command.index("-out") + 1]).write_bytes(b"smoke-output\n")
            return subprocess.CompletedProcess(command, 0, "", "")

        with mock.patch.object(release, "run_capture", side_effect=fake_run_capture):
            results = release.run_smokes(["losat"], REPO_ROOT, contract, threaded=True)
        self.assertEqual(len(results), 6)
        self.assertEqual(results[-1]["program"], "blastp")
        self.assertEqual(results[-1]["num_threads"], 2)
        self.assertEqual(sum(1 for r in results if r["num_threads"] == 2), 1)
        self.assertEqual(commands[-1][commands[-1].index("-num_threads") + 1], "2")

    def test_smoke_failures(self) -> None:
        contract = self.fake_contract(b"expected\n")
        for stdout, payload, pattern in [
            ("noise", b"expected\n", "console output"),
            ("", b"different\n", "hash mismatch"),
        ]:
            def fake_run_capture(command, cwd, *, clean_losat_environment=False, _o=stdout, _p=payload):
                del cwd, clean_losat_environment
                Path(command[command.index("-out") + 1]).write_bytes(_p)
                return subprocess.CompletedProcess(command, 0, _o, "")

            with mock.patch.object(release, "run_capture", side_effect=fake_run_capture):
                with self.assertRaisesRegex(release.ReleaseFailure, pattern):
                    release.run_smokes(["losat"], REPO_ROOT, contract)

        def no_output(command, cwd, *, clean_losat_environment=False):
            return subprocess.CompletedProcess(command, 0, "", "")

        with mock.patch.object(release, "run_capture", side_effect=no_output):
            with self.assertRaisesRegex(release.ReleaseFailure, "did not create"):
                release.run_smokes(["losat"], REPO_ROOT, contract)


def git(root: Path, *arguments: str) -> str:
    environment = dict(
        os.environ,
        GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@example.invalid",
        GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@example.invalid",
        GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_SYSTEM=os.devnull,
    )
    result = subprocess.run(
        ["git", "-c", "commit.gpgsign=false", *arguments],
        cwd=root, env=environment, capture_output=True, text=True, check=False,
    )
    if result.returncode:
        raise AssertionError(f"git {arguments} failed: {result.stderr}")
    return result.stdout.strip()


def put(root: Path, relative: str, text: str) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def commit(root: Path, message: str) -> str:
    git(root, "add", "--all", ".")
    git(root, "commit", "--quiet", "--allow-empty", "-m", message)
    return git(root, "rev-parse", "HEAD")


class FakeRepository:
    """A temporary git repository holding a copy of the contract and its fixtures.

    Commit 0 stands for the v0.1.0 ancestors, commit C holds the runtime inputs.
    """

    def __init__(self, root: Path) -> None:
        self.root = root
        git(root, "init", "--quiet", "-b", "main")
        put(root, "LICENSE", "license\n")
        put(root, "README.md", "readme\n")
        put(root, "LOSAT/Cargo.toml", '[package]\nname = "LOSAT"\nversion = "0.2.0"\n')
        self.ancestor = commit(root, "ancestor")
        self.contract = committed_contract()
        for key in release.V010_ANCESTOR_KEYS:
            self.contract["certification_lineage"]["v010_ancestors"][key] = self.ancestor
        for entry in self.contract["smoke"]:
            for key in ("query", "subject"):
                put(root, entry[key], ">s\nACGT\n")
        for item in self.contract["artifacts"]:
            if "runner_script" in item:
                put(root, item["runner_script"], "// runner\n")
        put(root, "LOSAT/src/main.rs", "fn main() {}\n")
        put(root, "LOSAT/tests/prepare_release_candidate_v020.py", "# tool\n")
        put(root, "LOSAT/tests/test_prepare_release_candidate_v020.py", "# test\n")
        put(root, "docs/release/notes.md", "notes\n")
        self.write_contract()
        self.certified = commit(root, "certified runtime")

    def write_contract(self) -> None:
        put(self.root, str(release.CONTRACT_PATH), json.dumps(self.contract, indent=2) + "\n")

    def fill(self) -> None:
        lineage = self.contract["certification_lineage"]
        lineage["candidate_sha"] = self.certified
        for key in release.REQUIRED_SUCCESS_KEYS:
            lineage["runs"][key].update(run_id=1000, conclusion="success")
        for index, run in enumerate(lineage["runs"]["G"]["runs"]):
            run.update(run_id=2000 + index, conclusion="success")
        for key in release.LOCAL_EVIDENCE_KEYS:
            lineage["local_evidence"][key]["evidence_sha256"] = "a" * 64
        self.write_contract()
        commit(self.root, "fill the contract")

    def validate(self, *, allow_unfilled: bool = False) -> None:
        contract = release.load_contract(self.root)
        release.validate_candidate(
            self.root, git(self.root, "rev-parse", "HEAD"), contract,
            allow_unfilled_lineage=allow_unfilled,
        )


class CandidateValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.repo = FakeRepository(Path(self.temp.name))

    def test_unfilled_candidate_sha_needs_the_flag(self) -> None:
        commit(self.repo.root, "head")
        with self.assertRaisesRegex(release.ReleaseFailure, "candidate_sha is not filled"):
            self.repo.validate()
        self.repo.validate(allow_unfilled=True)

    def test_filled_contract_passes_and_docs_only_change_after_c_passes(self) -> None:
        self.repo.fill()
        put(self.repo.root, "docs/release/v0.2.0.md", "go record\n")
        put(self.repo.root, "README.md", "readme, edited\n")
        put(self.repo.root, "changelog.d/x.md", "- entry\n")
        commit(self.repo.root, "release documents")
        self.repo.validate()

    def test_runtime_source_change_after_c_fails_even_with_the_flag(self) -> None:
        self.repo.fill()
        put(self.repo.root, "LOSAT/src/main.rs", "fn main() { println!(); }\n")
        commit(self.repo.root, "runtime change")
        for flag in (False, True):
            with self.assertRaisesRegex(release.ReleaseFailure, "changes a runtime input"):
                self.repo.validate(allow_unfilled=flag)

    def test_other_runtime_inputs_fail(self) -> None:
        self.repo.fill()
        for relative in (
            "LOSAT/tests/fasta/new.fa", "LOSAT/Cargo.lock", "LOSAT/build.rs",
            "LOSAT/.cargo/config.toml", "docs/product_decisions/PD-X.md",
        ):
            with self.subTest(relative):
                put(self.repo.root, relative, "changed\n")
                commit(self.repo.root, relative)
                with self.assertRaisesRegex(release.ReleaseFailure, "changes a runtime input"):
                    self.repo.validate()
                git(self.repo.root, "reset", "--quiet", "--hard", "HEAD~1")

    def test_release_tooling_change_is_excluded(self) -> None:
        self.repo.fill()
        put(self.repo.root, "LOSAT/tests/prepare_release_candidate_v020.py", "# tool, edited\n")
        put(self.repo.root, "LOSAT/tests/test_prepare_release_candidate_v020.py", "# test, edited\n")
        commit(self.repo.root, "tooling")
        self.repo.validate()

    def test_head_must_equal_the_candidate_sha(self) -> None:
        self.repo.fill()
        contract = release.load_contract(self.repo.root)
        with self.assertRaisesRegex(release.ReleaseFailure, "candidate SHA mismatch"):
            release.validate_candidate(self.repo.root, self.repo.certified, contract)
        with self.assertRaisesRegex(release.ReleaseFailure, "40 lowercase"):
            release.validate_candidate(self.repo.root, "abc", contract)

    def test_certified_commit_must_be_an_ancestor(self) -> None:
        git(self.repo.root, "switch", "--quiet", "-c", "side", self.repo.ancestor)
        put(self.repo.root, "side.txt", "side\n")
        side = commit(self.repo.root, "side")
        git(self.repo.root, "switch", "--quiet", "main")
        self.repo.certified = side
        self.repo.fill()
        with self.assertRaisesRegex(release.ReleaseFailure, "descend from the certified commit"):
            self.repo.validate()

    def test_v010_ancestors_must_be_ancestors(self) -> None:
        git(self.repo.root, "switch", "--quiet", "-c", "side", self.repo.ancestor)
        put(self.repo.root, "side.txt", "side\n")
        side = commit(self.repo.root, "side")
        git(self.repo.root, "switch", "--quiet", "main")
        self.repo.contract["certification_lineage"]["v010_ancestors"]["integrated_runtime"] = side
        self.repo.fill()
        for flag in (False, True):
            with self.assertRaisesRegex(release.ReleaseFailure, "v0.1.0 ancestor integrated_runtime"):
                self.repo.validate(allow_unfilled=flag)

    def test_package_version_must_match(self) -> None:
        put(self.repo.root, "LOSAT/Cargo.toml", '[package]\nname = "LOSAT"\nversion = "0.2.1"\n')
        # Cargo.toml is a runtime input: certify the wrong version so only the version rule fires.
        self.repo.certified = commit(self.repo.root, "version bump")
        self.repo.fill()
        with self.assertRaisesRegex(release.ReleaseFailure, "package version"):
            self.repo.validate()

    def test_runs_and_evidence_must_be_filled_and_green(self) -> None:
        self.repo.fill()
        lineage = self.repo.contract["certification_lineage"]
        mutations = [
            (lambda: lineage["runs"]["B"].update(run_id=None), "run B has no run_id"),
            (lambda: lineage["runs"]["C"].update(conclusion="failure"), "run C conclusion"),
            (lambda: lineage["runs"]["F"].update(conclusion=None), "run F conclusion"),
            (lambda: lineage["runs"]["H"].update(run_id=None), "run H has no run_id"),
            (lambda: lineage["runs"]["G"]["runs"][2].update(conclusion="cancelled"), r"G\[2\]"),
            (lambda: lineage["runs"]["G"]["runs"][0].update(run_id=None), r"G\[0\]"),
            (lambda: lineage["local_evidence"]["D"].update(evidence_sha256=None), "evidence D"),
        ]
        for mutate, pattern in mutations:
            with self.subTest(pattern):
                saved = deepcopy(self.repo.contract)
                mutate()
                self.repo.write_contract()
                commit(self.repo.root, pattern)
                with self.assertRaisesRegex(release.ReleaseFailure, pattern):
                    self.repo.validate()
                self.repo.validate(allow_unfilled=True)
                self.repo.contract.clear()
                self.repo.contract.update(saved)
                lineage = self.repo.contract["certification_lineage"]
                self.repo.write_contract()
                commit(self.repo.root, "restore " + pattern)

    def test_run_i_may_be_unset(self) -> None:
        self.repo.fill()
        self.assertIsNone(self.repo.contract["certification_lineage"]["runs"]["I"]["run_id"])
        self.repo.validate()

    def test_cli_validate_reports_a_failure(self) -> None:
        commit(self.repo.root, "head")
        head = git(self.repo.root, "rev-parse", "HEAD")
        argv = ["validate", "--repo-root", str(self.repo.root), "--candidate-sha", head]
        with mock.patch("sys.stderr", new_callable=io.StringIO) as err:
            self.assertEqual(release.main(argv), 1)
        self.assertIn("RC_HANDOFF_FAILED", err.getvalue())
        with mock.patch("sys.stdout", new_callable=io.StringIO) as out:
            self.assertEqual(release.main([*argv, "--allow-unfilled-lineage"]), 0)
        self.assertIn("RC_CANDIDATE_OK", out.getvalue())
        with mock.patch("sys.stdout", new_callable=io.StringIO) as out:
            self.assertEqual(release.main(["check-contract", "--repo-root", str(REPO_ROOT)]), 0)


class AssembleTests(unittest.TestCase):
    def test_native_assembly_records_metadata_and_flags(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        target = "x86_64-unknown-linux-gnu"
        spec = release.artifact_spec(contract, "native", target)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            binary = root / "losat"
            image = bytearray(64)
            image[:6] = b"\x7fELF\x02\x01"
            image[18:20] = (62).to_bytes(2, "little")
            binary.write_bytes(image)
            args = argparse.Namespace(
                repo_root=REPO_ROOT, candidate_sha="a" * 40, target=target, binary=binary,
                output_dir=root / "out", allow_unfilled_lineage=True,
            )
            with (
                mock.patch.object(release, "load_contract", return_value=contract),
                mock.patch.object(release, "validate_candidate") as validate,
                mock.patch.object(release, "capture_toolchain", return_value={}),
                mock.patch.object(release, "version_output", return_value="losat 0.2.0"),
                mock.patch.object(release, "run_smokes", return_value=[{"status": "PASS"}]),
            ):
                with contextlib.redirect_stdout(io.StringIO()):
                    release.assemble_binary(args, "native")
            self.assertTrue(validate.call_args.kwargs["allow_unfilled_lineage"])
            metadata = json.loads((root / "out" / f"{spec['filename']}.metadata.json").read_text())
            self.assertEqual(metadata["schema_version"], 2)
            self.assertEqual(metadata["binary"]["sha256"], release.sha256_path(binary))
            self.assertEqual(metadata["binary"]["architecture"], "x86_64")
            self.assertEqual(metadata["lineage_validation"],
                             {"allow_unfilled_lineage": True, "allow_node_difference": False})
            self.assertEqual(metadata["status"], "PASS")
            with tarfile.open(root / "out" / spec["filename"]) as archive:
                names = sorted(archive.getnames())
            root_name = spec["filename"].removesuffix(".tar.gz")
            self.assertEqual(
                names, sorted(f"{root_name}/{n}" for n in ("LICENSE", "README.md", "RELEASE-METADATA.json", "losat"))
            )

    def test_native_rejects_wrong_architecture_and_format(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            binary = root / "losat"
            image = bytearray(64)
            image[:6] = b"\x7fELF\x02\x01"
            image[18:20] = (183).to_bytes(2, "little")
            binary.write_bytes(image)

            def run(target: str) -> None:
                args = argparse.Namespace(
                    repo_root=REPO_ROOT, candidate_sha="a" * 40, target=target, binary=binary,
                    output_dir=root / "out", allow_unfilled_lineage=True,
                )
                with (
                    mock.patch.object(release, "load_contract", return_value=contract),
                    mock.patch.object(release, "validate_candidate"),
                    mock.patch.object(release, "capture_toolchain", return_value={}),
                ):
                    release.assemble_binary(args, "native")

            with self.assertRaisesRegex(release.ReleaseFailure, "architecture mismatch"):
                run("x86_64-unknown-linux-gnu")
            with self.assertRaisesRegex(release.ReleaseFailure, "format mismatch"):
                run("aarch64-apple-darwin")

    def test_wasm_runner_must_be_the_contract_runner(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            binary = root / "losat.wasm"
            binary.write_bytes(b"\0asm\x01\0\0\0")
            args = argparse.Namespace(
                repo_root=REPO_ROOT, candidate_sha="a" * 40, target="wasm32-wasip1-threads",
                binary=binary, node="node", runner=REPO_ROOT / "LOSAT/tests/run_losat_wasi.js",
                output_dir=root / "out", allow_unfilled_lineage=True, allow_node_difference=True,
            )
            with (
                mock.patch.object(release, "load_contract", return_value=contract),
                mock.patch.object(release, "validate_candidate"),
                mock.patch.object(release, "capture_wasm_toolchain", return_value={}),
            ):
                with self.assertRaisesRegex(release.ReleaseFailure, "runner_script"):
                    release.assemble_binary(args, "wasm")

    def test_wasm_node_difference_needs_the_flag(self) -> None:
        contract = release.load_contract(REPO_ROOT)

        def fake(command, cwd, *, clean_losat_environment=False):
            if command[0] == "rustc":
                return subprocess.CompletedProcess(command, 0, "rustc 1.92.0 (x)\nhost: x86_64-unknown-linux-gnu\n", "")
            if command[0] == "cargo":
                return subprocess.CompletedProcess(command, 0, "cargo 1.92.0 (x)\n", "")
            return subprocess.CompletedProcess(command, 0, "v99.0.0\n", "")

        environment = {key: "" for key in release.OUTPUT_ENV_KEYS}
        with mock.patch.dict(os.environ, environment), mock.patch.object(
            release._v010, "run_capture", side_effect=fake
        ), mock.patch.object(release, "run_capture", side_effect=fake):
            with self.assertRaisesRegex(release.ReleaseFailure, "Node differs"):
                release.capture_wasm_toolchain(
                    REPO_ROOT, "wasm32-wasip1", contract, "node", allow_node_difference=False
                )
            toolchain = release.capture_wasm_toolchain(
                REPO_ROOT, "wasm32-wasip1", contract, "node", allow_node_difference=True
            )
        self.assertEqual(toolchain["node_V"], "v99.0.0")


class AggregateTests(unittest.TestCase):
    def make_inputs(self, root: Path, contract: dict, candidate_sha: str, *, flags=None, skip=None) -> Path:
        inputs = root / "inputs"
        inputs.mkdir()
        for index, spec in enumerate(contract["artifacts"]):
            if spec["filename"] == skip:
                continue
            artifact = inputs / "nested" / spec["filename"]
            artifact.parent.mkdir(exist_ok=True)
            artifact.write_bytes(f"artifact-{index}\n".encode())
            digest = release.sha256_path(artifact)
            (inputs / f"{artifact.name}.sha256").write_text(f"{digest}  {artifact.name}\n", encoding="utf-8")
            release.write_json(
                inputs / f"{artifact.name}.metadata.json",
                {
                    "candidate_sha": candidate_sha,
                    "release": contract["release"],
                    "integration_base_sha": contract["integration_base_sha"],
                    "artifact": {
                        "kind": spec["kind"], "target": spec["target"], "filename": spec["filename"],
                        **({"sha256": digest} if spec["kind"] == "source" else {}),
                    },
                    **({"artifact_sha256": digest} if spec["kind"] != "source" else {}),
                    "certification_lineage": contract["certification_lineage"],
                    "binary_hash_policy": contract["binary_hash_policy"],
                    "publication": contract["publication"],
                    "lineage_validation": flags
                    or {"allow_unfilled_lineage": False, "allow_node_difference": False},
                    "status": "PASS",
                },
            )
        return inputs

    def run_aggregate(self, inputs: Path, outputs: Path, candidate_sha: str, contract: dict, **extra) -> None:
        args = argparse.Namespace(
            repo_root=REPO_ROOT, candidate_sha=candidate_sha, input_root=inputs, output_dir=outputs,
            workflow_run_id="123456", allow_unfilled_lineage=False, **extra,
        )
        with (
            mock.patch.object(release, "load_contract", return_value=contract),
            mock.patch.object(release, "validate_candidate"),
        ):
            with contextlib.redirect_stdout(io.StringIO()):
                release.aggregate(args)

    def test_requires_and_hashes_exactly_the_eight_artifacts(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        contract["certification_lineage"]["candidate_sha"] = "c" * 40
        candidate_sha = "a" * 40
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inputs = self.make_inputs(root, contract, candidate_sha)
            outputs = root / "outputs"
            self.run_aggregate(inputs, outputs, candidate_sha, contract)
            handoff = json.loads((outputs / "RC-HANDOFF.json").read_text(encoding="utf-8"))
            self.assertEqual(handoff["decision"], "RC_HANDOFF_READY")
            self.assertEqual(len(handoff["artifacts"]), 8)
            sums = (outputs / "SHA256SUMS").read_text().splitlines()
            self.assertEqual(len(sums), 8)
            for line in sums:
                digest, name = line.split("  ")
                self.assertEqual(digest, release.sha256_path(outputs / name))
            self.assertFalse(any(handoff["publication"].values()))
            markdown = (outputs / "RC-HANDOFF.md").read_text(encoding="utf-8")
            self.assertTrue(markdown.startswith("# LOSAT v0.2.0 release-candidate handoff"))
            self.assertIn(candidate_sha, markdown)
            self.assertIn("123456", markdown)
            self.assertIn("c" * 40, markdown)
            for _, _, filename in EXPECTED_ARTIFACTS:
                self.assertIn(filename, markdown)
            self.assertIn("tag_created: false", markdown)

    def test_missing_extra_and_tampered_inputs_fail(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        candidate_sha = "a" * 40
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inputs = self.make_inputs(root, contract, candidate_sha, skip="LOSAT-0.2.0-wasm32-wasip1.tar.gz")
            with self.assertRaisesRegex(release.ReleaseFailure, "exactly one downloaded"):
                self.run_aggregate(inputs, root / "o1", candidate_sha, contract)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inputs = self.make_inputs(root, contract, candidate_sha)
            (inputs / "LOSAT-0.2.0-extra.tar.gz").write_bytes(b"x")
            with self.assertRaisesRegex(release.ReleaseFailure, "unexpected artifact files"):
                self.run_aggregate(inputs, root / "o2", candidate_sha, contract)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inputs = self.make_inputs(root, contract, candidate_sha)
            (inputs / "nested/LOSAT-0.2.0.crate").write_bytes(b"tampered")
            with self.assertRaisesRegex(release.ReleaseFailure, "checksum sidecar mismatch"):
                self.run_aggregate(inputs, root / "o3", candidate_sha, contract)

    def test_dry_run_artifacts_are_not_releasable(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        candidate_sha = "a" * 40
        for flags in (
            {"allow_unfilled_lineage": True, "allow_node_difference": False},
            {"allow_unfilled_lineage": False, "allow_node_difference": True},
        ):
            with tempfile.TemporaryDirectory() as temp_dir:
                root = Path(temp_dir)
                inputs = self.make_inputs(root, contract, candidate_sha, flags=flags)
                with self.assertRaisesRegex(release.ReleaseFailure, "dry-run flag"):
                    self.run_aggregate(inputs, root / "out", candidate_sha, contract)

    def test_candidate_sha_mismatch_in_metadata_fails(self) -> None:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inputs = self.make_inputs(root, contract, "b" * 40)
            with self.assertRaisesRegex(release.ReleaseFailure, "identity/status mismatch"):
                self.run_aggregate(inputs, root / "out", "a" * 40, contract)


class CliTests(unittest.TestCase):
    def test_parser_has_no_flag_defaults_that_relax_the_gate(self) -> None:
        parser = release.build_parser()
        base = ["--repo-root", ".", "--candidate-sha", "a" * 40, "--output-dir", "o"]
        native = parser.parse_args(["native", *base, "--target", "t", "--binary", "b"])
        self.assertFalse(native.allow_unfilled_lineage)
        wasm = parser.parse_args(["wasm", *base, "--binary", "b", "--runner", "r"])
        self.assertFalse(wasm.allow_unfilled_lineage)
        self.assertFalse(wasm.allow_node_difference)
        self.assertEqual(wasm.target, "wasm32-wasip1")
        aggregate = parser.parse_args(
            ["aggregate", *base, "--input-root", "i", "--workflow-run-id", "1"]
        )
        self.assertFalse(aggregate.allow_unfilled_lineage)

    def test_module_docstring_forbids_the_flags_in_the_workflow(self) -> None:
        self.assertIn("NEVER PASS", release.__doc__)


if __name__ == "__main__":
    unittest.main()
