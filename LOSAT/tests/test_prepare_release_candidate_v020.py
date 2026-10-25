#!/usr/bin/env python3
"""Focused tests for the v0.2.0 exact-SHA release-candidate handoff."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import os
import subprocess
import tarfile
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from unittest import mock


MODULE_PATH = Path(__file__).with_name("prepare_release_candidate_v020.py")
SPEC = importlib.util.spec_from_file_location("prepare_release_candidate_v020", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
release = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release)

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_PATH = REPO_ROOT / ".github/workflows/release-readiness-v020.yml"
CONTRACT_REL = release.CONTRACT_PATH

TARGETS = {
    ("native", "x86_64-unknown-linux-gnu", "LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz"),
    ("native", "aarch64-unknown-linux-gnu", "LOSAT-0.2.0-aarch64-unknown-linux-gnu.tar.gz"),
    ("native", "x86_64-pc-windows-msvc", "LOSAT-0.2.0-x86_64-pc-windows-msvc.zip"),
    ("native", "aarch64-apple-darwin", "LOSAT-0.2.0-aarch64-apple-darwin.tar.gz"),
    ("native", "x86_64-apple-darwin", "LOSAT-0.2.0-x86_64-apple-darwin.tar.gz"),
    ("wasm", "wasm32-wasip1", "LOSAT-0.2.0-wasm32-wasip1.tar.gz"),
    ("wasm", "wasm32-wasip1-threads", "LOSAT-0.2.0-wasm32-wasip1-threads.tar.gz"),
    ("source", "source", "LOSAT-0.2.0.crate"),
}


def raw_contract() -> dict:
    return json.loads((REPO_ROOT / CONTRACT_REL).read_text(encoding="utf-8"))


def filled_contract(candidate_sha: str) -> dict:
    """The repository contract with every lineage field filled, as on the cut day."""
    contract = raw_contract()
    lineage = contract["certification_lineage"]
    lineage["candidate_sha"] = candidate_sha
    for label in ("B", "C", "F", "H"):
        lineage["runs"][label].update(run_id=1000 + ord(label), conclusion="success")
    for index, run in enumerate(lineage["runs"]["G"]["runs"]):
        run.update(run_id=2000 + index, conclusion="success")
    for label, item in lineage["local_evidence"].items():
        item["evidence_sha256"] = hashlib.sha256(label.encode()).hexdigest()
    for item in lineage["ws_evidence"]["J"].values():
        item["record"] = "docs/evidence/release_v020/ws/record.md"
    return contract


def write_contract(root: Path, contract: dict) -> None:
    path = root / CONTRACT_REL
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(contract), encoding="utf-8")


def git(root: Path, *args: str) -> str:
    environment = {
        **os.environ,
        "GIT_AUTHOR_NAME": "t",
        "GIT_AUTHOR_EMAIL": "t@example.invalid",
        "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@example.invalid",
    }
    result = subprocess.run(
        ["git", *args], cwd=root, env=environment, capture_output=True, text=True, check=True
    )
    return result.stdout.strip()


def commit_file(root: Path, relative: str, text: str) -> str:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    git(root, "add", relative)
    git(root, "commit", "-q", "-m", f"edit {relative}")
    return git(root, "rev-parse", "HEAD")


class TempRepo:
    """A throw-away git repository with the runtime-input layout of the contract."""

    def __init__(self) -> None:
        self._temp = tempfile.TemporaryDirectory()
        self.root = Path(self._temp.name)
        git(self.root, "init", "-q")
        commit_file(
            self.root,
            "LOSAT/Cargo.toml",
            '[package]\nname = "LOSAT"\nversion = "0.2.0"\n',
        )
        commit_file(self.root, "LOSAT/src/lib.rs", "// runtime\n")
        commit_file(self.root, "LOSAT/tests/fixture.txt", "fixture\n")
        commit_file(self.root, "README.md", "readme\n")
        self.base = git(self.root, "rev-parse", "HEAD")

    def close(self) -> None:
        self._temp.cleanup()

    def contract(self, candidate_sha: str) -> dict:
        contract = filled_contract(candidate_sha)
        for key in contract["certification_lineage"]["v010_ancestors"]:
            if key != "note":
                contract["certification_lineage"]["v010_ancestors"][key] = self.base
        return contract

    def head(self) -> str:
        return git(self.root, "rev-parse", "HEAD")


class ContractTests(unittest.TestCase):
    def test_contract_has_one_exact_artifact_per_supported_target(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        observed = {
            (item["kind"], item["target"], item["filename"]) for item in contract["artifacts"]
        }
        self.assertEqual(observed, TARGETS)
        self.assertEqual(len(contract["artifacts"]), 8)
        filenames = [item["filename"] for item in contract["artifacts"]]
        self.assertEqual(len(filenames), len(set(filenames)))
        self.assertFalse(any(contract["publication"].values()))
        self.assertFalse(
            contract["binary_hash_policy"]["cross_runner_binary_equality_required"]
        )
        binary_names = {
            item["target"]: item["binary_name"]
            for item in contract["artifacts"]
            if item["kind"] != "source"
        }
        self.assertEqual(binary_names["x86_64-pc-windows-msvc"], "losat.exe")
        self.assertEqual(binary_names["wasm32-wasip1"], "losat.wasm")
        self.assertEqual(binary_names["wasm32-wasip1-threads"], "losat.wasm")
        self.assertEqual(binary_names["aarch64-apple-darwin"], "losat")
        for item in contract["artifacts"]:
            historical = item["historical_certification_binary_sha256"]
            if historical is not None:
                self.assertRegex(historical, r"^[0-9a-f]{64}$")
            if item["kind"] == "wasm":
                self.assertTrue((REPO_ROOT / item["runner_script"]).is_file())

    def test_wasm_runners_are_serial_and_threaded(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        runners = {
            item["target"]: item["runner_script"]
            for item in contract["artifacts"]
            if item["kind"] == "wasm"
        }
        self.assertEqual(
            runners,
            {
                "wasm32-wasip1": "LOSAT/tests/run_losat_wasi.js",
                "wasm32-wasip1-threads": "LOSAT/tests/run_losat_wasi_threads.js",
            },
        )
        self.assertEqual(
            release.WASM_KIND_BY_TARGET,
            {
                "wasm32-wasip1": "serial-command",
                "wasm32-wasip1-threads": "threaded-command",
            },
        )

    def test_smoke_has_one_case_per_program_with_placeholders(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        smoke = contract["smoke"]
        self.assertEqual(len(smoke), 5)
        self.assertEqual(
            {case["program"] for case in smoke},
            {"blastn", "blastp", "blastx", "tblastn", "tblastx"},
        )
        for case in smoke:
            arguments = case["arguments"]
            self.assertEqual(arguments[0], case["program"])
            for placeholder in ("{query}", "{subject}", "{output}"):
                self.assertEqual(arguments.count(placeholder), 1)
            self.assertEqual(arguments[arguments.index("-num_threads") + 1], "1")
            self.assertRegex(case["expected_output_sha256"], r"^[0-9a-f]{64}$")
            self.assertTrue((REPO_ROOT / case["query"]).is_file(), case["query"])
            self.assertTrue((REPO_ROOT / case["subject"]).is_file(), case["subject"])
            # outfmt 6 or 8 only: outfmt 0 and 7 print the -subject path.
            self.assertIn(arguments[arguments.index("-outfmt") + 1], {"6", "8"})

    def test_contract_lists_rejected_tasks_and_toolchain(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        self.assertEqual(
            contract["public_cli"]["rejected_tasks"], ["blastp-short", "tblastn-fast"]
        )
        for case in contract["smoke"]:
            self.assertNotIn("blastp-short", case["arguments"])
            self.assertNotIn("tblastn-fast", case["arguments"])
        toolchain = contract["cert_toolchain"]
        self.assertEqual(toolchain["rust"], "1.92.0")
        self.assertEqual(toolchain["cargo"], "1.92.0")
        self.assertEqual(toolchain["node"], "24.21.0")
        self.assertEqual(toolchain["profile"], "release")
        self.assertTrue(toolchain["locked"])
        self.assertEqual(set(toolchain["output_affecting_environment"]), set(release.OUTPUT_ENV_KEYS))
        self.assertTrue(all(v == "" for v in toolchain["output_affecting_environment"].values()))

    def test_lineage_shape_is_fixed_but_unfilled_in_the_repository(self) -> None:
        lineage = release.load_contract(REPO_ROOT)["certification_lineage"]
        self.assertEqual(sorted(lineage["runs"]), ["B", "C", "F", "G", "H", "I"])
        self.assertEqual(lineage["runs"]["G"]["required_consecutive_green"], 3)
        self.assertEqual(len(lineage["runs"]["G"]["runs"]), 3)
        self.assertEqual(sorted(lineage["local_evidence"]), ["A", "D", "E"])
        self.assertEqual(lineage["candidate_sha"], "")
        self.assertTrue(lineage["rerun_required"])
        for name, sha in lineage["v010_ancestors"].items():
            if name != "note":
                self.assertRegex(sha, r"^[0-9a-f]{40}$")
        self.assertIn("LOSAT/src", lineage["runtime_input_paths"])
        self.assertIn(
            "LOSAT/tests/prepare_release_candidate_v020.py", lineage["runtime_input_excludes"]
        )

    def test_platform_claims_reference_contract_artifacts(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        claims = contract["platform_claims"]
        for claim in claims:
            index = claim["artifact_index"]
            if index is not None:
                self.assertTrue(1 <= index <= 7)
            self.assertTrue(claim["wording"])
        self.assertEqual(
            {claim["claim"] for claim in claims}, {"certified", "supported", "preview"}
        )

    def test_load_contract_rejects_changed_contracts(self) -> None:
        def mutate_release(c: dict) -> None:
            c["release"] = "v0.1.0"

        def mutate_version(c: dict) -> None:
            c["package_version"] = "0.1.0"

        def mutate_schema(c: dict) -> None:
            c["schema_version"] = 1

        def drop_artifact(c: dict) -> None:
            c["artifacts"].pop()

        def duplicate_filename(c: dict) -> None:
            c["artifacts"][1]["filename"] = c["artifacts"][0]["filename"]

        def swap_target(c: dict) -> None:
            c["artifacts"][1]["target"] = "x86_64-unknown-linux-gnu"

        def bad_historical(c: dict) -> None:
            c["artifacts"][0]["historical_certification_binary_sha256"] = "abc"

        def missing_runner(c: dict) -> None:
            del c["artifacts"][5]["runner_script"]

        def equality_required(c: dict) -> None:
            c["binary_hash_policy"]["cross_runner_binary_equality_required"] = True

        def drop_smoke(c: dict) -> None:
            c["smoke"].pop()

        def duplicate_program(c: dict) -> None:
            c["smoke"][1]["program"] = "blastn"
            c["smoke"][1]["arguments"][0] = "blastn"

        def missing_placeholder(c: dict) -> None:
            c["smoke"][2]["arguments"].remove("{output}")

        def many_threads(c: dict) -> None:
            arguments = c["smoke"][0]["arguments"]
            arguments[arguments.index("-num_threads") + 1] = "4"

        def published(c: dict) -> None:
            c["publication"]["tag_created"] = True

        def two_g_runs(c: dict) -> None:
            c["certification_lineage"]["runs"]["G"]["runs"].pop()

        cases = {
            "release": (mutate_release, "release/version"),
            "version": (mutate_version, "release/version"),
            "schema": (mutate_schema, "unsupported RC contract schema"),
            "count": (drop_artifact, "exactly eight artifacts"),
            "filename": (duplicate_filename, "unique and non-empty"),
            "target": (swap_target, "one per release target"),
            "historical": (bad_historical, "historical certification hash"),
            "runner": (missing_runner, "runner_script"),
            "equality": (equality_required, "cross-runner equality"),
            "smoke-count": (drop_smoke, "five smoke cases"),
            "smoke-program": (duplicate_program, "each of the five programs"),
            "placeholder": (missing_placeholder, r"\{output\}"),
            "threads": (many_threads, "-num_threads 1"),
            "publication": (published, "publication flags"),
            "g-runs": (two_g_runs, "three consecutive green"),
        }
        for name, (mutate, message) in cases.items():
            with self.subTest(name), tempfile.TemporaryDirectory() as temp_dir:
                contract = raw_contract()
                mutate(contract)
                write_contract(Path(temp_dir), contract)
                with self.assertRaisesRegex(release.ReleaseFailure, message):
                    release.load_contract(Path(temp_dir))

    def test_load_contract_accepts_null_historical_hash(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            contract = raw_contract()
            contract["artifacts"][0]["historical_certification_binary_sha256"] = None
            write_contract(Path(temp_dir), contract)
            release.load_contract(Path(temp_dir))


class LineageTests(unittest.TestCase):
    def setUp(self) -> None:
        self.repo = TempRepo()
        self.addCleanup(self.repo.close)
        self.certified = self.repo.head()

    def validate(self, contract: dict, **kwargs: object) -> bool:
        return release.validate_candidate(
            self.repo.root, self.repo.head(), contract, **kwargs
        )

    def test_accepts_filled_lineage_with_candidate_equal_to_head(self) -> None:
        self.assertTrue(self.validate(self.repo.contract(self.certified)))

    def test_accepts_ancestor_candidate_with_only_document_changes(self) -> None:
        contract = self.repo.contract(self.certified)
        commit_file(self.repo.root, "README.md", "new readme\n")
        commit_file(self.repo.root, "docs/release/v0.2.0.md", "go record\n")
        commit_file(
            self.repo.root,
            "LOSAT/tests/prepare_release_candidate_v020.py",
            "# tooling only\n",
        )
        commit_file(
            self.repo.root,
            "LOSAT/tests/test_prepare_release_candidate_v020.py",
            "# tooling only\n",
        )
        self.assertTrue(self.validate(contract))

    def test_rejects_runtime_input_diff(self) -> None:
        contract = self.repo.contract(self.certified)
        for relative in ("LOSAT/src/lib.rs", "LOSAT/tests/fixture.txt"):
            with self.subTest(relative):
                commit_file(self.repo.root, relative, f"changed after {relative}\n")
                with self.assertRaisesRegex(release.ReleaseFailure, "certification reuse is invalid"):
                    self.validate(contract)
                git(self.repo.root, "reset", "-q", "--hard", self.certified)

    def test_rejects_unfilled_candidate_sha(self) -> None:
        for value in ("", "abc", "A" * 40, None):
            with self.subTest(value=value):
                contract = self.repo.contract(self.certified)
                contract["certification_lineage"]["candidate_sha"] = value
                with self.assertRaisesRegex(release.ReleaseFailure, "contract not filled"):
                    self.validate(contract)

    def test_rejects_lineage_candidate_that_is_not_an_ancestor(self) -> None:
        contract = self.repo.contract(self.certified)
        git(self.repo.root, "checkout", "-q", "-b", "other", self.certified)
        other = commit_file(self.repo.root, "README.md", "other branch\n")
        git(self.repo.root, "checkout", "-q", "-")
        git(self.repo.root, "reset", "-q", "--hard", self.certified)
        commit_file(self.repo.root, "README.md", "main branch\n")
        contract["certification_lineage"]["candidate_sha"] = other
        with self.assertRaisesRegex(release.ReleaseFailure, "does not descend"):
            self.validate(contract)

    def test_rejects_head_mismatch_and_version_mismatch(self) -> None:
        contract = self.repo.contract(self.certified)
        with self.assertRaisesRegex(release.ReleaseFailure, "candidate SHA mismatch"):
            release.validate_candidate(self.repo.root, "a" * 40, contract)
        commit_file(
            self.repo.root, "LOSAT/Cargo.toml", '[package]\nname = "LOSAT"\nversion = "0.1.0"\n'
        )
        with self.assertRaisesRegex(release.ReleaseFailure, "package version"):
            self.validate(contract)

    def test_rejects_missing_v010_ancestor(self) -> None:
        contract = self.repo.contract(self.certified)
        contract["certification_lineage"]["v010_ancestors"]["integrated_runtime"] = "b" * 40
        with self.assertRaises(release.ReleaseFailure):
            self.validate(contract)

    def test_rejects_missing_or_failed_run_ids(self) -> None:
        def clear(label: str):
            def change(lineage: dict) -> None:
                lineage["runs"][label]["run_id"] = None
            return change

        def fail(label: str):
            def change(lineage: dict) -> None:
                lineage["runs"][label]["conclusion"] = "failure"
            return change

        def zero(lineage: dict) -> None:
            lineage["runs"]["B"]["run_id"] = 0

        def boolean(lineage: dict) -> None:
            lineage["runs"]["C"]["run_id"] = True

        def g_missing(lineage: dict) -> None:
            lineage["runs"]["G"]["runs"][2]["run_id"] = None

        def g_failed(lineage: dict) -> None:
            lineage["runs"]["G"]["runs"][1]["conclusion"] = "cancelled"

        changes = {
            "B missing": clear("B"),
            "C missing": clear("C"),
            "F missing": clear("F"),
            "H missing": clear("H"),
            "F failed": fail("F"),
            "zero": zero,
            "bool": boolean,
            "G missing": g_missing,
            "G failed": g_failed,
        }
        for name, change in changes.items():
            with self.subTest(name):
                contract = self.repo.contract(self.certified)
                change(contract["certification_lineage"])
                with self.assertRaisesRegex(release.ReleaseFailure, "run "):
                    self.validate(contract)

    def test_run_i_may_be_unrecorded(self) -> None:
        contract = self.repo.contract(self.certified)
        contract["certification_lineage"]["runs"]["I"].update(run_id=None, conclusion=None)
        self.assertTrue(self.validate(contract))

    def test_rejects_missing_evidence(self) -> None:
        contract = self.repo.contract(self.certified)
        contract["certification_lineage"]["local_evidence"]["D"]["evidence_sha256"] = None
        with self.assertRaisesRegex(release.ReleaseFailure, "local evidence D"):
            self.validate(contract)
        contract = self.repo.contract(self.certified)
        contract["certification_lineage"]["ws_evidence"]["J"]["speed_grid"]["record"] = None
        with self.assertRaisesRegex(release.ReleaseFailure, "ws evidence J speed_grid"):
            self.validate(contract)

    def test_unfilled_lineage_flag_skips_lineage_and_warns(self) -> None:
        contract = raw_contract()
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("GITHUB_ACTIONS", None)
            with mock.patch("sys.stderr", new_callable=io.StringIO) as stderr:
                checked = self.validate(contract, allow_unfilled_lineage=True)
        self.assertFalse(checked)
        self.assertIn("NOT checked", stderr.getvalue())

    def test_unfilled_lineage_flag_still_checks_head_and_version(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("GITHUB_ACTIONS", None)
            with self.assertRaisesRegex(release.ReleaseFailure, "candidate SHA mismatch"):
                release.validate_candidate(
                    self.repo.root, "c" * 40, raw_contract(), allow_unfilled_lineage=True
                )

    def test_unfilled_lineage_flag_is_refused_under_github_actions(self) -> None:
        with mock.patch.dict(os.environ, {"GITHUB_ACTIONS": "true"}):
            with self.assertRaisesRegex(release.ReleaseFailure, "GITHUB_ACTIONS"):
                self.validate(raw_contract(), allow_unfilled_lineage=True)

    def test_unfilled_repository_contract_fails_without_the_flag(self) -> None:
        with self.assertRaisesRegex(release.ReleaseFailure, "contract not filled"):
            self.validate(raw_contract())


class ArchiveTests(unittest.TestCase):
    def test_tar_and_zip_assembly_are_deterministic(self) -> None:
        entries = [
            ("LOSAT-0.2.0-test/losat", b"binary\0bytes", 0o755),
            ("LOSAT-0.2.0-test/README.md", b"readme\n", 0o644),
        ]
        self.assertEqual(
            hashlib.sha256(release.tar_bytes(entries)).digest(),
            hashlib.sha256(release.tar_bytes(reversed(entries))).digest(),
        )
        self.assertEqual(
            hashlib.sha256(release.zip_bytes(entries)).digest(),
            hashlib.sha256(release.zip_bytes(reversed(entries))).digest(),
        )

    def test_archive_extraction_rejects_path_escape(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            archive = root / "unsafe.tar.gz"
            with tarfile.open(archive, "w:gz") as handle:
                info = tarfile.TarInfo("LOSAT-0.2.0-test/../outside")
                data = b"unsafe"
                info.size = len(data)
                handle.addfile(info, io.BytesIO(data))
            with self.assertRaisesRegex(release.ReleaseFailure, "unsafe archive member"):
                release.extract_archive(archive, root / "extract", "LOSAT-0.2.0-test")

    def test_binary_architecture_recognizes_supported_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            for machine, expected in ((62, "x86_64"), (183, "arm64")):
                elf = bytearray(64)
                elf[:6] = b"\x7fELF\x02\x01"
                elf[18:20] = machine.to_bytes(2, "little")
                path = root / f"elf-{machine}"
                path.write_bytes(elf)
                self.assertEqual(release.binary_architecture(path), expected)

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

            macho[4:8] = (0x01000007).to_bytes(4, "little")
            macho_path.write_bytes(macho)
            self.assertEqual(release.binary_architecture(macho_path), "x86_64")

    def test_expected_architecture_covers_all_native_targets(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        observed = {
            item["target"]: release.expected_architecture(item["target"])
            for item in contract["artifacts"]
            if item["kind"] == "native"
        }
        self.assertEqual(observed["aarch64-unknown-linux-gnu"], "arm64")
        self.assertEqual(observed["aarch64-apple-darwin"], "arm64")
        self.assertEqual(observed["x86_64-pc-windows-msvc"], "x86_64")
        self.assertEqual(observed["x86_64-apple-darwin"], "x86_64")

    def test_source_crate_rejects_generated_content(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            crate = Path(temp_dir) / "LOSAT-0.2.0.crate"
            with tarfile.open(crate, "w:gz") as handle:
                for name in (
                    "LOSAT-0.2.0/Cargo.toml",
                    "LOSAT-0.2.0/Cargo.lock",
                    "LOSAT-0.2.0/README.md",
                    "LOSAT-0.2.0/target/debug/losat",
                ):
                    data = (
                        b'[package]\nname = "LOSAT"\nversion = "0.2.0"\nlicense = "MIT"\n'
                        if name.endswith("Cargo.toml")
                        else b"fixture"
                    )
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mode = 0o644
                    handle.addfile(info, io.BytesIO(data))
            with self.assertRaisesRegex(release.ReleaseFailure, "generated/scratch"):
                release.crate_members(crate, "LOSAT-0.2.0")


class SmokeTests(unittest.TestCase):
    def fake_runner(self, contract: dict, commands: list[list[str]]):
        by_program = {case["program"]: case for case in contract["smoke"]}

        def fake_run_capture(command, cwd, *, clean_losat_environment=False):
            del cwd
            assert clean_losat_environment
            commands.append(list(command))
            program = next(item for item in command if item in by_program)
            Path(command[command.index("-out") + 1]).write_bytes(
                f"output-of-{program}\n".encode()
            )
            return subprocess.CompletedProcess(command, 0, "", "")

        return fake_run_capture

    def contract_for_fake_outputs(self) -> dict:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        for case in contract["smoke"]:
            case["expected_output_sha256"] = hashlib.sha256(
                f"output-of-{case['program']}\n".encode()
            ).hexdigest()
        return contract

    def test_smoke_runs_all_five_programs_with_absolute_inputs(self) -> None:
        contract = self.contract_for_fake_outputs()
        commands: list[list[str]] = []
        with mock.patch.object(
            release, "run_capture", side_effect=self.fake_runner(contract, commands)
        ):
            result = release.run_smoke(["node", "runner.js", "losat.wasm"], REPO_ROOT, contract)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(
            [item["program"] for item in result["cases"]],
            [case["program"] for case in contract["smoke"]],
        )
        self.assertEqual(
            [item["case_id"] for item in result["cases"]],
            [case["case_id"] for case in contract["smoke"]],
        )
        self.assertEqual(len(commands), 5)
        for command, case in zip(commands, contract["smoke"]):
            self.assertEqual(command[:4], ["node", "runner.js", "losat.wasm", case["program"]])
            query = Path(command[command.index("-query") + 1])
            subject = Path(command[command.index("-subject") + 1])
            self.assertTrue(query.is_absolute() and subject.is_absolute())
            self.assertEqual(query, (REPO_ROOT / case["query"]).resolve())
            self.assertEqual(subject, (REPO_ROOT / case["subject"]).resolve())
            self.assertEqual(command[command.index("-num_threads") + 1], "1")
            self.assertNotIn("{output}", command)

    def test_wasm_command_prefix_uses_the_contract_runner(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        wasm = Path("/x/losat.wasm")
        for item in contract["artifacts"]:
            if item["kind"] != "wasm":
                continue
            runner = (REPO_ROOT / item["runner_script"]).resolve()
            prefix = release.wasm_command_prefix("node", runner, wasm)
            self.assertEqual(prefix, ["node", "--no-warnings", str(runner), str(wasm)])
            command = release.smoke_command(prefix, REPO_ROOT, contract["smoke"][0], Path("/o"))
            self.assertEqual(command[:5], [*prefix, "blastn"])
            self.assertEqual(command[-2:], ["-out", "/o"])
        self.assertTrue(
            release.wasm_command_prefix("node", Path("run_losat_wasi_threads.js"), wasm)[2]
            .endswith("run_losat_wasi_threads.js")
        )

    def test_native_smoke_command_has_no_wrapper(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        command = release.smoke_command(["/bin/losat"], REPO_ROOT, contract["smoke"][4], Path("/o"))
        self.assertEqual(command[:2], ["/bin/losat", "tblastx"])
        self.assertEqual(command[-2:], ["-out", "/o"])

    def test_smoke_rejects_hash_mismatch_and_console_output(self) -> None:
        contract = self.contract_for_fake_outputs()
        contract["smoke"][3]["expected_output_sha256"] = "0" * 64
        with mock.patch.object(
            release, "run_capture", side_effect=self.fake_runner(contract, [])
        ):
            with self.assertRaisesRegex(release.ReleaseFailure, "tblastn smoke hash mismatch"):
                release.run_smoke(["losat"], REPO_ROOT, contract)

        contract = self.contract_for_fake_outputs()

        def noisy(command, cwd, *, clean_losat_environment=False):
            del cwd, clean_losat_environment
            Path(command[command.index("-out") + 1]).write_bytes(b"x")
            return subprocess.CompletedProcess(command, 0, "", "warning\n")

        with mock.patch.object(release, "run_capture", side_effect=noisy):
            with self.assertRaisesRegex(release.ReleaseFailure, "unexpected console output"):
                release.run_smoke(["losat"], REPO_ROOT, contract)

    def test_smoke_requires_an_output_file_and_fixtures(self) -> None:
        contract = self.contract_for_fake_outputs()

        def silent(command, cwd, *, clean_losat_environment=False):
            del cwd, clean_losat_environment
            return subprocess.CompletedProcess(command, 0, "", "")

        with mock.patch.object(release, "run_capture", side_effect=silent):
            with self.assertRaisesRegex(release.ReleaseFailure, "did not create output"):
                release.run_smoke(["losat"], REPO_ROOT, contract)
        contract["smoke"][0]["query"] = "LOSAT/tests/fasta/does_not_exist.fa"
        with self.assertRaisesRegex(release.ReleaseFailure, "fixture is missing"):
            release.run_smoke(["losat"], REPO_ROOT, contract)

    def test_version_output_must_be_exactly_losat_0_2_0(self) -> None:
        contract = release.load_contract(REPO_ROOT)
        for stdout, stderr, ok in (
            ("losat 0.2.0\n", "", True),
            ("losat 0.2.0\r\n", "", True),
            ("losat 0.1.0\n", "", False),
            ("LOSAT 0.2.0\n", "", False),
            ("losat 0.2.0 extra\n", "", False),
            ("losat 0.2.0\n", "warn\n", False),
        ):
            with self.subTest(stdout=stdout, stderr=stderr):
                result = subprocess.CompletedProcess([], 0, stdout, stderr)
                with mock.patch.object(release, "run_capture", return_value=result):
                    if ok:
                        self.assertEqual(
                            release.version_output(["losat"], REPO_ROOT, contract), "losat 0.2.0"
                        )
                    else:
                        with self.assertRaises(release.ReleaseFailure):
                            release.version_output(["losat"], REPO_ROOT, contract)


class AssemblyTests(unittest.TestCase):
    def assemble(self, lineage_checked: bool, *, kind: str = "native") -> tuple[dict, Path]:
        contract = deepcopy(release.load_contract(REPO_ROOT))
        target = "x86_64-unknown-linux-gnu"
        spec = release.artifact_spec(contract, "native", target)
        spec["historical_certification_binary_sha256"] = "f" * 64
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        binary = root / "losat"
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
            allow_unfilled_lineage=not lineage_checked,
        )
        smoke = {"status": "PASS", "cases": []}
        with (
            mock.patch.object(release, "load_contract", return_value=contract),
            mock.patch.object(release, "validate_candidate", return_value=lineage_checked),
            mock.patch.object(release, "capture_toolchain", return_value={}),
            mock.patch.object(release, "version_output", return_value="losat 0.2.0"),
            mock.patch.object(release, "run_smoke", return_value=smoke),
        ):
            release.assemble_binary(args, kind)
        sidecar = json.loads(
            (output_dir / f"{spec['filename']}.metadata.json").read_text(encoding="utf-8")
        )
        return sidecar, output_dir / spec["filename"]

    def test_native_assembly_records_hashes_and_lineage_flag(self) -> None:
        sidecar, archive = self.assemble(True)
        self.assertTrue(sidecar["lineage_checked"])
        self.assertEqual(sidecar["candidate_sha"], "a" * 40)
        self.assertEqual(sidecar["release"], "v0.2.0")
        self.assertEqual(sidecar["binary"]["architecture"], "x86_64")
        self.assertFalse(sidecar["certification_binary_reference"]["same_bytes_in_this_build"])
        self.assertEqual(sidecar["certification_binary_reference"]["sha256"], "f" * 64)
        self.assertEqual(sidecar["artifact_sha256"], release.sha256_path(archive))
        with tarfile.open(archive) as handle:
            names = sorted(handle.getnames())
        root = archive.name.removesuffix(".tar.gz")
        self.assertEqual(
            names,
            [f"{root}/LICENSE", f"{root}/README.md", f"{root}/RELEASE-METADATA.json", f"{root}/losat"],
        )

    def test_dry_run_assembly_is_marked_unchecked(self) -> None:
        sidecar, _ = self.assemble(False)
        self.assertIs(sidecar["lineage_checked"], False)

    def test_assembly_is_byte_deterministic(self) -> None:
        _, first = self.assemble(True)
        _, second = self.assemble(True)
        self.assertEqual(release.sha256_path(first), release.sha256_path(second))

    def test_main_refuses_the_dry_run_flag_for_aggregate(self) -> None:
        with self.assertRaises(SystemExit), mock.patch("sys.stderr", new_callable=io.StringIO):
            release.build_parser().parse_args(
                [
                    "aggregate",
                    "--repo-root", ".",
                    "--candidate-sha", "a" * 40,
                    "--output-dir", "o",
                    "--input-root", "i",
                    "--workflow-run-id", "1",
                    "--allow-unfilled-lineage",
                ]
            )
        for command in ("native", "wasm", "source"):
            extra = {
                "native": ["--target", "x86_64-unknown-linux-gnu", "--binary", "b"],
                "wasm": ["--binary", "b"],
                "source": ["--crate", "c"],
            }[command]
            parsed = release.build_parser().parse_args(
                [command, "--repo-root", ".", "--candidate-sha", "a" * 40,
                 "--output-dir", "o", *extra, "--allow-unfilled-lineage"]
            )
            self.assertTrue(parsed.allow_unfilled_lineage)

    def test_main_reports_failures_with_exit_code_one(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir, mock.patch(
            "sys.stderr", new_callable=io.StringIO
        ):
            code = release.main(
                [
                    "native",
                    "--repo-root", temp_dir,
                    "--candidate-sha", "a" * 40,
                    "--output-dir", str(Path(temp_dir) / "o"),
                    "--target", "x86_64-unknown-linux-gnu",
                    "--binary", str(Path(temp_dir) / "missing"),
                ]
            )
        self.assertEqual(code, 1)


class AggregateTests(unittest.TestCase):
    candidate_sha = "a" * 40

    def populate(self, inputs: Path, contract: dict, *, lineage_checked: bool = True) -> None:
        inputs.mkdir()
        for index, spec in enumerate(contract["artifacts"]):
            artifact = inputs / spec["filename"]
            artifact.write_bytes(f"artifact-{index}\n".encode())
            digest = release.sha256_path(artifact)
            (inputs / f"{artifact.name}.sha256").write_text(
                f"{digest}  {artifact.name}\n", encoding="utf-8"
            )
            smoke = {
                "status": "PASS",
                "cases": [
                    {
                        "program": case["program"],
                        "case_id": case["case_id"],
                        "output_sha256": case["expected_output_sha256"],
                        "status": "PASS",
                    }
                    for case in contract["smoke"]
                ],
            }
            sidecar = {
                "candidate_sha": self.candidate_sha,
                "release": contract["release"],
                "lineage_checked": lineage_checked,
                "artifact": {
                    "kind": spec["kind"],
                    "target": spec["target"],
                    "filename": spec["filename"],
                    **({"sha256": digest} if spec["kind"] == "source" else {}),
                },
                "certification_lineage": contract["certification_lineage"],
                "binary_hash_policy": contract["binary_hash_policy"],
                "publication": contract["publication"],
                "status": "PASS",
            }
            if spec["kind"] == "source":
                sidecar["clean_install"] = {"smoke": smoke}
            else:
                sidecar["artifact_sha256"] = digest
                sidecar["extracted_artifact"] = {"smoke": smoke}
            release.write_json(inputs / f"{artifact.name}.metadata.json", sidecar)

    def run_aggregate(self, contract: dict, inputs: Path, outputs: Path) -> None:
        args = argparse.Namespace(
            repo_root=REPO_ROOT,
            candidate_sha=self.candidate_sha,
            input_root=inputs,
            output_dir=outputs,
            workflow_run_id="123456",
        )
        with (
            mock.patch.object(release, "load_contract", return_value=contract),
            mock.patch.object(release, "validate_candidate", return_value=True),
        ):
            release.aggregate(args)

    def test_aggregate_requires_and_hashes_the_exact_eight_files(self) -> None:
        contract = filled_contract("d" * 40)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.populate(root / "inputs", contract)
            self.run_aggregate(contract, root / "inputs", root / "outputs")
            outputs = root / "outputs"
            handoff = json.loads((outputs / "RC-HANDOFF.json").read_text(encoding="utf-8"))
            self.assertEqual(handoff["decision"], "RC_HANDOFF_READY")
            self.assertEqual(handoff["candidate_sha"], self.candidate_sha)
            self.assertEqual(handoff["certified_candidate_sha"], "d" * 40)
            self.assertEqual(len(handoff["artifacts"]), 8)
            self.assertEqual(len((outputs / "SHA256SUMS").read_text().splitlines()), 8)
            text = (outputs / "RC-HANDOFF.md").read_text(encoding="utf-8")
            self.assertTrue(text.startswith("# LOSAT v0.2.0 release-candidate handoff"))
            self.assertIn(self.candidate_sha, text)
            self.assertIn("d" * 40, text)
            self.assertIn("123456", text)
            self.assertIn("Rust `1.92.0`", text)
            self.assertIn("Node `24.21.0`", text)
            self.assertIn("2001, 2002", text.replace("`", ""))
            self.assertIn("No tag / release / upload performed", text)
            for claim in contract["platform_claims"]:
                self.assertIn(claim["platform"], text)
                self.assertIn(claim["wording"], text)

    def test_aggregate_refuses_a_dry_run_sidecar(self) -> None:
        contract = filled_contract("d" * 40)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.populate(root / "inputs", contract)
            first = contract["artifacts"][3]["filename"]
            path = root / "inputs" / f"{first}.metadata.json"
            sidecar = json.loads(path.read_text(encoding="utf-8"))
            sidecar["lineage_checked"] = False
            release.write_json(path, sidecar)
            with self.assertRaisesRegex(release.ReleaseFailure, "lineage_checked"):
                self.run_aggregate(contract, root / "inputs", root / "outputs")
            self.assertFalse((root / "outputs" / "RC-HANDOFF.md").exists())
            self.assertFalse((root / "outputs" / "RC-HANDOFF.json").exists())

    def test_aggregate_refuses_a_missing_lineage_flag(self) -> None:
        contract = filled_contract("d" * 40)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.populate(root / "inputs", contract)
            path = root / "inputs" / f"{contract['artifacts'][0]['filename']}.metadata.json"
            sidecar = json.loads(path.read_text(encoding="utf-8"))
            del sidecar["lineage_checked"]
            release.write_json(path, sidecar)
            with self.assertRaisesRegex(release.ReleaseFailure, "lineage_checked"):
                self.run_aggregate(contract, root / "inputs", root / "outputs")

    def test_aggregate_requires_the_exact_filename_set(self) -> None:
        contract = filled_contract("d" * 40)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.populate(root / "inputs", contract)
            (root / "inputs" / "LOSAT-0.2.0-extra.tar.gz").write_bytes(b"extra")
            with self.assertRaisesRegex(release.ReleaseFailure, "differs from the eight"):
                self.run_aggregate(contract, root / "inputs", root / "outputs")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.populate(root / "inputs", contract)
            (root / "inputs" / contract["artifacts"][6]["filename"]).unlink()
            with self.assertRaisesRegex(release.ReleaseFailure, "differs from the eight"):
                self.run_aggregate(contract, root / "inputs", root / "outputs")

    def test_aggregate_rejects_checksum_and_smoke_drift(self) -> None:
        contract = filled_contract("d" * 40)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.populate(root / "inputs", contract)
            name = contract["artifacts"][1]["filename"]
            (root / "inputs" / name).write_bytes(b"tampered")
            with self.assertRaisesRegex(release.ReleaseFailure, "checksum sidecar mismatch"):
                self.run_aggregate(contract, root / "inputs", root / "outputs")
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.populate(root / "inputs", contract)
            path = root / "inputs" / f"{contract['artifacts'][2]['filename']}.metadata.json"
            sidecar = json.loads(path.read_text(encoding="utf-8"))
            sidecar["extracted_artifact"]["smoke"]["cases"][0]["output_sha256"] = "0" * 64
            release.write_json(path, sidecar)
            with self.assertRaisesRegex(release.ReleaseFailure, "smoke record"):
                self.run_aggregate(contract, root / "inputs", root / "outputs")


class WorkflowTests(unittest.TestCase):
    @unittest.skipUnless(WORKFLOW_PATH.exists(), "PR 3 adds the workflow")
    def test_release_workflow_is_manual_read_only_and_non_publishing(self) -> None:
        workflow = WORKFLOW_PATH.read_text(encoding="utf-8")
        self.assertIn("workflow_dispatch:", workflow)
        self.assertIn("contents: read", workflow)
        self.assertNotIn("actions/create-release", workflow)
        self.assertNotIn("softprops/action-gh-release", workflow)
        self.assertNotIn("cargo publish\n", workflow)
        self.assertNotIn("git tag", workflow)
        self.assertNotIn("--allow-unfilled-lineage", workflow)
        for target in (
            "x86_64-unknown-linux-gnu",
            "aarch64-unknown-linux-gnu",
            "x86_64-pc-windows-msvc",
            "aarch64-apple-darwin",
            "x86_64-apple-darwin",
            "wasm32-wasip1",
            "wasm32-wasip1-threads",
        ):
            self.assertIn(target, workflow)
        self.assertIn("prepare_release_candidate_v020.py", workflow)


if __name__ == "__main__":
    unittest.main()
