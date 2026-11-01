#!/usr/bin/env python3
"""Tests for the v0.2.0 downloaded-archive verifier (argument, sidecar and path checks).

The smoke itself (version and five searches) is covered by the assembler's tests.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import sys
import tarfile
import tempfile
import shutil
import subprocess
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

import verify_release_archive_v020 as verify  # noqa: E402
from test_prepare_release_candidate_v020 import _FILL_MODULE  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
NATIVE_TARGET = "x86_64-unknown-linux-gnu"
NATIVE_ARCHIVE = "LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz"
ROOT = NATIVE_ARCHIVE.removesuffix(".tar.gz")


def write_tar(path: Path, members: dict[str, bytes]) -> None:
    with tarfile.open(path, "w:gz") as handle:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mode = 0o644
            handle.addfile(info, io.BytesIO(data))


def write_sidecar(archive: Path, checksum: str | None = None) -> None:
    value = checksum or hashlib.sha256(archive.read_bytes()).hexdigest()
    archive.with_name(f"{archive.name}.sha256").write_text(
        f"{value}  {archive.name}\n", encoding="utf-8"
    )


def run_main(argv: list[str]) -> tuple[int, str]:
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        code = verify.main(argv)
    return code, err.getvalue()


class VerifyReleaseArchiveTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temp = tempfile.TemporaryDirectory()
        self.addCleanup(self._temp.cleanup)
        self.temp = Path(self._temp.name)

    def argv(self, archive: Path, target: str = NATIVE_TARGET) -> list[str]:
        return ["--repo-root", str(REPO_ROOT), "--archive", str(archive), "--target", target]

    def test_missing_arguments_are_refused(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                verify.main(["--repo-root", str(REPO_ROOT)])
        self.assertEqual(caught.exception.code, 2)

    def test_unknown_target_is_refused(self) -> None:
        archive = self.temp / NATIVE_ARCHIVE
        write_tar(archive, {f"{ROOT}/losat": b"x"})
        write_sidecar(archive)
        code, message = run_main(self.argv(archive, "riscv64-unknown-linux-gnu"))
        self.assertEqual(code, 1)
        self.assertIn("not one artifact of the RC contract", message)

    def test_archive_name_must_match_the_contract(self) -> None:
        archive = self.temp / "renamed.tar.gz"
        write_tar(archive, {f"{ROOT}/losat": b"x"})
        write_sidecar(archive)
        code, message = run_main(self.argv(archive))
        self.assertEqual(code, 1)
        self.assertIn("differs from the contract", message)

    def test_missing_archive_is_refused(self) -> None:
        code, message = run_main(self.argv(self.temp / NATIVE_ARCHIVE))
        self.assertEqual(code, 1)
        self.assertIn("archive is missing", message)

    def test_missing_sidecar_is_refused(self) -> None:
        archive = self.temp / NATIVE_ARCHIVE
        write_tar(archive, {f"{ROOT}/losat": b"x"})
        code, message = run_main(self.argv(archive))
        self.assertEqual(code, 1)
        self.assertIn("cannot read checksum sidecar", message)

    def test_sidecar_mismatch_is_detected(self) -> None:
        archive = self.temp / NATIVE_ARCHIVE
        write_tar(archive, {f"{ROOT}/losat": b"x"})
        write_sidecar(archive, "0" * 64)
        code, message = run_main(self.argv(archive))
        self.assertEqual(code, 1)
        self.assertIn("differs from its sidecar", message)

    def test_malformed_sidecar_is_refused(self) -> None:
        archive = self.temp / NATIVE_ARCHIVE
        write_tar(archive, {f"{ROOT}/losat": b"x"})
        archive.with_name(f"{archive.name}.sha256").write_text("not a checksum\n")
        code, message = run_main(self.argv(archive))
        self.assertEqual(code, 1)
        self.assertIn("is not '<sha256>", message)

    def test_path_escape_is_refused(self) -> None:
        archive = self.temp / NATIVE_ARCHIVE
        write_tar(archive, {f"{ROOT}/losat": b"x", f"{ROOT}/../escape": b"y"})
        write_sidecar(archive)
        code, message = run_main(self.argv(archive))
        self.assertEqual(code, 1)
        self.assertIn("unsafe archive member", message)
        self.assertFalse((self.temp / "escape").exists())

    def test_member_outside_the_expected_root_is_refused(self) -> None:
        archive = self.temp / NATIVE_ARCHIVE
        write_tar(archive, {"other-root/losat": b"x"})
        write_sidecar(archive)
        code, message = run_main(self.argv(archive))
        self.assertEqual(code, 1)
        self.assertIn("outside expected root", message)

    def test_corrupt_archive_is_refused(self) -> None:
        archive = self.temp / NATIVE_ARCHIVE
        archive.write_bytes(b"not a gzip stream")
        write_sidecar(archive)
        code, message = run_main(self.argv(archive))
        self.assertEqual(code, 1)
        self.assertIn("cannot extract", message)

    def test_archive_without_the_binary_is_refused(self) -> None:
        archive = self.temp / NATIVE_ARCHIVE
        write_tar(archive, {f"{ROOT}/LICENSE": b"x"})
        write_sidecar(archive)
        code, message = run_main(self.argv(archive))
        self.assertEqual(code, 1)
        self.assertIn("has no losat", message)

    @unittest.skipUnless(shutil.which("node"), "needs Node.js")
    def test_threaded_archive_must_carry_a_verified_shared_memory_guard(self) -> None:
        target = "wasm32-wasip1-threads"
        contract = verify.rc.load_contract(REPO_ROOT)
        spec = next(item for item in contract["artifacts"] if item["target"] == target)
        root = spec["filename"].removesuffix(".tar.gz")
        raw, guarded = self.temp / "raw.wasm", self.temp / "guarded.wasm"
        raw.write_bytes(_FILL_MODULE)
        subprocess.run(
            ["node", str(REPO_ROOT / "LOSAT/tests/wasi_shared_memory.js"), "guard", str(raw), str(guarded)],
            check=True,
            capture_output=True,
        )
        archive = self.temp / spec["filename"]
        write_tar(archive, {f"{root}/{spec['binary_name']}": raw.read_bytes()})
        write_sidecar(archive)
        code, message = run_main(self.argv(archive, target))
        self.assertEqual(code, 1)
        self.assertIn("no verified shared-memory guard", message)
        write_tar(archive, {f"{root}/{spec['binary_name']}": guarded.read_bytes()})
        write_sidecar(archive)
        with (
            mock.patch.object(verify.rc, "version_output", return_value="losat 0.2.0"),
            mock.patch.object(verify.rc, "run_smoke", return_value={"status": "PASS", "cases": []}),
        ):
            result = verify.verify(REPO_ROOT, archive, target, "node", None)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["shared_memory_guard"]["fill"], 1)

    def test_source_crate_hash_and_layout_pass(self) -> None:
        archive = self.temp / "LOSAT-0.2.0.crate"
        write_tar(archive, {"LOSAT-0.2.0/Cargo.toml": b"[package]\n"})
        write_sidecar(archive)
        code, message = run_main(self.argv(archive, "source"))
        self.assertEqual((code, message), (0, ""))


if __name__ == "__main__":
    unittest.main()
