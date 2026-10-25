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
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import verify_release_archive_v020 as verify  # noqa: E402

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

    def test_source_crate_hash_and_layout_pass(self) -> None:
        archive = self.temp / "LOSAT-0.2.0.crate"
        write_tar(archive, {"LOSAT-0.2.0/Cargo.toml": b"[package]\n"})
        write_sidecar(archive)
        code, message = run_main(self.argv(archive, "source"))
        self.assertEqual((code, message), (0, ""))


if __name__ == "__main__":
    unittest.main()
