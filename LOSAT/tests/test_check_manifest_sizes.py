#!/usr/bin/env python3
"""Tests for the manifest size check (`check_manifest_sizes.py`).

No NCBI counterpart: the checker reads the size of LOSAT's own fixture lists.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_manifest_sizes.py")
SPEC = importlib.util.spec_from_file_location("manifest_sizes", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
sizes = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = sizes
SPEC.loader.exec_module(sizes)


def run(root: Path) -> tuple[int, str]:
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        code = sizes.main(["--root", str(root)])
    return code, out.getvalue()


def manifest_dir(root: Path) -> Path:
    directory = root / "LOSAT" / "tests" / "db_search_manifest"
    directory.mkdir(parents=True)
    return directory


class ManifestSizeTests(unittest.TestCase):
    def test_files_at_the_limit_pass(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = manifest_dir(Path(tmp))
            (directory / "blastp.tsv").write_bytes(b"x" * sizes.LIMIT)
            (directory / "blastn.tsv").write_bytes(b"row_id\n")
            code, text = run(Path(tmp))
        self.assertEqual(code, 0, text)

    def test_a_file_over_the_limit_fails_and_is_named(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = manifest_dir(Path(tmp))
            (directory / "blastp.tsv").write_bytes(b"x" * (sizes.LIMIT + 1))
            (directory / "blastn.tsv").write_bytes(b"row_id\n")
            code, text = run(Path(tmp))
        self.assertEqual(code, 1, text)
        self.assertIn("blastp.tsv", text)
        self.assertNotIn("blastn.tsv", text)

    def test_other_files_are_not_counted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = manifest_dir(Path(tmp))
            (directory / "README.md").write_bytes(b"x" * (sizes.LIMIT + 1))
            code, text = run(Path(tmp))
        self.assertEqual(code, 0, text)

    def test_a_missing_directory_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            code, text = run(Path(tmp))
        self.assertEqual(code, 2, text)


if __name__ == "__main__":
    unittest.main()
