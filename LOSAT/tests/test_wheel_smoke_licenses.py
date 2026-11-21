#!/usr/bin/env python3
"""Tests for the license-file check of python/tools/wheel_smoke.py (`license_failures`).

No NCBI counterpart: this is LOSAT's own packaging tooling. .github/workflows/python-wheels.yml runs
the check on a real wheel, but only for pull requests that touch python/**; these tests build a
fake wheel in a temporary directory, so they need neither maturin nor the network.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("wheel_smoke", REPO_ROOT / "python" / "tools" / "wheel_smoke.py")
assert SPEC is not None and SPEC.loader is not None
wheel_smoke = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = wheel_smoke
SPEC.loader.exec_module(wheel_smoke)

FILES = {name: f"text of {name}\n".encode() for name in wheel_smoke.LICENSE_FILES}


def metadata(declared: tuple[str, ...]) -> bytes:
    lines = ["Metadata-Version: 2.4", "Name: pylosat", "Version: 0.0.0", *(f"License-File: {n}" for n in declared)]
    return ("\n".join(lines) + "\n\nThe description, which may mention License-File: nothing\n").encode()


class LicenseFailuresTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temp = tempfile.TemporaryDirectory()
        self.addCleanup(self._temp.cleanup)
        self.repo = Path(self._temp.name) / "repo"
        (self.repo / "python").mkdir(parents=True)
        for name, data in FILES.items():
            (self.repo / "python" / name).write_bytes(data)

    def wheel(self, *, files: dict[str, bytes] | None = None, declared: tuple[str, ...] | None = None) -> Path:
        files = FILES if files is None else files
        declared = tuple(FILES) if declared is None else declared
        path = Path(self._temp.name) / "pylosat-0.0.0-py3-none-any.whl"
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("pylosat/__init__.py", "")
            archive.writestr("pylosat-0.0.0.dist-info/METADATA", metadata(declared))
            for name, data in files.items():
                archive.writestr(f"pylosat-0.0.0.dist-info/licenses/{name}", data)
        return path

    def test_a_wheel_with_the_files_of_python_passes(self) -> None:
        self.assertEqual(wheel_smoke.license_failures(self.wheel(), self.repo), [])

    def test_a_changed_file_is_reported(self) -> None:
        files = dict(FILES, THIRD_PARTY_LICENSES=b"stale\n")
        self.assertEqual(wheel_smoke.license_failures(self.wheel(files=files), self.repo),
                         ["pylosat-0.0.0.dist-info/licenses/THIRD_PARTY_LICENSES differs from python/THIRD_PARTY_LICENSES"])

    def test_a_missing_file_is_reported(self) -> None:
        files = {name: data for name, data in FILES.items() if name != "THIRD_PARTY_NOTICES"}
        self.assertEqual(wheel_smoke.license_failures(self.wheel(files=files), self.repo),
                         ["no .dist-info/licenses/THIRD_PARTY_NOTICES"])

    def test_a_file_that_metadata_does_not_name_is_reported(self) -> None:
        declared = ("LICENSE", "THIRD_PARTY_LICENSES")
        self.assertEqual(wheel_smoke.license_failures(self.wheel(declared=declared), self.repo),
                         ["METADATA has no License-File: THIRD_PARTY_NOTICES"])

    def test_only_the_header_of_metadata_counts(self) -> None:
        # "License-File:" in the long description (after the first empty line) is not a declaration.
        self.assertEqual(wheel_smoke.license_failures(self.wheel(declared=()), self.repo),
                         [f"METADATA has no License-File: {name}" for name in wheel_smoke.LICENSE_FILES])


if __name__ == "__main__":
    unittest.main()
