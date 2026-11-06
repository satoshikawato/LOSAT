#!/usr/bin/env python3
"""Tests for the public API check (`check_public_api.py`).

No NCBI counterpart: the checker compares the library's own surface with its listing.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


CHECKER_PATH = Path(__file__).with_name("check_public_api.py")
SPEC = importlib.util.spec_from_file_location("check_public_api", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
check = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = check
SPEC.loader.exec_module(check)

ROOT = Path(__file__).resolve().parents[2]
RUSTUP_LIST = (
    "1.92.0-x86_64-unknown-linux-gnu (active, default)\n"
    "1.97.0-x86_64-unknown-linux-gnu\n"
    f"{check.TOOLCHAIN}-x86_64-unknown-linux-gnu\n"
)


class ToolchainTests(unittest.TestCase):
    def test_names_drop_the_markers(self) -> None:
        names = check.toolchain_names(RUSTUP_LIST)
        self.assertEqual(names[0], "1.92.0-x86_64-unknown-linux-gnu")
        self.assertEqual(len(names), 3)

    def test_pinned_nightly_is_found_by_its_prefix(self) -> None:
        self.assertTrue(check.has_toolchain(check.toolchain_names(RUSTUP_LIST), check.TOOLCHAIN))

    def test_another_nightly_is_not_the_pinned_one(self) -> None:
        names = ["nightly-x86_64-unknown-linux-gnu", "nightly-2026-05-23-x86_64-unknown-linux-gnu"]
        self.assertFalse(check.has_toolchain(names, check.TOOLCHAIN))

    def test_tool_version(self) -> None:
        self.assertEqual(check.tool_version("cargo-public-api 0.52.0\n"), "0.52.0")
        self.assertIsNone(check.tool_version("something else 1.0\n"))


class ListingTests(unittest.TestCase):
    def test_equal_listings_have_no_diff(self) -> None:
        self.assertEqual(check.listing_diff("pub mod LOSAT\n", "pub mod LOSAT\n"), [])

    def test_a_removed_item_is_a_minus_line(self) -> None:
        diff = check.listing_diff(
            "pub const LOSAT::api::API_OWNED: &[&str]\npub const LOSAT::api::VERSION: &str\n",
            "pub const LOSAT::api::VERSION: &str\n",
        )
        self.assertIn("-pub const LOSAT::api::API_OWNED: &[&str]\n", diff)

    def test_command_pins_the_toolchain_and_leaves_out_blanket_impls(self) -> None:
        command = check.listing_command(Path("/repo"))
        self.assertEqual(command[:3], ["cargo", f"+{check.TOOLCHAIN}", "public-api"])
        self.assertIn(str(Path("/repo/LOSAT/Cargo.toml")), command)
        self.assertEqual(command.count("--simplified"), 1)

    def test_snapshot_exists_and_lists_the_api_module(self) -> None:
        snapshot = (ROOT / check.SNAPSHOT).read_text(encoding="utf-8")
        self.assertIn("pub mod LOSAT::api\n", snapshot)


class WorkflowPinTests(unittest.TestCase):
    """The versions in ci.yml are the ones this script names (a bump changes both)."""

    def setUp(self) -> None:
        self.workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")

    def test_ci_installs_the_pinned_tool_and_toolchain(self) -> None:
        self.assertIn(check.INSTALL_TOOL, self.workflow)
        self.assertIn(f"toolchain: {check.TOOLCHAIN}", self.workflow)

    def test_no_other_nightly_or_tool_version_is_left_in_ci(self) -> None:
        # A partial bump would leave one job (or its RUSTUP_TOOLCHAIN line) on another nightly.
        self.assertEqual(set(re.findall(r"nightly-\d{4}-\d{2}-\d{2}", self.workflow)), {check.TOOLCHAIN})
        self.assertEqual(set(re.findall(r"cargo-public-api@([0-9.]+)", self.workflow)), {check.TOOL_VERSION})

    def test_ci_runs_this_check(self) -> None:
        self.assertIn("python3 LOSAT/tests/check_public_api.py check", self.workflow)

    def test_semver_checks_use_the_same_toolchain(self) -> None:
        self.assertIn(f"cargo +{check.TOOLCHAIN} semver-checks", self.workflow)


class MainTests(unittest.TestCase):
    def test_a_missing_toolchain_is_exit_2_and_names_the_install_command(self) -> None:
        stderr = io.StringIO()
        with mock.patch.object(check.subprocess, "run") as run, contextlib.redirect_stderr(stderr):
            run.return_value = mock.Mock(stdout="1.92.0-x86_64-unknown-linux-gnu (default)\n")
            self.assertEqual(check.main(["check", "--root", str(ROOT)]), 2)
        self.assertIn(check.INSTALL_TOOLCHAIN, stderr.getvalue())

    def test_check_compares_with_the_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / check.SNAPSHOT).parent.mkdir(parents=True)
            (root / check.SNAPSHOT).write_text("pub mod LOSAT\n", encoding="utf-8")
            summary = root / "summary.md"
            # CI sets GITHUB_STEP_SUMMARY for every step: the test writes its fake diff to its own file.
            with mock.patch.dict(check.os.environ, {"GITHUB_STEP_SUMMARY": str(summary)}), mock.patch.object(
                check, "environment_error", return_value=None
            ), mock.patch.object(
                check, "current_listing", side_effect=["pub mod LOSAT\n", "pub mod LOSAT\npub mod LOSAT::x\n"]
            ), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(check.main(["check", "--root", tmp]), 0)
                self.assertFalse(summary.exists())
                self.assertEqual(check.main(["check", "--root", tmp]), 1)
            self.assertIn("+pub mod LOSAT::x\n", summary.read_text(encoding="utf-8"))

    def test_update_writes_the_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / check.SNAPSHOT).parent.mkdir(parents=True)
            with mock.patch.object(check, "environment_error", return_value=None), mock.patch.object(
                check, "current_listing", return_value="pub mod LOSAT\n"
            ), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(check.main(["update", "--root", tmp]), 0)
            self.assertEqual((root / check.SNAPSHOT).read_text(encoding="utf-8"), "pub mod LOSAT\n")


if __name__ == "__main__":
    unittest.main()
