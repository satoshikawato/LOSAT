#!/usr/bin/env python3
"""Focused tests for the pure-Rust runtime boundary checker.

NCBI references:
- ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
  qsort(..., ScoreCompareHSPs);
- ncbi-blast/c++/src/algo/blast/core/link_hsps.c:483-486
  qsort(..., s_RevCompareHSPsTbx);

The tests freeze LOSAT's ownership boundary around those external calls; they
do not reproduce or alter NCBI sorting behavior.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_pure_rust_runtime_boundary.py")
SPEC = importlib.util.spec_from_file_location("pure_rust_boundary", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
boundary = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = boundary
SPEC.loader.exec_module(boundary)


def empty_metadata() -> dict[str, object]:
    return {
        "packages": [{"id": "losat", "name": "LOSAT", "version": "0.1.0", "links": None}],
        "resolve": {"root": "losat", "nodes": [{"id": "losat", "deps": []}]},
    }


def write_fixture(root: Path, production: str, oracle: str = "") -> None:
    source_dir = root / "LOSAT" / "src"
    tests_dir = root / "LOSAT" / "tests"
    source_dir.mkdir(parents=True)
    tests_dir.mkdir(parents=True)
    (source_dir / "lib.rs").write_text(production, encoding="utf-8")
    (tests_dir / "oracle.rs").write_text(oracle, encoding="utf-8")
    (root / "LOSAT" / "Cargo.toml").write_text(
        '[package]\nname = "LOSAT"\nversion = "0.1.0"\nedition = "2021"\n',
        encoding="utf-8",
    )


class PureRustBoundaryCheckerTests(unittest.TestCase):
    def test_imported_abi_and_call_are_detected_but_rust_export_is_classified(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fixture(
                root,
                """
use std::process::{Command as ProcessRunner};

pub extern "C" fn losat_export() -> i32 { 7 }
unsafe extern "C" {
    fn host_sort();
}
fn invoke_host() {
    unsafe { host_sort(); }
}
fn launch_oracle_in_production() {
    let _ = ProcessRunner::new("blastn");
}
""",
            )
            report = boundary.audit_repository(root, metadata=empty_metadata())

        keys = {finding.key for finding in report.findings}
        self.assertIn(("rust.imported_abi", "LOSAT/src/lib.rs", "host_sort"), keys)
        self.assertIn(
            ("rust.imported_abi_call", "LOSAT/src/lib.rs", "host_sort@invoke_host"), keys
        )
        self.assertIn(
            (
                "rust.production_subprocess",
                "LOSAT/src/lib.rs",
                "Command::new@launch_oracle_in_production",
            ),
            keys,
        )
        self.assertIn(
            "losat_export", {observation.symbol for observation in report.observations}
        )
        self.assertNotIn(
            ("rust.imported_abi", "LOSAT/src/lib.rs", "losat_export"), keys
        )

    def test_comments_and_explicit_test_oracle_code_are_outside_production_findings(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fixture(
                root,
                """
// NCBI snippet only: unsafe extern "C" { fn qsort(); }
/* std::process::Command::new("blastn"); */
pub fn program_name() -> &'static str { "blastn" }
""",
                'std::process::Command::new("tblastx");\n',
            )
            report = boundary.audit_repository(root, metadata=empty_metadata())

        self.assertEqual(report.findings, [])

    def test_removed_qsort_adapter_names_and_paths_are_zero_state_findings(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fixture(root, "fn native_qsort_blastn_hsps_by() {}\n")
            obsolete_path = (
                root
                / "LOSAT"
                / "src"
                / "algorithm"
                / "tblastx"
                / "ncbi_qsort.rs"
            )
            obsolete_path.parent.mkdir(parents=True)
            obsolete_path.write_text(
                "pub fn pure_rust_placeholder() {}\n", encoding="utf-8"
            )
            report = boundary.audit_repository(root, metadata=empty_metadata())

        keys = {finding.key for finding in report.findings}
        self.assertIn(
            (
                "rust.obsolete_qsort_adapter",
                "LOSAT/src/lib.rs",
                "native_qsort_blastn_hsps_by",
            ),
            keys,
        )
        self.assertIn(
            (
                "path.obsolete_qsort_adapter",
                "LOSAT/src/algorithm/tblastx/ncbi_qsort.rs",
                "ncbi_qsort.rs",
            ),
            keys,
        )

    def test_project_build_script_is_a_review_signal(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fixture(root, "pub fn pure_rust() {}\n")
            (root / "LOSAT" / "build.rs").write_text(
                'fn main() { println!("cargo:rerun-if-changed=build.rs"); }\n',
                encoding="utf-8",
            )
            report = boundary.audit_repository(root, metadata=empty_metadata())

        self.assertIn(
            ("build.project_script", "LOSAT/build.rs", "build.rs"),
            {finding.key for finding in report.findings},
        )

    def test_renamed_native_dependency_is_detected_with_renamed_direct_call(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fixture(
                root,
                "fn invoke_host_qsort() { unsafe { host_c::qsort(); } }\n",
            )
            (root / "LOSAT" / "Cargo.toml").write_text(
                """[package]
name = "LOSAT"
version = "0.1.0"
edition = "2021"

[dependencies]
host_c = { package = "libc", version = "0.2" }
""",
                encoding="utf-8",
            )
            report = boundary.audit_repository(root, metadata=empty_metadata())

        self.assertIn(
            (
                "cargo.native_dependency",
                "LOSAT/Cargo.toml",
                "host_c->package=libc@dependencies",
            ),
            {finding.key for finding in report.findings},
        )

    def test_unknown_and_stale_entries_fail_exact_production_inventory(self) -> None:
        finding = boundary.Finding(
            "rust.imported_abi",
            "LOSAT/src/lib.rs",
            "host_sort",
            3,
            "production_native_boundary_candidate",
            "fixture",
        )
        exact = boundary.AllowlistEntry(
            finding.rule,
            finding.path,
            finding.symbol,
            boundary.TEMPORARY_DELEGATION,
        )
        self.assertTrue(boundary.evaluate_findings([finding], [exact]).passed)

        unknown = boundary.evaluate_findings([finding], [])
        self.assertFalse(unknown.passed)
        self.assertEqual(unknown.unexpected, [finding])

        stale_entry = boundary.AllowlistEntry(
            "rust.imported_abi",
            "LOSAT/src/removed.rs",
            "old_sort",
            boundary.TEMPORARY_DELEGATION,
        )
        stale = boundary.evaluate_findings(
            [finding], [exact, stale_entry]
        )
        self.assertFalse(stale.passed)
        self.assertEqual(stale.stale, [stale_entry])
        self.assertEqual(stale.invalid_allowlist, [])

        expanded_finding = boundary.Finding(
            "rust.imported_abi",
            "LOSAT/src/new_adapter.rs",
            "host_sort_2",
            7,
            "production_native_boundary_candidate",
            "fixture",
        )
        expanded_entry = boundary.AllowlistEntry(
            expanded_finding.rule,
            expanded_finding.path,
            expanded_finding.symbol,
            boundary.TEMPORARY_DELEGATION,
        )
        expanded = boundary.evaluate_findings(
            [finding, expanded_finding], [exact, expanded_entry]
        )
        self.assertTrue(expanded.passed)

        finalized = boundary.evaluate_findings(
            [finding], [exact], require_empty_allowlist=True
        )
        self.assertFalse(finalized.passed)
        self.assertEqual(finalized.unexpected, [finding])
        self.assertEqual(len(finalized.invalid_allowlist), 1)

    def test_cargo_links_metadata_uses_separate_exact_review_inventory(self) -> None:
        metadata = {
            "packages": [
                {"id": "losat", "name": "LOSAT", "version": "0.1.0", "links": None},
                {
                    "id": "native-sys",
                    "name": "native-sys",
                    "version": "1.2.3",
                    "links": "native",
                },
            ],
            "resolve": {
                "root": "losat",
                "nodes": [
                    {
                        "id": "losat",
                        "deps": [
                            {
                                "pkg": "native-sys",
                                "dep_kinds": [{"kind": None, "target": None}],
                            }
                        ],
                    },
                    {"id": "native-sys", "deps": []},
                ],
            },
        }
        findings = boundary._scan_cargo_links(metadata)
        self.assertEqual(len(findings), 1)
        finding = findings[0]
        self.assertEqual(finding.rule, "cargo.links_metadata")
        self.assertEqual(finding.classification, "review_required_metadata_not_delegation")
        self.assertFalse(
            boundary.evaluate_findings(
                findings,
                [],
                expected_classification=boundary.REVIEWED_NON_ALGORITHM,
            ).passed
        )

        reviewed = boundary.AllowlistEntry(
            finding.rule,
            finding.path,
            finding.symbol,
            boundary.REVIEWED_NON_ALGORITHM,
        )
        self.assertTrue(
            boundary.evaluate_findings(
                findings,
                [reviewed],
                expected_classification=boundary.REVIEWED_NON_ALGORITHM,
            ).passed
        )

        delegation = boundary.AllowlistEntry(
            finding.rule,
            finding.path,
            finding.symbol,
            boundary.TEMPORARY_DELEGATION,
        )
        self.assertFalse(
            boundary.evaluate_findings(
                findings,
                [delegation],
                expected_classification=boundary.REVIEWED_NON_ALGORITHM,
            ).passed
        )
        self.assertFalse(boundary.evaluate_findings(findings, [reviewed]).passed)

        changed = boundary.Finding(
            finding.rule,
            finding.path,
            "native-sys@1.2.4|links=native|kinds=normal",
            finding.line,
            finding.classification,
            finding.detail,
        )
        changed_evaluation = boundary.evaluate_findings(
            [changed],
            [reviewed],
            expected_classification=boundary.REVIEWED_NON_ALGORITHM,
        )
        self.assertFalse(changed_evaluation.passed)
        self.assertEqual(changed_evaluation.unexpected, [changed])
        self.assertEqual(changed_evaluation.stale, [reviewed])


if __name__ == "__main__":
    unittest.main()
