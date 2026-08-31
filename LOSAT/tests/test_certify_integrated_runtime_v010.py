#!/usr/bin/env python3
"""Focused contract tests for integrated native/serial-Wasm certification."""

from __future__ import annotations

from collections import Counter
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


# NCBI reference:
# ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-832
# The integrated gate delegates raw-output acceptance to the existing program
# authorities. These tests freeze orchestration counts and command identity;
# they do not add a second output classifier.
SCRIPT_PATH = Path(__file__).with_name("certify_integrated_runtime_v010.py")
SPEC = importlib.util.spec_from_file_location("integrated_runtime_cert", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
integrated = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = integrated
SPEC.loader.exec_module(integrated)


class IntegratedRuntimeCertificationTests(unittest.TestCase):
    def test_exact_matrix_and_repeatability_shape_is_frozen(self) -> None:
        self.assertEqual(integrated.EXPECTED_NATIVE_COUNTS, {"blastn": 14, "blastp": 9, "tblastx": 20})
        self.assertEqual(integrated.EXPECTED_WASM_COUNTS, {"blastn": 14, "blastp": 7, "tblastx": 20})
        self.assertEqual(len(integrated.BLASTP_THREAD_PAIRS), 2)
        self.assertEqual(len(integrated.REPEATABILITY_REPRESENTATIVES), 6)
        self.assertEqual(
            {row[2] for row in integrated.REPEATABILITY_REPRESENTATIVES},
            {
                "ordinary_exact",
                "equal_hsp",
                "representative",
                "linking_heavy_exact",
                "query_gencode_exact",
                "db_gencode_deviation",
            },
        )

    def test_existing_program_authorities_supply_the_exact_case_sets(self) -> None:
        repo_root = SCRIPT_PATH.parents[2]
        tests_dir = repo_root / "LOSAT" / "tests"
        blastn = integrated.load_authority(
            "test_integrated_blastn", tests_dir / "compare_blastn_parity.py"
        )
        blastp = integrated.load_authority(
            "test_integrated_blastp", tests_dir / "audit_blastp_v010.py"
        )
        tblastx = integrated.load_authority(
            "test_integrated_tblastx", tests_dir / "audit_tblastx_v010.py"
        )
        blastn_rows = blastn.read_manifest(tests_dir / "blastn_parity_manifest.tsv")
        blastp_cases = blastp.load_manifest(
            tests_dir / "blastp_v010_parity_manifest.tsv", repo_root
        )
        tblastx_cases = tblastx.load_manifest(
            tests_dir / "tblastx_v010_parity_manifest.tsv", repo_root
        )

        self.assertEqual(len(blastn_rows), 14)
        self.assertEqual(len(blastp_cases), 9)
        self.assertEqual(sum(case.num_threads == 1 for case in blastp_cases), 7)
        self.assertEqual(len(tblastx_cases), 20)
        self.assertEqual(
            Counter(case.contract for case in tblastx_cases),
            Counter({"parity": 14, "approved_db_gencode_deviation": 6}),
        )

    def test_wasm_command_preserves_losat_arguments_and_output_replacement(self) -> None:
        native = ["/tmp/LOSAT", "blastn", "-q", "/tmp/q.fa", "-o", "/tmp/a.out"]
        wasm = integrated.wasm_command(
            "node", Path("/repo/tests/run.js"), Path("/tmp/LOSAT.wasm"), native
        )
        self.assertEqual(
            wasm,
            [
                "node",
                "/repo/tests/run.js",
                "/tmp/LOSAT.wasm",
                "blastn",
                "-q",
                "/tmp/q.fa",
                "-o",
                "/tmp/a.out",
            ],
        )
        self.assertEqual(
            integrated.replace_output(wasm, "-o", Path("/tmp/run2.out"))[-1],
            "/tmp/run2.out",
        )

    def test_output_directory_must_be_fresh_and_under_tmp(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fresh = root / "fresh"
            self.assertEqual(integrated.validate_output_dir(fresh), fresh.resolve())
            (fresh / "existing.txt").write_text("evidence\n", encoding="utf-8")
            with self.assertRaises(integrated.CertificationFailure):
                integrated.validate_output_dir(fresh)
        with self.assertRaises(integrated.CertificationFailure):
            integrated.validate_output_dir(Path("/var/tmp/not-accepted"))

    def test_classification_counters_are_exact_not_minimums(self) -> None:
        integrated.assert_counter(
            Counter({"EXACT_TEXT": 14, "HSP_SET_DIFF": 6}),
            {"EXACT_TEXT": 14, "HSP_SET_DIFF": 6},
            "TBLASTX",
        )
        with self.assertRaises(integrated.CertificationFailure):
            integrated.assert_counter(
                Counter({"EXACT_TEXT": 15, "HSP_SET_DIFF": 5}),
                {"EXACT_TEXT": 14, "HSP_SET_DIFF": 6},
                "TBLASTX",
            )


if __name__ == "__main__":
    unittest.main()
