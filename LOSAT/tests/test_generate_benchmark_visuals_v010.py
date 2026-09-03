from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile
import unittest

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

import generate_benchmark_visuals_v010 as bench


class BenchmarkVisualizationTests(unittest.TestCase):
    def test_read_tsv_manifest_skips_comment_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.tsv"
            path.write_text(
                "# comment before header\n"
                "case_id\tvalue\n"
                "alpha\t1\n"
                "# inline comment\n"
                "beta\t2\n",
                encoding="utf-8",
            )
            rows = bench.read_tsv_manifest(path)
        self.assertEqual([row["case_id"] for row in rows], ["alpha", "beta"])

    def test_select_representative_cases_uses_declared_ids_only(self) -> None:
        cases = []
        for program in bench.PROGRAMS:
            for case_id in bench.REPRESENTATIVE_CASES[program]:
                cases.append(bench.CaseSpec(program, case_id, "parity", None))
            cases.append(bench.CaseSpec(program, f"{program}-extra", "parity", None))

        selected = bench.select_cases(cases, "representative")
        observed = [(case.program, case.case_id) for case in selected]
        expected = [
            (program, case_id)
            for program in bench.PROGRAMS
            for case_id in bench.REPRESENTATIVE_CASES[program]
        ]
        self.assertEqual(observed, expected)

    def test_select_representative_cases_fails_closed_when_authority_case_missing(self) -> None:
        cases = [
            bench.CaseSpec(program, case_id, "parity", None)
            for program in bench.PROGRAMS
            for case_id in bench.REPRESENTATIVE_CASES[program]
        ]
        cases.pop()
        with self.assertRaisesRegex(ValueError, "representative case missing"):
            bench.select_cases(cases, "representative")

    def test_parse_alignment_rows_ignores_outfmt7_comments_and_crlf(self) -> None:
        case = bench.CaseSpec("blastn", "example", "parity", None)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "output.tsv"
            path.write_bytes(
                b"# BLASTN 2.17.0+\r\n"
                b"q\ts\t99.5\t100\t0\t0\t1\t100\t5\t104\t1e-20\t80.0\r\n"
                b"# 1 hits found\r\n"
            )
            rows = bench.parse_alignment_rows(path, case, "NCBI BLAST+")

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].length, 100)
        self.assertEqual(rows[0].pident, 99.5)
        self.assertEqual(rows[0].fields[0:2], ("q", "s"))

    def test_write_alignment_data_preserves_metadata_and_standard_fields(self) -> None:
        row = bench.AlignmentRow(
            program="blastp",
            case_id="case",
            contract="parity",
            tool="LOSAT",
            fields=("q", "s", "90", "10", "1", "0", "1", "10", "2", "11", "1e-5", "20"),
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "rows.tsv"
            bench.write_alignment_data(path, [row])
            with path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.reader(handle, delimiter="\t"))

        self.assertEqual(rows[0][:4], ["program", "case_id", "contract", "tool"])
        self.assertEqual(rows[1][:4], ["blastp", "case", "parity", "LOSAT"])
        self.assertEqual(rows[1][4:], list(row.fields))

    def test_log_bins_are_positive_and_cover_observed_range(self) -> None:
        bins = bench.log_bins([10, 100, 1000], count=10)
        self.assertEqual(len(bins), 11)
        self.assertTrue(all(value > 0 for value in bins))
        self.assertLessEqual(bins[0], 10)
        self.assertGreaterEqual(bins[-1], 1000)


if __name__ == "__main__":
    unittest.main()
