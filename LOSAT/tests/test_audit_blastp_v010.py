from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT_PATH = Path(__file__).with_name("audit_blastp_v010.py")
SPEC = importlib.util.spec_from_file_location("audit_blastp_v010", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:39-45
# ```c++
# const char* kDfltArgTabularOutputFmt =
#     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
#     "evalue bitscore";
# ```
class BlastpAuditClassifierTests(unittest.TestCase):
    def classify(self, ncbi: str, losat: str) -> tuple[str, dict[str, object]]:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            ncbi_path = root / "ncbi.tsv"
            losat_path = root / "losat.tsv"
            ncbi_path.write_text(ncbi, encoding="utf-8")
            losat_path.write_text(losat, encoding="utf-8")
            return AUDIT.classify_outputs(ncbi_path, losat_path)

    def test_exact_text(self) -> None:
        row = "q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-05\t20.0\n"
        classification, metrics = self.classify(row, row)
        self.assertEqual(classification, "EXACT_TEXT")
        self.assertEqual(metrics["ncbi_rows"], 1)

    def test_exact_data_is_not_exact_text(self) -> None:
        row = "q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-05\t20.0"
        classification, _ = self.classify(row + "\n", row)
        self.assertEqual(classification, "EXACT_DATA")
        self.assertNotIn(classification, AUDIT.EXACT_CLASSIFICATIONS)

    def test_order_only(self) -> None:
        first = "q\ts1\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-05\t20.0\n"
        second = "q\ts2\t75.000\t4\t1\t0\t1\t4\t2\t5\t2e-04\t15.0\n"
        classification, _ = self.classify(first + second, second + first)
        self.assertEqual(classification, "ORDER_ONLY")

    def test_value_diff_with_same_hsp_key(self) -> None:
        ncbi = "q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-05\t20.0\n"
        losat = "q\ts\t75.000\t4\t1\t0\t1\t4\t1\t4\t2e-05\t19.0\n"
        classification, _ = self.classify(ncbi, losat)
        self.assertEqual(classification, "VALUE_DIFF")

    def test_missing_hsp_key(self) -> None:
        first = "q\ts1\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-05\t20.0\n"
        second = "q\ts2\t75.000\t4\t1\t0\t1\t4\t2\t5\t2e-04\t15.0\n"
        classification, metrics = self.classify(first + second, first)
        self.assertEqual(classification, "MISSING")
        self.assertEqual(metrics["missing_key_count"], 1)


class BlastpAuditManifestTests(unittest.TestCase):
    def test_manifest_is_the_nine_case_gbdraw_matrix(self) -> None:
        repo_root = SCRIPT_PATH.parents[2]
        cases = AUDIT.load_manifest(
            SCRIPT_PATH.with_name("blastp_v010_parity_manifest.tsv"),
            repo_root,
        )
        self.assertEqual(len(cases), 9)
        self.assertEqual(
            {profile for case in cases for profile in case.profile_ids.split(",")},
            {"P1", "P2", "P3"},
        )
        self.assertEqual({case.outfmt for case in cases}, {"6"})
        self.assertTrue(any(case.max_hsps_per_subject == 1 for case in cases))
        self.assertTrue(any(case.max_hsps_per_subject is None for case in cases))
        self.assertTrue(any(case.max_target_seqs == 1 for case in cases))
        self.assertTrue(any(case.max_target_seqs == 5 for case in cases))
        self.assertTrue(any(case.max_target_seqs is None for case in cases))
        self.assertEqual({case.num_threads for case in cases}, {1, 4})


if __name__ == "__main__":
    unittest.main()
