from __future__ import annotations

from collections import Counter
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT_PATH = Path(__file__).with_name("audit_tblastx_v010.py")
SPEC = importlib.util.spec_from_file_location("audit_tblastx_v010", SCRIPT_PATH)
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
class TblastxAuditClassifierTests(unittest.TestCase):
    def classify(self, ncbi: str, losat: str) -> tuple[str, dict[str, object]]:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            ncbi_path = root / "ncbi.tsv"
            losat_path = root / "losat.tsv"
            ncbi_path.write_text(ncbi, encoding="utf-8")
            losat_path.write_text(losat, encoding="utf-8")
            return AUDIT.classify_outputs(ncbi_path, losat_path)

    def test_exact_text(self) -> None:
        row = "q\ts\t100.000\t4\t0\t0\t1\t12\t1\t12\t1e-05\t20.0\n"
        classification, metrics = self.classify(row, row)
        self.assertEqual(classification, "EXACT_TEXT")
        self.assertEqual(metrics["ncbi_rows"], 1)

    def test_exact_data_is_not_exact_text(self) -> None:
        row = "q\ts\t100.000\t4\t0\t0\t1\t12\t1\t12\t1e-05\t20.0"
        classification, _ = self.classify(row + "\n", row)
        self.assertEqual(classification, "EXACT_DATA")
        self.assertNotIn(classification, AUDIT.EXACT_CLASSIFICATIONS)

    def test_order_only(self) -> None:
        first = "q\ts1\t100.000\t4\t0\t0\t1\t12\t1\t12\t1e-05\t20.0\n"
        second = "q\ts2\t75.000\t4\t1\t0\t1\t12\t4\t15\t2e-04\t15.0\n"
        classification, _ = self.classify(first + second, second + first)
        self.assertEqual(classification, "ORDER_ONLY")

    def test_value_diff_with_same_hsp_key(self) -> None:
        ncbi = "q\ts\t100.000\t4\t0\t0\t1\t12\t1\t12\t1e-05\t20.0\n"
        losat = "q\ts\t75.000\t4\t1\t0\t1\t12\t1\t12\t2e-05\t19.0\n"
        classification, _ = self.classify(ncbi, losat)
        self.assertEqual(classification, "VALUE_DIFF")

    def test_hsp_set_diff(self) -> None:
        first = "q\ts1\t100.000\t4\t0\t0\t1\t12\t1\t12\t1e-05\t20.0\n"
        second = "q\ts2\t75.000\t4\t1\t0\t1\t12\t4\t15\t2e-04\t15.0\n"
        classification, metrics = self.classify(first, second)
        self.assertEqual(classification, "HSP_SET_DIFF")
        self.assertEqual(metrics["missing_key_count"], 1)
        self.assertEqual(metrics["extra_key_count"], 1)


class TblastxAuditContractTests(unittest.TestCase):
    def make_case(self, contract: str) -> object:
        return AUDIT.Case(
            case_id="case",
            profile_ids="P1",
            contract=contract,
            query=Path("query.fa"),
            subject=Path("subject.fa"),
            query_gencode=1,
            db_gencode=1 if contract == AUDIT.PARITY_CONTRACT else 4,
            gencode_args="explicit",
            num_threads=1,
            outfmt="6",
            repeatability_required=True,
            coverage_note="test",
        )

    def test_parity_contract_rejects_sorted_only_equality(self) -> None:
        case = self.make_case(AUDIT.PARITY_CONTRACT)
        self.assertTrue(AUDIT.contract_accepts(case, "EXACT_TEXT"))
        self.assertFalse(AUDIT.contract_accepts(case, "ORDER_ONLY"))

    def test_deviation_contract_accepts_only_declared_hsp_set_difference(self) -> None:
        case = self.make_case(AUDIT.DEVIATION_CONTRACT)
        self.assertTrue(AUDIT.contract_accepts(case, "HSP_SET_DIFF"))
        self.assertFalse(AUDIT.contract_accepts(case, "VALUE_DIFF"))
        self.assertFalse(AUDIT.contract_accepts(case, "EXACT_TEXT"))

    def test_output_directory_must_remain_under_tmp(self) -> None:
        self.assertEqual(AUDIT.validate_output_dir(Path("/tmp/tblastx-audit")), Path("/tmp/tblastx-audit"))
        with self.assertRaises(ValueError):
            AUDIT.validate_output_dir(Path("/var/lib/losat-tblastx-audit"))


class TblastxAuditManifestTests(unittest.TestCase):
    def test_manifest_is_the_narrow_twenty_case_matrix(self) -> None:
        repo_root = SCRIPT_PATH.parents[2]
        cases = AUDIT.load_manifest(
            SCRIPT_PATH.with_name("tblastx_v010_parity_manifest.tsv"),
            repo_root,
        )
        self.assertEqual(len(cases), 20)
        self.assertEqual(
            Counter(case.contract for case in cases),
            Counter({AUDIT.PARITY_CONTRACT: 14, AUDIT.DEVIATION_CONTRACT: 6}),
        )
        self.assertEqual(
            {profile for case in cases for profile in case.profile_ids.split(",")},
            {"P1", "P2"},
        )
        self.assertEqual({case.outfmt for case in cases}, {"6"})
        self.assertEqual({case.num_threads for case in cases}, {1})
        self.assertEqual(sum(case.repeatability_required for case in cases), 15)
        self.assertEqual(sum(case.gencode_args == "implicit" for case in cases), 1)
        self.assertTrue(
            any(case.query_gencode == 4 and case.db_gencode == 1 for case in cases)
        )
        self.assertTrue(
            any(case.query_gencode == 1 and case.db_gencode == 4 for case in cases)
        )


if __name__ == "__main__":
    unittest.main()
