#!/usr/bin/env python3
"""Synthetic tests for the BLASTN v0.1.0 certification contract."""

from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("certify_blastn_v010.py")
SPEC = importlib.util.spec_from_file_location("certify_blastn_v010", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
certification = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = certification
SPEC.loader.exec_module(certification)

EXCEPTION_CASE = "Sakai.MG1655.megablast"


# NCBI reference: ncbi-blast-2.17.0+-src/c++/src/objtools/align_format/tabular.cpp:1100-1108
# ```c
# ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
#     if (iter != m_FieldsToShow.begin())
#         m_Ostream << m_FieldDelimiter;
#     x_PrintField(*iter);
# }
# m_Ostream << "\n";
# ```
def tabular_row(
    qstart: int,
    *,
    query: str = "query",
    pident: str = "90.000",
    length: int = 11,
    mismatch: int = 1,
    gapopen: int = 1,
    evalue: str = "1e-05",
    bitscore: str = "20.0",
) -> bytes:
    qend = qstart + 10
    return (
        f"{query}\tsubject\t{pident}\t{length}\t{mismatch}\t{gapopen}\t"
        f"{qstart}\t{qend}\t{qstart}\t{qend}\t{evalue}\t{bitscore}\n"
    ).encode()


def allowed_residual() -> tuple[bytes, bytes]:
    ncbi = [b"# exact header\n"]
    losat = [b"# exact header\n"]
    for index in range(5):
        qstart = 1 + index * 20
        ncbi.append(tabular_row(qstart))
        losat.append(
            tabular_row(
                qstart,
                pident="89.000" if index < 2 else "90.000",
                length=12 if index < 2 else 11,
                mismatch=2 if index < 2 else 1,
                gapopen=2,
            )
        )
    return b"".join(ncbi), b"".join(losat)


class CertificationFixture:
    def __init__(self, root: Path, case_ids: list[str], exception: bool = False) -> None:
        self.root = root
        self.manifest = root / "manifest.tsv"
        self.outputs = root / "paired"
        self.outputs.mkdir()
        self.exceptions = root / "exceptions.tsv"
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["case_id"])
            for case_id in case_ids:
                writer.writerow([case_id])
        with self.exceptions.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(certification.EXCEPTION_COLUMNS)
            if exception:
                writer.writerow(
                    [
                        EXCEPTION_CASE,
                        "pident,length,mismatch,gapopen",
                        "5",
                        "2",
                        "2",
                        "2",
                        "5",
                    ]
                )

    def write_pair(self, case_id: str, ncbi: bytes, losat: bytes) -> None:
        ncbi_path, losat_path = certification.paired_paths(self.outputs, case_id)
        ncbi_path.write_bytes(ncbi)
        losat_path.write_bytes(losat)

    def certify(self, *, exception: bool = False):
        expected_cases = (
            certification.EXPECTED_SOURCE_EXCEPTION_CASES
            if exception
            else frozenset()
        )
        return certification.certify_suite(
            self.manifest,
            self.outputs,
            self.exceptions,
            expected_manifest_case_count=len(certification.read_manifest_case_ids(self.manifest)),
            expected_exception_cases=expected_cases,
        )


class CertifyBlastnV010Tests(unittest.TestCase):
    def test_all_exact_suite_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), ["exact.one", "exact.two"])
            fixture.write_pair("exact.one", tabular_row(1), tabular_row(1))
            fixture.write_pair("exact.two", b"# exact\n", b"# exact\n")
            results = fixture.certify()
        self.assertEqual([result.classification for result in results], ["EXACT_TEXT"] * 2)

    def test_declared_source_underdetermined_residual_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), [EXCEPTION_CASE], exception=True)
            fixture.write_pair(EXCEPTION_CASE, *allowed_residual())
            result = fixture.certify(exception=True)[0]
        self.assertEqual(result.classification, "SOURCE_UNDETERMINED_ACCEPTED")
        self.assertEqual(result.differing_coord_keys, 5)
        self.assertEqual(
            result.field_diffs,
            {"gapopen": 5, "length": 2, "mismatch": 2, "pident": 2},
        )

    def test_undeclared_value_diff_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), ["exact"])
            fixture.write_pair("exact", tabular_row(1), tabular_row(1, pident="89.000"))
            with self.assertRaisesRegex(certification.CertificationError, "not byte-exact"):
                fixture.certify()

    def test_coordinate_key_difference_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), [EXCEPTION_CASE], exception=True)
            ncbi, losat = allowed_residual()
            losat = losat.replace(tabular_row(1, pident="89.000", length=12, mismatch=2, gapopen=2), tabular_row(2, pident="89.000", length=12, mismatch=2, gapopen=2), 1)
            fixture.write_pair(EXCEPTION_CASE, ncbi, losat)
            with self.assertRaisesRegex(certification.CertificationError, "unexpected qstart difference"):
                fixture.certify(exception=True)

    def test_row_count_difference_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), [EXCEPTION_CASE], exception=True)
            ncbi, losat = allowed_residual()
            fixture.write_pair(EXCEPTION_CASE, ncbi + tabular_row(201), losat)
            with self.assertRaisesRegex(certification.CertificationError, "physical line count differs"):
                fixture.certify(exception=True)

    def test_evalue_difference_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), [EXCEPTION_CASE], exception=True)
            ncbi, losat = allowed_residual()
            losat = losat.replace(b"1e-05", b"2e-05", 1)
            fixture.write_pair(EXCEPTION_CASE, ncbi, losat)
            with self.assertRaisesRegex(certification.CertificationError, "unexpected evalue difference"):
                fixture.certify(exception=True)

    def test_bitscore_difference_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), [EXCEPTION_CASE], exception=True)
            ncbi, losat = allowed_residual()
            losat = losat.replace(b"20.0\n", b"21.0\n", 1)
            fixture.write_pair(EXCEPTION_CASE, ncbi, losat)
            with self.assertRaisesRegex(certification.CertificationError, "unexpected bitscore difference"):
                fixture.certify(exception=True)

    def test_unexpected_field_difference_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), [EXCEPTION_CASE], exception=True)
            ncbi, losat = allowed_residual()
            losat = losat.replace(b"20.0\n", b"20.0\textra\n", 1)
            fixture.write_pair(EXCEPTION_CASE, ncbi, losat)
            with self.assertRaisesRegex(certification.CertificationError, "expected exactly 12 fields"):
                fixture.certify(exception=True)

    def test_expanded_differing_row_footprint_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), [EXCEPTION_CASE], exception=True)
            ncbi, losat = allowed_residual()
            fixture.write_pair(
                EXCEPTION_CASE,
                ncbi + tabular_row(201),
                losat + tabular_row(201, gapopen=2),
            )
            with self.assertRaisesRegex(certification.CertificationError, "footprint changed"):
                fixture.certify(exception=True)

    def test_missing_case_output_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = CertificationFixture(Path(tmp), ["missing"])
            with self.assertRaisesRegex(certification.CertificationError, "missing paired output"):
                fixture.certify()


if __name__ == "__main__":
    unittest.main()
