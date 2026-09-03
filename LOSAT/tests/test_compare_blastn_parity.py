#!/usr/bin/env python3
"""Focused, install-free tests for compare_blastn_parity.py."""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import io
import json
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("compare_blastn_parity.py")
SPEC = importlib.util.spec_from_file_location("compare_blastn_parity", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
comparison = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = comparison
SPEC.loader.exec_module(comparison)


# NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
# ```c
# ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
#     if (iter != m_FieldsToShow.begin())
#         m_Ostream << m_FieldDelimiter;
#     x_PrintField(*iter);
# }
# m_Ostream << "\n";
# ```
def tabular_row(qstart: int, *, pident: str = "100.000") -> str:
    qend = qstart + 10
    return f"query\tsubject\t{pident}\t11\t0\t0\t{qstart}\t{qend}\t{qstart}\t{qend}\t1e-05\t20.0\n"


# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-107,157-168,208
# ```c
# const string kArgQuery("query");
# const string kArgOutput("out");
# const string kArgSubject("subject");
# const string kTask("task");
# const string kArgNumThreads("num_threads");
# const string kArgEvalue("evalue");
# const string kArgWordSize("word_size");
# const string kArgDustFiltering("dust");
# const string kArgLookupTableMaskingOnly("soft_masking");
# const string kArgMaxHSPsPerSubject("max_hsps");
# ```
def manifest_row(case_id: str = "paired.case") -> dict[str, str]:
    return {
        "case_id": case_id,
        "task": "blastn",
        "reward": "2",
        "penalty": "-3",
        "gap_open": "5",
        "gap_extend": "2",
        "word_size": "11",
        "window_size": "0",
        "scan_range": "0",
        "dust": "true",
        "soft_masking": "true",
        "strand": "both",
        "evalue": "10",
        "max_target_seqs": "500",
        "max_hsps_per_subject": "0",
        "num_threads": "1",
        "outfmt": "6",
        "query": "query.fasta",
        "subject": "subject.fasta",
        "ncbi_out": "-",
        "losat_out": "-",
    }


class CompareBlastnParityTests(unittest.TestCase):
    # NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
    # The final newline and every preceding byte are written directly to the
    # output stream, so a byte-only gate must independently fail on comments.
    def test_fail_on_byte_diff_controls_exit_status_independently(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            temp_dir = Path(tmp)
            ncbi = temp_dir / "ncbi.out"
            losat = temp_dir / "losat.out"
            ncbi.write_text("# NCBI header\n" + tabular_row(1), encoding="utf-8")
            losat.write_text("# LOSAT header\n" + tabular_row(1), encoding="utf-8")

            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--ncbi",
                    str(ncbi),
                    "--losat",
                    str(losat),
                    "--fail-on-byte-diff",
                ],
                check=False,
                capture_output=True,
                text=True,
            )

        self.assertEqual(completed.returncode, 1, completed.stdout + completed.stderr)
        self.assertIn("raw_bytes: DIFFER", completed.stdout)
        self.assertIn("structured_equal: True", completed.stdout)

    # NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:828-832
    # ```c
    # ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
    #     tabinfo.SetFields(**itr, *m_Scope, &m_ScoringMatrix);
    #     tabinfo.Print();
    # }
    # ```
    def test_structured_equality_detects_row_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            temp_dir = Path(tmp)
            ncbi = temp_dir / "ncbi.out"
            losat = temp_dir / "losat.out"
            ncbi.write_text(tabular_row(1) + tabular_row(21), encoding="utf-8")
            losat.write_text(tabular_row(21) + tabular_row(1), encoding="utf-8")
            with (
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()),
            ):
                result = comparison.compare_case("order", ncbi, losat, 5)

        self.assertTrue(result["row_count_equal"])
        self.assertTrue(result["row_multiset_equal"])
        self.assertFalse(result["row_order_equal"])
        self.assertFalse(result["structured_equal"])
        self.assertFalse(result["ncbi_only"])
        self.assertFalse(result["losat_only"])

    # NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:828-832
    # Each CSeq_align_set entry is printed; coordinate-key multiplicity cannot
    # be discarded even when duplicate rows have identical values.
    def test_structured_equality_detects_duplicate_multiplicity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            temp_dir = Path(tmp)
            ncbi = temp_dir / "ncbi.out"
            losat = temp_dir / "losat.out"
            ncbi.write_text(tabular_row(1) + tabular_row(21) * 2, encoding="utf-8")
            losat.write_text(tabular_row(1) + tabular_row(21), encoding="utf-8")
            with (
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()),
            ):
                result = comparison.compare_case("duplicates", ncbi, losat, 5)

        self.assertFalse(result["row_count_equal"])
        self.assertFalse(result["row_multiset_equal"])
        self.assertFalse(result["structured_equal"])
        self.assertEqual(
            sum(count - 1 for count in result["ncbi_duplicates"].values()), 1
        )
        self.assertFalse(result["ncbi_only"])
        self.assertFalse(result["losat_only"])

    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-107
    # The fake executable only records the already-defined NCBI CLI boundary;
    # it is not an alternate BLAST implementation or an external dependency.
    def test_fresh_paired_mode_uses_same_cwd_and_manifest_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fake_executable = root / "fake_blast.py"
            fake_executable.write_text(
                textwrap.dedent(
                    """\
                    #!/usr/bin/env python3
                    import json
                    import os
                    import sys
                    from pathlib import Path

                    args = sys.argv[1:]
                    output_flag = "-out" if "-out" in args else "-o"
                    output = Path(args[args.index(output_flag) + 1])
                    output.write_bytes(b"query\\tsubject\\t100.000\\t11\\t0\\t0\\t1\\t11\\t1\\t11\\t1e-05\\t20.0\\n")
                    Path(str(output) + ".invocation.json").write_text(
                        json.dumps({"cwd": os.getcwd(), "args": args}), encoding="utf-8"
                    )
                    """
                ),
                encoding="utf-8",
            )
            fake_executable.chmod(0o755)
            manifest = root / "manifest.tsv"
            row = manifest_row()
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(row), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerow(row)
            output_dir = root / "paired-output"

            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--manifest",
                    str(manifest),
                    "--fresh-paired",
                    "--paired-output-dir",
                    str(output_dir),
                    "--ncbi-bin",
                    str(fake_executable),
                    "--losat-bin",
                    str(fake_executable),
                    "--fail-on-byte-diff",
                    "--fail-on-diff",
                ],
                cwd=root,
                check=False,
                capture_output=True,
                text=True,
            )

            self.assertEqual(
                completed.returncode, 0, completed.stdout + completed.stderr
            )
            ncbi_invocation = json.loads(
                (output_dir / "paired.case.ncbi.out.invocation.json").read_text(
                    encoding="utf-8"
                )
            )
            losat_invocation = json.loads(
                (output_dir / "paired.case.losat.out.invocation.json").read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(ncbi_invocation["cwd"], str(root))
        self.assertEqual(losat_invocation["cwd"], str(root))
        self.assertEqual(
            ncbi_invocation["args"][ncbi_invocation["args"].index("-query") + 1],
            row["query"],
        )
        self.assertEqual(
            losat_invocation["args"][losat_invocation["args"].index("-q") + 1],
            row["query"],
        )
        self.assertEqual(
            ncbi_invocation["args"][ncbi_invocation["args"].index("-subject") + 1],
            row["subject"],
        )
        self.assertEqual(
            losat_invocation["args"][losat_invocation["args"].index("-s") + 1],
            row["subject"],
        )
        self.assertIn("raw_bytes: exact", completed.stdout)

    # NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1275-1283,1322-1325
    # ```c
    # m_Ostream << "# " << num_hits << " hits found" << "\n";
    # m_Ostream << "# BLAST processed " << num_queries << " queries\n";
    # ```
    def test_permanent_manifest_covers_formats_no_hit_and_multi_query(self) -> None:
        manifest = SCRIPT.with_name("blastn_parity_manifest.tsv")
        fixture = SCRIPT.with_name("fasta") / "blastn_parity_compact.fasta"
        rows = [
            row
            for row in comparison.read_manifest(manifest)
            if row["case_id"].startswith("compact.")
        ]
        formats = {row["outfmt"] for row in rows}
        no_hit_rows = [row for row in rows if ".no_hit." in row["case_id"]]
        query_count = sum(
            1
            for line in fixture.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        )

        self.assertEqual(formats, {"6", "7"})
        self.assertTrue(no_hit_rows)
        self.assertTrue(all(row["ncbi_out"] == row["losat_out"] == "-" for row in rows))
        self.assertGreaterEqual(query_count, 2)


if __name__ == "__main__":
    unittest.main()
