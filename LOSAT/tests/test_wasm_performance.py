"""Gate negatives; no algorithm emulation or replacement expected output.

NCBI c++/src/objtools/align_format/tabular.cpp:1098-1108:
  x_PrintField(*iter); ... m_Ostream << "\n";
NCBI c++/src/algo/blast/format/blast_format.cpp:794-808:
  dbname = string("User specified sequence set (Input: ") + m_SubjectTag ...;
Synthetic subprocesses exercise the harness only, never certify a search.
"""

import hashlib
import os
from pathlib import Path
import sys
import tempfile
import unittest

import wasm_performance as perf
import summarize_wasm_performance as summary


class GateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)

    def gate(self, expected, observed):
        return perf.raw_gate(
            {
                "status": "PASS",
                "raw_output_sha256": hashlib.sha256(observed).hexdigest(),
            },
            hashlib.sha256(expected).hexdigest(),
        )

    def test_exporter_rejects_failure_after_successful_warmup(self):
        status, failed, median = summary.series_result(
            [{"status": "PASS"}, {"status": "TIMEOUT", "output": "failed.txt"}],
            {"warmup_count": 1, "timed_repetitions": 5},
        )
        self.assertEqual(status, "TIMEOUT")
        self.assertEqual(failed["output"], "failed.txt")
        self.assertIsNone(median)

    def test_exporter_requires_complete_series(self):
        status, _, median = summary.series_result(
            [{"status": "PASS"}],
            {"warmup_count": 1, "timed_repetitions": 5},
        )
        self.assertEqual(status, "NOT_RUN")
        self.assertIsNone(median)

    def test_exporter_median_excludes_warmup_and_diagnostic(self):
        rows = [
            {
                "status": "PASS",
                "repeat": i,
                "timed": 1 <= i < 6,
                "diagnostics_enabled": i == 6,
                "wall_seconds": seconds,
                "peak_rss_bytes": None,
            }
            for i, seconds in enumerate([999, 5, 1, 3, 4, 2, 888])
        ]
        status, _, median = summary.series_result(
            rows, {"warmup_count": 1, "timed_repetitions": 5}
        )
        self.assertEqual(status, "PASS")
        self.assertEqual(median["median_seconds"], 3)
        self.assertIsNone(median["peak_rss_bytes"])

    def test_wrong_hash(self):
        path = self.root / "artifact"
        path.write_bytes(b"wrong")
        with self.assertRaises(perf.GateFailure):
            perf.checked_file(path, "0" * 64)

    def test_missing_artifact_and_oracle(self):
        for name in ("artifact.wasm", "blastn"):
            with self.assertRaises(perf.GateFailure):
                perf.checked_file(self.root / name, "0" * 64)

    def test_header_path_is_not_normalized(self):
        self.assertEqual(
            self.gate(b"# Subject: /fixed/input\na\n", b"# Subject: /other/input\na\n"),
            "BASELINE_PARITY_FAIL",
        )

    def test_pairwise_header_is_not_normalized(self):
        self.assertEqual(
            self.gate(b"BLASTP 2.17.0+\nQuery= x\n", b"BLASTP 2.16.0+\nQuery= x\n"),
            "BASELINE_PARITY_FAIL",
        )

    def test_order_and_newlines_are_normative(self):
        for data in (b"b\na\n", b"a\r\nb\r\n"):
            self.assertEqual(self.gate(b"a\nb\n", data), "BASELINE_PARITY_FAIL")

    def test_serial_threaded_mislabel(self):
        serial = {"imports": [], "exports": [{"name": "_start"}]}
        threaded = {
            "imports": [{"module": "wasi", "name": "thread-spawn"}],
            "exports": [{"name": "_start"}, {"name": "wasi_thread_start"}],
        }
        with self.assertRaises(perf.GateFailure):
            perf.validate_kind("threaded", serial)
        with self.assertRaises(perf.GateFailure):
            perf.validate_kind("serial", threaded)
        perf.validate_kind("serial", serial)
        perf.validate_kind("threaded", threaded)

    def test_library_is_not_command(self):
        with self.assertRaises(perf.GateFailure):
            perf.validate_kind(
                "serial", {"imports": [], "exports": [{"name": "run_pair"}]}
            )

    def execute(self, script, timeout=5):
        return perf.execute(
            [sys.executable, "-c", script, "--out", "unused"],
            self.root,
            self.root / "run",
            os.environ.copy(),
            timeout,
        )

    def test_bad_exit_even_with_expected_output(self):
        result = self.execute(
            "import pathlib,sys;pathlib.Path(sys.argv[-1]).write_bytes(b'');sys.exit(7)"
        )
        self.assertEqual(result["status"], "BAD_EXIT")
        self.assertEqual(
            perf.raw_gate(result, hashlib.sha256(b"").hexdigest()), "BAD_EXIT"
        )

    def test_missing_output(self):
        self.assertEqual(self.execute("pass")["status"], "MISSING_OUTPUT")

    def test_timeout_never_zero_time_pass(self):
        result = self.execute("import time;time.sleep(5)", timeout=0.05)
        self.assertEqual(result["status"], "TIMEOUT")
        self.assertGreater(result["wall_seconds"], 0)

    def test_default_omission_is_preserved(self):
        catalog = perf.authority.load_catalog(
            perf.ROOT, perf.capture(["git", "rev-parse", "HEAD"])
        )
        case = catalog.blastp_cases["pairwise_default_serial"]
        _, argv = catalog.blastp_audit.build_commands(
            case, Path("oracle"), Path("LOSAT"), self.root
        )
        self.assertNotIn("--max-target-seqs", argv)

    def test_readonly_output_directory(self):
        (self.root / "run").mkdir()
        with self.assertRaises(FileExistsError):
            self.execute("pass")

    def test_sorted_diff_is_diagnostic_only(self):
        a, b = self.root / "a", self.root / "b"
        a.write_bytes(b"a\nb\n")
        b.write_bytes(b"b\na\n")
        perf.diagnostic_diff(a, b, self.root)
        self.assertEqual((self.root / "sorted-diagnostic.diff").read_text(), "")
        self.assertTrue((self.root / "ordered.diff").read_text())
        self.assertEqual(
            self.gate(a.read_bytes(), b.read_bytes()), "BASELINE_PARITY_FAIL"
        )

    def test_manifest_argv_tampering_is_rejected(self):
        # Construct the canonical portion without depending on locally built artifacts.
        head = perf.capture(["git", "rev-parse", "HEAD"])
        catalog = perf.authority.load_catalog(perf.ROOT, head)
        cases = []
        oracles = {p: Path(p) for p in perf.ORACLE_HASHES}
        for row in catalog.canonical_rows:
            ncbi, losat = perf.authority._base_commands(
                perf.ROOT,
                catalog,
                row["program"],
                row["case_id"],
                Path("LOSAT"),
                oracles,
                self.root,
            )
            cases.append(
                {
                    **row,
                    "serial_applicable": row["program"] != "blastp"
                    or catalog.blastp_cases[row["case_id"]].num_threads == 1,
                    "oracle_argv": ncbi,
                    "losat_argv": losat,
                }
            )
        manifest = {
            "schema_version": perf.SCHEMA,
            "cwd": str(perf.ROOT),
            "source": {"head": head},
            "cases": cases,
            "artifacts": {"native": {"path": "LOSAT"}},
            "oracle": {p: {"path": p} for p in oracles},
        }
        case = next(
            c for c in cases if c["case_id"] == "compact.multi_query.no_hit.outfmt7"
        )
        case["losat_argv"] += ["--evalue", "1e-180"]
        with self.assertRaisesRegex(
            perf.GateFailure, "ordered argv authority mismatch"
        ):
            perf.preflight(manifest)


if __name__ == "__main__":
    unittest.main()
