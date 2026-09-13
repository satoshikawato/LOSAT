"""Regression checks for benchmark data integrity, not BLAST algorithms."""
from contextlib import redirect_stdout
import io
import json
import sys
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import comparison_data as data
import plot_execution_time as timing
import plot_overall_trend as overall


# NCBI reference: c++/src/objtools/align_format/format_flags.cpp:38-40
# "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
# "evalue bitscore";
ROW = "q\ts\t100.000\t30\t0\t0\t1\t30\t1\t30\t1e-10\t50\n"


class SimpleComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "records"
        data.record_cli(["init", str(self.root)])
        self.patch = patch.object(data, "RESULT_DIR", self.root)
        self.patch.start()
        self.addCleanup(self.patch.stop)
        with patch.dict("os.environ", {"BENCHMARK_PROGRAMS": "tblastx", "BENCHMARK_CASE": ""}):
            self.cases = data.comparison_cases()[:2]

    def write(self, case, tools, row=ROW):
        for tool, path in data.result_paths(case).items():
            if tool in tools:
                path.parent.mkdir(parents=True, exist_ok=True)
                self.record(path, row)

    # NCBI reference: c++/src/objtools/align_format/tabular.cpp:1100-1108
    # x_PrintField(*iter); m_Ostream << "\n";
    # Synthetic records test evidence integrity, never generate search expectations.
    def record(self, path, row=ROW, log="real\t0.100\nexit_status=0\n"):
        path.with_suffix(".run.json").unlink(missing_ok=True)
        path.unlink(missing_ok=True)
        data.record_cli(["start", str(self.root), str(path.with_suffix("")), sys.executable, "-out", str(path)])
        path.write_text(row)
        path.with_suffix(".log").write_text(log)
        data.record_cli(["finish", str(self.root), str(path.with_suffix("")), "0"])


    # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-48
    # const string kArgQuery("query"); const string kArgSubject("subject");
    # The executable used for the oracle search supplies its own version identity.
    def test_selected_oracle_version_is_preserved(self):
        oracle = self.root / 'blastn'
        oracle.write_text(f'#!{sys.executable}\nimport sys\nassert sys.argv[1:] == ["-version"]\nprint("synthetic oracle version")\n')
        oracle.chmod(0o755)
        stem = self.root / 'oracle-result'
        data.record_cli(['start', str(self.root), str(stem), str(oracle), '-out', str(stem) + '.out'])
        record = json.loads(Path(str(stem) + '.run.json').read_text())
        self.assertEqual(record['oracle_version'], 'synthetic oracle version')
        self.assertEqual(record['files'][str(oracle.resolve())], data.sha256(oracle))

    # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
    # const string kArgQuery("query"); const string kArgSubject("subject");
    # The compilation profile changes only the selected TBLASTX host argv.
    def test_turbofan_profile_is_tblastx_only_and_can_be_disabled(self):
        with patch.dict(os.environ, {"NODE_ARGS_JSON": '["--trace-warnings"]'}, clear=True):
            expected = ["--trace-warnings", "--no-liftoff", "--no-wasm-tier-up"]
            self.assertEqual(data.node_args("tblastx"), expected)
            for program in (None, "blastn", "blastp", "megablast"):
                self.assertEqual(data.node_args(program), ["--trace-warnings"])
            with patch.dict(os.environ, {"NODE_TBLASTX_ARGS_JSON": "[]"}):
                self.assertEqual(data.node_args("tblastx"), ["--trace-warnings"])

    # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
    # const string kArgQuery("query"); const string kArgSubject("subject");
    # Reject malformed host options before any search rather than use a fallback.
    def test_tblastx_profile_rejects_malformed_argv(self):
        for value in ('{}', '[1]', '["bad\\u0000arg"]', 'not-json'):
            with self.subTest(value=value), patch.dict(os.environ, {"NODE_TBLASTX_ARGS_JSON": value}, clear=True):
                with self.assertRaises(ValueError):
                    data.node_args("tblastx")

    # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
    # const string kArgQuery("query"); const string kArgSubject("subject");
    # Persist exactly the same argument array that the shell consumes.
    def test_tblastx_profile_manifest_matches_nul_delimited_argv(self):
        directory = self.root / "profile"
        with patch.dict(os.environ, {"RUN_LOSAT_WASM": "0", "NODE_ARGS_JSON": '["--trace-warnings"]'}, clear=True):
            data.record_cli(["init", str(directory)])
        manifest = json.loads((directory / "run.json").read_text())
        flags = (directory / "node-tblastx-args.bin").read_bytes().split(b"\0")[:-1]
        self.assertEqual([x.decode() for x in flags], manifest["node_argv_by_program"]["tblastx"][1:])
        self.assertNotIn("--no-liftoff", manifest["node_argv_by_program"]["blastp"])
        self.assertNotIn("--no-liftoff", manifest["node_argv_by_program"]["blastn"])

    def test_one_thread_has_distinct_wasm_builds(self):
        with patch.object(data, "LOSAT_THREADS", 1):
            paths = data.result_paths(self.cases[0])
        self.assertEqual(len(paths), 4)
        self.assertEqual(len(set(paths.values())), 4)
        self.assertTrue(str(paths[data.WASM_SINGLE]).endswith(".wasm.out"))
        self.assertTrue(str(paths[data.WASM_MULTI]).endswith(".wasm.n1.out"))

    def test_nondefault_threads_find_matching_ncbi_and_native(self):
        with patch.object(data, "LOSAT_THREADS", 4):
            paths = data.result_paths(self.cases[0])
        self.assertTrue(str(paths["BLAST+"]).endswith(".n4.out"))
        self.assertTrue(str(paths[data.NATIVE_MULTI]).endswith(".n4.out"))

    def test_zero_hits_are_present_and_failed_outputs_are_absent(self):
        self.write(self.cases[0], ["BLAST+", data.WASM_SINGLE], row="")
        frames = data.load_case(self.cases[0])
        self.assertEqual(set(frames), {"BLAST+", data.WASM_SINGLE})
        self.assertTrue(all(frame.empty for frame in frames.values()))
        output = data.result_paths(self.cases[0])[data.WASM_SINGLE]
        output.with_suffix(".log").write_text("real\t0m0.100s\nexit_status=1\n")
        self.assertFalse(data.successful_output(output))
        self.assertIsNone(timing.parse_time(output.with_suffix(".log")))

    def test_malformed_output_fails_instead_of_dropping_rows(self):
        self.write(self.cases[0], ["BLAST+", data.WASM_SINGLE])
        self.record(data.result_paths(self.cases[0])[data.WASM_SINGLE], "q\ts\t100\n")
        with self.assertRaisesRegex(ValueError, "Expected 12"):
            data.load_case(self.cases[0])

    def test_wall_time_formats_and_explicit_failure(self):
        log = self.root / "time.log"
        for content, expected in [
            ("real\t2m3.456s\nexit_status=0\n", 123.456),
            ("real 1.25\nexit_status=0\n", 1.25),
            ("real 1.25\n", None),
            ("0.01user 0.00system 1:02:03.50elapsed\nexit_status=0\n", 3723.5),
            ("real\t0m0.001s\nexit_status=7\n", None),
            ("simple_benchmark=1\nreal\t0m0.001s\n", None),
        ]:
            log.write_text(content)
            self.assertEqual(timing.parse_time(log), expected)

    def test_overall_excludes_different_pair_sets(self):
        self.write(self.cases[0], ["BLAST+", data.NATIVE_SINGLE, data.WASM_SINGLE])
        self.write(self.cases[1], ["BLAST+", data.NATIVE_SINGLE])
        captured = []
        def inspect_hist(*args, **kwargs):
            captured.append(kwargs["data"])
        with patch.object(overall, "comparison_cases", return_value=self.cases), \
             patch.object(overall, "PLOT_DIR", self.root / "plots"), \
             patch.object(overall, "OUTPUT_IMAGE", self.root / "plots" / "overall.png"), \
             patch.object(overall.sns, "histplot", side_effect=inspect_hist), \
             patch.object(overall.sns, "scatterplot"), redirect_stdout(io.StringIO()) as output:
            overall.main()
        self.assertIn("1/2 common pairs", output.getvalue())
        self.assertEqual(len(captured), 2)
        for frame in captured:
            self.assertEqual(set(frame["Task"]), {self.cases[0]["name"]})
            self.assertEqual(len(frame), 3)

    def test_old_missing_tampered_and_cross_run_evidence_are_rejected(self):
        self.write(self.cases[0], ["BLAST+"])
        path = data.result_paths(self.cases[0])["BLAST+"]
        self.assertTrue(data.successful_output(path))
        for suffix in [".log", ".run.json"]:
            saved = path.with_suffix(suffix).read_bytes()
            path.with_suffix(suffix).unlink()
            self.assertFalse(data.successful_output(path))
            path.with_suffix(suffix).write_bytes(saved)
        path.write_text(ROW + ROW)
        self.assertFalse(data.successful_output(path))
        self.record(path)
        path.with_suffix(".log").write_text("real\t0.001\n")
        self.assertFalse(data.successful_output(path))
        self.record(path)
        record = json.loads(path.with_suffix(".run.json").read_text())
        record["run_id"] = "another-run"
        path.with_suffix(".run.json").write_text(json.dumps(record))
        self.assertFalse(data.successful_output(path))

    # NCBI reference: c++/src/objtools/align_format/tabular.cpp:1100-1108
    # x_PrintField(*iter); m_Ostream << "\n";
    # Numerically equivalent fields still violate the raw-output timing gate.
    def test_timing_plot_excludes_lexically_different_output(self):
        self.write(self.cases[0], ["BLAST+", data.WASM_SINGLE])
        self.record(data.result_paths(self.cases[0])[data.WASM_SINGLE], ROW.replace("1e-10", "1.0e-10"))
        with patch.object(timing, "comparison_cases", return_value=self.cases[:1]), \
             patch.object(timing, "PLOT_DIR", self.root / "plots"), \
             patch.object(timing, "OUTPUT_IMAGE", self.root / "plots/times.png"), \
             redirect_stdout(io.StringIO()):
            timing.main()
        import pandas as pd
        rows = pd.read_csv(self.root / "plots/execution_times.tsv", sep="\t")
        self.assertEqual(rows["Tool"].tolist(), ["BLAST+"])

    def test_runner_failure_timeout_and_missing_output_are_retained(self):
        binary = self.root / "fake-losat"
        env = dict(os.environ, BENCHMARK_PROGRAMS="blastp", BENCHMARK_CASE="WSSV.PajaWSV",
                   LOSAT_THREADS="4", LOSAT_BIN=str(binary), RUN_NATIVE="1", RUN_NCBI="0",
                   RUN_LOSAT_WASM="0", BENCHMARK_TIMEOUT="1")
        for number, (command, expected) in enumerate([("exit 7", 7), ("sleep 5", 124), ("exit 0", 125)]):
            with self.subTest(command=command):
                directory = self.root / f"run-{number}"
                binary.write_text(f"#!/bin/bash\n{command}\n"); binary.chmod(0o755)
                result = subprocess.run(["bash", str(data.SCRIPT_DIR / "run_comparison.sh")],
                    env={**env, "BENCHMARK_DIR":str(directory)}, capture_output=True, text=True, timeout=15)
                self.assertEqual(result.returncode, expected, result.stderr)
                output = directory / "losat_out/WSSV.PajaWSV.losatp.out"
                self.assertIn(f"exit_status={expected}", output.with_suffix(".log").read_text())
                self.assertFalse(data.successful_output(output))
                self.assertIsNone(timing.parse_time(output.with_suffix(".log")))
                self.assertEqual(json.loads(output.with_suffix(".run.json").read_text())["exit_status"], expected)
                # Retrying into the same directory cannot overwrite failed evidence or adopt old bytes.
                output.write_text(ROW)
                again = subprocess.run(["bash", str(data.SCRIPT_DIR / "run_comparison.sh")],
                    env={**env, "BENCHMARK_DIR":str(directory)}, capture_output=True, text=True, timeout=15)
                self.assertNotEqual(again.returncode, 0)
                self.assertEqual(output.read_text(), ROW)

    def test_runner_prints_each_time_before_starting_next_run(self):
        binary = self.root / "fake-losat"
        # A command's own 'real' line must not be mistaken for the Bash timer.
        binary.write_text("#!/bin/bash\nprintf 'real\\t99999.999\\n'\nfor ((i=1;i<=$#;i++)); do if [[ ${!i} == -out ]]; then j=$((i+1)); : > \"${!j}\"; fi; done\n")
        binary.chmod(0o755)
        env = dict(os.environ, BENCHMARK_DIR=str(self.root / "run"),
                   BENCHMARK_PROGRAMS="blastp", BENCHMARK_CASE="WSSV.PajaWSV",
                   LOSAT_THREADS="4", LOSAT_BIN=str(binary),
                   RUN_NATIVE="1", RUN_NCBI="0", RUN_LOSAT_WASM="0")
        result = subprocess.run(
            ["bash", str(data.SCRIPT_DIR / "run_comparison.sh")],
            env=env, capture_output=True, text=True, timeout=15,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for stem in ["WSSV.PajaWSV.losatp", "WSSV.PajaWSV.losatp.n4"]:
            log = (self.root / "run/losat_out" / f"{stem}.log").read_text()
            seconds = re.findall(r"^real\s+(\d+\.\d+)$", log, re.M)[-1]
            self.assertLess(float(seconds), 15)
            self.assertIn(f"Finished {stem}: {seconds} s", result.stdout)
            self.assertEqual(timing.parse_time(self.root / "run/losat_out" / f"{stem}.log"), float(seconds))
        self.assertLess(result.stdout.index("Finished WSSV.PajaWSV.losatp:"),
                        result.stdout.index("Running WSSV.PajaWSV.losatp.n4"))


if __name__ == "__main__":
    unittest.main()
