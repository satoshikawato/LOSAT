"""Regression checks for benchmark data integrity, not BLAST algorithms."""
from contextlib import redirect_stdout
import io
import json
import sys
import os
from pathlib import Path
import re
import shutil
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


# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:50,75
# const string kArgDb("db"); const string kArgNumThreads("num_threads");
# Both recorded DB thread variants share the n1 raw-output reference.
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
        with patch.dict(os.environ, {"RUN_LOSAT_WASM": "0", "RUN_LOSAT_WASM_THREADED": "0", "NODE_ARGS_JSON": '["--trace-warnings"]'}, clear=True):
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
        self.assertTrue(str(paths[data.NCBI_SINGLE]).endswith(".n1.out"))
        self.assertTrue(str(paths[data.NCBI_MULTI]).endswith(".n4.out"))
        self.assertTrue(str(paths[data.NATIVE_MULTI]).endswith(".n4.out"))

    def test_zero_hits_are_present_and_failed_outputs_are_absent(self):
        self.write(self.cases[0], [data.NCBI_SINGLE, data.WASM_SINGLE], row="")
        frames = data.load_case(self.cases[0])
        self.assertEqual(set(frames), {data.NCBI_SINGLE, data.WASM_SINGLE})
        self.assertTrue(all(frame.empty for frame in frames.values()))
        output = data.result_paths(self.cases[0])[data.WASM_SINGLE]
        output.with_suffix(".log").write_text("real\t0m0.100s\nexit_status=1\n")
        self.assertFalse(data.successful_output(output))
        self.assertIsNone(timing.parse_time(output.with_suffix(".log")))

    def test_malformed_output_fails_instead_of_dropping_rows(self):
        self.write(self.cases[0], [data.NCBI_SINGLE, data.WASM_SINGLE])
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
        self.write(self.cases[0], [data.NCBI_SINGLE, data.NATIVE_SINGLE, data.WASM_SINGLE])
        self.write(self.cases[1], [data.NCBI_SINGLE, data.NATIVE_SINGLE])
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
        self.write(self.cases[0], [data.NCBI_SINGLE])
        path = data.result_paths(self.cases[0])[data.NCBI_SINGLE]
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
        self.write(self.cases[0], [data.NCBI_SINGLE, data.NCBI_MULTI, data.WASM_SINGLE])
        for tool in (data.NCBI_MULTI, data.WASM_SINGLE):
            self.record(data.result_paths(self.cases[0])[tool], ROW.replace("1e-10", "1.0e-10"))
        with patch.object(timing, "comparison_cases", return_value=self.cases[:1]), \
             patch.object(timing, "PLOT_DIR", self.root / "plots"), \
             patch.object(timing, "OUTPUT_IMAGE", self.root / "plots/times.png"), \
             redirect_stdout(io.StringIO()):
            timing.main()
        import pandas as pd
        rows = pd.read_csv(self.root / "plots/execution_times.tsv", sep="\t")
        self.assertEqual(rows["Tool"].tolist(), [data.NCBI_SINGLE])

    # NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
    # CATCH_ALL(status) ... return status;
    # Preserve search failures and reject successful commands without output.
    def test_runner_failure_and_missing_output_are_retained(self):
        binary = self.root / "fake-losat"
        env = dict(os.environ, BENCHMARK_PROGRAMS="blastp", BENCHMARK_CASE="WSSV.PajaWSV",
                   LOSAT_THREADS="4", LOSAT_BIN=str(binary), RUN_NATIVE="1", RUN_NCBI="0",
                   RUN_LOSAT_WASM="0", RUN_LOSAT_WASM_THREADED="0")
        for number, (command, expected) in enumerate([("exit 7", 7), ("exit 0", 125)]):
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
                   RUN_NATIVE="1", RUN_NCBI="0", RUN_LOSAT_WASM="0", RUN_LOSAT_WASM_THREADED="0")
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

    # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:75
    # const string kArgNumThreads("num_threads");
    # A synthetic Node process tests harness routing only; actual artifact
    # validation and BLAST output remain covered by the real WASI gates.
    def test_threaded_only_and_serial_compatibility_are_independent(self):
        node = self.root / "fake-node"
        node.write_text(f'''#!{sys.executable}
import json, sys
from pathlib import Path
if sys.argv[1] == '-p':
    print('{{}}')
elif sys.argv[1] == '-':
    sys.stdin.read()
    serial, threaded, serial_enabled, threaded_enabled, inspector = sys.argv[2:]
    for filename, enabled in [(serial, serial_enabled), (threaded, threaded_enabled)]:
        if enabled == '1':
            assert Path(filename).is_file(), filename
else:
    args = sys.argv[1:]
    assert Path(args[0]).name in ['run_losat_wasi.js', 'run_losat_wasi_threads.js']
    assert Path(args[1]).is_file()
    Path(args[args.index('-out') + 1]).write_text({ROW!r})
''')
        node.chmod(0o755)
        for serial_flag, threaded_flag, n in [(0, 1, 4), (0, 1, 1), (1, 0, 4), (1, 1, 4), (0, 0, 4), ("", "", 4)]:
            serial_enabled = (str(serial_flag) or "0") == "1"
            threaded_enabled = (str(threaded_flag) or "1") == "1"
            with self.subTest(serial=serial_enabled, threaded=threaded_enabled, threads=n):
                label = f"selected-{serial_flag}-{threaded_flag}-{n}"
                directory = self.root / label
                serial = self.root / f"{label}-serial.wasm"
                threaded = self.root / f"{label}-threaded.wasm"
                if serial_enabled: serial.write_bytes(b"synthetic serial fixture")
                if threaded_enabled: threaded.write_bytes(b"synthetic threaded fixture")
                env = dict(os.environ, BENCHMARK_DIR=str(directory),
                           BENCHMARK_PROGRAMS="blastp", BENCHMARK_CASE="WSSV.PajaWSV",
                           RUN_NATIVE="0", RUN_NCBI="0", RUN_LOSAT_WASM=str(serial_flag),
                           RUN_LOSAT_WASM_THREADED=str(threaded_flag), LOSAT_THREADS=str(n),
                           BUILD_LOSAT_WASM="0", BUILD_LOSAT_WASM_THREADED="0",
                           LOSAT_WASM_BIN=str(serial), LOSAT_WASM_THREADED_BIN=str(threaded),
                           NODE_BIN=str(node), NODE_ARGS_JSON="[]")
                result = subprocess.run(["bash", str(data.SCRIPT_DIR / "run_comparison.sh")],
                                        env=env, capture_output=True, text=True, timeout=30)
                if not (serial_enabled or threaded_enabled):
                    self.assertEqual(result.returncode, 2, result.stderr)
                    self.assertIn("No runners selected", result.stderr)
                    continue
                self.assertEqual(result.returncode, 0, result.stderr)
                expected = ({"WSSV.PajaWSV.losatp.wasm.out"} if serial_enabled else set())
                if threaded_enabled:
                    expected.add("WSSV.PajaWSV.losatp.wasm.n1.out")
                    if n != 1: expected.add(f"WSSV.PajaWSV.losatp.wasm.n{n}.out")
                self.assertEqual({p.name for p in (directory / "losat_out").glob("*.out")}, expected)
                manifest = json.loads((directory / "run.json").read_text())
                self.assertEqual(manifest["enabled_runners"]["RUN_LOSAT_WASM"], bool(serial_enabled))
                self.assertEqual(manifest["enabled_runners"]["RUN_LOSAT_WASM_THREADED"], bool(threaded_enabled))
                self.assertIsNotNone(manifest["node_versions"])
                for output in (directory / "losat_out").glob("*.out"):
                    self.assertTrue(data.successful_output(output))


# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46,48,51,75
# const string kArgQuery("query"); const string kArgOutput("out");
# const string kArgSubject("subject"); const string kArgNumThreads("num_threads");
# Synthetic commands verify directory handoff only, not BLAST output parity.
class ComparisonDirectoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.scripts = self.root / "crate/tests"
        self.scripts.mkdir(parents=True)
        self.plots = ["plot_overall_trend.py", "plot_comparison.py", "plot_execution_time.py"]
        for name in ["comparison_data.py", "run_comparison.sh", *self.plots]:
            shutil.copy2(data.SCRIPT_DIR / name, self.scripts / name)
        for name in ["Cargo.toml", "Cargo.lock", "build.rs", ".cargo/config.toml"]:
            path = self.scripts.parent / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("synthetic metadata fixture\n")
        (self.scripts / "fasta").mkdir()
        (self.scripts / "fasta/query.faa").write_text(">q\nACDEFGHIKLMNPQRSTVWY\n")
        (self.scripts / "fasta/subject.faa").write_text(">s\nACDEFGHIKLMNPQRSTVWY\n")
        (self.scripts / "comparison_cases.tsv").write_text(
            "task\tquery\tsubject\tname\tlosat_stem\tncbi_stem\tquery_gencode\tdb_gencode\n"
            "blastp\tquery.faa\tsubject.faa\texample\texample.losatp\texample.blastp\t1\t1\n"
            "blastp\tquery.faa\tsubject.faa\texcluded\texcluded.losatp\texcluded.blastp\t1\t1\n"
        )
        self.binary = self.root / "synthetic-search"
        self.binary.write_text(f'''#!{sys.executable}
import sys
from pathlib import Path
if sys.argv[1:] == ['-version']:
    print('synthetic oracle version')
else:
    Path(sys.argv[sys.argv.index('-out') + 1]).write_text({ROW + ROW.replace('30', '60')!r})
''')
        self.binary.chmod(0o755)
        self.oracle = self.root / "blastp"
        shutil.copy2(self.binary, self.oracle)
        # NCBI reference: c++/src/app/blastdb/makeblastdb.cpp:236-247
        # arg_desc->SetConstraint(kArgDbType, &(*new CArgAllow_Strings, "nucl", "prot"));
        # arg_desc->AddFlag("parse_seqids", ...);
        # Synthetic DB files test harness routing without an NCBI dependency.
        self.makeblastdb = self.root / "makeblastdb"
        self.makeblastdb.write_text(f'''#!{sys.executable}
import json, sys
from pathlib import Path
args = sys.argv[1:]
if args == ['-version']:
    print('synthetic database builder version')
else:
    assert '-parse_seqids' in args
    dbtype = args[args.index('-dbtype') + 1]
    assert dbtype in ('nucl', 'prot')
    prefix = args[args.index('-out') + 1]
    path = Path(prefix + '.synthetic')
    assert not path.exists(), 'database must be built only once'
    path.write_text(json.dumps(args))
''')
        self.makeblastdb.chmod(0o755)
        self.env = {key: value for key, value in os.environ.items()
                    if not key.startswith(("BENCHMARK_", "LOSAT_", "LOSATP_", "RUN_", "BUILD_", "NODE_"))}
        self.env["MPLCONFIGDIR"] = str(self.root / "matplotlib")
        self.runner_env = dict(self.env, BENCHMARK_PROGRAMS="blastp", BENCHMARK_CASE="example",
                               LOSAT_THREADS="4", LOSAT_BIN=str(self.binary), BLASTP_BIN=str(self.oracle),
                               MAKEBLASTDB_BIN=str(self.makeblastdb),
                               RUN_NATIVE="1", RUN_NCBI="1", RUN_LOSAT_WASM="0", RUN_LOSAT_WASM_THREADED="0")

    def run_comparison(self, **settings):
        return subprocess.run(["bash", str(self.scripts / "run_comparison.sh")],
                              env={**self.runner_env, **settings}, cwd=self.root,
                              capture_output=True, text=True, timeout=30)

    def latest(self):
        marker = self.scripts / "benchmark-runs/latest.json"
        return Path(json.loads(marker.read_text())["directory"])

    def plot(self, name, **settings):
        return subprocess.run([sys.executable, str(self.scripts / name)],
                              env={**self.env, **settings}, cwd=self.root,
                              capture_output=True, text=True, timeout=45)

    def test_default_directory_and_settings_reach_all_three_plots(self):
        result = self.run_comparison()
        self.assertEqual(result.returncode, 0, result.stderr)
        directory = self.latest()
        self.assertEqual(directory.parent, self.scripts / "benchmark-runs")
        self.assertEqual(json.loads((directory / "run.json").read_text())["status"], "COMPLETE")
        for name in self.plots:
            result = self.plot(name)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(str(directory), result.stdout)
            self.assertNotIn("excluded", result.stdout)
        for name in ["overall_trend_comparison.png", "compare_example_BLASTP.png", "execution_time_comparison_all.png"]:
            self.assertGreater((directory / "plots" / name).stat().st_size, 0)
        table = (directory / "plots/execution_times.tsv").read_text()
        # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:75
        # const string kArgNumThreads("num_threads");
        self.assertIn("BLAST+ n1", table)
        self.assertIn("BLAST+ n4", table)
        self.assertIn("LOSAT native n4", table)
        self.assertNotIn("LOSAT native n8", table)

    # NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3223-3239,1052-1054
    # if (args.Exist(kArgSubject) && args[kArgSubject].HasValue() &&
    #     m_NumThreads != CThreadable::kMinNumThreads) { m_NumThreads = ...; }
    # opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
    # Exercise the shell boundary for every task, preserving translation options.
    def test_subject_targets_and_tblastx_database_preserve_both_thread_counts(self):
        (self.scripts / "fasta/query.fna").write_text(">q\nACGTACGT\n")
        (self.scripts / "fasta/subject.fna").write_text(">s\nACGTACGT\n")
        header = "task\tquery\tsubject\tname\tlosat_stem\tncbi_stem\tquery_gencode\tdb_gencode\n"
        tasks = ("tblastx", "blastn", "megablast", "blastp")
        rows = []
        for task in tasks:
            ext = "faa" if task == "blastp" else "fna"
            rows.append(f"{task}\tquery.{ext}\tsubject.{ext}\t{task}\t{task}\t{task}\t4\t4\n")
        rows.append("tblastx\tquery.fna\tsubject.fna\ttblastx_repeat\ttblastx_repeat\ttblastx_repeat\t4\t4\n")
        (self.scripts / "comparison_cases.tsv").write_text(header + "".join(rows))
        for n in (4, 1):
            with self.subTest(threads=n):
                directory = self.root / f"db-threads-{n}"
                result = self.run_comparison(
                    BENCHMARK_DIR=str(directory), BENCHMARK_PROGRAMS=",".join(tasks),
                    BENCHMARK_CASE="", LOSAT_THREADS=str(n), RUN_NATIVE="0",
                    BLASTN_BIN=str(self.oracle), TBLASTX_BIN=str(self.oracle))
                self.assertEqual(result.returncode, 0, result.stderr)
                records = list((directory / "blast_out").glob("*.run.json"))
                self.assertEqual(len(records), 5 if n == 1 else 10)
                for task in tasks:
                    for threads in sorted({1, n}):
                        record = json.loads((directory / f"blast_out/{task}.n{threads}.run.json").read_text())
                        args = record["ordered_argv"]
                        self.assertEqual(args[args.index("-num_threads") + 1], str(threads))
                        if task in {"blastn", "megablast"}:
                            self.assertEqual(args[args.index("-task") + 1], task)
                        if task == "tblastx":
                            self.assertNotIn("-subject", args)
                            db = Path(args[args.index("-db") + 1])
                            self.assertEqual(db.parent.name, "nucl")
                            self.assertIn(str(db) + ".synthetic", record["files"])
                            for option in ("-query_gencode", "-db_gencode"):
                                self.assertEqual(args[args.index(option) + 1], "4")
                        else:
                            self.assertNotIn("-db", args)
                            ext = "faa" if task == "blastp" else "fna"
                            self.assertEqual(args[args.index("-subject") + 1], str(self.scripts / f"fasta/subject.{ext}"))
                builds = list((directory / "blast_out/db").rglob("*.makeblastdb.json"))
                self.assertEqual(len(builds), 1)
                for build in builds:
                    record = json.loads(build.read_text())
                    args = record["ordered_argv"]
                    self.assertEqual(args[args.index("-dbtype") + 1], build.parent.name)
                    self.assertGreaterEqual(record["wall_seconds"], 0)
        # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:50-51
        # const string kArgDb("db"); const string kArgSubject("subject");
        # A subject-only selection neither requires nor invokes makeblastdb.
        directory = self.root / "subject-only"
        result = self.run_comparison(
            BENCHMARK_DIR=str(directory), BENCHMARK_PROGRAMS="blastn,megablast,blastp",
            BENCHMARK_CASE="", RUN_NATIVE="0", BLASTN_BIN=str(self.oracle),
            MAKEBLASTDB_BIN=str(self.root / "missing-makeblastdb"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((directory / "blast_out/db").exists())
        self.assertEqual(len(list((directory / "blast_out").glob("*.run.json"))), 6)

    # NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
    # CATCH_ALL(status) ... return status;
    # A required DB build must fail before any oracle search, without -subject fallback.
    def test_database_failure_stops_before_oracle_search(self):
        # NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1052-1054
        # if (m_Target == eDatabase && args[kArgDbGeneticCode] &&
        #     (program == eTblastn || program == eTblastx)) { ... }
        cases = self.scripts / "comparison_cases.tsv"
        cases.write_text(cases.read_text().replace("blastp\t", "tblastx\t"))
        self.makeblastdb.write_text("#!/bin/bash\nexit 9\n")
        result = self.run_comparison(RUN_NATIVE="0", BENCHMARK_PROGRAMS="tblastx", TBLASTX_BIN=str(self.oracle))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Database preparation failed", result.stderr)
        self.assertFalse(list((self.latest() / "blast_out").glob("*.run.json")))
        self.assertNotEqual(json.loads((self.latest() / "run.json").read_text())["status"], "COMPLETE")

    # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:50-51,75
    # const string kArgDb("db"); const string kArgSubject("subject");
    # const string kArgNumThreads("num_threads");
    def test_previous_protocols_are_not_relabeled_as_current_targets(self):
        result = self.run_comparison()
        self.assertEqual(result.returncode, 0, result.stderr)
        path = self.latest() / "run.json"
        manifest = json.loads(path.read_text())
        for schema in ("losat-simple-comparison-v2", "losat-simple-comparison-v3"):
            with self.subTest(schema=schema):
                manifest["schema"] = schema
                path.write_text(json.dumps(manifest))
                result = self.plot("plot_execution_time.py")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("previous target/thread protocol", result.stderr)

    # NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
    # CATCH_ALL(status) ... return status;
    # Both a failed search and failed metadata collection must block stale plots.
    def test_latest_failure_does_not_plot_an_older_success(self):
        first = self.run_comparison()
        self.assertEqual(first.returncode, 0, first.stderr)
        previous = self.latest()
        saved = (previous / "losat_out/example.losatp.out").read_bytes()
        self.binary.write_text("#!/bin/bash\nexit 7\n")
        failed = self.run_comparison()
        self.assertEqual(failed.returncode, 7, failed.stderr)
        self.assertNotEqual(self.latest(), previous)
        for name in self.plots:
            result = self.plot(name)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("did not complete", result.stderr)
        failed = self.run_comparison(RUN_LOSAT_WASM_THREADED="1", NODE_BIN=str(self.root / "missing-node"))
        self.assertNotEqual(failed.returncode, 0)
        self.assertFalse((self.latest() / "run.json").exists())
        result = self.plot("plot_execution_time.py")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No comparison run metadata", result.stderr)
        self.assertEqual((previous / "losat_out/example.losatp.out").read_bytes(), saved)
        self.assertFalse((previous / "plots").exists())

    def test_explicit_relative_directory_does_not_replace_default(self):
        result = self.run_comparison()
        self.assertEqual(result.returncode, 0, result.stderr)
        previous = self.latest()
        result = self.run_comparison(BENCHMARK_DIR="custom/run", LOSAT_THREADS="2")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.latest(), previous)
        result = self.plot("plot_execution_time.py", BENCHMARK_DIR="custom/run")
        self.assertEqual(result.returncode, 0, result.stderr)
        table = (self.scripts / "custom/run/plots/execution_times.tsv").read_text()
        self.assertIn("LOSAT native n2", table)
        result = self.plot("plot_execution_time.py", BENCHMARK_DIR="custom/missing")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No comparison run metadata", result.stderr)


if __name__ == "__main__":
    unittest.main()
