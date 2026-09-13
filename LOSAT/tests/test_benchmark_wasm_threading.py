"""Subprocess regressions for failed and incomplete benchmark evidence."""
import json
import io
from contextlib import redirect_stderr
from unittest.mock import patch
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import benchmark_wasm_threading as benchmark


# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1100-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Synthetic executables test orchestration, not the search or its oracle bytes.
class BenchmarkFailureTests(unittest.TestCase):
    def check_failure(self, phase):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            oracle = root / "oracle"
            oracle.mkdir()
            script = f'''#!{sys.executable}
import os,pathlib,sys
if '-version' in sys.argv:
    print('synthetic oracle 2.17.0');sys.exit(0)
pathlib.Path(sys.argv[sys.argv.index('-out')+1]).write_bytes(b'fixture\\n')
n=sys.argv[sys.argv.index('-num_threads')+1]
if os.environ.get('LOSAT_WASI_THREADS_DEBUG'):
    print('[losat-thread-pool] program=blastp requested_threads='+n+' pool_threads=0 caller_participates=false',file=sys.stderr)
'''
            for name in ["blastn", "blastp", "tblastx"]:
                path = oracle / name
                path.write_text(script)
                path.chmod(0o755)
            for version in ["baseline", "candidate"]:
                path = root / version / "native-command/release/LOSAT"
                path.parent.mkdir(parents=True)
                code = script
                if version == "candidate":
                    code += f"sys.exit(7 if bool(os.environ.get('LOSAT_WASI_THREADS_DEBUG')) == {phase == 'diagnostic'} else 0)\n"
                path.write_text(code)
                path.chmod(0o755)
            output = root / "results"
            command = [sys.executable, str(benchmark.TESTS / "benchmark_wasm_threading.py"),
                "--candidate-dir", str(root / "candidate"), "--baseline-dir", str(root / "baseline"),
                "--artifacts", str(root), "--baseline-runners", str(benchmark.TESTS),
                "--oracle-dir", str(oracle), "--output-dir", str(output),
                "--case", "WSSV.PajaWSV.losatp", "--kinds", "native", "--threads", "1",
                "--warmups", "0", "--repeats", "1", "--skip-reuse"]
            child = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(child.returncode, 0, child.stdout)
            samples = json.loads((output / "samples.json").read_text())
            failed = [x for x in samples if x.get("version") == "candidate" and x.get("phase") == phase]
            self.assertEqual(len(failed), 1, child.stderr)
            self.assertEqual(failed[0]["status"], "BAD_EXIT")
            self.assertEqual(failed[0]["exit_status"], 7)
            self.assertEqual(json.loads((output / "summary.json").read_text()), [])
            self.assertIn(json.loads((output / "run-status.json").read_text())["status"], {"FAILED", "PARTIAL"})
            if phase == "diagnostic":
                self.assertFalse(any(x["phase"] == "cold" for x in samples))
                self.assertTrue(json.loads((output / "excluded.json").read_text()))

    # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:178-188
    # (*thread)->Join(&result);
    # Non-finite limits must fail before starting a potentially unbounded child.
    def test_nonfinite_timeouts_rejected_before_run_creation(self):
        required = [value for name in ['candidate-dir', 'artifacts', 'baseline-dir',
                    'baseline-runners', 'oracle-dir', 'output-dir']
                    for value in ['--' + name, '/unused']]
        for option in ['--timeout', '--reuse-timeout']:
            for value in ['nan', 'inf', '-inf', '0']:
                with self.subTest(option=option, value=value), patch.object(sys, 'argv',
                        ['benchmark', *required, option + '=' + value]), \
                        patch.object(benchmark, 'run_benchmark') as run, redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as error:
                        benchmark.main()
                    self.assertEqual(error.exception.code, 2)
                    run.assert_not_called()

    def test_diagnostic_failure_retained_and_timings_skipped(self):
        self.check_failure("diagnostic")

    def test_cold_failure_retained_before_abort(self):
        self.check_failure("cold")


if __name__ == "__main__":
    unittest.main()
