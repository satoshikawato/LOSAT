"""gbdraw through pylosat and through the executable of the same commit: the same files.

gbdraw's own entry points run each scenario twice. With pylosat, starting any process fails.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import gbdraw.api  # noqa: F401
from gbdraw.circular import circular_main
from gbdraw.linear import linear_main

LOSAT_BIN = os.environ["LOSAT_BIN"]
OUT = Path(sys.argv[1])
ROOT = Path.cwd()
LAMBDA = ROOT / "tests/test_inputs/NC_001416.gb"
TUTORIAL = ROOT / "gbdraw/web/tutorial-data"
DE3 = TUTORIAL / "de3/NC_042057.1.gb"
HUMAN = TUTORIAL / "human-mitochondrion/HmmtDNA.gbk"
COMPARISON = TUTORIAL / "metazoan-mitochondria-comparison"
FASTAS = [COMPARISON / f"{accession}.fna" for accession in ("NC_002333.2", "NC_024511.2", "NC_001328.1")]
LINEAR = ["--gbk", str(LAMBDA), str(DE3), "--bitscore", "50", "--evalue", "0.01", "--identity", "0",
          "--alignment_length", "0", "--comparison_height", "120", "-f", "svg", "--save_session"]
CIRCULAR = ["--gbk", str(HUMAN), "--conservation_sequence", *map(str, FASTAS),
            "--bitscore", "50", "--evalue", "1e-5", "--identity", "40", "--alignment_length", "50",
            "-f", "svg", "--save_session"]

SCENARIOS = {
    "linear-losatn": (linear_main, [*LINEAR, "--losat", "losatn"]),
    "linear-losatn-blastn-task": (linear_main, [*LINEAR, "--losat", "losatn", "--losatn_task", "blastn"]),
    "linear-tlosatx": (linear_main, [*LINEAR, "--losat", "tlosatx"]),
    "linear-tlosatx-code-11": (linear_main, [*LINEAR, "--losat", "tlosatx", "--losat_gencode", "11", "11"]),
    "linear-losatp": (linear_main, [*LINEAR, "--losat", "losatp"]),
    "circular-tlosatx": (circular_main, [*CIRCULAR, "--losat", "tlosatx", "--losat_gencode", "2",
                                         "--conservation_losat_gencode", "2", "5", "5"]),
    "circular-losatn": (circular_main, [*CIRCULAR, "--losat", "losatn", "--losatn_task", "blastn"]),
}


def refuse(self, args, *rest, **options):
    raise AssertionError(f"a process was started: {args!r}")


real_init = subprocess.Popen.__init__


def run(name, main, argv, runtime):
    directory = OUT / name / runtime
    directory.mkdir(parents=True, exist_ok=True)
    extra = ["--losat_bin", LOSAT_BIN] if runtime == "executable" else []
    subprocess.Popen.__init__ = real_init if runtime == "executable" else refuse
    started = time.monotonic()
    try:
        main([*argv, *extra, "-o", str(directory / "figure")])
    finally:
        subprocess.Popen.__init__ = real_init
    return time.monotonic() - started, directory


def entries(directory):
    session = json.loads((directory / "figure.gbdraw-session.json").read_text())
    rows = session["losatCache"]["entries"]
    return rows, [{key: value for key, value in row.items() if key != "runtime"} for row in rows]


failures = 0
for name, (main, argv) in SCENARIOS.items():
    in_process_time, in_process = run(name, main, argv, "pylosat")
    executable_time, executable = run(name, main, argv, "executable")
    svg = (in_process / "figure.svg").read_bytes() == (executable / "figure.svg").read_bytes()
    rows, comparable = entries(in_process)
    other_rows, other_comparable = entries(executable)
    same = svg and comparable == other_comparable and len(rows) > 0
    failures += not same
    print(f"{name}: svg {svg}, {len(rows)} raw searches {comparable == other_comparable}, "
          f"{in_process_time:.2f} s in process, {executable_time:.2f} s with the executable, "
          f"runtime {rows[0]['runtime']['source']}/{rows[0]['runtime']['version']}", flush=True)
print("failures:", failures)
