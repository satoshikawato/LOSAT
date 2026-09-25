#!/usr/bin/env python3
"""Save fixed NCBI TBLASTN local cross-option Stage D comparisons."""
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
STAGE_C = HERE.parent.parent / "tlosan_stage_c"
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha

CASES = (("mode0_sumtrue", 0, "true"), ("mode2_sumfalse", 2, "false"))
D_PROBE = STAGE_D / "ncbi_d_call_trace.c"
K_PROBE = STAGE_D / "ncbi_kappa_traceback_trace.c"


def traced(command, plain, source, prefix, path):
    with tempfile.TemporaryDirectory(prefix="tlosan-d-option-") as tmp:
        run = probe_run(command, source, Path(tmp))
    selected = b"".join(line for line in run.stderr.splitlines(keepends=True) if line.startswith(prefix))
    ordinary = b"".join(line for line in run.stderr.splitlines(keepends=True) if not line.startswith(prefix))
    assert run.stdout == plain.stdout and ordinary == plain.stderr and selected
    path.write_bytes(selected)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_options.py NEW_OUTPUT_DIRECTORY")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    query = STAGE_C / "multi_query_20260924/query.faa"
    subject = STAGE_C / "multi_query_20260924/subjects.fna"
    manifest = [f"NCBI source commit: {SOURCE_COMMIT}", f"NCBI executable SHA256: {NCBI_SHA256}",
                f"D probe SHA256: {sha(D_PROBE)}", f"Kappa probe SHA256: {sha(K_PROBE)}",
                f"query SHA256: {sha(query)}", f"subject SHA256: {sha(subject)}"]
    for name, composition, sum_stats in CASES:
        command = [str(NCBI), "-task", "tblastn", "-query", str(query), "-subject", str(subject),
                   "-db_gencode", "1", "-num_threads", "1", "-comp_based_stats", str(composition),
                   "-sum_stats", sum_stats, "-outfmt",
                   "6 qseqid sseqid score bitscore evalue nident positive length mismatch gaps gapopen qstart qend sstart send sframe"]
        plain = subprocess.run(command, capture_output=True, check=True)
        (output / f"{name}.report.tsv").write_bytes(plain.stdout)
        (output / f"{name}.stderr").write_bytes(plain.stderr)
        traced(command, plain, D_PROBE, b"D_", output / f"{name}.calls.trace")
        if composition:
            traced(command, plain, K_PROBE, b"K_TRACE_", output / f"{name}.kappa.trace")
        manifest.append(f"{name} command: {command!r}")
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text("".join(f"{sha(path)}  {path.name}\n" for path in files))

if __name__ == "__main__":
    main()
