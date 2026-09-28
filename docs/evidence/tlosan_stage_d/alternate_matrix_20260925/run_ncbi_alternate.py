#!/usr/bin/env python3
"""Save pinned NCBI BLOSUM45/word-2 TBLASTN Stage D local trace."""
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
STAGE_C = HERE.parent.parent / "tlosan_stage_c" / "alternate_matrix_word2_20260924"
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha


def main():
    if len(sys.argv) != 2: raise SystemExit("usage: run_ncbi_alternate.py NEW_OUTPUT_DIRECTORY")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = STAGE_D / "ncbi_d_call_trace.c"
    assert subprocess.check_output(["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    query, subject = STAGE_C / "query.faa", STAGE_C / "subjects.fna"
    command = [str(NCBI), "-task", "tblastn", "-query", str(query), "-subject", str(subject),
               "-db_gencode", "1", "-matrix", "BLOSUM45", "-word_size", "2", "-threshold", "16",
               "-window_size", "60", "-gapopen", "14", "-gapextend", "2", "-evalue", "10000",
               "-num_threads", "1", "-comp_based_stats", "0", "-seg", "no", "-sum_stats", "false",
               "-outfmt", "6 qseqid sseqid score bitscore evalue nident positive length mismatch gaps gapopen qstart qend sstart send sframe"]
    plain = subprocess.run(command, capture_output=True, check=True)
    (output / "report.tsv").write_bytes(plain.stdout)
    (output / "stderr").write_bytes(plain.stderr)
    with tempfile.TemporaryDirectory(prefix="tlosan-d-alternate-") as tmp:
        run = probe_run(command, source, Path(tmp))
    selected = b"".join(line for line in run.stderr.splitlines(keepends=True) if line.startswith(b"D_"))
    ordinary = b"".join(line for line in run.stderr.splitlines(keepends=True) if not line.startswith(b"D_"))
    assert run.stdout == plain.stdout and ordinary == plain.stderr and selected
    (output / "calls.trace").write_bytes(selected)
    manifest = [f"NCBI source commit: {SOURCE_COMMIT}", f"NCBI executable SHA256: {NCBI_SHA256}",
                f"D probe SHA256: {sha(source)}", f"query SHA256: {sha(query)}",
                f"subject SHA256: {sha(subject)}", f"command: {command!r}"]
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    (output / "outputs.sha256").write_text("".join(f"{sha(path)}  {path.name}\n" for path in sorted(output.iterdir()) if path.is_file()))

if __name__ == "__main__": main()
