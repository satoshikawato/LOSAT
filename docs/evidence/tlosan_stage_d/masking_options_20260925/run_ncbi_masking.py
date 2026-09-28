#!/usr/bin/env python3
"""Save pinned NCBI Stage D local masking and composition option traces."""
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
STAGE_C = HERE.parent.parent / "tlosan_stage_c"
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import (  # noqa: E402
    NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha,
)
CASES = (
    ("seg_cross_traceback_20260924", True, False, False),
    ("seg_soft_traceback_20260924", True, True, False),
    ("lcase_query_20260924", False, False, True),
    ("lcase_soft_query_20260924", False, True, True),
    ("seg_lcase_overlap_20260924", True, False, True),
    ("seg_lcase_overlap_soft_20260924", True, True, True),
)
D_PROBE = STAGE_D / "ncbi_d_call_trace.c"
K_PROBE = STAGE_D / "ncbi_kappa_traceback_trace.c"


def run_probe(command, plain, source, prefix, output):
    with tempfile.TemporaryDirectory(prefix="tlosan-d-mask-") as tmp:
        traced = probe_run(command, source, Path(tmp))
    selected = b"".join(line for line in traced.stderr.splitlines(keepends=True)
                        if line.startswith(prefix))
    ordinary = b"".join(line for line in traced.stderr.splitlines(keepends=True)
                        if not line.startswith(prefix))
    assert traced.stdout == plain.stdout and ordinary == plain.stderr
    assert selected
    output.write_bytes(selected)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_masking.py NEW_OUTPUT_DIRECTORY")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    manifest = [
        f"NCBI source commit: {SOURCE_COMMIT}",
        f"NCBI executable SHA256: {NCBI_SHA256}",
        f"D probe SHA256: {sha(D_PROBE)}",
        f"Kappa probe SHA256: {sha(K_PROBE)}",
    ]
    for case, seg, soft, lowercase in CASES:
        query = STAGE_C / case / "query.faa"
        subject = STAGE_C / case / "subjects.fna"
        manifest.extend([f"{case} query SHA256: {sha(query)}",
                         f"{case} subject SHA256: {sha(subject)}"])
        for mode in (0, 2):
            stem = f"{case}.mode{mode}"
            command = [str(NCBI), "-task", "tblastn", "-query", str(query),
                       "-subject", str(subject), "-db_gencode", "1",
                       "-num_threads", "1", "-evalue", "10000",
                       "-max_target_seqs", "500", "-comp_based_stats", str(mode),
                       "-sum_stats", "false" if mode == 0 else "true",
                       "-seg", "12 2.2 2.5" if seg else "no",
                       "-outfmt", "6 qseqid sseqid score bitscore evalue nident positive "
                       "length mismatch gaps gapopen qstart qend sstart send sframe"]
            if soft: command.extend(["-soft_masking", "true"])
            if lowercase: command.append("-lcase_masking")
            plain = subprocess.run(command, capture_output=True, check=True)
            (output / f"{stem}.report.tsv").write_bytes(plain.stdout)
            (output / f"{stem}.stderr").write_bytes(plain.stderr)
            run_probe(command, plain, D_PROBE, b"D_", output / f"{stem}.calls.trace")
            if mode == 2:
                run_probe(command, plain, K_PROBE, b"K_TRACE_", output / f"{stem}.kappa.trace")
            manifest.append(f"{stem} command: {command!r}")
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files))


if __name__ == "__main__":
    main()
