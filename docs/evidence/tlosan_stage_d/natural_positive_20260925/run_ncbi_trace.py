#!/usr/bin/env python3
"""Save comparison-only pinned NCBI Stage D natural positive traces."""
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import (  # noqa: E402
    NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha,
)

PROBES = (
    (STAGE_D / "ncbi_parameter_trace.c", b"D_PARAM_", "params"),
    (STAGE_D / "ncbi_d_call_trace.c", b"D_", "calls"),
    (STAGE_D / "ncbi_kappa_traceback_trace.c", b"K_TRACE_", "kappa"),
)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_trace.py NEW_OUTPUT_DIRECTORY")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    manifest = [
        f"NCBI source commit: {SOURCE_COMMIT}",
        f"NCBI executable SHA256: {NCBI_SHA256}",
        f"query.faa SHA256: {sha(HERE / 'query.faa')}",
    ]
    for case, subject, evalue in (
        ("evalue_sort", "evalue_subjects.fna", "10000"),
        ("heap_replacement", "heap_subjects.fna", "10000"),
    ):
        manifest.append(f"{subject} SHA256: {sha(HERE / subject)}")
        command = [
            str(NCBI), "-task", "tblastn", "-query", str(HERE / "query.faa"),
            "-subject", str(HERE / subject), "-db_gencode", "1",
            "-num_threads", "1", "-seg", "no", "-evalue", evalue,
            "-max_target_seqs", "1",
            "-outfmt", "6 qseqid sseqid score bitscore evalue nident positive "
            "length mismatch gaps gapopen qstart qend sstart send sframe",
        ]
        manifest.append(f"{case} command: {command!r}")
        plain = subprocess.run(command, capture_output=True, check=True)
        (output / f"{case}.report.tsv").write_bytes(plain.stdout)
        (output / f"{case}.stderr").write_bytes(plain.stderr)
        for probe, prefix, label in PROBES:
            with tempfile.TemporaryDirectory(prefix="tlosan-d-positive-") as tmp:
                traced = probe_run(command, probe, Path(tmp))
            assert traced.stdout == plain.stdout
            selected = b"".join(line for line in traced.stderr.splitlines(keepends=True)
                                if line.startswith(prefix))
            ordinary = b"".join(line for line in traced.stderr.splitlines(keepends=True)
                                if not line.startswith(prefix))
            assert ordinary == plain.stderr, (case, label)
            assert selected, (case, label)
            (output / f"{case}.{label}.trace").write_bytes(selected)
            manifest.append(f"{label} probe SHA256: {sha(probe)}")
        kappa = (output / f"{case}.kappa.trace").read_text()
        if case == "evalue_sort":
            heap = [line.split("\t") for line in kappa.splitlines()
                    if line.startswith("K_TRACE_HEAP_HSP\t0\t")]
            retained = [line.split("\t") for line in kappa.splitlines()
                        if line.startswith("K_TRACE_RESULT_HSP\t0\t0\t")]
            assert [(int(row[3]), float(row[5])) for row in heap] == [
                (238, float(heap[0][5])), (189, float(heap[1][5])),
                (135, float(heap[2][5])),
            ]
            assert [int(row[4]) for row in retained] == [238, 135, 189]
            assert [float(row[5]) for row in retained] == sorted(
                float(row[5]) for row in retained)
        else:
            replacements = [line for line in kappa.splitlines()
                            if line.startswith("K_TRACE_HEAP_INSERT_RETURN\t")
                            and line.split("\t")[-1] != "-1"]
            assert replacements == [
                next(line for line in replacements if line.endswith("\t1"))
            ]
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files)
    )


if __name__ == "__main__":
    main()
