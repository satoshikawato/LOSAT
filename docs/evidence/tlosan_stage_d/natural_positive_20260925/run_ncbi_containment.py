#!/usr/bin/env python3
"""Save pinned comparison-only NCBI natural postredo containment evidence."""
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
    (HERE / "ncbi_containment_probe.c", b"CP_", "containment"),
    (STAGE_D / "ncbi_kappa_traceback_trace.c", b"K_TRACE_", "kappa"),
    (STAGE_D / "ncbi_d_call_trace.c", b"D_", "calls"),
)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_containment.py NEW_OUTPUT_DIRECTORY")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    command = [
        str(NCBI), "-task", "tblastn", "-query", str(HERE / "containment_query.faa"),
        "-subject", str(HERE / "containment_subject.fna"), "-db_gencode", "1",
        "-num_threads", "1", "-seg", "no", "-evalue", "10000",
        "-max_target_seqs", "1",
        "-outfmt", "6 qseqid sseqid score bitscore evalue nident positive "
        "length mismatch gaps gapopen qstart qend sstart send sframe",
    ]
    plain = subprocess.run(command, capture_output=True, check=True)
    (output / "report.tsv").write_bytes(plain.stdout)
    (output / "stderr").write_bytes(plain.stderr)
    manifest = [
        f"NCBI source commit: {SOURCE_COMMIT}",
        f"NCBI executable SHA256: {NCBI_SHA256}",
        f"containment_query.faa SHA256: {sha(HERE / 'containment_query.faa')}",
        f"containment_subject.fna SHA256: {sha(HERE / 'containment_subject.fna')}",
        f"Command: {command!r}",
    ]
    for source, prefix, label in PROBES:
        with tempfile.TemporaryDirectory(prefix="tlosan-d-containment-") as tmp:
            traced = probe_run(command, source, Path(tmp))
        selected = b"".join(line for line in traced.stderr.splitlines(keepends=True)
                            if line.startswith(prefix))
        ordinary = b"".join(line for line in traced.stderr.splitlines(keepends=True)
                            if not line.startswith(prefix))
        assert traced.stdout == plain.stdout
        assert ordinary == plain.stderr
        assert selected
        (output / f"{label}.trace").write_bytes(selected)
        manifest.append(f"{label} probe SHA256: {sha(source)}")
    containment = (output / "containment.trace").read_text().splitlines()
    assert [row for row in containment if row.startswith("CP_A\t")] == [
        "CP_A\t0\t2\t2\t0"
    ]
    assert [row for row in containment if row.startswith("CP_L\t")] == [
        "CP_L\t0\t1"
    ]
    kappa = (output / "kappa.trace").read_text().splitlines()
    assert sum(row.startswith("K_TRACE_RETURN\t") for row in kappa) == 2
    assert sum(row.startswith("K_TRACE_HEAP_HSP\t") for row in kappa) == 1
    calls = (output / "calls.trace").read_text().splitlines()
    assert any(row.startswith("D_LIST\t") and "\tlink_before\t0\t0\t1\t" in row
               for row in calls)
    assert len(plain.stdout.splitlines()) == 1
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files)
    )


if __name__ == "__main__":
    main()
