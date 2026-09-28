#!/usr/bin/env python3
"""Save pinned NCBI local-subject report fields as Stage D numeric evidence."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, sha


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_report_payload_trace.py NEW_OUTPUT_DIR")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    command = [
        str(NCBI), "-task", "tblastn",
        "-query", str(HERE / "query.faa"),
        "-subject", str(HERE / "subjects.fna"),
        "-db_gencode", "1", "-num_threads", "1",
        "-evalue", "10", "-max_target_seqs", "2",
        "-outfmt", "6 qseqid sseqid score nident positive length mismatch gaps gapopen qstart qend sstart send sframe",
    ]
    result = subprocess.run(command, capture_output=True, check=True)
    rows = result.stdout.decode().splitlines()
    assert len(rows) == 2
    assert [row.split("\t")[1] for row in rows] == ["weak_17", "weak_18"]
    (output / "report_fields.tsv").write_bytes(result.stdout)
    (output / "stderr.txt").write_bytes(result.stderr)
    broad_command = command.copy()
    broad_command[broad_command.index("-max_target_seqs") + 1] = "112"
    broad = subprocess.run(broad_command, capture_output=True, check=True)
    broad_rows = broad.stdout.decode().splitlines()
    assert len(broad_rows) == 19
    assert any(row.split("\t")[1] == "weak_16" and row.split("\t")[7] == "5"
               for row in broad_rows)
    (output / "report_fields_all.tsv").write_bytes(broad.stdout)
    (output / "stderr_all.txt").write_bytes(broad.stderr)
    (output / "manifest.txt").write_text(
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        f"Query SHA256: {sha(HERE / 'query.faa')}\n"
        f"Subject SHA256: {sha(HERE / 'subjects.fna')}\n"
        f"Command: {command!r}\n"
        f"Broad command: {broad_command!r}\n"
        "Comparison only; local -subject; code 1; one thread; Stage D numeric fields.\n"
    )
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files)
    )


if __name__ == "__main__":
    main()
