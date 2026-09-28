#!/usr/bin/env python3
"""Save pinned comparison-only NCBI initial/local parameter call states."""
from __future__ import annotations

# Pinned NCBI blast_setup.c:964-985,1001-1024 and blast_engine.c:1390-1445:
# initial hit parameters may be the inputs to search when TotLen is nonzero;
# the conditional OneSubjectUpdateParameters call is a separate call state.
import hashlib
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

STAGE_C = Path(__file__).resolve().parents[1] / "tlosan_stage_c"
sys.path.insert(0, str(STAGE_C))
from run_multi_query_trace import (  # noqa: E402
    NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha,
)

HERE = Path(__file__).resolve().parent
CASES = (
    "seg_hard_query_20260924_default",
    "seg_hard_query_20260924_control",
    "multi_query_20260924_default",
    "multi_query_20260924_control",
    "run_20260923_default",
    "run_20260923_control",
)
PROBES = (
    (HERE / "ncbi_parameter_trace.c", "D_PARAM_", "parameters"),
    (HERE / "ncbi_spouge_trace.c", "D_SPOUGE\t", "spouge"),
    (STAGE_C / "ncbi_context_cutoff_trace.c", "GAPPED_CONTEXT_CUTOFF\t", "hit_cutoffs"),
    (STAGE_C / "ncbi_wordfinder_context_cutoff_trace.c", "WORD_CONTEXT_CUTOFF\t", "word_cutoffs"),
)


def command_from_manifest(case: str) -> list[str]:
    manifest = (HERE / "run_20260924" / f"{case}.manifest.txt").read_text()
    line = manifest.split("Command: ", 1)[1].splitlines()[0]
    head, outfmt = line.split(" -outfmt ", 1)
    if " -comp_based_stats " in outfmt:
        outfmt, tail = outfmt.split(" -comp_based_stats ", 1)
        return shlex.split(head) + ["-outfmt", outfmt, "-comp_based_stats"] + shlex.split(tail)
    return shlex.split(head) + ["-outfmt", outfmt]


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_parameter_trace.py NEW_OUTPUT_DIR")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    manifest = [
        f"NCBI source commit: {SOURCE_COMMIT}",
        f"NCBI binary SHA256: {NCBI_SHA256}",
    ]
    for source, _, label in PROBES:
        manifest.append(f"{label} probe SHA256: {sha(source)}")
    for case in CASES:
        command = command_from_manifest(case)
        plain = subprocess.run(command, capture_output=True, check=True)
        assert plain.stdout == (HERE / "run_20260924" / f"{case}.out").read_bytes()
        manifest.append(f"{case} command: {command!r}")
        for source, prefix, label in PROBES:
            with tempfile.TemporaryDirectory(prefix="tlosan-d-parameter-") as tmp:
                traced = probe_run(command, source, Path(tmp))
            assert traced.stdout == plain.stdout
            lines = traced.stderr.splitlines(keepends=True)
            selected = [line for line in lines if line.startswith(prefix.encode())]
            if label == "parameters":
                assert b"".join(line for line in lines if not line.startswith(b"D_PARAM_")) == plain.stderr
                assert any(line.startswith(b"D_PARAM_HIT\t") for line in selected)
            assert selected
            (output / f"{case}.{label}.tsv").write_bytes(b"".join(selected))
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files)
    )


if __name__ == "__main__":
    main()
