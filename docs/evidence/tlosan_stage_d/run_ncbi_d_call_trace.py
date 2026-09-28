#!/usr/bin/env python3
"""Replay retained Stage C fixtures through comparison-only NCBI Stage D calls."""
from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
NCBI_SHA256 = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"
SOURCE_COMMIT = "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4"
FIELDS = "6 qseqid sseqid score bitscore evalue qstart qend sstart send sframe length"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(profile: str, fixture: Path, output: Path, probe: Path) -> None:
    command = [str(NCBI), "-task", "tblastn", "-query",
               str(fixture / "query.faa"), "-subject",
               str(fixture / "subjects.fna"), "-db_gencode", "1",
               "-num_threads", "1", "-outfmt", FIELDS]
    if profile == "control":
        command.extend(["-comp_based_stats", "0", "-sum_stats", "false"])
    elif profile != "default":
        raise ValueError(profile)
    plain = subprocess.run(command, capture_output=True, check=True)
    env = os.environ.copy()
    env["LD_PRELOAD"] = str(probe)
    traced = subprocess.run(command, env=env, capture_output=True, check=True)
    assert plain.stdout == traced.stdout, "probe changed NCBI stdout"
    non_probe = b"".join(line for line in traced.stderr.splitlines(keepends=True)
                         if not line.startswith((b"D_CALL\t", b"D_CONTEXT\t",
                                                 b"D_LIST\t", b"D_HSP\t",
                                                 b"D_RETURN\t", b"D_UPDATE_CONTEXT\t")))
    assert plain.stderr == non_probe, "probe changed NCBI diagnostic stderr"
    prefix = f"{fixture.name}_{profile}"
    (output / f"{prefix}.out").write_bytes(traced.stdout)
    (output / f"{prefix}.trace").write_bytes(traced.stderr)
    calls = [line for line in traced.stderr.decode().splitlines()
             if line.startswith("D_CALL\t")]
    (output / f"{prefix}.manifest.txt").write_text(
        f"Pinned NCBI source: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        f"Probe source SHA256: {sha(HERE / 'ncbi_d_call_trace.c')}\n"
        f"Query SHA256: {sha(fixture / 'query.faa')}\n"
        f"Subject SHA256: {sha(fixture / 'subjects.fna')}\n"
        f"Profile: {profile}\n"
        f"Command: {' '.join(command)}\n"
        f"Final output bytes unchanged by probe: yes\n"
        f"D calls: {len(calls)}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("fixture_dirs", nargs="+", type=Path)
    args = parser.parse_args()
    assert sha(NCBI) == NCBI_SHA256
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix="tlosan-d-probe-") as temp:
        probe = Path(temp) / "ncbi_d_call_trace.so"
        subprocess.run(["gcc", "-shared", "-fPIC", "-std=c11", "-O2",
                        "-Wall", "-Wextra", "-Werror", "-o", str(probe),
                        str(HERE / "ncbi_d_call_trace.c"), "-ldl"], check=True)
        for fixture_arg in args.fixture_dirs:
            fixture = fixture_arg.resolve()
            for profile in ("default", "control"):
                run(profile, fixture, output, probe)
    paths = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in paths))
    print(f"{len(args.fixture_dirs)} fixtures x 2 profiles; probe output parity passed")


if __name__ == "__main__":
    main()
