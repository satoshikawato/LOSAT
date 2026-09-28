#!/usr/bin/env python3
"""Stage G fixed-input TBLASTN command benchmark; correctness is gated separately."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[3]
QUERY_SOURCE = ROOT / "LOSAT/tests/fasta/AvCLPV.faa"
SUBJECT = ROOT / "LOSAT/tests/fasta/AvCLPV.fasta"
NATIVE = ROOT / "LOSAT/target/release/LOSAT"
SERIAL = ROOT / "LOSAT/target/serial-command/wasm32-wasip1/release/LOSAT.wasm"
THREADED = ROOT / "LOSAT/target/threaded-command/wasm32-wasip1-threads/release/LOSAT.wasm"
SERIAL_RUNNER = ROOT / "LOSAT/tests/run_losat_wasi.js"
THREADED_RUNNER = ROOT / "LOSAT/tests/run_losat_wasi_threads.js"
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
MAKEBLASTDB = Path("/home/kawato/micromamba/bin/makeblastdb")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def first_record(path: Path) -> bytes:
    lines = path.read_bytes().splitlines(keepends=True)
    assert lines and lines[0].startswith(b">")
    end = next((i for i, line in enumerate(lines[1:], 1) if line.startswith(b">")), len(lines))
    return b"".join(lines[:end])


def execute(command: list[str], out: Path, label: str, iteration: int) -> dict:
    # NCBI c++/src/app/blast/tblastn_app.cpp:290-301:
    # formatter.PrintOneResultSet writes the selected -out destination per run.
    # Remove the shared sink first so every checksum belongs to this command.
    output = out / "benchmark_result.out"
    output.unlink(missing_ok=True)
    started = time.perf_counter_ns()
    result = subprocess.run(command, cwd=ROOT, capture_output=True)
    elapsed = (time.perf_counter_ns() - started) / 1e9
    (out / f"{label}_{iteration}.stdout").write_bytes(result.stdout)
    (out / f"{label}_{iteration}.stderr").write_bytes(result.stderr)
    record = {
        "mode": label, "iteration": iteration,
        "timed": iteration > 0, "command": command, "elapsed_s": elapsed,
        "exit": result.returncode, "stdout_sha256": sha(out / f"{label}_{iteration}.stdout"),
        "stderr_sha256": sha(out / f"{label}_{iteration}.stderr"),
        "output_sha256": sha(output) if output.exists() else None,
        "output_bytes": output.stat().st_size if output.exists() else None,
    }
    if result.returncode != 0 or not output.exists():
        raise RuntimeError(record)
    return record


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    query = out / "first_AvCLPV_protein.faa"
    query.write_bytes(first_record(QUERY_SOURCE))
    assert all(p.is_file() for p in (SUBJECT, NATIVE, SERIAL, THREADED, NCBI, MAKEBLASTDB))
    db = out / "ncbi_db" / "AvCLPV"
    db.parent.mkdir()
    db_command = [str(MAKEBLASTDB), "-in", str(SUBJECT), "-dbtype", "nucl", "-out", str(db)]
    db_started = time.perf_counter_ns()
    built = subprocess.run(db_command, cwd=ROOT, capture_output=True)
    db_elapsed = (time.perf_counter_ns() - db_started) / 1e9
    (out / "makeblastdb.stdout").write_bytes(built.stdout)
    (out / "makeblastdb.stderr").write_bytes(built.stderr)
    assert built.returncode == 0
    output = out / "benchmark_result.out"
    local = ["tblastn", "-task", "tblastn", "-query", str(query), "-subject", str(SUBJECT),
             "-db_gencode", "1", "-outfmt", "6", "-out", str(output)]
    database = ["-task", "tblastn", "-query", str(query), "-db", str(db),
                "-db_gencode", "1", "-outfmt", "6", "-out", str(output)]
    node = ["node", "--no-warnings", "--experimental-wasi-unstable-preview1"]
    modes = {
        "losat_native_1": [str(NATIVE), *local, "-num_threads", "1"],
        "losat_native_4": [str(NATIVE), *local, "-num_threads", "4"],
        "losat_wasi_serial_1": [*node, str(SERIAL_RUNNER), str(SERIAL), *local, "-num_threads", "1"],
        "losat_wasi_threads_4": [*node, str(THREADED_RUNNER), str(THREADED), *local, "-num_threads", "4"],
        "ncbi_db_1": [str(NCBI), *database, "-num_threads", "1"],
        "ncbi_db_4": [str(NCBI), *database, "-num_threads", "4"],
    }
    environment = {
        "query_source": str(QUERY_SOURCE), "query_source_sha256": sha(QUERY_SOURCE),
        "extracted_query": str(query), "extracted_query_sha256": sha(query),
        "subject": str(SUBJECT), "subject_sha256": sha(SUBJECT),
        "outfmt": 6, "output_sink": str(output),
        "executable_sha256": {name: sha(path) for name, path in
                              (("losat_native", NATIVE), ("losat_serial_wasm", SERIAL),
                               ("losat_threaded_wasm", THREADED), ("ncbi_tblastn", NCBI),
                               ("makeblastdb", MAKEBLASTDB))},
        "ncbi_version": subprocess.check_output([str(NCBI), "-version"], text=True).strip(),
        "makeblastdb_version": subprocess.check_output([str(MAKEBLASTDB), "-version"], text=True).strip(),
        "db_build": {"command": db_command, "elapsed_s": db_elapsed,
                     "exit": built.returncode, "subject_sha256": sha(SUBJECT),
                     "stdout_sha256": sha(out / "makeblastdb.stdout"),
                     "stderr_sha256": sha(out / "makeblastdb.stderr")},
        "output_contract": "NCBI -db timing is separate from LOSAT local -subject parity and distribution outputs",
    }
    (out / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
    rows = []
    with (out / "measurements.jsonl").open("w") as log:
        for label, command in modes.items():
            hashes = []
            # NCBI c++/src/app/blast/tblastn_app.cpp:290-301:
            # lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
            # results = lcl_blast.Run();
            # formatter.PrintOneResultSet(**result, query);
            # Time search/report commands only. Database preparation is already complete.
            for iteration in range(4):
                record = execute(command, out, label, iteration)
                hashes.append(record["output_sha256"])
                rows.append(record)
                log.write(json.dumps(record, sort_keys=True) + "\n")
                log.flush()
            assert len(set(hashes)) == 1, (label, hashes)
            print(label, "warmup + 3 timed PASS", flush=True)
    summary = {}
    for label in modes:
        samples = [row["elapsed_s"] for row in rows if row["mode"] == label and row["timed"]]
        assert len(samples) == 3
        summary[label] = {"timed_samples_s": samples, "median_s": statistics.median(samples),
                          "min_s": min(samples), "max_s": max(samples),
                          "output_sha256": next(row["output_sha256"] for row in rows if row["mode"] == label)}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
