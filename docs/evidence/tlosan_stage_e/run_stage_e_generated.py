#!/usr/bin/env python3
"""Regenerate Stage D 15/30-MB local subjects and compare outfmt 0/6/7."""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

from run_stage_e_cli_matrix import ROOT, NCBI, LOSAT, SOURCE, PIN, NCBI_SHA, digest
import importlib.util
_generator = ROOT / "docs/evidence/tlosan_stage_d/extended_chunks_20260925/run_ncbi_extended_chunks.py"
_spec = importlib.util.spec_from_file_location("stage_d_generated", _generator)
assert _spec and _spec.loader
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
make_subject = _module.make_subject

STAGE_C = ROOT / "docs/evidence/tlosan_stage_c"
STAGE_D = ROOT / "docs/evidence/tlosan_stage_d"
CASES = ("long_chunk_20260924", "masked_chunk_boundary_20260924",
         "no_range_middle_20260924", "no_range_two_hits_20260924")


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"]).decode().strip() == PIN
    assert digest(NCBI.read_bytes()) == NCBI_SHA
    source = (STAGE_C / "run_20260923/subjects.fna").read_text()
    records = {block.splitlines()[0]: "".join(block.splitlines()[1:])
               for block in source.split(">") if block.strip()}
    insert = records["plus1"].encode()
    inputs = out / "generated_inputs"
    inputs.mkdir()
    rows = []
    commands = []
    for case in CASES:
        subject = inputs / f"{case}.fna"
        subject.write_bytes(make_subject(case, insert))
        query = STAGE_C / case / "query.faa"
        match = re.search(r"Input SHA256: query=([0-9a-f]+) subject=([0-9a-f]+)",
                          (STAGE_C / case / "manifest.txt").read_text())
        assert match and (digest(query.read_bytes()), digest(subject.read_bytes())) == match.groups(), case
        for mode in (0, 2):
            common = ["-task", "tblastn", "-query", str(query), "-subject", str(subject),
                      "-db_gencode", "1", "-num_threads", "1", "-evalue", "10000",
                      "-max_target_seqs", "500", "-comp_based_stats", str(mode),
                      "-sum_stats", "false" if mode == 0 else "true", "-seg", "no"]
            if case != "long_chunk_20260924":
                common.append("-lcase_masking")
            commands.append((f"{case}.mode{mode}", query, subject, common))
    # Stage D long_subject fixture reuses the saved unmasked subject with
    # three protein queries, including its invalid context.
    query = STAGE_C / "long_multi_query_20260924/query.faa"
    subject = inputs / "long_chunk_20260924.fna"
    manifest = (STAGE_D / "long_subject_20260925/run_20260925/manifest.txt").read_text()
    assert digest(query.read_bytes()) in manifest and digest(subject.read_bytes()) in manifest
    commands.append(("long_multi_query_20260924.mode0", query, subject,
                     ["-task", "tblastn", "-query", str(query), "-subject", str(subject),
                      "-db_gencode", "1", "-matrix", "BLOSUM62", "-word_size", "3",
                      "-threshold", "13", "-window_size", "40", "-gapopen", "11",
                      "-gapextend", "1", "-evalue", "10000", "-num_threads", "1",
                      "-comp_based_stats", "0", "-seg", "no", "-sum_stats", "false"]))
    (out / "environment.json").write_text(json.dumps({
        "ncbi_source_commit": PIN, "ncbi_executable_sha256": NCBI_SHA,
        "losat_executable_sha256": digest(LOSAT.read_bytes()), "cases": len(commands),
    }, indent=2) + "\n")
    for name, query, subject, common in commands:
        for fmt in ("0", "6", "7"):
            ncbi_cmd = [str(NCBI), *common, "-outfmt", fmt]
            losat_cmd = [str(LOSAT), "tblastn", *common, "-outfmt", fmt]
            expected = subprocess.run(ncbi_cmd, cwd=ROOT, capture_output=True)
            actual = subprocess.run(losat_cmd, cwd=ROOT, capture_output=True)
            equal = expected.returncode == actual.returncode == 0 and expected.stdout == actual.stdout
            record = {"case": name, "outfmt": fmt, "query_sha256": digest(query.read_bytes()),
                      "subject_sha256": digest(subject.read_bytes()), "ncbi_command": ncbi_cmd,
                      "losat_command": losat_cmd, "ncbi_exit": expected.returncode,
                      "losat_exit": actual.returncode,
                      "ncbi_stdout_sha256": digest(expected.stdout),
                      "losat_stdout_sha256": digest(actual.stdout),
                      "ncbi_stderr_sha256": digest(expected.stderr),
                      "losat_stderr_sha256": digest(actual.stderr),
                      "bytes": len(expected.stdout), "equal": equal}
            rows.append(record)
            if not equal:
                (out / "first_failure.json").write_text(json.dumps(record, indent=2) + "\n")
                (out / "first_failure.ncbi.out").write_bytes(expected.stdout)
                (out / "first_failure.losat.out").write_bytes(actual.stdout)
                (out / "first_failure.ncbi.err").write_bytes(expected.stderr)
                (out / "first_failure.losat.err").write_bytes(actual.stderr)
                offset = next((i for i, (a, b) in enumerate(zip(expected.stdout, actual.stdout)) if a != b),
                              min(len(expected.stdout), len(actual.stdout)))
                print(f"FAIL {name} fmt{fmt}: byte {offset}", flush=True)
                break
            print(f"PASS {name} fmt{fmt}: {len(expected.stdout)} bytes", flush=True)
        if rows and not rows[-1]["equal"]:
            break
    (out / "comparison.jsonl").write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))
    print(f"{sum(row['equal'] for row in rows)}/{len(rows)} byte-identical", flush=True)
    return 0 if len(rows) == len(commands) * 3 and all(row["equal"] for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
