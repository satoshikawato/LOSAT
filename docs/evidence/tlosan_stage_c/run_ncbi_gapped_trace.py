#!/usr/bin/env python3
"""Trace pinned NCBI's GetGappedScore for a Stage C local-subject fixture."""
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
FIELDS = "6 qseqid sseqid score qstart qend sstart send sframe qseq sseq"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--lcase-masking", action="store_true")
    args = parser.parse_args()
    fixture = args.fixture_dir.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert sha(NCBI) == NCBI_SHA256
    expected = fixture / ("lowercase_fields.out" if args.lcase_masking
                          else "raw_isolation_fields.out")
    assert expected.is_file(), expected
    with tempfile.TemporaryDirectory(prefix="tlosan-gapped-probe-") as tmp:
        probe = Path(tmp) / "ncbi_gapped_trace.so"
        subprocess.run(["gcc", "-shared", "-fPIC", "-std=c11", "-O2",
                        "-Wall", "-Wextra", "-o", str(probe),
                        str(HERE / "ncbi_gapped_trace.c"), "-ldl"], check=True)
        command = [str(NCBI), "-task", "tblastn", "-query",
                   str(fixture / "query.faa"), "-subject",
                   str(fixture / "subjects.fna"), "-db_gencode", "1",
                   "-matrix", "BLOSUM62", "-word_size", "3",
                   "-threshold", "13", "-window_size", "40", "-gapopen", "11",
                   "-gapextend", "1", "-evalue", "10000", "-num_threads", "1",
                   "-comp_based_stats", "0", "-seg", "no", "-sum_stats", "false",
                   "-outfmt", FIELDS]
        if args.lcase_masking:
            command.append("-lcase_masking")
        env = os.environ.copy()
        env["LD_PRELOAD"] = str(probe)
        result = subprocess.run(command, env=env, capture_output=True, check=True)
    assert result.stdout == expected.read_bytes(), "probe changed NCBI final output"
    (out / "ncbi_output.out").write_bytes(result.stdout)
    (out / "ncbi_gapped_trace.stderr").write_bytes(result.stderr)
    initial = (fixture / "ncbi_init_trace.tsv").read_text().splitlines()
    assert initial[0].split("\t") == [
        "call", "subject", "frame", "init_index", "q_seed", "s_seed",
        "q_start", "s_start", "length", "raw_score"]
    groups: list[tuple[str, str, list[list[str]]]] = []
    for line in initial[1:]:
        f = line.split("\t")
        key = (f[1], f[2])
        if not groups or groups[-1][:2] != key:
            groups.append((f[1], f[2], []))
        groups[-1][2].append(f)
    lines = result.stderr.decode().splitlines()
    inputs = [line.split("\t") for line in lines if line.startswith("GAPPED_INPUT\t")]
    inits = [line.split("\t") for line in lines if line.startswith("GAPPED_INIT\t")]
    outputs = [line.split("\t") for line in lines if line.startswith("GAPPED_OUTPUT\t")]
    hsps = [line.split("\t") for line in lines if line.startswith("GAPPED_HSP\t")]
    assert len(inputs) == len(outputs) == len(groups)
    assert len(inits) == len(initial) - 1
    input_rows = ["call\tsubject\tframe\tindex\tq_seed\ts_seed\tq_start\ts_start\tlength\traw_score\n"]
    for call, (subject, frame, expected_rows) in enumerate(groups):
        assert inputs[call] == ["GAPPED_INPUT", str(call), str(len(expected_rows))]
        actual_rows = [row for row in inits if int(row[1]) == call]
        assert len(actual_rows) == len(expected_rows)
        for index, (actual, expected_row) in enumerate(zip(actual_rows, expected_rows)):
            assert actual[2:] == [str(index), *expected_row[4:]], (actual, expected_row)
            input_rows.append("\t".join([str(call), subject, frame, *actual[2:]]) + "\n")
    (out / "gapped_input.tsv").write_text("".join(input_rows))
    output_rows = ["call\tsubject\tframe\tindex\traw_score\tcontext\tquery_frame\tq_start\tq_end\tq_gapped_start\tsubject_frame\ts_start\ts_end\ts_gapped_start\n"]
    for call, (subject, frame, _) in enumerate(groups):
        output = outputs[call]
        assert output[1:3] == [str(call), "0"], output
        actual_rows = [row for row in hsps if int(row[1]) == call]
        assert len(actual_rows) == int(output[3])
        for index, row in enumerate(actual_rows):
            assert row[2] == str(index)
            assert row[9] == frame
            output_rows.append("\t".join([str(call), subject, frame, *row[2:]]) + "\n")
    (out / "gapped_output.tsv").write_text("".join(output_rows))
    names = [line[1:].split()[0] for line in
             (fixture / "subjects.fna").read_text().splitlines()
             if line.startswith(">")]
    traceback_calls: list[dict] = []
    passes: dict[int, int] = {}
    for line in lines:
        row = line.split("\t")
        if row[0] == "TRACEBACK_INPUT":
            oid = int(row[1])
            assert 0 <= oid < len(names)
            passes[oid] = passes.get(oid, 0) + 1
            traceback_calls.append({
                "oid": oid, "pass": passes[oid], "before_count": int(row[2]),
                "before": [], "after": [],
            })
        elif row[0] == "TRACEBACK_IN_HSP":
            assert traceback_calls and int(row[1]) == traceback_calls[-1]["oid"]
            traceback_calls[-1]["before"].append(row[2:])
        elif row[0] == "TRACEBACK_OUTPUT":
            assert traceback_calls and int(row[1]) == traceback_calls[-1]["oid"]
            traceback_calls[-1].update(
                status=int(row[2]), after_count=int(row[3]), fence=int(row[4]))
        elif row[0] == "TRACEBACK_OUT_HSP":
            assert traceback_calls and int(row[1]) == traceback_calls[-1]["oid"]
            traceback_calls[-1]["after"].append(row[2:])
    assert traceback_calls, "traceback entry was not intercepted"
    tb_rows = [
        "call\toid\tsubject\tpass\tfence\tstatus\tcount_before\tcount_after\tindex\t"
        "before_raw\tbefore_frame\tbefore_q_start\tbefore_q_end\tbefore_s_start\tbefore_s_end\t"
        "after_raw\tafter_frame\tafter_q_start\tafter_q_end\tafter_s_start\tafter_s_end\n"
    ]
    for call, entry in enumerate(traceback_calls):
        assert entry["before_count"] == len(entry["before"])
        assert entry["after_count"] == len(entry["after"])
        assert entry["status"] == 0
        for index in range(max(entry["before_count"], entry["after_count"])):
            before = entry["before"][index] if index < len(entry["before"]) else ["-"] * 7
            after = entry["after"][index] if index < len(entry["after"]) else ["-"] * 7
            assert before[0] == str(index) or before[0] == "-"
            assert after[0] == str(index) or after[0] == "-"
            tb_rows.append("\t".join([
                str(call), str(entry["oid"]), names[entry["oid"]],
                str(entry["pass"]), str(entry["fence"]), str(entry["status"]),
                str(entry["before_count"]), str(entry["after_count"]), str(index),
                *before[1:], *after[1:]]) + "\n")
    (out / "traceback_hsps.tsv").write_text("".join(tb_rows))
    (out / "manifest.txt").write_text(
        "Comparison-only LD_PRELOAD NCBI trace; never used by LOSAT runtime/build.\n"
        "NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        f"Probe source SHA256: {sha(HERE / 'ncbi_gapped_trace.c')}\n"
        f"Fixture: {fixture}\n"
        f"Fixture output SHA256: {sha(expected)}\n"
        f"GetGappedScore calls: {len(groups)}; initial HSPs: {len(inits)}; output HSPs: {len(hsps)}\n"
        f"Traceback calls: {len(traceback_calls)}; fence retries: {sum(x['fence'] for x in traceback_calls)}\n"
        f"Final NCBI bytes unchanged by probe: yes\n"
        f"Command: {command!r}\n")
    paths = ["ncbi_output.out", "ncbi_gapped_trace.stderr",
             "gapped_input.tsv", "gapped_output.tsv", "traceback_hsps.tsv",
             "manifest.txt"]
    (out / "outputs.sha256").write_text(
        "".join(f"{sha(out / name)}  {name}\n" for name in paths))
    print(f"{len(groups)} gapped calls, {len(inits)} initial HSPs, {len(hsps)} output HSPs; final bytes unchanged")


if __name__ == "__main__":
    main()
