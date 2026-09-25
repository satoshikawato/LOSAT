#!/usr/bin/env python3
"""Comparison-only pinned NCBI TBLASTN outfmt 0/6/7 byte sweep."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import shlex
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
LOSAT = ROOT / "LOSAT/target/release/LOSAT"
SOURCE = Path("/mnt/c/Users/genom/GitHub/ncbi-blast")
PIN = "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4"
NCBI_SHA = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"
FAMILIES = (
    "masking_options_20260925/run_20260925",
    "remaining_local_20260925/run_20260925",
    "option_cross_20260925/run_20260925",
    "alternate_matrix_20260925/run_20260925",
    "natural_positive_20260925/run_20260925",
    "natural_positive_20260925/containment_run_20260925",
    "kappa_heap_rejection_20260925/natural_c_d_early_20260925",
    "kappa_heap_rejection_20260925/natural_mode2_20260925",
    "kappa_heap_rejection_20260925/report_payload_20260925",
    "kappa_heap_rejection_20260925/result_order_20260925",
    "kappa_heap_rejection_20260925/run_20260925",
    "kappa_mode2_20260924",
    "kappa_result_order_20260925",
    "kappa_composition_matrix_scores_20260924",
    "kappa_seg_20260924",
    "kappa_traceback_20260924",
    "erfc_20260924",
    "parameters_20260924",
    "seqsrc_callstate_20260924",
)



def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def normalize_path(value: str) -> str:
    low = "/mnt/c/users/genom/github/losat"
    if value.lower().startswith(low):
        return str(ROOT) + value[len(low):]
    return value


def read_commands():
    # Stage D manifests use "command:", "Command:", and
    # "case: ...; command [...]". The generated 15/30-MB subjects
    # are replayed by run_stage_e_generated.py.
    for family in FAMILIES:
        manifest = ROOT / "docs/evidence/tlosan_stage_d" / family / "manifest.txt"
        for index, line in enumerate(manifest.read_text().splitlines(), 1):
            marker = line.lower().find("command")
            if marker < 0:
                continue
            opening = line.find("[", marker)
            if opening < 0:
                continue
            try:
                command = ast.literal_eval(line[opening:])
            except (ValueError, SyntaxError):
                continue
            if not command or Path(command[0]).name != "tblastn":
                continue
            command = [normalize_path(str(item)) for item in command]
            if not all(arg in command for arg in ("-query", "-subject", "-outfmt")):
                continue
            query = Path(command[command.index("-query") + 1])
            subject = Path(command[command.index("-subject") + 1])
            assert query.is_file() and subject.is_file(), (manifest, index)
            format_index = command.index("-outfmt")
            del command[format_index:format_index + 2]
            name = line[:marker].strip(" :;").replace(" ", "_") or f"command_{index}"
            yield family, name, command
    # The earliest Stage D gate stores one shell-form command per fixture,
    # instead of Python-list commands in a shared manifest.txt.
    for family in ("run_20260924", "uneven_gap_run_20260924"):
        directory = ROOT / "docs/evidence/tlosan_stage_d" / family
        for manifest in sorted(directory.glob("*.manifest.txt")):
            line = next((line for line in manifest.read_text().splitlines()
                         if line.startswith("Command: ")), None)
            assert line is not None, manifest
            command = [normalize_path(item) for item in shlex.split(line.removeprefix("Command: "))]
            assert Path(command[0]).name == "tblastn", manifest
            assert Path(command[command.index("-query") + 1]).is_file(), manifest
            assert Path(command[command.index("-subject") + 1]).is_file(), manifest
            format_index = command.index("-outfmt")
            format_end = next((i for i in range(format_index + 2, len(command))
                               if command[i].startswith("-")), len(command))
            del command[format_index:format_end]
            yield family, manifest.name.removesuffix(".manifest.txt"), command


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    options = parser.parse_args()
    assert subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"]).decode().strip() == PIN
    assert digest(NCBI.read_bytes()) == NCBI_SHA
    assert LOSAT.is_file()
    out = options.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    commands = list(read_commands())
    (out / "environment.json").write_text(
        json.dumps({
            "ncbi_source_commit": PIN,
            "ncbi_executable": str(NCBI),
            "ncbi_executable_sha256": NCBI_SHA,
            "losat_executable": str(LOSAT),
            "losat_executable_sha256": digest(LOSAT.read_bytes()),
            "cases": len(commands),
        }, indent=2) + "\n"
    )
    for family, name, command in commands:
        query = Path(command[command.index("-query") + 1])
        subject = Path(command[command.index("-subject") + 1])
        for fmt in ("0", "6", "7"):
            ncbi_cmd = [str(NCBI), *command[1:], "-outfmt", fmt]
            losat_cmd = [str(LOSAT), "tblastn", *command[1:], "-outfmt", fmt]
            expected = subprocess.run(ncbi_cmd, cwd=ROOT, capture_output=True)
            actual = subprocess.run(losat_cmd, cwd=ROOT, capture_output=True)
            stem = f"{family.replace('/', '__')}__{name}__fmt{fmt}"
            record = {
                "case": stem,
                "query_sha256": digest(query.read_bytes()),
                "subject_sha256": digest(subject.read_bytes()),
                "ncbi_command": ncbi_cmd,
                "losat_command": losat_cmd,
                "ncbi_exit": expected.returncode,
                "losat_exit": actual.returncode,
                "ncbi_stdout_sha256": digest(expected.stdout),
                "losat_stdout_sha256": digest(actual.stdout),
                "ncbi_stderr_sha256": digest(expected.stderr),
                "losat_stderr_sha256": digest(actual.stderr),
                "bytes": len(expected.stdout),
                "equal": expected.returncode == 0 and actual.returncode == 0 and expected.stdout == actual.stdout,
            }
            rows.append(record)
            if not record["equal"]:
                (out / "first_failure.json").write_text(json.dumps(record, indent=2) + "\n")
                (out / "first_failure.ncbi.out").write_bytes(expected.stdout)
                (out / "first_failure.losat.out").write_bytes(actual.stdout)
                (out / "first_failure.ncbi.err").write_bytes(expected.stderr)
                (out / "first_failure.losat.err").write_bytes(actual.stderr)
                first = next((i for i, (a, b) in enumerate(zip(expected.stdout, actual.stdout)) if a != b),
                             min(len(expected.stdout), len(actual.stdout)))
                print(f"FAIL {stem}: first byte {first}, NCBI {len(expected.stdout)}, LOSAT {len(actual.stdout)}", flush=True)
                break
            print(f"PASS {stem}: {record['bytes']} bytes", flush=True)
        if rows and not rows[-1]["equal"]:
            break
    (out / "comparison.jsonl").write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))
    print(f"{sum(row['equal'] for row in rows)}/{len(rows)} byte-identical", flush=True)
    return 0 if len(rows) == 3 * len(commands) and all(row["equal"] for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
