#!/usr/bin/env python3
"""Capture LOSAT output hashes around the LOSAT Web engine refactors (S02-S04, SX).

The refactors must not change a single output byte. This script runs one LOSAT CLI
binary over the existing regression cases and records, per case, the SHA-256 of the
output and of the normalized stderr, and the exit code. Two captures (before and after
a change) are compared with `compare`.

Case sources (reused, not copied):
- BLASTN: LOSAT/tests/blastn_parity_manifest.tsv via compare_blastn_parity.
- BLASTP: LOSAT/tests/blastp_v010_parity_manifest.tsv via audit_blastp_v010, in the
  manifest format and additionally in outfmt 0, 7 and one custom tabular field list.
- TBLASTX: LOSAT/tests/tblastx_v010_parity_manifest.tsv via audit_tblastx_v010.
- TBLASTN: docs/evidence/tlosan_stage_g/matrix_162.jsonl (Stage G commands).
- BLASTX: a small set built from LOSAT/tests/fasta (outfmt 0/6/7, threads 1 and 4).

Where frozen expected hashes exist they are checked too: Gate A
(LOSAT/tests/platform_native_v010_canonical.tsv) for the BLASTN, BLASTP and TBLASTX
manifest cases, and the Stage G expected hashes for TBLASTN. Output formats 0 and 7
print the input paths, so inputs are passed with the path strings used when those
hashes were frozen: the Gate A lexical fixture root of certify_platform_native_v010
(HISTORICAL_LEXICAL_ROOT) and the absolute paths recorded by Stage G. Before a case
runs, every input read through such a path must be byte-identical to the file in this
checkout; otherwise the capture stops.

Usage:
  capture_outputs.py run --losat BIN --out DIR [--programs blastn,blastp,...] [--jobs N]
  capture_outputs.py compare BEFORE.tsv AFTER.tsv
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import dataclasses
import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TESTS = REPO / "LOSAT" / "tests"
sys.path.insert(0, str(TESTS))

import audit_blastp_v010  # noqa: E402
import audit_tblastx_v010  # noqa: E402
import certify_platform_native_v010 as certify  # noqa: E402
import compare_blastn_parity  # noqa: E402

PROGRAMS = ("blastn", "blastp", "tblastx", "tblastn", "blastx")
FIELDS = [
    "program", "case_id", "command", "exit", "output_sha256", "stderr_sha256",
    "expected_sha256", "expected_source", "matches_expected",
]
BLASTP_CUSTOM_FIELDS = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue "
    "bitscore score positive gaps ppos qframe sframe qseq sseq btop"
)
BLASTX_CASES = (
    ("mjenmv_self", "LOSAT/tests/fasta/MjeNMV.fasta", "LOSAT/tests/fasta/MjeNMV.faa"),
    ("small_vs_paja", "LOSAT/tests/fasta/small_test.fasta", "LOSAT/tests/fasta/PajaWSV.faa"),
)


@dataclasses.dataclass(frozen=True)
class Case:
    program: str
    case_id: str
    command: list[str]  # argv after the binary; "{OUT}" marks the output file, if any
    cwd: Path
    expected_sha256: str = ""
    expected_source: str = ""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def gate_a() -> dict[tuple[str, str], str]:
    rows = {}
    with open(TESTS / "platform_native_v010_canonical.tsv", newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        rows[(row["program"], row["case_id"])] = row["losat_sha256"]
    return rows


def strip_binary(command: list[str]) -> list[str]:
    return command[1:]


def require_same_file(lexical: str, local: Path) -> None:
    if sha256_bytes(Path(lexical).read_bytes()) != sha256_bytes(local.read_bytes()):
        raise SystemExit(f"input {lexical} differs from {local}")


def lexical_inputs(command: list[str], program: str) -> list[str]:
    """Rewrites absolute inputs under LOSAT/ to the Gate A lexical fixture root."""
    adapted = certify._adapt_fixture_arguments(command, program, REPO)
    for original, lexical in zip(command, adapted):
        if original != lexical:
            require_same_file(lexical, Path(original))
    return adapted


def blastn_cases(expected: dict[tuple[str, str], str]) -> list[Case]:
    rows = compare_blastn_parity.read_manifest(TESTS / "blastn_parity_manifest.tsv")
    cases = []
    for row in rows:
        absolute = certify._blastn_command_row(REPO, row)
        command = compare_blastn_parity.build_losat_command(absolute, Path("LOSAT"), "{OUT}")
        command = strip_binary(lexical_inputs(command, "blastn"))
        key = ("blastn", row["case_id"])
        cases.append(Case("blastn", row["case_id"], command, REPO,
                          expected.get(key, ""), "gate_a" if key in expected else ""))
    return cases


def blastp_cases(expected: dict[tuple[str, str], str]) -> list[Case]:
    manifest = audit_blastp_v010.load_manifest(TESTS / "blastp_v010_parity_manifest.tsv", REPO)
    cases = []
    for base in manifest:
        variants = [(base.case_id, base.outfmt), (f"{base.case_id}.fmt0", "0"), (f"{base.case_id}.fmt7", "7")]
        if base is manifest[0]:
            variants.append((f"{base.case_id}.fmt6custom", BLASTP_CUSTOM_FIELDS))
        for case_id, outfmt in variants:
            case = dataclasses.replace(base, outfmt=outfmt)
            _, losat = audit_blastp_v010.build_commands(case, Path("unused"), Path("LOSAT"), Path("{OUTDIR}"))
            command = strip_binary(lexical_inputs(losat, "blastp"))
            command[command.index("-out") + 1] = "{OUT}"
            key = ("blastp", case_id)
            cases.append(Case("blastp", case_id, command, REPO,
                              expected.get(key, ""), "gate_a" if key in expected else ""))
    return cases


def tblastx_cases(expected: dict[tuple[str, str], str]) -> list[Case]:
    manifest = audit_tblastx_v010.load_manifest(TESTS / "tblastx_v010_parity_manifest.tsv", REPO)
    cases = []
    for case in manifest:
        _, losat = audit_tblastx_v010.build_commands(case, Path("unused"), Path("LOSAT"), Path("unused"), Path("{OUT}"))
        key = ("tblastx", case.case_id)
        cases.append(Case("tblastx", case.case_id, strip_binary(lexical_inputs(losat, "tblastx")), REPO,
                          expected.get(key, ""), "gate_a" if key in expected else ""))
    return cases


def tblastn_cases() -> list[Case]:
    prefix = "/mnt/c/Users/genom/GitHub/LOSAT/"
    cases = []
    with open(REPO / "docs/evidence/tlosan_stage_g/matrix_162.jsonl") as handle:
        for line in handle:
            row = json.loads(line)
            command = row["losat_command"][1:]
            for arg in command:
                if arg.startswith(prefix):
                    require_same_file(arg, REPO / arg.removeprefix(prefix))
            case_id = f"{row['case']}.fmt{row['outfmt']}.n{row['threads']}.r{row['repetition']}"
            cases.append(Case("tblastn", case_id, command, REPO, row["expected_sha256"], "tlosan_stage_g"))
    return cases


def blastx_cases() -> list[Case]:
    cases = []
    for name, query, subject in BLASTX_CASES:
        for outfmt in ("0", "6", "7"):
            for threads in ("1", "4"):
                command = ["blastx", "-query", query, "-subject", subject, "-outfmt", outfmt,
                           "-num_threads", threads, "-out", "{OUT}"]
                cases.append(Case("blastx", f"{name}.fmt{outfmt}.n{threads}", command, REPO))
    return cases


def all_cases(programs: set[str]) -> list[Case]:
    expected = gate_a()
    builders = {
        "blastn": lambda: blastn_cases(expected),
        "blastp": lambda: blastp_cases(expected),
        "tblastx": lambda: tblastx_cases(expected),
        "tblastn": tblastn_cases,
        "blastx": blastx_cases,
    }
    return [case for program in PROGRAMS if program in programs for case in builders[program]()]


def run_case(case: Case, losat: Path, out_dir: Path) -> dict[str, str]:
    output = out_dir / "outputs" / case.program / f"{case.case_id}.out"
    output.parent.mkdir(parents=True, exist_ok=True)
    argv = [str(losat)] + [str(output) if arg == "{OUT}" else arg for arg in case.command]
    result = subprocess.run(argv, cwd=case.cwd, capture_output=True)
    uses_file = "{OUT}" in case.command
    data = output.read_bytes() if uses_file and output.exists() else (b"" if uses_file else result.stdout)
    if not uses_file:
        output.write_bytes(result.stdout)
    stderr = result.stderr.replace(str(output).encode(), b"{OUT}").replace(str(losat).encode(), b"{LOSAT}")
    digest = sha256_bytes(data)
    return {
        "program": case.program,
        "case_id": case.case_id,
        "command": json.dumps(case.command),
        "exit": str(result.returncode),
        "output_sha256": digest,
        "stderr_sha256": sha256_bytes(stderr),
        "expected_sha256": case.expected_sha256,
        "expected_source": case.expected_source,
        "matches_expected": "" if not case.expected_sha256 else str(digest == case.expected_sha256).lower(),
    }


def command_run(args: argparse.Namespace) -> int:
    losat = Path(args.losat).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    programs = set(args.programs.split(",")) if args.programs else set(PROGRAMS)
    cases = all_cases(programs)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(lambda case: run_case(case, losat, out_dir), cases))
    with open(out_dir / "hashes.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (out_dir / "binary.json").write_text(json.dumps({
        "losat": str(losat), "losat_sha256": sha256_bytes(losat.read_bytes()), "cases": len(rows),
    }, indent=2) + "\n")
    mismatched = [row for row in rows if row["matches_expected"] == "false"]
    for row in mismatched:
        print(f"expected-hash mismatch: {row['program']} {row['case_id']}", file=sys.stderr)
    print(f"{len(rows)} cases, {len(mismatched)} frozen-hash mismatches -> {out_dir / 'hashes.tsv'}")
    return 1 if mismatched else 0


def read_rows(path: str) -> dict[tuple[str, str], dict[str, str]]:
    with open(path, newline="") as handle:
        return {(row["program"], row["case_id"]): row for row in csv.DictReader(handle, delimiter="\t")}


def command_compare(args: argparse.Namespace) -> int:
    before, after = read_rows(args.before), read_rows(args.after)
    problems = [f"missing after: {key}" for key in before if key not in after]
    problems += [f"new after: {key}" for key in after if key not in before]
    for key in sorted(set(before) & set(after)):
        for field in ("command", "exit", "output_sha256", "stderr_sha256"):
            if before[key][field] != after[key][field]:
                problems.append(f"{field} differs: {key}")
    for line in problems:
        print(line)
    print(f"{len(set(before) & set(after))} shared cases, {len(problems)} differences")
    return 1 if problems else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    run = sub.add_parser("run")
    run.add_argument("--losat", required=True)
    run.add_argument("--out", required=True)
    run.add_argument("--programs", default="")
    run.add_argument("--jobs", type=int, default=6)
    compare = sub.add_parser("compare")
    compare.add_argument("before")
    compare.add_argument("after")
    args = parser.parse_args()
    return command_run(args) if args.action == "run" else command_compare(args)


if __name__ == "__main__":
    raise SystemExit(main())
