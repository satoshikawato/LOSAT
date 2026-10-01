#!/usr/bin/env python3
"""Frozen NCBI BLAST+ 2.17.0 BLASTN outputs for fast pull-request regressions.

The cases cover NCBI paths that the Gate A manifest and the outfmt 0 fixtures do
not reach: preliminary hit lists that overflow (more subjects with hits than
prelim_hitlist_size, ties broken by subject order), query batches, a split query
batch, IUPAC ambiguity, lowercase masking and non-default scores.

- generate: writes the inputs to LOSAT/tests/fixtures/blastn_regression/inputs/
  (deterministic; the files are committed, so later runs never regenerate them).
- freeze --oracle BLASTN: runs NCBI BLAST+ (comparison oracle only) from LOSAT/ and
  writes <case>.out, <case>.err and the hash columns of manifest.tsv.
- check --losat LOSAT: runs `LOSAT blastn <argv> <losat_extra>` from LOSAT/ and
  compares stdout, stderr and the exit status with the frozen files.

A case may set environment variables (`env`, space-separated KEY=VALUE) for both
programs; every other variable that changes NCBI's batches or report is unset.

The split case slices LOSAT/tests/fasta/EDL933.fna at run time into
LOSAT/target/blastn_regression/ (outfmt 6 prints no file names).

Usage:
  blastn_regression_fixtures.py generate
  blastn_regression_fixtures.py freeze --oracle /path/to/blastn
  blastn_regression_fixtures.py check --losat BIN [--jobs N] [--out TSV]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import os
import random
import shlex
import subprocess
import sys
from pathlib import Path

ENGINE = Path(__file__).resolve().parents[1]
FIXTURES = ENGINE / "tests/fixtures/blastn_regression"
INPUTS = FIXTURES / "inputs"
MANIFEST = FIXTURES / "manifest.tsv"
RUNTIME = ENGINE / "target/blastn_regression"
FASTA = ENGINE / "tests/fasta"
E2C_INPUTS = ENGINE.parent / "docs/evidence/losat_web_e2c/inputs"
FIELDS = ["case_id", "argv", "losat_extra", "exit", "stdout_sha256", "stdout_bytes", "stderr_sha256", "env"]
# NCBI reads these and each one changes the batches or the report.
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE", "ADAPTIVE_CBS")

I = "tests/fixtures/blastn_regression/inputs"
R = "target/blastn_regression"
P = f"-query {I}/q3k.fa -subject {I}/s600.fa"
T = f"-query {I}/q3k.fa -subject {I}/s600ties.fa"
B = f"-query {I}/mq.fa -subject {I}/s60k.fa"
S = f"-query {R}/edl933_2m.fa -subject {I}/sakai_60k.fa"
A = f"-query {I}/ambiguity_query_c.fa -subject {I}/ambiguity_subject.fa"
# (case_id, NCBI argv after `blastn`, extra LOSAT-only arguments[, environment])
_CASES = [
    ("prelim.default", f"{P} -outfmt 6", ""),
    ("prelim.task_blastn", f"{P} -task blastn -outfmt 6", ""),
    ("prelim.max1", f"{P} -max_target_seqs 1 -outfmt 6", ""),
    ("prelim.max3", f"{P} -max_target_seqs 3 -outfmt 6", ""),
    ("prelim.max6", f"{P} -max_target_seqs 6 -outfmt 6", ""),
    ("prelim.max3_besthit", f"{P} -task blastn -max_target_seqs 3 -subject_besthit -outfmt 6", ""),
    ("prelim.max5_maxhsps1", f"{P} -max_target_seqs 5 -max_hsps 1 -outfmt 6", ""),
    ("prelim.max2_evalue", f"{P} -max_target_seqs 2 -evalue 1e-5 -outfmt 6", ""),
    ("prelim.max3_fmt7", f"{P} -max_target_seqs 3 -outfmt 7", ""),
    ("prelim.max3_fmt0", f"{P} -max_target_seqs 3 -outfmt 0", ""),
    ("prelim.default_threads4", f"{P} -outfmt 6", "-num_threads 4"),
    ("prelim.scores_3_4", f"{P} -task blastn -reward 3 -penalty -4 -gapopen 10 -gapextend 3 -outfmt 6", ""),
    ("prelim.scores_1_1", f"{P} -task blastn -reward 1 -penalty -1 -gapopen 3 -gapextend 2 -outfmt 6", ""),
    ("ties.default", f"{T} -outfmt 6", ""),
    ("ties.max3", f"{T} -max_target_seqs 3 -outfmt 6", ""),
    ("ties.max1_task_blastn", f"{T} -task blastn -max_target_seqs 1 -outfmt 6", ""),
    ("batches.default", f"{B} -outfmt 6", ""),
    ("batches.task_blastn", f"{B} -task blastn -outfmt 6", ""),
    ("batches.word7", f"{B} -task blastn -word_size 7 -outfmt 6", ""),
    ("batches.besthit", f"{B} -subject_besthit -outfmt 6", ""),
    ("batches.fmt7", f"{B} -task blastn -outfmt 7", ""),
    ("batches.fmt0", f"{B} -outfmt 0", ""),
    ("batches.threads4", f"{B} -task blastn -outfmt 6", "-num_threads 4"),
    ("split.task_blastn", f"{S} -task blastn -outfmt 6", ""),
    ("split.besthit", f"{S} -task blastn -subject_besthit -outfmt 6", ""),
    ("ambiguity.default", f"{A} -outfmt 6", ""),
    ("ambiguity.task_blastn_word7", f"{A} -task blastn -word_size 7 -outfmt 6", ""),
    ("lcase.blastn", f"-query {I}/lcase_island_blastn_query.fa -subject {I}/lcase_island_blastn_subject.fa"
                     " -task blastn -lcase_masking -outfmt 6", ""),
    ("lcase.megablast", f"-query {I}/lcase_island_megablast_query.fa -subject {I}/lcase_island_megablast_subject.fa"
                        " -lcase_masking -outfmt 6", ""),
    # E2g T1: init hits of the two strands tied on score, subject start and length.
    ("pal.default", f"-query {I}/pal_query.fa -subject {I}/pal_subject.fa -outfmt 6", ""),
    ("pal.task_blastn", f"-query {I}/pal_query.fa -subject {I}/pal_subject.fa -task blastn -outfmt 6", ""),
    ("pal.word7_fmt0", f"-query {I}/pal_query.fa -subject {I}/pal_subject.fa -task blastn -word_size 7 -outfmt 0", ""),
    # E2g T11: showdefline.cpp kBits is "(bits)" when CTOOLKIT_COMPATIBLE is set (also empty).
    ("ctoolkit.fmt0", f"{P} -max_target_seqs 3 -outfmt 0", "", "CTOOLKIT_COMPATIBLE=1"),
    ("ctoolkit.empty_fmt0", f"{T} -task blastn -max_target_seqs 5 -outfmt 0", "", "CTOOLKIT_COMPATIBLE="),
]
CASES = [case if len(case) == 4 else (*case, "") for case in _CASES]


def genome(name: str) -> str:
    return "".join(line.strip() for line in (FASTA / name).read_text().splitlines()[1:])


def write_fasta(path: Path, records) -> None:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))


def mutate(rng: random.Random, seq: str, subst: float, indel: float = 0.0) -> str:
    out = []
    for base in seq:
        r = rng.random()
        if r < indel * 0.5:
            continue
        if r < indel:
            out += [base, rng.choice("ACGT")]
        elif rng.random() < subst:
            out.append(rng.choice([c for c in "ACGT" if c != base]))
        else:
            out.append(base)
    return "".join(out)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def command_generate(_args) -> int:
    INPUTS.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261002)
    edl, sakai = genome("EDL933.fna"), genome("Sakai.fna")
    q3k = edl[1_000_000:1_003_000]
    write_fasta(INPUTS / "q3k.fa", [("q3k EDL933 1000001-1003000", q3k)])
    records = []
    for index in range(600):
        length = rng.choice((60, 80, 100, 150, 250))
        start = rng.randrange(0, len(q3k) - length)
        window = mutate(rng, q3k[start:start + length], rng.choice((0.0, 0.0, 0.01, 0.03)))
        if rng.random() < 0.3:
            window = revcomp(window)
        left = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 200)))
        right = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 200)))
        records.append((f"s{index:03d}", left + window + right))
    write_fasta(INPUTS / "s600.fa", records)
    records = []
    for index in range(600):
        start = rng.randrange(0, len(q3k) - 100)
        records.append((f"t{index:03d}", q3k[start:start + 100]))
    write_fasta(INPUTS / "s600ties.fa", records)
    region = sakai[2_000_000:2_060_000]
    write_fasta(INPUTS / "s60k.fa", [("s60k Sakai 2000001-2060000 mutated", mutate(rng, region, 0.01, 0.005))])
    records = []
    for index in range(40):
        length = int(10 * (600 ** rng.random()))
        start = rng.randrange(0, len(region) - length)
        window = mutate(rng, region[start:start + length], rng.choice((0.0, 0.01, 0.03)))
        if rng.random() < 0.3:
            window = revcomp(window)
        records.append((f"mq{index:02d} len{length}", window))
    write_fasta(INPUTS / "mq.fa", records)
    write_fasta(INPUTS / "sakai_60k.fa", [("sakai_60k Sakai 500001-560000", sakai[500_000:560_000])])
    for name in ("ambiguity_query_c.fa", "ambiguity_subject.fa", "lcase_island_blastn_query.fa",
                 "lcase_island_blastn_subject.fa", "lcase_island_megablast_query.fa",
                 "lcase_island_megablast_subject.fa"):
        (INPUTS / name).write_bytes((E2C_INPUTS / name).read_bytes())
    generate_e2g_inputs(edl, sakai)
    return 0


def generate_e2g_inputs(edl: str, sakai: str) -> None:
    """Inputs added in E2g; a separate generator keeps the earlier inputs unchanged."""
    rng = random.Random(20261003)
    # T1: query X + revcomp(X); its minus-strand context is the same sequence, so both
    # strands give hits with equal score, subject start and length (context-local q_start).
    records, subject = [], []
    for index in range(40):
        start = rng.randrange(0, len(edl) - 400)
        x = edl[start:start + rng.choice((60, 90, 120, 150))]
        records.append((f"pal{index:02d} EDL933 {start + 1} X+revcomp(X)", x + revcomp(x)))
        subject.append("".join(rng.choice("ACGT") for _ in range(rng.randint(20, 120))))
        subject.append(mutate(rng, x, rng.choice((0.0, 0.02))))
    write_fasta(INPUTS / "pal_query.fa", records)
    write_fasta(INPUTS / "pal_subject.fa", [("pals EDL933 windows", "".join(subject))])


def prepare_runtime_inputs() -> None:
    RUNTIME.mkdir(parents=True, exist_ok=True)
    path = RUNTIME / "edl933_2m.fa"
    expected = ">edl933_2m EDL933 1-2000000\n" + genome("EDL933.fna")[:2_000_000] + "\n"
    if not path.exists() or path.read_text() != expected:
        path.write_text(expected)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def clean_env(case_env: str = "") -> dict[str, str]:
    env = {key: value for key, value in os.environ.items()
           if key not in REPORT_ENV and not key.startswith(("LOSAT_", "RAYON_"))}
    for item in shlex.split(case_env):
        key, _, value = item.partition("=")
        env[key] = value
    return env


def read_manifest() -> list[dict[str, str]]:
    lines = [line for line in MANIFEST.read_text().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def command_freeze(args) -> int:
    if any(key in os.environ for key in REPORT_ENV) or (Path.home() / ".ncbirc").exists():
        raise SystemExit(f"unset {REPORT_ENV} and remove ~/.ncbirc before freezing")
    prepare_runtime_inputs()
    oracle = str(Path(args.oracle).resolve())
    version = subprocess.run([oracle, "-version"], capture_output=True, text=True).stdout.strip().replace("\n", "; ")
    rows = []
    for case_id, argv, extra, case_env in CASES:
        result = subprocess.run([oracle, *shlex.split(argv)], cwd=ENGINE, capture_output=True, env=clean_env(case_env))
        (FIXTURES / f"{case_id}.out").write_bytes(result.stdout)
        err = FIXTURES / f"{case_id}.err"
        if result.stderr:
            err.write_bytes(result.stderr)
        elif err.exists():
            err.unlink()
        rows.append({"case_id": case_id, "argv": argv, "losat_extra": extra, "exit": str(result.returncode),
                     "stdout_sha256": sha256(result.stdout), "stdout_bytes": str(len(result.stdout)),
                     "stderr_sha256": sha256(result.stderr) if result.stderr else "", "env": case_env})
        print(f"{case_id}\texit {result.returncode}\t{result.stdout.count(b'\n')} lines\t{len(result.stderr)} stderr bytes", flush=True)
    with open(MANIFEST, "w", newline="") as handle:
        handle.write(f"# Frozen BLASTN outputs of {version} (comparison oracle only; see blastn_regression_fixtures.py).\n")
        handle.write("# Run from LOSAT/: env <env> blastn <argv> (LOSAT: LOSAT blastn <argv> <losat_extra>). Written by `freeze`; do not edit.\n")
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


def check_one(losat: str, row: dict[str, str]) -> tuple[str, str]:
    argv = [losat, "blastn", *shlex.split(row["argv"]), *shlex.split(row["losat_extra"])]
    result = subprocess.run(argv, cwd=ENGINE, capture_output=True, env=clean_env(row.get("env") or ""))
    expected_out = (FIXTURES / f"{row['case_id']}.out").read_bytes()
    err_path = FIXTURES / f"{row['case_id']}.err"
    expected_err = err_path.read_bytes() if err_path.exists() else b""
    problems = []
    if sha256(expected_out) != row["stdout_sha256"]:
        problems.append("frozen stdout file does not match the manifest")
    if str(result.returncode) != row["exit"]:
        problems.append(f"exit {result.returncode}, expected {row['exit']}")
    if result.stdout != expected_out:
        expected_lines, actual_lines = expected_out.split(b"\n"), result.stdout.split(b"\n")
        line = next((n for n, (a, b) in enumerate(zip(expected_lines, actual_lines), 1) if a != b),
                    min(len(expected_lines), len(actual_lines)))
        problems.append(f"stdout differs (expected {len(expected_lines) - 1} lines, got {len(actual_lines) - 1}; first difference at line {line})")
    if result.stderr != expected_err:
        problems.append(f"stderr differs: {result.stderr[:200]!r}")
    return row["case_id"], "; ".join(problems) or "same"


def command_check(args) -> int:
    prepare_runtime_inputs()
    losat = str(Path(args.losat).resolve())
    rows = read_manifest()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda row: check_one(losat, row), rows))
    lines = [f"{case_id}\t{verdict}" for case_id, verdict in results]
    if args.out:
        Path(args.out).write_text("case_id\tresult\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))
    differing = [case_id for case_id, verdict in results if verdict != "same"]
    print(f"{len(results)} cases, {len(differing)} differ")
    return 1 if differing else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("generate")
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--oracle", required=True)
    check = sub.add_parser("check")
    check.add_argument("--losat", required=True)
    check.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    check.add_argument("--out")
    args = parser.parse_args()
    return {"generate": command_generate, "freeze": command_freeze, "check": command_check}[args.action](args)


if __name__ == "__main__":
    sys.exit(main())
