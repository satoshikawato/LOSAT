#!/usr/bin/env python3
"""Sweep `-query_loc` and `-subject_loc` of BLASTN, BLASTP, TBLASTN and TBLASTX against NCBI
(Session S11, E2d; comparison only).

Every case runs NCBI BLAST+ and LOSAT from LOSAT/ on the same inputs with one range (or a
pair of ranges) in each requested output format, and is classified:

- same: both succeed with the same stdout and stderr;
- same-error: both fail with the same exit status, stdout and stderr (NCBI's error);
- losat-rejects: LOSAT fails with a message that names what LOSAT does not support (an
  explicit rejection: a range part that NCBI's `NStr::StringToInt` cannot convert, which
  stops NCBI with a CStringException that names its build's files, exit 255 (plan TD-15);
  an interval without letters, a range that starts just past a record's end, which NCBI
  searches as a sequence without data);
- arg-error: an argument that NCBI's argument parser rejects (USAGE, exit 1) and LOSAT's
  parser rejects (exit 2), approved exception 1 of PD-LOSAT-CLI-NONSEARCH-DIFFERENCES;
- DIFF: anything else, with the first difference.

The cases (`cases`): every spelling of a range (as NCBI's ParseSequenceRange and
StringToInt read it), for each role; ranges at the ends of the records (the start at
1, L-1, L, L+1 and L+2 of a record of L letters, the end before, at and past the record's
end); both ranges together; and the multi-record inputs of the E2d fixtures (records
shorter than the start, skipped or stopping the subject reader). Exits 1 when any case is
DIFF.

Usage: range_sweep.py --bin-dir DIR --losat LOSAT [--programs blastn,...] [--outfmt 0,6,7]
       [--jobs N] [--out TSV]
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ENGINE = REPO / "LOSAT"
F = "tests/fasta/outfmt0"
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE",
              "OVERLAP_CHUNK_SIZE", "ADAPTIVE_CBS", "BLASTDB", "PRE_FETCH_SEQS_LIMIT")
REJECTION_MARKERS = (b"not supported by LOSAT",)
# (query, subject, query record length, subject record length) of one-record inputs.
INPUTS = {
    "blastn": (f"{F}/e2d_n_q.fa", f"{F}/e2d_n_s_plus.fa", 8000, 7000),
    "blastp": (f"{F}/e2d_p_q.faa", f"{F}/e2d_p_sp75m.faa", 674, 871),
    "tblastn": (f"{F}/e2d_p_q.faa", f"{F}/e2d_n_sub8k.fa", 674, 8001),
    "tblastx": (f"{F}/e2d_x_q.fa", f"{F}/e2d_t_ts.fa", 768, 797),
}
# Multi-record inputs: (program, query, subject).
MULTI = [
    ("blastn", f"{F}/e2d_n_qm.fa", f"{F}/e2d_n_s_plus.fa"),
    ("blastn", f"{F}/e2d_n_qidx.fa", f"{F}/e2d_n_s_plus.fa"),
    ("blastn", f"{F}/e2d_n_qn.fa", f"{F}/e2d_n_sub_multi_t.fa"),
    ("blastn", f"{F}/e2d_n_qn.fa", f"{F}/e2d_n_sub_multi.fa"),
    ("blastp", f"{F}/e2d_p_qmulti.faa", f"{F}/e2d_p_s3.faa"),
    ("blastp", f"{F}/e2d_p_qp75.faa", f"{F}/e2d_p_sp_multi.faa"),
    ("blastp", f"{F}/e2d_p_qabshort.faa", f"{F}/e2d_p_s3.faa"),
    ("tblastn", f"{F}/e2d_p_qmulti.faa", f"{F}/e2d_n_sub_multi_t.fa"),
    ("tblastn", f"{F}/e2d_t_qabshort.faa", f"{F}/e2d_n_sub8k.fa"),
    ("tblastx", f"{F}/e2d_x_ab15.fa", f"{F}/e2d_n_s_plus.fa"),
    ("tblastx", f"{F}/e2d_n_qn.fa", f"{F}/e2d_n_sub_multi_t.fa"),
]
SPELLINGS = ["10-20", "1-2", "+1-5", "01-05", "1-+5", "1-2147483647", "1-1", "5-5", "20-10", "2-1",
             "0-10", "10-0", "-10", "10-", "10", "1--5", "1-2-3", "", "-", "--", "-5-10", " 1-5",
             "1- 5", "1-5 ", "a-5", "1-5a", "1.0-5", "0x10-0x20", "1,0-50", "0-a", "1-2147483648",
             "2147483648-1", "1-9999999999999999999"]
MULTI_RANGES = ["1-300", "100-400", "301-900", "502-1000", "600-5000", "700-900", "1001-2500",
                "2001-3000", "5000-25000", "100-300"]


def edge_ranges(length: int) -> list[str]:
    starts = sorted({1, max(1, length - 1), length, length + 1, length + 2})
    ranges = []
    for start in starts:
        for end in sorted({start + 1, length - 1, length, length + 1, 2147483647}):
            if end > start:
                ranges.append(f"{start}-{end}")
    return ranges


def cases(programs: list[str], outfmts: list[str]) -> list[tuple[str, str, list[str]]]:
    out = []
    for program in programs:
        query, subject, qlen, slen = INPUTS[program]
        base = ["-query", query, "-subject", subject]
        for outfmt in outfmts:
            for option in ("-query_loc", "-subject_loc"):
                for spelling in SPELLINGS:
                    out.append((program, f"spelling {option} {spelling!r}",
                                base + [option, spelling, "-outfmt", outfmt]))
            for role, length in (("-query_loc", qlen), ("-subject_loc", slen)):
                if length is None:
                    continue
                for value in edge_ranges(length):
                    out.append((program, f"edge {role} {value}", base + [role, value, "-outfmt", outfmt]))
            for qv, sv in (("2-1", "2-1"), ("10-20-30", "5-1"), ("1-50", "2-1"), ("2-1", "1-50"),
                           ("100-600", "50-500"), ("1-100000", "1-100000")):
                out.append((program, f"pair {qv} {sv}",
                            base + ["-query_loc", qv, "-subject_loc", sv, "-outfmt", outfmt]))
    for program, query, subject in MULTI:
        if program not in programs:
            continue
        for outfmt in outfmts:
            for value in MULTI_RANGES:
                for role in ("-query_loc", "-subject_loc"):
                    out.append((program, f"multi {Path(query).name} {Path(subject).name} {role} {value}",
                                ["-query", query, "-subject", subject, role, value, "-outfmt", outfmt]))
    return out


def clean_env() -> dict[str, str]:
    return {key: value for key, value in os.environ.items()
            if key not in REPORT_ENV and not key.startswith(("LOSAT_", "RAYON_"))}


def run(command: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=ENGINE, capture_output=True, env=clean_env(), timeout=1200,
                          stdin=subprocess.DEVNULL)


def first_difference(a: bytes, b: bytes) -> str:
    for number, (left, right) in enumerate(zip(a.split(b"\n"), b.split(b"\n")), 1):
        if left != right:
            return f"line {number}: {left[:80]!r} / {right[:80]!r}"
    return f"lengths {len(a)} / {len(b)}"


def classify(bin_dir: Path, losat: str, case) -> str:
    program, label, argv = case
    ncbi = run([str(bin_dir / program), *argv])
    ours = run([losat, program, *argv])
    if (ncbi.returncode, ncbi.stdout, ncbi.stderr) == (ours.returncode, ours.stdout, ours.stderr):
        kind = "same" if ncbi.returncode == 0 else "same-error"
        detail = ncbi.stderr.decode(errors="replace").strip().replace("\n", " | ")[:120]
    elif ours.returncode != 0 and any(marker in ours.stderr for marker in REJECTION_MARKERS):
        kind = "losat-rejects"
        detail = (f"NCBI exit {ncbi.returncode} {ncbi.stderr.decode(errors='replace').strip()[:80]!r}; "
                  f"LOSAT {ours.stderr.decode(errors='replace').strip()[:120]!r}").replace("\n", " | ")
    elif ncbi.returncode == 1 and b"USAGE" in ncbi.stdout + ncbi.stderr and ours.returncode == 2:
        kind = "arg-error"
        detail = ours.stderr.decode(errors="replace").strip().replace("\n", " | ")[:120]
    else:
        kind = "DIFF"
        parts = []
        if ncbi.returncode != ours.returncode:
            parts.append(f"exit {ncbi.returncode}/{ours.returncode}")
        if ncbi.stdout != ours.stdout:
            parts.append("stdout " + first_difference(ncbi.stdout, ours.stdout))
        if ncbi.stderr != ours.stderr:
            parts.append("stderr " + first_difference(ncbi.stderr, ours.stderr))
        detail = "; ".join(parts)
    return "\t".join([program, kind, label, " ".join(repr(a) if (" " in a or a == "") else a for a in argv), detail])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", required=True)
    parser.add_argument("--losat", required=True)
    parser.add_argument("--programs", default="blastn,blastp,tblastn,tblastx")
    parser.add_argument("--outfmt", default="0,6,7")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--out")
    args = parser.parse_args()
    bin_dir = Path(args.bin_dir).resolve()
    losat = str(Path(args.losat).resolve())
    all_cases = cases(args.programs.split(","), args.outfmt.split(","))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        lines = list(pool.map(lambda case: classify(bin_dir, losat, case), all_cases))
    header = "program\tclass\tcase\targv\tdetail"
    text = header + "\n" + "\n".join(lines) + "\n"
    if args.out:
        Path(args.out).write_text(text)
    else:
        sys.stdout.write(text)
    counts: dict[tuple[str, str], int] = {}
    for line in lines:
        program, kind = line.split("\t")[:2]
        counts[(program, kind)] = counts.get((program, kind), 0) + 1
    for (program, kind), count in sorted(counts.items()):
        print(f"# {program}\t{kind}\t{count}", file=sys.stderr)
    return 1 if any(line.split("\t")[1] == "DIFF" for line in lines) else 0


if __name__ == "__main__":
    sys.exit(main())
