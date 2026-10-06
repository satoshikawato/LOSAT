#!/usr/bin/env python3
"""Sweep the non-default options of BLASTP, TBLASTN and TBLASTX against NCBI (Session S08+; comparison only).

Every case runs NCBI BLAST+ and LOSAT on the same inputs (from LOSAT/, standard input
empty) with one set of options, in each requested output format, and is classified:

- same: both succeed with the same stdout and stderr;
- same-error: both fail with the same exit status, stdout and stderr (NCBI's rejection);
- losat-rejects: LOSAT fails with a message that names what LOSAT does not support (an
  explicit rejection: "not supported by LOSAT" or "which LOSAT does not reproduce");
- arg-error: an argument that NCBI's argument parser rejects (USAGE and a CArgException,
  exit 1) and LOSAT's parser rejects (exit 2, "error: "), approved exception 1 of
  PD-LOSAT-CLI-NONSEARCH-DIFFERENCES;
- exception-gencode: a TBLASTN or TBLASTX search with a non-default `-db_gencode`
  (AGENTS.md's approved exceptions) whose output equals NCBI's search of the subject as a
  BLAST database (`makeblastdb`, then `-db`; NCBI applies `-db_gencode` there) outside
  the database lines of the report, or, for TBLASTN with `--api`, NCBI's local-subject
  C++ API oracle at that code (gencode_api_check.py: the stdout of TLOSAN Stage E's
  oracle at the CLI's query batches, calibrated for outfmt 0, and the stderr of NCBI's
  CLI); for TBLASTN's ID 32, which NCBI's CLI rejects, LOSAT succeeds
  (`PD-TLOSAN-LOCAL-GENCODE-32`, checked by its own C++ API oracle);
- DIFF: anything else, with the first difference.

The options of each program are grouped in sets (`--sets`, default all): matrices x gap
costs (every pair of NCBI's table for the matrix, read from blast_stat.c, and pairs
outside it), word size x threshold, window, composition-based statistics (with
`-ungapped` and `-use_sw_tback`), SEG, e-values, integer spellings, tasks, genetic codes,
TBLASTN's intron length, X-drops, sum statistics and masking, TBLASTX's culling limit,
and the spellings of `-outfmt` (which carry their own format). Exits 1 when any case is
DIFF.

Usage: option_sweep.py --program blastp|tblastn|tblastx --bin-dir DIR --losat LOSAT
       [--ncbi-src C++DIR] [--outfmt 0,6,7] [--sets a,b] [--jobs N] [--work DIR] [--api ORACLE]
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2a"))
from check_losat import normalize_database_lines  # noqa: E402
from run_oracle import ENGINE  # noqa: E402

F = "tests/fasta/outfmt0"
INPUTS = {
    "blastp": ["-query", f"{F}/e2e_protein_query.faa", "-subject", f"{F}/e2e_protein_subject.faa"],
    "tblastn": ["-query", f"{F}/e2e_protein_query.faa", "-subject", f"{F}/e2e_tblastn_subject.fna"],
    "tblastx": ["-query", f"{F}/tblastx_multi_query.fasta", "-subject", f"{F}/tblastx_multi_subject.fasta"],
}
REJECTION_MARKERS = (b"not supported by LOSAT", b"which LOSAT does not reproduce")
MATRICES = ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250", "IDENTITY"]
TABLE_NAMES = {"BLOSUM45": "blosum45", "BLOSUM50": "blosum50", "BLOSUM62": "blosum62", "BLOSUM80": "blosum80",
               "BLOSUM90": "blosum90", "PAM30": "pam30", "PAM70": "pam70", "PAM250": "pam250",
               "IDENTITY": "prot_idenity"}
EXTRA_GAPS = [(1, 1), (0, 0), (11, 0), (0, 1), (32767, 32767), (100, 10), (-1, 1), (11, -1)]
GENCODES = [str(code) for code in range(0, 35)] + ["0x2", "02", "+1", "1.0", "x"]


def gap_tables(ncbi_src: Path) -> dict[str, list[tuple[int, int]]]:
    """The gap costs of NCBI's table for each matrix (blast_stat.c `<matrix>_values`); the
    INT2_MAX row (the ungapped values) is the pair 32767/32767."""
    text = (ncbi_src / "src/algo/blast/core/blast_stat.c").read_text()
    tables = {}
    for matrix, name in TABLE_NAMES.items():
        body = re.search(rf"static array_of_8 {name}_values\[[A-Z0-9_]+\] = \{{(.*?)\}};", text, re.S)
        if body is None:
            raise SystemExit(f"no table {name}_values in blast_stat.c")
        pairs = []
        for row in re.findall(r"\{([^{}]*)\}", body.group(1)):
            fields = [field.strip() for field in row.split(",")]
            value = lambda field: 32767 if "INT2_MAX" in field else int(float(field))  # noqa: E731
            pairs.append((value(fields[0]), value(fields[1])))
        tables[matrix] = pairs
    return tables


def cases(program: str, tables: dict[str, list[tuple[int, int]]]) -> dict[str, list[list[str]]]:
    """The option sets of the program: set name -> list of argument lists."""
    protein = program in ("blastp", "tblastn")
    sets: dict[str, list[list[str]]] = {}
    if protein:
        rows = []
        for matrix in MATRICES:
            rows.append(["-matrix", matrix])
            for gap_open, gap_extend in [*tables[matrix], *EXTRA_GAPS]:
                rows.append(["-matrix", matrix, "-gapopen", str(gap_open), "-gapextend", str(gap_extend)])
        rows += [["-matrix", name] for name in ("blosum62", "Blosum45", "pam30", "BLOSUM62_20", "FOO", "")]
        rows += [["-gapopen", "10"], ["-gapextend", "2"], ["-gapopen", "0x9", "-gapextend", "0x2"],
                 ["-gapopen", "010", "-gapextend", "1"], ["-gapopen", "+11", "-gapextend", "+1"]]
        sets["matrix"] = rows
        sets["cbs"] = [["-comp_based_stats", value] for value in
                       ("0", "1", "2", "3", "4", "D", "d", "F", "f", "T", "t", "2x", "", "yes", "0x1", "-1")]
        sets["cbs"] += [["-comp_based_stats", value, "-ungapped"] for value in ("0", "F", "1", "2", "3")]
        sets["cbs"] += [["-ungapped"]]
        sets["cbs"] += [["-matrix", matrix, "-comp_based_stats", value] for matrix in MATRICES
                        for value in ("0", "1", "3")]
        sets["cbs"] += [["-matrix", matrix, "-ungapped", "-comp_based_stats", "0"] for matrix in MATRICES]
        if program == "blastp":
            sets["cbs"] += [["-use_sw_tback", "-comp_based_stats", value] for value in ("0", "1", "2", "3")]
            sets["cbs"] += [["-use_sw_tback"], ["-use_sw_tback", "-matrix", "PAM30"]]
    words = {"blastp": ("1", "2", "3", "4", "5", "6", "7", "8"), "tblastn": ("1", "2", "3", "4", "5", "6", "7", "8"),
             "tblastx": ("1", "2", "3", "4", "5", "6")}[program]
    thresholds = ("0", "1", "11", "11.5", "13", "30", "999", "-1", "1e1", "0x10", "nan", "+inf")
    sets["word"] = [["-word_size", word] for word in (*words, "0x3", "03", "+3", "3.0")]
    sets["word"] += [["-threshold", value] for value in thresholds]
    sets["word"] += [["-word_size", word, "-threshold", value] for word in words[1:6]
                     for value in ("1", "11", "11.5", "21")]
    if protein:
        sets["word"] += [["-matrix", matrix, "-word_size", word] for matrix in ("PAM30", "BLOSUM45", "IDENTITY")
                         for word in ("2", "5", "6")]
    sets["window"] = [["-window_size", value] for value in ("0", "1", "2", "10", "40", "100", "-1", "0x28", "40.0")]
    sets["window"] += [["-window_size", "0", "-threshold", "15"], ["-window_size", "0", "-word_size", "2"]]
    sets["seg"] = [["-seg", value] for value in
                   ("yes", "no", "12 2.2 2.5", "10 1.8 2.1", "0 0 0", "-1 -1 -1", "12 2.2", "12 2.2 2.5 1",
                    "x y z", "12 2.5 2.2", "12 nan 2.5", "12 inf 2.5", "1e1 2.2 2.5", "12  2.2 2.5",
                    "0x10 2.2 2.5", "12 2.2 2.5 ", "YES", "No", "")]
    sets["evalue"] = [["-evalue", value] for value in
                      ("1e-5", "0.001", "1000", "1e3", "+inf", "inf", "nan", "-nan", "0", "-1", "0x10", "1e999",
                       ".5", "5.", "1E-5", " 1", "1e-300", "100000", "1e-50")]
    sets["ints"] = [["-max_target_seqs", value] for value in
                    ("1", "4", "5", "0x10", "010", "+5", "-1", "0", "2147483647", "1e1", "5.0", "0x")]
    sets["ints"] += [["-num_threads", value] for value in ("0", "0x1", "+1", "-1")]
    if program != "tblastx":
        sets["ints"] += [["-max_hsps", value] for value in ("1", "2", "0", "0x2", "-1", "+1")]
    else:
        sets["ints"] += [["-max_hsps", "1"]]
    sets["outfmt"] = [["-outfmt", value] for value in
                      ("6", "06", "+6", " 6", "6 ", "6 std", "7 std qseqid", "0x6", "6.0", "5", "13", "18",
                       "abc", "", " ", "0 std", "10", "17", "-1", "7 qseqid")]
    if program == "blastp":
        sets["task"] = [["-task", task] for task in ("blastp", "blastp-fast", "blastp-short", "BLASTP", "deltablast",
                                                    "blastp-fast ")]
        sets["task"] += [["-task", "blastp-short", "-word_size", value] for value in ("2", "3", "4")]
        sets["task"] += [["-task", "blastp-fast", "-word_size", value] for value in ("3", "5", "6")]
        sets["task"] += [["-task", "blastp-fast", "-matrix", "PAM30"], ["-task", "blastp-short", "-matrix", "BLOSUM62"]]
    if program == "tblastn":
        sets["task"] = [["-task", task] for task in ("tblastn", "tblastn-fast", "TBLASTN", "tblastx")]
        sets["task"] += [["-task", "tblastn-fast", "-word_size", value] for value in ("3", "5", "6")]
        sets["gencode"] = [["-db_gencode", code] for code in GENCODES]
        sets["intron"] = [["-max_intron_length", value] for value in ("0", "1", "50", "100", "1000", "-1", "0x10")]
        sets["intron"] += [["-max_intron_length", "100", "-sum_stats", "false"]]
        sets["xdrop"] = [[option, value] for option in ("-xdrop_gap", "-xdrop_gap_final", "-xdrop_ungap")
                         for value in ("0", "7", "15", "30", "1e9", "-1", "nan", "0x10")]
        sets["sumstats"] = [["-sum_stats", value] for value in ("true", "false", "T", "F", "yes", "no", "1", "0", "x")]
        sets["sumstats"] += [["-sum_stats", value, "-ungapped", "-comp_based_stats", "0"] for value in ("true", "false")]
        sets["mask"] = [["-lcase_masking"], ["-soft_masking", "true"], ["-soft_masking", "false"],
                        ["-soft_masking", "F"], ["-soft_masking", "x"], ["-lcase_masking", "-seg", "no"],
                        ["-lcase_masking", "-soft_masking", "false"]]
    if program == "tblastx":
        sets["gencode"] = [[option, code] for option in ("-query_gencode", "-db_gencode") for code in GENCODES]
        sets["culling"] = [["-culling_limit", value] for value in ("0", "1", "2", "5", "-1", "0x1", "+1")]
        sets["other"] = [[option] + values for option, values in
                         (("-sum_stats", ["false"]), ("-num_descriptions", ["5"]), ("-num_alignments", ["5"]),
                          ("-line_length", ["60"]), ("-sorthits", ["1"]), ("-sorthsps", ["1"]),
                          ("-lcase_masking", []), ("-strand", ["plus"]), ("-qcov_hsp_perc", ["50"]),
                          ("-xdrop_ungap", ["10"]), ("-searchsp", ["1000000"]), ("-dbsize", ["1000000"]))]
    return sets


def first_difference(left: bytes, right: bytes) -> str:
    for number, (a, b) in enumerate(zip(left.split(b"\n"), right.split(b"\n")), 1):
        if a != b:
            return f"line {number}: {a[:80]!r} / {b[:80]!r}"
    return f"{len(left)} / {len(right)} bytes"


def db_gencode(argv: list[str]) -> str | None:
    """The `-db_gencode` value of an argument list, or None."""
    return argv[argv.index("-db_gencode") + 1] if "-db_gencode" in argv else None


class Sweep:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.dbs: dict[str, str] = {}

    def database(self, subject: str) -> str:
        """The subject as a BLAST database (made once; NCBI's oracle of a non-default
        `-db_gencode`)."""
        if subject not in self.dbs:
            directory = self.args.work / "db" / Path(subject).stem
            shutil.rmtree(directory, ignore_errors=True)
            directory.mkdir(parents=True)
            subprocess.run([str(self.args.bin_dir / "makeblastdb"), "-in", subject, "-dbtype", "nucl",
                            "-out", str(directory / "subject")], cwd=ENGINE, capture_output=True, check=True)
            self.dbs[subject] = str(directory / "subject")
        return self.dbs[subject]

    def run(self, case: tuple[str, list[str], str]) -> list[str]:
        name, options, outfmt = case
        program = self.args.program
        argv = [*INPUTS[program], *options, *([] if "-outfmt" in options else ["-outfmt", outfmt])]
        try:
            ncbi = subprocess.run([str(self.args.bin_dir / program), *argv], cwd=ENGINE, capture_output=True,
                                  timeout=self.args.timeout, stdin=subprocess.DEVNULL)
            ours = subprocess.run([str(self.args.losat.resolve()), program, *argv], cwd=ENGINE, capture_output=True,
                                  timeout=self.args.timeout, stdin=subprocess.DEVNULL)
        except subprocess.TimeoutExpired:
            return [name, outfmt, " ".join(options), "timeout", "", ""]
        result = self.classify(program, argv, ncbi, ours)
        message = (ncbi.stderr if ncbi.returncode else ours.stderr).decode(errors="replace").strip().splitlines()
        message = next((line for line in message if "rror" in line), message[0] if message else "")
        lmessage = ours.stderr.decode(errors="replace").strip().splitlines()
        lmessage = next((line for line in lmessage if "rror" in line), lmessage[0] if lmessage else "")
        return [name, outfmt, " ".join(repr(o) if (" " in o or not o) else o for o in options), result,
                message[:160], lmessage[:160]]

    def classify(self, program: str, argv: list[str], ncbi: subprocess.CompletedProcess,
                 ours: subprocess.CompletedProcess) -> str:
        same_streams = ncbi.stdout == ours.stdout and ncbi.stderr == ours.stderr
        if ncbi.returncode == 0 and ours.returncode == 0 and same_streams:
            return "same"
        if ncbi.returncode and ncbi.returncode == ours.returncode and same_streams:
            return "same-error"
        if ours.returncode and any(marker in ours.stderr for marker in REJECTION_MARKERS):
            return "losat-rejects"
        if (ncbi.returncode == 1 and b"CArgException" in ncbi.stderr and ours.returncode == 2
                and ours.stderr.startswith(b"error: ") and ncbi.stdout == ours.stdout):
            return "arg-error"
        code = db_gencode(argv)
        if program in ("tblastn", "tblastx") and code not in (None, "1") and ours.returncode == 0:
            if program == "tblastn" and code == "32" and ncbi.returncode == 1:
                return "exception-gencode"
            index = argv.index("-subject")
            db_argv = [*argv[:index], "-db", self.database(argv[index + 1]), *argv[index + 2:]]
            oracle = subprocess.run([str(self.args.bin_dir / program), *db_argv], cwd=ENGINE, capture_output=True,
                                    stdin=subprocess.DEVNULL)
            if (oracle.returncode == 0 and oracle.stderr == ours.stderr
                    and normalize_database_lines(oracle.stdout) == normalize_database_lines(ours.stdout)):
                return "exception-gencode"
            if program == "tblastn" and self.args.api and ncbi.returncode == 0 and ncbi.stderr == ours.stderr:
                from gencode_api_check import api_output, ncbi_integer  # noqa: PLC0415
                fmt = argv[argv.index("-outfmt") + 1]
                if fmt in ("0", "6", "7") and api_output(self.args.api, argv, ncbi_integer(code), fmt) == ours.stdout:
                    return "exception-gencode"
        if ncbi.returncode == 0 and ours.returncode == 0:
            return "DIFF " + ("stdout " + first_difference(ncbi.stdout, ours.stdout) if ncbi.stdout != ours.stdout
                              else "stderr " + first_difference(ncbi.stderr, ours.stderr))
        return (f"DIFF exit {ncbi.returncode}/{ours.returncode}: "
                + ("stdout " + first_difference(ncbi.stdout, ours.stdout) if ncbi.stdout != ours.stdout
                   else "stderr " + first_difference(ncbi.stderr, ours.stderr)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--program", choices=sorted(INPUTS), required=True)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--ncbi-src", type=Path, default=Path("/mnt/c/Users/genom/GitHub/ncbi-blast/c++"))
    parser.add_argument("--outfmt", default="0,6,7")
    parser.add_argument("--sets", default="")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--work", type=Path, default=Path("/tmp/losat_e2e_option_sweep"))
    parser.add_argument("--api", type=Path, help="TLOSAN Stage E's local-subject API oracle (TBLASTN gencodes)")
    parser.add_argument("--list", action="store_true", help="print the cases and exit")
    args = parser.parse_args()
    args.work = args.work.resolve()
    sets = cases(args.program, gap_tables(args.ncbi_src))
    chosen = args.sets.split(",") if args.sets else list(sets)
    work = []
    for name in chosen:
        for options in sets[name]:
            formats = ["-"] if name == "outfmt" else args.outfmt.split(",")
            work += [(name, options, outfmt) for outfmt in formats]
    if args.list:
        for name, options, outfmt in work:
            print(f"{name}\t{outfmt}\t{options}")
        print(f"# cases={len(work)}")
        return 0
    sweep = Sweep(args)
    # The databases are made before the parallel runs.
    if args.program in ("tblastn", "tblastx") and any(db_gencode(options) not in (None, "1") for _, options, _ in work):
        sweep.database(INPUTS[args.program][3])
    with ThreadPoolExecutor(args.jobs) as pool:
        lines = list(pool.map(sweep.run, work))
    print("set\toutfmt\toptions\tresult\tncbi_message\tlosat_message")
    for line in lines:
        print("\t".join(line))
    counts: dict[str, int] = {}
    for line in lines:
        key = line[3].split()[0]
        counts[key] = counts.get(key, 0) + 1
    print(f"# program={args.program} cases={len(lines)} " + " ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    return 1 if counts.get("DIFF") or counts.get("timeout") else 0


if __name__ == "__main__":
    sys.exit(main())
