#!/usr/bin/env python3
"""Build the fixtures of masks whose order NCBI takes from memory addresses (LOSAT/tests/seqdb/mask_order).

When two Seq-ids of one sequence with the same packed number (`gi|5` and `gb|AY000005.1|`) have
masks on one strand, makeblastdb writes their ranges in the order of the memory addresses of their
handles, which LOSAT does not follow: it refuses such a build. That order shows only in the column
files of a build that succeeds. These cases are runs in which an error of makeblastdb that does not
depend on the order ends the build (the bounds check of the same sequence, also over a later
algorithm; a later sequence; a later input file), where LOSAT gives makeblastdb's outcome, and one
run in which nothing else fails (`tie`, LOSAT's refusal). They are the probes t1, t1b and t3 of the
second review of M4a-1 (2026-11-14).

NCBI BLAST+ 2.17.0 is the comparison oracle only. Each case runs, in its own directory of the work
directory, with BLAST_USAGE_REPORT=false:

  makeblastdb -in <files> -dbtype nucl -parse_seqids -title in.fna -mask_data <files> -out db
              -blastdb_version <4|5>

usage: make_mask_order_fixtures.py --ncbi-bin DIR --work DIR [--out LOSAT/tests/seqdb/mask_order]

Output files:
  in.fna, nodata.fna, *.asnt          the inputs
  cases.tsv                           one row per case: name, -in files (space separated), -mask_data
                                      files (comma separated)
  expected/<case>.v<4|5>.rc           makeblastdb's exit status
  expected/<case>.v<4|5>.stderr       its standard error
  expected/<case>.v<4|5>.stdout       its log, with the build time as <T>, the database path as <OUT>
                                      and the seconds of an "Adding sequences" line as <S>
  expected/<case>.v<4|5>.files        the files `db.*` left in the directory, one name per line
"""
import argparse
import os
import re
import shutil
import subprocess
from pathlib import Path

# The sequences of the review probes: `gi|5|gb|AY000005.1|` (200 nt; its two Seq-ids have the
# packed number 5) and `lcl|other` (200 nt).
IN_FNA = (
    ">gi|5|gb|AY000005.1| one\n"
    "AAAGCGGCACTTGTGAAGTGTTCCCCACGCCGCTTGGGTCTTCTGTGTTGTTCGCGTGGTGCTGAGACAAAGCACGCCATAAGGCCAAAAAAA"
    "GGCCCATACCAAGAGGTAGTAGTCTCAGAATCTTGCGGGTACAGACCCATCACCTAGACGGTGACATTCAACAAACCACATTGTCCTTAATC"
    "ATGAAGGGGATAAGC\n"
    ">lcl|other two\n"
    "ATATTTCAAGAGGACTCAGTTCGTAGAAAGTCAATATGGTCGGTTTTGTCCTGTAAAGCCTAAACGTCGTCGACTAGCGCCTCTGCTTATCT"
    "ATGTGTTGGACCTTAGTTCAATCTCATCGCTCATTGCTCAGATATGTGTAAGCTGCACTTTGCAGTAGATTCGTCTGAGGGGGTACTCAGAC"
    "TCGAAATGCGGAGTGC\n"
)
# A second input file whose only record has no sequence data: makeblastdb warns that it does not
# look like FASTA, then fails it with "No sequences added" after the first file was added.
NODATA_FNA = ">lcl|x\n\n"

GI = "gi 5"
GB = 'genbank { accession "AY000005", version 1 }'


def interval(start, stop, seq_id, strand="plus"):
    st = f"strand {strand}, " if strand else ""
    return f"int {{ from {start}, to {stop}, {st}id {seq_id} }}"


def mask_object(program, locs):
    body = ",\n      ".join(locs)
    return (f"Blast-db-mask-info ::= {{\n  algo-id 0,\n  algo-program {program},\n"
            f"  algo-options \"\",\n  masks {{\n    masks {{\n      {body}\n    }},\n"
            f"    more FALSE\n  }}\n}}\n")


MASKS = {
    # t1: the tied pair, one range past the end of the 200-nt sequence.
    "oob_same.asnt": mask_object("windowmasker", [interval(10, 20, GI), interval(300, 310, GB)]),
    # t1b: the tied pair in bounds; a second algorithm past the end.
    "tie_in.asnt": mask_object("windowmasker", [interval(10, 20, GI), interval(30, 40, GB)]),
    "oob_dust.asnt": mask_object("dust", [interval(150, 900, GI, None)]),
    # t3: the tied pair in bounds; the next sequence past its end.
    "oob_next.asnt": mask_object("windowmasker", [
        interval(10, 40, GI), interval(50, 60, GB), interval(150, 900, 'local str "other"', None)]),
    # The tied pair alone (NCBI builds, in the order of the addresses).
    "tie.asnt": mask_object("windowmasker", [interval(10, 40, GI), interval(50, 60, GB)]),
}

CASES = [
    # name, -in files, -mask_data files
    ("oob_same_sequence", "in.fna", "oob_same.asnt"),
    ("oob_later_algorithm", "in.fna", "tie_in.asnt,oob_dust.asnt"),
    ("oob_later_sequence", "in.fna", "oob_next.asnt"),
    ("no_sequences_later", "in.fna nodata.fna", "tie.asnt"),
    ("tie", "in.fna", "tie.asnt"),
]


def masked_log(text):
    """makeblastdb's log with what changes from run to run replaced (see the module doc)."""
    lines = []
    for line in text.split("\n"):
        if line.startswith("Building a new DB, current time: "):
            line = "Building a new DB, current time: <T>"
        elif line.startswith("New DB name:   "):
            line = "New DB name:   <OUT>"
        else:
            line = re.sub(r"( sequences in )\S+( seconds\.)$", r"\1<S>\2", line)
        lines.append(line)
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "mask_order"))
    args = ap.parse_args()
    bin_dir = Path(args.ncbi_bin)
    work = Path(args.work)
    out = Path(args.out)
    if work.exists():
        raise SystemExit(f"{work} exists; give a new work directory")
    work.mkdir(parents=True)
    (out / "expected").mkdir(parents=True, exist_ok=True)

    inputs = {"in.fna": IN_FNA, "nodata.fna": NODATA_FNA, **MASKS}
    for name, text in inputs.items():
        (out / name).write_text(text)
    env = dict(os.environ, BLAST_USAGE_REPORT="false")
    rows = []
    for name, in_files, mask_files in CASES:
        rows.append(f"{name}\t{in_files}\t{mask_files}\n")
        for version in (5, 4):
            run_dir = work / f"{name}.v{version}"
            run_dir.mkdir()
            for file in inputs:
                shutil.copyfile(out / file, run_dir / file)
            cmd = [bin_dir / "makeblastdb", "-in", in_files, "-dbtype", "nucl", "-parse_seqids",
                   "-title", "in.fna", "-mask_data", mask_files, "-out", "db",
                   "-blastdb_version", str(version)]
            result = subprocess.run([str(c) for c in cmd], cwd=run_dir, capture_output=True, env=env)
            stem = out / "expected" / f"{name}.v{version}"
            Path(f"{stem}.rc").write_text(f"{result.returncode}\n")
            Path(f"{stem}.stderr").write_bytes(result.stderr)
            Path(f"{stem}.stdout").write_text(masked_log(result.stdout.decode()))
            left = sorted(p.name for p in run_dir.iterdir() if p.name.startswith("db."))
            Path(f"{stem}.files").write_text("".join(f"{f}\n" for f in left))
    (out / "cases.tsv").write_text("".join(rows))


if __name__ == "__main__":
    main()
