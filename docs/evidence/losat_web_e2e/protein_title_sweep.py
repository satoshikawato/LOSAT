#!/usr/bin/env python3
"""Compare BLASTP outfmt 0 on protein subject titles with NCBI (Session S08+; comparison only).

NCBI writes the title of a protein subject with `CDeflineGenerator::GenerateDefline`
(`x_CleanAndCompress` with the protein rule: `. [` and `, [` become ` [`), in the
description table (with the TPA-like prefixes) and in the alignment heading. The subject
is the homolog of the first E2e query (`e2e_protein_subject.faa`), searched with that query,
under each defline:

- the 1023 deflines of one to five of `,;~` and space that do not start with a space (the
  deflines of docs/evidence/losat_web_e2e/title_sweep.py): NCBI's loop runs past the end of
  some of them and NCBI crashes; approved exception 2 of PD-LOSAT-NCBI-DEFECTS, extended to
  BLASTP by the maintainer in session SFd (DW-27, implemented in SFe): LOSAT stops the
  cleanup at the end of the string (until SFe LOSAT's BLASTP rejected them);
- the titles of the inventory (range RP, rows 18-20: trailing punctuation, spaces, brackets,
  prefixes, HTML character references, long titles).

Each is "same" (stdout, stderr and exit status), or "exception-2": NCBI dies of a signal,
LOSAT exits 0 with NCBI's stderr, and LOSAT's stdout equals NCBI's for the same subject with
a placeholder title of the length of LOSAT's title once the placeholder is replaced by the
title that NCBI's loop gives when it stops at the end of the string (the check of
title_sweep.py; the stand-in's file name has the length of the subject's). Any other outcome
is printed; exits 1 when any is.

Usage: protein_title_sweep.py --bin-dir DIR --losat LOSAT --work DIR [--jobs N]
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2g"))
from title_sweep import title  # noqa: E402  (NCBI's title stopped at the end of the string)

REPO = Path(__file__).resolve().parents[3]
INPUTS = REPO / "LOSAT/tests/fasta/outfmt0"
WORD_TITLES = [
    "", " ", "hypothetical protein", "hypothetical protein.", "hypothetical protein,",
    "hypothetical protein;", "hypothetical protein~", "hypothetical protein . , ; ~ ",
    "hypothetical  protein   double  space", "hypothetical protein [Homo sapiens]",
    "hypothetical protein [Homo sapiens].", "protein. [Homo sapiens]", "protein, [Homo sapiens]",
    "[Homo sapiens]", "x [a] [b]", "partial", "protein, partial", "&amp; protein", "a &lt; b",
    "a < b > c & d \"q\" 'r'", "MULTISPECIES: hypothetical protein", "TPA: hyp", "UNVERIFIED: hyp",
    "MAG: hyp", "MAG hyp", "TPA_exp: hyp", "TSA: hyp", "LOW QUALITY PROTEIN: hyp", "PREDICTED: hyp",
    "hyp (fragment)", "hyp ( fragment )", "a ,b", "a ;b", "a,,b", "a, ,b", "( x", "x )",
    "uncharacterized protein", "Hypothetical protein", "ABC", "A", ".", ",", "...", "(", ")", "!!!",
    "?", "<>", "&", "&&", "&#65;", "&nbsp;x", "x" * 100, "word " * 30,
    ",".join("abcdefghijklmnopqrstuvwxyz" * 2), "-".join("abcdefghijklmnopqrstuvwxyz" * 2),
    "word " * 13 + "end", "w" * 60 + " " + "v" * 20, "w" * 70 + " tail",
]


def deflines() -> list[str]:
    punct = ["".join(chars) for n in range(1, 6) for chars in itertools.product(",;~ ", repeat=n)
             if chars[0] != " "]
    return [f"s1 {title}" if title else "s1" for title in WORD_TITLES] + punct


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    query = (INPUTS / "e2e_protein_query.faa").read_text().split(">")[1]
    (args.work / "query.faa").write_text(">" + query)
    residues = (INPUTS / "e2e_protein_subject.faa").read_text().split(">")[1].split("\n", 1)[1]

    def run(case: tuple[int, str]) -> str:
        index, defline = case
        subject = args.work / f"s{index}.faa"
        subject.write_text(f">{defline}\n{residues}")
        argv = ["-query", str(args.work / "query.faa"), "-subject", str(subject), "-outfmt", "0"]
        ncbi = subprocess.run([str(args.bin_dir / "blastp"), *argv], capture_output=True, stdin=subprocess.DEVNULL)
        ours = subprocess.run([str(args.losat), "blastp", *argv], capture_output=True, stdin=subprocess.DEVNULL)
        if (ncbi.returncode, ncbi.stdout, ncbi.stderr) == (ours.returncode, ours.stdout, ours.stderr):
            return f"{index}\t{defline!r}\tsame"
        expected_title = title(defline)
        if ncbi.returncode < 0 and ours.returncode == 0 and expected_title:
            placeholder = "Z" * len(expected_title)
            stand_in = args.work / f"p{index}.faa"
            stand_in.write_text(f">{placeholder}\n{residues}")
            other = subprocess.run([str(args.bin_dir / "blastp"), *argv[:3], str(stand_in), *argv[4:]],
                                   capture_output=True, stdin=subprocess.DEVNULL)
            expected = (other.stdout.replace(placeholder.encode(), expected_title.encode())
                        .replace(str(stand_in).encode(), str(subject).encode()))
            if other.returncode == 0 and ours.stdout == expected and ours.stderr == other.stderr:
                return f"{index}\t{defline!r}\texception-2"
        result = f"DIFF ncbi={ncbi.returncode} losat={ours.returncode}"
        return f"{index}\t{defline!r}\t{result}"

    with ThreadPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(run, enumerate(deflines())))
    counts: dict[str, int] = {}
    for row in rows:
        print(row)
        key = row.split("\t")[2].split(" ")[0]
        counts[key] = counts.get(key, 0) + 1
    print(f"# {counts}")
    return 1 if counts.get("DIFF") else 0


if __name__ == "__main__":
    sys.exit(main())
