#!/usr/bin/env python3
"""Compare BLASTP outfmt 0 on protein subject titles with NCBI (Session S08+; comparison only).

NCBI writes the title of a protein subject with `CDeflineGenerator::GenerateDefline`
(`x_CleanAndCompress` with the protein rule: `. [` and `, [` become ` [`), in the
description table (with the TPA-like prefixes) and in the alignment heading. The subject
is the homolog of the first E2e query (`e2e_protein_subject.faa`), searched with that query,
under each defline:

- the 1023 deflines of one to five of `,;~` and space that do not start with a space (the
  deflines of docs/evidence/losat_web_e2e/title_sweep.py): NCBI's loop runs past the end of
  some of them and NCBI crashes; approved exception 2 of PD-LOSAT-NCBI-DEFECTS covers only
  the nucleotide subjects of BLASTN, TBLASTX and TBLASTN, so LOSAT's BLASTP rejects them;
- the titles of the inventory (range RP, rows 18-20: trailing punctuation, spaces, brackets,
  prefixes, HTML character references, long titles).

Each is "same" (stdout, stderr and exit status), or "losat-rejects": NCBI dies of a signal
or decodes an HTML character reference, and LOSAT exits non-zero naming what it does not
support. Any other outcome is printed; exits 1 when any is.

Usage: protein_title_sweep.py --bin-dir DIR --losat LOSAT --work DIR [--jobs N]
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

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
REJECTION = b"not supported by LOSAT's BLASTP"


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
            result = "same"
        elif ours.returncode and REJECTION in ours.stderr and (ncbi.returncode < 0 or b"&" in defline.encode()):
            result = "losat-rejects"
        else:
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
