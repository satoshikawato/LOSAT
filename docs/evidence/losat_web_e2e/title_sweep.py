#!/usr/bin/env python3
"""Compare TBLASTX and TBLASTN outfmt 0 on subject deflines of punctuation with NCBI (Session S08+; comparison only).

docs/evidence/losat_web_e2g/title_sweep.py (BLASTN) for the programs to which the
maintainer extended approved exception 2 of PD-LOSAT-NCBI-DEFECTS in session S08b (DW-17):
LOSAT stops the cleanup of the title at the end of the string where NCBI's
`x_CleanAndCompress` runs past it (NCBI crashes in outfmt 0). Every defline of one to five
of `,;~` and space that does not start with a space (1023 deflines) is the title of a
subject with hits: the first subject of the TBLASTX `many` fixture, searched with the
`many` query (TBLASTX) or with the frame +1 peptide of that query (TBLASTN, as
docs/evidence/losat_web_e2b/punct_defline.py makes it). Each is "same" (stdout, stderr and
exit status) or "exception-2": NCBI dies of a signal, LOSAT exits 0 with NCBI's stderr,
and LOSAT's stdout equals NCBI's for the same subject with a placeholder title of the
length of LOSAT's title once the placeholder is replaced by the title that NCBI's loop
gives when it stops at the end of the string. Any other outcome is printed; exits 1 when
any is.

Usage: title_sweep.py --bin-dir DIR --losat LOSAT --work DIR [--jobs N] [--programs tblastx,tblastn]
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2g"))
from title_sweep import title  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
MANY = REPO / "LOSAT/tests/fasta/outfmt0"
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"


def inputs(work: Path) -> tuple[str, dict[str, Path]]:
    """The subject's residues and each program's query (copied into `work`)."""
    residues = (MANY / "tblastx_many_subject.fasta").read_text().split(">")[1].split("\n", 1)[1]
    query_text = (MANY / "tblastx_many_query.fasta").read_text()
    query_nt = work / "query.fna"
    query_nt.write_text(query_text)
    seq = "".join(query_text.split("\n")[1:])
    index = {base: i for i, base in enumerate("TCAG")}
    peptide = "".join(CODE[index[seq[i]] * 16 + index[seq[i + 1]] * 4 + index[seq[i + 2]]]
                      for i in range(30, len(seq) - 32, 3)).replace("*", "")
    query_aa = work / "query.faa"
    query_aa.write_text(f">pep frame +1 of the many query\n{peptide}\n")
    return residues, {"tblastx": query_nt, "tblastn": query_aa}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--programs", default="tblastx,tblastn")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    residues, queries = inputs(args.work)
    titles = ["".join(chars) for length in range(1, 6) for chars in itertools.product(",;~ ", repeat=length)
              if chars[0] != " "]
    programs = args.programs.split(",")
    cases = [(program, index, defline) for program in programs for index, defline in enumerate(titles)]

    def run(case: tuple[str, int, str]) -> tuple[str, str, str, tuple[int, int] | None]:
        program, index, defline = case
        subject = args.work / f"{program}_title_{index}.fa"
        subject.write_text(f">{defline}\n{residues}")
        argv = ["-query", str(queries[program]), "-subject", str(subject)]
        try:
            ncbi = subprocess.run([str(args.bin_dir / program), *argv], capture_output=True, timeout=120)
            ours = subprocess.run([str(args.losat.resolve()), program, *argv], capture_output=True, timeout=120)
        except subprocess.TimeoutExpired:
            return program, defline, "timeout", None
        codes = (ncbi.returncode, ours.returncode)
        if (ncbi.returncode, ncbi.stdout, ncbi.stderr) == (ours.returncode, ours.stdout, ours.stderr):
            return program, defline, "same", codes
        expected_title = title(defline)
        if ncbi.returncode < 0 and ours.returncode == 0 and expected_title:
            placeholder = "Z" * len(expected_title)
            stand_in = args.work / f"{program}_title_{index}_p.fa"
            stand_in.write_text(f">{placeholder}\n{residues}")
            other = subprocess.run([str(args.bin_dir / program), "-query", str(queries[program]),
                                    "-subject", str(stand_in)], capture_output=True, timeout=120)
            expected = (other.stdout.replace(placeholder.encode(), expected_title.encode())
                        .replace(stand_in.name.encode(), subject.name.encode()))
            if other.returncode == 0 and ours.stdout == expected and ours.stderr == other.stderr:
                return program, defline, "exception-2", codes
        return program, defline, "different", codes

    counts: dict[tuple[str, str], int] = {}
    with ThreadPoolExecutor(args.jobs) as pool:
        for program, defline, outcome, codes in pool.map(run, cases):
            counts[(program, outcome)] = counts.get((program, outcome), 0) + 1
            if outcome not in ("same", "exception-2"):
                print(f"{outcome.upper()}\t{program}\t{defline!r}\t{codes}")
    unexpected = 0
    for program in programs:
        same = counts.get((program, "same"), 0)
        exception = counts.get((program, "exception-2"), 0)
        unexpected += len(titles) - same - exception
        print(f"# {program} deflines={len(titles)} same={same} exception-2={exception} "
              f"unexpected={len(titles) - same - exception}")
    return 1 if unexpected else 0


if __name__ == "__main__":
    sys.exit(main())
