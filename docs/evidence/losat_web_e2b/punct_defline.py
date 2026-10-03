#!/usr/bin/env python3
"""TBLASTX and TBLASTN on subject titles that NCBI's x_CleanAndCompress reads past (comparison only).

Approved exception 2 of PD-LOSAT-NCBI-DEFECTS (outfmt 0 titles made only of punctuation)
covers BLASTN, and the maintainer extended it to TBLASTX and TBLASTN in session S08b
(DW-17); since session S08+ LOSAT writes such a title stopped at the end of the string, as
BLASTN does (docs/evidence/losat_web_e2e/title_sweep.py checks 1023 deflines). The inputs:
three subjects of the TBLASTX `many` fixture with the deflines `, ,`, `x1 ok` and `;~ ;`
(the first and last reach the defect), the `many` query for TBLASTX and its frame +1
peptide for TBLASTN. Expected: NCBI dies of a signal in outfmt 0; LOSAT exits 0 with
NCBI's stderr and the stdout of NCBI for the same subjects with placeholder titles of the
lengths of LOSAT's titles once the placeholders are replaced by the titles that NCBI's loop
gives when it stops at the end of the string; outfmt 6 and 7 are the same as NCBI's.

Usage: punct_defline.py --bin-dir DIR --losat LOSAT --work DIR
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2g"))
from title_sweep import title  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
MANY = REPO / "LOSAT/tests/fasta/outfmt0"
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    records = [chunk.split("\n", 1)[1] for chunk in (MANY / "tblastx_many_subject.fasta").read_text().split(">")[1:4]]
    deflines = (", ,", "x1 ok", ";~ ;")
    subject = args.work / "subject.fna"
    subject.write_text("".join(f">{defline}\n{seq}" for defline, seq in zip(deflines, records)))
    # The stand-in subject: a placeholder of the length of each title that NCBI reads past.
    placeholders = {", ,": "Q" * len(title(", ,")), ";~ ;": "Z" * len(title(";~ ;"))}
    stand_in = args.work / "subject_p.fna"
    stand_in.write_text("".join(f">{placeholders.get(defline, defline)}\n{seq}"
                                for defline, seq in zip(deflines, records)))
    query_nt = MANY / "tblastx_many_query.fasta"
    seq = "".join(query_nt.read_text().split("\n")[1:])
    index = {base: i for i, base in enumerate("TCAG")}
    peptide = "".join(CODE[index[seq[i]] * 16 + index[seq[i + 1]] * 4 + index[seq[i + 2]]]
                      for i in range(30, len(seq) - 32, 3)).replace("*", "")
    query_aa = args.work / "query.faa"
    query_aa.write_text(f">pep frame +1 of the many query\n{peptide}\n")
    failures = 0
    print("program\toutfmt\tncbi_exit\tlosat_exit\tresult")
    for program, query in (("tblastx", query_nt), ("tblastn", query_aa)):
        for outfmt in ("0", "6", "7"):
            argv = ["-query", str(query), "-subject", str(subject), "-outfmt", outfmt]
            ncbi = subprocess.run([str(args.bin_dir / program), *argv], capture_output=True)
            losat = subprocess.run([str(args.losat.resolve()), program, *argv], capture_output=True)
            if outfmt == "0":
                other = subprocess.run([str(args.bin_dir / program), "-query", str(query), "-subject",
                                        str(stand_in), "-outfmt", outfmt], capture_output=True)
                expected = other.stdout.replace(stand_in.name.encode(), subject.name.encode())
                for defline, placeholder in placeholders.items():
                    expected = expected.replace(placeholder.encode(), title(defline).encode())
                crashed = ncbi.returncode < 0 or ncbi.returncode >= 128
                same = (losat.returncode, losat.stdout, losat.stderr) == (0, expected, other.stderr)
                result = "ncbi-crash, exception-2" if crashed and other.returncode == 0 and same else "UNEXPECTED"
            else:
                same = (ncbi.stdout, ncbi.stderr, ncbi.returncode) == (losat.stdout, losat.stderr, losat.returncode)
                result = "same" if same else "DIFF"
            failures += result in ("UNEXPECTED", "DIFF")
            print(f"{program}\t{outfmt}\t{ncbi.returncode}\t{losat.returncode}\t{result}")
    print(f"# unexpected={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
