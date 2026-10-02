#!/usr/bin/env python3
"""TBLASTX and TBLASTN on subject titles that NCBI's x_CleanAndCompress reads past (comparison only).

Approved exception 2 of PD-LOSAT-NCBI-DEFECTS (outfmt 0 titles made only of punctuation)
covers BLASTN. For TBLASTX and TBLASTN, LOSAT rejects such subjects in outfmt 0 until the
maintainer decides. The inputs: three subjects of the TBLASTX `many` fixture with the
deflines `, ,`, `x1 ok` and `;~ ;` (the first and last reach the defect), the `many`
query for TBLASTX and its frame +1 peptide for TBLASTN. Expected: NCBI dies of a signal in
outfmt 0; LOSAT exits 1 with the rejection; outfmt 6 and 7 are the same as NCBI's.

Usage: punct_defline.py --bin-dir DIR --losat LOSAT --work DIR
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

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
    subject = args.work / "subject.fna"
    subject.write_text("".join(f">{defline}\n{seq}" for defline, seq in zip((", ,", "x1 ok", ";~ ;"), records)))
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
                rejected = losat.returncode == 1 and b"reads past its end" in losat.stderr and not losat.stdout
                ok = ncbi.returncode < 0 or ncbi.returncode >= 128
                result = "ncbi-crash, losat-rejects" if ok and rejected else "UNEXPECTED"
            else:
                same = (ncbi.stdout, ncbi.stderr, ncbi.returncode) == (losat.stdout, losat.stderr, losat.returncode)
                result = "same" if same else "DIFF"
            failures += result in ("UNEXPECTED", "DIFF")
            print(f"{program}\t{outfmt}\t{ncbi.returncode}\t{losat.returncode}\t{result}")
    print(f"# unexpected={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
