#!/usr/bin/env python3
"""Write the inputs of the Session S08+ (E2e) fixtures under LOSAT/tests/fasta/outfmt0/.

Deterministic: the files are cut from committed fixtures, so a rerun writes the same bytes.

- punct_hits_subject.fna: the first three subjects of the TBLASTX `many` fixture with the
  deflines `, ,`, `x1 ok` and `;~ ;` (NCBI's x_CleanAndCompress reads past the end of the
  first and the last title and crashes in outfmt 0; approved exception 2 of
  PD-LOSAT-NCBI-DEFECTS), as docs/evidence/losat_web_e2b/punct_defline.py builds it.
- punct_hits_standin.fna: the same records with placeholders of the lengths of the titles
  that NCBI's loop gives when it stops at the end of the string (`, ` and `;~;`): the
  stand-in subject of the `approved_punct_title` rows of LOSAT/tests/outfmt0_manifest.tsv
  (docs/evidence/losat_web_e2a/run_oracle.py). Its path has the length of the subject's.
- punct_hits_query.faa: the frame +1 peptide of the `many` query (the TBLASTN query).

Usage: make_inputs.py
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "LOSAT/tests/fasta/outfmt0"
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
DEFLINES = (", ,", "x1 ok", ";~ ;")
STAND_INS = ("OO", "x1 ok", "UUU")


def main() -> None:
    records = [chunk.split("\n", 1)[1] for chunk in (OUT / "tblastx_many_subject.fasta").read_text().split(">")[1:4]]
    for name, deflines in (("punct_hits_subject.fna", DEFLINES), ("punct_hits_standin.fna", STAND_INS)):
        (OUT / name).write_text("".join(f">{defline}\n{seq}" for defline, seq in zip(deflines, records)))
    seq = "".join((OUT / "tblastx_many_query.fasta").read_text().split("\n")[1:])
    index = {base: i for i, base in enumerate("TCAG")}
    peptide = "".join(CODE[index[seq[i]] * 16 + index[seq[i + 1]] * 4 + index[seq[i + 2]]]
                      for i in range(30, len(seq) - 32, 3)).replace("*", "")
    (OUT / "punct_hits_query.faa").write_text(f">pep frame +1 of the many query\n{peptide}\n")


if __name__ == "__main__":
    main()
