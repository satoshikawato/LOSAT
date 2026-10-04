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
- e2e_protein_query.faa: seven MeenMJNV proteins with homologs in LvMJNV (identities from
  about 30 to 60 %, one weak pair) and the first AvCLPV protein (low complexity, for SEG);
  e2e_protein_subject.faa: the seven LvMJNV homologs and six other LvMJNV proteins (every
  17th record); the BLASTP inputs of option_sweep.py.
- e2e_tblastn_subject.fna: five windows of the LvMJNV genome that hold the TBLASTN hits of
  e2e_protein_query.faa, the second with 300 lowercase bases and the fourth with a run of
  40 N; the TBLASTN subject of option_sweep.py.
- e2e_many_query.faa, e2e_many_subject.faa, e2e_many_subject.fna: the first 200 residues
  of the first query and 300 variants of it (each residue replaced with probability 0.3,
  seeded), as proteins and as nucleotides (one codon per residue): more than 250 subjects
  with hits, for the numbers of descriptions (500) and alignments (250) of outfmt 0.
- e2e_titles_subject.faa: the first homolog of e2e_protein_subject.faa under protein titles
  that NCBI's GenerateDefline rewrites (trailing punctuation, double spaces, `. [`, `, [`,
  ` ,`, `( `, prefixes) and one ending in 60 letters (the CFastaReader title warning).
- e2e_amb_subject.fna: e2e_tblastn_subject.fna with 2 % of the bases replaced by IUPAC
  ambiguity codes (seeded): the displayed Sbjct rows of TBLASTN (B, Z, J, X).

Usage: make_inputs.py
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "LOSAT/tests/fasta/outfmt0"
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
DEFLINES = (", ,", "x1 ok", ";~ ;")
STAND_INS = ("OO", "x1 ok", "UUU")
FASTA = REPO / "LOSAT/tests/fasta"
QUERIES = ("BDT62569.1", "BDT62529.1", "BDT62562.1", "BDT62620.1", "BDT62567.1", "BDT62581.1", "BDT62533.1")
HOMOLOGS = ("BDT62125.1", "BDT62158.1", "BDT62172.1", "BDT62143.1", "BDT62170.1", "BDT62187.1", "BDT62149.1")
WINDOWS = ((93000, 98000), (143000, 152000), (160000, 165000), (199000, 206000), (242000, 247000))


def protein_records(path: Path) -> dict[str, str]:
    """The records of a FASTA file by ID, each with its defline and lines as written."""
    records: dict[str, str] = {}
    for chunk in path.read_text().split(">")[1:]:
        records[chunk.split(None, 1)[0]] = ">" + chunk
    return records


def sweep_inputs() -> None:
    meen = protein_records(FASTA / "MeenMJNV.faa")
    lv = protein_records(FASTA / "LvMJNV.faa")
    first_av = next(iter(protein_records(FASTA / "AvCLPV.faa").values()))
    (OUT / "e2e_protein_query.faa").write_text("".join(meen[key] for key in QUERIES[:6]) + first_av
                                               + meen[QUERIES[6]])
    others = [key for key in list(lv)[::17] if key not in HOMOLOGS][:6]
    (OUT / "e2e_protein_subject.faa").write_text("".join(lv[key] for key in (*HOMOLOGS[:6], *others, HOMOLOGS[6])))
    genome = "".join(line.strip() for line in (FASTA / "LvMJNV.fasta").read_text().splitlines()
                     if not line.startswith(">"))
    out = []
    for number, (start, end) in enumerate(WINDOWS, 1):
        seq = genome[start:end]
        if number == 2:
            seq = seq[:2500] + seq[2500:2800].lower() + seq[2800:]
        if number == 4:
            seq = seq[:3000] + "N" * 40 + seq[3040:]
        lines = "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70))
        out.append(f">LvMJNV_{start + 1}_{end} window {number} of LvMJNV\n{lines}\n")
    (OUT / "e2e_tblastn_subject.fna").write_text("".join(out))


def many_inputs() -> None:
    import random
    rng = random.Random(7)
    first = (OUT / "e2e_protein_query.faa").read_text().split(">")[1]
    base = first.split("\n", 1)[1].replace("\n", "")[:200]
    residues = "ACDEFGHIKLMNPQRSTVWY"
    codon = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
             "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT",
             "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}
    (OUT / "e2e_many_query.faa").write_text(f">q1\n{base}\n")
    variant = lambda: "".join(c if rng.random() > 0.3 else rng.choice(residues) for c in base)
    (OUT / "e2e_many_subject.faa").write_text("".join(f">s{i}\n{variant()}\n" for i in range(300)))
    (OUT / "e2e_many_subject.fna").write_text(
        "".join(f">n{i}\n" + "".join(codon.get(c, "NNN") for c in variant()) + "\n" for i in range(300)))


TITLES = ("hypothetical protein.", "hypothetical  protein   double  space", "protein. [Homo sapiens]",
          "protein, [Homo sapiens]", "a ,b ;c", "hyp ( fragment )", "TPA: hyp", "MULTISPECIES: hyp~",
          "x" + "y" * 60)


def titles_input() -> None:
    residues = (OUT / "e2e_protein_subject.faa").read_text().split(">")[1].split("\n", 1)[1]
    (OUT / "e2e_titles_subject.faa").write_text(
        "".join(f">t{i} {title}\n{residues}" for i, title in enumerate(TITLES)))


def amb_input() -> None:
    import random
    rng = random.Random(11)
    out = []
    for line in (OUT / "e2e_tblastn_subject.fna").read_text().split("\n"):
        if line.startswith(">") or not line:
            out.append(line)
            continue
        out.append("".join(rng.choice("NRYKMSWBDHV") if (c in "ACGT" and rng.random() < 0.02) else c
                           for c in line))
    (OUT / "e2e_amb_subject.fna").write_text("\n".join(out))


def main() -> None:
    records = [chunk.split("\n", 1)[1] for chunk in (OUT / "tblastx_many_subject.fasta").read_text().split(">")[1:4]]
    for name, deflines in (("punct_hits_subject.fna", DEFLINES), ("punct_hits_standin.fna", STAND_INS)):
        (OUT / name).write_text("".join(f">{defline}\n{seq}" for defline, seq in zip(deflines, records)))
    seq = "".join((OUT / "tblastx_many_query.fasta").read_text().split("\n")[1:])
    index = {base: i for i, base in enumerate("TCAG")}
    peptide = "".join(CODE[index[seq[i]] * 16 + index[seq[i + 1]] * 4 + index[seq[i + 2]]]
                      for i in range(30, len(seq) - 32, 3)).replace("*", "")
    (OUT / "punct_hits_query.faa").write_text(f">pep frame +1 of the many query\n{peptide}\n")
    sweep_inputs()
    many_inputs()
    titles_input()
    amb_input()


if __name__ == "__main__":
    main()
