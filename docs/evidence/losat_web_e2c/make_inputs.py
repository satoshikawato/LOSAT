#!/usr/bin/env python3
"""Generate the inputs that Session S07+ adds to LOSAT/tests/fasta/outfmt0/.

Every file comes from fixed random seeds or fixed text, so the bytes are always the same;
evidence.sha256 records them. The inputs of S06 and S07 come from
docs/evidence/losat_web_e2a/make_inputs.py, which this script does not change.

Usage: make_inputs.py [OUTDIR]   (default: LOSAT/tests/fasta/outfmt0 of this repository)
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
out = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "LOSAT/tests/fasta/outfmt0"
out.mkdir(parents=True, exist_ok=True)


def wrapped(defline: str, seq: str, width: int) -> str:
    return ">" + defline + "\n" + "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width))


def records(path: Path) -> list[tuple[str, str]]:
    chunks = [chunk.split("\n", 1) for chunk in path.read_text().split(">")[1:]]
    return [(head, body.replace("\n", "")) for head, body in chunks]


# iupac_query: the three multi_query records at three rates of IUPAC codes other than N.
# The ungapped Karlin block of each query depends on its composition, and so does the
# gapped block for gap costs beyond NCBI's tables (copied from the ungapped block).
rng = random.Random(7)
multi = records(out / "multi_query.fasta")
lines = []
for rate, codes in ((0.02, "RYKMSWBDHV"), (0.10, "RYKMSW"), (0.30, "RY")):
    for head, seq in multi:
        mutated = "".join(rng.choice(codes) if rng.random() < rate else base for base in seq)
        lines.append(wrapped(f"{head.split()[0]}_{int(rate * 100)} {head.split(' ', 1)[1]} with IUPAC codes",
                             mutated, 70))
(out / "iupac_query.fasta").write_text("".join(lines))


# rna_query: multi_query with U for T (upper and lower case) in its first two records; NCBI
# BLAST+ reads U as T.
lines = []
for number, (head, seq) in enumerate(multi):
    if number == 0:
        seq = seq.replace("T", "U")
    elif number == 1:
        seq = seq[:100] + seq[100:].replace("T", "u")
    lines.append(wrapped(head + (" as RNA" if number < 2 else ""), seq, 70))
(out / "rna_query.fasta").write_text("".join(lines))

# edge_short_invalid: a first query batch of 5000 residues (for the compact subjects),
# then two all-N queries of 70 residues together and a valid query. Later batches have at
# least 100 residues, so the N queries share a batch with the valid query: NCBI searches
# that batch and prints the normal footer for them.
compact = records(REPO / "LOSAT/tests/fasta/blastn_parity_compact.fasta")
alpha, beta = compact[0], compact[1]
(out / "edge_short_invalid.fasta").write_text(
    ">fill " + alpha[0].split()[0] + " repeated to 5000 residues\n" + (alpha[1] * 80)[:5000] + "\n"
    + ">shortN1 40 N\n" + "N" * 40 + "\n>shortN2 30 N\n" + "N" * 30 + "\n>" + beta[0] + "\n" + beta[1] + "\n")
