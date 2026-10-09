#!/usr/bin/env python3
"""Deterministic read-heavy inputs for the E2h V-PERF cases (not committed; regenerated from a seed).

Writes to $BUILD_ROOT/sfb-e2h/perf-inputs/ (or --out):

  q100k_300nt.fna        100000 nucleotide records of 300 nt (70 columns); every 100th record is a slice
                         of LC738884 with a few substitutions, so the search has hits
  genome5mb_1line.fna    one record, 5 Mb on ONE line (random, with three mutated copies of AP027152
                         inserted so the search has hits)
  genome5mb_80col.fna    the same record in 80-column lines
  protein_many.faa       20000 protein records of 300 aa (60 columns); every 50th is a slice of PajaWSV.faa
  inputs.sha256          SHA-256 of every file written

Usage: gen_perf_inputs.py [--out DIR] [--repo DIR] [--queries N] [--genome-mb N] [--proteins N]
Same seed, same bytes. The perf cases are in perf_cases.py (blastn-q100k, blastn-genome-1line,
blastn-genome-80col, blastp-many).
"""
import argparse
import hashlib
import os
import random
import sys
from pathlib import Path

SEED = 20261008
DNA = bytes.maketrans(bytes(range(256)), bytes(b"ACGT"[b & 3] for b in range(256)))
AA = b"ACDEFGHIKLMNPQRSTVWY"
PROT = bytes.maketrans(bytes(range(256)), bytes(AA[b % 20] for b in range(256)))


def read_fasta_seq(path: Path) -> bytes:
    return b"".join(line.strip() for line in path.read_bytes().split(b"\n") if line and not line.startswith(b">")).upper()


def mutate(rng: random.Random, seq: bytes, alphabet: bytes, rate: float) -> bytes:
    out = bytearray(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = alphabet[rng.randrange(len(alphabet))]
    return bytes(out)


def wrap(seq: bytes, width: int) -> bytes:
    return b"\n".join(seq[i:i + width] for i in range(0, len(seq), width)) + b"\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--repo", type=Path)
    ap.add_argument("--queries", type=int, default=100000)
    ap.add_argument("--genome-mb", type=int, default=5)
    ap.add_argument("--proteins", type=int, default=20000)
    args = ap.parse_args()
    repo = args.repo or Path(os.environ.get("WT") or Path(__file__).resolve().parents[4])
    fasta = repo / "LOSAT" / "tests" / "fasta"
    out = args.out or Path(os.environ.get("BUILD_ROOT", "/home/kawato/.cache/losat-work")) / "sfb-e2h" / "perf-inputs"
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    lc = read_fasta_seq(fasta / "LC738884.fasta")
    with (out / "q100k_300nt.fna").open("wb") as h:
        for n in range(args.queries):
            if n % 100 == 0:
                start = rng.randrange(0, len(lc) - 300)
                seq = mutate(rng, lc[start:start + 300], b"ACGT", 0.02)
            else:
                seq = rng.randbytes(300).translate(DNA)
            h.write(b">q%d\n" % n + wrap(seq, 70))

    ap_seq = read_fasta_seq(fasta / "AP027152.fasta")
    genome = bytearray(rng.randbytes(args.genome_mb * 1_000_000).translate(DNA))
    for copy in range(3):
        piece = mutate(rng, ap_seq, b"ACGT", 0.03)
        pos = rng.randrange(0, len(genome) - len(piece))
        genome[pos:pos + len(piece)] = piece
    genome = bytes(genome)
    (out / "genome5mb_1line.fna").write_bytes(b">genome5mb perf subject\n" + genome + b"\n")
    (out / "genome5mb_80col.fna").write_bytes(b">genome5mb perf subject\n" + wrap(genome, 80))

    pj = read_fasta_seq(fasta / "PajaWSV.faa")
    with (out / "protein_many.faa").open("wb") as h:
        for n in range(args.proteins):
            if n % 50 == 0 and len(pj) > 300:
                start = rng.randrange(0, len(pj) - 300)
                seq = mutate(rng, pj[start:start + 300], AA, 0.05)
            else:
                seq = rng.randbytes(300).translate(PROT)
            h.write(b">p%d\n" % n + wrap(seq, 60))

    names = ["q100k_300nt.fna", "genome5mb_1line.fna", "genome5mb_80col.fna", "protein_many.faa"]
    (out / "inputs.sha256").write_text("".join(f"{hashlib.sha256((out / n).read_bytes()).hexdigest()}  {n}\n" for n in names))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
