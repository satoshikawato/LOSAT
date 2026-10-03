#!/usr/bin/env python3
"""SD fixture discrimination: which template combinations of each input pair catch each mutant.

For the existing pair (dc_query x dc_subject) and the divergent pair (chosen_query x
chosen_subject), runs NCBI BLAST+ and each LOSAT binary on the 18 combinations
(word 11/12 x coding/optimal/coding_and_optimal x 16/18/21) and the megablast template row,
outfmt 6, and records whether stdout, stderr and exit status equal NCBI's.

Usage: compare.py BIN [BIN ...]  (writes compare.tsv and prints a summary)
"""
from __future__ import annotations

import concurrent.futures
import itertools
import os
import subprocess
import sys
from pathlib import Path

M = Path.home() / ".cache/losat-web-gui-target/sd-mutants"
REPO = Path("/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/tests/fasta/outfmt0")
DIV = Path.home() / ".cache/losat-web-gui-target/sd-disc-fixtures"
NCBI = "/home/kawato/micromamba/bin/blastn"
PAIRS = {"existing": (REPO / "dc_query.fasta", REPO / "dc_subject.fasta"),
         "divergent": (DIV / "chosen_query.fasta", DIV / "chosen_subject.fasta")}
COMBOS = [(f"w{w}_{t}_{l}", ["-task", "dc-megablast", "-word_size", str(w), "-template_type", t,
                             "-template_length", str(l)])
          for w, t, l in itertools.product((11, 12), ("coding", "optimal", "coding_and_optimal"), (16, 18, 21))]
COMBOS.append(("megablast_w11_coding_18", ["-task", "megablast", "-word_size", "11", "-template_type", "coding",
                                           "-template_length", "18"]))
ENV = {k: v for k, v in os.environ.items()
       if k not in ("BATCH_SIZE", "CHUNK_SIZE", "OVERLAP_CHUNK_SIZE", "BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE")}


def run(argv: list[str]) -> tuple[int, bytes, bytes]:
    result = subprocess.run(argv, capture_output=True, env=ENV, cwd=M)
    return result.returncode, result.stdout, result.stderr


def main() -> int:
    bins = sys.argv[1:]
    jobs = [(pair, name, args) for pair in PAIRS for name, args in COMBOS]
    def common(pair, args):
        q, s = PAIRS[pair]
        return ["-query", str(q), "-subject", str(s), *args, "-outfmt", "6"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        oracle = dict(zip([j[:2] for j in jobs], pool.map(lambda j: run([NCBI, *common(j[0], j[2])]), jobs)))
        rows = []
        for binary in bins:
            outs = list(pool.map(lambda j: run([binary, "blastn", *common(j[0], j[2])]), jobs))
            for job, out in zip(jobs, outs):
                rows.append((Path(binary).name, job[0], job[1], "same" if out == oracle[job[:2]] else "DIFF"))
    with open(M / "compare.tsv", "w") as handle:
        handle.write("binary\tpair\tcombination\tresult\n")
        handle.writelines("\t".join(r) + "\n" for r in rows)
    for binary in bins:
        name = Path(binary).name
        for pair in PAIRS:
            diff = [r[2] for r in rows if r[0] == name and r[1] == pair and r[3] == "DIFF"]
            print(f"{name}\t{pair}\t{len(diff)} of {len(COMBOS)} differ from NCBI\t{' '.join(diff)}")
    distinct = {pair: len({oracle[(pair, n)][1] for n, _ in COMBOS[:18]}) for pair in PAIRS}
    print(f"NCBI distinct outputs of the 18 combinations: {distinct}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
