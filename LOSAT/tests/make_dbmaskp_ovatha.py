#!/usr/bin/env python3
"""The database `ov_atha` of overlapping seg masks (line NP, series M, M4b-3b; fixture agent of ad-M4b-6, 2026-11-22):
40 proteins of the A. thaliana TAIR10 set (`ncbi-parity-20261019/M/m4b3/data/atha.faa`) whose `segmasker` (NCBI BLAST+
2.17.0 defaults, `-outfmt interval`) masks overlap, so that `makeblastdb -mask_data` of the maskinfo_asn1_bin output
holds inverted hard ranges (a mask that starts before the previous one ends).

Pick (deterministic, in the order of the sequence names): only proteins of 150-600 letters whose masks are sorted by
start, none reaching the last letter (a masked tail is another shape, `pe_tail`); then the first 10 with at least two
consecutive overlaps (mask i+1 starts at or before the last letter of mask i), the first 5 of the rest with a nested
mask (mask i+1 inside mask i), and the first of the rest with one overlap up to 40. `qova.faa` is 4 of the 40: the ones
with the most overlaps (the first by name among equals), `qova.fna` their back-translation (random codon per residue,
seed 20261122, the code of make_dbmaskp_ov.py). Headers keep the sequence name only.

    python3 make_dbmaskp_ovatha.py --atha PATH/atha.faa [--out DIR] [--ncbi-bin DIR]

writes DIR/ov_atha.faa, DIR/ov_athaw.faa (without the six proteins of WITHOUT), DIR/qova.faa, DIR/qova.fna, and with --ncbi-bin `segmasker -in ov_atha.faa -infmt fasta
-outfmt maskinfo_asn1_bin -out ov_atha.asnb` and `makeblastdb -in ov_atha.faa -dbtype prot -out ov_atha -title
"overlapping seg masks" -mask_data ov_atha.asnb` (BLAST_USAGE_REPORT=false) and removes ov_atha.pjs.
"""
import argparse
import os
import random
import subprocess
from pathlib import Path

from make_dbmaskp_ov import back_translate

# `ov_athaw` (decision 52): ov_atha without the six proteins whose masks leave diagonal entries that the next subject
# reads in the one-hit scan (-window_size 0), which LOSAT refuses before the search.
WITHOUT = ("AT1G02550.1", "AT1G02800.1", "AT1G04080.2", "AT1G04080.5", "AT1G06760.1", "AT1G09460.1")


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line[0] == ">":
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def intervals(ncbi_bin, faa, env):
    text = subprocess.run([str(Path(ncbi_bin) / "segmasker"), "-in", str(faa), "-infmt", "fasta", "-outfmt", "interval"],
                          env=env, check=True, capture_output=True, text=True).stdout
    iv, cur = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            iv[cur] = []
        else:
            a, b = line.split(" - ")
            iv[cur].append((int(a), int(b)))
    return iv


def shape(masks):
    """(overlaps, nested): consecutive pairs, inclusive intervals."""
    ov = sum(1 for p, q in zip(masks, masks[1:]) if q[0] <= p[1])
    ne = sum(1 for p, q in zip(masks, masks[1:]) if q[0] >= p[0] and q[1] <= p[1])
    return ov, ne


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atha", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "db_search" / "dbmaskp")
    ap.add_argument("--ncbi-bin", required=True)
    args = ap.parse_args()
    env = dict(os.environ, BLAST_USAGE_REPORT="false")
    seqs = read_fasta(args.atha)
    iv = intervals(args.ncbi_bin, args.atha, env)
    cand = {}
    for name in sorted(seqs):
        n = len(seqs[name])
        m = iv.get(name, [])
        if not 150 <= n <= 600 or m != sorted(m) or any(b < a for a, b in m) or (m and m[-1][1] >= n - 1):
            continue
        ov, ne = shape(m)
        if ov:
            cand[name] = (ov, ne)
    names = sorted(cand)
    pick = [n for n in names if cand[n][0] >= 2][:10]
    pick += [n for n in names if n not in pick and cand[n][1] >= 1][:5]
    pick += [n for n in names if n not in pick][:40 - len(pick)]
    pick.sort()
    assert len(pick) == 40, len(pick)
    out = args.out
    with open(out / "ov_atha.faa", "w") as f:
        for n in pick:
            s = seqs[n]
            f.write(f">{n}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    top = sorted(pick, key=lambda n: (-cand[n][0], n))[:4]
    top.sort()
    r = random.Random(20261122)
    with open(out / "qova.faa", "w") as fa, open(out / "qova.fna", "w") as fn:
        for n in top:
            s = seqs[n]
            fa.write(f">{n}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
            t = back_translate(r, s)
            fn.write(f">{n}\n" + "\n".join(t[i:i + 60] for i in range(0, len(t), 60)) + "\n")
    two = sum(1 for n in pick if cand[n][0] >= 2)
    nest = sum(1 for n in pick if cand[n][1] >= 1)
    print(f"picked 40: {two} with >= 2 consecutive overlaps, {nest} with a nested mask; queries {top}")
    exe = Path(args.ncbi_bin)
    with open(out / "ov_athaw.faa", "w") as f:
        for n in pick:
            if n not in WITHOUT:
                s = seqs[n]
                f.write(f">{n}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    for db, title in (("ov_atha", "overlapping seg masks"), ("ov_athaw", "overlapping seg masks without six")):
        subprocess.run([str(exe / "segmasker"), "-in", f"{db}.faa", "-infmt", "fasta", "-outfmt", "maskinfo_asn1_bin",
                        "-out", f"{db}.asnb"], cwd=out, env=env, check=True)
        subprocess.run([str(exe / "makeblastdb"), "-in", f"{db}.faa", "-dbtype", "prot", "-out", db, "-title", title,
                        "-mask_data", f"{db}.asnb"], cwd=out, env=env, check=True, stdout=subprocess.DEVNULL)
        (out / f"{db}.pjs").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
