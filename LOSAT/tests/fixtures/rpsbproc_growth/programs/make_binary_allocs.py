#!/usr/bin/env python3
# make_binary_allocs.py: for every row of motifs.tsv (label m01, m02, ... as runall_motifs.sh named the runs), the
# allocations of NCBI's rpsbproc 0.5.1 (mallog.so, >= 8 bytes) of the run minus those of the baseline run of the
# baseline motif of the same shape (motifs.tsv third column), as a signed count per size: "size:+n" more allocations,
# "size:-n" fewer. Allocations that do not come from the string or the vector of the motif check are removed from both
# runs first: the pattern's vector of positions (72 * 2^k bytes, 72 bytes a position) and the copy of the motif text
# (its length + 1 bytes when it is over 15 characters).
import collections, sys
def allocs(label, motif):
    c = collections.Counter()
    for l in open(f"rps/mal_{label}.log"):
        t = l.split()
        if t[0] not in ("m", "r"):
            continue
        n = int(t[1])
        if n % 72 == 0 and (n // 72) & (n // 72 - 1) == 0:
            continue
        c[n] += 1
    if len(motif) > 15:
        c[len(motif) + 1] -= 1
    return c
rows = [l.rstrip("\n").split("\t") for l in open("motifs.tsv")]
for i, (site, motif, basem) in enumerate(rows, 1):
    c = allocs("m%02d" % i, motif)
    b = allocs("b%02d" % i, basem)
    d = collections.Counter()
    for k in set(c) | set(b):
        v = c[k] - b[k]
        if v:
            d[k] = v
    print(site, motif, basem, " ".join(f"{k}:{v:+d}" for k, v in sorted(d.items())), sep="\t")
