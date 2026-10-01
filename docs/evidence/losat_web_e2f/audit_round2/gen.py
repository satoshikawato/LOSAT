#!/usr/bin/env python3
"""Deterministic generators for the round-2 hit-list audit (comparison only)."""
import random
from pathlib import Path
FASTA = Path("/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/tests/fasta")
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")

def genome(name):
    return "".join(l.strip() for l in (FASTA / name).read_text().splitlines()[1:])

def rc(s): return s.translate(COMP)[::-1]

def write(path, recs):
    Path(path).write_text("".join(f">{n}\n{s}\n" for n, s in recs))
    return str(path)

def rnd(rng, n): return "".join(rng.choice("ACGT") for _ in range(n))

def subst(rng, s, rate):
    if rate <= 0: return s
    out = list(s)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([c for c in "ACGT" if c != out[i]])
    return "".join(out)

def indel(rng, s, rate):
    out = []
    for c in s:
        r = rng.random()
        if r < rate * 0.5: continue
        if r < rate: out += [c, rng.choice("ACGT")]
        else: out.append(c)
    return "".join(out)

def subjects(kind, q, rng, n, lo=0, hi=None, lens=(60, 80, 100, 150, 250), flank="random"):
    """n subject records built from windows of q[lo:hi]."""
    hi = hi or len(q)
    recs = []
    base = []  # for dups
    def win(L):
        a = rng.randrange(lo, hi - L)
        return q[a:a + L]
    def fl():
        if flank == "none": return "", ""
        if flank == "fixed": return rnd(rng, 50), rnd(rng, 50)
        return rnd(rng, rng.randint(0, 200)), rnd(rng, rng.randint(0, 200))
    if kind == "dups":
        k = max(3, n // 40)
        base = [(win(rng.choice(lens)), rnd(rng, 40), rnd(rng, 40)) for _ in range(k)]
    for i in range(n):
        a, b = fl()
        if kind == "ties_exact":
            L = lens[0]
            s = a + win(L) + b
        elif kind == "ties_fixedlen":
            L = lens[0]
            s = win(L)  # identical length, exact
        elif kind == "classes":      # few score classes -> many ties
            s = a + win(rng.choice(lens)) + b
        elif kind == "dups":
            w, x, y = base[i % len(base)]
            s = x + w + y
        elif kind == "mixed":
            L = rng.choice(lens)
            s = a + subst(rng, win(L), rng.choice([0, 0, 0.02, 0.05, 0.08])) + b
        elif kind == "rc_mix":
            L = rng.choice(lens)
            w = subst(rng, win(L), rng.choice([0, 0, 0.03]))
            if rng.random() < 0.5: w = rc(w)
            s = a + w + b
        elif kind == "both_strands":
            L1, L2 = rng.sample(sorted(set(lens)) if len(set(lens)) > 1 else list(lens) * 2, 2)
            s = a + win(L1) + rnd(rng, 40) + rc(win(L2)) + b
        elif kind == "multi_hsp":
            parts = []
            for _ in range(rng.randint(2, 4)):
                w = win(rng.choice(lens))
                if rng.random() < 0.4: w = rc(w)
                parts.append(w)
            s = a + (rnd(rng, 60)).join(parts) + b
        elif kind == "indel":
            L = rng.choice(lens)
            w = indel(rng, subst(rng, win(L), 0.02), 0.03)
            if rng.random() < 0.5: w = rc(w)
            s = a + w + b
        elif kind == "same_window":   # all identical sequence
            if not base: base.append(win(lens[0]))
            s = base[0]
        else:
            raise ValueError(kind)
        recs.append((f"S{i}", s))
    return recs
