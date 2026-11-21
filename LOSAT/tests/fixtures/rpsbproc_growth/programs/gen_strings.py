#!/usr/bin/env python3
# Generates the operation lists for strprobe (seed fixed). Usage: gen_strings.py main|large|exc
import random, sys
M = (1 << 63) - 1
which = sys.argv[1]
rng = random.Random(20261120)
edges = set()
for v in (0, 1, 2, 3, 7, 8, 14, 15, 16, 17, 29, 30, 31, 32, 33, 59, 60, 61, 62, 63, 64, 100, 127, 128, 129, 200, 255, 256, 257):
    edges.add(v)
for k in range(9, 23):
    for d in (-2, -1, 0, 1, 2):
        edges.add((1 << k) + d)
for k in range(5, 22):  # 15 * 2^j neighbourhoods: capacities that doubling from 15 passes through
    for d in (-1, 0, 1):
        edges.add(15 * (1 << (k - 4)) + d)
edges = sorted(v for v in edges if v >= 0)
rows = []
if which == "main":
    for op in "RAP":
        for v in ([0] if op == "P" else edges):
            rows.append(f"{op}{v}" if op != "P" else "P")
    small = [0, 1, 14, 15, 16, 17, 29, 30, 31, 32, 45, 46, 60, 61, 100, 128, 255, 256, 1000, 1024]
    ops2 = [f"{o}{v}" for o in "RA" for v in small] + ["P", "C"]
    for a in ops2:
        for b in ops2:
            rows.append(f"{a} {b}")
    for a in ["A15", "A16", "R15", "R16", "A30", "A31", "R100", "A1000"]:
        for b in [f"{o}{v}" for o in "RA" for v in (0, 1, 15, 16, 29, 30, 31, 61, 100, 129, 2047, 4097)] + ["P"]:
            for c in ["P", "A1", "R1", "R400", "A400", "C P", "C A20"]:
                rows.append(f"{a} {b} {c}")
    for _ in range(4000):
        n = rng.randint(2, 16)
        ops = []
        for _ in range(n):
            r = rng.random()
            if r < 0.30:
                ops.append("P")
            elif r < 0.37:
                ops.append("C")
            else:
                v = rng.choice(edges) if rng.random() < 0.5 else rng.randint(0, 3000)
                ops.append(rng.choice("RAA") + str(v))
        rows.append(" ".join(ops))
elif which == "large":
    for k in range(22, 30):  # 4 MiB .. 512 MiB
        p = 1 << k
        for ops in (f"A{p-1}", f"A{p}", f"A{p+1}", f"R{p-1}", f"R{p}", f"R{p+1}", f"A{p-1} P", f"A{p} P", f"A{p-2} P P P",
                    f"R{p-1} P", f"R{p-1} A1 P", f"A{p//2} A{p//2}", f"A{p//2} A{p//2+1}", f"A{p//2} A{p//2} P C P", f"A{p//4*3} R{p}"):
            rows.append(ops)
    # 15 * 2^k capacities reached by appending one byte at a time from SSO (as `push_back`s do)
    rows.append("A15 " + " ".join(["P"] * 3))
    rows.append("A3000000 P P C A3000000 P R3000001 R6000001 P")
elif which == "host":   # ws only: 2^30 bytes or more at once (LOSAT refuses these as host dependent)
    for ops in ("A1073741823", "A1073741824", "A1073741822 P", "R1073741823", "A536870910 A1", "A536870911 A1 P",
                "A536870910 P", "A536870911 P", "A1000000000 P", "A2147483646 P", "A2147483647 P", "A1610612735 P"):
        rows.append(ops)
elif which == "exc":
    two47 = 1 << 47
    for ops in (
        f"A{M}", f"A{M-1}", f"A{M-2}", f"A1 A{M}", f"A1 A{M-1}", f"A1 A{M-2}", f"A16 A{M-16}", f"A16 A{M-15}", f"A16 A{M-17}",
        f"R{M}", f"R{M+1}", f"R{M-1}", f"R{(1<<64)-1}", f"A{(1<<64)-1}", f"A1 A{(1<<64)-1}", f"A{(1<<64)-16}", f"A16 A{(1<<64)-16}",
        f"R{1<<62}", f"A{1<<62}", f"A{(1<<62)-1}", f"A1 A{(1<<62)-1}", f"A1 A{(1<<62)-2}", f"R{(1<<62)-1}", f"R{(1<<62)+1}",
        f"R{two47-1}", f"R{two47}", f"A{two47-1}", f"A{two47}", f"A{two47+1}", f"A1 A{two47-2}", f"A1 A{two47-1}", f"A1 A{two47}",
        f"A16 A{two47}", f"A16 A{two47-16}", f"A16 A{two47-17}", f"A100 R{two47}", f"A100 A{two47}", "A100 P A100000 R9223372036854775808",
        f"A1 A{M} P", f"A1 A{M} A1", f"A15 A{M-15}", f"A15 A{M-14}", f"A15 A{M-16}", f"A100 A{M-100}", f"A100 A{M-99}",
    ):
        rows.append(ops)
print("\n".join(rows))
