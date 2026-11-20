#!/usr/bin/env python3
# Operation lists for vecprobe (seed fixed). Usage: gen_vectors.py main|large
import random, sys
which = sys.argv[1]
rng = random.Random(20261120)
edges = sorted({0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 127, 128, 129, 255, 256, 257, 1000, 1023, 1024, 1025, 4095, 4096, 4097})
rows = []
if which == "main":
    # every G n alone and in a ladder
    for n in edges:
        rows.append(f"G{n}")
    for a in edges:
        for b in edges:
            rows.append(f"G{a} G{b}")
    # copies: push a copy, then grow or shrink the vector, then take the copy back
    for a in edges:
        for b in edges:
            rows.append(f"G{a} S G{b} O")
            rows.append(f"G{a} S K O")
    for a in edges:
        for b in edges:
            rows.append(f"G{a} S N{b} O")
            rows.append(f"G{a} S G{b} S N{a//2} O O")
    for _ in range(6000):
        ops, size, stack = [], 0, []
        for _ in range(rng.randint(2, 24)):
            r = rng.random()
            if r < 0.38:
                n = rng.choice(edges) if rng.random() < 0.6 else rng.randint(0, 300)
                ops.append(f"G{n}"); size = n
            elif r < 0.60:
                ops.append("S"); stack.append(size)
            elif r < 0.82 and stack:
                ops.append("O"); size = stack.pop()
            elif r < 0.86:
                n = rng.randint(0, 40)
                ops.append(f"N{n}"); size = n
            elif r < 0.92:
                ops.append("K"); size = 0
            elif size > 0:
                ops.append("Z"); size -= 1
        rows.append(" ".join(ops))
elif which == "large":   # up to 2^26 records (512 MiB) in one vector
    for k in range(18, 27):
        p = 1 << k
        for ops in (f"G{p-1}", f"G{p}", f"G{p+1}", f"G{p} S", f"G{p-1} S G{p+1} O", f"G{p//2+1} S G{p} O", f"G{p} S G{p+1} S O O",
                    f"G{p//2} S G{p//2+1} S G{p} O O", f"G{p-1} S K O", f"G{p//2+1} S K G{p} O"):
            rows.append(ops)
print("\n".join(rows))
