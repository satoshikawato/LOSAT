#!/usr/bin/env python3
"""Writes erfc_double_rounding.tsv and spouge_double_rounding.tsv (manifest.txt): the rows where the
LOSAT value before this change (BASE) differs from the oracle (NCBI) and the value with the x87 model
(X87) equals it. Each input is the output of one driver on erfc_z.in or spouge_grid.in (gen_inputs.py):
lines "E <z hex> <bits hex> <value>" or "S <y> <m> <n> <db_length> <bits hex> <value>".
usage: select_rows.py NCBI BASE X87 OUT.tsv"""
import struct
import sys


def rows(path):
    return [line.split() for line in open(path)]


def main():
    ncbi, base, x87, out = sys.argv[1:5]
    n, b, x = rows(ncbi), rows(base), rows(x87)
    assert len(n) == len(b) == len(x)
    with open(out, "w") as f:
        for r, q, t in zip(n, b, x):
            assert r[:-2] == q[:-2] == t[:-2]
            if r[-2] == q[-2] or r[-2] != t[-2]:
                continue
            if r[0] == "E":
                z = struct.unpack("<d", struct.pack("<Q", int(r[1], 16)))[0]
                if abs(z) >= 26.0:  # the underflow range of the double-double exponential (manifest.txt)
                    continue
            f.write("\t".join(r[1:-1]) + "\n")


if __name__ == "__main__":
    main()
