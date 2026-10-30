#!/usr/bin/env python3
"""Writes the expected descendants of the `#[ignore]` test `real_file_equals_sqlite3` of
LOSAT/src/seqdb/tax_sqlite.rs (`LOSAT_TEST_TAX_SQLITE_EXPECTED`) from the real
`taxonomy4blast.sqlite3` (`LOSAT_TEST_TAX_SQLITE`), with Python's `sqlite3`.

Usage: tax_sqlite_ref.py TAXONOMY4BLAST_SQLITE OUT

OUT has one line `T: sorted descendants` per taxid: taxid 1 first (NCBI's query never ends
for it: taxid 1 is its own parent; the line is a walk of the whole tree that visits each
taxid once), then 12 fixed taxids (0 and -5 have none) and 200 drawn with seed 20261029, half
internal nodes and half leaves, each from NCBI's own query (tax4blastsqlite.cpp) run by SQLite.
The descendants leave out the taxid itself; taxids are read as NCBI's `Int4` (`get_int`).
"""
import random
import sqlite3
import sys

QUERY = """
WITH RECURSIVE descendants_from(taxid, parent) AS
    (SELECT ?1, NULL
     UNION ALL
     SELECT T.taxid, T.parent
       FROM descendants_from D, TaxidInfo T
       WHERE D.taxid = T.parent)
SELECT taxid from descendants_from;
"""
FIXED = [0, -5, 2, 561, 562, 543, 1224, 9606, 40674, 7742, 10090, 99999999]


def get_int(v):
    return ((v + (1 << 31)) % (1 << 32)) - (1 << 31)


def pick(con, n, seed):
    rnd = random.Random(seed)
    rows = con.execute("SELECT taxid FROM TaxidInfo").fetchall()
    parents = {r[0] for r in con.execute("SELECT DISTINCT parent FROM TaxidInfo WHERE parent IS NOT NULL")}
    inner = [r[0] for r in rows if r[0] in parents and r[0] != 1 and r[0] < 2**31]
    leaves = [r[0] for r in rows if r[0] not in parents and r[0] < 2**31]
    return rnd.sample(inner, n // 2) + rnd.sample(leaves, n - n // 2)


def tree1(con):
    kids = {}
    for t, p in con.execute("SELECT taxid, parent FROM TaxidInfo WHERE parent IS NOT NULL"):
        kids.setdefault(p, []).append(t)
    seen, todo = {1}, [1]
    while todo:
        for k in kids.get(todo.pop(), []):
            if k not in seen:
                seen.add(k)
                todo.append(k)
    return sorted(get_int(t) for t in seen if get_int(t) != 1)


def descendants(con, t):
    if t <= 0:
        return []
    return sorted(g for g in (get_int(r[0]) for r in con.execute(QUERY, (t,))) if g != t)


def main(db, out):
    con = sqlite3.connect("file:%s?mode=ro" % db, uri=True)
    lines = ["1: " + " ".join(map(str, tree1(con)))]
    for t in FIXED + pick(con, 200, 20261029):
        lines.append("%d: %s" % (t, " ".join(map(str, descendants(con, t)))))
    with open(out, "w") as f:
        f.write("".join(line.rstrip() + "\n" for line in lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
