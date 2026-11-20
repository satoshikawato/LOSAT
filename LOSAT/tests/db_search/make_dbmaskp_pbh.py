#!/usr/bin/env python3
"""Series NP, M4b PR 1, review H1 (2026-11-21): the inputs pbh of dbmaskp/ (README.md there). Finds blastx
-db_hard_mask -subject_besthit inputs where the filter after every chunk (NCBI, blast_engine.c:586-591) differs from one filter after the merge (LOSAT before the
fix). Subjects: [A][C] | hard mask | [B], where B (a close copy of a query stretch) holds A (a weaker copy of nearly
the same stretch) and A holds C (a weaker copy starting a little before A). Runs the old and the new LOSAT on rounds of
candidates; writes the differing subjects as pbh.faa / edge_pbh.asnt / qbh.fna into OUT and prints a summary.
usage: make_dbmaskp_pbh.py WORK NCBI_BIN OLD_LOSAT NEW_LOSAT OUT (OLD: a build that filters once, NEW: one that
filters after every chunk; OUT: dbmaskp)
"""
import os
import random
import subprocess
import sys

WORK, NB, OLD, NEW, OUT = sys.argv[1:6]
AA = "ACDEFGHIKLMNPQRSTVWY"
CODONS = {}
bases = "TCAG"
code = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(bases):
    for j, b in enumerate(bases):
        for k, c in enumerate(bases):
            CODONS.setdefault(code[16 * i + 4 * j + k], []).append(a + b + c)
FMT = "6 qseqid sseqid qstart qend sstart send score"
ENV = dict(os.environ, BLAST_USAGE_REPORT="false")


def rnd(r, n):
    return "".join(r.choice(AA) for _ in range(n))


def mutate(r, s, rate):
    return "".join(r.choice(AA) if r.random() < rate else a for a in s)


def mask_file(path, masks, text):
    lines = ["Blast-db-mask-info ::= {", "  algo-id 2,", "  algo-program seg,", f'  algo-options "{text}",',
             "  masks {", "    masks {"]
    for n, (frm, to) in enumerate(masks):
        lines += ["      packed-int {", "        {", f"          from {frm},", f"          to {to - 1},",
                  f"          id local id {n + 1}", "        }", "      }" + ("," if n + 1 < len(masks) else "")]
    lines += ["    },", "    more FALSE", "  }", "}"]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def build(d, name, seqs, masks):
    with open(f"{d}/{name}.faa", "w") as f:
        for n, s in enumerate(seqs):
            f.write(f">s{n + 1}\n{s}\n")
    mask_file(f"{d}/edge_{name}.asnt", masks, "besthit")
    subprocess.run([f"{NB}/makeblastdb", "-in", f"{name}.faa", "-dbtype", "prot", "-out", name, "-title",
                    "besthit across a mask", "-mask_data", f"edge_{name}.asnt"], cwd=d, env=ENV, check=True,
                   stdout=subprocess.DEVNULL)


def rows(binary, d, db, query, ncbi=False):
    cmd = ([f"{NB}/blastx"] if ncbi else [binary, "blastx"]) + [
        "-db", db, "-query", query, "-db_hard_mask", "21", "-subject_besthit", "-outfmt", FMT, "-num_threads", "1"]
    p = subprocess.run(cmd, cwd=d, env=ENV, capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f"{cmd[0]} rc {p.returncode}: {p.stderr[:500]}")
    by = {}
    for line in p.stdout.splitlines():
        by.setdefault(line.split("\t")[1], []).append(line)
    return p.stdout, by


def subject_index(sseqid):
    """The 0-based ordinal of a subject named `s<n>` (or `gnl|BL_ORD_ID|<oid>`)."""
    if sseqid.startswith("s") and sseqid[1:].isdigit():
        return int(sseqid[1:]) - 1
    return int(sseqid.rsplit("|", 1)[-1])


def main():
    os.makedirs(WORK, exist_ok=True)
    r = random.Random(20261121)
    q = rnd(r, 220)
    qdna = "".join(r.choice(CODONS[a]) for a in q)
    with open(f"{WORK}/qbh.fna", "w") as f:
        f.write(">qbh query with one domain, back-translated\n" + "\n".join(qdna[i:i + 60] for i in range(0, len(qdna), 60)) + "\n")
    for rnd_no in range(6):
        seqs, masks = [], []
        for _ in range(400):
            b0 = r.randint(50, 100)
            b1 = b0 + r.randint(80, 110)
            a0 = b0 - r.randint(0, 3)
            a1 = b1 + r.randint(-8, 3)
            c0 = a0 - r.randint(0, 4)
            c1 = c0 + r.randint(35, 80)
            A = mutate(r, q[a0:a1], r.uniform(0.15, 0.35))
            C = mutate(r, q[c0:c1], r.uniform(0.2, 0.4))
            B = mutate(r, q[b0:b1], r.uniform(0.0, 0.15))
            left = rnd(r, r.randint(5, 20)) + A + rnd(r, r.randint(5, 15)) + C + rnd(r, r.randint(4, 10))
            mask = rnd(r, r.randint(20, 40))
            seq = left + mask + rnd(r, r.randint(4, 10)) + B + rnd(r, r.randint(5, 20))
            seqs.append(seq)
            masks.append((len(left), len(left) + len(mask)))
        d = f"{WORK}/r{rnd_no}"
        os.makedirs(d, exist_ok=True)
        build(d, "cand", seqs, masks)
        _, old = rows(OLD, d, "cand", "../qbh.fna")
        _, new = rows(NEW, d, "cand", "../qbh.fna")
        diff = sorted((s for s in set(old) | set(new) if old.get(s) != new.get(s)), key=subject_index)
        print(f"round {rnd_no}: {len(seqs)} subjects, {len(old)} with rows, differing {len(diff)}", flush=True)
        if not diff:
            continue
        pick = [subject_index(s) for s in diff[:3]]
        # the fixture: the differing subjects alone, if they still differ there; else the round's whole database
        os.makedirs(OUT, exist_ok=True)
        build(OUT, "pbh", [seqs[i] for i in pick], [masks[i] for i in pick])
        os.replace(f"{WORK}/qbh.fna", f"{OUT}/qbh.fna")
        o, old2 = rows(OLD, OUT, "pbh", "qbh.fna")
        n, new2 = rows(NEW, OUT, "pbh", "qbh.fna")
        c, _ = rows(None, OUT, "pbh", "qbh.fna", ncbi=True)
        print(f"fixture pbh: subjects {[i + 1 for i in pick]} old==new {o == n} ncbi==new {c == n} ncbi==old {c == o}")
        print("--- old\n" + o + "--- new\n" + n + "--- ncbi\n" + c, flush=True)
        if o != n:
            print("BH FOUND")
            return
        print("the picked subjects alone do not differ; next round")
        os.replace(f"{OUT}/qbh.fna", f"{WORK}/qbh.fna")
        for ext in ("faa", "pdb", "phr", "pin", "pot", "psq", "ptf", "pto", "paa", "pab", "pac", "pjs"):
            for name in (f"{OUT}/pbh.{ext}",):
                if os.path.exists(name):
                    os.remove(name)
        os.remove(f"{OUT}/edge_pbh.asnt")
    print("BH NOT FOUND")


main()
