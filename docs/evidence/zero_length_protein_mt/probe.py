#!/usr/bin/env python3
"""NCBI BLAST+ 2.17.0 and LOSAT on a database holding a sequence of length 0, at 1 to 8 threads.

PD-LOSAT-NCBI-DEFECTS exception 36 (line RG-D2). Usage: python3 -I probe.py <ncbi_bin> <work_dir> [<losat>]
Builds the databases with NCBI makeblastdb, rewrites one sequence to length 0, and runs each search
at -num_threads 1, 2, 4, 8 (3 runs above 1). Writes results.tsv (one row per run) in work_dir.
SKIP_NCBI=1 runs LOSAT only.
"""
import hashlib
import os
import shutil
import struct
import subprocess
import sys

NCBI = sys.argv[1]
WORK = sys.argv[2]
LOSAT = sys.argv[3] if len(sys.argv) > 3 else None
HERE = os.path.dirname(os.path.abspath(__file__))
AA = "ACDEFGHIKLMNPQRSTVWY"
BIG = AA * 15
SMALL = AA * 15 + "Q" * 40 + AA * 3
SMALL2 = "W" * 10 + AA * 12 + "K" * 30
CODON = dict(A="GCT", C="TGT", D="GAT", E="GAA", F="TTT", G="GGT", H="CAT", I="ATT", K="AAA",
             L="CTG", M="ATG", N="AAT", P="CCG", Q="CAG", R="CGT", S="TCT", T="ACC", V="GTT",
             W="TGG", Y="TAT")


def dna(p):
    return "".join(CODON[c] for c in p)


def filler(i, n):
    # deterministic varied protein of length n, no strong hit to the query
    s = "".join("LKMNPGGSTWYRHHEDQQ"[(i * 7 + j * 13) % 18] for j in range(n))
    return s


def run(cmd, cwd, timeout=120):
    full = ["prlimit", "--as=8589934592", "--core=1", "timeout", str(timeout)] + cmd
    p = subprocess.run(full, cwd=cwd, capture_output=True)
    return p.returncode, p.stdout, p.stderr


def mkdb(dbdir, name, seqs, dbtype):
    fa = os.path.join(dbdir, name + ".fa")
    with open(fa, "w") as f:
        for sid, s in seqs:
            f.write(f">{sid}\n{s}\n")
    rc, out, err = run([os.path.join(NCBI, "makeblastdb"), "-in", fa, "-dbtype", dbtype, "-out", name,
                        "-title", name], dbdir)
    assert rc == 0, (name, out, err)


def parse_index(path, nucl):
    # version 5 index: version, seqtype, volume, title, lmdb name, date (each string = BE length + bytes),
    # then num(4 BE) total(8 LE) max(4 BE) and the offset arrays of (num+1) BE u32 each.
    b = open(path, "rb").read()
    ver = struct.unpack(">I", b[0:4])[0]
    assert ver == 5, ver
    p = 12
    for _ in range(3):
        n = struct.unpack(">I", b[p:p + 4])[0]
        p += 4 + n
    num = struct.unpack(">I", b[p:p + 4])[0]
    t = p + 16
    k = 3 if nucl else 2
    assert len(b) == t + k * 4 * (num + 1), (path, len(b), t, num)
    return b, t, num


def zero_prot(dbdir, name, k):
    pin = os.path.join(dbdir, name + ".pin")
    psq = os.path.join(dbdir, name + ".psq")
    b, t, num = parse_index(pin, False)
    hdr = list(struct.unpack(f">{num + 1}I", b[t:t + 4 * (num + 1)]))
    so = t + 4 * (num + 1)
    seq = list(struct.unpack(f">{num + 1}I", b[so:so + 4 * (num + 1)]))
    q = open(psq, "rb").read()
    parts = [q[seq[i]:seq[i + 1] - 1] for i in range(num)]
    parts[k] = b""
    out = bytearray(b"\x00")
    new = []
    for p in parts:
        new.append(len(out))
        out += p + b"\x00"
    new.append(len(out))
    total = sum(len(p) for p in parts)
    mx = max(len(p) for p in parts)
    nb = bytearray(b)
    nb[t - 12:t - 4] = struct.pack("<q", total)
    nb[t - 4:t] = struct.pack(">I", mx)
    nb[so:so + 4 * (num + 1)] = struct.pack(f">{num + 1}I", *new)
    open(pin, "wb").write(nb)
    open(psq, "wb").write(out)


def zero_nucl(dbdir, name, k):
    nin = os.path.join(dbdir, name + ".nin")
    nsq = os.path.join(dbdir, name + ".nsq")
    b, t, num = parse_index(nin, True)
    so = t + 4 * (num + 1)
    ao = so + 4 * (num + 1)
    seq = list(struct.unpack(f">{num + 1}I", b[so:so + 4 * (num + 1)]))
    amb = list(struct.unpack(f">{num + 1}I", b[ao:ao + 4 * (num + 1)]))
    q = open(nsq, "rb").read()
    lens = []
    packed = []
    ambs = []
    for i in range(num):
        pk = q[seq[i]:amb[i]]
        am = q[amb[i]:seq[i + 1]]
        if i == k:
            pk, am = b"\x00", b""
        packed.append(pk)
        ambs.append(am)
        lens.append((len(pk) - 1) * 4 + (pk[-1] & 3))
    out = bytearray(q[:seq[0]])
    ns, na = [], []
    for pk, am in zip(packed, ambs):
        ns.append(len(out))
        out += pk
        na.append(len(out))
        out += am
    ns.append(len(out))
    na.append(len(out))
    nb = bytearray(b)
    nb[t - 12:t - 4] = struct.pack("<q", sum(lens))
    nb[t - 4:t] = struct.pack(">I", max(lens))
    nb[so:so + 4 * (num + 1)] = struct.pack(f">{num + 1}I", *ns)
    nb[ao:ao + 4 * (num + 1)] = struct.pack(f">{num + 1}I", *na)
    open(nin, "wb").write(nb)
    open(nsq, "wb").write(out)
    return lens


def build(dbdir):
    os.makedirs(dbdir, exist_ok=True)
    many = [(f"f{i}", filler(i, 150 + (i * 37) % 250)) for i in range(300)]
    many[0] = ("small", SMALL)
    many[299] = ("small2", SMALL2)
    prot = {
        "db0": ([("big", BIG)], 0),
        "dbs": ([("small", SMALL)], None),
        "dbs2": ([("small2", SMALL2)], None),
        "s_first": ([("big", BIG), ("small", SMALL)], 0),
        "s_mid": ([("small", SMALL), ("big", BIG), ("small2", SMALL2)], 1),
        "s_last": ([("small", SMALL), ("big", BIG)], 1),
        "m_mid": (many, 150),
        "m_first": ([("z", BIG)] + many[1:], 0),
    }
    for name, (seqs, k) in prot.items():
        mkdb(dbdir, name, seqs, "prot")
        if k is not None:
            zero_prot(dbdir, name, k)
    for name, vols in {"a_0s": "db0 dbs", "a_s0": "dbs db0", "a_s0s": "dbs db0 dbs2"}.items():
        with open(os.path.join(dbdir, name + ".pal"), "w") as f:
            f.write(f"TITLE {name}\nDBLIST {vols}\n")
    nucl = {
        "n0": ([("big", dna(BIG))], 0),
        "ns": ([("small", dna(SMALL))], None),
        "ns_mid": ([("small", dna(SMALL)), ("big", dna(BIG)), ("small2", dna(SMALL2))], 1),
    }
    lens = {}
    for name, (seqs, k) in nucl.items():
        mkdb(dbdir, name, seqs, "nucl")
        if k is not None:
            lens[name] = zero_nucl(dbdir, name, k)
    with open(os.path.join(dbdir, "na_0s.nal"), "w") as f:
        f.write("TITLE na_0s\nDBLIST n0 ns\n")
    with open(os.path.join(dbdir, "q.fa"), "w") as f:
        f.write(f">q\n{AA * 3}\n")
    with open(os.path.join(dbdir, "qn.fa"), "w") as f:
        f.write(f">qn\n{dna(AA * 3)}\n")
    return lens


def main():
    dbdir = os.path.join(WORK, "db")
    shutil.rmtree(dbdir, ignore_errors=True)
    lens = build(dbdir)
    rows = []
    fmt = ["-outfmt", "6 std slen"]
    prot_dbs = ["a_0s", "a_s0", "a_s0s", "s_first", "s_mid", "s_last", "m_mid", "m_first"]
    progs = [
        ("blastp", ["blastp", "-query", "q.fa"]),
        ("blastp_ug", ["blastp", "-query", "q.fa", "-ungapped", "-comp_based_stats", "0"]),
        ("blastx", ["blastx", "-query", "qn.fa"]),
        ("psiblast", ["psiblast", "-query", "q.fa", "-num_iterations", "2"]),
    ]
    cases = [(p, db, a) for db in prot_dbs for p, a in progs]
    for db in ["na_0s", "ns_mid"]:
        cases += [("blastn_mega", db, ["blastn", "-query", "qn.fa"]),
                  ("blastn", db, ["blastn", "-task", "blastn", "-query", "qn.fa"]),
                  ("tblastn", db, ["tblastn", "-query", "q.fa"]),
                  ("tblastx", db, ["tblastx", "-query", "qn.fa"])]
    tools = ([] if os.environ.get("SKIP_NCBI") else [("ncbi", NCBI)]) + ([("losat", LOSAT)] if LOSAT else [])
    for tool, bindir in tools:
        for prog, db, args in cases:
            for th in (1, 2, 4, 8):
                for rep in range(1 if th == 1 else 3):
                    if tool == "ncbi":
                        cmd = [os.path.join(bindir, args[0])] + args[1:]
                    else:
                        cmd = [bindir] + args
                    cmd += ["-db", db, "-num_threads", str(th)] + fmt
                    rc, out, err = run(cmd, dbdir)
                    first = err.decode("utf-8", "replace").strip().split("\n")[0][:90] if err else ""
                    rows.append((tool, prog, db, th, rep, rc, out.count(b"\n"),
                                 hashlib.sha256(out).hexdigest()[:8], first))
                    print(*rows[-1], sep="\t", flush=True)
    with open(os.path.join(WORK, "results.tsv"), "w") as f:
        f.write(f"# nucleotide lengths after patch: {lens}\n")
        f.write("tool\tprog\tdb\tthreads\trep\trc\trows\tsha\tstderr1\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
