import random, sys
sys.path.insert(0, '.')
import gen
from gen import write, subjects, genome, rnd, rc, subst

W = "w"
edl, sak = genome("EDL933.fna"), genome("Sakai.fna")

def chunks(lengths, chunk_size, overlap=100):
    total = sum(lengths)
    count = total // (chunk_size - overlap)
    if count <= 1: return [(0, total)]
    size = (total + (count - 1) * overlap) // count
    if count < size - overlap: size += 1
    out, start = [], 0
    for k in range(count):
        end = start + size
        if end >= total or k + 1 == count: end = total
        out.append((start, end))
        start += size - overlap
        if start > total or end == total: break
    return out

cases = []
def add(name, cat, q, s, args):
    cases.append(dict(name=name, cat=cat, query=q, subject=s, args=list(args)))

# ---------- A. 560..2000 subjects, default hitlist ------------
qA = {}
def qfile(tag, seq): 
    return write(f"{W}/q_{tag}.fa", [(tag, seq)])
qs = {"e1": edl[100000:400000], "e2": edl[2000000:2300000], "s1": sak[800000:1100000], "e3": edl[4000000:4150000]}
qf = {k: qfile(k, v) for k, v in qs.items()}
kinds = ["classes", "ties_exact", "ties_fixedlen", "dups", "mixed", "rc_mix", "both_strands", "multi_hsp", "indel", "same_window"]
i = 0
for kind in kinds:
    for n, qk, task, fmt in ((560, "e1", "blastn", "6"), (1200, "e2", "megablast", "6"), (2000, "s1", "blastn", "6")):
        i += 1
        rng = random.Random(f"A-{kind}-{n}")
        lens = (60, 80, 100, 150, 250) if task == "blastn" else (45, 60, 100, 150, 250)
        if kind in ("ties_exact", "ties_fixedlen"): lens = (60 if task=="blastn" else 45,)
        sf = write(f"{W}/A_{kind}_{n}.fa", subjects(kind, qs[qk], rng, n, lens=lens, flank="fixed" if kind=="ties_exact" else "random"))
        add(f"A_{kind}_{n}_{task}", "A_big_default", qf[qk], sf, ["-task", task, "-outfmt", fmt])
# outfmt 7 / 0 and others on a few
for kind, n, fmt, task in (("classes", 700, "7", "blastn"), ("dups", 900, "0", "blastn"), ("mixed", 650, "0", "megablast"),
                           ("ties_exact", 600, "7", "megablast"), ("multi_hsp", 800, "7", "blastn"), ("rc_mix", 560, "0", "blastn")):
    rng = random.Random(f"A2-{kind}-{n}-{fmt}")
    lens = (60, 80, 100, 150, 250) if task == "blastn" else (45, 60, 100, 150, 250)
    if kind == "ties_exact": lens = (60 if task=="blastn" else 45,)
    sf = write(f"{W}/A2_{kind}_{n}_{fmt}.fa", subjects(kind, qs["e1"], rng, n, lens=lens, flank="fixed" if kind=="ties_exact" else "random"))
    add(f"A2_{kind}_{n}_{task}_o{fmt}", "A_big_default", qf["e1"], sf, ["-task", task, "-outfmt", fmt])
# boundary counts around 550 / 551 / 560 (exactly 549..552 hit subjects)
for n in (549, 550, 551, 552):
    rng = random.Random(f"A3-{n}")
    sf = write(f"{W}/A3_{n}.fa", subjects("classes", qs["e1"], rng, n, lens=(60, 70, 80)))
    add(f"A3_edge{n}", "A_big_default", qf["e1"], sf, ["-task", "blastn", "-outfmt", "6"])

# ---------- B. -max_target_seqs ----------
mts_vals = [1, 2, 3, 4, 5, 6, 10, 11]
for mts in mts_vals:
    for kind, n in (("classes", 40), ("ties_exact", 300), ("dups", 60)):
        rng = random.Random(f"B-{kind}-{n}-{mts}")
        lens = (60,) if kind == "ties_exact" else (60, 80, 100, 150, 250)
        sf = write(f"{W}/B_{kind}_{n}_{mts}.fa", subjects(kind, qs["e1"], rng, n, lens=lens, flank="fixed" if kind=="ties_exact" else "random"))
        add(f"B_mts{mts}_{kind}_{n}", "B_max_target_seqs", qf["e1"], sf, ["-task", "blastn", "-outfmt", "6", "-max_target_seqs", str(mts)])
for mts in (100, 550, 600):
    for kind, n in (("classes", 700), ("mixed", 1500)):
        rng = random.Random(f"B2-{kind}-{n}-{mts}")
        sf = write(f"{W}/B2_{kind}_{n}_{mts}.fa", subjects(kind, qs["e2"], rng, n))
        add(f"B2_mts{mts}_{kind}_{n}", "B_max_target_seqs", qf["e2"], sf, ["-task", "blastn", "-outfmt", "6", "-max_target_seqs", str(mts)])
for mts, fmt, task in ((3, "0", "blastn"), (5, "7", "megablast"), (1, "0", "megablast"), (2, "7", "blastn"), (11, "0", "blastn"), (10, "7", "blastn")):
    rng = random.Random(f"B3-{mts}-{fmt}")
    lens = (60, 80, 100, 150) if task == "blastn" else (45, 60, 100, 150)
    sf = write(f"{W}/B3_{mts}_{fmt}_{task}.fa", subjects("multi_hsp", qs["s1"], rng, 120, lens=lens))
    add(f"B3_mts{mts}_o{fmt}_{task}", "B_max_target_seqs", qf["s1"], sf, ["-task", task, "-outfmt", fmt, "-max_target_seqs", str(mts)])

# ---------- C. -subject_besthit ----------
for k, (kind, n, mts) in enumerate((("classes", 600, None), ("multi_hsp", 600, None), ("both_strands", 700, None), ("rc_mix", 560, None),
                                     ("multi_hsp", 80, 3), ("both_strands", 50, 1), ("classes", 40, 5), ("indel", 600, None),
                                     ("dups", 700, None), ("multi_hsp", 1000, 100), ("mixed", 650, None), ("same_window", 600, None))):
    rng = random.Random(f"C-{k}")
    sf = write(f"{W}/C_{k}.fa", subjects(kind, qs["e1"], rng, n))
    args = ["-task", "blastn", "-outfmt", "6", "-subject_besthit"] + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"C_besthit_{k}_{kind}_{n}_{mts}", "C_subject_besthit", qf["e1"], sf, args)

# ---------- D. -max_hsps ----------
for k, (kind, n, mh, mts, fmt) in enumerate((("multi_hsp", 600, 1, None, "6"), ("multi_hsp", 600, 2, None, "6"), ("multi_hsp", 600, 3, None, "6"),
                                          ("both_strands", 700, 1, None, "6"), ("multi_hsp", 60, 2, 3, "6"), ("multi_hsp", 60, 1, 1, "6"),
                                          ("indel", 600, 2, None, "7"), ("multi_hsp", 700, 3, 5, "6"), ("both_strands", 80, 1, 2, "0"),
                                          ("multi_hsp", 560, 2, None, "0"), ("classes", 600, 1, None, "6"), ("multi_hsp", 1000, 3, 100, "6"))):
    rng = random.Random(f"D-{k}")
    sf = write(f"{W}/D_{k}.fa", subjects(kind, qs["e2"], rng, n, lens=(60, 80, 100, 150, 250, 400)))
    args = ["-task", "blastn", "-outfmt", fmt, "-max_hsps", str(mh)] + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"D_maxhsps{mh}_{k}_{kind}_{n}_{mts}", "D_max_hsps", qf["e2"], sf, args)

# ---------- E. -num_threads 4 ----------
for k, (kind, n, mts, task) in enumerate((("classes", 600, None, "blastn"), ("ties_exact", 800, None, "blastn"), ("multi_hsp", 60, 3, "blastn"),
                                       ("dups", 700, None, "megablast"), ("mixed", 900, 5, "blastn"), ("rc_mix", 40, 1, "megablast"),
                                       ("both_strands", 1000, None, "blastn"), ("classes", 30, 2, "blastn"), ("same_window", 700, None, "blastn"),
                                       ("indel", 560, 10, "blastn"))):
    rng = random.Random(f"E-{k}")
    lens = (60, 80, 100, 150) if task == "blastn" else (45, 60, 100, 150)
    if kind == "ties_exact": lens = (60,)
    sf = write(f"{W}/E_{k}.fa", subjects(kind, qs["s1"], rng, n, lens=lens, flank="fixed" if kind=="ties_exact" else "random"))
    args = ["-task", task, "-outfmt", "6", "-num_threads", "4"] + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"E_threads_{k}_{kind}_{n}_{mts}_{task}", "E_num_threads", qf["s1"], sf, args)

# ---------- F. -evalue tiny / large ----------
for k, (kind, n, ev, mts) in enumerate((("mixed", 700, "1e-30", None), ("classes", 700, "1e-25", None), ("mixed", 700, "1e-60", None),
                                      ("classes", 800, "1000", None), ("classes", 800, "100000", None), ("ties_exact", 700, "10000", None),
                                      ("classes", 60, "1e-20", 3), ("mixed", 60, "1000", 4), ("dups", 700, "1e-10", None), ("indel", 900, "50", 5))):
    rng = random.Random(f"F-{k}")
    lens = (30, 40, 60, 100, 150, 250)
    if kind == "ties_exact": lens = (30,)
    sf = write(f"{W}/F_{k}.fa", subjects(kind, qs["e3"], rng, n, lens=lens, flank="fixed" if kind=="ties_exact" else "random"))
    args = ["-task", "blastn", "-outfmt", "6", "-evalue", ev] + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"F_evalue{ev}_{k}_{kind}_{n}_{mts}", "F_evalue", qf["e3"], sf, args)

# ---------- G. query split ----------
q2m = edl[:2_000_000]
write(f"{W}/q2m.fa", [("EDL933_2M", q2m)])
ch = chunks([2_000_000], 1_000_000)
bnd = ch[1][0]  # start of chunk 2
bend = ch[0][1]  # end of chunk 1
def spans(rng, n):
    out = []
    for i in range(n):
        L = rng.choice([300, 500, 800, 1500, 3000])
        a = bnd - rng.randint(50, L - 50)
        out.append((f"X{i}", subst(rng, q2m[a:a + L], rng.choice([0, 0.02]))))
    return out
for k, (nin, nspan, mts, extra) in enumerate((
        (560, 3, None, []), (560, 1, 1, []), (560, 20, 3, []), (100, 15, 5, []), (30, 3, 10, []), (30, 3, 11, []), (600, 5, None, []),
        (700, 40, None, ["-outfmt", "7"]), (80, 12, 4, []), (560, 3, None, ["-subject_besthit"]), (560, 8, 2, ["-subject_besthit"]),
        (560, 10, None, ["-max_hsps", "1"]), (120, 20, 3, ["-max_hsps", "2"]), (560, 3, None, ["-num_threads", "4"]), (60, 5, 2, ["-num_threads", "4"]),
        (900, 10, None, ["-outfmt", "0"]), (50, 4, 1, ["-outfmt", "0"]), (560, 30, 6, []))):
    rng = random.Random(f"G-{k}")
    # inside one chunk only: windows within first chunk (far from the boundary) and second chunk
    recs = []
    for i in range(nin):
        L = rng.choice([60, 80, 100, 150, 250, 400])
        lo, hi = (0, bnd - 5000) if i % 2 == 0 else (bnd + 5000, 2_000_000)
        a = rng.randrange(lo, hi - L)
        w = subst(rng, q2m[a:a + L], rng.choice([0, 0, 0.03]))
        if rng.random() < 0.3: w = rc(w)
        recs.append((f"I{i}", rnd(rng, rng.randint(0, 100)) + w + rnd(rng, rng.randint(0, 100))))
    recs += spans(rng, nspan)
    rng.shuffle(recs)
    recs = [(f"R{j}", s) for j, (_, s) in enumerate(recs)]
    sf = write(f"{W}/G_{k}.fa", recs)
    fmt = ["-outfmt", "6"] if "-outfmt" not in extra else []
    args = ["-task", "blastn"] + fmt + extra + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"G_split_{k}_in{nin}_span{nspan}_mts{mts}", "G_split", "w/q2m.fa", sf, args)
# spanning only / boundary-region tie-heavy: identical windows covering the boundary, many copies
for k, (n, mts, extra) in enumerate(((600, None, []), (600, 1, []), (60, 3, []), (560, 5, ["-outfmt", "7"]), (700, None, ["-task", "megablast"]))):
    rng = random.Random(f"G2-{k}")
    L = 400
    a = bnd - 200
    w = q2m[a:a + L]
    recs = [(f"B{i}", w if i % 2 else rc(w)) for i in range(n)]
    sf = write(f"{W}/G2_{k}.fa", recs)
    task = ["-task", "blastn"] if "-task" not in extra else []
    fmt = ["-outfmt", "6"] if "-outfmt" not in extra else []
    args = task + fmt + extra + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"G2_splitbnd_{k}_n{n}_mts{mts}", "G_split", "w/q2m.fa", sf, args)
# megablast split: 11 Mb query
both = edl + sak
write(f"{W}/qcat.fa", [("EDL933_Sakai", both)])
chm = chunks([len(both)], 5_000_000)
bm = chm[1][0]
for k, (nin, nspan, mts) in enumerate(((560, 3, None), (40, 6, 1), (60, 6, 3), (580, 12, 5))):
    rng = random.Random(f"G3-{k}")
    recs = []
    for i in range(nin):
        L = rng.choice([45, 60, 100, 150, 250, 400])
        lo, hi = (0, bm - 5000) if i % 2 == 0 else (bm + 5000, len(both))
        a = rng.randrange(lo, hi - L)
        w = both[a:a + L]
        if rng.random() < 0.3: w = rc(w)
        recs.append((f"I{i}", rnd(rng, rng.randint(0, 100)) + w + rnd(rng, rng.randint(0, 100))))
    for i in range(nspan):
        L = rng.choice([300, 800, 3000])
        a = bm - rng.randint(50, L - 50)
        recs.append((f"X{i}", both[a:a + L]))
    rng.shuffle(recs)
    sf = write(f"{W}/G3_{k}.fa", [(f"R{j}", s) for j, (_, s) in enumerate(recs)])
    args = ["-task", "megablast", "-outfmt", "6"] + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"G3_splitmega_{k}_in{nin}_span{nspan}_mts{mts}", "G_split", "w/qcat.fa", sf, args)
# whole EDL933 genome (5 chunks, blastn)
write(f"{W}/qedl.fa", [("EDL933", edl)])
che = chunks([len(edl)], 1_000_000)
for k, (nin, nspan, mts) in enumerate(((600, 8, None), (60, 8, 2))):
    rng = random.Random(f"G4-{k}")
    recs = []
    for i in range(nin):
        L = rng.choice([60, 80, 100, 150, 250])
        a = rng.randrange(0, len(edl) - L)
        w = edl[a:a + L]
        if rng.random() < 0.3: w = rc(w)
        recs.append((f"I{i}", w))
    for i in range(nspan):
        p = che[rng.randint(1, len(che) - 1)][0]
        L = rng.choice([300, 800, 3000]); a = p - rng.randint(50, L - 50)
        recs.append((f"X{i}", edl[a:a + L]))
    rng.shuffle(recs)
    sf = write(f"{W}/G4_{k}.fa", [(f"R{j}", s) for j, (_, s) in enumerate(recs)])
    args = ["-task", "blastn", "-outfmt", "6"] + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"G4_wholeEDL_{k}_in{nin}_span{nspan}_mts{mts}", "G_split", "w/qedl.fa", sf, args)

# ---------- H. multi-query ----------
def multiq(tag, seqs):
    return write(f"{W}/H_{tag}_q.fa", [(f"Q{i}", s) for i, s in enumerate(seqs)])
for k, (nq, qlen, mix, mts, extra, task) in enumerate((
        (4, 1000, (700, 30, 0, 600), None, [], "blastn"), (4, 1000, (700, 30, 0, 600), 1, [], "blastn"), (3, 1500, (650, 0, 40), 3, [], "blastn"),
        (6, 3000, (600, 20, 600, 0, 15, 580), None, [], "blastn"), (6, 3000, (600, 20, 600, 0, 15, 580), 5, ["-outfmt", "7"], "blastn"),
        (5, 800, (560, 560, 30, 12, 0), 2, [], "megablast"), (3, 1200, (900, 20, 700), None, ["-outfmt", "0"], "blastn"),
        (4, 1000, (700, 30, 5, 600), None, ["-subject_besthit"], "blastn"), (4, 1000, (700, 30, 5, 600), 4, ["-max_hsps", "1"], "blastn"),
        (3, 1500, (650, 600, 40), None, ["-num_threads", "4"], "blastn"), (8, 2500, (600, 5, 600, 0, 12, 3, 570, 0), 10, [], "blastn"),
        (2, 6000, (700, 600), 3, [], "blastn"), (3, 4000, (580, 600, 20), None, [], "megablast"))):
    rng = random.Random(f"H-{k}")
    starts = [rng.randrange(0, len(edl) - qlen) for _ in range(nq)]
    qseqs = [edl[a:a + qlen] for a in starts]
    qf_ = multiq(str(k), qseqs)
    recs = []
    for qi, cnt in enumerate(mix):
        lens = (60, 80, 100, 150, 250) if task == "blastn" else (45, 60, 100, 150, 250)
        for r in subjects("classes" if (qi + k) % 2 == 0 else "rc_mix", qseqs[qi], rng, cnt, lens=lens):
            recs.append((f"q{qi}_{r[0]}", r[1]))
    # some shared subjects hitting two queries
    for j in range(5):
        a, b = rng.sample(qseqs, 2) if nq >= 2 else (qseqs[0], qseqs[0])
        recs.append((f"shared{j}", a[100:260] + rnd(rng, 30) + rc(b[300:420])))
    rng.shuffle(recs)
    sf = write(f"{W}/H_{k}.fa", [(f"R{j}", s) for j, (_, s) in enumerate(recs)])
    args = ["-task", task] + (["-outfmt", "6"] if "-outfmt" not in extra else []) + extra + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"H_multiq_{k}_nq{nq}_mts{mts}_{task}", "H_multi_query", qf_, sf, args)

# ---------- I. scoring variants (best HSP != first HSP candidates) ----------
for k, sc in enumerate((["-reward", "1", "-penalty", "-3"], ["-reward", "1", "-penalty", "-1", "-gapopen", "5", "-gapextend", "2"],
                        ["-reward", "4", "-penalty", "-5"], ["-reward", "2", "-penalty", "-3", "-gapopen", "5", "-gapextend", "2"],
                        ["-word_size", "7"], ["-dust", "no"], ["-xdrop_gap", "10"], ["-ungapped"] )):
    rng = random.Random(f"I-{k}")
    sf = write(f"{W}/I_{k}.fa", subjects("multi_hsp", qs["e1"], rng, 600, lens=(60, 80, 100, 150, 250, 400)))
    add(f"I_scoring_{k}", "I_scoring", qf["e1"], sf, ["-task", "blastn", "-outfmt", "6"] + sc)
    rng = random.Random(f"I2-{k}")
    sf2 = write(f"{W}/I2_{k}.fa", subjects("both_strands", qs["e1"], rng, 40, lens=(60, 100, 150, 250, 400)))
    add(f"I2_scoring_{k}_mts2", "I_scoring", qf["e1"], sf2, ["-task", "blastn", "-outfmt", "6", "-max_target_seqs", "2"] + sc)

# ---------- K. e-values below 1e-180 compare equal in NCBI's heap (epsilon rule) ----------
for k, (n, mts, task, extra) in enumerate(((560, None, "blastn", []), (700, None, "blastn", []), (600, None, "megablast", []), (40, 1, "blastn", []),
                                        (40, 3, "blastn", []), (60, 5, "megablast", []), (800, 100, "blastn", []), (580, None, "blastn", ["-outfmt", "7"]),
                                        (30, 2, "blastn", ["-subject_besthit"]), (620, None, "blastn", ["-max_hsps", "1"]))):
    rng = random.Random(f"K-{k}")
    lens = (200, 250, 300, 340, 360, 380, 400, 450, 600, 900, 1500)
    sf = write(f"{W}/K_{k}.fa", subjects("rc_mix" if k % 2 else "classes", qs["e1"], rng, n, lens=lens))
    args = ["-task", task] + (["-outfmt", "6"] if "-outfmt" not in extra else []) + extra + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"K_eps_{k}_{n}_mts{mts}_{task}", "K_epsilon", qf["e1"], sf, args)
# ---------- J. hits near the e-value cutoff (prelim cut-off vs traceback) ----------
for k, (n, mts, ws, ev, kind) in enumerate(((700, None, 7, None, "classes"), (700, None, 7, None, "mixed"), (900, None, 8, None, "rc_mix"), (600, None, 7, "1", "classes"),
                                         (700, 3, 7, None, "mixed"), (50, 2, 7, None, "classes"), (800, None, 7, "100", "indel"), (700, 1, 7, "0.1", "mixed"),
                                         (1000, None, 11, "10", "classes"), (600, 5, 7, None, "rc_mix"))):
    rng = random.Random(f"J-{k}")
    lens = (14, 16, 18, 20, 22, 24, 26, 30)
    sf = write(f"{W}/J_{k}.fa", subjects(kind, qs["e3"], rng, n, lens=lens))
    args = ["-task", "blastn", "-outfmt", "6", "-word_size", str(ws)] + (["-evalue", ev] if ev else []) + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"J_cutoff_{k}_{n}_ws{ws}_ev{ev}_mts{mts}", "J_near_cutoff", qf["e3"], sf, args)
# ---------- G5. whole genome query vs 600 fragments of the other genome ----------
for k, (nfrag, mts, task, qk, extra) in enumerate(((600, None, "blastn", "qedl", []), (600, 3, "blastn", "qedl", []), (700, None, "megablast", "qedl", []),
                                                  (600, 1, "megablast", "qedl", []), (560, None, "megablast", "qedl", ["-subject_besthit"]),
                                                  (580, 5, "blastn", "qedl", ["-max_hsps", "2"]))):
    rng = random.Random(f"G5-{k}")
    recs = []
    for i in range(nfrag):
        L = rng.choice([800, 2000, 5000, 10000, 20000]); a = rng.randrange(0, len(sak) - L)
        w = sak[a:a + L]
        if rng.random() < 0.4: w = rc(w)
        recs.append((f"F{i}", w))
    sf = write(f"{W}/G5_{k}.fa", recs)
    args = ["-task", task, "-outfmt", "6"] + extra + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"G5_genome_{k}_{nfrag}_mts{mts}_{task}", "G_split", f"{W}/{qk}.fa", sf, args)

# ---------- P. subjects whose preliminary score is lower than their final score (gap of 20-40 letters beyond the
# preliminary X-drop): the preliminary rank differs from the final rank ----------
def gapsubj(rng, q, n, strong_frac=0.3):
    recs = []
    for i in range(n):
        r = rng.random()
        a = rng.randrange(0, len(q) - 8000)
        if r < strong_frac:           # clean long exact hit (strong in prelim)
            s = q[a:a + rng.choice([1500, 2500, 3500])]
        elif r < strong_frac + 0.35:  # insertion in the subject
            h1, h2 = rng.choice([800, 1200, 1800]), rng.choice([800, 1200, 1800])
            s = q[a:a + h1] + rnd(rng, rng.randint(20, 40)) + q[a + h1:a + h1 + h2]
        else:                         # deletion in the subject
            h1, h2 = rng.choice([800, 1200, 1800]), rng.choice([800, 1200, 1800])
            g = rng.randint(20, 40)
            s = q[a:a + h1] + q[a + h1 + g:a + h1 + g + h2]
        if rng.random() < 0.4: s = rc(s)
        recs.append((f"P{i}", s))
    return recs
for k, (n, mts, task, extra, frac) in enumerate(((12, 1, "blastn", [], 0.5), (12, 1, "blastn", [], 0.9), (14, 2, "blastn", [], 0.6), (20, 3, "blastn", [], 0.5),
        (30, 5, "blastn", [], 0.4), (15, 1, "megablast", [], 0.6), (25, 4, "megablast", [], 0.5), (600, None, "blastn", [], 0.3), (700, 100, "blastn", [], 0.2),
        (570, None, "megablast", [], 0.1), (16, 1, "blastn", ["-subject_besthit"], 0.6), (20, 3, "blastn", ["-max_hsps", "1"], 0.5),
        (20, 2, "blastn", ["-num_threads", "4"], 0.5), (14, 1, "blastn", ["-outfmt", "0"], 0.6), (18, 3, "blastn", ["-outfmt", "7"], 0.5),
        (13, 1, "blastn", ["-evalue", "1e-100"], 0.6), (20, 10, "blastn", [], 0.5), (30, 11, "blastn", [], 0.4))):
    rng = random.Random(f"P-{k}")
    sf = write(f"{W}/P_{k}.fa", gapsubj(rng, qs["e2"], n, frac))
    args = ["-task", task] + (["-outfmt", "6"] if "-outfmt" not in extra else []) + extra + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"P_gap_{k}_{n}_mts{mts}_{task}", "P_prelim_vs_final", qf["e2"], sf, args)
# multi-query with the same kind of subjects (each query has its own hit subjects)
for k, (nq, per, mts) in enumerate(((3, (14, 14, 0), 1), (4, (20, 600, 15, 20), None), (4, (20, 600, 15, 20), 3), (2, (570, 12), 2))):
    rng = random.Random(f"P2-{k}")
    qseqs = [edl[a:a + 20000] for a in [rng.randrange(0, len(edl) - 20000) for _ in range(nq)]]
    qf_ = write(f"{W}/P2_{k}_q.fa", [(f"Q{i}", s) for i, s in enumerate(qseqs)])
    recs = []
    for qi, cnt in enumerate(per):
        for r in gapsubj(rng, qseqs[qi], cnt, 0.4): recs.append((f"q{qi}_{r[0]}", r[1]))
    rng.shuffle(recs)
    sf = write(f"{W}/P2_{k}.fa", recs)
    add(f"P2_gapmultiq_{k}_nq{nq}_mts{mts}", "P_prelim_vs_final", qf_, sf, ["-task", "blastn", "-outfmt", "6"] + (["-max_target_seqs", str(mts)] if mts else []))
for k, (n, mts, task, frac) in enumerate(((40, 1, "blastn", 0.5), (60, 1, "blastn", 0.6), (60, 2, "blastn", 0.6), (100, 3, "blastn", 0.7), (100, 5, "blastn", 0.6),
                                          (40, 1, "megablast", 0.5), (80, 2, "megablast", 0.6), (100, 4, "blastn", 0.5), (50, 1, "blastn", 0.7), (70, 3, "megablast", 0.6))):
    rng = random.Random(f"P3-{k}")
    sf = write(f"{W}/P3_{k}.fa", gapsubj(rng, qs["e2"], n, frac))
    add(f"P3_gap_{k}_{n}_mts{mts}_{task}", "P_prelim_vs_final", qf["e2"], sf, ["-task", task, "-outfmt", "6", "-max_target_seqs", str(mts)])
# Q. duplicated gap/clean subjects: prelim ties decided by the oid rule while the final rank differs
for k, (ntypes, copies, mts, task, extra) in enumerate(((6, 10, 1, "blastn", []), (6, 20, 3, "blastn", []), (8, 10, 5, "blastn", []), (5, 120, None, "blastn", []),
                                                    (4, 160, 100, "blastn", []), (6, 12, 2, "megablast", []), (7, 90, None, "megablast", ["-outfmt", "7"]),
                                                    (6, 15, 4, "blastn", ["-outfmt", "0"]), (6, 100, None, "blastn", ["-max_hsps", "1"]), (5, 130, 6, "blastn", ["-subject_besthit"]))):
    rng = random.Random(f"Q-{k}")
    base = [r[1] for r in gapsubj(rng, qs["e2"], ntypes, 0.4)]
    recs = []
    for c in range(copies):
        for t, s in enumerate(base):
            recs.append((f"T{t}_{c}", s))
    rng.shuffle(recs)
    sf = write(f"{W}/Q_{k}.fa", recs)
    args = ["-task", task] + (["-outfmt", "6"] if "-outfmt" not in extra else []) + extra + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"Q_dupgap_{k}_{ntypes}x{copies}_mts{mts}_{task}", "Q_dup_gap", qf["e2"], sf, args)
# R. gap subjects (prelim rank != final rank) under besthit / threads / mts 6,10,11 / multi-query
for k, (n, mts, task, extra, frac) in enumerate(((60, 1, "blastn", ["-subject_besthit"], 0.6), (80, 3, "blastn", ["-subject_besthit"], 0.5), (600, None, "blastn", ["-subject_besthit"], 0.3),
        (60, 2, "megablast", ["-subject_besthit"], 0.6), (60, 1, "blastn", ["-num_threads", "4"], 0.6), (80, 5, "blastn", ["-num_threads", "4"], 0.6),
        (600, None, "blastn", ["-num_threads", "4"], 0.3), (60, 3, "megablast", ["-num_threads", "4"], 0.6), (60, 6, "blastn", [], 0.6), (60, 10, "blastn", [], 0.5),
        (60, 11, "blastn", [], 0.5), (60, 1, "blastn", ["-max_hsps", "1"], 0.6), (60, 2, "blastn", ["-max_hsps", "2"], 0.6), (80, 3, "blastn", ["-max_hsps", "3"], 0.6),
        (60, 3, "blastn", ["-evalue", "1e-300"], 0.6), (70, 2, "blastn", ["-evalue", "1e5"], 0.6))):
    rng = random.Random(f"R-{k}")
    sf = write(f"{W}/R_{k}.fa", gapsubj(rng, qs["e1"], n, frac))
    args = ["-task", task, "-outfmt", "6"] + extra + (["-max_target_seqs", str(mts)] if mts else [])
    add(f"R_gap_{k}_{n}_mts{mts}_{task}", "R_gap_options", qf["e1"], sf, args)
for k, (per, mts, extra, task) in enumerate((((40, 40, 5), 1, ["-subject_besthit"], "blastn"), ((60, 600, 8, 60), 3, ["-num_threads", "4"], "blastn"), ((60, 0, 70), 2, ["-max_hsps", "1"], "blastn"),
        ((50, 580), 5, [], "megablast"), ((50, 50, 50, 50, 50, 50), 1, [], "blastn"), ((80, 80, 80), 4, ["-outfmt", "7"], "blastn"))):
    rng = random.Random(f"R2-{k}")
    nq = len(per)
    qseqs = [edl[a:a + 20000] for a in [rng.randrange(0, len(edl) - 20000) for _ in range(nq)]]
    qf_ = write(f"{W}/R2_{k}_q.fa", [(f"Q{i}", s) for i, s in enumerate(qseqs)])
    recs = []
    for qi, cnt in enumerate(per):
        for r in gapsubj(rng, qseqs[qi], cnt, 0.4): recs.append((f"q{qi}_{r[0]}", r[1]))
    rng.shuffle(recs)
    sf = write(f"{W}/R2_{k}.fa", recs)
    args = ["-task", task] + (["-outfmt", "6"] if "-outfmt" not in extra else []) + extra + ["-max_target_seqs", str(mts)]
    add(f"R2_gapmultiq_{k}_nq{nq}_mts{mts}_{task}", "R_gap_options", qf_, sf, args)
# S. split query batches with gap subjects (prelim score != final score), subjects anywhere and across the boundary
def gapsubj_at(rng, q, n, centers, frac=0.5):
    recs = gapsubj(rng, q, n, frac)
    out = []
    for i, (nm, s) in enumerate(recs):
        if centers and i % 3 == 0:     # place around a chunk boundary: take the window centred there
            c = rng.choice(centers)
            h1 = rng.choice([800, 1200, 1800]); h2 = rng.choice([800, 1200, 1800]); g = rng.randint(20, 40)
            a = c - h1 + rng.randint(-100, 100)
            if i % 2: s = q[a:a + h1] + rnd(rng, g) + q[a + h1:a + h1 + h2]
            else: s = q[a:a + h1] + q[a + h1 + g:a + h1 + g + h2]
        out.append((nm, s))
    return out
for k, (n, mts, task, extra, frac) in enumerate(((40, 1, "blastn", [], 0.5), (60, 2, "blastn", [], 0.6), (80, 3, "blastn", [], 0.5), (100, 5, "blastn", [], 0.6),
        (600, None, "blastn", [], 0.3), (60, 1, "blastn", ["-subject_besthit"], 0.6), (60, 3, "blastn", ["-num_threads", "4"], 0.6), (60, 2, "blastn", ["-outfmt", "7"], 0.6),
        (60, 1, "blastn", ["-max_hsps", "1"], 0.6), (60, 4, "blastn", ["-outfmt", "0"], 0.6))):
    rng = random.Random(f"S-{k}")
    sf = write(f"{W}/S_{k}.fa", gapsubj_at(rng, q2m, n, [bnd, bend], frac))
    args = ["-task", task] + (["-outfmt", "6"] if "-outfmt" not in extra else []) + extra + ["-max_target_seqs", str(mts)] * (mts is not None)
    add(f"S_splitgap_{k}_{n}_mts{mts}_{task}", "S_split_gap", "w/q2m.fa", sf, args)
for k, (n, mts) in enumerate(((60, 1), (60, 3), (80, 5), (600, None))):
    rng = random.Random(f"S2-{k}")
    sf = write(f"{W}/S2_{k}.fa", gapsubj_at(rng, both, n, [bm, len(edl)], 0.5))
    add(f"S2_splitgapmega_{k}_{n}_mts{mts}", "S_split_gap", "w/qcat.fa", sf, ["-task", "megablast", "-outfmt", "6"] + ["-max_target_seqs", str(mts)] * (mts is not None))
# split batch followed/preceded by short queries (multi-query, a split batch plus others)
for k, (order, mts) in enumerate((("big_first", 2), ("big_last", 3), ("big_mid", None), ("big_first", 1))):
    rng = random.Random(f"S3-{k}")
    shorts = [edl[a:a + 20000] for a in [rng.randrange(2_100_000, len(edl) - 20000) for _ in range(3)]]
    qrecs = {"big_first": [q2m] + shorts, "big_last": shorts + [q2m], "big_mid": shorts[:1] + [q2m] + shorts[1:]}[order]
    qf_ = write(f"{W}/S3_{k}_q.fa", [(f"Q{i}", s) for i, s in enumerate(qrecs)])
    recs = []
    for qi, s in enumerate(qrecs):
        cnt = (60, 40, 600, 30)[qi % 4]
        centers = [bnd] if s is q2m else []
        for r in gapsubj_at(rng, s, cnt, centers, 0.5): recs.append((f"q{qi}_{r[0]}", r[1]))
    rng.shuffle(recs)
    sf = write(f"{W}/S3_{k}.fa", recs)
    add(f"S3_splitmultiq_{k}_{order}_mts{mts}", "S_split_gap", qf_, sf, ["-task", "blastn", "-outfmt", "6"] + ["-max_target_seqs", str(mts)] * (mts is not None))
# T. per-chunk overflow of the preliminary list (more than 550 subjects hit in EACH chunk) and the merge of the overflowed lists
def chunk_windows(rng, q, ranges, per, lens, mode, nspan=0, bounds=()):
    recs = []
    base = [q[rng.randrange(r[0], r[1] - max(lens)):][:lens[0]] for r in ranges for _ in range(8)]
    for ri, (lo, hi) in enumerate(ranges):
        for i in range(per):
            L = rng.choice(lens)
            if mode == "dups":
                w = base[ri * 8 + i % 8]
            else:
                a = rng.randrange(lo, hi - L); w = q[a:a + L]
            if rng.random() < 0.3: w = rc(w)
            recs.append((f"C{ri}_{i}", w if mode != "exact_fix" else w[:lens[0]]))
    for i in range(nspan):
        b = rng.choice(bounds); L = rng.choice([300, 800, 1500]); a = b - rng.randint(50, L - 50)
        recs.append((f"X{i}", q[a:a + L]))
    rng.shuffle(recs)
    return [(f"R{j}", s) for j, (_, s) in enumerate(recs)]
for k, (per, mode, mts, extra) in enumerate(((620, "mixed", None, []), (700, "dups", None, []), (650, "exact_fix", None, []), (600, "mixed", 3, []), (600, "exact_fix", 1, []),
                                             (600, "dups", 5, []), (620, "mixed", None, ["-outfmt", "7"]), (600, "mixed", 2, ["-subject_besthit"]), (580, "mixed", None, ["-max_hsps", "1"]))):
    rng = random.Random(f"T-{k}")
    lens = (60,) if mode == "exact_fix" else (50, 60, 80, 100, 150, 250)
    recs = chunk_windows(rng, q2m, [(0, bnd - 3000), (bnd + 3000, 2_000_000)], per, lens, mode, nspan=6, bounds=[bnd])
    sf = write(f"{W}/T_{k}.fa", recs)
    args = ["-task", "blastn"] + (["-outfmt", "6"] if "-outfmt" not in extra else []) + extra + ["-max_target_seqs", str(mts)] * (mts is not None)
    add(f"T_chunkover_{k}_{per}_{mode}_mts{mts}", "T_chunk_overflow", "w/q2m.fa", sf, args)
for k, (per, mode, mts) in enumerate(((620, "mixed", None), (600, "dups", 3), (600, "exact_fix", 1))):
    rng = random.Random(f"T2-{k}")
    lens = (45,) if mode == "exact_fix" else (45, 60, 100, 150, 250)
    recs = chunk_windows(rng, both, [(0, bm - 3000), (bm + 3000, len(both))], per, lens, mode, nspan=6, bounds=[bm])
    sf = write(f"{W}/T2_{k}.fa", recs)
    add(f"T2_chunkovermega_{k}_{per}_{mode}_mts{mts}", "T_chunk_overflow", "w/qcat.fa", sf, ["-task", "megablast", "-outfmt", "6"] + ["-max_target_seqs", str(mts)] * (mts is not None))
for k, (mode, mts) in enumerate((("mixed", None), ("exact_fix", 2))):
    rng = random.Random(f"T3-{k}")
    rngs = [(a, b) for a, b in [(c[0] + 3000, c[1] - 3000) for c in che]]
    lens = (60,) if mode == "exact_fix" else (50, 60, 80, 100, 150)
    recs = chunk_windows(rng, edl, rngs, 600, lens, mode, nspan=10, bounds=[c[0] for c in che[1:]])
    sf = write(f"{W}/T3_{k}.fa", recs)
    add(f"T3_chunkover5_{mode}_mts{mts}", "T_chunk_overflow", "w/qedl.fa", sf, ["-task", "blastn", "-outfmt", "6"] + ["-max_target_seqs", str(mts)] * (mts is not None))
# U. a whole-genome (chunked) subject together with many small subjects
for k, (n, mts, task) in enumerate(((600, None, "blastn"), (60, 3, "blastn"), (40, 1, "megablast"), (600, 100, "megablast"))):
    rng = random.Random(f"U-{k}")
    recs = [("BIG_sakai", sak)] + gapsubj(rng, qs["e1"], n, 0.4)
    rng.shuffle(recs)
    sf = write(f"{W}/U_{k}.fa", recs)
    add(f"U_bigsubject_{k}_{n}_mts{mts}_{task}", "U_big_subject", qf["e1"], sf, ["-task", task, "-outfmt", "6"] + ["-max_target_seqs", str(mts)] * (mts is not None))
# V. linear gap costs with greedy extension (megablast), a different Karlin block
for k, (sc, n, mts) in enumerate((("-reward 2 -penalty -3", 600, None), ("-reward 1 -penalty -2", 600, None), ("-reward 2 -penalty -7", 600, None),
                                  ("-reward 2 -penalty -3", 60, 2), ("-reward 1 -penalty -2", 60, 1), ("-reward 2 -penalty -7", 60, 3))):
    rng = random.Random(f"V-{k}")
    sf = write(f"{W}/V_{k}.fa", gapsubj(rng, qs["e1"], n, 0.4))
    add(f"V_linear_{k}_{n}_mts{mts}", "V_linear_gap", qf["e1"], sf, ["-task", "megablast", "-outfmt", "6", "-gapopen", "0", "-gapextend", "0", *sc.split()] + ["-max_target_seqs", str(mts)] * (mts is not None))
# W. -perc_identity (post-traceback drops), -dust no, other accepted options, with gap/mutated subjects
def mutsubj(rng, q, n):
    recs = []
    for i in range(n):
        a = rng.randrange(0, len(q) - 3000); L = rng.choice([300, 600, 1200, 2500])
        w = indel(rng, subst(rng, q[a:a + L], rng.choice([0.0, 0.03, 0.06, 0.1, 0.15])), rng.choice([0.0, 0.01, 0.02]))
        if rng.random() < 0.4: w = rc(w)
        recs.append((f"M{i}", rnd(rng, rng.randint(0, 100)) + w + rnd(rng, rng.randint(0, 100))))
    return recs
from gen import indel
for k, (opt, n, mts) in enumerate((("-perc_identity 99", 600, None), ("-perc_identity 90", 700, None), ("-perc_identity 80", 600, 100), ("-perc_identity 95", 60, 1),
                                   ("-perc_identity 90", 80, 2), ("-perc_identity 97", 80, 3), ("-perc_identity 85", 100, 5), ("-dust no", 600, None),
                                   ("-dust no", 60, 2), ("-gapopen 3 -gapextend 3", 600, None), ("-gapopen 3 -gapextend 3", 60, 3), ("-word_size 16", 600, None),
                                   ("-lcase_masking", 600, None), ("-perc_identity 92 -subject_besthit", 60, 2), ("-perc_identity 92 -max_hsps 1", 600, None))):
    rng = random.Random(f"W-{k}")
    sf = write(f"{W}/W_{k}.fa", mutsubj(rng, qs["e2"], n) + gapsubj(rng, qs["e2"], n // 4, 0.4))
    add(f"W_opt_{k}_{opt.replace(' ','_')}_{n}_mts{mts}", "W_other_options", qf["e2"], sf, ["-task", "blastn", "-outfmt", "6", *opt.split()] + ["-max_target_seqs", str(mts)] * (mts is not None))
