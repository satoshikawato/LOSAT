#!/usr/bin/env python3
"""Series NP, M4b PR 1, review H1 (2026-11-21): the inputs pbh of dbmaskp/ (README.md there). Finds blastx
-db_hard_mask -subject_besthit inputs where the filter after every chunk (NCBI, blast_engine.c:586-591) differs from
one filter after the merge (LOSAT before the fix). Subjects: [A][C] | hard mask | [B], where B (a close copy of a query
stretch) holds A (a weaker copy of nearly the same stretch) and A holds C (a weaker copy starting a little before A).
Draws rounds of 400 candidates (seed 20261121) and keeps the first three subjects that differ (and with
--ncbi-only the first three that differ in blastp with the protein of the query); writes them as
pbh.faa / edge_pbh.asnt / qbh.fna / qbh.faa (the protein the query back-translates) and the database pbh into OUT.

usage: make_dbmaskp_pbh.py WORK NCBI_BIN OLD_LOSAT NEW_LOSAT OUT
       make_dbmaskp_pbh.py --ncbi-only WORK NCBI_BIN OUT
The first form runs LOSAT before the fix (OLD: one filter after the merge) and after it (NEW: a filter after every
chunk) and keeps the subjects whose rows differ. The second (resume 2 of PR 1, 2026-11-21, while no LOSAT build was at
hand) runs NCBI BLAST+ alone: NCBI's rows with -subject_besthit against its rows without the option filtered once in
the order of their scores (the filter of blast_hits.c:2537-2580 on the query ranges of each frame, widened by 3
residues), and keeps the subjects where NCBI's filter deleted an alignment that the single filter keeps. Both forms
draw the same candidates.
"""
import os
import random
import subprocess
import sys

NCBI_ONLY = sys.argv[1:2] == ["--ncbi-only"]
if NCBI_ONLY:
    WORK, NB, OUT = sys.argv[2:5]
    OLD = NEW = None
else:
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
FMT_FRAME = "6 qseqid sseqid qstart qend sstart send score qframe"
ENV = dict(os.environ, BLAST_USAGE_REPORT="false")
# DEFAULT_SUBJECT_BESTHIT_PROT_MAX_RANGE_DIFF (blast_options.h): the widening of the query range.
RANGE_DIFF = 3


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


def rows(binary, d, db, query, ncbi=False, besthit=True, fmt=FMT, program="blastx"):
    cmd = ([f"{NB}/{program}"] if ncbi else [binary, program]) + [
        "-db", db, "-query", query, "-db_hard_mask", "21", "-outfmt", fmt, "-num_threads", "1"]
    if besthit:
        cmd.append("-subject_besthit")
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


def frame_range(line, query_length):
    """The frame and the query range in residues of the frame of a row of FMT_FRAME (blastx: positions in the
    nucleotide query; blastp: residues, frame 0)."""
    cells = line.split("\t")
    qstart, qend, frame = int(cells[2]), int(cells[3]), int(cells[7])
    if frame == 0:
        return frame, qstart - 1, qend
    shift = abs(frame) - 1
    if frame > 0:
        return frame, (qstart - 1 - shift) // 3, (qend - shift) // 3
    return frame, (query_length - qstart - shift) // 3, (query_length - qend + 1 - shift) // 3


def filtered_once(lines, query_length):
    """`lines` (one subject's rows) after one pass of the subject best-hit filter in the order of their scores."""
    order = sorted(range(len(lines)), key=lambda k: -int(lines[k].split("\t")[6]))
    hsps = [(lines[k], *frame_range(lines[k], query_length)) for k in order]
    kept = [True] * len(hsps)
    for i, (_, frame, offset, end) in enumerate(hsps):
        if not kept[i]:
            continue
        o, e = max(offset - RANGE_DIFF, 0), end + RANGE_DIFF
        for j in range(i + 1, len(hsps)):
            if kept[j] and hsps[j][1] == frame and hsps[j][2] >= o and hsps[j][3] <= e:
                kept[j] = False
    return sorted(h[0] for h, k in zip(hsps, kept) if k)


def ncbi_differing(d, db, query, query_length, program="blastx"):
    """The subjects where NCBI's rows with -subject_besthit hold fewer alignments than its rows without the option
    filtered once, in the order of the subjects."""
    _, plain = rows(None, d, db, query, ncbi=True, besthit=False, fmt=FMT_FRAME, program=program)
    _, best = rows(None, d, db, query, ncbi=True, besthit=True, fmt=FMT_FRAME, program=program)
    diff = []
    for s in sorted(set(plain), key=subject_index):
        once = filtered_once(plain[s], query_length)
        if len(best.get(s, [])) < len(once):
            diff.append(s)
    return diff, plain, best


def remove_fixture():
    for ext in ("faa", "pdb", "phr", "pin", "pot", "psq", "ptf", "pto", "paa", "pab", "pac", "pjs"):
        if os.path.exists(f"{OUT}/pbh.{ext}"):
            os.remove(f"{OUT}/pbh.{ext}")
    os.remove(f"{OUT}/edge_pbh.asnt")


def main():
    os.makedirs(WORK, exist_ok=True)
    r = random.Random(20261121)
    q = rnd(r, 220)
    qdna = "".join(r.choice(CODONS[a]) for a in q)
    with open(f"{WORK}/qbh.fna", "w") as f:
        f.write(">qbh query with one domain, back-translated\n" + "\n".join(qdna[i:i + 60] for i in range(0, len(qdna), 60)) + "\n")
    with open(f"{WORK}/qbh.faa", "w") as f:
        f.write(">qbhp query with one domain\n" + "\n".join(q[i:i + 60] for i in range(0, len(q), 60)) + "\n")
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
        if NCBI_ONLY:
            diff, _, _ = ncbi_differing(d, "cand", "../qbh.fna", len(qdna))
            # blastp with the protein of the query: three more subjects where its filter deletes more.
            diff_p, _, _ = ncbi_differing(d, "cand", "../qbh.faa", len(q), program="blastp")
            print(f"round {rnd_no}: {len(seqs)} subjects, NCBI's filter deleted more on {len(diff)} (blastx), "
                  f"{len(diff_p)} (blastp)", flush=True)
            diff = diff[:3] + [s for s in diff_p if s not in diff[:3]][:3]
        else:
            _, old = rows(OLD, d, "cand", "../qbh.fna")
            _, new = rows(NEW, d, "cand", "../qbh.fna")
            diff = sorted((s for s in set(old) | set(new) if old.get(s) != new.get(s)), key=subject_index)
            print(f"round {rnd_no}: {len(seqs)} subjects, {len(old)} with rows, differing {len(diff)}", flush=True)
        if not diff:
            continue
        pick = [subject_index(s) for s in (diff if NCBI_ONLY else diff[:3])]
        # the fixture: the differing subjects alone, if they still differ there
        os.makedirs(OUT, exist_ok=True)
        build(OUT, "pbh", [seqs[i] for i in pick], [masks[i] for i in pick])
        for name in ("qbh.fna", "qbh.faa"):
            os.replace(f"{WORK}/{name}", f"{OUT}/{name}")
        c, _ = rows(None, OUT, "pbh", "qbh.fna", ncbi=True)
        if NCBI_ONLY:
            still, plain, best = ncbi_differing(OUT, "pbh", "qbh.fna", len(qdna))
            bp_still, _, _ = ncbi_differing(OUT, "pbh", "qbh.faa", len(q), program="blastp")
            found = bool(still) and (bool(bp_still) or len(pick) <= 3)
            print(f"fixture pbh: subjects {[i + 1 for i in pick]} blastx differing {still} blastp differing {bp_still}")
            for s in sorted(plain, key=subject_index):
                print(f"--- {s} without -subject_besthit\n" + "\n".join(plain[s]))
                print(f"--- {s} NCBI -subject_besthit\n" + "\n".join(best.get(s, [])))
                print(f"--- {s} filtered once\n" + "\n".join(filtered_once(plain[s], len(qdna))), flush=True)
        else:
            o, _ = rows(OLD, OUT, "pbh", "qbh.fna")
            n, _ = rows(NEW, OUT, "pbh", "qbh.fna")
            found = o != n
            print(f"fixture pbh: subjects {[i + 1 for i in pick]} old==new {o == n} ncbi==new {c == n} ncbi==old {c == o}")
            print("--- old\n" + o + "--- new\n" + n + "--- ncbi\n" + c, flush=True)
        if found:
            print("BH FOUND")
            return
        print("the picked subjects alone do not differ; next round")
        for name in ("qbh.fna", "qbh.faa"):
            os.replace(f"{OUT}/{name}", f"{WORK}/{name}")
        remove_fixture()
    print("BH NOT FOUND")


main()
