#!/usr/bin/env python3
"""Protein databases with overlapping hard masks (line NP, series M, M4b-3b; written by the fixture agent of ad-M4b-6,
2026-11-22): the shapes of the probes of M4b-3 (`ncbi-parity-20261019/M/m4b3/probe/make_p1.py`, seed 20261121) and
M4b-3b (`.../M/m4b3b/probe/make_p2.py`, seed 20261122), whose letters and masks this script reproduces draw for draw.
NCBI BLAST+ 2.17.0 `-db_hard_mask` searches the letters of an inverted range (a mask that starts before the previous one
ends) through a chunk of negative length: the rows `db.pmask.{bp,bx}.hard.ov_*` freeze what it writes.

Each shape is one database of one subject (`ov_lk`: two, A then B); the letters are random without W but for the W runs
named. The queries are random (no W) + WWW + random: `qov1.faa` / `qov1.fna` (the back-translation, for BLASTX) of probe
1 (e3, e3c, g3, nest, big), `qov2.faa` / `qov2.fna` of probe 2 (s2, tl, r37, r36, nb, lk).
  ov_e3:   masks [40,60) [55,70): inverted chunk at 60 (length -5); WWW at 60 (inside the second mask) and at 70; the
           word at 70 pairs with the stale diagonal entry of the word at 60 in the two-hit scan (a left-only seed).
  ov_e3c:  control, masks [40,60) [60,70) (touching): no inverted chunk, same letters.
  ov_g3:   masks [40,60) [55,80) [75,90): two consecutive inverted chunks (60 and 80, length -5), WWW at 60 and 80.
  ov_nest: masks [40,100) [50,60) (a mask inside the previous one): inverted chunk of length -50 (> the window 40);
           WWW at 100 and at 64.
  ov_big:  masks [40,120) [60,130): overlap 60 > the window; WWW at 120, 130 and 150.
  ov_s2:   masks [40,60) [55,70) [72,90): inverted chunk at 60, then a chunk of 2 letters (70,72) < the word 3; WWW at 60
           and 70 (a seed saved in a chunk shorter than the word).
  ov_r37:  masks [40,80) [43,90): overlap 37 = window - word; WWW at 80 and 90.
  ov_r36:  masks [40,80) [44,90): overlap 36; WWW at 80 and 90.
  ov_nb:   masks [40,60) [50,70) [65,80) [100,120) [110,130): consecutive inverted chunks with different words
           (WCW at 60, WMW at 70, WWW at 80).
  ov_tl:   s2's letters, masks [65,240) then [55,64): a mask that runs to the end with a mask stored after it (fix
           round 1 of M4b-3b, review L1; probe M/m4b3b/probe/fix1), which LOSAT refuses before the search.
  ov_wpe:  subjects A then B (probe3, fix round 1, review H2): A's masks [10,60) [20,99) [30,40), whose one word at 99
           covers A's sentinel and B's first letter; searched, also for `qwpe.faa` (`W-W`: the FASTA reader drops `-`).
  ov_lk:   subjects A then B: A's masks [40,200) [60,238) (lengths 40, -140, 2), the query's first 33 letters ending at
           A's chunk end 40; B holds the same letters and the mask [230,235).
Masks are `Blast-db-mask-info` ASN.1 text with `local id N` (N the 1-based order of the subject; no -parse_seqids, as the
`edge_*.asnt`), algorithm seg; `makeblastdb -mask_data` gives them the algorithm number 21.
  ov_lkn:  the near input of ov_lk (PD-LOSAT-NCBI-DEFECTS exception 23, the Owner's answer Q5 = B of M4b, 2026-11-22):
           A, a spacer of 240 `X` (no word hit; at least the walk's margin, `near_db`), then B, with A's and B's
           masks, over which the alias `ov_lkn.pal` keeps ov_lk's number of sequences and total length: NCBI searches
           B with every diagonal entry of A out of reach, as LOSAT searches B of ov_lk where NCBI's one-thread start
           of B does not meet the entries that no sequence wrote (the two-hit scan).
  ov_lkx:  the other near input of ov_lk (`blank_db`): A's letters all `X` (it stores no entry), its masks and length
           the same, so NCBI's offsets are ov_lk's: NCBI searches B at ov_lk's one-thread offset with no entry of A, as
           LOSAT searches B of ov_lk at every thread count (in the one-hit scan B starts 98 below 0, on the entries
           that no sequence wrote: DECISIONS 55 of series M).
  ov_zt:   subjects A (`XXXXX`, no word hit) then B (240 letters, masks [10,100) [20,110): chunks 10, -80, 130, the
           one-hit walk 0, 10, -70; B[110,122) = the query's letters 18-29) of review 2 of M4b-3b (probe_zt, seed
           20261123): in the one-hit scan NCBI's one-thread search starts B at 5, and its walk meets the entries that
           no sequence wrote; `qzt.faa` / `qzt.fna` its query. `ov_ztl`: A of 100 `X` (B's walk stays above them);
           `ov_zto.pal`: an alias over ov_zt of its OID 2 alone (B is the first searched sequence).

    python3 make_dbmaskp_ov.py [--out DIR] [--ncbi-bin DIR]

writes DIR/ov_*.faa, DIR/ov_*.asnt, DIR/qov1.faa, DIR/qov1.fna, DIR/qov2.faa, DIR/qov2.fna, DIR/qwpe.faa, DIR/qzt.faa,
DIR/qzt.fna, DIR/ov_zto.pal (DIR defaults to db_search/dbmaskp next to this script), and with --ncbi-bin runs
`makeblastdb -dbtype prot -out DIR/ov_X -title "overlapping masks X" -mask_data ov_X.asnt` (BLAST_USAGE_REPORT=false)
and removes the .pjs files; ov_lkn as `near_db` says, ov_lkx as `blank_db` says.
"""
import argparse
import os
import random
import re
import subprocess
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVY"
BASES = "TCAG"
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
L = 240


def put(seq, pos, s):
    return seq[:pos] + s + seq[pos + len(s):]


def codons():
    table = {}
    for i, a in enumerate(BASES):
        for j, b in enumerate(BASES):
            for k, c in enumerate(BASES):
                table.setdefault(CODE[16 * i + 4 * j + k], []).append(a + b + c)
    return table


def back_translate(r, protein):
    c = codons()
    return "".join(r.choice(c[a]) for a in protein)


def asnt(masks_by_subject, options):
    """masks_by_subject: list of (ordinal, [(from, to_exclusive), ...]); one packed-int per subject."""
    blocks = []
    for ordinal, masks in masks_by_subject:
        ivs = ",\n".join(f"        {{\n          from {a},\n          to {b - 1},\n          id local id {ordinal}\n        }}"
                         for a, b in masks)
        blocks.append(f"      packed-int {{\n{ivs}\n      }}")
    body = ",\n".join(blocks)
    return ("Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program seg,\n  algo-options \"" + options + "\",\n"
            f"  masks {{\n    masks {{\n{body}\n    }},\n    more FALSE\n  }}\n}}\n")


def probe1():
    """make_p1.py: returns (query, back-translated query, {db: [(seq, masks)]})."""
    r = random.Random(20261121)

    def rnd(n):
        return "".join(r.choice(AA) for _ in range(n))
    dbs = {}
    base = rnd(L)
    s = put(put(base, 60, "WWW"), 70, "WWW")
    dbs["e3"] = [(s, [(40, 60), (55, 70)])]
    dbs["e3c"] = [(s, [(40, 60), (60, 70)])]
    s = put(put(rnd(L), 60, "WWW"), 80, "WWW")
    dbs["g3"] = [(s, [(40, 60), (55, 80), (75, 90)])]
    s = put(put(rnd(L), 100, "WWW"), 64, "WWW")
    dbs["nest"] = [(s, [(40, 100), (50, 60)])]
    s = put(put(put(rnd(L), 120, "WWW"), 130, "WWW"), 150, "WWW")
    dbs["big"] = [(s, [(40, 120), (60, 130)])]
    q = rnd(30) + "WWW" + rnd(30)
    return q, back_translate(r, q), dbs


def probe2():
    """make_p2.py (q, the shapes s2 r37 r36 nb and lk; its g3 is probe 1's)."""
    r = random.Random(20261122)

    def rnd(n):
        return "".join(r.choice(AA) for _ in range(n))
    q = rnd(30) + "WWW" + rnd(30)
    dbs = {}
    rnd(L)  # probe 2's own g3 (not frozen: the probe-1 g3 is ov_g3); keeps the draws of the others identical
    s = put(put(rnd(L), 60, "WWW"), 70, "WWW")
    dbs["s2"] = [(s, [(40, 60), (55, 70), (72, 90)])]
    # ov_tl: s2's letters, no draw of its own (the others' draws stay identical)
    dbs["tl"] = [(s, [(65, L), (55, 64)])]
    s = put(put(rnd(L), 80, "WWW"), 90, "WWW")
    dbs["r37"] = [(s, [(40, 80), (43, 90)])]
    s = put(put(rnd(L), 80, "WWW"), 90, "WWW")
    dbs["r36"] = [(s, [(40, 80), (44, 90)])]
    s = put(put(put(rnd(L), 60, "WCW"), 70, "WMW"), 80, "WWW")
    dbs["nb"] = [(s, [(40, 60), (50, 70), (65, 80), (100, 120), (110, 130)])]
    seg = q[:33]
    a = put(rnd(L), 7, seg)
    b = put(rnd(L), 7, seg)
    dbs["lk"] = [(a, [(40, 200), (60, 238)]), (b, [(230, 235)])]
    return q, back_translate(r, q), dbs


def probe3():
    """Fix round 1 of M4b-3b (review H2; seed 20261123): `wpe`, subjects A then B: A (100 letters) masks [10,60) [20,99)
    [30,40), whose one word at 99 covers A's sentinel and B's first letter (`W`, as A's letter 99); `WWW` at 45 and 85
    of A and at 120 of B. The query of probe 2 (no `-`) and `qwpe.faa`, probe 2's query with `W-W` after it."""
    r = random.Random(20261123)

    def rnd(n):
        return "".join(r.choice(AA) for _ in range(n))
    a = put(put(put(rnd(100), 45, "WWW"), 85, "WWW"), 99, "W")
    b = put(put(rnd(L), 0, "W"), 120, "WWW")
    return {"wpe": [(a, [(10, 60), (20, 99), (30, 40)]), (b, [])]}


def probe4():
    """Review 2 of M4b-3b, probe zt (`ncbi-parity-20261019/M/m4b3b/review/scripts/make_zt.py`, seed 20261123; its draws
    unchanged, then the back-translation): returns (query, back-translated query, {db: [(seq, masks)]})."""
    r = random.Random(20261123)

    def rnd(n):
        return "".join(r.choice(AA) for _ in range(n))
    q = rnd(20) + "WWHCMWYF" + rnd(20)
    seg = q[18:30]
    b = rnd(L)
    b = b[:110] + seg + b[110 + len(seg):]
    bm = [(10, 100), (20, 110)]
    dbs = {"zt": [("XXXXX", []), (b, bm)], "ztl": [("X" * 100, []), (b, bm)]}
    return q, back_translate(r, q), dbs


def walk_stores(length, masks, word, window):
    """The walk of a sequence's hard chunks in NCBI's word finder, as LOSAT's walk reads it
    (`algorithm::common::hard_chunk_walk::DiagonalWalk`, every one word taken to hit), relative to the offset before its
    first chunk: the lowest offset of a chunk that reads a word, the highest entry it can store without the flag and
    with it (two-hit; `None` for none), and the offset after its last chunk. The chunks are the ranges between the masks
    in their stored order (`core::subject_split::complement_masks`), each moving the offset by its length (negative for
    an inverted range) plus the window."""
    starts = [0] + [end for _, end in masks]
    ends = [start for start, _ in masks] + [length]
    at, lowest, unflagged, flagged = 0, None, None, None
    for left, right in zip(starts, ends):
        size = right - left
        if size >= word or window > 0:
            lowest = at if lowest is None else min(lowest, at)
            if window > 0:
                stored, stored_flag = at + max(size - word, 0), at + max(size - word + 1, 1)
                flagged = stored_flag if flagged is None else max(flagged, stored_flag)
            else:
                stored = at + size - word + 2
            unflagged = stored if unflagged is None else max(unflagged, stored)
        at += size + window
    return lowest, unflagged, flagged, at


def spacer_length(records, masks_of, before, scans):
    """The spacer of `near_db`: the database's average length, or longer where the walk needs it. NCBI's scan moves the
    offset past a spacer of `S` letters by `S` plus the window; a hit at offset `x` reads an entry `e` that an earlier
    sequence stored as it reads a cleared one where `x - e` is at least the window without the flag (two-hit; not below
    0 in the one-hit scan, window 0) and where `x` is not below `e` with the flag. So `S` must take the lowest chunk
    that reads a word of the sequence after the spacer (`walk_stores`) that far past every entry stored before it, for
    each scan `(word, window)` of `scans` (NCBI's one-thread walk from its new table at the window; a longer spacer
    only moves the later sequences further on)."""
    count, total = len(records), sum(len(s) for _, s in records)
    spacer = total // count
    for word, window in scans:
        offset, top, top_flag = window, None, None
        for i, (_, seq) in enumerate(records):
            lowest, stored, stored_flag, end = walk_stores(len(seq), masks_of.get(i, []), word, window)
            if i in before:
                if lowest is not None:
                    if top is not None:
                        spacer = max(spacer, top - offset - lowest)
                    if top_flag is not None:
                        spacer = max(spacer, top_flag - offset - lowest - window)
                offset += spacer + window
            if stored is not None:
                top = offset + stored if top is None else max(top, offset + stored)
            if stored_flag is not None:
                top_flag = offset + stored_flag if top_flag is None else max(top_flag, offset + stored_flag)
            offset += end
    return spacer


def blank_db(out, db, title, records, masks, blank, ncbi_bin):
    """The near input `{db}x` of a database whose diagonal entries reach across subjects (PD-LOSAT-NCBI-DEFECTS
    exception 23): `records` [(header, letters)] in order with the letters of the 0-based orders `blank` all `X` (no
    word of a query without `X` hits them, so they store no diagonal entry), `masks` its `Blast-db-mask-info` ASN.1
    text, unchanged: the lengths and masks are the database's, so NCBI's one-thread offsets and statistics are its
    own. Writes `{db}x.faa` and `{db}x.asnt`, and with `ncbi_bin` `makeblastdb -out {db}x`."""
    name = db + "x"
    with open(out / f"{name}.faa", "w") as f:
        for i, (header, seq) in enumerate(records):
            letters = "X" * len(seq) if i in blank else seq
            f.write(f">{header}\n" + "\n".join(letters[j:j + 60] for j in range(0, len(letters), 60)) + "\n")
    (out / f"{name}.asnt").write_text(masks)
    if ncbi_bin:
        env = dict(os.environ, BLAST_USAGE_REPORT="false")
        subprocess.run([str(Path(ncbi_bin) / "makeblastdb"), "-in", f"{name}.faa", "-dbtype", "prot", "-out", name,
                        "-title", title, "-mask_data", f"{name}.asnt"], cwd=out, env=env, check=True,
                       stdout=subprocess.DEVNULL)
        (out / f"{name}.pjs").unlink(missing_ok=True)


def near_db(out, db, title, records, masks, before, ncbi_bin, masks_of=None):
    """The near input `{db}n` of a database whose diagonal entries reach across subjects (PD-LOSAT-NCBI-DEFECTS exception
    23): `records` [(header, letters)] in order, `masks` its `Blast-db-mask-info` ASN.1 text (`local id N`, N the 1-based
    order), `masks_of` the same masks by 0-based order, `before` the 0-based orders that get a spacer before them (never
    0). Each spacer is a run of `X`, which no word of a query without `X` hits, at least as long as the walk's margin
    less the window (`spacer_length`: the database's average length, total // count, or longer), so that NCBI's scan
    moves the diagonal offset past it by its length and the window and searches the next subject with every entry of
    the ones before it out of reach. Writes `{db}n.faa` and `{db}n.asnt` (the masks with N renumbered), and with
    `ncbi_bin` `makeblastdb -out {db}nv` and the alias `{db}n.pal` over it, whose NSEQ, LENGTH, STATS_NSEQ and
    STATS_TOTLEN are the database's (the statistics and the reports read those, not the spacers)."""
    assert before and 0 not in before
    count, total = len(records), sum(len(s) for _, s in records)
    # The scans of the rows: two-hit at the window 40 and one-hit, word 3.
    spacer = "X" * spacer_length(records, masks_of or {}, before, [(3, 40), (3, 0)])
    rows, order = [], {}
    for i, (header, seq) in enumerate(records):
        if i in before:
            rows.append((f"xspacer{len(rows) + 1} hitless spacer", spacer))
        rows.append((header, seq))
        order[i + 1] = len(rows)
    name = db + "n"
    with open(out / f"{name}.faa", "w") as f:
        for header, seq in rows:
            f.write(f">{header}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    (out / f"{name}.asnt").write_text(re.sub(r"local id (\d+)", lambda m: f"local id {order[int(m.group(1))]}", masks))
    if ncbi_bin:
        env = dict(os.environ, BLAST_USAGE_REPORT="false")
        subprocess.run([str(Path(ncbi_bin) / "makeblastdb"), "-in", f"{name}.faa", "-dbtype", "prot", "-out", name + "v",
                        "-title", title, "-mask_data", f"{name}.asnt"], cwd=out, env=env, check=True,
                       stdout=subprocess.DEVNULL)
        (out / f"{name}v.pjs").unlink(missing_ok=True)
        (out / f"{name}.pal").write_text(f"TITLE {title}\nDBLIST {name}v\nNSEQ {count}\nLENGTH {total}\n"
                                         f"STATS_NSEQ {count}\nSTATS_TOTLEN {total}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "db_search" / "dbmaskp")
    ap.add_argument("--ncbi-bin", default="")
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    q2 = probe2()[0]
    (out / "qwpe.faa").write_text(f">qd query with a gap letter\n{q2}W-W\n")
    for tag, (q, qx, dbs) in (("qov1", probe1()), ("qov2", probe2()), (None, (None, None, probe3())),
                              ("qzt", probe4())):
        if tag is not None:
            (out / f"{tag}.faa").write_text(f">q query\n{q}\n")
            (out / f"{tag}.fna").write_text(f">qx query back-translated\n{qx}\n")
        for name, subjects in dbs.items():
            db = "ov_" + name
            with open(out / f"{db}.faa", "w") as f:
                for i, (seq, _) in enumerate(subjects):
                    sid = name if len(subjects) == 1 else f"{name}{'AB'[i]}"
                    f.write(f">{sid} probe\n{seq}\n")
            (out / f"{db}.asnt").write_text(asnt([(i + 1, m) for i, (_, m) in enumerate(subjects) if m],
                                                 "overlap " + name))
            if args.ncbi_bin:
                env = dict(os.environ, BLAST_USAGE_REPORT="false")
                subprocess.run([str(Path(args.ncbi_bin) / "makeblastdb"), "-in", f"{db}.faa", "-dbtype", "prot", "-out", db,
                                "-title", f"overlapping masks {name}", "-mask_data", f"{db}.asnt"],
                               cwd=out, env=env, check=True, stdout=subprocess.DEVNULL)
                (out / f"{db}.pjs").unlink(missing_ok=True)
    subjects = probe2()[2]["lk"]
    records = [(f"lk{'AB'[i]} probe", s) for i, (s, _) in enumerate(subjects)]
    lk_masks = asnt([(i + 1, m) for i, (_, m) in enumerate(subjects)], "overlap lk")
    near_db(out, "ov_lk", "overlapping masks lk", records, lk_masks, {1}, args.ncbi_bin,
            {i: m for i, (_, m) in enumerate(subjects)})
    blank_db(out, "ov_lk", "overlapping masks lk", records, lk_masks, {0}, args.ncbi_bin)
    # B of ov_zt alone (OID 2, 1-based in the alias file): the first searched sequence.
    (out / "ov_zto.pal").write_text("TITLE overlapping masks zt, OID 2\nDBLIST ov_zt\nFIRST_OID 2\nLAST_OID 2\n")


if __name__ == "__main__":
    main()
