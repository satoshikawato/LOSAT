#!/usr/bin/env python3
"""Build the inputs of the `db.x4b2.*` rows (line RG item X4b-2, BLASTX -use_sw_tback, 2026-11-21).

NCBI BLAST+ 2.17.0 makeblastdb is the comparison oracle only. The queries are X4's crafted proteins
(`tests/fasta/x4/sw_subjects.faa`) back-translated with one codon per amino acid, so that frame +1 gives the protein
(the reverse complement gives frame -1), and cuts of `tests/db_search/bslice_q.fna`; the subjects are X4's proteins.

Query files (`tests/fasta/x4b2/`):
- `w.fna`, `w_rc.fna`, `m.fna`: tks_W in frame +1 (BLASTX context 0) and -1 (context 3), tks_M in frame +1.
- `qok1_w12.fna`: qok1, tks_W and ten queries without a hit (`nohit`, CCCCGGGGAACC repeated): twelve queries in
  one batch, tks_W's frame +1 is context 6, below the batch's twelve queries (NCBI's heap test reads query 6's heap).
- `pre20_w12.fna`, `nohit_w12.fna`: the same with qok1's first 20 letters (`pre20`, one HSP on sok1) or a query
  without a hit first.
- `pre20_w.fna`: pre20 and tks_W, two queries (tks_W's context 6 is past the batch's two heaps).
- `pre20_pad_w.fna`: pre20, eleven queries without a hit and tks_W: with BATCH_SIZE=1644 (60 + 11 x 144 letters)
  the first batch has twelve queries and the second tks_W alone.
- `q12.fna`: the first 1000 bases of m_fam, bases 1001-2000 of m_mix (`bslice_q.fna`) and ten queries without a
  hit: every context of the hits is below twelve (NCBI's heap test reads a heap of the batch).
- `big12.fna`: m_fam, m_mix and ten queries without a hit, 13,440 letters: one batch with BATCH_SIZE=30000.
- `q240.fna`, `q120.fna`: the first 240 and 120 bases of m_fam (the `-db` rows of the cutoff score of 1, which NCBI
  redoes slowly: 1.4 to 3.2 seconds for q240 on ws, 0.1 for q120).
- `hf_w.fna` (fix round 1, review 1's S2): `A` + tks_W back-translated, so tks_W is in frame +2 (context 1, which
  reads query 1's heap), qok1 in frame +1 and ten queries without a hit: with `-max_target_seqs 1` qok1's match on sok1
  fills query 1's heap first, but tks_W's window on sfail has an e-value below the heap's 0.002, which a heap takes
  whatever it holds.
- `cut12.fna` (fix round 1, S2): `A` + tks_M's letters 22-33 back-translated (frame +2, context 1), qok1 and ten
  queries without a hit: with `-searchsp 10000000000` the window's e-value is above 0.002, and with
  `-max_target_seqs 1` query 1's full heap turns it away (indexed by the query, it would read query 0's empty heap).
- `pre20_pad.fna`, `pre20_pad_wrc.fna`, `pre20_pad_crash.fna`, `pre20_pad_m.fna` (fix round 1, S3): pre20 and the
  eleven queries without a hit (one batch of 1644 letters with BATCH_SIZE=1644), alone or followed by a second batch:
  tks_W in frame -1 (context 3 of one query: past the array), pre20 again (`pre20b`), tks_W and five more queries
  without a hit (`nohit12` to `nohit16`: tks_W's context 6 below the seven queries, the batch's first query keeps a
  match and the last match fails), or tks_M (the redo without end against X4's endless subjects).
- `empty_w.fna`, `ac_w.fna` (fix round 1, S3): a record without letters (it adds no letter, so tks_W joins its
  batch) or a two-base query (AC: with BATCH_SIZE=1 a first batch without a valid context, which NCBI reports without a
  search), then tks_W.
- `qok1.fna`, `qok1_rc.fna` (fix round 2): qok1 alone in frame +1 and in frame -1 (its hits against X4's subjects lie
  in other frames: contexts past the array of one heap).
- The padded files of approved exception 22 (fix round 2, the Owner's decision Q-HEAP=B, 2026-11-21): NCBI's heap test
  in e-value mode reads past its array of the batch's heaps (one per query) for a query context 6 x q + f at or above
  the batch's query count N; LOSAT tests an empty heap there. Each query batch of the as-is file followed by 5N queries
  without a hit (`pad001` on, CCCCGGGGAACC repeated) keeps every context of its queries inside NCBI's array, and the
  added heaps stay empty: NCBI's run of the padded file is the expectation of the as-is file.
  - `w_pad.fna`, `w_rc_pad.fna`, `qok1_pad.fna`, `qok1_rc_pad.fna`: one query and five pads, in one batch.
  - `pre20_pad_wrc_pad.fna`: `pre20_pad_wrc.fna` and five pads, which join its second batch with BATCH_SIZE=1644.
  - `q2.fna` (m_fam_1k and m_mix_1k, `q12.fna` without its ten queries without a hit) and `q2_b1_pad.fna` (with
    BATCH_SIZE=1720 for the as-is BATCH_SIZE=1): m_fam_1k and five pads (1720 letters), m_mix_1k and five pads (1720
    letters).
  - `big2.fna` (m_fam and m_mix, `big12.fna` without its ten queries without a hit: 10600 letters, one batch at the
    default 10002) and `big2_pad.fna` (with BATCH_SIZE=12040): m_fam, m_mix and ten pads (12040 letters). NCBI splits
    a BLASTX batch of more than one query into no chunks, so both runs search each batch whole.
  - `bh3_pad.fna` (the merge of dev's `-best_hit_*` for BLASTX, #1473): dev's `tests/fasta/outfmt0/bh3_q.fna` (dq0,
    dq1 and dq2, 4,655 letters, one batch at the default 10002) and fifteen pads (6,815 letters, still one batch).
  - No pair has a batch of queries without a hit only: NCBI warns about such a batch's queries (`Could not
    calculate ungapped Karlin-Altschul parameters`) by their numbers in the file, which the pads before them would
    shift.

Subject files: `ok_fail.faa` (sok1, then sfail: the failing match is the stream's last).

Databases (`tests/db_search/x4b2/`, protein, `-parse_seqids`):
- `fail`: sok1 (OID 0) and sfail (OID 1), each with an accession line (`dbj|AB30000N.1|`).
- `excl`: the same two and, at OID 2, `sx`: sfail's first half with 8 substitutions, whose only line has a ti
  (`gnl|ti|9003`). `neg.txt` lists only an id that no database holds: a negative seqid list with ids keeps sx and
  removes its line, and NCBI's fetch excludes it at the redo, after its early termination test; NCBI then sets the
  batch's status back to 0, which clears the failure of tks_W's match on sfail before it.

The databases are written with the creation date in their index files, so the frozen outputs need the committed
files, not a new run of this script.

usage: make_x4b2.py [--ncbi-bin DIR]   (without --ncbi-bin: the FASTA files only)"""
import os, subprocess, sys

CODON = dict(A="GCT", R="CGT", N="AAT", D="GAT", C="TGT", Q="CAA", E="GAA", G="GGT", H="CAT", I="ATT",
             L="CTG", K="AAA", M="ATG", F="TTT", P="CCG", S="TCT", T="ACT", W="TGG", Y="TAT", V="GTT")
COMP = str.maketrans("ACGT", "TGCA")
BODY = "QLILSHHLTMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEVIEGCTTVSSIK"
SOK1 = ("MATAPAALPFGKANAAACGNSNNTFDVLEGLVKRARVVLVFITALSIGLIVVAIMNLRMGMGNSAQGVPLLDEHGQPVLDGNDDPQIEEPNNTYPNKTRDYLI"
        "LMASILSFVIGAASSYTTRLK")
TKS_W, TKS_M = "W" + BODY, "M" + BODY
QOK1 = "W" + SOK1[1:]
PRE20 = QOK1[:20]
SFAIL = "D" + BODY + "D" + BODY
NOHIT = "CCCCGGGGAACC" * 12


def nt(protein):
    return "".join(CODON[a] for a in protein)


def rc(seq):
    return seq.translate(COMP)[::-1]


def weaker(protein):
    """sfail's first half with 8 substitutions (every tenth letter from the 6th becomes A, or G where it is A)."""
    letters = list(protein)
    for at in range(5, len(letters), 10):
        letters[at] = "G" if letters[at] == "A" else "A"
    return "".join(letters)


def fasta(records):
    out = []
    for name, seq in records:
        out.append(">" + name)
        out += [seq[i:i + 70] for i in range(0, len(seq), 70)]
    return "\n".join(out) + "\n"


def read_fasta(path):
    records, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name:
                records[name].append(line)
    return {k: "".join(v) for k, v in records.items()}


args = sys.argv[1:]
NCBI = args[1] if len(args) >= 2 and args[0] == "--ncbi-bin" else None
if args and NCBI is None:
    sys.exit(__doc__)
HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.normpath(os.path.join(HERE, "..", ".."))
FA = os.path.join(TESTS, "fasta", "x4b2")
slices = read_fasta(os.path.join(TESTS, "db_search", "bslice_q.fna"))
w = ("tks_W_nt tks_W in frame +1", nt(TKS_W))
pre20 = ("pre20_nt the first 20 letters of qok1 in frame +1", nt(PRE20))
nohit = [(f"nohit{i:02d} CCCCGGGGAACC repeated", NOHIT) for i in range(1, 17)]
pad = [(f"pad{i:03d} CCCCGGGGAACC repeated", NOHIT) for i in range(1, 16)]
bh3 = list(read_fasta(os.path.join(TESTS, "fasta", "outfmt0", "bh3_q.fna")).items())
qok1 = ("qok1_nt qok1 in frame +1", nt(QOK1))
qok1_rc = ("qok1_rc qok1 in frame -1", rc(nt(QOK1)))
w_rc = ("tks_W_rc tks_W in frame -1", rc(nt(TKS_W)))
m_fam_1k = ("m_fam_1k bases 1-1000 of m_fam", slices["m_fam"][:1000])
m_mix_1k = ("m_mix_1k bases 1001-2000 of m_mix", slices["m_mix"][1000:2000])
m_fam = ("m_fam MelaMJNV 26001-31000", slices["m_fam"])
m_mix = ("m_mix MelaMJNV 56001-61600", slices["m_mix"])
m = ("tks_M_nt tks_M in frame +1", nt(TKS_M))
queries = {
    "w.fna": [w],
    "w_rc.fna": [w_rc],
    "m.fna": [m],
    "qok1_w12.fna": [qok1, w] + nohit[1:11],
    "pre20_w12.fna": [pre20, w] + nohit[1:11],
    "nohit_w12.fna": [nohit[0], w] + nohit[1:11],
    "pre20_w.fna": [pre20, w],
    "pre20_pad_w.fna": [pre20] + nohit[:11] + [w],
    "hf_w.fna": [("f2_tks_W A + tks_W in frame +2", "A" + nt(TKS_W)), qok1] + nohit[1:11],
    "cut12.fna": [("f2_cut12 A + tks_M letters 22-33 in frame +2", "A" + nt(BODY[20:32])), qok1] + nohit[1:11],
    "pre20_pad.fna": [pre20] + nohit[:11],
    "pre20_pad_wrc.fna": [pre20] + nohit[:11] + [w_rc],
    "pre20_pad_crash.fna": [pre20] + nohit[:11] + [("pre20b_nt the first 20 letters of qok1 in frame +1, again", nt(PRE20)), w]
    + nohit[11:16],
    "pre20_pad_m.fna": [pre20] + nohit[:11] + [m],
    "empty_w.fna": [("empty a record without letters", ""), w],
    "ac_w.fna": [("ac two bases", "AC"), w],
    "q12.fna": [m_fam_1k, m_mix_1k] + nohit[1:11],
    "big12.fna": [m_fam, m_mix] + nohit[1:11],
    "qok1.fna": [qok1],
    "qok1_rc.fna": [qok1_rc],
    "w_pad.fna": [w] + pad[:5],
    "w_rc_pad.fna": [w_rc] + pad[:5],
    "qok1_pad.fna": [qok1] + pad[:5],
    "qok1_rc_pad.fna": [qok1_rc] + pad[:5],
    "pre20_pad_wrc_pad.fna": [pre20] + nohit[:11] + [w_rc] + pad[:5],
    "q2.fna": [m_fam_1k, m_mix_1k],
    "big2.fna": [m_fam, m_mix],
    "q2_b1_pad.fna": [m_fam_1k] + pad[:5] + [m_mix_1k] + pad[5:10],
    "big2_pad.fna": [m_fam, m_mix] + pad[:10],
    "bh3_pad.fna": bh3 + pad[:15],
    "q240.fna": [("m_fam_240 bases 1-240 of m_fam", slices["m_fam"][:240])],
    "q120.fna": [("m_fam_120 bases 1-120 of m_fam", slices["m_fam"][:120])],
    "ok_fail.faa": [("sok1 BDT63546.1 of PajaWSV", SOK1),
                    ("sfail tks_M with D first, twice: two HSPs with tks_W, W against D", SFAIL)],
}
os.makedirs(FA, exist_ok=True)
for name, records in queries.items():
    with open(os.path.join(FA, name), "w") as f:
        f.write(fasta(records))
assert sum(len(s) for _, s in queries["pre20_pad_w.fna"][:12]) == 1644


def batches(records, size):
    """NCBI's query batches (blast_input.cpp:140-171): a batch closes at the query that brings it to `size` letters."""
    out, current, letters = [], [], 0
    for name, seq in records:
        current.append(name.split()[0])
        letters += len(seq)
        if letters >= size:
            out.append(current)
            current, letters = [], 0
    return out + [current] if current else out


# Each padded batch is an as-is batch followed by pads (none for the last, whose queries have no hit), enough to keep
# every context of its queries with a hit inside the batch's array of heaps (5N pads keep all 6N contexts inside).
known = {**queries, "bh3_q.fna": bh3}
for as_is, as_size, padded, pad_size in [("q2.fna", 1, "q2_b1_pad.fna", 1720),
                                         ("big2.fna", 10002, "big2_pad.fna", 12040),
                                         ("bh3_q.fna", 10002, "bh3_pad.fna", 10002),
                                         ("pre20_pad_wrc.fna", 1644, "pre20_pad_wrc_pad.fna", 1644)]:
    own, oracle = batches(known[as_is], as_size), batches(known[padded], pad_size)
    stripped = [[name for name in batch if not name.startswith("pad")] for batch in oracle]
    hits = [batch for batch in own if not all(name.startswith("nohit") for name in batch)]
    assert stripped[:len(hits)] == hits, (as_is, own, oracle)
    for batch in oracle[:len(hits)]:
        last_hit = max(index for index, name in enumerate(batch) if not name.startswith(("nohit", "pad")))
        assert len(batch) >= 6 * (last_hit + 1), (padded, batch)
    assert sum(stripped[len(hits):], []) == sum(own[len(hits):], []), (as_is, own, oracle)
assert sum(len(s) for _, s in queries["pre20_pad_crash.fna"][12:]) <= 1644
if NCBI:
    base = [("dbj|AB300001.1| sok1 BDT63546.1 of PajaWSV", SOK1), ("dbj|AB300002.1| sfail tks_M with D first, twice", SFAIL)]
    dbs = {"fail": base, "excl": base + [("gnl|ti|9003 sx sfail's first half with 8 substitutions", weaker("D" + BODY))]}
    for name, records in dbs.items():
        fa = f"{HERE}/{name}.faa"
        with open(fa, "w") as f:
            f.write(fasta(records))
        subprocess.run([f"{NCBI}/makeblastdb", "-in", fa, "-dbtype", "prot", "-parse_seqids", "-out", f"{HERE}/{name}",
                        "-title", name, "-blastdb_version", "5"], check=True, stdout=subprocess.DEVNULL)
        os.remove(fa)
    with open(f"{HERE}/neg.txt", "w") as f:
        f.write("XX999999.1\n")
print("ok", sorted(queries), "databases" if NCBI else "no databases")
