"""Write the inputs of the rows of blastn -filtering_db (line RG item D4) from wm_query.fa.

Usage (from LOSAT/): python3 -I tests/db_search/make_rpt_d4.py
Reads tests/fixtures/blastn_regression/inputs/wm_query.fa and writes
- tests/db_search/rpt_d4.fna: the database of repeats (pieces of the queries, some changed, on
  either strand, and 500 random letters), made into the BLAST database tests/db_search/rpt_d4
  with makeblastdb (README.md);
- tests/fixtures/blastn_regression/inputs/rpt_query_lc.fa: wm_query.fa with more lower case,
  next to and over the repeats;
- tests/fixtures/blastn_regression/inputs/rpt_edge.fa: wm_mq13 in upper case with two lower-case
  runs that end 2 letters before and start 2 letters after the -query_loc intervals of the rows.
Fixed pseudo-random generator (seed 20261121); Python 3 only.
"""
import random


def read_fasta(path):
    records, name, seq = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                records[name] = "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        records[name] = "".join(seq)
    return records


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def mutate(rng, s, rate):
    out = []
    for c in s:
        if rng.random() < rate:
            out.append(rng.choice([b for b in "ACGT" if b != c.upper()]))
        else:
            out.append(c.upper())
    return "".join(out)


def indel(rng, s, every):
    """Insert or delete one letter every `every` letters (alternating)."""
    out, k = [], 0
    for i, c in enumerate(s):
        if i and i % every == 0:
            k += 1
            if k % 2:
                out.append(rng.choice("ACGT"))
                out.append(c)
            continue
        out.append(c)
    return "".join(out)


INPUTS = "tests/fixtures/blastn_regression/inputs"


def write_fasta(path, records):
    with open(path, "w") as out:
        for header, s in records:
            out.write(f">{header}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")


def headers(path):
    return {line[1:].split()[0]: line[1:].rstrip("\n") for line in open(path) if line.startswith(">")}


def lower(s, start, end):
    """`s` with its letters start..end (0-based, end excluded) in lower case."""
    return s[:start] + s[start:end].lower() + s[end:]


def main():
    q = read_fasta(f"{INPUTS}/wm_query.fa")
    names = headers(f"{INPUTS}/wm_query.fa")
    rng = random.Random(20261121)
    pieces = [
        # name, record, from, to, mutation rate, reverse complement, description
        ("rpt_a", "wm_mq00", 100, 400, 0.08, False, "wm_mq00 100-400, 8 % changed"),
        ("rpt_b", "wm_mq03", 500, 620, 0.05, True, "wm_mq03 500-620, 5 % changed, minus strand"),
        ("rpt_c", "wm_mq13", 1000, 1060, 0.0, False, "wm_mq13 1000-1060"),
        ("rpt_d", "wm_mq13", 1064, 1130, 0.0, False, "wm_mq13 1064-1130 (4 letters after rpt_c)"),
        ("rpt_e", "wm_mq13", 1200, 1260, 0.0, True, "wm_mq13 1200-1260, minus strand"),
        ("rpt_f", "wm_mq13", 1267, 1330, 0.0, False, "wm_mq13 1267-1330 (7 letters after rpt_e)"),
        ("rpt_g", "wm_mq35", 2000, 2400, 0.15, False, "wm_mq35 2000-2400, 15 % changed"),
        ("rpt_h", "wm_mq22", 300, 340, 0.0, False, "wm_mq22 300-340"),
        ("rpt_k", "wm_mq02", 0, 77, 0.0, False, "the whole wm_mq02"),
        ("rpt_m", "wm_mq04", 300, 700, 0.03, True, "wm_mq04 300-700, 3 % changed, an indel every 50, minus strand"),
        ("rpt_n", "wm_mq09", 1000, 1400, 0.03, False, "wm_mq09 1000-1400, 3 % changed, an indel every 50"),
    ]
    with open("tests/db_search/rpt_d4.fna", "w") as out:
        for name, rec, a, b, rate, rc, desc in pieces:
            s = mutate(rng, q[rec][a:b], rate)
            if name in ("rpt_m", "rpt_n"):
                s = indel(rng, s, 50)
            if rc:
                s = revcomp(s)
            out.write(f">{name} {desc}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")
        s = "".join(rng.choice("ACGT") for _ in range(500))
        out.write(">rpt_z 500 random letters\n")
        for i in range(0, len(s), 60):
            out.write(s[i:i + 60] + "\n")
    added = {"wm_mq00": (350, 450), "wm_mq13": (1100, 1180), "wm_mq35": (100, 160)}
    write_fasta(
        f"{INPUTS}/rpt_query_lc.fa",
        [(names[name], lower(s, *added[name]) if name in added else s) for name, s in q.items()],
    )
    edge = lower(lower(q["wm_mq13"].upper(), 1159, 1197), 1263, 1300)
    write_fasta(f"{INPUTS}/rpt_edge.fa", [(names["wm_mq13"], edge)])


main()
