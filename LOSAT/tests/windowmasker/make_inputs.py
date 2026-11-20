#!/usr/bin/env python3
"""Writes the FASTA, id-list, file-list and counts inputs of tests/windowmasker/in from
tests/winmask/fasta/dm2L_60k.fa (2L:1-60,000) and tests/winmask/counts/mg1655_u9.ascii.
Deterministic; run from anywhere: python3 -I make_inputs.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WINMASK = os.path.join(HERE, "..", "winmask")
IN = os.path.join(HERE, "in")


def letters(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">"))


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write(name, text, mode="w"):
    with open(os.path.join(IN, name), mode) as f:
        f.write(text)


def main():
    os.makedirs(IN, exist_ok=True)
    dm = letters(os.path.join(WINMASK, "fasta", "dm2L_60k.fa"))
    assert len(dm) == 60000, len(dm)
    # Records whose definition lines exercise the id lists (first words with and without bars,
    # a trailing period, a tab) and -parse_seqids (local, GenBank, Swiss-Prot, RefSeq, third-party).
    records = [
        ("rec1 first record", dm[0:4000]),
        ("lcl|rec2 local with bar", dm[4000:7000]),
        ("gb|AB000001.1| GenBank style", dm[7000:10500]),
        ("rec4.", dm[10500:13000]),
        ("sp|P12345|ABC_HUMAN protein-like id", dm[13000:15000]),
        ("ref|NC_000099.1|", dm[15000:18000]),
        ("tpg|BK000001.1| third party", dm[18000:20500]),
        ("rec8\twith a tab", dm[20500:22500]),
    ]
    write("multi.fa", "".join(f">{d}\n{wrap(s)}\n" for d, s in records))
    # Letters: lower case, runs of N, IUPAC codes, an all-N record, a record that starts with N.
    a = dm[30000:36000]
    low = a[:2000] + a[2000:4000].lower() + a[4000:]
    nrun = a[:1500] + "N" * 300 + a[1800:3000] + "n" * 7 + a[3007:]
    iu = list(a)
    for k, c in zip(range(50, 6000, 197), "RYKMSWBDHVNU*X." * 10):
        iu[k] = c
    iu = "".join(iu)
    write("letters.fa",
          ">t_lower some description\n" + wrap(low) + "\n"
          + ">t_nrun\n" + wrap(nrun) + "\n"
          + ">t_iupac\n" + wrap(iu) + "\n"
          + ">t_alln\n" + "N" * 100 + "\n"
          + ">t_nstart\n" + wrap("N" * 80 + a[:3000]) + "\n")
    # Short records (1 to 200 letters).
    write("short.fa", "".join(f">t_len{n}\n{dm[40000:40000 + n]}\n" for n in (1, 5, 9, 13, 20, 200)))
    # Records without letters: alone, in the middle, at the end.
    write("hdronly.fa", ">only\n")
    write("midempty.fa", f">m1\n{wrap(dm[0:3000])}\n>m2\n>m3\n{wrap(dm[3000:6000])}\n")
    write("lastempty.fa", f">l1\n{wrap(dm[0:3000])}\n>l2\n")
    write("empty.fa", "")
    write("crlf.fa", (">t_crlf with CRLF\r\n" + wrap(dm[44000:50000]).replace("\n", "\r\n") + "\r\n").encode(), "wb")
    write("gap.fa",
          ">t_gap\n" + wrap(dm[50000:53000]) + "\n>?100 [gap-type=within scaffold] [linkage-evidence=paired-ends]\n"
          + wrap(dm[53000:56000]) + "\n"
          + ">t_gap_plain\n" + wrap(dm[50000:52000]) + "\n>?50\n" + wrap(dm[52000:54000]) + "\n")
    # Id lists of CWinMaskConfig::FillIdList.
    write("ids_a.txt", "rec1\n>rec4\nAB000001.1\nlcl|rec2 extra words\n")
    write("ids_b.txt", "rec2\nP12345|ABC_HUMAN\n\nNC_000099.1|\n")
    write("ids_tab.txt", "rec8\tcomment\nBK000001.1\n")
    write("ids_crlf.txt", "rec1\r\nrec4\r\n".encode(), "wb")
    write("ids_none.txt", "nosuchid\n")
    write("ids_db.txt", "AB000001.1\nrec1\nNC_000099.1\n")
    # -mk_counts: lists of files for -fa_list (names relative to this directory; an empty line, a
    # CRLF line), an all-N genome, and the records of -checkdup: a part of 2L:5,004-56,000, and
    # 3,000 letters of a fixed pseudo-random sequence before 2L:10,001-60,000 (sequences of 50,000
    # letters or more, whose samples of 100 letters every 10,000 match those of dm2L_60k.fa
    # and of each other at offsets that are not multiples of 10,000).
    write("falist.txt", "../winmask/fasta/dm2L_3k.fa\n\nin/multi.fa\n")
    write("falist_crlf.txt", "../winmask/fasta/dm2L_3k.fa\r\n".encode(), "wb")
    write("alln.fa", ">alln\n" + "N" * 100 + "\n")
    state = 12345
    filler = []
    for _ in range(3000):
        state = (state * 1103515245 + 12345) % 2 ** 31
        filler.append("ACGT"[(state >> 16) & 3])
    write("cd.fa", f">cd_part\n{wrap(dm[5003:56000])}\n>cd_shift\n{wrap(''.join(filler) + dm[10000:60000])}\n")
    write("ids_cd.txt", "cd_part\n")
    # Counts files: garbage, empty, and the ascii counts of mg1655_u9 without their ##pct lines.
    write("garbage.counts", "hello\n")
    write("empty.counts", "")
    src = open(os.path.join(WINMASK, "counts", "mg1655_u9.ascii")).read().split("\n")
    write("nopct.ascii", "\n".join(line for line in src if not line.startswith("##pct:")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
