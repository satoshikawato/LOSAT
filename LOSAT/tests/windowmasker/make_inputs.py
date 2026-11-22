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
    write("ids_cr.txt", "rec1\rrec4\n".encode(), "wb")
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
    # A sequence whose first four samples of dm2L_60k.fa (offsets 0 to 30,000) follow each other
    # 9,995 letters apart from its start: the run's start in it falls 15 letters before 0, which
    # NCBI prints as a size_t that wraps (review of M5b, 2026-11-21).
    write("cd_wrap.fa", f">cd_wrap\n{wrap(dm[0:9995] + dm[10000:19995] + dm[20000:29995] + dm[30000:50100])}\n")
    # Id lists of -text_match F (CIdSet_SeqId: each id read as a Seq-id; the first six are those of
    # probe s1 of line NP M5, 2026-11-21): bare words, versions, case, forms, ordinal ids, ids that
    # CSeq_id cannot read (empty, a byte of 0x80 or more, too many parts, ...), a CRLF list.
    write("ids_unver.txt", "AB000001\nlcl|rec1\nrec1\nNC_000099\n")
    write("ids_ord.txt", "gnl|BL_ORD_ID|0\ngnl|BL_ORD_ID|3\n")
    write("ids_case.txt", "ab000001.1\nnc_000099.1\nREC1\nlcl|REC2\n")
    write("ids_forms.txt", "AB000001.2\nref|NC_000099.1|\ngb|AB000001.1|\nsp|P12345|ABC_HUMAN\nP12345\n"
          "BK000001.1\ntpg|BK000001.1|\n")
    write("ids_bad.txt", "gb|\nx y\n|||\nrec1\n")
    write("ids_lcl.txt", "lcl|rec2\nrec2\nrec4\nlcl|rec4.\n")
    write("ids_sp.txt", "P12345\n")
    write("ids_tpa.txt", "BK000001.1\n")
    write("ids_seqid_err.txt", b" x\n>\n> y\ngb||\ngb|A|B|C|D\ngi|abc\ngi|0\ncaf\xe9\nAB000001.1.2\n"
          b"ref|NC_000099.x|\npdb|\npat|US\ngnl|\ngnl|db\ntr|\nlcl|\n", "wb")
    write("ids_seqid_crlf.txt", b"AB000001.1\r\nrec1\r\nlcl|rec2\r\n", "wb")
    # Seq-id lists of the M5c review (2026-11-21): a list id that makes the handle's info of a
    # database id (another alternative, another number, another case), ordinal ids of a list of two
    # databases, and one id each for the database multi_ids.
    write("ids_emb2.txt", "emb|AB000002.1|\n")
    write("ids_case_mi.txt", "ab000001.1\nsp|P12345|abc_human\nREC1\nlcl|REC2\n")
    write("ids_ord3.txt", "gnl|BL_ORD_ID|3\n")
    write("ids_ord8.txt", "gnl|BL_ORD_ID|8\n")
    write("ids_gi5.txt", "5\n")
    write("ids_gi5_bar.txt", "gi|5\n")
    write("ids_ab5.txt", "AB000005.1\n")
    write("ids_ab6.txt", "AB000006.1\n")
    write("ids_aj7.txt", "AJ000007.1\n")
    # Accessions with a protein flag (P) whose division the accession guide does not give as
    # WGS, TSA or targeted: the parse's warning (error code 1306).
    write("ids_warn.txt", "ABCD12P12345\nAAAA01P000001\nlcl|rec2\n")
    # The database multi_ids (made with -parse_seqids by make_fixtures.py): a gi before a GenBank
    # id, two definition lines (^A), a local id before an EMBL id, a local id alone.
    write("multi_ids.fa",
          f">gi|5|gb|AB000005.1| gi first\n{wrap(dm[24000:27000])}\n"
          f">gb|AB000006.1| first line\x01emb|AJ000007.1| second line\n{wrap(dm[27000:30000])}\n"
          f">lcl|q|emb|AJ000008.1| local first\n{wrap(dm[36000:39000])}\n"
          f">lcl|r local alone\n{wrap(dm[39000:41000])}\n")
    # FASTA records whose Seq-id of the handle (GetAccessSeq_id_Handle: a gi, else a Textseq-id,
    # else the first) is not their first Seq-id.
    write("access.fa",
          f">lcl|a|gb|AB000001.1| GenBank after local\n{wrap(dm[0:4000])}\n"
          f">lcl|b|gi|7 gi after local\n{wrap(dm[4000:7000])}\n"
          f">gnl|db|e|ref|NC_000001.1| RefSeq after general\n{wrap(dm[7000:10500])}\n"
          f">emb|AJ000009.1|gb|AB000010.1| two Textseq-ids\n{wrap(dm[10500:13000])}\n"
          f">lcl|z|gb|AB000011.1|gi|9 gi last\n{wrap(dm[13000:15000])}\n")
    # The database wm_tx (made with -parse_seqids -taxid_map tx_taxid.txt by make_fixtures.py;
    # line RG item D5a-t, 2026-11-21): titles that sequence::GetTitle makes from the first source
    # (the taxids of ../seqdb/taxdb/t5: a lower case name, one with spaces around it, one with
    # entities, one starting "TPA: ", an empty name; a taxid it does not hold), a third-party id,
    # a patent id, a title of one comma, a title, and a title of dots.
    tx = [("lcl|t1", "", 100001), ("lcl|t2", ".", 100003), ("lcl|t3", "", 100004),
          ("lcl|t4", "", 100010), ("lcl|t5", "", 424242), ("tpg|BK000001.1|", "", 100001),
          ("tpg|BK000002.1|", "", 100005), ("pat|US|1234|5", "", 100003), ("lcl|t9", ",", 100003),
          ("lcl|t10", "titled", 100001), ("lcl|t11", "...", 100004)]
    write("tx.fa", "".join(f">{i}{' ' + t if t else ''}\n{wrap(dm[15000 + 800 * n:15800 + 800 * n])}\n"
                           for n, (i, t, _) in enumerate(tx)))
    write("tx_taxid.txt", "".join(f"{i}\t{taxid}\n" for i, _, taxid in tx))
    write("ids_tx.txt", "Uncultured\nNo\nTPA:\nSpaced\nSequence\n,\n")
    # Counts files: garbage, empty, and the ascii counts of mg1655_u9 without their ##pct lines.
    write("garbage.counts", "hello\n")
    write("empty.counts", "")
    src = open(os.path.join(WINMASK, "counts", "mg1655_u9.ascii")).read().split("\n")
    write("nopct.ascii", "\n".join(line for line in src if not line.startswith("##pct:")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
