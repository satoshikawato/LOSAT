#!/usr/bin/env python3
"""Builds tests/db_search/multidb/: the databases of the rows db.mdb.* (DI-8f: -db with several names).

Usage: make_multidb.py [--ordsl-only | --mixdb2-only | --mixtax-only | --loadord-only] MAKEBLASTDB   (writes multidb/*.fa, multidb/*.map and the databases below)

The records are those of taxfilter/tfp.faa (23 protein globins) and taxfilter/tfn.fna (16 globin cDNAs); the definition
lines (some records carry a second one) stay as they are, the taxids come from taxfilter/tfp.map and tfn.map.

- `pa` / `pb`: version 5 protein, `-parse_seqids`, `-taxid_map`; `pa` = records 1-12 of tfp.faa, `pb` = records 13-23 and
  `Q00001.1`, an exact copy of record 2's sequence (taxid 9606), so that two databases hold an identical sequence.
- `na` / `nb`: version 5 nucleotide, `-parse_seqids`, `-taxid_map`; `na` = records 1-8 of tfn.fna, `nb` = records 9-16 and
  `ENSTTEST0001.1`, a copy of record 1's sequence (taxid 9606).
- `nnp`: version 5 nucleotide without `-parse_seqids` (records of nb).
- `n4a` / `n4b`: version 4 nucleotide without `-parse_seqids` (records of na / nb).
- `p4a` / `p4b`: version 4 protein without `-parse_seqids` (records of pa / pb).
- `d1/g` and `d2/g`: version 5 `-parse_seqids` nucleotide of the records of na and nb, one base name in two directories.
- `et`: version 5 `-parse_seqids` nucleotide of the records of na with an empty title.
- `pgim`: version 5 protein, `-parse_seqids`, records 1-2 of tfp.faa whose only Seq-id is a giim (`gim|70001`, `gim|70002`;
  line DBID mixdb, review 1 H1: NCBI's prefetch looks a giim up in every LMDB file of a list).
- `al.nal`: an alias file (TITLE multidb alias of na, DBLIST na).
- `mix_u.nal`: an alias file of `na nnp` whose FIRST_OID 9 and LAST_OID 17 keep the sequences of nnp alone (line DBID mixdb).
- `sl.txt`: a seqid list of the first accession of na, nb, pa and pb.
- `ov_p` / `ov_n` (line DBID ordsl): `makeblastdb -input_type blastdb` of ../taxfilter/tfu_p and tfu_n (version 5, no
  parsed Seq-ids, four records) in two volumes of two records (`-max_file_sz` 400 and 500 bytes) under one LMDB file:
  makeblastdb stores `gnl|BL_ORD_ID|n` with n the OID in the database and keys the accession table with it, and NCBI's
  header reads add the volume's first OID (2) to the stored n of the second volume.
- `sl_ord0.txt`, `sl_ord1.txt`, `sl_ord2.txt`, `sl_ord9.txt`, `sl_ord13.txt`: the seqid list `gnl|BL_ORD_ID|n`;
  `sl_ord1.bsl` the binary list of `sl_ord1.txt` (`blastdb_aliastool -seqid_file_out`).
- `pxm` / `nxm` (line DBID ordsl, review 1): version 5, `-parse_seqids`, two records named to sort after `pa` / `na`: the
  letters of record 13 of tfp.faa (record 9 of tfn.fna) under two definition lines, `ref|NP_900003.1|` (`NM_`) and
  `gnl|BL_ORD_ID|7`, and a random sequence (seed 20261115) under `ref|NP_900005.1|gnl|BL_ORD_ID|0`, whose tag makeblastdb
  stores as the string "0" (a tag with a leading 0 is not made a number). `xq.faa`: that random protein as a query.
  Lists `sl_pxm7.txt` and `sl_pxm19.txt` (`NP_900003.1` and `gnl|BL_ORD_ID|7` / `19`: 7 plus the 12 sequences of pa),
  `sl_nxm7.txt`, `sl_ord7.txt`, `sl_ord15.txt` (7 plus the 8 of na), `sl_ord19.txt`, `sl_np5.txt` (`NP_900005.1`).
  `make_multidb.py --ordsl-only MAKEBLASTDB` builds only these and the ones above.
- `put` (line DBID mixdb2, review 1 H1): version 5 protein without `-parse_seqids`, the three proteins of
  ../ids/ids_q.faa (those of ../ids/noacc_v5) with taxid 9606 (`-taxid`), so that its loads look the taxonomy up in a
  list with `pa`.
- `mix2_q.faa`, `mix2_qx.fna`, `mix2_q.fna`, `mix2_qt.faa` (line DBID mixdb2): two queries each, the first of which
  hits only a database without the accession table (`-evalue 1e-10`), the second the parsed one (../taxfilter/q.faa or
  q.fna): `qp1` of ../ids/ids_q.faa (noacc_v5, put), the same back-translated (one codon per amino acid), the
  `FBtr0083245` record of nnp.fa (glob1), and the longest ORF of that record (153 aa).
  `make_multidb.py --mixdb2-only MAKEBLASTDB` builds only `put` and these.
- `nut` (line DBID mixtax): version 5 nucleotide without `-parse_seqids`, the records of nnp.fa with taxid 9606
  (`-taxid`): `put`'s counterpart for the nucleotide programs, whose loads look the taxonomy up in a list with `na`.
  `mix2_qc.faa` and `mix2_qcx.fna` (line DBID mixtax): the first query of mix2_q.faa and mix2_qx.fna, then one whose
  letters are those of ../taxfilter/q.faa followed by `qp1` (back-translated in the second), which hits both `pa` and
  `put`. `putl` and `nutl` (line DBID mixtax, review 1 L2): `qp2` of ../ids/ids_q.faa (back-translated in `nutl`) with
  taxid 99990, the last entry of ../../seqdb/taxdb/<variant>, without the table; `mix2_ql.faa`, `mix2_qlx.fna`,
  `mix2_ql.fna` and `mix2_qlt.faa`: the first query of mix2_q.faa, mix2_qx.fna, mix2_q.fna and mix2_qt.faa, then the
  letters of ../taxfilter/q.faa or q.fna followed by `qp2` (back-translated in the nucleotide files), which hits `pa` or
  `na`, `put` or `nut` and the volume `putl` or `nutl` that the first query does not.
  `make_multidb.py --mixtax-only MAKEBLASTDB` builds only `nut`, `putl`, `nutl` and these.
- `ldp` and `ldpn` (line DBID loadord): version 5, `-parse_seqids`, taxid 9606: `ref|NP_000002.1|` (`NM_`) with the letters
  of ../taxfilter/q.faa (q.fna), then `ref|NP_000900.1|` (`NM_`) with random letters (seed 20261116). Beside a database
  without the accession table, NCBI's prefetch loads the batch's Seq-align ids in the order of their Seq-id handles,
  `gnl|BL_ORD_ID|n` by n (0 last of the numbers) and the accession by its number 2, and stops at the accession's lookup:
  `lda` (`ldan`) holds `qp1` of ../ids/ids_q.faa (the glob1 record `FBtr0083245` of nnp.fa) alone, at OID 0 in a list with
  `ldp` (`ldpn`), which loads after the stop; `ldb` (`ldbn`) a random sequence and then that one, at OID 1, which loads
  before it; `ldu` (`ldun`) that one alone, sorted after `ldp` (`ldpn`), at OID 2: the same packed number 2, which NCBI
  orders by the address of the handles. Without `-parse_seqids`, taxid 9606 (`-taxid`), so that their loads open the
  taxonomy database. `ldg` (`ldgn`): `gi|1|ref|NP_000002.1|` (`NM_`), the letters of `ldp`'s first record under a gi whose
  load (handle 1, before the accession's 2) registers the accession in NCBI's data loader, so that NCBI never looks it up
  (review 1 of loadord, H1); `ldq` (`ldqn`): `lcl|Query_1` with those letters, which NCBI's scope finds among the queries
  (H2); both `-parse_seqids`. `ldo_qc.faa`: one query, the letters of ../taxfilter/q.faa followed by `qp1`, which hits both
  (psiblast's batch is one query). `make_multidb.py --loadord-only MAKEBLASTDB` builds only these.

`.pjs` and `.njs` (the metadata with the build's paths and dates; no search reads them) are dropped.
"""
import os
import random
import re
import subprocess
import sys

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "multidb")
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "taxfilter")


def makeblastdb(binary, *args):
    subprocess.run([binary, *args], check=True, cwd=HERE, stdout=subprocess.DEVNULL)


def read_fasta(name):
    return read_fasta_path(os.path.join(SRC, name))


def read_fasta_path(path):
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                records.append([line, []])
            else:
                records[-1][1].append(line)
    return records


def read_map(name):
    with open(os.path.join(SRC, name)) as f:
        return dict(line.split() for line in f if line.strip())


def accessions(header):
    """The first token of each definition line (control-A separated); `sp|P69905.2|` gives P69905.2."""
    found = []
    for part in header[1:].rstrip("\n").split("\x01"):
        token = part.split(" ", 1)[0]
        found.append(token.split("|")[1] if "|" in token else token)
    return found


def write_set(base, records, taxids):
    """Writes base.fa and, when `taxids` is not None, base.map (the map lines of the accessions in the records)."""
    with open(os.path.join(HERE, base + ".fa"), "w") as f:
        for header, lines in records:
            f.write(header)
            f.writelines(lines)
    if taxids is not None:
        with open(os.path.join(HERE, base + ".map"), "w") as f:
            for header, _ in records:
                for acc in accessions(header):
                    if acc in taxids:
                        f.write(f"{acc} {taxids[acc]}\n")


def make_pgim(binary, tfp):
    """Writes pgim.fa (records 1-2 of tfp.faa, each under a giim id before its title) and builds pgim."""
    records = [[f">gim|{70001 + i} " + header[1:].split(" ", 1)[1].split("\x01")[0].rstrip("\n") + "\n", lines]
               for i, (header, lines) in enumerate(tfp[:2])]
    write_set("pgim", records, None)
    makeblastdb(binary, "-in", "pgim.fa", "-dbtype", "prot", "-out", "pgim", "-parse_seqids", "-blastdb_version", "5",
                "-title", "multidb pgim")


def make_ordsl(binary):
    """Builds ov_p and ov_n and writes the ordinal-id seqid lists (line DBID ordsl)."""
    for kind, source, base, size in (("prot", "tfu_p", "ov_p", "400"), ("nucl", "tfu_n", "ov_n", "500")):
        makeblastdb(binary, "-in", os.path.join(SRC, source), "-input_type", "blastdb", "-dbtype", kind, "-out", base,
                    "-blastdb_version", "5", "-max_file_sz", size, "-title", "multidb " + base)
        os.remove(os.path.join(HERE, f"{base}.{kind[0]}js"))
    for n in (0, 1, 2, 9, 13):
        with open(os.path.join(HERE, f"sl_ord{n}.txt"), "w") as f:
            f.write(f"gnl|BL_ORD_ID|{n}\n")
    aliastool = os.path.join(os.path.dirname(binary), "blastdb_aliastool")
    subprocess.run([aliastool, "-seqid_file_in", "sl_ord1.txt", "-seqid_file_out", "sl_ord1.bsl"], check=True, cwd=HERE,
                   stdout=subprocess.DEVNULL)
    rng = random.Random(20261115)
    tfp, tfn = read_fasta("tfp.faa"), read_fasta("tfn.fna")
    for base, kind, source, prefix, alphabet, length in (("pxm", "prot", tfp[12], "NP", "ACDEFGHIKLMNPQRSTVWY", 120),
                                                          ("nxm", "nucl", tfn[8], "NM", "ACGT", 300)):
        letters = "".join(rng.choice(alphabet) for _ in range(length))
        random_lines = [letters[at:at + 60] + "\n" for at in range(0, length, 60)]
        records = [[f">ref|{prefix}_900003.1| two-line t3a\x01gnl|BL_ORD_ID|7 two-line t3b\n", source[1]],
                   [f">ref|{prefix}_900005.1|gnl|BL_ORD_ID|0 string tag t5\n", random_lines]]
        write_set(base, records, None)
        makeblastdb(binary, "-in", base + ".fa", "-dbtype", kind, "-out", base, "-parse_seqids", "-blastdb_version", "5",
                    "-title", "multidb " + base)
        os.remove(os.path.join(HERE, f"{base}.{kind[0]}js"))
        if kind == "prot":
            with open(os.path.join(HERE, "xq.faa"), "w") as f:
                f.write(">xq the random protein of pxm\n" + "".join(random_lines))
    lists = {"sl_pxm7.txt": "NP_900003.1\ngnl|BL_ORD_ID|7\n", "sl_pxm19.txt": "NP_900003.1\ngnl|BL_ORD_ID|19\n",
             "sl_nxm7.txt": "NM_900003.1\ngnl|BL_ORD_ID|7\n", "sl_ord7.txt": "gnl|BL_ORD_ID|7\n",
             "sl_ord15.txt": "gnl|BL_ORD_ID|15\n", "sl_ord19.txt": "gnl|BL_ORD_ID|19\n", "sl_np5.txt": "NP_900005.1\n"}
    for name, text in lists.items():
        with open(os.path.join(HERE, name), "w") as f:
            f.write(text)


def make_mixdb2(binary):
    """Builds put and writes the queries mix2_q*.f* (line DBID mixdb2)."""
    ids_q = os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa")
    makeblastdb(binary, "-in", ids_q, "-dbtype", "prot", "-out", "put", "-blastdb_version", "5", "-taxid", "9606",
                "-title", "multidb put")
    os.remove(os.path.join(HERE, "put.pjs"))

    def letters(path, index=0):
        header, lines = read_fasta_path(path)[index]
        return header.rstrip("\n"), "".join(line.strip() for line in lines)

    qp1 = letters(ids_q)[1]
    glob1 = next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(os.path.join(HERE, "nnp.fa"))
                 if header.split()[0] == ">FBtr0083245")
    bases = "TCAG"
    amino = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    code = {a + b + c: amino[16 * i + 4 * j + k] for i, a in enumerate(bases) for j, b in enumerate(bases)
            for k, c in enumerate(bases)}
    orf = ""
    for frame in range(3):
        protein = "".join(code.get(glob1[at:at + 3], "X") for at in range(frame, len(glob1) - 2, 3))
        for piece in protein.split("*"):
            start = piece.find("M")
            if start >= 0 and len(piece) - start > len(orf):
                orf = piece[start:]
    codon = {"A": "GCT", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG", "E": "GAA", "G": "GGT", "H": "CAC",
             "I": "ATC", "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTC", "P": "CCG", "S": "TCT", "T": "ACC", "W": "TGG",
             "Y": "TAC", "V": "GTT"}
    protein_q, nucleotide_q = letters(os.path.join(SRC, "q.faa")), letters(os.path.join(SRC, "q.fna"))
    for name, records in (
            ("mix2_q.faa", [(">mix2_u qp1 of ids/noacc_v5 (no accession table)", qp1), protein_q]),
            ("mix2_qx.fna", [(">mix2_u qp1 of ids/noacc_v5 back-translated", "".join(codon[a] for a in qp1)),
                             nucleotide_q]),
            ("mix2_q.fna", [(">mix2_u glob1 of multidb/nnp (no accession table)", glob1), nucleotide_q]),
            ("mix2_qt.faa", [(">mix2_u the glob1 ORF of multidb/nnp", orf), protein_q])):
        with open(os.path.join(HERE, name), "w") as f:
            for header, seq in records:
                f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")


def make_mixtax(binary):
    """Builds nut, putl and nutl and writes the queries mix2_qc*.f* and mix2_ql*.f* (line DBID mixtax)."""
    makeblastdb(binary, "-in", "nnp.fa", "-dbtype", "nucl", "-out", "nut", "-blastdb_version", "5", "-taxid", "9606",
                "-title", "multidb nut")
    os.remove(os.path.join(HERE, "nut.njs"))

    def records(name):
        out = []
        for header, lines in read_fasta_path(os.path.join(HERE, name)):
            out.append((header.rstrip("\n"), "".join(line.strip() for line in lines)))
        return out

    (u_header, qp1), _ = records("mix2_q.faa")
    (ux_header, qp1_x), _ = records("mix2_qx.fna")
    (un_header, glob1), _ = records("mix2_q.fna")
    (ut_header, orf), _ = records("mix2_qt.faa")
    globin = records(os.path.join(SRC, "q.faa"))[0][1]
    globin_x = records(os.path.join(SRC, "q.fna"))[0][1]
    qp2 = records(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa"))[1][1]
    codon = {"A": "GCT", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG", "E": "GAA", "G": "GGT", "H": "CAC",
             "I": "ATC", "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTC", "P": "CCG", "S": "TCT", "T": "ACC", "W": "TGG",
             "Y": "TAC", "V": "GTT"}
    qp2_x = "".join(codon[a] for a in qp2)
    for name, kind, header, seq in (("putl", "prot", ">qp2 unrelated protein", qp2),
                                    ("nutl", "nucl", ">qp2_x unrelated protein back-translated", qp2_x)):
        with open(os.path.join(HERE, name + ".fa"), "w") as f:
            f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")
        makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, "-blastdb_version", "5", "-taxid",
                    "99990", "-title", "multidb " + name)
        os.remove(os.path.join(HERE, name + (".pjs" if kind == "prot" else ".njs")))
    for name, pair in (
            ("mix2_qc.faa", [(u_header, qp1), (">mix2_c the letters of taxfilter/q.faa, then qp1", globin + qp1)]),
            ("mix2_qcx.fna", [(ux_header, qp1_x),
                              (">mix2_c the letters of taxfilter/q.fna, then qp1 back-translated", globin_x + qp1_x)]),
            ("mix2_ql.faa", [(u_header, qp1), (">mix2_l the letters of taxfilter/q.faa, then qp2", globin + qp2)]),
            ("mix2_qlx.fna", [(ux_header, qp1_x),
                              (">mix2_l the letters of taxfilter/q.fna, then qp2 back-translated", globin_x + qp2_x)]),
            ("mix2_ql.fna", [(un_header, glob1),
                             (">mix2_l the letters of taxfilter/q.fna, then qp2 back-translated", globin_x + qp2_x)]),
            ("mix2_qlt.faa", [(ut_header, orf), (">mix2_l the letters of taxfilter/q.faa, then qp2", globin + qp2)])):
        with open(os.path.join(HERE, name), "w") as f:
            for header, seq in pair:
                f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")


def make_loadord(binary):
    """Builds ldp, lda, ldb, ldu and their nucleotide counterparts ldpn, ldan, ldbn, ldun (line DBID loadord)."""
    rng = random.Random(20261116)

    def letters(path, pick=lambda header: True):
        return next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(path) if pick(header))

    def random_letters(alphabet, length):
        return "".join(rng.choice(alphabet) for _ in range(length))

    qp1 = letters(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa"))
    glob1 = letters(os.path.join(HERE, "nnp.fa"), lambda header: header.split()[0] == ">FBtr0083245")
    for kind, suffix, prefix, globin, unparsed, alphabet, length in (
            ("prot", "", "NP", letters(os.path.join(SRC, "q.faa")), qp1, "ACDEFGHIKLMNPQRSTVWY", 120),
            ("nucl", "n", "NM", letters(os.path.join(SRC, "q.fna")), glob1, "ACGT", 360)):
        sets = {
            "ldp": [(f">ref|{prefix}_000002.1| loadord packed number 2", globin),
                    (f">ref|{prefix}_000900.1| loadord random letters", random_letters(alphabet, length))],
            "lda": [(">loadord the unparsed hit at OID 0", unparsed)],
            "ldb": [(">loadord random letters", random_letters(alphabet, length)),
                    (">loadord the unparsed hit at OID 1", unparsed)],
            "ldu": [(">loadord the unparsed hit at OID 2", unparsed)],
            "ldg": [(f">gi|1|ref|{prefix}_000002.1| loadord the accession of ldp under a gi", globin)],
            "ldq": [(">lcl|Query_1 loadord a local id of the form of a generated query id", globin)],
        }
        for base, records in sets.items():
            name = base + suffix
            with open(os.path.join(HERE, name + ".fa"), "w") as f:
                for header, seq in records:
                    f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")
            parse = ["-parse_seqids"] if base in ("ldp", "ldg", "ldq") else []
            makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, *parse, "-blastdb_version", "5",
                        "-taxid", "9606", "-title", "multidb " + name)
            os.remove(os.path.join(HERE, name + (".pjs" if kind == "prot" else ".njs")))
    globin = letters(os.path.join(SRC, "q.faa"))
    with open(os.path.join(HERE, "ldo_qc.faa"), "w") as f:
        seq = globin + qp1
        f.write(">ldo_c the letters of taxfilter/q.faa, then qp1\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60))
                + "\n")


def main(binary):
    os.makedirs(os.path.join(HERE, "d1"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "d2"), exist_ok=True)
    tfp, tfn = read_fasta("tfp.faa"), read_fasta("tfn.fna")
    pmap, nmap = read_map("tfp.map"), read_map("tfn.map")
    pmap["Q00001.1"] = "9606"
    nmap["ENSTTEST0001.1"] = "9606"
    pcopy = [">sp|Q00001.1| copy of record 2 [tie]\n", tfp[1][1]]
    ncopy = [">ENSTTEST0001.1 copy of record 1\n", tfn[0][1]]
    sets = {
        "pa": (tfp[:12], pmap), "pb": (tfp[12:] + [pcopy], pmap),
        "na": (tfn[:8], nmap), "nb": (tfn[8:] + [ncopy], nmap),
    }
    for base, (records, taxids) in sets.items():
        write_set(base, records, taxids)
    write_set("nnp", sets["nb"][0], None)
    for base, src in (("n4a", "na"), ("n4b", "nb"), ("p4a", "pa"), ("p4b", "pb"), ("d1/g", "na"), ("d2/g", "nb"),
                      ("et", "na")):
        write_set(base, sets[src][0], None)
    for base, kind in (("pa", "prot"), ("pb", "prot"), ("na", "nucl"), ("nb", "nucl")):
        makeblastdb(binary, "-in", base + ".fa", "-dbtype", kind, "-out", base, "-parse_seqids", "-blastdb_version", "5",
                    "-taxid_map", base + ".map", "-title", "multidb " + base)
    makeblastdb(binary, "-in", "nnp.fa", "-dbtype", "nucl", "-out", "nnp", "-blastdb_version", "5", "-title", "multidb nnp")
    make_pgim(binary, tfp)
    for base, kind in (("n4a", "nucl"), ("n4b", "nucl"), ("p4a", "prot"), ("p4b", "prot")):
        makeblastdb(binary, "-in", base + ".fa", "-dbtype", kind, "-out", base, "-blastdb_version", "4",
                    "-title", "multidb " + base)
    for d in ("d1", "d2"):
        makeblastdb(binary, "-in", f"{d}/g.fa", "-dbtype", "nucl", "-out", f"{d}/g", "-parse_seqids",
                    "-blastdb_version", "5", "-title", f"multidb {d} g")
    makeblastdb(binary, "-in", "et.fa", "-dbtype", "nucl", "-out", "et", "-parse_seqids", "-blastdb_version", "5",
                "-title", "")
    for name in ("pa.pjs", "pb.pjs", "p4a.pjs", "p4b.pjs", "na.njs", "nb.njs", "nnp.njs", "n4a.njs", "n4b.njs",
                 "d1/g.njs", "d2/g.njs", "et.njs", "pgim.pjs"):
        path = os.path.join(HERE, name)
        if os.path.exists(path):
            os.remove(path)
    with open(os.path.join(HERE, "al.nal"), "w") as f:
        f.write("TITLE multidb alias of na\nDBLIST na\n")
    with open(os.path.join(HERE, "mix_u.nal"), "w") as f:
        f.write("TITLE multidb mix_u: the sequences of nnp in a list with na\nDBLIST na nnp\nFIRST_OID 9\nLAST_OID 17\n")
    firsts = [accessions(sets[b][0][0][0])[0] for b in ("na", "nb", "pa", "pb")]
    with open(os.path.join(HERE, "sl.txt"), "w") as f:
        f.write("".join(acc + "\n" for acc in firsts))


if __name__ == "__main__":
    if sys.argv[1] == "--ordsl-only":
        make_ordsl(sys.argv[2])
    elif sys.argv[1] == "--mixdb2-only":
        make_mixdb2(sys.argv[2])
    elif sys.argv[1] == "--mixtax-only":
        make_mixtax(sys.argv[2])
    elif sys.argv[1] == "--loadord-only":
        make_loadord(sys.argv[2])
    else:
        main(sys.argv[1])
        make_ordsl(sys.argv[1])
        make_mixdb2(sys.argv[1])
        make_mixtax(sys.argv[1])
        make_loadord(sys.argv[1])
