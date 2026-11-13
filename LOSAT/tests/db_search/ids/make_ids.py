#!/usr/bin/env python3
"""Build the databases with parsed Seq-ids of the `db.id.*` rows (line DBID: proteins ID2a, nucleotides ID2b).

NCBI BLAST+ 2.17.0 makeblastdb is the comparison oracle only: it writes the databases from the
records below, and every file it writes is kept (NCBI's formatters find a reported sequence
again through the ISAM or LMDB files, so the frozen NCBI outputs need them; LOSAT reads the
index, sequence and header files and checks that the LMDB file of version 5 exists).

usage: make_ids.py --ncbi-bin DIR [--out LOSAT/tests/db_search/ids]

Records (in this order): `LOSAT/tests/seqid/fasta/prot_ids.fa` (every Seq-id form makeblastdb
writes, Ctrl-A definition lines, GIs, a definition line without a title, `lcl|Subject_` ids),
`prot_extra.fa` (a GI in the second definition line only, 12 titles longer than 60 bytes, an
empty title), `variants.fa` without the records whose Seq-ids another record's lookup would
match (`pdb|1ABC|AA`, `pdb|1ABC|`, `sp|P99999|`, `ref|NP_000002|`: NCBI's formatters could find
the other sequence) and without three records that NCBI's formatters do not find in version 4
(`bbs|`, `gnl|PIG|`, `nat|`; version 5 finds `bbs|` and `nat|`), and `gen.faa` (records whose titles NCBI makes from their Seq-ids: `tpg`
and `tpe` ids with and without titles, `pdb` ids with and without a chain, a patent, RefSeq
NC_/NP_/WP_ accessions without titles, `lcl|Subject_` ids beside other ids). Taxids:
`LOSAT/tests/seqid/fasta/taxmap_prot.txt`.

Outputs: `prot_v4.*` (version 4, every record) and `prot_v5.*` (version 5, without the GI-only
record `gi|1414`, after which NCBI 2.17.0's makeblastdb aborts writing the LMDB file, and without
`lcl|123`, left out while LOSAT refused every local integer id of version 5; NCBI's formatters
find it by its LMDB key unless an earlier sequence holds the key, line DBID N2b).

`refuse_v5.*` (version 5): the records of `refuse.faa` with the taxids of `taxmap_refuse.txt`, each
the subject of its own query `refuse_<case>.faa` (its sequence; `refuse_nf.fna` a back-translation):
`nf` a taxid that the test taxonomy database (`tests/seqdb/taxdb/good`) does not know between two
that it does (NCBI's XML2 and JSON descriptions stop there), and the refusals of LOSAT: `adj` and
`tw` titles ending in `[Homo sapiens]` whose end NCBI rewrites (`x_AdjustProteinTitleSuffix`;
`adj` with a definition line without a title, `tw` a `lcl|Subject_` id), `long` an outfmt 0
description id longer than 65 bytes, `punct` a title that NCBI's `x_CleanAndCompress` reads past
the end of, `disp` a `lcl|777` before `gi|777` (NCBI's formatters find no sequence by the
Seq-align id `gi|777`). The queries hit their own record alone with `-evalue 1e-20` (NCBI's
outfmt 0 of a report that shows `punct` does not end).

Nucleotides (ID2b; every record holds the sequence of `LOSAT/tests/seqid/fasta/nucl_ids.fa`, which the
queries `qnucl.fa` and `qprot.fa` there hit): `nucl_v4.*` and `nucl_v5.*` from `nucl.fna` (the Seq-id forms
of `nucl_ids.fa`, with the GenBank, EMBL and DDBJ accessions each under prefixes of its own partner, so
that NCBI's process-wide Seq-id registry does not print one partner's type for another's (an approved
difference of line DBID, DECISIONS 5 Q2); titles that NCBI makes from the Seq-ids (`tpg` and RefSeq ids
without a title), a GI in the second definition line only, titles longer than 60 bytes, `lcl|Subject_`
beside an accession, a `gnl|ti|` integer tag) with the taxids of `taxmap_nucl.txt`; version 5 without the
GI-only record `gi|1111`. `nrefuse_v5.*` from `nrefuse.fna` (`lcl|777` before `gi|777`: NCBI's formatters
find no sequence by the Seq-align id `gi|777`). `njoined_a.*` and `njoined_b.*` (version 4) from
`njoined_a.fna` and `njoined_b.fna`, read through the alias `njoined.nal`: `gb|AF000040.1|` and
`gb|AF000040|` follow each other in the hits, and NCBI's XML and outfmt 0 take them for one sequence.
`nambig_v4.*` (version 4) from `nambig.fna` (the sequence of the first record of `nucl.fna`, and a copy of it
with runs of N, an R and a Y, so with ambiguity data; `tests/blastn_report_stream.rs` raises its
ambiguity count in a copy, which the traceback of a later batch refuses).

Seq-align ids that NCBI's lookup finds in another record (line DBID N1; the FASTA comes from
`make_found_fasta.py`): `pfound_v5.*` and `pfound_v4.*` from `pfound.faa`, five pairs of unrelated proteins, the
second record of each pair its first's sequence with three substitutions (the query `pfound_q<k>.faa` of pair k
hits both; `pfound_q1.fna` and `pfound_q4.fna` back-translations): `gb|AAA00001.1|` then `gb|AAA00001|` (NCBI
finds the first by `gb|AAA00001|`), `gb|BBB00002.1|` then `emb|BBB00002.1|` (the first by `emb|BBB00002.1|`),
`gb|CCC00003.1|` then `lcl|CCC00003` (no Bioseq by `lcl|CCC00003`), `gb|DDD00004|` then `gb|DDD00004.1|` and
`pdb|2XYZ|` then `pdb|2XYZ|AA` (each finds itself; their Seq-align ids match, so NCBI's XML and outfmt 0 join
them). `pdup_a.*` and `pdup_b.*` (version 5) from `pdup.faa`, the same record `gb|EEE00005.1|` (query
`pdup_q.faa`), read through the alias `pdup_ab.pal`: NCBI finds the first copy by the second's id, which shows the
same; `pdup_c.*` from `pdup_c.faa`, the record with another title, read with `pdup_a` through `pdup_ac.pal`.
`nfound_v5.*` and `nfound_v4.*` from `nfound.fna`: the sequence of `njoined_a.fna` as `gb|AF000050.1|`, then a
copy with three substitutions as `gb|AF000050|`.

Sequences whose only Seq-id is `gnl|BL_ORD_ID|n` (line DBID N2a; the FASTA comes from `make_ord_fasta.py`):
`pord_v5.*` and `pord_v4.*` from `pord.faa`, `pord_x.*` (version 5) from `pord_x.faa`, read after `pord_v5` through
the alias `pord_vx.pal` (which lists `pord_x` first: NCBI orders the volumes of an alias by their base names), and `nord_v5.*` and `nord_v4.*` from `nord.fna`.

Seq-align ids that NCBI's formatters look up like any other (line DBID N2b; the FASTA comes from `make_look_fasta.py`):
`plook_v5.*` and `plook_v4.*` from `plook.faa`, `plook_z.*` (version 5) from `plook_z.faa`, read after `plook_v5`
through the alias `plook_vz.pal`, and `nlook_v5.*` and `nlook_v4.*` from `nlook.fna`. A Seq-align id that another
sequence also holds (the same FASTA script): `pshare_v5.*` and `pshare_v4.*` from `pshare.faa`, `pshare_a.*` and
`pshare_b.*` (version 5) from `pshare_a.faa` and `pshare_b.faa`, read through the alias `pshare_ab.pal`, `pshare_c.*`
and `pshare_d.*` (version 5) through `pshare_cd.pal`, and `nshare_v5.*` from `nshare.fna`; the alias `pdup_abc.pal` reads
`pdup_a`, `pdup_b` and `pdup_c`.

Seq-ids that NCBI's scope matches weakly (line DBID N2c; the FASTA comes from `make_weak_fasta.py`): `pweak_v5.*` and
`pweak_v4.*` from `pweak.faa`, `pweak_r_v5.*` and `pweak_r_v4.*` from `pweak_r.faa` (review 1), `nweak_v5.*` and
`nweak_v4.*` from `nweak.fna`. `--weak-only` builds only these, so
that the other databases keep their bytes.

The same record with other taxids (line DBID N1b; `make_found_fasta.py`), all version 5: `pdup_t.*` from `pdup.faa`
with taxid 9606 and `pdup_u.*` without one, read through `pdup_at.pal` (`pdup_a` first: NCBI shows the copy without
a taxid) and `pdup_tu.pal` (`pdup_t` first); `pdup_n.*` and `pdup_o.*` from `pdup.faa`, both with taxid 99990 (the
last entry of `tests/seqdb/taxdb/*`, whose lookup writes the offset error of `offset_last_past_end`), through
`pdup_no.pal`; `pdupe_a.*` and `pdupe_t.*` from `pdupe.faa` (no title) without a taxid and with 9606, through
`pdupe_at.pal`; `ndup_a.*` from `ndup.fna` without a taxid and `ndup_t.*` with taxid 9606, through `ndup_at.nal`.
`--copytaxid-only` builds only these. Line DBID N1b L1 (titles that NCBI can make from the taxids): `ndupe_a.*` and
`ndupe_t.*` from `ndupe.fna` (no title) without a taxid and with 9606, through `ndupe_at.nal`; `nsubj_t.*` from
`nsubj.fna` (lcl|Subject_1, no title) and `nsubjg_t.*` from `nsubjg.fna` (lcl|Subject_3 beside gb|EEN00007.1|) with
9606; `psubj_a.*`, `psubj_t.*` from `psubj.faa` (lcl|Subject_1, no title)
and `psubjb_a.*`, `psubjb_t.*` from `psubjb.faa` (lcl|Subject_2, a title ending in ']'), through `psubj_at.pal` and
`psubjb_at.pal`. `--titlesrc-only` builds only these. Line DBID H3 (an empty title of a database without parsed
Seq-ids, which NCBI's outfmt 0 makes from the taxid): `nbare_t.*` from `nbare.fna` and `pbare_t.*` from `pbare.faa`
(`>` alone), without -parse_seqids, with taxid 9606. `--emptytitle-only` builds only these. Line DBID tabfields (a
table that reads nothing of the subject prints no id of its title): `ppunct_np.*` (version 4, without -parse_seqids)
from `ppunct.faa`, whose first title `, ,` is punctuation that NCBI's x_CleanAndCompress reads past the end of where
it makes the tabular id. `--tabfields-only` builds only these. Line DBID pfbatch (each batch's prefetch counts its
own Seq-aligns against a positive PRE_FETCH_SEQS_LIMIT): `nordtag_v5.*` from `nordtag.fna`, a sequence whose first id
`gnl|BL_ORD_ID|0` makeblastdb stores with a text tag (NCBI's lookup by that Seq-align id fails, as `pshare`'s) and
another sequence; the queries `nordtag_q.fna` and `nordtag_q.faa` hit one each (`pfb_aba.faa`, `pfb_aab.faa`,
`pfb_aab.fna` and `pfb_s.fna` place it in other batches, see tests/db_search/README.md; `pfb_rc_nambig.fna` and
`pfb_rc_ids_q.fna` are the reverse complements of `nambig.fna` and `ids_q.fna`, for the minus-strand alignments of
`nrefuse_v5`). `--pfbatch-only` builds only these. Line DBID pfwarn (NCBI's prefetch writes `Error pre-fetching
sequence data` where it fails): `ntag_t5.*` from `nordtag.fna`, `nlook_t4.*` (version 4) from `nlook.fna` and
`pshare_t5.*` from `pshare.faa`, every sequence with taxid 99990, whose lookup writes the offset error of
tests/seqdb/taxdb/offset_last_past_end; `pfw_weak_v5.*` from `pfw_weak.faa` (`sp||NAMEY`, `sp|P00042.1|NAMEY`, which
the query `pfw_q6.faa` matches exactly, `gb|FFF00099.1|` and `gnl|BL_ORD_ID|0` with a text tag, whose hit psiblast ranks
last); and `pfw_bad/ntag_t5.*`, a copy of `ntag_t5` whose accession table's leaf page (page 6 of the `.ndb`) has its
flags zeroed, so that the lookup of an accession fails on the file. `--pfwarn-only` builds only these. Line DBID
weakbatch (psiblast's rounds read the sequences before the cut for the next PSSM): `wkb_e3_v5.*` from `wkb_e3.faa`
(`gb|FFF00043.1|`, `gb|ZZZ00001.1|`, then `gb|FFF00043|`, which the lookup finds as the first) and `wkb_p1_v5.*` from
`wkb_p1.faa` (the same with `gb|FFF00043|` a weak hit above -inclusion_ethresh and `gb|WWW00001.1|` before it), for the
query `wkb_qy.faa`. `--weakbatch-only` builds only these. Line DBID pssmrows (psiblast's PSSM rows follow the
Seq-align ids): the databases `psr_*` from the FASTA of the same names (copies of one record, ids that `CSeq_id::Match`
and their controls, for the query `psr_q.faa`) and the aliases `psr_cp.pal`, `psr_cbxp.pal`, `psr_jx.pal`, `psr_cul.pal`,
`psr_cue.pal` and `psr_k.pal`, and the nucleotide database `psr_gimn` from `psr_gimn.fna` (queries `psr_qp.fna`,
`psr_qm.fna`). `--pssmrows-only` builds only these.
Line DBID pssmread (the letters psiblast's PSSM reads by a Seq-align id): `prd_av` (gb|AAA00001| of A, whose lookup
gives gb|AAA00001.1|, an unrelated sequence Y of the same length, first: the PSSM reads Y's letters), `prd_as1` (Y one
letter shorter than A's HSP: NCBI warns and the row is the query's) and its version 4 copy `prd_as1_4`, `prd_ad` (Y as
dbj|AAA00001.1|: read), `prd_ar` (Y as ref|AAA00001.1|, whose ids do not compare equal: the loader keeps none and NCBI
warns), `prd_t2` (A in two HSPs over a short Y: two warnings), from the FASTA of the same names, for the queries
`prd_q.faa`, `prd_q2.faa` (two queries), `prd_q3.faa` (three, the second with an O) and the PSSM query `prd_q.asn.gz`.
Review 1 of that line: `prd_u0y` (Y `gi|777|gb|YYY00001.1|`) and `prd_ua` (`lcl|777` before `gi|777`, so that NCBI's
formatters find no sequence by `gi|777` in its database, and B), read through the alias `prd_unf.pal` (`prd_u0y` first: NCBI reads Y's letters for the
unfetchable subject and does not warn); `prd_zunf0y` (Y `lcl|777` before `gi|777`) and `prd_zunfa` (`gi|777` and B), read
through `prd_zunf.pal` (the lookup gives Y, whose kept definition line does not hold `gi|777`: NCBI warns); `prd_otag`
(`gnl|BL_ORD_ID|0` stored with a text tag, which makes the lookup throw: NCBI warns and goes on) and its version 4 copy
`prd_otag_4`; `prd_wk5` (Y `ref|AAA00001.1|` kept by nothing, C `gb|AAA00001.2|` read first and matched weakly to A) and
`prd_ah4` (Y `lcl|AAA00001.1`, A `gb|AAA00001.1|`, Z `gb|ZZZ00001.1|` with `emb|AAA00001.1|` in its second definition
line: another sequence also holds A's id); `prd_txself` (A itself, `gb|AAA00001.1|`, with taxid 99990 by
`-taxid_map prd_txself.map`, B without; the PSSM's own load reads the taxonomy database).
Review 2 of that line: `prd_xg0y` (Y `gnl|DB|OTHER ^A gnl|DB|SEQA`) and `prd_xga` (A `gnl|DB|SEQA`, B) through
`prd_xg.pal` (NCBI's scope reads Y's second definition line, which holds A's id); `prd_tw0y` (Y
`gb|AAA00002.1|gb|AAA00001.1|`, taxid 99990 by `prd_tw0y.map`) and `prd_twa` (A `gb|AAA00001.1|`, B) through `prd_tw.pal`
(GetBioseq's filter keeps no line of Y: no taxid lookup); `prd_paha` (a copy of `prd_ua`) and `prd_pahc8` (C
`gi|888|gb|CCC00001.1| ^A gi|777`) through `prd_pah.pal` (another holder of `gi|777`; a list is read in the order of
its names); `prd_tu2y` (Y `gi|777|gb|YYY00001.1|` of 20 letters,
taxid 99990 by `prd_tu2y.map`) before `prd_ua` through `prd_tu.pal` (the prefetch reads Y for the unfetchable A); and
the nucleotide `prdn_u0y`/`prdn_ua` (`prdn_u.nal`) and `prdn_z0y`/`prdn_za` (`prdn_z.nal`), the same shapes as
`prd_unf`/`prd_zunf`, for the minus-strand query `prdn_qrc.fna` (BLASTN and TBLASTX prefetch the length of a
minus-strand Seq-align's subject); `prd_mb0ys` (Y `gb|BBB00001.1|`, 100 letters) and `prd_mbab` (A, B `gb|BBB00001|`,
C) through `prd_mb.pal` for `tests/psiblast_pssm_read_order.rs` (the last query's PSSM warning before the epilog).
Review 3 of that line: `prd_np` and its version 4 copy `prd_np_4` from `prd_np.faa` (A and B without `-parse_seqids`,
taxid 99990 by `-taxid 99990`: the PSSM reads them by `gnl|BL_ORD_ID|n` and looks their taxid up), and `prd_tc` (A
`gb|AAA00001.1|` taxid 99990 and B `gb|BBB00001.1|` taxid 150 by `prd_tc.map`, one slot of NCBI's taxonomy cache).
Review 4 of that line (the slot of NCBI's taxonomy cache that 99990 and 150 share, maps of the same names): `prd_m` (A
with the definition lines `gb|AAA00001.1|` taxid 99990, `gb|CCC00001.1|` 150 and `gb|DDD00001.1|` 99990: one load looks
99990 up twice around 150), `prd_m2` (A `gb|AAA00001.1|` 99990 and `gb|CCC00001.1|` 150, then B `gb|BBB00001.1|` 99990)
and `prd_ax` (A `gb|AAA00001.1|` 99990 and X `gb|XXX00001.1|` 150, a random sequence that no PSSM reads).
`--pssmread-only` builds only these.

Where the reports read their subjects one by one (line DBID taxrow, review 3 H2 of pssmread), version 5:
`txr_p` (protein, `txr_p.faa`) and `txr_n` (nucleotide, `txr_n.fna`, the proteins' back-translations), each A
`gb|AAA00001.1|` taxid 9606 (the best hit of `txr_q.faa` / `txr_q.fna`) and B `gb|BBB00001.1|` taxid 99990, whose lookup
writes the offset error of `tests/seqdb/taxdb/offset_last_past_end` (`txr.map`), for `tests/taxonomy_read_order.rs`.
`--taxrow-only` builds only these.

Copies of one record that follow each other in a hit list (line DBID N3; `make_join_fasta.py`), all version 5:
`pjoin_a.*`, `pjoin_b.*`, `pjoin_c.*` and `njoin_a.*`, `njoin_b.*`, `njoin_c.*` from the FASTA of the same names, read
through `pjoin_abc.pal` and `njoin_abc.nal`; `pjoin_bt.*` and `njoin_bt.*` from `pjoin_b.faa` and `njoin_b.fna` with
taxid 9606, read in place of the b databases through `pjoin_abt.pal` and `njoin_abt.nal` (copies with other taxids);
`pjoin_x.*` and `njoin_x.*` from the FASTA of the same names (the definition lines of r1 and r3 with other letters),
read after `pjoin_a` and `njoin_a` through `pjoin_ax.pal` and `njoin_ax.nal`; `njoin_ta.*` and `njoin_tb.*` from
`njoin_t.fna` with taxid 9606, through `njoin_tab.nal`. `--join-only` builds only these.
"""
import argparse
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEQID = HERE.parents[1] / "seqid" / "fasta"
DROPPED_VARIANTS = (
    "pdb|1ABC|AA",
    "pdb|1ABC|",
    "sp|P99999|",
    "ref|NP_000002|",
    "bbs|123",
    "gnl|PIG|123",
    "nat|NA000001.1|",
)


def records(path: Path) -> list[str]:
    return [">" + record for record in path.read_text().split(">")[1:]]


def build_weak(args: argparse.Namespace) -> None:
    for name, dbtype, extension in (("pweak", "prot", "faa"), ("pweak_r", "prot", "faa"), ("nweak", "nucl", "fna")):
        for version in (5, 4):
            subprocess.run(
                [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / f"{name}.{extension}"),
                 "-dbtype", dbtype, "-parse_seqids", "-blastdb_version", str(version),
                 "-title", f"{name}_v{version} (parsed Seq-ids)", "-out", str(args.out / f"{name}_v{version}")],
                check=True, stdout=subprocess.DEVNULL)


def build_copytaxid(args: argparse.Namespace) -> None:
    for name, fasta, dbtype, taxid in (("pdup_t", "pdup.faa", "prot", "9606"), ("pdup_u", "pdup.faa", "prot", None),
                                       ("pdup_n", "pdup.faa", "prot", "99990"), ("pdup_o", "pdup.faa", "prot", "99990"),
                                       ("pdupe_a", "pdupe.faa", "prot", None), ("pdupe_t", "pdupe.faa", "prot", "9606"),
                                       ("ndup_a", "ndup.fna", "nucl", None), ("ndup_t", "ndup.fna", "nucl", "9606")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta), "-dbtype", dbtype, "-parse_seqids",
             "-blastdb_version", "5", "-title", name, "-out", str(args.out / name)]
            + (["-taxid", taxid] if taxid else []),
            check=True, stdout=subprocess.DEVNULL)
    for alias, names, extension in (("pdup_at", "pdup_a pdup_t", "pal"), ("pdup_tu", "pdup_t pdup_u", "pal"),
                                    ("pdup_no", "pdup_n pdup_o", "pal"), ("pdupe_at", "pdupe_a pdupe_t", "pal"),
                                    ("ndup_at", "ndup_a ndup_t", "nal")):
        (args.out / f"{alias}.{extension}").write_text(
            f"#\n# Alias file of line DBID, N1b\n#\nTITLE {alias} (two databases, parsed Seq-ids)\n"
            f"DBLIST {names}\n")


def build_titlesrc(args: argparse.Namespace) -> None:
    for name, fasta, dbtype, taxid in (("ndupe_a", "ndupe.fna", "nucl", None), ("ndupe_t", "ndupe.fna", "nucl", "9606"),
                                       ("nsubj_t", "nsubj.fna", "nucl", "9606"),
                                       ("nsubjg_t", "nsubjg.fna", "nucl", "9606"),
                                       ("psubj_a", "psubj.faa", "prot", None), ("psubj_t", "psubj.faa", "prot", "9606"),
                                       ("psubjb_a", "psubjb.faa", "prot", None),
                                       ("psubjb_t", "psubjb.faa", "prot", "9606")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta), "-dbtype", dbtype, "-parse_seqids",
             "-blastdb_version", "5", "-title", name, "-out", str(args.out / name)]
            + (["-taxid", taxid] if taxid else []),
            check=True, stdout=subprocess.DEVNULL)
    for alias, names, extension in (("ndupe_at", "ndupe_a ndupe_t", "nal"), ("psubj_at", "psubj_a psubj_t", "pal"),
                                    ("psubjb_at", "psubjb_a psubjb_t", "pal")):
        (args.out / f"{alias}.{extension}").write_text(
            f"#\n# Alias file of line DBID, N1b L1\n#\nTITLE {alias} (two databases, parsed Seq-ids)\n"
            f"DBLIST {names}\n")


def build_emptytitle(args: argparse.Namespace) -> None:
    for name, fasta, dbtype in (("nbare_t", "nbare.fna", "nucl"), ("pbare_t", "pbare.faa", "prot")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta), "-dbtype", dbtype, "-blastdb_version", "5",
             "-taxid", "9606", "-title", name, "-out", str(args.out / name)],
            check=True, stdout=subprocess.DEVNULL)


def build_tabfields(args: argparse.Namespace) -> None:
    subprocess.run(
        [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / "ppunct.faa"), "-dbtype", "prot", "-blastdb_version", "4",
         "-title", "ppunct_np", "-out", str(args.out / "ppunct_np")],
        check=True, stdout=subprocess.DEVNULL)


def build_pfbatch(args: argparse.Namespace) -> None:
    subprocess.run(
        [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / "nordtag.fna"), "-dbtype", "nucl", "-parse_seqids",
         "-blastdb_version", "5", "-title", "nordtag_v5 (parsed Seq-ids)", "-out", str(args.out / "nordtag_v5")],
        check=True, stdout=subprocess.DEVNULL)


def build_pfwarn(args: argparse.Namespace) -> None:
    for name, fasta, dbtype, version in (("ntag_t5", "nordtag.fna", "nucl", "5"), ("nlook_t4", "nlook.fna", "nucl", "4"),
                                         ("pshare_t5", "pshare.faa", "prot", "5")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta), "-dbtype", dbtype, "-parse_seqids",
             "-blastdb_version", version, "-taxid", "99990", "-title", f"{name} (parsed Seq-ids, taxid 99990)",
             "-out", str(args.out / name)],
            check=True, stdout=subprocess.DEVNULL)
    subprocess.run(
        [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / "pfw_weak.faa"), "-dbtype", "prot", "-parse_seqids",
         "-blastdb_version", "5", "-title", "pfw_weak_v5 (parsed Seq-ids)", "-out", str(args.out / "pfw_weak_v5")],
        check=True, stdout=subprocess.DEVNULL)
    bad = args.out / "pfw_bad"
    bad.mkdir(exist_ok=True)
    for path in sorted(args.out.glob("ntag_t5.*")):
        data = bytearray(path.read_bytes())
        if path.suffix == ".ndb":
            # The flags of an LMDB page header (8-byte number, 2 bytes of padding, 2 of flags).
            flags = 6 * 4096 + 10
            assert data[flags] != 0, "page 6 of ntag_t5.ndb is no longer the accession table's leaf"
            data[flags] = 0
        (bad / path.name).write_bytes(bytes(data))


def build_weakbatch(args: argparse.Namespace) -> None:
    for name, fasta in (("wkb_e3_v5", "wkb_e3.faa"), ("wkb_p1_v5", "wkb_p1.faa")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta), "-dbtype", "prot", "-parse_seqids",
             "-blastdb_version", "5", "-title", f"{name} (parsed Seq-ids)", "-out", str(args.out / name)],
            check=True, stdout=subprocess.DEVNULL)


def build_pssmrows(args: argparse.Namespace) -> None:
    for name in ("psr_cpa", "psr_cpb", "psr_cbx", "psr_jm", "psr_jn", "psr_jx1", "psr_jx2", "psr_jctl", "psr_self",
                 "psr_cu", "psr_cl", "psr_ce", "psr_ka", "psr_kb", "psr_gim", "psr_gz", "psr_gimn", "psr_jv", "psr_vm"):
        nucleotides = name == "psr_gimn"
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / f"{name}.{'fna' if nucleotides else 'faa'}"),
             "-dbtype", "nucl" if nucleotides else "prot", "-parse_seqids",
             "-blastdb_version", "5", "-title", f"{name} (parsed Seq-ids)", "-out", str(args.out / name)],
            check=True, stdout=subprocess.DEVNULL)
    for alias, names in (("psr_cp", "psr_cpa psr_cpb"), ("psr_cbxp", "psr_cbx psr_cpb"), ("psr_jx", "psr_jx1 psr_jx2"),
                         ("psr_cul", "psr_cu psr_cl"), ("psr_cue", "psr_cu psr_ce"), ("psr_k", "psr_ka psr_kb")):
        (args.out / f"{alias}.pal").write_text(
            f"#\n# Alias file of line DBID, pssmrows\n#\nTITLE {alias} (two databases, parsed Seq-ids)\n"
            f"DBLIST {names}\n")


def build_pssmread(args: argparse.Namespace) -> None:
    for name, fasta, version in (("prd_av", "prd_av.faa", "5"), ("prd_as1", "prd_as1.faa", "5"),
                                 ("prd_ad", "prd_ad.faa", "5"), ("prd_ar", "prd_ar.faa", "5"),
                                 ("prd_t2", "prd_t2.faa", "5"), ("prd_as1_4", "prd_as1.faa", "4"),
                                 ("prd_u0y", "prd_u0y.faa", "5"), ("prd_ua", "prd_ua.faa", "5"),
                                 ("prd_zunf0y", "prd_zunf0y.faa", "5"), ("prd_zunfa", "prd_zunfa.faa", "5"),
                                 ("prd_otag", "prd_otag.faa", "5"), ("prd_otag_4", "prd_otag.faa", "4"),
                                 ("prd_wk5", "prd_wk5.faa", "5"), ("prd_ah4", "prd_ah4.faa", "5"),
                                 ("prd_txself", "prd_txself.faa", "5"),
                                 ("prd_xg0y", "prd_xg0y.faa", "5"), ("prd_xga", "prd_xga.faa", "5"),
                                 ("prd_tw0y", "prd_tw0y.faa", "5"), ("prd_twa", "prd_twa.faa", "5"),
                                 ("prd_paha", "prd_paha.faa", "5"), ("prd_pahc8", "prd_pahc8.faa", "5"),
                                 ("prd_tu2y", "prd_tu2y.faa", "5"), ("prd_mb0ys", "prd_mb0ys.faa", "5"),
                                 ("prd_mbab", "prd_mbab.faa", "5"),
                                 ("prdn_u0y", "prdn_u0y.fna", "5"), ("prdn_ua", "prdn_ua.fna", "5"),
                                 ("prdn_z0y", "prdn_z0y.fna", "5"), ("prdn_za", "prdn_za.fna", "5"),
                                 ("prd_tc", "prd_tc.faa", "5"), ("prd_m", "prd_m.faa", "5"),
                                 ("prd_m2", "prd_m2.faa", "5"), ("prd_ax", "prd_ax.faa", "5")):
        maps = {"prd_txself": "prd_txself.map", "prd_tw0y": "prd_tw0y.map", "prd_tu2y": "prd_tu2y.map",
                "prd_tc": "prd_tc.map", "prd_m": "prd_m.map", "prd_m2": "prd_m2.map", "prd_ax": "prd_ax.map"}
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta),
             "-dbtype", "nucl" if fasta.endswith(".fna") else "prot", "-parse_seqids",
             "-blastdb_version", version, "-title", f"{name} (parsed Seq-ids)", "-out", str(args.out / name)]
            + (["-taxid_map", str(HERE / maps[name])] if name in maps else []),
            check=True, stdout=subprocess.DEVNULL)
    # Without parsed Seq-ids: the PSSM reads by the ordinal id (review 3 of line DBID pssmread, H1).
    for name, version in (("prd_np", "5"), ("prd_np_4", "4")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / "prd_np.faa"), "-dbtype", "prot",
             "-taxid", "99990", "-blastdb_version", version, "-title", f"{name} (no parsed Seq-ids)",
             "-out", str(args.out / name)],
            check=True, stdout=subprocess.DEVNULL)
    for alias, names, extension in (("prd_unf", "prd_u0y prd_ua", "pal"), ("prd_zunf", "prd_zunf0y prd_zunfa", "pal"),
                                    ("prd_xg", "prd_xg0y prd_xga", "pal"), ("prd_tw", "prd_tw0y prd_twa", "pal"),
                                    ("prd_pah", "prd_paha prd_pahc8", "pal"), ("prd_tu", "prd_tu2y prd_ua", "pal"),
                                    ("prd_mb", "prd_mb0ys prd_mbab", "pal"),
                                    ("prdn_u", "prdn_u0y prdn_ua", "nal"), ("prdn_z", "prdn_z0y prdn_za", "nal")):
        (args.out / f"{alias}.{extension}").write_text(
            f"#\n# Alias file of line DBID, pssmread\n#\nTITLE {alias} (two databases, parsed Seq-ids)\n"
            f"DBLIST {names}\n")


def build_taxrow(args: argparse.Namespace) -> None:
    for name, fasta in (("txr_p", "txr_p.faa"), ("txr_n", "txr_n.fna")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta),
             "-dbtype", "nucl" if fasta.endswith(".fna") else "prot", "-parse_seqids",
             "-blastdb_version", "5", "-title", f"{name} (parsed Seq-ids)", "-out", str(args.out / name),
             "-taxid_map", str(HERE / "txr.map")],
            check=True, stdout=subprocess.DEVNULL)


def build_join(args: argparse.Namespace) -> None:
    for name, fasta, dbtype, taxid in (("pjoin_a", "pjoin_a.faa", "prot", None), ("pjoin_b", "pjoin_b.faa", "prot", None),
                                       ("pjoin_c", "pjoin_c.faa", "prot", None), ("pjoin_bt", "pjoin_b.faa", "prot", "9606"),
                                       ("njoin_a", "njoin_a.fna", "nucl", None), ("njoin_b", "njoin_b.fna", "nucl", None),
                                       ("njoin_c", "njoin_c.fna", "nucl", None), ("njoin_bt", "njoin_b.fna", "nucl", "9606"),
                                       ("pjoin_x", "pjoin_x.faa", "prot", None), ("njoin_x", "njoin_x.fna", "nucl", None),
                                       ("njoin_ta", "njoin_t.fna", "nucl", "9606"), ("njoin_tb", "njoin_t.fna", "nucl", "9606")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta), "-dbtype", dbtype, "-parse_seqids",
             "-blastdb_version", "5", "-title", name, "-out", str(args.out / name)]
            + (["-taxid", taxid] if taxid else []),
            check=True, stdout=subprocess.DEVNULL)
    for alias, names, extension in (("pjoin_abc", "pjoin_a pjoin_b pjoin_c", "pal"),
                                    ("pjoin_abt", "pjoin_a pjoin_bt pjoin_c", "pal"),
                                    ("njoin_abc", "njoin_a njoin_b njoin_c", "nal"),
                                    ("njoin_abt", "njoin_a njoin_bt njoin_c", "nal"),
                                    ("pjoin_ax", "pjoin_a pjoin_x", "pal"), ("njoin_ax", "njoin_a njoin_x", "nal"),
                                    ("njoin_tab", "njoin_ta njoin_tb", "nal")):
        count = {2: "two", 3: "three"}[len(names.split())]
        (args.out / f"{alias}.{extension}").write_text(
            f"#\n# Alias file of line DBID, N3\n#\nTITLE {alias} ({count} databases, parsed Seq-ids)\n"
            f"DBLIST {names}\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ncbi-bin", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=HERE)
    ap.add_argument("--weak-only", action="store_true", help="build only the databases of line DBID N2c")
    ap.add_argument("--copytaxid-only", action="store_true", help="build only the databases of line DBID N1b")
    ap.add_argument("--join-only", action="store_true", help="build only the databases of line DBID N3")
    ap.add_argument("--titlesrc-only", action="store_true", help="build only the databases of line DBID N1b L1")
    ap.add_argument("--emptytitle-only", action="store_true", help="build only the databases of line DBID H3")
    ap.add_argument("--tabfields-only", action="store_true", help="build only the database of line DBID tabfields")
    ap.add_argument("--pfbatch-only", action="store_true", help="build only the database of line DBID pfbatch")
    ap.add_argument("--pfwarn-only", action="store_true", help="build only the databases of line DBID pfwarn")
    ap.add_argument("--weakbatch-only", action="store_true", help="build only the databases of line DBID weakbatch")
    ap.add_argument("--pssmrows-only", action="store_true", help="build only the databases of line DBID pssmrows")
    ap.add_argument("--pssmread-only", action="store_true", help="build only the databases of line DBID pssmread")
    ap.add_argument("--taxrow-only", action="store_true", help="build only the databases of line DBID taxrow")
    args = ap.parse_args()
    if args.weak_only:
        build_weak(args)
        return 0
    if args.copytaxid_only:
        build_copytaxid(args)
        return 0
    if args.join_only:
        build_join(args)
        return 0
    if args.titlesrc_only:
        build_titlesrc(args)
        return 0
    if args.emptytitle_only:
        build_emptytitle(args)
        return 0
    if args.tabfields_only:
        build_tabfields(args)
        return 0
    if args.pfbatch_only:
        build_pfbatch(args)
        return 0
    if args.pfwarn_only:
        build_pfwarn(args)
        return 0
    if args.weakbatch_only:
        build_weakbatch(args)
        return 0
    if args.pssmrows_only:
        build_pssmrows(args)
        return 0
    if args.pssmread_only:
        build_pssmread(args)
        return 0
    if args.taxrow_only:
        build_taxrow(args)
        return 0
    variants = [r for r in records(SEQID / "variants.fa") if r[1:].split()[0] not in DROPPED_VARIANTS]
    all_records = (records(SEQID / "prot_ids.fa") + records(SEQID / "prot_extra.fa") + variants
                   + records(HERE / "gen.faa"))
    work = args.out / "work"
    work.mkdir(exist_ok=True)
    for version, kept in (
        (4, all_records),
        (5, [r for r in all_records if not r.startswith((">gi|1414 ", ">lcl|123 "))]),
    ):
        fasta = work / f"prot_v{version}.fa"
        fasta.write_text("".join(kept))
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(fasta), "-dbtype", "prot", "-parse_seqids",
             "-blastdb_version", str(version), "-title", f"prot_v{version} (parsed Seq-ids)",
             "-taxid_map", str(SEQID / "taxmap_prot.txt"), "-out", str(args.out / f"prot_v{version}")],
            check=True, stdout=subprocess.DEVNULL)
        fasta.unlink()
    work.rmdir()
    subprocess.run(
        [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / "refuse.faa"), "-dbtype", "prot",
         "-parse_seqids", "-blastdb_version", "5", "-title", "refuse_v5 (parsed Seq-ids)",
         "-taxid_map", str(HERE / "taxmap_refuse.txt"), "-out", str(args.out / "refuse_v5")],
        check=True, stdout=subprocess.DEVNULL)
    nucl = records(HERE / "nucl.fna")
    work.mkdir(exist_ok=True)
    for version, kept in ((4, nucl), (5, [r for r in nucl if not r.startswith(">gi|1111 ")])):
        fasta = work / f"nucl_v{version}.fna"
        fasta.write_text("".join(kept))
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(fasta), "-dbtype", "nucl", "-parse_seqids",
             "-blastdb_version", str(version), "-title", f"nucl_v{version} (parsed Seq-ids)",
             "-taxid_map", str(HERE / "taxmap_nucl.txt"), "-out", str(args.out / f"nucl_v{version}")],
            check=True, stdout=subprocess.DEVNULL)
        fasta.unlink()
    work.rmdir()
    for name, version, title in (("nrefuse_v5", 5, "nrefuse_v5 (parsed Seq-ids)"),
                                 ("njoined_a", 4, "njoined_a"), ("njoined_b", 4, "njoined_b"),
                                 ("nambig_v4", 4, "nambig_v4")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / (name.split("_v")[0] + ".fna")),
             "-dbtype", "nucl", "-parse_seqids", "-blastdb_version", str(version), "-title", title,
             "-out", str(args.out / name)],
            check=True, stdout=subprocess.DEVNULL)
    for name, dbtype, extension in (("pfound", "prot", "faa"), ("nfound", "nucl", "fna")):
        for version in (5, 4):
            subprocess.run(
                [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / f"{name}.{extension}"),
                 "-dbtype", dbtype, "-parse_seqids", "-blastdb_version", str(version),
                 "-title", f"{name}_v{version} (parsed Seq-ids)", "-out", str(args.out / f"{name}_v{version}")],
                check=True, stdout=subprocess.DEVNULL)
    for name, fasta in (("pdup_a", "pdup.faa"), ("pdup_b", "pdup.faa"), ("pdup_c", "pdup_c.faa")):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta), "-dbtype", "prot", "-parse_seqids",
             "-blastdb_version", "5", "-title", name, "-out", str(args.out / name)],
            check=True, stdout=subprocess.DEVNULL)
    for name, dbtype, extension in (("pord", "prot", "faa"), ("nord", "nucl", "fna")):
        for version in (5, 4):
            subprocess.run(
                [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / f"{name}.{extension}"),
                 "-dbtype", dbtype, "-parse_seqids", "-blastdb_version", str(version),
                 "-title", f"{name}_v{version} (parsed Seq-ids)", "-out", str(args.out / f"{name}_v{version}")],
                check=True, stdout=subprocess.DEVNULL)
    subprocess.run(
        [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / "pord_x.faa"), "-dbtype", "prot", "-parse_seqids",
         "-blastdb_version", "5", "-title", "pord_x", "-out", str(args.out / "pord_x")],
        check=True, stdout=subprocess.DEVNULL)
    for name, dbtype, extension in (("plook", "prot", "faa"), ("nlook", "nucl", "fna")):
        for version in (5, 4):
            subprocess.run(
                [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / f"{name}.{extension}"),
                 "-dbtype", dbtype, "-parse_seqids", "-blastdb_version", str(version),
                 "-title", f"{name}_v{version} (parsed Seq-ids)", "-out", str(args.out / f"{name}_v{version}")],
                check=True, stdout=subprocess.DEVNULL)
    subprocess.run(
        [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / "plook_z.faa"), "-dbtype", "prot", "-parse_seqids",
         "-blastdb_version", "5", "-title", "plook_z", "-out", str(args.out / "plook_z")],
        check=True, stdout=subprocess.DEVNULL)
    for name, fasta, dbtype, versions in (("pshare", "pshare.faa", "prot", (5, 4)),
                                          ("nshare", "nshare.fna", "nucl", (5,))):
        for version in versions:
            subprocess.run(
                [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / fasta), "-dbtype", dbtype, "-parse_seqids",
                 "-blastdb_version", str(version), "-title", f"{name}_v{version} (parsed Seq-ids)",
                 "-out", str(args.out / f"{name}_v{version}")],
                check=True, stdout=subprocess.DEVNULL)
    for name in ("pshare_a", "pshare_b", "pshare_c", "pshare_d"):
        subprocess.run(
            [str(args.ncbi_bin / "makeblastdb"), "-in", str(HERE / f"{name}.faa"), "-dbtype", "prot", "-parse_seqids",
             "-blastdb_version", "5", "-title", name, "-out", str(args.out / name)],
            check=True, stdout=subprocess.DEVNULL)
    for alias, names in (("pshare_ab", "pshare_a pshare_b"), ("pshare_cd", "pshare_c pshare_d"),
                         ("pdup_abc", "pdup_a pdup_b pdup_c")):
        count = "three" if alias == "pdup_abc" else "two"
        (args.out / f"{alias}.pal").write_text(
            f"#\n# Alias file of line DBID, N2b\n#\nTITLE {alias} ({count} databases, parsed Seq-ids)\n"
            f"DBLIST {names}\n")
    (args.out / "plook_vz.pal").write_text(
        "#\n# Alias file of line DBID, N2b\n#\nTITLE plook_vz (two databases, parsed Seq-ids)\n"
        "DBLIST plook_v5 plook_z\n")
    (args.out / "pord_vx.pal").write_text(
        "#\n# Alias file of line DBID, N2a\n#\nTITLE pord_vx (two databases, parsed Seq-ids)\n"
        "DBLIST pord_x pord_v5\n")
    for alias, second in (("pdup_ab", "pdup_b"), ("pdup_ac", "pdup_c")):
        (args.out / f"{alias}.pal").write_text(
            f"#\n# Alias file of line DBID, N1\n#\nTITLE {alias} (two databases, parsed Seq-ids)\n"
            f"DBLIST pdup_a {second}\n")
    (args.out / "njoined.nal").write_text(
        "#\n# Alias file of line DBID, ID2b\n#\nTITLE njoined (two volumes, parsed Seq-ids)\n"
        "DBLIST njoined_a njoined_b\n")
    build_weak(args)
    build_copytaxid(args)
    build_join(args)
    build_titlesrc(args)
    build_emptytitle(args)
    build_tabfields(args)
    build_pfbatch(args)
    build_pfwarn(args)
    build_weakbatch(args)
    build_pssmrows(args)
    build_pssmread(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
