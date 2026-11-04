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
(`>` alone), without -parse_seqids, with taxid 9606. `--emptytitle-only` builds only these.

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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
