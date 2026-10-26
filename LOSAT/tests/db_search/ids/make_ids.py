#!/usr/bin/env python3
"""Build the protein databases with parsed Seq-ids of the `db.id.*` rows (line DBID, ID2a).

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
the other sequence) and without the three records that no NCBI formatter finds (`bbs|`,
`gnl|PIG|`, `nat|`), and `gen.faa` (records whose titles NCBI makes from their Seq-ids: `tpg`
and `tpe` ids with and without titles, `pdb` ids with and without a chain, a patent, RefSeq
NC_/NP_/WP_ accessions without titles, `lcl|Subject_` ids beside other ids). Taxids:
`LOSAT/tests/seqid/fasta/taxmap_prot.txt`.

Outputs: `prot_v4.*` (version 4, every record) and `prot_v5.*` (version 5, without the GI-only
record `gi|1414`, after which NCBI 2.17.0's makeblastdb aborts writing the LMDB file, and without
`lcl|123`, a local integer id that NCBI's formatters do not find in an LMDB file: LOSAT refuses
it there, `report::db_subject::ParsedSubject::fetched`).

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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ncbi-bin", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=HERE)
    args = ap.parse_args()
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
