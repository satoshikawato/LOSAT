#!/usr/bin/env python3
"""Write the seqid list files of the rows db.sl.* (line DI-8b: -seqidlist, -negative_seqidlist).

Text lists hold the forms NCBI's reader takes (`SeqDB_ReadMemorySiList`: white-space separated
tokens, a leading `>` skipped, `#` to the end of the line a comment, CRLF lines) for the
databases tests/db_search/taxfilter/tfp (Swiss-Prot accessions), taxfilter/tfn (Ensembl / FlyBase
local ids in two volumes), ids/prot_v5 and ids/nucl_v5 (every Seq-id form; also read on the
version 4 databases ids/prot_v4 and ids/nucl_v4).

Binary lists (`.bsl`) are written here in the format of `blastdb_aliastool -seqid_file_out`
(`CSeqidlistRead`, c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:58-108): NUL, the
file size (u64 LE), the number of ids (u64), the title (u32 length + bytes), the creation date
(u8 length + bytes), the database length (u64; when not 0: the database's date, u8 length +
bytes, and its volume names, u32 length + bytes), then each id as a u8 length (0xFF: a u32
length follows) and its bytes. The ids are the processed forms (`s_ProcessSeqIDsForV5`: the
upper-cased accession with its version when given), sorted. `sl_tfp_badlen.bsl` records a
database length that is not tfp's: NCBI warns `Seqidlist file db info does not match input db`.

usage: make_seqidlist.py [--out LOSAT/tests/db_search/seqidlist]
"""
import argparse
import struct
from pathlib import Path

TEXT = {
    # tfp: human HBA (versioned), mouse (versionless), horse in FASTA form, two ids on a CRLF line, an unknown id
    "sl_tfp.txt": ">P69905.2\n# human; mouse without a version; horse in FASTA form; two per line\nP01942\nsp|P01935.1|\nQ9TS35 P63107.2\r\nUNKNOWN1\n",
    # the same ids without versions: a versionless id excludes the versioned stored id
    "sl_tfp_neg.txt": "P69905\nP01942.2\nP01935\nQ9TS35\nP63107\n",
    "sl_tfn.txt": "ENST00000320868.9\nENST00000251595.11\nENSDART00000113098.5\nFBtr0083245\n",
    # prot_v5: tpg, RefSeq and Swiss-Prot without versions, a pdb chain, local ids, a tpg beside a local id, a
    # general id as db:tag, a RefSeq without version, an unknown id, a gi (dropped)
    "sl_prot.txt": ">DAA10001.1\n# RefSeq and Swiss-Prot without versions\nNP_100006 P00012\n2XYZ_B\nmyprot9\nDAA10008.1\tmydb:prot10\nXP_000011\nUNKNOWN1\ngi|1414\n",
    # nucl_v5: a versionless RefSeq of an OID of three definition lines, a gb with its version, local ids, a
    # general id, a bare RefSeq, a versionless DDBJ, a gi (dropped), a versionless gb of a second definition line
    "sl_nucl.txt": "NM_000001\nAF000002.3 mynuc6\nmydb:nuc7\nNM_000008.2\nAB000009\ngi|1111\nSubject_30\nAF000020\n",
    "sl_gi.txt": "gi|654321\ngi|1111\n",
    "sl_comment.txt": "# a comment only\n\n",
    "sl_none.txt": "UNKNOWN1\nUNKNOWN2\n",
    "sl_empty.txt": "",
}


def binary(ids: list[bytes], title: bytes, db_length: int = 0, db_date: bytes = b"", db_vols: bytes = b"") -> bytes:
    body = struct.pack("<Q", len(ids)) + struct.pack("<I", len(title)) + title
    date = b"Oct 30, 2026  3:00 PM"
    body += struct.pack("<B", len(date)) + date + struct.pack("<Q", db_length)
    if db_length:
        body += struct.pack("<B", len(db_date)) + db_date + struct.pack("<I", len(db_vols)) + db_vols
    for sid in sorted(ids):
        body += (b"\xff" + struct.pack("<I", len(sid)) if len(sid) >= 0xFF else struct.pack("<B", len(sid))) + sid
    size = 1 + 8 + len(body)
    return b"\0" + struct.pack("<Q", size) + body


BINARY = {
    "sl_tfp.bsl": binary([b"P69905.2", b"P01942", b"P01935.1", b"Q9TS35", b"P63107.2"], b"sl_tfp"),
    "sl_tfp_badlen.bsl": binary([b"P69905.2", b"P01942", b"P01935.1", b"Q9TS35", b"P63107.2"], b"sl_tfp", 1, b"Oct 30, 2026  2:00 PM", b"tfp"),
    "sl_nucl.bsl": binary([b"NM_000001", b"AF000002.3", b"mynuc6", b"mydb:nuc7", b"NM_000008.2", b"AB000009", b"Subject_30", b"AF000020"], b"sl_nucl"),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, text in TEXT.items():
        (args.out / name).write_bytes(text.encode())
    for name, data in BINARY.items():
        (args.out / name).write_bytes(data)


if __name__ == "__main__":
    main()
