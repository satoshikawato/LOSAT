#!/usr/bin/env python3
"""Builds tests/db_search/paja_taxid.{phr,pin,psq}: the paja proteins with taxid 10090 in every
definition line, except BDT63547.1 with taxid 9606 (the only subject that the queries of
bpqb_q.faa hit in their second NCBI query batch with a different taxid: its hit comes after
the hits of the first batch, so a refusal of its taxonomy name would come after the first
block's bytes). OF-BPQB2 (outfmt 15 and 16 in the blastp block search).

Usage: make_paja_taxid.py MAKEBLASTDB OUTDIR   (OUTDIR/paja_taxid.* ; .pjs is dropped)

`makeblastdb -taxid 10090` writes the taxid as `a2 80 02 02 27 6a 00 00`; the record of
BDT63547.1 gets 9606 (`25 86`) in place of 10090 (`27 6a`): the same length.
"""
import os
import subprocess
import sys

TITLE = "PajaWSV proteins, taxid 10090 (BDT63547.1: 9606)"
SUBJECT = b"BDT63547.1"
OLD = bytes.fromhex("0202276a")
NEW = bytes.fromhex("02022586")


def main(makeblastdb, outdir):
    here = os.path.dirname(os.path.abspath(__file__))
    source = os.path.join(here, "..", "fasta", "PajaWSV.faa")
    out = os.path.join(outdir, "paja_taxid")
    subprocess.run(
        [makeblastdb, "-in", source, "-dbtype", "prot", "-out", out, "-title", TITLE,
         "-blastdb_version", "4", "-taxid", "10090"],
        check=True, stdout=subprocess.DEVNULL)
    os.remove(out + ".pjs")
    with open(out + ".phr", "rb") as f:
        data = f.read()
    at = data.index(SUBJECT)
    at = data.index(OLD, at)
    # the taxid of this definition line comes before the next record's title
    following = data.find(b"BDT6", at)
    assert following == -1 or following - at > len(OLD)
    data = data[:at] + NEW + data[at + len(OLD):]
    with open(out + ".phr", "wb") as f:
        f.write(data)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
