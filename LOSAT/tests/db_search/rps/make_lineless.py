#!/usr/bin/env python3
"""The profile databases with definition lines without a Seq-id ("lineless" lines; line RG, series D, rps-lineless,
2026-11-27), and their queries.

No NCBI counterpart: test inputs for build.sh (`build.sh - OUT_DIR lineless`). makeprofiledb writes one definition line
with a Seq-id per profile, so these are the committed headers of this directory rewritten: the `.phr` and the header
offsets of the `.pin` of rps_mini, rps_mini3v.00, rps_mini3v.02 and tx3v00 (inflated from their .gz files here), with
make_idlist.py's `drop_seqids` (its source run without its `main()` call, which would rebuild the idlist fixtures).
The other files of each database are the committed ones, so a row replaces only these two
(`rps-gz:rps_mini!phr=ln_a.phr.gz!pin=ln_a.pin.gz`). The rewrite reads only committed bytes: a rebuild gives the
committed files (after `gzip -9 -n`).

Usage: make_lineless.py OUT_DIR   (writes OUT_DIR/ln_*.phr, ln_*.pin and OUT_DIR/ln*.faa, ln*.fna)

rps_mini: OID 0 = gnl|CDD|101 (f101 TK-TMK ...), OID 1 = gnl|CDD|102 (f102 ribonucleotide reductase ...), OID 2 =
gnl|CDD|103 (f103 DNA polymerase ...), one line each. The variants:
  ln_a    OID 0: its line, then a copy without the Seq-id titled "f101x extra lineless"
  ln_a0   OID 0: its line, then a copy without the Seq-id and without a title
  ln_b    OID 1: a copy without the Seq-id titled "f102x lead lineless", then its line
  ln_b0   OID 1: a copy without the Seq-id and without a title, then its line
  ln_c    OID 2: its only line without its Seq-id (no Seq-align)
  ln_abc  the edits of ln_a, ln_b and ln_c in one database
rps_mini3v (three volumes of one profile each, searched by the alias rps_mini3v):
  ln3v_a00  rps_mini3v.00 (profile 101) edited as ln_a
  ln3v_c02  rps_mini3v.02 (profile 103) edited as ln_c
tx3v00 (rps_mini3v.00 made with makeprofiledb -taxid 9606), searched alone as rps_mini3v.00 with a taxonomy database:
  lntx_a  edited as ln_a
  lntx_b  its line after a copy without the Seq-id titled "f101y lead lineless"

Queries (the fragments of build.sh, ids q101 ...; `.fna`, the rpstblastn queries, back-translated with make_idlist.py's
codons: ln101, ln101d, ln102 and ln2b):
  ln101  q101 (profile 101)          ln102  q102 (profile 102)
  ln101d q101d: f101, 30 G, f101 (two HSPs on profile 101)
  ln31   q31: f103, 30 G, f101 (profiles 103 and 101)
  lnall  q101, q102, q103            ln2b   q102, then q101 (two query batches with BATCH_SIZE=1)
"""

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOSAT = HERE.parents[2]

# make_idlist.py's helpers, without its main() call (the file ends with a bare `main()`).
SOURCE = HERE.parent / "make_idlist.py"
text = SOURCE.read_text()
cut = text.rindex("\nmain()")
assert text[cut:].strip() == "main()", "unexpected end of make_idlist.py"
helpers = {"__file__": str(SOURCE), "__name__": "make_idlist_helpers"}
exec(compile(text[:cut] + "\n", str(SOURCE), "exec"), helpers)
drop_seqids, back_translate = helpers["drop_seqids"], helpers["back_translate"]

# (output name, committed volume, its .phr/.pin source names, drop_seqids arguments: lines, titles, copies)
EDIT_A = ({0: {1}}, {0: {1: b"f101x extra lineless"}}, {0: [(1, 0)]})
EDIT_B = ({1: {0}}, {1: {0: b"f102x lead lineless"}}, {1: [(0, 0)]})
EDIT_C = ({2: {0}}, None, None)
VARIANTS = [
    ("ln_a", "rps_mini", EDIT_A),
    ("ln_a0", "rps_mini", ({0: {1}}, {0: {1: None}}, {0: [(1, 0)]})),
    ("ln_b", "rps_mini", EDIT_B),
    ("ln_b0", "rps_mini", ({1: {0}}, {1: {0: None}}, {1: [(0, 0)]})),
    ("ln_c", "rps_mini", EDIT_C),
    ("ln_abc", "rps_mini", ({**EDIT_A[0], **EDIT_B[0], **EDIT_C[0]}, {**EDIT_A[1], **EDIT_B[1]},
                            {**EDIT_A[2], **EDIT_B[2]})),
    ("ln3v_a00", "rps_mini3v.00", EDIT_A),
    ("ln3v_c02", "rps_mini3v.02", ({0: {0}}, None, None)),
    ("lntx_a", "tx3v00", EDIT_A),
    ("lntx_b", "tx3v00", ({0: {0}}, {0: {0: b"f101y lead lineless"}}, {0: [(0, 0)]})),
]

# The fragments of build.sh: (CDD tag, WSSV record, first residue, last residue), 1-based.
FRAGMENTS = {"101": ("YP_009220601.1", 1, 190), "102": ("YP_009220529.1", 301, 480),
             "103": ("YP_009220649.1", 1101, 1290)}
LINKER = "G" * 30
NUCLEOTIDE_QUERIES = {"ln101", "ln101d", "ln102", "ln2b"}


def main():
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    for name, volume, (lines, titles, copies) in VARIANTS:
        for ext in ("phr", "pin"):
            (out / f"{name}.{ext}").write_bytes(gzip.decompress((HERE / f"{volume}.{ext}.gz").read_bytes()))
        drop_seqids(out / name, lines, titles=titles, copies=copies)
    records, record = {}, None
    for row in (LOSAT / "tests/fasta/WSSV.faa").read_text().splitlines():
        if row.startswith(">"):
            record = row[1:].split()[0]
            records[record] = ""
        else:
            records[record] += row
    frag = {tag: records[acc][first - 1:last] for tag, (acc, first, last) in FRAGMENTS.items()}
    queries = {
        "ln101": [("q101", frag["101"])],
        "ln102": [("q102", frag["102"])],
        "ln101d": [("q101d", frag["101"] + LINKER + frag["101"])],
        "ln31": [("q31", frag["103"] + LINKER + frag["101"])],
        "lnall": [("q101", frag["101"]), ("q102", frag["102"]), ("q103", frag["103"])],
        "ln2b": [("q102", frag["102"]), ("q101", frag["101"])],
    }
    for name, recs in queries.items():
        (out / f"{name}.faa").write_text("".join(f">{n}\n{s}\n" for n, s in recs))
        if name in NUCLEOTIDE_QUERIES:
            back_translate(out / f"{name}.faa", out / f"{name}.fna")


main()
