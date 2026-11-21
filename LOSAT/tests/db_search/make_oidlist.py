#!/usr/bin/env python3
"""Writes LOSAT/tests/db_search/oidlist/: alias files with the key OIDLIST over the databases of
taxfilter/ (tfp, 23 proteins; tfn, two volumes of 12 and 4 nucleotides), and their OID list files
(line RG, series D, D1b-1; remaining-gaps-20261112/D/SURVEY_D1b.md).

An OID list file is the big-endian index of its last OID, then one bit per OID, the first OID of a
byte in its most significant bit, set for an OID the list keeps (makeblastdb -oid_masks writes the
same form). The rows of the -db manifests that use these files are db.ol.*.
"""

from pathlib import Path

OUT = Path(__file__).resolve().parent / "oidlist"


def mask(last_oid: int, bits: bytes) -> bytes:
    return last_oid.to_bytes(4, "big") + bits


MASKS = {
    # tfp: OIDs 0, 2, 3, 6, 9, 10, 12, 15, 16, 17, 21, 22.
    "k.msk": mask(22, bytes([0b1011_0010, 0b0110_1001, 0b1100_0110])),
    # Each volume of tfn: its OIDs 0, 3, 5, 6, 8, 10 (tfn.01 has 4: 0 and 3).
    "n.msk": mask(11, bytes([0b1001_0110, 0b1010_0000])),
    # The last OID 3 and 4 bytes: the bits count up to the 32-bit word, every OID of tfp (probe H1).
    "h1.msk": mask(3, b"\xff\xff\xff\xff"),
    # The last OID 3 and 1 byte: the 2 bytes past the file read as zeros, OIDs 0-7 (probe H2).
    "h2.msk": mask(3, b"\xff"),
    # The last OID 0xFFFFFFFF: NCBI counts 0 OIDs, no byte is read (probe H3).
    "wrap.msk": mask(0xFFFF_FFFF, b"\xff\xff\xff"),
    # The header alone: NCBI's failed validation (probes G3, G8).
    "hdr.msk": mask(22, b""),
    # No byte: NCBI crashes (probe G5).
    "empty.msk": b"",
    # tfp: OID 2 alone, over the union of two ranges (p_cu, probe B1 of review A of D1b-1).
    "u2.msk": mask(22, bytes([0b0010_0000, 0, 0])),
}

# A path that is not a file: printed as it is (an absolute value is not joined to the alias file's
# directory), the same on every machine.
MISSING = "/losat-oidlist-fixture-missing/none.msk"

ALIASES = {
    "p_k.pal": "TITLE oidlist k\nDBLIST ../taxfilter/tfp\nOIDLIST k.msk\n",
    "p_knl.pal": "TITLE oidlist k walk\nDBLIST ../taxfilter/tfp\nOIDLIST k.msk\nNSEQ 99\nLENGTH 9999\n",
    "p_kn.pal": "TITLE oidlist k nseq\nDBLIST ../taxfilter/tfp\nOIDLIST k.msk\nNSEQ 99\n",
    "p_kr.pal": "TITLE oidlist k range\nDBLIST ../taxfilter/tfp\nOIDLIST k.msk\nFIRST_OID 3\nLAST_OID 15\n",
    "p_ku.pal": "TITLE oidlist k or range\nDBLIST p_k ../aliaskeys/p_range\n",
    "p_kpn.pal": "TITLE oidlist parent nseq\nDBLIST p_knl\nNSEQ 5\n",
    "p_h1.pal": "TITLE oidlist h1\nDBLIST ../taxfilter/tfp\nOIDLIST h1.msk\n",
    "p_h2.pal": "TITLE oidlist h2\nDBLIST ../taxfilter/tfp\nOIDLIST h2.msk\n",
    "p_wrap.pal": "TITLE oidlist wrap\nDBLIST ../taxfilter/tfp\nOIDLIST wrap.msk\n",
    "p_miss.pal": f"TITLE oidlist missing\nDBLIST ../taxfilter/tfp\nOIDLIST {MISSING}\n",
    "p_missl.pal": f"TITLE oidlist missing walk\nDBLIST ../taxfilter/tfp\nOIDLIST {MISSING}\nNSEQ 12\nLENGTH 1700\n",
    "p_dir.pal": "TITLE oidlist directory\nDBLIST ../taxfilter/tfp\nOIDLIST /\n",
    "p_hdr.pal": "TITLE oidlist header\nDBLIST ../taxfilter/tfp\nOIDLIST hdr.msk\n",
    "p_hdrl.pal": "TITLE oidlist header walk\nDBLIST ../taxfilter/tfp\nOIDLIST hdr.msk\nNSEQ 12\nLENGTH 1700\n",
    "p_empty.pal": "TITLE oidlist empty\nDBLIST ../taxfilter/tfp\nOIDLIST empty.msk\n",
    "n_n.nal": "TITLE oidlist n\nDBLIST ../taxfilter/tfn\nOIDLIST n.msk\n",
    "n_nnl.nal": "TITLE oidlist n walk\nDBLIST ../taxfilter/tfn\nOIDLIST n.msk\nNSEQ 30\nLENGTH 30000\n",
    "n_miss.nal": f"TITLE oidlist missing\nDBLIST ../taxfilter/tfn\nOIDLIST {MISSING}\n",
    "n_missl.nal": f"TITLE oidlist missing walk\nDBLIST ../taxfilter/tfn\nOIDLIST {MISSING}\nNSEQ 8\nLENGTH 5000\n",
    # The OIDs 1-3 and 7-9 of tfp, whose union is a bit array over OIDs 1-9: under an OID list
    # of OID 2 alone NCBI keeps OIDs 1 and 2 (review A of D1b-1, B1).
    "p_c1.pal": "TITLE oidlist c1\nDBLIST ../taxfilter/tfp\nFIRST_OID 2\nLAST_OID 4\n",
    "p_c2.pal": "TITLE oidlist c2\nDBLIST ../taxfilter/tfp\nFIRST_OID 8\nLAST_OID 10\n",
    "p_cu.pal": "TITLE oidlist union\nDBLIST p_c1 p_c2\nOIDLIST u2.msk\n",
    # A range across the two volumes of tfn with an OID list of each volume: OIDs 10 and 12.
    "n_nr.nal": "TITLE oidlist n range\nDBLIST ../taxfilter/tfn\nOIDLIST n.msk\nFIRST_OID 10\nLAST_OID 14\n",
}

# A text seqid list (two proteins of tfp) and a protein that no profile of rps_mini hits.
TEXT_FILES = {
    "ids.txt": "sp|P69905.2|\nsp|P01923.1|\n",
    "w36.faa": ">w36\n" + "W" * 36 + "\n",
}


def main() -> None:
    OUT.mkdir(exist_ok=True)
    (OUT / ".gitattributes").write_text(
        "# OID list and alias fixtures are compared byte for byte: no line-ending conversion.\n* -text\n")
    for name, data in MASKS.items():
        (OUT / name).write_bytes(data)
    for name, text in {**ALIASES, **TEXT_FILES}.items():
        (OUT / name).write_bytes(text.encode())


if __name__ == "__main__":
    main()
