#!/usr/bin/env python3
"""Builds LOSAT/tests/db_search/idlist/: databases with the alias keys GILIST, TILIST and SEQIDLIST
(line RG, series D, D1c-1; remaining-gaps-20261112/D/SURVEY_D1c.md), with NCBI BLAST+ 2.17.0's
makeblastdb (comparison-only tooling, as the fixtures' NCBI outputs; LOSAT never runs NCBI).

Usage: make_idlist.py NCBI_BIN   (writes into the idlist directory next to this script; g.faa, g.map,
n.fna, n.map, q.faa and nq.fna there are the committed sources, every other file is made again).

Sources: `g.faa` holds 14 proteins (OIDs 0-13), some with two or three definition lines (OID 0:
gi 101 NP_000101.1 and gi 102 P10102.1; OID 2: gi 104 and gi 105; OID 5: gi 108 AAA00108.1 and
Q30109.1 without a gi; OID 6 the local id loc7; OID 7 gi 110 with the pdb id 1ABC_A; OID 10 three
lines, gi 113, NP_000114.1 and gi 115), taxids from `g.map`; `n.fna` 6 nucleotides with gis 201-208,
taxids from `n.map`; `q.faa` and `nq.fna` one query per database sequence.

Databases (the `.pjs`/`.njs` metadata, with the build's paths and dates, are dropped):
- `g`: version 5 of g.faa with parsed Seq-ids and taxids (a GI ISAM, no string ISAM, no TI ISAM).
- `g4`: version 4 of g.faa (GI and string ISAMs).
- `gm`: version 5 of g.faa in two volumes (`gm.00` OIDs 0-7, `gm.01` OIDs 8-13) under `gm.pal`.
- `gu`: version 5 of g.faa without parsed Seq-ids (no ISAM).
- `gt`: `g4` with a TI ISAM (`gt.pti`, `gt.ptd`) written here, of the keys 1000 + 10 * OID; makeblastdb
  writes no TI ISAM. NCBI's lists find no OID through it (survey D1c section 2.4), while its single
  lookups (`-entry 'gnl|ti|1020'`) find OID 2.
- `n`: version 5 of n.fna, `n4`: version 4 of it.
- `tz`: version 4 of tz.faa (gi 7001, then the general ids ti|abc and ti|77).

The list files and alias files are written below (`LISTS`, `ALIASES`). The rows of the -db manifests
that use them are db.il.*; the blastdbcmd cases A272 on (D/fixtures/make_d1c_blastdbcmd.py).
"""

import shutil
import struct
import subprocess
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent / "idlist"
SOURCES = {"g.faa", "g.map", "n.fna", "n.map", "q.faa", "nq.fna", "tz.faa"}

# An absolute path that is not a file: printed as it is, the same on every machine.
MISSING = "/losat-idlist-fixture-missing/none.txt"


def gi_bin(gis):
    return b"\xff\xff\xff\xff" + struct.pack(">I", len(gis)) + b"".join(struct.pack(">I", g) for g in gis)


def ti_bin4(tis):
    return b"\xff\xff\xff\xfd" + struct.pack(">I", len(tis)) + b"".join(struct.pack(">i", t) for t in tis)


def ti_bin8(tis):
    return b"\xff\xff\xff\xfc" + struct.pack(">I", len(tis)) + b"".join(struct.pack(">q", t) for t in tis)


LISTS = {
    # The searches' own -seqidlist files beside the alias files' lists (rows db.il.bp.gi.sl*, sq4.sl).
    "u_np101.txt": b"NP_000101.1\n",
    "u_p10102.txt": b"P10102.1\n",
    "u_np103.txt": b"np_000103.1\n",
    # GIs of g: OIDs 0 (its second line), 2 (its second line), 5 (its first line), 12.
    "gi_a.txt": b"102\n105\n108\n117\n",
    "tz7001.txt": b"7001\n",
    "gi_a.bin": gi_bin([102, 105, 108, 117]),
    # OIDs 0 (its first line), 1, 2 (its first line and its second).
    "gi_b.txt": b"101\n103\n104\n105\n",
    # The last number without a separator after it is lost (117), carriage returns and spaces separate.
    "gi_nonl.txt": b"102\n105\n108\n117",
    "gi_cr.txt": b"102\r105\r108 117\r\n",
    # 2^32 + 102 wraps to 102; a 0 is no number.
    "gi_wrap.txt": b"4294967398\n0\n105\n",
    "gi_unk.txt": b"9991\n9992\n",
    # Errors of the file.
    "gi_hash.txt": b"#c\n102\n",
    "gi_alpha.txt": b"102\nabc\n",
    "gi_nl1.txt": b"\n102\n",
    "gi_empty.txt": b"",
    "gi_badcnt.bin": gi_bin([102, 105])[:4] + struct.pack(">I", 3) + gi_bin([102, 105])[8:],
    "gi_short.bin": b"\xff\xff\xff\xff\x00\x00\x00",
    # The GIs of n: OIDs 0 (its second line), 2, 4 (its first line).
    "gi_n.txt": b"202\n204\n206\n",
    # Seq-ids for the string ISAM of g4 (OIDs 0 by its second line, 1, 6, 13), in the forms the
    # ISAM's simplification finds (FASTA forms, lower case, a pdb id with its chain, a Swiss-Prot
    # name), and a list of no id.
    "sq_a.txt": b"NP_000103.1\nP10102.1\nloc7\nNP_000118.1\n",
    "sq_c.txt": b"ref|NP_000103.1|\nsp|P20106.1|N6_HUMAN\nlcl|loc7\npdb|1ABC|A\n",
    "sq_d.txt": b"np_000103.1\nn9_mouse\n# a comment\n>NP_000118\n",
    "sq_empty.txt": b"",
    "sq_n.txt": b"NM_000203.1\nab000207.1\n",
    # TIs for gt: OIDs 2 and 5 by its single lookups.
    "ti_a.txt": b"1020\n1050\n",
    "ti_a4.bin": ti_bin4([1020, 1050]),
    "ti_a8.bin": ti_bin8([1020, 1050]),
    "ti_bad.bin": gi_bin([1020]),
    # A taxid list (version 4 reads none).
    "tx_a.txt": b"9606\n",
    # An OID list of OIDs 1, 2 and 3 of g.
    "m123.msk": (13).to_bytes(4, "big") + bytes([0b0111_0000, 0]),
}

ALIASES = {
    # GILIST over g (version 5).
    "g_gi.pal": "TITLE gi a\nDBLIST g\nGILIST gi_a.txt\n",
    "g_gil.pal": "TITLE gi a walk\nDBLIST g\nGILIST gi_a.txt\nNSEQ 4\nLENGTH 472\n",
    "g_gib.pal": "TITLE gi a binary\nDBLIST g\nGILIST gi_a.bin\n",
    "g_nonl.pal": "TITLE gi no newline\nDBLIST g\nGILIST gi_nonl.txt\n",
    "g_cr.pal": "TITLE gi separators\nDBLIST g\nGILIST gi_cr.txt\n",
    "g_wrap.pal": "TITLE gi wrap\nDBLIST g\nGILIST gi_wrap.txt\n",
    "g_unk.pal": "TITLE gi unknown\nDBLIST g\nGILIST gi_unk.txt\n",
    "g_hash.pal": "TITLE gi hash\nDBLIST g\nGILIST gi_hash.txt\n",
    "g_hashl.pal": "TITLE gi hash walk\nDBLIST g\nGILIST gi_hash.txt\nNSEQ 4\nLENGTH 472\n",
    "g_alpha.pal": "TITLE gi alpha\nDBLIST g\nGILIST gi_alpha.txt\n",
    "g_nl1.pal": "TITLE gi newline first\nDBLIST g\nGILIST gi_nl1.txt\n",
    "g_empty.pal": "TITLE gi empty\nDBLIST g\nGILIST gi_empty.txt\n",
    "g_badcnt.pal": "TITLE gi bad count\nDBLIST g\nGILIST gi_badcnt.bin\n",
    "g_short.pal": "TITLE gi short\nDBLIST g\nGILIST gi_short.bin\n",
    "g_miss.pal": f"TITLE gi missing\nDBLIST g\nGILIST {MISSING}\n",
    "g_missl.pal": f"TITLE gi missing walk\nDBLIST g\nGILIST {MISSING}\nNSEQ 4\nLENGTH 472\n",
    # Not beside the alias file: no other directory is tried (the path names the alias file's).
    "g_rel.pal": "TITLE gi relative\nDBLIST g\nGILIST sub/none.txt\n",
    "g_multi.pal": "TITLE gi multiple\nDBLIST g\nGILIST gi_a.txt gi_b.txt\n",
    "g_multil.pal": "TITLE gi multiple walk\nDBLIST g\nGILIST gi_a.txt\tgi_b.txt\nNSEQ 4\nLENGTH 472\n",
    # With the node's range and OID list (ANDed), under another GILIST (ANDed; the lines of both
    # lists kept), beside another (ORed), beside the bare volume (every OID, the lines of the list).
    "g_gr.pal": "TITLE gi b range\nDBLIST g\nGILIST gi_b.txt\nFIRST_OID 3\nLAST_OID 4\n",
    "g_go.pal": "TITLE gi b oidlist\nDBLIST g\nGILIST gi_b.txt\nOIDLIST m123.msk\n",
    "g_gb.pal": "TITLE gi b\nDBLIST g\nGILIST gi_b.txt\n",
    "g_nest.pal": "TITLE gi b over gi a\nDBLIST g_gi\nGILIST gi_b.txt\n",
    "g_sib.pal": "TITLE gi a or b\nDBLIST g_gi g_gb\n",
    "g_bare.pal": "TITLE gi a or all\nDBLIST g_gi g\n",
    # TILIST over a database without a TI ISAM, and over gt (no OID found).
    "g_ti.pal": "TITLE ti over g\nDBLIST g\nTILIST ti_a.txt\n",
    "g_til.pal": "TITLE ti over g walk\nDBLIST g\nTILIST ti_a.txt\nNSEQ 2\nLENGTH 200\n",
    "gt_ti.pal": "TITLE ti a\nDBLIST gt\nTILIST ti_a.txt\n",
    "gt_ti4.pal": "TITLE ti a 4\nDBLIST gt\nTILIST ti_a4.bin\n",
    "gt_ti8.pal": "TITLE ti a 8\nDBLIST gt\nTILIST ti_a8.bin\n",
    "gt_tib.pal": "TITLE ti bad\nDBLIST gt\nTILIST ti_bad.bin\n",
    "gt_sib.pal": "TITLE ti or all\nDBLIST gt_ti gt\n",
    # GILIST over a database without parsed Seq-ids (no GI ISAM).
    "gu_gi.pal": "TITLE gi over gu\nDBLIST gu\nGILIST gi_a.txt\n",
    "gu_gil.pal": "TITLE gi over gu walk\nDBLIST gu\nGILIST gi_a.txt\nNSEQ 4\nLENGTH 472\n",
    # SEQIDLIST and TAXIDLIST over g4 (version 4: the string ISAM; TAXIDLIST read nowhere).
    "g4_sq.pal": "TITLE sq a\nDBLIST g4\nSEQIDLIST sq_a.txt\n",
    "g4_sqc.pal": "TITLE sq fasta forms\nDBLIST g4\nSEQIDLIST sq_c.txt\n",
    "g4_sqd.pal": "TITLE sq lower case\nDBLIST g4\nSEQIDLIST sq_d.txt\n",
    "g4_sqe.pal": "TITLE sq empty\nDBLIST g4\nSEQIDLIST sq_empty.txt\n",
    "g4_sqm.pal": f"TITLE sq missing\nDBLIST g4\nSEQIDLIST {MISSING}\n",
    "g4_gs.pal": "TITLE gi b and sq a\nDBLIST g4\nGILIST gi_b.txt\nSEQIDLIST sq_a.txt\n",
    "g4_gsu.pal": "TITLE gi a or sq a\nDBLIST g4_gi g4_sq\n",
    "g4_gi.pal": "TITLE gi a v4\nDBLIST g4\nGILIST gi_a.txt\n",
    "g4_tx.pal": "TITLE tx v4\nDBLIST g4\nTAXIDLIST tx_a.txt\n",
    "g4_txm.pal": f"TITLE tx v4 missing\nDBLIST g4\nTAXIDLIST {MISSING}\n",
    # tz: the general id ti|abc (OID 1), whose tag NCBI's filtered header read converts to a number and
    # throws on, read where the range beside the GI list keeps OIDs 1 and 2.
    "tz_g.pal": "TITLE tz g\nDBLIST tz\nGILIST tz7001.txt\n",
    "tz_r.pal": "TITLE tz r\nDBLIST tz\nFIRST_OID 2\nLAST_OID 3\n",
    "tz_rg.pal": "TITLE tz rg\nDBLIST tz_r tz_g\n",
    # Two volumes: the list over both, over one volume only, and each volume's own list.
    "gm_gi.pal": "TITLE gi a two volumes\nDBLIST gm\nGILIST gi_a.txt\n",
    "gm_v0.pal": "TITLE gi a volume 0\nDBLIST gm.00\nGILIST gi_a.txt\n",
    "gm_v1.pal": "TITLE gi a volume 1\nDBLIST gm.01\nGILIST gi_a.txt\n",
    "gm_v01.pal": "TITLE gi a each volume\nDBLIST gm_v0 gm_v1\n",
    "gm_mix.pal": "TITLE gi a volume 0 and volume 1\nDBLIST gm_v0 gm.01\n",
    "gm_fr.pal": "TITLE gi a range\nDBLIST gm\nGILIST gi_a.txt\nFIRST_OID 6\nLAST_OID 13\n",
    # Nucleotides.
    "n_gi.nal": "TITLE gi n\nDBLIST n\nGILIST gi_n.txt\n",
    "n4_sq.nal": "TITLE sq n v4\nDBLIST n4\nSEQIDLIST sq_n.txt\n",
    "n4_gi.nal": "TITLE gi n v4\nDBLIST n4\nGILIST gi_n.txt\n",
}


def run(argv):
    subprocess.run(argv, cwd=OUT, check=True, stdout=subprocess.DEVNULL, stdin=subprocess.DEVNULL)


def ti_isam(base: Path, num_oids: int):
    """A numeric ISAM of 4-byte keys in one page: key 1000 + 10 * OID for every OID."""
    rows = [(1000 + 10 * oid, oid) for oid in range(num_oids)]
    data = b"".join(struct.pack(">II", k, o) for k, o in rows)
    header = struct.pack(">9i", 1, 0, len(data), len(rows), 1, 256, 0, 0, 0)
    index = header + struct.pack(">II", rows[0][0], rows[0][1]) + struct.pack(">Ii", 0xFFFFFFFF, 0)
    base.with_suffix(".pti").write_bytes(index)
    base.with_suffix(".ptd").write_bytes(data)


def main():
    ncbi = Path(sys.argv[1])
    for path in OUT.iterdir():
        if path.name not in SOURCES:
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
    mk = str(ncbi / "makeblastdb")
    common = ["-parse_seqids", "-taxid_map"]
    run([mk, "-in", "g.faa", "-dbtype", "prot", "-out", "g", "-title", "g", *common, "g.map"])
    run([mk, "-in", "g.faa", "-dbtype", "prot", "-out", "g4", "-title", "g4", *common, "g.map",
         "-blastdb_version", "4"])
    run([mk, "-in", "g.faa", "-dbtype", "prot", "-out", "gm", "-title", "gm", *common, "g.map",
         "-max_file_sz", "1000B"])
    run([mk, "-in", "g.faa", "-dbtype", "prot", "-out", "gu", "-title", "gu"])
    run([mk, "-in", "g.faa", "-dbtype", "prot", "-out", "gt", "-title", "gt", *common, "g.map",
         "-blastdb_version", "4"])
    ti_isam(OUT / "gt", 14)
    run([mk, "-in", "n.fna", "-dbtype", "nucl", "-out", "n", "-title", "n", *common, "n.map"])
    run([mk, "-in", "tz.faa", "-dbtype", "prot", "-out", "tz", "-title", "tz", "-parse_seqids",
         "-blastdb_version", "4"])
    run([mk, "-in", "n.fna", "-dbtype", "nucl", "-out", "n4", "-title", "n4", *common, "n.map",
         "-blastdb_version", "4"])
    for path in list(OUT.glob("*.pjs")) + list(OUT.glob("*.njs")):
        path.unlink()
    for name, data in LISTS.items():
        (OUT / name).write_bytes(data)
    for name, text in ALIASES.items():
        (OUT / name).write_text(text)


main()
