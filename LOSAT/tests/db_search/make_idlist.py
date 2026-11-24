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
- `pz`: version 4 of pz.faa (gi 7001, then the general ids PIG|abc, xyz|77 and BL_ORD_ID|xyz: the tags
  of PIG and BL_ORD_ID are strings, which NCBI's simplification of the Seq-id throws on).
- `tq`: version 4 of tq.faa (`gnl|ti|abc|ref|NP_000009.1|`, then `ref|NP_000010.1|`).
- `cv0`, `cv1`: version 4 of cv0.faa (`lcl|cdd:5`, `ref|NP_000301.1|`) and of cv1.faa (`gnl|cdd|5`,
  `gnl|cdd|6`, `ref|NP_000302.1|`), the sequences of g.faa's OIDs 1-5, two volumes of one alias file.
- `gtx`: `gt` again, its TI data file 64 zero bytes longer than its index says (NCBI's TI ISAM is then
  not initialised: `x_InitSearch` fails).
- `gd`: `g4` again, its GI data file 64 zero bytes longer than its index says.
- `ln` (D1c-1b, review 2 R2 of D1c): version 5 of ln.faa (4 proteins, OIDs 1 and 3 with two lines) with
  taxids from ln.map, its headers rewritten so that the second line of OID 1, the only line of OID 2 and
  the first line of OID 3 have no Seq-id (makeblastdb writes none such; NCBI's lists erase such lines).
  `lnq.faa` holds the four sequences as queries, `lnq.fna` the same back-translated (BLASTX).
- `lx` (D1c-1b): version 5 of lx.faa (4 proteins, OID 3 with two lines) with the header taxids of lx.map
  and the taxid lookup files (`.pot`, `.ptf`, `.pto`) of a build with lx2.map: OIDs 1 and 2 are 9606 in
  their headers and 10090 in the lookup files, so a taxid list selects them and erases every line.
  `lxq.faa` / `lxq.fna` as for `ln`.
- `lt` (review 2 of D1c-1b, G1): version 5 of lt.faa in two volumes (`lt.00` OID 0, `lt.01` OID 1), taxids from
  lt.map: OID 1 holds `ref|NP_300002.1|` (10090) and `gnl|BL_ORD_ID|xyz` (9606), a string tag that NCBI's
  read of a volume after the first cannot make an OID (its exception), also where a list then erases that line.
  `ltq.faa` is OID 0's sequence.
- `lnn` (review 1 of D1c-1b, F1): version 5 of lnn.fna (3 nucleotides, taxids from lnn.map), the only line of
  OID 0 rewritten to have no Seq-id; OID 0 is `lnnq.fna` (the back-translation of `lnnq.faa`), OID 1 the same
  with every 37th base changed, OID 2 unrelated.
- `lr`, `ls` (D1c-1c, review 2 G2 of D1c-1b): version 5 of lr.faa and ls.faa in three volumes (one sequence each):
  OID 0 a filler, then the query `lrq.faa` (`ref|NP_700002.1| strong` in lr, `ref|NP_700001.1|` in ls) and a random
  1000-residue `gnl|BL_ORD_ID|xyz weak` (best E-value about 2.5), weak before strong in lr and after it in ls: under a
  negative seqid list NCBI's fetch of the weak subject cannot read its header (the string tag in a later volume), and
  the composition-based redo frees it unfetched once the strong hit fills the hit list. `lrn`, `lsn`: the same,
  back-translated (TBLASTN), `lrq.fna` the back-translated query (BLASTX). `lrqw.faa` (and `lrqw.fna`) holds that query
  and then the weak sequence `w`, whose own search fetches the weak subject: with `BATCH_SIZE=1` NCBI writes the first
  batch's report and then ends.
- `lm` (D1c-1c-b): `ln` with its first 30 letters made `Q` in every sequence and the SEG masks of `segmasker
  -parse_seqids` (algorithm 21) in the volume, the same lines rewritten without Seq-ids: blastdbcmd's masks of a
  sequence whose front line has none (NCBI's Seq-loc of the masks takes the null best id and throws).
- `pw`, `nw` (line DBID pfoff-7): version 4 of pw.faa and of nw.fna (its back-translation) with parsed Seq-ids: gi 7001
  (`ref|NP_7001.1|`), then the general ids PIG|abc (gi 7001's letters with 4 changes), ti|xyz, BL_ORD_ID|zz and xyz|77,
  whose lookups as Seq-align ids throw on the id itself where the tag is a string (PIG, BL_ORD_ID: `CObject_id::GetId()`)
  or not a number (ti); NCBI's prefetch of a batch fails there (`Error pre-fetching sequence data`) and its formats load
  the subjects one by one. `pwq.faa` (and `pwq.fna`, back-translated) holds a query without hits, then one query per
  sequence after the first (the PIG|abc query hits gi 7001 second).

The list files and alias files are written below (`LISTS`, `ALIASES`). The rows of the -db manifests
that use them are db.il.*; the blastdbcmd cases A272 on (D/fixtures/make_d1c_blastdbcmd.py).
"""

import shutil
import struct
import subprocess
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent / "idlist"
SOURCES = {"g.faa", "g.map", "n.fna", "n.map", "q.faa", "nq.fna", "tz.faa", "pz.faa", "tq.faa", "cv0.faa",
           "cv1.faa", "ln.faa", "ln.map", "lnq.faa", "lx.faa", "lx.map", "lx2.map", "lxq.faa", "lt.faa", "lt.map",
           "ltq.faa", "lnn.fna", "lnn.map", "lnnq.fna", "lnnq.faa", "lr.faa", "ls.faa", "lrq.faa", "lrqw.faa",
           "pw.faa", "nw.fna", "pwq.faa", "pwq.fna"}

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
    # loc7 (g4's OID 6): NCBI sorts an alias file's SEQIDLIST by its ids as written when the search has a
    # positive list, before it simplifies them, and its walk of sq_a.txt's simplified ids loses loc7.
    "u_loc7.txt": b"loc7\n",
    # A gi alone: the processing of a text list for a version 5 database drops it (no id is left).
    "u_gi.txt": b"gi|102\n",
    # tq's OID 0, whose line starts with the general id ti|abc.
    "u_np9.txt": b"NP_000009.1\n",
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
    # For cv0 and cv1: lcl|cdd:5 simplifies to cdd:5 for cv0, and again to gnl|cdd|5 for cv1 (each
    # volume's translation simplifies the list again); lcl|cdd:6 likewise.
    "cv_sq.txt": b"lcl|cdd:5\nlcl|cdd:6\nNP_000302.1\n",
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
    # lists kept), beside another (ORed), beside the bare volume (every OID and every line: the node
    # that names the volume reads no list below it, so NCBI neither reads g_gi's list nor attaches it).
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
    # pz: the general ids PIG|abc (OID 1) and BL_ORD_ID|xyz (OID 3), read where a range keeps them.
    "pz_g.pal": "TITLE pz g\nDBLIST pz\nGILIST tz7001.txt\n",
    "pz_r.pal": "TITLE pz r\nDBLIST pz\nFIRST_OID 2\nLAST_OID 3\n",
    "pz_rg.pal": "TITLE pz rg\nDBLIST pz_r pz_g\n",
    "pz_o.pal": "TITLE pz o\nDBLIST pz\nFIRST_OID 4\nLAST_OID 4\n",
    "pz_og.pal": "TITLE pz og\nDBLIST pz_o pz_g\n",
    # A SEQIDLIST over two version 4 volumes: each volume's translation simplifies the list again.
    "cv_sq.pal": "TITLE sq two volumes\nDBLIST cv0 cv1\nSEQIDLIST cv_sq.txt\n",
    # A TILIST over a TI ISAM that NCBI does not initialise, a GILIST over such a GI ISAM.
    "gtx_ti.pal": "TITLE ti a uninitialised\nDBLIST gtx\nTILIST ti_a.txt\n",
    "gd_gi.pal": "TITLE gi a uninitialised\nDBLIST gd\nGILIST gi_a.txt\n",
    # gu with every OID beside its GILIST under NSEQ and LENGTH (no GI ISAM).
    "gu_all.pal": "TITLE gu all\nDBLIST gu\n",
    "gu_alls.pal": "TITLE gu all or gi walk\nDBLIST gu_all gu_gil\n",
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


def ti_isam(base: Path, num_oids: int, pad: int = 0):
    """A numeric ISAM of 4-byte keys in one page: key 1000 + 10 * OID for every OID; `pad` zero bytes
    after the data file's rows, which its index does not count."""
    rows = [(1000 + 10 * oid, oid) for oid in range(num_oids)]
    data = b"".join(struct.pack(">II", k, o) for k, o in rows)
    header = struct.pack(">9i", 1, 0, len(data), len(rows), 1, 256, 0, 0, 0)
    index = header + struct.pack(">II", rows[0][0], rows[0][1]) + struct.pack(">Ii", 0xFFFFFFFF, 0)
    base.with_suffix(".pti").write_bytes(index)
    base.with_suffix(".ptd").write_bytes(data + bytes(pad))


def build_pf7(mk):
    """The databases of line DBID pfoff-7: pw, nw."""
    for name, source, kind in (("pw", "pw.faa", "prot"), ("nw", "nw.fna", "nucl")):
        run([mk, "-in", source, "-dbtype", kind, "-out", name, "-title", name, "-parse_seqids",
             "-blastdb_version", "4"])


def build_d1c_fix1(mk, common):
    """The databases of review 1 of D1c-1: pz, tq, cv0, cv1, gtx, gd."""
    for name in ("pz", "tq", "cv0", "cv1"):
        run([mk, "-in", f"{name}.faa", "-dbtype", "prot", "-out", name, "-title", name, "-parse_seqids",
             "-blastdb_version", "4"])
    run([mk, "-in", "g.faa", "-dbtype", "prot", "-out", "gtx", "-title", "gtx", *common, "g.map",
         "-blastdb_version", "4"])
    ti_isam(OUT / "gtx", 14, pad=64)
    run([mk, "-in", "g.faa", "-dbtype", "prot", "-out", "gd", "-title", "gd", *common, "g.map",
         "-blastdb_version", "4"])
    with open(OUT / "gd.pnd", "ab") as data:
        data.write(bytes(64))


def ber_parse(buf, pos):
    """One BER node of a header (`Blast-def-line-set`) at `pos`: (tag, bytes or child nodes), end."""
    tag, length = buf[pos], buf[pos + 1]
    pos += 2
    if length == 0x80:
        kids = []
        while buf[pos] or buf[pos + 1]:
            kid, pos = ber_parse(buf, pos)
            kids.append(kid)
        return (tag, kids), pos + 2
    if length & 0x80:
        width = length & 0x7F
        length = int.from_bytes(buf[pos:pos + width], "big")
        pos += width
    if tag & 0x20:
        end, kids = pos + length, []
        while pos < end:
            kid, pos = ber_parse(buf, pos)
            kids.append(kid)
        return (tag, kids), end
    return (tag, bytes(buf[pos:pos + length])), pos + length


def ber_encode(node):
    tag, value = node
    if isinstance(value, list):
        return bytes([tag, 0x80]) + b"".join(ber_encode(kid) for kid in value) + b"\0\0"
    if len(value) < 0x80:
        return bytes([tag, len(value)]) + value
    width = (len(value).bit_length() + 7) // 8
    return bytes([tag, 0x80 | width]) + len(value).to_bytes(width, "big") + value


def drop_seqids(base: Path, lines, kind="p"):
    """Empties the Seq-id list (member `seqid`, [1]) of the definition lines `lines` ({oid: {line}},
    0-based) of the version 5 volume `base` (.phr, or .nhr with `kind` "n"), with the header offsets of
    its index (.pin, .nin: the same layout up to them)."""
    pin = bytearray(base.with_suffix(f".{kind}in").read_bytes())
    phr = base.with_suffix(f".{kind}hr").read_bytes()
    at = 12  # version, sequence type, volume number
    for _ in range(3):  # title, LMDB name, date
        at += 4 + struct.unpack(">i", pin[at:at + 4])[0]
    num_oids = struct.unpack(">i", pin[at:at + 4])[0]
    at += 4 + 8 + 4  # OIDs, total length, longest
    offsets = [struct.unpack(">i", pin[at + 4 * i:at + 4 * i + 4])[0] for i in range(num_oids + 1)]
    out, starts = bytearray(), []
    for oid in range(num_oids):
        record = phr[offsets[oid]:offsets[oid + 1]]
        starts.append(len(out))
        if oid in lines:
            node, end = ber_parse(record, 0)
            assert end == len(record)
            for line in lines[oid]:
                members = node[1][line][1]
                for i, member in enumerate(members):
                    if member[0] == 0xA1:
                        members[i] = (0xA1, [(member[1][0][0], [])])
            record = ber_encode(node)
        out += record
    starts.append(len(out))
    for i, start in enumerate(starts):
        pin[at + 4 * i:at + 4 * i + 4] = struct.pack(">i", start)
    base.with_suffix(f".{kind}hr").write_bytes(bytes(out))
    base.with_suffix(f".{kind}in").write_bytes(bytes(pin))


def back_translate(faa: Path, fna: Path):
    """The proteins of `faa` as DNA, one fixed codon per letter (a BLASTX query of each)."""
    codon = dict(zip("ACDEFGHIKLMNPQRSTVWY", "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCG CAG CGT "
                                            "TCT ACC GTG TGG TAT".split()))
    out = []
    for line in faa.read_text().splitlines():
        out.append(line if line.startswith(">") else "".join(codon[letter] for letter in line))
    fna.write_text("\n".join(out) + "\n")


def build_d1c1b(mk, common):
    """The databases of D1c-1b: ln and lnn (lines without Seq-ids), lx (header taxids apart from the
    lookup's), lt (a string BL_ORD_ID tag in a later volume)."""
    run([mk, "-in", "ln.faa", "-dbtype", "prot", "-out", "ln", "-title", "ln", *common, "ln.map"])
    drop_seqids(OUT / "ln", {1: {1}, 2: {0}, 3: {0}})
    run([mk, "-in", "lx.faa", "-dbtype", "prot", "-out", "lx2", "-title", "lx", *common, "lx2.map"])
    run([mk, "-in", "lx.faa", "-dbtype", "prot", "-out", "lx", "-title", "lx", *common, "lx.map"])
    for ext in (".pot", ".ptf", ".pto"):
        shutil.move(OUT / f"lx2{ext}", OUT / f"lx{ext}")
    for path in OUT.glob("lx2.*"):
        path.unlink()
    back_translate(OUT / "lnq.faa", OUT / "lnq.fna")
    back_translate(OUT / "lxq.faa", OUT / "lxq.fna")
    run([mk, "-in", "lt.faa", "-dbtype", "prot", "-out", "lt", "-title", "lt", *common, "lt.map",
         "-max_file_sz", "1000B"])
    run([mk, "-in", "lnn.fna", "-dbtype", "nucl", "-out", "lnn", "-title", "lnn", *common, "lnn.map"])
    drop_seqids(OUT / "lnn", {0: {0}}, "n")


def build_d1c1c(mk):
    """The databases of D1c-1c: lr and ls (a later volume's string BL_ORD_ID tag on a weak subject), and their
    nucleotide forms lrn and lsn."""
    back_translate(OUT / "lrq.faa", OUT / "lrq.fna")
    back_translate(OUT / "lrqw.faa", OUT / "lrqw.fna")
    for name in ("lr", "ls"):
        run([mk, "-in", f"{name}.faa", "-dbtype", "prot", "-out", name, "-title", name, "-parse_seqids",
             "-max_file_sz", "1000B"])
        back_translate(OUT / f"{name}.faa", OUT / f"{name}.fna")
        run([mk, "-in", f"{name}.fna", "-dbtype", "nucl", "-out", f"{name}n", "-title", f"{name}n", "-parse_seqids",
             "-max_file_sz", "1000B"])
        (OUT / f"{name}.fna").unlink()


def build_d1c1cb(ncbi: Path):
    """The database of D1c-1c-b: lm (ln's lines without Seq-ids over a volume with SEG masks)."""
    lines = (OUT / "ln.faa").read_text().splitlines()
    (OUT / "lm.faa").write_text("\n".join(line if line.startswith(">") else "Q" * 30 + line[30:]
                                          for line in lines) + "\n")
    run([str(ncbi / "segmasker"), "-in", "lm.faa", "-infmt", "fasta", "-parse_seqids", "-outfmt",
         "maskinfo_asn1_bin", "-out", "lm.asnb"])
    run([str(ncbi / "makeblastdb"), "-in", "lm.faa", "-dbtype", "prot", "-out", "lm", "-title", "lm",
         "-parse_seqids", "-taxid_map", "ln.map", "-mask_data", "lm.asnb"])
    drop_seqids(OUT / "lm", {1: {1}, 2: {0}, 3: {0}})
    for name in ("lm.faa", "lm.asnb"):
        (OUT / name).unlink()


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
    build_d1c_fix1(mk, common)
    build_d1c1b(mk, common)
    build_d1c1c(mk)
    build_d1c1cb(ncbi)
    build_pf7(mk)
    for path in list(OUT.glob("*.pjs")) + list(OUT.glob("*.njs")):
        path.unlink()
    for name, data in LISTS.items():
        (OUT / name).write_bytes(data)
    for name, text in ALIASES.items():
        (OUT / name).write_text(text)


main()
