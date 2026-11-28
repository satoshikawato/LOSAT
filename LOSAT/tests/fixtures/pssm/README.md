# PSSM files read by `-in_pssm`

`d3a_fg_glimmer_overflow.asn`: a text file of two lines (`>x` and `A 1 2 99999999999 1.0`) for the rows
`db.ps.inpssm.fg.glimmer_overflow` and `db.tbn.inpssm.fg.glimmer_overflow`, which test how NCBI's format guess reads such a
file; origin not recorded (added in commit 0e9a920f9, 2026-11-15, RG-D3-1).

`m4b2_pl1.asn.gz`: the PSSM of the rows `db.tmask.tn.pssm.*` (TBLASTN `-in_pssm` with `-db_soft_mask` and
`-db_hard_mask`, series NP, M4b PR 2; `../db_search/README.md`, section `dbmaskt/`), made with NCBI BLAST+ 2.17.0 from
`../db_search/dbmaskt/pn1.faa` (`BLAST_USAGE_REPORT=false`):

```
makeblastdb -in pn1.faa -dbtype prot -out pn1db
psiblast -query pn1.faa -db pn1db -num_iterations 1 -out_pssm pn1.asn -save_pssm_after_last_round -num_threads 1
gzip -9 -n pn1.asn
```

`ps6c_other_iupacaa.asn.gz`: NCBI BLAST+ 2.17.0's `-out_pssm` of `psiblast -in_msa tests/fixtures/msa/ps6b_fa_mol_other.aln -db tests/db_search/paja` (a query whose Seq-inst molecule is `other`, `iupacaa` letters), gzipped; the PSSM query of row `db.ps.inpssm.other.0` (line PS, PS6c).

`ps9_<name>.asnb` (26 files; `ps9_<name>.asnb.gz`, `gzip -9 -n`, when the BER passes 20 KB: `tk_paja tkxbz_paja tk_ratios tk_zr
ts1_tslice ts4_tslice tf1_hba idq`): the binary (BER) twins of the text PSSMs of the rows `db.ps.inpssm.*` and
`db.tbn.inpssm.*`, the queries of the rows `db.ps.inpssm.bin.*` and `db.tbn.inpssm.bin.*` (each is the text row of the same
name without `.bin`, with the twin as its `-in_pssm` file; line RG, D3c-3a). Each was written by `datatool` 2.24.0 of the
bioconda BLAST 2.17.0 package (`datatool -m all_pssm.asn -v x.asn -e x.asnb`; `all_pssm.asn` is NCBI's ASN.1 modules from
`general` to `blastdb` with `scoremat.asn` of the NCBI C++ Toolkit 2.17.0, kept at
`remaining-gaps-20261112/D/d3/c3a/twins/all_pssm.asn`), so every length is indefinite. The generator is
`remaining-gaps-20261112/D/d3/c3a/tools/make_twins.py`; it also checks for each file that datatool's text print of the binary
equals its text print of the source. The sources, all text files of this tree: `tk_paja tkxbz_paja rnd26 two rnd_scaled
rnd_notprot rnd_noscores rnd_idO rnd_localidO rnd_noungap rnd_byrow rnd_gi rnd_ctrl tk_ratios tk_zr` from `ps4_<name>.asn.gz`
(`two` is the twins of that file's two objects one after the other); `rnd` from
`../db_search/db.ps.asn.rnd.6.asn` and `tkrnd` from `../db_search/db.ps.asn.tkrnd.6.asn`; `rnd_ident_utf8` from
`ps5_rnd_ident_utf8.asn.gz`; `rnd_kungneg` from `ps5c_rnd_kungneg.asn.gz`; `bad_iupacaa` and `lower_iupacaa` from
`ps_bad_iupacaa.asn.gz` and `ps_lower_iupacaa.asn.gz`; `ts1_tslice ts4_tslice tf1_hba idq one` from `ps7_<name>.asn.gz`.
`ps5_descr_empty.asn.gz` has no twin: datatool refuses that text. datatool's text reader replaces a byte outside 0x20 to 0x7E
by `#`, so for three twins, `ps9_rnd_ctrl.asnb`, `ps9_rnd_ident_utf8.asnb` and `ps9_bad_iupacaa.asnb`, the generator put the
source's raw bytes back after the conversion (same lengths; the title and identifier bytes `C3 A9` and the iupacaa byte `E9`).
The NCBI answers of the rows were frozen with the official NCBI BLAST+ 2.17.0 Linux build, `-num_threads 1`.

`ps9_ber_<case>.asnb` (107 files): hand-made BER inputs for the rows `db.ber.ps.<case>` (psiblast, 107 rows, `-comp_based_stats 0
-evalue 1000`, database `tests/db_search/paja`) and `db.ber.tbn.<case>` (tblastn, 8 rows: `base def title_utf8tag title_ff
real_pinf cut_title after_broken len89`, database `tests/db_search/tslice`), in `tests/db_search_manifest/psiblast.ber.tsv` and
`tblastn.ber.tsv`. `ps9_ber_base.asnb` (2,907 bytes) is `ps9_rnd.asnb` without the Pssm member `intermediateData` (its 840
frequency ratios); most cases are that file with a few bytes changed. The cases with the frequency ratios (`def_full`,
`real_list_15x`, `real_list_pinf`) start from the whole `ps9_rnd.asnb` (it is also the second object of `after_twice`), and
`store_ff` and `store_ctrl` from `ps9_bad_iupacaa.asnb`. They were made by byte surgery with `remaining-gaps-20261112/D/d3/c3a/tools/gen_defects.py` (it writes
the files here and a table of case, size and change to `D/d3/c3a/logs/gen_defects.tsv`; `ber.py` in the same folder
reads and writes the BER tree): definite lengths and long length forms (81 to 89, 80 on a primitive value), tags, members
out of order or missing, INTEGER, ENUMERATED, Uint4 and BOOLEAN of other sizes and values, REAL texts and binary forms,
strings with bytes outside the visible range, lengths past the end, truncations and bytes after the object. Two Uint4 cases,
`uint_ff` and `uint5_hi`, end psiblast by a signal (no row is possible) and are not kept here. The NCBI answers were frozen with
the official NCBI BLAST+ 2.17.0 Linux build, `-num_threads 1`, and eleven rows were repeated with the bioconda 2.17.0 build with
the same bytes. The probes of `makeblastdb` on the same inputs are in `remaining-gaps-20261112/D/d3/c3a/probe/ber/SUMMARY.tsv`.
The six `block_*` cases (the query's ncbistdaa OCTET STRING with a length past the end, or as a BIT STRING) are the inputs of
the review of D3c-3a (`remaining-gaps-20261112/D/d3/c3a/review_probe/`), frozen with the official build under an address limit
of 8 GiB.

`d3b_<case>.asn.gz` (24 text files), `d3b_<case>.asn` (7 text files whose NCBI answer is an error that names the file,
kept uncompressed for that reason) and `d3b_b_<case>.asnb` (7 BER files): the inputs of the rows `db.ptr.ps.<case>`
(psiblast, `-comp_based_stats 0 -evalue 1000`, database `tests/db_search/paja`) and `db.ptr.tbn.<case>` (tblastn,
`-evalue 1000`, database `tests/db_search/tslice`) in `tests/db_search_manifest/psiblast.ptr.tsv` and `tblastn.ptr.tsv`:
NCBI's pointer forms where its reader reads a member or an element held by reference (an OPTIONAL member of a class type,
an element of a list of references, and in a skip every such member): `NULL`, an object index `@n` and a value of a named
type `:Type value`, in text, and `05 00`, the object index tag `5e` and the long tag `7f` with a type name, in BER.
`d3b_base.asn.gz` is `../db_search/db.ps.asn.rnd.6.asn` without the Pssm member `intermediateData` (as `ps9_ber_base.asnb`);
the other text files are that file with one member or element changed (`dl_*`: a pub descriptor whose author list is the
delayed member `names`), and the BER files are `ps9_ber_base.asnb` with one member or element changed. They were made by
`remaining-gaps-20261112/D/d3/c3b/tools/gen_rows.py`, which also writes the rows; the probes behind them (official NCBI
BLAST+ 2.17.0 Linux build, `-num_threads 1`, an address limit of 8 GiB) are in `remaining-gaps-20261112/D/d3/c3b/probe/`.
The `losat-rejection` rows are LOSAT's refusals of a value of another type than the declared one and of an element stated
as `NULL` outside a skip (`docs/manual/src/refusals-edge-cases.md`).

`ps15b2_tk_other_stdaa.asn.gz`: `ps4_tk_paja.asn.gz` with its query's Seq-inst molecule `aa` changed to `other` (the first `mol aa,`; ncbistdaa letters), gzipped; the PSSM query of rows `db.ps15b2.p.othr.crash.*`, `.open.*`, `.prow.*`, `.loop.*` and `.past.*` (line PS, PS15b-2 fix 1).

`d3c_r_<case>.asn.gz` (6 text files): the inputs of the rows `db.ptr.ps.o11.r_<case>*` (line RG item D3c-3c) in
`tests/db_search_manifest/psiblast.ptr.tsv`: `../db_search/db.ps.asn.rnd.6.asn` with an element stated as `NULL` in a
list of its query that NCBI's object manager does not read and that holds no type the XML archive writer of LOSAT
cannot write yet (Seq-hist replaces ids, alone and second, and assembly; Seq-ext seg and delta; BioSource subtype), so
that NCBI's `-outfmt 11` stops at the element in its three streams (ASN.1 text, XML, binary). They were made by
`remaining-gaps-20261112/D/d3/c3c/tools/gen_probeA2.py` (official NCBI BLAST+ 2.17.0 Linux build, `-num_threads 1`,
probe `c3c/probeA2/`).

`ps7b_allx_ratios.asn`: a 10-residue all-X protein query with frequency ratios only (`intermediateData.freqRatios`, no `finalData`), hand-written for `db.ps7b.inpssm.allx_ratios.B45.6` and `db.ps7b.tbn.allx_ratios.B45.6` (line PS, item PS7b) to test NCBI's ungapped-before-gapped score-block error order in `PsiBlastComputePssmScores` when `-matrix BLOSUM45` does not support 11/1.

`ps_sw_tk_scores.asn.gz`: `ps4_tk_paja.asn.gz` without the Pssm member `intermediateData` (scores only, no frequency
ratios), gzipped without a name or time by `../db_search/pssw/make_pssw.py`; the PSSM query of rows `db.sw.scores.*`
(line PS-SW, PSIBLAST `-use_sw_tback` with a PSSM of scores only).

`ps_sw_hd0.asn.gz`: the PSSM of `../db_search/t5/hd_q0.faa` from NCBI psiblast (`-subject ../db_search/t5/hd.faa
-num_iterations 1 -save_pssm_after_last_round`), gzipped without a name or time by `../db_search/pssw/make_pssw.py
--fix2-only`; the PSSM query of rows `db.sw.utf8.*` and `db.sw.rej.hdr2.*` (line PS-SW fix round 2, PSIBLAST
`-use_sw_tback` on the databases `pssw/pfu8` and `pssw/pfu2t`).
