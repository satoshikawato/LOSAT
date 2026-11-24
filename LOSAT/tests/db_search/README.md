# Inputs of the `-db` fixtures

Used by `LOSAT/tests/db_search_fixtures.py` and `LOSAT/tests/db_search_manifest/` (rows that run
`losat blastp -db ...`, `losat blastn -db ...`, `losat tblastn -db ...`, `losat blastx -db ...` or `losat tblastx -db ...` against frozen NCBI BLAST+ 2.17.0
`blastp -db ...`, `blastn -db ...`, `tblastn -db ...`, `blastx -db ...` or `tblastx -db ... -num_threads 1` bytes; the manifest's `program` column says which). NCBI is the
comparison oracle only; the tests read the committed files and run no NCBI program.

FASTA files:

- `residues.faa`: the records of `subjects()` of the DB-3 residue experiment whose sequence has no `-` (17
  records: a plain copy, three- and four-letter replacements of O U * B Z J X, a beta chain with U and O, and a
  ten-letter record).
- `residues_gap.faa`: all records of `subjects()` (20 records; `s_dash`, `s_dashrun` and `s_mix` hold `-`).
- `residues_q.faa`: the records of `queries()` (five queries).
- `o_stored.faa` / `o_stored_q.faa`: an alpha chain with a tract of 24 Q and twelve O (`s_o`, `q_o`), the
  same with X (`s_x`), without them (`s_plain`, `q_plain`), and a ten-letter record. NCBI's composition redo
  reads the stored O (26), which is not the X the query holds for its own O, so `q_o` against `s_o` is not
  "near identical", the subject is SEG-filtered and scores 266 bits; against `s_x` it scores 307 (and
  `-subject` gives 307 for both).
- `oseed.faa` / `oseed_q.faa` (rows `db.oseed.*`, line NP item DBO): `oseed_o`, 20 V, `AOA`, 10 V, `AOA`, 20 V, and
  `oseed_x`, the same with X; the queries `q_o`, `q_x` and `q_ox` are the same with I for V and O, X, or O then X in the
  two words. A query word `AXA` (the query's O is read as X) scores 7 against itself with BLOSUM62, under the threshold
  11, so NCBI's lookup table holds it under its own letters only: it finds the stored `AXA` of `oseed_x` and not the
  stored `AOA` of `oseed_o`, which it would meet only as a neighbour (O scores as X). At `-threshold 7` (`db.oseed.t7.6`)
  it is a neighbour of `AOA` and both subjects are found; with words of 6 (`db.oseed.w6.6`) the compressed alphabet
  has no letter for O or X. `oseed_bx_q.fna` is `q_x` as DNA (ATT for I, GCT for A, NNN for X), the BLASTX query of
  `db.oseed.bx.6`.
- `len0_l0015` (rows `db.ung.len0.*`, line NP item O4b, PR (c)): the subject of `tests/fasta/outfmt0/len0_l0015_s.faa`
  (one protein of 189 letters), searched by the 3 proteins of `len0_l0015_q.faa` with `-ungapped -window_size 0
  -word_size 7 -threshold 9 -matrix PAM30 -seg yes`: query `q1` shows an ungapped HSP of length 0 at subject offset 23.
- `oreal.faa` / `oreal_q.faa` (row `db.oreal.psi.6`, line NP item DBO): four proteins of `LOSAT/tests/fasta/AP027132.faa`
  (records `a_BDV02837.1`, `a_BDV02907.1`, `a_BDV02938.1`, `a_BDV02956.1`) with 2 % of their letters replaced by O, X or U, and a protein of
  `AP027131.faa` (`BDV02481.1_x`) with 3 % replaced by X (a fixed pseudo-random generator; the set was cut down from 1,116
  such proteins while NCBI's result stayed apart from that of a search that reads a stored O as X). In PSIBLAST's
  iterations 2 and 3 the PSSM scores a stored O with its own score for O (`BLAST_SCORE_MIN` at a position that has
  aligned sequences, the matrix's score of X at one that has none), not with its score for X.
- `opssm.faa` (row `db.opssm.inpssm.6`, line NP item DBO): `opssm_o`, 40 random letters, 12 O and 40 random letters (a
  fixed pseudo-random generator), and `opssm_x`, the same with X. The row's `-in_pssm` query
  `LOSAT/tests/fixtures/pssm/dbo_qx_o9.asn.gz` is NCBI's `-out_pssm` of `psiblast -query` `q_x` of `oseed_q.faa`
  `-db oseed -num_iterations 2 -comp_based_stats 0`, with the score for O set to 9 at every position (the other scores
  and the frequency ratios as NCBI wrote them): NCBI finds the 12 stored O at every query offset, which a search that
  reads them as X (score -1) does not.
- `o_stored_bx_q.fna`: `s_x` back-translated with its twelve X as `NNN` codons, the BLASTX query of `db.o_stored.bx.8`
  (BLASTX reads the stored O of `s_o` as stored, 154 identities; `-subject o_stored.faa` reads it as X, 166, row
  `pyl.blastx.ostored.8` of `outfmt0_manifest.tsv`).
- `bx_title_stop_q.fna`: two records with a title that ends with 24 nucleotide letters; the first has a
  hyphen and an invalid letter (`7`) in its first data line, the second starts with a line that is not plausible
  data, the BLASTX query of `db.bxtitle.stop.6` (NCBI writes the first record's title warning after its data
  warnings, then stops; the second record's title warning is never posted).
- `bx_title_later_q.fna`: the same first record; the second has a data line `EF-` (no nucleotide letter: a hyphen
  and two invalid letters, with their warnings), then a line that is not plausible data, which NCBI still checks
  because the record has no residues yet (`db.bxtitle.later.6`).
- `bx_gap_first_q.fna`: a gap line `>?5` and 60 nucleotides, then the record `r1`, the BLASTX query of
  `db.bxgap.first.6` (with `-parse_deflines` NCBI stops at line 1: `Input doesn't start with a defline or
  comment`) and of the control `db.bxgap.first.nopd.6` (without it the gap line starts a record without a
  defline). `bx_gap_comment_q.fna` is the same after a comment, a blank line and a `#` line
  (`db.bxgap.comment.6`: NCBI names line 4), and `bx_gap_bad_q.fna` starts with `>?abc` instead
  (`db.bxgap.bad.6`: NCBI stops at line 1 before it parses the gap).
- `bx_prescan_data_q.fna`, `bx_prescan_gap_q.fna`: a data line, or the gap line `>?5`, before any defline,
  then the record `r1` whose title is not US-ASCII (`café x`), the BLASTX queries of `db.bxprescan.*`: with
  `-parse_deflines` NCBI stops at line 1 (`Input doesn't start with a defline or comment`) in every format,
  before any title is written (and before it fails on that one). `bx_prescan_badgi_q.fna`: the record `r0` with that title, then `>gi|abc`, which
  NCBI cannot read as a Seq-id under `-parse_deflines` (`Could not construct seq-id from 'gi|abc'`): in one batch
  the error comes before `r0` is written, with `BATCH_SIZE=1` after; with `BATCH_SIZE=0` (`db.bxprescan.bs0.*`, without
  `-parse_deflines`) NCBI reads no query and stops with `Empty CBlastQueryVector`. `bx_prescan_badgi_big_q.fna`: the same with
  10,200 nucleotides in `r0`, more than BLASTX's batch of 10,002 letters (`db.bxprescan.badgi.big.5`,
  `BATCH_SIZE=1000000`: one batch).
- The rows `db.bxgapline.*` search `bslice` with `LOSAT/tests/fasta/outfmt0/bx_gapline_q.fna`, the query of the
  `-subject` rows `gapline.blastx.*` of `outfmt0_manifest.tsv`: gap lines of bad sizes (`>?abc` before any
  defline, `>?0`, `>?unk`, `>?99999999999`, `>?-5` as a record's last line), text after the size and unknown,
  conflicting or unclosed modifiers, which NCBI reads as `ParseGapLine` does (a gap of 1 after `Bad gap size`,
  warnings for the rest).
- `p3.faa`: records p01, p03 and p08 of `LOSAT/tests/seqdb/prot.faa`, in that order.
- `d1_q.faa`: D1, a hit that `-db` keeps and `-subject` drops. Two proteins of `LOSAT/tests/fasta/WSSV.faa`
  (YP_009220551.1 and YP_009220614.1: for each, NCBI `blastp -db paja -outfmt 6` prints 2 rows and
  `blastp -subject tests/fasta/PajaWSV.faa -outfmt 6` prints 1), then `nohit`, a 72-residue low-complexity
  glycine/serine sequence for which NCBI `-db paja` reports no hit. For the whole file: 4 rows with `-db`, 2
  rows with `-subject` (outfmt 6).
- `gapmix.faa` (DI-12): six proteins of the CDS translations of Salmonella Typhi CT18 (AL513382.1, from the
  drop-in survey's Roary data, the first six of at least 120 residues without `-`) with `-`: one inserted in
  the middle (`ins1_pA`), one replacing a residue (`rep1_pA`), a run of three (`ins3_pB`), one at the
  end (`end1_pC`), two at the start (`lead2_pC`), at the start, in the middle and at the end (`sme_pD`), one every
  15 residues (`every15_pE`), none (`ctrl_pF`), and records of `-` only (`allgap5`, `onegap`).
  `gapmix_q.faa`: the six proteins without `-`; `gapmix_q.fna`: the same back-translated with one codon per
  letter (the BLASTX queries).
- `nohit.faa`: `GGGGG`, a query without hits in `nophr` (the database of `LOSAT/tests/seqdb/prot.faa`).
- `empty.faa`: an empty query file.
- `split_q.faa`: one query of 20,600 letters (a fixed pseudo-random filler with records p10, p05, p06 and p01 of
  `LOSAT/tests/seqdb/prot.faa` at 1001, 9851, 15001 and 20401), which NCBI searches in two chunks (CHUNK_SIZE
  10000, overlap 100: [0, 10000) and [9900, 20600)); rows `db.split.*` against `LOSAT/tests/seqdb/db/prot_v4`.
- `vol2.faa`: `A` (record p10 of `LOSAT/tests/seqdb/prot.faa`) and `B`, the same with every fifth letter
  replaced by `W` (`A` where it is `W`); `vol2_q.faa`: `A` as the query.
- `q3b_chains.faa` (rows `db.ps.q3b.*` of psiblast `-sum_stats true`, line RG item Q3b): eight subjects made of parts
  of `LOSAT/tests/fasta/psiblast/tk.faa` with some letters changed and random letters between them (a fixed
  pseudo-random generator), whose HSPs with that query form linked sets: two parts far apart (`split2`), close
  (`split2_close`) and in reverse order (`reversed2`), three parts (`split3`), overlapping parts (`overlap`), two weak
  parts (`weak2`), one part (`single`) and one part three times (`tandem`).
- `q3b_weak.faa` (rows `db.cnt.l0q3w.*` through `dbcounts/l0q3w.pal`, review 4 R4-S1 and review 5 R5-S1 and R5-N2 of
  D1a-3): five subjects, each 20 random letters, `LOSAT/tests/fasta/psiblast/tk.faa` 1-110 (the strong part), a barrier
  of 30 letters chosen to score lowest against the query letters they meet on the two parts' diagonals, a weak part a few
  letters after the strong part in the query and 30 after it in the subject (so the two HSPs link with a small gap) and a
  barrier tail. The weak parts: `weak10` tk.faa 118-127, `weak11` tk.faa 115-125, `gap7` tk.faa 118-124, `W`, tk.faa
  125-131 (two stretches that each score below the gap trigger), `gap7p` the same as `gap7`, `long8k` tk.faa 118-127 with
  D127E followed by 7,810 random letters (8,000 in all). `weak10`, `weak11` and `gap7` choose the barrier with BLOSUM62:
  with the matrix no gapped extension crosses it with an X-drop of 15 bits, while with the PSSM
  `tests/fixtures/pssm/ps4_tk_paja.asn.gz` the strong HSPs of `weak11` and `gap7` cross it. `gap7p` (a barrier of 36)
  and `long8k` choose it with BLOSUM62 and the PSSM, so that neither extension crosses it. Generator:
  `remaining-gaps-20261112/D/fixtures/make_q3b_weak.py` (Python 3 only; a fixed pseudo-random generator for the random
  letters).
- `q3b_long.faa` (rows `db.ps.q3b.long.*` and `db.oldfsc.blastp.split.6`): one protein of 21,000 random letters (a fixed
  pseudo-random generator) with `LOSAT/tests/fasta/psiblast/tk.faa` planted at 500, 9,800, 14,000 and 20,500 (cut at the
  end), a query that NCBI splits into chunks of 10,000 that overlap by 100.
- The database `x2g_c0054` (rows `db.x2g.*`, line RG item X2g) is the subject
  `LOSAT/tests/fasta/outfmt0/x2g_c0054_s.faa` (70 random letters, case c0054 of the line's survey; its query is
  `x2g_c0054_q.fna` there): with `-ungapped -window_size 0` and a low `-threshold`, the search reports one-hit
  ungapped HSPs of length 0. `x2g_f0616`, `x2g_f3147` and `x2g_f1144` are the subjects `x2g_f0616_s.faa`,
  `x2g_f3147_s.faa` and `x2g_f1144_s.faa` there (cases of the line's format survey: an HSP of length 0 in frame
  -2, in frame +1 at query offset 0, at subject offset 0), `x2g_ms` the two subjects `x2g_ms.faa` (with the two
  queries `x2g_mq.fna`, the second's HSP of length 0) and `x2g_rk` the two subjects `x2g_rk.faa` (the one
  ranked second holds the HSP of length 0 of `x2g_f0616_q.fna`); `x2g_rk_ver.faa` and `x2g_rk_lcl.faa` there are its
  letters under the ids `gb|AAA00001.1|` and `gb|AAA00001.2|`, and `lcl|AAA00001` and `gb|AAA00001.1|` (written by hand
  for the rows `db.html.bx.x2g.rk.ver.subj.pd.f0` and `db.html.bx.x2g.rk.lcl.subj.pd.f0`: two versions of one accession,
  and a local id before a versioned one). `x2g_s1`, `x2g_s1z` and `x2g_s1d` are the
  subjects `x2g_s1_s.faa`, `x2g_s1z_s.faa` and `x2g_s1d_s.faa` there (pairs of the line's probe s1c, built for
  one-hit HSPs of length 0: with `x2g_s1_q.fna` one at subject offset 1 and, without the subject's first letter,
  one at offset 0; with `x2g_s1d_q.fna` two, the subject written twice), whose every HSP has length 0.
  `x2g_p0127` is the subject `x2g_p0127_s.faa` there and `x2g_ord` the three subjects `x2g_ord.faa` (rows
  `db.x2g.bx.of.*`, line RG T3c after X2g; titles `x2g p0127 subject` and `x2g ord subjects`, `-blastdb_version 4`):
  `x2g_p0127_q.fna` / `x2g_p0127_s.faa` is the pair p0127 of the X2g-1 probe P1 (`p_gen.py p 150 20261116 36 300
  12 150`), whose HSP of length 0 starts at a codon of an HSP of the query's frame +1, where NCBI's
  query-anchored merge frees a segment it still reads; `x2g_ord.faa` holds p0127's subject, a protein that
  matches the first 48 nucleotides of P1's p0055 query in frame +1, and p0055's subject, and `x2g_ord_q.fna` and
  `x2g_ord_q2.fna` the two queries in both orders. `x2g_r11_001` and `x2g_r11_019` (rows `db.x2g.bx.of.r11_*`,
  line RG T3c fix 3; titles `x2g r11_001 subjects` and `x2g r11_019 subjects`) are the five subjects
  `x2g_r11_001_s.faa` and `x2g_r11_019_s.faa` there, with their queries `x2g_r11_001_q.fna` (one) and
  `x2g_r11_019_q.fna` (two): the sets r11_001 and r11_019 of the seed-11 fuzz of the review 3 of T3c
  (`gen_sets.py 11 30`), whose HSPs of length 0 with `-strand minus` make segments of length 0 of a query row,
  where NCBI's `FillUnalignedRegions` frees one that it still reads.

- The database `x2d_a0m` (rows `db.x2d.*`, line RG item X2d; title `x2d a0m subjects`, `-blastdb_version 4`) is the
  two subjects `LOSAT/tests/fasta/outfmt0/x2d_a0m_s.fna`: `x2d_a0_r0`, whose frame +1 ends with the word W V K and
  two letters, and `x2d_motif`, a match of the PSSM `LOSAT/tests/fixtures/pssm/x2d_d50.asn.gz`'s rows 18 to 29 (rows from 0).
  That PSSM scores its `-` column 50 at rows 13 to 15, so with `-window_size 0 -threshold 1` the one-hit alignment
  of length -1 on that word gets its gapped start on frame +1's NULLB; `x2d_b1m_s.fna` there is the same with
  `x2d_b1_r0` (the word and one letter: the start one past the NULLB). Made by the line's `gen_x2d.py` (the
  PSSM's rows of -6 but row 9 at -61, W +60 at row 10, V -60, K +1, S +4 at row 14, the motif MHYFREDQILTA at
  +25; subjects of random codons without another TGG or any CCA, so that no other frame holds a W).
  `x2d_m3_s.fna` there (rows `db.x2.tblastn.pssmpast.*`, line RG J2A, from its review 1) is the reverse complement of
  `AA` and `x2d_b1_r0`: its frame -3, the last frame of NCBI's translation buffer, ends with the word and one letter.
  `x2d_d50r.asn.gz` is `x2d_d50.asn.gz` with its `-` column at 50 on rows 16 to 20 too (title `x2d dash 50 rows
  13-20`; `gzip -9n` of the text PSSM), so that the window's NULLBs past frame -3's NULLB sum above 0.

- The database `xs_w12w11` (rows `db.xs.*`, line RG item X-short; title `xs w12 w11 subjects`, `-blastdb_version 4`) is
  `LOSAT/tests/fasta/outfmt0/xs_w12w11_s.fna`: `w12` (`TGGTGGTGGTTT`, 12 bases) and `w11` (`TGGTGGTGGTT`, 11 bases),
  whose frame +1 is W W W, the last three letters of the query `xs_q.faa` (`MHYWWW`). NCBI BLAST+ translates a subject's
  frames for its traceback only from 12 bases (blast_hits.c:1162-1164); with `-window_size 0` the one-hit HSP on `w11`
  reaches the traceback, which reads the translation never made (NCBI crashes, LOSAT refuses). `xs_*_s.fna` there:
  the same words and their reverse complements (`m11`, `m12`, frame -1) at 9, 11 and 12 bases, and `h11`/`h12`
  (`ATGCATTATTT`, M H Y in frame +1) for the PSSM `x2d_d50.asn.gz`, whose motif rows 18-20 score M H Y +25 each. `xs_long32_s.fna` is
  `w11` followed by 21 A, searched in `-subject_loc` 1-11 and 1-12 (NCBI measures the range).
  `xsb_*` (rows `db.xs.tblastn.bat.*`, line RG item X-short-bat): two query batches (`BATCH_SIZE=1`) whose second
  is refused there. `xsb_qa.faa` is `qa` (31 residues); `xsb_s.fna` is `pa`, its back-translation (93 bases), then
  `w11`: `qa` hits only `pa`. `xsb_q.faa` is `qa` then `q` (`MHYWWW`), `xsb_k.faa` is `qk` (`SRPP`, which hits only
  `w12` of `xs_w12w11`) then `q`, and `xsb_ql.faa` is `qa` then `ql` (`qa` ten times, 310 residues, which
  `CHUNK_SIZE=100 OVERLAP_CHUNK_SIZE=80` split into chunks that NCBI would split again, and
  `CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200` into a first wrapped chunk that covers it and chunks without a query, at
  which NCBI's `CQuerySplitter::Split` stops first), and `xsb_qg.faa` is `qa` then `qg` (`qa` eight times, 248
  residues, which `CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200` split into 2 chunks with gaps between them, row
  `db.xs.tblastn.bat.gap.6`). `xsb_qk.faa` is `qk` alone, row `db.xs.tblastn.bat.qk.6`, whose report is what NCBI writes
  before it crashes in row `db.xs.tblastn.bat.db.6`, as row `db.xs.tblastn.bat.qa.6` (`.qa.7`) is for rows `.6` (`.7`).
  Rows `db.xs.tblastn.bat.nosub.7` and `.nosubneg.7` search
  `xsb_ql.faa` in `taxfilter/tfn` with `-negative_seqidlist seqidlist/sl_tfn_all.txt` (every id: no subject left),
  where NCBI still splits `ql`'s batch and stops; `xsb_l.faa` is `ql` alone, one batch, for rows
  `db.xs.tblastn.bat.min1*.6` on `aliaskeys/n_min1` (`MIN_SEQ_LENGTH 1`), where the setup of the first chunk that NCBI
  searches fails (`The average subject length is too short`) unless a chunk stop comes first, and
  `db.xs.tblastn.bat.e0.6` on `aliaskeys/e_nF` (volumes without an OID), whose setup fails to read the length of OID 0
  (`Database memory map file error`), and `db.xs.tblastn.bat.nosq.6` on `arch_tbn_snucl_nonsq` (no sequence file), whose
  search fails to read the first sequence of chunk 0 (the same error).

- `x2f_c0384_q.fna` / `x2f_c0384_s.faa`, `x2f_c0326_q.fna` / `x2f_c0326_s.faa`: two random pairs of the line RG
  X2f survey (cases c0384 and c0326), rows `db.x2f.*`: with `-window_size 0` and a low `-threshold`, BLASTX's
  one-hit HSPs of negative length whose gapped start's window passes the end of the subject `s1` (18 letters)
  or of the query batch. `x2f_in.faa` is `s1` (as `x2f_s1`) followed by a filler protein `x2f_f1`, and
  `x2f_last.faa` the same two in the other order, for the databases `x2f_in` and `x2f_last`.
  `x2f_page.faa` is a filler of 4,075 `X` (`x2f_f1`, no hits), then `x2f_s1`, for the database `x2f_page`, whose
  `.psq` file is exactly 4,096 bytes (1 + 4,076 + 19): a window past `s1`'s NULLB leaves the file's last
  memory page. `x2f_pagepad` (line RG J2A, row `db.x2f.bx.pagepad.6`, the expectation of `db.x2f.bx.page.6`) is
  `x2f_page` with 16 zero bytes after its `.psq` file (4,112 bytes; `.pin` and `.phr` copied unchanged), so that
  NCBI's window past `s1` reads zeros inside its mapping. `x2f_c0346_q.fna` / `x2f_c0346_s.faa` (case c0346 of the
  survey, rows `db.x2f.bx.crash.*`; the database `x2f_c0346` of its subject, title `x2f c0346 subject`): windows past
  the query batch's buffer, where NCBI BLAST+ 2.17.0 exits with SIGSEGV natively; the rows are frozen from NCBI
  under valgrind.

- The databases `x2_negwin` and `x2_negtail` (rows `db.x2.tblastn.qpast.db.6` and `db.x2.tblastn.frame.db.6`, line
  RG J2A; titles `x2 negwin subject (VNV in frame 1)` and `x2 negtail subject (VNV at the end of frame -3)`,
  `-dbtype nucl -blastdb_version 4`, only `.nin .nsq .nhr` kept) are `LOSAT/tests/fasta/outfmt0/x2_negwin_s.fna` and
  `x2_negtail_s.fna`: the `-db` twins of the rows `db.x2.tblastn.qpast.6` (the query `x2_negpast_q.faa`, 15 letters,
  whose one-hit windows pass the query batch's buffer) and `db.x2.tblastn.frame.6` (windows past the subject's
  translation buffer).
- `LOSAT/tests/fasta/outfmt0/j2a_bp_q.faa` / `j2a_bp_s.faa` (row `db.j2a.bp.past.6`, line RG J2A, from its review 1):
  the round-2 audit's BLASTP case R2A2-2 (`docs/evidence/losat_web_e2e/audit/round2/a2_blastp.md`), a query of 37
  letters and two subjects of 20 and 14, whose one-hit windows pass the query batch's buffer and a subject's,
  where NCBI BLAST+ 2.17.0 reads past its buffers (natively it crashed in the audit's runs).

Databases (all `-blastdb_version 4`, no `-parse_seqids`; only `.pin .psq .phr` are kept, the `.pjs` file
makeblastdb also writes was deleted). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/residues.faa -dbtype prot -out tests/db_search/residues \
        -title "residues O U * B Z J X" -blastdb_version 4
    makeblastdb -in tests/db_search/residues_gap.faa -dbtype prot -out tests/db_search/residues_gap \
        -title "residues O U * - B Z J X" -blastdb_version 4
    makeblastdb -in tests/db_search/gapmix.faa -dbtype prot -out tests/db_search/gapmix \
        -title "gap character - at the start, in the middle, at the end, alone" -blastdb_version 4
    makeblastdb -in tests/fasta/PajaWSV.faa -dbtype prot -out tests/db_search/paja \
        -title "PajaWSV proteins" -blastdb_version 4
    makeblastdb -in tests/db_search/o_stored.faa -dbtype prot -out tests/db_search/o_stored \
        -title "stored O against the query's X" -blastdb_version 4
    makeblastdb -in tests/db_search/oseed.faa -dbtype prot -out tests/db_search/oseed \
        -title "stored O and X in words under the threshold" -blastdb_version 4
    makeblastdb -in tests/fasta/outfmt0/len0_l0015_s.faa -dbtype prot -out tests/db_search/len0_l0015 \
        -title "len0 l0015 subject" -blastdb_version 4
    makeblastdb -in tests/db_search/oreal.faa -dbtype prot -out tests/db_search/oreal \
        -title "four proteins with stored O, X and U" -blastdb_version 4
    makeblastdb -in tests/db_search/opssm.faa -dbtype prot -out tests/db_search/opssm \
        -title "12 stored O, and X, between random letters" -blastdb_version 4
    makeblastdb -in tests/db_search/vol2.faa -dbtype prot -out tests/db_search/vol2_nophr \
        -title "two volumes, the second without its header file" -blastdb_version 4 -max_file_sz 128B
    rm tests/db_search/vol2_nophr.01.phr
    makeblastdb -in tests/db_search/bxchunk.faa -dbtype prot -out tests/db_search/bxchunk \
        -title "proteins planted in bxchunk_q.fna" -blastdb_version 4
    makeblastdb -in tests/db_search/q3b_chains.faa -dbtype prot -out tests/db_search/q3b_chains \
        -title "subjects whose HSPs with tests/fasta/psiblast/tk.faa form linked sets" -blastdb_version 4
    makeblastdb -in tests/db_search/q3b_weak.faa -dbtype prot -out tests/db_search/q3b_weak \
        -title "subjects whose weak HSPs with tests/fasta/psiblast/tk.faa link with small gaps" -blastdb_version 4
    makeblastdb -in tests/db_search/x2f_in.faa -dbtype prot -out tests/db_search/x2f_in \
        -title "x2f_s1, then a filler" -blastdb_version 4
    makeblastdb -in tests/db_search/x2f_last.faa -dbtype prot -out tests/db_search/x2f_last \
        -title "a filler, then x2f_s1" -blastdb_version 4
    makeblastdb -in tests/db_search/x2f_page.faa -dbtype prot -out tests/db_search/x2f_page \
        -title "a filler of 4,075 letters, then x2f_s1" -blastdb_version 4
    for x in pin phr; do cp tests/db_search/x2f_page.$x tests/db_search/x2f_pagepad.$x; done
    { cat tests/db_search/x2f_page.psq; head -c 16 /dev/zero; } > tests/db_search/x2f_pagepad.psq
    makeblastdb -in tests/db_search/x2f_c0346_s.faa -dbtype prot -out tests/db_search/x2f_c0346 \
        -title "x2f c0346 subject" -blastdb_version 4
    makeblastdb -in tests/fasta/outfmt0/x2g_c0054_s.faa -dbtype prot -out tests/db_search/x2g_c0054 \
        -title "x2g c0054 subject" -blastdb_version 4
    for n in f0616:x2g_f0616_s.faa f3147:x2g_f3147_s.faa f1144:x2g_f1144_s.faa ms:x2g_ms.faa rk:x2g_rk.faa; do
        makeblastdb -in tests/fasta/outfmt0/${n#*:} -dbtype prot -out tests/db_search/x2g_${n%%:*} \
            -title "x2g ${n%%:*} subjects" -blastdb_version 4; done
    for n in s1 s1z s1d; do
        makeblastdb -in tests/fasta/outfmt0/x2g_${n}_s.faa -dbtype prot -out tests/db_search/x2g_$n \
            -title "x2g $n subject" -blastdb_version 4; done
    makeblastdb -in tests/fasta/outfmt0/x2g_p0127_s.faa -dbtype prot -out tests/db_search/x2g_p0127 \
        -title "x2g p0127 subject" -blastdb_version 4
    makeblastdb -in tests/fasta/outfmt0/x2g_ord.faa -dbtype prot -out tests/db_search/x2g_ord \
        -title "x2g ord subjects" -blastdb_version 4
    for n in r11_001 r11_019; do
        makeblastdb -in tests/fasta/outfmt0/x2g_${n}_s.faa -dbtype prot -out tests/db_search/x2g_$n \
            -title "x2g $n subjects" -blastdb_version 4; done

`wm_s60k` (rows `db.wm.*`, line NP M2: `-window_masker_db`) is the subject `s60k.fa` of
`LOSAT/tests/blastn_regression_fixtures.py`, with only `.nin .nsq .nhr` kept:

    makeblastdb -in tests/fixtures/blastn_regression/inputs/s60k.fa -dbtype nucl -out tests/db_search/wm_s60k \
        -title "s60k Sakai 2000001-2060000 mutated" -blastdb_version 4

`rpt_d4` (rows `db.rpt.*` and the `rpt.*` cases of `LOSAT/tests/blastn_regression_fixtures.py`, line RG item D4:
blastn `-filtering_db`) is a database of repeats: pieces of the queries of `wm_query.fa` (some with changed letters
or indels, some reverse-complemented, two 4 and 7 letters apart) and 500 random letters, written by
`tests/db_search/make_rpt_d4.py` (a fixed pseudo-random generator, Python 3 only; it also writes the queries
`rpt_query_lc.fa` and `rpt_edge.fa` of `tests/fixtures/blastn_regression/inputs/`), with only `.nin .nsq .nhr` kept:

    python3 -I tests/db_search/make_rpt_d4.py
    makeblastdb -in tests/db_search/rpt_d4.fna -dbtype nucl -out tests/db_search/rpt_d4 \
        -title "rpt_d4 repeats" -blastdb_version 4

`rpt_d4_badcount.nal` (row `db.rpt.badcount.6`) is an alias file of `rpt_d4` whose `STATS_NSEQ` (`2e1`) NCBI cannot
convert to a number, written by hand.

`txrand.fna` (rows `db.cnt.l0txr.*` through `dbcounts/l0txr.nal`, review 6 R6-M1 of D1a-3): two subjects of 12,000 random
nucleotides each, with two to five pieces of `tests/db_search/taxfilter/q.fna` (60-250 letters, 20-35 % of them changed,
half reverse-complemented) written over them (records `g2_3_12000` and `g5_2_12000` of review 6's generator; a fixed
pseudo-random generator: `remaining-gaps-20261112/D/fixtures/make_txrand.py`, Python 3 only), with only
`.nin .nsq .nhr` kept:

    makeblastdb -in tests/db_search/txrand.fna -dbtype nucl -out tests/db_search/txrand \
        -title "random subjects with mutated pieces of tests/db_search/taxfilter/q.fna" -blastdb_version 4

`nophr.pin` and `nophr.psq` are copies of `LOSAT/tests/seqdb/db/prot_v4.pin` and `.psq`: the database without
its header file. `vol2_nophr` (one sequence per volume, `vol2_nophr.pal`) lacks the header file of its second
volume. NCBI opens a volume's header file when it first reads a header of it, after the search and for the
final hit lists only, so these searches succeed when no reported subject is in such a volume.

Nucleotide inputs (BLASTN rows, `db.blastn.*`), written by `make_nucl_inputs.py` (run from `LOSAT/`; slices
of `tests/fasta/LC738875.fasta` and `AP027202.fasta`, two related viruses, mutations from a seeded generator):

- `nslice.fna`: six records for the database `nslice`: `a_l6k` (LC738875 bases 1-6000), `a_mut_8pct` (the same
  with 8 % substitutions and six short indels), `b_p6k` (AP027202 bases 1-6000), `c_amb` (LC738875 20001-23000
  with the codes R Y K M S W B D H V N, a block of 60 N and further runs of ambiguity letters), `d_mut_15pct`
  (LC738875 40001-44000 with 15 % substitutions), `e_p2k` (AP027202 30001-32500).
- `nslice_q.fna`: five queries: `q_a` (1500 bases of `a_l6k`), `q_c_amb` (1400 bases across the ambiguity codes of
  `c_amb`, 4 % substitutions), `q_d_minus` (reverse complement of `d_mut_15pct` bases, 12 % substitutions),
  `q_e_p` (1600 bases of `e_p2k`, 20 % substitutions) and `q_nohit` (900 random bases).
- `nslice_nohit.fna`: `q_nohit` alone.
- `nslice_qa.fna`: `q_a` alone (one query, for `-query_loc` and `-subject_loc`; rows `db.tab.bn.qloc.6`, `tab.bn.loc.6`).
- `nnophr.nin` and `nnophr.nsq` are copies of `LOSAT/tests/seqdb/db/nucl_v4.nin` and `.nsq`: the nucleotide
  database without its header file (NCBI's `blastn` fails with exit 2 when it opens the header file).

The database (version 4, no `-parse_seqids`; only `.nin .nsq .nhr` are kept, the `.njs` file makeblastdb also
writes was deleted). Run from `LOSAT/`:

    python3 tests/db_search/make_nucl_inputs.py
    makeblastdb -in tests/db_search/nslice.fna -dbtype nucl -out tests/db_search/nslice \
        -title "nucleotide slices of two viruses" -blastdb_version 4
    rm tests/db_search/nslice.njs

The nucleotide rows also use `tests/seqdb/db/nucl_v4`, `nucl_v5` and `nucl_mv` with queries cut from
`tests/seqdb/nucl.fna`.

TBLASTN `-max_intron_length` inputs (rows `db.x3.tn.*` and `x3.tblastn.*`, line RG item X3), written by
`make_x3_inputs.py` (run from `LOSAT/`; the generator is seeded):

- `x3ex_q.faa`: the query BDT61866.1 of `tslice_q.faa` (731 residues).
- `x3ex.fna`: `ex_plus`, six pieces of 45 residues of the query (25 % of their residues substituted, back-translated)
  with 91, 212, 450, 1201 and 4002 random bases between them and flanks of 300, and `ex_minus`, its reverse
  complement: which pieces NCBI links depends on the largest subject gap, which `-max_intron_length` sets.
  Database `x3ex`:

      makeblastdb -in tests/db_search/x3ex.fna -dbtype nucl -out tests/db_search/x3ex \
          -title "x3 exon-like subject" -blastdb_version 4

- `tests/fasta/outfmt0/x3_edl_q.faa`: three transposases of `tests/fasta/EDL933.fna` (ISAs1 family 272653..273789,
  ISEc26 276198..277370, ISEc8 346224..347762), queries of the rows on a subject longer than 5,000,000 bases.

Translated-search inputs (TBLASTN rows, `db.tblastn.*`), written by `make_tblastn_inputs.py` (run from `LOSAT/`;
slices of `tests/fasta/MelaMJNV.fasta`, `LvMJNV.fasta` and `AP027133.fasta`; the protein queries are proteins of
`PemoMJNVA.faa`, `LvMJNV.faa` and `AP027131.faa` that NCBI `tblastn -subject` finds in those genomes; mutations
and ambiguity codes come from a seeded generator):

- `tslice.fna`: five records for the database `tslice`: `m_plus` (MelaMJNV 36001-44000: plus-strand genes, a run
  of 15 ambiguity letters and an N block of 30 inside a gene, single codes every 97 bases in another),
  `m_minus` (MelaMJNV 59001-62000: a minus-strand gene with a run of the codes R Y K M S W B D H V N and single
  codes every 151 bases), `v_slice` (LvMJNV 21001-39000: genes on both strands, ambiguity runs on both),
  `v_mut_8pct` (LvMJNV 28501-30700 with 8 % substitutions and whole-codon indels) and `v_far_22pct` (the same
  region with 22 % substitutions).
- `tslice_q.faa`: seven queries (PemoMJNVA BDT61857.1, BDT61859.1, BDT61866.1; LvMJNV BDT62090.1, BDT62091.1,
  BDT62092.1, BDT62096.1). With the default options NCBI `tblastn -db tslice -outfmt 6` prints 37 rows and
  `tblastn -subject tslice.fna` prints 36 (D1: a database answers 3 as the shortest subject length, a FASTA
  set its real shortest; the extra row is BDT62092.1 against `m_plus`, E value 0.11).
- `tslice_q_lc.faa`: LvMJNV BDT62095.1 with residues 41-100 and 191-230 in lower case.
- `tslice_nohit.faa`: a random protein of 40 residues with no hit in `tslice` or `tgc4`;
  `tslice_nohit_masked.faa`: a 72-residue glycine/serine sequence that SEG masks completely (NCBI warns that it
  cannot calculate Karlin-Altschul parameters).
- `tgc4.fna`: three records for the database `tgc4`: `h_slice` (AP027133 106001-121000, genetic code 4 with an N
  block, a run of ambiguity letters and single codes), `h_mut_10pct` (AP027133 118901-120300, 10 % substitutions
  and indels) and `tiny` (8 bases, shorter than a translated search can use). `tgc4_q.faa`: five AP027131 proteins.
- `nnophr` (above) is also the database without a header file for `db.tblastn.no_nhr.nohit.0`.

The databases (version 4, no `-parse_seqids`; only `.nin .nsq .nhr` are kept, the `.njs` files were deleted).
Run from `LOSAT/`:

    python3 tests/db_search/make_tblastn_inputs.py
    makeblastdb -in tests/db_search/tslice.fna -dbtype nucl -out tests/db_search/tslice \
        -title "nucleotide slices of two shrimp viruses" -blastdb_version 4
    makeblastdb -in tests/db_search/tgc4.fna -dbtype nucl -out tests/db_search/tgc4 \
        -title "slices of a genetic code 4 genome" -blastdb_version 4
    makeblastdb -in tests/fasta/outfmt0/x2_negwin_s.fna -dbtype nucl -out tests/db_search/x2_negwin \
        -title "x2 negwin subject (VNV in frame 1)" -blastdb_version 4
    makeblastdb -in tests/fasta/outfmt0/x2_negtail_s.fna -dbtype nucl -out tests/db_search/x2_negtail \
        -title "x2 negtail subject (VNV at the end of frame -3)" -blastdb_version 4
    rm tests/db_search/tslice.njs tests/db_search/tgc4.njs

TBLASTN tblastn-fast inputs (rows `tfast.*` of `outfmt0_manifest.tsv` and `db.tfast.*`, line DI4: AMRFinderPlus runs
`tblastn -subject GENOME -query AMRProt.fa -task tblastn-fast -word_size 5 -window_size 15 -threshold 100 ...`),
written by `make_tfast_amr_inputs.py` from the nucleotide sequences of ABRicate's `ncbi` database (NCBI's AMR
reference genes, `db/ncbi/sequences`; the accessions are in the files):

- `tfast_amr.fna`: two contigs of seeded random spacer DNA with ten AMR genes: forward (TEM-1, CTX-M-15, Sul1,
  Tet(A)), reverse-complemented (KPC-2, NDM-1, QnrS1), with about 5 % substituted bases (OXA-48 forward, Erm(B)
  reverse-complemented), the first 60 % of VanA, and a run of 40 N. The database `tfast_amr` is these contigs.
- `tfast_amr_q.faa`: thirteen proteins (translations): TEM-116 and CTX-M-14 (other alleles of the families in the
  contigs), the eight others in the contigs, and SHV-12 (related to TEM), OXA-1 and APH(3')-Ia (not in them).

Run from `LOSAT/`:

    python3 tests/db_search/make_tfast_amr_inputs.py ABRICATE/db/ncbi/sequences tests/db_search
    mv tests/db_search/tfast_amr_s.fna tests/db_search/tfast_amr.fna
    makeblastdb -in tests/db_search/tfast_amr.fna -dbtype nucl -out tests/db_search/tfast_amr \
        -title "synthetic contigs with AMR genes" -blastdb_version 4
    rm tests/db_search/tfast_amr.njs

Translated-query inputs (BLASTX rows, `db.blastx.*`), written by `make_blastx_inputs.py` (run from `LOSAT/`;
slices of `tests/fasta/MelaMJNV.fasta`, proteins of `tests/fasta/PemoMJNVA.faa`; the decoy proteins and the query
without hits come from seeded generators):

- `bslice.faa`: 35 proteins for the database `bslice` (14,391 residues, shortest 75): the 11 PemoMJNVA proteins that
  NCBI `blastx -subject` finds in the two query slices (BDT61857.1, BDT61858.1, BDT61861.1, BDT61866.1, BDT61876.1,
  BDT61959.1 to BDT61964.1) and 24 more proteins of at most 1,000 residues.
- `bslice_q.fna`: two queries: `m_fam` (MelaMJNV 26001-31000: minus-strand genes, a family of five related proteins
  with several HSPs each) and `m_mix` (MelaMJNV 56001-61600: a plus-strand and a minus-strand gene). With the
  default options NCBI `blastx -db bslice -outfmt 6` prints 36 rows and `blastx -subject bslice.faa` prints 35
  (D1: a database answers 3 as the shortest subject length, a FASTA set its real shortest; the extra row is
  `m_fam` against BDT61959.1, E value 0.21). Query codes 1 and 11 give the same output (they differ in start codons
  only); code 4 changes the scores.
- `bslice_q_lc.fna`: `m_fam` bases 1-3000 with 301-700 and 1501-1800 in lower case (`-lcase_masking`).
- `bslice_nohit.fna`: 300 random bases with no hit in `bslice` (candidate 25 of 40 from `random.Random(20261019)`;
  the others all have weak hits at the default E value). `bslice_alln.fna`: 300 N (NCBI warns that it cannot
  calculate Karlin-Altschul parameters).
- `bxchunk_q.fna`: one query of 25,001 random bases (`random.Random(5)`, as the review of #702 made it) with a
  back-translated protein of 150 residues planted 200 bases before each of 3000, 5997, 6000, 9705, 10001, 12000,
  15000, 19410 and 20000 (the chunk borders of several `CHUNK_SIZE` values), and short lower-case runs at those
  points; `bxchunk.faa`: the nine planted proteins with a quarter of their residues changed, for the database
  `bxchunk` (rows `db.bxchunk.*`: BLASTX `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE`). `bxchunk_empty.fna`: a record of a
  defline only (NCBI reports the query splitter's errors before the empty query's `Sequence contains no data`).
- `bxempty_*.fna` (rows `db.bxempty.*`, against `bxchunk`, with `BATCH_SIZE=200000`: LOSAT searches such a batch in
  units of 10002 letters): `bxempty_le.fna` is the query of `bxchunk_q.fna` and an empty record (a defline only);
  `bxempty_lel.fna` adds a second query (its first 12,000 bases) and a second empty record; `bxempty_none.fna` is an
  empty record, 5,001 records of 2 bases (no context to search) and an empty record; `bxempty_none1.fna` the 5,001
  records and one empty record.
- `bnophr.pin` and `bnophr.psq` are copies of `bslice.pin` and `.psq`: the database without its header file, used
  with the query without hits (with a hit NCBI fails with exit 2).

The database (version 4, no `-parse_seqids`; only `.pin .psq .phr` are kept, the `.pjs` file was deleted).
Run from `LOSAT/`:

    python3 tests/db_search/make_blastx_inputs.py
    makeblastdb -in tests/db_search/bslice.faa -dbtype prot -out tests/db_search/bslice \
        -title "protein slice of a shrimp virus proteome" -blastdb_version 4
    rm tests/db_search/bslice.pjs
    cp tests/db_search/bslice.pin tests/db_search/bnophr.pin; cp tests/db_search/bslice.psq tests/db_search/bnophr.psq

Query-coverage inputs of the ungapped searches (row `db.blastx.ungapped_qcov60_mts1.6`, and the TBLASTX
`-subject` fixtures `qcov.tblastx.*` of `LOSAT/tests/outfmt0_manifest.tsv`), written by `make_qcov_inputs.py`
(run from `LOSAT/`; random letters from `random.Random(20261019)`):

- `qcovp_q.fna`: the query, 600 random bases.
- `qcovp.fna`: `s_part` (the query's first 330 bases unchanged between random flanks: 55 % of the query and the
  best e-value) and `s_full` (the whole query with 25 % substitutions: a weaker e-value).
- `qcovp.faa`: the translations of the same two (`p_part`, `p_full`; `p_full` has its own substitutions), the
  records of the database `qcovp`.
- With `-max_target_seqs 1 -qcov_hsp_perc 60` NCBI reports `s_full` (`p_full` with `blastx -ungapped
  -comp_based_stats 0`): an ungapped search removes the HSPs below the coverage before the preliminary hit list,
  which has the final size, is filled (`blast_engine.c:1538`). Without the coverage it reports `s_part` alone.

The database (version 4, no `-parse_seqids`; only `.pin .psq .phr` are kept, the `.pjs` file was deleted). Run from
`LOSAT/`:

    python3 tests/db_search/make_qcov_inputs.py
    makeblastdb -in tests/db_search/qcovp.faa -dbtype prot -out tests/db_search/qcovp \
        -title "query coverage, ungapped" -blastdb_version 4
    rm tests/db_search/qcovp.pjs

`bxbest` (row `db.blastx.best_hsp.0`) is the database of `tests/fasta/outfmt0/bx_best_hsp_subject.faa`, the Ensembl
R64-1-1 peptide YOL075C, searched with `tests/fasta/outfmt0/bx_best_hsp_query.fna`, the Ensembl BDGP6.54 cDNA
FBtr0100310 (both copied unchanged from the release files; the same pair is the `-subject` fixture `best_hsp.blastx` of
`LOSAT/tests/outfmt0_manifest.tsv`). The subject's first HSP by E-value is a sum-statistics pair of 119 and 45.4 bits
(Expect(2) = 8e-38) and its third has 131 bits (5e-35), so the description table shows 131 and 5e-35. Version 4, no
`-parse_seqids`, `.pjs` deleted; run from `LOSAT/`:

    makeblastdb -in tests/fasta/outfmt0/bx_best_hsp_subject.faa -dbtype prot -out tests/db_search/bxbest \
        -title "BLASTX best HSP: yeast YOL075C" -blastdb_version 4
    rm tests/db_search/bxbest.pjs

Nucleotide-versus-nucleotide translated inputs (TBLASTX rows, `db.tblastx.*`), written by `make_tblastx_inputs.py`
(run from `LOSAT/`; slices of `tests/fasta/MelaMJNV.fasta`, `LvMJNV.fasta` and `PemoMJNVA.fasta`; mutations and
ambiguity codes come from a seeded generator). TBLASTX is slow, so the database holds 14 kb and the queries 5 kb:

- `xslice.fna`: six records for the database `xslice`: `m_plus` (MelaMJNV 37001-41500: a plus-strand gene with a
  run of 15 ambiguity letters, an N block of 30 and single codes every 97 bases, all inside the region the
  queries are homologous to), `m_minus` (MelaMJNV 59301-61700: a minus-strand gene with a run of the codes
  R Y K M S W B D H V N and single codes every 151 bases), `v_mut_8pct` (LvMJNV 28501-30700 with 8 % substitutions
  and whole-codon indels), `v_far_22pct` (the same region with 22 % substitutions), `v_slice` (LvMJNV 37001-39000:
  genes on both strands with an ambiguity run, an N block of 12 and single codes) and `tiny` (8 bases). The
  codes are stored by makeblastdb as random bases, which the seed stage reads (D2), and as the real codes, which
  the ungapped extension and the alignments read.
- `xslice_q.fna`: four queries: `q_p_plus` (PemoMJNVA 10601-12200, 1600 bases of the plus-strand gene BDT61858.1,
  about 71 % identical to the MelaMJNV gene), `q_p_minus` (PemoMJNVA 32001-33500, 1500 bases of the minus-strand
  gene BDT61866.1), `q_v_a` (LvMJNV 28501-29500) and `q_v_b` (LvMJNV 37801-38600). With the default options NCBI
  `tblastx -db xslice -outfmt 6` prints 188 rows (`-db` and `-subject` agree: the search is ungapped, so
  there is no D1).
- `xslice_nohit.fna`: a random query of 45 bases without hits (every random query of 60 bases or more has weak hits
  at E value 10 in this database); `xslice_nohit_masked.fna`: 300 bases of CAG repeats that SEG masks completely
  (NCBI warns that it cannot calculate Karlin-Altschul parameters).
- `nnophr` (above) is also the database without a header file for `db.tblastx.no_nhr.nohit.0`.
- No row has a second query batch: NCBI's TBLASTX batch is 10,002 bases and a search of that size is too slow
  for the fixtures; the batches are covered by the `-subject` fixtures.

The database (version 4, no `-parse_seqids`; only `.nin .nsq .nhr` are kept, the `.njs` file was deleted). Run from
`LOSAT/`:

    python3 tests/db_search/make_tblastx_inputs.py
    makeblastdb -in tests/db_search/xslice.fna -dbtype nucl -out tests/db_search/xslice \
        -title "nucleotide slices of two shrimp viruses for TBLASTX" -blastdb_version 4
    rm tests/db_search/xslice.njs

`xslice_flip.{nin,nhr,nsq}` is `xslice` with the stored 2-bit base changed at every ambiguous position (126 positions;
the ambiguity data, the index and the header file are copies). makeblastdb stores at each ambiguous
letter the base that `-subject` also draws (`CRandom(length)`), so for `xslice` a `-db` search that translates the
stored bases and one that re-resolves the letters as `-subject` does give the same bytes. `xslice_flip` tells them
apart (rows `db.tblastx.flipped.6` and `.0`; D2: the preliminary stage translates the stored bases, the
re-evaluation and the reports read the stored codes). It is hand-edited, not written by makeblastdb: copy the three
files of `xslice` and rewrite the copy's `.nsq` in place with `flip_stored_bases.py` (the rule is in its docstring).
Run from `LOSAT/`:

    cp tests/db_search/xslice.nin tests/db_search/xslice_flip.nin
    cp tests/db_search/xslice.nhr tests/db_search/xslice_flip.nhr
    cp tests/db_search/xslice.nsq tests/db_search/xslice_flip.nsq
    python3 tests/db_search/flip_stored_bases.py tests/db_search/xslice.fna tests/db_search/xslice.nin \
        tests/db_search/xslice_flip.nsq tests/db_search/xslice_flip.nsq

NCBI `tblastx -db xslice_flip -query tests/db_search/xslice_q.fna -outfmt 6 -num_threads 1` prints 184 rows (188
for `xslice`). `db.tblastx.dbsize.0` and `.6` use `-dbsize 1000000000` (112 rows in outfmt 6): at that value a
search that feeds `-dbsize` to the link cutoffs gives 111 rows, at 100000000 it gave the same bytes.

`residues_gap` and `gapmix` hold stored `-` (NCBIstdaa 0), which NCBI searches as a residue that scores
`BLAST_SCORE_MIN` against every letter and its reports show as the letter of a gap (rows `db.residues_gap.*`,
`db.gap.*`). The other databases of the rows are the ones of
`LOSAT/tests/seqdb/db` (see `LOSAT/tests/seqdb/README.md`). A rebuild changes the dates stored in the index
files, so the frozen outputs (`Posted date:` lines) must then be frozen again:
`python3 LOSAT/tests/db_search_fixtures.py freeze --ncbi-bin DIR`. A freeze runs NCBI BLAST+ with its usage
report off (`BLAST_USAGE_REPORT=false` in each NCBI run's environment, `db_search_fixtures.ncbi_env`): without it NCBI
BLAST+ 2.17.0 sends each run's parameters to NCBI when it ends (`blast_usage_report.cpp:195-220`). LOSAT's runs
(`check`) do not get the variable.

OF5a (XML2 and JSON, outfmt 15 and 16): `x16_taxid.faa` (the two `-subject` records of the XML2 probe sprot_tax, `st1` and
`st2`) and `x16_taxid_q.faa` (its query), and `x16_taxid`, built from them with a taxid in each definition line (only `.pin
.psq .phr` are kept):

    makeblastdb -in tests/db_search/x16_taxid.faa -dbtype prot -out tests/db_search/x16_taxid \
        -title "x16 taxid 9606" -blastdb_version 4 -taxid 9606

NCBI writes the taxid in the `Hit` description of outfmt 13-16 (`HitDescr_taxid`); the rows `db.sa.x16_taxid.15` and
`db.sa.x16_taxid.16` freeze that, with no `taxdb.bti` on the search path.

OF5a part 1b (taxonomy names in outfmt 15 and 16): `x16_taxid_99990` and `x16_taxid_7777` are `x16_taxid` built with
`-taxid 99990` and `-taxid 7777` (same command); `x16_taxid_q2.faa` is `x16_taxid_q.faa` with a third query (`qp3`, hits
both subjects) so that a subject is a `Hit` of two queries. The rows `db.sa.x16_taxid.<case>.{15,16}` set
`BLASTDB=tests/seqdb/taxdb/<variant>` in `extra_args` (a `BLASTDB=DIR` token is the environment variable of both programs,
not an argument): `good`, `bad_magic`, `count_larger`, `no_btd`, `empty_fields`, `offset_last_past_end` and
`tests/seqdb/taxdb_xml2/no_9606` (`make_taxdb.py write` of `taxdb.tsv` without the row of 9606). The rows `.nonascii.*`
(`high_bytes_utf8`, the scientific name of 9606 has bytes >= 0x80) are frozen from NCBI like the others (K4h): its XML2
writer stops at the byte (outfmt 16 is cut there with a `Frame type= eFrameNamed` post and exit status 0, outfmt 14
leaves the side file cut and exits 255), its JSON writer copies the bytes (outfmt 15 and 13).

K4h (titles with bytes >= 0x80 or control bytes in the XML family): `xb/` holds the inputs of the rows `db.xb.*`, made by
`make_xb_inputs.py` (its docstring names the byte classes: `u2`, `u2end`, `l1x`, `c01first`, `ent1`, `x81`, `ovl3`; a query
title against a clean `-subject` or `-db`, a `-subject` title with a hit, a `-subject` title that no hit reaches). `xb_p` and
`xb_n` (only `.phr .pin .psq` and `.nhr .nin .nsq` kept) are the databases of the `-db` rows, made with
`makeblastdb -in xb/xb_sp_clean.faa -dbtype prot -out xb_p -title "xb proteins" -blastdb_version 4` and
`makeblastdb -in xb/xb_sn_clean.fna -dbtype nucl -out xb_n -title "xb nucleotides" -blastdb_version 4`. Every row is frozen from
NCBI BLAST+ 2.17.0 (`freeze --only db.xb.`), the trace lines of its exceptions naming NCBI's own build path, which LOSAT prints.
K4j (BLASTX titles that are not UTF-8, and every other byte class of K4h, in every report): `xj/` holds the inputs of the
rows `db.xj.*` of `blastx.tsv` and `rpstblastn.tsv`, made by `make_xj_inputs.py` (its docstring names the files: one query
title per byte class, `-subject` titles of every class, long titles, ASCII titles that `GenerateDefline` changes,
`-parse_deflines` titles with the byte in the first word, empty subjects, rpstblastn queries); the `-db` rows search `xb_p`.
Every row is frozen from NCBI BLAST+ 2.17.0 (`freeze --only db.xj.`). The titles that the rows `db.t7.*` (`t7/`) cover, a
first word `Subject_N` with a byte NCBI cannot decode and single `x81` titles, are not repeated here.
`xb/e2/` (K4h review 3, case E2: rows `db.xb.e2.*`) holds the BLASTP and BLASTN queries of `e2probe.py prep` (remaining-gaps
K/k4h/review3/e2): q0 without hits, q1 titled `q1 a\xc3\xa9b z` that hits only `multidb/put` (`nut`), q2 that hits `multidb/pa` (`na`)
and `put` (`nut`); `qmp` has the ASCII title `q1 ab z`, `qm1`'s q2 hits only `pa` (`na`), `qm2`'s only `put` (`nut`).

OF5a part 2 (XML2 and JSON of BLASTN): `x16_ntitles.fna` (three records: an empty definition line, a title with an
ampersand, angle brackets and quotes, and `.`) and `x16_ntaxid.fna` with its query `x16_ntaxid_q.fna` (`qn1` hits both
subjects, `qn2` none, `qn3` is a slice of `qn1`), and the databases built from them (only `.nhr .nin .nsq` are kept):

    makeblastdb -in tests/db_search/x16_ntitles.fna -dbtype nucl -out tests/db_search/x16_ntitles \
        -title "x16 ntitles" -blastdb_version 4
    makeblastdb -in tests/db_search/x16_ntaxid.fna -dbtype nucl -out tests/db_search/x16_ntaxid \
        -title "x16 ntaxid 9606" -blastdb_version 4 -taxid 9606
    makeblastdb -in tests/db_search/x16_ntaxid.fna -dbtype nucl -out tests/db_search/x16_ntaxid_99990 \
        -title "x16 ntaxid 99990" -blastdb_version 4 -taxid 99990

The rows `db.sa.x16_ntitles.{15,16}` and `db.sa.x16_ntaxid[.<case>].{15,16}` are the BLASTN counterparts of the BLASTP rows of
part 1 and 1b (`BLASTDB=DIR` tokens as there); `db.sa.x16_ntaxid.lastoff.16` freezes the messages that NCBI's BLASTN
writes (`Error: [blastn] ...`).

- `bpqb_q.faa`: 30 records of `tests/fasta/MeenMJNV.faa` of 280-560 letters (11,940 letters, so NCBI searches
  them in two query batches of its default size; LOSAT's block search takes one block at the default
  `LOSAT_X_BPQBLOCK` of 50,000 letters and two with a value of at most 10,000, which
  `tests/blastp_query_blocks_fmt.rs` runs against these frozen bytes), the first 30 of that length in file order. Rows `db.sa.bpqb.{5,8,9,12}`: `blastp -db paja -evalue 0.1
  -max_target_seqs 2 -max_hsps 1` (9 rows of outfmt 6: some queries have hits, some do not).
- `bpqb5.faa` (rows `db.cnt.l0paja.*`, review 7 R7-S2 of D1a-3): `bpqb_q.faa` five times, each record's id suffixed
  `_c0` to `_c4` in turn (`for c in 0 1 2 3 4; do sed "s/^>\([^ ]*\)/>\1_c$c/" bpqb_q.faa; done > bpqb5.faa`): 150
  records, 59,700 letters, six query batches of NCBI's default size, which LOSAT searches as two blocks at the default
  `LOSAT_X_BPQBLOCK` of 50,000 letters (batches 1 to 5, and 6).

## -culling_limit rows (line NP, O3b)

Rows `db.blastp.culling_*` and `db.tblastn.culling_*` search `paja` with `bpqb_q.faa` (BLASTP) and `tslice` with
`tslice_q.faa` (TBLASTN) with `-culling_limit`. The split rows `db.{blastp,tblastn}.culling_split_*` use two databases
of the `culling_dom_s.{faa,fna}` subjects in `tests/fasta/outfmt0/` (45 protein records; 25 nucleotide records, the
`-subject` rows `culling.{blastp,tblastn}.*` of `outfmt0_manifest.tsv` search the same files) and the queries
`culling_dom_q.faa` / `culling_domq_q.faa`; their `CHUNK_SIZE=300 OVERLAP_CHUNK_SIZE=50` tokens are environment
variables of both programs (a query is searched in chunks of 300 letters with an overlap of 50; each chunk is
culled by itself). Version 4, no `-parse_seqids`; only `.pin .psq .phr` / `.nin .nsq .nhr` are kept. Run from `LOSAT/`:

    makeblastdb -in tests/fasta/outfmt0/culling_dom_s.faa -dbtype prot -out tests/db_search/culldom \
        -title "culling dom" -blastdb_version 4
    makeblastdb -in tests/fasta/outfmt0/culling_dom_s.fna -dbtype nucl -out tests/db_search/culldomn \
        -title "culling dom nucl" -blastdb_version 4
    rm tests/db_search/culldom.pjs tests/db_search/culldomn.njs
Rows `db.{blastx,tblastn,tblastx}.bh*` (line NP, O6a-3: `-best_hit_overhang` and `-best_hit_score_edge`) search `bh3p`
(BLASTX) and `bh3n` (TBLASTN, TBLASTX), the subjects `bh3_s.faa` / `bh3_s.fna` of `tests/fasta/outfmt0/` (40 records of
overlapping, nested and repeated pieces of the three `bh3_q.faa` proteins; `bh3_s.fna` back-translated, every other one
reverse-complemented; generator `o6a3_gen.py 7` of line O6), with the queries `bh3_q.fna` (BLASTX, TBLASTX) and `bh3_q.faa`
(TBLASTN). `db.tblastn.bh_inpssm.6` reads `tests/fixtures/pssm/o6a3_q0_bh3p.asn.gz`, NCBI psiblast 2.17.0's PSSM
(`-num_iterations 2 -out_pssm`) of the first `bh3_q.faa` query against the `bh3_s.faa` proteins. Version 4, no
`-parse_seqids`; only `.pin .psq .phr` / `.nin .nsq .nhr` are kept. Run from `LOSAT/`:

    makeblastdb -in tests/fasta/outfmt0/bh3_s.faa -dbtype prot -out tests/db_search/bh3p \
        -title "best hit O6a-3 prot" -blastdb_version 4
    makeblastdb -in tests/fasta/outfmt0/bh3_s.fna -dbtype nucl -out tests/db_search/bh3n \
        -title "best hit O6a-3 nucl" -blastdb_version 4
    rm tests/db_search/bh3p.pjs tests/db_search/bh3n.njs
- OF-BPQB2 (outfmt 15 and 16 in the same block search): rows `db.sa.bpqb.15` and `db.sa.bpqb.16`, the same search as
  `db.sa.bpqb.8` (65 KB at most). `paja_taxid` is `paja` built with a taxid (`make_paja_taxid.py MAKEBLASTDB OUTDIR`: `makeblastdb
  -in tests/fasta/PajaWSV.faa -dbtype prot -blastdb_version 4 -taxid 10090`, then the taxid of BDT63547.1 patched to 9606 in
  the header file, the same length); BDT63547.1 is one of the two subjects that only the second NCBI query batch of `bpqb_q.faa`
  hits, so the row `db.sa.bpqb.taxdb.16` (`BLASTDB=tests/seqdb/taxdb/good`) writes `Mus musculus` in the first block and
  `Homo sapiens` in the second. `tests/seqdb/taxdb_xml2/last_10090_past_end` is `make_taxdb.py all` of a two-entry TSV (9606, 10090),
  variant `offset_last_past_end` (the last entry, 10090, has its offset past the end of `taxdb.btd`); the row
  `db.sa.bpqb.taxdb.lastoff.15` freezes NCBI's messages across the two query batches (the volume's lookup once, the `Hit`'s each time).

OF6 (SAM, outfmt 17, BLASTN): the rows `db.sa.sam.*` (stems `mega`, `dc`, `tb`, `mv`, `short`, `ambrc`, `mts1`, `maxhsps1`,
`cull1`, `qloc`, `nohit`, `b100`; suffix `.17p` = `-outfmt 17`, `.17q` = `17 SQ`, `.17r` = `17 SR`, `.17b` = `17 SQ SR`, `.5` = the
outfmt 5 sibling that `S/scripts/check_sam_layout.py` predicts the document from) search `nslice` with `nslice_q.fna`
(`mega`, `dc`, `tb`, `short`, the option rows and `b100`: `BATCH_SIZE=100`, one query per batch), `nucl_mv` with
`tests/seqdb/nucl.fna` (`mv`) and `sam_ambrc` with `tests/fasta/outfmt0/sam_ambig_q.fna` (`ambrc`: ambiguity letters in
the stored data of a minus-strand hit). `sam_ambrc` is built from `tests/fasta/outfmt0/sam_ambig_rc_s.fna` (only
`.nhr .nin .nsq` are kept):

    makeblastdb -in tests/fasta/outfmt0/sam_ambig_rc_s.fna -dbtype nucl -out tests/db_search/sam_ambrc \
        -title "sam ambiguity letters on a minus strand" -blastdb_version 4

The frozen documents are NCBI's bytes, including the `CL:` tag of the `@PG` line (the path of the NCBI binary that
froze them); the check removes that tag from both sides (`db_search_fixtures.py strip_sam_cl`).

OF5b part 1 (XML2 and JSON of BLASTX): `make_x16_blastx_inputs.py` writes the inputs (run from `LOSAT/`):
`x16_bx_taxid_q.fna` (three nucleotide queries for `x16_taxid`: a minus frame, a plus frame, one without hits) and
`x16_bx_taxid_q2.fna` (two queries of 10,003 bases, each its own NCBI query batch of 10,002 bases, that both hit
`x16_taxid_99990`), `tests/fasta/outfmt0/x16_bx_q_unsearched.fna` (a 10,003-base slice with lower-case stretches, then a
300-N query and a 2-base query: a batch NCBI does not search) and, in `tests/fixtures/outfmt16`, `qnx_bad.fna`,
`qnx_firstbad.fna`, `qnx_late.fna`. BLASTX does not read `BATCH_SIZE`, hence the padding to 10,003 bases. The rows
`db.x16.blastx.*` (bslice) and `db.x16_bx_taxid.*` (taxid and `BLASTDB=` variants, as in the BLASTP rows) are frozen from NCBI
BLAST+ 2.17.0; `tests/fixtures/outfmt16/db/sprot_tax` is `sprot.faa` built with `makeblastdb -dbtype prot -blastdb_version 4
-taxid 9606` (only `.phr .pin .psq` kept).

OF5b part 1 fix (DECISIONS 32, BLASTX `-db` 15/16 with a taxonomy database): `tests/fixtures/outfmt16/db/sprot_mixtax.pal` lists
`sprot_tax10090` (`sprot.faa`, `makeblastdb -dbtype prot -blastdb_version 4 -taxid 10090`) and `other_tax9606` (`other.faa`, one protein
nothing hits, `-taxid 9606`; only `.phr .pin .psq` kept). `mixtax_blastx_{good,badname}_{15,16}` are NCBI BLAST+ 2.17.0
`blastx -query qnx_late.fna -db db/sprot_mixtax -num_threads 1` with `BLASTDB=tests/seqdb/taxdb/good` and `.../high_bytes_utf8` (the name of
9606 has bytes >= 0x80; 9606 is in the database and never reported, so NCBI writes the same bytes as with `good`).

- `x16_ttaxid_q.faa`, `x16_ttaxid_qwarn.faa`, `x16_ttaxid_qwarn2.faa` (OF5b part 2, TBLASTN outfmt 15 and 16): protein queries
  translated from `sn1` (frame +1) and `sn2` (frame -3) of `x16_ntaxid.fna`, and an unrelated one without hits; the `qwarn` files
  add a title that ends in 60 amino-acid letters (the reader's warning) and the `qwarn2` order puts the query without hits
  first. They search the existing nucleotide databases `x16_ntaxid` (taxid 9606) and `x16_ntaxid_99990`: rows
  `db.sa.x16_ttaxid.*` (`BLASTDB=DIR` and `BATCH_SIZE=N` tokens as in the BLASTP and BLASTN rows; the messages on standard
  error name `tblastn`).

- `x16_xtaxid_q.fna`, `x16_xtaxid_qwarn.fna`, `x16_xtaxid_qwarn2.fna` (OF5b part 3, TBLASTX outfmt 15 and 16): nucleotide queries made
  from `sn1` (its first 150 bases, `xq1`) and `sn2` (the reverse complement of its first 150 bases, `xq2`) of `x16_ntaxid.fna`, and an
  unrelated one without hits (`xq3`); the `qwarn` files add a title that ends in 28 nucleotide letters (the reader's warning) to `xq2`
  and the `qwarn2` order puts the query without hits first. With `-evalue 1e-20` each query hits its own subject only. They search the
  existing databases `x16_ntaxid` (taxid 9606) and `x16_ntaxid_99990`: rows `db.sa.x16_xtaxid.*` (`BLASTDB=DIR` and `BATCH_SIZE=N`
  tokens as in the BLASTP, BLASTN and TBLASTN rows, the other tokens are options; the messages on standard error name `tblastx`).

## Databases with parsed Seq-ids (line DBID, ID2a)

`ids/`: rows `db.id.*`. `ids/make_ids.py` (its docstring has the usage and the records) writes
`ids/prot_v5` (version 5) and `ids/prot_v4` (version 4) with `makeblastdb -parse_seqids` from the records of
`LOSAT/tests/seqid/fasta` and `ids/gen.faa` (titles that NCBI makes from the Seq-ids), with every file
makeblastdb writes: NCBI's formatters find a reported sequence again through the ISAM or LMDB files, so the
frozen outputs need them (LOSAT reads the index, sequence and header files). Queries: `ids/ids_q.faa` (three
proteins: one hits every sequence, one hits a few weakly, one has no hit), `ids/ids_q.fna` (the same
back-translated, for BLASTX) and `ids/ids_qw.faa` (thirty W, which hits the records of
`LOSAT/tests/seqid/db/variants` that no NCBI formatter finds again). `ids/refuse_v5` (version 5, from
`ids/refuse.faa` and `ids/taxmap_refuse.txt`) holds one record per case of the rows `db.id.*.taxdb_nf.*` and
of LOSAT's refusals, each hit by its own query `ids/refuse_<case>.faa` (`.fna` for BLASTX) alone with
`-evalue 1e-20`: NCBI's outfmt 0 does not end on the records `long` and `punct`, so no other query may
report them.

Nucleotides (line DBID, ID2b; rows `db.id.blastn.*`, `db.id.tblastn.*`, `db.id.tblastx.*`): `make_ids.py` also
writes `ids/nucl_v5` and `ids/nucl_v4` from `ids/nucl.fna` and `ids/taxmap_nucl.txt` (the Seq-id forms of
`LOSAT/tests/seqid/fasta/nucl_ids.fa`, every record with its sequence, and records whose titles NCBI makes
from their Seq-ids), `ids/nrefuse_v5` from `ids/nrefuse.fna` (`lcl|777` before `gi|777`, which NCBI's
formatters do not find again) and the two volumes `ids/njoined_a` and `ids/njoined_b` of the alias
`ids/njoined.nal` (`gb|AF000040|` and `gb|AF000040.1|`, whose hits NCBI's XML and outfmt 0 join). Queries:
`LOSAT/tests/seqid/fasta/qnucl.fa` (BLASTN, TBLASTX) and `qprot.fa` (TBLASTN); the second query of each has
no hit. TBLASTN and TBLASTX run with `-max_hsps 1` to keep the frozen files small. `ids/nambig_v4` (from
`ids/nambig.fna`: a record with ambiguity data) has no row; `tests/blastn_report_stream.rs` searches it.
Rows `db.id.blastn.stream.q25.*` (line DBID, blastn's report stream of these databases): the 25 records of
`ids/nucl.fna` as queries in batches of `BATCH_SIZE=1500`, each batch with hits.

Seq-align ids by which NCBI's formatters find another sequence (line DBID N1; rows `db.id.*.found.*`):
`make_ids.py` writes, from the FASTA of `ids/make_found_fasta.py`, `ids/pfound_v5` and `ids/pfound_v4` from
`ids/pfound.faa` (five pairs of unrelated proteins, each hit by its query `ids/pfound_q<k>.faa`;
`ids/pfound_q1.fna` and `ids/pfound_q4.fna` for BLASTX), `ids/nfound_v5` and `ids/nfound_v4` from
`ids/nfound.fna` (`gb|AF000050.1|`, then `gb|AF000050|`; queries `qnucl.fa` and `qprot.fa`), and the aliases
`ids/pdup_ab.pal` (two databases with the same record `gb|EEE00005.1|`) and `ids/pdup_ac.pal` (the second
copy with another title). NCBI looks the Seq-align id of a reported sequence up again in the database and
takes the first sequence the lookup gives: by `gb|AAA00001|` it finds `gb|AAA00001.1|` before it, by
`emb|BBB00002.1|` the `gb|BBB00002.1|` before it, by `lcl|CCC00003` the `gb|CCC00003.1|` before it (whose ids
do not compare equal, so NCBI finds no sequence); `gb|DDD00004|` before `gb|DDD00004.1|` and `pdb|2XYZ|` before
`pdb|2XYZ|AA` each find their own. LOSAT refuses such a sequence in the formats that read its Bioseq, unless the
sequence found holds the same definition lines and letters (`pdup_ab`); `-negative_seqidlist
seqidlist/sl_pfound_neg.txt` leaves `gb|AAA00001.1|` out of the search, and NCBI's lookup with it.

The same record with other taxids (line DBID N1b; rows `db.id.*.found.taxid.*`): `make_ids.py --copytaxid-only`
writes `ids/pdup_t` (`pdup.faa` with taxid 9606), `ids/pdup_u` (without one), `ids/pdup_n` and `ids/pdup_o` (both
with taxid 99990, the last entry of the taxonomy databases of `../seqdb/taxdb`), `ids/ndup_a` and `ids/ndup_t`
(`ndup.fna`, the back-translation of `pdup.faa`, without a taxid and with 9606), `ids/pdupe_a` and `ids/pdupe_t`
(`pdupe.faa`, without a title, without a taxid and with 9606), and the aliases `ids/pdup_at.pal`, `ids/pdup_tu.pal`,
`ids/pdup_no.pal`, `ids/pdupe_at.pal` and `ids/ndup_at.nal` (queries `pdup_q.faa` and `ndup_q.fna`). NCBI shows for
both copies the one its lookup gives first, with its taxids: the same rows where no report writes taxids, which
LOSAT writes; LOSAT refuses the sequence in the formats of taxids and their names, for `pdupe_at` with the
taxonomy database (a title made from the taxids), and with a taxonomy database on `BLASTDB` whose texts depend on
the copy NCBI fetches (`bad_magic` for `pdup_at`, `offset_last_past_end` for the copies of 99990 in `pdup_no`).

Copies whose title NCBI can make from the taxids (line DBID N1b L1; rows `db.id.title.*`): `make_ids.py
--titlesrc-only` writes `ids/ndupe_a` and `ids/ndupe_t` (`ndupe.fna`, the letters of `ndup.fna` without a title,
without a taxid and with 9606, through `ids/ndupe_at.nal`), `ids/psubj_a`, `ids/psubj_t` (`psubj.faa`, lcl|Subject_1
without a title) and `ids/psubjb_a`, `ids/psubjb_t` (`psubjb.faa`, lcl|Subject_2 with a title ending in `]`), through
`ids/psubj_at.pal` and `ids/psubjb_at.pal`, and `ids/nsubj_t` (`nsubj.fna`, lcl|Subject_1 without a title, taxid
9606, one database). NCBI generates the title of the copy it fetches, from its taxids where the taxonomy database names them,
only in outfmt 0 and for the tabular word that replaces a `lcl|Subject_` id: LOSAT writes the other formats
(`pdupe_at` and `ndupe_at` with `tests/seqdb/taxdb/good`, `psubj_at`, whose word is `unnamed` either way) and refuses
those two where the word is printed (`psubjb_at`, whose word is taken after NCBI rewrites the end of the title). Of
`nsubj_t` alone with the taxonomy database LOSAT writes NCBI's word `Homo`, the organism's (line RG item T5; refused
before).

An empty title without parsed Seq-ids (line DBID H3; rows `db.id.bare.*`): `make_ids.py --emptytitle-only` writes
`ids/nbare_t` and `ids/pbare_t` (`nbare.fna` and `pbare.faa`, the letters of `ndup.fna` and `pdup.faa` under `>` alone,
without `-parse_seqids`, taxid 9606). The Bioseq that NCBI's formatters fetch has a source for the taxid, so outfmt 0
makes the title from it where the taxonomy database names the taxid (`tests/seqdb/taxdb/good`: `unnamed protein
product [Homo sapiens]`, `Homo sapiens`), as LOSAT makes it (line RG item T5; refused before), and without a
taxonomy database or with `no_btd` the title of a Bioseq without a name; the other formats show the stored title (XML)
or the first word of the made one (the tabular id of a protein is `unnamed` either way).

Titles made from the taxonomy database's names (line RG item T5; rows `db.t5.*`): `t5/make_t5.py` (run from `LOSAT/`
with the directory of NCBI's makeblastdb) writes the taxonomy database `../seqdb/taxdb/t5` (`t5.tsv`: `uncultured
bacterium`, `Bacillus sp.`, `  Spaced   name  `, `A &amp; B &lt;x&gt; &#65;`, `TPA: tpa organism`, `MAG uncultured
thing`, a name with the byte 0x81, which NCBI cannot decode, and an empty scientific name) and in `t5/` one version 4
database of one sequence without a title per taxid (`n_<taxid>`, `p_<taxid>`, the letters of `ids/ndup.fna` and
`ids/pdup.faa`; 77777 is not in the taxonomy database), without `-parse_seqids`, and `np_<taxid>`, `pp_<taxid>` with
it (`gb|AB0000nn.1|`), with the alias files `nall`, `pall` (every taxid), `nclean`, `pclean` (no name NCBI cannot
decode), `ntab` (the named taxids: an empty nucleotide title has no subject id) and `npall`, `ppall`. NCBI makes each
title from the name (`CDeflineGenerator`: a nucleotide's name without white space at its ends and its first letter
upper-cased, a protein's `unnamed protein product [name]`; the heading drops a `TPA:` or `MAG` prefix, the description
keeps it). The same databases serve `blastdbcmd -long_seqids` (`../blastdbcmd_discovery.rs`, cases `T501`-`T507`) and
`dustmasker` and `segmasker` (`-conffile conf/t5.ini`, rows `db_t5_*`). `t5/make_t5.py NCBI_BIN --mix-only` writes
the version 5 databases `t5/nmix` and `t5/pmix` (rows `db.t5.mix.*`, blastdbcmd cases `T508`-`T514`): with
`-parse_seqids` and a `-taxid_map` of `t5.tsv`'s taxids, one sequence per record, records of two definition lines
(taxids 100001 and 100002; the first without a taxid; the first untitled and the second titled; the other way round)
and one line of each Seq-id kind whose title NCBI makes from the source (`tpg`, `emb`, `gi` for `-show_gis`, `pdb`
without and with a chain, `lcl|Subject_1`, `gnl`; in `pmix` also `sp`, `ref|NP_`, `ref|WP_`), and the version 4
database `t5/psuffix` of one record whose first line, without a taxid, has a title ending in
`, partial [Bacillus sp.]` (blastdbcmd case `T512`; row `db.t5.psuffix.blastp.0`: the outfmt 0 heading of its untitled
second line, `AB362002.1 alpha [Bacillus sp.]`). The title is made from the Bioseq's first source, the first
positive taxid in `x_GetTaxDefline`'s order of the lines that the taxid list keeps (`-taxids 100002` and
`-negative_taxids 100001` with `-no_taxid_expansion`, which NCBI applies to a version 5 database only), with
`blastdbcmd -target_only` the target's line first and the others after it.

Protein titles ending in `]` (line RG item T5c; rows `db.t5.sfx.*`, `db.t5.sfs.*`, `db.t5.psuffix.*`, blastdbcmd cases
`T515`-`T518`, segmasker rows `db_t5c_*`): `t5/make_t5.py NCBI_BIN --suffix-only` writes `t5/psfx` (version 5,
`-parse_seqids`, `-taxid_map`: records of a titled line and an untitled one, whose outfmt 0 heading and outfmt 18 show
the titled line's title as `GenerateDefline` makes it, and `lcl|Subject_9nn` records, whose tabular ids are its first
word) and `t5/psfs_100002`, `psfs_77777`, `psfs_none` (version 4 without `-parse_seqids`, `-taxid`; the alias
`t5/psfs`: the tabular ids that replace `gnl|BL_ORD_ID|n`). NCBI's `x_AdjustProteinTitleSuffix` rewrites the end of
such a title from the Bioseq's first source, also without a name: a title ending in ` [<name>]` (any case) or ` [ ]`
loses it, then an organelle in parentheses of NCBI's list (`(mitochondrion)`, not `(virus)` or `(chromosome)`) and
`, partial`, and gets ` [<name>]` back where the taxonomy database names the taxid (`alpha, partial [Bacillus sp.]` ->
`alpha [Bacillus sp.]`; without the taxonomy database `zeta protein, partial [ ]` -> `zeta protein`, and the rows
`.nt.` show the rest unchanged); ` []` stays. The description (`fLeavePrefixSuffix`), outfmt 5, 13 to 16, `stitle`
and the stored headings show the title as stored.

Titles without a first word (fix 1 of line RG item T5c; rows `db.t5.word.*`): `t5/make_t5.py NCBI_BIN --word-only`
writes `t5/pword_77777` and `t5/pword_none` (version 4 without `-parse_seqids`, `-taxid 77777` and none; the alias
`t5/pword`) of the titles ` ` (one space), `TPA:`, `MAG: (plastid) [ ]`, `chi, partial [ ]` and `  MAG:  `. Where the
title that `GenerateDefline` makes has no first word (a title of spaces, a prefix alone, `MAG: (plastid) [ ]` cut to
`MAG:` by `x_AdjustProteinTitleSuffix` with the source of taxid 77777), NCBI keeps the sequence's `gnl|BL_ORD_ID|n` in
`sseqid` and `sallseqid` and writes `BL_ORD_ID:n` in `sacc`, `saccver` and `sallacc`; the nucleotide rows search
`t5/nclean` and `t5/nall`, whose empty titles of a taxid that the taxonomy database does not name (or without one)
give no word either.

Taxids in other encodings (fix 2 of line RG item T5c; rows `db.t5.taxid.*`, blastdbcmd cases `T519`-`T524`, segmasker
and dustmasker rows `db_t5c_taxid_*`): `t5/make_t5.py NCBI_BIN --taxid-only` writes the version 4 databases
`t5/ptx_<name>` and `t5/ntx_<name>` (without `-parse_seqids`, `-taxid 9606`, two records with the letters of
`ids/pdup.faa` and `ids/ndup.fna`) and rewrites the `taxid [2]` of the first record (and the header offsets of the
index): `app42` (the application INTEGER tag `0x42`) and `sx9` (nine content bytes, seven of them sign extension),
which NCBI's `ReadInt4` reads as 9606, and `big5` (2^31 in five bytes), on which NCBI fails (`eOverflow`: `TTaxId` is
an `int`) wherever it reads the header: exit 255 in the reports of the record (since line RG item T5d NCBI's rows
`db.t5.taxid.*.big5.*`, which LOSAT refused), `Error: [blastdbcmd] error while reading taxid` (exit 1). Probe:
`remaining-gaps-20261112/T/taxid_t5c.md`.

Header records that NCBI's read fails on (line RG item T5d; rows `db.t5d.*`, blastdbcmd cases `T523`, `T524`, `T527`, `T528`, `T541` and `T542` (`T525` and `T526` before the merge with item T5c), segmasker and
dustmasker rows `db_t5d_*`): `t5/make_t5.py NCBI_BIN --header-only` writes `t5/hd.faa` and `t5/hd.fna` (three records, the
nucleotides the codons of the proteins), the queries `t5/hd_q.*` (query n hits record n) and `t5/hd_q1.*` (record 1 alone), and
the version 4 databases `t5/phd_<name>` and `t5/nhd_<name>` (`-taxid 9606`) whose record of OID 1 is edited: `vstr` (the taxid's
tag a VisibleString), `extra` (data where the taxid's end is expected), `prim` (the taxid's tag primitive), `big5` (2^31), `title`
(the title's tag an OctetString), `cut` (the record without its last four bytes), `seqid` (`-parse_seqids`, `lcl|hd<n>`: the tag of
the local id's string); `title_z` edits OID 0, which `hd_q1.*` never reports. NCBI writes its exception after the outfmt 0
prolog and the reports of the batches before (`BATCH_SIZE=1` rows); `title_z` exits 0. Survey:
`remaining-gaps-20261112/T/SURVEY_T5d.md`. Records that NCBI reads (fix 1 of item T5d, `make_t5.py NCBI_BIN --header-only utf8
utf8id long longx`): `utf8` (the title a UTF8String; `utf8_z` at OID 0), `utf8id` (`-parse_seqids`, the local id's string a
UTF8String), read with NCBI's warning `(803.10) ... UTF8String data for VisibleString member <member>` once per run, where the
Seq-aligns of a batch first read the record (the `BATCH_SIZE=1` rows: in the second batch; `utf8_z` with `hd_q1.*`: none);
`long` (the title's member tag in the long form, `bf 00`), read with no message; `longx` (`bf 20`), NCBI's unexpected member. `eof`
(the record of OID 1 cut after its first byte): NCBI's `CEofException`, which the programs' `CATCH_ALL` takes as a
`CIOException` without a post, so the run ends there with exit status 0.
Probe: `remaining-gaps-20261112/T/fix1_t5d.md`. NCBI's other reads of a UTF8String header (fix 2 of item T5d, rows
`db.t5d.*.fx*`, blastdbcmd cases `T534`-`T537`): `make_t5.py NCBI_BIN --fetch-only` writes the version 5 databases
`pfx_utf8`, `pfx_utf8r`, `nfx_utf8`, `pfx_copy` and `nfx_copy` (`-parse_seqids`; hd0, which `hd_q0.*` matches, a weak
match with a UTF8String title (hd2, a chance alignment of E-value about 6, or a weaker copy of hd0), hd1), the negative seqid
lists `*fx_neg1.txt` and `*fx_neg2.txt`, and the version 4 database `pva_utf8` (`gb|HDA<n>.1|`, the UTF8String title at
OID 1) with `pva_batch.txt`: under a negative seqid list with ids NCBI's traceback and composition-based redo read the
header of every subject they fetch (the redo only those its early termination test leaves), and blastdbcmd's version check
of `HDA1.2` reads HDA1.1's. Probe: `remaining-gaps-20261112/T/fix2_t5d.md`. NCBI's scope's read of another sequence
that holds a reported sequence's Seq-align id (fix 3 of item T5d, rows `db.t5d.*.dup*`, blastdbcmd cases `T538`-`T540`):
`make_t5.py NCBI_BIN --dup-only` writes the version 5 databases `pds`, `pda` (the subject `gb|DUP00001.1|`, hd0),
`pdo_same` (that record, its title a UTF8String), `pdo_none` (`gb|DUP00001.1|`, hd1, a UTF8String title), their
nucleotide `nds`, `ndo_same`, `ndo_none`, the version 4 `pds4` and `pdo4_same`, and `dup_batch.txt`. NCBI orders a
list's volumes by path: `pdo_* pds` puts the other sequence first, whose header the data loader reads (`x_GetOid`)
wherever the scope resolves the subject's id, with the warning; `pda pdo_none` the subject first, whose other holder
is not read. Probe: `remaining-gaps-20261112/T/fix3_t5d.md`. Rows `db.t5d.*.dup.same.pf0.*` (fix 4 of item T5d): with
`PRE_FETCH_SEQS_LIMIT=0` the warning comes inside the report, where it reads the subject (a tabular row, outfmt 0's
description line, a JSON `Hit`; psiblast's PSSM before the first round's report where the report reads nothing),
`db.t5d.blastp.dup.same.15` at the JSON prefetch, `db.t5d.blastp.dup.none.magic.6e` (`BLASTDB=../seqdb/taxdb/bad_magic`)
before the taxonomy database's texts of opening of the same load. Probe: `remaining-gaps-20261112/T/fix4_t5d.md`.
Rows `db.t5d.*.dup5*` (fix 5 of item T5d): `make_t5.py NCBI_BIN --dup5-only` writes `pdo5_same`, `pds5`
(`gb|DUP00005.1|`, hd0, the copy's title a UTF8String), `ppig1` (`gnl|PIG|x001`, hd1: its lookup throws, and its
Seq-id handle, packed 1, loads before DUP00005's, packed 5), `ppig5` (`gnl|PIG|x005`, one packed number with
DUP00005's), `pdo5b_same`, `pdo5b_id`, `pds5b` (`gb|DUQ00005.1|`, hd2, the copy's title or accession a UTF8String),
`pgo5`, `pgs5` (`gi|5`, hd0, the copy's title a UTF8String), `pplain` (no `-parse_seqids`: a list's lookup of a text
id stops NCBI's prefetch), the nucleotide `ndo5_same`, `nds5`, `npig1`, and the queries `hd_q10.*` (q1, q0, q0b) and
`hd_q02.*` (q0, q2). A prefetch whose loads end at the throwing lookup writes `Error pre-fetching sequence data` and
no warning, which the next batch's prefetch writes; the reject rows are orders that NCBI takes from the addresses of
its Seq-id infos and a report that reads the throwing subject. Probe: `remaining-gaps-20261112/T/fix5_t5d.md`.

Rewritten ends that NCBI cannot decode (merge 4 of line RG item T5c with item T7; rows `db.t5.sfu.*`, `db.t5.sfux.*`):
`t5/make_t5.py NCBI_BIN --sfu-only` writes `t5/psfu_100009` (version 4 without `-parse_seqids`, `-taxid 100009`; the
alias `t5/psfu`) of `chi, partial [ ]`, `hypothetical protein [ ]` and `kappa, partial [Homo sapiens]`, and `t5/psfux`
(version 4, `-parse_seqids`, `-taxid_map`) of `lcl|Subject_901 chi, partial [ ]` (taxid 100009), a titled line `gb|AC400031.1|
omega, partial [ ]` followed by the untitled `gb|AD400031.1|` (taxid 100009) and `lcl|Subject_904 delta, partial [Homo
sapiens]` (taxid 100001). With `../seqdb/taxdb/t5`, `x_AdjustProteinTitleSuffix` writes the name of 100009 (`Bad \x81
name`) at the end of a ` [ ]` title, and `GenerateDefline` cannot decode the result: the tabular lookup fails (`Unknown`,
the Seq-align id after the `sall*` lists of the row before), the default fields stop at `saccver` (NCBI exit 255,
LOSAT's rejection) and the outfmt 0 heading of the untitled line is `Sequence with id AC400031 no longer exists in
database...alignment skipped`. Without a taxonomy database the ` [ ]` end goes (`chi`, `hypothetical`). makeblastdb
2.17.0 stores a title byte such as 0x81 as `#`, so a stored database title never carries such a byte itself.

Rewritten ends whose names hold bytes outside US-ASCII or control bytes (merge 5 of line RG item T5c with item K4h; rows
`db.t5.k5.*`, `db.t5.k5s.*`, `db.t5.k5x.*`): `t5/make_t5.py NCBI_BIN --k5-only` writes the taxonomy database
`../seqdb/taxdb/k5` (`k5.tsv`: `taxdb.tsv`'s rows and 5001 `Caf\xc3\xa9 sp.`, 5002 `Ctl\x07 name`, 5003 `Bad \x81 name`,
each with such bytes in all four names; 99990 stays the last entry), the version 4 databases `t5/pk5_<taxid>` (without
`-parse_seqids`, `-taxid` 99990, 5001, 5002, 77777 and 9606; the alias `t5/pk5`), `t5/pk5x` (`-taxid 5003`) and
`t5/pk5s` (`-parse_seqids`, `-taxid_map`: `lcl|Subject_9nn` records and titled lines followed by untitled ones) of titles
ending in `]`, and the queries `t5/k5_qa` and `t5/k5_qc` (`.faa`, `.fna` back-translated): q0 without hits, q1 with the
title `q1 a\xc3\xa9b z` (`qa`) or `q1 c\x07d z` (`qc`) hitting every record, q2 hitting only the records of group A (10
in `pk5`). The tabular word of `MAG: (plastid) [ ]` with 5001 is `[Caf\xc3\xa9`; XML2 writes a control byte of a
`sciname` as `#` with NCBI's 803.7 post and throws on the UTF-8 one and on 0x81 (outfmt 16 goes on); JSON writes them
(`\u0007`). With `../seqdb/taxdb/offset_last_past_end` and `-outfmt 5`, q1's `Iteration` is left out at its title and
the offset error of taxid 99990 comes at q2's load (`PRE_FETCH_SEQS_LIMIT=0`) or at its batch's prefetch (`BATCH_SIZE=1
PRE_FETCH_SEQS_LIMIT=10`). Probe: `remaining-gaps-20261112/T/merge5_t5c.md` (645 cases, NCBI 2.17.0 and LOSAT at 1
and 4 threads).

Subjects and queries whose title NCBI cannot decode (line RG item T7; rows `db.t7.*`): `t7/` holds the inputs of the
survey probe of T7 (its `gen.py`): `nq.fna`, `pq.faa` and their two-query forms `nq2x.fna`, `pq2x.faa`; `nqbad.fna`,
`pqbad.faa` (a query title `qbad x\x81y title`); subject sets `lone` (`Subject_9 bad x\x81y title`), `bg` and `gb`
(that title and a good one), `bb` (two such titles), `plain` (`bad x\x81y title`, no `Subject_` word), `lat`, `latp`
(cp1252 `caf\xe9`), `utf` (UTF-8 `\xc3\xa9`), `good` and `gen` (titles that `GenerateDefline` changes: `Subject_4.`,
`Subject_5 a &amp; b,`, `TPA: s7 third ~`), each as `.fna` and `.faa`. The alias files `t5/nt7_bf`, `t5/pt7_bf` (the
sequence of taxid 100009, whose name has 0x81, then 100001) and `t5/nt7_fb`, `t5/pt7_fb` (the other order) are
`make_t5.py`'s.
The rows `db.t7.b.*`, `db.t7.c.*`, `db.t7.p.*` and `db.t7.j.*` (the T7 fix of the review's probe
`T/t7/review`, its `mkin.py`) add `-parse_deflines` subject sets `lclS9` (`lcl|Subject_9 bad x\x81y title`),
`lclS9g` (that record, then `lcl|s2 good two`) and `pS9tpa` (`lcl|Subject_9 TPA: s7 third ~`), the queries
`qpd_local` (`qbad x\x81y title`), `qpd_gb` (`gb|AB000099.1| qbad ...`), `qpd_lcl` (`lcl|Query_1 qbad ...`), `qpd_gi`
(`gi|4504347|ref|NP_000549.1| qbad ...`) and `qpd_refgi` (`ref|NP_000549.1|gi|4504347| qbad ...`), each as `.fna` and
`.faa`, and the version 5 databases `t7/nlcls9`, `t7/plcls9` (one definition line `lcl|Subject_9`), `t7/ns9gb`,
`t7/ps9gb` (`lcl|Subject_9`, then `gb|AB000001.1|`) and `t7/ngbs9`, `t7/pgbs9` (the other order), each with an empty
title and the taxid 100009 of `../seqdb/taxdb/t5`, written by NCBI's `makeblastdb -parse_seqids -blastdb_version 5
-taxid_map` 2.17.0 from the letters of `t7/bg.fna` and `t7/bg.faa`.

A damaged taxonomy database in dustmasker, segmasker and blastdbcmd (line RG item T5b): `t5b/make_t5b.py` (run from
`LOSAT/` with the directory of NCBI's makeblastdb) writes the version 5 databases `t5b/nm`, `t5b/pm` (volume 1) and
`t5b/nmb`, `t5b/pmb` (volume 2) with `-parse_seqids` and a `-taxid_map` (the letters of `ids/ndup.fna` and
`ids/pdup.faa`), the alias files `t5b/nmv`, `t5b/pmv` over both, and `t5b/batch.txt` (the ids with taxid 99990).
Volume 1: two titled sequences without a taxid, titled and untitled 9606, untitled 99990, titled 150 (the slot of
99990 in NCBI's 256-slot taxonomy cache of a volume), titled 99990, one of two lines 99990 and 150, untitled 77777;
volume 2: untitled 99990, titled 9606. With `../seqdb/taxdb/<variant>` on the search path they serve the rows
`db_t5b_*` of `../dustmasker` and `../segmasker` (`-conffile conf/t5b_<variant>.ini`), the blastdbcmd cases
`T5B101`-`T5B504` (`../blastdbcmd_discovery.rs`) and `../taxdb_tool_texts.rs` (the merged streams, `-logfile`, and
the runs on `btd_empty` where NCBI crashes, which LOSAT refuses). `make_t5b.py --taxdb` (no NCBI program) writes the
taxonomy databases of fix 2, `t5b/taxdb/kingdom_past_page` and `t5b/taxdb/names_past_page`: the `taxdb.btd` of
`../seqdb/taxdb/good` (266 bytes) with an index whose record of 9606 reaches past NCBI's mapping of the file (the end
of its 4,096-byte page), after its third tab (68..5,000: its names inside the mapping, its kingdom not) or wholly
(4,096..5,000). NCBI reads memory past its mapping there; with them on the search path the databases serve the rows
`db_t5b_*_past_page_*` of the maskers, the blastdbcmd cases `T5B601`-`T5B608`, the rows `db.t5b.*` of
`../db_search_manifest/blastp.tsv` and `blastn.tsv`, and the refusals of `../taxdb_tool_texts.rs` (a name written from
past the mapping). Fix 3 adds `t5b/taxdb/zero_kingdom_past_page` (an entry of taxid 0 whose record runs past the
mapping after its third tab: the blastdbcmd cases `T5B609`-`T5B610`, `-tax_info` writes `N/A` for taxid 0) and
`t5b/taxdb/dup_last_btd_empty` (an empty `taxdb.btd` with good's index and its last taxid again after it: the search
for that taxid reads a record, so no lookup writes the offset error, and NCBI faults at the first lookup that reads a
record; the rows `db.t5b.de.*` of the six programs, LOSAT's refusals). The rows `db.t5b.ldo.*` of the six programs
(merge 2 with line DBID loadord) put `../seqdb/taxdb/btd_empty`, `t5b/taxdb/names_past_page` and
`t5b/taxdb/dup_last_btd_empty` on the search path of loadord's batch whose prefetch NCBI stops at a Seq-id lookup
(`multidb/ldl` or `ldb` with `ldp`, and their nucleotide copies): the loads before the stop write the offset error,
nothing, or make NCBI fault (LOSAT's refusal), then NCBI's error.

A title of punctuation without parsed Seq-ids (line DBID tabfields; rows `db.id.tf.*.punct.*`): `make_ids.py
--tabfields-only` writes `ids/ppunct_np` (version 4, without `-parse_seqids`) from `ids/ppunct.faa`, whose first title
`, ,` NCBI's `x_CleanAndCompress` reads past the end of where it makes the tabular id of a protein; the query
`ids/ppunct_q.faa` hits both records. A table that prints no subject id (`6 qseqid evalue`) does not make it, so LOSAT
writes NCBI's rows; `sseqid` is LOSAT's refusal.

Each batch's prefetch counts its own Seq-aligns (line DBID pfbatch; rows `db.id.pfb.*`): `make_ids.py --pfbatch-only`
writes `ids/nordtag_v5` (version 5, `-parse_seqids`) from `ids/nordtag.fna`, a sequence whose first id
`gnl|BL_ORD_ID|0` makeblastdb stores with a text tag, so NCBI's lookup by that Seq-align id fails (as in `pshare_v5`),
and another sequence. The queries `ids/nordtag_q.fna` and `ids/nordtag_q.faa` (qa, qb) hit one each; `ids/pfb_aba.faa`
(qa, qb, qa2: with `BATCH_SIZE=101` the batches [qa, qb] and [qa2]) and `ids/pfb_aab.faa`, `ids/pfb_aab.fna` (qa, qa2,
qb) put the sequence with the text tag in a batch of more Seq-aligns than the limit and in one within it, or only in the
first; `ids/pfb_s.fna` holds two pieces of it (two HSPs gapped, one Seq-align with `-ungapped`). `ids/pfb_rc_nambig.fna`
and `ids/pfb_rc_ids_q.fna`, the reverse complements of `ids/nambig.fna` and `ids/ids_q.fna`, align on the minus strand
of the query with `lcl|777` of `ids/nrefuse_v5`, whose Bioseq NCBI's prefetch does not get. With a positive
`PRE_FETCH_SEQS_LIMIT`, a batch whose Seq-aligns outnumber it reads nothing in its prefetch, so a table that reads
nothing of the subject writes NCBI's rows (NCBI writes no warning); a batch within it reads the subject, and where its
prefetch fails NCBI and LOSAT write `Error pre-fetching sequence data` (line DBID pfwarn).

The warning of a batch whose prefetch fails (line DBID pfwarn; rows `db.id.pfw.*`): `make_ids.py --pfwarn-only` writes
`ids/ntag_t5` (version 5) from `ids/nordtag.fna`, `ids/nlook_t4` (version 4) from `ids/nlook.fna` and `ids/pshare_t5`
(version 5) from `ids/pshare.faa`, every sequence with taxid 99990, whose lookup in
`tests/seqdb/taxdb/offset_last_past_end` writes the offset error. NCBI's prefetch writes the taxonomy database's texts of
the subjects it loads, then `Warning: [<program>] Error pre-fetching sequence data`, once for the batch, before the
batch's query warnings. A minus-strand Seq-align of `bbs|123` (no ISAM key in version 4, so no Bioseq) fails after every
load (`ids/pfw_look_fw_rc.fna`, `ids/pfw_look_rc_fw.fna`: `fw70` on the plus strand of `gb|AF000070.1|`, `rc123` on the
minus strand of `bbs|123`); the lookup of `gnl|BL_ORD_ID|0` with a text tag fails where the loads reach it, in the
order of NCBI's Seq-id handles (the ids its mapper packs first), which LOSAT does not model, so a batch with it whose
other subjects' loads write a text is LOSAT's refusal. In psiblast, where the report reads the sequences it shows
(`pfw_weak_v5` with `pfw_q6.faa`, `-max_target_seqs 3`: the text tag ranks last, beyond the cut), the failed prefetch
leaves `sp||NAMEY` unloaded and NCBI's formatter reads it one by one, matching `sp|P00042.1|NAMEY` weakly: LOSAT refuses
there, and where the PSSM (`-out_ascii_pssm`) holds such a sequence (NCBI drops its letters). `pfw_bad/ntag_t5` fails the
lookup of an accession on its file, which NCBI's prefetch rethrows (exit 2): LOSAT's refusal. With `-ungapped` BLASTN makes one Seq-align of a subject's HSPs, whose strand is its first
HSP's: `ids/pfw_fwL_rcS.fna` aligns with `bbs|123` on the plus strand first (240 letters) and on the minus strand (60),
`ids/pfw_fwS_rcL.fna` the other way round. `ids/pfw_ps_q21.faa` holds `pshare_q2.faa` and `pshare_q1.faa`.

Copies of one record that follow each other in a hit list (line DBID N3; rows `db.id.join.*`): `make_ids.py
--join-only` writes, from the FASTA of `ids/make_join_fasta.py`, `ids/pjoin_a`, `ids/pjoin_b`, `ids/pjoin_c` and
`ids/njoin_a`, `ids/njoin_b`, `ids/njoin_c` (version 5; r1 in all three, r3 with two definition lines and two
pieces of the query in a and b, r6 the letters of r1 under another id in b), `ids/pjoin_bt` and `ids/njoin_bt` (the b
records with taxid 9606), and the aliases `ids/pjoin_abc.pal`, `ids/pjoin_abt.pal`, `ids/njoin_abc.nal` and
`ids/njoin_abt.nal` (queries `pjoin_q.faa`, `njoin_q.fna` and `njoin_q2.fna`, the last with two records for two
query batches). NCBI orders tied hits by the later ordinal id first: r1 of c, r6, then r1 of b and a, which outfmt 0,
5 and 13 to 16 show in one hit (r3 of b and a too), as LOSAT does; `-num_descriptions`, `-num_alignments` and
psiblast's cut to `-max_target_seqs` count them as one subject, the tabular `qcovs` and `qcovus` take the group's
coverage. `ids/pjoin_x` and `ids/njoin_x` hold the definition lines of r1 and r3 with other letters (no copies,
refused; aliases `ids/pjoin_ax.pal`, `ids/njoin_ax.nal`); `ids/njoin_ta` and `ids/njoin_tb` (alias `ids/njoin_tab.nal`)
hold r10 with taxid 9606 and a second definition line `lcl|loc10`, which the outfmt 0 heading under `-taxids 9606
-no_taxid_expansion` hides where each copy has one HSP.

Sequences whose only Seq-id is `gnl|BL_ORD_ID|n` (line DBID N2a; rows `db.id.*.ord.*`): `make_ids.py` writes, from
the FASTA of `ids/make_ord_fasta.py`, `ids/pord_v5` and `ids/pord_v4` from `ids/pord.faa` (two accessions, then
`gnl|BL_ORD_ID|2`, its own number, `gnl|BL_ORD_ID|6`, the number of no sequence, and `gnl|BL_ORD_ID|1`, the
number of a sequence with an accession; each hit by its query `ids/pord_q<k>.faa`, `.fna` for BLASTX), `ids/pord_x`
from `ids/pord_x.faa`, read after `pord_v5` through the alias `ids/pord_vx.pal`, which lists it first (NCBI orders
the volumes of an alias by their base names; there `gnl|BL_ORD_ID|6` of `pord_v5` is the number of the second sequence of `pord_x`, whose
`gnl|BL_ORD_ID|1` is its own number in its volume), and `ids/nord_v5` and `ids/nord_v4` from `ids/nord.fna`
(`gnl|BL_ORD_ID|3` first, the number of no sequence, then `gnl|BL_ORD_ID|1` and an accession; queries
`ids/nord_q<k>.fna` and `.faa`, and `ids/nord_qall.fna` with the three in reverse order). NCBI looks n up
plus the first ordinal id of the volume, so its own number finds the sequence, which LOSAT reads as one without
parsed Seq-ids; another n finds no sequence or another one, which LOSAT refuses in the formats that read the
Bioseq. No record names `gnl|BL_ORD_ID|0`: makeblastdb stores it with a string tag.

Aliases that restrict a database with parsed Seq-ids by `FIRST_OID` and `LAST_OID` (line DBID ordnum; rows
`db.id.on.*`): `make_ids.py` writes `ids/pfound_r.pal` (OIDs 1 to 9 of `ids/pfound_v5`), `ids/refuse_r.pal` (OIDs 1
to 5 of `ids/refuse_v5`) and `ids/njoined_r.nal` (`nfound_v4`, `njoined_a` and `njoined_b`, OIDs 2 and 3 only: the
pair of `ids/njoined`), and the list `ids/pfound_nsl.txt` (AAA00001.1, OID 0 of `ids/pfound_v5`, for
`-negative_seqidlist`). The rank of a searched sequence is its OID less the OIDs left out before it; LOSAT's refusals
name the OID, as `blastdbcmd -outfmt %o` prints it (with `taxfilter/tfe_nr` and `taxfilter/tfe_pr`, OIDs 1 to 3).

Seq-align ids that NCBI's formatters look up like any other (line DBID N2b; rows `db.id.*.look.*`): `make_ids.py`
writes, from the FASTA of `ids/make_look_fasta.py`, `ids/plook_v5` and `ids/plook_v4` from `ids/plook.faa` (an
accession, `bbs|123`, `nat|NA000001.1|`, `lcl|777`, `gnl|PIG|124`, `bbs|300` then `lcl|300` and `lcl|301` then
`bbs|301` (one key each), `gnl|BL_ORD_ID|n` beside an accession with n its own number on two definition lines and on
one, the number of the PIG sequence and of no sequence, `bbs|500|gb|...|`, `lcl|778` beside an accession; each hit by
its query `ids/plook_q<k>.faa`, `.fna` for BLASTX), `ids/plook_z` from `ids/plook_z.faa`, read after `plook_v5` through
the alias `ids/plook_vz.pal` (its second sequence holds `gnl|BL_ORD_ID|1`, its own number in its volume), and
`ids/nlook_v5` and `ids/nlook_v4` from `ids/nlook.fna` (queries `ids/nlook_q<k>.fna` and `.faa`, `ids/nlook_qall.fna`
and `ids/nlook_qfound.fna` for the report streams). Version 5 finds `bbs|`, `nat|` and `lcl|<integer>` by their LMDB
keys unless an earlier sequence holds the key; version 4 finds neither `bbs|` nor `nat|`; no version finds `gnl|PIG|`.
A Seq-align id that another sequence also holds (review 1 of N2b): `ids/pshare_v5` and `ids/pshare_v4` from
`ids/pshare.faa` (`lcl|555`, then `gb|FFF00030.1|` with `lcl|555` in a second definition line and nearly the same
letters, so that the query `ids/pshare_q1.faa` hits both and NCBI's scope shows the second for the first; a
`gnl|BL_ORD_ID|0`, which makeblastdb stores with a text tag, query `ids/pshare_q2.faa`), `ids/pshare_a` and
`ids/pshare_b` (`bbs|123` in each) through the alias `ids/pshare_ab.pal` (query `ids/pshare_b_q0.faa`), and
`ids/nshare_v5` from `ids/nshare.fna` (query `ids/nshare_q1.fna`). Ids that differ in case alone (review 2 of N2b):
`ids/pshare.faa` 4 `lcl|abc` and 5 with `lcl|ABC` (query `ids/pshare_q5.faa`), and `ids/pshare_c` (`pdb|1ABC|A`,
`gnl|db|abc`) and `ids/pshare_d` (`pdb|1abc|A`, `gnl|DB|ABC`) through the alias `ids/pshare_cd.pal` (queries
`ids/pshare_d_q<k>.faa`); `ids/pdup_abc.pal` reads the copies `pdup_a` and `pdup_b` and the record of another title
`pdup_c`.

Seq-ids that NCBI's scope matches weakly (line DBID N2c; rows `db.id.*.weak.*`): `make_ids.py` (`--weak-only` builds
just these) writes, from the FASTA of `ids/make_weak_fasta.py`, `ids/pweak_v5` and `ids/pweak_v4` from `ids/pweak.faa`
(`sp|P00041.1|NAMEW`, `gb|FFF00040|`, `pir|PIRW1|`, `sp||NAMEX`, `gb|FFF00041.1|`, then three pairs S, T of nearly
equal letters: `sp||NAMEY` before `sp|P00042.1|NAMEY`, `gb|FFF00043|` before `gb|FFF00043.1|`, `sp|P00044.1|NAMEZ`
before `sp||NAMEZ`) and `ids/nweak_v5` and `ids/nweak_v4` from `ids/nweak.fna` (`gb|AF000081|` before
`gb|AF000081.1|`, then `gb|AF000082|`). NCBI's scope finds an accession without a version by the same one with a
version, an id of accession and name by the name-only id, a name-only id by an id of that name and an accession, and a
pir accession-only id by the name-only id; not the other way. The queries `ids/pweak_qa<k>.faa` (`.fna` for BLASTX,
`ids/nweak_qa2.fna`) carry such an id and are read with -parse_deflines (`qa4` the direction NCBI does not match);
`ids/pweak_q<k>.faa`, `ids/pweak_q8.fna`, `ids/nweak_q1.fna` and `.faa` hit record k (the second of a pair also hits
the first) and are run with `PRE_FETCH_SEQS_LIMIT=N` in `extra_args`, an environment variable of both programs (NCBI
reads the reported sequences one by one when N is 0 or below the batch's Seq-align count and then writes the names of
the sequence it read first). LOSAT refuses both cases, the second in a batch above N (line DBID weakbatch, rows
`db.id.wkb.*` below);
in a version 5 database LOSAT does not reproduce NCBI's lookup and refuses every Seq-align id with a name and every pir or prf id
without it, in version 4 the lookup of the id decides (a name alone with no other holder is written as NCBI does:
`db.id.blastp.weak.v4.pf.alone.6`, `v4.pf.none.6`). The `ncbi` rows are the searches without the variable or at 2147483647, the
version 4 sequences with no weak holder and outfmt 8. Review 1 of N2c: `ids/pweak_r_v5` and `ids/pweak_r_v4` from `ids/pweak_r.faa`,
four pairs S, T (T's first definition line holds its own id and the second id that S's matches weakly: `gb|FFF00051.1|`
with `sp||NAME_R` for `sp|P00051.1|name_r` (another letter case), `gb|FFF00053.1|` with `prf||NAME_S` for
`tpg|FFF00052.1|NAME_S`, `gb|FFF00054.1|` with `prf||PIRW3` for `pir|PIRW3|`, `gb|FFF00055.1|` with `emb|S12347.1|` for
`pir|S12347|`; queries `ids/pweak_r_q<k>.faa` are T's letters; rows `db.id.blastp.weak.r.*` and `.v4.r.*`: NCBI writes T's names
for S at `PRE_FETCH_SEQS_LIMIT=0` in both versions, and both its own rows without the variable).

Rows `db.id.wkb.*` (line DBID weakbatch): NCBI counts a positive `PRE_FETCH_SEQS_LIMIT` against each query batch's
Seq-aligns. A batch at or below it prefetches its sequences, each by its own id, and shows each pair of the rows above as
two sequences (the `ncbi` rows, in all six programs, versions 4 and 5, outfmt 0 and 5); a batch above it reads them one by
one and shows T for S (the rejection rows). `ids/wkb_q8x2.faa` and `ids/wkb_n1x2.fna` are `pweak_q8.faa` and
`nweak_q1.fna` twice, joined by 30 `X` and 40 `N` (two HSPs a subject, so 4 Seq-aligns in a gapped search), and
`ids/wkb_q0_q4x3.faa` is `pweak_q0.faa` followed by `pweak_q4.faa` three times (with `BATCH_SIZE=1`, a batch of 1 Seq-align
that prefetches beside one of 3 that does not). The psiblast rows `db.id.wkb.ps.*` run rounds (`-num_iterations 2`):
each round that makes the next round's PSSM reads the sequences before the cut by their Seq-align ids, so
`pir|S12347|` beyond `-max_target_seqs 1` (`pweak_r_v5`) is refused where the round's prefetch is off, and two
following hits whose Seq-align ids `Match` (`gb|FFF00043.1|` and `gb|FFF00043|` of `pweak_v5`), which NCBI puts in one
row of the PSSM's multiple alignment, make that PSSM as NCBI does (`db.id.wkb.ps.join.i2.6`, NCBI's bytes since line DBID
pssmrows; weakbatch refused it). `ids/wkb_e3_v5` (`make_ids.py
--weakbatch-only`, from `ids/wkb_e3.faa`) holds `gb|FFF00043.1|`, `gb|ZZZ00001.1|` and then `gb|FFF00043|`, which the
lookup finds as the first: beyond `-max_target_seqs 2` it is in round 1's PSSM, which reads the first one's letters
for it (`db.id.wkb.ps.e3.s1.7`, NCBI's bytes since line DBID pssmread; refused before); `ids/wkb_p1_v5` makes it a weak hit above `-inclusion_ethresh` behind
`gb|WWW00001.1|`, which the PSSM does not read (`db.id.wkb.ps.p1.i2.7`, NCBI's bytes); the query is `ids/wkb_qy.faa`.

Rows `db.id.psr.*` (line DBID, pssmrows; databases `ids/psr_*` from `ids/psr_*.faa` by `make_ids.py --pssmrows-only`, all
version 5, query `ids/psr_q.faa`, a random protein of 250 letters; `ids/psr_cpa` holds three variants of the query with 62, 72 and
82 substitutions): psiblast counts the rows of the PSSM's multiple alignment by Seq-align id handles, and a
Seq-align whose id `CSeq_id::Match`es the previous one's joins that row. `psr_cp.pal` (`psr_cpa`, then `psr_cpb`, a copy
of its first record) puts two copies of one record after each other: one row, as NCBI (`db.id.psr.cp.*`, round 2 and the
ASCII PSSM); `psr_cbxp.pal` puts another subject between the copies, where NCBI writes round 1's report and crashes
(SIGSEGV, rc 139) and LOSAT refuses after the same report (`db.id.psr.cbx.i2.6.reject`, stdout frozen from NCBI).
`psr_jm` (`gb|HHH00001|`, `gb|HHH00001.1|`, a third record) holds a Match pair whose second record fills only the cells
after the first's alignment (`db.id.psr.jm.*`: round 2, the ASCII and ASN.1 PSSMs, and `-inclusion_ethresh 1e-90`, which
leaves the second out); `psr_jn` the same with `sp||NAMEH` and `sp|P00009.1|NAMEH`; `psr_jx.pal` (`psr_jx1`, `psr_jx2`) the
pair in two databases (`-num_iterations 3`); `psr_jctl` the same letters under ids that never Match (the control: another
PSSM); `psr_self` a Match pair both identical to the query, purged as self hits (the empty row). `psr_cul.pal` (`gb|ABC00001.1|`,
then `gb|abc00001.1|`) and `psr_cue.pal` (`gb|ABC00001.1|`, then `emb|ABC00001.1|`) are handles that NCBI can take for one
sequence and LOSAT refuses (`db.id.psr.case.i2.6.reject`, `db.id.psr.gbemb.i2.6.reject`). `psr_k.pal` (`psr_ka`:
`gb|SSS00001.1|` identical to the query, `gb|AAA00001.1|`, `gb|ABC00001.1|`; `psr_kb`: `emb|ABC00001.1|`) puts that pair in
a search where neither is included in the PSSM nor shown, so NCBI's handles never compare them and LOSAT does not refuse
(`db.id.psr.k.i2.6`, `-inclusion_ethresh 1e-150 -max_target_seqs 1`). `psr_gim` holds `gim|123` (a giim id) in front of
`gb|BBB00001.1|` and `gb|CCC00001.1|`: NCBI's handles never take a giim id for one with an equal id, so NCBI writes
`Unknown` for it and does not converge in 4 rounds, and LOSAT refuses (`db.id.psr.gim.i4.6.reject`).
Review 2 (`ids/psr_gz.faa`: `gb|SSS00001.1|` identical to the query, `gb|AAA00001.1|`, `gim|123`; `ids/psr_gimn.fna`: `gim|456`
and `gb|NNN00001.1|`, queries `ids/psr_qp.fna` and `ids/psr_qm.fna`): `db.id.psr.last.k2.*` is `psr_k` with
`-inclusion_ethresh 1e-135 -max_target_seqs 2` (outfmt 6 and 0, the last row of the PSSM's multiple alignment) and
`db.id.psr.gz.r2m2.6` the same search of `psr_gz`. A giim id is refused wherever a report reads its sequence:
`db.id.psr.gim.bp6q`, `bnm6q` and `tbx6q` (`qseqid`, which reads nothing of the subject; the blastn row has NCBI's
`Error pre-fetching sequence data` warning in its `.err`) are NCBI's bytes, and `db.id.psr.gim.bp6.reject`,
`bnm6.reject` and `ni.reject` (psiblast with three rounds; `sseqid`, where NCBI writes `Unknown`) are LOSAT's refusals.
psiblast's PSSM reads the letters of the subjects it includes (below `-inclusion_ethresh`), not those of every reported
subject: at `-inclusion_ethresh 1e-140` the giim id of `psr_gim` is not included, so `db.id.psr.gim.r3q.6` (three
rounds) and `db.id.psr.gim.s1q` (round 1 and its ASCII PSSM) are NCBI's bytes. `psr_jv` (`gb|HHH00001.1|`, then
`gb|HHH00001|`) and `psr_vm` (`gb|AAA00001.1|`, then `gb|AAA00001|`) hold a sequence whose Seq-align id the lookup
does not give back as that sequence; it is not included at `-inclusion_ethresh 1e-100` (`psr_jv`) and `1e-150`
(`psr_vm`), so `db.id.psr.jv.r3q.6` and `db.id.psr.vm.r3q.6` (three rounds, `qseqid`) are NCBI's bytes. At the default
threshold `psr_jv` includes it: NCBI warns `Failed to retrieve sequence gb|HHH00001|` (its lookup gives the shorter
`gb|HHH00001.1|`) and the row, joined with `gb|HHH00001.1|`'s, becomes the query's; `db.id.psr.jv.r3q.i.6` is NCBI's bytes
(pssmread writes it).

Rows `db.id.prd.*` (line DBID, pssmread; databases `ids/prd_*` from `ids/prd_*.faa` by `make_ids.py --pssmread-only`, all
version 5 but `prd_as1_4`; queries `ids/prd_q.faa`, `ids/prd_q2.faa` (two queries), `ids/prd_q3.faa` (three, the second with an
`O`) and the PSSM query `ids/prd_q.asn.gz`, NCBI's `-out_pssm` of `prd_q.faa`): psiblast's PSSM reads the letters of each
included subject by its Seq-align id, so a record whose id the lookup gives back as another sequence Y (A is
`gb|AAA00001|`, Y `gb|AAA00001.1|`, listed first) has Y's letters in its row if A's HSP range fits in Y, and otherwise NCBI
warns `Failed to retrieve sequence gb|AAA00001|` once for each such HSP and the row is the query's. `prd_av` (Y of A's
length: read; `db.id.prd.av.s1`, the ASCII PSSM, and `.r3`, three rounds), `prd_ad` (Y as `dbj|AAA00001.1|`: read, no
warning; `db.id.prd.ad.s1`), `prd_as1` (Y one letter shorter than A's HSP: one warning; `db.id.prd.as1.*`: round 1's ASCII PSSM
(`s1`), two queries over three rounds (`r3q2`; each PSSM's warning comes after its round, before the next round's), two
queries and the last round's PSSM (`s2q2`; a warning at each PSSM), three queries (`s1q3`, `r2q3`; each query's PSSM warning
comes after its report, before the next query's warnings) and an `-in_pssm` search (`ip2`; the PSSM after round 1 warns)),
`prd_as1_4` (the same in a version 4 database; `db.id.prd.as1v4.s1`), `prd_ar` (Y as `ref|AAA00001.1|`, whose ids do not
compare equal to A's: the loader keeps none and NCBI warns; `db.id.prd.ar.r3`) and `prd_t2` (A in two HSPs over a short Y:
two warnings; `db.id.prd.t2.s1`). The `.err` files are frozen (the warning counts are 0 for `av` and `ad`, 1 for `s1`,
`as1v4.s1`, `ip2` and `ar.r3`, 2 for `r3q2` and `t2.s1`, 3 for `s1q3` and `r2q3`, 4 for `s2q2`).

Review 1 of pssmread (rows `db.id.prd.*` over `ids/prd_unf.pal`, `prd_zunf.pal`, `prd_otag`, `prd_otag_4`, `prd_wk5`, `prd_ah4`
and `prd_txself`, query `ids/prd_q.faa`, all psiblast against the same single query): `unf.s1`, `unf.r3` read an unfetchable
subject (`lcl|777` before `gi|777`) whose lookup gives Y (`prd_u0y`, listed first) first, so NCBI reads Y's letters and does
not warn; `zunf.s1`, `zunf.r3` do the same where Y's kept definition line is `lcl|777` and does not hold the id, so NCBI warns
`Failed to retrieve sequence gi|777`; `otag.s1`, `otag.pf0.s1` (`PRE_FETCH_SEQS_LIMIT=0`) and `otag4.s1` (version 4) store
`gnl|BL_ORD_ID|0` with a text tag, whose lookup throws: NCBI writes `Error pre-fetching sequence data` (prefetch on) and
`Failed to retrieve sequence gnl|BL_ORD_ID|0` and exits 0. `txself.pf0` (A itself with taxid 99990 by `prd_txself.map`,
`BLASTDB=tests/seqdb/taxdb/bad_magic` and the prefetch off: NCBI's PSSM read loads the Bioseq and writes the taxonomy
database's text after the report; a refusal until review 3 of the line). Two `losat-rejection` rows (rc 1, with the `pssm`
files kind: no PSSM file is written): `wk5.pf0.reject` (the prefetch off, a weakly matching sequence C that NCBI's scope has
read before: NCBI reads C's letters) and `ah4.reject` (another sequence that also holds A's id: NCBI's result depends on its
reads and memory layout).

Review 2 of pssmread (rows `db.id.prd2.*`; databases from `ids/prd_*.faa` and `ids/prdn_*.fna` by `make_ids.py
--pssmread-only`, read through `ids/prd_xg.pal`, `prd_tw.pal`, `prd_tu.pal`, `prd_pah.pal`, `prdn_u.nal` and
`prdn_z.nal`; the taxonomy databases of `tests/seqdb/taxdb`): `txself.tax.pf0.s1`, `txself.sid.pf0.r2` and
`txself.lastoff.pf0.0` read A (taxid 99990) by the query's report one by one (`staxids`, `sseqid`, outfmt 0, the prefetch
off) before the PSSM, so the taxonomy database's text (a wrong magic number, the offset error of taxid 99990) comes with
the report's read and the PSSM's read loads nothing; `xg.s1`, `xg.r3`: A `gnl|DB|SEQA`, whose lookup gives Y
(`gnl|DB|OTHER ^A gnl|DB|SEQA`) first: NCBI's reader keeps Y's second definition line, which holds the id, so the PSSM
reads Y's letters; `tw.pf0.s1`: the lookup gives Y (`gb|AAA00002.1|gb|AAA00001.1|`, taxid 99990), whose definition line
the reader's filter does not keep for A's `gb|AAA00001.1|` (`s_SeqDB_SeqIdIn` stops at the line's first comparable id,
`gb|AAA00002.1|`): `GetBioseq` throws before its taxid lookup, so NCBI writes
`Failed to retrieve sequence` and no taxonomy text; `tu.bad.s1` and `tu.lastoff.r1`: A `lcl|777 ^A gi|777` (unfetchable by
its own ids) after Y (`gi|777|gb|YYY00001.1|`, 20 letters, taxid 99990): the prefetch reads Y and looks its taxid up, and
the PSSM warns `Failed to retrieve sequence gi|777`; `tu.bad.pf0` (the same with the prefetch off: the PSSM's read loads Y
first, so NCBI writes the text, then the warning; a refusal until review 3); one `losat-rejection` row, `pah.reject` (C
`gi|888|gb|CCC00001.1| ^A gi|777` also holds A's `gi|777`: the refusal names database sequence 2). BLASTN (megablast) and TBLASTX with the
minus-strand query `ids/prdn_qrc.fna` (the prefetch asks the length of a minus-strand Seq-align's subject): `bn.u.rc`,
`tx.u.rc` (A `lcl|777 ^A gi|777` after Y `gi|777|gb|YYY00001.1|`: the scope gets Y's Bioseq, no warning) and `bn.z.rc`,
`tx.z.rc` (A `gi|777` after Y `lcl|777 ^A gi|777`, whose kept line `lcl|777` does not hold the id: `Error pre-fetching
sequence data`). `ids/prd_mb.pal` (`prd_mb0ys`: Y `gb|BBB00001.1|` of 100 letters; `prd_mbab`: A `gb|AAA00001.1|`, B
`gb|BBB00001|`, C) is `tests/psiblast_pssm_read_order.rs`'s: with `-max_target_seqs 1` no report reads B, whose PSSM read
fails, and the last query's warning comes before the epilog where both streams go to one file.

Where the reports read their subjects one by one (line DBID taxrow; `ids/txr_p` and `ids/txr_n` by `make_ids.py
--taxrow-only`, a decoy `gb|DDD00001.1|` without hits, A `gb|AAA00001.1|` taxid 9606, the best hit of `ids/txr_q.faa` /
`ids/txr_q.fna`, and B `gb|BBB00001.1|` taxid 99990; `ids/txr_q2.*` end with a query without hits): with `PRE_FETCH_SEQS_LIMIT=0` NCBI's formatter loads each subject when it formats it, so the offset error of
`tests/seqdb/taxdb/offset_last_past_end` comes between A's and B's tabular rows, description lines and BLAST XML `Hit`s
where both streams go to one file; the organism report (outfmt 18) reads each subject after the `Query=` block, before the
lookups of its taxids. `tests/taxonomy_read_order.rs` compares LOSAT's merged stream with NCBI 2.17.0's
(`tests/fixtures/taxonomy_read_order`) in the six programs, with `ids/prd_txself` for a plural name's open, and
`ids/txr_s` and `ids/txr_sn` (`make_ids.py --taxrow-only`; F `gb|FFF00001.1|` taxid 99990, a fragment of the query that
`-sorthits 3` puts before A) for the order of the pairwise report's reads (line DBID taxrow2): `bad_magic`'s open comes
before F's description line, where the sorted report reads first.

Review 3 of pssmread (rows `db.id.prd3.*`; `ids/prd_np` and its version 4 copy `ids/prd_np_4` from `ids/prd_np.faa`, A and B
without `-parse_seqids` and taxid 99990 by `-taxid 99990`, and `ids/prd_tc` from `ids/prd_tc.faa`, A taxid 99990 and B taxid
150 by `prd_tc.map`, all by `make_ids.py --pssmread-only`; the prefetch off and a report that reads nothing of the subjects):
psiblast's PSSM reads a subject without parsed Seq-ids by `gnl|BL_ORD_ID|n`, whose load looks its taxid up, so where that
is the run's first lookup NCBI writes the taxonomy database's text after the query's report, before the PSSM's warnings and
the next query's. `np.bad.pf0.s1` (a wrong magic number), `np4.lastoff.pf1.q3` (three queries, the offset error of taxid
99990 after the first query's report and before the second query's O warning), `np.count.pf0.r2` (two rounds, a count larger
than the entries: after round 1's report, before round 2's composition warning) and `np.bad.pf0.stop` (`-pseudocount -1`:
the text, then NCBI's PSSM error, exit 1) are NCBI's bytes; `tc.lastoff.pf0.q2.reject` is refused: each query's PSSM
reads B, then A; the second query reads B again, whose taxid takes the same slot of NCBI's 256-slot taxonomy cache of a
volume only where NCBI's scope loads B again, so A's read again writes the offset error or not as the scope decides (it
loads them again after 40 other included subjects, not after none).

Review 4 of pssmread (rows `db.id.prd4.*`; `ids/prd_m`, `ids/prd_m2` and `ids/prd_ax` from the FASTA and maps of the same
names by `make_ids.py --pssmread-only`, `BLASTDB=tests/seqdb/taxdb/offset_last_past_end`, the prefetch off and a report that
reads nothing of the subjects): NCBI's volume cache keeps one taxid per slot (`taxid & 255`), and each lookup of taxid 99990
that finds another taxid of its slot (150) there writes the offset error again, in the order of the PSSM's reads and of the
definition lines of each load. `m.lastoff.pf0.s1` (A's lines 99990, 150, 99990: twice), `m2.lastoff.pf0.s1` (A 99990 and 150,
then B 99990: twice), `ax.lastoff.pf0.q2` and `ax.lastoff.pf0.r2` (A 99990 and X 150, which no PSSM includes: once, the second
query's or round's read of A writes nothing whether NCBI's scope loads A again or not) are NCBI's bytes;
`m2.lastoff.pf0.q2.reject` is refused: the second query reads A again, which puts 150 in the slot only where NCBI's scope
loads A again, so B's read again writes the offset error or not as the scope decides.

## The taxonomy database's prefetch messages (line DBID, item T, review 1)

Rows `db.tm.bp.fast.*`, `db.tm.bp.2vol.*` and `db.tm.*.cache*`, all with `BLASTDB=tests/seqdb/taxdb/offset_last_past_end`
(the last entry, taxid 99990, begins past the end of `taxdb.btd`):

- `tm_fast_q289.faa`: 289 copies of `qp1` of `x16_taxid_q.faa`, the first with its first letter replaced by `O` (a warning
  of its own). With `-task blastp-fast` against `x16_taxid_99990` the search makes two query blocks; the offset error
  must come before that warning, at any thread count.
- `tm_2vol.pal` (`DBLIST tm_2vol_a tm_2vol_b`): two volumes that are copies of `x16_taxid_99990` (`.phr .pin .psq`); NCBI
  writes the offset error once per volume.
- `tm_cache.faa` (`st1`, `st2`, `st3`: the same protein), `tm_cache_q.faa` (its query), and the databases built from them
  with every file makeblastdb writes (version 5):

      makeblastdb -in tm_cache.faa -dbtype prot -out tm_cache150 -title "tm cache 150" \
          -parse_seqids -taxid_map tm_cache150.map
      makeblastdb -in tm_cache.faa -dbtype prot -out tm_cache151 -title "tm cache 151" \
          -parse_seqids -taxid_map tm_cache151.map

  (run in this directory). `tm_cache150` holds the taxids 99990, 150, 99990: 150 takes the slot of 99990 (`99990 & 255`) in
  NCBI's 256-slot taxonomy cache of a volume (`CSeqDBVol::m_TaxCache`), so NCBI writes the offset error again after a fetch
  of 150, as often as its order of fetches makes it (twice for the row's query, once for the same query with a version 4
  database or an absolute `BLASTDB`). NCBI's prefetch loads the batch's subjects in the order of their Seq-id handles; the
  mapper does not pack these local ids, so the order is that of the handles' addresses, and LOSAT refuses the run (line NP
  item O4b2-3). `tm_cache151` (151 instead of 150) is the control: once.
- `tm_cache_gi.faa` and `tm_cache_gi.map` (line NP item O4b2-3): 12 copies of the protein of `tm_cache.faa` with the ids
  `gi|7`, `gi|3`, `gi|12`, `gi|1`, `gi|9`, `gi|5`, `gi|11`, `gi|2`, `gi|8`, `gi|4`, `gi|10`, `gi|6` (in this OID order) and the
  taxids 99990 and 150 (the map), built with

      makeblastdb -in tm_cache_gi.faa -dbtype prot -out tm_cache_gi -title "tm cache gi" \
          -parse_seqids -taxid_map tm_cache_gi.map

  NCBI's mapper packs gi ids, so its prefetch loads the batch's subjects in the order of their gi numbers, where the taxids
  alternate from gi 1 to gi 10: the offset error comes six times (the order of the OIDs would give it once). Rows
  `db.tm.bp.cachegi.*`.

Rows `db.tm.psi.*` and `db.tm.bp.off0.*` (line DBID, item P1): psiblast's prefetch, and `PRE_FETCH_SEQS_LIMIT`, which turns
NCBI's prefetch of a query batch off when it is 0 or below the batch's number of Seq-aligns (one per HSP; the count of
outfmt 7's `# N hits found`): its formatters then fetch the subjects each query's report shows, after the query's own
warnings (where a tabular report has a field that reads them: not `6 qseqid evalue`).

- `tm_psi.faa` and `tm_psi.map`: eight proteins `gb|PTX00001.1|` to `gb|PTX00008.1|` (written by the line's probe script,
  kept outside the repository), six of them copies of `q1` with more and more letters changed, with the taxids
  9606, 10090, 99990, 562, 7777 and none, and two of `q2` (99990, 9606). With `-max_target_seqs 1` psiblast's traceback
  hit list of `q1` still holds `PTX00003.1` (99990), whose lookup writes the offset error of `offset_last_past_end`.
- `tm_psi_q.faa` (`q1`, `q2`, `q3` with weak hits only), `tm_psi_q1.faa` (`q1`), `tm_psi_qw.faa` (`q1` with an `O`
  first, a warning of its own, then `q2` and `q1`), `tm_psi_qb.faa` (`q3`, then `q1`; `psiblast_taxdb_blocks.rs`).
- The databases (version 5, every file makeblastdb writes; `tm_psi_u.faa` is `tm_psi.faa` with the titles
  `subject PTX0000N.1`):

      makeblastdb -in tm_psi.faa -dbtype prot -parse_seqids -taxid_map tm_psi.map -blastdb_version 5 \
          -title "tm psi" -out tm_psi
      makeblastdb -in tm_psi_u.faa -dbtype prot -taxid 9606 -blastdb_version 5 -title "tm psi u" -out tm_psi_u

  `tm_psi_u`, without parsed Seq-ids, is searched in blocks of query batches (`LOSAT_X_BPQBLOCK`) by
  `tests/psiblast_taxdb_blocks.rs`.
- Review 1 of P1: a subject beyond `-max_target_seqs 1` whose Seq-align id NCBI's prefetch finds nowhere, whose lookup
  fails, or that an earlier sequence holds. `tm_psi_rq.faa` (`q1`) is the query; `tm_psi_ra.faa` (`gb|PSX00001.1|`, the
  top hit), `tm_psi_rpig.faa` (`gnl|PIG|124`), `tm_psi_rz.faa` (the top hit and `gnl|BL_ORD_ID|0` with a second
  definition line, which makeblastdb stores with a text tag), `tm_psi_rhx.faa` (`gb|PSX00011.1|`, an unrelated holder of
  `lcl|555` in a second definition line), `tm_psi_rhy.faa` (`lcl|555`, the second hit), with the aliases
  `tm_psi_pig.pal` (`DBLIST tm_psi_ra0 tm_psi_rpig`), `tm_psi_held.pal` (`DBLIST tm_psi_rhx0 tm_psi_ra0 tm_psi_rhy`: the
  holder without a taxid, the subject with 99990) and `tm_psi_held2.pal` (`DBLIST tm_psi_rhx tm_psi_ra0 tm_psi_rhy0`: the
  holder with 99990, the subject without):

      makeblastdb -in tm_psi_ra.faa -dbtype prot -parse_seqids -blastdb_version 5 -title tm_psi_ra0 -out tm_psi_ra0
      makeblastdb -in tm_psi_rpig.faa -dbtype prot -parse_seqids -taxid 99990 -blastdb_version 5 -title tm_psi_rpig \
          -out tm_psi_rpig
      makeblastdb -in tm_psi_rz.faa -dbtype prot -parse_seqids -blastdb_version 5 -title tm_psi_rz -out tm_psi_rz
      makeblastdb -in tm_psi_rhx.faa -dbtype prot -parse_seqids -blastdb_version 5 -title tm_psi_rhx0 -out tm_psi_rhx0
      makeblastdb -in tm_psi_rhy.faa -dbtype prot -parse_seqids -taxid 99990 -blastdb_version 5 -title tm_psi_rhy \
          -out tm_psi_rhy
      makeblastdb -in tm_psi_rhx.faa -dbtype prot -parse_seqids -taxid 99990 -blastdb_version 5 -title tm_psi_rhx \
          -out tm_psi_rhx
      makeblastdb -in tm_psi_rhy.faa -dbtype prot -parse_seqids -blastdb_version 5 -title tm_psi_rhy0 -out tm_psi_rhy0
- Review 2 of P1: `tm_psi_c150`, one volume without parsed Seq-ids holding the taxids 99990, 150, 99990 (`st1` to `st3` of
  `tm_cache.faa`) and 9606 (`tm_psi_ra.faa`), made from four volumes of one taxid each (`makeblastdb -taxid`, version 5,
  `cn1` to `cn4`, not kept) listed in an alias `cnx` with `makeblastdb -in cnx -input_type blastdb -dbtype prot
  -blastdb_version 5 -out tm_psi_c150 -title tm_psi_c150`; `tm_psi_qbig.faa`: 365 copies of the letters of `tm_cache_q.faa`
  (`qp0` to `qp364`, more than one block of 50,000 letters), then `q1` of `tm_psi_rq.faa`.

## Taxid filtering (DI-8a)

Rows `db.tf.*`: `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and `-no_taxid_expansion` with `-db`, five
programs. `taxfilter/` holds the inputs and what `make_taxfilter.py MAKEBLASTDB TAXONOMY4BLAST_SQLITE` builds from them
(makeblastdb 2.17.0, version 5, `-parse_seqids` with a taxid map; `.pjs` and `.njs` dropped):

- `tfp.faa` and `tfp.map`: 23 Swiss-Prot globins of mammals, birds, reptiles, amphibians and fishes (the HBA_HUMAN hits of
  the DI-8 probe); the records of P69905.2, P01923.1, Q9TS35.2 and P63107.2 carry a second definition line (control-A) of a
  bird or another mammal. Database `tfp` (one volume).
- `tfn.fna` and `tfn.map`: 16 globin cDNAs of human (Ensembl GRCh38), zebrafish and fly; human HBA1 (ENST00000320868.9) carries
  a second definition line of a zebrafish transcript, zebrafish hbaa1 (ENSDART00000113098.5) one of a fly transcript.
  Database `tfn` in two volumes (`-max_file_sz 4KB`: OIDs 0-11 and 12-15) under one LMDB.
- `taxonomy4blast.sqlite3`: the 226 rows of NCBI's taxonomy file (the copy distributed with swissprot) for the taxids of both
  maps and all their ancestors, in a new file with the source's schema text. NCBI reads the table `TaxidInfo` and its index
  only, so a taxid's descendants within the databases are those of the whole file. The rows reach it by `BLASTDB=tests/db_search/taxfilter`;
  the `nosql` rows leave it out of reach (NCBI's warning, no restriction).
- `root/taxonomy4blast.sqlite3`: the same file without the rows of 7227 (fly) and 8237 (tuna). The rows `db.tf.bp.root.6` and
  `db.tf.tbn.root.6` freeze NCBI's `-taxids 131567` (the root's only child there) with it, and `tests/db_taxid_root.rs`
  checks that LOSAT's `-taxids 1` writes those bytes (taxid 1 is its own parent: NCBI's expansion of it does not end;
  DI-8 DECISIONS 10).
- `q.faa` (HBA_HUMAN) and `q.fna` (human HBA1 cDNA): the queries. `tl_crlf.txt` (40674 and 8782 with CRLF line ends, an
  empty line and blanks) and `tl_neg.txt` (9443): the taxid lists.
- `tfu_n`, `tfu_p` (rows `db.id.otx.*`, line DBID ordtax) and `tfe_n`, `tfe_p` (rows `db.id.ord.*`, line DBID ordname):
  databases without parsed Seq-ids of four records of `tfn.fna` and `tfp.faa`, taxids 7955, 9606, 7227 (protein 8237)
  and 9606 in OID order, under short titles (`tfu.fna`, `tfu.faa`) or the titles `zebrafish hbaa1 cDNA`, ` `, `  ` and
  the empty one (`tfe.fna`, `tfe.faa`), with the aliases `tfu_nr`, `tfu_pr`, `tfe_nr`, `tfe_pr` of OIDs 1 to 3
  (`make_taxfilter.py --unparsed-only MAKEBLASTDB`; its docstring says how).

The `v4` rows search `x16_taxid` and `x16_ntaxid` (version 4 databases with taxids, above).

## Seqid lists (DI-8b)

Rows `db.sa8.*`: the `ext` `use_this_seqid` of a positive list on the Seq-align formats and what it limits, on
`taxfilter/tfp`, `taxfilter/tfn` and `useseqid/`. `useseqid/` holds `p5m`, `p5mt` (protein, version 5, the records of
`p_m.fa`: a text id and a local id on one OID, a GI that is the Seq-align id, a local id first, a PRF name, PDB ids) and
`n5m` (nucleotide, `n_m.fa`), written by NCBI BLAST+ 2.17.0:
`makeblastdb -in p_m.fa -dbtype prot -out p5m -parse_seqids -blastdb_version 5 -taxid 9606 -title p5m`,
`makeblastdb -in p_m.fa -dbtype prot -out p5mt -parse_seqids -blastdb_version 5 -taxid_map p_m.map -title p5mt`,
`makeblastdb -in n_m.fa -dbtype nucl -out n5m -parse_seqids -blastdb_version 5 -taxid 562 -title n5m`;
`makeblastdb -in n_l.fa -dbtype nucl -out n5l -parse_seqids -blastdb_version 5 -taxid 562 -title n5l` (a local id and a
text id on one OID). `q.faa`, `q.fna`, `qx.fna` are the queries, `qnn.fna`, `qxx.fna` the queries with two copies apart (two
ungapped HSPs per subject), and `sl_m.txt` a seqid list.

Rows `db.sl.*`: `-seqidlist` and `-negative_seqidlist` with `-db`, five programs, on `taxfilter/tfp`, `taxfilter/tfn`
(two volumes), `ids/prot_v5`, `ids/nucl_v5` (version 5) and `ids/prot_v4`, `ids/nucl_v4` (version 4, string ISAM).
`seqidlist/` holds the lists, written by `seqidlist/make_seqidlist.py`:

- text lists in the forms NCBI's reader takes (white-space separated tokens, a leading `>` skipped, `#` to the end of the
  line a comment, CRLF lines): `sl_tfp.txt` (five Swiss-Prot ids of tfp, with and without versions, one in FASTA form, and an
  unknown id), `sl_tfp_neg.txt` (the same ids without versions: a versionless id excludes the versioned stored id),
  `sl_tfn.txt` (four local ids across tfn's two volumes), `sl_prot.txt` and `sl_nucl.txt` (every Seq-id form of
  `ids/prot_v5` and `ids/nucl_v5`: tpg, RefSeq and Swiss-Prot without versions, a pdb chain, local ids, a general id as
  `db:tag`, a versionless accession of an OID of several definition lines, a gi that NCBI drops), `sl_gi.txt` (gi ids only:
  NCBI's `Empty seqid list` warning on a version 5 database; as a positive list no restriction where one volume is reached
  through a path of alias files without a range that cuts it, else no OID: rows `db.gisl.*`; as a negative one no
  restriction; on a version 4 database no OID as a positive list),
  `sl_comment.txt` (no id: no list), `sl_none.txt` (ids no database has: NCBI searches an empty database, LOSAT refuses),
  `sl_empty.txt` (0 bytes: NCBI's `Failed to map seqidlist file ` error);
- binary lists in the format of `blastdb_aliastool -seqid_file_out` (NUL, file size, count, title, date, database length
  and, when it is not 0, the database's date and volume names, then the ids with a one-byte length): `sl_tfp.bsl` and
  `sl_nucl.bsl` (the processed ids of the text lists: no warning), `sl_tfp_badlen.bsl` (a database length that is not
  tfp's: NCBI's `Seqidlist file db info does not match input db` warning, the list is used).

The negative lists on `ids/prot_v5`, `ids/prot_v4`, `ids/nucl_v5` and `ids/nucl_v4` (rows `db.sl.bp.nforms.6`,
`db.sl.bp.v4neg.6`, `db.sl.bp.v4gi.6`, `db.sl.bn.nforms.6`, `db.sl.bn.v4neg.6`) are `losat-rejection` rows: those databases
hold sequences whose every definition line has only a gi or a ti, and a negative seqid list leaves such a sequence without a
line (`CSeqDBNegativeList::FindId` matches no gi or ti): NCBI searches it and drops its hits from the report without notice,
LOSAT refuses when the sequence is to be reported (`the seqid list keeps sequence N of the volume but none of its definition
lines`, DI-8b DECISIONS 17). A positive list that selects such a sequence is refused before the search.

`ids/noacc_v5` (three proteins of `ids/ids_q.faa`, `makeblastdb -dbtype prot -blastdb_version 5` without `-parse_seqids`,
NCBI 2.17.0) has no accession table: a seqid list on it is NCBI's `BLAST Database error: DB contains no accession info.`
(row `db.sl.bp.noacc.6`, rc 2).

A list file that is not on the search path (NCBI: `Error: NCBI C++ Exception: ... To be memory mapped the file must
exist: ''`, exit 255, with the path of NCBI's build tree) is checked by `src/cli_v2_tests.rs`, not by a row.

Rows `db.gisl.*` (line DBID, item gisl): `-seqidlist seqidlist/sl_gi.txt` (a positive list) on version 5 databases of one
and of two volumes, and on the alias files of `aliaskeys/` over one volume. NCBI's `x_ApplyUserGiList` clears every OID of
the list that lost every id, and the search finds an empty database (each query's `Filtering resulted in an empty
database.`), unless its OID bit set is the special all-set one, which a copy of one volume's special set makes (the `-db`
names reach one volume, and at least one path of alias files to it crosses no range that cuts it): the clearing writes into
a stale array, and the whole database is searched. So `taxfilter/tfp`, `taxfilter/tfn.00` (one volume of `tfn`),
`aliaskeys/p_r123` (the range of the whole volume), `p_twice`, `p_dupvol`, `p_chain` (an alias of `p_dupvol`), `p_above`
(the whole range above `p_dupvol`), `p_cover` and `p_coverr` (`p_r123` and `p_range` in both orders), `p_arrcov` (`p_union`,
a bit array, then `p_r123`) are searched whole; `taxfilter/tfn` (two volumes), the two databases `taxfilter/tfp
ids/prot_v5`, `aliaskeys/p_range`, `p_union`, `p_halves` (`p_lo` and `p_hi`, ranges 1-11 and 12-23 whose union is the
volume), `p_abovecut` (a range that cuts the volume above `p_dupvol`) and `n_vol1` are searched empty. The negative list on
`taxfilter/tfn` excludes nothing. `sl_gifull.txt` (`gi|1111|ref|NP_1.1|`) adds NCBI's two warnings `Ignoring extra parts
(synonyms?) in FASTA-style ID`.

Rows `db.ordsl.*` (line DBID, item ordsl): seqid lists that name `gnl|BL_ORD_ID|n` (`multidb/sl_ord*.txt`,
`sl_ord1.bsl`, `sl_pxm*.txt`, `sl_nxm7.txt`) where a line with such an id lies in a volume whose first OID is not 0: a
later database of the list (`taxfilter/tfu_p` behind `multidb/pa`, 12 sequences; `tfu_n` behind `multidb/na`, 8;
`multidb/pxm` and `nxm` likewise) or the second volume of `multidb/ov_p` and `ov_n`. The accession tables of
`makeblastdb -input_type blastdb` (`tfu_*`, `ov_*`) and of `-parse_seqids` FASTA ids (`pxm`, `nxm`) hold the stored
number, which the lookup reads, while NCBI's filtered header read tests the line as `gnl|BL_ORD_ID|<stored number plus the
volume's first OID>` (`CSeqDBVol::x_GetHdrAsn1` with `adjust_oids`). So a positive list can select a sequence whose every
line it then removes, and a negative list can remove a line of a sequence it did not exclude; NCBI leaves such a sequence's
hits out of the report without notice (a positive list's at the Seq-aligns, line DBID item lineless; a negative list's at
the fetches of the traceback, the redo and the ungapped re-evaluation, item lineless2), and so does LOSAT (`bp.list.6` on
`pa tfu_p`, `bp.listb.6` with the binary list, `bp.listneg.6` and `bn.listneg.6` negative, `bp.vol.6` the second volume of
`ov_p`, one row in each other program; LOSAT refused them before). Where a sequence keeps another line, the
test decides which lines NCBI writes (`bp.two.*`, `bn.two.*`: OID 12 of `pa pxm` and OID 8 of `na nxm` hold
`ref|NP_900003.1|` / `ref|NM_900003.1|` and `gnl|BL_ORD_ID|7`, read as 19 / 15; positive and negative lists of 7 and of
the sum, outfmt 6 with `sallseqid salltitles` and 16). `bp.tag.6` and `bp.tagq.6`: the line `gnl|BL_ORD_ID|0` of OID 13 of
`pa pxm`, which makeblastdb stores with the string tag "0", selected by `NP_900005.1`; NCBI reads the header only for a
sequence it reports and then throws (`GetId`, exit 255): a query that misses it (`-evalue 1e-5`) is NCBI's run, one that
hits it (`multidb/xq.faa`) LOSAT's refusal. The other `ncbi` rows: a list whose number no key holds (`bp.none.6`, an
empty database), and a sequence of the first volume, which keeps its line (`bp.vol0.0`, `bn.vol0.6`, `tx.vol0.6`).

`ll2/`: rows `db.ll2.*` (line DBID, item lineless2). `ll2/make_ll2.py` (its docstring has the usage and the records)
writes `ll2/prot_v4`, `prot_v5`, `nucl_v4` and `nucl_v5` with `makeblastdb -parse_seqids`: 14 subjects made from one query
with 0, 2, 4, ... substitutions (distinct ranks, one equal pair), of which the only line of r0 (the query itself), r2a
(the first of the pair) and r6 has a gi (version 4) or a ti (version 5). A `-negative_seqidlist` (`ll2/neg.txt`: r12's
accession, which the OID mask removes, and an id no database holds) keeps those three sequences and removes their lines;
NCBI's fetches with `check_oid_exclusion` exclude them, so they count in the statistics and take no place in the final
hit lists (`.m1.`: r1 is the best row; `db.ll2.tx.m10.6`, `db.ll2.bn.ug.m10.6` and `db.ll2.bn.ug.m1.6`: an ungapped
search of a nucleotide database never holds them in its preliminary hit list), except a protein subject of an ungapped
search (`db.ll2.bx.ug.*`), which takes its place before the cut and leaves at the Seq-aligns. `ll2/neggi.txt` (a gi) is a non-empty id set on version 4 and an empty list on version 5. `ll2/prot_v4bad` and
`nucl_v4bad` (`make_ll2.py --damaged-only`) damage the header of OID 8 (r7, which keeps its line): under `neg.txt`
NCBI's fetch reads it although r7 is not reported at `-max_target_seqs 3`, and the search fails (exit 255 in the
traceback, 3 in the ungapped re-evaluation); LOSAT refuses (the `losat-rejection` rows `db.ll2.*.bad.6`).
Queries: `ll2/q.faa`, `q.fna` (the query and its codons), `q2.faa`, `q2.fna` (the same and its first part).

`x4b/`: rows `db.x4b.dbfail.*`, `db.x4b.excl.*`, `db.x4b.cut0.db.*`, `db.x4b.sp.db.*`, `db.x4b.sp.two.db.*` and
`db.x4b.pssm.none.db.*` (line RG, item X4b-1 fix rounds 1 and 2: tblastn `-use_sw_tback`). `x4b/make_x4b.py` (its
docstring has the usage and the records) writes with NCBI BLAST+ 2.17.0 makeblastdb (version 5) `x4b/fail`
(`-parse_seqids`: X4's sok1 in frame +1 and sfail in frame -1, on which tks_W's Smith-Waterman redo fails), `x4b/excl`
(the same and, at OID 2, a weaker copy of sfail's first half whose only line has a ti) and `x4b/polycg`
(`tests/fasta/x4b/polycg.fna`, no hit for a PSSM). `x4b/neg.txt` lists an id that no database holds: as a
`-negative_seqidlist` it keeps OID 2 of `excl` without a line, so NCBI's redo excludes that subject at its fetch and
sets the batch's status back to 0 (`blast_kappa.c:3560-3575`), which clears the failure of the match on sfail before it
(exit 0 where `fail` ends with `Traceback failed`).

`x4excl/`: rows `db.x4excl.*` (line RG item X4-excl, BLASTP `-use_sw_tback`). `x4excl/make_x4excl.py` (its docstring has
the usage and the records) writes the version 5 databases `x4excl/fwd` (sexact, sfail, sok1 and sx_frag, whose only line
has a ti) and `x4excl/rev` (sx_strong, ti only, then sfail) with `makeblastdb -parse_seqids` from
`tests/fasta/rps/tks_w.faa` and `tests/fasta/x4/sw_subjects.faa`, and the lists `neg.txt` (an id no database holds) and
`neg_exact.txt` (also sexact's accession, which the OID mask removes). The negative list keeps the ti-only subjects
without a line, and NCBI's fetch of the Smith-Waterman redo excludes them after the early termination test of their
match, which sets the redo's status back to 0: a failed match before such a match does not fail the batch. On version 4
the list masks a gi-only sequence, which is then never searched, hence version 5 (as `ll2/prot_v5`).

## Alias files with OIDLIST (line RG, series D, D1b-1)

`oidlist/` (`make_oidlist.py`): alias files with `OIDLIST` over taxfilter's `tfp` (23 proteins) and `tfn` (two volumes of
12 and 4 nucleotides), and their OID list files (the big-endian index of the last OID, then one bit per OID, the most
significant bit first, set for a kept OID): `k.msk` keeps 12 OIDs of `tfp`, `n.msk` six OIDs of each volume of `tfn`
(read relative to each volume), `h1.msk` and `h2.msk` the bits past the header's last OID and past the file's end,
`wrap.msk` the last OID 0xFFFFFFFF (no OID), `hdr.msk` the header alone and `empty.msk` no byte. `p_k`, `p_kn` (with
`NSEQ`), `p_kr` (with `FIRST_OID` and `LAST_OID`), `p_ku` (`p_k` and `aliaskeys/p_range` over one volume) and `n_n` make
NCBI's constructor scan the kept OIDs; `p_knl`, `n_nnl` (with `NSEQ` and `LENGTH`) and `p_kpn` (`NSEQ` above `p_knl`)
do not. `p_miss`, `p_missl`, `n_miss` and `n_missl` name the absolute path of a file that does not exist (NCBI's error
prints it, the same on every machine), with and without `NSEQ` and `LENGTH`, and `p_dir` the directory `/`. The rows are
`db.ol.*` (every search program; `rps/olrps.pal` and `olwrps.pal` over `rps_mini` with `olall.msk`, `olno2.msk` or
`olonly2.msk` as `rps_mini.msk`); `segmasker` and `dustmasker` read them in their rows `db_oidlist_*`.

`x4b2/`: rows `db.x4b2.fail.w.db.*`, `db.x4b2.dbsize1.w.db.*`, `db.x4b2.rej.crash.db.*`, `db.x4b2.ctx.*.db.*`,
`db.x4b2.rej.ungapped.ac.b1.db.*`, `db.x4b2.excl.*` and `db.x4b2.past.excl.nolist.*` (line RG, item X4b-2: blastx
`-use_sw_tback`).
`x4b2/make_x4b2.py` (its docstring has the usage and the records) writes the queries `tests/fasta/x4b2/*.fna` (X4's
crafted proteins back-translated with one codon per letter, and cuts of `bslice_q.fna`; the `*_pad.fna` files are
query files whose batches are followed by queries without a hit, NCBI's runs of which are the expectations of the rows
`db.x4b2.past.*` and `db.x4b2.bh.false.6`, approved exception 23 of PD-LOSAT-NCBI-DEFECTS; `bh3_pad.fna` pads dev's
`tests/fasta/outfmt0/bh3_q.fna`), the subjects
`tests/fasta/x4b2/ok_fail.faa`, and with NCBI BLAST+ 2.17.0 makeblastdb (version 5, `-parse_seqids`) the protein
databases `x4b2/fail` (sok1 and sfail, on which tks_W's Smith-Waterman redo fails) and `x4b2/excl` (the same and, at
OID 2, a weaker copy of sfail's first half whose only line has a ti). `x4b2/neg.txt` lists an id that no database
holds: as a `-negative_seqidlist` it keeps OID 2 of `excl` without a line, so NCBI's redo excludes that subject at its
fetch and sets the batch's status back to 0, which clears the failure of the match on sfail before it. The databases
carry their creation date: the committed files, not a new run of the script, are the inputs of the frozen outputs.

## Alias files with OID_MASK_TYPE (line RG, series D, D1b-2)

`oidmask/` (`make_oidmask.sh`): `om.faa`, eight proteins whose definition lines carry predicted RefSeq ids (`XP_`) alone,
with a GenBank id, beside a normal RefSeq line (record 4, two lines) or none, and `omn.fna`, their back-translation with
the nucleotide ids (`XM_`, `XR_`, `NM_`, `NR_`). `omp` and `omp4` (versions 5 and 4) and `omn` (two volumes of four) are
`makeblastdb -parse_seqids -oid_masks 1` databases, `ompx`, `omp4x` and `omnx` the `blastdb_aliastool -oid_masks 1` aliases
over them (`OIDLIST` of the `.pxm`/`.nxm` mask and `OID_MASK_TYPE 1`). `tall.pal` and `ntall.nal` keep every OID
(`ones.msk`, `onesn.msk`) with the type 1, so that OIDs 1, 2 and 5 are left without a line; `t0` to `tfirst` write the
key's other values (0, 2, 3, 4294967295, +1, abc, 4294967296, -1, empty, the key alone, the key without `OIDLIST`, a bad
`FIRST_OID` before it, no `NSEQ` and `LENGTH`), `j1`, `j2`, `ompy` and `j5` its nesting (a range, a parent list without
the key, two siblings over one volume). `neg.txt` and `pos.txt` are the seqid lists of rows `db.om.bp.tall.neg.6` and
`.pos.6`. `rps/xp101.{phr,pin,pdb,pos}.gz` are the files of `rps_mini` built with profile 101 under the Seq-id
`ref|XP_000101.1|` (`rps/build.sh`; its other files are the same as `rps_mini`'s), searched through `rps/omrps.pal` and
`omwrps.pal` (`OID_MASK_TYPE 1`, the second with `NSEQ` and `LENGTH`) with `olall.msk` as `rps_mini.msk`.
`rps/xp3v01.{phr,pin,pdb,pos}.gz` are the files of volume `rps_mini3v.01` built with profile 102 under
`ref|XP_000102.1|`, searched through `rps/om3vrps.pal` (`DBLIST rps_mini3v`, `OIDLIST` of `olall.msk` as
`rps_mini3v.msk`, `OID_MASK_TYPE 1`): NCBI searches the three volumes bare, so the profile reaches the reports
(`db.om.rps.rej.x3v`, `db.om.rtbn.rej.x3v`) unless it is past the cut (`oidmask/qrps8.faa`, eight queries of
`tests/fasta/rps/q40.faa` whose first subject is another profile, rows `db.om.rps.x3vc.*` with `-max_target_seqs 1`).
`omb` (`make_oidmask.sh`) holds the byte 0x80 in the titles of two predicted lines, patched into its `.phr` after a
build with `#` there, behind `ombx.pal` (`OID_MASK_TYPE 1`): row `db.om.bp.rej.omb.6` and A271. `ompart.pal` (`nol.msk`:
OIDs 0, 3, 4, 6 and 7, each left with a line under the type) is the input of makeblastdb case `B_al_om`; `ompartn.nal`
(`noln.msk`: local OIDs 0 and 3 of each volume of `omn`, global 0, 3, 4 and 7), its nucleotide twin, that of `B_al_omn`;
`t0.pal`, `t2.pal`, `j2.pal` and `j5.pal` are copied in cases `B_al_t0`, `B_al_t2`, `B_al_j2` and `B_al_j5`. The rows are
`db.om.*` (every search program); `segmasker`, `dustmasker`, `windowmasker`, `makembindex` and `blastdbcmd` read the
databases in their rows `db_oidmask_*`, `oidmask_*` and A242 on.

## Alias files with GILIST, TILIST and SEQIDLIST (line RG, series D, D1c-1)

`idlist/` (`make_idlist.py`, whose docstring lists the sources and databases): `g` (version 5, 14 proteins, some with two
or three definition lines), `g4` (version 4), `gm` (two volumes), `gu` (no parsed Seq-ids, no ISAM), `gt` (`g4` with a TI
ISAM written by the script; makeblastdb writes none), `n` and `n4` (nucleotides), and the GI, TI, seqid and taxid list
files and alias files over them: the list alone, under `NSEQ` and `LENGTH` (read at the search's first use instead of
the constructor's scan), ANDed with a range or an `OIDLIST`, nested, two alias files over one volume, beside the bare
volume (sequences left with no line), over one or both volumes of `gm`, and the list files' forms and errors (binary,
no last separator, a GI past 2^32, `#`, empty, a bad count, missing at the absolute path
`/losat-idlist-fixture-missing/none.txt`, two files in one value). The rows are `db.il.*` (every search program;
`rps/ilmiss.pal` and `ilmissw.pal` over `rps_mini` with a missing GILIST) and blastdbcmd's cases A272-A366
(`tests/blastdbcmd_alias_keys.rs`). `u_np101.txt`, `u_p10102.txt` and `u_np103.txt` are the searches' own
`-seqidlist` files beside the alias files' lists.

## Ungapped nucleotide databases (line NP, O4a-1)

`ungapped_abr` (the 33 genes of `LOSAT/tests/fixtures/blastn_regression/inputs/culling_abricate_genes.fa`, as `ungapped_abr.fna`) is built like
the mlst databases (parsed Seq-ids and a hash index), `ungapped_ambrun` (the subject `ambrun_s.fa`, a stretch with 30 `N` and one `R` inside a match)
without, `ungapped_re` (the six subjects of `re_s.fa`: N every other base, a run of N, IUPAC codes, a minus-strand copy) without; rows
`db.blastn.ungapped*` (`blastn -ungapped`). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/ungapped_abr.fna -dbtype nucl -out tests/db_search/ungapped_abr \
        -parse_seqids -hash_index -title PubMLST -blastdb_version 4
    makeblastdb -in tests/db_search/ungapped_ambrun.fna -dbtype nucl -out tests/db_search/ungapped_ambrun \
        -title "ambiguity run" -blastdb_version 4
    makeblastdb -in tests/db_search/ungapped_re.fna -dbtype nucl -out tests/db_search/ungapped_re \
        -title "ungapped re-evaluation" -blastdb_version 4
    rm tests/db_search/ungapped_abr.njs tests/db_search/ungapped_ambrun.njs tests/db_search/ungapped_re.njs

`mejo` is `LOSAT/tests/fasta/MejoMJNV.fasta` (one sequence, 401,182 bases) as pyani ANIb builds its databases
(`makeblastdb -dbtype nucl -in FILE -title STEM -out OUT`, no parsed Seq-ids; version 4 for the repository), the
subject of rows `db.xdrop.*` (DI-6: the X-drops, `-window_size`, `-off_diagonal_range` and `-no_greedy`) and of some of
rows `db.minraw.*` (DI-6b: `-min_raw_gapped_score`; the queries
`tests/fasta/outfmt0/anib_frags*.fna` are 1020-nt fragments of `MelaMJNV.fasta`, pyani's fragment size). Run from `LOSAT/`:

    losat makeblastdb -dbtype nucl -in tests/fasta/MejoMJNV.fasta -title MejoMJNV -out tests/db_search/mejo -blastdb_version 4
    rm tests/db_search/mejo.njs

The rows `db.blastn.ungapped.{mb,bn}.{8,9,11,12}` (O4b) are NCBI's Dense-diag Seq-aligns of the same search; the rows `db.blastn.ungapped.mb.17` and `db.blastn.ungapped.bn.17b` (O4a-2) are its SAM (outfmt 17: one record per subject, its HSPs joined).

## Several databases in -db (DI-8f)

Rows `db.mdb.*`: `-db "a b"` searches the databases `a` and `b` as one virtual database, as NCBI's `CSeqDB` does (the `db` cell
of a row is one argument, so a cell with spaces is a list of names). `multidb/` holds the databases and what
`make_multidb.py MAKEBLASTDB` builds them from (makeblastdb 2.17.0; the records and taxids are those of `taxfilter/tfp.faa`,
`tfp.map`, `tfn.fna` and `tfn.map`, definition lines kept; `.pjs` and `.njs` dropped; the whole directory is about 0.8 MB):

- `pa` / `pb`: version 5 protein, `-parse_seqids`, taxid map; records 1-12 of `tfp.faa`, and records 13-23 plus `Q00001.1`, an
  exact copy of record 2's sequence (taxid 9606), so that two databases hold an identical sequence. `p4a` / `p4b`: the same
  records as version 4 without `-parse_seqids`.
- `na` / `nb`: version 5 nucleotide, `-parse_seqids`, taxid map; records 1-8 of `tfn.fna`, and records 9-16 plus
  `ENSTTEST0001.1`, a copy of record 1's sequence (taxid 9606). `n4a` / `n4b`: version 4 without `-parse_seqids`; `nnp`:
  the records of `nb` as version 5 without `-parse_seqids`.
- `d1/g` and `d2/g`: the records of `na` and `nb` under one base name in two directories. `et`: the records of `na` with an
  empty title. `al.nal`: an alias file (`TITLE multidb alias of na`, `DBLIST na`). `mix_u.nal` (line DBID mixdb, written by hand): an
  alias file of `na nnp` whose `FIRST_OID 9`, `LAST_OID 17` keep the sequences of `nnp` (without the accession table) alone. `sl.txt`: a seqid list of the first
  accession of `na`, `nb`, `pa` and `pb`.
- `ov_p` / `ov_n` (line DBID ordsl; `make_multidb.py --ordsl-only MAKEBLASTDB` builds only these): `-input_type blastdb` of
  `../taxfilter/tfu_p` and `tfu_n` in two volumes of two records each (`-max_file_sz` 400 and 500 bytes) under one LMDB file;
  `sl_ord0.txt`, `sl_ord1.txt`, `sl_ord2.txt`, `sl_ord9.txt`, `sl_ord13.txt`: the list `gnl|BL_ORD_ID|n`, and `sl_ord1.bsl`
  its binary form (rows `db.ordsl.*`).
- `pxm` / `nxm` (line DBID ordsl): version 5 `-parse_seqids` of two records: record 13 of `tfp.faa` (record 9 of `tfn.fna`)
  under the lines `ref|NP_900003.1|` (`NM_`) and `gnl|BL_ORD_ID|7`, and a random sequence under
  `ref|NP_900005.1|gnl|BL_ORD_ID|0` (stored with the string tag "0"); `xq.faa`, that random protein; the lists
  `sl_pxm7.txt`, `sl_pxm19.txt`, `sl_nxm7.txt`, `sl_ord7.txt`, `sl_ord15.txt`, `sl_ord19.txt`, `sl_np5.txt`.

The ncbi rows search `taxfilter/q.fna` (BLASTN, BLASTX, TBLASTX) or `taxfilter/q.faa` (BLASTP, TBLASTN, PSIBLAST) and cover:
two names in outfmt 0, 5, 6, 7, 8, 11, 15 and 16 (`bn`), the other order (`bn.rev`), a name twice (`dup`), Galaxy's call with
each name in double quotes (`galaxy`), extra spaces around the names (`space`), an alias file (`alias`), one base name in two
directories (`samename`), two version 4 databases (`v4`, `bp.v4`), a version 4 and a version 5 database (`v4v5`: NCBI exits 2,
`DB list contains both Version 4 and Version 5 dbs`), an empty title (`title`), `-taxids` and `-seqidlist` over both names
(`tx`, `sl`), the five other programs (`bp`, `bx`, `tbn`, `tbx`; `ps` searches the version 4 pair, because LOSAT's PSIBLAST
refuses databases with parsed Seq-ids). The three `losat-rejection` rows are LOSAT's refusals of a version 5 pair with and
without the accession table (`bn.mixed.6`: NCBI stops at the first hit of a parsed sequence), of an alias that lists the
volume of another name (`bn.aliasdup.0`), and of a negative seqid list (`multidb/sl_negpartial.txt`, FBtr0083245 of the
second volume) on `taxfilter/tfn.00`, the first of the two volumes of one LMDB file (`bn.negpartial.6`: NCBI clears bit -1 of
its OID mask and crashes).

## The environment around BLASTX's search and report (line DBID, N2d)

Rows `db.bxenv.*` search `bslice` with `bslice_q.fna` (BLASTX) with an environment variable in `extra_args`. NCBI converts
`PRE_FETCH_SEQS_LIMIT` with `NStr::StringToInt` before the report of each query batch (`blast_app_util.cpp:732-737`) and
stops with a `CStringException` (exit 255; with `-mt_mode 1` its nodes catch it, so NCBI writes nothing and exits 0) on a
value it cannot convert, an empty one too; the archive of outfmt 11 is written without that call (`blastx_app.cpp:283-288`).
`BL2SEQ_LEGACY` makes the subjects of `-subject` a search of each subject on its own (`blast_app_util.cpp:206-210`); a `-db`
search ignores it. LOSAT stops where NCBI's first batch converts the value, after the batch's search and before its report
(NCBI's exception, exit 255, after the outfmt 0 prolog; env parity H), and refuses `BL2SEQ_LEGACY` with `-subject`
before the search (a search whose report would not change, no hits or one subject record in outfmt 5, 6 or 10, is
refused too). The rows are `abc` in outfmt 6, 16 and 11, an empty value in outfmt 0, the integer 5 in outfmt 6, and
`BL2SEQ_LEGACY=1` in outfmt 0 and 6, all of kind ncbi. The `-subject` cases are
`environment_that_ncbi_reads_around_the_report_stops_where_ncbi_reads_it` in `tests/run_local_blastx.rs`.

## The same environment in the other programs (line DBID, N2e)

Rows `db.envshared.*` search `arch_sprot` (BLASTP, PSIBLAST), `arch_tbn_snucl` (TBLASTN), `arch_tbx_snucl` (TBLASTX) and
`arch_snucl` (BLASTN) with the queries of `tests/fasta/outfmt11/` and the environment of the section above. BLASTP, TBLASTN,
TBLASTX and BLASTN write the archive of outfmt 11 without the prefetch too (`blastp_app.cpp:280-285`,
`tblastn_app.cpp:294-299`, `tblastx_app.cpp:197-202`, `blastn_app.cpp:304-309`), so a value that NCBI cannot convert leaves it
as without the variable; psiblast's iterations fetch ahead before they write an archive (`psiblast_app.cpp:262-266`), which
LOSAT's PSIBLAST does not write. A `-db` search of each of them ignores `BL2SEQ_LEGACY`. NCBI's automatic split of a `-db`
search by queries (`-num_threads` 2 or more, a large query file, a small database) catches the exception as `-mt_mode 1` does. The rows are `abc` (BLASTN: an
empty value) in outfmt 11, `BL2SEQ_LEGACY=1` in outfmt 6 (TBLASTN: 0) and BLASTN's `abc` in outfmt 6 (NCBI's exception, exit
255), all of kind ncbi. The `-subject` cases of the five programs are `tests/search_environment.rs`.

## Values that NCBI cannot convert (env parity H)

Rows `db.envconv.*` search the databases of the two sections above (`bslice` for BLASTX, `rps_mini` with `-comp_based_stats 0`
for RPSBLAST) with a `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE`, `BATCH_SIZE` or `PRE_FETCH_SEQS_LIMIT` that `NStr::StringToInt`
(flags 0: a sign and digits, no spaces, within an `int`) cannot convert. NCBI stops with its `CStringException`
(`CATCH_ALL`: `Error: ` and its report, which names `ncbistr.cpp` in the build tree of the official Linux x64 release, exit
255) where it reads the value: `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` in the query splitter of the first batch, after the
outfmt 0 prolog (BLASTN reads `CHUNK_SIZE` first for its batch size mixer, after the prolog too); `BATCH_SIZE` before
anything (BLASTN after the prolog; PSIBLAST does not read it); `PRE_FETCH_SEQS_LIMIT` before the first batch's report,
except where every format is outfmt 11 (but in PSIBLAST). With outfmt 11 NCBI writes its error archive, 523 bytes for
every input. The probe of every program, mode and value is `env-parity-20261108/h_probe/REPORT.md` (outside the tree).
The `-subject` cases are `tests/search_environment.rs` and `tests/run_local_blastx.rs`. NCBI reads the titles of a
batch's reported sequences in its traceback, before the prefetch: with `nophr`, `bnophr` or `nnophr` (no header file)
and a query with hits it stops with the header file's error (exit 2) whatever `PRE_FETCH_SEQS_LIMIT` holds; the message
names the file's absolute path, so these cases are `a_header_file_that_ncbi_cannot_read_comes_before_the_prefetch` in
`tests/search_environment.rs`, not rows.

## TBLASTN `-query_loc` and the chunk masks (line NP, series SM)

Rows `db.tqloc.*` search `tqloc` with `tqloc_q.faa`, `-query_loc 5-3004 -lcase_masking` and
`CHUNK_SIZE=1000 OVERLAP_CHUNK_SIZE=100`. NCBI restricts the lower-case runs of the whole record to each chunk
(`split_query_cxx.cpp:278-281`), one letter longer, so a run that ends on the letter before the interval masks the
interval's first letter in the first chunk's preliminary search. The inputs come from `random.Random(20261104)`:

- `tqloc_q.faa`: `a` (`AC`, which the range skips), `p` (`MKVA` and 200 random residues) and `q` (`MKVa` and a random
  protein Q of 3,000 residues that starts with W). `q` is the batch's second query and its part starts the first chunk.
- `tqloc.fna`: `s14` and `s30`, the first 14 and 30 residues of Q between 50 random residues on each side,
  back-translated with one codon per residue. At `-evalue 1e-7` NCBI reports `q` against `s30` only; with `MKVA` in
  place of `MKVa` it reports `s14` too.

The database (version 4, no `-parse_seqids`; only `.nin .nsq .nhr` are kept). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/tqloc.fna -dbtype nucl -out tests/db_search/tqloc \
        -title "Q[0:14] and Q[0:30] of tqloc_q.faa between random letters, back-translated" -blastdb_version 4
    rm tests/db_search/tqloc.njs

The `-subject` rows `tqloc.tblastn.*` of `outfmt0_manifest.tsv` use `tests/fasta/outfmt0/tqloc_q.faa` (`MKVa`, the same
Q and 37,000 more random residues: 40,000 searched letters, which NCBI splits at the default sizes) and
`tests/fasta/outfmt0/tqloc_s.fna` (a copy of `tqloc.fna`), at `-evalue 1e-6`.

## BLASTX `-query_loc` and the chunk masks (line NP, series SM)

Rows `db.bxqloc.*` search `bxqloc` with `bxqloc_q.fna`, `-query_loc 10-9009 -lcase_masking` and
`CHUNK_SIZE=1002 OVERLAP_CHUNK_SIZE=99` (a BLASTX chunk is a multiple of 3). NCBI restricts the lower-case runs of the
whole record to each chunk (`split_query_cxx.cpp:278-281`), one letter longer, so a run that ends on the letter before
the interval masks the interval's first letter in the first chunk's preliminary search. NCBI splits a BLASTX query
batch only when it holds one query. The inputs come from `random.Random(20261104)` (Q and its back-translation, as in
the `db.tqloc.*` rows) and `random.Random(20261105)`:

- `bxqloc_q.fna`: `a` (`AC`, which the range skips), `p` (`ATGAAAGTT` and 600 random nucleotides) and `q`
  (`ATGAAAGTt` and Q of 3,000 residues, which starts with W, back-translated with one codon per residue: 9,000
  nucleotides). Without `BATCH_SIZE`, `p` and `q` share a query batch, which is not split (`db.bxqloc.two.6`: `s14` is
  reported); with `BATCH_SIZE=100` `q` is alone in its batch and split (`db.bxqloc.b100.6`: at `-evalue 1e-7` NCBI
  reports `q` against `s30` only; with `ATGAAAGTT` in place of `ATGAAAGTt` it reports `s14` too).
- `bxqloc.faa`: `s14` and `s30`, the first 14 and 30 residues of Q between 50 random residues on each side.

The database (version 4, no `-parse_seqids`; only `.pin .psq .phr` are kept). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/bxqloc.faa -dbtype prot -out tests/db_search/bxqloc \
        -title "Q[0:14] and Q[0:30] of bxqloc_q.fna's q between random letters" -blastdb_version 4
    rm tests/db_search/bxqloc.pjs

The `-subject` rows `bxqloc.blastx.*` of `outfmt0_manifest.tsv` use `tests/fasta/outfmt0/bxqloc_q.fna` (`ATGAAAGTt`,
the same 9,000 nucleotides and 11,000 more random ones: 20,000 searched nucleotides, which NCBI splits at the default
sizes, 10,002 and 297) and `tests/fasta/outfmt0/bxqloc_s.faa` (a copy of `bxqloc.faa`), at `-evalue 2e-7`.

## BLASTN `-query_loc` and the chunk masks (line NP, series SM)

Rows `db.bnqloc.*` search `bnqloc` with `bnqloc_q.fna`, `-query_loc 10-3009 -lcase_masking` and
`CHUNK_SIZE=1001 OVERLAP_CHUNK_SIZE=100` (NCBI BLASTN stops at 1000/100 with `Empty CBlastQueryVector`: its batch size
mixer's maximum is the chunk size less 1,000). NCBI restricts the lower-case runs of the whole record to each chunk
(`split_query_cxx.cpp:278-281`), one letter longer, so a run that holds the letter before the interval masks the
interval's first letter in the preliminary search of each chunk that holds it, when the query's part there is longer than
the interval's start (here 9 letters). Without `BATCH_SIZE` each record is a batch of its own; with `BATCH_SIZE=100000`
`p` and `q` share a split batch (four chunks of 1,126 letters) and `q`'s first letter lies 174 letters into the second
chunk. The inputs come from `random.Random(20261105)`:

- `bnqloc_q.fna`: `a` (`AC`, which the range skips), `p` (1,209 random nucleotides) and `q` (`ACGTACGTa` and B, 3,000
  random nucleotides).
- `bnqloc.fna`: `s11`, `s24` and `s40`, the first 11, 24 and 40 nucleotides of B between 50 random nucleotides on each
  side (the nucleotide after each differs from B's).

At `-task blastn -soft_masking false -evalue 1e-7` (`db.bnqloc.hard.*`, `db.bnqloc.two.6`) NCBI reports `q` against
`s40` only; with `ACGTACGTA` in place of `ACGTACGTa` it reports `s24` too. With the default `-soft_masking true` at
`-evalue 10` (`db.bnqloc.soft.6`) the masked letter leaves the lookup table, and `s11`, whose only seed starts there, is
not reported (with `ACGTACGTA` it is). Megablast at `-evalue 5e-17` (`db.bnqloc.mb.7`) reports no hit (with
`ACGTACGTA`, `s40`).

The database (version 4, no `-parse_seqids`; only `.nin .nsq .nhr` are kept). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/bnqloc.fna -dbtype nucl -out tests/db_search/bnqloc \
        -title "B[0:11], B[0:24] and B[0:40] of bnqloc_q.fna's q between random letters" -blastdb_version 4
    rm tests/db_search/bnqloc.njs

## `nident` alone and HSPs without identities (line DBID, L1)

Rows `db.nz.*`: with `PRE_FETCH_SEQS_LIMIT` turning NCBI's prefetch off, a tabular report whose only field of the
alignment block is `nident` (`6 qseqid nident`; in the translated programs also with `sframe` or `frames`) reads a
subject only on a row whose HSP has identities (the Seq-align's `num_ident` above 0, tabular.cpp:904), so a subject whose
shown HSPs have none is not looked up in the taxonomy database. The `pident` rows are the controls (every row reads the
subject). `-threshold 9` lets the seeds of alignments without identities through (a word of three such pairs scores at
most 9 in BLOSUM62).

- `nz_p.faa` (protein) and `nz_n.fna` (its proteins as codons): `sZ` (`M`, then `VY` 30 times, `CC` and `VY` 10 times),
  `sI` (`M` and B, 60 random residues), `sM` (`M`, `VY` 30 times, `W` 20 times and B) and `sN` (`M` and `VY` 20 times),
  with the taxids 9606, 562, 10090 and 99990 (`nz_p.map`, `nz_n.map`): only `sN` has the taxid whose lookup writes the
  offset error of `offset_last_past_end`. Queries `nz_q.faa` / `nz_q.fna`: `q1` (`IF` 30 times: aligned to the `VY`
  parts without identities), `q2` (`q1`, `D` 20 times and B), `q3` (B), `q4` (B, `D` 20 times and `q1`); `nz_q1.*` and
  `nz_q2.*` hold `q1` and `q2` alone. In BLASTX, TBLASTN and TBLASTX other frames align with identities too, hence
  `-max_hsps 1` and lower e-values in some rows.
- `nz_psi.faa`: psiblast's composition-based statistics drop the `IF`/`VY` alignments, so `sP` is `M` and `qP` (120
  residues of `IFKELSD`) with I/V, F/Y, K/R, E/Q, L/M, S/T and D/N exchanged, `sQ` is `M` and `qP`, `sR` the first half
  of `sP` and the second of `sQ`, with the taxids 99990, 9606 and 10090 (`nz_psi.map`); `nz_psi_q.faa` holds `qP`, its
  first half `qPh` and its second half `qP2`.
- `nz_pc.faa` / `nz_nc.fna` (`db.nz.*.c150.*`) are `nz_p` / `nz_n` with `sZ`'s taxid 150 (`nz_pc.map`, `nz_nc.map`),
  which NCBI's taxonomy cache (256 slots, taxid & 255) puts in the slot of 99990, the damaged last entry of
  `offset_last_past_end`: LOSAT refuses a run that may look both up (`TaxdbRun::cache_collision`), but not one whose
  rows look nothing up (`q1` alone with `nident`, a positive `PRE_FETCH_SEQS_LIMIT` below the batch's Seq-aligns).
  BLASTX writes each batch before it searches the next, so it refuses such a run at its first batch also where no row
  looks a tax id up (`db.nz.bx.c150.*`, as before line L1).
- Rows `db.nz.cq.*` (line DBID, L1b): TBLASTN and TBLASTX without the report stream make every batch's lookups before
  any report, so they ask for that refusal only where NCBI fetches a subject: not for fields that read none (`6 qseqid
  evalue`) where `PRE_FETCH_SEQS_LIMIT` turns the prefetch off, nor for a batch without hits, whose prefetch loads no
  sequence. `nz_qx.faa` / `nz_qx.fna` hold `qx` (`M` and 70 random residues, no hit at `-evalue 1e-3`) and
  `nz_qx2.*` `qx` then `q2`: with `BATCH_SIZE=1` the refusal of `q2`'s batch still comes before any byte. TBLASTX's
  report stream (outfmt 7 in two batches), which would write `qx`'s block before `q2`'s lookups, is off for a run that
  the taxonomy cache refuses, so `nz_qxy.fna` (`qx`, then `qy` without hits either) gets NCBI's bytes.
- Rows `db.nz.cq.bn.*` (line DBID, L1c): BLASTN asks for that refusal only where NCBI fetches a subject too: not for
  fields that read none (`6 qseqid evalue`) where a positive `PRE_FETCH_SEQS_LIMIT` below the batch's Seq-aligns turns
  the prefetch off, not for queries without hits whose reports would fetch (`PRE_FETCH_SEQS_LIMIT=0`, outfmt 6), and not
  for a batch without hits (`nz_qx.fna` in outfmt 0; `nz_qxy.fna` with `BATCH_SIZE=1` in outfmt 7, where the report
  stream writes each batch's reports). `nz_qx2.fna` with `BATCH_SIZE=1` is refused at `q2`'s batch, whose subjects
  (taxids 562 and 10090) NCBI fetches without a message (the superset kept since L1b); outfmt 6 writes nothing for
  `qx` before it, and `tests/blastn_report_stream.rs` checks that outfmt 0 and 7 write `qx`'s report first where the
  report stream acts (and only outfmt 0's prolog where it does not). The rows use `-task blastn -dust no -evalue 1e-3`
  (megablast and dc-megablast find no hit of these queries).
- Rows `db.nz.cq.x.*` (line DBID, L1d): the XML2 and JSON formats (outfmt 13 to 16), whose reports fetch every reported
  subject, ask for that refusal too, where a reported subject has a positive taxid. `nz_rep.fna` holds `a1` (the
  letters of `sN` as in `nz_n`), `b1` (`sZ`) and `c1` (`sI`) with the taxids 99990, 150 and 99990 (`nz_rep.map`), `d1`
  (the codons of 60 random residues) without a taxid and `e1` (another 60) with the taxid -5: one volume where NCBI
  really writes the offset error again (for `c1` after `b1`; `nz_nc`'s only 99990 is `sN`, whose entry NCBI keeps
  loaded, so it writes the message once there). `nz_rep_q.fna` holds `qN`, `qZ` and `qI`, the whole of `a1`, `b1` and
  `c1`, and `nz_rep_q.faa` their translations (`nz_repp.faa` holds those of `a1` to `d1`): NCBI writes the message 3 to
  5 times in outfmt 15 / 16, as in outfmt 6 with `BATCH_SIZE=1` (twice). LOSAT refused these runs before any byte
  until line NP item O4b2-6, which follows the slot through the XML2 and JSON reports' loads (rows `*.off0.*`,
  `db.nz.cq.x.bn.repo.*`, `db.nz.cq.x.b[npx].reppo.*`); a prefetch of `nz_rep`'s local ids, whose order NCBI takes from
  addresses, stays refused (`*.reject`). `nz_rep_qdq.fna` is `qD` then `nz_rep_q.fna` (one batch whose first `Hit`,
  `d1`, asks the taxonomy database nothing: `db.nz.cq.x.bn.repo.qd.16`), `nz_rep_qxn.faa` the proteins `qx` (of
  `nz_qx.faa`) then `qN` (BLASTP's `db.nz.cq.x.bp.reppo.b1.qx.16`).
  `nz_rep_qd.*` (`qD`, the start of `d1`) and `nz_rep_qe.fna` (`qE`, the start of `e1`) hit `d1` or `e1` alone: no
  volume looks a positive taxid up, so their `ncbi` rows are not refused. `nz_repo` and `nz_reppo` hold `a1` to `d1`
  (nucleotide, protein) without parsed Seq-ids: their taxid is that of the definition line, a branch of its own
  (`db.nz.cq.x.bn.repo.*`: NCBI writes the message 4 times with `BATCH_SIZE=1`). `nz_rep_qxn.fna` is `nz_qx.fna`'s `qx`
  then `qN`: BLASTX writes each batch's XML2 report before it searches the next, so on `nz_reppo` it held the
  document of a run with this refusal and dropped it at `qN`'s batch (no byte of `qx`'s; since O4b2-6 it writes NCBI's
  bytes, `db.nz.cq.x.bx.reppo.b1.qx.16`); `tm_cache_qx.fna` (`qx` then the three records of `x16_bx_taxid_q.fna`)
  still does so on `tm_cache150`, whose prefetch of local ids NCBI orders by addresses and whose parsed Seq-ids make
  BLASTX hold anyway. BLASTP's search in blocks of query batches takes the search of all queries for such a run
  (`tests/blastp_query_blocks.rs`, `nz_reppo`).
- The sequences come from the line's probe script (kept outside the repository). The databases (version 4, every file
  makeblastdb writes; `nz_pc` and `nz_nc` alike with their maps), run from `LOSAT/`:

      makeblastdb -in tests/db_search/nz_p.faa -dbtype prot -parse_seqids -taxid_map tests/db_search/nz_p.map \
          -blastdb_version 4 -title "nz p" -out tests/db_search/nz_p
      makeblastdb -in tests/db_search/nz_n.fna -dbtype nucl -parse_seqids -taxid_map tests/db_search/nz_n.map \
          -blastdb_version 4 -title "nz n" -out tests/db_search/nz_n
      makeblastdb -in tests/db_search/nz_psi.faa -dbtype prot -parse_seqids -taxid_map tests/db_search/nz_psi.map \
          -blastdb_version 4 -title "nz psi" -out tests/db_search/nz_psi
      makeblastdb -in tests/db_search/nz_rep.fna -dbtype nucl -parse_seqids -taxid_map tests/db_search/nz_rep.map \
          -blastdb_version 4 -title "nz rep" -out tests/db_search/nz_rep
      D=$PWD/tests/db_search; cd "$(mktemp -d)"   # nz_repo, nz_reppo: one database per record, then one volume
      for k in a1 b1 c1 d1; do  # a1, b1, c1, d1 of nz_rep.fna / nz_repp.faa, each in $k.fna / $k.faa
          t=""; case $k in a1|c1) t="-taxid 99990";; b1) t="-taxid 150";; esac
          makeblastdb -in $k.fna -dbtype nucl $t -blastdb_version 4 -out n$k
          makeblastdb -in $k.faa -dbtype prot $t -blastdb_version 4 -out p$k
      done
      makeblastdb -in "na1 nb1 nc1 nd1" -input_type blastdb -dbtype nucl -blastdb_version 4 -title "nz repo" \
          -out $D/nz_repo
      makeblastdb -in "pa1 pb1 pc1 pd1" -input_type blastdb -dbtype prot -blastdb_version 4 -title "nz reppo" \
          -out $D/nz_reppo

## `rps/`: the RPS-BLAST profile databases (line PS, series PR, PS11a, PS11b, PS11c-1 and PS12a)

- `rps_mini` (version 5), `rps_mini4` (version 4) and `rps_mini26` (version 4, a 26-column `.rps` with the 26-letter
  magic `0x1e16` in `.rps` and `.loo`, its `.freq` left at 28 columns, as NCBI's own test database `rpstest_le`):
  three profiles of WSSV protein fragments (`gnl|CDD|101` TK-TMK residues 1-190, `gnl|CDD|102` ribonucleotide
  reductase residues 301-480, `gnl|CDD|103` DNA polymerase residues 1101-1290, from `tests/fasta/WSSV.faa`), the
  frequency ratios of two psiblast rounds against `paja` with their scores computed by makeprofiledb at its default
  scale 100 (as CDD's). `rps/build.sh NCBI_BIN OUT_DIR` records how they were made with NCBI BLAST+ 2.17.0; makeprofiledb
  writes the date into the index files, so the committed `.gz` files are the fixtures (a rebuild differs).
  `rps_pam30.aux` is `rps_mini.aux` with `PAM30` and gap costs 9/1 (row `db.rps.rej.matrix`). `empty.aux` is an empty
  file (0 bytes), put in place of `rps_mini.aux` for NCBI's engine error of a `.aux` without a matrix name (row
  `db.rpstbn.emptyaux.6`).
- `rpscol` (PS15a-1): one profile (`makeprofiledb` of a PSSM whose query id is the local id `Query_1`, an entry that
  NCBI's scope takes for the first query of a run read without `-parse_deflines`; the PSSM is the one of
  `PR/followups/rps_lcl_query_title.md`'s fixture, `tests/fasta/rps/col_qp.faa` the protein of it, `col_qn.fna` its
  back-translation). Rows `db.rps.of11.col` and `db.rpstbn.of11.col` are NCBI's archive; `db.rps.rej.of5col` and
  `db.rpstbn.rej.of5col` are LOSAT's refusals of BLAST XML (NCBI's rpsblast shows the query's title for the entry,
  NCBI's rpstblastn stops with an exception after half a document). The files keep their original names
  (`rpscol`): the `.pin` and the `.aux` name them.
- `rpscoltax` (PS15a-2, review S1): `rpscol`'s profile made again with the taxid 9606: only its
  `.phr`, `.pin`, `.pdb`, `.pot` and `.ptf` are here (the `.pin` names the `.pdb`, so the files keep that name); a row's
  cell `rps-gz:rpscoltax!aux=rpscol.aux.gz!...` takes `rpscol`'s other files, which are the same. Rows `*.rej.of18coltax`
  and `*.rej.of20txcoltax`: NCBI shows the query's taxonomy (none) for the colliding entry.
- `rpstax` (PS15a-2, review S2): `rps_mini`'s three profiles with the taxids 7777 (`gnl|CDD|101`) and 9606 (`102`) and
  none (`103`), its `.phr`, `.pin`, `.pdb`, `.pot`, `.ptf` and `.pto` here and `rps_mini`'s other files by `!EXT=FILE`;
  rows `*.of*.tax` run with `BLASTDB=tests/seqdb/taxdb/good` (after the inflated directory in the search path), rows
  `*.taxn` without a taxonomy database.
- `rps_mini3v` (version 5, PS11b): the same three profiles as three volumes of one profile each (`rps_mini3v.00` to
  `.02`, each made by makeprofiledb from one profile, and the alias file `rps_mini3v.pal`, written by `build.sh`):
  the volumes' diagonal tables have 256 entries, so the hits of a long query differ from those of `rps_mini`'s one
  volume (the merge, the per-volume hit list size, outfmt 7's `# N hits found` before the cut, volume 0's K). A row
  `rps-gz:rps_mini3v!01.freq` leaves one volume's `.freq` out (`!00.freq`: the first, which NCBI answers with its
  warning and mode 0).
- `rps_wide` and `rps_wide_ug.loo.gz` (version 5, PS12m): one profile of 133,000 columns (`gnl|CDD|901`: the first
  fragment of `rps_mini`, 190 columns, tiled 700 times by `rps/build_wide.py tile`), made by makeprofiledb with
  `-threshold 20 -force`, so that its table holds the nine words that match exactly (700 hits each, one per tile) and
  the volume has 65 buckets of 2,048 rows: two groups of the scan of `losat rpsblast` and `rpstblastn` (the first
  of 64 buckets, `RPS_GROUP_BUCKETS`, the second of one bucket, rows 131,072 to 133,000), and each cell's hits lie in both. Its `.freq` is left out (rows with `-comp_based_stats 0`
  only) and `.rps.gz` is 113 KB (a 15 MB file of one tile's rows repeated). `rps_wide_ug.loo.gz` is its `.loo` with the
  hits of the nine cells not in ascending order (reversed, rotated by a third, first and last swapped; `build_wide.py
  ungroup`), which makeprofiledb never writes: rows `db.rps.ungrouped.6` and `db.rpstbn.ungrouped.6` put it in place
  of `rps_wide.loo` (`rps-gz:rps_wide!loo=rps_wide_ug.loo.gz`), `db.rps.wide.6` and `db.rpstbn.wide.6` search the
  table as it is; NCBI's bytes are the same for both (the order of a cell is invisible to its scan), LOSAT regroups
  the cells when it opens the volume. Queries `tests/fasta/rps/wide.faa` (the 190 residues of the fragment, which hit
  each tile) and `wide.fna` (their back-translation, `build_nq.py`, seed 21).
- `db_search_fixtures.py` runs a row whose `db` cell is `rps-gz:NAME[+NAME2][!EXT][!EXT=FILE][@ALIAS]` with the files
  `rps/NAME.*.gz` inflated into a fresh directory (without `NAME.EXT`, or with `rps/FILE` in its place, inflated when `FILE`
  ends in `.gz`), `BLASTDB` set to that directory and `-db "NAME NAME2"` (a name `DIR/NAME` has its files inflated into
  the subdirectory `DIR`: a copy of a volume under its own base name, line RG D1a-4 fix 2); with `@ALIAS` the alias file
  `rps/ALIAS.pal` is copied beside them and the run searches `-db ALIAS`. A row whose `covers` cell has `path: {db_dir}` keeps that directory as
  `{db_dir}` in the frozen bytes (NCBI's RPS errors print the volume's canonical path); the directory's path then has
  the `db_dir_len.txt` characters (the name's prefix is lengthened, `db_search_fixtures.py` `temp_dir_of_length`), because
  the wraps of an archive's message that holds it (`message "..."`) depend on its length, and `{db_dir}` stands for
  that length in `wording::tests::archive_messages_are_wrapped_again_as_ncbi_wraps_them`.
- `rps/uneq2.pal` (item PS15a-1, review S1 and S2; rows `db.rps.of11.udbs`, `db.rpstbn.of5.u2m`, `db.rpstbn.u2m.7`):
  `DBLIST "rps_mini3v.00" "rps_mini"`, two volumes of different size (one profile of 190 letters, three of 560) of which
  CDD:101 and CDD:102 are in both.
- Alias files (line RG, D1a-3 review 7 S1 and D1a-4), written by hand, each `TITLE t` and a `DBLIST` of `rps_mini`
  (one volume: `*rps.pal`), `rps_mini3v` (three volumes: `*rps3v.pal`) or `rps_mini3v.01 rps_mini3v.02` (two volumes:
  `l0rps2v.pal`) with one or two keys, as named: `l0` `LENGTH 0`, `l0n0` `LENGTH 0` and `NSEQ 0`, `n0` `NSEQ 0`, `nm5`
  `NSEQ -5`, `sn7` `STATS_NSEQ 7`, `snm5` `STATS_NSEQ -5`, `st5k` `STATS_TOTLEN 5000`, `st63` `STATS_TOTLEN
  9223372036854775809`, `min1` and `min100` `MIN_SEQ_LENGTH 1` and `100`, `l563n0` `LENGTH 563` and `NSEQ 0`, `fo11` and `fo22` `FIRST_OID` and `LAST_OID` 1
  and 2 (`fo11rps3v` over the three volumes), `fo3` `FIRST_OID 3` (`fo3rps`, `fo3rps3v`), `n1000`, `n4` and `n5` `NSEQ` 1000, 4
  and 5, `l400`, `l1000`, `l2p30`, `l2p31`, `l2p63` and `l2p63p1` `LENGTH` 400, 1000, 2^30, 2^31, 2^63 and 2^63 + 1 (the last
  two 0 and negative as NCBI's `Int8`), `nm300` and `nm562` `NSEQ` -300 and -562, `l400n200`, `l0n600` and `l2p63p1n600`
  `LENGTH` 400, 0 and 2^63 + 1 with `NSEQ` 200, 600 and 600. The rows `db.rps.{1v,2v,3v}.*`, `db.rpstbn.{1v,3v}.*`,
  `db.rpstbn.fo11.e5.6`, `db.{rps,rpstbn}.{empty,nolet}.*` and the rejection rows `db.rps.rej.{l0,nm300,l2p30,l2p31,over,n5.over,fo22,l563n0.oldfsc}`,
  `db.rpstbn.rej.{l0,nm300,fo22,l563n0.oldfsc}` search them (NCBI's search units of item D1a-4: one volume is searched with
  all its keys, its `LENGTH` plus `NSEQ` the concatenated length of the search, several volumes each bare under the
  database's length and count; the queries `tests/fasta/rps/lc403.faa`, `lcm2.faa`, `lc403.fna` and `over.faa` are records
  of `q40.faa` and `nq.fna` whose trees end below the rows, or whose extension runs past the last profile; fix round 1 of
  D1a-4: `q649.faa` is record YP_009220649.1 of `tests/fasta/WSSV.faa` (its residues 1101-1290 are `gnl|CDD|103`'s
  fragment), `n103.fna` those 190 residues back-translated with the first codon of each amino acid (A GCT, R CGT, N AAC, D
  GAC, C TGC, Q CAG, E GAA, G GGT, H CAC, I ATC, L CTG, K AAA, M ATG, F TTC, P CCG, S TCT, T ACC, W TGG, Y TAC, V GTT), and
  `nolet.fa` two records without letters (`>a`, `>b`); fix round 3 of D1a-4: `tests/fasta/rps/alln.fna` (`n1`, 60 N,
  and `n2`, 40 N; 110 bytes), a nucleotide run without a valid frame, and `allx.faa` give the rows
  `db.{rps,rpstbn}.1v.{l563n0,l2p31}.*` (with `OLD_FSC=1` or a count past `Int4`, a run without a valid query is
  reported as NCBI reports it) and `db.ak.bp.nseq0.oldfsc.allx.7` (BLASTP over `aliaskeys/p_nseq0`);
  `tests/fixtures/strategy/cases.tsv` copies some of them (`rpsblast.{3v,2v}.*.q40`, `rpstblastn.3v.*.nq.b2000`).
- Queries: `tests/fasta/WSSV.faa` (177 records, 10 query batches; its long records share diagonals of the 1,024-entry
  table of these databases), `tests/fasta/rps/q40.faa` (33 WSSV records, 11,434 letters, two query batches: the outfmt
  0 K of each batch's first query), `tests/fasta/rps/xoj.faa` (a record of X only, the TK-TMK record with an O at 30
  and a J at 60, part of the ribonucleotide reductase) and `tests/fasta/rps/mask.faa` (PS11c-1, the masks:
  `YP_009220478.1` of `WSSV.faa` as it is, whose acidic low-complexity run SEG masks, and `YP_009220601.1`, the TK-TMK
  record, as `tk_lc` with its residues 61-120 in lower case for `-lcase_masking`; lines of 60 letters, 819 bytes. Made
  from `WSSV.faa` by taking the two records, renaming the second and lowering the case of its letters 61-120), and
  `tests/fasta/rps/nident_gate.faa` (PS11c-1 review: record `BDT63522.1` of `tests/fasta/PajaWSV.faa`, lines of 60
  letters, 1,563 bytes, with 401 of its 1,473 letters in lower case at random: Python's `random.seed(7)`, then one
  `random.random() < 0.3` draw per letter over the first 30 records of `PajaWSV.faa` in file order, a letter lowered
  where it is true; with `-lcase_masking` some of its HSPs have identities on masked letters only).
- rpstblastn's queries (PS12a, rows `db.rpstbn.*`): `tests/fasta/rps/nq.fna`, `nedge.fna` and `odd.fna`, written by
  `tests/fasta/rps/build_nq.py` (seeded back-translations of `q40.faa`, the docstring there says how): `nq.fna`
  (12,440 bytes, two query batches of 10,000 letters; its first query `stopq` has stop codons only in frame +1, then 15
  records of alternating strands and flanks), `nedge.fna` (5,498 bytes: queries of 1 to 12 letters, cuts of 570 to 572
  letters and their reverse complements, a run of N, a lower-case run, ambiguity codes) and `odd.fna` (a 4-letter
  query, an empty one, one of N); `-query_loc 1-2` on `nedge.fna` gives a batch without a valid frame. The empty query
  file is `tests/db_search/empty.faa`.
- PS12a review A rows: `tests/fasta/rps/bound1.fna` (11,101 bytes: `n9999`, 9,999 letters of N, and `t2` of 2 letters,
  which fill rpstblastn's first query batch of 10,000 letters without a valid frame, then `v900`, 900 letters of a
  valid record in the next batch; lines of 60 letters) and `tests/fasta/rps/allx.faa` (42 bytes: two records of X
  only, 20 and 12 letters, a protein batch without a valid query), written by hand for the review's probes, and
  `tests/fasta/rps/tkcut570.fna` (the record `tkcut570` of `nedge.fna` alone, 589 bytes: rows `db.rpstbn.ll1.0`, `ll.0`,
  `sort.0`, `3v.0`, `mix2v.*`), and `tests/fasta/rps/nq4.fna` (the first four records of `nq.fna`, 2,922 bytes: rows
  `segsoft.0`, `pdsg.0`; review B S1 kept PS12's outfmt 0 rows small). Rows `db.rpstbn.mix2v.*` and `db.rps.mix2v.ds.0`
  search `rps_mini` (3 profiles) and `rps_tie.00` (1 profile) as two volumes (`rps-gz:rps_mini+rps_tie.00`).
- PS12f rows `db.rps.nopsq.later.*`: `tests/fasta/rps/b2v01.faa` (records `YP_009220493.1` and `YP_009220507.1` of
  `q40.faa`, in that order; 659 bytes) with `BATCH_SIZE=1` on `rps_mini3v` without `rps_mini3v.01.psq`: the first query
  hits volume 00 only and the second volume 01, so NCBI writes the first batch's report (without the epilog) and then
  fails the second batch's traceback (`BLAST Database error`, exit status 2). Rows `db.rps.noloo.inv.*` and
  `db.rps.3v.noloo1.inv.7` (review S1): `tests/fasta/rps/xfirst.faa` (`x1`, 20 X, then `YP_009220493.1` of `q40.faa`;
  386 bytes) with `BATCH_SIZE=1`, a first batch without a valid query, on a database without a volume's `.loo`: NCBI
  sets up every volume's RPS files in the first batch's set-up whatever its queries, so that batch fails (`BLAST
  engine error`, exit status 3).

## outfmt 11 titles where NCBI's search fails (line NP, series SM, DECISIONS 25)

Rows `db.archtitle.*.nosq.*` search databases without their sequence file:

- `arch_sprot_nopsq.pin` and `.phr` are copies of `arch_sprot.pin` and `.phr` (BLASTP, BLASTX); `arch_tbn_snucl_nonsq.nin`
  and `.nhr` are copies of `arch_tbn_snucl.nin` and `.nhr` (TBLASTN, BLASTN, TBLASTX). NCBI reads the sequence file first in
  the preliminary search of the first query batch with a valid context and fails that batch (`BLAST engine error: Database
  memory map file error`, exit 3); the batches before it are written.
- `aliaskeys/p_nopsq_short.pal` (line RG item X-bpgap, fix round 6) is an alias over `arch_sprot_nopsq` with statistics of
  2 letters in 3 sequences (`STATS_TOTLEN 2`, `STATS_NSEQ 3`): without the Gumbel block (`-searchsp`, `OLD_FSC`, an
  ungapped search) the setup's average subject length (0) fails before the search's read of the missing sequence file,
  and with it the read fails (rows `db.x4.chunk.avgdup.nps.*`). `../fasta/x4/chunk_gapvalid_tk.faa` is
  `chunk_gapvalid.faa` then `chunk_tk.faa` (row `db.x4.chunk.neg.gapvalid.tk.nopsq.7`).
- The queries `tests/fixtures/fasta_input/batch_title_unsearched.{prot.faa,nuc.fa}` and `batch_title_unsearched_next.*`
  hold `u1` and `u2`, 30 X (protein) or 90 N (nucleotide) only, which have no valid context, then the sequences of q2
  and q3 of `batch_title_first` as `q3` and `q4`; `café` is in the title of `u1` or of `q3`. With `BATCH_SIZE` 60 (180)
  the first batch is `{u1, u2}`.

A database without its header file fails the first batch whose hits read a definition line, with a message that holds the
file's absolute path; that case is `a_title_in_the_batch_that_a_missing_header_file_fails_is_not_refused` in
`tests/outfmt11_titles.rs` (a temporary database made by `LOSAT makeblastdb`).

## Databases with column files (line SW, SW8 follow-up 2)

Rows `db.col.*` search databases whose volumes have column files: a column is the files `.?Xa` (index), `.?Xb`
(big-endian data) and `.?Xc` (little-endian data) for a letter X, and `makeblastdb -mask_data` writes the masks as the
column `BlastDb/MaskData` (`.paa .pab .pac`, `.naa .nab .nac`). NCBI's searches open a volume's columns only through
`CSeqDBVol::x_OpenAllColumns`, which runs only for a masking algorithm named by `-db_soft_mask`, `-db_hard_mask` or an
imported search strategy (`-import_search_strategy`; LOSAT refuses all three: rows `db.col.rej.*` for the first two),
so a search without one writes the bytes of the same database without its columns, whatever the column files hold.
`coldb/` holds what `make_coldb.py NCBI_BIN` builds (makeblastdb, segmasker, dustmasker and windowmasker 2.17.0; the
records and taxids of `taxfilter/tfp.faa`, `tfp.map`, `tfn.fna` and `tfn.map`; `.pjs` and `.njs` dropped), and the
files below that it does not build:

- `p5m`: version 5 protein, `-parse_seqids`, taxid map, the seg masks of `segmasker -parse_seqids` (algorithm 21).
- `p4m`: version 4 protein without `-parse_seqids`, the seg masks of `segmasker`.
- `n5m`: version 5 nucleotide, `-parse_seqids`, taxid map, two volumes (`-max_file_sz 4KB`), each with its column, the
  masks of `dustmasker -parse_seqids` (11) and `windowmasker -parse_seqids` (30).
- `n4m`: version 4 nucleotide without `-parse_seqids`, the masks of `dustmasker`.
- `dmg/CASE/p4m` (twelve cases) and `dmg/CASE/n4m` (`idx_junk`, `data_cut`, `be_only`): copies of `p4m` and `n4m` whose
  column is damaged (`make_coldb.py` `DAMAGE`): the index cut, junk, empty or of format version 2, the data files cut or
  junk, one data file or both missing, the index missing, the column of a database of 6 sequences, the column twice under
  two letters. NCBI's column reader stops on most of them (`blastdbcmd -info`, or a search with `-db_soft_mask`); its
  searches without a database mask do not read them.
- `ucol/user-column-db.*` (not built by the script): every file of NCBI's
  `c++/src/algo/blast/unit_tests/seqdb_reader/data/user-column-db`, copied as it is (598d8ae6; six proteins with parsed Seq-ids and the user column `comedy`, its big-endian data file `.pab` alone; NCBI's
  tools cannot write a user column). `tests/seqdb/db` holds the same database without its ISAM files, on which NCBI's
  outfmt 6 and 7 stop with a null pointer exception, columns or not. `uq.faa` (calfumirin-1, XP_642131.1) and `uq.fna`
  (its back-translation by hand, one fixed codon per amino acid; neither built by the script) are the queries of its
  rows.

The other queries are `taxfilter/q.faa` (BLASTP, TBLASTN, PSIBLAST) and `taxfilter/q.fna` (BLASTN, BLASTX, TBLASTX).

## ADAPTIVE_CBS (line SW, env parity G2+G3)

- `acbs` (protein) and `acbs_n` (nucleotide, the same 1,700 subjects back-translated), with the queries
  `acbs_q.faa` and `acbs_q.fna`: written by `make_acbs_inputs.py` (its docstring gives the commands and the
  shape). The rows `db.acbs.*` run them with the environment variable `ADAPTIVE_CBS`: NCBI BLAST+ keeps 1,500
  preliminary subjects instead of 1,050 and cuts them at 1,100 (BLASTP, BLASTX with `-seg no`, TBLASTN with
  `-seg no`: 400 subjects reported, 350 without the variable; `-max_target_seqs 600`: 500, 550 without);
  PSIBLAST keeps the 1,500 and does not cut them. `db.acbs.rpsblast.*` run `rps/rps_mini` with the
  variable: a database of three profiles shows no difference, only that the variable is read (a local
  database of 1,700 single-sequence profiles of the same subjects gave rpsblast 351 -> 401 and 551 -> 501,
  too large to commit).
- `db.acbs.blastx.wrap*` and `db.acbs.blastx.crash*` search `bslice` with a `-max_target_seqs` whose
  preliminary size passes NCBI's 32-bit int: 2147483624 (NCBI keeps 2 subjects per query with
  composition-based statistics), 1073741824 without them (NCBI keeps 10), and 1073741824 with them (NCBI crashes). LOSAT computes the size
  without the overflow (approved exception 19 of PD-LOSAT-NCBI-DEFECTS) and keeps every hit: each row is
  compared with NCBI's run at 1073741798, its `ncbi-near` row `db.acbs.blastx.near*`.

## Alias-file keys (line RG, series D, D1a)

Rows `db.ak.*`: alias files of `aliaskeys/` over `taxfilter/tfp`, `taxfilter/tfn` (its volumes `tfn.00` and `tfn.01` too) and
`multidb/p4a` (PSIBLAST reads databases without parsed Seq-ids), written by hand for the keys NCBI's `CSeqDB` reads
(`remaining-gaps-20261112/D/fixtures/make_d1a_fixtures.py` writes them and the rows): `NSEQ` and `LENGTH` that are not the
sums of the volumes (`p_nseq`, `p_nseqhi`, `p_nseqonly`, `p_nseqplus`, `p_nseq0`: each query's warning `Filtering resulted
in an empty database.` and a search that goes on, `p_nseq4g`: NCBI's `int` of 4294967297 is 1, `n_nseq`, `n_nseq4g`,
`s_range`), `STATS_NSEQ` and `STATS_TOTLEN` (`p_stats`, `p_statsn`, `p_statsl`, `p_statsz`, `p_statsnest`, `p_statsover`, the
`-dbsize` row, `p_statsl` with `p_r123` in one `-db`, `n_stats`, `s_stats`), `MIN_SEQ_LENGTH` (`p_min`, `n_min1`: TBLASTN's
`The average subject length is too short`, exit 3, `n_min1000`), `MAX_SEQ_LENGTH abc` (never read), `VDBLIST` (`p_vdb`),
`FIRST_OID` / `LAST_OID` ranges (`p_range`, `p_nest` and `p_rin`: nested ranges ANDed, `p_union`: two alias files over one
volume ORed, `p_twice`: an alias file and the volume itself, `p_dupvol`: one volume twice in a `DBLIST`, `p_first0`, `p_last2`,
`p_empty`, `p_beyond`, `p_r312` and `p_r312s` with `-taxids` and `-negative_taxids`, `p_range` with `-seqidlist`, `n_range`,
`n_mv`: the global OIDs 9-13 of the two volumes, `n_vol1`: the second volume alone; `p_lo`, `p_hi`, `p_halves`, `p_cover`,
`p_coverr`, `p_arrcov`, `p_chain`, `p_above` and `p_abovecut` for the rows `db.gisl.*`). The ncbi rows were frozen with NCBI
BLAST+ 2.17.0 on the workstation (`env_blast`; two rows checked against the release build on the developer machine, the same
bytes); the `db.ak.bp.err.*` rows (NCBI's `CStringException` texts, exit 255, which carry the path of NCBI's build tree:
`p_bad_nseq`, `p_bad_first`, `p_bad_min`, `p_bad_stlen`, `p_bad_two` without and with `-taxids`) with the release build on the
developer machine.

Rows `db.ak.*.e0*` (D1a-2, `remaining-gaps-20261112/D/fixtures/make_d1a2_fixtures.py`): alias files `e_p*` over the
protein volume of no OID `../seqdb/db/empty` and `e_n*` over `aliaskeys/en` (a version 4 nucleotide index of no sequence)
with `NSEQ 1` and `LENGTH 100` (A), `LENGTH 100` (B), `STATS_TOTLEN 100` (C), `STATS_NSEQ 1` and `STATS_TOTLEN 100` (D),
`NSEQ 1` (F) and `NSEQ 1` and `STATS_TOTLEN 100` (G): NCBI reads the length of OID 0 only where the statistics' total is
not positive (F: `Database memory map file error`, PSIBLAST `OID not found`); otherwise it searches no sequence, with each
query's warning where `GetNumSeqs` is 0 (B, C, D); BLASTN and TBLASTX divide B's and C's 100 letters by 0 sequences
(SIGFPE: LOSAT's refusal). Frozen with the release build on the developer machine (2026-11-13; the same bytes as the
workstation's probes p4c_dp, p4c_dn).

Rows `db.ak.*.nseq0*`, `db.ak.*.e0len*`, `db.ak.*.e0stl*` (D1a-5 of line RG): a statistics count of 0 (`NSEQ 0`,
or `e_nB`'s and `e_nC`'s totals over a volume of no OID) where every query has no letters to search (`../fasta/rps/alln.fna`,
`allx.faa`): NCBI divides by the count only in the setup of a batch with a valid context, so it reports the run (exit 0);
the rejection rows give it a valid query, where NCBI crashes (SIGFPE) and LOSAT refuses: `.b1` with `BATCH_SIZE=1` and
`../fasta/rps/alln_q.fna` (`alln.fna`, then `taxfilter/q.fna`) or `allx_q649.faa` after the two reports NCBI writes
before its crash (`ncbi-near`), `nseq0ns` over `aliaskeys/n_nseq0ns` (`NSEQ 0` over `arch_tbn_snucl_nonsq`, which has no
sequence file) before the error of the missing file, which NCBI meets only after the division (TBLASTN with the Gumbel
block, which does not divide: that error). Frozen with NCBI BLAST+ 2.17.0 on the workstation (`env_blast`).

Rows `db.cnt.avg2p32.mb.*`, `db.cnt.avg2p32.mo00.split.6`, `db.cnt.avg2p32.mo52.split.6` and `db.ak.bn.rej.nseq0.mbsplit` (D1a-6 of line RG): a query batch split into chunks
(`CHUNK_SIZE=200`) over `dbcounts/avg2p32`, whose average length is too short, with options that LOSAT refuses for the
greedy extension (`-task megablast` with gap costs 40000 and 40000, or reward 0 and penalty 0): NCBI checks the subjects in
the setup of the first chunk with a valid context, before the greedy set-up, so it stops with `The average subject length is
too short` (exit 3), as an unsplit batch does (`.mbwhole.6`). `.b1.7` reports two queries without letters first (`BATCH_SIZE=1`,
`../fasta/rps/alln_q.fna`); `.nq` puts a query of 300 N (`../fasta/rps/n300_q.fna`: that query, then `taxfilter/q.fna`) before
it, so the first chunks have no valid context. The rejection row has `NSEQ 0`, where NCBI crashes in that setup (SIGFPE).
Rows `db.gr.bn.*` (the same item): `db.gr.bn.rej.mo00.b1` refuses matrix-only scoring with the greedy extension in the third
batch of `alln_q.fna` (`BATCH_SIZE=1`) after the reports of the two batches before, the 142 bytes NCBI writes before it
divides by zero there (SIGFPE); `db.gr.bn.alln_q.b1.7` is the near run NCBI finishes (`ncbi-near`).
Frozen with NCBI BLAST+ 2.17.0 on the developer machine and again on the workstation (`env_blast`), the same bytes.

Rows `db.ak.*.ord.*` (D1a-2): the order of a search database's errors, NCBI's (probes p4c_e*, p4c_f2): the list, the alias
files' counts, the formatter's reopening of each name, then the headers (before line DBID mixdb also the accession tables). `o_mix` and `o_mixs` are
`multidb/na` (with the accession table) and `multidb/nnp` (without) under `NSEQ abc` / `STATS_NSEQ abc`, with `-seqidlist`
(the list's `DB contains no accession info.` first) and `-taxids` (the alias file's error); `o_np` is `nnp` alone, `o_v4`
the version 4 `multidb/n4a` (`Taxonomy filtering is not supported in v4 BLAST dbs` first); `o_dh` and `o_dhs` are over
`dh`, a copy of `multidb/p4a` whose header file is 4,096 bytes of 0xff (the alias file's error, before any header is
read). Frozen with the release build on the developer machine (the `CStringException` texts carry its path).

Rows `db.oir.*` (D1a-2 part 2, `remaining-gaps-20261112/D/SURVEY_D1a2_A.md`, `remaining-gaps-20261112/D/fixtures/make_d1a2a_fixtures.py`
and `build_o.sh`): the probe cases `p4a_*` of that survey, each command of a search program a row (the row id ends with the
command's index in the probe case; the commands with `-num_threads 4` are left out, the harness runs each row at 1 and 4
threads; `-taxids` and `-negative_taxids` with `-no_taxid_expansion`). `oidrange/` holds the cases' alias files
(`<case>_<name>.pal|nal`, over `../taxfilter/tfp`, `tfn`, `tfn.00` and `o`) and `o`: `o.faa`, ten sequences of `tfp` as
`lcl|o0`..`lcl|o9`, made with NCBI BLAST+ 2.17.0 `makeblastdb -in o.faa -dbtype prot -out o -parse_seqids -blastdb_version 5
-taxid 10116 -title o`; NCBI sorts the volumes by name, so `o` holds the OIDs 0-9 and `tfp` 10-32 wherever both are in a
`-db`. `FIRST_OID` and `LAST_OID` of 2^31 or more (NCBI's negative `int`s: a negative first OID keeps no OID, a negative last
OID is no bound and ends the iteration range at 0, `FIRST_OID 2147483648` gives `INT_MAX`), and the iteration range of the only
alias file of a `-db` with another database or volume (a search reads the listed OIDs inside it; with `-taxids` or
`-seqidlist` every listed OID, and the statistics count NCBI's loop over the list, which can pass the first listed OID and the
OID at the range's end). Frozen with the release build on the developer machine (2026-11-13).

Rows `db.cnt.*` (D1a-3 part 1, `remaining-gaps-20261112/D/SURVEY_D1a3.md`, `remaining-gaps-20261112/D/fixtures/make_d1a3_fixtures.py`):
`dbcounts/` holds alias files `<name>.pal` over `../taxfilter/tfp` and `<name>.nal` over `../taxfilter/tfn`, each with
`TITLE t` and one key: `neg1` `LENGTH 9223372036854775809`, `negm1` `LENGTH 18446744073709551615`, `l2p63` `LENGTH
9223372036854775808` (negative as NCBI's `Int8`), `l4295` `LENGTH 4294967296000000000` (NCBI's commas lose a digit),
`l4294` `LENGTH 4294967295999999999`, `l5e18` `LENGTH 5000000000000000000`, `nsm5` `NSEQ -5`, `ns2p31` `NSEQ 2147483648`,
`nsm1` `NSEQ 4294967295` (-1 as NCBI's `int`), `ns2e9` `NSEQ 2000000000`, and `avg2p32` (`NSEQ 23` / `16` with `LENGTH`
23 x 2^32 / 16 x 2^32: an average whose `Int4` is 0). Frozen with the release build on the developer machine (2026-11-14).
`l0b.pal` is `LENGTH 0` over `../bslice` (`DBLIST ../bslice`; rows `db.cnt.l0b.*`: `-subject_besthit` removes HSPs of
`bslice_q.fna` there, which it does not over `tfp`). `x300r.faa` is a query of 300 `X` (no countable letter: NCBI's
Karlin-Altschul setup fails on it, "Could not calculate ungapped Karlin-Altschul parameters") followed by a 142-letter
random protein (`random.Random(1002)` choices of `ACDEFGHIKLMNPQRSTVWY`), for the split-query rows
`db.cnt.*.bp_c0.split.fsc.inv`: with `CHUNK_SIZE=80` a chunk's first context is a part of the random query while the
batch's first query is the invalid one (review 2 R2-M1 of D1a-3). Frozen with the local NCBI BLAST+ 2.17.0 release build
(2026-11-14). `l0q3b.pal` is `LENGTH 0` over `../q3b_chains` (`DBLIST ../q3b_chains`; rows `db.cnt.l0q3b.*`: psiblast
`-sum_stats true` with the statistics of each subject, whose HSPs with `tests/fasta/psiblast/tk.faa` form linked sets:
the rounds with the matrix guard the linking's small-gap cutoff at each subject's length, as round 2 of the
`-num_iterations 2` row does with the PSSM, the `-in_pssm` rows the linking's search space at each subject's length; review
3 R3-M1 and review 5 R5-S1 of D1a-3). At their `-evalue` of 100 and 1000 the hit-saving cutoff is
already below the sum statistics' cutoff, so its lowering does not act there (review 4 R4-S1). `l0q3w.pal` is `LENGTH 0`
over `../q3b_weak` (rows `db.cnt.l0q3w.*`, `-evalue 1e-6`, the hit-saving cutoff of each subject lowered to the sum
statistics' cutoff at the subject's length: with the matrix the weak parts of `weak10` and `weak11` score below the
hit-saving cutoff and above the lowered one, so only that lowering keeps them for the preliminary stage, which links each
to its strong part with a small gap, and a build without the lowering prints their strong parts only; with the PSSM
`weak10`'s weak part does so, while the PSSM's extension joins the weak parts of `weak11` and `gap7` to their strong HSPs;
`long8k`'s weak part (43 with the matrix, 42 with the PSSM) scores below the cutoff lowered at its length of 8,000 and above
one lowered at the setup's subject length, which a build that lowers at that length prints; the weak parts of `gap7`
(with the matrix) and `gap7p` (with the PSSM), whose two stretches score below the gap trigger, are never extended, as the
initial words' cutoff reads the hit-saving cutoff before the lowering: a build whose words read the lowered cutoff prints
them; `-xdrop_gap_final 15` keeps the traceback from joining the parts). Frozen with the local NCBI BLAST+ 2.17.0 release
build (2026-11-14 and 2026-11-15). `l0txr.nal` is `LENGTH 0` over `../txrand` (rows `db.cnt.l0txr.*`: tblastx at
`-num_threads 1`, whose frame loop saves the two-hit extensions against each subject's word cutoffs, not the setup's, which
are 1 for a database length of 0; review 6 R6-M1 of D1a-3; `db.cnt.l0txr.tx.w0`: the one-hit word finder of `-window_size 0`
the same, review 8 R8-S1). Frozen with the local NCBI BLAST+ 2.17.0 release build
(2026-11-15). `l0paja.pal` is `LENGTH 0` over `../paja` (`DBLIST ../paja`; rows `db.cnt.l0paja.*` with `bpqb5.faa` and
`-evalue 0.01`: blastp's preliminary stage in blocks (outfmt 10) and its block search (outfmt 6, blastp and psiblast's one
round), each block's call with the statistics of each subject; a build whose blocks take the setup's statistics prints 130
rows instead of 95; review 7 R7-S2 of D1a-3). Frozen with the local NCBI BLAST+ 2.17.0 release build (2026-11-15). `l2p63b.pal` is `LENGTH
9223372036854775808` over `../bslice` (row `db.cnt.l2p63b.bx_c0.16`: a database without parsed Seq-ids or taxonomy, so that
only the report sums can make BLASTX's library hold its XML2 start before the refusal; review 3 R3-S1 of D1a-3).

## psiblast's id options (line RG, series Q, Q1)

- `q1_sq1` (protein, version 5): `makeblastdb -in tests/fasta/psiblast/q1_subj_q1.faa -dbtype prot -parse_seqids
  -out tests/db_search/q1_sq1` (NCBI BLAST+ 2.17.0, run from `LOSAT/`): its first sequence is `lcl|Query_1`, the
  Seq-id of an `-in_pssm` query, which NCBI's `-parse_deflines` then shows as the PSSM's query (rows
  `db.ps.q1.*q1db*`).
- `tests/fixtures/pssm/q1_tfp.asn.gz`: NCBI BLAST+ 2.17.0 `psiblast -query tests/db_search/taxfilter/q.faa -db
  tests/db_search/taxfilter/tfp -save_pssm_after_last_round -out_pssm` (gzip -n -9), the `-in_pssm` of the rows
  `db.ps.q1.in.*`.
- `tests/fixtures/pssm/q1_tfp_id5.asn.gz`: the same PSSM with its query's `local str "Query_1"` replaced by `local id 5`
  by hand (rows `db.ps.q1.in.id5.*`).
- `tests/fixtures/pssm/ps15b_rnd_id0.asn.gz` and `ps15b_rnd_idmax.asn.gz`: `tests/fixtures/db_search/db.ps.asn.rnd.6.asn`
  with its query's `local str "Query_1"` replaced by `local id 0` and `local id 2147483647` by hand (gzip -n -9; rows
  `db.ps15b.fix.id0.*` and `db.ps15b.fix.idmax.*`: the Seq-aligns keep the integer Object-id).

Query-anchored inputs (rows `db.of1.bp.*`, line RG T3a: BLASTP `-outfmt 1` to `4`), in `of1/`, from the
survey's probes (`T/t3/scripts/gen.py`, a seeded generator, and three hand-made records):

- `qp.fa`: three queries: `qa` (150 letters), `qnone` (no hits with `sp`: two weak ones), `qc`; `qa.fa` and `qc.fa`
  hold one of them each. `qgi.fa`: `qa` under `gi|7777|ref|XP_000009.1|` (with `-parse_deflines` its rows are
  labelled `XP_000009.1`, with `-show_gis` too `7777`).
- `sp.fa`: six subjects of `qa` and `qc`: `s_multi` (two HSPs with `qc` whose rows come in the order of their
  longest segments, not of their E-values), `s_ins` (insertions and a deletion), `s_ins2` (an insertion at the
  same place), `s_twice` (two HSPs on the same query range), `s_far` (weak), `s_c`. `sp_ids.fa`: the same with
  Seq-ids of each kind (`ref`, `sp`, `gi|...|ref`, `gnl`, `lcl`, `pdb`); `s_multi.fa`: `s_multi` alone (with
  `qc` its two HSPs leave the query stretch 11-37 between them, which the query row fills).
- `q100.fa`, `s100.fa`: one identical pair of 100 letters (`-line_length 99`: the start column is two digits
  wide, from the largest 0-based position 99).
- `s_twice_lc.fa`: `s_twice` with the lower-case runs 85-92 (in the ranges of both HSPs with `qa`) and
  120-125 (in the second's only). With `-lcase_masking` NCBI maps each run to the columns of the first row
  of the subject whose letters it meets and puts those columns in lower case in every row of the subject
  (the rows' own letters there are other positions of the subject):

  NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4300-4329
- `qseg.fa`, `sseg.fa`: `qa` with a run of 20 Q at 76-95, which `-seg yes` masks; `sseg` is the same with every
  11th letter changed, and `s_ins`.
- `sp_dupid.fa` (two subjects `lcl|dup1`): the `-parse_deflines` refusal `db.of1.bp.dupid.reject`. A query with a
  subject's id is refused in the query's batch, after the prolog as in outfmt 0, which a rejection row cannot
  hold: `tests/parse_deflines_subject.rs` checks it in formats 0 to 4.
- `self1.fa` and `selfseg3.fa` (the T3a review's probes, `T/t3/review_t3a/in/`), each searched against itself
  with `-parse_deflines`, so that the subject is the query's own Bioseq. `self1` is a record A, a spacer, A with
  GG inserted, a spacer, A with a quarter of its letters changed and PW inserted before its last 10, and a
  spacer: with `-query_loc 1-48 -subject_loc 49-223` its HSP with the third copy holds the longest segment
  and takes the subject's first row from the HSP with the second copy, which NCBI then shows below it
  (`db.of1.bp.self1.*`, `alnmerger.cpp:499-501`). `selfseg3` is 30 letters, 80 Q and 40 letters: `-seg yes`
  splits its self diagonal into two ungapped HSPs, the second the longer, where NCBI aborts
  (`db.of1.bp.selfseg3.reject`, PD-LOSAT-NCBI-DEFECTS 1.23).
- `deflineonly.fa` (a record without letters) and `allx_q3.fa` (an all-X query, then a query): runs that stop
  at their first query batch, which write outfmt 0's prolog (`db.of1.bp.defl.*`, `db.of1.bp.emptydb.*`, the
  empty database `tests/seqdb/db/empty`).
- `lo24.fa`, `lp17.fa`, `lp20.fa`, `lp58.fa` (review 2 of T3a, `T/t3/review2_t3a/in2/`, `in3/`), each searched
  against itself with `-parse_deflines` and a `-query_loc` and `-subject_loc` that keep one copy as the query:
  the HSP of the longest segment follows a subject insertion, so the subject row that keeps its other segments
  ends with that insertion, after the alignment's last column, where NCBI writes its `\` and `|` onto the
  strings' terminators (`db.of1.bp.lo24.f1`, `lp17.f2`, `lp58.f2`; `lp20.ll20.f1` with `-line_length 20`).
- `qid` (the database of `qid.fa`) and `qid_q1.fa`, `qid_q2.fa` (review 2 of T3a, D2): database sequences with the
  ids that queries read without `-parse_deflines` get, `lcl|Query_1` (the letters of `qid_q1.fa`'s query) and
  `lcl|query_2` (those of `qa`), and `other1`. NCBI's scope finds a query of the batch for such a sequence and
  shows the query in its place (and aborts on its merge with `-seg yes`): `db.of1.bp.qid*.reject`; with
  `BATCH_SIZE=1` the second query is in the next batch and NCBI shows the database sequence
  (`db.of1.bp.qid2.e1.b1.f2`).
- `dispid` (the database of `dispid.fa`), `dispref` (of `dispref.fa`) and the queries `dispid_q.fa` (three
  records without ids), `dispid_q2.fa` (the second of them alone) and `dispid_pq.fa` (`lcl|abc`, `lcl|def`, read
  with `-parse_deflines`) (review 3 of T3a, D1; `T/t3/review3_t3a/`): `gi|777|lcl|Query_1` and `gi|779|lcl|abc`
  (the first query's letters with a tenth changed), whose Seq-align id is the GI and whose description line
  shows the display id `lcl|Query_1` or `lcl|abc` (`GetDisplayIds`, the `WorstRank` best id); NCBI looks the
  title up again by that id and shows the query's title there (`db.of1.bp.dispid.*.reject`, `psiblast`'s
  `db.of1.ps.dispid.pd.f0.reject`), not past `-num_descriptions` nor its `-sorthits` order
  (`db.of1.bp.dispid.nd0.f3`, `db.of1.bp.dispid.q2.s3.f1`). `ref|NP_000777.1|lcl|Query_1` shows the display id
  `ref|NP_000777.1|`, so NCBI shows the database sequence (`db.of1.bp.dispref.f1`).
- The rows `db.of1.bp.otx.*` (review 3 of T3a, D2) search `../taxfilter/tfu_pr` (an alias of `FIRST_OID 2` and
  `LAST_OID 4`) and `../taxfilter/tfu_p` with `-taxids 9606`, without parsed Seq-ids: the rows are labelled with
  each sequence's OID, the tag of `gnl|BL_ORD_ID|<OID>`.

The databases (version 4; `.pjs` deleted). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/of1/sp.fa -dbtype prot -out tests/db_search/of1/sp \
        -title "subjects of the query-anchored rows" -blastdb_version 4
    makeblastdb -in tests/db_search/of1/sp_ids.fa -dbtype prot -out tests/db_search/of1/sp_ids \
        -title "subjects of the query-anchored rows with Seq-ids" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/qid.fa -dbtype prot -out tests/db_search/of1/qid \
        -title "a database sequence under the id of a query" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/dispid.fa -dbtype prot -out tests/db_search/of1/dispid \
        -title "database sequences whose display id is a query id" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/dispref.fa -dbtype prot -out tests/db_search/of1/dispref \
        -title "a database sequence whose display id is not its local id" -parse_seqids -blastdb_version 4
    rm tests/db_search/of1/*.pjs

BLASTN and TBLASTN inputs of the query-anchored rows (rows `db.of1.bn.*` and `db.of1.tn.*`, line RG T3b), in
`of1/`, from the T3b probes (`T/t3/scripts_t3b/gen3b.py`, a seeded generator; `tnq1.pssm` from NCBI's
`psiblast -query tnq.fa(first record) -db <six mutated copies> -num_iterations 2 -out_pssm`):

- `bnq.fa`: two nucleotide queries, `bnq` (300 letters with a run of 25 T that DUST masks) and `bnq2`; `bns.fa`:
  `s_pmo` (a plus-strand and a minus-strand HSP over the query's 61-160), `s_mins` and `s_mins2` (minus-strand
  copies with insertions of 4 and 2 letters at one query column), `s_mdel` (a minus-strand deletion), `s_plus`
  (a plus-strand insertion at the same column); `bnslc.fa`: a subject whose first 40 letters are lower case
  (`-lcase_masking`: the mask shows on the columns of its first row, also in its minus-strand row);
  `bnids.fa`: the subjects of `bns.fa` under `gi|...|gb|`, `ref|`, `lcl|` and `gnl|` ids.
- `bself.fa` (copies of a 80-letter stretch: the second with a 5-letter insertion, the third ungapped and
  mutated; `-query_loc 1-80 -subject_loc 81-325 -parse_deflines`, NCBI aborts) and `bselfrc.fa` (a record with an
  inverted repeat and a lower-case run; with itself under `-parse_deflines` its rows of both strands).
- `bnqid.fa` (`lcl|Query_1`) and `bnqgi.fa` (`gi|777|lcl|Query_1`): database sequences named as the first query
  read without `-parse_deflines`, which BLASTN refuses in outfmt 1 to 4 after its prolog
  (`tests/parse_deflines_subject.rs`).
- `tnq.fa`: two protein queries (`tnq` with 16 Q that SEG masks); `tns.fa`: `t_pm` (a plus and a minus frame over
  the query's 61-80), `t_fr` (two plus frames over 46-55), `t_ins`, `t_mins` (insertions on a plus and a minus
  frame), `t_mdel`, `t_code` (codons that genetic code 2 translates otherwise: AGA, AGG, ATA, TGA);
  `tnsids.fa`: four of them under parsed ids; `tnqid.fa` (`lcl|Query_1`) and `tnqgi.fa` (`gi|778|lcl|Query_1`):
  the refusals `db.of1.tn.qid.reject` and `db.of1.tn.qgi.reject` (NCBI's merge reads the protein query as the
  sequence's nucleotides; its description line shows the query's title).
- `bnsr.nal`: an alias of the OIDs 2 to 4 of `bns` (rows labelled by OID).

The databases (version 4; `.njs` deleted). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/of1/bns.fa -dbtype nucl -out tests/db_search/of1/bns \
        -title "nucleotide subjects of the query-anchored rows" -blastdb_version 4
    makeblastdb -in tests/db_search/of1/bnids.fa -dbtype nucl -out tests/db_search/of1/bnids \
        -title "nucleotide subjects of the query-anchored rows with Seq-ids" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/bnqid.fa -dbtype nucl -out tests/db_search/of1/bnqid \
        -title "a nucleotide database sequence under the id of a query" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/bnqgi.fa -dbtype nucl -out tests/db_search/of1/bnqgi \
        -title "a nucleotide database sequence whose display id is a query id" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/tns.fa -dbtype nucl -out tests/db_search/of1/tns \
        -title "translated subjects of the query-anchored rows" -blastdb_version 4
    makeblastdb -in tests/db_search/of1/tnsids.fa -dbtype nucl -out tests/db_search/of1/tnsids \
        -title "translated subjects of the query-anchored rows with Seq-ids" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/tnqid.fa -dbtype nucl -out tests/db_search/of1/tnqid \
        -title "a translated database sequence under the id of a query" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/tnqgi.fa -dbtype nucl -out tests/db_search/of1/tnqgi \
        -title "a translated database sequence whose display id is a query id" -parse_seqids -blastdb_version 4
    rm tests/db_search/of1/*.njs

BLASTX and TBLASTX inputs of the query-anchored rows (rows `db.of1.bx.*` and `db.of1.tx.*`, line RG T3c), in `of1/`,
from the T3c probes (`T/t3/scripts_t3c/gen3c.py`, a seeded generator that places coding pieces in the frames it names):

- `bxq.fa`: one nucleotide query (1229 letters) with pieces in all six frames (+1 with isoleucines coded ATA, which
  genetic code 2 reads as M; +2 with 16 Q that SEG masks; +3, -1, -2, -3); `bxqlc.fa`: the same with two lower-case
  runs (in the +1 piece and in the -2 piece); `bxqminus.fa`: a query with HSPs in one minus frame (-3) only;
  `bxqm.fa`: `bxq`, a query without hits, a query of one frame and `bxqminus`.
- `bxs.fa`: seven proteins: `bx_a` (frames +1 and -2), `bx_ins` and `bx_ins2` (insertions of 5 and 3 residues at one
  query column of frame +2), `bx_rep` (two copies over frame +3), `bx_mins` (an insertion and a deletion on frame -1),
  `bx_m6` (frame -3), `bx_far` (weak); `bxsids.fa`: the same under `gi|...|gb|`, `ref|`, `lcl|`, `gnl|`, `sp|` and
  `gi|...|emb|` ids; `bxslc.fa`: the same with lower-case runs (`-lcase_masking`: shown on every row of the subject);
  `bxqid.fa` (`lcl|Query_1`) and `bxqgi.fa` (`gi|779|lcl|Query_1`): database sequences named as the first query read
  without `-parse_deflines`, which BLASTX refuses after its prolog (`tests/parse_deflines_subject.rs`).
- `txq.fa`: one nucleotide query (571 letters) with pieces in frames +1 (isoleucines coded ATA), +2 (12 Q), -1 and -3;
  `txqlc.fa`: the same with a lower-case run in the -1 piece; `txqm.fa`: `txq`, a query without hits and one of one
  frame.
- `txs.fa`: `tx_pm` (a strong plus copy and a weaker minus copy of the +1 piece: in the query's frame +1 the minus
  HSP comes later and NCBI leaves it out of formats 1 to 4), `tx_mp` (a strong minus copy and two weaker plus copies of
  the +2 piece: the plus HSPs are left out), `tx_c` (the -1 piece on the subject's plus strand and the -3 piece on its
  minus strand), `tx_two` (two plus frames over the +1 piece), `tx_far` (weak); `txsids.fa`: the same under parsed
  ids; `txqid.fa` (`lcl|Query_1`) and `txqgi.fa` (`gi|781|lcl|Query_1`): refused as `bxqid`, `bxqgi`;
  `txself.fa`: one record with repeats on both strands, read as a query and a subject under `-parse_deflines`
  (`tests/parse_deflines_subject.rs`). From the review of T3c (`T/t3/review_t3c/in/t_self*.fa`), rows
  `db.of1.tx.self.*` (line RG T3c fix 1: TBLASTX's sections of a query record that is also a subject, which NCBI's
  `CAlnMix` folds into the query's row): `txselfq.fa` (`lcl|selfx`, 453 letters), `txselfq2.fa` (`lcl|other1`
  before `lcl|selfx`), `txselfs.fa` (`lcl|ydup`, two copies of `selfx`, then `selfx`), `txselfs2.fa` (`lcl|ycopy`,
  the letters of `selfx`, then `selfx`: equal hits) and `txselfs3.fa` (`selfx`, then `ydup`). `txrep2.fa`
  (`lcl|rep2`, two copies of a 300-letter coding piece; `T/t3/fix1_t3c/in/rep2.fa`) and `txrep2b.fa` (`lcl|other1`
  before it): with `-query_loc 1-300` or a second query, and `-subject_loc 301-600`, NCBI's merge throws on its
  sections (`Start iterator already exists for seq2`); `txmg0.fa` (`lcl|mg0`, a tandem record `A+A`) and
  `txmg0s.fa` (`mg0`, then `lcl|mx0`, which holds `A`; review 2 of T3c fix 1, `T/t3/review2_t3c/in/mg0_*.fa`):
  with `-subject_loc 154-303 -seg no` NCBI's `Build` throws (`Infinite loop detected`). Rows `db.of1.tx.throw.*`
  (T3c fix 2: NCBI's report up to the query's description table, then its error, exit 255) and
  `tests/parse_deflines_subject.rs`.
- `bxsr.pal` and `txsr.nal`: aliases of the OIDs 2 to 4 of `bxs` and `txs`.

The databases (version 4; `.pjs` and `.njs` deleted). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/of1/bxs.fa -dbtype prot -out tests/db_search/of1/bxs \
        -title "protein subjects of the BLASTX query-anchored rows" -blastdb_version 4
    makeblastdb -in tests/db_search/of1/bxsids.fa -dbtype prot -out tests/db_search/of1/bxsids \
        -title "protein subjects of the BLASTX query-anchored rows with Seq-ids" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/bxqid.fa -dbtype prot -out tests/db_search/of1/bxqid \
        -title "a protein database sequence under the id of a query" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/bxqgi.fa -dbtype prot -out tests/db_search/of1/bxqgi \
        -title "a protein database sequence whose display id is a query id" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/txs.fa -dbtype nucl -out tests/db_search/of1/txs \
        -title "nucleotide subjects of the TBLASTX query-anchored rows" -blastdb_version 4
    makeblastdb -in tests/db_search/of1/txsids.fa -dbtype nucl -out tests/db_search/of1/txsids \
        -title "nucleotide subjects of the TBLASTX query-anchored rows with Seq-ids" -parse_seqids \
        -blastdb_version 4
    makeblastdb -in tests/db_search/of1/txqid.fa -dbtype nucl -out tests/db_search/of1/txqid \
        -title "a nucleotide database sequence under the id of a TBLASTX query" -parse_seqids -blastdb_version 4
    makeblastdb -in tests/db_search/of1/txqgi.fa -dbtype nucl -out tests/db_search/of1/txqgi \
        -title "a nucleotide database sequence whose display id is a TBLASTX query id" -parse_seqids \
        -blastdb_version 4
    rm tests/db_search/of1/*.pjs tests/db_search/of1/*.njs
## `dbmask/`: database masks of BLASTN searches (line NP, series M, M4a-2a)

Rows `db.dbmask.*` (`-db_soft_mask`, `-db_hard_mask`). The inputs are the probe inputs of the M4a-2 survey
(`ncbi-parity-20261019/M/m4a2/probe/in`, written by its `scripts/gen.py`, seed 20261114): `subj.fna` (five subjects:
`s1` with a 100-letter CA repeat at 1000 and a 24-letter poly-A at 2000, `s2` without low complexity, `s3` with
repeats at both ends, `s4` all `GA`, `s5` random), the queries `q.fna` (all five), `q1.fna` (over the CA repeat of
`s1`), `q4.fna` (over its poly-A), `q5.fna` (over `s5`), `q56.fna` (`q5` and the `s4` query of the edge cases) and the
seqid list `ids_s1.txt`. Made with NCBI BLAST+ 2.17.0 from `dbmask/` (`BLAST_USAGE_REPORT=false`):

```
dustmasker -in subj.fna -infmt fasta -outfmt maskinfo_asn1_bin -out dust.asnb
dustmasker -in subj.fna -infmt fasta -level 40 -window 32 -outfmt maskinfo_asn1_bin -out dust40.asnb
windowmasker -mk_counts -in subj.fna -infmt fasta -out wm.counts
windowmasker -ustat wm.counts -in subj.fna -infmt fasta -outfmt maskinfo_asn1_bin -out wm.asnb
dustmasker -in subj.fna -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out dust_p.asnb
makeblastdb -in subj.fna -dbtype nucl -out a5 -title "mask probe v5" -mask_data dust.asnb,dust40.asnb,wm.asnb
makeblastdb -in subj.fna -dbtype nucl -out a4 -title "mask probe v4" -blastdb_version 4 -mask_data dust.asnb,dust40.asnb,wm.asnb
makeblastdb -in subj.fna -dbtype nucl -out b5 -title "dust40 only" -mask_data dust40.asnb
makeblastdb -in subj.fna -dbtype nucl -out k5 -title "mask ids v5" -mask_data dust.asnb,wm.asnb -mask_id mydust,mywm
makeblastdb -in subj.fna -dbtype nucl -out nm5 -title "no mask v5"
makeblastdb -in subj.fna -dbtype nucl -out p5 -title "parsed v5" -parse_seqids -mask_data dust_p.asnb
makeblastdb -in subj.fna -dbtype nucl -out xb4 -title "big-endian mask data only" -blastdb_version 4 -mask_data dust.asnb
rm xb4.nac
blastdb_aliastool -dblist "a5 nm5" -dbtype nucl -out al_a_nm -title "a5 and nm5"
blastdb_aliastool -dblist "a5 b5" -dbtype nucl -out al_a_b -title "a5 and b5"
blastdb_aliastool -dblist "b5 a5" -dbtype nucl -out al_b_a -title "b5 and a5"
sed 's/^>s/>t/' subj.fna > t.fna
dustmasker -in t.fna -infmt fasta -level 40 -window 32 -parse_seqids -outfmt maskinfo_asn1_bin -out t_dust40.asnb
makeblastdb -in t.fna -dbtype nucl -out zz5 -title "zz dust40" -parse_seqids -mask_data t_dust40.asnb
rm t.fna
```

(`wm.counts` is not kept.) `a5` and `a4` have the algorithms 11 (dust, defaults), 12 (dust `window=32; level=40;
linker=1`) and 30 (windowmasker, default options); `b5` has the second dust as its 11; `k5` has 100 (`mydust`) and 101
(`mywm`); `p5` has the dust of `-parse_seqids` ids as 11; `nm5` has no mask data; `xb4` is a mask column without its
little-endian `.nac`, which LOSAT refuses for a search (NCBI reads its ranges byte-reversed, M/DECISIONS 9); `zz5`
(review 1 of M4a-2a, F2) holds the same subjects as `t1`-`t5` with the second dust alone (its 11), so that in the list
`p5 zz5` the algorithm 11 is `p5`'s dust and `zz5` lacks it: NCBI reads no mask of `zz5` when `-seqidlist`
(`ids_s1s2s3.txt`) searches `p5`'s sequences only, and LOSAT refuses the search only where it would read one. The seqid
lists `ids_s1.txt`, `ids_s1s2s3.txt` and `ids_s3s5.txt` hold those local ids, one per line. The rows
also use `../seqdb/db/mask_nucl_colon_id` (`-mask_id c:id`, stored with 0x01 for the colon),
`../seqdb/db/mask_nucl_mv` (two volumes, dust 11 and windowmasker 30) and `../fasta/outfmt11/arch_empty.fna`
(an empty query).

The edge databases hold hand-made masks of `s5` (local id 5, 1200 letters) and `s4` (local id 4, 300 letters) in
ASN.1 text (`edge_<case>.asnt`, the survey's `scripts/mkedge.py` form: `algo-program dust`, `algo-options "edge
<case>"`, one `packed-int` per sequence, `from`/`to` inclusive): `unsorted` (600-699, 100-199), `overlap` (100-299,
200-399), `adjacent` (100-199, 200-299), `ends` (0-99, 1100-1199), `tiny` (301-301, 500-502, 703-705), `nested`
(100-499, 200-299), `many` (25 masks of 5 letters every 37 from 150), `head` (0-99), `tail` (1000-1199), `whole5`
(0-1199) and `whole4` (`s4` 0-299, `s5` 400-449). makeblastdb stores each sequence's masks in the order given, so the
searches meet unsorted and overlapping lists:

```
makeblastdb -in subj.fna -dbtype nucl -out edge1 -title "edge masks 1" -blastdb_version 4 \
  -mask_data edge_unsorted.asnt,edge_overlap.asnt,edge_adjacent.asnt,edge_ends.asnt,edge_tiny.asnt
makeblastdb -in subj.fna -dbtype nucl -out edge2 -title "edge masks 2" -blastdb_version 4 \
  -mask_data edge_nested.asnt,edge_many.asnt,edge_head.asnt,edge_tail.asnt,edge_whole5.asnt,edge_whole4.asnt
```

`edge1` has the algorithms 11 to 15 in that order, `edge2` 11 to 16.

Review 2 of M4a-2a added damaged copies of `a5` and alias files over `p5 zz5`:

```
mkdir dmg dmg3 && cp a5.* dmg/ && cp a5.* dmg3/
truncate -s 360 dmg/a5.nac dmg/a5.nab
for f in dmg3/a5.nac dmg3/a5.nab; do printf '\x00\x00\x10\x00' | dd of=$f bs=1 seek=8 conv=notrunc; done
```

`dmg/a5` has mask data files cut inside the blob of `s1` (a range past the end of the data file: NCBI's
`SEQDB_FILE_ASSERT`), `dmg3/a5` a first blob whose algorithm 11 states 4096 pairs (`CBlastDbBlob::x_ReadRaw: hit end
of data`); NCBI's search reports both as `BLAST engine error: Database memory map file error`, exit status 3. The alias
files `zr_t2.nal`, `zr_p5.nal` and `zr_t1.nal` give `DBLIST p5 zz5` with the OID ranges 7-7 (`t2` alone, no mask
data), 3-5 (`s3`-`s5`) and 1-6 (`s1`-`s5` and `t1`, which has mask data that `zz5`'s volume reads without the
algorithm 11 of the list): NCBI reads masks of `zz5` only in the last.

M4a-2b (`-db_hard_mask`, line NP) added `edge3` (s5's one-letter mask 1001, `edge_short.asnt`: a last mask shorter
than the two letters before its end that the last chunk starts with, which LOSAT refuses) and four alias files that
leave out the sequences whose mask runs to their end (`a5`'s and `k5`'s `s4`, `mask_nucl_mv`'s `m04`, OID 3): `zh_a5.nal`
and `zh_k5.nal` (OIDs 1-3, `s1`-`s3`), `zh_mv_a.nal` (`../../seqdb/db/mask_nucl_mv`, OIDs 1-3: `m01`-`m03`) and
`zh_mv_b.nal` (OIDs 5-8: `m05`-`m08`, the second volume):

```
makeblastdb -in subj.fna -dbtype nucl -out edge3 -title "edge masks 3" -blastdb_version 4 -mask_data edge_short.asnt
```

`edge_short.asnt` is `edge_tail.asnt` with `algo-options "edge short"` and `from 1001, to 1001`.

The review fixes of M4a-2b (ad-M4a-6) added `edge4` and `edgew`. `edge4` holds two masks of `s5` with an empty mask
at 0 first (`from 0, to 4294967295`, which makeblastdb stores as `[0, 0)`): algorithm 11 (`edge_e0tail.asnt`, then
1000-1199: the chunk after the empty mask starts at 0 and the next mask runs to the end, refused; review F1) and 12
(`edge_e0mid.asnt`, then 500-599: searched). `edgew` holds four sequences of 3 to 6 letters (`wb.fna`, `w3`-`w6`) whose
mask runs to their end (`edge_wb.asnt`: 1-2, 2-3, 3-4 and 3-5), searched alone through `zw_s3.nal`, `zw_s4.nal` and
`zw_s5.nal` (OIDs 1, 2 and 3) with the query `qw.fna` at words 4, 5 and 7 (review F2: NCBI scans every BLASTN
subject, and its scan of the whole sequence reads a word when the sequence has at least the word's letters):

```
makeblastdb -in subj.fna -dbtype nucl -out edge4 -title "edge masks 4" -blastdb_version 4 -mask_data edge_e0tail.asnt,edge_e0mid.asnt
makeblastdb -in wb.fna -dbtype nucl -out edgew -title "edge masks word" -blastdb_version 4 -mask_data edge_wb.asnt
```

The second review fixes of M4a-2b (ad-M4a-6, 2026-11-15) added `edge5`: last masks shorter than the letters before their
end that the final chunk starts with, refused only where NCBI's merge after that chunk can take an alignment of an earlier
chunk by the overlap it reads past its list of ranges (review 2 F1). Algorithm 11 (`edge_ov.asnt`: `s5`'s one-letter
masks 300 and 302; the chunk between them scans no word) and 12 (`edge_em.asnt`: 700-701, then a mask of no letters at
702, `from 702, to 701`) are searched; 13 (`edge_ovr.asnt`: 301 and 302; the chunk 0-301 ends after the final chunk's
offset 300) is refused; 14 (`edge_fw.asnt`: `s4`'s mask 289; the final chunk 288-300 is scanned from 290) is searched at
words of 11 and 28 letters and refused at blastn-short's 7:

```
makeblastdb -in subj.fna -dbtype nucl -out edge5 -title "edge masks 5" -blastdb_version 4 -mask_data edge_ov.asnt,edge_em.asnt,edge_ovr.asnt,edge_fw.asnt
```

M4a-2c (line NP, 2026-11-15; PD-LOSAT-NCBI-DEFECTS exception 13: with `-db_hard_mask`, a sequence whose first mask
runs to its end is searched before the mask at a scanned word of 11 or more) added four databases. `edge6` holds the
near inputs of the rows on which NCBI BLAST+ 2.17.0 crashes: algorithm 11 (`edge_whole4n.asnt`) is `edge2`'s 16 with
`s4`'s mask one letter shorter (0-298), 12 (`edge_e0tailn.asnt`) is `edge4`'s 11 with the mask after the empty one one
letter shorter (1000-1198). `mvn` holds `../../seqdb/mask_nucl.fna` with `mask_nucl_mv`'s dust (algorithm 11, from
`../../seqdb/expected/mask_nucl_mv.masks.tsv`) and `m04`'s mask one letter shorter (`mvn_dust.asnt`: 300-388, one
volume; a search's numbers and rows are those of the two volumes). `edge7` (`tl.fna`, queries `q7.fna`) holds four
sequences of the M4a-2c probes with tail masks (`edge_tl.asnt`, algorithm 11): `x1`-`x3` of 916, 326 and 788 letters
whose last 1, 2 and 3 letters are masked (queries that match to the end: the traceback reads the whole mask as `B`, so
no alignment reaches the last letter, while in the near input of `edge_tln.asnt`, algorithm 12, each mask one letter
shorter and `x1` unmasked, the alignments reach it) and `x4`, a short subject (74 letters, mask 55-73). `edge8`
(`ts.fna`, queries `q8.fna`) holds `y1` (37 letters, mask 12-36) and `y2` (1877 letters, mask 1682-1876) whose
queries hold only 12 and 17 letters before the masks, where NCBI's scan past the chunk finds words across the masks'
start (`edge_ts.asnt`, algorithm 11) that its normal chunks of the near input (`edge_tsn.asnt`, 12) do not: the rows
of `edge8`'s 11 expect the near input's output (`ncbi-near:` in the manifest). The masks of `edge7` and `edge8` name
the sequences by number (`local id N`), as the other edge files do:

```
makeblastdb -in subj.fna -dbtype nucl -out edge6 -title "edge masks 6" -blastdb_version 4 -mask_data edge_whole4n.asnt,edge_e0tailn.asnt
makeblastdb -in tl.fna -dbtype nucl -out edge7 -title "edge masks 7" -blastdb_version 4 -mask_data edge_tl.asnt,edge_tln.asnt
makeblastdb -in ts.fna -dbtype nucl -out edge8 -title "edge masks 8" -blastdb_version 4 -mask_data edge_ts.asnt,edge_tsn.asnt
makeblastdb -in ../../seqdb/mask_nucl.fna -dbtype nucl -out mvn -title "mask_nucl_mv dust near" -blastdb_version 4 -mask_data mvn_dust.asnt
```

The O1 follow-up of M4a-2c (2026-11-15, `M/m4a2c/O1_SURVEY.md` of line NP) added `edge9` and `edge10`, from the survey's
probes `co` and `rp`. `edge9` (`tc.fna`, query `q9.fna`) holds `c1` (620 letters, mask 20-619 in `edge_tc.asnt`, algorithm
11; 20-618 in `edge_tcn.asnt`, 12) and then `c2` (no mask), which holds the query's first 300 letters on the diagonal of
`c1`'s masked copy: NCBI's search of `c1`'s tail drops `c2`'s row in the following subject (always at one thread, in some
runs at two and four). `edge10` (`tr.fna`, query `q10.fna`) holds `r1` (900 letters) with two copies of the query inside the
mask 20-899 (`edge_tr.asnt`, 11; 20-898 in `edge_trn.asnt`, 12), on which NCBI's interval tree runs out of memory and NCBI
crashes at words of 11 and more. The rows of their algorithm 11 expect NCBI's output of algorithm 12 (`ncbi-near:`):

```
makeblastdb -in tc.fna -dbtype nucl -out edge9 -title "edge masks 9" -blastdb_version 4 -mask_data edge_tc.asnt,edge_tcn.asnt
makeblastdb -in tr.fna -dbtype nucl -out edge10 -title "edge masks 10" -blastdb_version 4 -mask_data edge_tr.asnt,edge_trn.asnt
```

## `dbmaskp/`: database masks of BLASTP and BLASTX searches (line NP, series M, M4b)

Rows `db.pmask.*` (`blastp.pmask.tsv`, `blastx.pmask.tsv`: `-db_soft_mask`, `-db_hard_mask`). The inputs are the probe
inputs of the M4b survey (`ncbi-parity-20261019/M/m4b/probe/in`, written by its `scripts/gen.py`, seed 20261120):
`prot.faa` (seven proteins: `p1` with a Q block and an SG block, `p2` random, `p3` with an EK head and a K tail, `p4`
all GS, `p5` random for the hand masks, `p6` a homolog of `p1`, `p7` with a PA block), `proth.faa` (the same without
`p4`, whose seg mask covers it whole), the protein queries `qp.faa` (`qp_all.fa`) and the nucleotide queries `qx.fna`
(`bx_all.fa`), and the hand masks of `p5` `edge_*.asnt` (`pedge_*.asnt`, ASN.1 text, `local id 5`, algorithm seg).
Made with NCBI BLAST+ 2.17.0 from `dbmaskp/` (`BLAST_USAGE_REPORT=false`; the script is
`ncbi-parity-20261019/M/m4b/ws/fixtures.sh`):

```
segmasker -in prot.faa -infmt fasta -outfmt maskinfo_asn1_bin -out seg.asnb
segmasker -in prot.faa -infmt fasta -outfmt maskinfo_asn1_bin -window 8 -locut 1.8 -hicut 2.1 -out seg8.asnb
segmasker -in proth.faa -infmt fasta -outfmt maskinfo_asn1_bin -out segh.asnb
segmasker -in proth.faa -infmt fasta -outfmt maskinfo_asn1_bin -window 8 -locut 1.8 -hicut 2.1 -out segh8.asnb
makeblastdb -in prot.faa -dbtype prot -out ps5 -title "seg masks v5" -mask_data seg.asnb,seg8.asnb
makeblastdb -in proth.faa -dbtype prot -out ph5 -title "seg masks v5 without p4" -mask_data segh.asnb,segh8.asnb
makeblastdb -in proth.faa -dbtype prot -out ph4 -title "seg masks v4 without p4" -blastdb_version 4 -mask_data segh.asnb,segh8.asnb
makeblastdb -in proth.faa -dbtype prot -out phmv -title "seg masks, volumes" -max_file_sz 600B -mask_data segh.asnb
for c in mid head adjacent gap2 tiny neartail two_tail tail unsorted; do
  makeblastdb -in prot.faa -dbtype prot -out pe_$c -title "edge $c" -mask_data edge_$c.asnt
done
makeblastdb -in prot.faa -dbtype prot -parse_seqids -out pl3a -title "parsed ids, p5 masked in the middle" -mask_data l3_mid.asnt
makeblastdb -in prot.faa -dbtype prot -parse_seqids -out pl3b -title "parsed ids, p5 masked at the head" -mask_data l3_head.asnt
blastdb_aliastool -dblist "ph5 phmv" -dbtype prot -out pal -title "alias of ph5 and phmv"
python3 ../make_dbmaskp_pbh.py --ncbi-only WORK $NCBI_BIN .   # pbh.faa, edge_pbh.asnt, qbh.fna, qbh.faa, pbh
makeblastdb -in pshort.faa -dbtype prot -out pshort -title "short pieces between masks" -mask_data edge_short.asnt
python3 - <<'PY'   # prot_p4x.faa, prot_p5x.faa: prot.faa with p4[119], p5[399] stored as X (the near inputs' freed letters)
def put(src, dst, name, pos):
    recs, cur = [], None
    for line in open(src):
        if line.startswith('>'):
            cur = [line, '']; recs.append(cur)
        else:
            cur[1] += line.strip()
    with open(dst, 'w') as f:
        for head, seq in recs:
            if head[1:].split()[0] == name:
                seq = seq[:pos] + 'X' + seq[pos + 1:]
            f.write(head + '\n'.join(seq[i:i + 60] for i in range(0, len(seq), 60)) + '\n')
put('prot.faa', 'prot_p4x.faa', 'p4', 119)
put('prot.faa', 'prot_p5x.faa', 'p5', 399)
PY
makeblastdb -in prot_p4x.faa -dbtype prot -out ps5n -title "seg masks v5, p4 one shorter" -mask_data near_ps5.asnt
makeblastdb -in prot_p5x.faa -dbtype prot -out pe_tailn -title "edge tail, one shorter" -mask_data edge_tailn.asnt
makeblastdb -in prot.faa -dbtype prot -out pe_e0tail -title "edge e0tail" -mask_data edge_e0tail.asnt
makeblastdb -in prot_p5x.faa -dbtype prot -out pe_e0tailn -title "edge e0tail, one shorter" -mask_data edge_e0tailn.asnt
rm -f *.pjs
```

The review of PR 1 added `pl3a` and `pl3b` (parsed Seq-ids; `l3_mid.asnt` and `l3_head.asnt` are `edge_mid.asnt` and
`edge_head.asnt` with `local str "p5"` and one algorithm text, so that the two databases share the Seq-id `p5` under
different masks), the alias file `pal`, and `pbh` with the blastx query `qbh.fna` and the blastp query `qbh.faa` (the
protein that `qbh.fna` back-translates): subjects that hold a weaker copy of a query stretch, a still weaker one starting
a little before it, a hard-masked stretch, and a close copy of the stretch, on which `-subject_besthit` after every
chunk (NCBI) and once after the merge differ. `make_dbmaskp_pbh.py` (seed 20261121) draws rounds of 400 such subjects.
With `--ncbi-only` (resume 2 of PR 1, run with NCBI 2.17.0 alone) it keeps the first three of a round where NCBI's rows
with `-subject_besthit` hold fewer alignments than its rows without the option filtered once, and the first three where
blastp's do (round 0: blastx s1, s2, s4, s6 and blastp s4, s6 of the six differ so); with OLD and NEW (LOSAT before
and after the fix) it keeps the first three where those two builds differ.

`ps5`, `ph5` and `ph4` have the algorithms 21 (seg, defaults) and 22 (seg `-window 8 -locut 1.8 -hicut 2.1`); `phmv`
(several volumes) and the `pe_*` databases have 21 only. `ps5` keeps `p4`, masked whole, and `pe_tail` has a masked
tail: with `-db_hard_mask` both are searched over the part before the mask (the Owner's answer Q1 = B of M4b,
PD-LOSAT-NCBI-DEFECTS exception 21, M4b-3), and their rows hold NCBI's run of the near input where NCBI's own run
crashes or differs from it: `ps5n` (`near_ps5.asnt`: `ps5`'s algorithm 21 written out by `segmasker -outfmt interval`,
with `p4`'s mask one letter shorter, [0, 119)) and `pe_tailn` (`edge_tailn.asnt`: [300, 399)). `pe_e0tail`
(`edge_e0tail.asnt`: an empty mask at 0, `from 0, to 4294967295`, which makeblastdb stores as [0, 0), then [300, 400))
is searched the same way, as BLASTN's `edge4`; its near input is `pe_e0tailn` ([0, 0) and [300, 399)). The near
inputs store the letter that the shorter mask frees as `X` (`prot_p4x.faa`, `prot_p5x.faa`), so that it stays masked in
the alignments and the composition-based statistics, as LOSAT reads the masked letters. `ph5`, `ph4` and
`phmv` are the other hard-mask databases; `pe_unsorted` (masks out of order) is refused.

`pshort` (resume 2 of PR 1, 2026-11-21; `pshort.faa`, `edge_short.asnt`, `qshort.faa`, `qshort.fna`, written by
`ncbi-parity-20261019/M/m4b/probe3/make.py`, seed 20261121) holds four copies of one random protein without `W` but for
a run of six `W` at 60-65, under hand masks (algorithm 21, local ids 1 to 4): `s1` [30, 60) and [62, 90), which leave the
piece `WW` before a masked `W`; `s2` [30, 60) and [63, 90), the piece `WWW`; `s3` [30, 60) and [60, 90), touching masks;
`s4` [1, 10), a first piece of one letter. The query `qshort.faa` (and its back-translation `qshort.fna`) holds one
`WWW`. NCBI's one-hit scan (`-window_size 0`) aligns that `WWW` with `s2`'s piece only: it reads no word of a hard piece
shorter than the word (aa_ungapped.c:762-766, no raise of the range's end), where a word read at the piece's start
would have given `s1` the same alignment; with `-db_soft_mask` it reads no word of `s4` at all, whose first scan range
is shorter than the word. The rows `*.w0.*` of `pe_tiny`, `pe_gap2` and `pe_adjacent` search the hand masks of `p5`
with pieces of two and no letters the same way. The rows `*.ungapped.*` of BLASTP (restored after dev's BLASTP
`-ungapped`, #1420, and `ps5`/`ph5` with seg masks) search each hard piece's ungapped hits and merge them as the
gapped pieces are merged.

`ov_*` and `ov_atha` (M4b-3b, 2026-11-22; overlapping hard masks, rows `db.pmask.{bp,bx}.hard.ov_*`): databases whose masks
overlap, so that `-db_hard_mask` searches an inverted range (a chunk of negative length: a mask that starts before the
previous one ends). `../../make_dbmaskp_ov.py` (seeds 20261121 and 20261122, the letters and masks of the probes of
`ncbi-parity-20261019/M/m4b3/probe/make_p1.py` and `M/m4b3b/probe/make_p2.py`) writes `ov_<shape>.faa` and `ov_<shape>.asnt`
(ASN.1 text, algorithm seg, `local id N`, N the order of the subject) and the queries `qov1.faa`/`qov1.fna` (probe 1: `e3`,
`e3c`, `g3`, `nest`, `big`) and `qov2.faa`/`qov2.fna` (probe 2: `s2`, `tl`, `r37`, `r36`, `nb`, `lk`; a random protein with one `WWW`
and its back-translation for BLASTX); each shape is one subject of 240 random letters without `W` but for the `W` runs named
in the script's header (`ov_lk`: two subjects), and the database is made as the `pe_*`:

```
python3 ../../make_dbmaskp_ov.py --ncbi-bin /home/kawato/micromamba/bin   # = for each shape X:
makeblastdb -in ov_X.faa -dbtype prot -out ov_X -title "overlapping masks X" -mask_data ov_X.asnt && rm -f ov_X.pjs
python3 ../../make_dbmaskp_ovatha.py --atha ATHA_FAA --ncbi-bin /home/kawato/micromamba/bin   # ov_atha, qova.faa, qova.fna
```

`ov_e3` (masks [40,60) [55,70), `WWW` at 60 and 70: NCBI's two-hit scan pairs the word at 70 with the stale diagonal entry of
the word at 60) and its control `ov_e3c` (touching masks [60,70), the same letters), `ov_g3` (two consecutive inverted chunks;
NCBI's gapped two-hit BLASTP search of it never ends: the seed of the chunk of length -5 goes into an interval tree whose
subject bounds that length sets, which splits its nodes without end, so those runs are rejection rows), `ov_nest` (a mask
inside the previous one), `ov_big` (overlap 60 > the window 40), `ov_s2` (an inverted chunk, then a chunk of 2 letters shorter
than the word: the seed of its one word, read from the stored letters past the chunk, extended with the chunk's length; NCBI
writes the same rows in every run), `ov_tl` (`ov_s2`'s letters with the masks [65,240) then [55,64): a mask that runs to the
end after the chunk from the start, with a mask stored after it, which LOSAT refuses before the search), `ov_wpe` (A, of 100
letters with masks [10,60) [20,99) [30,40), whose one word at 99 covers A's sentinel and B's first letter, then B; searched
for `qov2.faa` and for `qwpe.faa`, `qov2.faa`'s protein with `W-W`, whose `-` the FASTA reader drops), `ov_r37` and
`ov_r36` (overlaps 37 = window - word and 36), `ov_nb` (consecutive inverted chunks with different words) and `ov_lk` (two
subjects, A with masks [40,200) [60,238): diagonal entries that outlive A reach B, 98 positions past A's end, 137 in the
one-hit scan (`-window_size 0`); NCBI's rows of B depend on the thread order). LOSAT searches B with A's entries out of reach
(PD-LOSAT-NCBI-DEFECTS exception 29, the Owner's answer Q5 = B of M4b, 2026-11-22; it refused `ov_lk` before). At the default
NCBI's rows are those of the near input `ov_lkn`; the rows at `-evalue 1000` and `-ungapped` expect NCBI's run of `ov_lkn`
(`ncbi-near:`, rows `.hard.near.ov_lkn.*`): `make_dbmaskp_ov.py` (`near_db`) writes `ov_lkn.faa` (A, a spacer of 240 `X`, which
no word of these queries hits, as they hold no `X`, then B; at least as long as the walk's margin less the window, `spacer_length`, here the database's average
length) and `ov_lkn.asnt` (A's and B's masks, `local id` renumbered), and the alias `ov_lkn.pal` over the volume `ov_lknv` keeps
`ov_lk`'s NSEQ, LENGTH, STATS_NSEQ and STATS_TOTLEN, so that NCBI's statistics are `ov_lk`'s. At `-window_size 0` NCBI's
one-thread search starts B 98 below 0, on the entries that no sequence wrote, as LOSAT's does (fix round 2 of M4b-3b,
DECISIONS 55 of series M), so those rows expect NCBI's own bytes; the near input `ov_lkx` (`blank_db`: A's letters all `X`, the
lengths and masks the same, so NCBI's offsets too) gives B's rows the same (rows `.hard.near.ov_lkx.21.w0.6`):

```
makeblastdb -in ov_lkn.faa -dbtype prot -out ov_lknv -title "overlapping masks lk" -mask_data ov_lkn.asnt && rm -f ov_lknv.pjs
printf 'TITLE overlapping masks lk\nDBLIST ov_lknv\nNSEQ 2\nLENGTH 480\nSTATS_NSEQ 2\nSTATS_TOTLEN 480\n' > ov_lkn.pal
makeblastdb -in ov_lkx.faa -dbtype prot -out ov_lkx -title "overlapping masks lk" -mask_data ov_lkx.asnt && rm -f ov_lkx.pjs
```

`ov_zt` and `ov_ztl` (fix round 2 of M4b-3b, review 2's probe `M/m4b3b/review/probe_zt`, seed 20261123; rows
`db.pmask.{bp,bx}.hard.ov_zt*`): A (`ov_zt`: 5 `X`; `ov_ztl`: 100), which stores no entry, then B (240 letters, masks [10,100)
[20,110): chunks 10, -80 and 130, a one-hit walk 70 below its start); the query `qzt.faa` (`qzt.fna` its back-translation)
holds B's letters 110-121 at 18. NCBI's one-thread search starts B at 5 in `ov_zt`, where its walk meets the entries that no
sequence wrote, and at 100 in `ov_ztl`, where it does not; `ov_zto.pal` is an alias over `ov_zt` of its OID 2 alone (FIRST_OID 2,
LAST_OID 2), where B is the first searched sequence. Made as the other shapes (`make_dbmaskp_ov.py`, the alias by its
`main`).

`ov_th10`, `ov_thm20`, `ov_th10f` and `ov_thm20f` (fix round 3 of M4b-3b, review 3's probe `M/m4b3b/review/probe_r3`,
`probe5` of `make_dbmaskp_ov.py`, no draw of its own; rows `db.pmask.{bp,bx}.hard.ov_th*`): A (430 `X`, or 400 in the `thm20`
forms; masks [10,300) [20,310)), which stores no entry and whose two-hit walk goes back, then B (`ov_zt`'s B with the query's
letters 18-29 at 3, or in the `f` forms its letters 20-24, one run of three word hits; no masks). NCBI's one-thread search
starts B on a new table at 10 (`thm20`: -20), within the window of the entries that no sequence wrote (below them at -20).
The rows of `qzt.faa` / `qzt.fna` (one query) are NCBI's one-thread bytes, and LOSAT gives them at every thread count
(NCBI's own runs at 4 threads can differ: `ov_thm20f` gave 4 rows instead of 10); `q1q.faa` / `q1q.fna` (q1, 60 letters,
then `qzt`'s query) are two NCBI query batches at `BATCH_SIZE=50` and one by default, and q's rows differ between the two in
NCBI's runs (BLASTP: 5 and 4 rows of q in `ov_th10`), as NCBI searches each batch alone, its queries counted from the batch's
start: the rows `.q1q.` (and `.q1q.b50.`) check both batchings, and `LOSAT/tests/blastp_prelim_blocks.rs` runs them with
LOSAT_X_BPQBLOCK 0 and 1 too.

The other batch rows of BLASTP (fix round 4 of M4b-3b, review 4's L2, `probe6` of `make_dbmaskp_ov.py`): `q1q.faa` on
`ov_s2` and `ov_nb` (sequences kept apart because their walks read back; q's rows 1 and 3, 2 and 4 in one batch and in two),
`ov_zt` and `ov_ztl` with `-window_size 0` (the one-hit rings; q1's rows 15 and 14, and 22 in two batches of `ov_ztl`);
`q4.faa` (q1, q, q1b, qb: three batches at `BATCH_SIZE=50`, of 60, 108 and 48 letters) on `ov_th10`, `ov_s2` and `ov_zt`
(`-window_size 0`), which LOSAT_X_BPQBLOCK=120 searches in two blocks, the first of two batches; and `qbig10000.faa` (a
random query of 10,000 letters without `W`, seed 20261125, then q) on `ov_th10`, two batches at BLASTP's default batch
size (q's rows 5, as in a batch of its own). All are NCBI's one-thread bytes; `blastp_prelim_blocks.rs` runs every
`.q1q.`, `.q4.` and `.qbig10000.` row with LOSAT_X_BPQBLOCK unset, 0, 1 and 120 at 1 and 4 threads and reads from
LOSAT_X_REACH that the blocks ran exactly where an arm plans two or more.

The batch rows of BLASTX (line NP item M4b-3b2, decision 57 of series M; `qfc.fna` of `make_dbmaskp_ov.py`): 10,002 `A`,
whose frames SEG masks whole (no word), then the query `cand` (75 letters, `M/m4b3b2/scripts/find_fx.py` candidate 1 of
seed 1, back-translated with seed 20261125), on `ov_th10`: one NCBI query batch at `BATCH_SIZE=30000` (rows
`db.pmask.bx.hard.ov_th10.qfc.b30000.*`: 18 rows of `cand`, 20 with `-ungapped`), two at the default batch size
(`.qfc.21.e1000.6`: 19, `cand` first in its batch; the `.err` file holds NCBI's warning for the filler), and the search
without masks (`db.pmask.bx.none.ov_th10.qfc.b30000.*`: 5 rows). A word of `cand`'s first frame near
B's start reads the entries of B's new table at `cand`'s offset in the batch, so BLASTX, which searches a batch in units
of 10,002 letters, searches the whole batch as one unit where the walk keeps a sequence apart
(`blastx::query_setup::unit_ranges`); the units wrote 19 and 21 rows at `BATCH_SIZE=30000` before.

`ov_m1` (fix round 4 of M4b-3b, review 4's L1 and its database `M/m4b3b/review/probe_r4/m1`, `probe6`, no draw of its
own; rows `db.pmask.{bp,bx}.hard.ov_m1.*`, `qzt.faa` / `qzt.fna`): 600 subjects `s0` to `s599`, each one `X` but `s62`
(5,000 `X`, masks [10,4990) [20,4995): its walk goes about 4,955 back), `s63` (P, 50 `X` without masks, the last of the
first 64 subjects, which NCBI's one-thread search starts on a new table below the entries that no sequence wrote), `s64`
(5,000 `X`) and `ov_zt`'s B without masks at `s128`, `s192`, ..., `s512`, the first subject of each later group of 64.
BLASTX at two threads or more in waves of 64 subjects per worker (LOSAT_X_BXSUBJFLOW=0) or in tasks of 64 subjects
(LOSAT_X_BXTASKRES=0) searches one of those B right after P on P's table, so its rows are NCBI's only where P's search
moves the offset past those entries; `LOSAT/tests/blastx_subject_flow.rs` runs the two BLASTX rows so at 2 and 4
threads (before that move, CI's build of 5fcc98420 lost the rows of one to four B there). Made as the other shapes
(`makeblastdb -in ov_m1.faa -dbtype prot -out ov_m1 -title "overlapping masks m1" -mask_data ov_m1.asnt`).

6 of `ov_atha`'s 40 proteins (AT1G02550.1, AT1G02800.1, AT1G04080.2, AT1G04080.5, AT1G06760.1, AT1G09460.1) scan below
their starts or leave diagonal entries past the next subject's start in the one-hit scan, which LOSAT refused before
exception 29 (decision 52); NCBI's rows of `ov_atha` at `-window_size 0` equal those of its near input with a spacer of `X`
at each of the six boundaries (at 1 and 4 threads, BLASTP and BLASTX; `ncbi-parity-20261019/M/m4b3b/Q5B_PROGRESS.md`), so
they are `ncbi` rows. `ov_athaw` (`make_dbmaskp_ovatha.py` too: the 34 other proteins, `segmasker` and `makeblastdb` as for
`ov_atha`) has the one-hit rows `.ov_athaw.21.w0.6`. NCBI BLAST+ 2.17.0 gave the same stdout bytes at 1 thread twice and at 4 threads for every
frozen row. `ov_atha` holds 40 proteins of the A. thaliana TAIR10 set (150-600 letters) whose
`segmasker` masks (`-outfmt maskinfo_asn1_bin`, algorithm 21) overlap: the first ten by name with two or more consecutive
overlaps, the first five of the rest with a nested mask and then the first with an overlap up to 40, none with a mask that
reaches the last letter (`make_dbmaskp_ovatha.py` has the rule); `qova.faa` is the four with the most overlaps and `qova.fna`
their back-translation (seed 20261122). The rows run the default two-hit search, `-evalue 1000`, `-window_size 0` and
`-ungapped -comp_based_stats 0 -evalue 1000`, outfmt 6 with the long field list, and `ov_atha` also outfmt 0 at the default.

## `dbmaskt/`: database masks of TBLASTN and TBLASTX searches (line NP, series M, M4b PR 2)

Rows `db.tmask.*` (`tblastn.tmask.tsv`, `tblastx.tmask.tsv`: `-db_soft_mask`, `-db_hard_mask`). The inputs are written by
`ncbi-parity-20261019/M/m4b/scripts2/gen_fix2.py` (seed 20261121): `nucl.fna` (four subjects: `t1`, 1,500 bases with a
60-base CA repeat at 600 and a 30-base poly-A at 1100; `t2`, `t1`'s bases 300-1199 with 8 % substitutions; `t3`, 6,000
bases with an 80-base CA repeat at 2500 and a 50-base GT repeat at 4800, long enough for the traceback's partial fetch;
`t4`, 1,200 random bases for the hand masks), `nuclf.fna` (the same and `t5`, 300 bases of `GA`, which dust masks whole),
the TBLASTN queries `qn.faa` (translations, with 10 % substitutions, of stretches that cross the masks in plus and minus
frames) and the TBLASTX queries `qx.fna` (stretches with 8 % substitutions, both strands), and the hand masks of `t4`
(`hedge_<case>.asnt`, `local id 4`, from/to inclusive): `mid0` (300-389), `mid1` (301-390), `mid2` (302-391), `adjacent`
(600-699, 700-799), `overlap` (600-749, 700-799), `mid` (300-389, 900-959), `unsorted` (600-699, 300-389), `head`
(0-99), `tail` (1000-1199) and `near` (1000-1197). Made with NCBI BLAST+ 2.17.0 from `dbmaskt/` (`BLAST_USAGE_REPORT=false`;
the script is `ncbi-parity-20261019/M/m4b/ws2/fixtures.sh`):

```
dustmasker -in nucl.fna -infmt fasta -outfmt maskinfo_asn1_bin -out dust.asnb
dustmasker -in nucl.fna -infmt fasta -level 10 -outfmt maskinfo_asn1_bin -out dust10.asnb
dustmasker -in nucl.fna -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out dust_p.asnb
dustmasker -in nuclf.fna -infmt fasta -outfmt maskinfo_asn1_bin -out dustf.asnb
makeblastdb -in nucl.fna -dbtype nucl -out tm5 -title "dust masks v5" -mask_data dust.asnb,dust10.asnb
makeblastdb -in nucl.fna -dbtype nucl -out tb5 -title "dust level 10 only" -mask_data dust10.asnb
makeblastdb -in nucl.fna -dbtype nucl -out tm4 -title "dust masks v4" -blastdb_version 4 -mask_data dust.asnb,dust10.asnb
makeblastdb -in nucl.fna -dbtype nucl -out tp5 -title "dust masks, parsed ids" -parse_seqids -mask_data dust_p.asnb
makeblastdb -in nucl.fna -dbtype nucl -out tmv -title "dust masks, volumes" -blastdb_version 4 -max_file_sz 1600B -mask_data dust.asnb
makeblastdb -in nuclf.fna -dbtype nucl -out tf4 -title "dust masks with t5" -blastdb_version 4 -mask_data dustf.asnb
blastdb_aliastool -dblist "tm5 tb5" -dbtype nucl -out tal -title "tm5 and tb5"
makeblastdb -in nucl.fna -dbtype nucl -out th1 -title "hand masks 1" -blastdb_version 4 \
  -mask_data hedge_mid0.asnt,hedge_mid1.asnt,hedge_mid2.asnt,hedge_adjacent.asnt,hedge_overlap.asnt,hedge_mid.asnt
makeblastdb -in nucl.fna -dbtype nucl -out th2 -title "hand masks 2" -blastdb_version 4 \
  -mask_data hedge_unsorted.asnt,hedge_head.asnt,hedge_tail.asnt,hedge_near.asnt
sed '164s/....$/NNNN/' nucl.fna > nucl_t4n.fna   # t4's bases 1196-1199, which hedge_tailn frees, stored as N
makeblastdb -in nucl_t4n.fna -dbtype nucl -out th2tn -title "hand mask tail, 4 bases shorter" -blastdb_version 4 \
  -mask_data hedge_tailn.asnt
rm -f *.njs
```

`tm5`, `tm4` and `tp5` have the algorithms 11 (dust, defaults) and 12 (dust `level=10`; `tp5` has 11 only); `tb5` has
the second dust alone, as its 11, which the alias `tal` and the list `tm5 tb5` number 12 as `tm5` does; `tmv` (three
volumes: `t1` and `t2`, `t3`, `t4`) and `tf4` have 11. `th1` has the hand masks `mid0`, `mid1`, `mid2`, `adjacent`,
`overlap` and `mid` as 11 to 16, `th2` `unsorted`, `head`, `tail` and `near` as 11 to 14. A mask is mapped to the six
reading frames by its positions divided by 3 (the minus frames from frame -1's length); with `-db_hard_mask` a
sequence whose masks leave NCBI one part of a frame to search, from the frame's start, that stops before the frame's
end (`tf4`'s `t5`, `head` for frame -1, `tail` for frame +1, `near` for frames +2 and +3) is searched over that part
alone (the Owner's answers Q1 = B and Q2 = B of M4b, PD-LOSAT-NCBI-DEFECTS exception 21, M4b-3). TBLASTX's rows of
`tail` and `near`, where NCBI writes an HSP once per seed past the part, hold NCBI's run of the near input `th2tn`
(`hedge_tailn.asnt`: `t4` [1000, 1196), the fewest bases short of the end after which frames +1 to +3 each keep a
letter after the mask; `nucl_t4n.fna` stores the 4 freed bases as `N`, so that they stay masked as LOSAT reads them);
the other rows are NCBI's own, which equal its run of the near input there. Masks out of
order (`unsorted`, `overlap`) are refused, as for the other programs.

The review rows of PR 2 (2026-11-21; `db.tmask.*.te.*`, `*.tl.*`, `*.th1.*`, `*.tm5.*.w0.*`, `*.3f`, `*.1a`) read
inputs written by `ncbi-parity-20261019/M/m4b/scripts3/gen_fix3.py` (seed 20261124): `edge.fna` (`e1`, 1,200 bases whose
bases 901-1020 are 40 sense codons of frame +2; `e2`, 1,202 bases, 400 sense codons and `AT`; `e3`, 1,200 bases, and
`e4`, 1,201, the minus strands of 400 sense codons; `e5` to `e8`, 10, 11, 16 and 17 random bases, around the shortest
subjects searched: 11 bases at word size 3, the default of both programs, 17 at TBLASTN's `-word_size 5`), `long.fna`
(`l1`, 7,000 bases with coding stretches at 3000-3599, frame +1, and 5000-5449, frame -1; `l2`, 3,000 bases: longer than
the average, so the traceback fetches parts of `l1`), the TBLASTX queries `qxe.fna`, the TBLASTN queries `qne.faa`, the
PSSM query `pn1.faa` (the PSSM is `../fixtures/pssm/m4b2_pl1.asn.gz`) and the hand masks `hedge3_<case>.asnt` (from/to inclusive below as `[start, end)`):

| `te` | Mask | Case |
|---|---|---|
| 11 `f1` | `e1` [600, 902) | covers the first base of frame +2's residue 300 (review F1) |
| 12 `n47` | `e2` [4, 7) | searched |
| 13 `n1200` | `e2` [1200, 1202) | NCBI crashes (the minus frames start with an empty part); searched as a frame without a mask, the masked bases staying masked (PD-LOSAT-NCBI-DEFECTS exception 27, the Owner's answers Q4 = B and Q6 = B of M4b, M4b-3c; near input `te13n`) |
| 14 `h01`, 15 `h12` | `e3` [0, 1), [1, 2) | NCBI crashes (the plus frames start with an empty part); searched as a frame without a mask, the masked bases staying masked (PD-LOSAT-NCBI-DEFECTS exception 27, the Owner's answers Q4 = B and Q6 = B of M4b, M4b-3c; near input `te14n`, `te15n`) |
| 16 `h23`, 17 `h13` | `e3` [2, 3), [1, 3) | frame -1's part before the mask searched alone (no near input keeps it; NCBI's own rows) |
| 18 `h46`, 19 `h46b` | `e3` [4, 6), `e4` [4, 6) | frames -2 and -3, frame -3: the part before the mask searched alone, as 16 |
| 31 `two` | `e4` [0, 3), [1200, 1201) | NCBI crashes (frame -1 starts with an empty part, then a part before the mask [0, 3)); searched as a frame without a mask, the masked bases staying masked (PD-LOSAT-NCBI-DEFECTS exception 27, the Owner's answers Q4 = B and Q6 = B of M4b, M4b-3c; near input `te31n`), frame -1's part searched alone as 16 |
| 32 `m10` to 35 `m17` | the last 3 bases of `e5` to `e8` | the plus frames' parts before the mask searched alone where the subject is searched (from 11 bases at word size 3, 17 at TBLASTN's `-word_size 5`: row `te.34.w5`), else not searched |
| 36 `s01`, 37 `s1200` | `e3` [0, 1), [600, 700); `e2` [300, 400), [1200, 1202) | searched (near misses of 14 and 13) |

`tl` has 11 (`long`): `l1` [3270, 3330) and [5201, 5251). The two databases were made with:

```
makeblastdb -in edge.fna -dbtype nucl -out te -title "edge masks" -blastdb_version 4 -mask_data hedge3_f1.asnt,\
hedge3_n47.asnt,hedge3_n1200.asnt,hedge3_h01.asnt,hedge3_h12.asnt,hedge3_h23.asnt,hedge3_h13.asnt,hedge3_h46.asnt,\
hedge3_h46b.asnt,hedge3_two.asnt,hedge3_m10.asnt,hedge3_m11.asnt,hedge3_m16.asnt,hedge3_m17.asnt,hedge3_s01.asnt,\
hedge3_s1200.asnt
makeblastdb -in long.fna -dbtype nucl -out tl -title "long masks" -blastdb_version 4 -mask_data hedge3_long.asnt
rm -f *.njs
```

makeblastdb numbers nine algorithms per masking program (dust 11-19), so the masks after the ninth name windowmasker
(31-37); the program only names the algorithm.

The rows of M4b-3c (2026-11-23; `db.tmask.*.hard.te.{13,14,15,31}.*`, `*.hard.tq4.*` and their near rows
`*.hard.near.te{13,14,15,31}n.*`, `*.hard.near.tq4n.*`) read inputs written by `../make_dbmaskt_q4.py` from `edge.fna`:
a mask that leaves a reading frame an empty first part makes NCBI crash in every run (PD-LOSAT-NCBI-DEFECTS exception
27, the Owner's answers Q4 = B of M4b of 2026-11-22 and Q6 = B of 2026-11-23), so each row's expectation is NCBI's
run of its near input, the same database with that short mask dropped and its bases stored as `N` (not the literal
mask removal, which keeps the real bases and gives other rows). The near input is equal to LOSAT's reading only as
long as NCBI's preliminary stage reads no masked residue: makeblastdb stores an `N` as a random base (different from
the real one at every masked position here), which the near input's seed scan and early extensions read, where
LOSAT and NCBI's own runs with a far mask read the real base. The near inputs' bases (`edge_n13.fna`: `e2`'s bases 1200 and 1201; `edge_n14.fna`: `e3`'s base 0; `edge_n15.fna`: `e3`'s
base 1; `edge_n31.fna`: `e4`'s base 1200, its mask [0, 3) kept), with no other mask on a searched sequence (`hedge4_n13`,
`hedge4_n14`, `hedge4_n15`: `e5`'s last 3 bases, a subject no search reads; `hedge4_n31`: `e4` [0, 3)). `tq4` holds
`e3`'s empty mask at 0 (`from 0, to 4294967295`) and [1100, 1200) (`hedge4_q4`); its near input `tq4n` the tail 4 bases
shorter, [1100, 1196), with `e3`'s bases 1196-1199 stored as `N` (`edge_nq4.fna`, `hedge4_nq4`: exception 21's near
input). The TBLASTN queries `qneq4.faa` are 60 residues at each frame end that holds a masked base, the TBLASTX query
`qxeq4.fna` (`x4t`) `e4`'s bases 1001-1200. Each database has the one algorithm 11:

```
python3 ../make_dbmaskt_q4.py . $NCBI_BIN   # writes the files above, then for each of te13n te14n te15n te31n tq4 tq4n:
makeblastdb -in edge_n13.fna -dbtype nucl -out te13n -title "near te 13" -blastdb_version 4 -mask_data hedge4_n13.asnt
makeblastdb -in edge_n14.fna -dbtype nucl -out te14n -title "near te 14" -blastdb_version 4 -mask_data hedge4_n14.asnt
makeblastdb -in edge_n15.fna -dbtype nucl -out te15n -title "near te 15" -blastdb_version 4 -mask_data hedge4_n15.asnt
makeblastdb -in edge_n31.fna -dbtype nucl -out te31n -title "near te 31" -blastdb_version 4 -mask_data hedge4_n31.asnt
makeblastdb -in edge.fna -dbtype nucl -out tq4 -title "empty mask at 0 and a tail" -blastdb_version 4 -mask_data hedge4_q4.asnt
makeblastdb -in edge_nq4.fna -dbtype nucl -out tq4n -title "near tq4" -blastdb_version 4 -mask_data hedge4_nq4.asnt
rm -f *.njs
```

## A database title outside the visible ASCII range (M4a-1b part 5, M4a-1b2)

`title_bad` (version 4, nucleotide): the records of `tests/fasta/outfmt11/arch_s.fna` and `arch_tbn_s.fna`, built by NCBI
BLAST+ 2.17.0 `makeblastdb -dbtype nucl -blastdb_version 4 -title "café db <DEL>"` (the title holds the bytes 0xC3 0xA9 and
0x7F). The rows `db.dbtitle.*` (BLASTN and TBLASTN, outfmt 8 and 9) write the title into each query's Seq-annot (`desc` user
object `Blast Database Title`) with `#` for each such byte, and NCBI posts `Error: (803.7) [blastn]
Seq-annot.desc.E.user.data.E.label.str` / `Bad char [0xC3] in string at line N` (outfmt 9: `at byte N`) / the title for each
byte of each query's document.

`seqannot/` (M4a-1b2): `tdel` (version 5, `-parse_seqids`: `subj_del.fa`, one sequence of 16,800 nucleotides, 84 copies of a
200-nucleotide unit, whose id `lcl|s<DEL>1` makeblastdb stores as `s#1`) and `pbad` (version 4, `pdel.faa`, the protein of
`subj_del.fa`'s first frame, 5,600 letters), both built by NCBI BLAST+ 2.17.0 makeblastdb with `-title "café db <DEL>"`
(`-dbtype nucl -parse_seqids -blastdb_version 5`, `-dbtype prot -blastdb_version 4`); `q2n_del.fa` holds the unit and its
reverse complement as `lcl|q<DEL>1` and `lcl|q<DEL>2`, and `qp.faa` the protein of the unit's first frame (stops as `X`).
The rows `db.dbtitle.blastn.qdel.8` and `db.dbtitle.blastx.qdel.9` write the title and, with `-parse_deflines`, the query ids
with `#` and NCBI's 803.7 message for each byte; fasta_input's `blastp.seqannot_delid_s.*` searches `qp.faa` against
`pdel.faa` (`-subject -parse_deflines`, its id `lcl|p<DEL>1`).

A `pdb` Seq-id whose chain is an INTEGER outside 0-255 without a chain-id (line RG item K4j; rows `db.pdbchain.*`):
`pdbchain/make_pdbchain.py` (run from `LOSAT/` with the directory of NCBI's makeblastdb) writes the version 4 database
`pdbchain/pdbchain` of five sequences without a title, each with a `gb|AB60000<n>.1|` Seq-id and a `pdb` one with the
chain 321, -191, 383, 65 or none, by makeblastdb from `gb|AB60000<n>.1|pdb|1AB<n>|<letter>` and a rewrite of the stored
deflines and of the chain letter in the keys of `.psd`. NCBI's `x_SetTitleFromPDB` reads the chain as
`(unsigned char) m_PDBChain`: `Chain A` for 321, -191 and 65. The gb Seq-id lets NCBI load the sequence for its report
(with the pdb one alone it prints `Unknown` and `no longer exists in database`).
- `tn_amb_q.faa` / `tn_amb_s.fna` (rows `db.tn.amb.*`, line O4c): a protein of 80 random letters with an X, a B, a Z
  and a W, and its back-translation with an `RAT` codon (D or N) at the X, the B and a D, an `SAR` codon (E or Q) at
  the Z and `TGA` (a stop in code 1) at the W, between `CC` and `GG`; then (review X2 of line O4c) 32 more letters,
  written by hand without repeats so that SEG keeps them and the alignments run through them, with an I, an L, a J and
  an X against an `MTT` codon (I or L) and an A and an X against `NNN`. TBLASTN's search reads the ambiguous codons as
  X, its reports as B, Z and J (`NNN` stays X in both): the identities and positives of the reports differ from the
  Seq-aligns'.
- `tn_ac_q.faa` (rows `db.tn.rej.ungapped.sw.ac.b1.*` and `db.tn.ungapped.ac.b1.nosw.7`, line O4c): a protein of W, C
  and H alone, from which NCBI computes no ungapped Karlin-Altschul parameters (a query batch without a valid
  context, which NCBI reports with a warning and does not search), then the five G of `nohit.faa`.

`x1bp/`: rows `db.x1bp.*` (line RG, item X1b-p). `d3m`, `rs` and `d3` are `makeblastdb -dbtype nucl` (NCBI BLAST+
2.17.0, run from `LOSAT/`) of `tests/fasta/outfmt0/x1bp_d3m_s.fna`, `x1bo_R_short_s.fna` and `x1bo_d3_s.fna`: subjects
on which `-ungapped -xdrop_ungap 1e10` (and megablast for `rs`) keeps an HSP outside its query context. Fix 2 (review 1
M1 and S1): `qpa` and `r3` the same of `x1bp_qpA_s.fna` and `x1bo_r3_s.fna`; `d3t` of `x1bo_d3_s.fna` with
`-taxid 9606 -blastdb_version 4`, for the prefetch of row `db.x1bp.d3t.db.6.rej`, whose
`BLASTDB=tests/seqdb/taxdb/bad_magic` makes NCBI's prefetch write a message.

## Where some inputs come from (provenance)

- `tm_psi_rhx.faa`, `tm_psi_rhy.faa` (commit eb27ebee4, 2026-11-04): the FASTA of the databases above; origin of the letters
  not recorded.
- `ids/plook_q*.faa`, `ids/plook_q*.fna`, `ids/plook_z_q0.faa`, `ids/nlook_q*.faa`, `ids/nlook_q*.fna`: written by
  `ids/make_look_fasta.py` (random proteins, seed 20261103; commit d191a9d3a, 2026-11-02).
- `ids/pord_q*.faa`, `ids/pord_q*.fna`, `ids/pord_x_q0.faa`, `ids/nord_q*.faa`, `ids/nord_q*.fna`: written by
  `ids/make_ord_fasta.py` (random proteins, seed 20261102; commit 99074fc9f, 2026-11-02).
- `ids/prdn_u0y.fna`: input of `ids/make_ids.py --pssmread-only` for the database `prdn_u0y`; origin of the letters not
  recorded (commit a18dd4682, 2026-11-11).
- `multidb/ld*.fa`, `multidb/l5b.fa`, `multidb/l6*.fa`, `multidb/nutl.fa`, `multidb/putl.fa`, `multidb/pfz*.fa`: written by
  `make_multidb.py` (its docstring tells, per mode, which records they hold and which come from `ids/ids_q.faa`, from
  `multidb/nnp.fa` or from seeded random letters); commits 2104f0950, dfa9763cd and db2e5a338 (2026-11-16), 6b1b2a43d,
  bfa554af4 and dc934f255 (2026-11-17), 1a529b829 (2026-11-18).
- `xpost/s60*.faa`: written by `make_xpost_inputs.py` (commit be806c7d4, 2026-11-16).
- `xpost/m2b_*_qm.*`: made by agent4's `m2/probe_m2b.py`, and `xpost/db/{nST,pA,pB,pR3}/*.fa` by the probes of K4h review 2,
  as `make_xpost_inputs.py`'s docstring says (commit 135255760, 2026-11-16); neither script is in this repository, and the
  origin of the letters is not recorded.

## `html/`: inputs of the `-html` rows (line RG, item T4a2)

Rows `db.html.*` (`../db_search_manifest/<program>.html.tsv`): `-html` in formats 0 to 4 and 18 with `-db` (formats 1 to 4: line RG T4b), and the registry
and environment cases of NCBI's HTML formatter. The inputs are those of the survey T4 of line RG
(`remaining-gaps-20261112/T/t4`, written by its `scripts/gen.py` and `in2/gen2.py`), copied byte for byte so that the
frozen bytes can be compared with the survey's runs. The rows also use `of18/ps_ids`, `of18/ns_ids`, `of18/ps_plain`,
`of18/ps_plain_8890` and `of18/pq.faa` (the T2 ids databases), `dbmask/edge1`, `dbmask/edge2` and `dbmask/q5.fna`, and
the taxonomy databases `../seqdb/taxdb/good` and `../seqdb/taxdb/of18bl`, and the queries and the `-subject` file of
`../fasta/html` (the inputs of the `html.*` rows of `../outfmt0_manifest.tsv`; its README.md): `pq_spec.faa`,
`nq_spec.fna`, `pq_ids.faa` (`-parse_deflines`) and `ps_spec.faa` (the `-subject` of the registry rows, `Subject_N`).

- `pq_nohit.faa` (`pq2`, a query without hits) and `pq_nohitfirst.faa` (that query, then `q1` of `pq_spec.faa`).
- `pmix_ids` (protein) and `nmix_ids` (nucleotide), version 5 databases of the id kinds that `of18/ps_ids` and `ns_ids`
  lack: `../fasta/html/pmix_ids.faa` holds `pir||A12345`, `prf||1234567A`, `pat|US|5000000|1`, `tpg|DAA00001.1|`,
  `gi|1000|gb|AAA00001.1|`, `gnl|TI|12345`, `sp|Q99999.1|NAME_HUMAN` with a second defline `emb|CAA99999.1|` (joined by
  `\x01`), a sequence with nine deflines `ref|XP_000901.1|` to `ref|XP_000909.1|` (more than NCBI's eight shown: `6 more
  sequence titles`), `lcl|loc_x` and `gi|999`; `nmix_ids.fna` holds `gi|2000`, `gi|2001|emb|CAA00001.1|`,
  `tpe|BN000001.1|`, `gnl|TI|54321`, `pat|EP|1000000|2` and `gb|AY000001.1|` with a second defline `ref|NM_000099.1|`.
  Built by NCBI BLAST+ 2.17.0 in the survey's directory, so their titles are the `-in` text
  (`BLAST_USAGE_REPORT=false`):

```
makeblastdb -in in/pmix_ids.faa -dbtype prot -parse_seqids -out db/pmix_ids
makeblastdb -in in/nmix_ids.fna -dbtype nucl -parse_seqids -out db/nmix_ids
```

- `rc/<name>/`: the run directories of the `cwd=` rows and the `$NCBI` directories of the `NCBI=` rows, each with the
  file that NCBI reads there: `.ncbirc` of `p_http` (`[BLASTFMTUTIL]` `PROTOCOL=http:`), `p_https` (`https:`), `p_empty`
  (`PROTOCOL=` empty), `lower` (`[blastfmtutil]` `protocol=http:`), `taxbr` (`[TAX_BROWSER]` `BLASTFMTUTIL=...`),
  `localid` (`[LOCAL_ID]` `TOOL_URL=` with the tags `seq_id`, `db_name`, `taxid`, `protocol` and `title`), `feat`
  (`[FEATURE_INFO]` `FEATURE_FILE`, `FEATURE_FILE_INDEX`), `plain_only` (`[BLAST]` `DATA_LOADERS=blastdb`), `ncbidir` (a
  `$NCBI` directory: `[LOCAL_ID]` `TOOL_URL`), and the malformed files, whose first line `this line has no equals sign`
  NCBI's registry refuses: `bad/.ncbirc` (the run directory's), `ncbibad/.ncbirc` (a `$NCBI` directory's),
  `fmtbad/bad.rc` (the file that `FMTCFG=bad.rc` names in the run directory) and `ncbifmt/bad.rc` (the same in a `$NCBI`
  directory); `fmtcfg/other.rc` is a valid file named by `FMTCFG=other.rc`. Two files that NCBI's registry reads with a
  post on stderr at each read (line RG T4a2 fix 1, B1: the references', the description table's and the alignments'
  reads of `-html`): `rep/.ncbirc` (`[BLASTFMTUTIL]` with `PROTOCOL=https:` twice, `Found multiple ... settings`) and
  `bom/.ncbirc` (the one byte 0xEF, the start of a UTF-8 byte order mark, `Error reading the registry after line 1`,
  survey T2's `bom` file).

## Titles with bytes outside the visible ASCII range (line RG-KB item 9(b))

`tchar/` (rows `db.tchar.*`): `tchar/make_tchar.py` (run from `LOSAT/` with the directory of NCBI's makeblastdb,
and the stems to write, or none for all) writes the version 4 databases `tc_p`, `tc_b_p`, `tc_big_p`, `tc_big2_p`,
`tc_x_p`, `tc_dmg_p`, `tc_ps_p`, `tc_u8_p` (proteins) and `tc_n`, `tc_b_n`, `tc_big_n`, `tc_big2_n`, `tc_x_n`,
`tc_dmg_n`, `tc_ps_n`, `tc_u8_n` (their back-translations), the aliases `tc_alias_p.pal` and `tc_alias_n.nal`
(`tc_*` then `tc_b_*`), `tc_bigal_p.pal` and `tc_bigal_n.nal` (`tc_big_*` then `tc_big2_*`) and the queries
`tc_q1`, `tc_q2`, `tc_qslot` and `tc_qu8` (`.faa`, `.fna`). makeblastdb stores each byte outside ' '..'~' of a FASTA
title as `#`; the script rewrites those stored titles to raw bytes (0x00, 0x01, 0x1B, 0x7F, 0x80, 0xC3 0xA9, 0xE9,
0xFF; titles of 133, 300 and 1,100 bytes). `tc_big_*` holds 300 sequences, of which the OIDs 0, 1, 2, 3, 256 and 257
have such bytes and hit the queries with OID 258 (clean): 0 and 256, 1 and 257, 2 and 258 share a slot of NCBI's
cache of 256 header records per volume. `tc_qslot` reads 0 and 256 in both orders, then 2, 258 and 2 again, then 3
twice (`-evalue 1e-10` keeps the random records out). `tc_big2_*` holds 300 sequences, OID 215 with 0xE9, which
`tc_bigal_*` reads between the reads of `tc_big_*`'s OID 3 (215 is 3's slot by the alias's OIDs, 515, not by the
volume's). `tc_dmg_*` holds OID 0 with 0xE9 and OID 1 whose title tag is rewritten to 0x1B (NCBI's error, exit 255,
after the post of OID 0); `tc_ps_*` has parsed Seq-ids (`-parse_seqids`); `tc_u8_*` holds a UTF8String title,
which `tc_qu8`'s first query reads (NCBI's warning 803.10), and titles with 0xE9 for the second query.
NCBI BLAST+ 2.17.0 reads a reported subject's definition line where it makes the Seq-aligns of its query batch, and
posts `Error: (803.7) [<program>] Blast-def-line-set.E.title` / `Bad char [0x..] in string at byte N` / the title, for
each such byte in order, at each read that the cache does not answer (N is the position after the whole title in its
record; the bytes before the posted one are already `#`). Outfmt 8, 9, 11 and 12 read the definition line there only
(`PRE_FETCH_SEQS_LIMIT=2147483647` makes 8, 9 and 12 prefetch); the other formats read it again, which LOSAT refuses
(the `losat-rejection` rows), except in a batch that fails at a header read, which reaches no formatter (`dmg.6`).
