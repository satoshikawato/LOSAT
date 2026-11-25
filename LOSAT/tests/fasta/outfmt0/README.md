# Inputs of the `outfmt0_manifest.tsv` rows: where the files below come from

These files are inputs of rows of `LOSAT/tests/outfmt0_manifest.tsv` (some also of `db_search_manifest` rows); the frozen
outputs in `../../fixtures/outfmt0/` are NCBI BLAST+ 2.17.0's. "Origin not recorded" means that no generator or source of the
letters is kept in the repository; the commit is the one that added the files.

- `culling_dom_s_rev.faa`, `culling_dom_s_rev.fna`: the records of `culling_dom_s.faa` and `culling_dom_s.fna` in reverse
  order (checked: the same records, ids stay with their records); rows `culling.blastp.dom_rev_cl2` and
  `culling.tblastn.dom_rev_cl1` (commit 6986c11d1, 2026-10-22).
- `blastp_mt_redo_query_rev.faa`: the two records of `blastp_mt_redo_query.faa` in reverse order, row
  `mt.blastp.redo_last.6` (commit 76a2137ee, 2026-10-12). The defline names the proteins: MTB_ONCMY (UniProtKB P68501) and
  PUR4_HUMAN (UniProtKB O15067); UniProt data are CC BY 4.0.
- `sam_long_q.fna`, `sam_long_s.fna`: a pair of 30,000 bases that differ in one base, rows `sa.sam.long.17p` and
  `sa.sam.long.17b` (commit b1543abf2, 2026-10-22); origin of the letters not recorded.
- `x2c_d7k_*`, `x2c_e7k_*`, `x2c_g48_*`: written by `gen_x2c.py` (kinds d, e, f, g; the script is not in the repository),
  inputs of `tests/tblastx_query_batch.rs` (commit 232ae2d70, 2026-11-16).
- `x1bp_cov_q.faa`, `x1bp_cov_s.faa`, `x1bp_cov_s.fna`: the record `q_plain` of
  `../../db_search/residues_q.faa` (142 residues), its residues 11-130, and its residues back-translated between two
  random flanks of 300 bases (checked: the 426 bases from base 301 translate to `q_plain`); rows `db.x1bp.cov.bp.6` and
  `db.x1bp.cov.tbn.6` (line RG X1b-p, review 2's probe D, 2026-11-23); the flanks' generator is not in the repository.
- `x2b_r7_*`, `x2b_r11_*`: pairs of rows `x2b.tblastx.*` (`-window_size 0`, HSPs of length 0 or less); origin not recorded
  (commit f28d5c5a3, 2026-11-14).
- `xtl2_wdw.faa` (W D W), `xtl2_wdw_q.fna` (`TGGGATTGG`, W D W in frame +1) and `xtl2_s11pad.fna` (an 11-base subject
  `TGGTGGTGGTA`, W W W in frame +1, then 100 copies of `TAA`): written for rows `xtl2.*` (line RG item X-short-tbx-link-2,
  2026-11-25): with PAM30 and `-threshold 11` the word W D W against W W W is the only one that reaches 11 and gives an HSP of
  length -1; the 300 bases after the subject find nothing and make NCBI's linking cutoff finite.
- `x2g_c0037_*` (commit 03b352b90) and `x2g_p0002_q.fna`, `x2g_p0086_*` (commit a8e16bfda), rows `x2g.blastx.*`: pairs of
  the line RG X2g survey (the commit messages: "rand2 cases", "random pairs of probe P1"); the generator is not in the
  repository (2026-11-14).
- `q3c_*`: inputs of rows `bx.q3c.*` (BLASTX `-ungapped`, linked sets of HSPs); origin not recorded (commit fdeabff90,
  2026-11-14).
- `sm3_*`: inputs of rows `sm3.*` (BLASTN `-soft_masking`); origin not recorded (commit 6314d8a9b, 2026-10-31).
- `e2d_*`: inputs of rows `e2d.*` (`-query_loc`, `-subject_loc`); origin not recorded (commit d67175bd3, 2026-10-05).
- `e2e_fast_lone*`: inputs of rows `e2e.blastp.fast_lone_hsp*` (`-task blastp-fast`); origin not recorded (commit 3b84ce6c4,
  2026-10-05).
- `x16_title_*`, `x16n_emptymid_s.fna`, `x16t_*`, `x16x_emptymid_s.fna`: inputs of rows `sa.x16.*` (outfmt 16); origin not
  recorded (commits 40aec1fea, 0748a144c, 6892d1e0e, f06344780, 2026-10-21 and 2026-10-22).
- `bx_title_q.fna` (commit fea5852cf, 2026-11-04), `bx_gapline_*` (commits 7452c1274 and 290933e7b, 2026-11-07),
  `seg_bx_subj_q.fna` (commit 364e28742, 2026-11-02), `x2_tbx_lc_q.fna` (commit c3c55d244, 2026-11-13) and
  `sml_bp_mixlc.faa` (commit 23ce32af3, 2026-10-31): inputs of the rows `title.blastx.6`, `gapline.blastx.*`,
  `seg.blastx.subj.6`, `x2.tblastx.lc.6` and `sml.blastp.mixlc.0`; origin not recorded.
- `twidth_q.faa`, `twidth_s.fna`: written by `docs/evidence/losat_web_e2a/make_inputs.py` (commit 63ddf7380, 2026-09-29),
  row `width.tblastn`.
