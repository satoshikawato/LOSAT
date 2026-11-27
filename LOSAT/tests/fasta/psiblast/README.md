# Inputs of the psiblast fixtures

Used by the `db.ps.*` rows of `LOSAT/tests/db_search_manifest/psiblast.tsv` (against `tests/db_search/paja`) and
the `ps.subj.*` rows of `LOSAT/tests/outfmt0_manifest.tsv` (against `-subject tests/fasta/PajaWSV.faa`).
The oracle is NCBI BLAST+ 2.17.0 `psiblast ... -num_threads 1`.

- `tk.faa`: YP_009220601.1 of `tests/fasta/WSSV.faa` (TK-TMK, 398 letters).
- `two.faa`: YP_009220601.1 and YP_009220615.1 of `tests/fasta/WSSV.faa`.
- `pol.faa`: YP_009220649.1 of `tests/fasta/WSSV.faa` (DNA polymerase, 2,351 letters).
- `rnd.faa`: a random 30-letter protein without hits.
- `tk_rnd.faa`: `tk.faa`, then `rnd.faa` (PS2: the second query writes no PSSM).
- `tk_xbz.faa`: `tk.faa` with X at positions 50-55, B Z U J at 100-103 and `*` at 201 (1-based), its
  definition line ending in ` with X B Z U J *`, 60 letters a line (PS2: the PSSM's X columns and
  printed letters).
- `tk_o.faa`: `tk.faa` with O at positions 30 and 121 (1-based), its definition line ending in ` with O`, 60
  letters a line (PS2: an O of the query is X in the search's buffer and in the PSSM, with NCBI's warning).
- `tk_x_pol.faa`: `tk.faa`, a record of 20 X (no valid Karlin-Altschul parameters), `pol.faa`: psiblast
  stops at the second query with `BLAST engine error: Could not calculate ungapped Karlin-Altschul
  parameters ...`.
- `tk_lc_pol.faa`: `tk.faa`, a record of 60 A, `pol.faa`: with `-seg yes` the second query is masked
  entirely and stops the run the same way.
- `tk_e_pol.faa`: `tk.faa`, a record without residues, `pol.faa`: the empty record and the next query form
  one batch of psiblast (`CBlastInput(&*fasta, 1)`), which `CPsiBlast` refuses (`PSI-BLAST only accepts one
  protein sequence as query`).
- `e_only.faa`: a record without residues alone, on which NCBI 2.17.0 psiblast crashes (LOSAT rejects it).
- `x1800.faa`: a record of 1,800 X (no valid Karlin-Altschul parameters). With `CHUNK_SIZE=1001` it is
  not split (1,800 < 2 x (1,001 - 100)) and psiblast stops with the BLAST engine error; with
  `CHUNK_SIZE=1000` NCBI 2.17.0 splits it and crashes (LOSAT rejects it).
- `notitle.faa`, `quote.faa`, `long.faa`, `lc.faa`, `tab.faa`, `utf8.faa` (PS3, `-out_pssm`): `tk.faa`'s letters
  under other definition lines, for the title of the PSSM's query Bioseq: none (`>` alone: `title ""`), quotes
  (doubled in the ASN.1 string), a 685-byte title (wrapped at word ends), `tk tab<TAB>here` (the reader keeps
  `tk tab`), `tk café title` in UTF-8 (LOSAT refuses it with `-out_pssm`; NCBI writes `#` for each byte with an
  error message); `lc.faa` has letters 61-120 in lower case (the PSSM's query is the search's buffer, upper case).
- `ps4b_q2.faa` (PS4b, design question Q2): a record without residues (`>e1 empty record`), then BDT63550.1,
  BDT63556.1 and BDT63563.1 of `tests/fasta/PajaWSV.faa`: the `Subject sequence contains no data` warning comes
  before the PSSM query's warning line.
- `ps4b_stale_s1.faa` (PS4b review N3): `s_bz` (24 `BZ`, no true amino acid) and `s_good` (`tk.faa`'s 398 letters; the query
  of `tests/fixtures/pssm/ps4_stale_q` is three TK stretches with 40-letter `DE` runs between them), for NCBI's adjusted PSSM kept across matches: `s_good`
  is redone first, then `s_bz`'s first HSP is not redone (its adjustment fails) and its later HSPs are redone on the
  PSSM `s_good` left.
- `ps5_*.faa` (PS5a, psiblast's rounds; the queries of the step-1 cases, `PS/scripts/ncbi_ps5_check.py`):
  `ps5_pk.faa` YP_009220615.1 (putative protein kinase; with `-inclusion_ethresh 5` five rounds without
  convergence, a sixth that converges), `ps5_rr1.faa` YP_009220529.1, `ps5_exhaust.faa` YP_009220499.1 and
  `ps5_conv6.faa` YP_009220492.1 (with `-inclusion_ethresh 10`: six rounds without convergence, and convergence on
  the sixth), `ps5_nohit2.faa` YP_009220505.1 and `ps5_nohit2b.faa` YP_009220545.1 (round 1's hits are all above
  0.002: round 2 searches the PSSM of the query alone and has no hits), all of `tests/fasta/WSSV.faa`;
  `ps5_three.faa`: `tk.faa`, `ps5_nohit2.faa`, `qnone` (36 G, no hits) and `rnd.faa` (the rounds start again for
  each query; outfmt 7 counts the rounds); `ps5_subj_one.faa` (BDT63550.1) and `ps5_subj_two.faa` (BDT63550.1 and
  BDT63556.1) of `tests/fasta/PajaWSV.faa`, subject sets of one and two records (the two-part table with `-subject`).
- `ps15b_ploop_rnd.faa` (PS15b-1, psiblast's per-round archive and XML in rows `db.ps15b.*` of
  `tests/db_search_manifest/psiblast.ps15b.tsv`): letters 20-35 of `tk.faa` (its P-loop), then `rnd.faa`. With
  `-num_iterations 3 -inclusion_ethresh 0.005 -evalue 1` the first query takes two rounds (the second converges) and
  `rnd` one, so that the archive of three rounds and two queries stays small: a PSSM round's request carries a
  PSSM of 28 columns a letter.
- `ps15b_lcmix.faa` and `ps15b_pd_pol.faa` (review fixes of PS15b-1, rows `db.ps15b.fix.*`): `tk.faa` with letters
  31-70 and 151-175 in lower case (`-lcase_masking`: round 2 aligns across them, where psiblast's XML shows the
  query's letters, not X); `pol.faa` under the made-up defline `>sp|P99998.1|POL_WSSV DNA polymerase test`, searched
  against itself with `-parse_deflines -seg yes` (psiblast's outfmt 0 draws no query masks, so the record that is
  both a query and a subject is written, where blastp refuses it).
- `ps15b_lcmix_pd.faa` (second review fixes of PS15b-1, rows `db.ps15b.fix.self.lcm.*`): `ps15b_lcmix.faa` under the made-up
  defline `>sp|Q99997.1|TKLC_TEST lowercase stretches`, searched against itself with `-parse_deflines -lcase_masking`
  (psiblast's outfmt 0 draws no query masks, so the record that is both a query and a subject is written).
- `ps5c_big.faa` (PS5c): YP_009220597.1 of `tests/fasta/WSSV.faa` (capsid protein, 6,077 letters), searched
  against itself with `-searchsp 9223372036854775807 -evalue 1e-300`: a hit cutoff whose `K × searchsp / E`
  is past the double (`ncbi_cutoffs::cutoff_score_from_evalue`).
- `ps5c_qx.faa` (PS5c): one query of 30 X, which gives no Karlin-Altschul parameters (NCBI does not search its
  batch); `ps5c_qx_tk.faa`: `ps5c_qx.faa`, then `tk.faa`. With `BATCH_SIZE=1` and a search that fails in its
  setup (here the empty database `tests/seqdb/db/empty`), NCBI reports the first batch, then the error of the
  second (blastp rows `db.empty.unsearched.*`).
- `ps7_mixed.fna` (PS7, tblastn `-in_pssm` with `-subject`, rows `ps7.subj.*` of `tests/outfmt0_manifest.tsv`): the
  first record of `tests/db_search/tslice.fna`, an empty record (`Subject sequence contains no data`) and the second
  record in lower case.
- `q1_*.faa` (line RG, series Q, Q1: psiblast's filters, `-show_gis`, `-parse_deflines` and `-lcase_masking` in
  its rounds; rows `db.ps.q1.*` and `ps.q1.*`): `q1_lc.faa` is `tests/db_search/taxfilter/q.faa` (P69905.2, 142
  letters) with letters 30-60 in lower case; `q1_tfp_lc.faa` is `tests/db_search/taxfilter/tfp.faa` with letters
  20-70 of records 2, 5 and 9 in lower case; `q1_subj_ids.faa` holds the first 8 records of `tfp.faa` under the ids
  `gi|4504347|ref|NP_000549.1|`, `sp|P01942.2|HBA_MOUSE`, `lcl|my_seq1`, `gnl|db|tag7`, `P01966.2`, `pdb|1ABC|A`,
  `gb|AAA12345.1|`, `emb|CAA01234.1|` and the first title of each; `q1_gi.faa`, `q1_pdb.faa`, `q1_lclnum.faa`,
  `q1_gnlid.faa`, `q1_word.faa` and `q1_ref.faa` are the letters of `q.faa` under `gi|4504347|ref|NP_000549.1|`,
  `pdb|1A3N|A`, `lcl|123`, `gnl|mydb|42`, `myquery` and `ref|NP_000549.1|` (with `-parse_deflines` the PSSM's query
  takes that Seq-id); `q1_subj_q1.faa` is `q1_subj_ids.faa` after a first record `lcl|Query_1` (the letters of
  `q.faa` with W at position 71), the local id of an `-in_pssm` query (`tests/db_search/q1_sq1` is its database).
- `q4a_skip.faa` (line RG, series Q, Q4a): a record of 10 letters, which `-query_loc 20-300` skips, before the letters
  of `tk.faa`; `q4a_xprefix.faa`: 400 X, then the letters of `tk.faa` (with `-query_loc 401-798` the later rounds read
  the PSSM's X columns, whose frequency ratios are 0).
- `q4b_lc1_200.faa` (line RG, series Q, Q4b): the letters of `tk.faa` with letters 1-200 in lower case; with
  `-lcase_masking` round 1 masks them, and NCBI's later rounds index them again (round 1's filtering takes the
  lower-case mask out of the query's sequence block), so the later rounds find alignments in them.
  `q4b_inv.faa`: 20 lower-case letters, then `WCYHWCYHWC`; `q4b_inv_s.faa`: two subjects with `WCYHWCYHWC`. With
  `-lcase_masking -seg yes` round 1 masks the lower case only; round 2's SEG on those X letters masks the island too,
  so round 2 has no Karlin-Altschul parameters (NCBI 2.17.0 psiblast crashes there, exit 139).
- `q4c_tk_pol.faa` (line RG, series Q, Q4c): `tk.faa` then `pol.faa`; with `CHUNK_SIZE=1000` NCBI splits only the second
  query (rows `db.ps.it.split.two.*`).
- `ps8_*.faa` (line PS, PS8a fix round 1; the `db.phi.*` rows of `psiblast.phi.tsv`, with the patterns
  `tests/fixtures/phi/c2c.pat` `C-x(2)-C`, `x3.pat` `x(3)` and `brknum.pat` `[C1]-x-C`, a digit in a bracket):
  `ps8_q3.faa`: YP_009220477.1, YP_009220540.1 and YP_009220504.1 of `tests/fasta/WSSV.faa` (a search strategy of
  three PHI queries); `ps8_e1_q3.faa`: an empty record `e1`, then `ps8_q3.faa` (psiblast's batch of one query takes
  the empty record with the next protein); `ps8_emptylast.faa`: YP_009220477.1, then an empty last record
  (`Sequence contains no data`, exit 3); `ps8_x10_p30.faa` and `ps8_p30_x10.faa`: 10 `X` and the first 30 letters of
  `ps8_c2c.faa`, in both orders (an X-only query is searched); `ps8_s_dneq.faa`: a 63-letter subject rich in D, N,
  E and Q; `ps8_c2c.faa`: YP_009220489.1 of `tests/fasta/WSSV.faa`; `ps8_s_lc.faa`: its letters with the first 60 in
  lower case (`-lcase_masking` masks a subject's lower case); `ps8_sid_q.faa` and `ps8_sid_s.faa`: a query and a
  subject set whose first record has the query's Seq-id `gb|ABC12345.1|` (under `-parse_deflines` NCBI shows them as
  one sequence; LOSAT's rule refuses it). For the `db.phi.long.*` rows (PS8a dev merge): `ps8_q24k.faa`, one
  24,334-letter record made of BDT63125.1, BDT61687.1, BDT63580.1 and BDV50172.1 (the wsv360-like proteins of the
  `tests/fasta/` virus sets, longest first), long enough for blastp's query split, which NCBI never applies to
  PHI-BLAST; `ps8_x21000p.faa`: 21,000 `X`, then the letters of `ps8_c2c.faa`, in one record; `ps8_s_bdt63125.faa`:
  BDT63125.1 alone. For the `db.phi.cull.*` and `db.phi.bh.*` rows (PR #1654 R2): `ps8_cull_q.faa` (two-domain
  synthetic query) and `ps8_cull_s.faa` (three synthetic subjects, with pattern `tests/fixtures/phi/cwh.pat`
  `C-W-x(2)-H-C`), exercising preliminary culling and best-hit HSP drops under `-max_target_seqs 1`. For the
  `db.phi.varlen.*` rows (PR #1654 R3): `ps8_varlen_q1.faa` / `ps8_varlen_s1.faa` and `ps8_varlen_q2.faa` /
  `ps8_varlen_s2.faa` with `tests/fixtures/phi/cx12c.pat` (`C-x(1,2)-C`), exercising variable-length pattern
  hits starting at the same offset (`CACC`) under `-culling_limit` and `-best_hit_*`.
- `ps15b2_othr_*.faa` (line PS, PS15b-2 fix 1): subjects of the molecule-`other` PSSM
  `tests/fixtures/pssm/ps15b2_tk_other_stdaa.asn.gz` (rows `db.ps15b2.p.othr.*`), from `tk.faa`: `del100` without
  letter 101 (a query letter against a gap of the subject on a query row of its own: NCBI 2.17.0's merge writes
  through a null row, SIGSEGV), `cover` that subject and letters 2-300 (another subject's block opens that row: NCBI
  answers), `t2a_s1` letters 1-121 with a `W` after letter 61 (a later pass finds no row of its HSP), `t2a` that
  subject and letters 101-192 with a `W` after letter 132 (with `-evalue 1e-10`, the two HSPs order their shared
  segments in two ways), `tk195` letters 195-313 (at `-num_iterations 2`, round 2's two HSPs leave a segment's
  start behind its row's end when NCBI orders the segments: NCBI reads past the row's end, `CObject is corrupted`).
