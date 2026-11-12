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
- `ps5c_big.faa` (PS5c): YP_009220597.1 of `tests/fasta/WSSV.faa` (capsid protein, 6,077 letters), searched
  against itself with `-searchsp 9223372036854775807 -evalue 1e-300`: a hit cutoff whose `K × searchsp / E`
  is past the double (`ncbi_cutoffs::cutoff_score_from_evalue`).
- `ps5c_qx.faa` (PS5c): one query of 30 X, which gives no Karlin-Altschul parameters (NCBI does not search its
  batch); `ps5c_qx_tk.faa`: `ps5c_qx.faa`, then `tk.faa`. With `BATCH_SIZE=1` and a search that fails in its
  setup (here the empty database `tests/seqdb/db/empty`), NCBI reports the first batch, then the error of the
  second (blastp rows `db.empty.unsearched.*`).
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
