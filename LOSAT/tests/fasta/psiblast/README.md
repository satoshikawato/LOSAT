# Inputs of the psiblast fixtures

Used by the `db.ps.*` rows of `LOSAT/tests/db_search_manifest.tsv` (against `tests/db_search/paja`) and
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
