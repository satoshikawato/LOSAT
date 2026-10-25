# Inputs of the psiblast fixtures

Used by the `db.ps.*` rows of `LOSAT/tests/db_search_manifest.tsv` (against `tests/db_search/paja`) and
the `ps.subj.*` rows of `LOSAT/tests/outfmt0_manifest.tsv` (against `-subject tests/fasta/PajaWSV.faa`).
The oracle is NCBI BLAST+ 2.17.0 `psiblast ... -num_threads 1`.

- `tk.faa`: YP_009220601.1 of `tests/fasta/WSSV.faa` (TK-TMK, 398 letters).
- `two.faa`: YP_009220601.1 and YP_009220615.1 of `tests/fasta/WSSV.faa`.
- `pol.faa`: YP_009220649.1 of `tests/fasta/WSSV.faa` (DNA polymerase, 2,351 letters).
- `rnd.faa`: a random 30-letter protein without hits.
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
