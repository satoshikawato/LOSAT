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
