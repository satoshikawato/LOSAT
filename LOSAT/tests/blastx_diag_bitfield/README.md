# blastx_diag_bitfield

NCBI BLAST+ 2.17.0 (`blastx: 2.17.0+`, bioconda build of 2026-08-22) outputs for
`tests/blastx_diag_bitfield.rs`, made 2026-10-17 on the LOSAT workstation (one thread):

```
blastx -query query.fna -subject subjects.faa -window_size W -outfmt 6 -num_threads 1   # ncbi.subject.wW.out6
makeblastdb -in subjects.faa -dbtype prot -out s2db
blastx -query query.fna -db s2db -window_size W -outfmt 6 -num_threads 1               # ncbi.db.wW.out6
```

for W = 40 (the default) and W = 536870000.

- `query.fna`: Ensembl BDGP6.54 cDNA FBtr0078009 (Dbp21E2), title cut to the identifier.
- `subjects.faa`: Ensembl R64-1-1 proteins YLR106C (MDN1, 4,910 residues) and then YBR142W
  (MAK5), titles cut to the identifiers.

At W = 536870000 the diagonal offset of the second subject is 2 x 536870000 + 4,910 =
1,073,744,910, past 2^30 = 1,073,741,824, so every `s + diag_offset` of YBR142W wraps in the
31-bit `last_hit` of NCBI's `DiagStruct` (blast_extend.h:57-60) and its two-hit test never lets
a seed through: NCBI reports no YBR142W row (3 rows with `-subject`, 5 with `-db`), while the
6 / 7 rows of W = 40 include three YBR142W rows. YBR142W alone (offset W) gives the same 6 rows
at both windows.
