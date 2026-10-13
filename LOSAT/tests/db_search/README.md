# Inputs of the `-db` fixtures

Used by `LOSAT/tests/db_search_fixtures.py` and `LOSAT/tests/db_search_manifest.tsv` (rows that run
`LOSAT blastp -db ...` against frozen NCBI BLAST+ 2.17.0 `blastp -db ... -num_threads 1` bytes). NCBI is the
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
- `p3.faa`: records p01, p03 and p08 of `LOSAT/tests/seqdb/prot.faa`, in that order.
- `d1_q.faa`: D1, a hit that `-db` keeps and `-subject` drops. Two proteins of `LOSAT/tests/fasta/WSSV.faa`
  (YP_009220551.1 and YP_009220614.1: for each, NCBI `blastp -db paja -outfmt 6` prints 2 rows and
  `blastp -subject tests/fasta/PajaWSV.faa -outfmt 6` prints 1), then `nohit`, a 72-residue low-complexity
  glycine/serine sequence for which NCBI `-db paja` reports no hit. For the whole file: 4 rows with `-db`, 2
  rows with `-subject` (outfmt 6).

Databases (all `-blastdb_version 4`, no `-parse_seqids`; only `.pin .psq .phr` are kept, the `.pjs` file
makeblastdb also writes was deleted). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/residues.faa -dbtype prot -out tests/db_search/residues \
        -title "residues O U * B Z J X" -blastdb_version 4
    makeblastdb -in tests/db_search/residues_gap.faa -dbtype prot -out tests/db_search/residues_gap \
        -title "residues O U * - B Z J X" -blastdb_version 4
    makeblastdb -in tests/fasta/PajaWSV.faa -dbtype prot -out tests/db_search/paja \
        -title "PajaWSV proteins" -blastdb_version 4
    makeblastdb -in tests/db_search/o_stored.faa -dbtype prot -out tests/db_search/o_stored \
        -title "stored O against the query's X" -blastdb_version 4

`residues_gap` is the database LOSAT rejects (a stored `-`). The other databases of the rows are the ones of
`LOSAT/tests/seqdb/db` (see `LOSAT/tests/seqdb/README.md`). A rebuild changes the dates stored in the index
files, so the frozen outputs (`Posted date:` lines) must then be frozen again:
`python3 LOSAT/tests/db_search_fixtures.py freeze --ncbi-bin DIR`.
