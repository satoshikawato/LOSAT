# Inputs of the `-db` fixtures

Used by `LOSAT/tests/db_search_fixtures.py` and `LOSAT/tests/db_search_manifest.tsv` (rows that run
`LOSAT blastp -db ...`, `LOSAT blastn -db ...` or `LOSAT blastx -db ...` against frozen NCBI BLAST+ 2.17.0
`blastp -db ...`, `blastn -db ...` or `blastx -db ... -num_threads 1` bytes; the manifest's `program` column says which). NCBI is the
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
- `nohit.faa`: `GGGGG`, a query without hits in `nophr` (the database of `LOSAT/tests/seqdb/prot.faa`).
- `empty.faa`: an empty query file.
- `split_q.faa`: one query of 20,600 letters (a fixed pseudo-random filler with records p10, p05, p06 and p01 of
  `LOSAT/tests/seqdb/prot.faa` at 1001, 9851, 15001 and 20401), which NCBI searches in two chunks (CHUNK_SIZE
  10000, overlap 100: [0, 10000) and [9900, 20600)); rows `db.split.*` against `LOSAT/tests/seqdb/db/prot_v4`.
- `vol2.faa`: `A` (record p10 of `LOSAT/tests/seqdb/prot.faa`) and `B`, the same with every fifth letter
  replaced by `W` (`A` where it is `W`); `vol2_q.faa`: `A` as the query.

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
    makeblastdb -in tests/db_search/vol2.faa -dbtype prot -out tests/db_search/vol2_nophr \
        -title "two volumes, the second without its header file" -blastdb_version 4 -max_file_sz 128B
    rm tests/db_search/vol2_nophr.01.phr

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
- `bnophr.pin` and `bnophr.psq` are copies of `bslice.pin` and `.psq`: the database without its header file, used
  with the query without hits (with a hit NCBI fails with exit 2).

The database (version 4, no `-parse_seqids`; only `.pin .psq .phr` are kept, the `.pjs` file was deleted).
Run from `LOSAT/`:

    python3 tests/db_search/make_blastx_inputs.py
    makeblastdb -in tests/db_search/bslice.faa -dbtype prot -out tests/db_search/bslice \
        -title "protein slice of a shrimp virus proteome" -blastdb_version 4
    rm tests/db_search/bslice.pjs
    cp tests/db_search/bslice.pin tests/db_search/bnophr.pin; cp tests/db_search/bslice.psq tests/db_search/bnophr.psq

`residues_gap` is the database LOSAT rejects (a stored `-`). The other databases of the rows are the ones of
`LOSAT/tests/seqdb/db` (see `LOSAT/tests/seqdb/README.md`). A rebuild changes the dates stored in the index
files, so the frozen outputs (`Posted date:` lines) must then be frozen again:
`python3 LOSAT/tests/db_search_fixtures.py freeze --ncbi-bin DIR`.
