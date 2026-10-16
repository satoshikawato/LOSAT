# Inputs of the `-db` fixtures

Used by `LOSAT/tests/db_search_fixtures.py` and `LOSAT/tests/db_search_manifest.tsv` (rows that run
`LOSAT blastp -db ...`, `LOSAT blastn -db ...`, `LOSAT tblastn -db ...` or `LOSAT tblastx -db ...` against frozen NCBI BLAST+ 2.17.0
`blastp -db ...`, `blastn -db ...`, `tblastn -db ...` or `tblastx -db ... -num_threads 1` bytes; the manifest's `program` column says which). NCBI is the
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
    rm tests/db_search/tslice.njs tests/db_search/tgc4.njs

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

`residues_gap` is the database LOSAT rejects (a stored `-`). The other databases of the rows are the ones of
`LOSAT/tests/seqdb/db` (see `LOSAT/tests/seqdb/README.md`). A rebuild changes the dates stored in the index
files, so the frozen outputs (`Posted date:` lines) must then be frozen again:
`python3 LOSAT/tests/db_search_fixtures.py freeze --ncbi-bin DIR`.
