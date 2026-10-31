# Inputs of the `-db` fixtures

Used by `LOSAT/tests/db_search_fixtures.py` and `LOSAT/tests/db_search_manifest.tsv` (rows that run
`losat blastp -db ...`, `losat blastn -db ...`, `losat tblastn -db ...`, `losat blastx -db ...` or `losat tblastx -db ...` against frozen NCBI BLAST+ 2.17.0
`blastp -db ...`, `blastn -db ...`, `tblastn -db ...`, `blastx -db ...` or `tblastx -db ... -num_threads 1` bytes; the manifest's `program` column says which). NCBI is the
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
- `gapmix.faa` (DI-12): six proteins of the CDS translations of Salmonella Typhi CT18 (AL513382.1, from the
  drop-in survey's Roary data, the first six of at least 120 residues without `-`) with `-`: one inserted in
  the middle (`ins1_pA`), one replacing a residue (`rep1_pA`), a run of three (`ins3_pB`), one at the
  end (`end1_pC`), two at the start (`lead2_pC`), at the start, in the middle and at the end (`sme_pD`), one every
  15 residues (`every15_pE`), none (`ctrl_pF`), and records of `-` only (`allgap5`, `onegap`).
  `gapmix_q.faa`: the six proteins without `-`; `gapmix_q.fna`: the same back-translated with one codon per
  letter (the BLASTX queries).
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
    makeblastdb -in tests/db_search/gapmix.faa -dbtype prot -out tests/db_search/gapmix \
        -title "gap character - at the start, in the middle, at the end, alone" -blastdb_version 4
    makeblastdb -in tests/fasta/PajaWSV.faa -dbtype prot -out tests/db_search/paja \
        -title "PajaWSV proteins" -blastdb_version 4
    makeblastdb -in tests/db_search/o_stored.faa -dbtype prot -out tests/db_search/o_stored \
        -title "stored O against the query's X" -blastdb_version 4
    makeblastdb -in tests/db_search/vol2.faa -dbtype prot -out tests/db_search/vol2_nophr \
        -title "two volumes, the second without its header file" -blastdb_version 4 -max_file_sz 128B
    rm tests/db_search/vol2_nophr.01.phr

`wm_s60k` (rows `db.wm.*`, line NP M2: `-window_masker_db`) is the subject `s60k.fa` of
`LOSAT/tests/blastn_regression_fixtures.py`, with only `.nin .nsq .nhr` kept:

    makeblastdb -in tests/fixtures/blastn_regression/inputs/s60k.fa -dbtype nucl -out tests/db_search/wm_s60k \
        -title "s60k Sakai 2000001-2060000 mutated" -blastdb_version 4

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
- `nslice_qa.fna`: `q_a` alone (one query, for `-query_loc` and `-subject_loc`; rows `db.tab.bn.qloc.6`, `tab.bn.loc.6`).
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

TBLASTN tblastn-fast inputs (rows `tfast.*` of `outfmt0_manifest.tsv` and `db.tfast.*`, line DI4: AMRFinderPlus runs
`tblastn -subject GENOME -query AMRProt.fa -task tblastn-fast -word_size 5 -window_size 15 -threshold 100 ...`),
written by `make_tfast_amr_inputs.py` from the nucleotide sequences of ABRicate's `ncbi` database (NCBI's AMR
reference genes, `db/ncbi/sequences`; the accessions are in the files):

- `tfast_amr.fna`: two contigs of seeded random spacer DNA with ten AMR genes: forward (TEM-1, CTX-M-15, Sul1,
  Tet(A)), reverse-complemented (KPC-2, NDM-1, QnrS1), with about 5 % substituted bases (OXA-48 forward, Erm(B)
  reverse-complemented), the first 60 % of VanA, and a run of 40 N. The database `tfast_amr` is these contigs.
- `tfast_amr_q.faa`: thirteen proteins (translations): TEM-116 and CTX-M-14 (other alleles of the families in the
  contigs), the eight others in the contigs, and SHV-12 (related to TEM), OXA-1 and APH(3')-Ia (not in them).

Run from `LOSAT/`:

    python3 tests/db_search/make_tfast_amr_inputs.py ABRICATE/db/ncbi/sequences tests/db_search
    mv tests/db_search/tfast_amr_s.fna tests/db_search/tfast_amr.fna
    makeblastdb -in tests/db_search/tfast_amr.fna -dbtype nucl -out tests/db_search/tfast_amr \
        -title "synthetic contigs with AMR genes" -blastdb_version 4
    rm tests/db_search/tfast_amr.njs

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

Query-coverage inputs of the ungapped searches (row `db.blastx.ungapped_qcov60_mts1.6`, and the TBLASTX
`-subject` fixtures `qcov.tblastx.*` of `LOSAT/tests/outfmt0_manifest.tsv`), written by `make_qcov_inputs.py`
(run from `LOSAT/`; random letters from `random.Random(20261019)`):

- `qcovp_q.fna`: the query, 600 random bases.
- `qcovp.fna`: `s_part` (the query's first 330 bases unchanged between random flanks: 55 % of the query and the
  best e-value) and `s_full` (the whole query with 25 % substitutions: a weaker e-value).
- `qcovp.faa`: the translations of the same two (`p_part`, `p_full`; `p_full` has its own substitutions), the
  records of the database `qcovp`.
- With `-max_target_seqs 1 -qcov_hsp_perc 60` NCBI reports `s_full` (`p_full` with `blastx -ungapped
  -comp_based_stats 0`): an ungapped search removes the HSPs below the coverage before the preliminary hit list,
  which has the final size, is filled (`blast_engine.c:1538`). Without the coverage it reports `s_part` alone.

The database (version 4, no `-parse_seqids`; only `.pin .psq .phr` are kept, the `.pjs` file was deleted). Run from
`LOSAT/`:

    python3 tests/db_search/make_qcov_inputs.py
    makeblastdb -in tests/db_search/qcovp.faa -dbtype prot -out tests/db_search/qcovp \
        -title "query coverage, ungapped" -blastdb_version 4
    rm tests/db_search/qcovp.pjs

`bxbest` (row `db.blastx.best_hsp.0`) is the database of `tests/fasta/outfmt0/bx_best_hsp_subject.faa`, the Ensembl
R64-1-1 peptide YOL075C, searched with `tests/fasta/outfmt0/bx_best_hsp_query.fna`, the Ensembl BDGP6.54 cDNA
FBtr0100310 (both copied unchanged from the release files; the same pair is the `-subject` fixture `best_hsp.blastx` of
`LOSAT/tests/outfmt0_manifest.tsv`). The subject's first HSP by E-value is a sum-statistics pair of 119 and 45.4 bits
(Expect(2) = 8e-38) and its third has 131 bits (5e-35), so the description table shows 131 and 5e-35. Version 4, no
`-parse_seqids`, `.pjs` deleted; run from `LOSAT/`:

    makeblastdb -in tests/fasta/outfmt0/bx_best_hsp_subject.faa -dbtype prot -out tests/db_search/bxbest \
        -title "BLASTX best HSP: yeast YOL075C" -blastdb_version 4
    rm tests/db_search/bxbest.pjs

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

`residues_gap` and `gapmix` hold stored `-` (NCBIstdaa 0), which NCBI searches as a residue that scores
`BLAST_SCORE_MIN` against every letter and its reports show as the letter of a gap (rows `db.residues_gap.*`,
`db.gap.*`). The other databases of the rows are the ones of
`LOSAT/tests/seqdb/db` (see `LOSAT/tests/seqdb/README.md`). A rebuild changes the dates stored in the index
files, so the frozen outputs (`Posted date:` lines) must then be frozen again:
`python3 LOSAT/tests/db_search_fixtures.py freeze --ncbi-bin DIR`.

OF5a (XML2 and JSON, outfmt 15 and 16): `x16_taxid.faa` (the two `-subject` records of the XML2 probe sprot_tax, `st1` and
`st2`) and `x16_taxid_q.faa` (its query), and `x16_taxid`, built from them with a taxid in each definition line (only `.pin
.psq .phr` are kept):

    makeblastdb -in tests/db_search/x16_taxid.faa -dbtype prot -out tests/db_search/x16_taxid \
        -title "x16 taxid 9606" -blastdb_version 4 -taxid 9606

NCBI writes the taxid in the `Hit` description of outfmt 13-16 (`HitDescr_taxid`); the rows `db.sa.x16_taxid.15` and
`db.sa.x16_taxid.16` freeze that, with no `taxdb.bti` on the search path.

OF5a part 1b (taxonomy names in outfmt 15 and 16): `x16_taxid_99990` and `x16_taxid_7777` are `x16_taxid` built with
`-taxid 99990` and `-taxid 7777` (same command); `x16_taxid_q2.faa` is `x16_taxid_q.faa` with a third query (`qp3`, hits
both subjects) so that a subject is a `Hit` of two queries. The rows `db.sa.x16_taxid.<case>.{15,16}` set
`BLASTDB=tests/seqdb/taxdb/<variant>` in `extra_args` (a `BLASTDB=DIR` token is the environment variable of both programs,
not an argument): `good`, `bad_magic`, `count_larger`, `no_btd`, `empty_fields`, `offset_last_past_end` and
`tests/seqdb/taxdb_xml2/no_9606` (`make_taxdb.py write` of `taxdb.tsv` without the row of 9606). The rows `.nonascii.*`
(`high_bytes_utf8`) are LOSAT refusals: NCBI aborts its XML2 writer there.

OF5a part 2 (XML2 and JSON of BLASTN): `x16_ntitles.fna` (three records: an empty definition line, a title with an
ampersand, angle brackets and quotes, and `.`) and `x16_ntaxid.fna` with its query `x16_ntaxid_q.fna` (`qn1` hits both
subjects, `qn2` none, `qn3` is a slice of `qn1`), and the databases built from them (only `.nhr .nin .nsq` are kept):

    makeblastdb -in tests/db_search/x16_ntitles.fna -dbtype nucl -out tests/db_search/x16_ntitles \
        -title "x16 ntitles" -blastdb_version 4
    makeblastdb -in tests/db_search/x16_ntaxid.fna -dbtype nucl -out tests/db_search/x16_ntaxid \
        -title "x16 ntaxid 9606" -blastdb_version 4 -taxid 9606
    makeblastdb -in tests/db_search/x16_ntaxid.fna -dbtype nucl -out tests/db_search/x16_ntaxid_99990 \
        -title "x16 ntaxid 99990" -blastdb_version 4 -taxid 99990

The rows `db.sa.x16_ntitles.{15,16}` and `db.sa.x16_ntaxid[.<case>].{15,16}` are the BLASTN counterparts of the BLASTP rows of
part 1 and 1b (`BLASTDB=DIR` tokens as there); `db.sa.x16_ntaxid.lastoff.16` freezes the messages that NCBI's BLASTN
writes (`Error: [blastn] ...`).

- `bpqb_q.faa`: 30 records of `tests/fasta/MeenMJNV.faa` of 280-560 letters (11,940 letters, so NCBI searches
  them in two query batches of its default size; LOSAT's block search takes one block at the default
  `LOSAT_X_BPQBLOCK` of 50,000 letters and two with a value of at most 10,000, which
  `tests/blastp_query_blocks_fmt.rs` runs against these frozen bytes), the first 30 of that length in file order. Rows `db.sa.bpqb.{5,8,9,12}`: `blastp -db paja -evalue 0.1
  -max_target_seqs 2 -max_hsps 1` (9 rows of outfmt 6: some queries have hits, some do not).

## -culling_limit rows (line NP, O3b)

Rows `db.blastp.culling_*` and `db.tblastn.culling_*` search `paja` with `bpqb_q.faa` (BLASTP) and `tslice` with
`tslice_q.faa` (TBLASTN) with `-culling_limit`. The split rows `db.{blastp,tblastn}.culling_split_*` use two databases
of the `culling_dom_s.{faa,fna}` subjects in `tests/fasta/outfmt0/` (45 protein records; 25 nucleotide records, the
`-subject` rows `culling.{blastp,tblastn}.*` of `outfmt0_manifest.tsv` search the same files) and the queries
`culling_dom_q.faa` / `culling_domq_q.faa`; their `CHUNK_SIZE=300 OVERLAP_CHUNK_SIZE=50` tokens are environment
variables of both programs (a query is searched in chunks of 300 letters with an overlap of 50; each chunk is
culled by itself). Version 4, no `-parse_seqids`; only `.pin .psq .phr` / `.nin .nsq .nhr` are kept. Run from `LOSAT/`:

    makeblastdb -in tests/fasta/outfmt0/culling_dom_s.faa -dbtype prot -out tests/db_search/culldom \
        -title "culling dom" -blastdb_version 4
    makeblastdb -in tests/fasta/outfmt0/culling_dom_s.fna -dbtype nucl -out tests/db_search/culldomn \
        -title "culling dom nucl" -blastdb_version 4
    rm tests/db_search/culldom.pjs tests/db_search/culldomn.njs
- OF-BPQB2 (outfmt 15 and 16 in the same block search): rows `db.sa.bpqb.15` and `db.sa.bpqb.16`, the same search as
  `db.sa.bpqb.8` (65 KB at most). `paja_taxid` is `paja` built with a taxid (`make_paja_taxid.py MAKEBLASTDB OUTDIR`: `makeblastdb
  -in tests/fasta/PajaWSV.faa -dbtype prot -blastdb_version 4 -taxid 10090`, then the taxid of BDT63547.1 patched to 9606 in
  the header file, the same length); BDT63547.1 is one of the two subjects that only the second NCBI query batch of `bpqb_q.faa`
  hits, so the row `db.sa.bpqb.taxdb.16` (`BLASTDB=tests/seqdb/taxdb/good`) writes `Mus musculus` in the first block and
  `Homo sapiens` in the second. `tests/seqdb/taxdb_xml2/last_10090_past_end` is `make_taxdb.py all` of a two-entry TSV (9606, 10090),
  variant `offset_last_past_end` (the last entry, 10090, has its offset past the end of `taxdb.btd`); the row
  `db.sa.bpqb.taxdb.lastoff.15` freezes NCBI's messages across the two query batches (the volume's lookup once, the `Hit`'s each time).

OF6 (SAM, outfmt 17, BLASTN): the rows `db.sa.sam.*` (stems `mega`, `dc`, `tb`, `mv`, `short`, `ambrc`, `mts1`, `maxhsps1`,
`cull1`, `qloc`, `nohit`, `b100`; suffix `.17p` = `-outfmt 17`, `.17q` = `17 SQ`, `.17r` = `17 SR`, `.17b` = `17 SQ SR`, `.5` = the
outfmt 5 sibling that `S/scripts/check_sam_layout.py` predicts the document from) search `nslice` with `nslice_q.fna`
(`mega`, `dc`, `tb`, `short`, the option rows and `b100`: `BATCH_SIZE=100`, one query per batch), `nucl_mv` with
`tests/seqdb/nucl.fna` (`mv`) and `sam_ambrc` with `tests/fasta/outfmt0/sam_ambig_q.fna` (`ambrc`: ambiguity letters in
the stored data of a minus-strand hit). `sam_ambrc` is built from `tests/fasta/outfmt0/sam_ambig_rc_s.fna` (only
`.nhr .nin .nsq` are kept):

    makeblastdb -in tests/fasta/outfmt0/sam_ambig_rc_s.fna -dbtype nucl -out tests/db_search/sam_ambrc \
        -title "sam ambiguity letters on a minus strand" -blastdb_version 4

The frozen documents are NCBI's bytes, including the `CL:` tag of the `@PG` line (the path of the NCBI binary that
froze them); the check removes that tag from both sides (`db_search_fixtures.py strip_sam_cl`).

OF5b part 1 (XML2 and JSON of BLASTX): `make_x16_blastx_inputs.py` writes the inputs (run from `LOSAT/`):
`x16_bx_taxid_q.fna` (three nucleotide queries for `x16_taxid`: a minus frame, a plus frame, one without hits) and
`x16_bx_taxid_q2.fna` (two queries of 10,003 bases, each its own NCBI query batch of 10,002 bases, that both hit
`x16_taxid_99990`), `tests/fasta/outfmt0/x16_bx_q_unsearched.fna` (a 10,003-base slice with lower-case stretches, then a
300-N query and a 2-base query: a batch NCBI does not search) and, in `tests/fixtures/outfmt16`, `qnx_bad.fna`,
`qnx_firstbad.fna`, `qnx_late.fna`. BLASTX does not read `BATCH_SIZE`, hence the padding to 10,003 bases. The rows
`db.x16.blastx.*` (bslice) and `db.x16_bx_taxid.*` (taxid and `BLASTDB=` variants, as in the BLASTP rows) are frozen from NCBI
BLAST+ 2.17.0; `tests/fixtures/outfmt16/db/sprot_tax` is `sprot.faa` built with `makeblastdb -dbtype prot -blastdb_version 4
-taxid 9606` (only `.phr .pin .psq` kept).

OF5b part 1 fix (DECISIONS 32, BLASTX `-db` 15/16 with a taxonomy database): `tests/fixtures/outfmt16/db/sprot_mixtax.pal` lists
`sprot_tax10090` (`sprot.faa`, `makeblastdb -dbtype prot -blastdb_version 4 -taxid 10090`) and `other_tax9606` (`other.faa`, one protein
nothing hits, `-taxid 9606`; only `.phr .pin .psq` kept). `mixtax_blastx_{good,badname}_{15,16}` are NCBI BLAST+ 2.17.0
`blastx -query qnx_late.fna -db db/sprot_mixtax -num_threads 1` with `BLASTDB=tests/seqdb/taxdb/good` and `.../high_bytes_utf8` (the name of
9606 has bytes >= 0x80; 9606 is in the database and never reported, so NCBI writes the same bytes as with `good`).

- `x16_ttaxid_q.faa`, `x16_ttaxid_qwarn.faa`, `x16_ttaxid_qwarn2.faa` (OF5b part 2, TBLASTN outfmt 15 and 16): protein queries
  translated from `sn1` (frame +1) and `sn2` (frame -3) of `x16_ntaxid.fna`, and an unrelated one without hits; the `qwarn` files
  add a title that ends in 60 amino-acid letters (the reader's warning) and the `qwarn2` order puts the query without hits
  first. They search the existing nucleotide databases `x16_ntaxid` (taxid 9606) and `x16_ntaxid_99990`: rows
  `db.sa.x16_ttaxid.*` (`BLASTDB=DIR` and `BATCH_SIZE=N` tokens as in the BLASTP and BLASTN rows; the messages on standard
  error name `tblastn`).

- `x16_xtaxid_q.fna`, `x16_xtaxid_qwarn.fna`, `x16_xtaxid_qwarn2.fna` (OF5b part 3, TBLASTX outfmt 15 and 16): nucleotide queries made
  from `sn1` (its first 150 bases, `xq1`) and `sn2` (the reverse complement of its first 150 bases, `xq2`) of `x16_ntaxid.fna`, and an
  unrelated one without hits (`xq3`); the `qwarn` files add a title that ends in 28 nucleotide letters (the reader's warning) to `xq2`
  and the `qwarn2` order puts the query without hits first. With `-evalue 1e-20` each query hits its own subject only. They search the
  existing databases `x16_ntaxid` (taxid 9606) and `x16_ntaxid_99990`: rows `db.sa.x16_xtaxid.*` (`BLASTDB=DIR` and `BATCH_SIZE=N`
  tokens as in the BLASTP, BLASTN and TBLASTN rows, the other tokens are options; the messages on standard error name `tblastx`).

## Databases with parsed Seq-ids (line DBID, ID2a)

`ids/`: rows `db.id.*`. `ids/make_ids.py` (its docstring has the usage and the records) writes
`ids/prot_v5` (version 5) and `ids/prot_v4` (version 4) with `makeblastdb -parse_seqids` from the records of
`LOSAT/tests/seqid/fasta` and `ids/gen.faa` (titles that NCBI makes from the Seq-ids), with every file
makeblastdb writes: NCBI's formatters find a reported sequence again through the ISAM or LMDB files, so the
frozen outputs need them (LOSAT reads the index, sequence and header files). Queries: `ids/ids_q.faa` (three
proteins: one hits every sequence, one hits a few weakly, one has no hit), `ids/ids_q.fna` (the same
back-translated, for BLASTX) and `ids/ids_qw.faa` (thirty W, which hits the records of
`LOSAT/tests/seqid/db/variants` that no NCBI formatter finds again). `ids/refuse_v5` (version 5, from
`ids/refuse.faa` and `ids/taxmap_refuse.txt`) holds one record per case of the rows `db.id.*.taxdb_nf.*` and
of LOSAT's refusals, each hit by its own query `ids/refuse_<case>.faa` (`.fna` for BLASTX) alone with
`-evalue 1e-20`: NCBI's outfmt 0 does not end on the records `long` and `punct`, so no other query may
report them.

Nucleotides (line DBID, ID2b; rows `db.id.blastn.*`, `db.id.tblastn.*`, `db.id.tblastx.*`): `make_ids.py` also
writes `ids/nucl_v5` and `ids/nucl_v4` from `ids/nucl.fna` and `ids/taxmap_nucl.txt` (the Seq-id forms of
`LOSAT/tests/seqid/fasta/nucl_ids.fa`, every record with its sequence, and records whose titles NCBI makes
from their Seq-ids), `ids/nrefuse_v5` from `ids/nrefuse.fna` (`lcl|777` before `gi|777`, which NCBI's
formatters do not find again) and the two volumes `ids/njoined_a` and `ids/njoined_b` of the alias
`ids/njoined.nal` (`gb|AF000040|` and `gb|AF000040.1|`, whose hits NCBI's XML and outfmt 0 join). Queries:
`LOSAT/tests/seqid/fasta/qnucl.fa` (BLASTN, TBLASTX) and `qprot.fa` (TBLASTN); the second query of each has
no hit. TBLASTN and TBLASTX run with `-max_hsps 1` to keep the frozen files small.

## The taxonomy database's prefetch messages (line DBID, item T, review 1)

Rows `db.tm.bp.fast.*`, `db.tm.bp.2vol.*` and `db.tm.*.cache*`, all with `BLASTDB=tests/seqdb/taxdb/offset_last_past_end`
(the last entry, taxid 99990, begins past the end of `taxdb.btd`):

- `tm_fast_q289.faa`: 289 copies of `qp1` of `x16_taxid_q.faa`, the first with its first letter replaced by `O` (a warning
  of its own). With `-task blastp-fast` against `x16_taxid_99990` the search makes two query blocks; the offset error
  must come before that warning, at any thread count.
- `tm_2vol.pal` (`DBLIST tm_2vol_a tm_2vol_b`): two volumes that are copies of `x16_taxid_99990` (`.phr .pin .psq`); NCBI
  writes the offset error once per volume.
- `tm_cache.faa` (`st1`, `st2`, `st3`: the same protein), `tm_cache_q.faa` (its query), and the databases built from them
  with every file makeblastdb writes (version 5):

      makeblastdb -in tm_cache.faa -dbtype prot -out tm_cache150 -title "tm cache 150" \
          -parse_seqids -taxid_map tm_cache150.map
      makeblastdb -in tm_cache.faa -dbtype prot -out tm_cache151 -title "tm cache 151" \
          -parse_seqids -taxid_map tm_cache151.map

  (run in this directory). `tm_cache150` holds the taxids 99990, 150, 99990: 150 takes the slot of 99990 (`99990 & 255`) in
  NCBI's 256-slot taxonomy cache of a volume (`CSeqDBVol::m_TaxCache`), so NCBI writes the offset error again after a fetch
  of 150, as often as its order of fetches makes it (twice for the row's query, once for the same query with a version 4
  database or an absolute `BLASTDB`); LOSAT refuses the run. `tm_cache151` (151 instead of 150) is the control: once.

## Taxid filtering (DI-8a)

Rows `db.tf.*`: `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and `-no_taxid_expansion` with `-db`, five
programs. `taxfilter/` holds the inputs and what `make_taxfilter.py MAKEBLASTDB TAXONOMY4BLAST_SQLITE` builds from them
(makeblastdb 2.17.0, version 5, `-parse_seqids` with a taxid map; `.pjs` and `.njs` dropped):

- `tfp.faa` and `tfp.map`: 23 Swiss-Prot globins of mammals, birds, reptiles, amphibians and fishes (the HBA_HUMAN hits of
  the DI-8 probe); the records of P69905.2, P01923.1, Q9TS35.2 and P63107.2 carry a second definition line (control-A) of a
  bird or another mammal. Database `tfp` (one volume).
- `tfn.fna` and `tfn.map`: 16 globin cDNAs of human (Ensembl GRCh38), zebrafish and fly; human HBA1 (ENST00000320868.9) carries
  a second definition line of a zebrafish transcript, zebrafish hbaa1 (ENSDART00000113098.5) one of a fly transcript.
  Database `tfn` in two volumes (`-max_file_sz 4KB`: OIDs 0-11 and 12-15) under one LMDB.
- `taxonomy4blast.sqlite3`: the 226 rows of NCBI's taxonomy file (the copy distributed with swissprot) for the taxids of both
  maps and all their ancestors, in a new file with the source's schema text. NCBI reads the table `TaxidInfo` and its index
  only, so a taxid's descendants within the databases are those of the whole file. The rows reach it by `BLASTDB=tests/db_search/taxfilter`;
  the `nosql` rows leave it out of reach (NCBI's warning, no restriction).
- `root/taxonomy4blast.sqlite3`: the same file without the rows of 7227 (fly) and 8237 (tuna). The rows `db.tf.bp.root.6` and
  `db.tf.tbn.root.6` freeze NCBI's `-taxids 131567` (the root's only child there) with it, and `tests/db_taxid_root.rs`
  checks that LOSAT's `-taxids 1` writes those bytes (taxid 1 is its own parent: NCBI's expansion of it does not end;
  DI-8 DECISIONS 10).
- `q.faa` (HBA_HUMAN) and `q.fna` (human HBA1 cDNA): the queries. `tl_crlf.txt` (40674 and 8782 with CRLF line ends, an
  empty line and blanks) and `tl_neg.txt` (9443): the taxid lists.

The `v4` rows search `x16_taxid` and `x16_ntaxid` (version 4 databases with taxids, above).

## Seqid lists (DI-8b)

Rows `db.sa8.*`: the `ext` `use_this_seqid` of a positive list on the Seq-align formats and what it limits, on
`taxfilter/tfp`, `taxfilter/tfn` and `useseqid/`. `useseqid/` holds `p5m`, `p5mt` (protein, version 5, the records of
`p_m.fa`: a text id and a local id on one OID, a GI that is the Seq-align id, a local id first, a PRF name, PDB ids) and
`n5m` (nucleotide, `n_m.fa`), written by NCBI BLAST+ 2.17.0:
`makeblastdb -in p_m.fa -dbtype prot -out p5m -parse_seqids -blastdb_version 5 -taxid 9606 -title p5m`,
`makeblastdb -in p_m.fa -dbtype prot -out p5mt -parse_seqids -blastdb_version 5 -taxid_map p_m.map -title p5mt`,
`makeblastdb -in n_m.fa -dbtype nucl -out n5m -parse_seqids -blastdb_version 5 -taxid 562 -title n5m`;
`q.faa`, `q.fna`, `qx.fna` are the queries and `sl_m.txt` a seqid list.

Rows `db.sl.*`: `-seqidlist` and `-negative_seqidlist` with `-db`, five programs, on `taxfilter/tfp`, `taxfilter/tfn`
(two volumes), `ids/prot_v5`, `ids/nucl_v5` (version 5) and `ids/prot_v4`, `ids/nucl_v4` (version 4, string ISAM).
`seqidlist/` holds the lists, written by `seqidlist/make_seqidlist.py`:

- text lists in the forms NCBI's reader takes (white-space separated tokens, a leading `>` skipped, `#` to the end of the
  line a comment, CRLF lines): `sl_tfp.txt` (five Swiss-Prot ids of tfp, with and without versions, one in FASTA form, and an
  unknown id), `sl_tfp_neg.txt` (the same ids without versions: a versionless id excludes the versioned stored id),
  `sl_tfn.txt` (four local ids across tfn's two volumes), `sl_prot.txt` and `sl_nucl.txt` (every Seq-id form of
  `ids/prot_v5` and `ids/nucl_v5`: tpg, RefSeq and Swiss-Prot without versions, a pdb chain, local ids, a general id as
  `db:tag`, a versionless accession of an OID of several definition lines, a gi that NCBI drops), `sl_gi.txt` (gi ids only:
  NCBI's `Empty seqid list` warning and no restriction on a version 5 database, no OID on a version 4 one),
  `sl_comment.txt` (no id: no list), `sl_none.txt` (ids no database has: NCBI searches an empty database, LOSAT refuses),
  `sl_empty.txt` (0 bytes: NCBI's `Failed to map seqidlist file ` error);
- binary lists in the format of `blastdb_aliastool -seqid_file_out` (NUL, file size, count, title, date, database length
  and, when it is not 0, the database's date and volume names, then the ids with a one-byte length): `sl_tfp.bsl` and
  `sl_nucl.bsl` (the processed ids of the text lists: no warning), `sl_tfp_badlen.bsl` (a database length that is not
  tfp's: NCBI's `Seqidlist file db info does not match input db` warning, the list is used).

The negative lists on `ids/prot_v5`, `ids/prot_v4`, `ids/nucl_v5` and `ids/nucl_v4` (rows `db.sl.bp.nforms.6`,
`db.sl.bp.v4neg.6`, `db.sl.bp.v4gi.6`, `db.sl.bn.nforms.6`, `db.sl.bn.v4neg.6`) are `losat-rejection` rows: those databases
hold sequences whose every definition line has only a gi or a ti, and a negative seqid list leaves such a sequence without a
line (`CSeqDBNegativeList::FindId` matches no gi or ti): NCBI searches it and drops its hits from the report without notice,
LOSAT refuses when the sequence is to be reported (`the seqid list keeps sequence N of the volume but none of its definition
lines`, DI-8b DECISIONS 17). A positive list that selects such a sequence is refused before the search.

`ids/noacc_v5` (three proteins of `ids/ids_q.faa`, `makeblastdb -dbtype prot -blastdb_version 5` without `-parse_seqids`,
NCBI 2.17.0) has no accession table: a seqid list on it is NCBI's `BLAST Database error: DB contains no accession info.`
(row `db.sl.bp.noacc.6`, rc 2).

A list file that is not on the search path (NCBI: `Error: NCBI C++ Exception: ... To be memory mapped the file must
exist: ''`, exit 255, with the path of NCBI's build tree) is checked by `tests/cli_v2.rs`, not by a row.

## Ungapped nucleotide databases (line NP, O4a-1)

`ungapped_abr` (the 33 genes of `LOSAT/tests/fixtures/blastn_regression/inputs/culling_abricate_genes.fa`, as `ungapped_abr.fna`) is built like
the mlst databases (parsed Seq-ids and a hash index), `ungapped_ambrun` (the subject `ambrun_s.fa`, a stretch with 30 `N` and one `R` inside a match)
without, `ungapped_re` (the six subjects of `re_s.fa`: N every other base, a run of N, IUPAC codes, a minus-strand copy) without; rows
`db.blastn.ungapped*` (`blastn -ungapped`). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/ungapped_abr.fna -dbtype nucl -out tests/db_search/ungapped_abr \
        -parse_seqids -hash_index -title PubMLST -blastdb_version 4
    makeblastdb -in tests/db_search/ungapped_ambrun.fna -dbtype nucl -out tests/db_search/ungapped_ambrun \
        -title "ambiguity run" -blastdb_version 4
    makeblastdb -in tests/db_search/ungapped_re.fna -dbtype nucl -out tests/db_search/ungapped_re \
        -title "ungapped re-evaluation" -blastdb_version 4
    rm tests/db_search/ungapped_abr.njs tests/db_search/ungapped_ambrun.njs tests/db_search/ungapped_re.njs

`mejo` is `LOSAT/tests/fasta/MejoMJNV.fasta` (one sequence, 401,182 bases) as pyani ANIb builds its databases
(`makeblastdb -dbtype nucl -in FILE -title STEM -out OUT`, no parsed Seq-ids; version 4 for the repository), the
subject of rows `db.xdrop.*` (DI-6: the X-drops, `-window_size`, `-off_diagonal_range` and `-no_greedy`) and of some of
rows `db.minraw.*` (DI-6b: `-min_raw_gapped_score`; the queries
`tests/fasta/outfmt0/anib_frags*.fna` are 1020-nt fragments of `MelaMJNV.fasta`, pyani's fragment size). Run from `LOSAT/`:

    losat makeblastdb -dbtype nucl -in tests/fasta/MejoMJNV.fasta -title MejoMJNV -out tests/db_search/mejo -blastdb_version 4
    rm tests/db_search/mejo.njs

The rows `db.blastn.ungapped_rej.*` (outfmt 8, 11, 17 with `-ungapped`) are refusals until O4a-2.
