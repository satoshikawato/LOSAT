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
- `o_stored_bx_q.fna`: `s_x` back-translated with its twelve X as `NNN` codons, the BLASTX query of `db.o_stored.bx.8`
  (BLASTX reads the stored O of `s_o` as stored, 154 identities; `-subject o_stored.faa` reads it as X, 166, row
  `pyl.blastx.ostored.8` of `outfmt0_manifest.tsv`).
- `bx_title_stop_q.fna`: two records with a title that ends with 24 nucleotide letters; the first has a
  hyphen and an invalid letter (`7`) in its first data line, the second starts with a line that is not plausible
  data, the BLASTX query of `db.bxtitle.stop.6` (NCBI writes the first record's title warning after its data
  warnings, then stops; the second record's title warning is never posted).
- `bx_title_later_q.fna`: the same first record; the second has a data line `EF-` (no nucleotide letter: a hyphen
  and two invalid letters, with their warnings), then a line that is not plausible data, which NCBI still checks
  because the record has no residues yet (`db.bxtitle.later.6`).
- `bx_gap_first_q.fna`: a gap line `>?5` and 60 nucleotides, then the record `r1`, the BLASTX query of
  `db.bxgap.first.6` (with `-parse_deflines` NCBI stops at line 1: `Input doesn't start with a defline or
  comment`) and of the control `db.bxgap.first.nopd.6` (without it the gap line starts a record without a
  defline). `bx_gap_comment_q.fna` is the same after a comment, a blank line and a `#` line
  (`db.bxgap.comment.6`: NCBI names line 4), and `bx_gap_bad_q.fna` starts with `>?abc` instead
  (`db.bxgap.bad.6`: NCBI stops at line 1 before it parses the gap).
- `bx_prescan_data_q.fna`, `bx_prescan_gap_q.fna`: a data line, or the gap line `>?5`, before any defline,
  then the record `r1` whose title is not US-ASCII (`café x`), the BLASTX queries of `db.bxprescan.*`: with
  `-parse_deflines` NCBI stops at line 1 (`Input doesn't start with a defline or comment`) in every format,
  before any title is written. `bx_prescan_badgi_q.fna`: the record `r0` with that title, then `>gi|abc`, which
  NCBI cannot read as a Seq-id under `-parse_deflines` (`Could not construct seq-id from 'gi|abc'`): in one batch
  the error comes before `r0` is written, with `BATCH_SIZE=1` after; with `BATCH_SIZE=0` (`db.bxprescan.bs0.*`, without
  `-parse_deflines`) NCBI reads no query and stops with `Empty CBlastQueryVector`. `bx_prescan_badgi_big_q.fna`: the same with
  10,200 nucleotides in `r0`, more than BLASTX's batch of 10,002 letters (`db.bxprescan.badgi.big.5`,
  `BATCH_SIZE=1000000`: one batch).
- The rows `db.bxgapline.*` search `bslice` with `LOSAT/tests/fasta/outfmt0/bx_gapline_q.fna`, the query of the
  `-subject` rows `gapline.blastx.*` of `outfmt0_manifest.tsv`: gap lines of bad sizes (`>?abc` before any
  defline, `>?0`, `>?unk`, `>?99999999999`, `>?-5` as a record's last line), text after the size and unknown,
  conflicting or unclosed modifiers, which NCBI reads as `ParseGapLine` does (a gap of 1 after `Bad gap size`,
  warnings for the rest).
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
    makeblastdb -in tests/db_search/bxchunk.faa -dbtype prot -out tests/db_search/bxchunk \
        -title "proteins planted in bxchunk_q.fna" -blastdb_version 4

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
- `bxchunk_q.fna`: one query of 25,001 random bases (`random.Random(5)`, as the review of #702 made it) with a
  back-translated protein of 150 residues planted 200 bases before each of 3000, 5997, 6000, 9705, 10001, 12000,
  15000, 19410 and 20000 (the chunk borders of several `CHUNK_SIZE` values), and short lower-case runs at those
  points; `bxchunk.faa`: the nine planted proteins with a quarter of their residues changed, for the database
  `bxchunk` (rows `db.bxchunk.*`: BLASTX `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE`). `bxchunk_empty.fna`: a record of a
  defline only (NCBI reports the query splitter's errors before the empty query's `Sequence contains no data`).
- `bxempty_*.fna` (rows `db.bxempty.*`, against `bxchunk`, with `BATCH_SIZE=200000`: LOSAT searches such a batch in
  units of 10002 letters): `bxempty_le.fna` is the query of `bxchunk_q.fna` and an empty record (a defline only);
  `bxempty_lel.fna` adds a second query (its first 12,000 bases) and a second empty record; `bxempty_none.fna` is an
  empty record, 5,001 records of 2 bases (no context to search) and an empty record; `bxempty_none1.fna` the 5,001
  records and one empty record.
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
no hit. TBLASTN and TBLASTX run with `-max_hsps 1` to keep the frozen files small. `ids/nambig_v4` (from
`ids/nambig.fna`: a record with ambiguity data) has no row; `tests/blastn_report_stream.rs` searches it.
Rows `db.id.blastn.stream.q25.*` (line DBID, blastn's report stream of these databases): the 25 records of
`ids/nucl.fna` as queries in batches of `BATCH_SIZE=1500`, each batch with hits.

Seq-align ids by which NCBI's formatters find another sequence (line DBID N1; rows `db.id.*.found.*`):
`make_ids.py` writes, from the FASTA of `ids/make_found_fasta.py`, `ids/pfound_v5` and `ids/pfound_v4` from
`ids/pfound.faa` (five pairs of unrelated proteins, each hit by its query `ids/pfound_q<k>.faa`;
`ids/pfound_q1.fna` and `ids/pfound_q4.fna` for BLASTX), `ids/nfound_v5` and `ids/nfound_v4` from
`ids/nfound.fna` (`gb|AF000050.1|`, then `gb|AF000050|`; queries `qnucl.fa` and `qprot.fa`), and the aliases
`ids/pdup_ab.pal` (two databases with the same record `gb|EEE00005.1|`) and `ids/pdup_ac.pal` (the second
copy with another title). NCBI looks the Seq-align id of a reported sequence up again in the database and
takes the first sequence the lookup gives: by `gb|AAA00001|` it finds `gb|AAA00001.1|` before it, by
`emb|BBB00002.1|` the `gb|BBB00002.1|` before it, by `lcl|CCC00003` the `gb|CCC00003.1|` before it (whose ids
do not compare equal, so NCBI finds no sequence); `gb|DDD00004|` before `gb|DDD00004.1|` and `pdb|2XYZ|` before
`pdb|2XYZ|AA` each find their own. LOSAT refuses such a sequence in the formats that read its Bioseq, unless the
sequence found holds the same definition lines and letters (`pdup_ab`); `-negative_seqidlist
seqidlist/sl_pfound_neg.txt` leaves `gb|AAA00001.1|` out of the search, and NCBI's lookup with it.

The same record with other taxids (line DBID N1b; rows `db.id.*.found.taxid.*`): `make_ids.py --copytaxid-only`
writes `ids/pdup_t` (`pdup.faa` with taxid 9606), `ids/pdup_u` (without one), `ids/pdup_n` and `ids/pdup_o` (both
with taxid 99990, the last entry of the taxonomy databases of `../seqdb/taxdb`), `ids/ndup_a` and `ids/ndup_t`
(`ndup.fna`, the back-translation of `pdup.faa`, without a taxid and with 9606), `ids/pdupe_a` and `ids/pdupe_t`
(`pdupe.faa`, without a title, without a taxid and with 9606), and the aliases `ids/pdup_at.pal`, `ids/pdup_tu.pal`,
`ids/pdup_no.pal`, `ids/pdupe_at.pal` and `ids/ndup_at.nal` (queries `pdup_q.faa` and `ndup_q.fna`). NCBI shows for
both copies the one its lookup gives first, with its taxids: the same rows where no report writes taxids, which
LOSAT writes; LOSAT refuses the sequence in the formats of taxids and their names, for `pdupe_at` with the
taxonomy database (a title made from the taxids), and with a taxonomy database on `BLASTDB` whose texts depend on
the copy NCBI fetches (`bad_magic` for `pdup_at`, `offset_last_past_end` for the copies of 99990 in `pdup_no`).

Copies whose title NCBI can make from the taxids (line DBID N1b L1; rows `db.id.title.*`): `make_ids.py
--titlesrc-only` writes `ids/ndupe_a` and `ids/ndupe_t` (`ndupe.fna`, the letters of `ndup.fna` without a title,
without a taxid and with 9606, through `ids/ndupe_at.nal`), `ids/psubj_a`, `ids/psubj_t` (`psubj.faa`, lcl|Subject_1
without a title) and `ids/psubjb_a`, `ids/psubjb_t` (`psubjb.faa`, lcl|Subject_2 with a title ending in `]`), through
`ids/psubj_at.pal` and `ids/psubjb_at.pal`, and `ids/nsubj_t` (`nsubj.fna`, lcl|Subject_1 without a title, taxid
9606, one database). NCBI generates the title of the copy it fetches, from its taxids where the taxonomy database names them,
only in outfmt 0 and for the tabular word that replaces a `lcl|Subject_` id: LOSAT writes the other formats
(`pdupe_at` and `ndupe_at` with `tests/seqdb/taxdb/good`, `psubj_at`, whose word is `unnamed` either way) and refuses
those two where the word is printed (`psubjb_at`, whose word is taken after NCBI rewrites the end of the title), also
for `nsubj_t` alone with the taxonomy database (NCBI's word is `Homo`, the organism's; not with `6 salltitles`, which
prints no id).

An empty title without parsed Seq-ids (line DBID H3; rows `db.id.bare.*`): `make_ids.py --emptytitle-only` writes
`ids/nbare_t` and `ids/pbare_t` (`nbare.fna` and `pbare.faa`, the letters of `ndup.fna` and `pdup.faa` under `>` alone,
without `-parse_seqids`, taxid 9606). The Bioseq that NCBI's formatters fetch has a source for the taxid, so outfmt 0
makes the title from it where the taxonomy database names the taxid (`tests/seqdb/taxdb/good`: `unnamed protein
product [Homo sapiens]`, `Homo sapiens`): LOSAT refuses that outfmt 0 and writes it without a taxonomy database or with
`no_btd`, and the other formats (XML shows the stored title, the tabular id of a protein is `unnamed` either way).

A title of punctuation without parsed Seq-ids (line DBID tabfields; rows `db.id.tf.*.punct.*`): `make_ids.py
--tabfields-only` writes `ids/ppunct_np` (version 4, without `-parse_seqids`) from `ids/ppunct.faa`, whose first title
`, ,` NCBI's `x_CleanAndCompress` reads past the end of where it makes the tabular id of a protein; the query
`ids/ppunct_q.faa` hits both records. A table that prints no subject id (`6 qseqid evalue`) does not make it, so LOSAT
writes NCBI's rows; `sseqid` is LOSAT's refusal.

Each batch's prefetch counts its own Seq-aligns (line DBID pfbatch; rows `db.id.pfb.*`): `make_ids.py --pfbatch-only`
writes `ids/nordtag_v5` (version 5, `-parse_seqids`) from `ids/nordtag.fna`, a sequence whose first id
`gnl|BL_ORD_ID|0` makeblastdb stores with a text tag, so NCBI's lookup by that Seq-align id fails (as in `pshare_v5`),
and another sequence. The queries `ids/nordtag_q.fna` and `ids/nordtag_q.faa` (qa, qb) hit one each; `ids/pfb_aba.faa`
(qa, qb, qa2: with `BATCH_SIZE=101` the batches [qa, qb] and [qa2]) and `ids/pfb_aab.faa`, `ids/pfb_aab.fna` (qa, qa2,
qb) put the sequence with the text tag in a batch of more Seq-aligns than the limit and in one within it, or only in the
first; `ids/pfb_s.fna` holds two pieces of it (two HSPs gapped, one Seq-align with `-ungapped`). `ids/pfb_rc_nambig.fna`
and `ids/pfb_rc_ids_q.fna`, the reverse complements of `ids/nambig.fna` and `ids/ids_q.fna`, align on the minus strand
of the query with `lcl|777` of `ids/nrefuse_v5`, whose Bioseq NCBI's prefetch does not get. With a positive
`PRE_FETCH_SEQS_LIMIT`, a batch whose Seq-aligns outnumber it reads nothing in its prefetch, so a table that reads
nothing of the subject writes NCBI's rows (NCBI writes no warning); a batch within it reads the subject, and where its
prefetch fails NCBI and LOSAT write `Error pre-fetching sequence data` (line DBID pfwarn).

The warning of a batch whose prefetch fails (line DBID pfwarn; rows `db.id.pfw.*`): `make_ids.py --pfwarn-only` writes
`ids/ntag_t5` (version 5) from `ids/nordtag.fna`, `ids/nlook_t4` (version 4) from `ids/nlook.fna` and `ids/pshare_t5`
(version 5) from `ids/pshare.faa`, every sequence with taxid 99990, whose lookup in
`tests/seqdb/taxdb/offset_last_past_end` writes the offset error. NCBI's prefetch writes the taxonomy database's texts of
the subjects it loads, then `Warning: [<program>] Error pre-fetching sequence data`, once for the batch, before the
batch's query warnings. A minus-strand Seq-align of `bbs|123` (no ISAM key in version 4, so no Bioseq) fails after every
load (`ids/pfw_look_fw_rc.fna`, `ids/pfw_look_rc_fw.fna`: `fw70` on the plus strand of `gb|AF000070.1|`, `rc123` on the
minus strand of `bbs|123`); the lookup of `gnl|BL_ORD_ID|0` with a text tag fails where the loads, in an order of
NCBI's Seq-id handles that LOSAT does not reproduce, reach it, so a batch with it whose other subjects' loads write a
text is LOSAT's refusal. With `-ungapped` BLASTN makes one Seq-align of a subject's HSPs, whose strand is its first
HSP's: `ids/pfw_fwL_rcS.fna` aligns with `bbs|123` on the plus strand first (240 letters) and on the minus strand (60),
`ids/pfw_fwS_rcL.fna` the other way round. `ids/pfw_ps_q21.faa` holds `pshare_q2.faa` and `pshare_q1.faa`.

Copies of one record that follow each other in a hit list (line DBID N3; rows `db.id.join.*`): `make_ids.py
--join-only` writes, from the FASTA of `ids/make_join_fasta.py`, `ids/pjoin_a`, `ids/pjoin_b`, `ids/pjoin_c` and
`ids/njoin_a`, `ids/njoin_b`, `ids/njoin_c` (version 5; r1 in all three, r3 with two definition lines and two
pieces of the query in a and b, r6 the letters of r1 under another id in b), `ids/pjoin_bt` and `ids/njoin_bt` (the b
records with taxid 9606), and the aliases `ids/pjoin_abc.pal`, `ids/pjoin_abt.pal`, `ids/njoin_abc.nal` and
`ids/njoin_abt.nal` (queries `pjoin_q.faa`, `njoin_q.fna` and `njoin_q2.fna`, the last with two records for two
query batches). NCBI orders tied hits by the later ordinal id first: r1 of c, r6, then r1 of b and a, which outfmt 0,
5 and 13 to 16 show in one hit (r3 of b and a too), as LOSAT does; `-num_descriptions`, `-num_alignments` and
psiblast's cut to `-max_target_seqs` count them as one subject, the tabular `qcovs` and `qcovus` take the group's
coverage. `ids/pjoin_x` and `ids/njoin_x` hold the definition lines of r1 and r3 with other letters (no copies,
refused; aliases `ids/pjoin_ax.pal`, `ids/njoin_ax.nal`); `ids/njoin_ta` and `ids/njoin_tb` (alias `ids/njoin_tab.nal`)
hold r10 with taxid 9606 and a second definition line `lcl|loc10`, which the outfmt 0 heading under `-taxids 9606
-no_taxid_expansion` hides where each copy has one HSP.

Sequences whose only Seq-id is `gnl|BL_ORD_ID|n` (line DBID N2a; rows `db.id.*.ord.*`): `make_ids.py` writes, from
the FASTA of `ids/make_ord_fasta.py`, `ids/pord_v5` and `ids/pord_v4` from `ids/pord.faa` (two accessions, then
`gnl|BL_ORD_ID|2`, its own number, `gnl|BL_ORD_ID|6`, the number of no sequence, and `gnl|BL_ORD_ID|1`, the
number of a sequence with an accession; each hit by its query `ids/pord_q<k>.faa`, `.fna` for BLASTX), `ids/pord_x`
from `ids/pord_x.faa`, read after `pord_v5` through the alias `ids/pord_vx.pal`, which lists it first (NCBI orders
the volumes of an alias by their base names; there `gnl|BL_ORD_ID|6` of `pord_v5` is the number of the second sequence of `pord_x`, whose
`gnl|BL_ORD_ID|1` is its own number in its volume), and `ids/nord_v5` and `ids/nord_v4` from `ids/nord.fna`
(`gnl|BL_ORD_ID|3` first, the number of no sequence, then `gnl|BL_ORD_ID|1` and an accession; queries
`ids/nord_q<k>.fna` and `.faa`, and `ids/nord_qall.fna` with the three in reverse order). NCBI looks n up
plus the first ordinal id of the volume, so its own number finds the sequence, which LOSAT reads as one without
parsed Seq-ids; another n finds no sequence or another one, which LOSAT refuses in the formats that read the
Bioseq. No record names `gnl|BL_ORD_ID|0`: makeblastdb stores it with a string tag.

Seq-align ids that NCBI's formatters look up like any other (line DBID N2b; rows `db.id.*.look.*`): `make_ids.py`
writes, from the FASTA of `ids/make_look_fasta.py`, `ids/plook_v5` and `ids/plook_v4` from `ids/plook.faa` (an
accession, `bbs|123`, `nat|NA000001.1|`, `lcl|777`, `gnl|PIG|124`, `bbs|300` then `lcl|300` and `lcl|301` then
`bbs|301` (one key each), `gnl|BL_ORD_ID|n` beside an accession with n its own number on two definition lines and on
one, the number of the PIG sequence and of no sequence, `bbs|500|gb|...|`, `lcl|778` beside an accession; each hit by
its query `ids/plook_q<k>.faa`, `.fna` for BLASTX), `ids/plook_z` from `ids/plook_z.faa`, read after `plook_v5` through
the alias `ids/plook_vz.pal` (its second sequence holds `gnl|BL_ORD_ID|1`, its own number in its volume), and
`ids/nlook_v5` and `ids/nlook_v4` from `ids/nlook.fna` (queries `ids/nlook_q<k>.fna` and `.faa`, `ids/nlook_qall.fna`
and `ids/nlook_qfound.fna` for the report streams). Version 5 finds `bbs|`, `nat|` and `lcl|<integer>` by their LMDB
keys unless an earlier sequence holds the key; version 4 finds neither `bbs|` nor `nat|`; no version finds `gnl|PIG|`.
A Seq-align id that another sequence also holds (review 1 of N2b): `ids/pshare_v5` and `ids/pshare_v4` from
`ids/pshare.faa` (`lcl|555`, then `gb|FFF00030.1|` with `lcl|555` in a second definition line and nearly the same
letters, so that the query `ids/pshare_q1.faa` hits both and NCBI's scope shows the second for the first; a
`gnl|BL_ORD_ID|0`, which makeblastdb stores with a text tag, query `ids/pshare_q2.faa`), `ids/pshare_a` and
`ids/pshare_b` (`bbs|123` in each) through the alias `ids/pshare_ab.pal` (query `ids/pshare_b_q0.faa`), and
`ids/nshare_v5` from `ids/nshare.fna` (query `ids/nshare_q1.fna`). Ids that differ in case alone (review 2 of N2b):
`ids/pshare.faa` 4 `lcl|abc` and 5 with `lcl|ABC` (query `ids/pshare_q5.faa`), and `ids/pshare_c` (`pdb|1ABC|A`,
`gnl|db|abc`) and `ids/pshare_d` (`pdb|1abc|A`, `gnl|DB|ABC`) through the alias `ids/pshare_cd.pal` (queries
`ids/pshare_d_q<k>.faa`); `ids/pdup_abc.pal` reads the copies `pdup_a` and `pdup_b` and the record of another title
`pdup_c`.

Seq-ids that NCBI's scope matches weakly (line DBID N2c; rows `db.id.*.weak.*`): `make_ids.py` (`--weak-only` builds
just these) writes, from the FASTA of `ids/make_weak_fasta.py`, `ids/pweak_v5` and `ids/pweak_v4` from `ids/pweak.faa`
(`sp|P00041.1|NAMEW`, `gb|FFF00040|`, `pir|PIRW1|`, `sp||NAMEX`, `gb|FFF00041.1|`, then three pairs S, T of nearly
equal letters: `sp||NAMEY` before `sp|P00042.1|NAMEY`, `gb|FFF00043|` before `gb|FFF00043.1|`, `sp|P00044.1|NAMEZ`
before `sp||NAMEZ`) and `ids/nweak_v5` and `ids/nweak_v4` from `ids/nweak.fna` (`gb|AF000081|` before
`gb|AF000081.1|`, then `gb|AF000082|`). NCBI's scope finds an accession without a version by the same one with a
version, an id of accession and name by the name-only id, a name-only id by an id of that name and an accession, and a
pir accession-only id by the name-only id; not the other way. The queries `ids/pweak_qa<k>.faa` (`.fna` for BLASTX,
`ids/nweak_qa2.fna`) carry such an id and are read with -parse_deflines (`qa4` the direction NCBI does not match);
`ids/pweak_q<k>.faa`, `ids/pweak_q8.fna`, `ids/nweak_q1.fna` and `.faa` hit record k (the second of a pair also hits
the first) and are run with `PRE_FETCH_SEQS_LIMIT=N` in `extra_args`, an environment variable of both programs (NCBI
reads the reported sequences one by one when N is 0 or below the batch's Seq-align count and then writes the names of
the sequence it read first). LOSAT refuses both cases and every N but 2147483647 (it does not count the Seq-aligns);
in a version 5 database LOSAT does not reproduce NCBI's lookup and refuses every Seq-align id with a name and every pir or prf id
without it, in version 4 the lookup of the id decides (a name alone with no other holder is written as NCBI does:
`db.id.blastp.weak.v4.pf.alone.6`, `v4.pf.none.6`). The `ncbi` rows are the searches without the variable or at 2147483647, the
version 4 sequences with no weak holder and outfmt 8. Review 1 of N2c: `ids/pweak_r_v5` and `ids/pweak_r_v4` from `ids/pweak_r.faa`,
four pairs S, T (T's first definition line holds its own id and the second id that S's matches weakly: `gb|FFF00051.1|`
with `sp||NAME_R` for `sp|P00051.1|name_r` (another letter case), `gb|FFF00053.1|` with `prf||NAME_S` for
`tpg|FFF00052.1|NAME_S`, `gb|FFF00054.1|` with `prf||PIRW3` for `pir|PIRW3|`, `gb|FFF00055.1|` with `emb|S12347.1|` for
`pir|S12347|`; queries `ids/pweak_r_q<k>.faa` are T's letters; rows `db.id.blastp.weak.r.*` and `.v4.r.*`: NCBI writes T's names
for S at `PRE_FETCH_SEQS_LIMIT=0` in both versions, and both its own rows without the variable).

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

Rows `db.tm.psi.*` and `db.tm.bp.off0.*` (line DBID, item P1): psiblast's prefetch, and `PRE_FETCH_SEQS_LIMIT`, which turns
NCBI's prefetch of a query batch off when it is 0 or below the batch's number of Seq-aligns (one per HSP; the count of
outfmt 7's `# N hits found`): its formatters then fetch the subjects each query's report shows, after the query's own
warnings (where a tabular report has a field that reads them: not `6 qseqid evalue`).

- `tm_psi.faa` and `tm_psi.map`: eight proteins `gb|PTX00001.1|` to `gb|PTX00008.1|` (written by the line's probe script,
  kept outside the repository), six of them copies of `q1` with more and more letters changed, with the taxids
  9606, 10090, 99990, 562, 7777 and none, and two of `q2` (99990, 9606). With `-max_target_seqs 1` psiblast's traceback
  hit list of `q1` still holds `PTX00003.1` (99990), whose lookup writes the offset error of `offset_last_past_end`.
- `tm_psi_q.faa` (`q1`, `q2`, `q3` with weak hits only), `tm_psi_q1.faa` (`q1`), `tm_psi_qw.faa` (`q1` with an `O`
  first, a warning of its own, then `q2` and `q1`), `tm_psi_qb.faa` (`q3`, then `q1`; `psiblast_taxdb_blocks.rs`).
- The databases (version 5, every file makeblastdb writes; `tm_psi_u.faa` is `tm_psi.faa` with the titles
  `subject PTX0000N.1`):

      makeblastdb -in tm_psi.faa -dbtype prot -parse_seqids -taxid_map tm_psi.map -blastdb_version 5 \
          -title "tm psi" -out tm_psi
      makeblastdb -in tm_psi_u.faa -dbtype prot -taxid 9606 -blastdb_version 5 -title "tm psi u" -out tm_psi_u

  `tm_psi_u`, without parsed Seq-ids, is searched in blocks of query batches (`LOSAT_X_BPQBLOCK`) by
  `tests/psiblast_taxdb_blocks.rs`.
- Review 1 of P1: a subject beyond `-max_target_seqs 1` whose Seq-align id NCBI's prefetch finds nowhere, whose lookup
  fails, or that an earlier sequence holds. `tm_psi_rq.faa` (`q1`) is the query; `tm_psi_ra.faa` (`gb|PSX00001.1|`, the
  top hit), `tm_psi_rpig.faa` (`gnl|PIG|124`), `tm_psi_rz.faa` (the top hit and `gnl|BL_ORD_ID|0` with a second
  definition line, which makeblastdb stores with a text tag), `tm_psi_rhx.faa` (`gb|PSX00011.1|`, an unrelated holder of
  `lcl|555` in a second definition line), `tm_psi_rhy.faa` (`lcl|555`, the second hit), with the aliases
  `tm_psi_pig.pal` (`DBLIST tm_psi_ra0 tm_psi_rpig`), `tm_psi_held.pal` (`DBLIST tm_psi_rhx0 tm_psi_ra0 tm_psi_rhy`: the
  holder without a taxid, the subject with 99990) and `tm_psi_held2.pal` (`DBLIST tm_psi_rhx tm_psi_ra0 tm_psi_rhy0`: the
  holder with 99990, the subject without):

      makeblastdb -in tm_psi_ra.faa -dbtype prot -parse_seqids -blastdb_version 5 -title tm_psi_ra0 -out tm_psi_ra0
      makeblastdb -in tm_psi_rpig.faa -dbtype prot -parse_seqids -taxid 99990 -blastdb_version 5 -title tm_psi_rpig \
          -out tm_psi_rpig
      makeblastdb -in tm_psi_rz.faa -dbtype prot -parse_seqids -blastdb_version 5 -title tm_psi_rz -out tm_psi_rz
      makeblastdb -in tm_psi_rhx.faa -dbtype prot -parse_seqids -blastdb_version 5 -title tm_psi_rhx0 -out tm_psi_rhx0
      makeblastdb -in tm_psi_rhy.faa -dbtype prot -parse_seqids -taxid 99990 -blastdb_version 5 -title tm_psi_rhy \
          -out tm_psi_rhy
      makeblastdb -in tm_psi_rhx.faa -dbtype prot -parse_seqids -taxid 99990 -blastdb_version 5 -title tm_psi_rhx \
          -out tm_psi_rhx
      makeblastdb -in tm_psi_rhy.faa -dbtype prot -parse_seqids -blastdb_version 5 -title tm_psi_rhy0 -out tm_psi_rhy0
- Review 2 of P1: `tm_psi_c150`, one volume without parsed Seq-ids holding the taxids 99990, 150, 99990 (`st1` to `st3` of
  `tm_cache.faa`) and 9606 (`tm_psi_ra.faa`), made from four volumes of one taxid each (`makeblastdb -taxid`, version 5,
  `cn1` to `cn4`, not kept) listed in an alias `cnx` with `makeblastdb -in cnx -input_type blastdb -dbtype prot
  -blastdb_version 5 -out tm_psi_c150 -title tm_psi_c150`; `tm_psi_qbig.faa`: 365 copies of the letters of `tm_cache_q.faa`
  (`qp0` to `qp364`, more than one block of 50,000 letters), then `q1` of `tm_psi_rq.faa`.

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
`makeblastdb -in n_l.fa -dbtype nucl -out n5l -parse_seqids -blastdb_version 5 -taxid 562 -title n5l` (a local id and a
text id on one OID). `q.faa`, `q.fna`, `qx.fna` are the queries, `qnn.fna`, `qxx.fna` the queries with two copies apart (two
ungapped HSPs per subject), and `sl_m.txt` a seqid list.

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

## Several databases in -db (DI-8f)

Rows `db.mdb.*`: `-db "a b"` searches the databases `a` and `b` as one virtual database, as NCBI's `CSeqDB` does (the `db` cell
of a row is one argument, so a cell with spaces is a list of names). `multidb/` holds the databases and what
`make_multidb.py MAKEBLASTDB` builds them from (makeblastdb 2.17.0; the records and taxids are those of `taxfilter/tfp.faa`,
`tfp.map`, `tfn.fna` and `tfn.map`, definition lines kept; `.pjs` and `.njs` dropped; the whole directory is about 0.5 MB):

- `pa` / `pb`: version 5 protein, `-parse_seqids`, taxid map; records 1-12 of `tfp.faa`, and records 13-23 plus `Q00001.1`, an
  exact copy of record 2's sequence (taxid 9606), so that two databases hold an identical sequence. `p4a` / `p4b`: the same
  records as version 4 without `-parse_seqids`.
- `na` / `nb`: version 5 nucleotide, `-parse_seqids`, taxid map; records 1-8 of `tfn.fna`, and records 9-16 plus
  `ENSTTEST0001.1`, a copy of record 1's sequence (taxid 9606). `n4a` / `n4b`: version 4 without `-parse_seqids`; `nnp`:
  the records of `nb` as version 5 without `-parse_seqids`.
- `d1/g` and `d2/g`: the records of `na` and `nb` under one base name in two directories. `et`: the records of `na` with an
  empty title. `al.nal`: an alias file (`TITLE multidb alias of na`, `DBLIST na`). `sl.txt`: a seqid list of the first
  accession of `na`, `nb`, `pa` and `pb`.

The ncbi rows search `taxfilter/q.fna` (BLASTN, BLASTX, TBLASTX) or `taxfilter/q.faa` (BLASTP, TBLASTN, PSIBLAST) and cover:
two names in outfmt 0, 5, 6, 7, 8, 11, 15 and 16 (`bn`), the other order (`bn.rev`), a name twice (`dup`), Galaxy's call with
each name in double quotes (`galaxy`), extra spaces around the names (`space`), an alias file (`alias`), one base name in two
directories (`samename`), two version 4 databases (`v4`, `bp.v4`), a version 4 and a version 5 database (`v4v5`: NCBI exits 2,
`DB list contains both Version 4 and Version 5 dbs`), an empty title (`title`), `-taxids` and `-seqidlist` over both names
(`tx`, `sl`), the five other programs (`bp`, `bx`, `tbn`, `tbx`; `ps` searches the version 4 pair, because LOSAT's PSIBLAST
refuses databases with parsed Seq-ids). The three `losat-rejection` rows are LOSAT's refusals of a version 5 pair with and
without the accession table (`bn.mixed.6`: NCBI stops at the first hit of a parsed sequence), of an alias that lists the
volume of another name (`bn.aliasdup.0`), and of a negative seqid list (`multidb/sl_negpartial.txt`, FBtr0083245 of the
second volume) on `taxfilter/tfn.00`, the first of the two volumes of one LMDB file (`bn.negpartial.6`: NCBI clears bit -1 of
its OID mask and crashes).

## The environment around BLASTX's search and report (line DBID, N2d)

Rows `db.bxenv.*` search `bslice` with `bslice_q.fna` (BLASTX) with an environment variable in `extra_args`. NCBI converts
`PRE_FETCH_SEQS_LIMIT` with `NStr::StringToInt` before the report of each query batch (`blast_app_util.cpp:732-737`) and
stops with a `CStringException` (exit 255; with `-mt_mode 1` its nodes catch it, so NCBI writes nothing and exits 0) on a
value it cannot convert, an empty one too; the archive of outfmt 11 is written without that call (`blastx_app.cpp:283-288`).
`BL2SEQ_LEGACY` makes the subjects of `-subject` a search of each subject on its own (`blast_app_util.cpp:206-210`); a `-db`
search ignores it. LOSAT refuses before the search: an error that NCBI raises in the first batch's search before the
prefetch (a query record without letters, `-query_loc` past every record: `BLAST engine error`, exit 3) and a `-subject`
search with `BL2SEQ_LEGACY` whose report would not change (no hits, or one subject record in outfmt 5, 6 or 10) are
refused too. The three `losat-rejection` rows are `abc` in outfmt 6
and 16 and an empty value in outfmt 0; the ncbi rows are `abc` in outfmt 11, the integer 5 in outfmt 6, and `BL2SEQ_LEGACY=1`
in outfmt 0 and 6. The `-subject` cases are `environment_that_ncbi_reads_around_the_report_is_rejected_where_it_reads_it`
in `tests/run_local_blastx.rs`.

## The same environment in the other programs (line DBID, N2e)

Rows `db.envshared.*` search `arch_sprot` (BLASTP, PSIBLAST), `arch_tbn_snucl` (TBLASTN), `arch_tbx_snucl` (TBLASTX) and
`arch_snucl` (BLASTN) with the queries of `tests/fasta/outfmt11/` and the environment of the section above. BLASTP, TBLASTN,
TBLASTX and BLASTN write the archive of outfmt 11 without the prefetch too (`blastp_app.cpp:280-285`,
`tblastn_app.cpp:294-299`, `tblastx_app.cpp:197-202`, `blastn_app.cpp:304-309`), so a value that NCBI cannot convert leaves it
as without the variable; psiblast's iterations fetch ahead before they write an archive (`psiblast_app.cpp:262-266`), which
LOSAT's PSIBLAST does not write. A `-db` search of each of them ignores `BL2SEQ_LEGACY`. NCBI's automatic split of a `-db`
search by queries (`-num_threads` 2 or more, a large query file, a small database) catches the exception as `-mt_mode 1` does. The ncbi rows are `abc` (BLASTN: an
empty value) in outfmt 11 and `BL2SEQ_LEGACY=1` in outfmt 6 (TBLASTN: 0); the `losat-rejection` row is BLASTN's `abc` in
outfmt 6. The `-subject` cases of the five programs are `tests/search_environment.rs`.

## TBLASTN `-query_loc` and the chunk masks (line NP, series SM)

Rows `db.tqloc.*` search `tqloc` with `tqloc_q.faa`, `-query_loc 5-3004 -lcase_masking` and
`CHUNK_SIZE=1000 OVERLAP_CHUNK_SIZE=100`. NCBI restricts the lower-case runs of the whole record to each chunk
(`split_query_cxx.cpp:278-281`), one letter longer, so a run that ends on the letter before the interval masks the
interval's first letter in the first chunk's preliminary search. The inputs come from `random.Random(20261104)`:

- `tqloc_q.faa`: `a` (`AC`, which the range skips), `p` (`MKVA` and 200 random residues) and `q` (`MKVa` and a random
  protein Q of 3,000 residues that starts with W). `q` is the batch's second query and its part starts the first chunk.
- `tqloc.fna`: `s14` and `s30`, the first 14 and 30 residues of Q between 50 random residues on each side,
  back-translated with one codon per residue. At `-evalue 1e-7` NCBI reports `q` against `s30` only; with `MKVA` in
  place of `MKVa` it reports `s14` too.

The database (version 4, no `-parse_seqids`; only `.nin .nsq .nhr` are kept). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/tqloc.fna -dbtype nucl -out tests/db_search/tqloc \
        -title "Q[0:14] and Q[0:30] of tqloc_q.faa between random letters, back-translated" -blastdb_version 4
    rm tests/db_search/tqloc.njs

The `-subject` rows `tqloc.tblastn.*` of `outfmt0_manifest.tsv` use `tests/fasta/outfmt0/tqloc_q.faa` (`MKVa`, the same
Q and 37,000 more random residues: 40,000 searched letters, which NCBI splits at the default sizes) and
`tests/fasta/outfmt0/tqloc_s.fna` (a copy of `tqloc.fna`), at `-evalue 1e-6`.

## BLASTX `-query_loc` and the chunk masks (line NP, series SM)

Rows `db.bxqloc.*` search `bxqloc` with `bxqloc_q.fna`, `-query_loc 10-9009 -lcase_masking` and
`CHUNK_SIZE=1002 OVERLAP_CHUNK_SIZE=99` (a BLASTX chunk is a multiple of 3). NCBI restricts the lower-case runs of the
whole record to each chunk (`split_query_cxx.cpp:278-281`), one letter longer, so a run that ends on the letter before
the interval masks the interval's first letter in the first chunk's preliminary search. NCBI splits a BLASTX query
batch only when it holds one query. The inputs come from `random.Random(20261104)` (Q and its back-translation, as in
the `db.tqloc.*` rows) and `random.Random(20261105)`:

- `bxqloc_q.fna`: `a` (`AC`, which the range skips), `p` (`ATGAAAGTT` and 600 random nucleotides) and `q`
  (`ATGAAAGTt` and Q of 3,000 residues, which starts with W, back-translated with one codon per residue: 9,000
  nucleotides). Without `BATCH_SIZE`, `p` and `q` share a query batch, which is not split (`db.bxqloc.two.6`: `s14` is
  reported); with `BATCH_SIZE=100` `q` is alone in its batch and split (`db.bxqloc.b100.6`: at `-evalue 1e-7` NCBI
  reports `q` against `s30` only; with `ATGAAAGTT` in place of `ATGAAAGTt` it reports `s14` too).
- `bxqloc.faa`: `s14` and `s30`, the first 14 and 30 residues of Q between 50 random residues on each side.

The database (version 4, no `-parse_seqids`; only `.pin .psq .phr` are kept). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/bxqloc.faa -dbtype prot -out tests/db_search/bxqloc \
        -title "Q[0:14] and Q[0:30] of bxqloc_q.fna's q between random letters" -blastdb_version 4
    rm tests/db_search/bxqloc.pjs

The `-subject` rows `bxqloc.blastx.*` of `outfmt0_manifest.tsv` use `tests/fasta/outfmt0/bxqloc_q.fna` (`ATGAAAGTt`,
the same 9,000 nucleotides and 11,000 more random ones: 20,000 searched nucleotides, which NCBI splits at the default
sizes, 10,002 and 297) and `tests/fasta/outfmt0/bxqloc_s.faa` (a copy of `bxqloc.faa`), at `-evalue 2e-7`.

## BLASTN `-query_loc` and the chunk masks (line NP, series SM)

Rows `db.bnqloc.*` search `bnqloc` with `bnqloc_q.fna`, `-query_loc 10-3009 -lcase_masking` and
`CHUNK_SIZE=1001 OVERLAP_CHUNK_SIZE=100` (NCBI BLASTN stops at 1000/100 with `Empty CBlastQueryVector`: its batch size
mixer's maximum is the chunk size less 1,000). NCBI restricts the lower-case runs of the whole record to each chunk
(`split_query_cxx.cpp:278-281`), one letter longer, so a run that holds the letter before the interval masks the
interval's first letter in the preliminary search of each chunk that holds it, when the query's part there is longer than
the interval's start (here 9 letters). Without `BATCH_SIZE` each record is a batch of its own; with `BATCH_SIZE=100000`
`p` and `q` share a split batch (four chunks of 1,126 letters) and `q`'s first letter lies 174 letters into the second
chunk. The inputs come from `random.Random(20261105)`:

- `bnqloc_q.fna`: `a` (`AC`, which the range skips), `p` (1,209 random nucleotides) and `q` (`ACGTACGTa` and B, 3,000
  random nucleotides).
- `bnqloc.fna`: `s11`, `s24` and `s40`, the first 11, 24 and 40 nucleotides of B between 50 random nucleotides on each
  side (the nucleotide after each differs from B's).

At `-task blastn -soft_masking false -evalue 1e-7` (`db.bnqloc.hard.*`, `db.bnqloc.two.6`) NCBI reports `q` against
`s40` only; with `ACGTACGTA` in place of `ACGTACGTa` it reports `s24` too. With the default `-soft_masking true` at
`-evalue 10` (`db.bnqloc.soft.6`) the masked letter leaves the lookup table, and `s11`, whose only seed starts there, is
not reported (with `ACGTACGTA` it is). Megablast at `-evalue 5e-17` (`db.bnqloc.mb.7`) reports no hit (with
`ACGTACGTA`, `s40`).

The database (version 4, no `-parse_seqids`; only `.nin .nsq .nhr` are kept). Run from `LOSAT/`:

    makeblastdb -in tests/db_search/bnqloc.fna -dbtype nucl -out tests/db_search/bnqloc \
        -title "B[0:11], B[0:24] and B[0:40] of bnqloc_q.fna's q between random letters" -blastdb_version 4
    rm tests/db_search/bnqloc.njs

## `nident` alone and HSPs without identities (line DBID, L1)

Rows `db.nz.*`: with `PRE_FETCH_SEQS_LIMIT` turning NCBI's prefetch off, a tabular report whose only field of the
alignment block is `nident` (`6 qseqid nident`; in the translated programs also with `sframe` or `frames`) reads a
subject only on a row whose HSP has identities (the Seq-align's `num_ident` above 0, tabular.cpp:904), so a subject whose
shown HSPs have none is not looked up in the taxonomy database. The `pident` rows are the controls (every row reads the
subject). `-threshold 9` lets the seeds of alignments without identities through (a word of three such pairs scores at
most 9 in BLOSUM62).

- `nz_p.faa` (protein) and `nz_n.fna` (its proteins as codons): `sZ` (`M`, then `VY` 30 times, `CC` and `VY` 10 times),
  `sI` (`M` and B, 60 random residues), `sM` (`M`, `VY` 30 times, `W` 20 times and B) and `sN` (`M` and `VY` 20 times),
  with the taxids 9606, 562, 10090 and 99990 (`nz_p.map`, `nz_n.map`): only `sN` has the taxid whose lookup writes the
  offset error of `offset_last_past_end`. Queries `nz_q.faa` / `nz_q.fna`: `q1` (`IF` 30 times: aligned to the `VY`
  parts without identities), `q2` (`q1`, `D` 20 times and B), `q3` (B), `q4` (B, `D` 20 times and `q1`); `nz_q1.*` and
  `nz_q2.*` hold `q1` and `q2` alone. In BLASTX, TBLASTN and TBLASTX other frames align with identities too, hence
  `-max_hsps 1` and lower e-values in some rows.
- `nz_psi.faa`: psiblast's composition-based statistics drop the `IF`/`VY` alignments, so `sP` is `M` and `qP` (120
  residues of `IFKELSD`) with I/V, F/Y, K/R, E/Q, L/M, S/T and D/N exchanged, `sQ` is `M` and `qP`, `sR` the first half
  of `sP` and the second of `sQ`, with the taxids 99990, 9606 and 10090 (`nz_psi.map`); `nz_psi_q.faa` holds `qP`, its
  first half `qPh` and its second half `qP2`.
- `nz_pc.faa` / `nz_nc.fna` (`db.nz.*.c150.*`) are `nz_p` / `nz_n` with `sZ`'s taxid 150 (`nz_pc.map`, `nz_nc.map`),
  which NCBI's taxonomy cache (256 slots, taxid & 255) puts in the slot of 99990, the damaged last entry of
  `offset_last_past_end`: LOSAT refuses a run that may look both up (`TaxdbRun::cache_collision`), but not one whose
  rows look nothing up (`q1` alone with `nident`, a positive `PRE_FETCH_SEQS_LIMIT` below the batch's Seq-aligns).
  BLASTX writes each batch before it searches the next, so it refuses such a run at its first batch also where no row
  looks a tax id up (`db.nz.bx.c150.*`, as before line L1).
- Rows `db.nz.cq.*` (line DBID, L1b): TBLASTN and TBLASTX without the report stream make every batch's lookups before
  any report, so they ask for that refusal only where NCBI fetches a subject: not for fields that read none (`6 qseqid
  evalue`) where `PRE_FETCH_SEQS_LIMIT` turns the prefetch off, nor for a batch without hits, whose prefetch loads no
  sequence. `nz_qx.faa` / `nz_qx.fna` hold `qx` (`M` and 70 random residues, no hit at `-evalue 1e-3`) and
  `nz_qx2.*` `qx` then `q2`: with `BATCH_SIZE=1` the refusal of `q2`'s batch still comes before any byte. TBLASTX's
  report stream (outfmt 7 in two batches), which would write `qx`'s block before `q2`'s lookups, is off for a run that
  the taxonomy cache refuses, so `nz_qxy.fna` (`qx`, then `qy` without hits either) gets NCBI's bytes.
- Rows `db.nz.cq.bn.*` (line DBID, L1c): BLASTN asks for that refusal only where NCBI fetches a subject too: not for
  fields that read none (`6 qseqid evalue`) where a positive `PRE_FETCH_SEQS_LIMIT` below the batch's Seq-aligns turns
  the prefetch off, not for queries without hits whose reports would fetch (`PRE_FETCH_SEQS_LIMIT=0`, outfmt 6), and not
  for a batch without hits (`nz_qx.fna` in outfmt 0; `nz_qxy.fna` with `BATCH_SIZE=1` in outfmt 7, where the report
  stream writes each batch's reports). `nz_qx2.fna` with `BATCH_SIZE=1` is refused at `q2`'s batch, whose subjects
  (taxids 562 and 10090) NCBI fetches without a message (the superset kept since L1b); outfmt 6 writes nothing for
  `qx` before it, and `tests/blastn_report_stream.rs` checks that outfmt 0 and 7 write `qx`'s report first where the
  report stream acts (and only outfmt 0's prolog where it does not). The rows use `-task blastn -dust no -evalue 1e-3`
  (megablast and dc-megablast find no hit of these queries).
- Rows `db.nz.cq.x.*` (line DBID, L1d): the XML2 and JSON formats (outfmt 13 to 16), whose reports fetch every reported
  subject, ask for that refusal too, where a reported subject has a positive taxid. `nz_rep.fna` holds `a1` (the
  letters of `sN` as in `nz_n`), `b1` (`sZ`) and `c1` (`sI`) with the taxids 99990, 150 and 99990 (`nz_rep.map`), `d1`
  (the codons of 60 random residues) without a taxid and `e1` (another 60) with the taxid -5: one volume where NCBI
  really writes the offset error again (for `c1` after `b1`; `nz_nc`'s only 99990 is `sN`, whose entry NCBI keeps
  loaded, so it writes the message once there). `nz_rep_q.fna` holds `qN`, `qZ` and `qI`, the whole of `a1`, `b1` and
  `c1`, and `nz_rep_q.faa` their translations (`nz_repp.faa` holds those of `a1` to `d1`): NCBI writes the message 3 to
  5 times in outfmt 15 / 16 and LOSAT refuses before any byte, as in outfmt 6 with `BATCH_SIZE=1` (twice).
  `nz_rep_qd.*` (`qD`, the start of `d1`) and `nz_rep_qe.fna` (`qE`, the start of `e1`) hit `d1` or `e1` alone: no
  volume looks a positive taxid up, so their `ncbi` rows are not refused. `nz_repo` and `nz_reppo` hold `a1` to `d1`
  (nucleotide, protein) without parsed Seq-ids: their taxid is that of the definition line, a branch of its own
  (`db.nz.cq.x.bn.repo.*`: NCBI writes the message 4 times with `BATCH_SIZE=1`). `nz_rep_qxn.fna` is `nz_qx.fna`'s `qx`
  then `qN`: BLASTX writes each batch's XML2 report before it searches the next, so on `nz_reppo` it holds the
  document of a run with this refusal and drops it at `qN`'s batch (no byte of `qx`'s); `tm_cache_qx.fna` (`qx` then
  the three records of `x16_bx_taxid_q.fna`) does the same on `tm_cache150`, whose parsed Seq-ids make BLASTX hold
  anyway. BLASTP's search in blocks of query batches takes the search of all queries for such a run
  (`tests/blastp_query_blocks.rs`, `nz_reppo`).
- The sequences come from the line's probe script (kept outside the repository). The databases (version 4, every file
  makeblastdb writes; `nz_pc` and `nz_nc` alike with their maps), run from `LOSAT/`:

      makeblastdb -in tests/db_search/nz_p.faa -dbtype prot -parse_seqids -taxid_map tests/db_search/nz_p.map \
          -blastdb_version 4 -title "nz p" -out tests/db_search/nz_p
      makeblastdb -in tests/db_search/nz_n.fna -dbtype nucl -parse_seqids -taxid_map tests/db_search/nz_n.map \
          -blastdb_version 4 -title "nz n" -out tests/db_search/nz_n
      makeblastdb -in tests/db_search/nz_psi.faa -dbtype prot -parse_seqids -taxid_map tests/db_search/nz_psi.map \
          -blastdb_version 4 -title "nz psi" -out tests/db_search/nz_psi
      makeblastdb -in tests/db_search/nz_rep.fna -dbtype nucl -parse_seqids -taxid_map tests/db_search/nz_rep.map \
          -blastdb_version 4 -title "nz rep" -out tests/db_search/nz_rep
      D=$PWD/tests/db_search; cd "$(mktemp -d)"   # nz_repo, nz_reppo: one database per record, then one volume
      for k in a1 b1 c1 d1; do  # a1, b1, c1, d1 of nz_rep.fna / nz_repp.faa, each in $k.fna / $k.faa
          t=""; case $k in a1|c1) t="-taxid 99990";; b1) t="-taxid 150";; esac
          makeblastdb -in $k.fna -dbtype nucl $t -blastdb_version 4 -out n$k
          makeblastdb -in $k.faa -dbtype prot $t -blastdb_version 4 -out p$k
      done
      makeblastdb -in "na1 nb1 nc1 nd1" -input_type blastdb -dbtype nucl -blastdb_version 4 -title "nz repo" \
          -out $D/nz_repo
      makeblastdb -in "pa1 pb1 pc1 pd1" -input_type blastdb -dbtype prot -blastdb_version 4 -title "nz reppo" \
          -out $D/nz_reppo
