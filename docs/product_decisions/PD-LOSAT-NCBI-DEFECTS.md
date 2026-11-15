# Product Decision: NCBI BLAST+ behaviour that is a defect

- Decision ID: `PD-LOSAT-NCBI-DEFECTS`
- Version: 1.32
- Date: 2026-10-02 (1.0); 1.1 the same day (the three confirmations below, Session S07+++b);
  1.2 2026-10-03 (exception 2 for TBLASTX and TBLASTN, Session S08b); 1.3 2026-10-05
  (BLASTP, TBLASTN and TBLASTX of stage E2e: exception 3 and the items below, Session S08+b);
  1.4 2026-10-22 (exception 4, the culling limit near `INT4_MAX`, NCBI parity line O3b);
  1.5 2026-10-30 (exception 5, the expansion of a taxid whose ancestors form a cycle, line DI-8a);
  1.6 2026-10-30 (exception 6, an accession stored under more than about 500 sequences of a
  version 5 database, and the seqid-list items of rules 2 and 3, line DI-8b);
  1.7 2026-10-31 (exception 6 also for the id lookups of `losat blastdbcmd`, line BDC's BDC2b);
  1.8 2026-11-01 (exception 5 also for `losat blastdbcmd -taxids`, line BDC's BDC3a);
  1.9 2026-11-01 (a version 5 database list that mixes databases with and without accession
  tables, and a negative seqid list that excludes a sequence of a volume the database does
  not hold, rule 3, line DI-8f);
  1.10 2026-11-02 (the type of a `gb`, `emb` or `dbj` id, exception 2 for the tabular first
  word of an `lcl|Subject_N` id, and the reported sequences of a database with parsed Seq-ids
  that are refused, line DBID);
  1.11 2026-11-04 (the makeblastdb exception: a `-in` file of several bzip2 streams or zstd
  frames, line GZ; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07);
  1.12 2026-11-12 (exception 8, a TBLASTX cutoff quotient past the range of a double, line PS);
  1.13 2026-11-12 (exception 9, the `losat rpsbproc` exception: PSSM-IDs of 2^31 and more,
  line PS item PS16a; decided 2026-11-11, numbered after 1.12);
  1.14 2026-11-12 (exception 10, the `losat rpsbproc` exception: a motif check that reads past
  the stored positions, line PS item PS16b; decided 2026-11-11, numbered after 1.13);
  1.15 2026-11-12 (the explicit rejections of `losat rpsbproc`, rule 3, line PS items PS16a and
  PS16b);
  1.16 2026-11-13 (the explicit rejections of `losat rpsbproc` for nucleotide queries, rule 3,
  line PS item PS16c);
  1.17 2026-11-13 (the explicit rejections of `losat rpsbproc` in the XML input, rule 3, line PS
  item PS16d);
  1.18 2026-11-14 (the explicit rejections of TBLASTN's one-hit gapped start past NCBI's buffers
  and of the endless split it can cause, rule 3, line RG item X2);
  1.19 2026-11-14 (exception 11, a REAL whose exponent overflows NCBI's `int`; the XML input of
  `losat rpsbproc` reads a REAL of more than 18 significant digits as NCBI does, which 1.17
  refused; line SW item SW8 B28);
  1.20 2026-11-14 (the explicit rejections of `-outfmt 18`, rule 3, line RG items T2 and T2c);
  1.21 2026-11-14 (TBLASTX's one-hit HSPs of length 0 or less ported, rule 4, and the explicit
  rejection of their merge across subject chunks with a sum of lengths of 0, rule 3, line RG item
  X2b);
  1.22 2026-11-14 (the explicit rejections of BLASTX's one-hit gapped start past NCBI's buffers,
  rule 3, line RG item X2f);
  1.23 2026-11-14 (the explicit rejection of BLASTP's query-anchored section that NCBI's merge
  aborts on, rule 3, line RG item T3a; Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07);
  1.24 2026-11-14 (the rule-2 rejection of a query's id widened to the display id and to generated
  ids, line RG item T3a; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07,
  within the class the maintainer accepted with 1.10 on 2026-11-02);
  1.25 2026-11-14 (a version 5 database list that mixes databases with and without accession
  tables is searched; the rejection of 1.9 keeps only a run whose reports look a parsed Seq-id up
  in a file without the table, rule 3, line DBID item mixdb, Owner-delegated);
  1.26 2026-11-14 (exception 12, the `losat rpsbproc` exception: coordinate arithmetic past NCBI's
  `int`, protein and nucleotide queries and the XML input alike; line PS item PS16-Q1, the Owner's
  Q1 = B with the Owner-delegated additions of line PS `DECISIONS.md` 58 and 59; `MapSegSet`'s segment
  end reproduced, rule 2);
  1.27 2026-11-15 (the rule-2 item of a seqid list of gi ids only narrowed to what NCBI's runs
  show: a positive list keeps every sequence only where NCBI's OID bit set is its special
  all-set one, and none elsewhere, line DBID item gisl, Owner-delegated);
  1.28 2026-11-15 (the explicit rejection of 1.23 and the rule-2 rejection of a query's id of 1.24
  in the query-anchored views of BLASTN and TBLASTN, line RG item T3b; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07, within the classes of 1.23 and 1.10);
  1.29 2026-11-15 (the rejection of a sequence that a seqid list keeps while matching none of its
  lines, of 1.6, also where a list's `gnl|BL_ORD_ID|n` reaches a line in a volume whose first OID is
  not 0, which NCBI's header read numbers differently from the lookup; line DBID item ordsl,
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07, within the class of 1.6);
  1.30 2026-11-15 (the explicit rejection of RPSBLAST's Smith-Waterman redo that NCBI repeats
  without end, rule 3, line PS item PS11c-3; Owner-delegated, line PS `DECISIONS.md` 63);
  1.31 2026-11-15 (the stop of a batch's prefetch at a Seq-id that a version 5 list without the
  accession table in one of its databases cannot look up is written as NCBI writes it; the
  refusal of 1.25 narrowed to the prefetch off and to the orders of loads that LOSAT does not
  model; line DBID item mixdb2, Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07);
  1.32 2026-11-15 (the explicit rejections of psiblast and tblastn `-in_pssm` where NCBI's object
  manager or format guess prints addresses, reads outside its data, overflows a signed integer or
  depends on the machine's memory, rule 3, line RG item D3-1; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; no rule-1 exception: the Owner's questions of D3-1
  stay open)
- Status: Accepted by the maintainer on 2026-10-02, in Session S07+++b (E2g), on the
  NCBI BLAST+ 2.17.0 behaviours that the BLASTN inventory
  (`docs/evidence/losat_web_e2g/INVENTORY.tsv`; moved to a release asset on 2026-10-23:
  [ARCHIVED.md](../evidence/losat_web_e2g/ARCHIVED.md)) and the independent audits found to be
  defects. Plan decision DW-15 in [`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md).

## Scope

NCBI BLAST+ 2.17.0 behaviour that is a defect rather than a design: a crash, an
exception that only a debug-build assertion was meant to prevent, a read past the end of
a buffer, or an integer that wraps. Everything else stays under the root
[`AGENTS.md`](../../AGENTS.md) bit-perfect rule. Versions 1.0 and 1.1 list the BLASTN
items; version 1.2 adds TBLASTX and TBLASTN; version 1.3 the BLASTP, TBLASTN and TBLASTX
items of their non-default options (E2e); BLASTX adds its own in its session under the same
rule.

## Rule

1. **NCBI fails, and LOSAT can give a valid result:** an approved exception. The result
   is valid when it equals NCBI's output for a nearby input or configuration that does
   not reach the defect (evidence under `docs/evidence/`).
2. **NCBI gives a deterministic result, even one that looks wrong:** LOSAT reproduces it
   byte for byte. This is not an exception. The maintainer may keep an explicit rejection
   instead for a setting without a practical use whose reproduction is costly (version
   1.1: the two items so marked below).
3. **NCBI fails, and no valid result can be defined or checked:** an explicit rejection
   ("... is not supported by LOSAT's BLASTN"). This is not an exception.
4. **NCBI accepts an input deterministically that LOSAT rejected:** LOSAT ports it.

## Approved exceptions (BLASTN)

1. **A query chunk that NCBI would split again.** With `CHUNK_SIZE` and
   `OVERLAP_CHUNK_SIZE` such that a query chunk is long enough to be split again (an
   overlap close to the chunk size), NCBI skips the chunk's lookup table and stops with a
   `CCoreException` (null pointer) that names its build's source files, exit 3
   (`split_query_aux_priv.cpp:190-201`, `blast_aux_priv.cpp:206-208`; only a debug build
   asserts against it). LOSAT searches each chunk once. Evidence
   `docs/evidence/losat_web_e2g/resplit/`: 42 configurations (3 queries of 30 to 60 kb,
   blastn and megablast, 7 chunk and overlap pairs); NCBI exit 3 in all; LOSAT's output
   equals NCBI's unsplit output in 35, and NCBI's output at the largest overlap that does
   not split a chunk again in 74 of 84 runs (outfmt 6 and 0). The 10 others differ only
   by the HSPs at chunk boundaries that NCBI's own splits also lose. The round-3 audit
   (`~/.cache/losat-web-gui-target/e2g-audit/r3b/`): 61 of 65 sampled rows equal NCBI at
   that largest overlap byte for byte; in batches of several queries with weak HSPs
   (E-values near 1 to 10) LOSAT's output differs from NCBI's at the non-re-splitting
   overlaps by about as many weak HSPs as NCBI's outputs at two such overlaps differ from
   each other. Fixtures `env.resplit_*` (expected output from NCBI at that largest overlap,
   `oracle_env`).
   With a negative `CHUNK_SIZE` above a negative `OVERLAP_CHUNK_SIZE` the same failure is
   an explicit rejection (below).
2. **outfmt 0 titles made only of punctuation.** For a subject defline of commas,
   semicolons, tildes and spaces that ends in a run of spaces and separators (such as
   `, ,`), NCBI's `x_CleanAndCompress` (`create_defline.cpp:219-312`) lets its `size_t`
   count of the remaining letters wrap, reads past the end of the string and crashes when
   it writes the title of such a subject with hits (SIGSEGV). LOSAT stops the cleanup at
   the end of the string (`, ,` gives `, `). Evidence
   `docs/evidence/losat_web_e2g/punct_defline/`: 8 such deflines on subjects with hits,
   megablast and blastn, with and without `-max_target_seqs 3`; NCBI crashes in all 4
   runs; LOSAT's report equals NCBI's report for the same subjects with placeholder
   deflines once each placeholder is replaced by LOSAT's title. Subjects with such
   deflines and no hits, where NCBI runs, match NCBI (fixtures `punct.nohit_*`).

## Approved exceptions (TBLASTX and TBLASTN)

Version 1.2, accepted by the maintainer on 2026-10-03 in Session S08b (E2b), plan DW-17.

2. **outfmt 0 titles made only of punctuation** (BLASTN exception 2, extended). NCBI
   tblastx and tblastn build the outfmt 0 title of a nucleotide subject with the same
   `CDeflineGenerator` and `x_CleanAndCompress` (`create_defline.cpp:219-312`) and crash
   (SIGSEGV) on the same deflines when such a subject has hits; outfmt 6 and 7 do not build
   the title, and there NCBI and LOSAT agree. Evidence `docs/evidence/losat_web_e2b/`:
   `punct_defline.py` (deflines `, ,` and `;~ ;`: NCBI dies of the signal in outfmt 0,
   outfmt 6 and 7 equal), and the S08 audit (c), rounds 1 and 2: over all 15624 TBLASTX
   deflines of up to 6 and all 3124 TBLASTN deflines of up to 5 of `,;~ a` and space, NCBI
   crashes on 272 and 66, exactly those that LOSAT's test of such a title
   (`report/defline.rs` `ncbi_nucleotide_title_reads_past_end`) finds; the other reports are
   byte-identical. The approved result: LOSAT stops the cleanup at the end of the string, as
   for BLASTN. Implemented in Session S08+ (E2e; until then LOSAT rejected such subjects
   with hits in outfmt 0 explicitly) and checked as BLASTN's was: over all 1023 deflines of
   one to five of `,;~` and space, each TBLASTX and TBLASTN report equals NCBI's (957) or
   equals NCBI's report for the same subject with a placeholder defline once the
   placeholder is replaced by LOSAT's title (66, where NCBI crashes),
   `docs/evidence/losat_web_e2e/title_sweep.py`; fixtures `punct.tblastx` and
   `punct.tblastn` of `LOSAT/tests/outfmt0_manifest.tsv` (contract `approved_punct_title`).
   This records the implementation; the decision of version 1.2 is unchanged.

## BLASTP, TBLASTN and TBLASTX non-default options (version 1.3)

Version 1.3, accepted by the maintainer on 2026-10-05 in Session S08+b (E2e), plan DW-19:
the decisions D11 to D15 of `docs/evidence/losat_web_e2e/AUTHORITY.md` §M as recommended.

3. **Approved exception (BLASTP): a one-hit gapped start that reads past a sequence.**
   With `-window_size 0` (one-hit) and a low `-threshold`, NCBI's ungapped extension can
   give a negative `Int4` length that `BlastGetStartForGappedAlignment` receives as
   `Uint4` (`aa_ungapped.c:1054,1083`, `blast_gapalign.c:3393-3437`): it scores only the
   first window of 11 letters, and when that window passes the end of the query or the
   subject it reads past the sequence buffer and indexes the matrix with the stray byte;
   some such searches crash (SIGSEGV, exit 139; 21 of 5,000 random short searches with
   `-threshold` 1 to 3 in the round-2 audit). LOSAT reads the letters past the sequence
   as the sentinel (`NULLB`, `BLAST_SCORE_MIN`), as for the letter right after it. NCBI run
   under `valgrind -q`, which reports the invalid reads and does not crash, gives LOSAT's
   output byte for byte in all 23 crashing cases of the audit's second and third passes
   (`docs/evidence/losat_web_e2e/audit/round2/a2_blastp.md`, `a3_blastp.md`). Searches
   whose window stays inside the sequences and its sentinel match NCBI (fixture
   `e2e.blastp.one_hit_negative_width`).

Explicit rejections (rule 3):

- `-evalue` infinite or of DBL_MAX or more (BLASTP, TBLASTN; decision D12): a hit list
  whose HSPs all drop keeps the best e-value DBL_MAX, which passes `best_evalue <=
  expect_value` (`blast_kappa.c:409,3687`), and the empty list is read
  (`blast_hits.c:3266`): NCBI crashes on some inputs, and LOSAT cannot tell in advance
  which. `1.7976931348623156e308` and below run and match NCBI. TBLASTX runs these values
  and matches NCBI.
- Query-split settings that NCBI fails on (BLASTP, TBLASTN; decision D13): a
  `CHUNK_SIZE`/`OVERLAP_CHUNK_SIZE` pair with which a chunk would be split again, and a
  negative `CHUNK_SIZE` that splits a batch. Unlike BLASTN's exception 1, the maintainer
  kept the rejection for these programs.
- Settings on which NCBI crashes (decision D3): an IDENTITY word size of 5 or 6,
  `-use_sw_tback -ungapped -comp_based_stats 0`, compressed-lookup thresholds that
  exhaust its overflow bank, and `-max_target_seqs` values that make the preliminary hit
  list size not positive. A query length plus `-window_size` in (2^30, 2^31 - 1], where
  NCBI does not finish (decision D8).

Kept as explicit rejections although NCBI's result is deterministic (rule 2's
maintainer clause; decision D11): a query length plus `-window_size` above 2^31 - 1, where
NCBI's `Int4` diagonal table wraps to one cell and no hit is found (BLASTP, TBLASTX;
together with D8, every sum above 2^30 is rejected).

Reproduced as NCBI (rule 2): the out-of-range `double` to `Int4` conversions of
`-threshold` and the X-drops (INT_MIN, decision D4); `-max_target_seqs` of 2147483624 or
more in BLASTP, whose preliminary hit list size wraps to 2 to 48. TBLASTX's culling limit
plus 3 wrapping in `Int4` (decision D1) was reproduced until version 1.4 (exception 4).

## Approved exception (every program with `-culling_limit`, version 1.4)

Version 1.4, decided by the maintainer on 2026-10-22 ("広い型になおしましょう", NCBI parity
line `DECISIONS.md` 17), under the rule of 2026-10-10 that LOSAT does not reproduce a signed
`Int4` overflow (undefined behaviour in C), as it already adds TBLASTX's chain sums in `i64`.

4. **A culling limit near `INT4_MAX`.** The preliminary culling writer adds 3 to a limit
   above 1 (`setup_factory.cpp:338-340`) and every HSP's merit counts down from it
   (`hspfilter_culling.c:123-189`), both in `Int4`. From 2147483645 on the sum overflows:
   NCBI 2.17.0's binary wraps it, so 2147483645 starts every merit at `INT_MIN` and culls
   nothing, and 2147483646 and 2147483647 start negative and act as a limit of 1 in the
   preliminary stage (the pipe after the traceback uses the limit itself and culls
   nothing). LOSAT computes the sum and the merits in `i64` (BLASTN, BLASTX, TBLASTX and the
   programs that add the option later), so these limits cull nothing, as 2147483644 does in
   NCBI. Evidence: NCBI's output at 2147483644 equals its output at 2147483645 and differs at
   2147483646 and 2147483647 (BLASTN and TBLASTX fixtures `culling.near_int4_max_*` of
   `LOSAT/tests/{blastn,tblastx}_regression_fixtures.py`, whose oracle runs 2147483644 for
   the last two, `ORACLE_ARGV`). Every other limit gives NCBI's bytes.

## Approved exception (`-taxids` and the other taxid lists with `-db`, version 1.5)

Version 1.5, decided under the maintainer's standing delegation on 2026-10-29 after the
maintainer asked whether such a list could be filtered instead of refused (line DI-8
`DECISIONS.md` 10, replacing its decision 9 to refuse it); its PR (DI-8a) asks the
maintainer to accept it.

5. **A taxid whose ancestors form a cycle in `taxonomy4blast.sqlite3`.** The root, taxid 1,
   is its own parent in NCBI's file, and NCBI expands a taxid of a list to its descendants
   with `WITH RECURSIVE ... UNION ALL` (`tax4blastsqlite.cpp`, `GetLeafNodeTaxids`), which
   does not end for a start whose descendants reach the cycle: NCBI 2.17.0's `blastp -db
   swissprot -taxids 1` writes nothing in 120 s. LOSAT visits each taxid once
   (`TaxonomySqlite::descendants`), so the root expands to every taxid of the file and the
   search ends; the other taxids expand as in NCBI. Evidence: with
   `LOSAT/tests/db_search/taxfilter/root/taxonomy4blast.sqlite3` (whose root has the one
   child 131567), LOSAT's `-taxids 1` writes NCBI's bytes of `-taxids 131567`
   (`LOSAT/tests/db_taxid_root.rs`, rows `db.tf.bp.root.6` and `db.tf.tbn.root.6` of
   `LOSAT/tests/db_search_manifest/`, BLASTP and TBLASTN) and fewer than the whole
   database, whose two taxids outside the file are not descendants.
   Version 1.8 (line BDC's BDC3a, the same standing delegation): `losat blastdbcmd -taxids`
   and `-taxidlist` expand their taxids with the same walk. NCBI's `blastdbcmd -db
   18S_fungal_sequences -taxids 1` writes nothing in 120 s (line BDC's probe J120, and K13 on
   a small taxonomy); LOSAT writes NCBI's bytes of the same run with the explicit list of
   every taxid that taxid 1 reaches and `-no_taxid_expansion` (K01, K25:
   `dropin-survey-20261029/BDC/probe/make_expected.py`, `LOSAT/tests/blastdbcmd_discovery.rs`
   case `K13`).

## Approved exception (`-seqidlist` and `-negative_seqidlist` with `-db`, version 1.6)

Version 1.6, decided under the maintainer's standing delegation on 2026-10-29 (line DI-8
`DECISIONS.md` 3, of the kind the maintainer approved for the survey on 2026-10-29: NCBI
defects are fixed and presented as differences); its PR (DI-8b) asks the maintainer to accept it.

6. **An accession stored under more than about 500 sequences of a version 5 database.**
   LMDB keeps the OIDs of a duplicated accession as duplicates of one key of `acc2oid`, in a
   sub-database once they are many, and NCBI's lookup (`CSeqDBLMDB::GetOids`,
   `seqdb_lmdb.cpp:360-400`, `MDB_SET` then `MDB_GET_CURRENT`) fails on such a key: a seqid
   list naming the accession stops with `BLAST Database error: Accessions to Oids lookup
   error in <db>`, exit 2 (NCBI 2.17.0 on a database with 501 copies of `DUP00001`, positive
   and negative lists, probe S03 and S62 of `dropin-survey-20261029/DI8`). LOSAT reads the
   first duplicate as for any other key (`LmdbEntry::accessions_to_oids`, `first_oid`), as
   NCBI does for an accession of fewer sequences: a positive list keeps that first OID, a
   negative list judges it, and the other ids of the list are unaffected. `blastdbcmd -entry`
   and `-entry_batch` stop on such a key the same way (`Error: [blastdbcmd] Accessions to Oids
   lookup error in <db>`, exit 1, probes H107-H109 of `dropin-survey-20261029/BDC/probe`): LOSAT
   answers with the OIDs of the key in LMDB's order, what NCBI gives for a key of fewer
   sequences (H110), and with the first of them in the bulk lookup of `-entry_batch`
   (`LmdbEntry::oids`, `IdLookup`; `LOSAT/tests/blastdbcmd_retrieve.rs`).

## Approved difference (`gb`, `emb` and `dbj` ids, version 1.10)

The maintainer accepted this on 2026-10-21 for the tabular `sseqid` of a database with parsed
Seq-ids (line DBID `DECISIONS.md` 5). NCBI's result is deterministic, so rule 2 would otherwise
reproduce it. The other places the same table reaches (below) follow under the maintainer's
standing delegation, which the maintainer accepted on 2026-11-02 with version 1.10.

NCBI's Seq-id table keeps these three types in one tree for the whole process
(`c++/src/objects/seq/seq_id_tree.cpp:255-259`, `CSeq_id_GB_Tree` `:1707-1719`). While it lives,
an entry holds the ids of one prefix, digit count and version (`CSeq_id_Textseq_Info::ParseAcc`,
`:825-880`) with the type of the id that made it (`FindOrCreate`, `:1161-1177`;
`GetPackedSeqId`, `:938-968`), so NCBI can print `gb|AF000002.1|` for a sequence stored as
`emb|AF000002.1|`. LOSAT prints the stored type. That is NCBI's own output with
`OBJECTS_PACK_TEXTID=0`, which turns the packing off (the probes compared them in every format
in BLASTN and BLASTP). LOSAT does not read this variable, and it refuses the registry entry
`[OBJECTS] PACK_TEXTID`, in a file or in the environment. `tpg`, `tpe`, `tpd`, `ref`, `sp` and
`gpp` ids keep their type in NCBI, and so do an id of another version or digit count, an id
with a locus name and an accession with a prefix longer than 4 characters.

Where it shows, in all programs that take the option: the tabular `sseqid` with `-db` (the default
fields, `sacc`, `saccver`, `sallseqid` and the other formats agree), `qseqid` with
`-parse_deflines`, and with `-subject -parse_deflines` also `sseqid`, `sallseqid` and the
Seq-align id of outfmt 5, 8, 9, 11, 12 and 13 to 16. It needs definition lines that give one key two of these types: in a copy of NCBI's
core_nt, the 57,267,555 such ids of its 145,503,554 definition lines have no such key (a scan of
2026-11-02), so the outputs agree there. Probes: `blastdb-ids-20261021/R/id9` of the line's
records.

## Exception 2 on a database with parsed Seq-ids (version 1.10)

Decided under the maintainer's standing delegation (line DBID `DECISIONS.md` 10) and accepted by
the maintainer on 2026-11-02 with version 1.10. For a nucleotide record `lcl|Subject_N , ,` NCBI's
tabular id is the first word of the generated title, so BLASTN, TBLASTN and TBLASTX crash in
outfmt 6 as well as in outfmt 0 (SIGSEGV). LOSAT stops at the end of the string and writes `,`,
which is what NCBI writes for the nearby record `lcl|Subject_N , x`. Other ids with such a title
do not build it in the tabular formats, and NCBI and LOSAT agree there, as exception 2 says.

## Approved exception (makeblastdb, a compressed `-in` file of several streams, version 1.11)

7. **A `-in` file of several bzip2 streams or zstd frames.** makeblastdb reads an `-in` name
   ending in `.gz`, `.bz2` or `.zst` through `CDecompressIStream` (`makeblastdb.cpp:1037-1050`).
   Its gzip stream reads every member of a concatenated file (`fGZip` holds
   `fAllowConcatenatedGZip`), but its bzip2 and zstd streams stop after the first stream or
   frame and drop the rest without a message: `cat a.fa.bz2 b.fa.bz2` and `cat a.fa.zst
   b.fa.zst` each build a database of the first file's sequence only (probe 2026-11-04,
   `docs/evidence/compressed_input/README.md`; pbzip2 and pzstd write such files for any input
   above their block size). A valid result exists and can be checked: the database makeblastdb
   builds from the same bytes in one stream, or from the plain file. LOSAT reads every stream
   and frame (`blastinput::compressed`), so `losat makeblastdb -in parts.fna.bz2` equals
   `makeblastdb -in parts.fna`; a single-stream file builds the same files as makeblastdb
   (`LOSAT/tests/cli_smoke.rs`, `makeblastdb_reads_compressed_inputs_as_makeblastdb_does`).
   Line GZ, DECISIONS 9 (Owner-delegated).

## Approved exception (TBLASTX, a cutoff quotient past the range of a double, version 1.12)

Accepted by the maintainer on 2026-11-12 (line PS `DECISIONS.md` 36).

8. **A cutoff whose quotient `K * searchsp / E` is past the range of a double.** NCBI computes a
   cutoff score as `(Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ))`
   (`blast_stat.c:4059-4061`, with `E` at least `kSmallFloat`). When the quotient overflows to
   infinity, the conversion to `Int4` is undefined in C; on x86-64 it gives `INT_MIN`, so
   `BLAST_Cutoffs` keeps a cutoff of 1, and TBLASTX's word cutoff (the smallest of the update
   cutoff, the gap trigger and `cutoff_score_max`) becomes 1, so every ungapped HSP reaches the
   sum statistics. TBLASTX asks this cutoff for an E-value of `1e-300` (`CUTOFF_E_TBLASTX`) over
   `MIN(subject length, query length) * subject length` (`blast_parameters.c:348-363`), so with
   the default options the quotient overflows once K times that search space passes about 1.8e11
   (the largest double times `kSmallFloat`), which ordinary genome pairs reach. On E. coli K-12
   MG1655 positions 1 to 1,200,000 against the 5.5 Mb E. coli O157:H7 Sakai genome (outfmt 6,
   one thread), NCBI BLAST+ 2.17.0 stopped with SIGSEGV in `BLAST_LinkHsps` after 12 min 58 s
   and 43.8 GB under a 48 GB address-space limit, and
   wrote nothing in 4 hours under a 300 GB limit (91.8 GB resident after 1.5 hours). With an
   extreme `-searchsp` and `-evalue` on a small input NCBI does write rows: 320 for a
   6,077-codon sequence against itself with `-searchsp 9223372036854775807 -evalue 1e-300`.
   A valid result exists and can be checked: LOSAT keeps the finite cutoff of the mathematically
   equal `ln K + ln searchsp - ln E` (`algorithm::tblastx::ncbi_cutoffs`), and on a nearby pair
   whose quotient stays finite (the same 1.2 Mb query against a 2 Mb subject made of 30 kb of
   the Sakai genome inside N) NCBI and LOSAT write the same 457 rows byte for byte. On the
   genome pair above LOSAT writes 88,562 rows in 2 minutes, and on the small input 262 rows. The
   maintainer's rule of 2026-10-10 (no reproduction of NCBI's undefined behaviour or integer
   overflow) applies, also on these ordinary inputs.

## Approved exception (`losat rpsbproc`, PSSM-IDs of 2^31 and more, version 1.13)

Version 1.13, decided under the maintainer's standing delegation on 2026-11-11 (line PS
`DECISIONS.md` 27, addendum of the PS16a review), under the rule of 2026-10-10 that LOSAT does
not reproduce a signed `Int4` overflow (undefined behaviour in C++), as exception 4 widens the
culling limit; its PR (PS16a) asks the maintainer to accept it, as version 1.5's did.

9. **A PSSM-ID of 2^31 or more in rpsbproc's data.** rpsbproc 0.5.1 orders its maps of domains
   and superfamilies by the `int` difference of two ids (`TSortByPssmId`, `(int)v1.m_uiPssmId -
   (int)v2.m_uiPssmId`, `common/basedata.hpp:796-807`; `TSortByClusterId`, `:809-820`), used as
   `m_comp(*p1, *p2) < 0` by its `std::sort` and by the sign of the difference in its binary
   search (`common/ptrmap.hpp:63-70`, `:109-127`). Two ids whose `int` values are 2^31 or more
   apart overflow it. NCBI's 2025-05-13 binary (GCC, no `-fwrapv`) sorts by the signed values
   (the compiler folds `a - b < 0` to `a < b`) but searches with the wrapped difference: with
   the fixture's `cddid.tbl` plus eight domain rows of PSSM-IDs between 2^31 - 2.5e6 and 2^31 +
   1.3e6 (and in runs 2 and 3 a row of PSSM-ID 2^31), it prints the hit of PSSM-ID 488084 (run
   1) and of 2148253328 (run 3) but not that of 2^31 (run 2; PS16a review oracle A2); with 1 to
   30 rows of PSSM-IDs between 2^31 - 6e8 and 2^32 that no archive hits, it prints the output of
   the data without them (oracle A4-1 cases `big3`, `big86` and `big88`). LOSAT compares the signed values
   in the sort and in the search (`rpsbproc::data`, `int_order`), so it finds every row, the
   one of 2^31 too; elsewhere its maps are NCBI's. The valid result is the one NCBI gives for the
   row under a PSSM-ID that does not overflow the difference (only the PSSM-ID column, which
   prints 0 from `INVALIDPSSMID` on, would differ). No CDD release has a PSSM-ID of 2^31 or
   more (v3.22's largest is 533,962; the pseudo PSSM-IDs of accessions are about 1e8).
   Unit test `pssm_ids_from_2_to_the_31_are_found` (`LOSAT/src/rpsbproc/data_tests.rs`).

## Approved exception (`losat rpsbproc`, a motif check past the stored positions, version 1.14)

Version 1.14, decided by the Owner on 2026-11-11 (line PS `DECISIONS.md` 25, Q2 = A, answer
「すべて推奨通り」; numbered after 1.13); its PR (PS16b) asks the maintainer to accept it as the PR
of exception 9 (PS16a) did.

10. **`MotifCheck` reads `st_motif_pos` past its end.** rpsbproc 0.5.1 keeps the positions of
    the first motif-bearing site that it checks in a function-static `vector<SeqPos_t>
    st_motif_pos` (`common/basedata.cpp:671-676`) and compares the mapped positions of every
    later motif-bearing site, of any domain and query of the process, with it: for each of the
    `mtfLen` positions of the motif's minimal map, `rMappedRes[mapIdx].ori > st_motif_pos[mtfIdx]`
    (`:704`). When the motif's minimal map is longer than the stored list (`mtfIdx >=
    st_motif_pos.size()`) the read is outside the vector (C++ undefined behavior: heap bytes),
    so NCBI's result depends on the heap. LOSAT treats the check as failed (returns 1, the site
    is dropped), as the loop does when the mapped residues are exhausted (`:714-715`). The list
    itself is NCBI's process-wide state, reproduced per run (rule 2). `MotifState::motif_check`
    (`LOSAT/src/rpsbproc/sites.rs`); unit test `a_minimal_map_longer_than_the_stored_positions_fails_the_check`.
    The PR gives the number of occurrences on E. coli K-12 MG1655 x CDD v3.22.

## Approved exception (a REAL whose exponent overflows NCBI's `int`, version 1.19)

Version 1.19, decided under the maintainer's standing delegation on 2026-11-14 (line SW item
SW8 B28, `sw8/DECISIONS.md` 58), under the rule of 2026-10-10 that LOSAT does not reproduce a
signed `Int4` overflow (undefined behaviour in C++), as exceptions 4 and 9 widen theirs; its PR
asks the maintainer to accept it, as the PRs of versions 1.5 and 1.13 did.

11. **A REAL whose exponent overflows `int` in `NStr::StringToDoublePosix`.** NCBI's own decimal
    reader (`corelib/ncbistr.cpp:985-1278`), which reads a REAL of the ASN.1 text form (a PSSM of
    `-in_pssm` in psiblast and tblastn, an archive of rpsbproc) and of the XML form (rpsbproc's
    BlastOutput input), keeps the position of the point and the exponent in `int`
    (`dot_position - digits`, `exponent - expvalue` and `exponent + expvalue`, `:1113` and
    `:1155`); a text whose exponent passes `int`'s range, such as `1.55e-2147483647`, overflows
    it. NCBI's corelib wraps it (to 2147483647 there) and gives the infinity, `DBL_MAX` with
    `fDecimalPosixFinite`, the flag of the word form `1.55e-2147483647` of both readers. LOSAT
    computes the exponent in 64 bits (`ncbi_str::double::string_to_double_posix`), so the
    number is read with its exact exponent: `DBL_MIN` with `fDecimalPosixFinite`, 0 without.
    Every exponent that `int` holds gives NCBI's result: of 900,184 texts (184 edge cases,
    mantissas of up to 5,700 bytes, 900,000 random) read by NCBI's corelib (598d8ae6) and by
    LOSAT, only this one differs (`sw8/nfix/probe/x87`). NCBI's writers print an exponent of a
    few hundred at most. Unit tests `refusals_wait_for_the_whole_object`
    (`LOSAT/src/rpsbproc/archive_tests.rs`, the archive) and `numbers` (`LOSAT/src/report/serial/xml_reader_tests.rs`,
    the XML input).

## Approved exception (`losat rpsbproc`, coordinate arithmetic past NCBI's `int`, version 1.26)

Version 1.26. The Owner decided on 2026-11-14 (line PS `DECISIONS.md` 57, Q1 = B, words 「Q1=B」)
the question of PS16c's PR (#1081) and of PS16d's review: a hit whose coordinate overflows NCBI's
`int` (`m_iFrom + COORDSBASE` and `m_iTo + COORDSBASE` of the DOMAINS and SUPERFAMILIES lines, with
`Calculate`'s overlap sets of the hit), for protein and nucleotide queries alike, and the XML input's
coordinates (`- COORDSBASE` of an HSP's coordinates and what follows from it), are widened under one
approved exception, one rule for one expression, under the Owner's rule of 2026-10-10 that LOSAT does
not reproduce a signed `Int4` overflow (undefined behaviour in C++), as exception 4 widens the culling
limit. The other `int` overflows of the same coordinate arithmetic, which the question did not name,
were added to this exception after a code check by an Owner-delegated decision (line PS
`DECISIONS.md` 58, standing instruction 2026-09-29, restated 2026-10-07): the site's `xres +
COORDSBASE` and the motif's `GetLeft() + COORDSBASE`, `MinusPr2PlusNA`'s `pr - 1`, the set algebra of
`CSegSet`, `CalcMasterGaps`, `GetTranslatedPosMap`'s `++c`, a `Std-seg`'s `t - f + 1` and `m_iTo + 1`,
and the data files' `SetData` `seg.to + 1` and `GetMotifResPos` `++c` (rejected there: leaving them
refused, two rules for one kind of overflow); and, by line PS `DECISIONS.md` 59 (Owner-delegated,
after the delta review), `++alnResIdx` and `ParseAlignStrings`' `++start` of an `Hsp_qseq` of 2^31
letters or more (rejected there: refusing such an `Hsp_qseq`, by its length or by the memory rule of
`DECISIONS.md` 40, which measures what NCBI builds from a number or pattern in the input, not the
strings the input carries). The PR (PS16-Q1) asks the maintainer to accept version
1.26 with both parts. It moves the coordinate overflows that versions 1.16 and 1.17 rejected (rule 3)
here, and names the widening that PS16a made for protein queries.

12. **rpsbproc's `int` coordinate arithmetic past `int`'s range.** rpsbproc 0.5.1 keeps positions in
    `SeqPos_t` (`int`) and adds to them in `int`. A hit that starts or ends at base 2^31 - 1, a site or
    motif mapped there, a `Std-seg` subject interval of 2^31 positions or more, or a coordinate of
    -2147483648 in a BlastOutput XML document overflows that arithmetic. NCBI's two 0.5.1 builds wrap
    it: a hit of rpstblastn that ends at 2^31 - 1 prints `-2147483648` as its end (PS16c review oracle
    `plusmax`, PS16-Q1 oracles `ntplusmax`, `ntminusmax`), and so does an ANK motif of sd00045 that
    starts there (`ntmotifmax`); an `Hsp_query-from` of -2147483648 throws `basic_string::replace` at
    position 2147483647 (`xqfmin`). LOSAT computes every such `int` operation as its mathematical
    value (`i64`) and keeps the value where NCBI keeps the `int`; what NCBI defines after it stays
    NCBI's: a conversion to `unsigned int` is the value modulo 2^32, `unsigned int` arithmetic wraps,
    a `size_t` is the value modulo 2^64, a print prints the value. A value that then reads a string
    outside its range is refused as before (rule 3). The expressions, for protein and nucleotide
    queries alike:
    - `m_iFrom + COORDSBASE` and `m_iTo + COORDSBASE` of the DOMAINS and SUPERFAMILIES lines, the
      residues `+ COORDSBASE` of the SITES and MOTIFS lines, and a nucleotide query's `xres +
      COORDSBASE` of a site's codon and `GetLeft() + COORDSBASE` of a motif (one owner,
      `rpsbproc::printed_position`);
    - the set algebra of `CSegSet` (`from - 1`, `to + 1`, the total residues, the `ori_from` sums),
      `MapSegSet`'s `diff` and `end + 1`, `CalcMasterGaps`' gap start and length, `GetTranslatedPosMap`'s
      `++c` and `ori_from + inc` (PS16a's widening, `rpsbproc::segset::Seg`; `Calculate`'s overlap
      sets of a hit that ends at 2^31 - 1);
    - `ParseAlignSegs` of a `Std-seg`: a subject interval's `t - f + 1` (stored in the `unsigned int`
      lengths) and a minus-strand hit's `m_iTo + 1` (divided as an `unsigned int`); `MinusPr2PlusNA`'s
      `pr - 1` at `INT_MIN`;
    - the XML input's `- COORDSBASE` of an `Hsp_query-from`, `Hsp_query-to` or `Hsp_hit-from` of
      -2147483648, "see if NA"'s `m_iTo - m_iFrom + 1`, `ParseAlignStrings`' `++qf`, `--qf`, `++hf`
      and `--hf` (`m_vecSStarts` holds the values: `rpsbproc::domains::SubjectStarts`), and the
      alignment strings' positions `++alnResIdx` and `++start` in an `Hsp_qseq` of 2^31 letters or more
      (compared with the `unsigned int` length, used as `size_t` indices and in the `unsigned int`
      `start - segLen`; the segment lengths are `unsigned int` already);
    - the data files' site positions (a hand-edited `cddannot.dat` with a position of 2^31 - 1):
      `SetData`'s `sorted[i] > seg.to + 1` (`common/impl/segset_impl.hpp:25-55`,
      `rpsbproc::data::site_data_segs`) and `GetMotifResPos`'s `for (SeqPos_t c = v->from; c <= v->to;
      ++c)` (`common/basedata.cpp:622`, `rpsbproc::sites::MotifState::motif_check`), whose `++c` past
      2^31 - 1 would never end the loop in a wrap; LOSAT computed both in `i64` before version 1.26.

    Evidence: for each kind, an input at the limit run by both 0.5.1 builds (the same bytes in all
    33 jobs) and a nearby input whose coordinates are moved away from the limit by a known shift;
    LOSAT's output of the limit equals NCBI's output of the nearby input up to that shift, and LOSAT's
    output of the nearby input equals NCBI's (task folder `PR/tmp/impl16q1/oracle/INDEX.md`; the
    fixture rows of `use` `widened` in `LOSAT/tests/fixtures/rpsbproc/CASES.tsv`, each with its nearby
    `compare` row, README.md there):
    - `ntplusmax` (nt67's first hit ends at 2^31 - 1; `-t doms -m std`, `-t doms -m rep -f`, `-t both
      -m std`): NCBI prints `-2147483648`, LOSAT 2147483648; the hit 3 bases shorter
      (`ntplusmax3`) prints 2147483645 in NCBI, nothing else differs.
    - `ntminusmax` (nt67's first minus-strand hit ends at 2^31 - 1): `alignedLen` is 715,827,882
      codons as for the hit one base shorter (`ntminusmax1`, which NCBI prints with end 2147483647);
      LOSAT prints 2147483648, NCBI `-2147483648`.
    - `ntslen` (a subject interval from 2^32 - 5 to 2^31 - 2, `t - f + 1` = 2^31 + 4): NCBI's and
      LOSAT's bytes are the same, and the same as for the interval 5 positions shorter (`ntslen1`, no
      overflow).
    - `ntmotifmax` (nt68's sd00045 hit moved so that an ANK motif starts at base 2^31 - 1, with the
      oracle data `ntmot` that keeps its motifs only): NCBI prints `-2147483648`, LOSAT 2147483648;
      the hit 3 bases lower (`ntmotifmax3`) prints 2147483645 and every other number 3 lower.
    - `xqfmin` (an `Hsp_query-from` of -2147483648): NCBI's `replace` throws at position 2147483647,
      LOSAT's at 18446744071562067967, one below NCBI's for -2147483647 (`xqfmin1`).
    - `xqtmin` (`Hsp_query-to`): NCBI keeps the range from 4 to 2^31 - 1 (its wrap) and prints the hit
      in `-m rep`; LOSAT's range ends at -2^31 - 1, below its start, as NCBI's for -2147483647
      (`xqtmin1`): the hit is no representative and prints `-2147483648` in `-m full`, one below
      `xqtmin1`'s.
    - `xhfmin` (`Hsp_hit-from`): NCBI's wrapped domain start 2^31 - 1 marks the hit `N`; LOSAT's
      -2^31 - 1 does not, as NCBI for -2147483647 (`xhfmin1`).
    - `xhfmax` (`++hf` past 2^31 - 1 after a query gap) and `xspan` (a span of 2^31 for a query of
      length 0): NCBI's bytes, LOSAT's and the nearby input's (`xhfmax5`, 5 lower; `xspan1`, a span of
      2^31 - 1) are the same.
    - Nucleotide queries of the XML input (the PS16-Q1 review's oracle `nx`: `x_ntpair.xml` with the
      first minus-frame or plus-frame HSP's `Hsp_query-from`, `Hsp_query-to` or `Hsp_hit-from` at
      -2147483648, 12 limit jobs, through `(m_iTo / 3) << 2`, `PlusNA2PlusPr`, `PlusNA2MinusPr` and `hf`
      on both strands): LOSAT's output equals NCBI's for -2147483647, byte for byte or with the one
      printed coordinate one lower; fixture rows `x.h_ntmqfmin.rep` (NCBI drops the hit its wrapped
      range reclassifies; LOSAT prints `-2147483648` where NCBI prints `-2147483647` for the nearby
      input) and `x.h_ntphfmin.rep` (NCBI marks the hit `N`, LOSAT and the nearby input do not).
    The fixture test `widened_rows_are_their_nearby_rows_moved_back` (`LOSAT/tests/rpsbproc_fixtures.rs`)
    checks every `widened` row against NCBI's bytes of its nearby row and the shift of
    `widened/PAIRS.tsv`.
    `pr - 1` of `INT_MIN` and a site's `xres + 1` at 2^31 - 1 are not reached end to end: a site
    reads its residue in the frame's translation first, which no query below 2^31 bases holds there
    (refused, rule 3), and a minus-strand motif reaches residue `INT_MIN` only through a segment of
    about 7e8 residues; their values are unit tested against the residue next to them (the codon 3
    bases away).
    Unit tests (`LOSAT/src/rpsbproc/`): `mod.rs` `printed_positions_past_int`, `query.rs`
    `refusals_of_translated_hits`, `domains_tests.rs` `std_seg_forms`,
    `a_minus_strand_hit_that_ends_at_int_max`, `ends_that_overflow_int`, `frames_tests.rs`
    `frame_conversions`, `nucleotide_tests.rs` `codons_and_their_overflow`,
    `motif_ends_and_their_overflow`, `offline_tests.rs` `coordinates_of_int_min`,
    `span_that_overflows_int`, `domain_positions_that_overflow_int`, `string_reads_at_int_positions`,
    `see_if_na_counts_and_refuses`.
    An `Hsp_qseq` of 2^31 letters or more (PS16-Q1 fix round 2 oracle `q31`, ws): the fixture's `x_one.xml`
    with L leading query gaps in `Hsp_qseq` and L letters in `Hsp_hseq` before the HSP's 70 columns, which
    `++alnResIdx` trims, so the HSP's segments and residues are those of the 70 columns whatever L is. For
    L = 2^31 - 128 (the nearby input, no overflow), LOSAT's stdout, stderr and exit code equal both NCBI 0.5.1
    builds'. For L = 2^31 + 8 (the limit), LOSAT, in a build with overflow checks, prints the same bytes
    except the line that echoes the input file's name. NCBI's binary prints the same there: its compiled
    loop shows no wrap in this output, which the undefined behaviour allows.

## Reproduced as NCBI (rule 2)

- `losat rpsbproc`, the end of an aligned segment when a site is mapped (`MapSegSet`'s `SeqPos_t end
  = m_vecSStarts[idx] + m_vecLens[idx] - 1`, `common/basedata.cpp:1180-1250`): an `unsigned int` sum
  stored in an `int`, which wraps for a segment that passes model position 2^31 - 1, so the sites
  beyond it are not mapped. Version 1.26 (PS16-Q1) reproduces it; PS16a had added it in `i64` and
  printed a SITES line that NCBI does not (the PS16-Q1 review's oracle `mse`: `three.asn` with the
  subject starts of its 199908 hits and that model's site positions moved by K; at K = 2147483586 both
  0.5.1 builds and LOSAT print no SITES line, at K = 2147483576 all print it). Fixture rows
  `mse86.both`, `mse76.both`; unit test `mapped_segment_ends_wrap_as_unsigned_int`
  (`LOSAT/src/rpsbproc/domains_tests.rs`).

- A text seqid list whose ids all go in NCBI's rewriting for a version 5 database (gi ids
  only): NCBI prints `Empty seqid list`. A negative list excludes nothing. A positive list
  clears every OID (`x_ApplyUserGiList`): no sequence is searched, and each query gets
  `Filtering resulted in an empty database.`, except where NCBI's OID bit set is its special
  all-set one: the clearing writes into a stale array, and the whole database is searched.
  That is the case when the `-db` names reach one volume and at least one path of alias files
  to it crosses no `FIRST_OID`-`LAST_OID` range that cuts the volume (two ranges whose union is
  the volume do not count; two or more volumes never do) (line DI-8 `DECISIONS.md` 18; line
  DBID item gisl, 2026-11-15).

- `CHUNK_SIZE=1000` without `BATCH_SIZE`: the batch size becomes 0, NCBI prints the
  outfmt 0 prolog and then `BLAST engine error: Empty CBlastQueryVector`, exit 3.
- `-max_target_seqs` from 2^30 to 2^31-51: the preliminary hit list size
  (`MIN(MAX(2*h,10),h+50)` in `Int4`) wraps and becomes 10.
- Defects without an effect on the output that LOSAT follows: the negative `Int4`
  offset of `BlastGetStartForGappedAlignmentNucl`, the out-of-range read of
  `hard_ranges[1]` for the last subject chunk, the out-of-range `double` to `Int4`
  conversion (INT_MIN), and the catch in `BlastSetupPreliminarySearchEx` that never
  matches.

## Explicit rejections (rule 3)

- A scoring system without Karlin-Altschul values when an invalid query (such as all N)
  comes first in its batch: NCBI goes on without a gapped block and crashes (SIGSEGV).
- `-max_target_seqs` above 2^31-51: the preliminary hit list size becomes negative and
  NCBI crashes.
- A reward of 32767 or a penalty of -32768 (after NCBI's 16-bit conversion):
  `BLAST_SCORE_MAX`/`BLAST_SCORE_MIN` fall outside NCBI's score range; a reward of 32767
  is counted past the end of the frequency array (invalid queries or a crash), a penalty
  of -32768 makes every query invalid. Version 1.1 (maintainer, 2026-10-02): a reward of
  32768 or more, which NCBI's 16-bit field wraps to 0 or below (NCBI: every query invalid,
  no hits, a deterministic result), stays an explicit rejection too.
- Subjects of 2^31 letters or more in total: NCBI's 32-bit total wraps.
- megablast gap costs above 32767: NCBI's 32-bit greedy distance arithmetic can wrap.
- A negative `CHUNK_SIZE` above a negative `OVERLAP_CHUNK_SIZE` that splits a query batch
  (`SplitQuery_CalculateNumChunks` gives more than one chunk): NCBI's `size_t` chunk ranges
  wrap; NCBI stops as in exception 1 where a chunk would be split again, and otherwise
  searches chunk ranges with gaps between them (hits are lost). Version 1.1 (maintainer,
  2026-10-02): the second case, although NCBI's result is deterministic, stays an explicit
  rejection (a setting without a practical use; reproducing it would need NCBI's wrapped
  `size_t` arithmetic through the chunk ranges, masks and merge). Where such a pair does
  not split the batch (for example `-1` above `-2147483648`), NCBI searches as without the
  variables and so does LOSAT (fixtures `env.negative_pair_*`).

- A `CHUNK_SIZE`/`OVERLAP_CHUNK_SIZE` pair whose chunk ranges leave a chunk without a
  query (an overlap close to or above the chunk size): NCBI's `CQuerySplitter::Split`
  (`split_query_cxx.cpp:872-880`) or the chunk's search stops with a null-pointer
  `CCoreException`; LOSAT rejects "a query chunk without a query" (round-3 audit: 27
  oracle cases, the reason holds in all).

- A sequence that a seqid list keeps while matching none of its definition lines (a pdb
  molecule named without its chain, a local id in other letters, a Swiss-Prot name on a
  version 4 database, a sequence whose every line has only a gi or a ti under a negative
  list; version 1.29: a list's `gnl|BL_ORD_ID|n` that the accession table resolves by the
  stored number, of `makeblastdb -input_type blastdb` or of a FASTA id `gnl|BL_ORD_ID|n`, to a
  sequence in a volume whose first OID is not 0 (a later volume, or a later database of the
  list), whose lines NCBI's header read tests with the volume's first OID added
  (`seqdbvol.cpp:2222-2236`), and a negative list's line removed there the same way; where such
  a sequence keeps another line, LOSAT writes the lines NCBI keeps): NCBI's
  filtered header leaves no line, `GetBioseq` returns null
  (`seqdbvol.cpp:1370-1376`) and the sequence's hits leave every report without notice
  while its length counts in the statistics. LOSAT refuses the search with the sequence
  named (`the seqid list keeps sequence N of the volume but none of its definition lines, and
  NCBI leaves its hits out of the report without notice`): a positive list before the
  search, a negative list when such a sequence is to be reported (version 1.6, line DI-8
  `DECISIONS.md` 17 and 22). The result is deterministic, but it loses hits without a
  word, which no report can show; a user who wants NCBI's report names the chain or the
  stored id.

- A version 5 database list (several names in `-db`, or the `DBLIST` of an alias file) that
  mixes databases with and without the accession table (`acc2oid`) of `makeblastdb
  -parse_seqids`: NCBI searches it, and before it writes a batch's reports (every format but
  the Seq-annots and the archive) its BlastDb data loader looks each hit's Seq-id up in every
  LMDB file of the list (`BlastFormatter_PreFetchSequenceData`,
  `CSeqDBLMDBSet::AccessionToOids`); a hit whose Seq-id it looks up as a string (an accession, a
  name, a giim; not a gi, not `gnl|BL_ORD_ID|n`) then stops the run with `BLAST Database error:
  DB contains no accession info.` (exit 2), after the outfmt 0 header and the earlier batches. Where `PRE_FETCH_SEQS_LIMIT` turns the prefetch off, the reports write
  `Unknown` for such a hit, and then outfmt 0 to 4, 13, 14 and 17 stop with that error, outfmt 6
  crashes (`CCoreException`, exit 255), the JSON and XML2 documents 15 and 16 end with its trace
  (exit 0), and 5, 7, 10 and 18 go on (exit 0).
  Without such a hit (a `gnl|BL_ORD_ID|n` of a database without parsed Seq-ids, or of one made
  with `-input_type blastdb`, needs no accession table), and for the Seq-annots and the archive,
  NCBI runs to the end (the Seq-annots too stop where `PRE_FETCH_SEQS_LIMIT=2147483647` prefetches
  for them); `-seqidlist` and `-negative_seqidlist` of accessions stop it up front with the same
  error. LOSAT searches such a list as NCBI does (line DBID item mixdb, probe of 4,080 runs:
  6 programs, 20 formats, the prefetch on and off, 17 lists and alias files) and, where the
  batch's prefetch is on, stops as NCBI does: the prologs and the reports of the batches before
  (without the epilog, in every format that prefetches, the JSON and XML2 documents and side
  files too), the batch's reader messages and the error (line DBID item mixdb2, probe of 4,800
  runs; psiblast's PSSM drops such a sequence's letters with `Failed to retrieve sequence`, as
  NCBI's). It refuses a run where the reports read such a Seq-id with the prefetch off (`... whose
  lookup in the database fails (DB contains no accession info.; NCBI BLAST+ 2.17.0's prefetch
  stops with a BLAST Database error, and where the prefetch is off its formatters fail to read
  it)`), and a stopping batch whose other sequences' loads write texts of the
  taxonomy database or whose other lookups throw on their ids (NCBI's prefetch loads them in the
  order of its Seq-id handles, which decides what it writes before the error).
  Versions 1.9 to 1.24 refused every such list before the search (line DI-8 `DECISIONS.md` 31),
  versions 1.25 to 1.30 every run whose reports or prefetch looked such a Seq-id up.

- A negative seqid list that excludes a sequence (every id of its `.nos`/`.pos` record
  listed) of a volume that the database does not hold, while the database holds other
  volumes of the same LMDB file (`-db m.00` of a database in several volumes, or `-db "d1/m
  d2/m"` whose volumes interleave): `x_AdjustOidsOffset` turns that sequence's OID into
  `kSeqDBEntryNotFound`, and NCBI clears bit -1 of its OID mask and crashes (SIGSEGV, exit
  139; review of line DI-8f, fixture row `db.mdb.bn.negpartial.6`). LOSAT dropped the OID
  and searched; it now refuses the list (`the negative seqid list excludes sequence N of the
  LMDB file, in a volume that the database does not hold`; version 1.9). A positive list on
  such a database is NCBI's and LOSAT's alike.

- On a database with parsed Seq-ids, where NCBI stops (a reported sequence that its formatters
  cannot fetch, in the default tabular fields), pads the line without end (an outfmt 0
  description id of 65 bytes or more) or crashes (a protein title of punctuation): version 1.10,
  the table in [`PD-LOSAT-DATABASE-SEARCH`](PD-LOSAT-DATABASE-SEARCH.md).

- `-outfmt 18`, the organism report (version 1.20, 2026-11-14, line RG item T2c; all five
  programs): a query batch that NCBI does not search (no query of it has a valid context, "Could
  not calculate ungapped Karlin-Altschul parameters"): NCBI leaves the Seq-align-set of each of its
  queries null (`local_blast.cpp:204-205`), and `x_PrintTaxReport` reads the query's one
  (`blast_format.cpp:893`, `:898-899`) and stops with a null pointer
  exception and a stack trace of its build, exit status 255; LOSAT writes NCBI's bytes up to the
  query's `Length=` line and then the error. Also (line RG item T2, 2026-11-13): a taxonomy name
  with `<@bl_taxid@>` (NCBI prints the uninitialised `STaxInfo::blNameTaxid`, `taxFormat.hpp:142`)
  and a `./.ncbirc` with a syntax error (NCBI's report stops at it with exit status 255).

- `losat rpsbproc` (version 1.15, 2026-11-12, line PS items PS16a and PS16b, line PS
  `DECISIONS.md` 37 and 40, Owner-delegated): where NCBI's rpsbproc 0.5.1 reads memory outside
  its strings or its result depends on the host's memory, LOSAT refuses (exit 1) and names the
  data line, site or position; the report keeps the lines before it, and no part of a line.
  - A `cdtrack.txt` line on which NCBI's loader reads past the line's terminating 0 (a line that
    ends without a root or after a field it skips; PS16a).
  - A motif match that reads beyond the residues of a site, where NCBI's answer depends on the
    bytes it reads there: LOSAT runs the match once with every such read failing and once with
    every such read passing, and gives the answer where the two agree (NCBI's search tries every
    count, so that answer holds for any bytes); where they differ, NCBI's answer depends on its
    stack's bytes (`Hx(0,10)W` on `HKY`). Also refused where the two agree, if a read may leave
    NCBI's buffer of the residues (the capacity of `minMap`; past it a read may fault, by the
    memory layout; `DECISIONS.md` 40).
  - A motif whose check makes NCBI hold 1 GiB or more at once in allocations below 2^47 bytes:
    the minimal map's string with its old buffer while it grows, or the match records of
    `curr.flex` with the copies of them on NCBI's stack and the string of the residues. NCBI's
    allocations then succeed or throw `std::bad_alloc` by the host's memory and overcommit: the
    oracle `H-{<}(0,67108863)-{<}(0,67108863)-K` (1.5 GiB of records in allocations of 512 MB)
    and `Hx(536870910)K` (a map of 512 MB, then 1 GiB) exit 0 under a 4 GiB address space and
    with `std::bad_alloc` under 1.5 GiB in both 0.5.1 builds. A repeat that does not advance
    (`<*`) grows NCBI's records until the memory runs out. One allocation of 2^47 bytes or more
    throws `std::bad_alloc` on every host, and a map longer than the string's `max_size()`
    (2^63 - 1 in NCBI's C++20 builds) throws `basic_string::_M_replace_aux`: there LOSAT ends as
    NCBI does (exit 3), and below 1 GiB it gives NCBI's result.
  - A site or motif position outside the query's residues other than their end, also a query
    without residues: NCBI reads memory outside the string (its two 0.5.1 builds print different
    bytes for a query without `seq-data`).
  - Nucleotide queries (rpstblastn's archives; version 1.16, 2026-11-13, line PS item PS16c,
    `DECISIONS.md` 37 and 40): a first `Std-seg` with fewer than 2 ids, or a `Std-seg` without
    locations or whose only location is not `null` (`ParseAlignSegs` reads past them; PS16c oracle
    `noids1`: segmentation fault in both 0.5.1 builds; a `Std-seg` whose only location is `null`
    is skipped, as NCBI's `seqlocs[0]->IsNull() || ...` skips it: review oracle `nullone`), a
    minus-strand alignment with a subject location other than `int` or `pnt` (its loop reads
    `m_vecLens` past its end); a query whose hits carry sites or motifs and whose bases, or the
    IUPACna letters of their reverse complement, are fewer than 3 (`TranslateAll` erases the first
    byte of an empty string, and libstdc++ moves `SIZE_MAX` bytes: oracle `short2`, segmentation
    fault); a site at a residue outside its reading frame's translation other than its end
    (`operator[]` past the string; oracle `natoprot` and the review's `otherstrand`, where the two
    builds print different bytes). Version 1.26 moved the arithmetic of the coordinates that
    overflows NCBI's `int` (a subject interval's length, `m_iFrom + 1` and `m_iTo + 1` of a hit that
    starts or ends at 2^31 - 1, `pr - 1` of a residue at `INT_MIN`, `xres + 1` of a site's or motif's
    base) to exception 12. Every case is made by hand; none occurs in the rpstblastn archives of
    PS16c probe 1 (891 jobs).
  - The BlastOutput XML input (`-outfmt 5` of rpsblast and rpstblastn; version 1.17, 2026-11-13,
    line PS item PS16d, `DECISIONS.md` 37 and 40): an `Hsp_hseq` shorter than its `Hsp_qseq` by
    more than one letter (`ParseBlastOutput` reads `hseq` at the indices of `qseq`, past its
    terminating 0), an alignment that writes a residue other than NUL at or past the end of NCBI's
    string of the query's frame (`operator[]`; PS16d oracle `n_qlen2`, a nucleotide query of 2 bases
    whose hits overwrite the string's terminator in both 0.5.1 builds; a NUL written at the end is
    the defined write of the terminator and is reproduced), a query whose strings NCBI holds as 1
    GiB or more at once (`NACommit`'s six frames, `PRCommit`'s residues and their copy, in
    `unsigned int` lengths: a nucleotide query of 0 or 1 base gives frames of more than a billion
    residues; oracles `n_qlen1`, `n_noqlen` and `p_qlenneg` end with `basic_string::replace` or
    `std::bad_alloc` by the address space), an `Hsp_qseq` of 2^31 letters or more that "see if NA"
    reads past its end (its `t0 = t1 + 1` stores a `size_t` of 2^31 or more in an `int`: a negative
    index) or over which that loop starts over and never ends (2^32 - 1 letters or more; version 1.26,
    which widens the rest of such an `Hsp_qseq`, exception 12). (A REAL of more than 18
    significant digits, refused in 1.17, is read as NCBI reads it since 1.19.) Of NCBI's XML
    reader, LOSAT refuses the constructs that NCBI's BLAST XML writers never produce and that NCBI
    reads with state LOSAT does not port: attributes (skipped, with NCBI's handling of `xmlns` and
    `nil`), namespace prefixes (stripped), control characters in character data (replaced by `#`
    with a warning that quotes NCBI's stack of the read). Version 1.26 moved the coordinate arithmetic
    of an HSP that overflows NCBI's `int` (`- COORDSBASE` of an `Hsp_query-from`, `Hsp_query-to` or
    `Hsp_hit-from` of `-2147483648`; `m_iTo - m_iFrom + 1` of "see if NA"; `++hf` and `--hf` of
    `ParseAlignStrings`) to exception 12. Every case is made by hand; none occurs in the 2,310 jobs
    of PS16d probe 1 (NCBI's rpsblast and rpstblastn 2.17.0 XML).
  - NaN E-values and bit scores of hits, archives and XML alike (version 1.17, the PS16d review
    S3-3; it narrows PS16a's refusal of every NaN): NCBI's `HUGE_VAL/HUGE_VAL` (sign set, printed
    `-nan`) is kept and the hits are sorted as libstdc++'s `std::sort` sorts them with
    `TDomAlignOrderByEValue`, which a NaN makes no strict weak order. Refused only where one of
    that sort's unguarded loops would read outside its vector (the comparison stays irreflexive and
    asymmetric with a NaN, so the partition's sentinels hold; the final insertion's walk can leave
    the vector only after the heap sort of the leftmost part, when it reaches the depth limit). The
    permutations of 86 NaN inputs equal libstdc++ 13's (unit vectors); the fixture rows
    `x.h_*nan*`, `nan_hit`, `nan_threshold` equal NCBI.

- TBLASTN with `-window_size 0` (version 1.18, 2026-11-14, line RG item X2, Owner-delegated): with
  a low `-threshold`, NCBI's one-hit ungapped extension can give a negative `Int4` length, which
  `BlastGetStartForGappedAlignment` reads as `Uint4` (`aa_ungapped.c:1054,1083`,
  `blast_gapalign.c:3393-3413`): it scores the first window of 11 letters from the HSP's start,
  which may pass the end of its context and its frame. Inside NCBI's buffers the window reads the
  NULLB and then the next contexts' or frames' letters, and LOSAT reads the same bytes. LOSAT
  refuses (exit 1) where the window reads outside them:
  - past the byte after the last context's NULLB of the query batch (`calloc(buflen+1)` with
    `buflen` from `QueryInfo_GetSeqBufLen`, `blast_setup_cxx.cpp:500-501`,
    `blast_query_info.c:246-250`);
  - past the NULLB after frame -3 of the subject's translation (`BLAST_GetAllTranslations` writes
    frames 1, 2, 3, -1, -2, -3 into `2*(nucl_length+1)+2` bytes and leaves the last one unset,
    `blast_util.c:1058,1063,1098-1101`);
  - with `-in_pssm`, past the PSSM's last position (`pssm->data` holds one row per position).

  It also refuses a result of such a gapped start that NCBI does not finish: two gapped starts
  past the end of their chunk of the frame (the window of a chunk that is not the frame's last),
  whose HSPs start past the preliminary interval tree's `s_max` (the chunk's length + 1,
  `blast_gapalign.c:3829-3830`) in one subject subtree. `BlastIntervalTreeAddHSP` sends both
  right at every node, splits the node `[s_max, s_max]` into `[s_max + 1, s_max]`, whose middle is
  `s_max` again, and repeats that split while its node pool doubles (`blast_itree.c:57-112,
  603-793`), so the search never returns. A gapped start at the chunk's last letter or past it can
  also give an alignment of score 0 whose right extension did not run (its start one past its end,
  kept by a cutoff of 0: `blast_gapalign.c:4290-4295`, `blast_stat.c:5236-5266`); such a range
  never straddles a node's middle, and two of them split a node `[p, p]` the same way. LOSAT stops
  where that first inverted node would be made; ranges in order inside the tree's bounds never
  make one, and only this word finder's HSPs are out of order or leave them.

  BLASTP's approved exception 3 (version 1.3), which reads such letters as the sentinel, is not
  extended to TBLASTN.

- BLASTX with `-window_size 0` (version 1.22, 2026-11-14, line RG item X2f, Owner-delegated): with
  a low `-threshold`, NCBI's one-hit ungapped extension can give a negative `Int4` length
  (`aa_ungapped.c:1054,1083`), which `BlastGetStartForGappedAlignment` reads as `Uint4`: it scores
  the first window of 11 letters from the HSP's start (`blast_gapalign.c:3394-3437`). A window that
  passes the end of its context or of its subject holds that NULLB, which scores
  `BLAST_SCORE_MIN` (`blast_stat.c:1559-1579`), so NCBI returns the HSP's start whatever bytes
  follow, and so does LOSAT. LOSAT refuses (exit 1) where the window reads outside NCBI's buffers:
  - past the byte after the last context's NULLB of the query batch (`calloc(buflen+1)` with
    `buflen` from `QueryInfo_GetSeqBufLen`, `blast_setup_cxx.cpp:500-501`,
    `blast_query_info.c:246-250`); a batch that LOSAT searches in parts (with `BATCH_SIZE` above
    10002) is one buffer, and each chunk of a split query has its own;
  - past the NULLB after a `-subject` protein (`malloc(len+2)`, `blast_setup_cxx.cpp:885-893,
    1236-1239`);
  - past the last memory page of a database volume's sequence file. NCBI maps the file whole from
    offset 0 and reads the sequences in place (`seqdbvol.cpp:1802`, `seqdbatlas.hpp:700-722`,
    `seqdbatlas.cpp:114-127`, `ncbifile.cpp:5636-5638`), and the system fills the rest of the
    file's last page with zeros, which NCBI reads as NULLBs, as LOSAT does. So only a window from
    the volume's last sequence in a file that ends within 10 bytes of a 4,096-byte page boundary
    reaches past it.

  NCBI indexes its matrix with the bytes it reads there. In 400 random searches with this setting
  and `-subject` (one to three queries of 24 to 300 bases, one to three subjects of 12 to 150
  letters, word sizes 2 to 5, `-threshold` 1 to 4), NCBI crashed (SIGSEGV) in 6, and in 2 its report
  lost a line, in one of them only in some runs; 115 of the 366 gapped searches read past a buffer.
  Under valgrind the 13 such searches checked (the 8 above among them) finish, with the output of
  reading those bytes as NULLBs. With `-db`, a window past the volume's last sequence into the rest
  of the page gave the same output in five runs and no invalid read under valgrind; one past the
  page of a 4,096-byte file made valgrind report an invalid read and crash. BLASTP's
  approved exception 3 (version 1.3), which reads such letters as the sentinel, is not extended to
  BLASTX.

- TBLASTX with `-window_size 0` (version 1.21, 2026-11-14, line RG item X2b, Owner-delegated): with
  a low `-threshold`, NCBI's one-hit word finder saves ungapped HSPs of length 0 or less
  (`aa_ungapped.c:1054,1083`), whose `Int4` ends lie at or before their starts. A subject frame
  of more than 5,000,000 letters is searched in chunks that overlap, and `Blast_HSPListsMerge`
  merges an HSP of the later chunk into one of the earlier chunk that contains its start or end
  (`s_BlastMergeTwoHSPs`, `blast_hits.c:1488-1531`) with the density `(score1 + score2) /
  (length1 + length2)`. With such an HSP the lengths can add up to 0 (an HSP of length 0 that
  both chunks find, or HSPs of lengths 2 and -2 or 1 and -1): the density is `+inf`, and the
  conversion `(int)` of its product, which sets the merged score, is undefined (C11 6.3.1.4;
  x86-64 gives `INT_MIN`, which the `MAX` with the old score then ignores). LOSAT refuses the
  search (exit 1) where that merge would happen. Every other merge of such HSPs is NCBI's (`Int4`
  lengths, a density that may be negative or raised by a negative length, an HSP of length 0
  that absorbs another and becomes positive).

- BLASTP and BLASTN `-outfmt 1` to `4` (the query-anchored views; version 1.23, 2026-11-14, line RG
  item T3a, Owner-delegated; BLASTN version 1.28, 2026-11-15, line RG item T3b, probes
  `bn.t3b.bself*`): a query read with `-parse_deflines` that is also a subject (one record, the
  same definition line and letters) is the query's own Bioseq for NCBI's `CAlnMix`, whose `Add`
  presets that subject's row to the HSP added first (`alnseq.cpp:165-186`). The first match that
  `x_Merge` places on the query (the longest aligned segment) sets that row to its own HSP
  (`alnmerger.cpp:499-501`) and goes to an extra row; when that HSP has no other segment, the
  row keeps no start, `BuildRows` leaves it and its extra rows out (`alnseq.cpp:266-283`), and
  `FillUnalignedRegions` indexes its vectors with their row index of -1 (`alnsegments.cpp:286-338`):
  NCBI writes before the vectors' buffers and aborts (`double free or corruption (out)`, exit 134;
  review probes `selfseg3`, a record whose SEG mask splits its self diagonal, and `self2`). LOSAT
  runs the same merge and refuses (exit 1) where it reaches that row, before the first byte of
  the reports of the search that holds the query: the report's first byte, or in the block
  search of many queries the first byte of the query's block, after the reports of the blocks
  before (as outfmt 0's refusal of such a record's masks). NCBI's output before the abort is what
  its stream had flushed, which depends on its buffer (to a file or a pipe, the version line of
  the prolog), so no prefix of it is NCBI's. BLASTN writes its prolog before the search, as NCBI
  does, so its refusal follows the prolog. The other sections of such a record are NCBI's
  (fixture rows `db.of1.bp.self1.*`, `db.of1.bn.selfrc.*`).

- RPSBLAST `-use_sw_tback` with composition-based statistics (version 1.30, 2026-11-15, line PS
  item PS11c-3, Owner-delegated, line PS `DECISIONS.md` 63): a profile with two or more
  preliminary HSPs when the e-value of a score of 0 passes `-evalue` (the query's search space
  times the profile's K below it: a large `-evalue` or a small `-searchsp`) and the profile's first
  position, after the composition adjustment, scores the query's first letter (after SEG) above 0.
  `Blast_RedoOneMatchSmithWaterman`'s loop (`redo_alignment.c:1438-1529`) then makes every
  positive alignment and comes to a round of score 0 whose alignment is that one cell, at (0, 0)
  (`smith_waterman.c:276-341, 437-446`, `blast_kappa.c:844-879`); its range changes no score the
  next scan sees, so every round is the same, and the loop leaves only when an allocation fails
  (after memory runs out, or by the OOM killer): no checkable result. LOSAT refuses (exit 1, "a
  profile whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI BLAST+ 2.17.0 ... is not
  supported by LOSAT's RPSBLAST") after the profile's first Smith-Waterman scan (the composition
  adjustment and the significance test of that scan's score come first), where its other refusals
  of a profile in the RPS traceback come. Where that cell scores 0 or less, NCBI's redo fails at the round of score 0
  and drops the profile, which LOSAT reproduces (rule 2; fixture rows `db.rps.rej.swendless` and
  `db.rps.mini.sw.zero.6`).

- psiblast and tblastn `-in_pssm` (version 1.32, 2026-11-15, line RG item D3-1, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07): LOSAT reads the PSSM's ASN.1 text as
  NCBI's reader does (every member of every type it reaches) and indexes the query's annotations
  as NCBI's object manager does when the PSSM's Bioseq joins the scope
  (`scope->AddTopLevelSeqEntry(pssm->SetQuery())`, `psiblast_app.cpp:537`, `tblastn_app.cpp:226`;
  `CSeq_annot_Info::x_InitAnnotList`). It refuses (exit 1, "a PSSM query whose annot has ... is
  not supported by LOSAT's PSIBLAST", or TBLASTN) where NCBI's output holds memory addresses,
  where NCBI reads outside its data or overflows a signed integer, and where its outcome depends
  on the machine's memory:
  - an empty location of a feature, graph, alignment, location set or Seq-table before NCBI's
    diagnostic stream fails: `Empty region in 0x...` (1204.1 to 1204.4, 1204.7, 1204.8,
    `seq_annot_info.cpp:761, 825, 898`) prints a heap address;
  - a Seq-table row without a Seq-id or a location: NCBI's `CCoreException` NULL-pointer report
    holds a stack of addresses;
  - two Seq-ids of one Seq-table row location that get different messages and that no gi,
    accession or general id number orders: NCBI writes them in the order of its `CSeq_id_Handle`s,
    which is then the address of their info (`seq_id_handle.hpp:207`);
  - a Packed-seg with fewer `present`, `starts`, `ids` or `lens` than its `dim` and `numseg`
    (NCBI reads past its vectors);
  - a signed overflow that NCBI's build runs without a sign of it (undefined behaviour): the sums
    of a Seq-table `int-delta` column, an `int-scaled` value past `Int8`, a column `field-id` below
    -2147483638 (its difference with the column group's base);
  - a Dense-seg of no segments with `starts` or `strands`: NCBI divides by `numseg` after its
    warnings (SIGFPE, exit status 136);
  - a feature Seq-table of more than 10,000,000 rows, a Std-seg `dim` above 1,000,000: NCBI's
    allocations succeed or fail by the machine's memory (10,000,000 rows take NCBI 1.4 GB,
    100,000,000 rows end in `BLAST ran out of memory`);
  - two `giim` Seq-ids of one number in one location after an earlier annotation, a Seq-table
    column or a Std-seg looked that number up, or with `-subject`: NCBI's `giim` tree
    (`seq_id_tree.cpp:2926-2964`) takes them for one Seq-id or two by the `giim` handles alive at
    that moment, which LOSAT follows only within one annotation (a conservative refusal; the port
    of those handles' lives stays open);
  - BitMagic data of a Seq-table that NCBI decodes (`bit-bvector` values, a `bit-set-bvector`'s
    sparse index): NCBI decodes well-formed streams deterministically and reads past or crashes on
    some malformed ones; refused until the Owner answers question Q-BM of D3-1 (port the decode or
    keep refusing). A `bit-bvector`'s size, which NCBI reads without a decode, is NCBI's.

  The `-in_pssm` file's format guess (NCBI's `CFormatGuess`, ported) is refused where one of its
  tests reads memory that the file does not decide: `TestFormatGvf`'s ninth column of a GVF line
  of eight (`tokens[8]`), `TestFormatPhrapAce`'s tokens of an `AS` line of fewer than three
  (`values[1]`, `values[2]`; NCBI 2.17.0 crashes on both, exit status 139) and `TestFormatFasta`'s
  search for `|` past a test buffer of 10 bytes or fewer.

  None of these is a rule-1 exception: D3-1's questions to the Owner, Q1 (the Dense-seg division by
  zero), Q2 (the memory bounds), Q3 (the signed overflows) and Q-BM, are open, and their cases stay
  refused until the Owner answers. Not rule-3 items but refusals until LOSAT writes NCBI's ASN.1
  text of the values it reads (line RG item D3-2): the annotation's text that NCBI prints after its
  `Failed to parse location of` message (1204.6), and with `-out_pssm`, `-export_search_strategy`
  or psiblast `-outfmt 11` a PSSM with members that LOSAT's PSSM writer does not write yet; nor
  the XML and JSON reports of an untitled query whose annotations hold features (NCBI's
  definition line reads them; not ported). The rest is NCBI's (rules 2 and 4): the index's
  messages (`Invalid 'lens' size in denseg`, the Seq-table column messages, `Failed to parse
  location of`, after which NCBI's later messages print their headers only) and exceptions
  (`Error: NCBI C++ Exception:`, `Error: vector::reserve` of a Dense-seg or Seq-table reserve), the
  reader's format errors, number-range and member-tag reports, its `#` for a byte outside the
  visible ASCII range and the deprecated members of a `Variation-ref`; where the index's
  exception ends the run, each output format's bytes (tblastn's `-outfmt 11` error archive with
  every message saved whole with its severity, psiblast's `-outfmt 10` and `20`), as with the
  reader's errors (fixture rows `db.ps.inpssm.d3.*`, `db.tbn.inpssm.d3.*`, `db.ps.inpssm.fg.*`,
  `db.tbn.inpssm.fg.*`).

## Explicit rejections of deterministic results (rule 2, version 1.10)

On a database with parsed Seq-ids NCBI gives these results without failing, and they look wrong
or need NCBI's object manager to reproduce. LOSAT refuses them instead (line DBID
`DECISIONS.md` 9, 20 and 21, line DI-7 for the query's id), under the maintainer's standing
delegation, which the maintainer accepted on 2026-11-02 with version 1.10. Each is a row of the table
in [`PD-LOSAT-DATABASE-SEARCH`](PD-LOSAT-DATABASE-SEARCH.md), which also lists refusals that rule 4
would port:

- A reported sequence whose id NCBI's formatters do not find, in the formats and field lists
  where NCBI exits 0: it writes `Unknown` or the id alone.
- A reported sequence whose id NCBI's lookup finds as another sequence first: NCBI writes that
  sequence's names, in every program but BLASTN also its identities (for the same record with
  other tax ids, its tax ids).
- Following hits whose ids NCBI takes for one sequence, other than copies of one record (which
  LOSAT shows in one hit as NCBI does): NCBI shows them as one hit, writes their tabular `qcovs`
  and `qcovus` as one subject's, and PSIBLAST's cut to `-max_target_seqs` counts them as one
  subject.
- A reported sequence with the id, or the display id of its description line (outfmt 0 to 4 of
  BLASTP, PSIBLAST and RPSBLAST, 1 to 4 of BLASTN and TBLASTN since version 1.28), of a query read
  with `-parse_deflines` in the query batch that reports it (with outfmt 5, in any batch), and in
  outfmt 1 to 4 of BLASTP, BLASTN and TBLASTN with the id or display id `lcl|Query_N` of a query
  read without it (TBLASTN's merge then reads the protein query's letters as nucleotides, probe
  `tn.t3b.d.qid.f1`): NCBI shows the query's title and letters for the
  id (its second lookup, by the display id, the query's title in that line).
- A stored title or Seq-id with a byte outside ' ' to '~': NCBI prints `#` and an error message
  each time it reads the definition line.
- Titles NCBI takes from the tax ids, where a report shows them (outfmt 0, a printed tabular id
  that replaces `lcl|Subject_N`): a title whose end `x_AdjustProteinTitleSuffix` rewrites (an
  empty title LOSAT makes from the tax id's name as NCBI does since line RG item T5).

## Ported (rule 4)

- TBLASTX `-window_size 0` with a low `-threshold` (version 1.21, 2026-11-14, line RG item X2b):
  the one-hit word finder's ungapped HSPs of length 0 or less (`aa_ungapped.c:1054,1083`). NCBI
  keeps their `Int4` ends, sorts them by score, links them with the other HSPs in the first
  `BLAST_LinkHsps` (where one whose score is above the large-gap cutoff, or with small gaps on
  above the small-gap cutoff, can join a chain and so change the chain's E-value and the list's
  order) and deletes them in its re-evaluation (a
  length of 0 or less scores 0, below every cutoff), so they are never reported. The X2 build
  refused these searches; LOSAT now gives NCBI's output (the rule-3 merge case above aside).
- `-evalue` with a signed infinity or NaN (`+inf`, `-nan`, `+nan(1)`, `1e999`): NCBI's
  `CArg_Double` reads them with `strtod` and its option check (`<= 0`) lets them pass;
  NCBI searches as with the largest e-value (the cutoff stays 1, no HSP is reaped). LOSAT
  rejected them before Session S07+++b.
