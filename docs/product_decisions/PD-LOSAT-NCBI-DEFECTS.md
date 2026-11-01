# Product Decision: NCBI BLAST+ behaviour that is a defect

- Decision ID: `PD-LOSAT-NCBI-DEFECTS`
- Version: 1.9
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
  not hold, rule 3, line DI-8f)
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
   `LOSAT/tests/db_search_manifest.tsv`, BLASTP and TBLASTN) and fewer than the whole
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

## Reproduced as NCBI (rule 2)

- A text seqid list whose ids all go in NCBI's rewriting for a version 5 database (gi ids
  only): NCBI prints `Empty seqid list` and searches the whole database (its code would
  keep no OID; the run keeps all), for a positive and for a negative list (line DI-8
  `DECISIONS.md` 18).

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
  list): NCBI's filtered header leaves no line, `GetBioseq` returns null
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
  the archive) its BlastDb data loader looks each hit's Seq-id up in every LMDB file of the
  list (`BlastFormatter_PreFetchSequenceData`, `CSeqDBLMDBSet::AccessionToOids`); a hit with
  a parsed Seq-id then stops the run with `BLAST Database error: DB contains no accession
  info.` (exit 2), after the outfmt 0 header and the earlier batches. Without such a hit, and
  for outfmt 11, NCBI runs to the end; a taxid list stops it up front. LOSAT refuses the list
  before the search (`the databases mix version 5 files with the accession table of
  makeblastdb -parse_seqids ... and without it`; version 1.9, line DI-8 `DECISIONS.md` 31,
  probe DI-8f F42, P14, H09, H15-H19): whether NCBI fails depends on the hits, and
  reproducing that would need the lookup before every batch's reports in every program.

- A negative seqid list that excludes a sequence (every id of its `.nos`/`.pos` record
  listed) of a volume that the database does not hold, while the database holds other
  volumes of the same LMDB file (`-db m.00` of a database in several volumes, or `-db "d1/m
  d2/m"` whose volumes interleave): `x_AdjustOidsOffset` turns that sequence's OID into
  `kSeqDBEntryNotFound`, and NCBI clears bit -1 of its OID mask and crashes (SIGSEGV, exit
  139; review of line DI-8f, fixture row `db.mdb.bn.negpartial.6`). LOSAT dropped the OID
  and searched; it now refuses the list (`the negative seqid list excludes sequence N of the
  LMDB file, in a volume that the database does not hold`; version 1.9). A positive list on
  such a database is NCBI's and LOSAT's alike.

## Ported (rule 4)

- `-evalue` with a signed infinity or NaN (`+inf`, `-nan`, `+nan(1)`, `1e999`): NCBI's
  `CArg_Double` reads them with `strtod` and its option check (`<= 0`) lets them pass;
  NCBI searches as with the largest e-value (the cutoff stays 1, no HSP is reaped). LOSAT
  rejected them before Session S07+++b.
