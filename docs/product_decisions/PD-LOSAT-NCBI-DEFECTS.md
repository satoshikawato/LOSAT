# Product Decision: NCBI BLAST+ behaviour that is a defect

- Decision ID: `PD-LOSAT-NCBI-DEFECTS`
- Version: 1.51
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
  1.32 2026-11-15 (the sequence that a positive seqid list keeps while matching none of its lines,
  of 1.6 and 1.29, is no longer an explicit rejection: NCBI's result is deterministic, and LOSAT
  reproduces it, rule 2; a negative list's stays refused; line DBID item lineless, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07);
  1.33 2026-11-15 (BLASTX's ungapped one-hit HSPs of negative length and of length 0 ported, rule 4,
  line RG item X2g, with the refusals of decisions X2g-4 and X2g-5 of that line and the tables of
  X2g-5's one alignment ported by X2g-7: Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07);
  1.34 2026-11-16 (the sequence that a negative seqid list keeps while removing all its lines, of
  1.6, 1.29 and 1.32, is no longer an explicit rejection: NCBI's fetches exclude it, and LOSAT
  reproduces them, rule 2; line DBID item lineless2, Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07);
  1.35 2026-11-16 (exception 13, BLASTN `-db_hard_mask` on a sequence whose first mask runs to its
  end, at a scanned word of 11 or more; the Owner's decision of 2026-11-15 「M4a Q1=B」, line NP
  item M4a-2c);
  1.36 2026-11-16 (the refusal of 1.31 of a stopping batch whose other sequences' loads write texts
  of the taxonomy database narrowed in the JSON and XML2 formats, outfmt 13 to 16, to the texts
  after the lookups of the batches before, as in the other formats: a stop after a batch that
  opened the taxonomy database is written as NCBI writes it; line DBID item mixtax,
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07);
  1.37 2026-11-16 (a subject whose title, or for an empty database title its name in the taxonomy
  database, NCBI's `GenerateDefline` cannot decode: the tabular rows and XML values of NCBI's
  catches reproduced, rule 2, and the rows at which NCBI stops or prints an uninitialised taxid
  refused, rule 3; line RG item T7, Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07);
  1.38 2026-11-16 (the end of a database protein's title that `x_AdjustProteinTitleSuffix`
  rewrites from the tax ids, which version 1.10 listed among the explicit rejections of
  deterministic results of rule 2, moved to the results ported under rule 4, line RG item T5c;
  a name written there that `GenerateDefline` cannot decode fails the lookup as in 1.37;
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; written 2026-11-14 as
  1.26);
  1.39 2026-11-16 (the query-anchored views of BLASTX and TBLASTX, line RG item T3c: TBLASTX's
  strand drop, the spread of a subject's masks and TBLASTX's merge of a query record that is also a
  subject into the query's sequence reproduced, with the two exceptions that this merge can throw
  (the report up to the query's description table, then NCBI's error, exit 255; 2026-11-16, T3c
  fix 2, in place of fix 1's explicit rejection), rule 2; and the rule-2 rejection of a query's id
  of 1.24 and 1.28 in BLASTX and TBLASTX; BLASTX's HSPs of length 0 of 1.33 in these views
  reproduced, rule 2, with the explicit rejection of a section whose merge NCBI frees a segment in,
  rule 3 (2026-11-16: the split of T3c after X2g, and `FillUnalignedRegions`'s free of T3c fix 3);
  Owner-delegated, standing instruction 2026-09-29,
  restated 2026-10-07, within rule 2, rule 3 as 1.23 and the class of 1.10);
  1.40 2026-11-17 (TBLASTX's ungapped extension past the end of a query frame into the next frame
  of the query batch: reproduced where it scores no higher there, rule 2, and an explicit
  rejection where it does, rule 3, the Owner's decision of 2026-11-16 on line RG's question
  Q-X2c; line RG item X2c, the reproduction Owner-delegated, standing instruction 2026-09-29,
  restated 2026-10-07; the TBLASTX one-hit merge refusal of 1.21 now also comes after the reports
  of the batches before);
  1.41 2026-11-17 (that refusal narrowed to the loads whose order NCBI takes from the addresses
  of its Seq-id handles: the texts of the loads before the stopping lookup in the order of the
  handles are written as NCBI writes them, and of several stopping ids the first that NCBI looks
  up ends the loads; refused, a stopping id that NCBI may never look up, as its scope or its data
  loader, which keeps the ids of what the run loaded, finds it first; a directory in place of a
  version 4 ISAM file is no ISAM, as NCBI's; line DBID item loadord, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07);
  1.42 2026-11-17 (a version 4 lookup returns at the first volume that holds the id, as NCBI's, so
  an ISAM file of a later volume that does not open stops none of its lookups; refused, where
  LOSAT's list of the other holders fails on such a file; line DBID item v4late, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07);
  1.43 2026-11-17 (the explicit rejections of a taxonomy record read past NCBI's mapping of
  `taxdb.btd`, rule 3, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07, and
  exception 14, LOSAT's output where NCBI may fault there, which the Owner decided on 2026-11-15;
  line RG item T5b; written 2026-11-15 as 1.26 with exception 13);
  1.44 2026-11-17 (TBLASTN `-in_pssm -window_size 0`, a one-hit gapped start on its frame's NULLB
  or after it, which a PSSM whose '-' column scores above BLAST_SCORE_MIN allows: its preliminary
  gapped score reproduced, rule 2, and the explicit rejection of its traceback after the NULLB,
  rule 3; line RG item X2d, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07);
  1.45 2026-11-17 (psiblast and tblastn `-in_pssm`, line RG item D3-1: exception 15, a PSSM whose
  values nest deeper than NCBI's stack holds, read to any depth, rule 1, decided by the Owner on
  2026-11-17 (Q-DEPTH = B); the explicit rejections where NCBI's object
  manager or format guess prints addresses, reads outside its data, overflows a signed integer or
  depends on the machine's memory, rule 3, and of well-formed BitMagic data and of a `-subject`
  duplicate-id message that may spell a subject's Seq-id, rule 2, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; the Owner decided D3-1's questions Q1 to Q3 and
  Q-BM on 2026-11-16 (「推奨通り」: Q1 to Q3 refused, Q-BM to port BitMagic); written 2026-11-15,
  numbered after 1.44);
  1.46 2026-11-17 (TBLASTN, the explicit rejection of the traceback of an alignment on a subject of
  fewer than 12 bases, whose frame NCBI never translates and reads from address 1 on, rule 3; line RG
  item X-short, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07);
  1.47 2026-11-17 (outfmt 13 to 16 of a reported sequence whose lookup in a version 5 list fails
  where the batch's prefetch is off are written as NCBI writes them: 15 and 16 post the exception
  and leave the query's report out, 13 and 14 stop with `BLAST Database error`; the taxonomy
  database's texts at the `Hit`s that NCBI loads, none of the throwing `Hit` or of those after it
  in that report; the refusal of 1.31 narrowed to the other formats; line DBID item pfoff, review 1
  fix included, Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07);
  1.48 2026-11-17 (approved exception 3 extended from BLASTP to TBLASTN and BLASTX: the window of a
  one-hit gapped start reads the letters past NCBI's buffers as the sentinel, and the explicit
  rejections of 1.18 and 1.22 there are removed, with a narrower rule-3 rejection of a TBLASTN
  `-in_pssm` gapped start past the NULLB after frame -3; decided by the Owner on 2026-11-17, line RG
  DECISIONS.md OQ-1117, J2 = A, words 「A 広げる (推奨)」; line RG item J2A; numbered after 1.47);
  1.49 2026-11-17 (psiblast and tblastn `-in_pssm`, line RG item D3-2: exception 15 also for NCBI's
  ASN.1 text writer, which crashes writing back with `-out_pssm` the nests it read; the refusals of
  the annotation's text after `Failed to parse location of` (1204.6) and of the PSSM members that
  LOSAT did not write removed; the rule-3 rejection of Seq-ids in the order of addresses widened to
  a feature's and a Seq-table's location; the PSSM's most negative Int4 and Int8 values written as
  their digits; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered
  after 1.48);
  1.50 2026-11-18 (outfmt 6, 7, 10 and 20 of such a sequence are written as NCBI writes them: a
  field list of the subject block writes the row that its `catch` leaves and goes on, one that reads
  the sequence after that block stops with `BLAST Database error`, `nident` alone only at a row with
  identities; the organism report (outfmt 18) leaves such a sequence's HSPs out and goes on; the
  refusal of 1.31 narrowed to the other formats and to the first row of a query's report with
  `sacc`, `saccver` or `staxid`, where NCBI stops or prints an uninitialised number as in 1.37;
  line DBID item pfoff-2, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07;
  numbered after 1.49);
  1.51 2026-11-18 (outfmt 0 to 4 of such a sequence are written as NCBI writes them: its
  description line is the Seq-align id and `Unknown`; outfmt 0 stops with `BLAST Database error`
  where the first alignment shown is of such a sequence under a Seq-align id that is not local,
  else writes the skipped line for each of its HSPs; outfmt 1 to 4 stop where any alignment shown
  is; the description lines load the sequences that they show, in psiblast's rounds after the
  first those of the two-part table of the first `-num_descriptions` HSPs; the refusal of 1.31
  no longer covers outfmt 0 to 4; line DBID item pfoff-3, Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.50)
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

3. **Approved exception (BLASTP, TBLASTN and BLASTX): a one-hit gapped start that reads past a
   sequence.**
   With `-window_size 0` (one-hit) and a low `-threshold`, NCBI's ungapped extension can
   give a negative `Int4` length that `BlastGetStartForGappedAlignment` receives as
   `Uint4` (`aa_ungapped.c:1054,1083`, `blast_gapalign.c:3393-3437`): it scores only the
   first window of 11 letters, and when that window passes the end of the query or the
   subject it reads past the sequence buffer and indexes the matrix with the stray byte;
   some such searches crash (SIGSEGV, exit 139; 21 of 5,000 random short searches with
   `-threshold` 1 to 3 in the round-2 audit). LOSAT reads the letters past the sequence
   as the sentinel (`NULLB`, `BLAST_SCORE_MIN`), as for the letter right after it. NCBI run
   under `valgrind -q`, which reports the invalid reads and did not crash in those searches, gives LOSAT's
   output byte for byte in all 23 crashing cases of the audit's second and third passes
   (`docs/evidence/losat_web_e2e/audit/round2/a2_blastp.md`, `a3_blastp.md`). Searches
   whose window stays inside the sequences and its sentinel match NCBI (fixture
   `e2e.blastp.one_hit_negative_width`).

   TBLASTN and BLASTX (version 1.48, decided by the Owner on 2026-11-17, line RG item J2A): their
   one-hit word finders give the same negative lengths (`aa_ungapped.c:1054,1083`), and the same
   window reads past NCBI's buffers: past the byte after the query batch's last NULLB
   (`calloc(buflen+1)`, `blast_setup_cxx.cpp:500-501`); for TBLASTN past the NULLB after frame -3
   of the subject's translation, the last byte NCBI's buffer sets (`2*(nucl_length+1)+2` bytes,
   `blast_util.c:1058,1063,1098-1101`); for BLASTX past a `-subject` protein's `malloc(len+2)`
   (`blast_setup_cxx.cpp:885-893`) and past the last memory page of a database volume's sequence
   file (`seqdbatlas.cpp:114-127`). LOSAT reads every letter past those buffers as the sentinel, as
   for BLASTP. A window that passes its context or subject already holds that NULLB, which scores
   `BLAST_SCORE_MIN` with every matrix (`blast_stat.c:1559-1574`), so inside NCBI's buffers the
   start is the HSP's whatever follows: BLASTX and TBLASTN's query side read NULLBs from there on,
   as BLASTP does. TBLASTN's subject side reads the bytes NCBI reads (the next frames' letters,
   which a hand-made PSSM's '-' column can score: the TBLASTN item of version 1.18 below). In 450
   random TBLASTN searches with this setting (300 with `-subject`, 150 with `-db`; one to three
   queries of 12 to 120 letters, one to three subjects of 36 to 360 bases, word sizes 2 to 5,
   `-threshold` 1 to 4), 129 read past a buffer (106 first past the query batch, 23 past the
   subject's translation); NCBI crashed (SIGSEGV) in 23 of them, in each of three runs, and no
   output varied between runs. NCBI under valgrind, which reports the invalid reads in
   `BlastGetStartForGappedAlignment`, finishes all 129 with LOSAT's output byte for byte; the 106
   that NCBI finishes natively and the 321 searches inside the buffers equal NCBI's native output.
   For BLASTX, in the 400 random searches of the BLASTX item of version 1.22 below, LOSAT equals
   NCBI under valgrind in the 115 that read past a buffer (6 of them crash natively) and NCBI's
   native output in the other 285; natively NCBI's report differs from LOSAT's in 2 of the 109 past
   ones it finishes. LOSAT's output is the same at 1 and 4 threads in all 850. Fixtures
   `db.x2.tblastn.qpast.*`, `db.x2.tblastn.frame.*`, `db.x2.tblastn.pssmpast.6` (`-in_pssm`: the
   window reads the byte after the NULLB after frame -3, which NCBI's buffer never sets) and
   `db.x2f.bx.*` (`db.x2f.bx.page.6` is an ncbi-near row: it takes the output of NCBI's native run on the same
   volume with 16 zero bytes after its last page, `db.x2f.bx.pagepad.6`, because valgrind itself faults on it) of
   `LOSAT/tests/db_search_manifest/` (`-subject` and `-db`, the query and the
   subject side), and BLASTP's `db.j2a.bp.past.6` (the round-2 audit's case R2A2-2, on which NCBI
   crashed natively in the audit), frozen from NCBI under valgrind (except `db.x2f.bx.page.6`) where NCBI reads past its
   buffers. Not covered, and still refused: with TBLASTN `-in_pssm`, a window past the PSSM's last
   position (`pssm->data` holds one row per position and no sentinel row: the read is past a
   matrix, not past a sequence), and a gapped start that a PSSM whose '-' column scores above
   `BLAST_SCORE_MIN` puts past the NULLB after frame -3, the last byte the translation buffer sets,
   where NCBI's preliminary gapped alignment then reads the byte after it, which the buffer
   allocates and never sets, or past it (the TBLASTN item of version 1.18 below; fixture
   `db.x2.tblastn.pssmpast.rej`).

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

## Approved exception (BLASTN `-db_hard_mask`, a mask that runs to the end, version 1.35)

Version 1.35. The Owner decided on 2026-11-15 (line NP `DECISIONS.md` 45, words 「M4a Q1=B」) question
Q1 of M4a-2b (PR #1173): option B, the search over the range before the mask at a scanned word of 11
or more, blastn-short keeping the refusal. Not taken: option A, refusing the shape for every task
(it refuses searches that NCBI ends normally, such as a transcript database whose poly-A tails are
masked). Line NP item M4a-2c.

13. **`-db_hard_mask`: a chunk from the start of a sequence whose mask runs to its end.** With
    `-db_hard_mask` NCBI cuts each database sequence at its masks into chunks
    (`s_GetNextSubjectChunk`, `blast_engine.c:221-276`). A sequence whose first mask with letters runs
    to its end, with at most 5,000,000 letters before it and no other mask before it than masks of no
    letters at 0 (masked whole, or a masked tail such as a poly-A tail), gets one chunk from 0 that ends
    at the mask, and for that chunk the test of the "no chunking" branch (`offset == 0 && residual == 0
    && next == full_range.right`, `blast_engine.c:264-268`) holds: its `seq_ranges` are the soft
    ranges, which `-db_hard_mask` leaves as the full range (`blast_engine.c:162`), while
    `subject->length` is the chunk's: NCBI uses two lengths for one subject. The preliminary stage
    reads the subject's stored letters in the mapped volume file, unmasked (`seqsrc_seqdb.cpp:326-330,
    375`, `seqdbvol.cpp:1817`: masks reach the engine only as ranges), and the scan
    (`na_ungapped.c:1647-1668`) and the word checks (`s_BlastNaExtend`, `s_TypeOfWord`, bounded by the
    scan's end, `na_ungapped.c:1123, 508-554, 1683`) read them from 0 to the full length: the letters
    under the mask are read as they are stored, nothing outside the sequence (valgrind finds no error in
    such runs). The extensions take the chunk's length (`na_ungapped.c:165, 213-214`,
    `blast_gapalign.c:2784, 4021`), and so does the interval tree of the gapped HSPs, `[0, a + 1]`
    (`blast_gapalign.c:3830-3831`). Three effects follow:
    - a word that crosses the mask's start, or lies inside the mask, is found and extended up to the
      mask's start: NCBI reports a row that its normal chunks of the near input (the mask one letter
      shorter) do not find, where the alignment has fewer letters before the mask than a word;
    - the diagonal table keeps those words' ends, up to the full length past its offset
      (`na_ungapped.c:757, 768`), while `Blast_ExtendWordExit` moves the offset on by the chunk's
      length and the window only (`blast_extend.c:172`, called at `na_ungapped.c:1686`), so the
      subject that the same thread searches next loses its seeds on those diagonals (`na_ungapped.c:672`).
      Each thread has its own table (`blast_engine.c:1002, 1377`) and takes the next database
      sequences as it asks for them (`seqdbimpl.cpp:301-395`, `prelim_stage.cpp:142-188`), so at two
      threads and more which subject loses its hits changes from run to run;
    - two preliminary HSPs past the tree's end with the same query range are never separated: the
      tree's node `[l, r]` with `l > r` keeps its middle at `r` (`blast_itree.c:95, 107`), each insertion
      splits again and doubles the node pool (`blast_itree.c:70-79`) until the allocation fails; the
      search then stops (`blast_gapalign.c:4087-4090`, `blast_engine.c:1469-1478`) and frees the
      subject's block (`blast_engine.c:1595`), whose sequence pointer, marked allocated although it
      points into the mapped volume (`seqsrc_seqdb.cpp:392`, `blast_util.c:225-227`), is passed to
      `free()`: NCBI BLAST+ 2.17.0 ends with `free(): invalid pointer`, `munmap_chunk(): invalid
      pointer` (SIGABRT) or SIGSEGV, at blastn-short's word 7 and at words of 11 and more alike.
    At a scanned word (the word size, or the discontiguous template's length)
    of 11 letters or more (megablast, blastn, dc-megablast) LOSAT scans that chunk as NCBI scans every
    other hard chunk, from its start to its end: the range before the mask, which NCBI's normal chunks
    also search in the near input whose mask is one letter shorter. Its extensions and its diagonal
    state take that one length, so the subjects searched after it are searched as without it, at every
    thread count, and the search ends normally (exit 0). The traceback reads every masked
    letter as `B`, as NCBI's traceback does (`s_SeqDBMaskSequence`, `seqdbvol.cpp:1518-1585`). Below
    a scanned word of 11 (blastn-short's 7) LOSAT refuses such a database before it writes anything,
    as it did for every word before version 1.35 (rule 3). `blastinput::db_mask::
    TAIL_MASK_SEARCHED_WORD`, `algorithm::blastn::blast_engine::run::next_hard_chunk`.
    Evidence (NCBI BLAST+ 2.17.0, each run capped; inputs, scripts and results in the line's task
    folder, `M/m4a2b/q1` and `M/m4a2c`):
    - The question's corpus: 300 random subjects of 300 to 2,000 letters, each with one tail mask of
      12 letters to two thirds of the subject (200 of random letters, 100 with a masked poly-A tail),
      and a query across the mask's start that runs to the subject's end. At words of 11 and more
      NCBI's output of the tail mask, its output of the near input and LOSAT's are the same bytes in
      all 900 runs. At word 7, NCBI's output of the tail mask differs from the near input's in 29 of
      300 runs and NCBI crashes in 71.
    - Tails of 1 to 4 letters (200 subjects, queries that match to the subject's end): LOSAT's
      output is NCBI's own in all 600 runs at words of 11 and more. The near input's alignment reaches
      the subject's last letter, which is not masked there, in 50, 97 and 97 of the 200 runs
      (megablast, blastn, dc-megablast); NCBI's run of the tail mask and LOSAT's read it as `B`.
    - Alignments with few letters before the mask (200 subjects of 30 to 150 letters whose chunk
      before the mask holds 0 to all but 12 of their letters, queries over the whole subject; and 200
      subjects of 300 to 2,000 letters whose queries hold 4 to 60 letters before the mask and 60 in
      it): NCBI's output of the tail mask equals LOSAT's, the near input's, in 115, 177 and 116, and
      in 113, 178 and 105 of the 200 runs (megablast, blastn, dc-megablast). In the others NCBI adds a
      row from a word across the mask's start, checked against the masked letters and cut at the
      mask's start, often shorter than the word (12 letters at megablast's 28), which its normal
      chunks of the near input do not find; LOSAT gives the near input's rows there.
    - The row follows the stored letters under the mask: with those letters complemented, NCBI's
      extra row goes away (the O1 survey of M4a-2c, `M/m4a2c/O1_SURVEY.md` §3.1).
    - The next subject: a tail-masked subject X (mask after 20 letters) followed by an unmasked Y that
      holds the query's first 300 letters on the same diagonal (fixture `edge9`): NCBI drops Y's
      300-letter 100 % row (E-value 5.6e-162) in 20 of 20 runs at `-num_threads 1`, 3 of 20 at 2 and 5
      of 20 at 4 threads (two different outputs at 2 and at 4); with the near mask Y's row is in all 40
      runs, and so it is with Y alone, Y before X, or Y's copy on the next diagonal. LOSAT prints Y's row,
      the near input's output, in 180 of 180 runs (megablast, blastn, dc-megablast at 1, 2 and 4
      threads, 20 runs each).
    - A crash at words of 11 and more: two copies of the query inside a tail mask (fixture `edge10`):
      NCBI SIGSEGV (exit 139) in megablast, blastn and dc-megablast; one copy, or the near mask, exit 0
      and no row. LOSAT: exit 0 and the near input's output (no row) at 1 and 4 threads. Under
      valgrind the same run shows the node pool doubling, "BLAST ran out of memory" and an invalid
      `free()` of an address in the mapped `.nsq`, the mechanism of blastn-short's 71 crashes too.
    - The fixture databases: NCBI crashes in each of 10 runs of edge2's algorithm 16 (s4 masked
      whole) and edge4's 11 (the mask after an empty mask at 0) in megablast, blastn and
      dc-megablast, and of `mask_nucl_mv`'s 11 (m04's tail) in megablast and dc-megablast; LOSAT
      writes NCBI's bytes of the near input there.
    Fixture rows `db.dbmask.hard.tail.*` and `db.dbmask.hard.whole.*` of
    `LOSAT/tests/db_search_manifest/blastn.tsv`: NCBI's bytes, or with `ncbi-near:` those of NCBI's
    run of the near input (the rows `db.dbmask.hard.near.*`; `LOSAT/tests/db_search_fixtures.py`), `edge9`'s and
    `edge10`'s (`db.dbmask.hard.tail.edge9.6`, `db.dbmask.hard.tail.edge10.6`) among them;
    `db.dbmask.hard.rej.tail.*.short.6` and `db.dbmask.hard.rej.whole.*.short.6`, the refusal at word 7.
    Unit tests `the_hard_mask_shapes_follow_ncbi_s_chunks` and
    `a_hard_masked_search_refuses_the_masks_that_ncbi_s_chunks_do_not_handle`
    (`LOSAT/src/blastinput/db_mask_tests.rs`), `hard_masked_chunks_follow_ncbi_s_get_next_subject_chunk`
    and `hard_mask_refusals_are_the_merges_that_ncbi_decides_past_its_list`
    (`LOSAT/src/algorithm/blastn/blast_engine/run.rs`).

## Approved exception (psiblast and tblastn `-in_pssm`, values nested deeper than NCBI's stack holds, version 1.45; its writer 1.49)

Version 1.45, line RG item D3-1. The Owner decided on 2026-11-17 (line RG `DECISIONS.md` OQ-1117,
words 「B 深さに関係なく読む (推奨)」) question Q-DEPTH of D3-1 (PR #1262): option B, reading every depth and
writing NCBI's output with a stack large enough. Not taken: option A, refusing the depths at which
NCBI may crash, from NCBI's stack bytes per nested value computed for the run (D3-1 fix round 6):
those bytes differ with the path of each recursion and with the types open together, so such a rule
refuses depths that NCBI reads at every start and reads depths at which NCBI crashes at every start.

15. **A PSSM whose values nest deeper than NCBI's stack holds.** NCBI's reader of the `-in_pssm`
    ASN.1 text (`CObjectIStreamAsn`) reads a value whose type holds values of its own type by
    recursion: `ReadObject` and `SkipObject` of each value (`include/serial/impl/objistr.inl:93-96`,
    `:99`), entered from each member and variant (`src/serial/member.cpp:710`,
    `src/serial/variant.cpp:471`), with the class, choice and container functions between them, take
    some hundreds to about 2,100 bytes of the main thread's stack for each value open. The PSSM's types nest in themselves through Bioseq-sets
    (`seq-set`), Seq-locs (`mix`, `equiv`), User-objects and User-fields (`data object`, `data
    objects`, `data fields`), Seq-aligns (`segs disc`), Pubs (`equiv`), Variation-refs (`data set
    { variations }`, `consequence { variation }`), SeqTable-multi-data (`int-scaled`, `int-delta`,
    `real-scaled`) and OrgNames (`name hybrid`). With the default stack limit of 8 MiB NCBI BLAST+
    2.17.0 overflows its stack from some thousands of levels and ends with SIGSEGV (exit status 139,
    no output): psiblast with the most stack it gets (a minimal environment, no random start) crashes
    from 2,976 User-objects nested through `data objects`, 5,082 Bioseq-sets, 6,888 Seq-loc mixes
    and 9,519 SeqTable-multi-data through `int-scaled`, for example. The depth is not the input's: the
    stack limit, the lengths of the command line and of the environment (which the kernel lays out
    at the stack's end) and a random start of 0 to 8,191 bytes at each run move it, and each
    recursion path and each combination of the types open takes its own bytes, so one PSSM is read
    in one run and crashes in the next. With a stack large enough (`ulimit -s unlimited`) NCBI reads
    every such PSSM and writes its normal output: the search, or for a query Seq-entry set its
    `CSeq_entry::GetSeq()` exception. LOSAT keeps the values it reads, and the object manager's walks
    over them, on the heap, so it reads every depth on every target as far as memory allows and
    writes NCBI's output with a stack large enough. NCBI's writer of the same values
    (`CObjectOStreamAsn`: `WriteObject`, include/serial/impl/objostr.inl:106-110, from each member
    and variant, src/serial/member.cpp:934-970, src/serial/variant.cpp:577-651) recurses the same
    way when psiblast writes the PSSM it read back with `-out_pssm` (or the search strategy's request
    holds it): with 8 MiB NCBI 2.17.0 read and then crashed writing the `-out_pssm` file of
    a descriptor's User-object nested 4,000 deep and a Seq-table column's data nested 9,600 deep
    through `int-delta` (exit status 139, no file, no output); LOSAT writes it from a list on the heap (line RG item D3-2,
    `asn_types::write`), as NCBI writes it with a stack large enough. No PSSM that psiblast writes
    from a FASTA query nests this way.
    Evidence: `docs/evidence/pssm_nest_depth/` (for each kind of nest, NCBI's run with 8 MiB, a
    crash, and with no stack limit, the same bytes as LOSAT's; the generator of the inputs). Fixture
    rows `db.ps.inpssm.d3.om.r6.nest_*`, `db.ps.inpssm.d3.om.r7.nest_*`, `db.tbn.inpssm.d3.om.r6.nest_*`
    and `db.tbn.inpssm.d3.om.r7.nest_*` of `LOSAT/tests/db_search_manifest/`, with their `.outpssm`
    twins (those marked `ncbi-stack: unlimited` frozen from NCBI's run with no stack limit,
    `LOSAT/tests/db_search_fixtures.py`). Unit tests `every_recursive_type_nests_to_any_depth`,
    `a_nest_of_any_depth_is_read_without_recursion` (`LOSAT/src/report/serial/asn_types/read_tests.rs`),
    `a_nest_of_any_depth_is_written_without_recursion` (`write_tests.rs`),
    `nests_of_any_depth_are_walked_without_recursion` and `nested_multi_data_are_read_without_recursion`
    (`LOSAT/src/algorithm/pssm_search/annot_index/tests.rs`, `table_tests.rs`).

## Reproduced as NCBI (rule 2)

- TBLASTX `-outfmt 1` to `4` (the query-anchored views; version 1.39, 2026-11-15, line RG item T3c,
  DECISIONS T3-6): NCBI adds TBLASTX's HSPs of each query frame with `CAlnMix::fForceTranslation`
  (`eTranslateNucToNucAlignment`, showalign.cpp:2058-2063), whose `CAlnMixMatches::Add` throws
  "Unable to mix strands when forcing translation!" for an HSP whose subject strand differs from
  the sign of the subject's strand score so far, the strand of its first HSP in that frame
  (`alnmatch.cpp:191-205`); `DisplaySeqalign` catches it and leaves the HSP out of the section
  (`showalign.cpp:2073-2076`), while the description table and outfmt 0 show it. The throw comes
  after `m_DsCnt++`, so the HSPs after it keep their numbers. LOSAT leaves the same HSPs out
  (`report::query_anchored::aln_mix`; probes `tx.s.f1`, `tx.t3c.s.f1`; fixture rows
  `db.of1.tx.s.f*`). BLASTX's HSPs of a minus query frame are added with `fNegativeStrand`, which
  `Add` reads as `fForceTranslation` (`alnmix.hpp:77,102`): their query rows run from high to low, and
  none is left out (a protein row has no strand).
- TBLASTX `-outfmt 1` to `4` of a query record that is also a `-subject` record under
  `-parse_deflines` (version 1.39, line RG item T3c fix 1): without `fPreserveRows`
  (`fForceTranslation`, showalign.cpp:2058-2063) `CAlnMixSequences::Add` makes no other row for the
  record's subject row (alnseq.cpp:173-196), so both rows of its HSPs with itself are the query's one
  `CAlnMixSeq`. No sequence then holds as many Dense-segs as the mix (alnmerger.cpp:108-127: no
  single refseq), `x_SetSeqFrame` of the second row compares it with the first (alnmerger.cpp:
  1272-1276) and takes the query's own row, and `x_SecondRowFits` gives that row the first such HSP
  and every later one a new extra row of the query's sequence. The first one's subject positions are
  the query's own starts on the diagonal (`this position already exists`, alnmerger.cpp:753-771), so
  the record's alignment with itself in each query frame disappears into the query row, and its
  other alignments with itself are rows labelled with the record above the other subjects' rows.
  NCBI writes that section with exit 0, from the inputs alone; LOSAT writes the same
  (`report::query_anchored::aln_mix`; the review of T3c's cases `x.self.*`, `t.self.*`,
  `v.self.b1.f1`, the probes `tx.t3c.self.*`, the 13 rows `pd.self.tblastx.*` of
  `tests/outfmt0_manifest.tsv` in formats 1 and 3; fixture rows `db.of1.tx.self.*`).
- The same TBLASTX sections where that merge throws (version 1.39, 2026-11-16, line RG item T3c
  fix 2, in place of fix 1's explicit rejection). The first HSP of the record with itself that the
  merge takes (the longest, ties in their order, an HSP of the other strand included) goes to the
  query's own row, and its walk can create a start whose segment already has that row
  (`CAlnMixMerger::x_Merge(): Internal error: Start iterator already exists for seq2.`,
  alnmerger.cpp:779-794, a repeat whose copies `-query_loc` and `-subject_loc` keep apart, a
  palindrome) or join segments of the query's row into a loop that `CAlnMixSegments::Build` meets
  (`CAlnMixSegments::Build(): Internal error: Infinite loop detected.`, alnsegments.cpp:162-169, a
  tandem repeat). Only `Add` is in `DisplaySeqalign`'s `try`, not `Merge` (showalign.cpp:2060-2083):
  the exception leaves the formatter before the section's first byte, and the application's
  `CATCH_ALL` reports `Error: NCBI C++ Exception` and exits with 255, its stream flushed. So NCBI's
  output is deterministic: the reports of the queries before (the batches before and the batch's
  queries before), the query's heading and description table, no epilog; standard error the
  warnings so far, then the exception (the same bytes to a file, a pipe and `-out`, with 400 queries
  before and with `BATCH_SIZE=500`: review 2 of T3c fix 1, `T/t3/review2_t3c/probes/buf/`). LOSAT
  runs the same merge and writes the same (`report::query_anchored::write_section` fails with
  NCBI's exception, `cli::NativeError::of_writer`; the review's 936 cases where NCBI throws, its
  108 `mg*` cases of both exceptions, the probes `T/t3/fix1_t3c/cases_fix1.tsv`; fixture rows
  `db.of1.tx.throw.*`).
- `-outfmt 1` to `4` of BLASTP, BLASTN and BLASTX (version 1.39, DECISIONS T3-6): a subject's masks
  (the `-lcase_masking` runs of a `-subject` record, BLASTN's `-db_soft_mask`) are placed on the columns of the subject's first row that
  meets each and shown in lower case on those columns of every row of the subject
  (`x_FillLocList`, showalign.cpp:4291-4338; `x_OutputSeq`, :2486-2523), also on its rows of other
  HSPs, where the letters are not masked (probes `bn.lc2.f2`, `bn.t3b.lc.*`, `bx.t3c.slc.*`).
- BLASTX `-outfmt 1` to `4` with the HSPs of length 0 of `-ungapped -window_size 0` (version 1.39,
  2026-11-16, line RG item T3c after X2g, 1.33): the description table and the sorts stop the report
  as in outfmt 0 (`Row is empty`, the Dense-seg's exception). Past the table, the section makes a
  Dense-seg of every aligned HSP before its first byte and outside its `try`
  (showalign.cpp:2002-2009): `-subject`'s mapped Std-seg throws there (`Seq_align.cpp:801`, exit
  255). The Dense-seg of empty rows as made (`-db`, and `-subject`'s that the mapper leaves) has one
  segment of length 0 and no widths (`Seq_align.cpp:774-934`), so `DisplaySeqalign` puts it in the
  group of the query's frame +1 whatever its frame (showalign.cpp:2010-2041) and adds it with
  `fPreserveRows`. Its match is merged last (score 0): at the query's own frame row of that group it
  is a segment of length 0 that shows nothing (a frame +1 HSP at query offset 0 starts the frame +1
  row there, and the unaligned letters fill up to the next HSP), at another nucleotide a query row
  of that frame (`x_SetSeqFrame`) without a subject row (`x_CreateDenseg` drops the empty one).
  The anchored views (1, 2) show no such row; the flat views (3, 4) translate it as a row after row
  0 with BLASTX's database genetic code 0 (`CBlastOptionsLocal::SetProgram`), which throws
  `Unable to find genetic code number 0` inside the group's `try`: the group shows its heading
  `For reading frame +1` alone. A section of such HSPs only is an alignment of no columns
  (`(int)GetAlnStop()` is -1): its blank line alone. LOSAT writes the same (probes
  `T/t3/x2g/`: 1,880 cases of the X2g fixtures' inputs and of the X2g-1 probe P1's pairs, 1,662 of
  P1's 277 pairs with such HSPs, 18 of two queries; fixture rows `db.x2g.bx.of.*`).

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

- A sequence that a positive seqid list (`-seqidlist`) keeps while matching none of its
  definition lines (a pdb molecule named without its chain, a local id in other letters, a
  Swiss-Prot name on a version 4 database, a list's `gnl|BL_ORD_ID|n` that the accession table
  resolves by the stored number to a sequence in a volume whose first OID is not 0, whose lines
  NCBI's header read tests with the volume's first OID added, `seqdbvol.cpp:2222-2236`, of 1.29):
  the sequence is searched and its length counts in the statistics, but NCBI's filtered header
  leaves it no Seq-id, and `BlastHitList2SeqAlign_OMF` makes no Seq-align of its hits
  (`blast_seqalign.cpp:1582-1586`; the post is at the `Info` level and not printed). The hits leave
  every report, the prefetch, the taxonomy lookups and psiblast's PSSM without a word, after they
  took their places in the hit list before the cut to `-max_target_seqs` (psiblast's Seq-aligns are
  made of the longer preliminary hit list, so its cut comes after: they take none there). Versions
  1.6 and 1.29 made this an explicit rejection; version 1.32 reproduces it in the six search
  programs, each of which drops such hit lists where NCBI makes the Seq-aligns
  (`SeqDb::keeps_line`, `algorithm::common::database::SeqAligned`; line DBID item lineless, probe of
  1,375 runs: 6 programs, 15 lists of the five kinds, 20 formats, `-max_target_seqs` 1 and 3, a
  query without hits). A read of such a sequence's lines still refuses, and no report reaches it.

- A sequence that a negative seqid list (`-negative_seqidlist`) keeps while removing every one
  of its definition lines (a sequence whose every line has only a gi or a ti, which
  `CSeqDBNegativeList::FindId` does not match by type; a `BL_ORD_ID` line of a volume whose first
  OID is not 0, of 1.29): NCBI's id set of the list is not blank (`CSeqDBImpl::x_InitIdSet`,
  seqdbimpl.cpp:1909-1946; on a version 5 database, the ids left after the gis are dropped), so
  the subject fetches with `check_oid_exclusion` (`seqsrc_seqdb.cpp:302-313`) exclude it: the
  gapped traceback's (blast_traceback.c:1567-1577), the composition-based redo's
  (blast_kappa.c:1378, 3560-3575) and the ungapped re-evaluation's of a nucleotide subject
  (blast_hits.c:2651-2665). The sequence is searched and counts in the statistics; a gapped
  search's preliminary hit list holds it, and it takes no place in the final one; an ungapped
  search of a nucleotide database frees its HSPs before the preliminary hit list. A protein
  subject of an ungapped search (blastx `-ungapped`) is fetched by none of them and leaves at the
  Seq-aligns, as a positive list's does. Versions 1.6, 1.29 and 1.32 made this an explicit
  rejection; version 1.34 reproduces it in the six search programs, each of which removes such
  subjects where NCBI fetches them (`SeqDb::excludes_lineless_at_fetch`,
  `algorithm::common::database::FetchExcluded`; line DBID item lineless2). Under such a list
  NCBI's fetch reads the header of every subject it fetches, so a header that cannot be read
  (a damaged record, a `BL_ORD_ID` line with a string tag in a later volume) ends NCBI's search
  there (exit 255 in the traceback, `BLAST engine error` in the ungapped re-evaluation) even
  where the subject is not reported; LOSAT refuses such a search at the same place. The
  composition-based redo fetches only the matches that its early termination test leaves, and
  LOSAT checks every subject of the preliminary hit lists, so it may refuse such a search where
  NCBI ends it without reading that header.

- TBLASTX, an ungapped extension that reaches the end of its query frame (version 1.40,
  2026-11-17, line RG item X2c, Owner-delegated): `s_BlastAaExtendLeft` and
  `s_BlastAaExtendRight` count the letters of the query batch, not of the frame (`n = MIN(s_off,
  q_off)` and `n = MIN(subject->length - s_off, query->length - q_off)` with offsets into the
  batch, `aa_ungapped.c:841,897`). A walk that reaches the end of its frame without a break scores
  the NULLB between it and the next frame (`BLAST_SCORE_MIN`, `blast_stat.c:1559-1576`; one NULLB
  between two frames with letters, `blast_setup_cxx.cpp:65-88`) and, while that does not break it,
  the next frame's letters, up to the batch's first letter or the end of its last frame. With a
  raw X-drop of 32,768 or less (with BLOSUM62 below about 15,000 bits of `-xdrop_ungap`) the
  NULLB always ends the walk with the values of a walk that stops at the frame's end. With a
  larger one, a right walk that reaches the frame's end above 32,768 goes on; where it scores no
  higher there, only its stop moves (`*s_last_off = s_off + i`), and the diagonal's `last_hit`
  then covers the next frame's hits on that diagonal (`aa_ungapped.c:593-597`), which NCBI does
  not extend: a pair of letters that a long similar pair just before the frame's end runs into is
  not reported. A left walk that scores no higher past the NULLB changes nothing. (The one-hit
  extension, `-window_size 0`, would start its right walk on the next frame's first letter after
  a word at a frame's end whose running score never rises above 0, `aa_ungapped.c:1030-1031,
  1073-1076`; such a word seeds only as an exact match of a subject word of the same letters,
  `blast_aalookup.c:498-504`, and must start with `X`, the only letter that scores below 0 against
  itself, while the seed stage's subject frames are translated from 2-bit bases,
  `blast_engine.c:1403-1405`, and hold no `X`: no such seed is found; probes with query frames
  ending in `XXX`, as `-subject` and from a database, gave the same output in NCBI BLAST+ 2.17.0
  and LOSAT.) LOSAT's walks go on into the batch in the same way
  (`tblastx::extension::QueryBatch`), in the reference frame loop and in the units; where a unit
  of a query batch split into several context groups sees a right walk go on past a NULLB, the
  subject's units run again in group order on one table of the batch's diagonals, so the hits it
  covers in another group's frames are left out as in NCBI's one table, at every `-num_threads`
  (`x_pair_par.rs`, module doc item 3). Fixture rows `x2c.tblastx.*` of
  `tests/outfmt0_manifest.tsv`; the unit tests of
  `LOSAT/src/algorithm/tblastx/extension/query_batch_tests.rs` compare the two-hit and one-hit
  extensions with NCBI's functions on the batch's buffer.
- `CHUNK_SIZE=1000` without `BATCH_SIZE`: the batch size becomes 0, NCBI prints the
  outfmt 0 prolog and then `BLAST engine error: Empty CBlastQueryVector`, exit 3.
- `-max_target_seqs` from 2^30 to 2^31-51: the preliminary hit list size
  (`MIN(MAX(2*h,10),h+50)` in `Int4`) wraps and becomes 10.
- A subject whose `GenerateDefline` throws (version 1.37, 2026-11-16, line RG item T7): a title,
  or for a database sequence with an empty stored title the scientific name from which NCBI makes
  it (since version 1.38 also the name that `x_AdjustProteinTitleSuffix` writes at the end of a
  protein title ending in ` [ ]`, line RG item T5c), in no encoding that `CUtf8::GuessEncoding` recognises (not UTF-8, with a byte 0x81, 0x8D,
  0x8F, 0x90 or 0x9D; `NStr::HtmlDecode` throws, `ncbistr.cpp:4528-4534`). The tabular formats
  ask it where `CShowBlastDefline::GetSeqIdList` replaces an id (`showdefline.cpp:218-242`): a
  `-subject` record whose title's first word starts with `Subject_`, or read with
  `-parse_deflines` with the id `lcl|Subject_N` (NCBI keeps the parsed local id and generates the
  title from the rest of the defline: `s7` for `>lcl|Subject_9 TPA: s7 third ~`) (`SetSubjectId`,
  for `sseqid`, `sgi`, `sacc` or `saccver`), and a database sequence without parsed Seq-ids or
  with an `lcl|Subject_N` id (`SetSubjectId` where the subject Bioseq's ids hold it, and
  `x_SetSubjectIds` for `sallseqid`, `sallgi` or `sallacc`, which pushes the ids of the
  definition lines before the first that holds it). The catch of `SetFields`
  (`tabular.cpp:862-875`) adds the Seq-align id to the `sall*` lists of the rows before
  (`SetSubjectId(list)`, `tabular.cpp:1374-1378`; after `x_SetSubjectIds` cleared them and pushed
  those lines where that threw) and returns: the row prints `Unknown` for `sseqid`, 0 for `sgi`, the scores, the query's
  ids, `qcovs`, `qcovus` and `qcovhsp`, `qlen` where a query id field or `qcovhsp` is asked (else
  0), `slen` 0 (where only `x_SetSubjectIds` threw, the length and the subject's own `sseqid`,
  `sgi`, `sacc`, `saccver`), the zeros of `x_ResetFields`
  (`pident` 0.000, `ppos` 0.00, `frames` 0/0, `sstrand` N/A, empty letters and BTOP), and the titles
  and taxonomy of the last row of the query's (or psiblast round's) report that found its subject
  (`N/A` before one). In XML the catches write the Seq-id: outfmt 5 `Hit_id` `lcl|Subject_N`,
  `Hit_def` and `Hit_accession` `Unknown`, `Hit_len` 0 (`blastxml_format.cpp:398-407`), a query's
  Seq-loc id (`lcl|Query_N`; with `-parse_deflines` the best of its ids, `lcl|qbad`,
  `gb|AB000099.1|`, `gi|4504347` of `gi|4504347|ref|NP_000549.1|`) and `No definition line`
  (`:550-566`); XML2 and JSON a `HitDescr` with the id alone and the sequence's length
  (`blastxml2_format.cpp:410-415`), a query's id without a title (`:584-596`). Outfmt 0, 11 and 18
  were NCBI's before (`Unknown`, the line in place of each alignment), but for a `-parse_deflines`
  record with an `lcl|Subject_N` id, whose line names the Seq-align id (`Subject_9`, LOSAT wrote
  `Subject_1`). BLASTX reads such a title as its bytes, as the other programs do. Each stream of
  LOSAT Web is the CLI's run of its format: where the default tabular formats stop (rule 3 below)
  and the argv asks for another format, streams 6 and 7 fail alone and the others are written in
  the same search; with the argv's own `-outfmt 6` or `7` the run fails, as the CLI's
  (`docs/web/abi_v2.md` §7). So does every error of one format's report, at a check before its
  rows or in its middle (`-sorthits` with a description total above 99999 bits in outfmt 0, an
  NCBI error of outfmt 0 alone such as BLASTX's `Row is empty`), and a stream that fails alone
  after the host received a chunk of it (BLASTX's later query batch) is discarded by the host at
  the run's signal (`docs/web/abi_v2.md` §5, line RG T7-W).
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

- A version 5 database list (several names in `-db`, or the `DBLIST` of an alias file) that
  mixes databases with and without the accession table (`acc2oid`) of `makeblastdb
  -parse_seqids`: NCBI searches it, and before it writes a batch's reports (every format but
  the Seq-annots and the archive) its BlastDb data loader looks each hit's Seq-id up in every
  LMDB file of the list (`BlastFormatter_PreFetchSequenceData`,
  `CSeqDBLMDBSet::AccessionToOids`); a hit whose Seq-id it looks up as a string (an accession, a
  name, a giim; not a gi, not `gnl|BL_ORD_ID|n`) then stops the run with `BLAST Database error:
  DB contains no accession info.` (exit 2), after the outfmt 0 header and the earlier batches. Where `PRE_FETCH_SEQS_LIMIT` turns the prefetch off, the reports write
  `Unknown` for such a hit, and then outfmt 0 to 4 and 17 stop with that error, outfmt 6
  crashes (`CCoreException`, exit 255), and 5, 7, 10 and 18 go on (exit 0). The JSON and XML2
  reports throw at that `Hit` (`s_SeqAlignToXMLHit`: `GetBioseqHandle`, then `GetLength` outside
  its catch): 15 and 16 catch it for the query (`BlastJSON_FormatReport`,
  `BlastXML2_FormatReport`), post `Error: [<program>] ` with the exception's text and trace, leave
  the query's report out (its `,` separator stays) and go on with the next query (exit 0); 13 and
  14 leave the query's file empty and stop with that error (exit 2).
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
  NCBI's). Before the error come the taxonomy database's texts of the loads that the prefetch
  makes before the stopping lookup (since version 1.41, line DBID item loadord, probe of 184
  runs, all six programs): it loads a batch's Seq-align ids, up to 200 of them, in the order of
  their Seq-id handles, those that its mapper packs (gi, an accession of a prefix and digits
  without a name or release, a `general` id with an integer tag or a string tag of at least
  three digits; a tag of 0 or less, such as `gnl|BL_ORD_ID|0`, after every positive number) by
  their number and the others after them, and ends at the first lookup that fails. It refuses a
  run where the reports read such a Seq-id with the prefetch off (`... whose lookup in the
  database fails (DB contains no accession info.; NCBI BLAST+ 2.17.0's prefetch stops with a
  BLAST Database error, and where the prefetch is off its formatters fail to read it)`) but
  writes NCBI's results where every report of the run is outfmt 13 to 16 (since version 1.47,
  line DBID item pfoff: the texts of the loads before the throwing `Hit`, none of its own or of the
  `Hit`s after it in that report, and a text that NCBI writes once at its next load, a later `Hit` or
  a later batch's prefetch) or outfmt 6, 7, 10, 18 and 20 (since version 1.50, line DBID item
  pfoff-2: the row that the `catch` of the subject block leaves, or the stop at a row that reads the
  sequence after that block; the organism report leaves its HSPs out; a throwing row loads nothing,
  and the names and tax ids are looked up only for the sequences that a report loads; refused at
  the first row of a query's report with `sacc`, `saccver` or `staxid`, below) or outfmt 0 to 4
  (since version 1.51, line DBID item pfoff-3: the description line of the Seq-align id and
  `Unknown`, then the stop at outfmt 0's first alignment that loads such a sequence under a
  Seq-align id that is not local, or at any alignment of outfmt 1 to 4 that loads one, or outfmt
  0's skipped line for each of its HSPs; the description lines load the sequences that they show,
  in psiblast's rounds after the first those of the two-part table of the first
  `-num_descriptions` HSPs, then the alignments theirs up to the stop), a
  stopping batch where a load whose order to the stop NCBI takes from the addresses of its
  handles (the same packed number as the stopping id's, or neither packed, which changes from
  one command line to another; a Seq-align id of two sequences; every load of a batch of more
  than 200 Seq-align ids, which NCBI sorts by `CSortableSeq_id` and loads in pieces of 100)
  writes a text of the taxonomy database (after the lookups of the batches before: their
  prefetch, or in the JSON and XML2 formats their reports, since version 1.36), one where
  another lookup that throws on its id may come before the stop, and one whose stopping ids NCBI
  may all never look up: another sequence of the database that holds an id of its Seq-id handle
  may be loaded first (NCBI's data loader keeps the ids of every sequence it has loaded and
  gives that sequence for them): `gi|1|ref|NP_000002.1|` before `ref|NP_000002.1|` in the
  batch's order, or a sequence that the run loaded before (in a batch, a query or a psiblast
  round before: every subject of its Seq-aligns where its prefetch is on, psiblast's whole
  traceback beyond `-max_target_seqs` too, else, with `PRE_FETCH_SEQS_LIMIT`, the subjects that
  the rows of the reports NCBI writes read, of the hits that they show after the formatter's cuts
  to the hitlist size (`PruneSeqalign`: to `-max_target_seqs`; in the pairwise report, outfmt 0,
  of its description lines to `-num_descriptions` and of its alignments to `-num_alignments`,
  whose subjects it loads both, review 7 of item loadord; none for `6 qseqid evalue`, none of a
  BLAST XML report whose iteration NCBI's XML stream cut at a query title or id outside US-ASCII,
  review 6 of item loadord), and those that psiblast's PSSM includes) where the lookup of its
  Seq-align id gives it first
  (the loader loads the first sequence that the lookup gives: of two holders of one gi, the
  first; review 4 of item loadord), or the first sequence that a lookup of an id of such a
  sequence gives; not a sequence that no batch loaded (NCBI then stops), nor a copy with
  the same Seq-align id, whose load is the same lookup; and, in a call of the library's query
  batches after the first, which has no record of the calls before, any such sequence. Or a query that NCBI's scope holds at the prefetch has the stopping id, or one
  that it matches weakly: the records that NCBI's reader read for the batch, a record that
  `-query_loc` drops included (it keeps its number), by their generated ids (`lcl|Query_N`, the
  numbers of the whole run, also in psiblast's rounds) or with `-parse_deflines` by the ids of
  their deflines alone, or the query of a PSSM (`-in_pssm`, `lcl|myprot`); with outfmt 5, whose
  report keeps the scope, also those of the batches before (reviews 1 and 2 of item loadord:
  NCBI then searches to the end, exit 0, or blastx ends with its CSeqalign exception). Where a
  batch has several stopping ids, NCBI's loads end at the first lookup that it makes, in the
  order of the handles (review 3 of item loadord): a stopping id before it that NCBI may skip so
  does not refuse the batch by itself, but the loads from the first stopping id to the one that
  NCBI surely looks up, and the sequence that the skipped lookup is given instead, are loads that
  NCBI makes only where it skips the lookups before them, which LOSAT does not know (refused
  where one writes a text), and the stops that the loads may end at must fail with one message.
  A volume of a version 4 database whose string ISAM (`.psi`/`.psd`, `.nsi`/`.nsd`) or PIG ISAM
  (`.ppi`/`.ppd`) has a file that is not a regular file (a directory) has no such ISAM, as NCBI's
  `CSeqDBIsam::IndexExists` (`CFile::Exists`, a regular file) says: the lookups find nothing in
  that volume and go on to the next (review 4 of item loadord: NCBI exit 0 where a later volume
  holds the id, 255 where none does; LOSAT stopped with NCBI's `Error: Could not open input
  file`, which NCBI writes for a regular file that does not open). Where its index or data file
  is empty, NCBI crashes where its lookup reads it (SIGSEGV), which LOSAT refuses. NCBI's
  lookup (`SeqidToOid`, `SeqidToOids` without `multi`) returns at the first OID it finds, in the
  order of the volumes, so an ISAM file of a later volume that does not open (a regular file of
  mode 000) stops a lookup only of an id that no volume before it holds in the search's list
  (review 3 of item loadord, H3, and item v4late: NCBI exit 0 where LOSAT stopped with `Error:
  Could not open input file`); LOSAT, which also lists the other sequences that hold the id,
  refuses where that list fails on such a file and a report shows what it can change
  (`PD-LOSAT-DATABASE-SEARCH`).
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

- TBLASTN with `-window_size 0` (version 1.18, 2026-11-14, line RG item X2, Owner-delegated;
  version 1.48): with a low `-threshold`, NCBI's one-hit ungapped extension can give a negative
  `Int4` length, which `BlastGetStartForGappedAlignment` reads as `Uint4`
  (`aa_ungapped.c:1054,1083`, `blast_gapalign.c:3393-3413`): it scores the first window of 11
  letters from the HSP's start, which may pass the end of its context and its frame. Inside NCBI's
  buffers the window reads the NULLB and then the next contexts' or frames' letters. LOSAT reads
  the subject's bytes as NCBI does (the next frames' letters, which a PSSM's '-' column can score)
  and the query's context followed by NULLBs: the matrix scores the context's NULLB
  `BLAST_SCORE_MIN`, so what follows it never changes the start, and a PSSM reads its rows, not the
  query. Past the buffers LOSAT reads NULLBs (approved exception 3, version 1.48; 1.18 refused
  these searches): past the byte after the last context's NULLB of the query batch
  (`calloc(buflen+1)` with `buflen` from `QueryInfo_GetSeqBufLen`, `blast_setup_cxx.cpp:500-501`,
  `blast_query_info.c:246-250`) and past the NULLB after frame -3 of the subject's translation
  (`BLAST_GetAllTranslations` writes frames 1, 2, 3, -1, -2, -3 into `2*(nucl_length+1)+2` bytes
  and leaves the last one unset, `blast_util.c:1058,1063,1098-1101`). LOSAT refuses (exit 1):
  - with `-in_pssm`, a window past the PSSM's last position (`pssm->data` holds one row per
    position and no sentinel row);
  - with `-in_pssm` and a PSSM whose '-' column scores above BLAST_SCORE_MIN (below), a gapped
    start that the window's sum puts past the NULLB after frame -3, where NCBI's preliminary gapped
    alignment reads the byte after that NULLB, which the buffer allocates and never sets, or past
    the buffer (version 1.48; exception 3 covers the window's letters, not that alignment's;
    fixture `db.x2.tblastn.pssmpast.rej`).

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

  With `-in_pssm` and a PSSM whose '-' column scores above BLAST_SCORE_MIN (a hand-made one: the
  matrices score it BLAST_SCORE_MIN, `blast_stat.c:1559-1562`, and so does NCBI's PSSM engine, the
  gap's background probability being 0, `blast_psi_priv.c:2435-2436`), the window can sum above 0
  over the NULLB after a subject frame's last letter, and the gapped start then lies on that NULLB
  or one or two letters after it (version 1.44, 2026-11-17, line RG item X2d, Owner-delegated).
  NCBI's preliminary gapped alignment reads back from the start over the six frames' buffer (the
  next frame's letters, the NULLB scored with the '-' column, the frame's letters; no right
  extension, `blast_gapalign.c:4250-4347`), and LOSAT reads the same bytes; an HSP that NCBI drops
  before its traceback leaves the search as NCBI's. LOSAT refuses where such an HSP reaches the
  traceback with its start after the NULLB: `BLAST_GappedAlignmentWithTraceback` reads the subject
  at the start first (`blast_gapalign.c:4609`, `:569-575`), a byte past the fence of the partial
  translation (`blast_hits.c:1208,1225-1227`) and, on the retry after that fence, past the NULLB of
  the whole frame's translation (`subject.offset` -1: `blast_traceback.c:424-428`,
  `blast_hits.c:1175-1177,1189,1208`, no fence; the NULLB of `blast_util.c:453`) or past its
  allocation. A start on the NULLB is defined there (the fence, then the whole frame's NULLB); LOSAT
  does not reproduce the reports of that alignment's column of the NULLB yet and refuses it where
  they would show it, as a gap of its own, not under this rule.

  Version 1.48 extends BLASTP's approved exception 3, which reads the letters past a sequence as the
  sentinel, to TBLASTN's window (above).

- BLASTX with `-window_size 0` (version 1.22, 2026-11-14, line RG item X2f, Owner-delegated;
  version 1.48): with a low `-threshold`, NCBI's one-hit ungapped extension can give a negative
  `Int4` length (`aa_ungapped.c:1054,1083`), which `BlastGetStartForGappedAlignment` reads as
  `Uint4`: it scores the first window of 11 letters from the HSP's start
  (`blast_gapalign.c:3394-3437`). A window that
  passes the end of its context or of its subject holds that NULLB, which scores
  `BLAST_SCORE_MIN` (`blast_stat.c:1559-1579`), so NCBI returns the HSP's start whatever bytes
  follow, and so does LOSAT. Where the window reads outside NCBI's buffers LOSAT reads NULLBs
  (approved exception 3, version 1.48; 1.22 refused these searches, exit 1):
  - past the byte after the last context's NULLB of the query batch (`calloc(buflen+1)` with
    `buflen` from `QueryInfo_GetSeqBufLen`, `blast_setup_cxx.cpp:500-501`,
    `blast_query_info.c:246-250`);
  - past the NULLB after a `-subject` protein (`malloc(len+2)`, `blast_setup_cxx.cpp:885-893,
    1236-1239`);
  - past the last memory page of a database volume's sequence file. NCBI maps the file whole from
    offset 0 and reads the sequences in place (`seqdbvol.cpp:1802`, `seqdbatlas.hpp:700-722`,
    `seqdbatlas.cpp:114-127`, `ncbifile.cpp:5642-5643`), and the system fills the rest of the
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
  page of a 4,096-byte file made valgrind report an invalid read and crash; there LOSAT's output is
  NCBI's on the same volume with 16 zero bytes after that page (fixture `db.x2f.bx.page.6`).
  With version 1.48 LOSAT runs all 400 searches: it equals NCBI under valgrind in the 115 that read
  past a buffer (valgrind reports invalid reads in each) and NCBI's native output in the other 285
  (34 of them `-ungapped`, which both refuse with composition statistics); natively NCBI's report
  differs from LOSAT's in 2 of the 109 past ones it finishes.

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
  search (exit 1) where that merge would happen, after the reports of the query batches before
  (since version 1.40, as the refusal of an extension past a query frame's end; before, they were
  written only where the reports were streamed). Every other merge of such HSPs is NCBI's (`Int4`
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

- BLASTX `-outfmt 1` to `4` with `-ungapped -window_size 0` (version 1.39, 2026-11-16, line RG item
  T3c after X2g and T3c fix 3, Owner-delegated, rule 3 as 1.23): an HSP of length 0 among the
  aligned subjects whose segment of length 0 NCBI's merge of the section frees while it still reads
  it. A row's map of starts (`CAlnMixStarts`) holds its segments by `CRef`, `m_Segments` by raw
  pointer (alnsegments.hpp:67,80,89), so a write into a map that replaces a segment's last `CRef`
  frees it. That happens in two places:
  - The split. An HSP of length 0 whose query row starts at the first nucleotide of a codon of an
    HSP of the query's frame +1: its match (merged with that frame's HSPs, above) starts at or
    inside a segment of letters of that row, and `x_Merge` "splits" the segment at the match's
    length of 0 (alnmerger.cpp:580-614): `SetStarts()[start + 0] = seg` replaces each row's
    reference to the segment, whose rows are all on the query's strand, so the last one frees it
    while `ITERATE` still walks its `m_StartIts`. SIGSEGV in every such probe (37 of the 1,880
    cases in each of the four views, 98 of the grid's 277 searches, the second query of
    `x2g_ord_q2.fna`).
  - `FillUnalignedRegions` (alnsegments.cpp:286-334). An HSP of length 0 whose segment its query
    row's next segment does not follow at once (unaligned letters between them): the fill of that
    stretch starts at the segment's own start (`new_start = prev_start + prev_len * width`, with a
    length of 0), so `row->SetStarts()[new_start] = seg` (:323) replaces the row's reference to the
    segment. It is the last one unless the segment is the last that `x_Merge` made (its local
    `CRef<CAlnMixSegment> seg`, alnmerger.cpp:161, alive until `x_CreateDenseg` has returned,
    :822-830). Besides `seg`, only the rows' maps hold a segment then, and none holds one of length
    0 elsewhere: on the first row, whose segments `x_MinimizeGaps` never takes
    (alnsegments.cpp:254-256), such a segment has one start (the walk of a match's second row does
    not run for a length of 0, alnmerger.cpp:752); the first frame group's other query rows hold
    segments of length 0 only, which `x_MinimizeGaps` makes into one segment per run that nothing
    on those rows follows, so no fill starts at one. (A map that still holds a replaced segment at
    another start of its row keeps it alive: where a run of `x_MinimizeGaps` ends at a segment that
    shares a row with an earlier one, each new segment is bound at both starts,
    alnsegments.cpp:478-480, :526, and a later fill can write at the earlier one. That happens to
    segments of letters, as in an alignment with an insertion next to a deletion, and frees
    nothing.) `x_CreateDenseg` then reads the freed segment through `m_Segments`
    (alnmerger.cpp:1196-1210). Official NCBI BLAST+
    2.17.0 exits 139 (SIGSEGV), and valgrind finds the invalid reads of a block "free'd by
    `CAlnMixSegments::FillUnalignedRegions()`" (review 3 of T3c: the 8 runs of its seed-11 fuzz
    with `-strand minus -num_descriptions 0`, cases `r11_001` and `r11_019`); where the segment is
    still held, or no fill starts at it, NCBI exits 0 and valgrind is clean (`r21_173`, `r21_322`,
    `r21_381`, review 4's `r41_079` on the default strand and the x2g case f3147 at query offset 0,
    held by `seg`; `r21_470`, `r11_003`).

  No checkable result in either place. LOSAT runs the same merge and refuses (exit 1) where it
  reaches that split or replaces that last reference
  (`report::query_anchored::aln_mix::MixError::NcbiFreesSegment`, the predicate
  `aln_mix::Mix::held`), before the first byte of the reports of the search's call (as 1.23),
  unless the report of a query before it ends with NCBI's exception, where NCBI stops first
  (`pairwise::blastx_anchored_stop`; fixture rows `db.x2g.bx.of.p0127.nd0.1`,
  `db.x2g.bx.of.ord.n1a5.1`, `db.x2g.bx.of.ord2.n1a5.3`, `db.x2g.bx.of.r11_001.nd0m.1`,
  `db.x2g.bx.of.r11_019.nd0m.3`).

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

- The tabular `sacc` or `saccver` of a row whose subject lookup fails as in rule 2 above (version
  1.37, 2026-11-16, line RG item T7, Owner-delegated): outfmt 6, 7 and 10 (and 20) with them in the
  field list, which the default fields of 6, 7 and 10 have (`saccver`). NCBI BLAST+ 2.17.0 takes the
  label of `FindBestChoice` of the empty `m_SubjectId`, a null pointer (`tabular.cpp:174-185`;
  `CCoreException::eNullPtr`, `ncbiobj.cpp:1010`), and stops with exit status 255 and a stack trace
  of its build after the row's earlier fields. Also `staxid` of such a row when it is the first row
  of its query's (or round's) report: `m_SubjectTaxId` is printed before any row set it
  (`tabular.cpp:262-268`; the member is not initialised, a different number in each run). LOSAT
  refuses (exit 1) before the query batch's report, where it refused every such subject before
  1.37, so its output is a prefix of NCBI's. A failed lookup does not depend on the other
  subjects: one alone, the first, the last or two give the same rows. Since version 1.50 (line
  DBID item pfoff-2) also the first row of a query's report whose subject's load throws where the
  batch's prefetch is off (`DB contains no accession info.`, a version 5 list of databases with
  and without accession tables): the `catch` of the subject block leaves `m_SubjectId` empty
  before any row of the report loaded its subject, and `m_SubjectTaxId` unset; a later such row
  prints what the last loaded row set.
- psiblast and tblastn `-in_pssm` (version 1.45, 2026-11-16, line RG item D3-1, Owner-delegated,
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
    `seq_annot_info.cpp:761, 825, 898, 971, 1057, 1138`) prints a heap address;
  - a Seq-table row without a Seq-id or a location: NCBI's `CCoreException` NULL-pointer report
    holds a stack of addresses;
  - two Seq-ids of one location (a Seq-table row's, and since line RG item D3-2 also a feature's
    or a Seq-table's, where one of them has an empty region and the other NCBI's `Failed to parse
    location of` message) that get different messages and that no
    gi, accession or general id number orders: NCBI writes them in the order of its
    `CSeq_id_Handle`s, which is then the address of their info (`seq_id_handle.hpp:207`);
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
  - BitMagic data of a Seq-table that NCBI decodes (`bit-bvector` values, a `bit-set-bvector`'s
    sparse index) in a stream that NCBI reads past or crashes on (well-formed data, which NCBI
    decodes deterministically, is a rule-2 rejection below). A `bit-bvector`'s size, which NCBI
    reads without a decode, is NCBI's;

  The `-in_pssm` file's format guess (NCBI's `CFormatGuess`, ported) is refused where one of its
  tests reads memory that the file does not decide: `TestFormatGvf`'s ninth column of a GVF line
  of eight (`tokens[8]`), `TestFormatPhrapAce`'s tokens of an `AS` line of fewer than three
  (`values[1]`, `values[2]`; NCBI 2.17.0 crashes on both, exit status 139) and `TestFormatFasta`'s
  search for `|` past a test buffer of 10 bytes or fewer.

  None of these is a rule-1 exception: D3-1's questions to the Owner, Q1 (the Dense-seg division by
  zero), Q2 (the memory bounds), Q3 (the signed overflows) and Q-BM, were decided
  by the Owner on 2026-11-16 (「推奨通り」): Q1 to Q3 refuse their cases, as built, and Q-BM ports BitMagic
  in a later item (refused until then, below). Not a rule-3 item but a refusal: the XML and JSON
  reports of an untitled query whose annotations hold features (NCBI's definition line reads them;
  not ported). Since line RG item D3-2 LOSAT writes NCBI's ASN.1 text of the values it reads: the
  annotation's text that NCBI prints after its `Failed to parse location of` message (1204.6) is
  written, and with `-out_pssm` and `-export_search_strategy` every member of the PSSM. A value of
  -2147483648 (`Int4`) or -9223372036854775808 (`Int8`) that the PSSM holds is written as its
  digits: NCBI's `COStreamBuffer::PutInt4` and `PutInt8` negate it with a signed overflow, and its
  build writes `-8` and `-0`: a result of undefined behaviour, which LOSAT does not reproduce under
  the Owner's rule of 2026-10-10 and 2026-10-22 (NCBI's signed overflows widened), as line DI-4's
  decision 13 writes the archive's `WordThreshold` through the same `PutInt4` (Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07). The rest is
  NCBI's (rules 2 and 4): the index's messages (`Invalid 'lens' size in denseg`, the Seq-table
  column messages, `Failed to parse location of`, after which NCBI's later messages print their
  headers only and an empty region the object's ASN.1 text as far as NCBI's 4096-byte output
  buffer writes it to the failed stream, exit status 0) and exceptions
  (`Error: NCBI C++ Exception:`, `Error: vector::reserve` of a Dense-seg or Seq-table reserve), the
  reader's format errors, number-range and member-tag reports, its `#` for a byte outside the
  visible ASCII range and the deprecated members of a `Variation-ref`; where the index's
  exception ends the run, each output format's bytes (tblastn's `-outfmt 11` error archive with
  every message saved whole with its severity, psiblast's `-outfmt 10` and `20`), whatever NCBI
  reads only later (the PSSM's values and the alphabet of its letters, the `-subject` records'
  Seq-ids under `-parse_deflines`, the search's options), and no PSSM file of `-out_pssm` and an
  empty file of `-export_search_strategy`, as with the reader's errors; a query Seq-entry set,
  whose annotations NCBI's object manager takes (the members' first as it makes them, the set's own
  first as it indexes them) and at which NCBI then stops (`CSeq_entry::GetSeq()`), and a Seq-id
  twice in the query's entry (the attach's `duplicate Bioseq id` exception, before the index, for
  every two Seq-ids that NCBI's Seq-id handles take as one, with NCBI's spelling of each; with
  `-subject`, where a subject record's spelling may show in the message, a rule-2 rejection below)
  (fixture rows `db.ps.inpssm.d3.*`, `db.tbn.inpssm.d3.*`, `db.ps.inpssm.fg.*`, `db.tbn.inpssm.fg.*`).
  Values nested deeper than NCBI's stack holds are read to any depth: exception 15 (rule 1, above).

- TBLASTX, an ungapped extension that scores higher past the end of its query frame, in the frame
  next to it in the query batch (version 1.40, 2026-11-17, line RG item X2c; the Owner's decision
  of 2026-11-16 20:00 JST on line RG's question Q-X2c, 「ごめんAで」: refuse, the results of the query
  batches before are written; line RG `DECISIONS.md` X2c-1 and X2c-Q; the walks are those of the
  rule-2 item above): the HSP then spans two frames. NCBI keeps it in the frame of its seed with
  query offsets outside that frame (`s_AdjustInitialHSPOffsets`, `blast_gapalign.c:2383-2392`: a
  start before the frame fails its `ASSERT(init_hsp->ungapped_data == NULL ||
  init_hsp->ungapped_data->q_start >= 0)`, compiled out of release builds; an end past the frame is
  not tested), re-evaluates it on the batch (`Blast_HSPListReevaluateUngapped`,
  `blast_hits.c:2609-2725`), and its formatters take its coordinates from those offsets
  (`blast_seqalign.cpp:1344-1376`). One that starts in the frame before a minus frame is written
  with query coordinates past the query and no identities (`0.000 7040 7040 0 42245 21126 1
  21120` for a query of 21,122 bases, input `x2c_e7k`), one that ends in the frame after its own
  with an end past the query (`14404 43206` of 28,803, `x2c_g48`), and one of a plus frame that
  starts in the frame before is left out of every report (`x2c_d7k`; where the Seq-align layer
  drops it was not found), while it takes its place among the HSPs and its extension covers the
  hits it passed. LOSAT refuses the search (exit 1, "an ungapped extension that passes the end of
  a query frame into the frame next to it in the query batch and scores higher there ... is not
  supported by LOSAT's TBLASTX") where such an extension first scores higher: the reports of the
  query batches before are written, then the failing batch's messages and the error, as after an
  NCBI error of a batch, in every format and whether or not the reports are streamed, the same at
  every `-num_threads` (`tests/tblastx_query_batch.rs`: 1, 4, 8 and 16 threads, the one-hit word
  finder, and a run of two query batches). It takes a raw X-drop above 32,768 (with BLOSUM62 about
  15,000 bits of `-xdrop_ungap` or more).

- TBLASTN, an alignment on a subject (or its `-subject_loc` range: NCBI's `subject_blk->length` is
  the searched range's) of fewer than 12 bases that reaches the traceback (version 1.46,
  2026-11-17, line RG item X-short, Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07). `Blast_HSPGetTargetTranslation` translates a frame only when `start || stop <
  length/3 - 3` (`blast_hits.c:1161-1163`), false from the start for 11 bases or fewer, so it
  returns the translation it never made (`calloc`'d to NULL, `blast_util.c:1284`) plus 1
  (`:1233-1234`), and every caller reads from that address on next: the traceback's gapped start
  and `ALIGN_EX` (`blast_traceback.c:429-440,509-512`), the identities after the composition-based
  redo (`blast_kappa.c:527-528`, which run only for an alignment left by the redo's E-value purge,
  `:3677-3700`), the ungapped re-evaluation (`blast_hits.c:2700-2705`; LOSAT's TBLASTN refuses
  `-ungapped`). NCBI BLAST+ 2.17.0 dies with SIGSEGV in every such run of ours (ws and the local
  official binary, 2026-11-17: 11 bases in frames +1 and -1, `-subject`, `-subject_loc` and `-db`,
  a matrix query at `-comp_based_stats` 0 to 3, `-in_pssm`, `-outfmt` 0, 5 and 6, a second query, a
  12-base subject before or after it; gdb: `s_Blast_HSPGetNumIdentitiesAndPositives` under
  `Blast_RedoAlignmentCore_MT`, `ALIGN_EX` under `Blast_TracebackFromHSPList`; an independent
  review's 850 random runs and about 150 targeted ones found no other outcome). LOSAT refuses at
  that read (`TargetTranslation::get`, which `x_whole_window` defers to for such a subject), where
  dev panicked. A subject that leaves no alignment to that point (cut by a cutoff, the redo's purge
  or the hit list; at the default word size, 9 or 10 bases in the probes; with `-word_size 2` they
  reach it too) searches as NCBI's, and so does every run at a two-hit window in our probes (no
  frame of such a subject, three letters at most, holds two words that do not overlap); 12 bases
  are translated (rows `db.xs.tblastn.*`). With `-subject` and a `-db_gencode` other than 1, LOSAT
  translates in the selected code (the approved deviation PD-TLOSAN-LOCAL-GENCODE-32, whose
  reference is NCBI's `-db` search) and refuses where NCBI's `-db` search crashes, although NCBI's
  `-subject` search, translating in code 1, may find nothing there. LOSAT's TBLASTN refusals write
  no report, also where NCBI writes the reports of query batches before its crash (`BATCH_SIZE`).

## Explicit rejections of a taxonomy record past NCBI's mapping (rule 3, version 1.43)

`CSeqDBTaxInfo::GetTaxNames` maps `taxdb.btd` and reads a record from its offset to the next
entry's (the last entry's: to the end of the file) without checking it against the file's size
(`seqdbtax.cpp:258-298`). Past the end of the file NCBI reads the zero bytes of its mapping's last
4,096-byte page, which LOSAT reproduces, and after that page memory that does not belong to the
file. LOSAT refuses (exit 1, one line naming its program) only:

- every lookup that reads a record of an empty `taxdb.btd`, which NCBI maps as a null pointer: it
  faults for certain (SIGSEGV in every program; the loads of dustmasker, segmasker, blastdbcmd's
  FASTA and the searches' prefetch included);
- a name that the output writes and NCBI reads past the mapping: the super kingdom of a record
  that reaches past it (`blastdbcmd -outfmt %K`, `-tax_info` with `%K`, `sskingdom` and
  `sskingdoms`), and any name of a record whose third tab is past it (whether NCBI finds names at
  all depends on that memory): a title made from the scientific name, `%S %L %B %C %N`, the
  tabular name fields, the `sciname` of outfmt 13 to 16, outfmt 18. Line RG item T5b's review
  measured NCBI 2.17.0 writing 4,835 bytes and more of other mappings as the kingdom of such a
  record (`-tax_info -outfmt "%T|%K|%B"`), the same in two runs.

## Approved exception (a taxonomy record past NCBI's mapping of `taxdb.btd`, version 1.43)

Version 1.43. The Owner decided on 2026-11-15 (line RG `DECISIONS.md` T5b-6, Q1 = A, words
「RG-T Q1=A」) the question of line RG item T5b's fix 2: LOSAT keeps writing NCBI's completed output
here, as rule 1's approved exception, against rule 3's rejection of every lookup of such a record,
which would also refuse runs that NCBI completes with the same bytes (44 runs of the item's probes).
Numbered after exception 13 (version 1.35).

14. **A record that NCBI reads past its mapping where the output writes nothing of the bytes
    there.** Where NCBI reads such a record past its mapping and the output writes nothing of the
    bytes there (the lookups of `GetBioseq`, whose names the maskers' and blastdbcmd's output do
    not write, a stored title, the names before a third tab inside the mapping, the searches'
    prefetch), LOSAT writes what NCBI writes when it completes, which is the result of the same run
    with a `taxdb.btd` that holds those bytes (rule 1's valid result). Whether NCBI completes
    depends on its memory map, not on the input: in the review of line RG item T5b, a record of a
    266-byte `taxdb.btd` that ends at 5,000 bytes completed in `dustmasker`, `segmasker` and
    `blastdbcmd`, one that ends at 9,000 or 20,000 bytes in `blastdbcmd -tax_info`, and one that
    ends at 1 MiB faulted (SIGSEGV) in `dustmasker` and `blastdbcmd`. LOSAT writes the output in
    all of them, as in the 20 of the item's 115 search runs and 2 of its 32 tool runs where NCBI
    faulted.
    This was LOSAT's behaviour for the maskers' and blastdbcmd's loads before line RG item T5b
    (dev), and is now its behaviour in the searches as well. What the output writes from past the
    mapping stays refused (the rejections above).

## Explicit rejections of deterministic results (rule 2, version 1.10)

On a database with parsed Seq-ids NCBI gives these results without failing, and they look wrong
or need NCBI's object manager to reproduce. LOSAT refuses them instead (line DBID
`DECISIONS.md` 9, 20 and 21, line DI-7 for the query's id), under the maintainer's standing
delegation, which the maintainer accepted on 2026-11-02 with version 1.10. Each is a row of the table
in [`PD-LOSAT-DATABASE-SEARCH`](PD-LOSAT-DATABASE-SEARCH.md), which also lists refusals that rule 4
would port. The last two entries (version 1.45, line RG item D3-1) are not of these: they concern a
PSSM query's annotations and Seq-ids, not a database; they are Owner-delegated (standing instruction
2026-09-29, restated 2026-10-07), not under the acceptance of 2026-11-02; they have no row in
PD-LOSAT-DATABASE-SEARCH; and the first is a refusal until a port lands, which the Owner decided on
2026-11-16 (「推奨通り」, Q-BM):

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
  BLASTP, PSIBLAST and RPSBLAST, 1 to 4 of BLASTN and TBLASTN since version 1.28, of BLASTX and
  TBLASTX since version 1.39), of a query read with `-parse_deflines` in the query batch that
  reports it (with outfmt 5, in any batch), and in outfmt 1 to 4 of BLASTP, BLASTN, BLASTX, TBLASTN
  and TBLASTX with the id or display id `lcl|Query_N` of a query read without it (TBLASTN's merge
  then reads the protein query's letters as nucleotides, probe `tn.t3b.d.qid.f1`; BLASTX's and
  TBLASTX's merges of the query's rows with the sequence's throw, exit 255, probes
  `bx.t3c.d.qid.f1`, `tx.t3c.d.qid.f1`): NCBI shows the query's title and letters for the
  id (its second lookup, by the display id, the query's title in that line).
- A stored title or Seq-id with a byte outside ' ' to '~': NCBI prints `#` and an error message
  each time it reads the definition line.
- psiblast and tblastn `-in_pssm` (version 1.45, line RG item D3-1, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07): well-formed BitMagic data of a Seq-table of the
  PSSM query's annotations that NCBI's object manager decodes (`bit-bvector` values, a
  `bit-set-bvector`'s sparse index). NCBI's decode is deterministic; LOSAT refuses it ("a PSSM
  query whose annot has ... is not supported") until that port lands: the Owner decided on 2026-11-16 (「推奨通り」, Q-BM) to
  port BitMagic in a later item. Streams that NCBI reads past or crashes on are rule-3
  rejections (above).
- psiblast and tblastn `-in_pssm` with `-subject` (version 1.45, line RG item D3-1,
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07): a Seq-id twice in the
  PSSM query's entry whose `duplicate Bioseq id` message may spell a subject record's Seq-id (NCBI
  reads the subjects into its scope first, so the handle's spelling can be the subject's: a
  GenBank, EMBL or DDBJ accession alone, a name or a patent with a letter, more than 64 letters in
  a local or general id or an accession). NCBI's message is deterministic; LOSAT refuses it ("a
  PSSM query entry with a Seq-id twice in a search of -subject records ... is not supported")
  because its step runs when the PSSM is read, before the subjects, until it takes the subjects'
  Seq-ids. With `-db`, and for the Seq-ids whose message is the entry's own, LOSAT
  writes NCBI's message (above).

## Ported (rule 4)

- TBLASTX `-window_size 0` with a low `-threshold` (version 1.21, 2026-11-14, line RG item X2b):
  the one-hit word finder's ungapped HSPs of length 0 or less (`aa_ungapped.c:1054,1083`). NCBI
  keeps their `Int4` ends, sorts them by score, links them with the other HSPs in the first
  `BLAST_LinkHsps` (where one whose score is above the large-gap cutoff, or with small gaps on
  above the small-gap cutoff, can join a chain and so change the chain's E-value and the list's
  order) and deletes them in its re-evaluation (a
  length of 0 or less scores 0, below every cutoff), so they are never reported. The X2 build
  refused these searches; LOSAT now gives NCBI's output (the rule-3 merge case above aside).
- BLASTX `-ungapped -window_size 0` with a low `-threshold` (version 1.33, 2026-11-15, line RG item
  X2g): the same HSPs of the one-hit word finder. NCBI links them with the others in the first
  `BLAST_LinkHsps` and reaps them by their E-values; `Blast_HSPListReevaluateUngapped` does not
  rescore them (the subject is a protein) and deletes those of negative length in `Blast_HSPTest`
  (`blast_hits.c:993-1000,2719`: the alignment length is below `min_hit_length`, which is 0), so
  they are never reported. `Blast_HSPListSubjectBestHit` (`-subject_besthit`) sees them too
  (`blast_hits.c:2569`). LOSAT panicked on these searches and now gives NCBI's output. No merge
  divides by a sum of lengths of 0 here: NCBI never splits the query of an ungapped search
  (`split_query_cxx.cpp:60-61`), and LOSAT refuses subjects that NCBI would search in chunks.
  Those of length 0 are kept and reported; their Std-seg rows end one before they start
  (`x_UngappedHSPToStdSeg`). LOSAT now writes them as NCBI's formatters do. With `-db` the
  Seq-aligns stay as made: tabular rows of empty ranges (a subject whose only HSP it is shows the
  query's Seq-align id and fields at 0: `x_SetQueryCovSeqalign` throws on the empty row and
  `SetFields` leaves before its alignment block), Seq-aligns whose empty rows end at `TSeqPos`
  4294967295, `CSeq_align::GetSeqRange`'s `Row is empty` (`Seq_align.cpp:227`, exit status 255)
  where `qcovs` and outfmt 0's description table read the rows, and a block without rows past
  that table, and the same exception where outfmt 0's `-sorthits` and `-sorthsps` sorts read them
  (`SortSeqalignForSortableFormat`: `-sorthits 4`, `-sorthsps 2` and `4`). With `-subject` NCBI's
  Seq-loc mapper rewrites the query row to one nucleotide with a fuzz
  (`seq_align_mapper_base.cpp:1538-1559`): `CreateDensegFromStdseg`'s exception
  (`Seq_align.cpp:801`) where the tabular rows, outfmt 0's alignments and sorts (`-sorthits 3`,
  `-sorthsps 3`; `-sorthsps 4` gives 227 as with `-db`) and the XML make a Dense-seg of it, that
  nucleotide in `qcovs`. BLAST XML, XML2 and JSON catch the exception (exit status 0, the message
  on stderr, the query's report cut or left out), outfmt 13 and 14 do not (an empty side file,
  exit status 255). A subject whose every HSP has length 0 has a Seq-align without a row of
  letters: with `-db` NCBI's prefetch finds its subject row empty and writes `Error pre-fetching
  sequence data` (every format but 8, 9, 11 and 12; 8, 9 and 12 too with
  `PRE_FETCH_SEQS_LIMIT=2147483647`; the archive is never prefetched); LOSAT writes the same. With
  `-subject` the mapper finds no sequence type in such a subject and takes both rows as
  nucleotides; the stop of an empty row's segment, of the unknown length -1, is `start + m_Len - 1`
  in the unsigned `TSeqPos` (`seq_align_mapper_base.cpp:1067-1068`), which wraps only for an HSP
  at the first or second letter of the searched range (counted from a `-subject_loc` interval's
  start). Further in, the mapper keeps the segment's length -1 and its scores, the Seq-align stays
  as made, and `-subject` gives `-db`'s output, the prefetch's warning too; LOSAT writes it with
  `-db`'s port. Where it wraps, the mapper maps the subject row to the end of its interval: one HSP
  becomes an alignment without scores
  of the rest of the subject's letters against as many query nucleotides (from the start of the
  empty query row up, or on a minus strand ending two below it, as the mapper moves the other rows
  of a segment of unknown length; an E-value of 0.0, a bit score of -1.0, a score of -1; identities
  where a nucleotide letter equals an amino acid letter, positives by the matrix; past the query's
  last nucleotide the gap character; a minus-strand row that would start before the query's first
  nucleotide wraps, and NCBI then writes the caught row of the query's Seq-align id, `Row is empty`
  for `qcovs` or the Dense-seg's exception). LOSAT writes the tabular rows (outfmt 6, 7, 10, 20) of
  that alignment as NCBI does (decision X2g-7), reading the query record's letters outside a
  `-query_loc` interval where the row does (the cut keeps the longest subject interval's letters and
  one more on each side). LOSAT refuses the rest of `-subject`'s mapped alignments. Decision X2g-4:
  an HSP of length 0 at the first letter of the searched range of a subject with HSPs of letters, which
  the mapper cuts into two alignments without scores (the mapper knows such a subject as a protein
  and counts its start in nucleotides, three per letter, so its stop wraps at the first letter
  only). Decision X2g-5: a subject whose every HSP has length 0, one of them at the first or second
  letter of the searched range, in every format but the tables and outfmt 18, and in the tables
  where it has two or more (more rows, exceptions and out-of-memory stops). NCBI's
  output there is partly not defined (on the 4 GiB limit of the line's probes `BLAST ran out of
  memory`, exit status 4, for a segment of about 2^32 letters), and the rest comes from alignments
  that the mapper makes by accident.
- The end of a database protein's title ending in `]` that `x_AdjustProteinTitleSuffix` rewrites
  from the Bioseq's first tax id (with or without a taxonomy database), where a report shows the
  generated title (the outfmt 0 heading of a definition line without a title, outfmt 18, a tabular
  id that the title's first word replaces; version 1.38): ported by line RG item T5c on
  2026-11-14; version 1.10 listed it above as a deterministic result that LOSAT refused. Where the
  name that it writes makes the title one that `GenerateDefline` cannot decode (`chi, partial [ ]`
  with the taxid 100009 of `tests/seqdb/taxdb/t5`), the lookup fails as in rule 2's item of version
  1.37: `Unknown` in the tabular ids and the alignment skipped in the outfmt 0 heading (rows
  `db.t5.sfu.*`, `db.t5.sfux.*`).
- `-evalue` with a signed infinity or NaN (`+inf`, `-nan`, `+nan(1)`, `1e999`): NCBI's
  `CArg_Double` reads them with `strtod` and its option check (`<= 0`) lets them pass;
  NCBI searches as with the largest e-value (the cutoff stays 1, no HSP is reaped). LOSAT
  rejected them before Session S07+++b.
