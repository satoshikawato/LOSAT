# Product Decision: NCBI BLAST+ behaviour that is a defect

- Decision ID: `PD-LOSAT-NCBI-DEFECTS`
- Version: 1.78
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
  1.51 2026-11-18 (TBLASTX, the explicit rejection of the ungapped re-evaluation of HSPs on a subject
  of fewer than 12 bases, whose frame NCBI never translates and reads from address 1 on, rule 3, and
  names X-short-bat open in 1.46; line RG item X-short-tbx, Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.50, line DBID item pfoff-2);
  1.52 2026-11-18 (makeblastdb `-mask_data`, line NP item M4a-1c: the explicit rejection of a mask
  file's string whose length passes the end of the file and is 1 GiB to 2^47 bytes, which NCBI's
  reader reserves before it reads and which succeeds or fails by the host's memory, and of an
  element of a mask's list stated as `NULL` that makeblastdb dereferences, whose report holds the
  run's memory addresses, rule 3; Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07; numbered after 1.51, line RG item X-short-tbx);
  1.53 2026-11-18 (BLASTP `-use_sw_tback`, line RG item X4: the explicit rejection of 1.30 of a
  Smith-Waterman redo that NCBI repeats without end extended from RPSBLAST to BLASTP, rule 3,
  within the class of 1.30; it and BLASTP's rejections of a query chunk split again or without a
  query (decision D13) come after the reports of the query batches before, the form of 1.40 (the
  Owner's decision of 2026-11-16 on line RG's question Q-X2c); Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; numbered after 1.52, line NP item M4a-1c; the
  same day, X4 fix round 3: also the explicit rejection of a failed batch whose first query holds
  a redone match, on which NCBI crashes, rule 3; X4 fix round 4: a negative `CHUNK_SIZE` that
  leaves a chunk without a query also after the reports of the batches before; 2026-11-20, review
  3 of X4: so does the rejection that version 1.1 kept of the negative pair whose chunks leave
  gaps);
  1.54 2026-11-20 (BLASTP unified P-values, `-comp_based_stats 1u`, `2u` and `3u`, line RG item X5:
  the e-value adjustment with the values of the next query of the batch reproduced, rule 2;
  exception 16 for the last query of a batch, whose values NCBI reads past its array of query
  contexts: the Owner's decision of 2026-11-20 on question Q-X5 (「Q1=A」, Q1 of that day's batch of
  questions = Q-X5); the rule-2 rejection of a batch in which a query is followed by a query
  without a search space, Owner-delegated (standing instruction 2026-09-29, restated 2026-10-07)
  and part of the option A that the Owner chose; numbered after 1.53, line RG item X4);
  1.55 2026-11-20 (BLASTN, line RG item X6: reward 0 and penalty 0, NCBI's matrix-only scoring
  of `-task rmblastn`, now runs; the explicit rejection of it with the greedy extension of
  `-task megablast`, where NCBI divides by zero or does not finish, rule 3; Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.54, line RG item X5);
  1.56 2026-11-17 (RPSBLAST and RPSTBLASTN: a reported hit of a profile that an alias file's
  `FIRST_OID` and `LAST_OID` leave out, which NCBI's report does not find, refused in every format,
  within the classes of 1.10, rule 2, where NCBI writes `Unknown` and exits 0, and of 1.37, rule 3,
  where it stops with a null pointer exception: reported sequences only; and the explicit
  rejections, rule 3, of the concatenated length of an RPS search unit at which NCBI's interval tree
  splits a node without end, its diagonal table never stops growing, it frees the profiles twice,
  or an ungapped extension reads past the row after the profiles or that row past the `.rps`
  file's last memory page, precedent 1.30; line RG item D1a-4, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; fix round 1 of D1a-4: the reported profiles only,
  the refusals at the batch where NCBI's engine meets them, the row after the profiles read as
  NCBI reads it, DECISIONS D1a4-9; fix round 2: following profiles of one Seq-id count as one
  subject in the cut, and a reported profile out of the range is written as the included copy of
  its record that NCBI's report finds by its Seq-id; fix round 3: that copy has the same letters
  too, else it is refused; numbered after 1.46 at the merge of dev, then after 1.48 and
  1.49 at fix round 4's two merges of dev, then after 1.50 at fix round 7's merge of dev, then
  after 1.51 and after 1.52 at the two merges of dev before the pull request, then after 1.53, 1.54
  and 1.55 at the last three merges of dev);
  1.57 2026-11-17 (approved exception 17, the Owner's decision of 2026-11-17 on line RG's question
  D1a-3 Q1, 「A 正しい数を書く (推奨)」: a database length or `-dbsize` of exactly -2^63, which NCBI's
  serial writer prints as `-0`, written as `-9223372036854775808`; the rule-3 refusal of that one
  value removed; numbered after the D1a-4 entry at each of its merges of dev; exception 16 until
  dev's exception 16 of line RG item X5, version 1.54, came first);
  1.58 2026-11-20 (psiblast and tblastn `-in_pssm`, line RG item Q-BM: the BitMagic data of a
  Seq-table of the PSSM query's annotations decoded as NCBI's object manager decodes it, with
  NCBI's exception `BM::Invalid serialization format`; the rule-2 rejection of well-formed BitMagic
  data removed and the rule-3 rejection narrowed to the streams that NCBI reads past, runs with
  undefined behaviour or a write past BitMagic's blocks or buffers, or decodes from memory it did
  not write; porting the data and refusing its broken streams decided by the Owner on 2026-11-16
  (「推奨通り」, Q-BM). A GAP block that is not one is decoded, and a row of the index that reaches
  it answered, where NCBI's scalar, SSE2 and SSE4.2 builds of BitMagic agree, and refused under
  rule 3 only where they differ, as NCBI's result then depends on the machine it was built for
  (Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07); numbered after 1.57, line
  RG's approved exception 17);
  1.59 2026-11-21 (BLASTP `-use_sw_tback`, line RG item X4-excl: a match whose subject NCBI's
  fetch excludes sets the batch's status back to 0 unless NCBI's early termination skips it, so a
  failed Smith-Waterman redo before it neither fails the batch nor crashes NCBI, reproduced, rule 2;
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.58);
  1.60 2026-11-21 (`losat makembindex`, line RG item D5a: the explicit rejections, rule 3, of a
  `-legacy false` index whose `Save` NCBI never finishes (a stride above `-ws_hint` − `-nmer` + 1)
  and of a stride of 32768 or more, whose offsets NCBI computes by shifting an `int` past its
  width; the explicit rejection, rule 2, of an index volume of more than 4 GiB of packed
  sequences, whose sizes NCBI writes wrapped at 32 bits; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07);
  1.61 2026-11-21 (`losat makembindex`, line RG item D5a-t: the explicit rejection, rule 3, of a
  database sequence without letters, whose null Bioseq NCBI dereferences; Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07);
  1.62 2026-11-21 (makeblastdb, an index file that a failed write left empty or cut short: the
  outcomes NCBI reaches through its mapping of the file reproduced, rule 2, and the explicit
  rejections, rule 3, of an index file of no bytes, which NCBI maps as a null pointer and crashes
  on, and of a header field that reaches past the file's last page, where NCBI's result depends
  on its memory map; line LIB-O item O1e, PR #1433; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.61);
  1.63 2026-11-21 (alias files with `OIDLIST`, line RG item D1b-1: the explicit rejections of an empty
  OID list file, of bits past the memory page after it and of a file of its header alone, rule 3, and
  the two RPSBLAST and RPSTBLASTN rejections of D1a-4 extended to the profiles `OIDLIST` leaves out
  and to the totals of its scan; the explicit rejection of an alias file's `FIRST_OID`-`LAST_OID`
  range outside the bit array of the ranges below it, where NCBI clears bits outside that array
  (before it a crash, SIGSEGV, probes D1a_P and D1a_R of review A of D1b-1; past it writes over
  NCBI's heap, probes RB9 of review B), rule 3; Owner-delegated, standing instruction 2026-09-29,
  restated 2026-10-07);
  1.64 2026-11-21 (approved exception 18, the numbers of a `secondary-accession` modifier past
  `int` in psiblast `-in_msa`, line PS item PS6b, PR #1440; proposed under the standing instruction
  2026-09-29, restated 2026-10-07, on the Owner's rules of 2026-10-10 「NCBI のバグを再現する必要はない」
  and 2026-10-22 「広い型になおしましょう」, and accepted by the maintainer on 2026-11-21 by marking
  PR #1440 ready for review, as the PR asked; numbered after 1.63);
  1.65 2026-11-21 (`losat blastn -use_index true`, line RG item D5b-1: the explicit rejections,
  rule 3, of an empty first MegaBLAST index volume, whose header NCBI reads through a null
  pointer, of a new-style index volume whose sequence count NCBI reads from words its stream never
  set, and of an `-index_name` range of volumes whose loop NCBI never ends; Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07);
  1.66 2026-11-21 (approved exception 19, hit list sizes past `Int4`, NCBI parity line O1c, PR
  #1448: the preliminary hit list size, and the best-hit writer's twice that size, computed without
  NCBI's `Int4` overflow, which ends the rule-2
  reproduction of its wrapped sizes and the rule-3 rejection of the sizes on which NCBI crashes;
  proposed under the standing instruction 2026-09-29, restated 2026-10-07, on the Owner's rules of
  2026-10-10 「NCBI のバグを再現する必要はない」 and 2026-10-22 「広い型になおしましょう」, and accepted by the
  maintainer on 2026-11-21 by marking PR #1448 ready for review, as the PR asked; numbered after 1.65);
  1.67 2026-11-20 (the approved difference of `losat rpsbproc`'s sorts: Rust's standard sort, so hits
  equal in E-value and bit score and data rows of one key keep their input order, and a NaN
  E-value or bit score sorts after every number; the explicit rejection of 1.17 of hits whose NaN
  values would make NCBI's sort read outside its vector removed; decided by the Owner on
  2026-11-20, line LIB `DECISIONS.md` 33; the place of a NaN Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; line LIB-S item S1; numbered after 1.66; first 1.54, then 1.56, 1.58, 1.62 and 1.66 at
  the merges of dev that brought line RG item X5's 1.54, then line RG item X6's 1.55 and line RG
  item D1a-4's 1.56 and 1.57, then line RG's 1.58 to 1.61, then the 1.62 to 1.65 of lines LIB-O,
  RG and PS, then line NP's 1.66);
  1.68 2026-11-21 (approved exception 20, TBLASTX's linking of a query batch whose average query
  length is 0 against one subject, where NCBI crashes: LOSAT's sums in 64 bits, whose rows no run
  that NCBI completes gives; decided by the Owner on 2026-11-21 (「O2b-Q1=A」, line NP item O2b's
  question Q1, line NP `DECISIONS.md` 42); line NP item O2b; numbered after 1.67);
  1.69 2026-11-21 (approved exception 21, the Owner's answers of 2026-11-21 to questions Q1 and Q2
  of line NP's M4b, 「すべて推奨通りで」: with `-db_hard_mask`, BLASTP, BLASTX, TBLASTN and TBLASTX
  search a sequence whose first mask runs to its end (in TBLASTN and TBLASTX, a reading frame whose
  first part runs into such a mask) over the part before the mask, as NCBI's normal chunks search the
  near input, the masked letters staying masked; the refusals of those shapes removed; line NP item
  M4b-3; numbered after 1.68);
  1.70 2026-11-21 (alias files with `OID_MASK_TYPE`, line RG item D1b-2: the explicit rejection, rule 3,
  of a Seq-id on which NCBI's accession guide posts a warning under the key, which NCBI posts again at
  each read of the header that its cache does not answer; the rejection of 1.61 extended to an old
  style `losat makembindex` index of a sequence whose every definition line the key erases;
  dustmasker's `CUnassignedMember` error at such a sequence reproduced, rule 2; the rejection of a
  definition line with bytes NCBI reads as '#' extended to the lines that the key or a list removes,
  as NCBI's decoder writes its error for them before they are removed; the rejection of 1.56
  extended to a reported RPS profile of a database of several volumes whose lines the key removes,
  in the classes of 1.10 and 1.37; Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07; numbered after 1.69; first 1.65, then 1.66, 1.68 and 1.69 at the merges of dev that brought
  line RG item D5b-1's 1.65, then line NP's 1.66 and line LIB-S's 1.67, then line NP's 1.68, then
  line NP's 1.69)
  1.71 2026-11-22 (outfmt 0 to 4 of such a sequence are written as NCBI writes them: its
  description line is the Seq-align id and `Unknown`; outfmt 0 stops with `BLAST Database error`
  where the first alignment shown is of such a sequence under a Seq-align id that is not local,
  else writes the skipped line for each of its HSPs; outfmt 1 to 4 stop where any alignment shown
  is; the description lines load the sequences that they show, in psiblast's rounds after the
  first those of the two-part table of the first `-num_descriptions` HSPs (at least the first:
  `-num_descriptions 0` shows its sequence too; psiblast ignores `-sorthits` and `-sorthsps`, as
  NCBI's does); a report that stops the run ends it there, and no later batch is prefetched; the
  refusal of 1.31 no longer covers outfmt 0 to 4, but in a call of LOSAT's library that writes
  them with a format of another writer it covers the whole call, also where the call searches
  once for each hit list size; line DBID item pfoff-3, Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07);
  1.72 2026-11-22 (`losat windowmasker`, line NP items M5a and M5b: NCBI's defined results of
  `-mk_counts`, `-convert` and `-checkdup` reproduced, rule 2, and the explicit rejections, rule 3,
  of the masking and counts inputs whose result NCBI leaves undefined, takes from the C runtime or
  the host's memory, or fetches from NCBI's servers; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.71; first 1.68, then 1.70 and 1.72 at the merges
  of dev that brought line NP items O2b's 1.68 and M4b-3's 1.69, then line RG item D1b-2's 1.70 and line
  DBID item pfoff-3's 1.71);
  1.73 2026-11-21 (TBLASTN `-use_sw_tback`, line RG item X4b-1: the explicit rejections of 1.53,
  a Smith-Waterman redo that NCBI repeats without end and a failed batch whose first query holds
  a redone match, on which NCBI's clean-up crashes, extended from BLASTP to TBLASTN, rule 3,
  within the class of 1.30, after the reports of the query batches before; a subject excluded
  at the redo's fetch clears the batch's status as in NCBI; `-in_pssm -use_sw_tback`, on which
  NCBI fails, crashes, does not end or writes rows that change from run to run, refused at the
  first match the redo reaches instead of before the search; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; numbered after 1.55, line RG item X6, then after
  1.56 and 1.57, line RG item D1a-4, after 1.58, line RG item Q-BM, and after 1.59 to 1.72
  at the merges of dev);
  1.74 2026-11-22 (`losat blastn -use_index true` with a MegaBLAST index of only some database
  volumes, line RG item D5b-3: the results that NCBI keeps of an index volume before a volume
  without an index, read by later query batches and chunks, reproduced, rule 2, at every thread
  count as NCBI's one thread keeps them (the Owner's rule of 2026-10-13); the explicit rejection,
  rule 3, of a later search whose scan saves initial hits into the list that NCBI then frees;
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.73; first
  1.72, then 1.73 and 1.74 at the merges of dev that brought line NP items M5a and M5b's 1.72, then
  line RG item X4b-1's 1.73);
  1.75 2026-11-22 (approved exception 22, the Owner's answer of 2026-11-22 「Q8=B」 to line NP item M5's
  question Q10: `losat windowmasker` posts NCBI's `bad id` error for each line of an `-ids` or
  `-exclude_ids` list matched as text that gives an empty id, and skips the id, where NCBI crashes; the
  explicit rejection of 1.72 removed; line NP item M5; numbered after 1.74; first 1.73, then 1.74 and 1.75
  at the merges of dev that brought line RG item X4b-1's 1.73, then line RG item D5b-3's 1.74);
  1.76 2026-11-21 (BLASTX `-use_sw_tback`, line RG item X4b-2: the explicit rejections of 1.53
  and 1.73, a Smith-Waterman redo that NCBI repeats without end and a failed batch whose first
  query holds a redone match, extended to BLASTX, rule 3, within the class of 1.30, after the
  reports of the query batches before; approved exception 23, the Owner's decision of 2026-11-21
  on line RG's question Q-HEAP, 「Q-HEAP=B」: the heap test of NCBI's redo in e-value mode that
  reads past its array of the batch's heaps tests an empty heap (NCBI's run as is with an empty
  heap's answer there, which NCBI's run of the batch followed by queries without a hit checks
  where it searches the batch alike); a new explicit rejection, rule 3: `-ungapped
  -comp_based_stats 0 -use_sw_tback`, on which NCBI crashes at every query batch it searches;
  otherwise Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after
  1.73, line RG item X4b-1, after 1.74, line RG item D5b-3, and after 1.75, line NP item M5, at
  the merges of dev, and after exception 22);
  1.77 2026-11-22 (psiblast and tblastn `-in_pssm` of a PSSM in ASN.1 binary, line RG item D3c-3a:
  the explicit rejections, rule 3, of a binary PSSM's string or byte block whose length passes the
  end of the file and is 1 GiB to 2^47 bytes, the class of 1.52, and of a query whose letters are
  missing or whose number is not its Seq-inst's `length`, on which NCBI crashes or writes past its
  buffer for most lengths and which it answers for one less than the letters, an answer LOSAT does
  not yet give; and a refusal that is not a
  rule-3 item, of a PSSM's Karlin-Altschul values that are not finite numbers where LOSAT does not
  yet compute NCBI's output of them; Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07; numbered after 1.76; first 1.71, then 1.72, 1.74, 1.76 and 1.77 at the merges of dev that brought line DBID item
  pfoff-3's 1.71, line NP items M5a and M5b's 1.72 and line RG item X4b-1's 1.73, then line RG item D5b-3's 1.74 and
  line NP item M5's 1.75, then line RG item X4b-2's 1.76);
  1.78 2026-11-22 (approved exception 24, BLASTP and BLASTX `-db_hard_mask` with overlapping masks
  after which NCBI's word hits reach the table of diagonals of the sequences searched before or
  after a sequence: each such sequence searched on its own, the Owner's answer Q5 = B of M4b of
  2026-11-22 22:34 JST, in place of the explicit rejection `StaleDiagonals` of this version's first
  draft, at NCBI's one-thread offset on a new table where NCBI's one-thread search meets the
  entries that no sequence wrote, decision 55 of line NP's series M, 2026-11-23, Owner-delegated,
  and in BLASTP on NCBI's table of each query batch where its hits can meet entries that its own
  chunks did not write;
  BLASTP and BLASTX `-db_hard_mask` with masks that overlap, nest or are stored out
  of order: searched as NCBI searches them, rule 4, also the seed that NCBI's two-hit scan saves from
  the one word of a part between overlapping masks shorter than the word; the explicit rejections,
  rule 3, of such a part whose seed's alignment NCBI's interval tree splits without
  end (`short_chunk_tree_refusal`), of a part between masks shorter than the word whose one word NCBI
  reads past the sequence's end, only for a BLASTP query that holds residue 0, which no query read
  from FASTA does (`WordPastEnd`), of a mask that runs to the end after the part from the start with
  masks stored after it (`TailMaskWithLaterMasks`, within exception 21's class) and of a mask that
  ends before it starts (`Inverted`); not yet ported, as the range longer than 5,000,000 letters:
  masks out of order whose parts NCBI merges across an overlap (`PartsOutOfOrder`) and such a short
  part whose alignments meet another part's in that merge (`short_chunk_merge_refusal`), line NP item
  M4b-3d; TBLASTN and TBLASTX unchanged; line NP item M4b-3b, the Owner's decision of 2026-11-21
  (「すべて推奨通りで」, M4b, overlap = A); Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07, for the rejections; numbered after 1.77, line RG item D3c-3a, as the next free numbers
  on dev when it merged dev (agreed with line NP item M5 on 2026-11-23); first 1.72, then 1.73 at the
  merge of dev that brought line NP items M5a and M5b's 1.72, then 1.74 after dev's 1.73, then 1.76 and
  exception 23 after dev's 1.74, then 1.78 and exception 24 after dev's 1.75 to 1.77 of line NP item M5
  (exception 22), line RG item X4b-2 (exception 23) and line RG item D3c-3a)
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
  kept the rejection for these programs. BLASTP's rejections of a chunk split again or of a
  chunk without a query come after the reports of the query batches before, where NCBI's
  exception stops the search of the batch (version 1.53, the form of 1.40; psiblast, which
  searches one query at a time, stops in the query's round 1 after the rounds of the queries
  before, as NCBI does: review 3 of X4 found that NCBI splits that round too). A chunk without a
  query is refused first: NCBI's `CQuerySplitter::Split` makes every chunk's query factory before
  it sets a chunk up and stops at such a chunk's missing query vector
  (`split_query_cxx.cpp:234,876-879`, `CCoreException`), before any chunk is split again (as with
  a positive pair whose evened chunk size falls below the overlap, whose chunk 0 covers the
  batch). A negative `CHUNK_SIZE` is one of these where its wrapped evened size passes the
  batch's length (`split_query_cxx.cpp:145-171`): chunk 0 covers the batch and the chunks after
  it hold no query, so NCBI stops at the split after the batches before (version 1.53, X4 fix
  round 4; review 2 of X4 found NCBI's exit status 3 with the reports of the batches before,
  rows `db.x4.chunk.neg.*`). A negative pair whose wrapped ranges leave gaps between chunks
  that all hold a query and are not split again is not a failure of NCBI's: NCBI searches on
  with the hits in the gaps lost (exit status 0), and the maintainer kept that setting rejected
  (version 1.1, rule 2). BLASTP's rejection comes where the batch's search would start, after
  the reports of the query batches before (version 1.53, review 3 of X4: the form of 1.40 for a
  rejection in the middle of a run; row `db.x4.chunk.neg.gap.6`, whose bytes are NCBI's report
  of the batch before, row `db.x4.chunk.gap_anchor.6`; rows `db.x4.chunk.empty.*` and
  `db.x4.chunk.again.*` hold NCBI's bytes of the other two stops).
- Settings on which NCBI crashes (decision D3): an IDENTITY word size of 5 or 6,
  `-use_sw_tback -ungapped -comp_based_stats 0`, compressed-lookup thresholds that
  exhaust its overflow bank. (`-max_target_seqs` values that make NCBI's preliminary hit
  list size not positive were rejected before version 1.66, approved exception 19.) A query length plus `-window_size` in (2^30, 2^31 - 1], where
  NCBI does not finish (decision D8).

Kept as explicit rejections although NCBI's result is deterministic (rule 2's
maintainer clause; decision D11): a query length plus `-window_size` above 2^31 - 1, where
NCBI's `Int4` diagonal table wraps to one cell and no hit is found (BLASTP, TBLASTX;
together with D8, every sum above 2^30 is rejected).

Reproduced as NCBI (rule 2): the out-of-range `double` to `Int4` conversions of
`-threshold` and the X-drops (INT_MIN, decision D4); `-max_target_seqs` of 2147483624 or
more in BLASTP, whose preliminary hit list size wraps to 2 to 48, before version 1.66 (approved
exception 19 computes it without the overflow). TBLASTX's culling limit
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

## Approved exception 16 (BLASTP unified P-values, the last query of a query batch, version 1.54)

Version 1.54, line RG item X5 (line RG `DECISIONS.md` X5-1 to X5-4). The Owner decided on 2026-11-20
(words 「Q1=A Q2=B Q3=A Q4=A Q5=A Q6=B Q7=A Q8=OK」, whose Q1 was this question) question Q-X5: option A,
running unified P-values with NCBI's bytes for a query that has a next query in its batch (rule 2
below), the query's own values for the last query of a batch (this exception) and the refusal of a
batch whose next query has no search space (rule 2's rejections). Until this version LOSAT refused every request for
unified P-values (NCBI parity line `DECISIONS.md` 23 a, 2026-10-25; `OWNER_DECISIONS.md` row 23).
Not taken: keeping that refusal for every query (a search of one query is the commonest use, and
all but the last query of a batch have NCBI's bytes, rule 2 below); the query's own values for
every query (NCBI's bytes for a query with a next query are deterministic, rule 2).

16. **Unified P-values: the last query of a query batch.** A `u` after the mode of
    `-comp_based_stats` (`1u`, `2u`, `3u`; BLASTP only, `blast_args.cpp:880-884`) makes NCBI's
    composition-based redo combine the P-value of each alignment with the P-value of the pair's
    compositions (`s_AdjustEvaluesForComposition`, `blast_kappa.c:134-210`). The step scales the
    e-value with the query's length, length adjustment and effective search space, which it reads
    from `queryInfo->contexts[context_index]` (`:147-158`). `Blast_RedoAlignmentCore_MT` sets
    `context_index = query_index * numFrames` (`blast_kappa.c:3543`), its loop over the query's
    frames adds one per frame (`:3590-3592`; BLASTP has one frame), and the value after the loop
    goes to `s_HitlistEvaluateAndPurge` (`:3675-3683`): the context of the next query of the batch.
    For the last query of a batch, and so for every search of one query, that is the element after
    the array's end (`calloc(last_context + 1, 32)`, `blast_query_info.c:77-79`): the three values
    are heap bytes. In our runs of NCBI BLAST+ 2.17.0 (the evidence below) the report of one query
    alone took 7 forms in 40 runs at 1 and 4 threads, and that of a query whose hit list fills 3
    forms in 10 runs. LOSAT adjusts the last query of a batch with that
    query's own three values, which the step's parameter is documented to hold ("info about this
    query context", `blast_kappa.c:126-127`), and redoes every match of it. That result is
    NCBI's: in NCBI's search of the batch's queries with a copy of the last query after it in the
    same batch, the next context of the query holds the copy's values, which are the query's own,
    and NCBI's report of the query is deterministic (the same bytes in every run, at 1, 4 and 8
    threads, and with a third query after the copy). LOSAT's report of the last query of a batch
    equals NCBI's report of that query in that search (outfmt 0, 5, 6 and 7 are frozen), where
    the copy stays in the query's batch (a batch closes at `BATCH_SIZE` letters, 10,000 by
    default, so the search with a copy may need a larger value) and NCBI searches that batch in
    the same chunks (a batch above `CHUNK_SIZE` is split, the copy moves the split, and the split
    can change the preliminary hits; the frozen rows search one chunk). With `-use_sw_tback` the
    failure of a batch's traceback follows the batch's last redone match (version 1.53, X4): in
    NCBI's search with a copy that match is the copy's, which NCBI may end early where the copy's
    hit list fills; where it does not, LOSAT's batch fails where NCBI's search with the copy does.
    The redo's early termination (`BlastCompo_EarlyTermination`, `redo_alignment.c:1560-1580`)
    belongs to the same result. NCBI redoes the matches in the order of their queries
    (`BlastHSPStreamRead`, `blast_hspstream.c:288-319`), and the test needs the heap of every
    query of the batch to hold `-max_target_seqs` matches with an e-value of at most 0.002
    (`inclusion_ethresh`, `blast_kappa.c:3074-3077`, `compo_heap.c:405-409`); the heaps of the
    queries after the one in work are still empty, so the test can end only matches of the
    batch's last query (those with a preliminary e-value above 0.01). In the search with a copy
    the query is not the last one and NCBI redoes every match of it. So LOSAT ends no match early
    in a search with unified P-values (`ResolvedBlastpArgs::redone_match_heap`); the queries
    before the last one are never ended early in NCBI either.
    Evidence: `docs/evidence/unified_p_last_query/` (NCBI's runs: the variation of the last
    query, the stable reports with a copy, LOSAT's reports against them). Fixture rows
    `db.x5u.*` of `LOSAT/tests/db_search_manifest/blastp.x5u.tsv`: the rows with `ncbi-copies:`
    hold NCBI's whole output of the search with copies, and the check derives the expected bytes
    from it by dropping the copies' reports (`LOSAT/tests/db_search_fixtures.py`). Unit tests
    `evalue_adjust_contexts_are_the_next_querys_of_the_batch_or_the_querys_own` and
    `test_adjust_evalues_for_composition_matches_ncbi_sequence_scaling`
    (`LOSAT/src/algorithm/blastp/kappa.rs`), `a_heap_without_early_end_is_never_filled_to_its_cutoff`
    (`LOSAT/src/algorithm/blastp/hsp.rs`).

## Approved exception 21 (BLASTP, BLASTX, TBLASTN and TBLASTX `-db_hard_mask`, a mask that runs to the end, version 1.69)

The Owner decided on 2026-11-21 (words 「すべて推奨通りで」, line NP `M/m4b/OWNER_QUESTIONS.md`) questions
Q1 and Q2 of M4b: option B of each, the search over the part before the mask as a normal chunk for
the four programs (Q1), and each HSP written once where NCBI's TBLASTX (and TBLASTN with
`-ungapped`) writes it once per seed past that part (Q2, the same fix in the same place). Not taken:
option A, refusing such databases (what LOSAT did from M4b PR 1 and PR 2 until this version): segmasker
masks reach the C-terminus of many proteins (189 of the 48,321 proteins of A. thaliana TAIR10 have
the shape with segmasker's defaults), so it refuses whole proteomes that NCBI often searches without
a crash. Its version is in the list at the top; elsewhere it is cited by its number and name,
"approved exception 21 (a protein or translated mask that runs to the end)". Line NP item M4b-3.

21. **`-db_hard_mask` in the protein and translated programs: a part from the start that a mask
    running to the end follows.** NCBI cuts each database sequence at its masks into chunks
    (`s_GetNextSubjectChunk`, `blast_engine.c:221-311`); the protein programs' chunks start at their
    range's left (no residual). TBLASTN and TBLASTX map each stretch between the masks to every
    reading frame by its positions divided by 3 (the minus frames from frame -1's length) and cut
    each frame as a protein subject (`blast_engine.c:805-830`). A sequence (or a frame) whose first
    chunk, from 0 and with at most 5,000,000 letters, is followed by a mask that runs to its end
    (masked whole or a masked tail; in BLASTP and BLASTX also after empty masks at 0; in a frame,
    masks at or near either end of the sequence) takes the "no chunking" branch (`blast_engine.c:264-268`), as in exception 13:
    its `seq_ranges` are the full range (`:162`) while `subject->length` is the chunk's. The protein
    scan reads the ranges (`masksubj.inl:43-59`, `blast_aascan.c:75`) and so the stored letters under
    the mask, while the right extensions (`aa_ungapped.c:841`), the diagonal offset
    (`Blast_ExtendWordExit`, `blast_extend.c:172`) and the interval tree of the gapped HSPs
    (`blast_gapalign.c:3830`) take the chunk's length. NCBI BLAST+ 2.17.0 then:
    - crashes on some such sequences: the preliminary search stops and frees the subject block, whose
      sequence points into the mapped volume (SIGABRT `munmap_chunk(): invalid pointer` or `free():
      invalid pointer`, or SIGSEGV; gdb: `BlastSequenceBlkFree <- BLAST_PreliminarySearchEngine`), at
      every thread count; after empty masks at 0 the chunk before follows an allocated chunk and the
      branch frees memory that NCBI did not allocate, in every run;
    - gives rows that differ from the near input's on others: hits lost, gained, or extended into
      the mask;
    - in TBLASTX (and TBLASTN with `-ungapped`), extends each seed read past the part to the left
      only and writes the same ungapped HSP, which ends where the mask starts, once per such seed.
    LOSAT scans that part as every other hard chunk, over its own letters (`[0, length)`; the chunk
    of a whole frame over the frame's, `[0, frame length)`), with the extensions, the diagonal state
    and the HSP tree of that one length: the chunks that NCBI searches in the near input, the same
    database with the mask shortened until a letter follows the part (one letter in BLASTP and
    BLASTX; in TBLASTN and TBLASTX the fewest bases after which every frame of the sequence keeps a
    letter after its part: 2 to 6 bases at an end in the fixtures). The masked letters stay masked:
    the traceback and the composition-based statistics read them as `X` (the masked bases as any
    base in the translated programs), as NCBI's do. So LOSAT's result is NCBI's result for the near
    input with the freed letters stored as `X` (`N`); the near input with those letters unmasked
    gives other rows where its alignments reach them (a tail mask of 1 or 2 letters) or under
    composition-based statistics (the default of BLASTP, BLASTX and TBLASTN). Where a mask of one or
    two bases at an end of the sequence leaves the frame no part after its first part (no near
    input keeps that part), the result is the part's alone, which NCBI's own run gives in the
    fixtures (its scan reads the masked letter, the extensions stop at the part's end). A mask that
    leaves a frame an empty first part (`FrameStartMask`: an empty first part of a frame, then a part
    that NCBI searches as the whole frame; in TBLASTN and TBLASTX this includes an empty mask at 0
    before a tail mask, which BLASTP and BLASTX search under Q1 as BLASTN does under exception 13)
    stays refused (rule 3) until the Owner's question Q4 of M4b is answered, as do a mask that ends
    past its sequence and in BLASTP and BLASTX a range longer than 5,000,000 letters. Masks out of
    order and overlapping masks stayed refused in this version; from version 1.78 BLASTP and BLASTX
    search them (rule 4, and approved exception 24) except the rejections of 1.78 (the end of the
    Explicit rejections list),
    and TBLASTN and TBLASTX still refuse them. `blastinput::db_mask::HardScan`,
    `HardMaskShape::of_scan`, `algorithm::common::translated_chunks::translated_db_mask_chunks`,
    `algorithm::blastp::blast_engine::db_masked::hard_chunks`.
    Evidence (NCBI BLAST+ 2.17.0, each run capped; line NP's task folder `M/m4b` (SURVEY.md §4.1,
    §4.2, NOTES_M4B2.md) and `M/m4b3` (`near/`):
    - BLASTP `-db_hard_mask 21` on `ps5` (p4, 120 GS letters masked whole by segmasker): SIGABRT in
      every run (3 of 3 at 1 thread, at 4 threads, cbs 0 to 3, blastp-fast, versions 4 and 5, parsed
      ids, volumes, alias); BLASTX on a hand mask of a whole sequence: SIGABRT `free(): invalid
      pointer`, the same stack. BLASTP and BLASTX on `pe_e0tail` (an empty mask at 0, then p5's
      [300, 400)): SIGABRT `munmap_chunk(): invalid pointer` at 1 and 4 threads.
    - Output without a crash: p5's tail [300, 400) against its near input [300, 399): the query
      p5[250, 400) finds 1-50/251-300 (score 272) in the near input and nothing with the tail mask
      (BLASTP, likewise BLASTX and TBLASTN). A corpus of 300 random subjects, 100 with one tail mask,
      and 40 queries across the mask's start: no crash, the same bytes at 1 and 4 threads and across
      runs, and the tail mask's rows differ from the near input's in 48 of 1,577 query-subject pairs
      of BLASTP (cbs 0), 5 of 157 of TBLASTN and 6 of 435 of TBLASTX.
    - TBLASTX duplicates: the corpus' 1,749 rows hold 579 distinct ones (one row 25 times); the hand
      tail mask 81 rows of 59 distinct (one row 23 times); TBLASTN `-ungapped` 52 of 31 distinct; the
      same bytes at 1 and 4 threads and across runs. In the fixtures, `th2`'s tail [1000, 1200) and its
      mask ending 2 bases before the end give 26,168 and 27,632 bytes where the near input
      [1000, 1196) gives 24,848.
    - The fixtures where NCBI ends normally and gives the near input's bytes (BLASTX `ps5`; TBLASTN
      and TBLASTX `tf4`'s t5 masked whole, `th2`'s head [0, 100);
      TBLASTN `th2`'s tail, `te`'s masks of the last three bases of subjects of 11 to 17 bases): LOSAT's
      rows are NCBI's own there.
    Fixture rows `db.pmask.*.hard.tail.*`, `db.tmask.*.hard.tf4.11.6`, `db.tmask.*.hard.th2.1[234].6`
    and `db.tmask.*.hard.te.{16,17,18,19,33,34,35}.6` of `LOSAT/tests/db_search_manifest/{blastp,blastx}.pmask.tsv`
    and `{tblastn,tblastx}.tmask.tsv`: NCBI's bytes, or with `ncbi-near:` those of NCBI's run of the
    near input (the rows `db.pmask.*.hard.near.*`, `db.tmask.tx.hard.near.th2tn.11.6`; databases
    `ps5n`, `pe_tailn`, `pe_e0tailn`, `th2tn`, built from copies of the FASTA files with the freed
    letters stored as `X` or `N`: `prot_p4x.faa`, `prot_p5x.faa`, `nucl_t4n.fna`). With the near
    input's own letters BLASTP and BLASTX give other E-values at the default composition-based
    statistics (`pe_tail`), and a 60-letter protein whose last 2 letters are masked aligns 1-58 with
    the tail mask (NCBI's own run, LOSAT the same) but 1-60 with the near input's mask of one letter.
    TBLASTN with `-window_size 0` on `th2`'s tail (`db.tmask.tn.hard.th2.13.w0.6`): NCBI's own
    bytes, the same as its run of the near input (with `-ungapped -comp_based_stats 0` too, 4,951
    bytes, no row repeated; LOSAT refuses TBLASTN `-ungapped` for every database, so no row);
    BLASTP with `-window_size 0` on `pe_tail` (`db.pmask.bp.hard.tail.pe_tail.w0.6`): the near
    input's. Unit tests `protein_shapes_of_the_hard_masks` and
    `translated_hard_masks_are_refused_where_a_frame_restarts_at_its_start`
    (`LOSAT/src/blastinput/db_mask_tests.rs`),
    `a_whole_frame_scans_to_the_frame_s_end_and_a_tail_masked_one_its_chunk`
    (`LOSAT/src/algorithm/common/translated_chunks_tests.rs`).

## Approved exception 24 (BLASTP and BLASTX `-db_hard_mask`, diagonal entries that reach another sequence, version 1.78)

The Owner decided on 2026-11-22 at 22:34 JST (words 「Q1=B Q2=B Q3=A Q4=A Q5=C Q6=B Q7=B Q8=B Q9=A
Q10=A Q11=A Q12=A Q13=OK」 to the question batch of that day, whose Q7 is line NP's question Q5 of
M4b, `M/m4b/OWNER_QUESTIONS.md`) option B: 「配列ごとに独立に探す（前の配列の印が届かない所まで
目盛りを進めます。今の判定の式を置き換えるだけで、速さは変わりません）」, with the merit
「こうした DB も探せて、スレッド数によらず同じ結果になります」 and the demerit 「前の配列の印が次の配列に
届く場合だけ、NCBI の 1 スレッドと bytes が違います…新しい例外になります」. Not taken: option A,
refusing a database with such a sequence before the search (the explicit rejection `StaleDiagonals`
of the first draft of version 1.78), which refuses whole seg-masked proteomes; option C, writing
NCBI's 1-thread rows at every thread count, which carries the entries of one sequence into the
next on purpose and makes the search depend on the order of the sequences. Its version is in the
list at the top; elsewhere it is cited by its number and name, "approved exception 24 (diagonal
entries that reach another sequence)". Line NP item M4b-3b.

23. **`-db_hard_mask` in BLASTP and BLASTX: overlapping masks after which NCBI's word hits reach
    the table of diagonals of another sequence.** NCBI's protein scan keeps one table of diagonals
    per thread (made at `blast_extend.c:40-66`, cleared at `:87-107`) and moves its offset after
    each chunk of a sequence by the chunk's length plus the window (`Blast_ExtendWordExit`,
    `blast_extend.c:162-173`, which clears the table once the offset reaches INT4_MAX / 4). A chunk
    between overlapping masks has a negative length (`blast_engine.c:221-276`), so the offset moves
    back and a later chunk's hits (`aa_ungapped.c:519-551` two-hit, `:771-802` one-hit) meet the
    entries of earlier chunks: of the same sequence (searched as NCBI searches them, rule 4) or,
    where the walk of the chunks reaches below the sequence's start or leaves entries past the next
    sequence's start, of another sequence. Which sequence a thread searched before depends on how
    the threads share out the database, so NCBI's rows for such a sequence depend on its neighbours
    at 1 thread and change from run to run with several threads.
    LOSAT searches each such sequence with every entry of the other sequences out of reach, at
    NCBI's one-thread place on the table where that place changes a result. Before the search it
    walks NCBI's one-thread offsets over the searched sequences in the search's order (the
    sequences that the lists and the alias files' ranges include, in the order of their ids; NCBI
    does not search one of no letters, `blast_engine.c:1429-1432`; each chunk's exit moves the
    offset or, at INT4_MAX / 4, clears the table: `algorithm::common::hard_chunk_walk::NcbiWalk`,
    with `-db_hard_mask` only, in the pass that looks for LOSAT's refusals). A sequence whose walk
    there meets the entries that no sequence wrote since NCBI's thread made or cleared its table (0
    at first, the window less than 0 after a clear: a hit of its lowest chunk that reads a word lies
    less than the window past them, or below them in the one-hit scan), or inside which NCBI's walk
    clears the table, starts on such a table at NCBI's one-thread offset (`TableStart::New`, only
    the searching thread's table): the first searched sequence where its walk goes below its start,
    a sequence searched within a few hundred letters of the database's start or of a clear, or one
    after walks that went back (`ov_lk`'s B in the one-hit scan, 98 below 0, with or without masks
    of its own). Every other sequence meets those entries nowhere at NCBI's offset. Before such a
    sequence the search moves the offset up by the walk's margin below its start, and after a
    sequence by the margin of its entries past the next sequence's start (`DiagonalWalk::start` and
    `DiagonalWalk::after`: the one walk that `algorithm::blastp::blast_engine::db_masked` and
    `algorithm::blastx::search_db_masked` read and apply to their tables, `TableStart::apply` and
    `TableEnd::apply`); a sequence with no such margin keeps its offset (NCBI's until the first
    move). Where a move would bring the offset to INT4_MAX / 4 or past it, or a chunk before a
    sequence's last would exit at INT4_MAX / 4 or past it, LOSAT clears the table first, so no
    clear falls inside a sequence that does not start on NCBI's table. From the window a walk
    reaches INT4_MAX / 4 (536,870,911) only for a protein of more than 536 million letters or about
    12 million masks at the window 40 (each chunk adds its length plus the window; the chunks with
    words are disjoint parts of the sequence); the longest walk of A. thaliana TAIR10 reaches 5,324
    positions (of the UniRef50 sample 77,483). The search without `-db_hard_mask` is unchanged. The
    kernels count every one word of a chunk past a sequence's end as a hit (as for a query with the
    gap letter), so they also move the offset for a word that no query read from FASTA hits, which
    changes no row. In BLASTX, which keeps a worker's table across subjects that are not consecutive,
    a sequence without masks that starts on a new table moves the offset after it past the entries
    that no sequence wrote as well as past its own (`DiagonalWalk::plain`), so the next sequence that
    the worker's table searches, not always the next one in the database, meets none of them; BLASTP
    searches such a sequence on a table of its own and leaves the worker's table as it was (review 3
    of M4b-3b, M1).
    NCBI searches each query batch of a query file alone (`blastp_app.cpp:216, 257`; 10,000 letters
    a batch in BLASTP, `BATCH_SIZE` overrides it): its table is the batch's own, of the power of two
    just at least the batch's length plus the window (`blast_extend.c:46-65`), and a hit's diagonal
    and the context test (unsigned, `aa_ungapped.c:566-573`) count the query offsets from the batch's
    start. BLASTP searches several batches at once (all queries on one table, or the blocks of whole
    batches of LOSAT_X_BPQBLOCK), which gives each batch NCBI's result where every entry within
    reach of a hit was written by the same chunk's scan. A sequence whose hits can meet other
    entries, one that starts on a new table or whose walk reads back (a chunk within reach of an
    entry that an earlier chunk can store), is searched on NCBI's table of each query batch instead
    (`algorithm::common::query_batch_tables::QueryBatchTables`, `DiagonalWalk::apart`: each batch's
    ring, the diagonal and the context test in the batch's offsets, one offset for all, by the
    reference scan), the second at the window plus its margin below its start, where no entry
    reaches its hits; the thread's table stays as the sequence before it left it (review 3 of
    M4b-3b, H1). The scans of the other sequences, and every search without `-db_hard_mask`, are
    unchanged. BLASTX searches each reader batch alone, as NCBI searches it with the default batch
    size.
    So the rows of each sequence are NCBI's one-thread rows of it, for the same query file and
    batches, in the same database with the letters of every other sequence replaced by letters that
    no query word hits (the lengths and masks, and so NCBI's offsets and statistics, stay the
    database's), the same at every thread count: they do not
    depend on the entries of other sequences or on `-num_threads`, and like NCBI's one-thread rows
    they can depend on the lengths of the sequences searched before it, through NCBI's offsets
    (`ov_zt` below: B's 7 rows after a sequence of 5 letters, 13 after one of 100, 6 as the first).
    They differ from NCBI's one-thread rows only where an entry of another sequence changes a
    result. Where NCBI's one-thread start of the sequence does not meet the entries that no sequence
    wrote, they are also NCBI's rows for the near input with a spacer of `X` before each sequence
    that such an entry reaches, at least as long as the walk's margin less the window (the margin
    past the end of the sequence before it plus the margin below the start of the next;
    `LOSAT/tests/make_dbmaskp_ov.py`, `spacer_length`, takes the database's average length where
    that is longer), under an alias with the database's NSEQ, LENGTH, STATS_NSEQ and STATS_TOTLEN
    (the statistics read those, `seqsrc_seqdb.cpp:160-205`).
    Not taken: a new table at the window before every sequence, whose entries a hit below the start
    meets where NCBI's one-thread search meets none, so it would differ from NCBI's 1-thread rows
    also where no entry of another sequence reaches; and the first draft of this exception's
    search, which started only the first searched sequence on NCBI's new table and moved every
    other past the entries that no sequence wrote as past another sequence's (review 2 of M4b-3b,
    M1): its rows differ from NCBI's one-thread rows where no entry of another sequence exists
    (`ov_zt`: 13 rows where NCBI's one-thread run writes 7). The start at NCBI's one-thread offset is
    decision 55 of line NP's series M (2026-11-23, Owner-delegated, standing instruction
    2026-09-29, restated 2026-10-07), within the Owner's option, which moves the offset only "where
    the previous marks do not reach", and its demerit, the bytes differing from NCBI's 1-thread run
    only where an entry of another sequence reaches a result.
    Evidence (NCBI BLAST+ 2.17.0 under a 4 GB limit and 60 s; line NP's task folder `M/m4b3b`,
    `Q5B_PROGRESS.md`, `q5b/ncbi`, `q5b/ncbi_near.sh`):
    - `ov_lk` (A with masks [40,200) [60,238), then B with [230,235)): A's word hits reach 98
      positions past its end, 137 in the one-hit scan. Its near input `ov_lkn` (A, 240 `X`, then B,
      the masks the same, an alias with `ov_lk`'s statistics) keeps A's rows and gives the same bytes
      at 1 and 4 threads. NCBI's own runs of `ov_lk`: at the default the near input's rows; at
      `-evalue 1000` one more row of B in BLASTP (two in BLASTX) at 1 thread and the near input's
      rows at 4; at `-window_size 0` BLASTP misses B's identical row 1-33/8-40 (E 3.27e-23) and the
      row 6-13/86-93 at 1 and 4 threads, BLASTX misses B's two rows at 1 thread and has them at 4;
      with `-ungapped -comp_based_stats 0 -evalue 1000` six more rows of B in BLASTP and nine in
      BLASTX, where two of B's rows also take other E-values. At `-window_size 0` NCBI's one-thread search starts B
      98 below 0 (A's walk goes back by 140), so B's hits in its first 98 letters meet the entries
      that no sequence wrote and are skipped there too; the near input `ov_lkx` (A's letters all
      `X`, the lengths and masks the same) gives the same rows of B as `ov_lk` at 1 and 4 threads
      (BLASTP 1, BLASTX 10), so LOSAT, which starts B on such a table at -98, writes NCBI's
      one-thread bytes of `ov_lk` there (fix round 2 of M4b-3b).
    - `ov_atha` (40 TAIR10 proteins whose seg masks overlap): 6 of them reach another sequence in
      the one-hit scan (3 below their start, 3 past their end); its near input, a spacer of 357 `X`
      at each of the six boundaries, gives NCBI's own bytes of `ov_atha` in every run (BLASTP and
      BLASTX; the default, `-evalue 1000`, `-window_size 0`, `-ungapped`; 1 and 4 threads).
    - `lk400` (200 pairs of A then B, BLASTP `-window_size 0`, `M/m4b3b/NOTES_M4B3B.md`): at 1 thread
      two runs gave B's identical row in 0 of 200 pairs and A's in 200 of 200; at 4 threads three runs
      gave B's in 5, 6 and 3 of 200 and A's in 197, 196 and 200, every run another md5.
    - Counts with segmasker 2.17.0's defaults (each database in its order, `M/m4b3b/scripts/stale_stats.py`,
      the model of `DiagonalWalk::start` and `after`): A. thaliana TAIR10 118 of its 35,642 proteins
      with masks at the default window (each past its end), 578 with `-window_size 0` (318 below their
      start, 293 past their end), 118 with `-task blastp-short` and 116 at word 6; a UniRef50 sample
      (its first 60,000 sequences) 242 of its 59,484 with masks, 1,438 with `-window_size 0` (839 below
      their start, 736 past their end; 2 more than the first draft's `StaleDiagonals`, whose criterion
      did not count the entries of the sequences before a walk that goes back further than its margin
      below the start), 242 with `-task blastp-short` and 229 at word 6. Both databases are searched.
    - `ov_zt` (review 2 of M4b-3b, probe `zt`, `M/m4b3b/review/REVIEW_M4B3B_R2.md`; A of 5 `X`,
      which stores no entry, then B with masks [10,100) [20,110), whose one-hit walk goes 70 below
      its start; BLASTP `-window_size 0 -evalue 1000 -comp_based_stats 0`): NCBI's one-thread search
      starts B at 5, where its walk meets the entries that no sequence wrote: 7 rows, where the near
      input with a spacer of 300 `X` before B has 13 (B's 40.8-bit row 19-35/111-127, E 8.03e-11,
      and 5 more at B's letters 124-154 are missing); a 4-thread run gave 6 in review 2 (without
      the row 29-31/177-179) and 7 in the runs of fix round 2. With A of 100 `X` (`ov_ztl`)
      NCBI starts B at 100 and writes 13 rows, and as the first sequence (`ov_zto`, an alias of OID
      2) 6; BLASTX 27, 48 and 23. LOSAT starts B of `ov_zt` and `ov_zto` on a new table at NCBI's
      offset (5, 0) and moves past the entries in `ov_ztl`: NCBI's one-thread bytes in each.
    - `ov_th10` (review 3 of M4b-3b, probe `r3`, `M/m4b3b/review/REVIEW_M4B3B_R3.md`; A of 430 `X`
      with masks [10,300) [20,310), whose two-hit walk goes back, then B without masks, which NCBI's
      one-thread search starts on a new table at 10; BLASTP `-evalue 1000 -comp_based_stats 0`):
      q alone 5 rows; queries q1 (60 letters) then q in two batches (`BATCH_SIZE=50`) the same 5 rows
      of q, in one batch 4 (q's row 16-17/27-28, E 82, whose diagonal reads the new table's entries
      only where q is the first query of its batch, is gone). `ov_thm20` (B at -20): 7, 7 and 5;
      `ov_th10f` (another B): 7, 7 and 6; with A of 480 `X` (B at 60, no new table) 4, 4 and 4.
      BLASTX: q's 8 rows in one batch and in two (`BATCH_SIZE=50`). LOSAT: NCBI's one-thread
      bytes in each, with one table of all queries, LOSAT_X_BPQBLOCK's blocks and its default.
    - A. thaliana TAIR10 (48,321 proteins, segmasker's defaults, algorithm 21) for 55 queries, the
      proteins that such a protein's entries reach or that reach below their start (BLASTX: the first
      20, back-translated; ws, `M/m4b3b/q5b/tair`, `q5b/ws/q5b_threads.sh`): LOSAT wrote the same bytes at
      1, 4 and 16 threads in every case; at the default window LOSAT's bytes are NCBI's (8,716 BLASTP
      rows, 1,185 BLASTX rows); with `-window_size 0` NCBI gives the same bytes at 1 and 4 threads and
      LOSAT 3 more BLASTP rows (AT1G47340.1 and AT1G50880.1 against AT2G17310.1 from its first letters,
      E 0.012 and 0.001; AT2G27380.1 against itself at 750-761, E 1.1) and 1 more BLASTX row (the first
      of these), whose one-hit seeds NCBI skips at entries of the protein searched before; no row of
      NCBI's is missing. (Measured before fix rounds 2 and 3 of M4b-3b, whose start at NCBI's
      one-thread offset and NCBI's tables of each query batch can change only a sequence that NCBI's
      walk starts within reach of the entries that no sequence wrote, or, in a BLASTP search of
      several query batches, one whose walk reads back.)
    Fixture rows `db.pmask.{bp,bx}.hard.ov_lk.*` (NCBI's bytes at the default and at
    `-window_size 0`; with `ncbi-near:` those of NCBI's run of the near input, rows
    `db.pmask.{bp,bx}.hard.near.ov_lkn.*`, otherwise), `db.pmask.{bp,bx}.hard.near.ov_lkx.21.w0.6`,
    `db.pmask.{bp,bx}.hard.ov_zt.*`, `.ov_ztl.*`, `.ov_zto.*`, `.ov_th10.*`, `.ov_thm20.*`,
    `.ov_th10f.*` and `.ov_thm20f.*` (NCBI's one-thread bytes; the `.q1q.` rows of two queries in
    one batch and, `BATCH_SIZE=50`, in two, also with LOSAT_X_BPQBLOCK 0 and 1 in
    `LOSAT/tests/blastp_prelim_blocks.rs`) and
    `db.pmask.{bp,bx}.hard.ov_atha.21.w0.6` (NCBI's bytes) of
    `LOSAT/tests/db_search_manifest/{blastp,blastx}.pmask.tsv`, checked at 1 and 4 threads; the
    databases `ov_lkn`, `ov_lkx`, `ov_zt`, `ov_ztl`, the alias `ov_zto.pal`, `ov_th10`, `ov_thm20`,
    `ov_th10f` and `ov_thm20f` of `LOSAT/tests/db_search/dbmaskp/` (`make_dbmaskp_ov.py`: `near_db`,
    `blank_db`, `probe4`, `probe5`). Unit tests `each_sequence_is_searched_on_its_own`,
    `a_new_table_holds_entries_below_the_start`,
    `ncbi_walk_names_the_sequences_that_meet_unwritten_entries`, `ncbi_walk_follows_the_clears`,
    `hard_starts_by_ordinal`, `starts_and_ends_apply_to_the_table`,
    `a_plain_sequence_on_a_new_table_ends_past_the_unwritten_entries` and
    `sequences_whose_hits_meet_other_entries_are_kept_apart`
    (`LOSAT/src/algorithm/common/hard_chunk_walk_tests.rs`), the tests of
    `LOSAT/src/algorithm/common/query_batch_tables_tests.rs`, `hard_mask_starts_and_ends_move_the_table`
    (`LOSAT/src/algorithm/blastp/blast_engine_tests.rs`, `LOSAT/src/algorithm/blastx/seed_tests.rs`).

## Approved difference (`losat rpsbproc`, the order of tied hits, version 1.67)

Decided by the Owner on 2026-11-20 (line LIB `DECISIONS.md` 33; the Owner's words are in
[`OWNER_DECISIONS.md`](../decisions/OWNER_DECISIONS.md)): `losat rpsbproc` sorts with Rust's
standard sort. NCBI's order is deterministic for its binary, so rule 2 would otherwise reproduce
it. The place of a NaN (below) is Owner-delegated (2026-11-20, standing instruction 2026-09-29,
restated 2026-10-07).

rpsbproc 0.5.1 sorts with `std::sort` in five places whose order the output reads: the hits of a
query, twice, by `TDomAlignOrderByEValue` (the E-value ascending, then the bit score descending,
`common/cdalignproc_base.hpp:62-70`), before it picks the representative, standard and specific hits
(`common/cdalignproc_base.cpp:197-270`) and the structural domains (`:305-342`); and the rows of
its CDD data by PSSM-ID, by accession and by cluster (`CPtrMap::Reset`, `common/ptrmap.hpp:167-179`).
A sixth map, the clusters by PSSM-ID (`m_pssmid2fam.Reset`, `common/basedata.cpp:2585` and `:2688`),
is sorted but never searched, so its order does not reach the output; LOSAT does not build it.
The C++ standard does not keep the order of elements that compare equal in `std::sort`; NCBI's
Linux builds leave them in the order that their C++ library's sort gives. LOSAT sorts with
Rust's `slice::sort_unstable_by`, which sorts in place, by NCBI's comparisons made a total order
and then by the element's place in the input, which gives the order of a stable sort
(`rpsbproc::data::sort_ties_in_input_order`):

- Hits equal in E-value and bit score (the numbers compared as C++ compares doubles: -0.0 and 0.0
  are equal) keep the order of the archive or XML. Data rows of one key keep the order of their
  file; no CDD v3.22 file has two rows of one key.
- A hit whose E-value or bit score is NaN (printed `-nan`) sorts after every hit with a number in
  that key (Owner-delegated, 2026-11-20): a NaN E-value after every E-value that is a number, a NaN
  bit score last among the hits of one E-value; two NaNs are equal. NCBI's comparison orders a NaN
  against no number (no strict weak ordering), so its place for such a hit depends on the steps of its library's sort. Of the six
  fixtures whose sorts hold NaN keys, the place after every number leaves 2 equal to NCBI and 46
  lines different in the other 4, the place before every number also 2, with 59 lines.

Where it shows: in the DOMAINS lines of `-m std` and `-m full` within one query, tied hits change
order, and which of two tied hits becomes `Specific` can change (`-m std` prints the other hit's
line, `-m full` swaps their hit types). In a query with a NaN hit, the order of its hits, which hit
is representative, `Specific` or its superfamily's, and the SITES and MOTIFS lines that follow can
change. Line LIB measured the stable sort on 2026-11-20 (dev e851ce5a7, the same order of numbers
as here): 6 of the 408 fixture outputs and 18 of 74 real-size outputs (E. coli K-12 MG1655 proteins
with rpsblast, `-evalue` 0.01 and 10, and its genome in 929 pieces of 5 kb with rpstblastn, against
CDD v3.22, as archives and as XML) differed, by 1 to 6 lines each, every one from the first sort of
hits; in 7 of them one query (`Query_4017` of the `-evalue` 10 search) printed smart00370 in place
of smart00369 as `Specific` (both E-value 0.759289 and bit score 28.0864). No real-size output holds
a NaN. With the NaN rule, 9 of the 408 fixture outputs differ: `tie1.std`, `tie1.std.f`, `tie1.full`,
`tie1.full.f` and `x.tie1.doms.full` by 2 lines each (two tied hits change places), and
`x.h_tnan3.full`, `x.h_tnan3.both`, `x.h_ntnanthird.full` and `x.h_ntnanthird.both`, inputs with NaN
values made by hand (fixture README).

Fixtures: the rows of `use` `tied` in `LOSAT/tests/fixtures/rpsbproc/CASES.tsv` hold LOSAT's bytes
(`tied/`), and the test `tied_rows_differ_from_ncbi_only_in_the_order_of_tied_hits`
(`LOSAT/tests/rpsbproc_fixtures.rs`) compares them with NCBI's frozen bytes (`oracle/`): standard
error and exit status are equal, and the standard output differs only inside query blocks, either in
runs of DOMAINS lines that share one query, E-value and bit score with the same hit types, or in a
query with a NaN E-value or bit score, with the number of such blocks and lines that
`tied/CLASSES.tsv` lists. Code: `rpsbproc::domains::e_value_order`, the sorts of
`rpsbproc::data` (`CddData`'s lookups) and `rpsbproc::data::sort_ties_in_input_order`.

## Reproduced as NCBI (rule 2)

- BLASTP unified P-values, `-comp_based_stats 1u`, `2u` and `3u` (version 1.54, 2026-11-20, line RG
  item X5, DECISIONS X5-1): for a query that has a next query in its query batch, NCBI adjusts the
  e-values with the next query's length, length adjustment and effective search space
  (`blast_kappa.c:3543`, `:3590-3592`, `:3675-3683`; exception 16 above describes the index). The
  result is deterministic: NCBI's `-num_threads 1` writes the same bytes in every run, and they
  follow the next query (a next query of another length changes the e-values and the rows that
  pass `-evalue`). LOSAT reproduces them byte for byte from the batch's queries at every thread
  count (NCBI with more threads may form its batches another way, `-mt_mode 1`)
  (`kappa::evalue_adjust_contexts`); a next query that is not valid (no residue that can be scored)
  gives a length adjustment of 0 and the search space of `-searchsp`, as NCBI's
  `BLAST_CalcEffLengths` leaves them (`blast_setup.c:773-784`, `:846-847`). So the report of a
  query depends on the query after it and on `BATCH_SIZE`; only the last query of a batch is an
  exception (16). In PSIBLAST, BLASTX, TBLASTN, RPSBLAST and RPSTBLASTN NCBI reads the mode's
  first character and ignores the `u` (`blast_args.cpp:880-884`), and LOSAT does the same; the
  refusal list had named PSIBLAST by mistake.
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
- makeblastdb, an index file that a failed write left cut short (version 1.62, 2026-11-21, line LIB-O
  item O1e, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07). NCBI opens the
  database it wrote and reads the index through its mapping, where bytes the file lacks read as
  zeros to the end of its last page. LOSAT reproduces these outcomes: "Not a valid version 4 or 5
  database.", "requested sequence type does not match DB.", a header-less database and, for a date
  that was cut, the exception of CTime (exit status 255, no metadata file). The rule-3 refusals of
  an index file of no bytes and of a header field past the last page are in the list below.
- Defects without an effect on the output that LOSAT follows: the negative `Int4`
  offset of `BlastGetStartForGappedAlignmentNucl`, the out-of-range read of
  `hard_ranges[1]` for the last subject chunk, the out-of-range `double` to `Int4`
  conversion (INT_MIN), and the catch in `BlastSetupPreliminarySearchEx` that never
  matches.
- `losat blastn -use_index true` with a MegaBLAST index of only some database volumes (version
  1.74, 2026-11-22, line RG item D5b-3, Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07). NCBI keeps an index volume's seeds while a search thread may still need them and
  counts a volume down as a thread passes it, but a thread that goes on to a database volume
  without an index counts none down (`blast_dbindex.cpp:812-857`). So the last index volume that
  the search reaches before each such volume that it reaches (all of them, or those with a
  sequence in the database's OID list) keeps the seeds of the first search for the rest of the
  run, and a later query batch or query chunk takes them as the first search's word finder left
  them: every list moved out and empty (`blast_extend.c:248-259`). There it finds no seed: a
  sequence whose list the first search had is read and searched without seeds, the others are
  skipped, and the hits of the later queries in that volume are lost (fixture rows `m_part_b1`,
  `m_part_one` and `m_part_c2000` of `tests/blastn_mb_index`; NCBI's `-num_threads 1` writes the
  same bytes in every run). LOSAT reproduces this with NCBI's counts in one thread's order
  (`mb_index_seeds::plan_search`). A first search, and an index of every volume, are not
  affected. With more threads, which volumes keep their results depends on which thread passes
  which volume, and NCBI's result can differ from its one-thread result (with `-num_threads 4`,
  row `m_mixed1_t4` crashed, exit status 139, in 8 of 8 runs, where its `-num_threads 1` search
  completes); LOSAT gives the one-thread result at every thread count, by the Owner's rule of
  2026-10-13 (「DB-3: NCBI の 1 スレッド。LOSAT はどのスレッド数でも同じ出力」; rows `m_mixed1_t4`
  and `m_part_b1_t4` are compared with NCBI's `-num_threads 1` files).

- `losat windowmasker -mk_counts`, `-convert` and `-checkdup` (version 1.72, 2026-11-22, line NP
  item M5b, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07): NCBI's defined
  but surprising results are written byte for byte (fixture rows of `LOSAT/tests/windowmasker`):
  - `-smem` passes through `atoi` into a `Uint2` (`seq_masker_ostat_factory.cpp:63-72`): 65536 is
    0 megabytes, which gives the insufficient-memory error (exit status 3) after an empty output
    file; 70000 is 4464; -1 is 65535.
  - The `-t_*_pct` options apply only when `-out` names a file (`win_mask_gen_counts.cpp:216-219`).
  - Binary counts written to the standard output start with an extra 0 word
    (`seq_masker_ostat_bin.cpp:62-67`), which NCBI's own reader then rejects (`wrong length`).
  - An `-out` that cannot be opened gives exit status 0, no file and no message (NCBI's writer
    writes to a failed `ofstream`).
  - Counts above 511 in the optimized formats run into the check bits (`Uint2` counts, 9-bit
    fields, `seq_masker_ostat_opt.cpp:82-99, 311-325`); counts above 65,535 keep their low 16 bits.
  - `-t_high_pct 0` gives `t_high` 0 and every count 0 (the ascii file then fails to load, the
    obinary one masks nothing); `-t_high_pct 100` and `-t_low_pct 100` keep the counts NCBI keeps.
  - A negative `-genome_size` wraps as a `Uint8` and gives units of 16 letters.
  - The FASTA reader's warnings are written once per read of the input, so their number follows
    `-mem` (one read per prefix and pass).
  - `-convert` of an ascii counts file whose unit size line is 0 (or a multiple of 256) to `ascii`
    or `binary`; a `##pct:` line whose value is missing, or whose word does not start a number,
    reads as 0, as NCBI's `istringstream` stores it.
  - A `-checkdup` interval whose start falls below 0 is printed as NCBI's `size_t` wrap.

## Explicit rejections (rule 3)

- A scoring system without Karlin-Altschul values when an invalid query (such as all N)
  comes first in its batch: NCBI goes on without a gapped block and crashes (SIGSEGV).
- A reward of 32767 or a penalty of -32768 (after NCBI's 16-bit conversion):
  `BLAST_SCORE_MAX`/`BLAST_SCORE_MIN` fall outside NCBI's score range; a reward of 32767
  is counted past the end of the frequency array (invalid queries or a crash), a penalty
  of -32768 makes every query invalid. Version 1.1 (maintainer, 2026-10-02): a reward of
  32768 or more, which NCBI's 16-bit field wraps to 0 or below (NCBI: every query invalid,
  no hits, a deterministic result), stays an explicit rejection too.
- Subjects of 2^31 letters or more in total: NCBI's 32-bit total wraps.
- megablast gap costs above 32767: NCBI's 32-bit greedy distance arithmetic can wrap.
- Reward 0 and penalty 0 (matrix-only scoring, the scores of `-task rmblastn`) with the greedy
  extension of `-task megablast` (line RG X6, 2026-11-20, Owner-delegated): with gap costs 0
  and 0 (megablast's own) NCBI divides by zero in the greedy set-up (`blast_gapalign.c:209`,
  SIGFPE, exit status 136, also with `-ungapped`); with other gap costs the greedy alignment
  converts an infinite quotient to `Int4` (`greedy_align.c:981`, undefined behaviour) and does
  not finish on pairs of 60-100 kb. LOSAT refuses it where it checks the greedy gap costs (a
  batch with a valid query), so a batch of invalid queries ends with NCBI's warning as in NCBI.
  It differs from NCBI where NCBI finishes with other gap costs because no seed reaches the
  greedy alignment (small inputs). `-no_greedy` and `-ungapped` with other gap costs run.
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
  (since version 1.71, line DBID item pfoff-3: the description line of the Seq-align id and
  `Unknown`, then the stop at outfmt 0's first alignment that loads such a sequence under a
  Seq-align id that is not local, or at any alignment of outfmt 1 to 4 that loads one, or outfmt
  0's skipped line for each of its HSPs; the description lines load the sequences that they show,
  in psiblast's rounds after the first those of the two-part table of the first
  `-num_descriptions` HSPs, at least the first, then the alignments theirs up to the stop, after
  which no later batch is prefetched), a
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
  - NaN E-values and bit scores of hits, archives and XML alike: no longer refused since version
    1.67 (NCBI's `HUGE_VAL/HUGE_VAL`, sign set, printed `-nan`, is kept and sorted after every
    number; the approved difference of version 1.67). Version 1.17 (the PS16d review S3-3) had
    refused the hits whose sort would read outside NCBI's vector.

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

- BLASTP `-use_sw_tback` with composition-based statistics (version 1.53, 2026-11-18, line RG item
  X4, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07, within the class of
  1.30; P2: no checkable result exists, P3: one round loop decides it for both programs,
  `core::composition_adjustment::smith_waterman::smith_waterman_redo_window`): a match (a subject's
  preliminary HSPs for one query) with two or more preliminary HSPs when the e-value of a score of
  0 passes `-evalue` (the query's effective search space times the redo's K below it: a large
  `-evalue` or a small `-searchsp`) and the first query letter, after the composition adjustment,
  scores the subject's first letter (after SEG) above 0. The loop of 1.30
  (`Blast_RedoOneMatchSmithWaterman`, `redo_alignment.c:1438-1529`) is the same function with the
  rows of the adjusted square matrix (`matrix[query[i]]`), so it never ends either. LOSAT refuses
  (exit 1, "a match whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI BLAST+ 2.17.0
  ... is not supported by LOSAT's BLASTP") at that match's first significant Smith-Waterman round,
  after the reports of the query batches before (without the epilog; for the first batch the
  outfmt 0 prolog) and the match's batch's reader messages and search strategy request, as NCBI's
  error of a batch's search (`Traceback failed`, below) comes: NCBI prints each batch's reports
  before it searches the next (`blastp_app.cpp:257-276`) and then never ends in the match's batch,
  the form of 1.40 (the Owner's decision of 2026-11-16 on line RG's question Q-X2c: "the results
  of the query batches before are written"). Where that cell scores 0 or less, NCBI's redo of the
  match fails: a failed match that is not the batch's last redone one is dropped, and a failed last
  one fails the batch's traceback (`BLAST engine error: Traceback failed`, exit status 3, after the
  reports of the batches before), which LOSAT reproduces (rule 2; fixture rows `db.x4.rej.swendless`,
  `db.x4.sw.zero.6` and `db.x4.sw.fail.*`; `LOSAT/tests/blastp_sw_tback.rs`). Before that
  round NCBI makes every positive alignment of the window, which takes it about a minute for a
  subject of 1,200 letters against an 80-letter query; LOSAT gives the same outcome at the first
  round. A failed batch whose first query holds a redone match crashes NCBI instead (rule 3,
  2026-11-18, X4 fix round 3): its clean-up frees the arrays of every query's heap of redone
  matches and then pops the first query's heap, whose count it left, through the freed array
  (`blast_kappa.c:3796-3802`, `compo_heap.c:444-457`; SIGSEGV in `BlastCompo_HeapPop` called
  from `Blast_RedoAlignmentCore_MT`, with `-subject` and `-db`), and LOSAT refuses that batch
  (exit 1, "a query batch whose last match fails its Smith-Waterman redo (-use_sw_tback) after a
  match of the batch's first query was redone and kept, with which NCBI BLAST+ 2.17.0's blastp
  crashes, is not supported by LOSAT's BLASTP") after what NCBI formats before that step, as
  above (fixture rows `db.x4.sw.fail.c.*` and `db.x4.sw.fail.h.*`, whose outfmt 0 and 7 stdout is
  NCBI's of a near input that fails at the same step, rows `db.x4.sw.fail.f.0` and
  `db.x4.sw.fail.h_anchor.*`; `db.x4.sw.fail.f.6` and `db.x4.sw.fail.g.6`, whose first query holds
  none, end in `Traceback failed`). The batch's status is that of the last match the redo reaches: a
  match whose subject NCBI's fetch excludes (a subject that a `-negative_seqidlist` keeps without
  any of its lines) sets it back to 0, unless the redo skips that match by early termination: by
  then every query of the batch already holds `-max_target_seqs` alignments whose worst e-value is
  0.002 or less, and the match's preliminary e-value is above 0.01 (`blast_kappa.c:3525-3575`,
  `redo_alignment.c:1560-1580`). A failure before such a match therefore neither fails the batch
  nor crashes NCBI, and LOSAT reproduces that (rule 2; version 1.59; 2026-11-21, line RG item
  X4-excl, rows `db.x4excl.*`).

- TBLASTN `-use_sw_tback` (version 1.73, 2026-11-21, line RG item X4b-1, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07, within the class of 1.30 and the rules of 1.53; P2,
  P3: the same match body and round loop as BLASTP's,
  `core::composition_adjustment::smith_waterman_match`): TBLASTN's redo of a match holds the
  windows of all its frames (`s_WindowsFromTranslatedAligns`), sorted by subject position
  (`s_SubjectCompareWindows`: the subject range's begin, end and frame, then the query range's;
  windows of different frames interleave), in one call
  (`Blast_RedoOneMatchSmithWaterman`), with composition-based statistics 0 to 3. A window of two or
  more preliminary HSPs in which a score of 0 is significant (in e-value mode, `-sum_stats false`
  or a `-max_intron_length` of 1 to 4, when the query's effective search space times K is below
  `-evalue`; with sum statistics, the default, where the redo's cutoff score is 0 or below, which
  a large `-evalue` gives with `-subject` and with `-db`: rows `db.x4b.sp.*`, `-evalue 1e4`; at
  more than one thread LOSAT redoes a batch's matches in parallel where the batch has two or more
  lists, rows `db.x4b.sp.two.*`, and a batch of one list in the one-thread loop; `-dbsize 1` does
  not, rows `db.x4b.cut0.*`) makes NCBI's round loop reach the score of 0: where the query's first
  letter scores the window's first subject letter (after SEG) above 0, by the scores of the matrix
  that the composition-based statistics give (tks_M's subjects never end with 2 and fail with 0 and
  1, rows `db.x4b.rej.endless.6`, `db.x4b.end.c*` and `db.x4b.sp.end.*`), NCBI never ends, and LOSAT
  refuses (exit 1, "a match whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI
  BLAST+ 2.17.0 ... is not supported by LOSAT's TBLASTN"; probes `tbn.fail.*`, `tbn.endless.*` of
  line RG X4b's survey, NCBI still running at 120 s); otherwise the match's redo fails and frees the
  alignments of its earlier windows too, a failed last match fails the batch's traceback (`BLAST
  engine error: Traceback failed`, exit status 3, NCBI probe `tbn.fail.a.rc`), and a failed last
  match after the batch's first query kept a redone match crashes NCBI's clean-up as in 1.53, which
  LOSAT refuses (exit 1, "... with which NCBI BLAST+ 2.17.0's tblastn crashes, is not supported by
  LOSAT's TBLASTN"). Both refusals come after the reports of the query batches before, as in 1.53.
  A match whose subject NCBI's fetch excludes (`check_oid_exclusion`: a subject that a
  `-negative_seqidlist` keeps without a line, `seqsrc_seqdb.cpp:306-314`) and that passes the
  early termination test sets the batch's status back to 0 (`blast_kappa.c:3560-3575`), so a
  failure before it neither fails the batch nor crashes its clean-up: LOSAT keeps such a list in
  the redo's stream as a marker that clears the status (`FetchMarks`,
  `SwBatchStatus::excluded_at_fetch`; rows `db.x4b.excl.*`). The status is per thread in NCBI
  (with `-num_threads 8` the same search ends in `Traceback failed`); LOSAT gives the outcome of
  `-num_threads 1` at every thread count.
  With `-in_pssm`, NCBI's Smith-Waterman redo of a PSSM takes the PSSM's rows by the window's
  positions, not the query's, and reads memory it never wrote (`smith_waterman.c:86-90`,
  `blast_kappa.c:866-873`): depending on the input it fails its traceback (exit 3), throws
  (`Mapping resulted in an empty alignment`, exit 255), crashes (SIGABRT `free(): invalid
  pointer` or `invalid size`, SIGSEGV, with `-subject` and with `-db`), does not end (still
  running at 120 s with `-db`), or ends with exit 0 and rows that change from run to run (three
  runs of the same command gave three outputs; review 2 of X4b-1, probes 1 and 3, 20 PSSM and
  subject pairs). None of these is a result LOSAT can give, so LOSAT, which does not port that
  redo, refuses at the first match it reaches, after the prolog (exit 1, "-use_sw_tback with
  -in_pssm on a search with hits, on which NCBI BLAST+ 2.17.0's tblastn fails, crashes, does not
  end or writes rows that differ from run to run, is not supported by LOSAT's TBLASTN"), where it
  refused the option before the search, and a search without hits ends as NCBI's (exit 0; rows
  `db.x4b.pssm.*`, `db.tbn.inpssm.rej.swtback`).
  Fixture rows `db.x4b.fail.*`, `db.x4b.dbfail.*` and `db.x4b.rej.*`;
  `LOSAT/tests/tblastn_sw_tback.rs`; crafted inputs `LOSAT/tests/fasta/x4b` (X4's proteins
  back-translated) and the databases `LOSAT/tests/db_search/x4b` (`make_x4b.py`). With `-db`, sum
  statistics and the default `-evalue` the cutoff score of the linking is NCBI's floor of 1
  (`BLAST_Cutoffs` starts from a cutoff of 1, SeqDB's shortest length 10): the redo makes every
  local alignment down to that score, slow but finite: LOSAT gives NCBI's output (rows
  `db.x4b.cut.*`).

- BLASTX `-use_sw_tback` (version 1.76, 2026-11-21, line RG item X4b-2, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; the endless redo and the failed batch within the
  class of 1.30 and the rules of 1.53 and 1.73; a new rule-3 rejection, line decision X4b-0 (c);
  approved exception 23, the Owner's decision Q-HEAP=B of 2026-11-21;
  P2, P3: the same match body and round loop as BLASTP's and TBLASTN's,
  `core::composition_adjustment::smith_waterman_match`, with the windows of the query's frames):
  NCBI redoes each frame of a match on its own, with the Karlin block of the frame's query context
  and one array of heaps of redone matches per batch, one heap per query
  (`blast_kappa.c:3601-3622`), with composition-based statistics 0 to 3 (0: scale 1 and the
  unadjusted matrix). The endless redo and the failed batch behave as in 1.73: LOSAT refuses a match
  whose redo NCBI never ends (exit 1, "a match whose Smith-Waterman redo (-use_sw_tback) does not
  end in NCBI BLAST+ 2.17.0 ... is not supported by LOSAT's BLASTX"; tks_M back-translated against
  X4's endless subjects, rows `db.x4b2.rej.endless.*`; with `-dbsize 1` and sum statistics the
  redo's cutoff score is 0, and on the first 240 bases of `bslice_q.fna`'s m_fam NCBI aborts in
  malloc's heap check in one run and still runs at 60 s in others), copies `Traceback failed` (exit
  3) where the batch's last redone match fails (`blast_kappa.c:3775-3779`; rows `db.x4b2.fail.*`,
  `db.x4b2.dbsize1.*`) and refuses the batch where its first query holds a redone match, on whose
  clean-up NCBI crashes (`blast_kappa.c:3795-3802`; exit 1, "... with which NCBI BLAST+ 2.17.0's
  blastx crashes, is not supported by LOSAT's BLASTX"; rows `db.x4b2.rej.crash.*`). A subject that
  NCBI's fetch excludes clears the batch's status as in 1.73 (rows `db.x4b2.excl.*`). New in BLASTX:
  in e-value mode (`-sum_stats false`, or a `-max_intron_length` of 1 to 4, with which NCBI turns
  sum statistics off in a gapped search, `blast_parameters.c:793-805`) NCBI tests a window against
  the heap that its query context indexes in the batch's array of heaps
  (`redo_alignment.c:1465-1474`), and a BLASTX query has six contexts: the context of frame f of the
  batch's query q is 6 x q + f, so the test reads another query's heap while the context is below
  the batch's query count, which LOSAT reproduces (rule 2; rows `db.x4b2.ev.*` and `db.x4b2.ctx.*`,
  where a full heap of another query turns a window away), and past the array at or above it, into
  NCBI's other allocations, where LOSAT tests an empty heap (approved exception 23 below; rows
  `db.x4b2.past.*`, each NCBI's run of the batch followed by queries without a hit). Also new:
  `-ungapped -comp_based_stats 0 -use_sw_tback` (the only composition mode that lets `-ungapped` reach a
  search) makes NCBI's traceback run the redo's set-up, which reads the gapped Karlin blocks that an
  ungapped search does not make (`blast_traceback.c:1486-1498`, `blast_kappa.c:3092-3095`,
  `2073-2074`), and NCBI crashes (SIGSEGV) at every query batch it searches, with hits or without,
  with `-subject` and `-db`; LOSAT refuses at the traceback of the first batch it searches, after
  the reports of the batches without a valid context before it, which NCBI reports with a warning
  and does not search (rows `db.x4b2.rej.ungapped.ac.b1.*`) (exit 1, "-use_sw_tback with -ungapped,
  on which NCBI BLAST+ 2.17.0's blastx crashes at every query batch it searches, is not supported by
  LOSAT's BLASTX"; rows `db.x4b2.rej.ungapped.*`). With `-use_sw_tback` LOSAT searches a query batch
  of NCBI's reader as one unit (the batch's status, its first query's heap and the contexts that
  index the heaps are per batch), and gives the outcome of `-num_threads 1` at every thread count.
  Fixture rows `db.x4b2.*`; `LOSAT/tests/blastx_sw_tback.rs`; inputs `LOSAT/tests/fasta/x4b2` and
  the databases `LOSAT/tests/db_search/x4b2` (`make_x4b2.py`).

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
- makeblastdb `-mask_data` (version 1.52, 2026-11-18, line NP item M4a-1c, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07): a string of a binary mask file whose
  length passes the end of the file and is 1 GiB to 2^47 bytes. NCBI's reader reserves the
  length before it reads the bytes (`CIStreamBuffer::GetChars`, `util/strbuffer.cpp:496-523`):
  from 2^47 bytes, which no x86-64 host can allocate, it throws `std::bad_alloc` ("BLAST ran out
  of memory", exit status 4, which LOSAT writes), past the largest `std::string` (2^63 - 1 bytes)
  `std::length_error` ("Error: basic_string::_M_create", exit status 255, which LOSAT writes), and
  below 1 GiB its reservation succeeds and the read ends at the end of the file (NCBI's
  end-of-file report, which LOSAT writes); between them the reservation succeeds or fails by the
  host's memory and overcommit, the bound of line PS's DECISIONS 40 that `losat rpsbproc` uses.
  LOSAT refuses it (exit 1, "losat makeblastdb: mask data that LOSAT does not read: a string whose
  length ... passes the end of the file ..."), before the build. Fixtures
  `LOSAT/tests/seqdb/mask_parse/b_length_*.asnb` (`make_mask_parse_fixtures.py`).
  An element of a mask's list of `CRef`s (`masks`, a `mix`, an `equiv`, a `packed-int`) stated as
  `NULL` (`05 00` in binary), which NCBI's reader takes and leaves as an empty `CRef`
  (`CObjectIStream::ReadPointer`, `serial/objistr.cpp:1133-1213`), is explicitly rejected under
  rule 3 where makeblastdb dereferences it: the `GetId` of each mask in `x_ProcessMaskData`
  (`**iter`) and the walk of `CSeq_loc::CheckId` (`**ii`, `(*li)->CheckId`,
  `objects/seqloc/Seq_loc.cpp:3541-3620`) stop with a `CCoreException` whose report prints a stack
  trace with the run's library addresses and the path of the installed binary. LOSAT refuses it
  there (exit 1, "losat makeblastdb: mask data that LOSAT does not read: an element of ... stated
  as NULL ..."): after the masks before it, after any error of NCBI's reader in the rest of its
  object, and only where the walk reaches it (after two Seq-ids or a `feat` the walk stops first,
  and LOSAT writes NCBI's "Cannot get masked range Seq-id"). Unit tests in
  `LOSAT/src/writedb/mask_info/tests.rs`; the fixtures `h_pointer_*` and `b_pointer_*` of
  `make_mask_parse_fixtures.py` hold the cases NCBI completes or fails on without that report.
- makeblastdb, an index file that a failed write left empty or cut short (version 1.62, 2026-11-21,
  line LIB-O item O1e, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; the
  outcomes that NCBI reaches and LOSAT reproduces, rule 2, are in the list above). An index file of
  no bytes is mapped as a null pointer, and NCBI crashes with SIGSEGV (76 of the 1,592 cases of
  O1a's sweep); LOSAT refuses it (exit 1, "NCBI crashes on it"). A header field that reaches past
  the last page of the file reads whatever memory follows, so what NCBI does there depends on its
  memory; LOSAT refuses that too ("depends on its memory map"). Tests:
  `LOSAT/tests/makeblastdb_failed_index.rs` (fixture `fsize_rows.txt`) and
  `LOSAT/src/seqdb/damaged_files_tests.rs`.
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
    sparse index; since version 1.58, line RG item Q-BM, LOSAT decodes it as NCBI does) in a
    stream that NCBI reads past (empty or cut data, a big-endian array whose length moves NCBI's
    position past the data), that dereferences NCBI's null reference vector (an XOR block), that
    runs undefined behaviour or writes past BitMagic's blocks or buffers (an Elias gamma code of
    32 or more zeros, `bit_scan_reverse32(0)` in the interpolative decoder, an interpolative value,
    a run or a bit interval past its block, a run of a GAP block that ends before it starts where
    `or_bit_block` shifts by 32 bits or writes past the block, an index past BitMagic's arrays, a
    GAP block of length 0, a window of no words), that reads entries of BitMagic's
    work arrays that this decoding has not written (an empty little-endian id list that no earlier
    list wrote the first entry for, run lengths that no array wrote); line RG D3-1 `bm/SPEC.md`
    §3.3 and §5.2 and line RG item Q-BM. A stream is refused at the first such event; up to it
    LOSAT decodes as NCBI does, NCBI's exception included. A GAP block whose ends do not increase
    up to 65535 stays in the vector as NCBI keeps it; a row of the index that reaches it (its bit,
    or a count of set bits up to a row in or past it) is refused where NCBI's answer differs
    between its builds of BitMagic, as it then depends on the machine NCBI was built for:
    BitMagic's search of a GAP block is SSE2 or SSE4.2 code on x86-64 (BitMagic 8.0.1:
    `bmsse2.h:1638-1642`, `bmsse4.h:2260-2264`; NCBI's configure adds `-msse4.2` there) and scalar
    code elsewhere (BitMagic 8.0.1: `bmfunc.h:1835-1955`), which find a block's bits the same way
    only where its ends increase, and the SSE builds count such a block modulo 2^16 per lane. LOSAT
    runs each build's search and count and answers where they agree (a count past several such
    blocks by each build's own sum); rows that do not reach the block are answered. A bit that the
    decoding sets in such a block (`gap_set_value`) may leave each build a different block: LOSAT
    then decodes the data once for each build and answers each row where the three builds' vectors
    agree; the data is refused only where the builds' decodings end differently (one throws and
    another does not, or the size or the bytes read differ) or where one of them reaches an event
    above (line RG item Q-BM, Owner-delegated, standing instruction 2026-09-29, restated
    2026-10-07). A `bit-bvector`'s size, which NCBI reads without a decode, is NCBI's;

  The `-in_pssm` file's format guess (NCBI's `CFormatGuess`, ported) is refused where one of its
  tests reads memory that the file does not decide: `TestFormatGvf`'s ninth column of a GVF line
  of eight (`tokens[8]`), `TestFormatPhrapAce`'s tokens of an `AS` line of fewer than three
  (`values[1]`, `values[2]`; NCBI 2.17.0 crashes on both, exit status 139) and `TestFormatFasta`'s
  search for `|` past a test buffer of 10 bytes or fewer.

  Since version 1.77 (line RG item D3c-3a, Owner-delegated, standing instruction 2026-09-29,
  restated 2026-10-07) LOSAT also reads a PSSM in ASN.1 binary as NCBI's reader does
  (`CObjectIStreamAsnBinary`, through the same type table and frames as the text), and refuses
  (exit 1) two more cases there. A string or a byte block (OCTET STRING or BIT STRING) whose length
  passes the end of the file and is 1 GiB to 2^47 bytes: NCBI's reader reserves the length before
  it reads the bytes (a string's in `CIStreamBuffer::GetChars`, `util/strbuffer.cpp:496-523`,
  when the bytes are not in its buffer; a byte block's, whose length the binary stream knows, in
  its `vector<char>` reader, `serial/stdtypes.cpp:1538-1548`), which succeeds or fails by the
  host's memory and overcommit (NCBI's probe host: `Unrecognized format for PSSM`, exit status 1);
  below 1 GiB LOSAT writes NCBI's end-of-file answer, from 2^47 bytes `BLAST ran out of memory`
  (exit status 4), the bounds of 1.52, and over 2^63 - 1 bytes NCBI's `std::length_error`
  (`Error: basic_string::_M_create` or `Error: vector::reserve`, exit status 255); a byte block's
  reserve that throws is posted first (`Error: (803.4) [psiblast] unended byte block`), as in
  NCBI (rows `db.ber.ps.block_*`). A query whose letters are missing or whose number is not the
  `length` of its Seq-inst ("a PSSM query whose letters are missing or of another length than its
  Seq-inst is not supported by LOSAT's PSIBLAST", or TBLASTN; refused in the text form before):
  for a query of 30 letters NCBI 2.17.0's psiblast ends by a signal with the lengths 31, 60 and
  255 (SIGSEGV), 0, 1, 2, 10 and 4,294,967,295 (SIGABRT), and writes past its buffer with 28
  (valgrind: an invalid write) though it exits 0. With 29, one less than the letters, NCBI answers
  (exit 0, clean under valgrind) in ASN.1 text and binary alike: it searches the 29 letters with
  the PSSM, and its E-values are those of a query of 29 letters. LOSAT does not yet give that
  answer (the PSSM of 30 columns over a query of 29 letters belongs to line PS: hand-off
  `D/d3/c3a/HANDOFF_PS_nonfinite_lambda.md` of line RG), so it is refused too. Not a rule-3 item but
  a refusal: with `-out_pssm`, author names (`Auth-list` `names`) of a binary PSSM, which NCBI keeps
  in a delay buffer and reads again only while it writes them, where that read posts a message or
  fails (NCBI reports it in the middle of the file it writes; not ported). The warning of a
  UTF8String tag (803.10) is not one of them: NCBI posts it once per process, when it first reads
  the file, so the second read posts nothing and the names are written. Also not a rule-3 item
  but a refusal, of the text form too: a gapped `lambda` of `PLUS-INFINITY`, a gapped `lambda`,
  `kappa` or `h` of `NOT-A-NUMBER`, and psiblast's `lambdaUngapped` of `PLUS-INFINITY` ("a PSSM
  whose gapped Karlin-Altschul parameter lambda is infinite (LOSAT does not yet compute the hits,
  E-values and bit scores that NCBI BLAST+ 2.17.0 computes from it) ..."). NCBI goes on with them
  (E-values 0.0, bit scores `inf` or -1.0); LOSAT's search wrote NaN or found other hits, and its
  arithmetic of these values belongs to line PS (hand-off `D/d3/c3a/HANDOFF_PS_nonfinite_lambda.md`
  of line RG), which lifts the refusal when it computes NCBI's output. `MINUS-INFINITY` stays in the
  refusal of a negative gapped value above. Tests:
  `LOSAT/src/algorithm/pssm_search/in_pssm_tests.rs`,
  `LOSAT/src/report/serial/asn_types/read_ber_tests.rs` and `write_tests.rs`; fixture rows
  `db.ber.ps.*` and `db.ber.tbn.*`.

  None of these is a rule-1 exception: D3-1's questions to the Owner, Q1 (the Dense-seg division by
  zero), Q2 (the memory bounds), Q3 (the signed overflows) and Q-BM, were decided
  by the Owner on 2026-11-16 (「推奨通り」): Q1 to Q3 refuse their cases, as built, and Q-BM ports BitMagic
  (line RG item Q-BM, version 1.58). Not a rule-3 item but a refusal: the XML and JSON
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
  length/3 - 3` (`blast_hits.c:1162-1164`), false from the start for 11 bases or fewer, so it
  returns the translation it never made (`calloc`'d to NULL, `blast_util.c:1284`) plus 1
  (`:1234-1235`), and every caller reads from that address on next: the traceback's gapped start
  and `ALIGN_EX` (`blast_traceback.c:429-440,509-512`), the identities after the composition-based
  redo (`blast_kappa.c:527-528`, which run only for an alignment left by the redo's E-value purge,
  `:3677-3700`), the ungapped re-evaluation (`blast_hits.c:2700-2705`; LOSAT's TBLASTN refuses
  `-ungapped`; TBLASTX reaches it in every search, version 1.51 below). NCBI BLAST+ 2.17.0 dies
  with SIGSEGV in every such run of ours (ws and the local
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
  no report, also where NCBI has written the reports of the query batches before when it crashes
  (but for what its output buffer still holds; 20,000 query letters a batch by default, or
  `BATCH_SIZE`; open for this refusal, line RG item X-short-bat, version 1.51 below).

- TBLASTX, HSPs on a subject (or its `-subject_loc` range) of fewer than 12 bases that reach the
  ungapped re-evaluation (version 1.51, 2026-11-17, line RG item X-short-tbx, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07). TBLASTX is always ungapped: after the
  first linking call and the preliminary E-value (`s_BlastSearchEngineCore`,
  `blast_engine.c:870-899`), `BLAST_PreliminarySearchEngine` re-evaluates the HSPs of every subject
  that has some (`Blast_HSPListReevaluateUngapped`, `blast_engine.c:1493`,
  `blast_hits.c:2609-2737`). It fetches the subject, then takes the translation of each HSP's frame
  from `Blast_HSPGetTargetTranslation` (`:2700`), which for 11 bases or fewer is the translation it
  never made plus 1 (the TBLASTN item above), and `Blast_HSPReevaluateWithAmbiguitiesUngapped`
  reads the HSP's letters from there (`:692-700`). NCBI BLAST+ 2.17.0 dies there with SIGSEGV in
  every such run of ours whose HSPs pass its first linking call (ws, 2026-11-17: 11 bases in all six
  frames, `-subject`, `-subject_loc` and
  `-db` of format versions 4 and 5, `-word_size` 2 (8 to 10 bases too) and 3, `-threshold 1`,
  `-matrix PAM30`, the default SEG, `-outfmt` 0, 5, 6 and 7, `-max_target_seqs 1`,
  `-culling_limit 1`, `-num_threads 4`, ambiguous bases, several subjects, a 12-base subject before
  or after; gdb: `Blast_HSPReevaluateWithAmbiguitiesUngapped` under
  `Blast_HSPListReevaluateUngapped` under `BLAST_PreliminarySearchEngine`; the short queries at the
  end of this item die before, in that linking call). LOSAT refuses there
  (`SubjectEmit::refuse_untranslated_frame`, the predicate and the message shared with the TBLASTN
  item, `core::blast_hits`), where LOSAT wrote rows before this version. A subject that leaves no
  HSP to that point searches as NCBI's: 9 or 10 bases at the default word size, which NCBI does not
  search (fewer than `3 * word_size + 2` bases, `blast_engine.c:1318-1330,1429`), the default
  two-hit window (no frame of three letters holds two words that do not overlap), no HSP out of the
  ungapped extension (`-evalue 1e-30`: the gap trigger caps the cutoff above the HSP's score), and
  HSPs that the preliminary E-value removes (`-evalue 1e-4`: the empty list is freed,
  `blast_engine.c:899,903-904`, and not re-evaluated, `:1480`); so does a subject that the
  re-evaluation's fetch excludes (a `-negative_seqidlist` that names each of its ids, which leaves it
  no definition line: NCBI frees its HSPs before
  it reads the subject's letters or translates them, `blast_hits.c:2659-2665`, and LOSAT drops them
  at that step, `SubjectEmit::reevaluation_fetch_keeps`); 12 bases are translated (rows
  `db.xs.tblastx.*`). Under such a list a subject header that this fetch cannot read (line DBID
  item lineless2, above) ends NCBI's search with a `BLAST engine error` (exit 3) after it has
  printed the query batches before; LOSAT's TBLASTX refusal there now comes after the same reports,
  in the form below (`subject_loop::x_fetch_keeps`; before this version LOSAT wrote none of them).
  The other programs' fetch refusals still come before any report: open (line RG item
  X-short-bat). The refusal is the error of the earliest such subject at every
  `-num_threads`, and the run writes what it writes for the X2c refusal above, as after an NCBI
  error of a batch: the reports of the query batches before, whole and without the epilog, or,
  when the first batch fails, its prologs (`-outfmt 0`'s, and the first lines of the XML2 and JSON
  documents and of `-outfmt 20`), then the failing batch's messages and the error; nothing of the
  failing batch's reports (`tests/tblastx_short_subject.rs`). NCBI has formatted the same when it
  dies, and its crash loses what its output buffer still holds: in our `BATCH_SIZE=1` probes its
  `-outfmt 6` output was as long as the first query's own report (11,318 bytes), its `-outfmt 0`
  output ended a few bytes into that report's footer (50,001 bytes), and in one batch it wrote
  nothing with `-outfmt 6` and the first 18 bytes of the prolog with `-outfmt 0`. Those cuts are
  NCBI's buffer, not a result: LOSAT writes the reports and prologs whole. LOSAT's TBLASTN refusal
  of the same read (the item above) writes no report, although NCBI's TBLASTN searches and prints
  each query batch before it reads the next (`tblastn_app.cpp:275-304`, as TBLASTX's
  `tblastx_app.cpp:176-207`) and so has formatted the reports of the batches before when it dies
  (line RG X-short review 1: with `BATCH_SIZE=1`, NCBI wrote the first query's 43-byte row before
  its SIGSEGV, LOSAT nothing); that difference is open (line RG item X-short-bat). With `-subject`
  and a `-db_gencode` other than 1, LOSAT translates in the selected code
  (PD-TBLASTX-LOCAL-SUBJECT-GENCODE) and refuses where those HSPs reach the re-evaluation. TBLASTX
  translates the query whole before the search, so a short query reaches no such read; but in our
  probes a query of 10 or 11 bases against a subject of 12 to 14 bases, a query of 11 or 12 bases
  against a subject of 11 bases, and 9 bases against 12 or 13 at `-word_size 2` make NCBI die
  earlier, in its first linking call (`BLAST_LinkHsps` under `s_BlastSearchEngineCore`, gdb on ws,
  2026-11-17 and 2026-11-18). Against 12 to 14 bases LOSAT writes rows; against 11 bases LOSAT
  refuses with this item's message, at a later step than NCBI's. Both are open (line RG item
  X-short-tbx-link).

- RPSBLAST and RPSTBLASTN over an alias file whose `FIRST_OID` and `LAST_OID` leave out profiles of
  its volumes (2026-11-17, line RG item D1a-4, Owner-delegated), or whose `OIDLIST` does (line RG
  D1b-1, 2026-11-21: the same form, probes R5 and R7 of D1b): NCBI's RPS search
  ignores the range (its lookup table holds every profile of a volume, `blast_aalookup.c:182-184`;
  the volumes of several are searched bare), while its formatter opens the `-db` name with its keys,
  whose `SeqidToOids` keeps only the OIDs of the range (`seqdbimpl.cpp:1414-1418`): a hit of a
  profile outside it is not found. Outfmt 0 then writes `Unknown` and "Sequence with id CDD:N no
  longer exists in database...alignment skipped", the tabular formats `Unknown`, `N/A` and a length
  of 0 (the class of rule 2's version 1.10, exit status 0), and `sacc`, `saccver` and the default
  fields of outfmt 6, 7 and 10 stop with a null pointer exception, exit status 255 (the class of
  1.37; probes p6_1v_fo22 and p6b_fo22 of line RG D1a-4). LOSAT refuses (exit 1, "a hit of the
  profile of OID N outside the RPS-BLAST database's alias files' FIRST_OID, LAST_OID and OIDLIST ... is not
  supported by LOSAT's RPSBLAST", one text for both programs) where a query's reported profiles
  hold such a profile: the profiles of the first `-max_target_seqs` subjects of its hit list (the
  hit list size, `PruneSeqalign` with `m_HitlistSize`, `blast_format.cpp:810-813`, which counts
  following profiles of one Seq-id as one subject, as with a volume given twice under two paths:
  fix round 2 of D1a-4; outfmt 0's descriptions and alignments lie within it). NCBI's report looks
  a reported profile up by its Seq-id among the profiles the range includes: where an included copy
  of its record holds it (a volume given twice under two paths, one past the range), NCBI shows that
  copy and exits 0, and LOSAT writes the same (probe dup6 of D1a-4 fix round 2); an included profile
  of another record or of other letters that holds it is refused too ("... whose Seq-id NCBI BLAST+
  2.17.0's report finds in another record or other letters ...": NCBI shows that profile's letters
  and length with the HSPs of the one out of the range; fix round 3, review C N4). RPSBLAST refuses before any output (its search of every query batch
  comes before the reports); RPSTBLASTN after the reports of the query batches before. With several
  volumes a profile outside the range past that cut of the merged lists, which no report of NCBI's
  looks up, is NCBI's run: its HSPs count in outfmt 7's `# N hits found` (and in
  `PRE_FETCH_SEQS_LIMIT`'s count), as NCBI counts the whole merged set (probes r1_3v_fo3 and
  r1_3v_fo11 of review A of D1a-4). A run whose reported hits are all in the range is NCBI's
  (fixture rows `db.rps.rej.fo22`, `db.rpstbn.rej.fo22`, `db.rpstbn.fo11.e5.6`, `db.rps.1v.fo3.*`,
  `db.rps.3v.fo3.*`, `db.rpstbn.3v.fo11.*`; for `OIDLIST` `db.ol.rps.rej.no2w`, `db.ol.rps.rej.only2w`,
  `db.ol.rtbn.rej.no2w` and `db.ol.rps.only2w.mts1.6`).
- RPSBLAST and RPSTBLASTN over an alias file whose `OID_MASK_TYPE` removes a definition line of a
  profile, in a database of several volumes (version 1.70, 2026-11-21, line RG item D1b-2,
  Owner-delegated, the form above): NCBI searches the volumes bare (`rpsblast_local.cpp:380-399`),
  so its Seq-aligns keep a profile that the key leaves without a line under its stored Seq-id,
  while its formatter, over the `-db` name with the key, does not find it (probe xp3v of D1b-2,
  profile `ref|XP_000102.1|` in the second of three volumes): outfmt 0 writes `XP_000102.1 Unknown`
  and "Sequence with id XP_000102 no longer exists in database...alignment skipped" (exit status
  0, the class of 1.10), and the default fields of outfmt 6, 7 and 10 stop with a null pointer
  exception, exit status 255 (the class of 1.37). LOSAT refuses (exit 1, "a hit of the profile of
  OID N whose definition lines the RPS-BLAST database's alias files' OID_MASK_TYPE removes, in a
  database of several volumes ... is not supported by LOSAT's RPSBLAST") where such a profile is
  among a query's reported profiles (the same cut as above, the Seq-ids of the stored headers);
  past the cut it is NCBI's run, its HSPs counted in outfmt 7 (rows `db.om.rps.x3vc.*`). With one
  volume NCBI searches the `-db` name, whose Seq-aligns leave such a profile out, and LOSAT does
  the same (rows `db.om.rps.x.*`, `db.om.rtbn.x.*`; refusal rows `db.om.rps.rej.x3v`,
  `db.om.rtbn.rej.x3v`).

- RPSBLAST and RPSTBLASTN over one volume whose alias file's `LENGTH` and `NSEQ` make the
  concatenated length `Lc = (Int4)(LENGTH + NSEQ)` shorter than the profiles (2026-11-17, line RG
  item D1a-4 part B, Owner-delegated; precedent 1.30): NCBI's gapped stage roots
  its interval tree at `[0, Lc + 1]` (`blast_gapalign.c:3830-3831`) while the HSPs of later profiles
  lie past it. A region (the profile's offset plus the HSP's range) that starts past `R = Lc + 1`
  goes right at every node `[l, R]` of the right spine, whose middles `(l + R) / 2` are at most `R`,
  and becomes the leaf at its end; a second one in the same tree (one tree per query, per frame for
  RPSTBLASTN) meets that leaf, and `BlastIntervalTreeAddHSP` (`blast_itree.c:603-793`) splits
  `[R, R]` into `[R + 1, R]`, whose middle is `R` again, without end, doubling its node pool until
  an allocation fails (`BLAST ran out of memory`, exit status 4 under an 8 GiB limit, after outfmt
  0's prolog; the OOM killer without it). A region that crosses `R` straddles a middle of the spine
  and ends in a subject tree; below `R = 0` every middle of the spine is below the regions' starts.
  LOSAT refuses (exit 1, "... where two gapped alignments of a query start past position ... is not
  supported by LOSAT's RPSBLAST") at the insertion that would make the node, after NCBI's end-point
  checks (which may drop the HSP or the leaf it meets): `BlastIntervalTree::try_add_hsp`, the stop
  of TBLASTN's endless split of 1.18; RPSBLAST before any output, RPSTBLASTN after the reports of the
  query batches before (NCBI writes the batches before, then fails). With at most one such
  region in each tree NCBI ends, and LOSAT gives its bytes: the truncated extensions to the right
  and the diagonal table of `Lc` (fixture rows `db.rps.rej.l0`, `db.rpstbn.rej.l0`,
  `db.rps.rej.nm300`, `db.rpstbn.rej.nm300`, and the one-volume rows that end).
  An alias file's `OIDLIST` without `NSEQ` and `LENGTH` above it makes NCBI's constructor scan the
  profiles its list keeps, and the scan's total and count are `LENGTH` and `NSEQ` here (line RG D1b-1,
  2026-11-21, probe R1 of D1b: `BLAST ran out of memory`, exit 4; LOSAT's text names the list; rows
  `db.ol.rps.rej.no2`, `db.ol.rtbn.rej.no2`).

- RPSBLAST and RPSTBLASTN over one volume whose concatenated length `Lc` is longer than the profiles
  (an alias file's `LENGTH` and `NSEQ`, or `-db "a x"` with an alias `a` of the volume `x`, which counts
  it twice; 2026-11-17, line RG item D1a-4 part B, Owner-delegated): the ungapped
  extension to the right runs to `MIN(subject->length - s_off, Lc - q_off)` (`aa_ungapped.c:948`),
  past the rows. The last profile's sentinel row (`-BLAST_SCORE_MAX`, `makeprofiledb.cpp:1166-1167`)
  stops it unless the X-drop is above that score; past it NCBI reads the row pointer after the last
  row (`rps_pssm[rows]`, the last of the `rows + 1` pointers, `blast_aalookup.c:185-190`), which
  points to the bytes after the rows in the `.rps` file's mapping, and then past its pointer array.
  At `Lc = rows + 1` the extension reads that one row: the file ends there (makeprofiledb writes no
  byte after the rows), so NCBI reads the zeros of the file's last memory page, the row scores 0
  and the extension ends at its bound (`NSEQ 4` over the 563 rows of `rps_mini`, whose `.rps` of
  63,080 bytes leaves 2,456 zero bytes in its last 4 KiB page: NCBI exit 0, 431 bytes, probe g9#0 of
  review B of D1a-4); where the file ends on a page boundary that read faults (SIGBUS). LOSAT reads
  that row as NCBI does (the file's bytes after the rows, then the page's zeros: DECISIONS D1a4-9,
  rule 2, the precedent of DECISIONS X2f-4 for `.psq`) and refuses (exit 1, "... whose ungapped
  extension of a query runs past the row after the last profile, or reads that row past the .rps
  file's last memory page ... is not supported by LOSAT's RPSBLAST") an extension that reads its
  bytes past the last page (rule 3) or goes on past it, where `Lc` is the rows plus 2 or more: with
  `-xdrop_ungap 100000` and a query that matches the profiles to their end NCBI crashes (SIGSEGV,
  exit status -11; probes p9_over of line RG D1a-4 and g9#1, `NSEQ 5`); RPSBLAST before any output,
  RPSTBLASTN after the reports of the query batches before. Elsewhere `Lc` above the rows is NCBI's
  search (its diagonal table of `Lc`; fixture rows `db.rps.rej.over`, `db.rps.rej.n5.over`,
  `db.rps.1v.n4.over.6`).

- RPSBLAST and RPSTBLASTN with a concatenated length `Lc` whose diagonal table never stops growing
  (2026-11-17, line RG item D1a-4 part B, Owner-delegated): `s_BlastDiagTableNew`
  (`blast_extend.c:52-57`) doubles an `Int4` from 1 while it is below `Lc` plus the word window; past
  a sum of 2^30 it doubles to 2^31 (`INT_MIN`) and then to 0 and stays there (an alias file's
  `LENGTH` 2^30 over one volume: no end, timeout); a sum above `INT4_MAX` (`LENGTH` 2^31 - 4 and
  `NSEQ` 3) overflows NCBI's `Int4` addition. LOSAT refuses (exit 1) a unit whose `Lc` plus the
  window is above 2^30, for every unit (also a bare volume of that many rows; the volume test of
  earlier versions refused a sum of exactly 2^30, where NCBI's loop ends; fixture row
  `db.rps.rej.l2p30`), where NCBI's engine meets it: in the unit's search of the first query batch
  that NCBI searches, after the batch's query set-up (an empty query file ends with NCBI's `Query is
  Empty!` and exit 0, a first batch without letters with NCBI's engine error and exit 3: rows
  `db.rps.empty.l2p31`, `db.rps.nolet.l2p30`, `db.rpstbn.nolet.l2p30`; review B of D1a-4, MF1).

- RPSBLAST and RPSTBLASTN over one volume whose `LENGTH` and `NSEQ` sum above `INT4_MAX`
  (2026-11-17, line RG item D1a-4 part B, Owner-delegated): `s_RPSPreliminarySearchEngine`
  returns -3 (`blast_engine.c:1131-1135`) after `RPSPsiMatrixAttach` gave the score block the
  lookup table's profiles, and the score block's free frees them again (`double free or corruption`,
  SIGABRT, exit status -6 after outfmt 0's prolog; `LENGTH` 2^31). LOSAT refuses (exit 1; fixture row
  `db.rps.rej.l2p31`) where NCBI's engine meets it, as the diagonal table's refusal above.

- `losat makembindex -legacy false` with a `-stride` above `-ws_hint` − `-nmer` + 1 (version 1.60,
  2026-11-21, line RG item D5a, Owner-delegated): the list's `mult_` is 0 and the loop over its
  divisor classes in `COffsetList::Save` (`dbindex_factory.cpp:1222-1252`) starts below 1 and
  never ends (probes `-stride 20` with the defaults and `-stride 30 -nmer 8`, random input and an
  all-N record: timeout, a partial `.00.idx` of 262,188 bytes). LOSAT refuses (exit 1) at the
  first sequence of an index volume; an input without sequences still gives NCBI's
  `removed (empty)` and exit 0.
- `losat makembindex -legacy false` with a `-stride` of 32768 or more (same version and line):
  `GetMinOffset` (`dbindex.cpp:446-447`) shifts an `int` by 32 bits or more, which C++ leaves
  undefined; the official x86 build finishes with what its shift instruction gives (exit 0, the
  same file on a rerun). Computed in a wider type, the offsets no longer fit the index's 32-bit
  words, so no valid index can be defined. LOSAT refuses (exit 1) at the first sequence of an
  index volume.
- `losat makembindex` of a database with a sequence without letters (version 1.61, 2026-11-21, line
  RG item D5a-t, Owner-delegated): `CSeqDB::GetBioseq` gives a null Bioseq for it
  (`seqdbvol.cpp:1380-1386`) and `CSequenceIStreamBlastDB::next` dereferences it
  (`sequence_istream_bdb.cpp:105-107`): `CCoreException::eNullPtr`, exit status 3, with a stack
  trace of the run's memory addresses (probe `remaining-gaps-20261112/D/d5/probe/zero`, a database
  whose first sequence's offsets were patched to 0 letters; `-idmap` true or false). LOSAT refuses
  (exit 1) at that sequence. windowmasker skips such a sequence, as NCBI's, and makes no title or
  taxonomy lookup for it.
- `losat makembindex -old_style_index true -iformat blastdb` of an alias file whose `OID_MASK_TYPE`
  erases every definition line of a sequence (version 1.70, 2026-11-21, line RG item D1b-2,
  Owner-delegated): `CSeqDBVol::GetBioseq` gives no Bioseq for it either (`seqdbvol.cpp:1362-1366`),
  dereferenced as above (probe `mk_old_ntall` of D1b-2: exit status 3, a stack trace). LOSAT refuses
  (exit 1) at that sequence; a new style index reads no alias key, as NCBI's. NCBI's dustmasker stops
  at such a sequence with `CUnassignedMember::eGet` `Seq_entry.set` (exit status 3), which LOSAT
  writes as NCBI does (rule 2); segmasker and the windowmasker masking step skip it, as NCBI's.
- A sequence of a volume under an alias file's `OID_MASK_TYPE` with bit 0 set, one of whose Seq-ids
  NCBI's accession guide posts a warning for (same version and line, Owner-delegated): NCBI
  identifies the ids of the lines at each read of the header that its cache of filtered headers
  does not answer (`seqdbvol.cpp:1977-1988`, 2006-2058) and posts the warning each time, so their
  number follows that cache. LOSAT refuses (exit 1) where the sequence's lines are first read.

- An alias file's `OIDLIST` file (line RG D1b-1, 2026-11-21, Owner-delegated; DECISIONS D1b1-1 of line RG) that
  is empty: NCBI maps nothing and reads the 4-byte header through a null pointer (SIGSEGV, probe G5 of
  D1b); or whose bits NCBI copies for a volume past the memory page after the file's end (a volume
  of more OIDs than the file's page holds bits for, under a header that does not stop the copy): NCBI
  reads memory that is not the file's. A read inside that page reads the page's zeros, as NCBI does
  (probes H2 and H7, not refused). A file of 1 to 4 bytes fails NCBI's assertion that the file is
  longer than its header, whose text names a line of NCBI's source: LOSAT's words, refused where NCBI's constructor
  reads the file (its scan; exit 2 in NCBI), and NCBI's `Database memory map file error` (exit 3) at
  a search's first read under `NSEQ` and `LENGTH`, as for the column files' assertions
  (`seqdb::column`). Rows `db.ol.bp.rej.empty`, `db.ol.bp.rej.hdr` and `db.ol.bp.hdrl.6`.

- An alias file's `FIRST_OID`-`LAST_OID` range that lies outside the bit array of the OIDs that the
  ranges of the alias files below it keep (line RG D1b-1, 2026-11-21, Owner-delegated; DECISIONS
  D1b1-10 of line RG): NCBI's intersection clears the array's bits outside the range at the OID less
  the array's first, so a range that ends before the array's first OID wraps that index below 0 and
  writes far outside the array (SIGSEGV where the OID list is built: `blastdbcmd -entry all`, probes
  D1a_P and D1a_R of review A of D1b-1; `-info` builds no list and works), and one that starts past
  the array's bytes writes past them, over NCBI's heap: NCBI aborts (`free(): invalid pointer`, a
  failed assertion of glibc's allocator) or crashes (SIGSEGV, a reference counter's overflow) where
  the clearing reaches 24 to 2,500 bytes from the array's first (probes RB9_T200 to RB9_T19990,
  RB9_U5000 and RB9_U19990 of review B of D1b-1), and keeps no OID with exit status 0 where it
  overwrote nothing it uses next (RB_pc, RB9_T21, RB9_U300; DECISIONS D1b1-13): a result of NCBI's
  allocator, not of its reader, refused as the reads past an OID list file's page above. Inside the array's last byte the
  clearing is harmless and no OID is kept (probe D1a_Q, case A237). LOSAT refuses where the list is
  first used (cases A238 and A241 of `blastdbcmd_alias_keys.rs`, `alias::OidSet::clears_outside`).

- `losat blastn -use_index true` of a megablast search (version 1.65, 2026-11-21, line RG item
  D5b-1, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07), where NCBI's
  loader of the MegaBLAST database index meets one of these. LOSAT refuses each (exit 1) before
  the search:
  - The first index volume (`<name>.NN.idx`, new or old style) is a file of 0 bytes. NCBI maps it
    as no memory and reads its header through a null pointer when it checks the word size
    (SIGSEGV, exit status 139, also without read access; fixture row `used_old_empty` of
    `tests/blastn_mb_index`). A first volume that is not a regular file, does not open or has a
    word size field above 15 gives NCBI's warning and its search without the index (rows
    `failread_*`); one shorter than its header reads zeros past its end, as NCBI's mapping does.
  - A new-style index volume of 36 to 39 bytes, or one NCBI cannot open: NCBI's count of its
    sequences reads a word that its stream never set, an uninitialized value (row
    `undefined_idx38`, where NCBI printed 42,008,582). Shorter volumes give NCBI's exception
    (`at start oid`, `at end oid`), reproduced.
  - An `-index_name` range whose last volume converts to the largest `unsigned long` (a last
    volume of -1, such as `oidx,,0,-1`), or whose last volume is past `LONG_MAX` while its first
    is not: NCBI's loop over the volume numbers never ends (row `old_endless`). Other ranges,
    negative ones too, end and are reproduced.
- `losat blastn -use_index true` with a MegaBLAST index of only some database volumes, in a later
  query batch or query chunk that reads the results NCBI keeps (version 1.74, 2026-11-22, line RG
  item D5b-3, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; the kept
  results are reproduced, rule 2 above), where NCBI's search thread holds the empty list of the
  kept results as its list of initial hits, of no memory. The scan of a database volume without an
  index then saves its first initial hit with a reallocation of 0 bytes, which glibc answers with a
  block of one hit, and its second with another, which frees the list (`blast_extend.c:325-356`);
  NCBI then writes to and frees freed memory: `double free detected` (exit status 134, rows
  `m_part_two` and `m_part_c1200`), a crash (139, row `m_mixed2`), or a result read from freed
  memory. LOSAT refuses (exit 1) a search whose scan saves a second initial hit into such a list,
  at that batch or chunk after the output of the batches before (`IndexSeeds::check`), at every
  thread count by the list of NCBI's one thread; one hit, or none, is searched as NCBI does (row
  `m_part_one`).

- `losat windowmasker` (version 1.72, 2026-11-22, line NP items M5a and M5b, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07). LOSAT refuses each (exit status 1, the
  manual's `refusals.md`) where NCBI meets it:
  - `-mk_counts` of an input without a unit of unambiguous letters (all `N`, or an empty input
    with `-unit`): NCBI writes the percentages of the file as its C runtime prints 0/0 (`-nan` on
    x86-64 Linux), a text that depends on the runtime.
  - `-mk_counts -infmt seqids`: NCBI fetches each sequence from NCBI's servers, as masking with
    `-infmt seqids` is refused.
  - `-mk_counts` or `-convert` to `oascii` or `obinary` with units of 1 to 3 letters and an
    `-smem` whose low 16 bits are 0: NCBI's `Uint1` key size wraps below 0 to 253
    (`seq_masker_ostat_opt.cpp:222-242`), and NCBI then allocates 8 GiB and shifts by 253 bits
    (undefined behaviour in C++). Units of 4 letters or more give NCBI's insufficient-memory
    error, which LOSAT writes.
  - `-mk_counts` whose counters (4^suffix `Uint4`) or cache bit array (4^unit/32 words) cannot be
    allocated on the host: NCBI's result then depends on the memory of the computer
    (`std::bad_alloc`, or a file without the cache bit array after a warning).
  - `-convert` of an ascii counts file with units of 17 letters or more: the converter's reverse
    complement shifts the 32-bit unit past its width (undefined behaviour), as masking with such a
    file is refused.
  - `-convert` of an ascii counts file whose unit size line is 0 (or a multiple of 256) to `oascii`
    or `obinary`: NCBI's search for the size of its hash table never ends
    (`seq_masker_ostat_opt.cpp:217-223`).
  - A counts file whose `##pct:` line holds a number that C++ libraries read differently (`inf`,
    `nan`, a hexadecimal number, a number with letters after it, a number past the range of a
    double), with `-convert` or a `-t_*_pct` option; and a `binary` or `obinary` file whose
    `##pct:` line names 4294967295 or more as its largest count, where NCBI's table has no room for
    it and its loop reads past the table.
  - An optimized ascii counts file whose hash-parameter line leaves a value unread, which NCBI then
    checks uninitialised; an `obinary` file with units of 1 or 2 letters, whose empty array NCBI
    reads past; a FASTA record or database sequence of 4,294,967,296 letters or more, whose
    positions NCBI counts in 32 bits.

- BLASTP and BLASTX `-db_hard_mask` on a database whose masks overlap, nest or are stored out of
  order (version 1.78, 2026-11-22, line NP item M4b-3b, Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07). Such a database is searched as NCBI searches it (rule 4): NCBI
  cuts each sequence into parts at the complement of its masks in the order stored, a part
  between overlapping masks "ends before it starts" and holds no letters, and its table of
  diagonals moves back by it (`blast_engine.c:221-276`); segmasker writes such masks (A. thaliana
  TAIR10, UniRef50). A part shorter than the word gets one word from NCBI's two-hit scan, read
  from the stored letters at the part's start (`aa_ungapped.c:495-505`); the seeds it saves are
  extended with the part's length as `subject->length` (`blast_gapalign.c:4246-4295`: the left
  extension reads the letters up to the gapped start, the right one runs only before the part's
  last letter) and kept in an interval tree whose subject bounds that length sets
  (`blast_itree.c:510-795`), or are the part's HSPs in an ungapped search; LOSAT does the same
  (`BlastpPrelimKernel::one_word_gapped`, `blastx::search::one_word_hits`). A sequence whose
  chunks' word hits reach the table of diagonals of another sequence is searched on its own
  (approved exception 24). Fixture rows `db.pmask.bp.hard.ov_*` and
  `db.pmask.bx.hard.ov_*` of `LOSAT/tests/db_search_manifest/{blastp,blastx}.pmask.tsv`, databases
  `ov_*` of `LOSAT/tests/db_search/dbmaskp/` (`ov_s2`: the seed of a part of 2 letters, NCBI's rows
  the same at 1 thread twice and at 4 threads; `ov_nest` and `ov_big` at `-window_size 0`: hits
  below the start of the only sequence). LOSAT refuses (exit 1) the shapes below, one refused
  sequence refusing the whole database. Evidence: NCBI BLAST+ 2.17.0 runs of probe 2 of line NP's
  `M/m4b3b/NOTES_M4B3B.md` and of `M/m4b3b/probe/fix1`, each capped.
  - `algorithm::common::hard_chunk_scan::short_chunk_tree_refusal`, during the search: a part
    between overlapping masks shorter than the word whose seed's gapped alignment lies outside the
    subject bounds of the tree that the part's length (negative) sets, where NCBI's insertion splits
    the tree's nodes without end (`BlastIntervalTree::try_add_hsp`, an exact port that reports
    where the insertion never ends). Evidence: `g3` (consecutive inverted parts, `WWW` at both
    offsets), BLASTP two-hit gapped: NCBI does not end (rc -11 at 60 s under a 4 GB limit); with
    `-ungapped` it ends, and LOSAT gives its rows; BLASTX's search of `g3` ends with no row, and so
    does LOSAT's.
  - `HardMaskShape::WordPastEnd`, before the search, BLASTP only and only for a query that holds
    residue 0 (the gap letter `-`, NCBIstdaa 0): a part between masks shorter than the word whose one
    word NCBI reads past the sequence's end, from the next sequence's letters, where its hits can
    change the rows. No other query can hit such a word: the matrices score residue 0 at
    `BLAST_SCORE_MIN` (`blast_stat.c:1559-1578`), so no neighbouring word holds it, and the exact words
    that the lookup table always adds are the query's (`blast_aalookup.c:503-509`); BLASTX's
    translation gives no residue 0. No query read from FASTA holds it: NCBI's `CFastaReader` drops
    `-` with the warning "Hyphens are invalid and will be ignored", as LOSAT's reader does, so this
    rejection does not occur with NCBI's query input. Evidence: `ov_wpe` (`M/m4b3b/probe/fix1/wpe`):
    NCBI's rows of A then B are B's alone, the same in every run, also for a query file with `W-W`
    (rows `db.pmask.bp.hard.ov_wpe.*`). Counts for a query that would hold residue 0: TAIR10 36, the
    UniRef50 sample 78.
  - `HardMaskShape::TailMaskWithLaterMasks`, before the search: a mask that runs to the end after the
    part from the start, with masks stored after it. NCBI's chunk loop ends at that mask and searches
    the part past its end in the no-chunking branch (exception 21's class) without the later masks,
    and applies them to the letters of its alignments; no near input defines the result. Evidence
    (`M/m4b3b/probe/fix1`, ov_s2's letters): the mask [65, 240) alone gives one row (q 31-32,
    s 62-63), with [55, 64) stored after it none, with [100, 120) inside it the same row; with
    `-window_size 0` none of the three ends in 30 s (one BLASTX run of the first SIGSEGV). Counts:
    TAIR10 0, the UniRef50 sample 0.
  - Not yet ported (LOSAT does not do NCBI's merge across a chunk overlap yet, as with a range longer
    than 5,000,000 letters; line NP item M4b-3d): `HardMaskShape::PartsOutOfOrder`, before the
    search: masks stored out of order whose part after the later-stored mask starts before the end of
    a part before the earlier one, both able to hold a word: NCBI joins the two parts' alignments by
    `Blast_HSPListsMerge` across that overlap (as the rejection of a last mask shorter than the
    letters before its end, `ShortLastMask`); masks out of order whose parts do not meet so are
    searched. And `hard_chunk_scan::short_chunk_merge_refusal`, during the search: a part shorter
    than the word whose alignments from its one word meet those of another part in that merge.
    Counts of `PartsOutOfOrder`: TAIR10 0, the UniRef50 sample 0.
  - `HardMaskShape::Inverted`, before the search: a mask that ends before it starts, which makeblastdb
    never writes.
  TBLASTN and TBLASTX refuse masks that overlap or are stored out of order (`Unsorted`), as before
  (line NP item M4b-3c); BLASTN is unchanged.

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

## Approved exception (a database length or `-dbsize` of exactly -2^63 in NCBI's serial output, line RG item D1a-4)

The Owner decided on 2026-11-17 (line RG `DECISIONS.md` OQ-1117, D1a-3's Q1 = A,
words 「A 正しい数を書く (推奨)」) the question of line RG items D1a-3 and D1a-4: LOSAT writes the
value's digits, as it writes the `Int4` `INT_MIN` of an archive as `-2147483648` where NCBI writes
`-8` (line DI4, decision 13), against the rule-3 refusal of that one value that LOSAT made until
then. Numbered after exception 16 (version 1.54, line RG item X5; written as exception 15, then 16,
and renumbered at the merges of dev). Its version is in the list at the top; elsewhere it is cited by its number and
name, "approved exception 17 (Int8 -2^63 written as digits)", and D1a-4's other entry as "the
D1a-4 entry of PD-LOSAT-NCBI-DEFECTS", so that a merge of dev that takes their version numbers
moves only that list (line RG `DECISIONS.md` D1a4-13).

17. **A database length or `-dbsize` of exactly -2^63 that NCBI's serial writer prints.** `COStreamBuffer::PutInt8`
    (`util/strbuffer.cpp:1002-1010`) negates a negative value before it writes the digits; for
    INT64_MIN that is a signed overflow, and NCBI BLAST+ 2.17.0's build prints `-0`. NCBI writes
    such a value where an alias file's `LENGTH` (or `STATS_TOTLEN` over several RPS volumes) is
    exactly 9223372036854775808, NCBI's Int8 -2^63: BLAST XML's `<Statistics_db-len>-0</Statistics_db-len>`
    and XML2's `<db-len>-0</db-len>` (outfmt 5, 14 and 16; probes p5_neg_2p63 and p5b_neg_fmt of
    line RG D1a-3), the `DbLength` of RPSBLAST's and RPSTBLASTN's search strategy requests after the
    first (`big-integer -0`; probe r2_3v_l2p63 of D1a-4's review A), and with `-dbsize
    -9223372036854775808` the `DbLength` of every program's strategy and of the BLAST archive of ASN.1
    text or XML (probe x_dbsize_min of D1a-4's fix round 1). Its JSON writer prints
    `-9223372036854775808`, the binary writer the value's bytes. LOSAT writes `-9223372036854775808`
    in every one of them (`report::serial` `big_int`); every other byte of these runs is NCBI's
    (probe minus0 of D1a-4's fix round 2: the nine commands that LOSAT refused, each NCBI's output
    with that one value written as its digits). The writers' unit tests hold the value, as for
    `INT_MIN` (no fixture row kind holds NCBI's bytes with one value written otherwise).

## Approved exception (psiblast `-in_msa`, a `secondary-accession` number past `int`, line PS item PS6b, version 1.64)

Proposed under the standing instruction of 2026-09-29 (restated 2026-10-07), under the Owner's rules
of 2026-10-10 「NCBI のバグを再現する必要はない」 and 2026-10-22 「広い型になおしましょう」 and rule 1 above;
accepted by the maintainer on 2026-11-21 by marking PR #1440 ready for review, as the PR asked.

18. **A number of a `secondary-accession` value that does not fit `int`.** `SSeqIdRange`
    (`objects/seqloc/Seq_id.cpp:3600-3640`; `int start; int stop;` in `Seq_id.hpp`) reads the digits
    with `start = start * 10 + digit`, a signed overflow from the tenth digit on (undefined behavior
    in C++). NCBI BLAST+ 2.17.0's build wraps: `AB9999999999` gives the accession `AB1410065407` in the
    GB-block's `extra-accessions` and `hist.replaces` (rc 0), and `AB2147483648` or a range that
    crosses 2^31 gives a negative number that `CSeq_id` refuses, `Malformatted ID AB-2147483648`
    (rc 255; probes `M.secondary_num_big`, `M.secondary_num_over`, `M.secondary_range_over` of
    the PS6b review). LOSAT reads the numbers in 64 bits and writes each accession with the digits as
    written, zero-padded to the first number's width as NCBI pads every number that fits (NCBI 2.17.0
    `AB1000000000-AB1000000002`, `AB0999999999-AB1000000001` and `AB1000000000`: the same bytes as
    LOSAT's, probes `S.a`, `S.b`, `S.c` of the fix round); the accession then goes through the same
    Seq-id parse as every other, which decides by its own rule whether the run goes on. A
    number of more than 18 digits (past 64 bits) and a range that runs backwards (NCBI's loop does not
    end) stay refused (rule 3). The range's size has no limit: NCBI expands any range into as many
    Seq-ids, finishes up to the memory it has (2e7 ids in 115 s and 6.9 GB), and LOSAT writes the ids
    as it expands them, in less memory than NCBI.

## Approved exception (hit list sizes past `Int4`, NCBI parity line O1c, version 1.66)

Proposed under the standing instruction of 2026-09-29 (restated 2026-10-07), under the Owner's rules
of 2026-10-10 「NCBI のバグを再現する必要はない」 and 2026-10-22 「広い型になおしましょう」 and rule 1 above;
accepted by the maintainer on 2026-11-21 by marking PR #1448 ready for review, as the PR asked.

19. **A hit list size whose preliminary size, or twice it with `-best_hit_*`, passes `Int4`.** `GetPrelimHitlistSize`
    (`core/blast_hits.c:44-70`) computes the preliminary hit list size in `Int4`:
    `MIN(MAX(2*h,10),h+50)` for a gapped search without composition-based statistics, `2*h+50`
    above 500 with them (from 1000 with `ADAPTIVE_CBS`), a signed overflow (undefined behavior in C)
    from h = 2^30 without them and from h = 2^30 - 25 (1073741799) with them. NCBI BLAST+ 2.17.0's
    build wraps: without composition-based statistics a size from 2^30 to 2^31 - 51 keeps 10
    subjects per query and a larger one crashes (SIGSEGV); with them a size from 1073741799 to
    2147483623 crashes and one from 2147483624 keeps 2 to 48 (probes of line O1c, `-max_target_seqs`
    and BLASTN's `-num_descriptions` with outfmt 0, in BLASTN, BLASTP, TBLASTN, BLASTX, PSIBLAST,
    RPSBLAST and RPSTBLASTN). The size is the one of every collector of the preliminary stage
    (`hspfilter_collector.c`, `hspfilter_culling.c`, `hspfilter_besthit.c`, whose `2 *
    prelim_hitlist_size` LOSAT computes in 64 bits too), and PSIBLAST's traceback takes it as its hit
    list size (`traceback_stage.cpp:214-223`). LOSAT computes the size in 64 bits, once for every
    program (`blastn::hsp::get_prelim_hitlist_size`; on a 32-bit target a size past `usize` is its
    largest value, which no hit list reaches, since NCBI's OIDs are `Int4`), so each such size keeps
    every hit, as NCBI's 1073741798 (2^30 - 26) does: NCBI's output at 1073741798 equals LOSAT's at
    these sizes (the rows with the contract `approved_hitlist_wrap:1073741798` of
    `LOSAT/tests/outfmt0_manifest.tsv` and the `ncbi-near` rows `db.acbs.blastx.*` and
    `db.hitlist.*` of `LOSAT/tests/db_search_manifest/`). Ungapped searches and TBLASTX keep the
    size itself, in NCBI as in LOSAT. With `-best_hit_overhang` or `-best_hit_score_edge` the
    best-hit writer holds `2 * prelim_hitlist_size` HSPs of a query before it rebuilds its list
    (`max_hsps`, `core/hspfilter_besthit.c:79`), also in `Int4`, which overflows from a smaller size:
    1073741774 for a gapped search without composition-based statistics, 536870887 with them
    (`2*(2*h+50)`) and 2^30 for an ungapped search. NCBI BLAST+ 2.17.0 then rebuilds the list from
    its first HSP on and keeps HSPs that the best-hit test removes at the size below (line O1c's
    probe: BLASTN 69 lines where 1073741773 gives 56, ungapped BLASTN 69 where 1073741823 gives 57,
    BLASTP and PSIBLAST 11 where 536870886 gives 10; line O6a-3's probe: BLASTX 85 where 536870886
    gives 72, TBLASTN 62 where 536870886 gives 49, TBLASTX 252 where 1073741823 gives 214; each a
    superset). LOSAT computes that count in
    64 bits too (`core/hspfilter_besthit.rs`, since line O6a), so with these options such a size
    gives NCBI's output at 1073741773, 536870886 or 1073741823 (the rows `hitlist.*.bh*` with those
    contracts). Every smaller size gives NCBI's bytes. Before version 1.66 LOSAT
    reproduced the wrapped sizes (rule 2) and rejected those on which NCBI crashes (rule 3).

## Approved exception (TBLASTX, a query batch whose average query length is 0, against one subject, version 1.68)

Decided by the Owner on 2026-11-21 (「O2b-Q1=A」, line NP item O2b's question Q1, line NP
`DECISIONS.md` 42), under the Owner's rule of 2026-10-10 「Int4は気にせず続けてください」. Rule 1
asks for a run that NCBI completes to check the result; none gives LOSAT's rows here, and the Owner
chose this on the question's statement that no NCBI result can check them. The last paragraph of
the item was found after the decision, for this record, and was put to the Owner in the pull
request that added it.

20. **A TBLASTX query batch whose average query length is 0, linked against one subject.** The
    linking of TBLASTX's sum statistics takes one query length per query batch: the end of the last
    context in the batch's concatenated query (its letters and one separator after each context
    that has letters) divided by the number of contexts, all six of every query
    (`blast_parameters.c:1023-1026`). With `-strand plus` or `minus` the contexts of the other
    strand have no letters, so a batch of one 1,200-nt query and 300 queries of 3 nt averages 0, as
    one with 1,300 such queries does with both strands. NCBI's expected length is then
    `(Int4)BLAST_Nint(log(0) / H)`, the low 32 bits of the `long` of minus infinity, 0 on x86-64
    (`blast_parameters.c:1036-1037`). When the subjects together hold more protein letters than
    the subject, NCBI completes, and LOSAT gives its bytes (rule 2; the cases `strand.plus_avg0`
    and `strand.minus_avg0` of `LOSAT/tests/tblastx_regression_fixtures.py`). With one subject
    (`-subject` of one sequence, or a database of one), the database's length divided by 3 is the
    subject's, so the linking's `y_variable` and `x_variable` are 0 and its big-gap cutoff is
    `(Int4)floor(log(0) / Lambda) + 1`, `INT_MIN + 1` (`blast_parameters.c:1046-1076`). The
    linking's `Int4` sums with `cutoff[1]` (`link_hsps.c:868`, `942-943`) then overflow, and NCBI
    BLAST+ 2.17.0 stops with SIGSEGV in `BLAST_LinkHsps` (exit status 139, no output; with
    `-subject` and with a one-sequence `-db`). LOSAT computes the same cutoff and adds those sums
    in 64 bits. With a subject of up to 60,002 bases the linking keeps big gaps only, so each pass
    takes the chain of the most HSPs and, among those, of the highest score: the limit of NCBI's
    choice as the cutoff falls without bound. With a longer subject the search space passes
    8 x 50^2 and each pass chooses between that chain and the small-gap chain by their e-values,
    as NCBI does (`blast_parameters.c:1062-1070`). The default linking kernel and the literal port
    of NCBI's (`LOSAT_LINKING_LEGACY`) give the same rows at 1, 4 and 8 threads, and `-subject` and
    `-db` give the same rows. `LOSAT/tests/tblastx_avg0_one_subject.rs` keeps them against the
    3,000-nt subject `Sa0mut`: 15 rows with `-strand plus`, 25 with `minus`, 40 with both strands
    and 1,300 queries of 3 nt; and against `Sa0mut` followed by 57,003 N (60,003 nt), 13 rows with
    `-strand minus`.

    Found after the decision: NCBI divides the database's and the subject's lengths by 3 before it
    compares them, so it takes the same branch, and crashes, when the other subjects add fewer than
    3 bases in all (more exactly, when the subject's length modulo 3 and their lengths add up to
    less than 3). With `Sa0mut` and a second subject of 2 nt, `-strand plus`, NCBI exits with
    status 139 and LOSAT writes the 15 alignments above with the e-values of two subjects, which the
    same test keeps. NCBI's nearest runs that complete add a second subject of 12 bases (4 protein
    letters) of N, which has no hits, with `-dbsize 3000` (the cases
    `strand.plus_avg0_pad_subject`, `strand.minus_avg0_pad_subject` and
    `strand.both_avg0_pad_subject`, where LOSAT gives NCBI's bytes): their big-gap cutoff is -4. A
    second subject of 3 N gives the lowest cutoff that a completing run reaches, -8, and the same
    bytes. They write 15, 20 and 35 rows, each of whose alignments is one of LOSAT's rows above, but
    LOSAT's linking joins 5 more weak HSPs (15 to 18 bits) into sets with `-strand minus` and with
    both strands, so the sets and their e-values differ: with `-searchsp 1000000` in both runs, 6
    of NCBI's 9 rows, 3 of 17 and 9 of 26 are LOSAT's byte for byte. Until PR #1401 LOSAT raised an
    average of 0 to 1 and wrote the rows of other cutoffs here.

## Approved exception 22 (`losat windowmasker`, an empty id in a list matched as text, version 1.75)

Decided by the Owner on 2026-11-22 22:3x JST: the answer 「Q8=B」 (the whole answer 「Q1=B Q2=B Q3=A
Q4=A Q5=C Q6=B Q7=B Q8=B Q9=A Q10=A Q11=A Q12=A Q13=OK」) to question Q8 of that day's batch, line NP
item M5's question Q10, whose option B read 「同じ意味の警告を LOSAT の文言で出し、その id を飛ばして続ける」
(post a warning of the same meaning in LOSAT's words, skip that id and go on). Not taken: option A,
the explicit rejection of such lists (rule 3) that LOSAT made from version 1.72 until this version,
which stopped LOSAT on an input that NCBI's own error says it ignores. Owner-delegated (standing
instruction 2026-09-29, restated 2026-10-07): the post is NCBI's line through LOSAT's one mapping of
NCBI's words (`crate::wording`, as windowmasker's `can't understand id` post of `-text_match false`),
at NCBI's severity `Error`; not taken: a sentence of LOSAT's own or the severity `Warning`, which
would give the same skip two texts in one program and a line that NCBI never writes. Line NP item
M5.

22. **`losat windowmasker -ids` or `-exclude_ids`, matched as text: a line that gives an empty id.**
    `CWinMaskConfig::FillIdList` takes from each line that is not empty the text from its start
    (after a `>` at 0) to its first space or tab (`win_mask_config.cpp:467-482`), and
    `CIdSet_TextMatch::insert` stores it by its number of `|`-words (`win_mask_util.cpp:75-123`). An
    id of no words, which is the empty id (a line that starts with a space or a tab, or whose first
    word, up to a space, a tab or the line's end, is `>`), `|` (a first word `|` or `>|`), `>` (a
    first word `>>`) or `>|` (a first word `>>|`), makes `insert` post `Error:
    CWinMaskConfig::CIdSet_TextMatch::insert(): bad id: <id>: ignoring` and then store the id at
    `nword_sets_[nwords - 1]` with `nwords` 0: the `Uint4` index wraps to 4294967295, far past the
    end of the vector, and for the empty id the check of its last letter reads before the string
    (undefined behaviour in C++). NCBI BLAST+ 2.17.0 then stops with SIGSEGV (exit status 139, no output), in
    masking (`-ustat`) and in `-mk_counts`, for either list. A list matched as Seq-ids (`-text_match
    false`, `CIdSet_SeqId`) is not concerned: NCBI's parse posts such an id as one it cannot read
    and skips it, as LOSAT's does. LOSAT posts NCBI's error for each such id as it reads the list,
    where NCBI posts it, and skips the id, as the error says. The resize
    to the id's number of words that comes before NCBI's store changes nothing for 0 words, so the
    list is the one that NCBI reads from the file without those lines, the near input, which NCBI
    completes: LOSAT's output and exit status are NCBI's for the near input, and its standard error
    is NCBI's for the near input after those errors, in the order of their lines
    (`windowmasker::ids::EmptyId`).
    Evidence (NCBI BLAST+ 2.17.0 on ws, 2026-11-22; line NP's task folder `M/m5`, `SURVEY.md` §6.4
    and question Q10, probes b3 `ids_gt_mask` and `ids_pipe_mask` of 2026-11-20): fixture rows
    `e_ids_<form>` of `LOSAT/tests/windowmasker/cases.tsv`, whose lists `in/ids_e_<form>.txt` are
    `in/ids_a.txt` with one such line of each form after its second line (`>`, `> y`, a tab, a space,
    `|`, `>|`, `>>`, `>>|`), and `e_excl_gt` and `e_k_gt` (`-exclude_ids`, `-mk_counts`): NCBI's crash
    as frozen (exit status 139, no output, its one error). `e_ids_all`, `e_excl_all` and `e_k_all`
    hold the eight lines between the lines of `in/ids_a.txt`, `e_ids_only` the eight lines alone;
    their near rows `i_ids_a`, `i_excl_a`, `k_ids` and `i_ids_emptyfile` (the empty list
    `in/ids_empty.txt`) are NCBI's runs that complete. `LOSAT/tests/masker_fixtures.rs`
    (`EMPTY_ID_ROWS`) expects of each row the near row's bytes after the errors that NCBI posts in the
    rows of one line, through the command line and the library; the test
    `the_empty_id_rows_are_ncbis_crash_and_near_input` checks that each row of one line is NCBI's
    crash and that each row's list is its near row's with those lines added. Unit test
    `ids_of_no_words_are_skipped_with_ncbis_error` (`LOSAT/src/windowmasker/ids_tests.rs`).

## Approved exception 23 (BLASTX `-use_sw_tback`, a heap test past the batch's heaps in e-value mode, version 1.76)

Version 1.76, line RG item X4b-2 (line RG `DECISIONS.md` X4b-9). On 2026-11-21 (21:1x JST) the
Owner answered line RG's question Q-HEAP with option B (「O2b-Q1=A Q-HEAP=B Q-CAP=A」): a heap test
that NCBI makes past its array of the batch's heaps tests an empty heap. Until this version
LOSAT refused such a batch under rule 3 (line RG `DECISIONS.md` X4b-0 (a), 2026-11-20). Not taken:
keeping that refusal (option A; with one query in a batch every hit outside frame +1 reaches the
test, so most searches with hits in this mode were refused). Numbered after exception 22 (version
1.75).

23. **BLASTX `-use_sw_tback` in e-value mode: the heap test past the batch's heaps.** In e-value
    mode (`-sum_stats false`, or a `-max_intron_length` of 1 to 4, `blast_parameters.c:793-805`)
    NCBI's Smith-Waterman redo tests the first significant window of each query context of a match
    against the heap of redone matches that the context indexes (`redo_alignment.c:1465-1474`), in an array of one heap per query
    of the batch (`blast_kappa.c:3282-3296`). A BLASTX query has six contexts, 6 x q + f for the
    batch's query q and its frame f (0 to 5 for frames +1, +2, +3, -1, -2, -3). Below the batch's
    query count the index names another query's heap, which LOSAT reproduces (rule 2, the 1.76
    entry of the rule-3 list above); at or above it the index lies past the array's end, in NCBI's
    other allocations. With one query in a batch (a search of a single query, as most are, the
    last query of a file after a full batch, or a query that fills BLASTX's batch of 10,002 letters
    by itself) that is every hit outside frame +1. There NCBI BLAST+ 2.17.0
    crashes (SIGSEGV; at every thread count above 1 in our runs), fails the batch (`Traceback
    failed`, exit status 3), or ends with rows that depend on what lies there: in the runs of the
    evidence below one query with a hit in frame -1 wrote 2 rows in one run and 3 in others, and
    a query that NCBI crashed on wrote 6 rows with `MALLOC_PERTURB_=165`. LOSAT tests such a
    window against an empty heap, which takes every window: an empty heap holds fewer matches than
    its hit list size (`compo_heap.c:252-275`). One function owns the lookup of a heap by context,
    `CompoHeap::would_insert_at` (`LOSAT/src/algorithm/blastx/results.rs`). That is NCBI's run as
    is with an empty heap's answer there, and NCBI's own runs check it: in NCBI's search of the
    same batch followed by 5 x N queries without a hit (N the batch's query
    count; each `CCCCGGGGAACC` twelve times), every context of the batch's queries indexes a heap
    in the array, and the heap of a query without a hit is empty. The added queries write no row
    in outfmt 6, and no byte in a batch that fails. In all 9 cases of the evidence in which NCBI's
    padded run ends (exit status 0 or 3), LOSAT's run of the input as is equals it, at 1, 4 and 8
    threads; where NCBI's run as is ends with the same rows in every run (exit status 3 twice,
    exit status 0 once), the padded run equals it. The redo's early termination
    (`redo_alignment.c:1560-1580`, `compo_heap.c:405-409`) reads only the batch's own heaps, inside
    the array, and LOSAT keeps it; the padded run's empty heaps would turn it off, so the padded run
    is the check where the batch's heaps do not all fill (every row of the evidence has the default
    hit list size, 500, which none fills). Nor is a batch of one query that NCBI splits for its
    preliminary search (19,410 letters or more: `split_query_aux_priv.cpp:73-95` splits BLASTX only
    for a batch of one query, in chunks of 10,002 letters that overlap by 297, `local_blast.cpp:54-99`)
    split with the pads, so the padded run checks only shorter queries (the longest of the evidence
    has 5,600 letters); LOSAT splits such a query as NCBI does as is. Where the padded batch reaches the clean-up crash of
    1.30's class (its last redone match fails after its first query kept one), NCBI crashes with the
    pads too, and LOSAT refuses the batch (rule 3, the 1.76 entry above).
    Evidence: the `ncbi` rows `db.x4b2.past.*.pad.*` and `db.x4b2.bh.false.pad.6` of
    `LOSAT/tests/db_search_manifest/blastx.x4b2.tsv` hold NCBI's runs of the padded inputs
    (`LOSAT/tests/fasta/x4b2/*_pad.fna`, `LOSAT/tests/db_search/x4b2/make_x4b2.py`), and the rows
    `db.x4b2.past.*` without `.pad` and `db.x4b2.bh.false.6` take them as the expectation of the
    inputs as is (`ncbi-near:`,
    `LOSAT/tests/db_search_fixtures.py`): one query with hits in frames -1, -2 and -3, with
    `-subject` and `-db`, two one-query batches of `BATCH_SIZE=1`, NCBI's default batch of two
    queries of 5,000 and 5,600 letters, batches that fail, and the three queries of the best-hit
    inputs with `-best_hit_overhang` against `-db bh3p`, on which NCBI crashes as is. Unit tests
    `a_context_past_the_batch_tests_an_empty_heap` (`LOSAT/src/algorithm/blastx/results_boundary_tests.rs`) and
    `a_heap_test_past_the_batch_reads_an_empty_heap` (`LOSAT/tests/blastx_sw_tback.rs`).

## Explicit rejections of deterministic results (rule 2, version 1.10)

On a database with parsed Seq-ids NCBI gives these results without failing, and they look wrong
or need NCBI's object manager to reproduce. LOSAT refuses them instead (line DBID
`DECISIONS.md` 9, 20 and 21, line DI-7 for the query's id), under the maintainer's standing
delegation, which the maintainer accepted on 2026-11-02 with version 1.10. Each is a row of the table
in [`PD-LOSAT-DATABASE-SEARCH`](PD-LOSAT-DATABASE-SEARCH.md), which also lists refusals that rule 4
would port. The last three entries are not of these. The last (version 1.60, line RG item D5a)
and the one before it (version 1.54, line RG item X5) say so themselves. The one before those (version 1.45, line RG item D3-1) is not either: it concerns a PSSM
query's Seq-ids, not a database; it is Owner-delegated (standing instruction 2026-09-29, restated
2026-10-07), not under the acceptance of 2026-11-02; and it has no row in PD-LOSAT-DATABASE-SEARCH.
Version 1.45's rejection of well-formed BitMagic data was removed in version 1.58: line RG item
Q-BM decodes it as NCBI does, as the Owner decided on 2026-11-16 (「推奨通り」, Q-BM):

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
- BLASTP unified P-values (version 1.54, 2026-11-20, line RG item X5, DECISIONS X5-2,
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; not under the acceptance
  of 2026-11-02 and without a row in PD-LOSAT-DATABASE-SEARCH): a query batch in which a query
  with a residue to score is followed by a query without a search space, that is, a query that is
  not valid (all `X`, empty, or masked whole) when `-searchsp` is not given or is 0. NCBI's next
  context then holds an effective search space of 0 (`blast_setup.c:773-784`, `:846-847`), the
  scale of the adjustment is infinite, and NCBI reports every e-value of the query before it as
  0 (case F of `docs/evidence/unified_p_last_query/`: all 805 rows of the query, `-evalue 1000`);
  an alignment whose e-value was 0 before the step gets a NaN, which NCBI's sorts and cut-offs
  then compare. NCBI's bytes are deterministic, but the result has no use, and the order of a list with
  NaNs follows the C library's sort. LOSAT refuses the batch (exit 1, "unified P-values ... are
  not supported by LOSAT's BLASTP for query ..., which is followed in its query batch by a query
  without a search space ...", naming both queries), whether or not the query has a match. It
  finds the batch once the queries are set up, before the redo, and the error comes after the
  reports of the query batches before (the form of 1.40) and before any byte of the batch (a
  search that could meet it is left to the search of all queries, which writes those reports,
  not to the search in blocks of queries). With a `-searchsp` above 0 the next context's search space is that value
  and LOSAT writes NCBI's bytes (rule 2 above).

- `losat makembindex` with an index volume of more than 4 GiB of packed sequences (version 1.60,
  2026-11-21, line RG item D5a, Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07; not of the database class above): NCBI keeps the sequences' start positions, the
  volume's total and the size words in its 32-bit `TWord` and exits 0 with an index whose recorded
  sizes wrap and no longer match its data. That needs a `-volsize` of about 18,000 or more and more
  than 16 Gbp of sequences in one index volume; LOSAT refuses (exit 1) at the sequence that passes
  4 GiB.

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
- psiblast and tblastn `-in_pssm`, the BitMagic data of a Seq-table of the PSSM query's
  annotations (version 1.58, 2026-11-20, line RG item Q-BM; the Owner decided on 2026-11-16
  (「推奨通り」, Q-BM) to port it): `bit-bvector` values and a `bit-set-bvector`'s sparse index are
  decoded where NCBI's object manager decodes them, at the first value read
  (`CBVector_data::GetBitVector`: `bm::deserialize` of BitMagic 8.0.1 (Apache-2.0, bundled in NCBI BLAST+) in
  NCBI's 32-bit build, ported in LOSAT with its notice and licence text in
  THIRD_PARTY_NOTICES, section 10 (BitMagic)), with
  the tokens in NCBI's order and NCBI's checks: its `std::logic_error` `BM::Invalid serialization
  format` ends the run with `Error: BM::Invalid serialization format` and exit status 255 where
  the index reads the value, and is the `NCBI_CATCH` message of the `Sorted, max length` and
  `Seq-table location` columns (`Error: [<program>] [Bad '<column>' column] Exception: ...`),
  after which the run goes on. A stream that reaches one of the events of the rule-3 rejection
  (above) is refused at that event; a GAP block that is not one is decoded, and its rows answered,
  where NCBI's scalar, SSE2 and SSE4.2 builds agree and refused (rule 3) where they differ.
  Version 1.45 refused all such data.
