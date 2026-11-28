# Product Decision: NCBI BLAST+ behaviour that is a defect

- Decision ID: `PD-LOSAT-NCBI-DEFECTS`
- Version: 1.127
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
  line NP's 1.69);
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
  buffer for most lengths and which it answers for one less than the letters, an answer that LOSAT gives
  since line PS item PS7c; Karlin-Altschul values that are not finite numbers are searched as NCBI
  does since PS7c; Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07; numbered after 1.76; first 1.71, then 1.72, 1.74, 1.76 and 1.77 at the merges of dev that brought line DBID item
  pfoff-3's 1.71, line NP items M5a and M5b's 1.72 and line RG item X4b-1's 1.73, then line RG item D5b-3's 1.74 and
  line NP item M5's 1.75, then line RG item X4b-2's 1.76);
  1.78 2026-11-23 (approved exception 24, the Owner's answer of 2026-11-22 「Q5=C」 to line RG's
  question Q-HARDIDX: `losat blastn -use_index true` with `-db_hard_mask` searches the MegaBLAST
  index with the masked letters matching nothing and each seed in the hard range that holds it,
  where NCBI BLAST+ 2.17.0 takes for each hard chunk the index's list of the chunk with the same
  number; the refusal that LOSAT made from line RG item D5b-2 until this version removed; line
  RG item D5b-5; first 1.77, then 1.78 at the merge of dev that brought line RG item D3c-3a's 1.77);
  1.79 2026-11-23 (TBLASTN `-ungapped` runs, line NP item O4c: the explicit rejection of 1.46
  extended to the ungapped re-evaluation of an alignment on a subject of fewer than 12 bases,
  which reads the same translation that NCBI never makes, rule 3; the explicit rejection of 1.76,
  `-ungapped -comp_based_stats 0 -use_sw_tback`, on which NCBI crashes at every query batch it
  searches, extended from BLASTX to TBLASTN, rule 3; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.78; first 1.72, then 1.78 and 1.79 at the merges of
  dev that brought 1.72 to 1.77 (line NP items M5a, M5b and M5, line RG items X4b-1, D5b-3, X4b-2
  and D3c-3a), then line RG item D5b-5's 1.78);
  1.80 2026-11-21 (RPSTBLASTN `-use_sw_tback`, line PS item PS12c: the explicit rejection of 1.30
  of a Smith-Waterman redo that NCBI repeats without end extended from RPSBLAST to RPSTBLASTN, rule
  3, within the class of 1.30; Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07; numbered after 1.79; first 1.68, then 1.70, 1.72, 1.74, 1.77, 1.79 and 1.80 at the merges
  of dev that brought line NP's 1.68 and 1.69, lines RG's 1.70 and DBID's 1.71, line NP's 1.72
  and line RG's 1.73, then 1.74 to 1.76, then 1.77 and 1.78, then 1.79; with it, under
  rule 2, a Smith-Waterman alignment that starts with a gap, refused until line PS item PS12d);
  1.81 2026-11-23 (makeblastdb `-gi_mask`, line NP item M4c-1: exception 25, the alias file of a
  database of one volume, which names the volume where NCBI's names the database itself, the
  Owner's decision of 2026-11-20 (Q6=B); the explicit rejection, rule 2, of two masking
  algorithms that would write one GI-mask data file at the same time, whose bytes depend on how
  the C++ runtime buffers two streams on one file, the Owner's decision of 2026-11-21 23:07 JST
  (question Q-GIDUP, 「推奨通り」); the GI-mask writer's other defects reproduced under rule 2, two
  algorithms of one name whose files are written one after the other included; the form of
  exception 25 Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered
  after 1.80, line PS item PS12c; first 1.60, then 1.61, 1.63, 1.68, 1.70, 1.72, 1.74, 1.76 and 1.81 at the
  merges of dev that brought line RG's 1.60, then line RG's 1.61 and line LIB-O's 1.62, then the
  1.63 to 1.67 of lines RG, PS, NP and LIB-S, then line NP's 1.68 and 1.69, then line RG's 1.70
  and line DBID's 1.71, then line NP's 1.72 and line RG's 1.73, then line RG's 1.74 and line
  NP's 1.75, then the 1.76 to 1.80 of lines RG, NP and PS; the exception numbered 18 at first,
  then 20, 22 and 23);
  1.82 2026-11-22 (approved exception 26, makeprofiledb's word table past `Int4`, psiblast line item
  PS10a: the scaled word score threshold and the neighbouring-word score sums of `-dbtype rps` computed
  in 64 bits; proposed under the standing instruction 2026-09-29, restated 2026-10-07, on the Owner's
  rules of 2026-10-10 「NCBI のバグを再現する必要はない」 and 2026-10-22 「広い型になおしましょう」, and accepted
  by the maintainer on 2026-11-23 by marking PR #1497 ready for review, as the PR asked; and the explicit
  rejections, rule 3, of `losat makeprofiledb`, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.81; first
  1.72, then 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.80, 1.81 and 1.82 at the merges of dev that brought line NP
  items M5a and M5b's 1.72, then line RG item X4b-1's 1.73, line RG item D5b-3's 1.74, line NP item M5's 1.75, line
  RG item X4b-2's 1.76, line RG item D3c-3a's 1.77, line RG item D5b-5's 1.78, line NP item O4c's 1.79, line PS
  item PS12c's 1.80 and line NP item M4c-1's 1.81, and after exception 25; the exception numbered 22 at first, then 24, 25 and 26);
  1.83 2026-11-23 (approved exception 3 extended to BLASTP and BLASTX with `-db_hard_mask`: NCBI's
  traceback reads one letter before and one letter past the hard-masked copy of a subject that has
  masks, a copy that holds no sentinel byte, and ends with a crash (SIGSEGV) where the byte past the
  copy happens to equal the fence byte twice; LOSAT reads the bytes past both ends of that copy as
  the NULLBs of the mapped volume and writes the rows of NCBI's runs that do not meet the fence byte;
  it changes no code and no output of LOSAT, the behaviour being on dev already; the no-effect entry
  of `hard_ranges[1]` generalised to `hard_ranges[n]` after the last chunk; decided by the Owner on
  2026-11-23 13:49 JST, question Q-ALIGNEX, words 「推奨通り」, option A (recommended); line NP item
  M4b-AX; numbered after 1.82);
  1.84 2026-11-23 (the rule-3 rejection of a query letter of `losat makeprofiledb` narrowed, psiblast line
  item PS10c: NCBI's conversion gives every byte of the query a code, 255 for one that is no residue, so a
  lower-case IUPACaa query is built as NCBI builds it, and only a code past the 28 residues is refused,
  where NCBI's engine counts it, after the errors of the lengths and the matrix that NCBI gives first;
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.83; first 1.79, then
  1.83 at the merge of dev that brought 1.78 to 1.82 (psiblast line item PS10a, line RG item D5b-5, line NP items
  O4c and M4c-1, psiblast line item PS12c), then 1.84 at the merge of dev that brought line NP item M4b-AX's 1.83);
  1.85 2026-11-22 (alias files with `GILIST`, `TILIST` and `SEQIDLIST`, line RG item D1c-1: an
  explicit rejection of a `SEQIDLIST` token that NCBI warns about, kept under rule 2 by the
  maintainer (Owner's answer 「Q-D1c-R2=A」, 2026-11-22 22:4x JST); numbered after 1.71; first
  1.72, then 1.74, 1.78, 1.81, 1.82, 1.83 and 1.85 at the merges of dev that brought 1.72 and 1.73,
  then 1.74 to 1.77, then 1.78 to 1.80, then 1.81, then 1.82, then 1.83 and 1.84);
  1.86 2026-11-23 (two rule-3 rejections of `losat makeprofiledb`, psiblast line item PS10d: a
  query in a nucleotide coding whose data is shorter than NCBI's writer reads for its length, and
  in a profile with frequency ratios only a query whose `Seq-data` is a gap; a query in `iupacna`,
  `ncbi2na` or `ncbi4na` of enough bytes is built as NCBI builds it, and in a profile with
  frequency ratios only the query of every other coding gets NCBI's conversion: `Unknown
  conversion.` for a nucleotide coding, NCBI8aa's bytes as its codes, the unknown exception for
  NCBIpna and NCBIpaa; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07;
  numbered after 1.85; first 1.80, then 1.84 at the merge of psiblast line item PS10c, which brought
  dev's 1.80 to 1.82 and PS10c's 1.83, then 1.85 at the merge of dev that brought line NP item
  M4b-AX's 1.83 and PS10c's 1.84, then 1.86 at the merge of dev that brought line RG item D1c-1's 1.85);
  1.87 2026-11-23 (BLASTP `-ungapped -comp_based_stats 0 -use_sw_tback`, line NP item O4b2-2: the explicit
  rejection of decision D3 comes at the traceback of the first query batch that NCBI searches, after the
  reports of the batches before, the form of 1.76 and 1.79, instead of before the first byte; a run that
  searches no batch ends as NCBI's does; no new exception, rule 3; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.86; first 1.83, then 1.84 at the merge of dev that brought line NP
  item M4b-AX's 1.83, then 1.86 at the merge of dev that brought 1.84 and 1.85, then 1.87 at the merge of dev
  that brought psiblast line item PS10d's 1.86);
  1.88 2026-11-23 (approved exception 27, the Owner's answer of 2026-11-22 22:34 JST to question Q4
  of line NP's M4b (question Q6 of that day's batch, 「Q1=B Q2=B Q3=A Q4=A Q5=C Q6=B Q7=B Q8=B Q9=A
  Q10=A Q11=A Q12=A Q13=OK」): with `-db_hard_mask`, TBLASTN and TBLASTX search a reading frame whose
  masks leave it an empty first part (a mask within the first 2 bases of a sequence or within the 1
  or 2 bases after its last whole codon, or an empty mask at 0 before a tail mask), on which NCBI
  crashes in every run, as 「mask されていないものとして普通に探す」; the refusal of that shape
  removed; how the 1 or 2 masked bases are read (masked in every frame): the Owner's answer
  「推奨通り」 of 2026-11-23 13:49 JST, on the board, to the report of line NP's session ad-M4b-8 that
  asked Q-ALIGNEX and M4b-3c's question Q6 (recommended A and B), which is Q6 = B; line NP item
  M4b-3c; numbered after 1.87 (1.83 is line NP item M4b-AX, 1.84 psiblast line item PS10c, 1.85 line RG
  item D1c-1, 1.86 psiblast line item PS10d, 1.87 line NP item O4b2-2); first 1.77, then 1.85 at the merge of dev that brought its 1.75 to
  1.82 (the branch had assumed 1.76 for line NP item M4b-3b; dev's 1.76 and 1.77 are line RG's X4b-2 and
  D3c-3a), then 1.86 at the merge of dev that brought line RG item D1c-1's 1.85, then 1.87 at the merge of dev
  that brought PS10d's 1.86, then 1.88 at the merge of dev that brought O4b2-2's 1.87; the exception numbered 24 at first, then 27 after dev's 24 to 26);
  1.89 2026-11-23 (TBLASTX and TBLASTN `-db_soft_mask` with masks stored out of order on a frame
  of more than 5,000,000 residues, line NP item M4b-T2: the scan of a chunk's later range that
  starts before the chunk, and of the ranges after a first range raised to one word on the
  frame's last letter, ported, rule 4; a TBLASTN gapped search of a chunk that saves an initial
  alignment starting before the chunk, on which NCBI crashes or starts its gapped alignment
  outside the chunk, an explicit rejection, rule 3; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.88, line NP item M4b-3c; first 1.81, then 1.82 at the merge of dev
  that brought 1.81, then 1.84 after dev's 1.82 (line PS item PS10a) and M4b-AX's 1.83, then 1.86 after dev's 1.84
  (psiblast line item PS10c), then 1.87 after dev's 1.85 (line RG item D1c-1), then 1.88 after dev's 1.86
  (psiblast line item PS10d), then 1.89 after dev's 1.87 (line NP item O4b2-2) and M4b-3c's 1.88);
  1.90 2026-11-16 (the read of a header record by NCBI's object stream reproduced, rule 2: a
  definite length in the long form ends the value one byte per length byte after the first too
  early, and the bytes after the set are not read; line RG item T5d, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; written 2026-11-16 as 1.39, renumbered at the merges
  with dev to 1.44, 1.46, 1.47, 1.50, 1.51, 1.52, 1.68, 1.70, 1.71, 1.72, 1.73, 1.74, 1.75, 1.78, 1.81,
  1.82, 1.83, 1.85, 1.86, 1.89 and then 1.90, numbered after 1.89);
  1.91 2026-11-23 (RPSTBLASTN `-use_sw_tback`, line PS item PS12d: the Smith-Waterman alignment
  that starts with a gap in the profile's row, refused by 1.80, ported under rule 2; a start with
  a gap in the query's row, which no search was found to give, stays refused; Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.90; first 1.82, then 1.83
  at the merge of dev that brought psiblast line item PS10a's 1.82, then 1.89 at the merge of dev that
  brought 1.83 to 1.88 (line NP items M4b-AX, O4b2-2 and M4b-3c, psiblast line items PS10c and PS10d, line RG
  item D1c-1), then 1.91 at the merge of dev that brought line NP item M4b-T2's 1.89 and line RG
  item T5d's 1.90);
  1.92 2026-11-24 (BLASTP `-ungapped` under a negative seqid list, line NP item O4b2-4: the
  sequence that the list keeps while removing all its lines, of 1.34, stays in its hit list
  through the traceback's cut to the hit list size, its culling and best hit pipes and the prune
  and leaves at its Seq-align, as BLASTX's does, where 1.34 removed it before the traceback; its
  header is read only at the Seq-aligns (NCBI's UTF8String warning and a damaged header's error
  there); rule 2; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07);
  1.93 2026-11-24 (under a negative seqid list, a subject whose header NCBI's fetch cannot read
  ends the search only where the composition-based or Smith-Waterman redo fetches it, rule 2: its
  list stays in the redo's stream, line RG item D1c-1c; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; first 1.92, then 1.93 at the merge of dev that brought line NP
  item O4b2-4's 1.92);
  1.94 2026-11-23 (the readers of GI-indexed masks, line NP item M4c-2: NCBI's results on such
  databases reproduced, rule 2, in blastdbcmd, the five search programs, makeblastdb's copy and
  makembindex, with the search programs' GI-indexed masks on more than one thread, where NCBI
  stops or crashes in some runs, giving NCBI's `-num_threads 1` output as every `-db` search does
  (PD-LOSAT-DATABASE-SEARCH, the Owner's decision DB-3 of 2026-10-13); the `.nog` read as NCBI's
  mapping reads it, and the explicit rejections, rule 3, of an empty `.nog`, of a `.nog` entry
  past its mapping and of a `MASKLIST` name that is not UTF-8 where the masks are read (review
  round of 2026-11-24; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07);
  numbered after 1.93; first
  1.64, then 1.82 at the merge of line NP item M4c-1's final head, which brought dev's 1.64 to 1.80
  and its own 1.81, then 1.85 at the merge of dev that brought 1.82 to 1.84 (psiblast line items PS10a
  and PS10c, line NP item M4b-AX), then 1.91 at the merge of dev that brought 1.85 to 1.90 (line RG
  items D1c-1 and T5d, psiblast line item PS10d, line NP items O4b2-2, M4b-3c and M4b-T2), then
  1.92 at the merge of dev that brought psiblast line item PS12d's 1.91, then 1.93 at the merge of dev
  that brought line NP item O4b2-4's 1.92, then 1.94 at the merge of dev e89f6237c that brought
  line RG item D1c-1c's 1.93);
  1.95 2026-11-22 (approved exception 29, BLASTP and BLASTX `-db_hard_mask` with overlapping masks
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
  M4b-3d (both ported in version 1.122, leaving only the final-chunk overlap read past the list of
  ranges as an explicit rejection under rule 3); TBLASTN and TBLASTX unchanged; line NP item M4b-3b, the Owner's decision of 2026-11-21
  (「すべて推奨通りで」, M4b, overlap = A); Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07, for the rejections; numbered after 1.94 as the next free numbers
  on dev when it merged dev (agreed with line NP item M5 on 2026-11-23); first 1.72, then 1.73 at the
  merge of dev that brought line NP items M5a and M5b's 1.72, then 1.74 after dev's 1.73, then 1.76 and
  exception 23 after dev's 1.74, then 1.78 and exception 24 after dev's 1.75 to 1.77 of line NP item M5
  (exception 22), line RG item X4b-2 (exception 23) and line RG item D3c-3a, then 1.82 and exception
  26 after dev's 1.78 to 1.81 of line RG item D5b-5 (exception 24), line NP item O4c, line PS item
  PS12c and line NP item M4c-1 (exception 25), then 1.95 and exception 29 after dev's 1.82 to 1.94 (exceptions 26 and 27; exception 28 held by
  line NP item M4a's mask depth pull request));
  1.96 2026-11-24 (BLAST XML, outfmt 5, of such a sequence is written as NCBI writes it: its
  `Hit` is the Seq-align id's FASTA form, `Unknown` and length 0, each of its HSPs its number,
  scores and E-value, the coordinates 0, the Seq-align's identities (-1 where it has none) and no
  letters, and the report goes on (exit 0); the refusal of 1.31 no longer covers outfmt 5 (outfmt
  17 alone, and in a call of LOSAT's library outfmt 5 with outfmt 0 to 4); review 1 of the item:
  such a sequence is refused, in every format that reads it (outfmt 0 to 4, 5, 6, 7, 10, 13 to 16,
  18 and 20), where the run loaded before it, in NCBI's order, another sequence that holds an id
  NCBI's scope takes for its Seq-align id, which NCBI writes for it without the lookup (H1, H2;
  this changes what 1.47, 1.50 and 1.71 wrote for outfmt 13 to 16, 6, 7, 10, 18, 20 and 0 to 4
  there: silent differences, and for outfmt 10 a wrong refusal), and a call of LOSAT's library is
  refused where NCBI's runs of its formats stop at different places after such a sequence (M1:
  a format that ends the run there with one whose run goes on, or two that end at different
  rows); review 2 of the item: NCBI gives a held sequence for such a Seq-align id from its data
  loader, by the id's own type, accession and version, for the whole run, and from its scope, by
  another type of the accession and version or a weak match (an accession without a version and
  its versions, a name), until the end of the query batch (of the run in outfmt 5), and LOSAT
  keeps the two apart (H1: a weak match is refused where outfmt 0, 5, 6, 7 and 18 wrote
  `Unknown`; M1: a sequence of another type loaded in an earlier batch, which this version
  refused at first, leaves it thrown outside outfmt 5, as NCBI does); such a refusal also takes
  534 runs of the review's probe whose output LOSAT wrote as NCBI does before, where NCBI's
  output with the holder equals the one without it (M2); and a call of the library is also
  refused where a format that reads nothing of the sequence, or one whose run goes on and posts
  nothing, comes with one whose run ends there, a pairwise report, or outfmt 15 or 16, which post
  the error on standard error (M3); review 3 of the item: NCBI's scope holds a sequence by the ids
  of the definition line that its load made its Bioseq of (the line of its own Seq-align id: its
  first line with a gi, else its first line), its data loader by the ids of every line (M1: a
  holder whose other line holds an id that the scope takes for the Seq-align id or matches weakly,
  `gi|7052|ref|NP_100001.1|` after `gi|7051|ref|XP_900001.1|`, leaves it thrown, as NCBI does,
  where this version refused it at first), and with `-html` outfmt 0 ends the run at the first
  description line of such a sequence (exit 2), so a holder after it is not refused (M2); at the
  merge of dev that brought line PS item PS15b-1 (psiblast's outfmt 5, 8, 9 and 11), psiblast's
  outfmt 5 of such a sequence is written as BLASTP's, and its PSSM, which reads each sequence that
  it includes by its Seq-align id after the round's reports, in the order of the hit list, refuses
  such a sequence where NCBI gives a held one for it (outfmt 8 wrote NCBI's `Failed to retrieve
  sequence` for a read that NCBI's scope held: probe of the merge, 6 runs); review 5 of the item:
  whether a query batch prefetches has one owner, of the limit, the batch's Seq-align count, the
  program and the first format (NCBI prefetches for outfmt 8, 9, 11 and 12 only at
  `PRE_FETCH_SEQS_LIMIT=2147483647`, and psiblast alone before its archive), so psiblast's outfmt
  8, 9 and 11 read one by one in NCBI's default environment and their PSSM reads are refused where
  NCBI gives a held sequence for them (H2 and H3: dev's #1538 took those batches for prefetched
  and wrote the warning, or `Unknown`, where NCBI reads the holder's letters, 548 runs of outfmt 8
  and 548 of outfmt 11 in the review's probe), psiblast's outfmt 11 at that limit stops where its
  prefetch throws (exit 2, H4), and the PSSM's loads write the taxonomy database's texts whatever
  the first format (outfmt 8, 9 and 11 with rounds wrote none where NCBI writes the error of a
  database with a wrong magic number), which a BLAST archive saves with the severity of the
  line's header (`code 2` for `Error: `); a held PSSM read is refused only where NCBI's read takes
  the holder's letters, where the subject range of an HSP that the PSSM includes ends within
  them (a shorter holder's read throws, and NCBI writes `Failed to retrieve sequence`, as LOSAT
  does: H1, M4), because LOSAT does not write a holder's letters into the PSSM (line DBID's next
  item, M3); and the Seq-annot and archive reports, which load nothing, are not loads before a
  PSSM read (M1); line
  DBID item pfoff-4, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered
  after 1.95; first 1.75, then 1.77, 1.79, 1.86, 1.90, 1.92, 1.93, 1.94, 1.95 and 1.96 at the merges of dev that brought
  line RG item D5b-3's 1.74, line NP item M5's 1.75 and line RG item X4b-2's 1.76, then line RG
  items D3c-3a's 1.77 and D5b-5's 1.78, then line NP item O4c's 1.79, line PS item PS12c's 1.80,
  line NP item M4c-1's 1.81, psiblast line item PS10a's 1.82, line NP item M4b-AX's 1.83, psiblast
  line item PS10c's 1.84 and line RG item D1c-1's 1.85, then psiblast line item PS10d's 1.86, line
  NP item O4b2-2's 1.87, line NP item M4b-3c's 1.88 and line NP item M4b-T2's 1.89, then line RG
  item T5d's 1.90 and line PS item PS12d's 1.91, then line NP item O4b2-4's 1.92, then line RG item
  D1c-1c's 1.93, then line NP item M4c-2's 1.94, then line NP item M4b-3b's 1.95);
  1.97 2026-11-18 (TBLASTN, line RG item X-short-bat: the rejections of a step where NCBI's search
  of a query batch crashes, does not finish or stops come after the reports of the query
  batches before, as NCBI's crash does and as TBLASTX's rejections do (the form of 1.40, the
  Owner's decision of 2026-11-16 on line RG's question Q-X2c), also D13's where NCBI stops (also
  on a database left without a subject, whose batch LOSAT reported before) and, with `-in_pssm`,
  that of a query whose letters give no Karlin-Altschul parameters, which came before the run's
  first byte; and moves D13's rejection of a negative `CHUNK_SIZE` whose wrapped ranges leave gaps
  between them (NCBI searches on) to the same place, where the search of the batch's first chunk
  with a valid context would read its subjects (after that chunk's setup checks and the reads of
  the database that fail before it; without a subject in hand nothing is lost and the batch is
  reported); Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07;
  numbered after 1.96, line NP item M4b-3b (1.95) and line DBID item pfoff-4 (1.96), at the merge of dev of 2026-11-24 that brought them;
  first 1.52, then 1.53 at the merge of dev of 2026-11-18, then 1.72 and 1.74 at the merges of dev
  of 2026-11-22, then 1.76, 1.77, 1.82, 1.86, 1.89 and 1.91 at the merges of dev of 2026-11-23, then 1.93 at
  the merge of dev of 2026-11-24 that brought 1.92, then 1.94 at the merge of dev of 2026-11-24
  that brought 1.93, then 1.95 at the merge of dev of 2026-11-24
  that brought 1.94, then 1.97 at the merge of dev of 2026-11-24
  that brought 1.95 and 1.96);
  1.98 2026-11-23 (BLASTP and psiblast, line RG item X-bpgap: D13's rejection of a negative
  `CHUNK_SIZE` whose wrapped ranges leave gaps between them moved to TBLASTN's step of 1.97, where
  the search of the batch's first chunk with a valid context would read subjects in hand (a batch
  without a subject in hand or without such a chunk is reported, as NCBI 2.17.0 reports it), and
  the stops of a split batch's chunks before the errors of the database and of the subject length
  that NCBI meets later, through one owner of the steps of a split batch that BLASTP, psiblast and
  TBLASTN share; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered
  after 1.97; first 1.85, then 1.87 at the merge of dev that brought 1.82 to 1.85, then 1.90 at the
  merge of dev that brought 1.86 to 1.88, then 1.92 at the merge of dev that brought 1.89 and 1.90,
  then 1.93 at the merge of dev that brought 1.91, line PS item PS12d, then 1.94 at the
  merge of dev of 2026-11-24 that brought 1.92, line NP item O4b2-4, then 1.95 at the merge of dev
  of 2026-11-24 that brought 1.93, line RG item D1c-1c, then 1.96 at the merge of dev
  of 2026-11-24 that brought 1.94, line NP item M4c-2, then 1.98 at the merge of dev
  of 2026-11-24 that brought 1.96, line NP item M4b-3b (1.95) and line DBID item pfoff-4 (1.96));
  1.99 2026-11-24 (approved exception 20 also for a short query against one short subject:
  TBLASTX's, TBLASTN's `-ungapped` and BLASTX's `-ungapped` linking against one subject whose
  expected length is 0 or less, where NCBI's big-gap cutoff is `INT_MIN + 1` as for a batch
  averaging 0 and NCBI crashes in its first linking call; LOSAT's sums in 64 bits, and the cutoff
  of a negative expected length taken as NCBI's `INT_MIN + 1` (LOSAT computed 1 there); its rows
  equal NCBI's own output for the nearest search it completes in every probed case; the Owner's
  answer of 2026-11-24 to question Q-XT1 of line RG item X-short-tbx-link, 「Q-XT1=A」; and the
  same for several subjects on the route of a negative expected length, the Owner's answer of
  2026-11-24 to question Q-XT2, 「A」; numbered after 1.98; first 1.97);
  1.100 2026-11-24 (FASTA records past NCBI's 32-bit `TSeqPos`, line RG item KC-11d: rule 3 for a
  nucleotide piece that NCBI's reader packs into 0 bytes and for `makembindex -old_style_index
  true`'s chunk end past 2^32; `losat windowmasker` reads records of 4,294,967,296 letters or more
  as NCBI's reader does instead of refusing them; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.99; first 1.94, then 1.97 at the merge of dev that brought
  line NP item M4c-2's 1.94, line NP item M4b-3b's 1.95 and line DBID item pfoff-4's 1.96, then 1.99 at
  the merge of dev that brought line RG items X-short-bat's 1.97 and X-bpgap's 1.98, then 1.100 at
  the merge of dev that brought line RG item X-short-tbx-link's 1.99);
  1.101 2026-11-24 (approved exception 30, the Owner's answer of 2026-11-24 12:44 JST 「Q-MIXIDX=C」 to
  line RG's question Q-MIXIDX: `losat blastn -use_index true` with a MegaBLAST index of only some
  database volumes searches each index volume that a query batch or query chunk reaches with that
  search's query, where NCBI BLAST+ 2.17.0 gives later searches the first search's used results of
  some volumes; the item of rule 2 and the explicit rejection of rule 3 of 1.74 removed; line RG
  item D5c; first 1.95, then 1.99 at the merge of dev of 2026-11-24 that brought 1.95 to 1.98, line
  NP item M4b-3b, line DBID item pfoff-4, line RG items X-short-bat and X-bpgap, then 1.100 at the
  merge of dev of 2026-11-24 that brought 1.99, line RG item X-short-tbx-link, then 1.101 at the
  merge of dev of 2026-11-24 that brought 1.100, line RG item KC-11d);
  1.102 2026-11-23 (rule-3 rejections of RPSBLAST and RPSTBLASTN at the traceback of a profile whose
  `.psq` letters differ from its `.rps` columns, psiblast line item PS10e: such a profile is searched
  with its `.pin` length, and a query batch whose hit makes NCBI's traceback leave the profile's data
  is refused; the refusal at the volume's open narrowed to offsets that do not ascend; these
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; and one rule-2 rejection by
  the maintainer's decision of 2026-11-24 (the answer 「Q-PS10e=A」, psiblast `DECISIONS.md` 103),
  RPSTBLASTN without composition-based statistics on an alignment that starts past such a profile's
  letters; numbered after 1.101; first 1.87, then 1.90 at the merge of dev
  that brought 1.86 to 1.89 (psiblast line item PS10d, line NP items O4b2-2, M4b-3c and M4b-T2), then 1.92
  at the merge of dev that brought line RG item T5d's 1.90 and line PS item PS12d's 1.91, then 1.93 at the merge of dev that
  brought line NP item O4b2-4's 1.92, then 1.97 at the merge of dev that brought line RG item
  D1c-1c's 1.93, line NP items M4c-2's 1.94 and M4b-3b's 1.95 and line DBID item pfoff-4's 1.96, then 1.99 at the merge of dev that brought line RG items
  X-short-bat's 1.97 and X-bpgap's 1.98, then 1.101 at the merge of dev that brought line RG items
  X-short-tbx-link's 1.99 and KC-11d's 1.100, then 1.102 at the merge of dev that brought line RG
  item D5c's 1.101);
  1.103 2026-11-24 (a `MASKLIST` name that is not UTF-8 is read by its bytes, as NCBI's
  `std::string` is: the rejection of 1.94 is withdrawn, and blastdbcmd, the search programs,
  makeblastdb's copy and makembindex read such GI-indexed masks as NCBI's do; line RG item K10a,
  round 2 of its review, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.102;
  first 1.97, then 1.99 at the merge of dev of 2026-11-24 that brought 1.97 and 1.98, then 1.100 at
  the merge of dev of 2026-11-24 that brought 1.99, line RG item X-short-tbx-link, then 1.101 at the
  merge of dev of 2026-11-25 that brought 1.100, line RG item KC-11d, then 1.102 at the merge of
  dev of 2026-11-25 that brought 1.101, line RG item D5c, then 1.103 at the merge of
  dev of 2026-11-25 that brought 1.102, psiblast line item PS10e);
  1.104 2026-11-25 (TBLASTN and TBLASTX `-db_soft_mask` with masks stored out of order on a frame of
  more than 5,000,000 residues, line NP item M4b-T3: a TBLASTN gapped search gives NCBI's rows for
  an initial alignment of at most 11 letters that starts before its chunk with its middle in the
  chunk, and refuses only those on which NCBI crashes, rule 3 narrowed; NCBI's one-thread diagonal
  offset and its clears taken from the database order for such sequences at every thread count,
  rule 4; line NP decision 58, Q1 = B and Q2 = B, 2026-11-23; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; numbered after 1.103);
  1.105 2026-11-25 (rule-3 rejections of `losat makeprofiledb` with ASN.1 text and `-binary` PSSM
  inputs and binary definition lines, psiblast line item PS10b: a query whose `iupacna` data holds
  a byte >= 0x80 in a pair position before any final odd overhang byte, where NCBI's
  `x_ConvertIupacnaTo4na` indexes its static table with a signed `char` and reads before the table;
  a `Seq-id` element stated as `NULL` inside binary definition lines (`ASN1_BlastDefLine`) when
  writing a version 5 database with no prior sequence written yet and no active taxid map (including
  `-taxid` or an empty `-taxid_map`), where NCBI's taxonomy writer aborts on close while
  `-taxid_map` or an `org` descriptor, `-blastdb_version 4`, `-index false` on the first slot, and
  version 5 after an earlier profile match NCBI; and `:OtherType` pointers in both ASN.1 text and
  `-binary` profiles via `Names::PROFILE` as well as binary definition lines, plus host-dependent
  1 GiB..2^47 lengths, whose NCBI behaviour depends on process registration history or host memory;
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.104 at
  the merge of dev that brought line NP item M4b-T3's 1.104; first 1.104);
  1.106 2026-11-25 (TBLASTN and TBLASTX `-db_hard_mask` on masks that overlap or nest in the order
  of their start, line NP item M4b-3c part 2: searched as NCBI searches them, rule 4, with NCBI's
  one-thread table where a sequence meets entries that no sequence wrote; the rule-3 rejections
  narrowed to masks out of the order of their start, NCBI's never-ending gapped search, a frame
  searched whole past overlapping masks and approved exception 29 extended to TBLASTN and TBLASTX
  (Owner's decision of 2026-11-25 17:54 JST on question Q-3C2-STALE = B of line NP's series M,
  decision 77, chosen option 「B: 例外 29 を tblastn/tblastx にも広げ、配列ごとに独立探索する（4 検索統一）」);
  series M decision 46, 2026-11-21, Owner's answer A on overlapping masks; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07, for the shape of the change; numbered after 1.105 at
  the merge of dev that brought psiblast line item PS10b's 1.105; first 1.105);
  1.107 2026-11-23 (makeblastdb `-mask_data`, line NP item M (b): exception 28, a mask file whose
  Seq-locs nest deeper than NCBI's stack holds, read to any depth, rule 1, as approved exception 15
  reads a PSSM; approved by the Owner's answer 「Q-MDEPTH=B」 of 2026-11-24 to question Q-MDEPTH; the
  refusal of nests deeper than 400 object frames removed; numbered after 1.106, the last on dev
  22ac798e5 (line NP item M4b-3c part 2); first 1.77 (exception
  24), then 1.81 (exception 25) at the merge of dev that brought 1.76 to 1.80 (line RG items X4b-2,
  D3c-3a and D5b-5, line NP item O4c, line PS item PS12c), then 1.92 (exception 28) at the merge of
  dev that brought 1.81 to 1.91 and exceptions 25 to 27, then 1.93 at the merge of dev that
  brought line NP item O4b2-4's 1.92, then 1.94, then 1.95 at the merge of dev that brought line RG
  item D1c-1c's 1.93, then 1.99 at the merge of dev 4626ddd5d that brought 1.94 to 1.98 (line NP
  items M4c-2 and M4b-3b, line DBID item pfoff-4, line RG items X-short-bat and X-bpgap), then
  1.104 at the merge of dev c2e1e4c36 that brought 1.99 to 1.103 (line RG items X-short-tbx-link, KC-11d,
  D5c and K10a, psiblast line item PS10e), then 1.105 at the merge of dev 163323c56 that brought
  1.104 (line NP item M4b-T3), then 1.106 at the merge of dev a5595ed11 that brought 1.105
  (psiblast line item PS10b), then 1.107 at the merge of dev 22ac798e5 that brought 1.106 (line NP
  item M4b-3c part 2));
  1.108 2026-11-24 (BLASTP, line RG item X-short-bat-2: the rejection of a compressed lookup table
  whose neighbouring words fill NCBI's overflow banks (decision D3) comes in the search of the query
  batch whose table fills them (of a split batch, its chunk's table), after the reports of the
  query batches before, as TBLASTN's of 1.97 (it came before the run's first byte); the three
  builders of the table name the first of NCBI's tables that fills them; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; numbered after 1.107; first 1.99, then 1.100 at
  the merge of dev of 2026-11-24 that brought 1.99, line RG item X-short-tbx-link, then 1.103 at
  the merge of dev of 2026-11-25 that brought 1.100 to 1.102, line RG items KC-11d and D5c and
  psiblast line item PS10e, then 1.104 at the merge of dev of 2026-11-25 that brought 1.103, line
  RG item K10a, then 1.106 at the merge of dev of 2026-11-25 that brought 1.104 and 1.105, line NP
  item M4b-T3 and psiblast line item PS10b, then 1.107 at the merge of dev of 2026-11-25 that
  brought 1.106, line NP item M4b-3c part 2, then 1.108 at the merge of dev of 2026-11-26 that
  brought 1.107, line NP item M (b));
  1.109 2026-11-25 (version 5 alias `SEQIDLIST` combined with a search `-seqidlist` or
  `-negative_seqidlist`, line RG item D1c-2: explicit rejection, rule 3, when the alias
  `SEQIDLIST` or the search's `-seqidlist` holds two differently spelled Seq-ids that NCBI's
  `s_CompareSeqId` treats as equivalent, or two unparseable Seq-ids whose `s_CompareSeqId`
  exceptions differ, or a Seq-id that posts a warning alongside an unparseable Seq-id, where
  NCBI's `std::sort` violates strict weak ordering or depends on the C++ library's sort steps;
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.108 at
  the merge of dev that brought 1.108 (line RG item X-short-bat-2); first 1.105, then 1.108);
  1.110 2026-11-24 (BLASTN and BLASTX, line RG item X-short-bat-2: their rejections inside the search
  of a query batch come after the reports of the query batches before, in NCBI's order: a query
  chunk without a query (in BLASTN also a negative `CHUNK_SIZE` whose first range covers the batch,
  in BLASTX every negative pair that splits a batch) and, in BLASTX, a chunk split again, at NCBI's
  split, before the errors of the database that NCBI meets in a chunk's setup; BLASTN's rejection
  of a negative pair whose chunks leave gaps (version 1.1) at TBLASTN's step of 1.97, where the
  search of the batch's first chunk with a valid context would read subjects in hand whose
  sequences can be read; BLASTX's rejection of an overlap that is not a multiple of 3 only with
  such subjects (a batch without them is reported, or fails, as NCBI's is); the fetch rejections
  under a negative seqid list of 1.93 and BLASTX's rejection of the compressed-lookup thresholds
  that exhaust NCBI's overflow banks (decision D3, LOSAT's rejection with exit status 1, no
  longer reported as an engine error) after the same reports; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; numbered after 1.109; first 1.99, then 1.109 at the merge of dev that brought 1.99 to 1.108,
  then 1.110 at the merge of dev of 2026-11-26 that brought 1.109, line RG item D1c-2);
  1.111 2026-11-24 (approved exception 31, the Spouge cutoff search past `Int4`, psiblast line item
  PS7c: psiblast and tblastn `-in_pssm` search for the score of the hit-saving cutoff
  (`BLAST_SpougeEtoS`) without the signed overflow of its `Int4` doubling and bisection; the Owner's
  answer of 2026-11-24, 12:44 JST, 「Q-PS7c-SPOUGE=B」; the explicit rejection of these Lambdas that the
  item's fix round 1 had built, never on dev, removed; numbered after 1.110; first 1.99, then 1.100 at
  the merge of dev of 2026-11-24 that brought 1.99, line RG item X-short-tbx-link, then 1.104 at the
  merge of dev of 2026-11-25 that brought 1.100, line RG item KC-11d, 1.101, line RG item D5c (which
  took the number 30 for its exception), 1.102, psiblast line item PS10e, and 1.103, line RG item
  K10a, then 1.105 at the merge of dev of 2026-11-25 that brought 1.104, line NP item M4b-T3, then
  1.106 at the merge of dev of 2026-11-25 that brought 1.105, psiblast line item PS10b, then 1.107 at
  the merge of dev of 2026-11-25 that brought 1.106, line NP item M4b-3c part 2, then 1.108 at the
  merge of dev of 2026-11-26 that brought 1.107, line NP item M (b), then 1.109 at the merge of dev of
  2026-11-26 that brought 1.108, line RG item X-short-bat-2, then 1.110 at the merge of dev of
  2026-11-26 that brought 1.109, line RG item D1c-2, then 1.111 at the merge of dev of
  2026-11-26 that brought 1.110, line RG item X-short-bat-2; this item's exception was 30 until
  the 1.104 merge);
  1.112 2026-11-25 (`losat makeblastdb -input_type asn1_txt|asn1_bin`, line RG item D3c-4: the
  explicit rejections, rule 3, of a nucleotide ASN.1 sequence whose `Seq-data` (`iupacna`,
  `ncbi2na` or `ncbi4na`) is shorter than `Seq-inst.length` asks for (including a `-hash_index`
  read of 1 to 3 packed tail bases past the data, class of 1.86), of `NULL` elements inside
  `Bioseq-set.seq-set`, `Bioseq.id`, `Seq-descr`, `Delta-ext`, `Blast-def-line-set`, or
  `Blast-def-line.seqid` and an empty `Delta-ext` list (class of 1.45 and 1.52), and of a
  punctuation-only Bioseq title that NCBI's `x_CleanAndCompress` reads past the end of (class of
  1.2 and 1.37); Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered
  after 1.111 at the merge of dev that brought line RG item X-short-bat-2's 1.110 (#1647) and
  psiblast line item PS7c's 1.111; first 1.105, then 1.106, 1.107, 1.108, 1.109 and 1.110);
  1.113 2026-11-26 (rule-3 rejections of `losat makeprofiledb -dbtype delta`, psiblast line item
  PS14: profiles stored by row (`byRow TRUE`) or with `numRows < 28` and non-empty
  `weightedResFreqsPerPos`, where NCBI's `x_CheckDelta` and `x_UpdateDelta` copy into an empty
  `std::list::begin()`; an empty `numIndeptObsr` list, where NCBI dereferences `end()`; a
  `weightedResFreqsPerPos` list whose length is not a multiple of 28 when `len / 28 ==
  numIndeptObsr.len()` and every full column of 28 sums within 0.0001 of 1.0, where `x_ValidateCd`
  increments past `end()`; `weightedResFreqsPerPos` or `numIndeptObsr` values whose scaled
  conversion to `Uint4` is out of range or `NaN`; and a missing or unwritable `TMPDIR` where
  NCBI's `CTmpFile` stops with a `mkstemp()` message containing its process ID; Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.112 at the merge of dev
  245b7e2e3; first 1.109, then 1.110);
  1.114 2026-11-25 (a Seq-align id whose lookup throws on the id itself, a string tag of `PIG` or
  `BL_ORD_ID` (`gnl|BL_ORD_ID|0`, which makeblastdb stores with a text tag) or a `ti` tag that is
  not a number, is written as NCBI writes it with the prefetch on or off: NCBI's prefetch of the
  batch catches the exception, writes `Error pre-fetching sequence data` and ends its loads, and
  the formats then load the batch's sequences one by one, as with the prefetch off (1.47, 1.50,
  1.71 and 1.96), with `Error: ` and the exception's text and exit status 255 where a format stops
  at it; the rejection of 1.10 stays for such a sequence's tabular row that prints `sacc` or
  `saccver` (NCBI's null pointer exception, with addresses of its binary), outfmt 17 and a holder
  that NCBI's scope or data loader gives for the id, and the batch's other sequences are refused as
  where the prefetch is off (a Seq-align id that an id of another sequence of the database matches
  weakly, or that a version 5 database cannot list; review 1 of item pf7, M1); line DBID item pfoff-7, Owner-delegated, standing instruction 2026-09-29,
  restated 2026-10-07; numbered after 1.113 (psiblast line item PS14); first 1.101, then 1.103 at the
  merge of dev that brought line RG item KC-11d's 1.100, exception 30's 1.101 and psiblast line
  item PS10e's 1.102, then 1.104 at the merge of dev that brought line RG item K10a's 1.103, then
  1.106 at the merge of dev that brought line NP item M4b-T3's 1.104 and psiblast line item PS10b's
  1.105, then 1.107 at the merge of dev that brought line NP item M4b-3c part 2's 1.106, then 1.108
  at the merge of dev that brought line NP item M (b)'s 1.107, then 1.109 at the merge of dev that
  brought line RG item X-short-bat-2's 1.108, then 1.113 at the merge of dev that brought 1.109 to
  1.112, then 1.114 at the merge of dev that brought psiblast line item PS14's 1.113);
  1.115 2026-11-26 (approved exception 32, megablast greedy gap costs above 32767 and up to
  `Int4` maximum, line RG item G1: the rule-3 explicit rejection of megablast gap costs above
  32767 is removed; affine greedy gap-cost and distance arithmetic are widened to 64 bits over
  non-empty distances so every gap cost accepted by the nucleotide Karlin-Altschul tables runs
  without NCBI's signed 32-bit overflow or dense-table memory exhaustion; approved by the Owner
  on 2026-11-26 at 15:41 JST (「承認 A」, choosing option A to accept megablast gap costs above
  32767 as an approved exception under the Owner's rules of 2026-10-10 「NCBI のバグを再現する
  必要はない」 and 2026-10-22 「広い型になおしましょう」); numbered after 1.114 at the merge of
  dev 42678f399 that brought line DBID item pfoff-7's 1.114; first 1.105, then 1.106, 1.107,
  1.108, 1.113 and 1.114; exception 31 until the 1.113 merge);
  1.116 2026-11-26 (BLASTN query chunks, line RG item X-chunk: unifies BLASTN's split-batch steps on
  `SplitBatchSteps` (`without_resplit_stop()`), keeps 64-bit wrapping `size_t` chunk-range
  arithmetic in `query_split.rs` across all three empty-chunk shapes (including exact-overlap
  `L <= overlap` and unwrapped negative pairs reaching or passing L before the last chunk), checks
  split chunks in duplicate-id batch setup with batch validity, shared chunk registry loads
  (`WindowMasker::load_posts`) and chunk masks, and clarifies in approved exception 1 that
  re-splitting chunks apply only when every chunk holds a query; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; numbered after 1.115 at the merge of dev that
  brought line RG item G1's 1.115; first 1.105, then 1.106, 1.107, 1.108, 1.113, 1.114 and 1.115);
  1.117 2026-11-26 (all programs, command line and library API, line WQA item D23, Owner decision
  `Q-D23=A` / `web-qa-20261119/DECISIONS.md` 49 of 2026-11-24: a `-logfile` naming an input file
  read by the run — including `-query`, `-subject`, `.ncbirc` / `ncbi.ini`, `<program>.ini`,
  `-conffile`, `-in`, `-ustat`, `-in_pssm`, `-in_msa`, `-gilist`, `-negative_gilist`, `-ipglist`,
  `-negative_ipglist`, `-seqidlist`, `-negative_seqidlist`, `-taxidlist`, `-negative_taxidlist`,
  `-window_masker_db`, `-window_masker_taxid` counts, `-import_search_strategy`, `-mask_data`,
  `-taxid_map`, `-entry_batch`, `-ipg_batch`, `-ids`, `-exclude_ids`, `-i`, `rpsbproc` data files,
  or a BLAST database alias, alias list, GI-mask, MegaBLAST index, or component file of `-db` /
  `-filtering_db`, via relative paths, symlinks, hard links, or missing inputs that `-logfile`
  would create — is refused before opening or creating the file, exit status 1, leaving existing
  inputs unchanged, matching LOSAT Web's D19 rule; numbered after 1.116 at the merge of dev
  f71a5cb85 that brought 1.116 (#1619); first 1.110, then 1.113, then 1.114, then 1.115, then 1.116);
  1.118 2026-11-24 (approved exception 33, the Owner's decision of 2026-11-20 on line RG's question
  F2-Q1, 「Q1=C」: an imported strategy's subjects searched once at any `-num_threads`; line RG item
  F2, `-import_search_strategy`; numbered after 1.117; first 1.70, taking exception 22 at the merge
  of dev 540759a2f, then 1.71 at the merge of dev 194053478 that brought line RG item D1b-2's 1.70,
  then 1.72 at the merge of dev 1b6e05a0d that brought line DBID item pfoff-3's 1.71, then 1.74 at
  the merge of dev 2e9653a71 that brought line NP items M5a and M5b's 1.72 and line RG item X4b-1's
  1.73, then 1.75 at the merge of dev d7641b99f that brought line RG item D5b-3's 1.74, then 1.81 and
  exception 25 at the merge of dev 34cc301c4 that brought 1.75 to 1.80 and exceptions 22 to 24
  (line NP items M5 and O4c, line RG items X4b-2, D3c-3a and D5b-5, line PS item PS12c), then 1.82
  and exception 26 at the merge of dev 328728a05 that brought line NP item M4c-1's 1.81 and
  exception 25; then 1.83 and exception 27 at the merge of dev 50474994c that brought psiblast line
  item PS10a's 1.82 and exception 26; then 1.86 at the merge of dev 719709924 that
  brought line NP item M4b-AX's 1.83, psiblast line item PS10c's 1.84 and line RG item D1c-1's 1.85;
  then 1.92 and exception 28 at the merge of dev ad35fcaee that brought 1.86 to 1.91 and exception 27
  (psiblast line item PS10d, line NP items O4b2-2, M4b-3c and M4b-T2, line RG item T5d, line PS item
  PS12d); then 1.93 at the merge of dev 2cd5a1906 that brought line NP item O4b2-4's 1.92; then 1.94
  at the merge of dev 73242b104 that brought line RG item D1c-1c's 1.93; then 1.95 at the merge of dev
  bb64581ea that brought line NP item M4c-2's 1.94; then 2.00 at the merge of dev 4626ddd5d that
  brought 1.95 to 1.98 and exception 29: line NP item M4b-3b, line DBID item pfoff-4, line RG items
  X-short-bat and X-bpgap, with 1.99 left to line RG item D5's pull request; then 1.104 at the merge
  of dev 518d403cc that brought 1.99 to 1.103 and exception 30 (line RG items X-short-tbx-link,
  KC-11d, D5c and K10a, psiblast line item PS10e); then 1.106 at the merge of dev f25e2276a that
  brought 1.104 and 1.105 (line NP item M4b-T3, psiblast line item PS10b); then 1.108 and exception 31
  at the merge of dev 5c7e93761 that brought 1.106 and 1.107 and exception 28 (line NP items M4b-3c
  part 2 and M4a-mdepth); then 1.109 at the merge of dev df03b20e2 that brought 1.108 (line RG item
  X-short-bat-2); then 1.113 and exception 32 at the merge of dev 245b7e2e3 that brought 1.109 to
  1.112 and exception 31 (line RG items D1c-2, X-short-bat-2, D3c-4 and psiblast line item PS7c);
  then 1.114 at the merge of dev fe1452e21 that brought 1.113 (psiblast line item PS14); then 1.115
  at the merge of dev 42678f399 that brought 1.114 (line DBID item pfoff-7); then 1.116 and
  exception 33 at the merge of dev 9758dbd59 that brought 1.115 and exception 32 (line RG item G1);
  then 1.117 at the merge of dev f71a5cb85 that brought 1.116 (line RG item X-chunk); then 1.118
  at the merge of dev 351661f13 that brought 1.117 (line WQA item D23));
  1.119 2026-11-24 (approved exception 34, the Owner's answer of 2026-11-24 00:04 JST, on the board, to
  line RG's question Q-F2-28 (recommended A), 「Q-F2-28=A」: a negative `WindowSize` and a
  `LongestIntronLength` of an imported strategy whose arithmetic passes `Int4`, and the one-thread
  result of a negative window; line RG item F2 fix round 3, review 3's R3-B2; numbered after 1.118;
  first 1.87 (then proposed as exception 28 and pending the Owner), then 1.93 and exception 29 at the
  merge of dev ad35fcaee that brought 1.86 to 1.91 and exception 27; then 1.94 at the merge of dev
  2cd5a1906 that brought line NP item O4b2-4's 1.92; then 1.95 at the merge of dev 73242b104 that
  brought line RG item D1c-1c's 1.93; then 1.96 at the merge of dev bb64581ea that brought line NP
  item M4c-2's 1.94; then 2.01 and exception 31 at the merge of dev 4626ddd5d that brought 1.95 to
  1.98 and exception 29 (line NP item M4b-3b), with 1.99 and exception 30 left to line RG item D5's
  pull request; then 1.105 at the merge of dev 518d403cc that brought 1.99 to 1.103 and exception
  30 (line RG items X-short-tbx-link, KC-11d, D5c and K10a, psiblast line item PS10e); then 1.107
  at the merge of dev f25e2276a that brought 1.104 and 1.105 (line NP item M4b-T3, psiblast line
  item PS10b); then 1.109 and exception 32 at the merge of dev 5c7e93761 that brought 1.106 and
  1.107 and exception 28 (line NP items M4b-3c part 2 and M4a-mdepth); then 1.110 at the merge of dev
  df03b20e2 that brought 1.108 (line RG item X-short-bat-2); then 1.114 and exception 33 at the merge
  of dev 245b7e2e3 that brought 1.109 to 1.112 and exception 31 (line RG items D1c-2, X-short-bat-2,
  D3c-4 and psiblast line item PS7c); then 1.115 at the merge of dev fe1452e21 that brought 1.113
  (psiblast line item PS14); then 1.116 at the merge of dev 42678f399 that brought 1.114 (line DBID
  item pfoff-7); then 1.117 and exception 34 at the merge of dev 9758dbd59 that brought 1.115 and
  exception 32 (line RG item G1); then 1.118 at the merge of dev f71a5cb85 that brought 1.116 (line
  RG item X-chunk); then 1.119 at the merge of dev 351661f13 that brought 1.117 (line WQA item D23));
  1.120 2026-11-26 (the rule-2 explicit rejection of a stored definition-line title with a byte
  outside `' '` to `'~'`, of 1.10, narrowed to the cases that read the definition line outside
  a query batch's Seq-aligns: formats other than `-outfmt` 8, 9, 11 and 12, or 8, 9 and 12 with
  `PRE_FETCH_SEQS_LIMIT=2147483647`; `psiblast`, `rpsblast` and `rpstblastn` (except a query
  batch that ends at a header NCBI cannot read); searches restricted
  by a seqid, GI or taxid list, alias files with a list or an OID mask, or a volume with an OID
  mask type; a query batch whose headers also give NCBI's warning of a UTF8String; and a stored
  Seq-id with such a byte. Otherwise BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX write NCBI's
  output and `Error: (803.7)` messages, rule 2; and in any format of all eight search programs
  (including `psiblast`, `rpsblast` and `rpstblastn`), a query batch that ends at a header NCBI
  cannot read writes NCBI's `Error: (803.7)` messages for the titles read before it, then NCBI's
  error; line RG-KB item 9(b), Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.119 at the merge of dev bf0e63c3d; first 1.107,
  then 1.113, then 1.114, then 1.115, then 1.116, then 1.117, then 1.118);
  1.121 2026-11-24 (a database sequence whose Seq-align id is a pdb id whose handle has another
  chain, which NCBI's data loader serves no sequence for: the empty Bioseq handle reproduced in
  every report, rule 2, and the explicit rejections, rule 3, of its tabular `sacc` and `saccver`
  in the class of 1.37, of its tabular sequence-reading rows without a subject-block field (where
  NCBI stops with `CAlnException::eInvalidSeqId`, exit 255) or minus-strand translated-frame rows
  with an alignment-block field and no subject-block field (where NCBI stops with
  `CCoreException::eNullPtr`, exit 255), of its printed `-show_gis` description line,
  whose GI NCBI never sets, and of NCBI's throw where another id of the sequence holds the
  handle's id but no definition line does as NCBI's load reads the ids (where another definition
  line holds it, the reports show that line's title as NCBI's, review 2); line RG-KB item 12, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; first 1.92, then 1.93 at the merge of dev that
  brought line NP item O4b2-4's 1.92, then 1.94 at the merge of dev that brought line RG item
  D1c-1c's 1.93, then 1.95 for line NP item M4c-2's 1.94, which dev took after that merge, then
  1.96 at the merge of dev 2a5889938 that brought line NP item M4b-3b's 1.95, then 1.97 at the
  merge of dev 3d8c77b16 that brought line DBID item pfoff-4's 1.96, then 1.99 at the
  merge of dev 4626ddd5d that brought line RG item X-short-bat's 1.97 and line RG item X-bpgap's
  1.98, then 1.100 at the merge of dev 0dbeba38b that brought line RG item X-short-tbx-link's
  1.99, then 1.103 at the merge of dev fbf312db0 that brought line RG items KC-11d's 1.100 and
  D5c's 1.101 and psiblast line item PS10e's 1.102, then 1.104 for line RG item K10a's 1.103, which
  dev took after that merge, then 1.106 at the merge of dev 17c46261c that brought line RG item
  K10a's 1.103, line NP item M4b-T3's 1.104 and psiblast line item PS10b's 1.105, then 1.107 at
  the merge of dev 22ac798e5 that brought line NP item M4b-3c part 2's 1.106, then 1.108 at the
  merge of dev f2be0e496 that brought line NP item M (b)'s 1.107, then 1.109 at the merge of dev
  10f846fa3 that brought line RG item X-short-bat-2's 1.108, then 1.111 at the merge of dev
  924bf3be4 that brought line RG item D1c-2's 1.109 and line RG item X-short-bat-2's 1.110, then
  1.113 at the merge of dev 245b7e2e3 that brought psiblast line item PS7c's 1.111 and line RG
  item D3c-4's 1.112, then 1.114 at the merge of dev f17536f02 that brought psiblast line item
  PS14's 1.113, then 1.115 at the merge of dev 42678f399 that brought line DBID item pfoff-7's
  1.114, then 1.116 at the merge of dev c2544a7dd that brought line RG item G1's 1.115, then
  1.117 at the merge of dev 00d811a41 that brought line RG item X-chunk's 1.116, then 1.118 at
  the merge of dev 351661f13 that brought line WQA item D23's 1.117, then 1.120 at the merge of
  dev bf0e63c3d that brought line RG item F2's 1.118 and 1.119, then 1.121 at the merge of dev
  9f87bd3a2 that brought line RG-KB item 9(b)'s 1.120);
  1.122 2026-11-27 (BLASTP and BLASTX `-db_hard_mask`, line NP item M4b-3d: `PartsOutOfOrder`, hard
  ranges longer than 5,000,000 residues (`LongRange`), and in-bounds `short_chunk_merge_refusal`
  ported under rule 4 via the shared `core::blast_hits::merge_subject_chunk_hsps` helper; explicit
  rejection, rule 3, of `HardMaskShape::LastRangeOverlap` (before the search) and
  `last_range_overlap_merge_refusal` (during the search for short `OneWord` chunks) where the final
  chunk's out-of-bounds overlap read (`blast_engine.c:579-580`, `hard_ranges[num_hard_ranges].left`)
  would decide a merge with an alignment from an earlier hard range (same B-d boundary as BLASTN's
  `ShortLastMask`: merges within a final range longer than 5,000,000 letters keep
  `DBSEQ_CHUNK_OVERLAP = 100`); version 1.21's explicit rejection, rule 3, of a chunk-overlap merge
  whose query lengths sum to 0 or whose density-scaled score leaves `Int4`
  (`chunk_merge_undefined_refusal`, typed `ShortChunkSeedRefusal` naming the subject OID and, for
  TBLASTN, reading frame) extended to BLASTP and BLASTX `-db_hard_mask`, TBLASTN subject-chunk
  merges (both frames longer than 5,000,000 residues and `-db_hard_mask`), and BLASTN subject-chunk
  and query-chunk merges (scaled score leaving `Int4`); Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07; numbered after 1.121 at the merge of dev cb7044346 that brought
  line RG-KB item 12's 1.121; first 1.113, then 1.114, then 1.115, then 1.116, then 1.117,
  then 1.118, then 1.120, then 1.121);
  1.123 2026-11-24 (the tabular rows of 1.37, 1.50, 1.114 and 1.121, and `sallacc` of a database
  sequence with a definition line without a Seq-id, at which NCBI stops with its null pointer or
  alignment exception or prints an uninitialised taxid: NCBI's bytes before the stop are written,
  then the refusal of rule 3, in every program, through one owner of the rows' stops; line RG-KB
  item 18b, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after
  1.122; first 1.100, then 1.101 at the merge of dev that brought line RG item X-short-tbx-link's
  1.99 and line RG item KC-11d's 1.100, then 1.104 at the merge of dev of 2026-11-25 that brought
  line RG item D5c's 1.101, psiblast line item PS10e's 1.102 and line RG item K10a's 1.103, then
  1.108 at the merge of dev of 2026-11-26 that brought line NP item M4b-T3's 1.104, psiblast line
  item PS10b's 1.105, line NP item M4b-3c part 2's 1.106 and line NP item M (b)'s 1.107, then
  1.109 at the merge of dev of 2026-11-26 that brought line RG item X-short-bat-2's 1.108, then
  1.113 at the merge of dev of 2026-11-26 that brought line RG item D1c-2's 1.109, line RG item
  X-short-bat-2's 1.110, psiblast line item PS7c's 1.111 and line RG item D3c-4's 1.112, then
  1.114 at the merge of dev of 2026-11-26 that brought psiblast line item PS14's 1.113, then 1.115
  at the merge of dev of 2026-11-26 that brought line DBID item pfoff-7's 1.114, then 1.116 at the
  merge of dev of 2026-11-26 that brought line RG item G1's 1.115, then 1.117 at the merge of dev
  of 2026-11-26 that brought line RG item X-chunk's 1.116, then 1.118 at the merge of dev of
  2026-11-26 that brought line WQA item D23's 1.117, then 1.121 at the merge of dev of 2026-11-27
  that brought line RG item F2's 1.118 and 1.119 and line RG-KB item 9(b)'s 1.120, then 1.122 at
  the merge of dev of 2026-11-27 that brought line RG-KB item 12's 1.121, then 1.123 at the merge
  of dev of 2026-11-27 that brought line NP item M4b-3d's 1.122);
  1.124 2026-11-26 (approved exception 35, `losat blastdbcheck` alias count checks and
  `-ends` OID arithmetic widened beyond NCBI's signed 32-bit/64-bit integer overflow,
  line RG-D2 items D2b-2..5; Owner 2026-11-25 17:54
  JST, Q-D2b-OV=A; numbered after 1.123 at the merge of dev c66bea41b that brought 1.123
  (line RG-KB item 18b); first 1.109, then 1.110, then 1.113, then 1.114, then 1.115, then 1.116,
  then 1.117, then 1.118, then 1.120, then 1.121, then 1.122, then 1.123; exception 31 until 1.111 took 31,
  exception 32 until 1.115 took 32, and exception 33 until 1.118/1.119 took 33 and 34);
  1.125 2026-11-27 (approved exception 36, a protein database with a sequence of length 0
  searched at more than one thread, line RG-D2, found in review R4 of PR #1660; Owner 2026-11-27
  09:25 JST, Q-DBA0-MT=KEEP; first 1.124);
  1.126 2026-11-24 (BLASTP, psiblast, RPSBLAST and TBLASTN, line RG item T4b-f: the prolog of
  `-outfmt` 0 to 4 is written before the first query batch, where NCBI writes it (before its first
  search strategy request), as BLASTN's, BLASTX's, TBLASTX's and RPSTBLASTN's, so the explicit
  rejections that come at a query batch, in the search or at the report follow the prolog instead
  of coming before the first byte, as NCBI's stops there follow it (rows whose stdout is NCBI's
  prolog of a search of the same database, `db.pro.*` and `db.id.pro.*`); no new exception, rule
  3; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07; numbered after 1.125;
  first 1.99, then 1.100 at the merge of dev of 2026-11-24 that brought 1.99, line RG item
  X-short-tbx-link, then 1.102 at the merge of dev of 2026-11-24 that brought 1.100, line RG item
  KC-11d, and 1.101, line RG item D5c, then 1.106 at the merge of dev of 2026-11-25 that brought
  1.102 to 1.105, psiblast line items PS10e and PS10b, line RG item K10a and line NP item M4b-T3,
  then 1.107 at the merge of dev of 2026-11-25 that brought 1.106, line NP item M4b-3c part 2, then
  1.108 at the merge of dev of 2026-11-26 that brought 1.107, line NP item M (b), then 1.109 at the
  merge of dev of 2026-11-26 that brought 1.108, line RG item X-short-bat-2, then 1.114 at the
  merge of dev of 2026-11-26 that brought 1.109 to 1.113, line RG items D1c-2, X-short-bat-2 and
  D3c-4 and psiblast line items PS7c and PS14, then 1.115 at the merge of dev of 2026-11-26 that
  brought 1.114, line DBID item pfoff-7, then 1.118 at the merge of dev of 2026-11-26 that brought
  1.115 to 1.117, line RG items G1 and X-chunk and line WQA item D23, then 1.120 at the merge of dev
  of 2026-11-26 that brought 1.118 and 1.119, line RG item F2, then 1.121 at the merge of dev of
  2026-11-27 that brought 1.120, line RG-KB item 9(b), then 1.122 at the merge of dev of 2026-11-27
  that brought 1.121, line RG-KB item 12, then 1.123 at the merge of dev of 2026-11-27 that brought
  1.122, line NP item M4b-3d, then
  1.124 at the merge of dev of 2026-11-27 that brought 1.123, line RG-KB item 18b, then 1.126 at the merge
  of dev of 2026-11-27 that brought 1.124 and 1.125, line RG-D2 items D2b-2..5 and exception 36);
  1.127 2026-11-28 (`-evalue` NaN and the infinities, line RG item X-evalue-nan: the item of rule
  4 now holds in every program, BLASTX and RPSTBLASTN search NaN, Spouge's cutoff search, the
  reaps and the redo's test of the best e-value take NCBI's comparisons, and the explicit
  rejection of an infinite `-evalue` (decision D12) is limited to the composition-based redo of a
  database search, as DBL_MAX's is (RPSBLAST and RPSTBLASTN run both, and their redo hands lists
  back raw as NCBI's does); no new exception, rules 3 and 4; Owner-delegated, standing instruction 2026-09-29,
  restated 2026-10-07)
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
   `OVERLAP_CHUNK_SIZE` such that every chunk holds a query and a query chunk is long enough to
   be split again (an overlap close to the chunk size, such as `CHUNK_SIZE` 101 to 150 with the
   default overlap 100 and a query length of at least 101), NCBI skips the chunk's lookup table
   and stops with a `CCoreException` (null pointer) that names its build's source files, exit 3
   (`split_query_aux_priv.cpp:190-201`, `blast_aux_priv.cpp:206-208`,
   `prelim_search_runner.hpp:77`; only a debug build asserts against it; at `-num_threads` > 1
   with `-db` worker threads log the exception and exit 0 without that chunk's hits). Where any
   chunk holds no query (for example `CHUNK_SIZE` 101 to 150 with overlap 100 and length at most
   100), NCBI stops earlier in `CQuerySplitter::Split`, which is a rule-3 explicit rejection
   below. LOSAT searches each non-empty re-splitting chunk once. Evidence
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
   sequence; and (BLASTP and BLASTX with `-db_hard_mask`) the traceback's read past a hard-masked
   copy of a subject.**
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

   BLASTP and BLASTX with `-db_hard_mask` (version 1.83, decided by the Owner on 2026-11-23 13:49 JST,
   line NP item M4b-AX; question Q-ALIGNEX, recommended option A, the Owner's words 「推奨通り」):
   NCBI's traceback fetches a masked copy of each subject that has masks
   (`seqsrc_seqdb.cpp:326-330, 346-353`), allocated at the sequence's length with no sentinel byte
   for a protein (`seqdbvol.cpp:1633-1635`, `:1488`), where the mapped volume holds a NULLB before
   and after every sequence. `ALIGN_EX` reads one letter past each end of its subject when a row
   of its dynamic program reaches that end (`blast_gapalign.c:481-492, 555-578, 671-675`): the byte
   before the copy (the allocator's header, 0 with glibc) and the byte after it (stale heap data).
   That byte feeds only the fence test (`:571-576`) and a score that no cell uses. Where it equals
   the fence byte (`FENCE_SENTRY`, 201; `blast_util.h:364`), the traceback is rolled back and
   retried once (`blast_traceback.c:516-520, 701-709, 1651-1688`); where the retry reads that byte
   again, NCBI formats preliminary HSPs that have no edit script and crashes (SIGSEGV,
   `blast_seqalign.cpp:676-677`), sometimes after part of its rows. Whether it does depends on the
   heap's history: the batch size, the layout of the query file, the process environment. With no
   mask, with `-db_soft_mask`, or for a subject without a mask, NCBI reads the mapped volume's NULLBs
   and does not meet this. LOSAT reads the bytes past both ends of a hard-masked subject as those
   NULLBs, as this exception reads the letters past a sequence, and writes the rows of NCBI's runs
   that do not meet the fence byte; its code and its output do not change with version 1.83.
   TBLASTN and TBLASTX are not affected: their subject is a translation buffer with NULLBs
   (`blast_util.c:438-453`). Evidence: with `MALLOC_PERTURB_=201` (glibc fills freed memory with that
   byte) NCBI crashes (exit 139) on the overlapping-mask database (BLASTP self-search and BLASTX) and
   on the database whose masks do not overlap (BLASTX); without it NCBI ends with exit 0 and
   the bytes that LOSAT writes: BLASTX against the overlapping-mask database with `-window_size 0`,
   23,074 rows at every batch size; BLASTP of that database against itself, 3,695 rows; BLASTP of
   the other, 311 rows (`cmp`-equal). Under `-db_soft_mask` and without a mask NCBI does not crash
   with `MALLOC_PERTURB_=201`. The kernel log shows 15 crashes of this form, at `s_BlastHSP2SeqAlign`
   reading the missing edit script. So NCBI's output is either LOSAT's bytes or a crash, never other
   rows: the value read feeds only the fence test and a score that is never used.

Explicit rejections (rule 3):

- `-evalue` infinite or of DBL_MAX or more with the composition-based redo (`-comp_based_stats`
  1 to 3; BLASTP, PSIBLAST and TBLASTN, and an infinite value in BLASTX; decision D12): a hit list whose HSPs all drop keeps the best e-value DBL_MAX, which passes
  `best_evalue <= expect_value` (`blast_kappa.c:409,3687`), and the empty list is read
  (`blast_hits.c:3266`): NCBI crashes on some inputs, and LOSAT cannot tell in advance
  which. `1.7976931348623156e308` and below run and match NCBI. Without the redo NCBI never
  puts an empty list into the hit list, and LOSAT runs these values: DBL_MAX since line O4c's
  review D1 (2026-11-21), an infinite value since line RG item X-evalue-nan (version 1.127; it was
  refused in every mode, BLASTX's and RPSTBLASTN's too, though NCBI's cutoff of it is the cutoff of
  any e-value that no score reaches: `BLAST_Cutoffs` keeps 1, Spouge's search gives 0). TBLASTX
  and BLASTN run these values and match NCBI, and so do RPSBLAST and RPSTBLASTN with the redo since
  version 1.127: RPS-BLAST's traceback gives the redo no sequence source, gets every list back and
  drops a list without HSPs before the hit list (`blast_traceback.c:1238-1275`), so NCBI does not
  crash there (they were refused before).
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
  (version 1.1, rule 2). BLASTP's rejection comes after the reports of the query batches before
  (version 1.53, review 3 of X4: the form of 1.40 for a rejection in the middle of a run; row
  `db.x4.chunk.neg.gap.6`, whose bytes are NCBI's report of the batch before, row
  `db.x4.chunk.gap_anchor.6`; rows `db.x4.chunk.empty.*` and `db.x4.chunk.again.*` hold NCBI's
  bytes of the other two stops), at TBLASTN's step since version 1.98 (below).
  TBLASTN refuses where NCBI stops in the batch's search (a chunk without a query, also where the
  first range of a negative `CHUNK_SIZE` covers the batch, at which `CQuerySplitter::Split` stops
  before any chunk is searched, and a chunk split again, at which the search of that chunk stops
  (at more than one thread NCBI logs that null pointer, reports the batch without that chunk's
  hits and exits 0; LOSAT refuses at every thread count, as NCBI does at `-num_threads 1`); also on
  a database left without a subject: a seqid list that keeps none, an alias over volumes without
  sequences) after the reports of the batches before, in NCBI's order with the setup of each
  chunk's search and its read of the database (version 1.97), and so does the rejection of a
  negative `CHUNK_SIZE` whose wrapped ranges leave gaps between them, which NCBI searches on: where
  the search of the batch's first chunk with a valid context would read its subjects, after that
  chunk's setup checks and the reads of the database that fail before it, with subjects in hand
  (a batch searched without a subject is reported as NCBI reports it: nothing is lost) (version
  1.97; before, before the run's first byte; row `db.xs.tblastn.bat.gap.6`, whose bytes are
  NCBI's report of the batch before, row `db.xs.tblastn.bat.qa.6`). Since version 1.98 BLASTP
  and psiblast take the same steps, through the one owner of the steps of a split batch that
  TBLASTN's search uses (line RG item X-bpgap, NCBI 2.17.0 probes): the rejection of the gaps comes
  where the search of the batch's first chunk with a valid context would read subjects in hand,
  and a batch searched without a subject in hand (a seqid list that keeps none) or whose chunks
  have no valid context is reported as NCBI reports it (rows `db.x4.chunk.neg.gap.nosub.*`,
  `db.x4.chunk.neg.gapvalid.*`); the stops of the chunks come before the errors of the database
  and of the subject length that NCBI meets later, in the setup of a chunk or its search (NCBI
  stops at the null pointer first); a batch of invalid queries only, which NCBI does not search,
  is not split and meets no stop; and a split batch whose chunks have no valid context ends with
  NCBI's `Traceback failed` (exit status 3) where the traceback's own setup of the batch fails the
  subject-length check, and reads nothing where the database's sequence file cannot be read (exit
  status 0), in BLASTP, psiblast and TBLASTN. A refused batch of two queries with one id that NCBI splits
  gives NCBI's error where its chunks reach the setup's failure, and the refusal where a chunk
  stops first.
  BLASTX (version 1.110, line RG item X-short-bat-2, rows `db.xs.blastx.bat.*`) refuses a chunk
  split again and a chunk without a query at NCBI's split, after the reports of the query batches
  before and before the database's errors of a chunk's setup; every negative pair that splits a
  BLASTX batch leaves such chunks (BLASTX evens the chunk size to a multiple of 3 and does not
  adjust it again, so chunk 0 covers the batch). Its rejection of an overlap that is not a
  multiple of 3, whose chunks NCBI searches out of the query's reading frames, comes after the same
  reports, only where the batch has subjects in hand whose sequences can be read.
- Settings on which NCBI crashes (decision D3): an IDENTITY word size of 5 or 6,
  `-use_sw_tback -ungapped -comp_based_stats 0` (since version 1.87 at the first query batch NCBI searches), compressed-lookup thresholds that
  exhaust its overflow bank (in BLASTP since version 1.108, as in TBLASTN since 1.97, in the query
  batch whose table exhausts it, after the reports of the batches before; in BLASTX since version
  1.110, after the reports of the query batches before, exit status 1). (`-max_target_seqs` values that make NCBI's preliminary hit
  list size not positive were rejected before version 1.66, approved exception 19.) A query length plus `-window_size` in (2^30, 2^31 - 1], where
  NCBI does not finish (decision D8).

Kept as explicit rejections although NCBI's result is deterministic (rule 2's
maintainer clause; decision D11): a query length plus `-window_size` above 2^31 - 1, where
NCBI's `Int4` diagonal table wraps to one cell and no hit is found (BLASTP, TBLASTX;
together with D8, every sum above 2^30 is rejected); and (version 1.117, line WQA item D23,
Owner decision `Q-D23=A` of 2026-11-24) a `-logfile` pointing to any input file read by the run
(all programs, command line and library API, matching LOSAT Web's D19 rule): NCBI BLAST+ 2.17.0
opens `-logfile` in append mode before loading configuration or reading inputs, which causes
`tblastx -query q.fa -logfile q.fa` to read its own warnings back in an endless loop (growing
`q.fa` to 10 GiB in 5 minutes in probe `s2-query`), corrupts `.ncbirc` on its second read
(`s2-ncbirc-self`), and appends diagnostics into `-subject` and other inputs; LOSAT refuses any
such `-logfile` before opening or creating the file (exit status 1, existing inputs unchanged).

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
   equal `ln K + ln searchsp - ln E` (`stats::cutoffs`), and on a nearby pair
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

## Approved exception 28 (makeblastdb `-mask_data`, Seq-locs nested deeper than NCBI's stack holds, version 1.107)

Version 1.107, line NP item M (b). The Owner answered on 2026-11-24 (09:13 JST, through the session
board, record ce1e4751-4972-4b38-bd21-36be4d8fbd98) 「Q-MDEPTH=B」. The question Q-MDEPTH was
「入れ子の深いマスクファイルを makeblastdb -mask_data でどこまで読むか」 (how deep should makeblastdb
`-mask_data` read a mask file whose Seq-locs nest?), and its option B read 「B: 深さに関係なく読む（推奨）」
(B: read at any depth (recommended)), the form of approved exception 15. Not taken: option A, refusing
the depths at which NCBI may crash: NCBI's crash depth moves by some dozen levels from run to run and differs
between text and binary files, so any limit refuses files that NCBI reads at every start or reads
files on which NCBI crashes at every start. Until this version LOSAT refused nests deeper than 400
object frames, 97 levels of `mix` and more, which NCBI reads at every start.

28. **A mask file whose Seq-locs nest deeper than NCBI's stack holds.** makeblastdb's reader of the
    `-mask_data` files (`CObjectIStreamAsn` for text, `CObjectIStreamAsnBinary` for binary) reads
    a Seq-loc, the one type of a mask file that holds values of its own type (`mix` and `equiv`,
    and a value written as `:Seq-loc` in NCBI's pointer syntax), by recursion: each level opens the
    choice, its variant, the list and the list's element (`ReadChoiceSimple`,
    `src/serial/objistr.cpp:1482-1500`; `ReadContainer`, `src/serial/objistrasn.cpp:1245-1274` in
    text and `src/serial/objistrasnb.cpp:955-1000` in binary), and an element written as `:Seq-loc`
    the named frame of `ReadPointer` (`src/serial/objistr.cpp:1173-1198`), about 1.3 KB of the main
    thread's stack for each level. With the default stack limit of 8 MiB
    NCBI BLAST+ 2.17.0 overflows its stack and ends with SIGSEGV (exit status 139; the log's first
    lines, nothing on standard error, no database): in a minimal environment, 5 runs at each
    depth, the bioconda build read every nest of up to 6,449 levels and crashed in every run from
    6,474 levels in binary files, and read up to 6,610 levels and crashed from 6,637 in text files,
    crashing in some runs only at the depths between; NCBI's Linux x86-64 release crashes deeper
    (binary files read up to 6,784 levels, crashes from 6,809; text files 6,872 and 6,899); `mix`,
    `equiv`, the two alternating and the pointer syntax alike. As in
    exception 15 the depth is not the input's: the stack limit, the lengths of the command line and
    of the environment and the stack's random start at each run move it. With no stack limit
    (`ulimit -s unlimited`) NCBI reads every depth tried, to 100,000 levels (170 to 239 MB of
    memory, by build and kind of nest), and writes its normal database, or for a file damaged deep
    inside its report, which names every frame of the nest (at 100,000 levels about 25 s, for a
    report of 2 MB; a number too large for `int` in a text
    file, whose report gives each frame a line of its own, took 245 s at 5,000 levels and 4,933 s
    and 253 MB at 20,000, for a report of 17 MB). LOSAT keeps the levels and the frames on the heap
    (`TypedReader::nest` of `report::serial::asn_types::typed`, which `writedb::mask_info::reader`
    reads a Seq-loc with), so it reads every depth as far as memory allows and writes NCBI's output with a stack large enough, its reports in time linear
    in the depth. No mask file that dustmasker, segmasker, windowmasker or convert2blastmask writes
    nests this way (their locations are `int` and `packed-int`).
    Evidence: `docs/evidence/mask_nest_depth/` (for each kind of nest and depth, NCBI's runs with
    8 MiB and with no stack limit, LOSAT's runs and the generator of the inputs). Fixture rows
    `n_*` of `LOSAT/tests/seqdb/mask_parse/expected.tsv`, frozen from NCBI's run with no stack
    limit (`LOSAT/tests/seqdb/make_mask_parse_fixtures.py`): every kind of nest, damaged innermost
    in `mix`, `equiv`, alternating and pointer nests, and at a middle level, at each step that
    begins or ends a level. Unit tests `a_nest_of_a_million_levels_is_read_without_recursion`,
    `every_nesting_path_is_read_at_depth`,
    `a_null_element_deep_inside_is_found_where_ncbi_dereferences_it`,
    `an_error_deep_inside_leaves_every_frame_of_the_nest`,
    `an_error_deep_inside_other_nests_leaves_every_frame`,
    `an_error_at_a_middle_level_leaves_the_frames_outside_it` (each fixture's damage at 1,000,001
    levels gives its report with every level added) and
    `a_number_too_large_takes_a_line_for_each_frame`
    (`LOSAT/src/writedb/mask_info/depth_tests.rs`).

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
    leaves a frame an empty first part (an empty first part of a frame, then a part that NCBI
    searches as the whole frame; in TBLASTN and TBLASTX this includes an empty mask at 0 before a
    tail mask, which BLASTP and BLASTX search under Q1 as BLASTN does under exception 13), refused
    (rule 3) until version 1.88, is searched under approved exception 27 (a translated mask that
    leaves a frame an empty first part). A mask that ends past its sequence and in BLASTP and BLASTX
    a range longer than 5,000,000 letters stay refused (rule 3). Masks out of order and overlapping
    masks stayed refused in this version; from version 1.95 BLASTP and BLASTX search them (rule 4,
    and approved exception 29) except the rejections of 1.95 (the end of the Explicit rejections
    list), and TBLASTN and TBLASTX still refuse them. `blastinput::db_mask::HardScan`,
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
    bytes, no row repeated; LOSAT's the same since TBLASTN `-ungapped` runs, line NP item O4c,
    row `db.tmask.tn.hard.th2.13.ung.6`);
    BLASTP with `-window_size 0` on `pe_tail` (`db.pmask.bp.hard.tail.pe_tail.w0.6`): the near
    input's. Unit tests `protein_shapes_of_the_hard_masks` and
    `translated_hard_masks_in_order_are_searched_also_where_a_frame_restarts_at_its_start`
    (`LOSAT/src/blastinput/db_mask_tests.rs`),
    `a_whole_frame_scans_to_the_frame_s_end_and_a_tail_masked_one_its_chunk`
    (`LOSAT/src/algorithm/common/translated_chunks_tests.rs`).

## Approved exception 27 (TBLASTN and TBLASTX `-db_hard_mask`, a mask that leaves a reading frame an empty first part, version 1.88)

The Owner decided on 2026-11-22 at 22:34 JST (words 「Q1=B Q2=B Q3=A Q4=A Q5=C Q6=B Q7=B Q8=B Q9=A
Q10=A Q11=A Q12=A Q13=OK」; question Q6 of that batch is question Q4 of line NP's M4b,
`M/m4b/OWNER_QUESTIONS.md`) option B, which the batch worded 「B（推奨）. その frame は mask されて
いないものとして普通に探す（その短い mask を除いた時の NCBI の結果と同じ規則）」 (the frame is
searched as a frame without a mask, by the rule of NCBI's result with that short mask removed): in
TBLASTN and TBLASTX with `-db_hard_mask`, a mask within the first 2 bases of a sequence, or within the
1 or 2 bases after its last whole codon, on which NCBI BLAST+ 2.17.0 crashes in every run, is no
longer refused. Not taken: option A, refusing such databases (what LOSAT did from M4b PR 2 until
this version), which leaves every database with a short mask at a sequence end unsearchable. With
Q1 = B (approved exception 21) the same NCBI mechanism, an empty mask at 0 before a tail mask, is
already searched in BLASTP and BLASTX, so the four programs now agree.

The words of option B can be read two ways for the 1 or 2 masked bases themselves, which are also
letters of the frames of the other strand (e3's [0, 1) is the first residue of frame +1 and the last
of frame -1). Read literally, 「その短い mask を除いた」 is a database with the mask removed and the
real bases kept (the near input of a reading A); read as the Owner's earlier answers Q1 = B and
Q2 = B (exception 21) read masked letters, the bases stay masked in every frame and only the
crashing frame's search is the unmasked one (reading B; a reading C would give the real bases to the
crashing frame alone). The Owner answered that follow-up question (§ Q6 of
`M/m4b/OWNER_QUESTIONS.md`, whose option B was "the crashing frame is searched as without the mask
(the whole frame), and the 1 or 2 masked bases are read as masked in every frame") on 2026-11-23 at
13:49 JST on the board (`tools/session-board/sent.jsonl`, uuid
ed3e0671-abb5-4c21-b3bd-eb22f24af44f) with 「推奨通り」, to the report of line NP's session ad-M4b-8,
which asked Q-ALIGNEX and M4b-3c's question Q6 (recommended A and B and gave the answer forms
「Q-ALIGNEX=A Q6=B」 and 「すべて推奨通り」): Q6 = B, which is what this version implements. Its version
is in the list at the top; elsewhere it is cited by its number and name, "approved exception 27 (a
translated mask that leaves a frame an empty first part)". Line NP item M4b-3c.

27. **`-db_hard_mask` in TBLASTN and TBLASTX: a frame whose masks leave it an empty first part.**
    NCBI maps each stretch between a sequence's masks to every reading frame by its positions divided
    by 3, the minus frames from frame -1's length (`blast_engine.c:805-830`), and cuts each frame as a
    protein subject (`s_GetNextSubjectChunk`, `blast_engine.c:221-311`). A mask within the first 2
    bases of a sequence maps to an empty first part [0, 0) in all three plus frames, a mask
    within the 1 or 2 bases after the last whole codon to an empty first part in all three minus frames
    (the mapping has no frame offset: [1, 2) of 1,200 bases gives an empty first part, (0, 0), in frames
    +1, +2 and +3), and an empty mask at 0 to one in the plus frames. The empty part takes the branch that
    allocates the subject's ranges (`s_AllocateSeqRange`, `blast_engine.c:184-199`); the next part,
    from 0, then takes the "no chunking" branch (`blast_engine.c:264-269`), which points the ranges at
    the backed-up full range, and NCBI searches the frame whole (or, where a mask runs to the end of
    the frame after it, as in exception 21). Its restore then frees that pointer, which it did not
    allocate (`s_RestoreSubject`, `blast_engine.c:200-215`): SIGABRT `free(): invalid pointer` (exit
    134) or SIGSEGV (139), in every run, at 1 and 4 threads. The empty part adds nothing: its scan
    reads at most one word at 0, and the diagonal offset that the word leaves (`Blast_ExtendWordExit`,
    `blast_extend.c:162-173`) is a window or more behind every hit of the next part on that diagonal
    (`aa_ungapped.c:495-500, 525-545`). LOSAT searches the frame as its parts come, the whole frame
    after the empty part, as every other hard chunk, and the masked bases stay masked: the traceback,
    the composition-based statistics and the re-scoring of TBLASTX read them as any base, in every
    frame, as for every other mask (a fourfold-degenerate codon keeps its letter). So LOSAT gives what
    NCBI computes before its free(). That is what NCBI's own runs give where another mask of the
    sequence avoids the crash: in 12 of 12 probes (6 forms, TBLASTN and TBLASTX), the rows that the
    short mask moves in NCBI's own run with a far mask added are the rows that storing its bases as
    `N` moves in the near input, and nothing outside the form's sequence moves. The near input of
    this reading (B) is the database with that short mask's bases stored as `N` and the mask dropped;
    the near input of reading A is the literal database with the mask removed and the real bases
    kept, which gives other rows: each alignment through such a residue is one residue longer
    (TBLASTN on e3 with [0, 1): `2 60 4 180` with the mask, `1 60 1 180` with the real base; TBLASTX
    on e2 with [1200, 1202): frame -1's alignment starts at base 1199 with the mask, at 1202 with
    the real bases). In the 32 fixture rows LOSAT gives the `N`-stored near input's bytes, and CI
    passes at 1 and 4 threads, but the two are equal only as long as NCBI's preliminary stage does not
    read a masked residue: that stage translates the stored 2-bit bases (`blast_engine.c:760, 772`),
    makeblastdb stores an `N` as a random base (`writedb_convert.cpp:361-392`; in the fixtures the
    stand-in differs from the real base at every masked position), so in the near input the seed scan
    and the early extensions read that stand-in, where LOSAT, like NCBI's own computation and its
    far-mask runs, reads the real base (`X` appears only in the traceback). Any preliminary step that
    reads such a residue (the seed scan, the ungapped extension, the preliminary gapped score)
    decides which HSPs reach the traceback and from where, so LOSAT could differ from the near input
    in a case that no probe showed: in the crashing frame itself (its masked residue lies inside the
    chunk from 0) as well as in a frame that the empty part does not touch. A second item is the
    end-sentinel word of a frame one residue shorter than its strand's first (frames -2 and -3 of
    1,200 bases), which LOSAT scans as NCBI's own chunks do. Masks out of order,
    overlapping masks and a mask that ends past its sequence stay refused (rule 3).
    `blastinput::db_mask::HardScan`, `HardMaskShape::of_frames`,
    `algorithm::common::translated_chunks::translated_db_mask_chunks`.
    Evidence (NCBI BLAST+ 2.17.0, each run capped; line NP's task folder `M/m4b3c`, `Q4B_PROGRESS.md`,
    `Q4B_BLOCKED`, `probe/`):
    - The 7 forms (e3 [0, 1), [1, 2), [0, 2); e2 [1200, 1202), [1201, 1202); e4 [0, 3) with [1200,
      1201); e3's empty mask at 0 with [1100, 1200)) crash in both programs at 1 and 4 threads (28 of
      28 runs: exit 134 or 139); the near inputs, the proxies with a far mask and the controls end
      normally, and each of the 50 probe inputs gives the same bytes at 1 and 4 threads.
    - Fixture rows `db.tmask.{tn,tx}.hard.te.{13,14,15,31}.6`, `db.tmask.tn.hard.te.{13,14,15,31}.q4.6`,
      `db.tmask.tx.hard.te.31.q4.6`, `db.tmask.{tn,tx}.hard.tq4.11.6` and `db.tmask.tn.hard.tq4.11.q4.6`
      of `LOSAT/tests/db_search_manifest/{tblastn,tblastx}.tmask.tsv` (LOSAT on the database with the
      mask; the refusal rows `db.tmask.*.hard.rej.te.{13,14,15,31}.6` before), with `ncbi-near:` the
      bytes of NCBI's run of the near input (the rows `db.tmask.*.hard.near.te{13,14,15,31}n.*` and
      `db.tmask.*.hard.near.tq4n.*`; databases `te13n`, `te14n`, `te15n`, `te31n` and `tq4n`, and
      `tq4` with the empty mask, built by `LOSAT/tests/db_search/make_dbmaskt_q4.py` from copies of
      `edge.fna` with the bases stored as `N`; `tq4n`'s tail is also shortened to [1100, 1196) as
      exception 21's near input; the queries `qneq4.faa` and `qxeq4.fna` align through the masked
      residues). Unit tests `translated_hard_masks_in_order_are_searched_also_where_a_frame_restarts_at_its_start`
      (`LOSAT/src/blastinput/db_mask_tests.rs`) and `a_frame_after_an_empty_first_part_scans_the_whole_frame`
      (`LOSAT/src/algorithm/common/translated_chunks_tests.rs`).

## Approved exception 29 (BLASTP, BLASTX, TBLASTN and TBLASTX `-db_hard_mask`, diagonal entries that reach another sequence, versions 1.95 and 1.106)

The Owner decided on 2026-11-22 at 22:34 JST (words 「Q1=B Q2=B Q3=A Q4=A Q5=C Q6=B Q7=B Q8=B Q9=A
Q10=A Q11=A Q12=A Q13=OK」 to the question batch of that day, whose Q7 is line NP's question Q5 of
M4b, `M/m4b/OWNER_QUESTIONS.md`) option B: 「配列ごとに独立に探す（前の配列の印が届かない所まで
目盛りを進めます。今の判定の式を置き換えるだけで、速さは変わりません）」, with the merit
「こうした DB も探せて、スレッド数によらず同じ結果になります」 and the demerit 「前の配列の印が次の配列に
届く場合だけ、NCBI の 1 スレッドと bytes が違います…新しい例外になります」. On 2026-11-25 at
17:54 JST (line NP's series M decision 77, question Q-3C2-STALE = B, `qbatch-20261125/ANSWERS.md`,
chosen option 「B: 例外 29 を tblastn/tblastx にも広げ、配列ごとに独立探索する（4 検索統一）」), the
Owner extended this exception to TBLASTN and TBLASTX (version 1.106, line NP item M4b-3c part 2).
Not taken: option A,
refusing a database with such a sequence before the search (the explicit rejection `StaleDiagonals`
of the first draft of version 1.95), which refuses whole seg-masked proteomes (and in TBLASTN and
TBLASTX would refuse databases that BLASTP and BLASTX search); option C, writing
NCBI's 1-thread rows at every thread count, which carries the entries of one sequence into the
next on purpose and makes the search depend on the order of the sequences. Its version is in the
list at the top; elsewhere it is cited by its number and name, "approved exception 29 (diagonal
entries that reach another sequence)". Line NP items M4b-3b and M4b-3c part 2.

26. **`-db_hard_mask` in BLASTP, BLASTX, TBLASTN and TBLASTX: overlapping masks after which NCBI's word hits reach
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
    unchanged. BLASTX searches a query batch in units that end at the first query that brings them
    to 10,002 letters (BLASTX's default batch), which give the batch's rows where the walk keeps no
    sequence apart; where it keeps one apart, BLASTX searches the whole batch as one unit, as NCBI
    does (`blastx::query_setup::unit_ranges`, `HardStarts::keeps_apart`; item M4b-3b2, decision 57
    of line NP's series M, 2026-11-23, Owner-delegated), since the rows of a query there depend on
    its place in the batch.
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
      the one-hit scan (3 below their start, 3 past their end); for the query `qova` (`qova.faa`,
      `qova.fna`) its near input, a spacer of 357 `X` at each of the six boundaries, gives NCBI's own
      bytes of `ov_atha` in every run (BLASTP and BLASTX; the default, `-evalue 1000`,
      `-window_size 0`, `-ungapped`; 1 and 4 threads).
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
      BLASTX with `qfc.fna` (10,002 `A`, no word, then a query of 75 letters; item M4b-3b2): one
      batch at `BATCH_SIZE=30000` 18 rows of the query (`-ungapped` 20), where the query first in
      its batch (the default batch) has 19 (21): LOSAT's units of 10,002 letters wrote 19 and 21
      at 1 and 4 threads, the whole batch as one unit NCBI's bytes.
    - A. thaliana TAIR10 (48,321 proteins, segmasker's defaults, algorithm 21) for 55 queries, the
      proteins that such a protein's entries reach or that reach below their start (BLASTX: the first
      20, back-translated; ws, `M/m4b3b/q5b/tair`, `q5b/ws/q5b_threads.sh`): LOSAT wrote the same bytes at
      1, 4 and 16 threads in every case; at the default window LOSAT's bytes are NCBI's (8,716 BLASTP
      rows, 1,185 BLASTX rows); with `-window_size 0` NCBI gives the same bytes at 1 and 4 threads and
      LOSAT 3 more BLASTP rows (AT1G47340.1 and AT1G50880.1 against AT2G17310.1 from its first letters,
      E 0.012 and 0.001; AT2G27380.1 against itself at 750-761, E 1.1) and 1 more BLASTX row (the first
      of these), whose one-hit seeds NCBI skips at entries of the protein searched before; no row of
      NCBI's is missing. (Measured before fix rounds 2 and 3 of M4b-3b and again after fix round 7,
      on 0c6e07281 with `q5b/ws/q5b_threads2.sh`, 2026-11-23: the same rows and bytes.)
    Fixture rows `db.pmask.{bp,bx}.hard.ov_lk.*` (NCBI's bytes at the default and at
    `-window_size 0`; with `ncbi-near:` those of NCBI's run of the near input, rows
    `db.pmask.{bp,bx}.hard.near.ov_lkn.*`, otherwise), `db.pmask.{bp,bx}.hard.near.ov_lkx.21.w0.6`,
    `db.pmask.{bp,bx}.hard.ov_zt.*`, `.ov_ztl.*`, `.ov_zto.*`, `.ov_th10.*`, `.ov_thm20.*`,
    `.ov_th10f.*` and `.ov_thm20f.*` (NCBI's one-thread bytes; the `.q1q.` rows of two queries in
    one batch and, `BATCH_SIZE=50`, in two, with those of `ov_s2`, `ov_nb`, `ov_zt` and `ov_ztl`,
    the `.q4.` rows of four queries in three batches and `.ov_th10.qbig10000.` of two batches at
    BLASTP's default size, also with LOSAT_X_BPQBLOCK 0, 1 and 120 in
    `LOSAT/tests/blastp_prelim_blocks.rs`), `db.pmask.{bp,bx}.hard.ov_m1.*` (NCBI's one-thread bytes;
    BLASTX also in waves and in tasks of 64 subjects at 2 and 4 threads,
    `LOSAT/tests/blastx_subject_flow.rs`: the subject after a sequence on a new table in the same
    worker, review 4 of M4b-3b, L1) and
    `db.pmask.{bp,bx}.hard.ov_atha.21.w0.6` (NCBI's bytes) of
    `LOSAT/tests/db_search_manifest/{blastp,blastx}.pmask.tsv`, checked at 1 and 4 threads; the
    databases `ov_lkn`, `ov_lkx`, `ov_zt`, `ov_ztl`, the alias `ov_zto.pal`, `ov_th10`, `ov_thm20`,
    `ov_th10f`, `ov_thm20f` and `ov_m1` of `LOSAT/tests/db_search/dbmaskp/` (`make_dbmaskp_ov.py`:
    `near_db`, `blank_db`, `probe4`, `probe5`, `probe6`). Unit tests `each_sequence_is_searched_on_its_own`,
    `a_new_table_holds_entries_below_the_start`,
    `ncbi_walk_names_the_sequences_that_meet_unwritten_entries`, `ncbi_walk_follows_the_clears`,
    `hard_starts_by_ordinal`, `starts_and_ends_apply_to_the_table`,
    `a_plain_sequence_on_a_new_table_ends_past_the_unwritten_entries` and
    `sequences_whose_hits_meet_other_entries_are_kept_apart`
    (`LOSAT/src/algorithm/common/hard_chunk_walk_tests.rs`), the tests of
    `LOSAT/src/algorithm/common/query_batch_tables_tests.rs`, `hard_mask_starts_and_ends_move_the_table`
    (`LOSAT/src/algorithm/blastp/blast_engine_tests.rs`, `LOSAT/src/algorithm/blastx/seed_tests.rs`).

    For TBLASTN and TBLASTX (version 1.106, line NP item M4b-3c part 2), where a database has a
    sequence whose stored masks overlap or nest in the order of their start and leave a reading
    frame a negative chunk, `blastinput::db_mask::translated_walk` runs `NcbiWalk::translated_masked`
    over all six frames of every searched sequence in database order (`HardStarts`): a sequence
    whose walk meets entries that no sequence wrote or clears inside it starts on a table of its own
    at NCBI's one-thread offset (`NcbiStart`), and any other walked sequence starts on a table of
    its own at `window + below_start` with every other sequence's entries out of reach; all six
    frames of such a sequence run in order on that one table (`tblastn/search_gapped.rs`,
    `tblastx/blast_engine/frame_chunks.rs`). Databases without overlapping masks never run the walk.
    Evidence (NCBI BLAST+ 2.17.0 under a 4 GB limit and 30 s; series M `M/m4b3c2/ncbi`):
    - `lk` (A with masks [40,200) [60,238) in frame +1 residues, then B with [230,235)) at 1 thread:
      A's one-hit walk moves 628 positions below its start and leaves entries 627 past its end,
      skipping 25 (`tblastn -window_size 0`), 30 (`tblastn -ungapped -window_size 0`) and 132
      (`tblastx -window_size 0`) of B's rows compared with B searched first (`lkBA`, same database
      size); all but 1 row in each search come from B's start at -588 on unwritten entries (`f2b`'s
      B, after F2 which moves -588 while storing no entries: 42 / 44 / 215 rows vs `lk`'s B:
      41 / 43 / 214), so A's own entries skip 1 row of B in each one-hit search.
    - `lkBA` (B then A) at 1 thread: B's entries skip 4 / 4 / 5 of A's rows (`lkBA`'s A: 6 / 8 / 55
      vs `pa`'s A, after P which moves +952 storing no entries so no entry reaches A: 10 / 12 / 60).
    - `lk400` (200 copies of A then B): one-hit rows vary with threads and across runs (`tblastn`
      ungapped one-hit 1,800 rows at 1 thread vs 1,812 at 4 threads, 2 distinct outputs in 3 runs;
      `tblastx` one-hit 7,586 vs 7,616); two-hit `lk`, `lkBA` and `lk400` give the same rows at 1
      and 4 threads (400 / 2,200 in `lk400`).
    - `fb` (F of 722 bases moves offset to -1410 storing no entries, then B): 0 rows at 1 thread vs
      66 / 73 / 346 at 4 threads; `f2b`: 42 / 44 / 215 at 1 thread vs 66 / 73 / 346 at 4 threads;
      `pa`: 10 / 12 / 60 at 1 thread (4 / 17 in a 4-thread run). LOSAT writes NCBI's 1-thread bytes
      for `fb`, `f2b` and `pa` at every thread count, and on `lk`, `lkBA` and `lk400` writes each
      sequence's independent rows identically at 1, 2, 4 and 8 threads (`LOSAT/tests/t_db_hard_overlap.rs`,
      `LOSAT/src/algorithm/common/hard_chunk_walk_tests.rs`).

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

## Approved exception 25 (makeblastdb `-gi_mask`, a database of one volume, version 1.81)

Version 1.81, line NP item M4c-1. The Owner answered on 2026-11-20 「Q1=A Q2=B Q3=A Q4=A Q5=A Q6=B Q7=A
Q8=OK」. Q6 was the question Q-GIMASK 「makeblastdb -gi_mask が 1 volume の DB を作る時、NCBI と同じにしますか」
(when `makeblastdb -gi_mask` builds a database of one volume, should LOSAT do as NCBI does?), and its option
B read 「LOSAT だけが 1 volume でも読める DB を書く。NCBI と違う file になり、新しい例外が 1 つ増える。」 (LOSAT
alone writes a database of one volume that can be read; the files differ from NCBI's, and one more exception
is added). Not taken: option A, NCBI's files and its exit status 2 under rule 2. The form below (which bytes
differ, and how) is Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07. Not taken there:
naming the volume `<name>.00` as a database of several volumes does (the name of every volume file would
differ from NCBI's), and keeping NCBI's alias file while only LOSAT's reader takes it for the volume (no
file would differ, and NCBI's programs still could not open the database).

25. **The alias file of a database of one volume built with `-gi_mask`.** With `-gi_mask` NCBI's
    writer gives every database an alias file (`writedb_impl.cpp:378-380`); its `DBLIST` is the list
    of the volumes' file names when there are several, and the database's own name `m_Dbname`, an
    absolute path, when there is one (`CWriteDB_Impl::x_MakeAlias`, `writedb_impl.cpp:403-416`).
    NCBI's reader takes a `DBLIST` entry with a directory for another alias file or volume to look
    up; this one names the alias file itself, so every reader of the database stops with `BLAST
    Database error: Illegal configuration: DB alias files are mutually recursive.`
    (`CSeqDBAliasNode::x_ExpandAliases`, `seqdbalias.cpp:735-813`). makeblastdb opens the new
    database for its metadata file (`makeblastdb.cpp:1296-1317`) and ends there with that message
    and exit status 2, after every other file is written; no NCBI program can open the database.
    LOSAT writes the database's file name without a directory in that line
    (the file name of `m_Dbname`, the form of the list of several volumes), which NCBI's reader
    takes for the volume of the alias file's own name (`seqdbalias.cpp:773-797`): the database
    opens, and the build goes on as a build of several volumes does, with the metadata file and
    exit status 0. Every other file, the messages and the log are NCBI's. The metadata file is the
    one NCBI's own code gives for that database (`CSeqDB::GetDBMetaData`, `seqdb.cpp:1610-1749`;
    NCBI BLAST+ 2.17.0 `blastdbcmd -metadata` on NCBI's database with that one line changed): it
    lists the alias file twice, as the alias file and as a file of the volume (the extensions of
    `SeqDB_GetFileExtensions` begin with the alias file's, `seqdbcommon.cpp:2629`), and counts its
    bytes twice. An alias file that NCBI's makeblastdb left from such a build is not LOSAT's: the
    removal of the old database cannot read it (`seqdb.cpp:1542, 1573-1574`), a build with
    `-gi_mask` writes its own alias file over it, and a build without `-gi_mask` keeps it and ends
    with the message above, as NCBI's does (rule 2).
    Evidence: `docs/evidence/gi_mask_single_volume/` (NCBI BLAST+ 2.17.0 opens such databases,
    version 5 and 4, nucleotide and protein, `-out` in a subdirectory, the default GI-mask name:
    `blastdbcmd -info`, `%m0`, `-mask_sequence_with 0`, `-metadata`, `blastn -db_soft_mask 0` and
    `-db_hard_mask 0`, `blastp -db_soft_mask 0`; the masks are those of the mask column but for the
    sequences without a GI, as for a database of several volumes; from another directory NCBI does
    not find the GI-mask files, as for a database of several volumes). Fixtures: the rows of one
    volume of `LOSAT/tests/seqdb/gi_mask` (`make_gi_mask_fixtures.py`, `single_volume_exception`:
    NCBI's run with that one line changed, taken as the successful run where NCBI ended on the
    message above), among them `x01_*` (a rebuild over LOSAT's alias file) and `x10_*` (builds
    beside NCBI's alias file). Tests `LOSAT/tests/makeblastdb_gi_mask.rs`.

## Reproduced as NCBI (rule 2)

- RPSTBLASTN `-use_sw_tback`: a Smith-Waterman alignment that starts with a gap (version 1.80,
  2026-11-22, line PS item PS12c, review B1; ported in version 1.91, 2026-11-23, line PS item
  PS12d; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07): with forbidden
  ranges NCBI's redo can give a profile's second alignment that starts with a gap in the profile's
  row (rows `db.rpstbn.sw.nest.*`: `nest.faa`, the profile's source twice, nested,
  back-translated, on both strands and as two queries of one batch). Its Std-seg's first segment
  has no location in that row, so `CSeq_align::GetSeqRange` takes no id from it and throws "Row
  seqids not consistent" on the row (`Seq_align.cpp:171-191`) wherever a report reads it. The
  output is deterministic: the prefetch, which reads the profile's row of every Seq-align, warns
  (every format but 8, 9, 11 and 12; not with `PRE_FETCH_SEQS_LIMIT=0`), outfmt 0 to 4 stop after
  the query's `Length=` (the description table reads both rows) with exit status 255, 13 and 14
  too, BLAST XML, JSON and XML2 stop their document with the exception (exit status 0), and the
  tables, the Seq-aligns and outfmt 0 to 4 without the description table are the alignments as
  made, the leading gap shown as end gaps (spaces) in outfmt 0. LOSAT writes these reports as
  NCBI's: one predicate of the rows whose range `GetSeqRange` cannot take (an empty row, a row
  that starts with a gap; `report::serial::seqalign::std_seg_row_range`) stops each reader of a
  row where NCBI's stops (BLASTX's owner of these rows, which RPSTBLASTN's reports share). An
  alignment that starts with a gap in the query's row was not found (series PR item PS12d's grids
  of profiles and queries); LOSAT refuses it, with the query batch that holds the profile, as the
  endless redo of 1.80 (below). rpsblast's Dense-seg has a location in every segment; the X-drop
  redo's alignments start with a substitution.

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
- The header record of a database sequence as NCBI's object stream reads it
  (`CObjectIStreamAsnBinary`, `x_GetHdrAsn1`), version 1.90 (line RG item T5d): a constructed
  value with a definite length in the long form ends where NCBI's `ExpectIndefiniteLength` puts
  it, at the position after the first length byte plus the length (objistrasnb.inl:296 takes the
  stream position before `ReadLengthLong` reads the other length bytes):
  one byte per further length byte before the value's end, so that NCBI 2.17.0 fails on such a
  record with `EndOfContentsByte expected` (probes c30 and c31 of
  `remaining-gaps-20261112/T/t5d/impl/probe1`); makeblastdb writes indefinite lengths. The bytes
  after the set are not read (`EndOfRead`; probe c25). `seqid::ber::decode_full_set` reads as
  NCBI does, and the ordinal-id decoder of `seqdb::header` fails wherever it does; unit tests
  `header_reads_fail_as_ncbis_object_stream_fails` and
  `definite_lengths_read_as_ncbis_reader_reads_them` (`LOSAT/src/seqid/tests.rs`).

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
  subject of an ungapped search (blastx and blastp `-ungapped`) is fetched by none of them and
  leaves at the Seq-aligns, as a positive list's does: it keeps its place in the hit list through
  the traceback's cut to the hit list size, its culling and best hit pipes (its HSPs cull those of
  other subjects) and the prune (`blast_traceback.c:1705-1777`, `blast_seqalign.cpp:1582-1586`;
  version 1.92 for BLASTP, rows `db.ll2.bp.ug.*`, `db.ll2.bx.ug.*`). Versions 1.6, 1.29 and 1.32 made this an explicit
  rejection; version 1.34 reproduces it in the six search programs, each of which removes such
  subjects where NCBI fetches them (`SeqDb::excludes_lineless_at_fetch`,
  `algorithm::common::database::FetchExcluded`; line DBID item lineless2). Under such a list
  NCBI's fetch reads the header of every subject it fetches, so a header that cannot be read
  (a damaged record, a `BL_ORD_ID` line with a string tag in a later volume) ends NCBI's search
  there (exit 255 in the traceback, `BLAST engine error` in the ungapped re-evaluation) even
  where the subject is not reported; LOSAT refuses such a search at the same place. An ungapped
  search of a protein database reads such a header at the Seq-aligns only, where the subject is
  in the final hit list: there NCBI's Seq-align fails (exit 255), and LOSAT's report gives the
  same error (rule 2). The composition-based redo (and the Smith-Waterman redo) fetches only the
  matches that its early termination test leaves: since version 1.93 such a subject's list stays
  in the redo's stream (also one whose lines LOSAT refuses where it reads them) and LOSAT refuses
  the search only where NCBI's redo fetches it, after the reports of the query batches before (as
  NCBI's search ends in that batch), so a search whose redo frees that match unfetched writes
  NCBI's output (BLASTP, PSIBLAST, BLASTX, TBLASTN; line RG item D1c-1c, rows `db.il.*.lr*`,
  `db.il.*.ls*`).

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
- A database sequence whose Seq-align id is a pdb id that its Seq-id handle holds with another
  chain (version 1.121, 2026-11-24, line RG-KB item 12, Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07). The handle normalizes a pdb id (`CSeq_id_PDB_Info::Normalize`,
  `seq_id_tree.cpp:3381-3464`): where it rebuilds the id (a release date that normalizes, a
  lower-case letter in the molecule, or a chain and chain-id that disagree) the chain becomes the
  first byte of the chain-id, or the chain's low byte, read back as a signed char, and `GetSeqId`
  (`:3293-3346`) gives back the stored id with that chain. NCBI's data loader looks the Seq-align
  id up as that Seq-id and serves the sequence only where one of its stored ids compares equal to
  it (`bdbloader.cpp:386-431`); `CPDB_seq_id::Match` compares the chains where both are set
  (`PDB_seq_id.cpp:53-67`). So a stored chain outside -128 to 127 without a chain-id (321, -191,
  383), or a chain-id whose first letter is not the chain with a rebuilt id (`1ab1`, chain 66 and
  chain-id `A`), leaves the report's Bioseq handle empty, without an exception, where no other
  definition line of the sequence holds the handle's id. Outfmt 0 writes the description line of
  the `catch` without its links (`x_FillDeflineAndId` returns at once, `showdefline.cpp:388-394`):
  `1AB1_A Unknown`, and `Sequence with id 1AB1_A no longer exists in database...alignment
  skipped` for each alignment, no heading counted in the `-html` ranks after it
  (`showalign.cpp:1952-1981`); outfmt 1 to 4 leave the subject out of the merged alignment and
  write a line feed alone where every subject is left out (`:1997-2003`, `:2086-2110`); XML the
  `Hit` of the `catch` (outfmt 5 `Unknown`, `Hit_len` 0, `blastxml_format.cpp:399-407`; XML2 and
  JSON the `HitDescr` of the id alone and the sequence's own length,
  `blastxml2_format.cpp:410-415`, as `GetSequenceLength` skips the comparison,
  `bdbloader.cpp:312-319`) and the `Hsp` of the `catch` (the scores, the Seq-align's identities,
  the coordinates 0 and no letters, `blastxml_format.cpp:142-165`, `blastxml2_format.cpp:130-153`);
  the tabular formats the rows of a failed lookup above (`SetSubjectId` clears the subject id and
  throws, or `GetBioseqLength` does, `tabular.cpp:833-875`); outfmt 18 no row
  (`taxFormat.cpp:996-998`); psiblast's PSSM reads none of its letters. A release date
  normalizes where it is a Date-std of year 1 to 4095, month 1 to 15, day 1 to 31, hour 0 to 30,
  minute and second 0 to 62 and no season (`x_NormalizeDate`, `:3381-3405`). The same records with
  a chain of 65 or none, a date of month 16 or a season, or a normalized chain that is the stored
  one (chain 65 and chain-id `A` with a date of month 13), are served as before. LOSAT finds it from
  the decoded Seq-ids alone (`ParsedSubject::scope_finds_none`, the handle's id
  `seqid::mask_handle::pdb_handle_seq_id`) and writes these outputs in BLASTP, BLASTX, TBLASTN,
  TBLASTX, BLASTN and psiblast (rows `db.pdbchain.only.*` and `db.pdbchain.norm.*`: a lower-case
  molecule, a date of month 13 and a chain of 200 without a chain-id, and the dates and chain
  that load). Where another id of the same sequence compares equal to the handle's id, NCBI's load
  makes the Bioseq of the first definition line that holds it as `s_SeqDB_SeqIdIn` reads the ids,
  which stops at the first id of the line that does not match (`seqdbvol.cpp:1328-1359`): another
  definition line than the Seq-align id's (probes `e15_twoline` of review 1 of item 12 and `e15nt`
  of review 2) makes the Bioseq of the reports of that line, and the reports read it in two ways.
  Those that read the Bioseq's line take that line: its title is NCBI's description line's (outfmt
  0 to 4), the generated title of the organism report (outfmt 18) and the generated title of a
  line without one in the alignment's heading (`showalign.cpp:2250-2255`, `taxFormat.cpp:915-926`),
  and its ids are the tabular `sseqid`, `sacc` and `saccver`, the labels of outfmt 1 to 4, the
  SAM name and the `-html` link. Those that read the Bioseq's whole set of lines list the
  sequence's own lines: the alignment's heading, the tabular `stitle`, `salltitles` and
  `sallseqid`, and the XML `Hit` (outfmt 5, 15 and 16; `tabular.cpp:415-439`,
  `showdefline.cpp:248-293`; probes `e15g` and `e15tu` of review 3). LOSAT takes the first group
  from that Bioseq's line (`ParsedSubject::target`, the line of `DeflineSet::bioseq_defline` of
  the handle's id; review 2, N3, where review 1 had refused the description line and the
  organism report) and the second from the stored set, as NCBI does; no definition line (one
  line with the chain 321 before the chain 65, `e15b_oneline`) makes NCBI's load throw, in every
  format that reads the Bioseq or prefetches it (outfmt 0, 1, 2, 5, 6, 7, 10, 15, 16 and 18 of the
  probe), which LOSAT refuses
  (rule 3 below), while the Seq-annots and the archive (outfmt 8, 9, 11, 12) read nothing of it
  and are NCBI's bytes. A handle whose id another sequence's definition line holds, where NCBI's
  loader serves that sequence (`share1` of the probe), stays the refusal of a Seq-align id that
  another sequence holds.
- makeblastdb, an index file that a failed write left cut short (version 1.62, 2026-11-21, line LIB-O
  item O1e, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07). NCBI opens the
  database it wrote and reads the index through its mapping, where bytes the file lacks read as
  zeros to the end of its last page. LOSAT reproduces these outcomes: "Not a valid version 4 or 5
  database.", "requested sequence type does not match DB.", a header-less database and, for a date
  that was cut, the exception of CTime (exit status 255, no metadata file). The rule-3 refusals of
  an index file of no bytes and of a header field past the last page are in the list below.
- makeblastdb `-gi_mask` (version 1.81, 2026-11-23, line NP item M4c-1): the GI-mask writer
  (`seqdb_writer/writedb_gimask.cpp`, `writedb_impl.cpp:338-342, 405-438, 1453-1456`) is ported
  with its defects, which give deterministic files, messages and exit statuses, but for the alias
  file of a database of one volume (approved exception 25). An alias file of that kind that NCBI's
  makeblastdb left in the directory is kept by a later build into the same `-out` (the removal of
  the old database cannot read it and ignores the error, `seqdb.cpp:1542, 1573-1574`); a build
  without `-gi_mask` then ends at its metadata step with `BLAST Database error: Illegal
  configuration: DB alias files are mutually recursive.`, exit status 2, as NCBI's does. The masks
  of a sequence without a GI are dropped before their ranges are checked (a range past such a
  sequence's end is no error), and an algorithm that reached no sequence with a GI writes no files
  and is left out of `MASKLIST`. A data volume takes its first record
  whatever its size, and when the first record does not fit `-max_file_sz` the volume `.00` is
  never created. A `-gi_mask_name` in a directory that does not exist writes nothing there and is
  still listed in `MASKLIST`. GIs from 2^31 to 2^32 - 1 are written as the low 32 bits
  (NCBI converts a GI to `int` with `GI_TO`). The size of a record in `CanFit` (`writedb_gimask.hpp:96-99`: twice the
  number of ranges plus one, times four, which NCBI computes in `int` and which overflows past
  268,435,455 ranges) is computed without overflow. Two algorithms whose GI-mask files are the
  same files are written as NCBI writes them where their files are written one after the other
  (the rule-2 rejection of version 1.81 below covers the others). Fixtures `LOSAT/tests/seqdb/gi_mask`
  (`make_gi_mask_fixtures.py`), tests `LOSAT/tests/makeblastdb_gi_mask.rs` and
  `LOSAT/src/writedb/gi_mask_tests.rs`.
- The readers of GI-indexed masks (version 1.94, 2026-11-23, line NP item M4c-2; review round
  2026-11-24): `blastdbcmd` (`-info`, `-outfmt %m`, `-mask_sequence_with`), the search programs
  `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx` (`-db_soft_mask`, `-db_hard_mask`),
  `makeblastdb -input_type blastdb` (the masks copied into the copy's mask column) and
  `makembindex` (`-db_mask`, `-show_filters`) read the masks that an alias file's `MASKLIST` names
  as NCBI's reader does
  (`seqdbalias.cpp:117-139`, `seqdbgimask.cpp`), with its deterministic defects: the masks apply
  only when `-db` is one name whose own alias file has `MASKLIST` (in an alias file over the
  database or with several names they are passed over without a message); the mask files are
  looked for from the current directory and `BLASTDB`, not from the database's directory
  (`Could not open gi-mask index files.` from elsewhere); a sequence without a GI takes the masks
  of GI 4294967295 when that GI is in the files (`(Uint4)` of `INVALID_GI`); with equal GIs the
  record is the one NCBI's binary search reaches; a data volume that was never written gives
  NCBI's memory-map text with the empty name (`To be memory mapped the file must exist: ''`);
  the algorithm ids and names are those of `MASKLIST`, with NCBI's names table in the errors of
  an unknown id or name; makembindex takes an id or a `MASKLIST` name, and its new-style index
  opens each volume alone, which names no masks, so on a database of several volumes it stops
  with `unrecognized filter algorithm id`, while an old-style index and a database of one volume
  read them (`sequence_istream_bdb.cpp:62-96`, `mkindex_app.cpp:268-296`). The `.nog` of a volume
  (the GI of each sequence), which NCBI reads through its mapping with no check of its length
  (`seqdbvol.cpp:52-71`), reads as that mapping does: bytes past the end of the file read as zeros
  to the end of its last page, so a file cut inside its header gives a count of 0 and every
  sequence takes the first GI of its definition lines, as without the file; a cut-off entry is GI
  0, which no mask has; a count below the volume's sends the sequences after it to their
  definition lines. A protein in FASTA with such masks is approved difference 5 of
  PD-LOSAT-BLASTDBCMD-DEVIATIONS, as with mask columns. blastdbcmd, makeblastdb and makembindex
  read the masks on one thread (they have no `-num_threads`). The five search programs read the
  masks of the sequences their threads fetch (`seqsrc_seqdb.cpp:316-321`,
  `seqdbgimask.cpp:65-112`), and with `-db_soft_mask` or `-db_hard_mask` and more than one thread
  NCBI stops in some runs with `BLAST engine error: Database memory map file error` (exit status
  3) or crashes (SIGSEGV, and at 8 threads glibc's double-free or unaligned-chunk abort): in the
  fixtures' masked searches, blastn failed in 4 to 18 of 60 runs at 2, 4 and 8 threads, and in the
  probe of 2026-11-24 (review S1) tblastn in 3 of 60 runs at 4 threads and 10 of 30 at 8, blastp
  in 3 and 2, blastx in 2 and 2, tblastx in 1 and 10, blastn in 8 and 4; the other runs gave the
  one-thread bytes, and with one thread every run of each program (20 of 20), as with mask columns
  at any thread count (60 of 60), gave the same bytes. LOSAT writes NCBI's `-num_threads 1` output
  at every thread count in all five, as every `-db` search does (PD-LOSAT-DATABASE-SEARCH, DB-3).
  The explicit rejections of an empty `.nog` and of a `.nog` entry past its mapping are in the
  rule-3 list; a `MASKLIST` name that is not UTF-8 is read by its bytes since version 1.103
  (`LOSAT/tests/non_utf8_names`, cases `K10_*_gimask_*`). Fixtures
  `LOSAT/tests/seqdb/gi_mask_read` (`make_gi_mask_read_fixtures.py`, NCBI BLAST+ 2.17.0: the rows
  `bdc_wide_*` of a data file past 64 KiB, `bdc_sv4n*` and `bn_sv4n*` of the `.nog` cut or edited,
  `mbi_*` of makembindex), tests `LOSAT/tests/gi_mask_read.rs` (each search row again with
  `-num_threads 4`) and `LOSAT/src/seqdb/gi_mask_tests.rs`.
- Defects without an effect on the output that LOSAT follows: the negative `Int4`
  offset of `BlastGetStartForGappedAlignmentNucl`, the out-of-range read of
  `hard_ranges[n]` after the last chunk (version 1.83; version 1.82 and before named `[1]`; for a
  masked subject `n` is the number of hard ranges, a read 4 bytes past the array of
  `CSeqDB`; the overlap it chooses only merges HSPs that intersect across a chunk boundary, and
  hard chunks are disjoint), and the catch in `BlastSetupPreliminarySearchEx` that never
  matches.
- The out-of-range `double` to `Int4` conversion (`INT_MIN`) of BLASTN's `BLAST_Cutoffs`
  (`BlastKarlinEtoS_simple`, `blast_stat.c:4041-4064`), listed until 2026-11-28 among the defects
  without an effect, does change the output. With an extreme `-searchsp` and `-evalue` (rows
  `db.etos.bn.int4min.ss.e300`: `-searchsp 9223372036854775807 -evalue 1e-300`) the quotient
  `K * searchsp / E` is past the range of a double, NCBI's x86-64 build converts it to `INT_MIN`,
  the hit-saving cutoff stays 1 and so does the word cutoff, every seed is extended, and NCBI
  BLAST+ 2.17.0 writes an alignment that the finite value of exception 8 would not find (one row
  against none; with `-evalue 1e-200` the quotient is finite and NCBI writes none,
  `db.etos.bn.finite.ss.e200`). LOSAT reproduces it (`stats::cutoffs::CutoffConversion::BLASTN`),
  as decision D4 reproduces the conversions of `-threshold` and the X-drops. Whether it stays so
  is the maintainer's open question Q-ETOS8 (line RG item X-etos-owner), with the other programs'
  conversion of the same quotient.
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
  (`SplitQuery_CalculateNumChunks` gives more than one chunk): NCBI's 64-bit `size_t` chunk ranges
  wrap (`split_query_cxx.cpp:145-171`). Where the wrapped chunk size covers the batch, or where
  `chunk_start` reaches or passes the batch's end before the last chunk (for example `-100/-200`
  on 310 bases, `-133/-200` on 200 bases, or `-132/-200` on 201 bases), later chunks hold no query
  and NCBI stops in `CQuerySplitter::Split` with a null-pointer `CCoreException` (exit 3), handled
  by the empty-chunk rejection below (a negative pair never produces a non-empty chunk that would
  be split again). Otherwise every chunk holds a query with gaps between them (for example
  `-100/-200` on 248 bases, `-1/-200` on 930 bases): NCBI searches on and loses hits in the gaps.
  Version 1.1 (maintainer, 2026-10-02): that gap case, although NCBI's result is deterministic,
  stays an explicit rejection (a setting without a practical use; reproducing it would need NCBI's
  wrapped `size_t` arithmetic through the chunk ranges, masks and merge). Where such a pair does
  not split the batch (for example `-1` above `-2147483648`), NCBI searches as without the
  variables and so does LOSAT (fixtures `env.negative_pair_*`). Since versions 1.110 (line RG item
  X-short-bat-2, rows `db.xs.blastn.bat.*`) and 1.116 (line RG item X-chunk, rows `db.xc.*`), in
  BLASTN, as in TBLASTN (1.97) and BLASTP (1.98), the gap rejection comes after earlier batches'
  reports and prologs, at the first chunk with a valid context with subjects in hand, after that
  chunk's masking, setup checks (`e_nF`, average subject length) and database sequence-file read
  (`nonsq`); without subjects in hand (empty alias `e_nA`, seqid list leaving no subject) or when
  no chunk has a valid context (`gapvalid`, whose valid letters lie only in the gap, which reaches
  traceback setup: `Traceback failed` on short subjects, `Database memory map file error` on OID 0,
  or exit 0 with no hits), the batch is reported as NCBI reports it.

- A `CHUNK_SIZE`/`OVERLAP_CHUNK_SIZE` pair whose chunk ranges leave a chunk without a
  query: a positive pair with evened chunk size below or equal to overlap at batch length `<=`
  overlap (for example `140/100` on 90 bases, `150/100` on 100 bases, `100/95` on 31 bases), a
  negative pair whose wrapped chunk size covers the batch (`-100/-200` on 310 bases), or an
  unwrapped negative pair whose `chunk_start` reaches or passes the batch's end before the last
  chunk (`-133/-200` on 200 bases, `-132/-200` on 201 bases). NCBI's `CQuerySplitter::Split`
  (`split_query_cxx.cpp:234-249,872-880`) stops with a null-pointer `CCoreException` (exit 3) at
  every thread count before any chunk is masked, set up, or searched, after earlier batches'
  reports and prologs (in `-outfmt` 0, 15, 16 and 20 on batch 1, and all formats on batch 2, with
  `-subject` outfmt 15/16 holding back the previous query). Since versions 1.110 (line RG item
  X-short-bat-2, rows `db.xs.blastn.bat.empty.*`, `db.xs.blastn.bat.ord.empty.*`) and 1.116 (line
  RG item X-chunk, rows `db.xc.*` of `blastn.xc.tsv`, `tests/blastn_chunk_stops.rs`; round-3 audit:
  27 oracle cases, plus 461 official 2.17.0 probe runs across all four tasks, `-db` and `-subject`,
  1 and 4 threads), LOSAT rejects "a query chunk without a query" at that exact step in BLASTN,
  also with no subject in hand.

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
  which no later batch is prefetched) or outfmt 5 (since version 1.96, line DBID item pfoff-4: each
  `Hit` of such a sequence is the Seq-align id's FASTA form, `Unknown` and length 0, with HSPs of
  its scores, coordinates 0, the Seq-align's identities and no letters, and the report goes on);
  in each of these formats, and in psiblast's PSSM read (in every format where the query's batch
  reads one by one: its prefetch turned off by `PRE_FETCH_SEQS_LIMIT` or, in outfmt 8, 9 and 11
  without `PRE_FETCH_SEQS_LIMIT=2147483647`, by the format; after the round's reports and the reads
  before it in the PSSM, in the order of the hit list; where the subject range of an HSP that the
  PSSM includes ends within the holder's letters, which NCBI's read takes, while a shorter holder's
  read throws and NCBI writes `Failed to retrieve sequence`; refused because LOSAT does not write
  a holder's letters into the PSSM; since version 1.96, at the merge of line PS item PS15b-1 and
  review 5 of item pfoff-4, H1 to H3), it refuses such a sequence where the run loaded before it,
  in NCBI's order and while NCBI keeps it, another sequence of the database whose own lookup
  does not fail
  and that NCBI gives for its Seq-align id without the lookup (listed in the LMDB files that have
  the accession table), from one of two places: NCBI's data loader keeps the ids of every
  definition line of every sequence that it loaded for the whole run and gives that sequence for
  an id of the same type, accession and version, in any letter case (`gi|5001|ref|NP_100001.1|`
  for `ref|NP_100001.1|`, also for `ref|np_100001.1|`), and NCBI's scope keeps the sequences that it loaded, each with the ids of
  the definition line that its load made its Bioseq of (the line of its own Seq-align id: its
  first line with a gi, else its first line; since version 1.96, review 3 of item pfoff-4, M1),
  until the end of the query batch (in outfmt 5 until the end of the run) and gives one of them
  also for an id of another type of the same accession and version (`emb|NP_100001.1|` for `ref|NP_100001.1|`) or for a weak match (the versions of an
  accession without one, `ref|NP_100001.1|` for `ref|NP_100001|`; an id of the same name,
  `gb|AAA09999.1|LOCUSN` for `gb||LOCUSN`; in a version 5 database LOSAT cannot list the weak
  matches of a name or of a `pir` or `prf` id and takes any sequence loaded in that time); before
  it is a report before (a query, a batch, a psiblast round, a call of the library's batches
  before, of which LOSAT has no record and takes any such sequence, each call a batch of its
  own), or a place before it in its report (a row, a `Hit`, a description line before the
  alignments; of two sequences with equal scores NCBI's reports show the later one in the
  database first, and NCBI orders the databases of a `-db` list by name), with the words of
  another sequence that holds the Seq-align id (`... which database sequence N also holds (where
  both are reported, NCBI BLAST+ 2.17.0 can write that sequence's names and letters for it)`) or
  of holders that LOSAT cannot list: NCBI writes that sequence's record for it (since version
  1.96, review 1 of item pfoff-4, H1 and H2: `gi|5001|ref|NP_100001.1|` loaded in the first
  query's report for `ref|NP_100001.1|` in the second's, also in outfmt 0, 6, 7, 10 and 18,
  where versions 1.47 to 1.76 wrote `Unknown` or refused outfmt 10 for a crash that NCBI does not
  have; a holder that loads after the throw leaves it thrown, and a report after a table's first
  row that ends NCBI's run at the null pointer of `sacc` or `saccver` is not reached; review 2 of
  item pfoff-4, H1 and M1: the weak matches, which LOSAT wrote as `Unknown` in outfmt 0, 5, 6, 7
  and 18, and the scope's holders only while NCBI keeps them, so a sequence of another type or a
  weak match loaded in an earlier batch leaves it thrown outside outfmt 5, as NCBI's run does;
  review 3 of item pfoff-4, M2: with `-html` outfmt 0 ends NCBI's run at the first description
  line of such a sequence, exit 2, so a holder after that line is not reached). A lookup that
  does not fail still checks every definition line of the other sequences, for the scope's keys
  too, so LOSAT refuses more often than NCBI holds (review 4 of item pfoff-4: 800 runs of its
  probe, no hold missed); NCBI's two places for those checks are line DBID's next item.
  Where the holder comes first in its report, NCBI's output can equal what LOSAT wrote before
  (a table's row of a thrown lookup repeats the previous row's names, here the holder's), and
  those runs are refused too (review 2 of item pfoff-4, M2: 534 runs of the review's
  probe that LOSAT wrote as NCBI does before; writing them would need the holder's record, which
  LOSAT does not write). It refuses a call of
  LOSAT's library whose formats NCBI's runs of each end apart after such a sequence, with
  another exit status or another standard error, where the call has one result and one standard
  error (outfmt 0, outfmt 1 to 4, the formats that end the run there, at their rows, outfmt 15
  and 16, which post the error on standard error, and the others, those that read nothing of the
  sequence among them; since version 1.96, review 1 of item pfoff-4, M1: outfmt 5 with
  `6 qseqid slen`; review 2, M3: outfmt 8 with outfmt 0, outfmt 16 with outfmt 6), a
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

  Version 1.83 extends approved exception 3 once more, to the traceback's read past a hard-masked
  copy of a subject in BLASTP and BLASTX (the BLASTP and BLASTX paragraph of exception 3 above).

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
  that absorbs another and becomes positive). Version 1.122 (2026-11-27, line NP item M4b-3d,
  Owner-delegated, rule 3) extends the same rule-3 refusal (`chunk_merge_undefined_refusal`, typed
  `ShortChunkSeedRefusal` via the shared `core::blast_hits::merge_two_chunk_hsps` /
  `merge_subject_chunk_hsps`) to BLASTP and BLASTX `-db_hard_mask` `-ungapped -window_size 0`
  subject-chunk merges, TBLASTN `-ungapped -window_size 0` subject-chunk merges (both frames longer
  than 5,000,000 residues and `-db_hard_mask`, naming OID and frame), and BLASTN subject-chunk and
  query-chunk merges (`merge_prelim_hit_lists`) whose density-scaled merged score leaves `Int4`
  (naming OID, where NCBI's `(Int4)` float conversion at `blast_hits.c:1526` is undefined).

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
  the reports of the search that holds the query: after the prolog, or in the block search of
  many queries before the first byte of the query's block, after the reports of the blocks
  before (as outfmt 0's refusal of such a record's masks). NCBI's output before the abort is what
  its stream had flushed, which depends on its buffer (to a file or a pipe, the version line of
  the prolog), so no prefix of it is NCBI's. BLASTN and, since version 1.126, BLASTP write the
  prolog before the search, as NCBI does, so the refusal follows the prolog. The other sections of such a record are NCBI's
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
  LOSAT's BLASTX"; rows `db.x4b2.rej.ungapped.*`). TBLASTN `-ungapped` (version 1.79, line NP item
  O4c) crashes there in the same way, with `-matrix IDENTITY` too, before the redo's frequency ratios
  that fail a gapped batch (NCBI BLAST+ 2.17.0 from Bioconda, 2026-11-23: hits and none, `-subject` and
  `-db`, a first batch without a valid context reported before the crash); LOSAT refuses it as BLASTX
  does (exit 1, "-use_sw_tback with -ungapped, on which NCBI BLAST+ 2.17.0's tblastn crashes at every
  query batch it searches, is not supported by LOSAT's TBLASTN"; one producer of the text,
  `smith_waterman_match::sw_ungapped_error`; rows `db.tn.rej.ungapped.sw.*`, anchor
  `db.tn.ungapped.ac.b1.nosw.7`). BLASTP refuses it in the same way since version 1.87 (line NP item
  O4b2-2; before, before the first byte: decision D3 above), with the same text for `blastp`: NCBI BLAST+
  2.17.0 from Bioconda, 2026-11-23, crashes at the first query batch it searches and reports the batches
  before it, whose queries give no Karlin-Altschul parameters (rows `db.ung.sw.*`). With `-use_sw_tback` LOSAT's BLASTX searches a query batch
  of NCBI's reader as one unit (the batch's status, its first query's heap and the contexts that
  index the heaps are per batch), and gives the outcome of `-num_threads 1` at every thread count.
  Fixture rows `db.x4b2.*`; `LOSAT/tests/blastx_sw_tback.rs`; inputs `LOSAT/tests/fasta/x4b2` and
  the databases `LOSAT/tests/db_search/x4b2` (`make_x4b2.py`).

- RPSTBLASTN `-use_sw_tback` with composition-based statistics (version 1.80, 2026-11-21, line PS
  item PS12c, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07, within the
  class of 1.30; P2: no checkable result exists, P3: the round loop of 1.30 and 1.53 decides it,
  `core::composition_adjustment::smith_waterman::smith_waterman_redo_window`): a profile with two
  or more preliminary HSPs in a window of a query frame (`s_WindowsFromTranslatedAligns`) when the
  e-value of a score of 0 passes `-evalue` (the search space of the query's context 0 times the
  profile's K below it: a large `-evalue`, a small `-searchsp`, or `-strand minus`, whose invalid
  frame +1 leaves that search space 0, so that every score passes) and the profile's first
  position, after the window's composition adjustment, scores the window's first letter (after
  SEG) above 0. NCBI runs the loop of 1.30 (`Blast_RedoOneMatchSmithWaterman`,
  `redo_alignment.c:1438-1529`) on each window's letters, so it never ends either: NCBI BLAST+
  2.17.0's rpstblastn did not end within 120 and 300 seconds on the inputs of rows
  `db.rpstbn.rej.swendless` and `db.rpstbn.rej.swendless.minus`, whose searches without
  `-use_sw_tback` (`tks_m.fna`, and `nq.fna` on both strands) took 0.1 to 0.2 seconds. LOSAT
  refuses (exit 1, "a profile whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI BLAST+ 2.17.0 ... is not supported by LOSAT's RPSTBLASTN")
  with the query batch that holds the profile, after the reports of the batches before it (in outfmt
  0 after the prolog, which NCBI writes before it searches its first batch; the form of 1.40 and of
  BLASTP's 1.53; rows `db.rpstbn.rej.swendless.0` and, on a later batch, the probe of PS12c's
  review S2: the same bytes at 1, 4 and 8 threads). Where that cell scores 0 or less, NCBI's redo fails
  and drops the profile, which LOSAT reproduces (rule 2; rows `db.rpstbn.sw.*`, among them the
  `-strand minus` inputs on which NCBI ends, `db.rpstbn.sw.tksm.minus.6`,
  `db.rpstbn.sw.tksw.minus.6` and `db.rpstbn.sw.s5e.minus.*`).

- The tabular `sacc` or `saccver` of a row whose subject lookup fails as in rule 2 above (version
  1.37, 2026-11-16, line RG item T7, Owner-delegated): outfmt 6, 7 and 10 (and 20) with them in the
  field list, which the default fields of 6, 7 and 10 have (`saccver`). NCBI BLAST+ 2.17.0 takes the
  label of `FindBestChoice` of the empty `m_SubjectId`, a null pointer (`tabular.cpp:174-185`;
  `CCoreException::eNullPtr`, `ncbiobj.cpp:1010`), and stops with exit status 255 and a stack trace
  of its build after the row's earlier fields. Also `staxid` of such a row when it is the first row
  of its query's (or round's) report: `m_SubjectTaxId` is printed before any row set it
  (`tabular.cpp:262-268`; the member is not initialised, a different number in each run), and
  `sallacc` of a database sequence with a definition line without a Seq-id (`tabular.cpp:257-264`,
  same null pointer after earlier lines' accessions and `;`). Since
  version 1.123 (line RG-KB item 18b) LOSAT writes NCBI's bytes up to the stop, in every program:
  the rows before it and the row's fields before that field, each after its delimiter (`Print`
  writes the delimiter before each field, `tabular.cpp:1100-1108`; NCBI's stream is flushed as
  the formatter is destroyed), and then refuses (exit 1): NCBI's stack trace names the paths of
  its build and install and the addresses of its run, and the uninitialised number differs from
  run to run, so no program can write either. NCBI's run ends at the null pointer, so no later
  row, query or batch is written; after the uninitialised number NCBI goes on (exit 0). Before
  1.123 LOSAT refused before the query batch's report, and before 1.37 every such subject; its
  output was a prefix of NCBI's then too (probe of line RG-KB item 18b: blastp, psiblast, blastx,
  tblastn, tblastx and blastn). A failed lookup does not depend on the other
  subjects: one alone, the first, the last or two give the same rows. Since version 1.50 (line
  DBID item pfoff-2) also the first row of a query's report whose subject's load throws where the
  batch's prefetch is off (`DB contains no accession info.`, a version 5 list of databases with
  and without accession tables): the `catch` of the subject block leaves `m_SubjectId` empty
  before any row of the report loaded its subject, and `m_SubjectTaxId` unset; a later such row
  prints what the last loaded row set. Since version 1.121 (line RG-KB item 12) also the rows of a
  subject whose pdb Seq-align id NCBI's data loader serves no sequence for (rule 2 above): exit
  status 255 after the row's earlier fields for `sacc`/`saccver`, and also exit status 255 in
  `CAlnVec::GetBioseqHandle` (`CAlnException::eInvalidSeqId`) when no field of the subject block
  (`reads_subject_bioseq`) is requested alongside a sequence-reading field (`qseq`, `sseq`, `btop`,
  `positive`, `ppos`; outside BLASTN also `pident`, `nident`, `mismatch`), or in
  `CAlignFormatUtil::GetFrame` (`CCoreException::eNullPtr`, `ncbiobj.cpp:1010`) in TBLASTN and
  TBLASTX when a minus-strand subject alignment is reported with `sframe`/`frames`, an
  alignment-block field, and no subject-block field (plus-strand rows compute `(start % 3) + 1`
  without touching the handle); since version 1.123 through the same row-stop path (earlier rows
  and the stopped row's fields before the stop are emitted before LOSAT's refusal).
- `-show_gis` with outfmt 0 to 4 over a database sequence whose pdb Seq-align id NCBI's data loader
  serves no sequence for (version 1.121, 2026-11-24, line RG-KB item 12, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07): the description info of the empty Bioseq handle
  keeps the GI it was allocated with (`x_GetDeflineInfo` sets no GI before `x_FillDeflineAndId`
  returns, `showdefline.cpp:1576-1589`, `:388-394`), which the description line prints when it is
  positive (`:886-892`): three runs of one search gave no GI, `gi|300995824|` and `gi|395411296|`.
  LOSAT refuses (exit 1, "... whose outfmt 0 to 4 description line with -show_gis prints a GI that
  NCBI BLAST+ 2.17.0 never sets ...") where a report prints the description line of such a
  subject (within `-num_descriptions`, `common::database::refuse_described_lines`, as the refusal
  of a description id longer than 65 bytes; row `db.pdbchain.only.blastp.0gs`); with
  `-num_descriptions 0`, or where the subject's rank is past it, the output is NCBI's (row
  `db.pdbchain.only.blastp.0g`).
- A database sequence whose pdb Seq-align id's handle another id of the same sequence holds (rule
  2 above; version 1.121, 2026-11-24, line RG-KB item 12, reviews 1 and 2, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07). Where another definition line than the
  Seq-align id's holds it, the reports show that line's title as NCBI's (rule 2 above; review 1
  refused the description line and the organism report, review 2 writes them). Where no definition
  line holds it as NCBI's load reads the ids, NCBI stops with "BLAST Database error: Error: oid
  headers do not contain target gi/seq_id." (exit status 2, in the same bytes in each run, after
  the prolog of outfmt 0 to 4, 15 and 16, before any byte of outfmt 5): LOSAT refuses (exit 1,
  "... its load throws: BLAST Database error: Error: oid headers do not contain target
  gi/seq_id.") in the formats that read the Bioseq or prefetch it
  (`ParsedSubject::handle_load_refusal`), not NCBI's error and exit status: NCBI's throw comes
  from the batch's prefetch, whose stop LOSAT makes in line DBID's `prefetch_stop`, which item 12
  leaves to that line (open point of review 1, R6; line DBID's item pfoff-9 after item 12). Item
  12 keys a pdb id by its handle in NCBI's scope (`defline_set::scope_keys`), so a query
  `pdb|1AB1|A` read with `-parse_deflines` holds the Seq-align id of a sequence whose pdb id has
  the chain 321 and no chain-id: the refusal of a sequence with the id of a query
  (`common::database::refuse_id_of_a_query`, PD-LOSAT-DATABASE-SEARCH) covers it, since review 2
  (N1) only in the formats in which NCBI's scope looks the sequence up by its Seq-align id (a
  report that reads its Bioseq, or the batch's prefetch; psiblast every format where the call
  makes a PSSM from the round, review 3, F1). The Seq-annots
  and the archive (outfmt 8, 9, 11, 12) write NCBI's bytes, as on dev. NCBI shows the query's
  Bioseq for the sequence in the reports that read it (`6 sseqid stitle` of probe `c321`:
  `pdb|1AB1|A` and `N/A`), and since decision KB-12-1 (2026-11-25) so does LOSAT, from one owner
  of the run's query Bioseqs (`report::query_bioseq`) that the tables and the organism report
  read, for every sequence with the id of a query in BLASTP, BLASTN, BLASTX, TBLASTN, TBLASTX
  and a psiblast call that makes no PSSM (PD-LOSAT-DATABASE-SEARCH): the query's ids,
  `stitle` `N/A`, the query's length as `slen`, no taxonomy, and the query's line in the organism
  report. The fields and formats that read the found Bioseq's letters keep the refusal.
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
  Seq-inst (NCBI BLAST+ 2.17.0 writes outside its query buffer or crashes unless there is one letter more than the Seq-inst's length) is not supported by LOSAT's PSIBLAST", or TBLASTN; refused in the text form before):
  for a query of 30 letters NCBI 2.17.0's psiblast ends by a signal with the lengths 31, 60 and
  255 (SIGSEGV), 0, 1, 2, 10 and 4,294,967,295 (SIGABRT), and writes past its buffer with 28
  (valgrind: an invalid write) though it exits 0. With 29, one less than the letters, NCBI answers
  (exit 0, clean under valgrind) in ASN.1 text and binary alike: it searches the 29 letters with
  the PSSM, and its E-values are those of a query of 29 letters. LOSAT gives that answer since
  line PS item PS7c (the query is searched as its first 29 letters), and refuses the other lengths. Not a rule-3 item but
  a refusal: with `-out_pssm`, author names (`Auth-list` `names`) of a binary PSSM, which NCBI keeps
  in a delay buffer and reads again only while it writes them, where that read posts a message or
  fails (NCBI reports it in the middle of the file it writes; not ported). The warning of a
  UTF8String tag (803.10) is not one of them: NCBI posts it once per process, when it first reads
  the file, so the second read posts nothing and the names are written. Also not a rule-3 item
  but a refusal, of the text form too: a gapped `lambda` of `PLUS-INFINITY`, a gapped `lambda`,
  `kappa` or `h` of `NOT-A-NUMBER`, a negative gapped value, and psiblast's `lambdaUngapped` of
  `PLUS-INFINITY` were refused here while LOSAT did not compute NCBI's output of them. Since line PS
  item PS7c LOSAT does: NCBI goes on with these values (E-values 0.0, bit scores `inf` or -1.0) and
  so does LOSAT, so none of them is refused for what it is. A gapped `lambda` below 0 (a finite one)
  or so small that NCBI's search for the score of the hit-saving cutoff (`BLAST_SpougeEtoS`) doubles
  or adds its `Int4` scores past 2^31 - 1 is searched without that signed overflow: approved
  exception 31 (the Spouge cutoff search past `Int4`). Tests:
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
  `:3677-3700`), the ungapped re-evaluation (`blast_hits.c:2700-2705`; TBLASTN with `-ungapped`,
  refused there since version 1.79; TBLASTX reaches it in every search, version 1.51 below). NCBI BLAST+ 2.17.0 dies
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
  `-subject` search, translating in code 1, may find nothing there. A refusal in a later query batch
  (`BATCH_SIZE`, or queries whose letters fill a batch) comes after the reports of the batches
  before, as NCBI's crash does (`tblastn_app.cpp:275-304` prints each batch before it searches the
  next): LOSAT writes them whole (with `-subject`, outfmt 15 and 16 hold each query's report until
  the next query's results arrive, as NCBI's do, so the last query before the refused batch is left
  out), for the first batch the prologs, and leaves the epilog out (version 1.97, 2026-11-18, line
  RG item X-short-bat, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07;
  `tests/tblastn_batch_refusal.rs`, rows `db.xs.tblastn.bat.*`). TBLASTN's other rejections of a
  step where NCBI's search of a batch crashes, does not finish or stops take the same form,
  TBLASTX's (version 1.40): two `-window_size 0` one-hit alignments at the end of their subject
  chunk or past it, which NCBI's interval tree splits without end; the overflow banks of the seed
  search; a query chunk split again or holding no query (the latter also where a negative
  `CHUNK_SIZE`'s first wrapped range covers the batch; both also on a database left without a subject: a
  seqid list that keeps none, an alias over volumes without sequences; D13); and with `-in_pssm` (one batch: the prologs) a
  query whose letters give no Karlin-Altschul parameters (one that `-seg` masks completely: NCBI's
  `PsiBlastSetupScoreBlock` reads the `kbp_psi` block that `Blast_ScoreBlkKbpUngappedCalc` did not
  make, in the setup of its preliminary search, `psiblast_aux_priv.cpp:90-93`; refused before the
  run's first byte until this version, line PS DECISIONS 35 Q5), a window past the PSSM (1.18), a
  gapped start past the subject's translation (1.48) and a traceback that starts past its frame's
  NULLB (1.44). LOSAT writes what NCBI has formatted before the failing step; what NCBI's crash
  loses from its output buffer is not imitated (NCBI's file lacks the end of the text that buffer
  held: 21 bytes of the `-in_pssm` prolog reach it in the probe of 2026-11-18; decision S1 (a) of
  line RG item X-short-bat, Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07). The refusal of a negative `CHUNK_SIZE` whose wrapped ranges leave gaps between
  them, which NCBI searches on (D13), comes after the same reports, where the search of the
  batch's first chunk with a valid context would read its subjects (after that chunk's setup
  checks and the reads of the database that fail before it). Before the run's
  first byte stays the refusal of a reported column of that NULLB (1.44), whose alignment NCBI
  reports.

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
  X-short-bat-2). The refusal is the error of the earliest such subject at every
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
  of the same read (the item above) comes after the same reports (version 1.97, line RG item
  X-short-bat): NCBI's TBLASTN searches and prints each query batch before it reads the next
  (`tblastn_app.cpp:275-304`, as TBLASTX's `tblastx_app.cpp:176-207`) and so has formatted the
  reports of the batches before when it dies (line RG X-short review 1: with `BATCH_SIZE=1`, NCBI
  wrote the first query's 43-byte row before its SIGSEGV, LOSAT nothing before version 1.97). With `-subject`
  and a `-db_gencode` other than 1, LOSAT translates in the selected code
  (PD-TBLASTX-LOCAL-SUBJECT-GENCODE) and refuses where those HSPs reach the re-evaluation. TBLASTX
  translates the query whole before the search, so a short query reaches no such read; but a
  short query against a short subject (the probes of exception 20's paragraph on them) makes NCBI
  die earlier, in its first linking call (`BLAST_LinkHsps` under `s_BlastSearchEngineCore`, gdb on ws,
  2026-11-17 and 2026-11-18): their big-gap cutoff is `INT_MIN + 1`, approved exception 20 (since
  version 1.99, a short query against one short subject). LOSAT links them with its sums in 64
  bits; against 12 bases or more it writes rows, and against fewer than 12 bases it refuses with
  this item's message at the re-evaluation that follows the linking.

- TBLASTN with `-db_soft_mask`, a gapped search of a chunk that saves an initial alignment starting
  before the chunk (version 1.89, 2026-11-23, line NP item M4b-T2, Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; narrowed in version 1.104, line NP item M4b-T3).
  Masks stored out of order on a frame of more than 5,000,000 residues let the scan find words
  before the chunk (the rule-4 item of version 1.89); the one-hit extension of such a word
  (`-window_size 0`), or a two-hit extension whose word lies wholly before the chunk, gives an
  initial alignment that starts there. No gapped alignment found before holds it: those start in
  the chunk (`blast_itree.c:822-831`). NCBI's gapped start of an alignment of at most 11 letters
  (`HSP_MAX_WINDOW`) is its middle, read from no letter (`blast_gapalign.c:3401-3404`), and the
  subject's start moves with it (`blast_gapalign.c:3942-3945`). Where that start lies in the
  chunk, NCBI's preliminary gapped alignment runs from it on the chunk's letters and its traceback
  on the frame: NCBI BLAST+ 2.17.0 gives the same rows every run, at 1 and 2 threads (probes of
  2026-11-24: alignments of 11 letters from the chunk's offsets -1, -4 and -5, `-window_size 0`,
  with and without `-in_pssm`), and LOSAT gives them (version 1.104, line NP item M4b-T3). Where
  the alignment is longer than 11 letters (its length read as `Uint4`, so also a negative one),
  NCBI reads the subject at the alignment's start taken as `Uint4`, about 4 GiB past the frame
  (`blast_gapalign.c:3405-3413`; gdb: `BlastGetStartForGappedAlignment`), and where the middle lies
  before the chunk, its right extension reads the subject from that start taken as `Uint4`
  (`blast_gapalign.c:4300-4333`; gdb: `Blast_SemiGappedAlign`): NCBI BLAST+ 2.17.0 ends with
  SIGSEGV every run (with and without `-in_pssm`, at 1 and 2 threads; the database `before` of
  `tests/tblastx_db_soft_lead.rs` and the probes of 2026-11-24). LOSAT refuses exactly these (exit
  1, "TBLASTN: database sequence <OID>, frame <frame>: an initial alignment of the chunk from
  residue <n> starts before the chunk (...)"; one producer,
  `search_gapped::x_gapped_start_before_chunk`; the manual's refusals page), on the initial
  alignments of NCBI's one-thread search: where NCBI's thread clears its table inside the
  sequence, the words before a later chunk are not extended and NCBI completes the search, and
  LOSAT gives that result at every thread count (the rule-4 item; PD-LOSAT-DATABASE-SEARCH, DB-3).
  Before version 1.104 LOSAT refused every such alignment. An initial alignment that starts in the
  chunk from a word before it (a seed at the chunk's offset -1 or -2) is not refused: its gapped
  start is the seed's `Uint4` moved as NCBI's (`blast_gapalign.c:3942-3945`), in the chunk. An
  ungapped search (`-ungapped`) keeps these alignments as NCBI does (rule 4).

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
  and length with the HSPs of the one out of the range; fix round 3, review C N4). RPSBLAST refuses before any report (its search of every query batch
  comes before the reports; with `-outfmt` 0 to 4 after the prolog, which it writes before the
  search since version 1.126); RPSTBLASTN after the reports of the query batches before. With several
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
  of TBLASTN's endless split of 1.18; RPSBLAST before any report (with `-outfmt` 0 to 4 after the prolog,
  version 1.126), RPSTBLASTN after the reports of the query batches before (NCBI writes the batches
  before, then fails). With at most one such
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
  exit status -11; probes p9_over of line RG D1a-4 and g9#1, `NSEQ 5`); RPSBLAST before any report
  (with `-outfmt` 0 to 4 after the prolog, version 1.126), RPSTBLASTN after the reports of the query
  batches before. Elsewhere `Lc` above the rows is NCBI's
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

- An alias file's `SEQIDLIST` over a version 4 database (version 1.85, 2026-11-22, line RG item
  D1c-1; rule 2, the rejection kept by the maintainer: Owner's answer 「Q-D1c-R2=A」, 2026-11-22
  22:4x JST), in a shape whose
  NCBI result is deterministic but has no known use: a token of the list that NCBI warns about and
  drops (a token of form feeds and vertical tabs alone: `WARNING:  <token> is not a valid seqid
  string.` when NCBI reads the list, at its constructor's scan or at the search's first read under
  `NSEQ` and `LENGTH`). LOSAT refuses (exit 1) where NCBI reads the list. A list that two volumes
  or more translate, with an id that a second simplification changes, is read in NCBI's order:
  the ids simplified again for each later volume and not sorted again (rows `db.il.bp.cv.*`).

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
    reads past. (A FASTA record of 4,294,967,296 letters or more, refused here until version 1.100,
    is read as NCBI's reader reads it, with its length modulo 2^32; see the packing item below.)
- `losat makeprofiledb -dbtype rps` (version 1.82, 2026-11-22, psiblast line item PS10a, Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07): where NCBI's makeprofiledb 2.17.0 crashes,
  aborts, never ends or reads memory that is not its data. LOSAT refuses each (exit 1, "... is not
  supported by LOSAT's MAKEPROFILEDB") where NCBI's run meets it; the rows are those of the line's
  oracle (`PR/ps10/oracle/cases.tsv`) and of its review.
  - `-max_smp_vol` 0, or a negative value that gives 0 volumes: NCBI divides by the number of
    volumes (`makeprofiledb.cpp:1336-1337`; SIGFPE, rows V05 and RN03).
  - `-in` naming a directory: NCBI's loop over the list's lines (`makeprofiledb.cpp:585`) never ends
    (row O16).
  - A `-taxid_map` file with a line of white space only: NCBI crashes (SIGSEGV, row RX06).
  - A profile stored by row (`byRow TRUE`) with more than 28 positions: NCBI's `.freq` file reads the
    frequency ratios past their matrix of 28 rows (rows S_byrow and S_fbyrow); a profile with
    frequency ratios only whose list is shorter than the profile: NCBI iterates its list past the end
    (row S_ffrs).
  - A frequency ratio whose value times 10^9 does not fit a 64-bit integer (row S_frhuge), a word
    score threshold that is not a number after the scaling, a `-scale` whose ceiling does not fit
    `int`: NCBI's conversions of them have no defined result.
  - A profile with frequency ratios only whose query has a residue that is none of the 28 NCBIstdaa
    codes (version 1.84, PS10c): NCBIeaa other than a capital letter, `*` or `-`, IUPACaa other than a
    letter of either case (NCBI's `CSeqConvert` tables give them 255), or an NCBIstdaa code of 28 or
    more. NCBI's engine counts the query's composition in an array of the 28 codes
    (`blast_stat.c:1962-1971,1984-2007`) and writes outside it; what follows depends on the code
    (PS10c probes, `PR/ps10c/probe*`: SIGSEGV for every such byte of 42 to 90, ten bytes outside and
    the codes 28, 29, 64, 128 and 255, malloc's abort on its corrupted heap for 31, a database for
    32, with or without a later error such as a negative frequency ratio). A query of the wrong
    length and a matrix that NCBI's engine does not know are NCBI's errors first (fixture rows
    QC02-QC05), and a lower-case IUPACaa query is built (row QC01).
  - A query in `iupacna` coding holding a byte `>= 0x80` in a pair position before any final odd
    overhang byte (version 1.105, PS10b; reachable in ASN.1 binary where `StringStore` keeps raw
    bytes without `#` replacement): `WriteDB_IupacnaToBinary` (`writedb_convert.cpp:590-604`)
    converts all input bytes through `x_ConvertIupacnaTo4na` (`sequtil_convert_imp.cpp:313-339`),
    whose pair loop indexes `sm_IupacnaTo4na` with a signed `char` and reads before the static table
    in `.rodata` (undefined behaviour; fixture inputs `l_QR08.lst`-`l_QR10.lst`, scored and
    frequency-ratio-only alike). Bytes `< 0x80` (row `B19`), the single odd overhang byte which
    casts to `Uint1` (row `B18`), and `mol na` (which stops at NCBI's molecule error first, row
    `B20`) are built or reported as in NCBI.
  - Binary definition lines (`ASN1_BlastDefLine`) holding a `Seq-id` element stated as `NULL`
    (version 1.105, PS10b; fixture input `l_QR12.lst` and rows `RB17`-`RB29`): carried to NCBI's
    exact failing point rather than refused up front. With an active taxid map (a non-empty
    `-taxid_map` or an earlier/current `org` descriptor without `-taxid`) `FixTaxId` stops before
    writing with `Attempt to access NULL pointer.` (exit 255, rows `RB19`, `RB20`, `RB27`); with
    `-index false` when the `NULL` element is the first `Seq-id` of the first definition line,
    replacing it with `gnl|BL_ORD_ID|<oid>` removes it and builds the database (exit 0, rows `RB24`,
    `RB25`); otherwise BER serialization in `WriteDb` posts `(802.5)` and fails at the exact BER
    stream byte offset (`byte N: cannot write unassigned member NULL element while writing container `,
    exit 255 in `-blastdb_version 4` and in version 5 when an earlier profile has already written an
    OID; rows `RB17`, `RB18`, `RB21`-`RB23`, `RB26`, `RB28`, `RB29`). Only in version 5 when no OID
    has been written yet (no `-taxid_map` or an empty map, including `-taxid`), where
    `~CWriteDB_TaxID` throws `No tax info for any oid` out of its destructor and NCBI aborts
    (exit 134), or when a later input error unwinds through `~CWriteDB_Impl`'s `Close` (v4 or v5,
    exit 134), does LOSAT refuse via the writer-destructor-abort refusal below (fixture inputs
    `l_QR12.lst`, `l_QR15.lst`). Damaged binary definition lines, replaced visible characters (803.7), `NULL`
    `Blast-def-line` elements (`Attempt to access NULL pointer.`), and huge lengths (`Unknown
    exception`) match NCBI through the shared table-driven BER reader (fixture rows `RB04`, `RB04m`,
    `RB04b`, `RB08`-`RB16`).
  - `:OtherType` pointers and host-dependent 1 GiB..2^47 string or octet lengths in `makeprofiledb`
    profiles (both ASN.1 text and `-binary`, read through `Names::PROFILE` as `makeblastdb` uses
    `Names::MASK`) and binary definition lines (version 1.105, PS10b; fixture inputs `l_QR11.lst`,
    `l_QR13.lst`, `l_QR14.lst`): because `makeprofiledb` prints NCBI's `CSerialException` text
    directly, whether `:OtherType` prints `class not found` or `incompatible member type` depends
    on which ASN.1 classes the process has registered so far (for example, on text `:Seq-interval`,
    NCBI 2.17.0 prints `class not found: Seq-interval`, whereas LOSAT before 1.105 (`Names::PSSM`)
    erroneously printed `incompatible member type`), so both text and binary `:OtherType` pointers
    are refused under one rule.
  - A query in a nucleotide coding with fewer bytes than NCBI's writer reads for its length (version
    1.86, PS10d): `ncbi2na` of fewer than `ceil(length / 4)` bytes, `ncbi4na` of fewer than
    `ceil(length / 2)`, `iupacna` whose NCBI4na form and its terminating zero hold fewer than
    `ceil(length / 2)` bytes (four letters or more short of an even length, three of an odd one).
    NCBI's `x_CookSequence` converts the kept sequence when it is published
    (`writedb_convert.cpp:444-528,569-610`; at the next query, at the close, or at the close after an
    error, `writedb_impl.cpp:100-108,158-180,324-331`), and its checks of the size are debug-only: it
    reads past the data (fixture inputs `l_QR01.lst`-`l_QR05.lst`: `.psq` bytes that differ between
    two runs, SIGSEGV for empty `ncbi2na`). The data of
    enough bytes (more than the length needs, unused bits set, iupacna up to three letters short of an
    even length and two of an odd one, whose read ends at the string's zero or before it) is NCBI's
    (fixture rows QD01-QD40).
  - A profile with frequency ratios only whose query's `Seq-data` is a `gap` (version 1.86, PS10d):
    NCBI's conversion of the query to NCBIstdaa looks its coding up in a table of the `CSeqUtil`
    codings that ends before the gap's number (`seqport_util.cpp:1947-1959,2032-2035`) and reads
    past it. The other codings that the conversion meets are NCBI's (fixture rows QD41-QD46:
    `Unknown conversion.` for NCBI8na, NCBI8aa's bytes as the query's codes, the unknown exception
    of NCBIpna and NCBIpaa, which the table has no coding for; then the writer's error for the
    coding it cannot store, posted or reported where NCBI's is).
  - A profile with frequency ratios only whose matrix has no frequency ratios, or whose query has no
    ungapped Karlin block (NCBI's engine reads `kbp_psi[0]`, a null pointer).
  - `-scale 1` on frequency ratios whose scaling search does not converge: NCBI's search never ends
    (review rows RF_BLOSUM45_s1, RF_BLOSUM50_s1, RF_BLOSUM80_s1, RF_tiny1).
  - A query that has only GI ids between queries with other ids, in a version 5 database with
    `-index true`: NCBI's LMDB writer aborts (rows I01g and I01x; makeblastdb's refusal of the same
    shape).
  - A query without Seq-ids and without descriptors, in a version 5 database with `-index true`:
    NCBI's taxonomy writer aborts as it closes after the error (row RI03); with `-index false`
    (version 4 or 5), and for a query whose first binary definition line has no Seq-id, NCBI writes the
    `gnl|BL_ORD_ID` id into the first element of an empty list (`writedb_impl.cpp:865-874`; exit 0, rows
    RI03x, RI04 and RB02x). In a version 4 database with `-index true` NCBI's duplicate check stops
    with `Duplicate seq_ids are found:` (exit 255, rows RI02 and RI02b), which LOSAT writes as NCBI does
    (rule 2).
  - An error after which the database's writer fails again as it closes: NCBI's destructor closes the
    writer during the unwinding and the second exception ends the process (`makeprofiledb.cpp:334-377`;
    `terminate`, exit 134). Rows RZ07 (an empty `-out` in version 5 with an input file after the first
    missing), RB03 and RB03b (a binary definition line set without lines in version 5) and RI03.
  - A profile stored by row without scores whose query length is 4,294,967,295 with `numColumns -1`:
    the `int` length is -1, NCBI's matrix has no row (`makeprofiledb.cpp:144-156`) and its extra row is
    written at index -1 (`makeprofiledb.cpp:1106-1167`; SIGSEGV, row RL05). The other shapes of a length
    past `int` are written as NCBI writes them (rule 2: rows RL01 to RL04, RL06 and RL07).
  - A volume whose offsets or positions pass 2^31, or whose word table passes `Int4` bytes: NCBI's
    counters of them are `Int4`, and the files hold 32-bit fields, so no valid file can be defined.
  - With `-dbtype delta` (version 1.113, 2026-11-26, psiblast line item PS14, Owner-delegated,
    standing instruction 2026-09-29, restated 2026-10-07): profiles stored by row (`byRow TRUE`) or
    with `numRows < 28` and non-empty `weightedResFreqsPerPos` (`makeprofiledb.cpp:1610-1647,
    1763-1800`, where NCBI copies into an empty `std::list::begin()`, overwriting the list header
    node; rows `QR16`, `QR17`); an empty `numIndeptObsr` list (`makeprofiledb.cpp:1659, 1803`,
    where NCBI dereferences `end()`; row `QR18`); a `weightedResFreqsPerPos` list whose length is
    not a multiple of 28 when `len / 28 == numIndeptObsr.len()` and every full column of 28 sums
    within 0.0001 of 1.0 (`makeprofiledb.cpp:1720-1740`, where `x_ValidateCd` increments past
    `end()`; row `QR19`; when an earlier full column fails the sum check, NCBI returns false
    deterministically before reaching the partial tail, reproduced in rows `D33`-`D35`);
    `weightedResFreqsPerPos` or `numIndeptObsr` values whose scaled conversion `(val * 1000)` or
    `((val + 1) * 1000)` to `Uint4` is out of range or `NaN` (`makeprofiledb.cpp:1674, 1688`; rows
    `QR20`-`QR22`); and a missing or unwritable temporary directory (`TMPDIR`) when `x_UpdateDelta`
    creates its two temporary files (`makeprofiledb.cpp:1572-1575`, `ncbifile.cpp:6580-6599`, where
    NCBI throws `mkstemp() failed for '...'` with the process and thread IDs in the message, which
    changes from run to run). Delta 32-bit offset wrapping (`CurrFreqOffset += FreqBuff.size()`,
    `CurrObsrOffset += ObsrBuff.size()`, `Uint4 num++`) is deterministic in NCBI and reproduced
    with wrapping arithmetic.
- The readers of GI-indexed masks, a damaged `.nog` (version 1.94, 2026-11-24, line NP item
  M4c-2, review S4 and N8; the rejection of a `MASKLIST` name that is not UTF-8 withdrawn in 1.103; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07; the outcomes of a `.nog` cut short that NCBI
  reaches and LOSAT reproduces, rule 2, are in the list above). An empty `.nog` is mapped as a
  null pointer: NCBI's blastdbcmd and blastn crash with SIGSEGV where a sequence's GI is read
  (`-outfmt %m`, `-db_soft_mask`; 5 of 5 and 3 of 3 runs; `-info` reads no GI and completes).
  LOSAT refuses it there (exit 1, "(NCBI crashes on it)"). An entry past the last 4 KiB page of
  the file's mapping reads whatever memory follows it: in the probe a `.nog` cut to 100 bytes for
  1,100 sequences, NCBI's blastdbcmd completed, the 84 sequences past the page taking GIs that
  matched no mask, which depends on its memory map; LOSAT refuses it ("depend on its memory map").
  Rejected: reproducing the bytes past the page (a model of NCBI's memory map). Tests:
  `LOSAT/src/seqdb/gi_mask_tests.rs` (`the_nog_reads_as_ncbi_s_mapping_of_it`).

- BLASTP, BLASTX, TBLASTN and TBLASTX `-db_hard_mask` on a database whose masks overlap, nest or
  are stored out of order (version 1.95, 2026-11-22, line NP item M4b-3b, Owner-delegated, standing instruction
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
  (approved exception 29). Fixture rows `db.pmask.bp.hard.ov_*` and
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
  - `HardMaskShape::LastRangeOverlap` (before the search) and
    `hard_chunk_scan::last_range_overlap_merge_refusal` (during the search, for short `OneWord`
    chunks) (version 1.122, 2026-11-27, line NP item M4b-3d, Owner-delegated, rule 3): version 1.122
    ported NCBI's `Blast_HSPListsMerge` across overlapping subject chunks for BLASTP and BLASTX
    `-db_hard_mask` (ranges longer than 5,000,000 letters and out-of-order masks whose parts overlap
    when NCBI reads the merge overlap in bounds, replacing version 1.95's unported `PartsOutOfOrder`
    and in-bounds `short_chunk_merge_refusal`). On the final chunk of the final hard range
    (`overlap_past_list`), NCBI reads `backup.hard_ranges[num_hard_ranges].left` one slot past its
    list (`blast_engine.c:579-580`), which holds memory left by a prior sequence on the same thread.
    Evidence (`po_last`, masks `[250, 280)` then `[100, 130)`, `blastp -ungapped`): with a preceding
    sequence whose `hard_ranges[3].left == 130`, official NCBI BLAST+ 2.17.0 at 1 thread writes 50
    rows (`0e5f2d460c84`), while the control database and 4-thread run write 33 rows
    (`f4b4cf155cea`; `POISON_DIFF=True`, `T1_T4_DIFF=True`). Matching BLASTN's decision B-d
    (`ShortLastMask`, `PD-LOSAT-DATABASE-SEARCH.md`), LOSAT refuses before the search where the
    final chunk scans a word and starts before the end of an earlier chunk from another hard range
    (`LastRangeOverlap`, fixture rows `db.pmask.{bp,bx}.hard.rej.pe_unsorted.*`), and during the
    search where an HSP from an earlier hard range (such as a short `OneWord` chunk) meets an HSP of
    the final chunk across that out-of-bounds overlap (`last_range_overlap_merge_refusal`); merges
    within a final hard range longer than 5,000,000 letters keep `DBSEQ_CHUNK_OVERLAP = 100`.
    Out-of-order merges in middle ranges and long-range merges are searched byte-for-byte
    (`db.pmask.{bp,bx}.hard.pe_unsorted_mid.*`, `LOSAT/tests/protein_db_hard_long.rs`).
    `hard_chunk_scan::chunk_merge_undefined_refusal` refuses a chunk-overlap merge whose query
    lengths sum to 0 or whose scaled score leaves `Int4` (version 1.21 extended in version 1.122 to
    BLASTP and BLASTX `-db_hard_mask`, TBLASTN subject-chunk merges, and BLASTN subject/query-chunk
    merges). Counts of `LastRangeOverlap`: TAIR10 0, the UniRef50 sample 0.
  - `HardMaskShape::Inverted`, before the search: a mask that ends before it starts, which makeblastdb
    never writes.
  TBLASTN and TBLASTX (version 1.106, line NP item M4b-3c part 2) search masks that overlap or nest
  in the order of their start as NCBI searches them: each strand maps the hard ranges once
  (`blast_engine.c:805-830`), every frame's chunk moves the one table of diagonals by its `Int4`
  length and the window, a chunk shorter than the word gets its one word in the two-hit scan, and
  all six frames of such a sequence run in order on one table, at NCBI's one-thread offset where
  the sequence meets entries that no sequence wrote. They refuse masks out of the order of their
  start (`Unsorted`): NCBI's traceback masks the slices of its partial fetch of a long nucleotide
  sequence (`seqsrc_seqdb.cpp:233-242`) by a scan that stops at the first mask that starts past the
  slice (`seqdbvol.cpp:1518-1538`), exact for masks in start order and not otherwise. They refuse a
  sequence with overlapping masks one of whose frames NCBI searches whole (the no-chunking branch
  of exceptions 21 and 27), `short_chunk_tree_refusal` (evidence `g3`: NCBI tblastn's gapped search
  never ends, memory to the 4 GB cap within 4 s) and nothing else (approved exception 29). BLASTN is
  unchanged.
- RPSBLAST and RPSTBLASTN on a profile whose `.psq` letters (its `.pin` length `L_pin`) differ from
  its `.rps` columns `L_rps` (version 1.102, 2026-11-23, psiblast line item PS10e, line PS
  `DECISIONS.md` 97, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07).
  makeprofiledb writes such a profile for a scored `.smp` whose `seq-data` is nucleotide-coded or has
  more or fewer letters than its length (it reads no letter of a scored profile and checks only the
  length). NCBI's preliminary stage reads the profile's `.rps` rows only; its traceback
  (`s_RPSComputeTraceback`, `blast_traceback.c:1046-1316`) takes the profile's `.psq` letters as its
  query and `L_pin` as the profile's length. LOSAT searches such a profile as NCBI does, with `L_pin`
  in the traceback (the rows from the profile's offset, the composition, the words, the redo's
  matrix and query, the e-value length; fixture rows `db.rps.odd.*` and `db.rpstbn.odd.*`), and
  refuses the query batch (exit 1, "a hit of the profile of OID N, whose consensus has
  L_pin letters in .psq for L_rps columns of scores in .rps (makeprofiledb writes it for a scored
  profile whose seq-data is nucleotide-coded or has more or fewer letters than its length), where
  NCBI BLAST+ 2.17.0's traceback READ, is not supported by LOSAT's RPSBLAST", one text for both
  programs; `LOSAT/src/algorithm/rps_search/profile.rs`) at the first read where NCBI's traceback
  leaves the profile's data, in NCBI's order (survey `PR/ps10e/SURVEY_ps10e.md` of the line's
  records, its probes named here):
  - A, the rows from the profile's offset: `L_pin` `.freq` row pointers with composition-based
    statistics, `L_pin` `.rps` row pointers for RPSBLAST's rescale without them, of NCBI's `R + 1`
    (`blast_traceback.c:968-982`; entry `R`, the row after the rows, read in the file's mapping): the
    offset plus `L_pin` above `R + 1`, or a column of row `R` past the file's last page that NCBI
    reads (the rescale and the `.freq` rows read all of its columns; RPSTBLASTN without
    composition-based statistics reads the mapped row `R` one letter at a time, so only a letter its
    alignment reads there counts) (probes A_P23x3L_q3, A_P23x3L_q3_c0 and A_P23x3L_nq3, SIGSEGV;
    rows `db.rps.rej.odd.p23x3l.q3`, `.p23x3l.q3c0`, `db.rpstbn.rej.odd.p23x3l.nq3`).
    READ is "reads past its rows".
  - B, the composition with composition-based statistics (`Blast_ReadAaComposition`,
    `composition_adjustment.c:1162-1199`): a letter of 28 or more, such as the packed bytes of a
    nucleotide-coded profile (B_SC_q2, SIGSEGV; rows `db.rps.rej.odd.sc.q2`, `.sc.q1q2`, `.sc2v.q2`,
    `.hl.q12.m1`, `db.rpstbn.rej.odd.sc.nq2`). READ is "reads letter B past its composition tables".
  - C, the near-identical test of the redo (`s_TestNearIdentical`, `blast_kappa.c:1259-1340`, on the
    copy of `s_SequenceGetRange`, `:1671-1708`): a read past the copy's `L_pin + 2` bytes or past its
    `L_pin - 7` words, or, below 8 letters, of the null word table at an offset other than 0
    (`:1161-1166`). T_AA1_q1, Q_P16_q1, Q_P10_q1 and Q_P10_nq1 exit 0 with invalid reads under
    valgrind; Q_P6_q1 ends with SIGSEGV (rows `db.rps.rej.odd.aa1.q1`, `.p16.q1`, `.p10.q1`, `.p6.q1`,
    `db.rpstbn.rej.odd.p10.nq1`). READ is "reads past its copy of the consensus".
  - D, an alignment that starts at or past `L_pin` in NCBI's heap matrix of `L_pin` rows: the redo's
    (`s_RedoOneAlignment`, `blast_kappa.c:1924-1946`), and RPSBLAST's rescaled matrix without
    composition-based statistics (T_AA1_q1_c0, H_HL_q12_m6c0: SIGSEGV; rows
    `db.rps.rej.odd.aa1.q1c0`, `.hl.q12.m6c0`). READ is "starts an alignment past the consensus".
  - E, the identities and positives of an alignment traced without composition-based statistics
    (`blast_traceback.c:589`, `blast_hits.c:796-812`): a profile letter of 28 or more in a
    substituted column reads past the 28 rows of the score matrix (E_NA1_qi11_c0 and
    E_NA1_nqi11_c0 of the line's `PR/ps10e/probe2`: c1 nucleotide-coded, 6 letters 143 18 159 6 210
    210 for 22 columns, and c1's first 11 letters or their back-translation; the alignment starts at
    column 5, below `L_pin`; SIGSEGV natively and under valgrind, invalid read in
    `s_Blast_HSPGetNumIdentitiesAndPositives`; rows `db.rps.rej.odd.na1.qi11c0`,
    `db.rpstbn.rej.odd.na1.nqi11c0`). Until those rows no fixture row reached E: in N_SC_nq2_c0
    NCBI's run ends at E after a start past `L_pin`, where LOSAT refuses first (the rule-2 item
    below), and with 12 letters or more the start is at or past column 6 (`L_pin`), read D
    (E_NA1_qi12_c0 to E_NA1_qi14_c0). READ is "reads letter B past its score matrix".

  A profile whose `L_pin` is `L_rps` never meets them, and a hit that NCBI's preliminary e-value or
  hit list leaves out of the traceback meets none (rows `db.rps.odd.sc.q2.e14`,
  `db.rps.odd.hl.q12.m5c0`). The refusal at the volume's open that compared each profile's rows with
  its `.pin` length (from series PR item PS11b) is narrowed to its part that only a hand-edited
  `.rps` meets, offsets that do not ascend (a profile without its sentinel row, which makeprofiledb
  writes for every profile: "a profile database whose .rps profile N has no rows, not even its
  sentinel row (its offsets A and B do not ascend) is not supported"); LOSAT refused such a file
  before too. NCBI crashes on one when a query hits its profiles (N2_SB0eq_qcdd, N2_SB0desc_qcdd and
  their `-comp_based_stats 0` twins of `PR/ps10e/probe2`: SB0 with profile 1's second offset equal
  to its first or 3 below it, SIGSEGV natively and under valgrind, invalid reads in
  `s_ExtendRight` and `ALIGN_EX`). The preliminary tree's refusal of a concatenated length below the
  rows (line RG item D1a-4) names this further cause: profiles whose `.pin` lengths are below their
  `.rps` rows (N_AA1s_nqt2_c0, NCBI exit 4; row `db.rpstbn.rej.odd.aa1s.nqt2c0`).

- A nucleotide FASTA record that NCBI's `CFastaReader` packs into a segment of 4294967293 to
  4294967295 bases of NCBI2na, or of 4294967295 bases of NCBI4na, in `losat makeblastdb -dbtype
  nucl`, `dustmasker`, `windowmasker`, `makembindex` and `convert2blastmask` (version 1.100,
  2026-11-24, line RG item KC-11d, Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07): `AssembleSeq` packs the record's pieces with `CDelta_ext::AddAndSplit`
  (`fasta.cpp:1456-1458,1504`), whose `CDelta_ext_PackTarget::NewSegment` sizes a segment as
  `(length + 3) / 4` or `(length + 1) / 2` bytes in a `TSeqPos` (`Delta_ext.cpp:163-164`); that
  wraps to 0 and the conversion writes past it. Every program ends with SIGSEGV (exit status 139;
  probes `remaining-gaps-20261112/K/kc/p11d/p2` s3n and s7n, through the `N` of `>?` gap lines,
  and `p3` b4n and b5n, 4294967295 letters of `ACGT` and of `ACGTN`, in dustmasker and
  makeblastdb). LOSAT refuses (exit 1, "... is not supported by LOSAT's <program>") at that
  record, after the substitution's `std::out_of_range`, which NCBI throws first. A protein is
  stored unpacked and meets no such segment.
- `losat makembindex -old_style_index true` of a FASTA record of 4,294,914,101 letters or more,
  counting the letters of its `>?` gap lines and modulo 2^32 (same version, line and delegation):
  `CSubjectMap_Factory_TBase::AddSequenceChunk` takes a chunk's end as
  `(TSeqPos)(chunk_start + chunk_size_)` (`dbindex_factory.cpp:724-725`), which wraps for a chunk
  that starts 4,289,967,296 letters or more into the sequence (the chunks are 5,000,000 letters,
  4,999,900 apart); the chunk's one valid interval is then 5,000,000 letters from its start, past
  the sequence's end, and the index reads them (probes `p3` l93n and l94n, 4294967293 and 4294967294
  letters: SIGSEGV after `creating ....00.idx...`, no file left). LOSAT refuses (exit 1) at that
  sequence.
- A version 5 alias `SEQIDLIST` combined with a search `-seqidlist` or `-negative_seqidlist`, when
  the alias `SEQIDLIST` or the search's `-seqidlist` holds two differently spelled Seq-ids that
  NCBI's `s_CompareSeqId` treats as equivalent, or two unparseable Seq-ids whose `s_CompareSeqId`
  exceptions differ, or a Seq-id that posts a warning alongside an unparseable Seq-id (version
  1.109, 2026-11-25, line RG item D1c-2, Owner-delegated, standing instruction 2026-09-29,
  restated 2026-10-07): `s_ProcessSeqIdFilters` (`seqdboidlist.cpp:578-649`) sorts the search's
  `-seqidlist` (`:592`) and each alias `SEQIDLIST` (`:611`) with `s_CompareSeqId` (`:565-576`),
  which parses both ids inside the comparator (`CSeq_id` with `fParse_AnyRaw | fParse_ValidLocal`),
  returns `false` whenever `CSeq_id::Match` succeeds, and falls back to raw byte order otherwise.
  Because `CSeq_id::Match` is not transitive (for example, `NP_000103` matches both `NP_000103.1`
  and `NP_000103.2`, while `NP_000103.1` and `NP_000103.2` do not match; similarly `LOC7` and
  `loc7` bracket unrelated ids in ASCII order), `s_CompareSeqId` is not a strict weak ordering.
  NCBI's `std::sort` therefore invokes undefined ordering behaviour governed by its C++ library's
  sort steps, which changes which spelling `std::set_intersection` (`:615`) or
  `std::set_difference` (`:626`) keeps for case-sensitive LMDB `AccessionsToOids` (observed
  against NCBI BLAST+ 2.17.0 on database `g` with `NP_000103.2` / `NP_000103` / `NP_000103.1`
  against `-seqidlist NP_000103` or `NP_000103.1`: NCBI matches `NP_000103.1`, whereas another
  sort leaves `NP_000103.2` or nothing). For the same reason, when a list of two or more distinct
  entries contains unparseable ids with different exception messages (observed against NCBI BLAST+
  2.17.0 on `sp||` / `tr||` / `NP_000103.1`, where `libstdc++` compares entry 2 before entry 1 and
  throws `tr||`'s `release unreviewed` instead of `sp||`'s `release reviewed`, and on a 0-entry
  `.bsl` with spaces in its date whose header splits into `Malformatted ID ` and
  `Malformatted ID PM`) or a warning-posting id alongside an unparseable id, which exception or
  warning NCBI prints depends on the C++ library's sort steps. Copying a C++ library's sort
  implementation is barred by `AGENTS.md` rule 4. LOSAT parses all unique entries of each sorted
  list up front (`seqdb::alias::validate_compare_seqids`) and, in order: (1) refuses if any entry
  posts a Seq-id warning, (2) refuses if unparseable entries raise two or more distinct exception
  texts, (3) raises NCBI's `CSeqIdException` when all unparseable entries raise the same exception
  text, and (4) refuses (exit status 1) if any two distinct spellings match under `CSeq_id::Match`;
  lists without such pairs have `s_CompareSeqId` equal to raw bytewise `<` and sort identically in
  every sort. A search `-negative_seqidlist` is sorted by raw `<` in NCBI (`seqdboidlist.cpp:597`),
  so equivalent pairs inside `-negative_seqidlist` alone are not refused (equivalent pairs inside
  the alias `SEQIDLIST` remain refused, since `:611` still sorts it with `s_CompareSeqId`), and an
  alias `SEQIDLIST` without a search seqid list is not sorted and is not refused. Fixture rows
  `db.il.bp.v5.sq.amb*.rej.6`, `db.il.bp.v5.sq.uamb.rej.6`, `db.il.bp.v5.sq.bad*.rej.6`,
  `db.il.bp.v5.sqbe.sl.rej.6` (refused), `db.il.bp.v5.sq.slunp.6`, `db.il.bp.v5.sqbe1.sl.6`,
  `db.il.bp.v5.sq.uns.*`, `db.il.bp.v5.sq.amb1.plain.6`, and `db.il.bp.v5.sq.nslamb.6` of
  `LOSAT/tests/db_search_manifest/blastp.il.tsv`; unit test
  `v5_rejects_non_strict_weak_equivalent_seqid_pairs` in `LOSAT/src/seqdb/alias_tests.rs`.
- `losat makeblastdb -input_type asn1_txt` and `-input_type asn1_bin` (version 1.112, 2026-11-25,
  line RG item D3c-4, Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07):
  (a) a nucleotide `Bioseq` whose `iupacna`, `ncbi2na` or `ncbi4na` `Seq-data` holds fewer bytes or
  letters than `Seq-inst.length` asks for (including a `-hash_index` `CSeqVector` read of 1 to 3
  packed tail bases past the end of `ncbi2na`/`ncbi4na` data), where `WriteDB_Ncbi4naToBinary`,
  `WriteDB_Ncbi2naToBinary` and `copy_2bit_any`/`copy_4bit_any` read past the buffer (class of
  1.86; fixture `llong.asn`); (b) a list element stated as `NULL` (`05 00` in BER) in
  `Bioseq-set.seq-set`, `Bioseq.id`, `Seq-descr`, `Delta-ext`, `Blast-def-line-set` or
  `Blast-def-line.seqid`, or an empty `Delta-ext` list, where NCBI dereferences a null pointer or
  `dext.Get().front()` (`build_db.cpp:309`) and crashes (class of 1.45 and 1.52); and (c) a
  punctuation-only Bioseq title that `CDeflineGenerator::x_CleanAndCompress` reads past the end of
  and crashes on (class of 1.2 and 1.37). LOSAT refuses each with exit status 1 and erases the
  partial build (`erase_after_refusal`).

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

    **A short query against one short subject (version 1.99, line RG item X-short-tbx-link;
    the Owner's answer of 2026-11-24 to question Q-XT1, 「Q-XT1=A」).** The same cutoff comes from another
    input: NCBI's expected length (`BLAST_Nint(log(K * q * s) / H)`, with the average query length
    `q` and the subject's length `s` in protein letters) is 0 when `K q s` is near 1, and negative
    below about `e^{-H/2}`; the subject's length less the expected length is then the subject's or
    longer, and where the database's length in protein letters is not above it
    (`blast_parameters.c:1046`) `y_variable` is the logarithm of the subject's length over it. At
    an expected length of 0 that takes one subject (other subjects of fewer than 3 bases in all
    count as none, as above): `y_variable` and `x_variable` are 0 and the logarithm -inf, as above.
    Below 0, `y_variable` and `x_variable` are negative and the logarithm is NaN, which NCBI's
    x86-64 conversion also takes to `INT_MIN` (`blast_parameters.c:1036-1072`); there other
    subjects whose length in all, in protein letters, is not above minus the expected length keep
    the database short enough (the several-subject case below). On aarch64, whose conversion
    instruction takes a NaN to 0, NCBI's build would link with a cutoff of 1 there (not run; LOSAT
    follows the x86-64 build, `ncbi_int4_from_double`). In both the big-gap cutoff is
    `INT_MIN + 1`, the linking's
    `maxscore` (`-cutoff[1]`) is `INT_MAX`, no sum that overflows `Int4` reaches it, and NCBI reads
    the chain it never found (`link_hsps.c:694,868,931-932`; gdb on 2.17.0: a read at address 0x2a
    in `BLAST_LinkHsps` under `s_BlastSearchEngineCore`). This takes short inputs. In our probes
    with one subject (TBLASTX and BLASTX with `-window_size 0`): TBLASTX queries of 9 to 11 bases
    against subjects of 11 to 14 bases and of 12 to 14 bases against 11 bases; at `-word_size 2
    -threshold 1`, 6 to 8 bases against 8 to 20 bases, 9 to 11 against 8 to 14, 12 to 14 against 8
    to 11 and 15 or 16 against 8; TBLASTN `-ungapped` `WWWW` against 12 to 14 bases at `-word_size
    2 -threshold 1`; BLASTX `-ungapped` at `-word_size 2 -threshold 1`, queries of 6 to 8 bases
    against 2 to 6 residues, 9 to 11 bases against 2 to 4, 12 to 14 against 2 or 3 and 15 or 16
    against 2 (with the two-hit word finder no HSP of these reached BLASTX's linking). NCBI dies in
    every probed run whose subject list reaches the linking with an HSP (89 of 418 TBLASTX runs of
    `TGG` repeats, 42 of 600 runs of random short sequences, the TBLASTN one and 32 of 132 BLASTX
    runs). LOSAT links with the same cutoff and its sums in 64 bits, as above; before this version
    it took Rust's conversion of the NaN, 0, a cutoff of 1. NCBI's nearest completed run is the
    subject with a second one that has no hits (12 N, or 300 N where 12 leave the database's length
    in protein letters under the subject's length less the expected length; 4 X for BLASTX), which
    makes the database longer and the cutoff finite, and `-dbsize` of the subject's length. In
    every probed run it wrote LOSAT's one-subject bytes: all 45 TBLASTX runs above with subjects of
    12 bases or more, all 42 random ones, the TBLASTN one and all 32 BLASTX ones. That is not a
    property of the inputs: with a finite cutoff NCBI's linking can choose other chains where a
    frame offers competing ones, as for the batches averaging 0; the frames of these probes, of a
    few letters, offer almost none. The cases `link.q11_s12`, `link.q8_s12_ws2` and
    `link.rq76_rs76_ws2` of `tblastx_regression_fixtures.py` (NCBI's arguments in its
    `ORACLE_ARGV`), `tests/tblastn_ungapped_short_link.rs` and
    `tests/blastx_ungapped_short_link.rs` keep them; the unit tests of
    `core::link_hsps::calculate_link_hsp_cutoffs` keep the cutoff of both routes. With a TBLASTX
    subject of fewer than 12 bases the next step of such a search is the re-evaluation that the
    TBLASTX item of version 1.51 refuses (rule 3), so LOSAT refuses there, after the linking.

    *Several subjects (found after the decision of Q-XT1; the Owner's answer of 2026-11-24 to
    question Q-XT2 of line RG item X-short-tbx-link, 「A」, its option A: exception 20 covers every
    ungapped linking whose big-gap cutoff is NCBI's `INT_MIN + 1`, whatever the number of
    subjects).* On the NaN route other subjects whose length in all is not
    above minus the expected length leave the cutoff at `INT_MIN + 1`, and NCBI dies the same way:
    10 of the 42 random TBLASTX runs above died with the subject and a second subject of 12 N.
    LOSAT computes the same cutoff and writes rows there (as it did before this version, with its
    cutoff of 1).

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

## Approved exception 24 (BLASTN `-use_index true` with `-db_hard_mask`, the seeds in their hard ranges, version 1.78)

Version 1.78, line RG item D5b-5. The Owner decided on 2026-11-22 22:3x JST: the answer 「Q5=C」
(the whole answer 「Q1=B Q2=B Q3=A Q4=A Q5=C Q6=B Q7=B Q8=B Q9=A Q10=A Q11=A Q12=A Q13=OK」) to
question Q5 of that day's batch, line RG's question Q-HARDIDX, whose option C read 「種を正しい範囲に置く
（NCBI が出さない正しい結果、どの thread 数でも同じ）」 (place the seeds in their right ranges: a correct
result that NCBI does not give, the same at every thread count). Not taken: option A, NCBI's
one-thread result at every thread count through a pass over the sequences in OID order that
decides which list each range takes (the line's recommendation before the answer), and option B,
the refusal (`-use_index true with a MegaBLAST database index (...) and -db_hard_mask`) that LOSAT
made from item D5b-2 (PR #1472) until this version. Owner-delegated (standing instruction
2026-09-29, restated 2026-10-07; P2, P3): what the right range of a seed means, below (in the
index's search a masked letter matches no query letter, and two hits of a diagonal with a masked
letter between them do not pair). Not taken: NCBI's lists cut at the masked letters after the
search (the item's first version), under which a masked letter could still qualify a two-hit
seed (6 of 400 random databases). Also delegated: the masks of sequences that an OID list leaves
out of the search are not read.

24. **BLASTN `-use_index true` with `-db_hard_mask`.** NCBI's indexed word finder
    (`MB_IndexedWordFinder`, `na_ungapped.c:1697-1782`) takes the seeds of `GetResults` of the
    sequence and `subject->chunk` (`blast_dbindex.cpp:869-875`): the list of the index's chunk
    whose number is the number of the chunk being searched. The index numbers each sequence's
    chunks over its whole range (5,000,000 letters that overlap by 100), while with
    `-db_hard_mask` each hard range starts chunks of its own (`s_GetNextSubjectChunk`,
    `blast_engine.c:221-262`). So NCBI extends the seeds of a masked sequence at the offsets of
    another range's chunk, and a sequence with more hard ranges than index chunks takes the lists
    of the sequences after it, which the word finder moves out, so that what each sequence finds
    depends on which thread reaches it first: in the survey's probe g_f (dustmasker's masks,
    `-db_hard_mask 11`) NCBI wrote 5 rows at `-num_threads 1`, 7 or 8 rows from run to run at 4
    threads, where its search without the index wrote 12. LOSAT searches the index with each letter of the
    `-db_hard_mask` masks matching no query letter (`dbindex::search`): a root whose Nmer holds a
    masked letter is no root, the extensions (`ExtendLeft`, `ExtendRight`) stop before a masked
    letter, and the two-hit search takes no earlier hit of a diagonal with a masked letter between
    it and the seed. (NCBI's search reads the index's letters, masked or not, and its two-hit
    search can qualify a seed by masked letters.) So each seed lies in one unmasked stretch, and
    no masked letter decides one. Each seed then goes to the hard range that holds its left end,
    into the last chunk of the range that starts at or before it (the index's own rule where
    chunks overlap), in the order of the lists and of each list (`mb_index_seeds::HardSeeds`).
    The refusals of the search without the masks still apply: a chunk that is a later sequence's
    in the index, and a seed at or past its chunk's last letter. The word finder, the ungapped
    and gapped stages and the traceback are those of every hard chunk. A local id of the index
    (its group of sequences) without masked letters is searched as without `-db_hard_mask`. A
    sequence without masked letters keeps NCBI's seeds, unless one of its exact runs reaches past
    its end into masked letters of the next sequence of its local id (NCBI's index search reads on
    there). The result depends neither on `-num_threads` nor on the masked letters. The masks of a
    sequence that an OID list leaves out of the search are not read: its letters count as
    unmasked, and only a run past the end of the sequence before it reads them. The refusals of
    `HardMaskShape` and approved exception 13 (a mask that runs to the end) apply as they do
    without the index. Evidence: the rows `plain NAME` of `LOSAT/tests/blastn_mb_index/cases.tsv` (`hm_i*`,
    `sm_i_hard`; database hm of `gen_d5b5.py`: masks inside, at the start and at the end of
    sequences and hard ranges of 1 and 2 letters; nmer 8, 9 with stride 3, and 11; word sizes 16
    and 28; `-dust no`; 1 and 4 threads; an old-style index; database hmr of `gen_d5b5_r.py`, hm with every masked
    letter replaced at random), whose output is NCBI's one-hit search without the index with the
    same `-db_hard_mask` (hm: 32 rows, where NCBI's search with the index writes 5). The two-hit
    search (`-window_size`) of NCBI's index pairs other hits than its scan even without masks, so
    no NCBI search gives its result: the rows `expect hm_i_win` (`-window_size 20` at 1 and 4
    threads, and on hmr) give LOSAT's 31 rows, NCBI's 28 rows without the index and 3 of the
    index's pairs (`the_two_hit_index_search_adds_the_index_s_pairs`: the row of qhm7 in hm4 is in
    NCBI's own indexed search without the mask, `hm_iu_win`; the other two are runs that a mask
    ends, paired within the window); and NCBI's own two-hit search with the index of the database
    of hm's unmasked stretches, each its own sequence (`hm_win_stretches.out`), gives each of
    LOSAT's rows but for three that NCBI's search without the index and LOSAT's carry across a
    7-letter mask (the traceback reads the whole sequence with the masked letters as B, where the
    stretches end). Random comparison `remaining-gaps-20261112/D/d5/probe/d5b5/rand_d5b5.py` (400
    databases of dustmasker's masks with queries across the masked stretches; NCBI 2.17.0): the
    308 one-hit searches give NCBI's rows without the index; of the 92 two-hit searches, 61 give
    them, and 31 are of databases where NCBI's index and scan pair other hits even without the
    mask. For those 31, NCBI's indexed search of the database of unmasked stretches
    (`st_oracle.py`, `st_explain.py`) gives LOSAT's rows in 17; in 12 the only differences are
    alignments carried across a short mask as above, and in the last 2 (indexes that store every
    19th and every 2nd position, `-legacy false`) the stretch database moves the stored positions,
    after which NCBI's two-hit rule pairs other hits (traced: a root 7 letters later drops the earlier
    hit). No run depended on `-num_threads` (1 and 4) or on the masked letters (the same database
    with every masked letter replaced at random). Unit tests `a_masked_letter_ends_a_run_and_pairs_no_hits` and
    `the_masks_are_of_their_own_sequence` (`LOSAT/src/dbindex/search_tests.rs`), and the placement
    tests of `LOSAT/src/algorithm/blastn/blast_engine/mb_index_seeds_tests.rs`.

## Approved exception 30 (BLASTN `-use_index true` with a MegaBLAST index of only some database volumes, every query batch searched with the index, version 1.101)

Version 1.101, line RG item D5c. The Owner decided on 2026-11-24 12:44 JST: the answer 「Q-MIXIDX=C」 (the
whole answer 「Q-AGVI-CAP=A Q-FM-R2=B Q-D23=A Q-KD-1=A Q-MIXIDX=C Q-GBHELP=A Q-DL-TIE=B Q-DL-QLOC=B
Q-PS7c-SPOUGE=B Q-PS10e=A」) to line RG's question Q-MIXIDX (the query batches and query chunks after
the first of a search with an index of only some database volumes), whose option C read
「どのバッチでも種を失わない正しい検索を行う（11-22 の `Q-HARDIDX=C` と同一方針。PD 例外更新＋`RG-D5` 再開で 1 PR）」
(a correct search that loses no seeds in any batch, the course of the answer C to Q-HARDIDX of
2026-11-22; this record updated and line RG item D5 resumed, in one PR). Not taken: option A, the
form of versions 1.74 to 1.94: NCBI's kept results reproduced (rule 2) at every thread count as its
one thread keeps them, and the explicit rejection (rule 3) of a later search whose scan saves initial
hits into the list that NCBI frees.

30. **BLASTN `-use_index true` with a MegaBLAST index of only some database volumes, after the first
    search.** NCBI keeps an index volume's seeds while a search thread may still need them and counts
    a volume down as a thread passes it, but a thread that goes on to a database volume without an
    index counts none down (`blast_dbindex.cpp:812-857`). So the last index volume that the search
    reaches before each such volume that it reaches (all of them, or those with a sequence in the
    database's OID list) keeps the seeds of the first search for the rest of the run, and each later
    query batch or query chunk takes them as the first search's word finder left them: every list
    moved out and empty (`blast_extend.c:248-259`). There it finds no seed: a sequence whose list the
    first search had is read and searched without seeds, the others are skipped, and the later
    queries' hits in that volume are lost. NCBI's search thread then holds the empty list as its list
    of initial hits, of no memory: the scan of a database volume without an index saves its first
    initial hit with a reallocation of 0 bytes, which glibc answers with a block of one hit, and its
    second with another, which frees the list (`blast_extend.c:325-356`); NCBI then writes to and
    frees freed memory: `double free detected` (exit status 134), a crash (139), or a result read
    from freed memory. With more threads, which volumes keep their results depends on which thread
    passes which volume (NCBI's `-num_threads 4` search of the fixture row `m_mixed1_t4` crashed, exit
    status 139, in 8 of 8 runs, where its `-num_threads 1` search completes). LOSAT searches each index
    volume that a search reaches with that search's query, the first search and the later ones alike
    (`mb_index_seeds::IndexSeeds`), as NCBI's first search does: no seed is lost, no such search is
    refused, and the result does not depend on `-num_threads`. A first search, and an index of every
    volume, keep NCBI's result. Evidence: the rows `batch NAME` of
    `LOSAT/tests/blastn_mb_index/cases.tsv` (`m_mixed1`, `m_mixed1_t4`, `m_mixed2`, `m_mixed1_win` of
    `-window_size 40`, `m_part_b1`, `m_part_b1_t4`, `m_part_two`) give NCBI's own search with the index
    of the same queries as one query batch (`BATCH_SIZE=1000000000`, rows `NAME_1b`), a first search:
    759 rows of `m_mixed1`, where NCBI's search writes 755 (and crashes at 4 threads); 759 of
    `m_mixed2`, where it crashes after 35; 95 of the two-hit `m_mixed1_win`, where it writes 92; 5 of
    `m_part_two`, where it stops after 3 (134). The rows `plain NAME` of query chunks (`m_part_c1200`
    and `m_part_c2000`, `CHUNK_SIZE` 1200 and 2000), whose later chunks no NCBI search with the index
    searches anew, give NCBI's search without the index (5 rows; NCBI's search with the index writes
    4, or none and stops with 134). Each of NCBI's own outputs of these rows (`ncbi/<row>.*`) is that
    whole list less hits in the volume that keeps the results, or a beginning of it before a crash
    (`ncbi_loses_seeds_after_the_first_search` of `LOSAT/tests/blastn_mb_index_fixtures.rs`).
    Random comparison `remaining-gaps-20261112/D/d5/probe/d5c/rand_d5c.py` (596 random databases, 537
    of them of 2 to 9 volumes with 1 to all but one of their `.shd` files removed; queries in several
    query batches (`BATCH_SIZE` 1, 300, 1,000 or 2,500, or NCBI's own) in 405 runs and in query chunks
    (`CHUNK_SIZE` 600 to 2,500, one-hit) in 132; 121 two-hit searches; nmer 8 to 12, 163 indexes of
    `-legacy false` with a stride; 1 and 4 threads): each LOSAT run gave the oracle's rows (NCBI's
    search with the index as one query batch, or for query chunks its search without the index);
    NCBI's own run of the same command lost rows in 163 and stopped with 134 or 139 in 25. Two-hit
    searches in query chunks, which no NCBI search of such a database gives anew, were compared where
    every volume without an index comes before the index volumes, so that NCBI keeps no results
    (`rand_d5c2.py`, 267 databases of 2 to 9 volumes, `CHUNK_SIZE` 600 to 2,500, 1 and 4 threads):
    each LOSAT run gave NCBI's output.

## Approved exception 26 (makeprofiledb's word table past `Int4`, psiblast line item PS10a, version 1.82)

Version 1.82, proposed under the standing instruction 2026-09-29, restated 2026-10-07, on the Owner's
rules of 2026-10-10 「NCBI のバグを再現する必要はない」 and 2026-10-22 「広い型になおしましょう」; accepted by the
maintainer on 2026-11-23 by marking PR #1497 ready for review, as the PR asked. Its version is in the list at
the top; elsewhere it is cited by its number and name, "approved exception 26 (makeprofiledb's word table
past `Int4`)". Psiblast line item PS10a.

25. **A word score threshold that, times the profile's scaling factor, does not fit `Int4`, and
    neighbouring-word score sums past `Int4`.** makeprofiledb indexes the words of every profile of a
    volume with `s_AddPSSMNeighboringWords` (`blast_aalookup.c:623-748`), whose threshold is the
    double `scale_factor * wordScoreThreshold` stored in an `Int4` (`makeprofiledb.cpp:928-937`,
    `1465-1469`) and whose recursion adds the profile's `Int4` scores in an `Int4`. A double past
    `Int4` converted to `Int4` is undefined behavior; NCBI BLAST+ 2.17.0's x86-64 build gives INT_MIN,
    so every word of the profile is a hit (oracle rows V25, S_cwstbig, S_cwstbig1: a `.loo` of
    5,178,168 or 1,841,556 bytes where the threshold's value gives 524,344). Scores large enough that
    a word's sum passes `Int4` (row S_scedge) wrap (signed overflow, undefined) and keep or drop the
    word by the wrapped value. LOSAT keeps the threshold and the sums in 64 bits
    (`algorithm::pssm_search::lookup`, `makeprofiledb::rps_files::word_threshold`), so a word is kept
    by its score: the table NCBI builds with the threshold's value. A threshold that is not a number
    after the scaling stays refused (rule 3). The recursion is the one psiblast and tblastn
    `-in_pssm` use; their PSSM scores stay within `BLAST_SCORE_MIN..=BLAST_SCORE_MAX` and their
    threshold is an `Int4`, so their sums never leave `Int4` and their output does not change.

## Approved exception 31 (psiblast and tblastn `-in_pssm`, the Spouge cutoff search past `Int4`, version 1.111)

Decided by the Owner on 2026-11-24 at 12:44 JST, 「Q-PS7c-SPOUGE=B」 (one answer of the batch 「Q-AGVI-CAP=A
Q-FM-R2=B Q-D23=A Q-KD-1=A Q-MIXIDX=C Q-GBHELP=A Q-DL-TIE=B Q-DL-QLOC=B Q-PS7c-SPOUGE=B Q-PS10e=A」; psiblast line
`DECISIONS.md` 103 and 104). The option chosen reads 「Owner 10-10 の「符号付きオーバーフローは広げる」方針（承認例外 17 の
culling と同じ形）で数学的に正しい値へ広げ、PD 例外として記録する」: under the Owner's rule of 2026-10-10 that a signed
overflow is widened, in the form of the culling limit's exception 4 (line NP entry 17 of `docs/decisions/OWNER_DECISIONS.md`, 「広い型になおしましょう」),
the overflowing values are computed as their mathematical values and recorded here as an exception. Not taken: option
A, the explicit rejection at the overflow that line PS item PS7c's fix round 1 had built (never on dev). Its version
is in the list at the top; elsewhere it is cited by its number and name, "approved exception 31 (the Spouge cutoff
search past `Int4`)". Psiblast line item PS7c, fix round 1 (review B1).

31. **A gapped Lambda for which NCBI's search for the hit-saving cutoff passes `Int4`.** With the Gumbel block of
    the search's matrix, the hit-saving cutoff is the score of the preliminary e-value that `BLAST_SpougeEtoS`
    searches for (`blast_parameters.c:931-945`): from `b = MAX((int)(log(db_length / e) / Lambda), 2)` it doubles
    `b` while the e-value of `b` is above `e`, then bisects (`c = (a + b) / 2`) down to the score `a`
    (`blast_stat.c:5236-5269`), all in `Int4`. The gapped Lambda of a matrix keeps the search far below 2^30 at
    every e-value; a PSSM file's gapped Lambda (psiblast and tblastn `-in_pssm`, `kbp_gap_psi`) need not. The cutoff
    is about t / Lambda, where t (7.4 for the PSSM of 30 letters of the rows below, at `-evalue 10` against
    `tests/db_search/paja`) depends on the e-value, the lengths and the Gumbel block, so it passes 2^31 - 1 below
    about 3.4e-9 for that PSSM; with a Lambda below 0 the e-value grows with the score and no score reaches `e`.
    Then `b *= 2` or `a + b` overflows (signed overflow, undefined behaviour). NCBI BLAST+ 2.17.0's x86-64 build
    (under gdb): at the doubling with a small positive Lambda (1e-9, 1e-300) it wraps `b` to `INT_MIN`, then to 0,
    and never ends; at the doubling with a negative Lambda (-0.5) it wraps `b` to `INT_MIN`, whose e-value is 0,
    and returns 2^30 (1,073,741,824) as the cutoff; at the bisection (psiblast 5e-9 to 7e-9, tblastn 4e-9 and 5e-9
    here) it wraps the sum, takes it as the lower bound and returns a negative cutoff (psiblast 6e-9: 733,828,724 +
    1,467,657,448, the cutoff -1,046,740,562; tblastn 4e-9: 876,829,330 + 1,753,658,660, -832,239,653), with which
    every ungapped hit starts a gapped extension and every gapped HSP passes the cutoff. LOSAT takes NCBI's steps
    with `b`, `a + b` and the scores in 64 bits (`stats::spouge::blast_spouge_etos`, from NCBI's first bound, its
    cast of the double as `core::blast_util::ncbi_int4_from_double`): a search inside `Int4` gives NCBI's score, a
    bisection past it its score in 64 bits (psiblast 6e-9: 1,231,258,856; tblastn 4e-9: 1,689,538,156), and a
    search whose lower bound `a`, a score whose e-value is above `e` that no step lowers, passes 2^31 - 1 ends
    there with `i32::MAX`, the cutoff of NCBI's context that saves nothing (`blast_parameters.c:920`): the cutoff is
    then above every `Int4` score (psiblast 1e-9: about 7.4e9; 1e-300: about 4.7e300; -0.5: none). So the doubling
    runs at most 31 times and the bisection at most 33, for every Lambda. With the PSSM of the rows below,
    whose alignments score far below these cutoffs, each report has no hits.
    Where NCBI ends, its wrapped cutoff keeps no hit either and LOSAT writes its bytes (rows
    `db.ps7c.ps.lambda_6e-9`, `db.ps7c.ps.lambda_neg`, `db.ps7c.tbn.lambda_4e-9` and `db.ps7c.tbn.lambda_neg`, and
    their outfmt 0 rows `.0`). Where NCBI never ends (1e-9), LOSAT writes the bytes of NCBI's run of the near
    Lambda 8e-9, whose search stays inside `Int4` (the cutoffs 923,444,142 and 844,769,078) and which reports no
    hits either, with X-drops of 1 bit that NCBI's kernels size (rows `db.ps7c.ps.lambda_1e-9` and
    `db.ps7c.tbn.lambda_1e-9`, outfmt 7, which prints neither; in outfmt 0 LOSAT's report is the 8e-9 report with
    the footer's Lambda of the file). NCBI's runs of the bisection and of the negative Lambda gap every ungapped
    hit; LOSAT's keep NCBI's cutoff of the words (the gap trigger) and spend no time there.
    `LOSAT/src/stats/spouge_tests.rs` keeps the cutoffs: NCBI's arguments of these rows (gdb) and, for gapped
    Lambdas from 1e-300 to 1e300, below 0, 0, infinite and not a number, NCBI's search in 128-bit integers.

## Approved exception 32 (megablast greedy gap costs past `Int4`, line RG item G1, version 1.115)

Version 1.115, decided by the Owner on 2026-11-26 at 15:41 JST (「承認 A」, choosing option A:
accept megablast gap costs above 32767 up to 2147483647 as an approved exception, following the
Owner's rules of 2026-10-10 「NCBI のバグを再現する必要はない」 and 2026-10-22
「広い型になおしましょう」; not taken: option C, refusing values above NCBI's
input-length-dependent 32-bit overflow boundary). Its version is in the list at the top;
elsewhere it is cited by its number and name, "approved exception 32 (megablast greedy gap costs
past `Int4`)". Line RG item G1.

32. **Megablast affine greedy alignment with gap costs whose 32-bit setup or distance arithmetic
    overflows or exhausts dense-table memory.** `s_BlastGreedyAlignMemAlloc`
    (`core/blast_gapalign.c:171-267`) and `BLAST_AffineGreedyAlign`
    (`core/greedy_align.c:755-1249`) compute doubled gap costs (for an odd reward),
    `GE_cost = gap_extend + reward / 2`, `max_d * GE_cost`, `gap_open + GE_cost`,
    `2 * (max_d + 1 + max_cost)`, and `d * score_common_factor` in signed `Int4`, and allocate
    dense arrays sized by unreduced `max_d` and `max_cost`. When those products overflow `Int4`
    (undefined behaviour in C) or exceed available memory, NCBI BLAST+ 2.17.0 fails in setup
    with `BLAST engine error: Unknown error code -1` (exit status 3, from a sign-extended
    negative `calloc` count or host-RAM exhaustion: for example `-gapopen 1000000 -gapextend
    1000000` on `taxfilter/tfn`, where `max_d_1 = 560` overflows `2 * (max_d + 1 + max_cost)` at
    955286, while 955285 needs 73 GB and exits 0 at 128 GB; or 535799 on a 4 kb subject where
    535798 exits 0 at 128 GB), crashes mid-search with `free(): invalid pointer` or SIGSEGV (for
    example `40000/2` at 8 GB, `2000000000/2`, `1073741824/2`, `2147483647/2147483647`), or wraps
    `2 * gap_open` to a negative penalty (`2147483647/2` with odd reward wraps to -2 and alters
    hits). LOSAT widens the affine greedy gap-cost, GCD, X-drop offset and distance arithmetic
    to 64 bits (`algorithm::blastn::alignment::greedy`) and steps only over non-empty distances
    without dense tables, so every `-gapopen` and `-gapextend` accepted by the nucleotide
    Karlin-Altschul tables (including values above 32767 up to 2147483647) completes in low
    memory. Wherever NCBI completes without `Int4` overflow or memory exhaustion (such as
    `32768/32768`, `40000/40000`, `100000/100000`, `40000/2` at 64 GB, `2/40000`, high-GCD
    `-reward 1000 -penalty -2000 -gapopen 40000 -gapextend 2000`, `955285/955285` on `tfn` at
    128 GB, `1000000/1000000` on a 600 bp pair at 64 GB), LOSAT gives NCBI's exact bytes (with
    LOSAT wording). Where NCBI overflows `Int4`, LOSAT writes the alignments of NCBI's nearby
    non-overflowing run (the `ncbi-near` rows `db.g1.*` of `LOSAT/tests/db_search_manifest/blastn.g1.tsv`
    and the `g1.*` cases of `LOSAT/tests/blastn_regression_fixtures.py`). Before version 1.115
    LOSAT refused megablast gap costs above 32767 (rule 3).

## Approved exception 33 (an imported strategy's subjects at more than one thread, version 1.118)

Version 1.118, line RG item F2 (DECISIONS F2-Q1). The Owner
decided on 2026-11-20 at 20:05 JST, with the words 「Q1=C」 after asking 「重複すること自体バグでは？」 ("is the
duplication itself not a bug?"), that a search strategy whose subjects are sequences, imported with
`-num_threads` above 1, gives the output of `-num_threads 1` in every program. Not taken: NCBI's copies
(option A, rule 2 as written: the copies come from a check that reads the wrong argument, so the same
alignment printed N times is not a result of the search), and a refusal of the combination (option B:
the search itself is valid and its `-num_threads 1` result is deterministic).

33. **An imported strategy's subjects at more than one thread.** NCBI decides whether a search takes
    `-num_threads` from the command line's `-subject` only (`CMTArgs::x_ExtractAlgorithmOptions`, NCBI
    reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3203-3236): a strategy that carries its subjects
    (`s_ImportSearchStrategy`, NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:532-628)
    leaves the threads on, and each thread searches every subject. In our runs of NCBI BLAST+ 2.17.0
    `blastp`, `blastx` and `tblastn` at 2 to 8 threads, and `blastn` at 2 and 3, print each alignment
    once per thread (probes `I.<program>.sub.nt3.r1` to `.r3`, `I.<program>.sub.nt4`,
    `I.blastp.sub.nt4f6`, `I.blastn.sub.nt2.*`); `blastn` at 4 and 8 threads and `tblastx` at 2, 3, 4
    and 8 stop with `malloc(): invalid size (unsorted)` and a core dump (SIGABRT; `I.blastn.sub.nt4`,
    `I.blastn.sub.nt8.*`, `I.tblastx.sub.nt2.*` to `I.tblastx.sub.nt8.*`). The copies are rule 2's deterministic NCBI output, which this exception
    replaces by the Owner's decision; the aborts are rule 1. LOSAT searches a strategy's subjects
    once, with the command line's `-num_threads`, as it searches the records of `-subject` (its
    output is the same at every thread count, the Owner's rule of 2026-10-13): the engine's import
    (`import::apply`, which the command line, the library API and LOSAT Web share, also where one
    program runs another's strategy) keeps the threads, with no warning and no refusal, and
    the output is NCBI's `-num_threads 1` output, frozen. On a terminal, `losat blastn` draws its
    progress bar for them as for the records of `-subject` (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES
    exception 11). `-mt_mode 1` gives the same bytes (its
    output is that of one thread, as for every search). NCBI's warning `Number of threads was reduced
    to N to match the number of available CPUs` is not printed, as for every search (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES). A strategy whose subject is a
    database, and the command line's `-subject` or `-db` over a strategy, keep their usual threads.
    Evidence: the probes of line RG item F2 named above (NCBI BLAST+ 2.17.0).
    Fixture rows: `LOSAT/tests/fixtures/strategy/import/import_cases.tsv` holds the `-num_threads 1`
    rows `I.<program>.sub.alone`; the test
    `a_strategys_subjects_are_searched_once_at_any_thread_count`
    (`LOSAT/tests/import_search_strategy.rs`) runs them at 2, 4 and 8 threads against those bytes for
    the five programs (a search of the subjects once per thread fails it), and
    `a_strategys_subjects_keep_the_threads_of_the_command_line`
    (`LOSAT/src/report/archive/read/strategy_tests.rs`) checks that the import keeps `-num_threads 4`.

## Approved exception 34 (a strategy's negative window and intron length past `Int4`, and a negative window at more than one thread, version 1.119)

Version 1.119, line RG item F2 fix round 3 (review 3, R3-B2; DECISIONS F-Q28). The Owner decided on
2026-11-24 at 00:04 JST, on the board, with the words 「Q-F2-28=A」 (option A of line RG's question
Q-F2-28, the recommended one; proposed 2026-11-23 under the standing instruction of 2026-09-29, restated
2026-10-07), under the Owner's rules of 2026-10-10 「NCBI のバグを再現する必要はない」 and 2026-10-22
「広い型になおしましょう」 and rule 1 above. Not taken: NCBI's wrapped values (option B: undefined behavior in
C, whatever NCBI's build gives, and NCBI's result at several threads, which changes from run to run,
below, cannot be reproduced), and a refusal of the values at the import (option C: the searches are valid
and their values of a nearby input are deterministic).

34. **A negative `WindowSize` or `LongestIntronLength` of an imported strategy.** The command line takes
    0 or more for both (`-window_size`, `-max_intron_length`); a hand-edited strategy can carry any `int`,
    and NCBI runs it with no check.
    - *Window, arithmetic.* A window of 0 or less runs the one-hit word finder over NCBI's diagonal
      table, whose state passes from one word finder call to the next: the table's offset grows by the
      call's subject length plus the window after each call and is cleared at `INT4_MAX / 4`
      (`Blast_ExtendWordExit`, NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-185),
      and a word is extended where its subject offset passes the diagonal's last hit less the offset
      (aa_ungapped.c:715-810, na_ungapped.c:640-945). From a window of about -2^29 down, the offset of a
      search passes the range of `Int4` within a few calls, an overflow of signed values (undefined
      behavior). NCBI BLAST+ 2.17.0's build wraps: `tblastx` with the strategy `tblastx.db.ss` and
      `WindowSize` -536870912 reports 11 alignments, -2^30 14, -1.5 x 10^9 43, -2147483647 and -2^31 27,
      where every window from -1000 to -10^6 (defined arithmetic) reports none (probes
      `F/f2/tmp/wneg_tblastx/out_<w>.tsv` of line RG item F2). LOSAT computes these sums in 64 bits, their
      mathematical value, and keeps NCBI's comparisons, its reset and the 31 bits of a diagonal's last hit
      (a defined truncation): it reports no alignment for those windows, as NCBI does for -10^6. Where
      NCBI's arithmetic stays in `Int4` (the windows of the frozen rows: -1, -5, -40, -1000, -10^5,
      -10^6), LOSAT's bytes are NCBI's.
    - *Window, threads.* NCBI's threads take the subjects (or database chunks) by timing and each keeps
      its own table, so the state that one search passes to the next depends on the schedule: in our runs
      of NCBI BLAST+ 2.17.0 with `WindowSize` -40 at 4 and 8 threads (3 runs each), `blastp` gave 2
      different outputs in 6 runs, `tblastn` 3 and `tblastx` 3, and those of `tblastn` and `tblastx` all
      differ from their `-num_threads 1` output. LOSAT gives NCBI's `-num_threads 1` output at every
      thread count: the word finder calls of a negative window run in NCBI's one-thread order on one
      table (rows `S3.<program>.window_neg40.nt4`, which hold NCBI's `-num_threads 1` bytes).
    - *Intron length.* `blastx` and `tblastn` with sum statistics subtract 2 from the length and divide
      by 3 (blast_parameters.c:788-815): for -2^31 and -2^31 + 1 the subtraction passes `Int4`, which
      NCBI's build wraps to a large positive gap. LOSAT subtracts in 64 bits: the gap is negative, so a
      gapped search links no HSPs, as NCBI does for -2^31 + 2 and every other negative length
      (`core::link_hsps::translated_longest_intron`).
    Evidence: the probes named above; rows `S3.*.window_*` of
    `LOSAT/tests/fixtures/strategy/import/import_cases.tsv` (NCBI's bytes where its arithmetic is
    defined), the unit tests of the negative window in the five engines and
    `the_longest_intron_of_a_translated_search_follows_ncbis_branches`
    (`LOSAT/src/core/link_hsps_tests.rs`).

## Approved exception 35 (`blastdbcheck` alias count checks and `-ends` OID arithmetic past signed overflow, version 1.124)

Decided by the Owner on 2026-11-25 at 17:54 JST, 「Q-3C2-STALE=B Q-D2b-OV=A Q-KB9A-SPEED=A Q-AGVI-CAP2=A
Q-AGVI-MODEL=A Q-CITE-ALL=B」 (`qbatch-20261125/ANSWERS.md` 2; `remaining-gaps-20261112/DECISIONS.md` OQ-1125).
The chosen option A reads 「A: 広い型で正確に計算し誤判定を出さない（PD に新例外）」. Not taken: (B) reproducing
NCBI's signed wraparound; (C) refusing overflowing inputs.

35. **`blastdbcheck` alias count checks and `-ends` OID arithmetic past signed 32-bit or
    64-bit overflow (Owner 2026-11-25 17:54 JST, Q-D2b-OV=A, version 1.124).** NCBI BLAST+ 2.17.0
    computes `blastdbcheck`'s alias-file count guards (`FIRST_OID`, `LAST_OID`, `MAXOID`, `NSEQ`,
    `LENGTH`, `OIDLIST` combinations in `blastdbcheck.cpp:782-805`) and `-ends` OID bounds
    (`blastdbcheck.cpp:666-682`) in signed 32-bit `int` and 64-bit `Int8`, where
    products such as `(FIRST_OID + 1) * (LAST_OID + 1)` (e.g. `65535 * 65535`) or `num_oids - n`
    for `-ends -2147483648` overflow signed integer arithmetic (undefined behavior in C++): NCBI's
    x86-64 build wraps them to negative or zero values, emitting false alias errors (rows `a_ov1`,
    `a_ov2`, and `-dir .` rows `a_dir` and `a_dir_v1` that visit them) or wrapping `-ends -2147483648`
    back into the valid OID range (row `ends_intmin`, passing exit 0 while nearby `-ends -2147483638`
    jumps to `oid 2147483647` and fails with `OID not found`). LOSAT computes the alias count
    expressions in 128-bit signed integers (`blastdbcheck::alias::check_values`) and
    `-ends` OID arithmetic in 64-bit signed integers (`blastdbcheck::volume`), avoiding
    signed overflow and producing exact results (`a_ov1` and `a_ov2` pass with exit 0 while keeping
    non-overflow warnings such as `.msk` extension warnings; `a_dir` and `a_dir_v1` drop those two
    false errors and report 19 alias errors and 37 total errors instead of 21 and 39; `ends_intmin`
    jumps past the volume to `oid 2147483657` and fails with `OID not found`, matching nearby
    `-ends -2147483638`).

## Approved exception 36 (a protein database with a sequence of length 0 at more than one thread, version 1.125)

Decided by the Owner on 2026-11-27 at 09:25 JST, in the Antigravity (Gemini) session that recorded it
(`gemini-collab-20261122/DECISIONS.md` 48; `remaining-gaps-20261112/DECISIONS.md` OQ-1127), with the words
「dev 既存の「長さ 0 protein DB の MT abort」の扱い（後続項目）：有効な出力をそのまま残す」 ("keep the valid output as it
is"). Not taken: a refusal (rule 3) of the searches where NCBI aborts.

36. **A protein database with a sequence of length 0, searched at more than one thread (Owner
    2026-11-27 09:25 JST, Q-DBA0-MT=KEEP, version 1.125).** At `-num_threads` 2 or more NCBI BLAST+
    2.17.0 reads the database through one buffer per thread. `GetSequence` checks a sequence out of the
    thread's buffer (`CSeqDBImpl::x_GetSeqBuffer`, NCBI reference (598d8ae6):
    c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:667-681). For a length of 0 the search's reader
    returns an error without giving the sequence back (`s_SeqDbGetSequence`, NCBI reference (598d8ae6):
    c++/src/algo/blast/api/seqsrc_seqdb.cpp:378), the engine goes on to the next sequence (NCBI
    reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1425-1427), and the thread's next fill
    of its buffer throws `Sequence not returned.` (NCBI reference (598d8ae6):
    c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:656-661). In our runs `blastp` (gapped and
    `-ungapped`), `blastx` and `psiblast` at 2, 4 and 8 threads stopped in each of 3 runs on each of 8
    databases: the sequence of length 0 first, in the middle or last of one volume of 2, 3 or 300
    sequences, or in the first, the last or the middle volume of an alias of 2 or 3 volumes. Each printed
    `BLAST engine error: Database memory map file error` (`psiblast`: `BLAST Database error: Sequence
    not returned.`), then ended on the uncaught `CSeqDBException` with SIGABRT (exit 134) and no output.
    At 1 thread each search finished with exit 0. LOSAT skips a sequence of length 0 at every thread
    count, as NCBI does at 1 thread (`seqdb::Volume::read_protein_search_block`), and prints NCBI's
    `-num_threads 1` output at any thread count, with no warning and no refusal. Evidence:
    [`docs/evidence/zero_length_protein_mt/`](../evidence/zero_length_protein_mt/README.md) (the probe
    and its rows, NCBI BLAST+ 2.17.0). Fixture rows: `LOSAT/tests/fixtures/zero_length_mt/` holds
    NCBI's `-num_threads 1` output of `blastp`, `blastp -ungapped`, `blastx`, `blastx -ungapped` and
    `psiblast` on `dba0` (a volume whose one sequence has length 0, then a volume of one sequence); the
    test `a_protein_of_length_0_keeps_the_one_thread_output_at_any_thread_count`
    (`LOSAT/tests/protein_db_32bit_lengths.rs`) runs them at 1, 2 and 4 threads against those bytes.

## Explicit rejections of deterministic results (rule 2, version 1.10)

On a database with parsed Seq-ids NCBI gives these results without failing, and they look wrong
or need NCBI's object manager to reproduce. LOSAT refuses them instead (line DBID
`DECISIONS.md` 9, 20 and 21, line DI-7 for the query's id), under the maintainer's standing
delegation, which the maintainer accepted on 2026-11-02 with version 1.10. Each is a row of the table
in [`PD-LOSAT-DATABASE-SEARCH`](PD-LOSAT-DATABASE-SEARCH.md), which also lists refusals that rule 4
would port. The last five entries are not of these. The last (version 1.102, psiblast line item PS10e), the one
before it (version 1.81, line NP item M4c-1), the one before that (version 1.60, line RG item D5a)
and the one before those (version 1.54, line RG item X5) say so themselves. The one before them (version 1.45, line RG item D3-1) is not either: it concerns a PSSM
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
- A stored Seq-id with a byte outside ' ' to '~', and a stored title with such a byte in the
  cases narrowed by version 1.120 (line RG-KB item 9(b), Owner-delegated: `psiblast`, `rpsblast`
  and `rpstblastn`; a format other than `-outfmt` 8, 9, 11 and 12, or 8, 9 and 12 with
  `PRE_FETCH_SEQS_LIMIT=2147483647`; a search restricted by a seqid, GI or taxid list, an alias
  file with a list or an OID mask, or a volume with an OID mask type; a query batch whose headers
  also give NCBI's warning of a UTF8String; except a query batch that ends at a header NCBI cannot
  read, which writes NCBI's `Error: (803.7)` messages for the titles read before it, then NCBI's
  error, in any format of all eight search programs): NCBI prints `#` and an error message each
  time it reads the definition line. Otherwise BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX write
  NCBI's output and `Error: (803.7)` messages.
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

- makeblastdb `-gi_mask`, two masking algorithms that would write one GI-mask data file at the
  same time (version 1.81, 2026-11-23, line NP item M4c-1, the Owner's decision; not of the
  database class above). Two algorithms whose names
  lead to the same files (one `-gi_mask_name` twice, two names of one path such as `x` and `./x`,
  or a name that reaches the same directory through a link; on Windows and macOS also two names
  that differ only in case) each get a `CWriteDB_GiMask` (`writedb_impl.cpp:1453-1456`) with its
  own streams on the same paths (`writedb_gimask.cpp:46-57`). NCBI exits 0. Where the files are
  written one after the other, its bytes follow from the order of the writes, and LOSAT writes
  them as NCBI does: an algorithm that reaches no sequence with a GI creates no file
  (`writedb_gimask.cpp:92-96, 290-295`), the offset and index files are each created, written
  and closed in turn at the close (`:109-119`), and a data volume that one algorithm closed before
  the other created it is created anew (`:70-80`). Where a data volume of one algorithm would be
  created while the other's of the same number is open, NCBI writes both algorithms' records
  through two streams into one file, and its bytes then depend on how the C++ runtime buffers two
  streams on one file (survey D6, probe `m04_dupname`: the files hold the second algorithm's
  masks); LOSAT does not reproduce that, and its model of the runtime's file streams
  (`writedb/ofstream.rs`) is under licence audit (2026-11-21). LOSAT refuses the record that
  would create that volume (exit 1, "losat makeblastdb: -gi_mask: two masking algorithms would
  write one GI-mask data file at the same time ..."). The build ends there as an error of NCBI's
  `SetMaskData` ends it (case `x03gi`): the volume files are removed, and the GI-mask files of
  the other algorithm and the alias file stay. The Owner answered 「推奨通り」 on 2026-11-21 at
  23:07 JST to question Q-GIDUP, whose recommended option A was to keep this rejection; option
  B, to reproduce the C++ runtime's buffering, was not taken. Unit tests in `LOSAT/src/writedb/gi_mask_tests.rs`, cases `w25`, `w25b` and `m04`
  (refused) and `x11_1`, `x11_2` and `x12` (written in turn) of `LOSAT/tests/seqdb/gi_mask`
  (`make_gi_mask_fixtures.py`).

- RPSTBLASTN without composition-based statistics, an alignment that starts at or past the `.psq`
  letters of a profile whose letters are fewer than its `.rps` columns (version 1.102, 2026-11-23,
  psiblast line item PS10e, the rule-3 entry of such profiles above; not of the database class
  above; decided by the maintainer on 2026-11-24, the answer 「Q-PS10e=A」 at 12:44 JST, psiblast
  `DECISIONS.md` 103). NCBI's rows
  there are the mapped `.rps` rows from the profile's offset (`blast_traceback.c:1194-1198`), so its
  run ends with exit 0 and a clean valgrind unless the alignment then reads a profile letter of 28
  or more (read E of that entry), but its gapped traceback skips the extension to the
  right and ends the HSP before its start (`blast_gapalign.c:4618-4651`), counts the identities on
  the `.psq` bytes after the profile's letters, and prints HSPs whose columns, scores and identities
  disagree (probes N_AA1_nq1_c0: profile columns 1 to 10 of a 6-letter profile; N_AA1_nqt2_c0: a bit
  score of 0.0 with an e-value of 3.38e-11, columns 26 to 31 of a 22-column profile; N_P6d_nq1_c0:
  the archive's `num_ident` 9 where the tabular formats count 6). LOSAT refuses it with the
  text of read D of that entry, "starts an alignment past the consensus" (rows
  `db.rpstbn.rej.odd.aa1.nq1c0` and `db.rpstbn.rej.odd.sc.nq2c0`; in the second NCBI's run goes on
  to read E of that entry and ends with SIGSEGV, probe N_SC_nq2_c0, as in E_NA1_nqi12_c0 to
  E_NA1_nqi14_c0 of the line's `PR/ps10e/probe2`), as its own arm
  (`ProfileRead::StartInMapping` of `LOSAT/src/algorithm/rps_search/profile.rs`). The other option,
  which the maintainer rejected, to print NCBI's bytes, would have changed besides that arm: the profile's letters (NCBI's traceback reads the
  mapped `.psq` bytes after its `L_pin` letters: the profile's trailing 0, the next profile's
  letters, past the file's end the zeros of its last page; LOSAT holds the `L_pin` letters); the
  traceback's identities and positives on those bytes (`final_hsp` of
  `LOSAT/src/algorithm/pssm_search/traceback0.rs` and the test of read E); the reports' identities
  on them (`formatter_identities` of `LOSAT/src/algorithm/rps_search/traceback.rs`); the
  archive's count of the traceback beside the formatter's recount, which differ there (9 and 6 in
  N_P6d_nq1_c0); and the reports' rows of the profile past its `L_pin` letters, which NCBI's
  formatter prints as `-` (sseq `IVLARI----` in N_P6d_nq1_c0, `-----B----` in N_AA1_nq1_c0) where
  LOSAT's renderer indexes the letters (`display_alignment` of
  `LOSAT/src/algorithm/blastx/report.rs`: sseq, the pairwise rows, btop, the XML and JSON hit
  sequences). A run such as N_SC_nq2_c0, whose alignment then reads a letter of 28 or more,
  stays refused, at read E. LOSAT's gapped traceback already skips the extension to the right as
  NCBI's does. This list is what reading the code found; no such HSP has been run through LOSAT's
  later stages (the e-values and bit scores, the order of the hits, the reports).

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
  NCBI searches as with the largest e-value (`BLAST_Cutoffs` keeps the cutoff 1, Spouge's search
  of the cutoff (`blast_stat.c:5236-5269`) gives 0, no HSP is reaped: a reap drops an HSP whose
  e-value is above the limit). NCBI's tests that keep what is at most `-evalue` are false for NaN:
  the composition-based redo (`blast_kappa.c:3687`) and the Smith-Waterman redo keep no match
  there in a database search, while RPS-BLAST's traceback, whose redo has no sequence source, gets
  every list back as the redo made it (`blast_kappa.c:3687-3733`): the scores stay scaled by 32
  (`0.0 bits (880)`), the bit scores 0 and the identities uncounted (the Seq-align's `num_ident`
  is 0; the reports count them from the letters). LOSAT rejected them before Session S07+++b. Since line RG item X-evalue-nan (version
  1.127) every program gives NCBI's output for NaN: until then BLASTX and RPSTBLASTN refused it,
  the search of Spouge's cutoff raised it to the smallest normal double (BLASTP, PSIBLAST, RPSBLAST,
  TBLASTN and BLASTX with a Gumbel block, that is without `-searchsp`: no hit), and TBLASTN,
  RPSBLAST, PHI-BLAST and the redo of PSSM searches reaped with a test of at most `-evalue` (no
  hit), while TBLASTN's and BLASTX's redo had no test of the best e-value and RPS-BLAST's redo
  normalized every list. An infinite value runs
  without the composition-based redo (rows `db.evnan.*`); with it, the explicit rejection of rule 3
  (decision D12, above). `-inf` and `-1e999` stop at NCBI's option check in every program.
- TBLASTX and TBLASTN `-db_soft_mask` with masks stored out of order on a subject frame of more than
  5,000,000 residues (version 1.89, 2026-11-23, line NP item M4b-T2; Owner-delegated, standing
  instruction 2026-09-29, restated 2026-10-07). NCBI moves a chunk's ranges to the chunk's offsets
  and raises only the first one's left to 0 (`blast_engine.c:286-307`), so a later range of masks
  stored out of order can start before the chunk (NCBI cuts each translated frame into chunks of
  5,000,000 residues, `blast_engine.c:246-250,805-812`: subjects of more than about 15 Mbp). NCBI's
  scan reads the frame's letters there and keeps the hits at the chunk's offsets, below 0 as `Uint4`
  (`blast_aascan.c:76-77,110`). The two-hit test reads them as `Int4` (`aa_ungapped.c:516-547`), the
  extensions read the frame's letters at them and walk left to the chunk's first letter and right to
  its last (`aa_ungapped.c:1089-1158`), and the one-hit extension gives an alignment that starts
  before the chunk (`aa_ungapped.c:1020-1083`). A first range raised to one word on the frame's last
  letter (`aa_ungapped.c:499-500`) holds the NULLB after the frame, which no word holds, and the
  scan goes on with the later ranges (`masksubj.inl:50-58`). The HSPs and seeds keep the chunk's
  offsets as `Int4`, below 0 before the chunk, until the chunk's offset moves them
  (`blast_hits.c:150-189,3037-3051`), and the sort of a chunk's HSPs reads them so
  (`blast_hits.c:1330-1353`). LOSAT read past its slice there (TBLASTX; a debug build panicked) or
  missed those words, and ended the chunk's scan at such a first range (TBLASTX). The probes of
  2026-11-23 (NCBI BLAST+ 2.17.0; five databases of 15,300,000-base subjects with three or six masks
  stored out of order, both strands, leads of 1, 2 and 3,899,900 letters, a range across the chunks'
  boundary, alignments before the chunk, across its start, and from a seed at its offset -1; 1 and 2
  threads, each run twice, the same bytes) gave rows found only from the letters before a chunk for
  TBLASTX and TBLASTN, default, `-window_size 0`, `-ungapped -comp_based_stats 0` (also with
  `-window_size 0`) and `-in_pssm`. LOSAT gives NCBI's output there: the chunks of a frame (TBLASTN:
  the units of a subject) whose scan reads letters before a chunk run in order on one diagonal
  table, as NCBI's thread runs them, on a table made for the subject (`tests/tblastx_db_soft_lead.rs`,
  NCBI's rows at 1 and 4 threads). Where NCBI's thread starts its table again between two chunks of
  such a frame (its offset past INT4_MAX / 4, `blast_extend.c:166-172`), NCBI's result depends on
  the subjects that thread searched before; LOSAT gives NCBI's `-num_threads 1` result there at
  every thread count (PD-LOSAT-DATABASE-SEARCH, DB-3; version 1.104, line NP item M4b-T3, line NP
  decision 58, Q2 = B, 2026-11-23, Owner-delegated, standing instruction 2026-09-29, restated
  2026-10-07). At the first such sequence of a search it walks NCBI's one-thread offset over the
  searched sequences in the database's order (`blast_engine.c:1377`, one table per thread: a
  sequence of fewer than 3 × word + 2 bases is skipped, `blast_engine.c:1318-1331,1368-1369,1429`;
  each of the six frames of every other one is cut into chunks of 5,000,000 residues, a chunk
  without a range searched by no word finder, `blast_engine.c:478-483`; the exit of each chunk's
  word finder moves the offset by the chunk's length and the window, or clears the table at
  INT4_MAX / 4, `blast_extend.c:87-107,162-173`), and searches each such sequence's chunks in order
  on a table at its offset whose entries are a new or a cleared table's (0, or the window below 0),
  so that its clears fall where NCBI's do. The entries of the sequences before it are never within
  reach of its hits. The probes of 2026-11-24 (databases of 284 Mbp in which the clear falls after
  the first chunk of the long sequence's frame +1, or, one base shorter, after its second chunk)
  gave rows at one thread that differ from NCBI's at 2 and 4 threads, and a gapped search that NCBI
  completes at one thread and crashes on at 2 and 4; LOSAT gives the one-thread rows at 1, 2, 4 and
  8 threads. The TBLASTN gapped search of an initial alignment that starts before its chunk is
  refused only where NCBI crashes (rule 3, above).
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
