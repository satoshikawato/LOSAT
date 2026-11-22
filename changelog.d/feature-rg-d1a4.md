- `rpsblast` and `rpstblastn` search databases whose alias files give `FIRST_OID`, `LAST_OID`,
  `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH`, or a `LENGTH` or `NSEQ` of 0, negative or of 2^63
  or more, as NCBI BLAST+ 2.17.0 does (LOSAT refused them): a database of one volume with all its
  keys, one of several volumes each volume with its own length, count and shortest length under the
  database's length and count where these are above 0; the footers' search space is the first
  volume's and the search strategy's later requests hold the database length and the last volume's
  search space. Over one volume the `LENGTH` and `NSEQ` also set the concatenated length that NCBI's
  search covers (their sum): the diagonal table, the ungapped extensions to the right and the gapped
  stage's interval tree follow it, shorter or longer than the profiles. The preliminary hit cutoff
  without the Gumbel block reads NCBI's preliminary search space (it differs from the traceback's
  over one volume of a `LENGTH` of 2^63 or more, and with `-dbsize` over one volume of a total of 0 or less), and a total length of 0 or less keeps the order of
  the hits after `-max_hsps` and `-qcov_hsp_perc`, as NCBI's traceback does. `rpstblastn` now reads
  `STATS_NSEQ` and `STATS_TOTLEN` and each volume's shortest length, and refuses `NSEQ 0` over one
  volume with `-searchsp` or `OLD_FSC` where a query is searched, where NCBI BLAST+ crashes (a run
  of queries without letters it searches is written as NCBI's, also by `blastp` and `psiblast`
  over `NSEQ 0`, which refused it); `rpsblast` searches that `NSEQ 0`
  over several volumes, where NCBI does; both fail as NCBI where the average subject length of one
  volume is below one letter without the Gumbel block. Refused: a reported hit (within
  `-max_target_seqs`) of a profile that `FIRST_OID` and `LAST_OID` leave out (NCBI's report does
  not find it; one past that cut of several volumes is searched and counted as NCBI's), a
  concatenated length above 2,147,483,647 (NCBI frees the profiles twice) or whose diagonal table
  passes 2^30 entries (NCBI's loop never ends), both where NCBI's search meets them (an empty query
  file still ends as NCBI's), a query whose two gapped alignments start past a concatenated length
  shorter than the profiles plus one (NCBI's interval tree splits without end), an ungapped
  extension past the row after the last profile of a concatenated length longer than them (NCBI
  reads past its score matrix; the profiles plus one is searched as NCBI). Each search's diagonal tables are no longer than its profiles and longest query
  need (NCBI's length or less, the same hits): no table of 4 GiB per thread near 2^30.
- `rpstblastn` over several volumes that hold one profile (a volume given twice under two paths)
  writes both volumes' hits of it within `-max_target_seqs`, as NCBI BLAST+ 2.17.0 and `rpsblast`
  do: NCBI's cut counts following hits of one Seq-id as one subject (LOSAT wrote one volume's).
  With `FIRST_OID` and `LAST_OID`, such a profile that the range leaves out but whose Seq-id an
  included copy of its record (the same definition lines and letters) holds is written as that
  copy (NCBI's report finds it by its Seq-id), in both programs; LOSAT dropped it past the cut.
- A database length of exactly 2^63 (an alias file's `LENGTH` 9223372036854775808, NCBI's 64-bit
  -9223372036854775808) and `-dbsize -9223372036854775808` are written as `-9223372036854775808` in
  BLAST XML and XML2 (`-outfmt` 5, 14, 16), the search strategy and the BLAST archive, where NCBI
  BLAST+ 2.17.0's writer prints `-0` (the Owner's decision of 2026-11-17; LOSAT refused them).
