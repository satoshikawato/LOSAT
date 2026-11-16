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
  over one volume of a `LENGTH` of 2^63 or more), and a total length of 0 or less keeps the order of
  the hits after `-max_hsps` and `-qcov_hsp_perc`, as NCBI's traceback does. `rpstblastn` now reads
  `STATS_NSEQ` and `STATS_TOTLEN` and each volume's shortest length, and refuses `NSEQ 0` over one
  volume with `-searchsp` or `OLD_FSC`, where NCBI BLAST+ crashes; `rpsblast` searches that `NSEQ 0`
  over several volumes, where NCBI does. Refused: a hit of a profile that `FIRST_OID` and
  `LAST_OID` leave out (NCBI's report does not find it), a concatenated length above 2,147,483,647
  (NCBI frees the profiles twice) or whose diagonal table passes 2^30 entries (NCBI's loop never
  ends), a query whose two gapped alignments start past a concatenated length shorter than the
  profiles (NCBI's interval tree splits without end), an ungapped extension past the last profile of
  a concatenated length longer than them (NCBI reads past its score matrix), and with
  `-export_search_strategy` a database length of exactly 2^63 over several volumes.
