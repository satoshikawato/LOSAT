- `rpsblast` and `rpstblastn` search databases whose alias files give `FIRST_OID`, `LAST_OID`,
  `STATS_NSEQ`, `STATS_TOTLEN`, `MIN_SEQ_LENGTH`, or a `LENGTH` or `NSEQ` of 0, negative or of 2^63
  or more, as NCBI BLAST+ 2.17.0 does (LOSAT refused them): a database of one volume with all its
  keys, one of several volumes each volume with its own length, count and shortest length under the
  database's length and count where these are above 0; the footers' search space is the first
  volume's and the search strategy's later requests hold the database length and the last volume's
  search space. `rpstblastn` now reads `STATS_NSEQ` and `STATS_TOTLEN` and each volume's shortest
  length, and refuses `NSEQ 0` over one volume with `-searchsp` or `OLD_FSC`, where NCBI BLAST+
  crashes; `rpsblast` searches that `NSEQ 0` over several volumes, where NCBI does. Refused: a hit of
  a profile that `FIRST_OID` and `LAST_OID` leave out (NCBI's report does not find it), a `LENGTH`
  and `NSEQ` over one volume whose sum is below the profiles' concatenated length (for now) or above
  2,147,483,647, or whose diagonal table passes 2^30 entries, and with `-export_search_strategy` a
  database length of exactly 2^63 over several volumes.
