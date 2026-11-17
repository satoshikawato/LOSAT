- BLASTN: with an ungapped X-drop of 2^30 or more (`-xdrop_ungap 1e10`, which NCBI BLAST+ turns into
  INT_MIN), the debug and test builds no longer abort, `-ungapped` searches no longer crash, and
  searches with a word size of 4 to 8 (such as `-task blastn-short`) give NCBI BLAST+'s hits. Where
  NCBI's ungapped walk runs past two sentinels and its 32-bit sum wraps, the saved hit starts left of
  its query context; LOSAT now keeps NCBI's signed start and length for it (in the sort, the
  containment test, the greedy seed and the `-ungapped` re-evaluation, which reads it from NCBI's
  start), and the word sizes 4 to 8 extend on the whole query block as NCBI does. Default searches
  are unchanged.
- BLASTN `-ungapped`: an alignment whose re-evaluation keeps its part in an earlier query context,
  which NCBI BLAST+ keeps with a query start left of its context, goes through the later steps as
  NCBI's does (its e-value, percent identity, query coverage, `-max_hsps`, `-max_target_seqs` and
  `-culling_limit` ranks), so a search whose limits drop it gives NCBI BLAST+'s report; where it
  reaches the report, on which NCBI BLAST+ stops with an exception (plus strand) or shows query
  coordinates outside the query (minus strand, or from a database), it is refused (a result LOSAT
  does not reproduce yet; README, What LOSAT refuses).
