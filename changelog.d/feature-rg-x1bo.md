- BLASTN: with a raw ungapped X-drop of 2^30 or more (`-xdrop_ungap 1e10`, which NCBI BLAST+
  turns into INT_MIN), the debug and test builds no longer abort, `-ungapped` searches no longer
  crash, and searches with a word size of 4 to 8 (such as `-task blastn-short`) give NCBI BLAST+'s
  hits. Where NCBI's ungapped walk runs past two sentinels and its 32-bit sum wraps, the saved hit
  passes an end of its query context; LOSAT now keeps NCBI's signed query offsets for it and for
  the alignments made from it (in the sort, the containment test, the `-ungapped` re-evaluation,
  and the megablast greedy alignment, which runs from NCBI's seed on the whole query block, on
  both sides of the context), and the word sizes 4 to 8 extend on the whole query block as NCBI
  does. So where the greedy alignment from such a hit comes back into the context, a megablast
  search gives NCBI BLAST+'s alignment instead of other alignments, also where the query batch is
  split into chunks (a batch of more than 1,000,000 letters, 5,000,000 for megablast, or a smaller
  `CHUNK_SIZE`). Default searches give the same output in every case checked; a hit that the
  approximate two-hit extension keeps below the cutoff and that starts left of its context, which
  a default search can reach (the debug build stopped at an assertion there), now sorts and aligns
  as NCBI BLAST+'s does.
- BLASTN: an alignment that NCBI BLAST+ keeps outside its query context (an `-ungapped`
  alignment whose re-evaluation keeps its part there, or a megablast alignment whose greedy
  extension stays there) goes through the later steps as NCBI's does (its e-value, percent
  identity, query coverage, `-max_hsps`, `-max_target_seqs`, `-culling_limit` and
  `-subject_besthit` ranks), so a search whose limits drop it gives NCBI BLAST+'s report; where it
  reaches the report, on which NCBI BLAST+ stops with an exception (query coordinates that leave
  the query on the left) or shows query coordinates past the query's end, it is refused, on both
  sides of the context (a result LOSAT does not reproduce yet; README, What LOSAT refuses). This
  includes a megablast alignment on a subject of 90,000 letters or more whose traceback starts
  more than 3,000 letters outside its context, which LOSAT left out without a message.
