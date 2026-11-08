- Engine: the last two steps of the BLASTN traceback of a subject (the purges and the
  re-evaluation after the traceback loop, then the score sort, the final interval tree and the
  e-values) move out of `search_query_batch` into their own module. No output changes.
