- Engine: the statistics a TBLASTX subject reads, the end of each subject (its final HSPs to the
  hit lists), the subject's linking parameters and the stages after its re-evaluation are made by
  their own functions, moved out of `search_query_batch` into the modules of the subject loop, the
  cutoffs, the linking calls and the report hits. No output changes.
