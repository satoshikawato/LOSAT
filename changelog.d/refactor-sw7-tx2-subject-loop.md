- Engine: the seed tables, subject loop and hit collection of a TBLASTX query batch move out of
  `search_query_batch`. No output changes; `LOSAT_TIMING` drops its `identity_calc` line (always 0).
