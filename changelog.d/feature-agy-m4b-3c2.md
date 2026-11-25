- `tblastn` and `tblastx` with `-db_hard_mask` now search databases whose masks overlap or nest when
  they are stored in the order of their start, with NCBI BLAST+ 2.17.0's `-num_threads 1` rows at
  every `-num_threads` (and searching any sequence whose word hits would reach another sequence's
  diagonal entries on its own without those entries, as `blastp` and `blastx` do under
  `PD-LOSAT-NCBI-DEFECTS` approved exception 29); they refused every such database before. Masks
  stored out of the order of their start are still refused, as are the rare sequences listed in the
  manual's refusals page (line NP, M4b-3c).
