- Engine: a BLASTN query batch's preliminary-stage search spaces, its cutoffs and the statistics
  of one subject are made by their own functions, moved out of `search_query_batch` into the
  module of the batch's other set-up steps. No output changes.
