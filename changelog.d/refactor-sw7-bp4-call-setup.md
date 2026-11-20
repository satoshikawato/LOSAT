- Engine: the setup of a BLASTP call before its search (the reports that its formats ask
  for, the queries' IDs, titles and lengths, the query frames, the encoded subjects and
  their titles, the scoring values and X-drops, and the query batches, the lookup table
  and the query contexts) moves out of `run_call` into its own module. No output changes.
