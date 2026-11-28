- `tblastn` and `tblastx` `-db`: a nucleotide database sequence of negative length (an ambiguity
  offset at or below its start, inside the sequence file) is searched as NCBI BLAST+ 2.17.0
  searches it, instead of the whole search being refused: at `-num_threads 1` (and with `tblastn`'s
  query split, `-mt_mode 1` or its automatic choice) it is skipped, and where the threads split the
  database the search stops before it, as `blastn`, `blastp` and `blastx` do, also before a later
  volume whose sequence file is empty. Ambiguity data that cannot be read (`amb > next`, or fewer
  entries than their count) in a sequence the search reads remain refused.
