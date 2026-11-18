- `blastn` with a `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` that leave a query chunk without a query (an evened chunk size
  below the overlap, or a negative `CHUNK_SIZE` above a negative `OVERLAP_CHUNK_SIZE` whose first piece covers the
  query batch) now stops where NCBI BLAST+ 2.17.0 does, before any chunk is searched, and writes the results of the
  query batches before it first in every format (before, only the formats written batch by batch had them, and the
  negative pair was refused before any output). The overlap is no longer an unsigned subtraction in debug builds.
  Searches without such a pair are unchanged.
