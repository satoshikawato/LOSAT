- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where `PRE_FETCH_SEQS_LIMIT`
  turns a query batch's prefetch off and a query's report reads a parsed Seq-id that a database of the list cannot
  look up: the pairwise report (`-outfmt 0`) and the query-anchored reports (`-outfmt 1` to `4`) now give the results
  of NCBI BLAST+ instead of LOSAT's refusal, in BLASTP, psiblast (`-outfmt 0`, each round's report), BLASTN (every
  task), BLASTX, TBLASTN and TBLASTX. The description line of such a sequence shows its Seq-align id and `Unknown`
  (`NP_100001.1 Unknown`, with `-show_gis` `ref|NP_100001.1|Unknown`) with its scores. `-outfmt 0` stops after the
  query's description lines with `LOSAT database error: DB contains no accession info.` (exit status 2, NCBI
  BLAST+'s `BLAST Database error:` in LOSAT's words) where the
  first alignment it shows is of such a sequence under a parsed Seq-id (not a local one), and otherwise writes
  `Sequence with id NP_100001 no longer exists in database...alignment skipped` for each of its alignments and goes
  on; `-outfmt 1` to `4` stop there where any alignment they show is of such a sequence. The taxonomy database's texts
  are those of the sequences that the report reads before it stops (in psiblast's rounds after the first, those of
  the two-part description table, which counts the first `-num_descriptions` HSPs, not sequences, and shows at least
  the first; psiblast ignores `-sorthits` and `-sorthsps` there, as NCBI BLAST+ does). A run that such a report stops
  ends there, as NCBI BLAST+'s does: no later query batch is prefetched. LOSAT's library still refuses such a sequence
  in a call that writes `-outfmt 0` to `4` with another format that reads it and stops at another place (`-outfmt 0`
  with `1` to `4`, or either with a table, XML or JSON format), for the whole call, also where `-num_descriptions` or
  `-num_alignments` make it search once for each hit list size.
