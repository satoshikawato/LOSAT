- Database searches: a sequence with a definition line without a Seq-id beside lines with Seq-ids (a
  hand-edited database; makeblastdb writes none) is reported as NCBI BLAST+ 2.17.0 reports it instead of
  the run being refused: `Unknown` in the description line where NCBI's walk of the lines meets that line
  (up to 4096 bytes), the skipped id in place of each alignment's heading, outfmt 5's `Unknown` hit where
  that line has a title after the first titled line, the XML2 and JSON descriptions up to that line, and
  `Unknown`, 0 and `N/A` in the tabular fields. LOSAT refuses `sallacc` of such a sequence (NCBI stops with
  a null pointer, exit 255), and such a sequence in the organism report (outfmt 18) and in `rpsblast` and
  `rpstblastn`.
- Database searches: where LOSAT refuses a tabular row at which NCBI BLAST+ stops or prints an uninitialised
  number (`sallacc` above; `sacc`, `saccver` or `staxid` of a sequence whose lookup fails), the reports of
  the query batches before it are written first, as NCBI writes them, and then the messages that NCBI
  writes of the failing batch before that row (its reader messages, and for its queries up to the row's
  their warnings, the taxonomy database warning of a taxonomy name field and, with the prefetch off, the
  messages of reading the subjects up to the row); `blastp`, `psiblast`, `tblastn` and `tblastx` wrote no
  report before the refusal, nor did `blastn` where it writes its reports after the search (with
  `-dbsize` near NCBI's 64-bit bound or a negative `CHUNK_SIZE`), and no program wrote all of those
  messages.
- Database searches: the outfmt 0 65-byte description-id length check measures the Seq-align id written
  by NCBI's catch path whenever generating the description line throws (a definition line without a
  Seq-id, or an undecodable title or source name) for both protein and nucleotide databases.
