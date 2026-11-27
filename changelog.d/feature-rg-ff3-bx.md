- `losat blast_formatter` reads BLAST archives in ASN.1 binary and XML too (`-outfmt 11` with
  `ARCHIVE_FORMAT=bin` or `xml`, from LOSAT or NCBI BLAST+), which it refused before: the reports,
  messages and exit statuses are NCBI BLAST+ 2.17.0's, the same as for the archive's ASN.1 text, also
  for archives of several query batches, from the standard input, and for a cut archive or bytes
  after it.
