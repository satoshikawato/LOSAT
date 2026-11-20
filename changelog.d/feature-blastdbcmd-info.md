- `losat blastdbcmd`: NCBI BLAST+ 2.17.0 blastdbcmd's `-info` (with `-db`, `-dbtype`, `-exact_length`, `-out`,
  `-logfile`), `-h`, `-help` and `-version`, byte-identical to NCBI (the text, the database errors with exit status 2,
  and NCBI's usage report for every argument error). Its command line is NCBI's own. Sequence retrieval (`-entry`,
  `-entry_batch`, `-outfmt`), `-list`, the taxonomy options, `-metadata`, a database with the masks of
  `makeblastdb -mask_data` and a `-db` value that names several databases are refused.
