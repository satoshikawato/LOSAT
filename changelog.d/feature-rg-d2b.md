- New program `blastdbcheck`, NCBI BLAST+ 2.17.0's check of BLAST databases: one database (`-db`) or every database
  of a directory (`-dir`, `-recursive`), its metadata, files, sequences, ISAM indexes, alias files and taxonomy ids
  (`-must_have_taxids`), every sequence (`-full`), every N-th (`-stride`), a random sample (`-random`, by default
  200) or both ends (`-ends`), at `-verbosity` 0 to 4 and to `-logfile`. It writes NCBI's report with NCBI's exit
  statuses. `-cdd_delta` is refused until it is ported, and so are cut database files that NCBI reads past their end
  and alias files whose counts overflow NCBI's arithmetic where the overflow decides the check.
