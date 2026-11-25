- BLAST database alias files with `SEQIDLIST` and `TAXIDLIST` over version 5 databases are now
  read through LMDB in all search programs, `blastdbcmd` and `blastdb_aliastool` (removing the
  version 5 list refusal; only `MEMB_BIT` remains refused). Alias `SEQIDLIST` (text and binary
  `.bsl`) and `TAXIDLIST` filter OIDs and definition lines as in NCBI BLAST+ 2.17.0, combine with
  Pass 1 alias filters and command-line `-seqidlist`, `-negative_seqidlist`, `-taxids` and
  `-negative_taxids` filters, and preserve NCBI's eager and lazy error messages and exit statuses.
