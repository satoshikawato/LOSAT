- `losat blastdbcmd -ipg` / `-ipg_batch` write the sequence of a PIG whose OID an alias file's
  `FIRST_OID` / `LAST_OID` range leaves out, as NCBI BLAST+ 2.17.0 does (LOSAT refused it); `%o` is
  that OID.
- An alias file's `LENGTH` or `STATS_TOTLEN` over volumes without sequences is searched as NCBI
  BLAST+ does (no hit, each query's `Filtering resulted in an empty database.` where the alias file
  gives no sequence count; LOSAT refused it). BLASTN and TBLASTX still refuse it without `NSEQ` or
  `STATS_NSEQ`, where NCBI BLAST+ divides by 0 sequences and crashes.
- The search programs report a database's errors in NCBI's order: a list's error and an alias file's
  `NSEQ`/`STATS_NSEQ`/`LENGTH`/`MIN_SEQ_LENGTH` error come before a damaged first header.
- Alias files whose `FIRST_OID` or `LAST_OID` is 2^31 or more are read as NCBI BLAST+ 2.17.0 reads
  them, as negative numbers (LOSAT refused them), and the `FIRST_OID`-`LAST_OID` range of the only
  alias file of a `-db` that also names another database is searched as NCBI does (LOSAT refused
  it): the searches read the sequences inside that range, a list's search reads every listed
  sequence with NCBI's statistics, and `losat blastdbcmd -entry all`, `-entry`, `-entry_batch`,
  `-taxids`, `-tax_info` and `-metadata` give the sequences and taxids NCBI's blastdbcmd gives.
