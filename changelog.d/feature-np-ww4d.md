- LOSAT Web: a local BLAST database whose path holds a space (a folder such as `my dbs`) is now searched.
  Search's Database (`-db`), BLASTN's Repeats database (`-filtering_db`) and WindowMasker's `-infmt
  blastdb` input write its name between double quotes (`-db '"my dbs/ndb"'`), as the command line reads
  it; before, the search looked for a database named `my` and failed. Run details shows the command
  with the name quoted, so it works when pasted. A database whose name holds a double quote, which
  the command line cannot read in these options, is refused where it is chosen, with the reason.
