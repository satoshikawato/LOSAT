- LOSAT Web: a local BLAST database whose path holds a space (a folder such as `my dbs`) is now
  searched. Search's Database (`-db`), BLASTN's Repeats database (`-filtering_db`) and the `-infmt
  blastdb` input of WindowMasker and of the Masking tab's DustMasker and SegMasker write its name
  between double quotes (`-db '"my dbs/ndb"'`), as the command line reads it; before, the search looked
  for a database named `my` and failed. Run details shows the command with the name quoted, so it works
  when pasted. A database whose name holds a double quote, which the command line cannot read in these
  options, is refused where it is chosen, with the reason. The refusal stays in the input block, the
  check and Run LOSAT for as long as that database is chosen or its files are in the run. With
  WindowMasker's `-infmt seqids` or `-fa_list` and the Masking tools' FASTA input, whose input is one
  file's name, the files of the run do not fill in the input name. The Database picker reads the
  database list of an alias file as BLAST+ does, so a multi-volume database whose volume names hold a
  wide space (U+3000) is one choice, not several.
