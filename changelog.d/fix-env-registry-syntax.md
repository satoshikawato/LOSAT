- A configuration file (`.ncbirc`, `<program>.ini`) is read as NCBI BLAST+ 2.17.0 reads it, escapes, quotes,
  continued lines and byte order mark included, and one that NCBI BLAST+ reports as malformed no longer stops
  LOSAT: a malformed `.ncbirc` gets NCBI's `Critical:` message for each of its two reads (the first before
  `-version`) and the program goes on with the entries above the bad line; a malformed `<program>.ini` ends
  the program with NCBI's messages and exit status 2, `-h` included. With `-logfile`, blastdbcmd and
  makeblastdb write the messages where NCBI does (makeblastdb `-version -logfile F` now makes `F`, as NCBI
  does). A file found that cannot be opened is no file, as in NCBI, instead of stopping LOSAT. On Windows
  the system-wide file is `ncbi.ini`, as in NCBI. LOSAT still stops for a file in UTF-16 or a malformed line
  that is not UTF-8, and its library refuses a malformed file.
