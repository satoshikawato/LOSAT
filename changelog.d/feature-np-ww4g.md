- LOSAT Web reads a `-db` database's files in far fewer pieces: each engine worker keeps the last 8 MiB
  that it read from the user's files, so the index's offsets and the sequences that a search fetches
  again are not read again, and each place of a file read in order (a database's index is read in two
  places in turn; a protein database's letters are read in blocks that begin one byte before the last
  one ended) is read ahead. `blastp` of 5 proteins against 600,000 proteins in Chromium read the files
  95 times instead of 333 (1.6 s instead of 2.6 s), and 607 instead of 1,330 times with 4 threads.
