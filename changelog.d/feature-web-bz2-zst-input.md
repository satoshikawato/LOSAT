- Web: FASTA files named `.zst` opened in "Or, upload file" are inflated in the browser, as `.gz`
  files are, and kept under their name without the suffix; every zstd frame of a concatenated file
  is read.
