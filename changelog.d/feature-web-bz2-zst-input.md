- Web: FASTA files named `.bz2` or `.zst` opened in "Or, upload file" are inflated in the browser,
  as `.gz` files are, and kept under their name without the suffix; every zstd frame and every
  bzip2 stream of a concatenated file is read.
