- Web: FASTA files named `.bz2` or `.zst` opened in "Or, upload file" are inflated in the browser,
  as `.gz` files are, and kept under their name without the suffix; every bzip2 stream and zstd
  frame of a concatenated file is read.
