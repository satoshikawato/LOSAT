- Input files named `.bz2` and `.zst` are read through bzip2 and zstd, as `.gz` files are through
  gzip: `-query` and `-subject` of every search program and `makeblastdb -in`, concatenated streams
  and frames (pbzip2, pzstd) included; makeblastdb's refusal of `.bz2` / `.zst` names is gone. NCBI's
  makeblastdb reads only the first bzip2 stream or zstd frame of such a file (`PD-LOSAT-NCBI-DEFECTS`
  exception 7). Pure-Rust decoders (bzip2 with libbz2-rs-sys, ruzstd).
