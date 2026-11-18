- LOSAT Web and the WebAssembly (wasm32) builds now do what the command line does where they differed:
  a BLASTX batch whose lookup table would need over 268 million entries is built by the reference
  builder (same output) instead of stopping with "capacity overflow"; rpsbproc prints NCBI's text
  for a `cdtrack.txt` line that ends early (`__pos (which is 18446744073709551615)`, not
  `4294967295`); a `-window_masker_taxid` run without a counts file reports NCBI's "could not open "
  where WASI used to open its root directory; and on WASI hosts `-export_search_strategy X -out X`
  and `-out_pssm X -out X` are refused as on Linux (a symbolic or hard link to the `-out` file is
  not seen there, because Rust's standard library has no stable file identity on WASI).
