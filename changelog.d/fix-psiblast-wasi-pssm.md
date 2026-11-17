- WASI command builds (`wasm32-wasip1`, `wasm32-wasip1-threads`): psiblast with a PSSM as the query
  (`-in_pssm FILE`) now runs, as natively, instead of stopping with "reading the query from standard
  input ('-') on this platform is not supported by LOSAT's PSIBLAST" (NCBI reads no query then; LOSAT
  had opened the default `-query -`, whose standard input WASI cannot copy), and its own errors come in
  NCBI's order; `-out_pssm` (the ASN.1 PSSM file) is written instead of the command panicking
  (`std::process::id()`, which WASI does not have, named the staging file). The serial and threaded
  command checks (`LOSAT/tests/check_wasm_threading.py`) now run every psiblast row of
  `outfmt0_manifest.tsv` through each WASI command.
