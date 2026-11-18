- WASI (`wasm32-wasip1`, `wasm32-wasip1-threads`): a value in psiblast's ASN.1 PSSM file (`-out_pssm`
  with `-out_ascii_pssm`) can differ from the other platforms' by 1 in its last (15th) digit, because
  WebAssembly's `log()` has no FMA; NCBI BLAST+ 2.17.0 writes the same digit on an x86-64 CPU without
  FMA. PD-LOSAT-PLATFORM-LIBM 1.1 records this for WASI (Owner, 2026-11-18), as 1.0 did for macOS. The
  WASI checks (`LOSAT/tests/check_wasm_threading.py`) now run every psiblast `-db` row through each WASI
  command, with that record's rule on the PSSM files only (`db_search_fixtures.py check
  --libm-last-digit`); the rule, in `xinclude_fixtures.py` and in the macOS `-outfmt 11` test, accepts
  only a 15-digit mantissa 1 away.
