- `exp`, `log`, `log10` and `pow` now give the bits of glibc's x86-64 FMA versions, which NCBI BLAST+
  2.17.0 used for the reference output, on every platform LOSAT builds for (`LOSAT/src/utils/glibc_fma_libm.rs`).
  On WebAssembly (WASI and the web page) and on x86-64 CPUs without FMA, a value printed with all 15
  digits no longer differs from Linux by 1 in its last digit: the E values of `-outfmt 11` and
  `-outfmt 8` (for example `tblastx`) and the values of
  psiblast's ASN.1 PSSM. Tests: the WASI checks (`LOSAT/tests/check_wasm_threading.py`) run every
  psiblast `-db` row through each WASI command and compare the output files, the ASN.1 PSSM included,
  byte for byte with the frozen bytes of NCBI BLAST+ 2.17.0; `clippy` refuses the transcendental
  methods of `f64` outside that module.
