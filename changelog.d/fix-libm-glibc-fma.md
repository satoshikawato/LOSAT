- `exp`, `log`, `log10` and `pow` now give the bits of glibc's x86-64 FMA versions, which NCBI BLAST+
  2.17.0 used for the reference output, on every platform LOSAT builds for (`LOSAT/src/utils/glibc_fma_libm.rs`).
  On WebAssembly (WASI and the web page) and on x86-64 CPUs without FMA, a value printed with all 15
  digits no longer differs from Linux by 1 in its last digit: the E values of `-outfmt 11` and
  `-outfmt 8` (for example `tblastx`, `blastp` with composition-based statistics) and the values of
  psiblast's ASN.1 PSSM. Tests: the WASI checks compare psiblast's PSSM files byte for byte; `clippy`
  refuses the transcendental methods of `f64` outside that module.
