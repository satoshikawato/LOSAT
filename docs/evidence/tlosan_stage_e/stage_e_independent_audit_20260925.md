# Independent read-only TBLASTN Stage E audit — 2026-09-25

Auditor: `ncbi_parity_auditor` subagent, independent read-only review. Verdict: **SUPPORTED** for the declared native Linux x86_64 single-thread local `-subject` fixtures and supported option boundary, against fixed NCBI source `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4` and BLAST+ 2.17.0+.

- Final LOSAT release binary SHA-256: `da85b0a87a0bc92fcfa1bcc92fb8aa76716424286468a60a92994d2ddfca5e87`.
- Reviewed complete outfmt 0/6/7 byte comparisons: physical fixtures **219/219**, generated fixtures **27/27**, all 27 subject genetic codes **81/81**; code-1 CLI/API same-input calibration **81/81**.
- Independently reran uneven-gap outfmt 0, natural heap-order outfmt 7, and code-32 outfmt 0 cases; checked their commands and recorded output SHA-256.
- Confirmed Stage D oracle replay **130/130**, checksum gate **615** files, focused tests **106/106**, and explicit unsupported-path rejection **14/14**.
- Checked NCBI source and Rust timing, values, ordering, masks, frames, and absence of NCBI runtime/build/fallback dependencies. Found no remaining concrete mismatch.

The auditor's earlier findings were corrected before the verdict: final binary SHA in every comparison row, inclusion of eight per-fixture manifests including uneven-gap, and a pinned NCBI source snippet above the `Expect(n)` link-count branch. The verdict does not cover Wasm, multithreading, performance, or arbitrary untested option combinations.
