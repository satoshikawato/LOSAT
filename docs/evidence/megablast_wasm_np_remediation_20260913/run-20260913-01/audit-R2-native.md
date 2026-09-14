# Independent read-only R2 audit

The ncbi_parity_auditor checked the C++ zero-initialized options, traditional megablast defaults (no SetXDropoff), blastn and dc-megablast defaults (20 bits), and BlastInitialWordParametersUpdate (zero initialization selects the current word cutoff). The task-specific Rust branch is supported. No blocking source finding.

The auditor directly compared all three B1 native n1/outfmt6 raw pairs: NZ 38,401 bytes, EDL933/Sakai 473,032 bytes, Sakai/MG1655 536,518 bytes. Each is exact and each NCBI output matches the corresponding B0 oracle hash. B1 native executable: ae494b43e50ad15ae004d23ded13b325cfe8c3bc3b036928ca3302d2d2adf922.

At this pass, native n8, Wasm, task/format/boundary regressions, frozen Gate A/B and performance adoption were not audited. This pass does not certify those gates. No build, test, or source edit was delegated to the auditor.
