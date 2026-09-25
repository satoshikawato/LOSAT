# Independent read-only audit — Stage G, 2026-09-26

Auditor: `ncbi_parity_auditor` custom agent. Verdict: **SUPPORTED for the declared, fixture-scoped Stage G TBLASTN parity contract and fixed-build absolute performance record**. The auditor's two final conditions were regeneration/check of `evidence.sha256` after final edits, and committing the audited source and evidence.

The auditor independently checked the pinned NCBI source and BLAST+ binary, accepted `PD-TLOSAN-LOCAL-GENCODE-32`, exact 162 unique mandatory rows, expanded native/real/no-hit/batch-boundary/Wasm rows, raw code-32 0/6/7 samples, final build/test/fmt/clippy/negative gates, BLASTN 14/14, BLASTP 9/9, focused TBLASTX 12/12, and recomputed six benchmark medians/ranges. It confirmed that NCBI's prebuilt `-db` benchmark output is separate from LOSAT local `-subject` distribution output.

The verdict does **not** cover arbitrary untested inputs/options, a LOSAT-vs-NCBI speedup, release-wide parity, or the separate broad TBLASTX audit. A later non-code-only difference in any relevant fixture reopens the gate.
