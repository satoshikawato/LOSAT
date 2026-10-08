---
name: verify-ncbi-parity-and-speed
description: Implement, diagnose, or verify LOSAT behavior against the NCBI BLAST source and executable oracle while preserving performance. Use for LOSAT TBLASTX, BLASTN, BLASTP, TBLASTN, BLASTX, native or Wasm parity, output differences, hit-count discrepancies, coordinate or score mismatches, NCBI source ports, benchmarks, SIMD, threading, or release parity evidence.
---

# Verify NCBI parity and speed (Claude Code entry)

The skill body is shared with Codex and lives in
`.agents/skills/verify-ncbi-parity-and-speed/`. Read
`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md` now and follow it, with its
`references/` files.

Claude Code names for what it mentions:

- NCBI source tree: `$NCBI_SRC`; NCBI BLAST+ executables: `$NCBI_BIN` (see `CLAUDE.local.md`).
- The `ncbi_parity_auditor` custom agent: the `losat-reviewer` agent, asked for an independent
  audit of one angle.
- Running the oracle (concurrency, lock, work directories): skill `losat-oracle-runs`.
- Which checks to run for a change: skill `losat-gates`.
