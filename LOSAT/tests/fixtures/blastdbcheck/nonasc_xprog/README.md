# Non-ascending index offset cross-program fixtures (PR #1660 R2)

Frozen outputs of official NCBI BLAST+ 2.17.0 (`blastdbcmd`, `blastp`, `blastn`, `blastdbcheck`)
on databases with non-ascending index offset arrays, tested in `tests/blastdbcheck_fixtures.rs`
(`nonasc_cross_program_parity`) under `[profile.test]` (integer overflow checks enabled).

## Variants
- `d_p4_nonasc_seq`: `p4` with `p4.pin` offset 106 poked to `01f4`, setting
  `seq_offsets = [1, 62, 500, 179, ...]` (OID 1 spans 437 residues across the rest of `.psq`
  with zero-padded tail; OID 2 has signed length `-322`).
  - `nonasc_seq_dbcmd_fasta`, `nonasc_seq_dbcmd_fmt`, `nonasc_seq_blastp_t1`: LOSAT matches
    NCBI BLAST+ 2.17.0 byte-for-byte.
  - `nonasc_seq_blastp_t4`: frozen NCBI `-num_threads 4` output is kept for provenance. In NCBI,
    `CSeqDBImpl::x_FillSeqBuffer` (`seqdbimpl.cpp:699-707`) returns an empty chunk at OID 2
    (`length = -322 < 0`), causing `GetNextOIDChunk` to stop early and drop OIDs 3..8. LOSAT
    matches NCBI `-num_threads 1` (`nonasc_seq_blastp_t1`) at both `-num_threads 1` and
    `-num_threads 4` (cross-thread invariance, without debug-build subtraction underflow).
- `d_p4_nonasc_hdr`: `p4` with `p4.pin` offset 68 poked to `00000000`, setting
  `hdr_offsets = [0, 79, 167, 0, 334, ...]` (OID 2 has descending header range `[167, 0)`;
  OID 3 has range `[0, 334)` starting at OID 0's header).
  - `nonasc_hdr_dbcmd_fasta`, `nonasc_hdr_dbcmd_fmt`, `nonasc_hdr_blastp_q2_t1` (hits OID 2,
    silent `CEofException` end), `nonasc_hdr_blastp_q4_t1` (hits OID 4): LOSAT matches NCBI
    byte-for-byte at `-num_threads 1` and `-num_threads 4`.
  - `nonasc_hdr_blastp_q3_t1`: frozen NCBI output is kept for provenance (OID 3's header decodes
    OID 0's `gi|101` / `NP_000101.1`, so NCBI formats OID 3 with OID 0's letters at 0.000%
    identity). LOSAT refuses OID 3 under the existing parsed-id lookup-mismatch refusal.
- `n4_amb_eq_start`: `n4` with `n4.nin` offset 100 poked to `00000020`, setting
  `amb[1] = seq[1] = 32` (breaking `start < amb`). NCBI `x_GetSequence` (`seqdbvol.cpp:1833-1835`)
  indexes `(*buffer)[-1]` out of bounds before OID 1, and `blastn -num_threads 4` gives a
  different hit set from `-num_threads 1`. LOSAT refuses at the shared `seqdb` read site across
  `blastdbcmd`, `blastn` (`-num_threads 1` and `4`), and `blastdbcheck`.
- `n4_amb_gt_next`: `n4` with `n4.nin` offset 100 poked to `0000003b`, setting
  `amb[1] = seq[2] + 4 = 59` (breaking `amb <= next`). NCBI reads 4 bytes past OID 1 into OID 2
  in `GetSeqLengthExact`/`x_GetSequence`, and `blastn -num_threads 4` aborts with `SIGABRT`
  (`rc=-6`, uncaught `ncbi::CSeqDBException`). LOSAT refuses at the shared `seqdb` read site
  across `blastdbcmd`, `blastn` (`-num_threads 1` and `4`), and `blastdbcheck`.
