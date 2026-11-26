# Non-ascending index offset cross-program fixtures

Frozen outputs of official NCBI BLAST+ 2.17.0 (`blastdbcmd`, `blastp`, `blastx`, `psiblast`,
`blastn`, `blastdbcheck`) on databases with non-ascending index offset arrays, tested in
`tests/blastdbcheck_fixtures.rs` (`nonasc_cross_program_parity`) under `[profile.test]`
(integer overflow checks enabled).

## Variants
- `d_p4_nonasc_seq`: `p4` with `p4.pin` offset 106 poked to `01f4`, setting
  `seq_offsets = [1, 62, 500, 179, ...]` (OID 1 spans 437 residues across the rest of `.psq`
  with zero-padded tail; OID 2 has signed length `-322`).
  - `nonasc_seq_dbcmd_fasta`, `nonasc_seq_dbcmd_fmt`, `nonasc_seq_blastp_t1`,
    `nonasc_seq_blastp_t2`, `nonasc_seq_blastp_t4`, `nonasc_seq_blastp_mt1_t4`,
    `nonasc_seq_blastx_t1`, `nonasc_seq_blastx_t4`, `nonasc_seq_blastx_mt1_t4`,
    `nonasc_seq_psiblast_t1`, `nonasc_seq_psiblast_t4`: LOSAT matches NCBI BLAST+ 2.17.0
    byte-for-byte. In multithreaded DB-split mode (`-num_threads >= 2` without `-mt_mode 1`
    query-split), NCBI's `CSeqDBImpl::x_FillSeqBuffer` (`seqdbimpl.cpp:699-707`) returns an empty
    chunk at OID 2 (`length = -322 < 0`), stopping the scan before OIDs 3..8; `-num_threads 1`
    and `-mt_mode 1 -num_threads 4` scan all OIDs.
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
  `amb[1] = seq[1] = 32` (NCBI reads file byte 31 inside its mapping, yielding signed length
  `-4`). `n4_ambeq_dbcmd`, `n4_ambeq_dbcmd_fmt`, `n4_ambeq_blastn_t1`, `n4_ambeq_blastn_t4`
  (stops at OID 1 in DB-split mode), `n4_ambeq_blastn_mt1_t4`, `n4_ambeq_check`, and
  `n4_ambeq_check_v3` all match NCBI BLAST+ 2.17.0 byte-for-byte.
- `n4_amb_gt_next`: `n4` with `n4.nin` offset 100 poked to `0000003b`, setting
  `amb[1] = seq[2] + 4 = 59` (length 107 inside the mapping, then `x_GetAmbChar` throws
  `File error: could not get ambiguity data.`). `n4_ambgt_dbcmd`, `n4_ambgt_dbcmd_fmt`,
  `n4_ambgt_blastn_t1`, `n4_ambgt_check`, and `n4_ambgt_check_v3` match NCBI BLAST+ 2.17.0
  byte-for-byte. `n4_ambgt_blastn_t4` (and `-mt_mode 1 -num_threads 4`) aborts in NCBI
  (`SIGABRT`, `rc=-6`) and is refused by LOSAT.
- `n4_amb0_eq_0`: `n4` with `n4.nin` offset 96 poked to `00000000`, setting `amb[0] = 0`
  (NCBI reads byte `-1` before its sequence file mapping and crashes with `SIGSEGV`). Refused
  across `blastdbcmd`, `blastn` (`-num_threads 1` and `4`), and `blastdbcheck`.
