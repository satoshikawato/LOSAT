# Non-ascending index offset cross-program fixtures

Frozen outputs of official NCBI BLAST+ 2.17.0 (`blastdbcmd`, `blastp`, `blastx`, `psiblast`,
`blastn`, `blastdbcheck`) on databases with non-ascending index offset arrays, tested in
`tests/blastdbcheck_fixtures.rs` (`nonasc_cross_program_parity`) under `[profile.test]`
(integer overflow checks enabled).

How the rows of `n6db`, `mv3db`, `p2db` and `fdb` were made: the test of each row sets up its
run directory (the database copied, the poke, the alias, the list) and names the program and its
arguments; NCBI's program ran with the same arguments in a directory set up the same way, with
the environment of the test's LOSAT run (`PATH=/usr/bin:/bin`, `TZ=UTC`, `LC_ALL=C`, `HOME` and
`NCBI_CONFIG_PATH` an empty directory, `BLAST_USAGE_REPORT=false`), and its standard output,
standard error and exit status (a signal as `-N`) are the row's `.out`, `.err` and `.rc`.

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
  - `nonasc_seq_blastp_neggi2_t4`, `nonasc_seq_blastp_neggi2_t1`, `nonasc_seq_blastx_neggi2_t4`,
    `nonasc_seq_psiblast_neggi2_t4` (`-negative_gilist` of GIs 103 and 105: OIDs 2 and 4
    excluded, so the included OIDs form three blocks, 0..2, 3 and 5..8): LOSAT matches NCBI
    BLAST+ 2.17.0 byte-for-byte. At 4 threads the excluded OID 2 before the block of OID 3 stops
    the scan, and the block of OIDs 5..8 after it, whose own excluded OID 4 does not stop, is
    read but dropped (frozen 2026-11-27 for review R5 F3 of PR #1660).
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
- `n4_amb_past_end`: `n4` with `n4.nin` offset 108 poked to `10000000`, setting `amb[3] = 2^28`
  (OID 3's last byte lies 256 MiB past the 114-byte `.nsq`; NCBI reads it and crashes with
  `SIGSEGV`, `rc=-11`). `n4_ambpast_blastn_neggi_t1`, `n4_ambpast_blastn_neggi_mt1_t4`
  (`-negative_gilist` of GI 304, OID 3 excluded: NCBI never reads OID 3 at 1 thread or with
  `-mt_mode 1`), `n4_ambpast_blastn_negall_t1` (every OID excluded: `Filtering resulted in an
  empty database.`), `n4_ambpast_blastn_range3_t1` (`-db n4r`, an alias file whose only filter is
  `FIRST_OID 1`, `LAST_OID 3`: OIDs 0..2, iteration range end 3) and `n4_ambpast_blastn_range2_t4`
  (`-db n4r2`, `LAST_OID 2`: NCBI reads OID 2 at the range's end, not OID 3) match NCBI BLAST+
  2.17.0 byte-for-byte. `n4_ambpast_blastn_t1`, `n4_ambpast_blastn_t4`,
  `n4_ambpast_blastn_neggi_t4`, `n4_ambpast_blastn_negall_t4` and `n4_ambpast_blastn_range3_t4`
  crash in NCBI and are refused by LOSAT with the same message (`byte range
  268435455..268435456 passes the end of the file`): in multithreaded DB-split mode NCBI's
  `x_FillSeqBuffer` (`seqdbimpl.cpp:684-710`) reads the excluded OIDs after the last included
  one, up to the iteration range's end, every OID of a list that includes none, and the OID at
  the range's end when it is inside a volume (it reads the next OID before it compares it with
  the end) (frozen 2026-11-27 for review R5 F1 of PR #1660).
- `n6db` (PR #1660 review R6; test `nonasc_gap_stop_before_crash_parity`): six sequences made
  with NCBI `makeblastdb` (v4, GIs 601..606), `n6.nin` poked so that OID 1 has `amb == start`
  (signed length `-4`) and OID 5 has `amb = 2^28` (last byte far past the file).
  `n6_<row>_t<T>.{out,err,rc}` are NCBI BLAST+ 2.17.0 `blastn -db n6 -outfmt 6` at
  `-num_threads` 1, 2 and 4, without a list (`n6_t<T>`) or with the `-negative_gilist` of GIs
  602, 604 and 606 (`neg135`), 604 and 606 (`neg35`) or 606 (`neg5`): `t1` crashes (`rc=-11`,
  refused by LOSAT); at 2 and 4 threads only OID 0's rows (the scan stops at OID 1); the other
  rows complete and LOSAT matches them byte-for-byte (frozen 2026-11-27, RG-D2 tail item 4; the
  2-thread rows and `neg35_t1` frozen 2026-11-28 for review R4 finding 4).
- `mv3db` (RG-D2 tail item 3, 2026-11-27): three volumes `v0` (5 OIDs), `v1` (6), `v2` (4), GIs
  601..615, made with NCBI `makeblastdb`. `mv3e_rows.tsv` and `mv3e_o*.out`, `mv3e_e6.err` (test
  `nonasc_empty_later_volume_sequence_file_parity`): `v1.nsq` of 0 bytes (`empty`), also OID 2 of
  `v0` with signed length `-4` (`neg`), also the first sequence offset of `v1` 0 (`zero`,
  `negzero`); aliases `n`, `r_head`, `r_tail`, `r_pre`; lists none/GI lists; 1, 2 and 4 threads.
  NCBI's read of an OID of the empty file that starts at 0 gives the length -1, which stops the
  multithreaded scan and which one thread skips, and it crashes (`SIGSEGV`) at any other start.
  The rows of the lists `posv1first` (`-gilist` of GI 606, the first OID of `v1`) and
  `posv1first2` (606 and 607), kinds `zero` and `empty` (frozen 2026-11-28 for review R4 finding
  3(a)): with `zero`, `posv1first` completes without rows at every thread count, and
  `posv1first2` crashes at one thread (OID 6 is read) and completes at 2 and 4. LOSAT matches all
  204 rows (59 refused, 145 identical). Each crash is NCBI's
  first read of a sequence of `v1`, inside the iteration range of its row (`r_head`'s range ends
  at OID 7, inside `v1`, after OIDs 5 and 6; with a `-gilist` or `-negative_gilist` NCBI's range is
  the whole database), so no row depends on NCBI's reads past the end of a range.
- `p2db` (RG-D2 tail R2, review finding 8, 2026-11-28): two protein volumes `pv0` (GIs
  7001-7003) and `pv1` (7004-7006), made with NCBI `makeblastdb`; the test
  `nonasc_protein_scan_stop_parity` adds the alias `p` (`DBLIST pv0 pv1`) and a GI list `gl.txt`.
  `p2_<case>_t<T>.{out,err,rc}`: NCBI `blastp -db p -outfmt 6 -gilist gl.txt` at 1 and 2 threads
  (`p2x_`: `blastx` with the query of the codons of the protein query). A protein of length -1
  has its end (the next start) set to its start. `gapstop`: OID 1 of `pv0` of length -1, left out
  by the list before the block of `pv1`; with `pv1.psq` empty or missing, NCBI's 2-thread scan
  stops before `pv1` and completes. `v1gap`: the negative OID is OID 1 of `pv1`, after its first
  OID (a missing file is NCBI's memory map error). `emptyv1_skipfirst`: no stop before the empty
  `pv1.psq` (`SIGSEGV`). Review R4 finding 1 (frozen 2026-11-28, blastp and
  blastx): `trail_missingv1` (list 7001, `pv1.psq` missing): the 2-thread scan reads the excluded
  OIDs after OID 0 and stops with NCBI's memory map error (`rc=3`) where one thread completes;
  `trail_neg_missingv1` (and OID 1 of `pv0` of length -1, which stops that scan first) and
  `trail_emptyv1` (`pv1.psq` empty: a protein's length comes from the index) complete;
  `none_missingv1` (list 9999, no OID included): one thread warns that the filtering left an
  empty database, two threads give the memory map error without the warning;
  `none_neg_missingv1`: the warning at both. Review R4 finding 2 (frozen 2026-11-28):
  `v1firstneg_emptyv1` (no list; OID 0 of `pv1` of length -1, `pv1.psq` empty): the 2-thread scan
  stops at that OID (`pv0`'s rows), one thread skips it and crashes at OID 4 (`SIGSEGV`);
  `v1firstneg_emptyv1_list` (list 7001 and 7004): nothing is read from the empty file (blastp:
  OID 0's row; blastx: `BLAST engine error: Unknown error code 1`, `rc=3`, at both thread
  counts). LOSAT matches every row (the crashes refused).
- `fdb` (RG-D2 tail R4 item 5, 2026-11-28; test `nonasc_filtering_db_repeat_search_runs_one_thread`):
  NCBI `makeblastdb` v4 databases `main` (`s0` holds R,
  the first 150 bases of the query `fdb_q.fna`; `s1` holds U, its last 150), and the repeat
  databases `rep` (`r2` holds R) and `repa` (`r1` holds R). `fdb_<case>_t<T>.{out,err,rc}`: NCBI
  `blastn -db main -filtering_db rep|repa -outfmt 6` at 1 and 2 threads. `clean`: R is masked, only
  `s1` is reported. `neg`: OID 1 of `rep` with `amb == start` (signed length -3): NCBI's search of
  the repeats runs one thread, reads past it and still masks R (the same rows). `ambnext`: OID 1 of
  `repa` with its ambiguity offset 4 bytes past the next start: the search of the repeats cannot
  read its ambiguity data, and NCBI stops with `BLAST engine error: Warning: File error: could not
  get ambiguity data. ` (exit status 3) at both thread counts.
