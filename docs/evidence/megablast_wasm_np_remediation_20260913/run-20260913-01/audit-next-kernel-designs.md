# Independent preimplementation kernel design checks

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), 2026-09-13. Read-only; no implementation, builds, measurements, or speed claim.

## BLASTN fixed bands

B1 `gapped.rs` SHA256 `85cfb33aebd5dbbf7e7460284bf51e55625e656793747ebf8a7d03df950513c4`. Reviewed the passive `prepare_n_band_candidate.py` before execution.

The BLASTN row allocator resizes to at least row_capacity. This covers the entire current band for positive gap extension; with zero extension, n+3-first_b_index also covers b_size<=n+1. Thus the old trace-row bounds check is always true for valid state. Capturing band_start independently from the X-drop cursor retains every ascending visited cell; band_index equals the original trace index. Fence stops before writing exactly as before, with the same row_end, final row size, reverse access and terminal sentinel. Existing unwritten row bytes from scratch reuse are preserved.

Only the first raw score pointer in each trace function is removed; the second pointer after DP reallocation remains for band growth and sentinel assignment. The mutable slices end before that allocation. The packed loop keeps dp_cells and all arithmetic/ties unchanged. The script originally replaced NCBI comment text along with Rust accesses; this was corrected to skip comment lines. Its packed-loop reference is now 3146–3197.

Boundary evidence should cover both actual trace functions/directions, first/interior fence, contracting/expanding band, zero gap extension, n+1 sentinel, nonmultiple-of-four packed lengths and offsets. Full BLASTN row comparisons require **the same scratch history in B1 and candidate**: fresh and reused rows need not have equal unwritten bytes. Serial Wasm formerly used a separate Vec access path, so it must be tested alongside native/threaded.

## BLASTP initialized trace prefix

B1 `gapalign.rs` SHA256 `d70563703948c8d6ae19b3cced19350d626e50f06f5fbd8e8c3e4bbcacf0ec96`.

After row clear, reserve(band_len) guarantees the actual band capacity independently of the existing relative-reserve shortfall. Pair exact equal-length DP and spare-capacity slices. Every non-fence iteration must write script, including X-drop failure. orig_b_index and band_start remain equal; moving the X-drop cursor makes no holes. Normal completion initializes band_len bytes, fence initializes only the preceding prefix, and an empty band initializes zero. Existing row_end_b_index - band_start is therefore the exact initialized length.

Expose that length once after normal loop completion and before the outer fence/empty-band break. Do not expose it from a panic Drop guard. The two vectors do not alias, and the spare-capacity borrow must finish before set_len/growth. Retain the existing final used_cells +1 **resize**, not an unsafe extension into uninitialized bytes. NCBI allocator accounting is at blast_gapalign.c:669.

Existing BLASTP fresh/reused tests cover matrix variants, reverse, zero gap extension, changing band and fence, comparing used row bytes and offsets. They are useful but cannot alone prove B1 equivalence if both executions share the same bug; raw oracle comparisons remain required. The current cost diagnostics' 2,094,514,564 / 3,240,532,528 visited trace cells and saved output hashes were independently verified. These counts do not identify the helper's wall-time share.
