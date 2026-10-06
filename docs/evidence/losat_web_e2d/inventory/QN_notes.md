# QN notes: the query range for nucleotide queries (blastn, tblastx)

All NCBI citations are relative to `c++/`, commit 598d8ae6; oracle = NCBI BLAST+ 2.17.0 `blastn`/`tblastx` from `/home/kawato/micromamba/bin`. Scratch files (inputs, outputs `ev_*.{cmd,out,err,rc}`, helper scripts `mk.py`, `r.sh`, `eq.py`, `fr.py`) are in `scratch_QN/` and its sub-directories `edge/ lc/ dust/ bn/ tx/ split/`.

## 0. Method: the "cut and shift" model, and where NCBI deviates from it
The decisive question for the porting session: is a ranged search = a search of the record cut to the interval, with qstart/qend shifted by from-1? `scratch_QN/eq.py` runs NCBI with `-query_loc F-T` and, separately, NCBI on a copy of every record cut to `[F-1, min(T,len))` (no `-query_loc`), shifts the cut run's qstart/qend by F-1 and compares all columns of `-outfmt "6 std qlen ..."` (qlen normalised, because qlen = range as typed, RF). Result:
- blastn, megablast and `-task blastn`, query plus-strand hits and query minus-strand hits (both `q1.fa` and its reverse complement `q1rc.fa`), ranges 1-8000, 1001-6000, 1500-6003, 3000-9000 (end past record end 8000), 2400-3300, `-strand plus|minus`, three records with one range: byte-identical to the cut run on every row (coordinates, % identity, e-value, bit score, order). Evidence: `bn/bneq_*_{ranged,cut}.out`.
- tblastx, ranges with start offsets 0, 1, 2 mod 3 and lengths 5000, 4999, 4998, 5001, 5002: identical rows (same coordinates, e-values, bit scores, same order) EXCEPT the `qframe` column, which is record-relative (section 4). Evidence: `tx/txeq_*`, `tx/fr.py` output.
So: the engine (lookup, extension, statistics, cutoffs, search space) sees exactly the interval. The deviations from the cut model that matter are listed in section 7 (D1-D6).

## 1. Call path followed (NCBI, bl2seq `-query Q -subject S`, blastn / tblastx)
1. `blastn_app.cpp:116-170` `CBlastnApp::Run` -> `x_RunMTBySplitDB` (bl2seq has no `-db`, so the by-query MT branch is never taken; `-mt_mode` is rejected by LOSAT) / `tblastx_app.cpp:100-215`. Order: options (`SetOptions`; `-query_loc` is parsed by `CQueryOptionsArgs::ExtractAlgorithmOptions`, blast_args.cpp:1996 `ParseSequenceRange`, 1-based -> 0-based, `from>=1,to>=1,from<to` else "Invalid specification of query location (...)": LA's range), subject read, `CBlastInputSourceConfig iconfig(dlconfig, strand, lcase, believe_defline, query_opts->GetRange())`, `IsIStreamEmpty` ("Query is Empty!"), `CBlastFastaInputSource`, `CBlastInput` (batch size), formatter (`SetQueryRange` = the range AS TYPED, RF), prolog, then per batch: `input.GetNextSeqBatch` -> `CObjMgr_QueryFactory` -> `CLocalBlast::Run` -> `PrintOneResultSet`.
2. `CBlastInput::GetNextSeqBatch` (blast_input.cpp:135-170) -> `CBlastFastaInputSource::GetNextSequence` (blast_fasta_input.cpp:484-505) -> `x_FastaToSeqLoc` (377-466): `SaveMask` (lowercase, whole record), `ReadOneSeq` (whole record, invalid letters removed with a warning), type checks, strand of the interval, range checks (433-453), `SetFrom(from0)`, `SetTo((to0>0 && to0<seqlen) ? to0 : seqlen-1)`, id. The batch loop swallows every `std::exception` of `GetNextSequence` (line 153) and adds `sequence::GetLength(id)` = FULL record length to `size_read` (160-163).
3. `CObjMgr_QueryFactory(CBlastQueryVector&)` ("Empty CBlastQueryVector" if the batch is empty), `CObjMgr_LocalQueryData`, `CBlastQuerySourceOM`:
   - `GetLength(i)` = `sequence::GetLength(seqloc)` = interval length (blast_objmgr_tools.cpp:324).
   - `SetupQueryInfo_OMF` (blast_setup_cxx.cpp:153-281): contexts of interval length (blastn 2, tblastx 6 with protein lengths of the interval).
   - `SafeSetupQueries` -> `SetupQueries_OMF` (485-658): per query `s_GetRestrictedBlastSeqLocs` (masks restricted to the interval and shifted, `BlastSeqLoc_RestrictToInterval` blast_setup.c:1030), `GetBlastSequence` -> `CBlastSeqVectorOM(seqloc)` (letters of the interval, strand of the interval), translation of the interval for tblastx, `s_AddMask`; per-query exceptions ("Sequence contains no data") become warnings; `BlastSetup_Validate` throws when no context is valid.
   - masks: `x_CalculateMasks` (blast_objmgr_tools.cpp:176-228) -> `Blast_FindDustFilterLoc(CBlastQueryVector&)` (dust_filter.cpp:166-186): DUST on `CSeqVector(*GetQuerySeqLoc(i))` (interval letters), result mapped to record coordinates by a `CSeq_loc_Mapper(Whole -> interval)`, merged with the lowercase masks (record coordinates) and stored with `SetMaskedRegions`. SEG (tblastx) is core: `BlastSetUp_GetFilteringLocations` on the translated frame buffers of the interval.
   - `CLocalBlast::Run`: query splitting `CQuerySplitter` (query_data.cpp `GetSumOfSequenceLengths` = sum of INTERVAL lengths; `split_query_cxx.cpp`), lookup table / diag table / statistics from `BlastQueryInfo` (all interval-based), traceback, `BLAST_Results2CSeqAlign` + `RemapToQueryLoc` (blast_seqalign.cpp:1520: `OffsetRow(query_row, interval.from)`), `CSetupFactory::CreateScoreBlock` -> `Blast_GetSeqLocInfoVector` (setup_factory.cpp:190) for the displayed masks.
4. Formatting (RF): qstart/qend are record positions, `qlen` = typed range length, outfmt 0 "Length=" = record length, frames record-relative (section 4 of this file).

## 2. Coordinate conventions at each step
| step | coordinates |
|---|---|
| `-query_loc a-b` typed | 1-based closed, applies to every query record |
| `ParseSequenceRange` | from0=a-1, to0=b-1, 0-based closed |
| interval (`CSeq_interval`) | [from0, min(to0, len-1)] 0-based closed, record coordinates, nucleotide letters; strand both (or plus/minus with -strand) |
| `BlastQueryInfo`, query buffers, HSPs inside the engine | RELATIVE to the interval start (context offsets, plus strand = interval, minus = revcomp of interval; tblastx frames start at the interval's first letter / last letter) |
| lowercase masks (CFastaReader) / DUST-after-mapper masks in the query vector | record coordinates, 0-based closed, frames +1/-1 |
| mask restriction at setup (`BlastSeqLoc_RestrictToInterval`) | record -> interval: left=max(0,left-from0), right=min(right,to0)-from0 |
| HSPs after `RemapToQueryLoc` | record coordinates (interval offset added) |
| query split chunks | concatenated INTERVAL-length coordinates (query_ranges built from interval lengths); chunk Seq-locs are record coordinates (offset from0 added); chunk masks: see row 23 |

## 3. Oracle observations (command, evidence)
All runs from `scratch_QN`; unset BLASTDB/BATCH_SIZE/CHUNK_SIZE except where `X_CHUNK_SIZE` is stated (the helper `r.sh` unsets and re-exports them). `eq.py PROG Q S FROM TO OUTFMT [args]` = ranged run vs cut run (see section 0).
- Inputs: `q1.fa` = LC738874 233001-241000 (8000), `s_plus.fa` = LC738875 264001-271000 (7000, hits on the minus strand of q1; `q1rc.fa` gives query-minus hits), `qm.fa`/`qidx.fa` multi-record, `lq.fa/ls.fa` (lowercase 601-800), `dust/q60.fa` (AT*30 at 301-360), `tx/sq.fa` (CAG*40 at 301-420), `split/qlc.fa`, `qplain.fa`, `q4lc.fa`, `qdust.fa` (EDL933 first 3.5 Mb with lowercase / low-complexity inserts), `tx/big4.fa` (291 kb).
- Cut model: 20 blastn runs (megablast/blastn, both strands, 5 ranges), 8 `-strand`, 4 `-subject_besthit/-max_hsps/dc-megablast/blastn-short`, 2 multi-record, 12 lowercase, 12 DUST, 8 `-strand`+DUST, N/R ambiguity inside/outside; tblastx 8 ranges (start offsets 0,1,2 mod 3, lengths mod 3), 9 SEG, 8 lowercase, gencode 4/11/2, `-num_threads 4`, 2 split ranges: all identical to the cut run except qframe (section 4), the split-mask defect (row 23) and tblastx outfmt 0 lowercase (row 27, UNSURE). Files: `bn/bneq_*`, `lc/lceq_*`, `dust/dueq_*`, `dust/duxs_*`, `tx/txeq_*`, `tx/sgeq_*`, `tx/lcx_*`, `tx/gceq*`, `amb/amb*`.
- DUST independence of outside letters: `dust/dust_range_hardmask.txt` (hard-masked query, AT*30 at 301-360): -query_loc 331-660/321-660/311-660/340-660 -> HSP from 361 (masked), 351-660 -> HSP from 351 (10 run letters inside: not masked), 1-340 -> HSP 1-299 (masked), 1-310 -> unmasked hits inside the run. The run of 30 letters inside a 700-letter query IS masked when it sits at the interval start (DUST window edge effect) but NOT inside the full record (30 letters, not masked in `dust/` first sweep: threshold between 50 and 60 letters in the middle). So DUST depends only on the interval letters.
- Lowercase straddling the start: `lc/lcv_*` (hard masking): -query_loc 701-1400 with -lcase_masking -> HSP 801-1400 only (mask clipped to 701-800 and shifted); without -lcase_masking the whole range aligns; 601-800 fully lowercase -> no HSP and the invalid-query warning of an unranged fully masked query.
- Batches: `batch/ev_b_AB_*` / `ev_b_BA_*` and `tx/ev_bt_*` (full-length batching, rows 5-6); `edge/*`, `batch/ev_skip2_*` (rows 4, 7, 8).
- Split: `split/ev_sp_*`, `ev_m4_*`, `speq_*`, `spdu_*` (rows 21-24); tx split `tx/txsplit_*` (row 25).
- LOSAT-before: `-query_loc` -> exit 2 "error: the NCBI BLAST+ option -query_loc is not supported by LOSAT's BLASTN" (`bn/ev_L_blastn_qloc.*`; same for TBLASTX); an empty record -> "Error: query record 1 (empty) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's BLASTN" (`edge/emp_cmp.txt`).

## 4. Where NCBI uses the full record, the interval, or the typed range
| quantity | source |
|---|---|
| query batch composition (`size_read`) | FULL record length (`sequence::GetLength(id)`) -> rows 5, 6 |
| `GetSumOfSequenceLengths` (split decision, chunk sizes) | INTERVAL length (empty/invalid queries 0) -> row 21 |
| context lengths, buffers, lookup table, diag table, X-drop, cutoffs, length adjustment, effective search space, Karlin-Altschul blocks, avg query length (tblastx) | INTERVAL (verified through identical e-values and "Effective search space used") |
| lowercase masks | FULL record at read, restricted to the interval at setup |
| DUST | INTERVAL letters |
| SEG | INTERVAL translation |
| invalid-letter warning of the FASTA reader | FULL record (even when the invalid letter is outside the range) |
| `Query_<n>` of warnings | position of the record in the file (skipped records counted) |
| outfmt 0 "Length=" | FULL record (RF) |
| tabular `qlen` | `m_QueryRange.GetLength()` = typed range length, even past the record end (RF) |
| qstart/qend | record positions (interval offset added) |
| tblastx qframe / "Frame =" | RECORD-relative: plus f -> ((from0+f-1) mod 3)+1; minus -k -> -(((e+k-1) mod 3)+1) with e = record length - min(typed end, record length) (verified, 14 runs, 0 mismatches) |
| displayed soft-masks (blastn) | interval masks + from0 (pure shift, verified); tblastx: see row 27 (UNSURE) |

## 5. Errors and warnings with exit status (all observed)
| situation (blastn and tblastx alike) | text | rc |
|---|---|---|
| every record of the first batch has from0 > length | stdout: outfmt 0 prolog only; stderr `BLAST engine error: Empty CBlastQueryVector` | 3 |
| every record of a later batch skipped | earlier batches fully printed, no epilog; same stderr | 3 |
| record with from0 > length among others | silent, record vanishes | 0 |
| empty interval (from0 == length, or empty record with start 1) in a batch with a valid query | stderr `Warning: [blastn] Query_<n> <title>: Sequence contains no data ` (trailing space), empty result block | 0 |
| every query of a batch has an empty interval | `BLAST engine error: Warning: Sequence contains no data ` repeated once per query, no separator, then newline | 3 (timing: before any result of that batch; earlier batches already printed) |
| 1-letter interval alone in a tblastx batch | `Warning: [tblastx] Query_1 tiny: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options ` | 0 |
| end past the record end | nothing (silently clipped) | 0 |
| fully hard-masked interval | invalid-query Karlin warning as for an unranged record | 0 |
The exceptions of `ParseSequenceRange` / `CInputException eInvalidRange "Invalid sequence range"` (to<from) are unreachable from the command line (the parser rejects them first) -> LA.

## 6. Surprises
1. **Records with start > length are dropped silently** (the batch reader swallows the `CInputException`), with all the consequences of rows 4-7 (the first batch empty -> "Empty CBlastQueryVector" exit 3 after the outfmt 0 prolog; later batch empty -> error after earlier output).
2. **Batches use the full record length**: a 300 kb record with range 1-1000 is its own batch; an invalid (empty-interval) query alone in a batch changes a warning into an exit-3 error.
3. **NCBI defect in the split path** (row 23): masks restricted with an interval-relative window although they are in record coordinates; reproduced in the oracle (lowercase masks lost in a split, ranged query with start > 1). Policy: reproduce deterministic defects.
4. **tblastx qframe is record-relative**, coordinates +from0; outfmt 0 "Frame =" likewise.
5. **tblastx outfmt 0 lowercase (SEG masks) of ranged queries is not the shift of the cut run** (row 27, UNSURE): tblastx runs `-seg` by default, so this touches every ranged tblastx outfmt 0 run (outfmt 6/7 are fine).
6. DUST mask extents depend only on interval letters, so a low-complexity run cut by the range start/end is judged on the part inside (identical to a record edge).
7. Query-id numbering in warnings counts skipped records.

## 7. Deviations from the cut model (D1-D6) = what the session must reproduce beyond "cut and shift"
- D1 batches by FULL length and skipped records (rows 4-7, 28).
- D2 empty interval handling (row 8: reject unless empty records are ported).
- D3 split-path mask defect (row 23).
- D4 report side: qframe/frames record-relative (row 29), displayed masks + from0 (row 26), tblastx SEG display (row 27 UNSURE), qlen/Length (RF).
- D5 `Query_<n>` numbering (row 28).
- D6 `lengths` vectors: three different ones are needed in LOSAT's blastn `run_in_pool`/`search_query_batch`: full lengths (batch composition), interval lengths (split, contexts, subject_besthit query_len), record lengths for "Length=" (RF).

## 8. Where LOSAT must apply the range (summary for the port)
BLASTN: in `run.rs` `search()`/`run_in_pool`: after `check_residues`, `with_u_as_t` and the empty-record checks (full records), compute per record `from0`/interval; skip records with from0 > len (keep original input index); build `lengths_full` (batching), cut `query_records` to the interval (keep case) for `search_query_batch`/`search_query_chunks`/DUST/lowercase; keep `from0` per query for the report (+from0 on qstart/qend, masks) and for the split-mask defect (record-coordinate masks = cut masks + from0, restrict with window [pf, pt-1) then [pf+d, pt-1+d] as in row 23). TBLASTX: same in `run_impl.rs` `search()`/`run_in_pool` (cut before `generate_frames`; `lengths` full for `next_query_batch_end`; frames relabelled record-relative in the writers).
After that cut the existing code behaves like NCBI for contexts, masks (non-split), DUST, SEG, statistics, e-values, lookup tables (everything in rows 11-22 is faithful). The `blastn/input.rs` rejection of empty records stays; the empty interval needs the same decision.

## 9. Decisions for the session
1. Empty interval (row 8): keep rejected (cheapest, consistent with empty records) or port both together (M).
2. Whether to reproduce the split-mask defect (row 23, rare, port:S, recommended per the deterministic-defect policy).
3. Row 27 (tblastx outfmt 0 SEG lowercase with a range) needs a targeted investigation with RF before byte parity can be claimed; fallback: reject `-query_loc` for TBLASTX outfmt 0 while seg is on (not decided here).
4. `-query_loc` with `-outfmt 0` for BLASTN depends on RF's rows for "Length=", Query coordinates and displayed soft masks (+from0).
