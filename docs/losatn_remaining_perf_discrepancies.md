# LOSATN Remaining Performance Discrepancies (LOSAT vs NCBI BLASTN)

## Scope
- Focus: blastn/megablast engine performance differences that remain after recent optimizations.
- Source of truth: NCBI BLAST C/C++ code only.
- Goal: identify discrepancies that likely keep LOSATN slower than BLASTN.

## Findings (Remaining)

### 1) Subject sequences are loaded/encoded on the fly vs SeqDB pre-encoded sequences
- LOSAT preloads subject FASTA records for `use_parallel` or `limit_lookup`, then encodes each subject into BLASTNA + packed ncbi2na in `process_subject`.
- NCBI streams subjects from `SeqSrc` and stores pre-encoded sequences in `BLAST_SequenceBlk` via SeqDB.
- Impact: extra I/O pass and per-subject encoding overhead, especially for large DBs.
- Recheck: encoding uses `IUPACNA_TO_BLASTNA` + `& 3` packing with remainder in the last byte; behavior matches NCBI, so this is a pure performance gap.
- Status: open — closing the gap requires SeqDB/BLAST DB ingestion or a persistent pre-encoded cache; the streaming path still does a full metadata scan before processing.
- LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:3658`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:4342`, `LOSAT/src/algorithm/blastn/coordination.rs:522`, `LOSAT/src/core/blast_encoding.rs:429`.
- NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:1154`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85`.

### 2) Lowercase masking is a separate pass vs inline parsing
- LOSAT scans sequences to find lowercase masks after reading (`collect_lowercase_masks`), for queries and subjects.
- NCBI's `CFastaReader` tracks masks inline during parsing (`x_OpenMask` / `x_CloseMask`).
- Impact: extra full-sequence scan and mask Vec allocation.
- LOSAT refs: `LOSAT/src/algorithm/blastn/coordination.rs:561`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:4307`.
- NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/readers/fasta.cpp:856`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/readers/fasta.cpp:1079`.

## Recently addressed (no longer discrepancies)
- Parallel subject chunk processing now streams bounded batches in offset order (no full `SubjectChunk` materialization or per-chunk sort before merge); LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:7075`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:7135`; NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:478-536`.
- Query lookup build now uses pre-encoded BLASTNA buffers (no ASCII->2-bit conversion in lookup builders); LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:3923`, `LOSAT/src/algorithm/blastn/lookup.rs:23`, `LOSAT/src/algorithm/blastn/lookup.rs:953`; NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:498-604`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:90-120`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:878-893`.
- limit_lookup PV fill now iterates unmasked ranges (no per-k-mer mask checks); LOSAT refs: `LOSAT/src/algorithm/blastn/lookup.rs:249`, `LOSAT/src/algorithm/blastn/lookup.rs:505`, `LOSAT/src/algorithm/blastn/coordination.rs:722`; NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:856-889`.
- Masked chunk `seq_ranges` now build into per-thread scratch for both sequential and parallel paths (no per-chunk Vec allocations); LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:180`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:4568`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:2900`; NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:184-198`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:298-307`.
- HashMap lookup replaced with array-backed `NaLookupTable` (thick_backbone + overflow + PV), and offset buffer sizing now includes `longest_chain`; LOSAT refs: `LOSAT/src/algorithm/blastn/lookup.rs:119`, `LOSAT/src/algorithm/blastn/lookup.rs:1160`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:3982`; NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:109`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:442`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:41`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255`.
- Offset-pair buffers now use fixed-length arrays with an index (no `push`/`drain` in scan loop); LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:3030`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:5850`; NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:991-1041`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1482-1534`.
- Scan callback preselected via `select_mb_scan_kind` and passed into the scan loop (no per-range `mb_scan_kind` dispatch); LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:2103`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:5815`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs:6023`; NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2602-2677`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1651-1667`, `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2994-3006`.
- Buffered output for outfmt0/7 (`LOSAT/src/report/pairwise.rs`, `LOSAT/src/report/outfmt6.rs`).
- Reuse packed subject encoding when `limit_lookup` preloads subjects (`LOSAT/src/algorithm/blastn/lookup.rs`, `LOSAT/src/algorithm/blastn/blast_engine/run.rs`).

## Notes
- Items above focus on performance differences; parity behavior should remain unchanged if aligned to NCBI.
- Recheck after the NaLookup swap: no additional missing steps found in the non-direct lookup construction beyond the remaining items.
- Recheck after fixed offset-pair buffers: no remaining `Vec<OffsetPair>` push/drain usage in the blastn scan path.
- Recheck after scan callback preselection: scan kind selection now happens once per subject scan, mirroring `s_MBChooseScanSubject`/`BlastChooseNucleotideScanSubjectAny`.
- Recheck after parallel chunk streaming: merge order still follows `s_GetNextSubjectChunk` offsets; no additional discrepancies found.
- Recheck: query lookup now uses pre-encoded BLASTNA; items 1-2 remain open; no other discrepancies found.
- If you want, I can attach a timing profile plan to validate which items dominate your workload.
