# SEG mask investigation findings
Fri Oct  2 21:20:56 JST 2026

## Log
### 0. Setup
- LOSAT SEG: LOSAT/src/utils/seg.rs (core/blast_seg.rs re-exports it). Prior note (nisland/FINDINGS.md sec. 5): divergence at ctx0 q aa 54-113 (NCBI score 148) vs 54-84 (LOSAT 99).

### 1. NCBI masks via LD_PRELOAD shim (works: SeqBufferSeg/BlastSeqLocCombine interposable)
- shim: segmask/shim/segshim.c (+.so) intercepts SeqBufferSeg (dumps input bytes in ncbistdaa, output locs) and BlastSeqLocCombine. Run: `LD_PRELOAD=segmask/shim/segshim.so tblastx -query q_only.fa -subject s_noamb.fa -outfmt 6 2> shim.err`
- NCBI SeqBufferSeg is called 6 times (ctx0..5), each len=177 aa (frames of the 533-nt query), window 12, default params. NCBI per-context masks (0-based inclusive aa coordinates):
  ctx0: 33-52 96-105 113-142 ; ctx1: 33-50 78-92 96-106 113-141 ; ctx2: 31-55 95-104 ; ctx3: 35-64 67-93 121-144 ; ctx4: 35-64 67-83 126-144 ; ctx5: 34-63 71-80 125-143
  BlastSeqLocCombine changes nothing here.

### 2. LOSAT masks (instrumented copy inst/, env LOSAT_DUMP_SEG=1) vs NCBI, repro q_only.fa
- LOSAT (0-based inclusive aa): ctx0(+1): 33-52 **84-95** 96-105 113-142 ; ctx1: 33-50 78-92 96-106 113-141 ; ctx2: 31-55 95-104 ; ctx3..5 identical to NCBI.
- NCBI ctx0: 33-52 96-105 113-142. ONLY difference: LOSAT has an extra mask 84-95 (AAKLAARSAAVA, 12-mer with 7 A: entropy 1.947 <= locut) in frame +1. All other 5 frames identical.
- Literal Python port of blast_seg.c (py/ncbiseg.py) reproduces NCBI's masks for all 6 contexts exactly (incl. the trace below).
- Trace of NCBI s_SegSeq on ctx0: trigger i=87 (loi 86, hii 142) trim (81,148)->(113,142); i+upset-1=93 < leftend=113 => recursion on left piece [81..112] (offset 81). Inside it: trigger i=6 trim (0,31)->(15,24) => seg 96-105, again i+upset-1=12 < 15 => nested recursion on [81..95] which finds seg (3,14)+81 = 84-95. The nested call returns leftsegs=[84-95] (1 elem), pushed; then 96-105 pushed => the first-level leftsegs list = [96-105, 84-95] (2 elems).

### 3. ROOT CAUSE (proved)
NCBI s_SegSeq (blast_seg.c:2086-2101) recurses on the left part of a trimmed segment with a fresh list `SSeg *leftsegs = NULL` and then links it with `leftsegs->next = *segs; *segs = leftsegs;` (lines 2096-2100). That assignment overwrites the `next` pointer of the head of leftsegs, so when the recursive call found MORE THAN ONE segment, only the head (most recently found) survives; the rest of leftsegs is lost (leaked). LOSAT's seg.rs:811-857 passes the shared `segs` Vec into the recursive call, so all segments of the recursion are kept.
In the repro: ctx0 trigger i=87 -> trim (81,148)->(113,142) -> recursion on [81..112] finds [96-105, 84-95] (2 segs: nested recursion on [81..95] found 84-95 and the level itself found 96-105). NCBI keeps only the head 96-105, drops 84-95. LOSAT keeps both => extra mask 84-95 (AAKLAARSAAVA) in frame +1 => lookup words overlapping 84-95 are not indexed/seeded (query masked as X) so the seed that extended to 54-113 in NCBI is lost; LOSAT's HSP is 54-84 (99).
Python literal port: dropbug=True == NCBI on all 6 ctx; dropbug=False == LOSAT's masks.

### 4. Fix applied in private copy fixclean/ (diff: proposed_fix.diff, 23 added lines incl. comments)
- Repro (q_only x s_noamb, and rnd2/q1898 x s1898), tblastx default SEG, outfmt 6/0/7: cur-vs-NCBI diff lines 45/195/45 (and 45/182/45) -> fix-vs-NCBI 0/0/0.
- Direct SEG test (examples/segdump.rs calling SegMasker::mask_sequence on the exact bytes NCBI's SeqBufferSeg received, via blastp -seg ... with the shim; 3000 random protein queries with low-complexity runs, 4 SEG param sets "12 2.2 2.5","10 1.8 2.1","20 2.0 3.0","15 1.5 1.8"): orig differs from NCBI in 1/12000, fix 0/12000. The 1 case is a plain protein (blastp reproducer, p-index in corpus/prot3000.fa).

### 5. BLASTP reproducer (plain protein, found by the direct SEG scan above)
- `repro_blastp/q173.fa` (458 aa, protein p173 of corpus/prot3000.fa; NCBI SEG `-seg yes` masks 31-120 127-166 236-255 261-300 332-361 395-410, LOSAT additionally 224-235 KQPKKQPKKQPK) vs `repro_blastp/s_a.fa`: `blastp -query q173.fa -subject s_a.fa -seg yes -outfmt 6`: NCBI `217-226 31-40 evalue 0.006 bitscore 21.2`; current LOSAT `217-224, 0.003, 21.9`; fix == NCBI (also s_b, s_c; s_d outfmt 0). outfmt 0: only remaining diff after fix is the 'Method: Composition-based stats' vs 'Compositional matrix adjust' label (pre-existing, unrelated to SEG; cur has it too).

### 6. Direct SEG function test, large (NCBI SeqBufferSeg via shim vs LOSAT SegMasker orig/fix on the identical byte sequences)
- 100000 protein SEG calls (generators py/genprot.py, py/genprot2.py: runs of homopolymers/short repeats/biased composition, adjacent LC pieces, X/*/B/Z/U/O/J; params "12 2.2 2.5" x40000, "10 1.8 2.1" x20000, "20 2.0 3.0" x20000, "8 1.5 2.0" x20000): orig differs from NCBI in 93 (0.093%), fix in 0. Earlier runs (12000 + 6834 calls): orig 13, fix 0. No other SEG divergence was found in these 119k calls (no X/stop/alphabet/trim/merge difference).

### 7. End-to-end vs NCBI 2.17.0 (cur = native/release/LOSAT [before], fix = fixclean/target/release/LOSAT [after]); cases/ built by py/build_cases.py, runner cmpt.sh / runcases.sh, raw lines res_*.txt, table summary_cases.txt
Pool "bug" = 97 protein sequences (default SEG params) on which orig SegMasker != NCBI SeqBufferSeg (from the 100k direct scan), x4 variants each (subject = mutated copy/window, random strand/flanks) = 388 cases/program. Pool "ctl" = 400 random LC-rich proteins (py/genprot2.py seeds 300000+), 1 case each. blastp/tblastn use `-seg yes`; blastx/tblastx default SEG.
| prog | set | fmt | runs | cur!=NCBI | fix!=NCBI |
| blastp | bug | 6 | 388 | 344 | 0 |
| blastp | bug | 0 | 388 | 388 | 0 (*) |
| blastp | ctl | 6 | 400 | 2 | 0 |
| blastp | ctl | 0 | 400 | 188 -> 2 (*) | 187 -> 0 (*) |
| tblastn | bug | 6 / 0 | 388 | 337 / 388 | 0 / 0 |
| tblastn | ctl | 6 / 0 | 400 | 1 / 2 | 0 / 0 |
| blastx | bug | 6 / 0 | 388 | 344 / 343 | 0 / 0 |
| blastx | ctl | 6 / 0 | 400 | 2 / 2 | 0 / 0 |
| tblastx | bug | 6 / 0 | 388 | 355 / 355 | 0 / 0 |
| tblastx | ctl | 6 / 0 | 400 | 5 / 5 | 0 / 0 |
(*) blastp -outfmt 0 has two UNRELATED pre-existing differences present in cur and fix alike: LOSAT prints SEG-masked query residues uppercase where NCBI prints lowercase, and the 'Method:' label (Composition-based stats vs Compositional matrix adjust) differs in some alignments. Comparing case- and Method-insensitively: bug cur!=NCBI 388 -> fix 0; ctl cur 2 -> fix 0.
fix!=NCBI is 0 everywhere else, so the fix never moves an output away from NCBI (fix!=cur only where cur!=NCBI).
- tblastx rnd2 corpus (nisland/gen2.py seeds 1-5000, DNA-level LC inserts, N/IUPAC): default SEG, outfmt 6 and 0: 5000 runs each: cur!=NCBI 1 (seed 1898) / 1, fix!=NCBI 0 / 0.
- Larger LC-rich random corpus ("big": 3000 random proteins per program, 2/3 genprot2 + 1/3 genprot, back-translated/mutated per program as above; cases/<prog>_big): 
| blastp | 6 | cur!=NCBI 6 | fix 0 | ; blastp fmt0 case/Method-insensitive: see (*) ; tblastn fmt6: 5 -> 0 ; tblastn fmt0: 41 -> 35 (the 35 are an UNRELATED pre-existing footer difference: Lambda/K/H line 0.317 vs 0.318, present in cur too; the 6 SEG ones fixed) ; blastx fmt6/0: 6/6 -> 0/0 ; tblastx fmt6/0: 8/8 -> 0/0.
- Unit tests in fixclean (private copy with diff): `cargo test --release --lib`: 655 passed, 0 failed, 1 ignored (incl. 37 seg-related tests that compare against NCBI segmasker/core SEG results: utils::seg::tests::*ncbi*, kappa/redo_alignment seg tests, tblastn hard_seg tests).

### 8. PROPOSED CHANGE (not applied to any repository; applied only in segmask/fixclean/LOSAT). File segmask/proposed_fix.diff, against LOSAT/src/utils/seg.rs (working tree LOSAT-web-gui, 2026-09-28):
```diff
--- /mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/src/utils/seg.rs	2026-09-28 23:42:23.515604400 +0900
+++ src/utils/seg.rs	2026-10-02 21:26:41.251447182 +0900
@@ -840,7 +840,28 @@
                         if lend < leftend {
                             let rend = leftend - 1;
                             if rend < seq.len() && lend <= rend {
-                                seg_seq(masker, &seq[lend..=rend], offset + lend, segs);
+                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2086-2101
+                                // ```c
+                                // SSeg *leftsegs = (SSeg*) NULL;
+                                // status = s_SegSeq(leftseq, sparamsp, &leftsegs, offset+lend);
+                                // ...
+                                // /* prepend here, order will be restored in s_SegToSeqLoc */
+                                // if (leftsegs!=NULL)
+                                // {
+                                //    leftsegs->next = *segs;
+                                //    *segs = leftsegs;
+                                // }
+                                // ```
+                                // `leftsegs->next = *segs` overwrites the link of the head of
+                                // the recursive result, so when the recursion found more than
+                                // one segment only its head (the most recently found one)
+                                // survives; the remaining segments are dropped (leaked in C).
+                                // This is deterministic NCBI behavior and is reproduced.
+                                let mut left_segs: Vec<(usize, usize)> = Vec::new();
+                                seg_seq(masker, &seq[lend..=rend], offset + lend, &mut left_segs);
+                                if let Some(&head) = left_segs.first() {
+                                    segs.insert(0, head);
+                                }
                             }
                         }
                     }
```
Only `seg_seq`'s recursive call inside `SegMasker::mask_sequence` changes (seg.rs:843). core/blast_seg.rs re-exports this, so BLASTP query/subject (kappa/redo_alignment), BLASTX, TBLASTN (query + subject), TBLASTX all get the same NCBI-faithful behavior with one edit.
NCBI citations: c++/src/algo/blast/core/blast_seg.c:2030 s_SegSeq; 2086-2101 recursion (2087 `SSeg *leftsegs = NULL`, 2088 recursive call, 2095-2100 `if (leftsegs!=NULL) { leftsegs->next = *segs; *segs = leftsegs; }`); caller BlastSetUp_Filter blast_filter.c:1122-1160 (SeqBufferSeg, overlaps=TRUE) via s_GetFilteringLocationsForOneContext blast_filter.c:1217-1256 (per context/frame, mask in protein coords, BlastSeqLocCombine(…,0)).
LOSAT citations: LOSAT/src/utils/seg.rs:811 `fn seg_seq`, 836-846 recursion (shared `segs` Vec -> keeps all left segments), 856 `segs.insert(0, …)`.

### 9. Other findings (out of scope, pre-existing, present in cur and fix)
- blastp -outfmt 0 with `-seg yes`: LOSAT does not print SEG-masked query residues in lowercase in the alignment (NCBI does). tblastx/tblastn/blastx print them (lowercase masks matched in all fmt-0 tests).
- blastp -outfmt 0 'Method:' label (Composition-based stats vs Compositional matrix adjust) differs in some alignments (e.g. repro_blastp/s_a).
- tblastn -outfmt 0 footer (Lambda/K/H line e.g. 0.317 vs 0.318) differs in ~35/3000 LC-rich random cases.
- Gates run with the fix binary: tests/tblastx_regression_fixtures.py check: 51 cases, 0 differ; docs/evidence/losat_web_e2a/check_losat.py --programs blastp,tblastn,tblastx: 28 fixtures, differing=[] (2 approved code4 exceptions); cargo test --release --lib 655 passed.
