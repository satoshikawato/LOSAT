# shortq FINDINGS (append-only log)

## [1] Preliminary (instrumented copy of 0d533ba76 in shortq/src/LOSAT, binary shortq/target/release/LOSAT)
- Repro: in/pq2_n3.fa (>a ACG, then q2) vs q2 alone, -subject in/ss.fa -outfmt 6.
- Instrumented values for batch [3nt, q2(3000nt)] (LOSAT_SHORTQ_DBG=1):
  nq=2 avg_query_length=500 (q2 alone: 1000); contexts listed = 8 (frame +1 and -1 of 'a', 6 of q2) instead of NCBI's 12;
  min lambda = 0.2779 (frame +1 of 'a' = T, 1 aa) vs q2 alone 0.3162; cutoff_score_min=11 vs 33;
  link cutoffs small/big: 22/37 (s1) vs 33/37 (q2 alone).
- NCBI also treats the 1-aa context as valid: `tblastx -query a3.fa -outfmt 0` footer: Lambda 0.278 K 0.0810 H 0.188, eff. search space 2500
  == LOSAT's values for ctx 0 (lambda=0.277864592 k=0.081043369 h=0.188476204, eff_sp=2500).
- ABLATION: overriding avg_query_length (->1000), smallest-Lambda block (exclude 1-aa contexts) and cutoff_score_min (->33)
  ALL to the q2-alone values does NOT fix LOSAT's output  => the batch-wide cutoff values are NOT the cause.
- ROOT CAUSE FOUND (not a batch-wide scalar): context NUMBERING. LOSAT's `contexts` vector is compact (generate_frames skips frames
  with i+3 > seq_len, translation.rs:98-134), but sum_stats_linking/linking.rs:180-182 `translated_context_group = hit.ctx_idx / 3`
  mirrors NCBI's `context/(NUM_FRAMES/2)` (link_hsps.c:343-349, 510-531) which assumes 6 contexts per query (blast_setup_cxx.cpp:139,
  SetupQueryInfo_OMF always creates kNumContexts=6 per query incl. zero-length contexts).
  With a 3-nt query first (2 contexts: +1,-1) q2's frames get ctx_idx 2..7 -> groups {+1 | +2,+3,-1 | -2,-3}: wrong strand grouping
  => HSPs of q2 frames +1 / +2 / -1 are linked/not linked in the wrong sets.
- Fix test (private copy): group = q_idx*2 + (q_frame<0)  -> LOSAT == NCBI == q2 alone for pq2_n3.

## [2] FINAL: root cause, evidence, proposed change  (everything below reproducible from shortq/w/*.sh)

### Root cause
NOT a batch-wide scalar. The grouping of HSPs into "query frame groups" for sum-statistics linking uses `ctx_idx / 3`, which
equals NCBI's `context/(NUM_FRAMES/2)` only if LOSAT's context vector holds 6 contexts per query like NCBI's. It does not for
queries of 1-4 nt: `generate_frames` drops every frame with `i + 3 > seq_len`
(LOSAT translation.rs:105 / :119): 1-2 nt -> 0 contexts, 3 nt -> 2 (+1,-1), 4 nt -> 4 (+1,+2,-1,-2), >=5 nt -> 6 (also all-N >=5 nt:
the frames exist, only is_valid=false).  `prepare_lookup_query` (lookup/backbone.rs:352-412) numbers contexts by position in
that compact list, UngappedHit.ctx_idx is that position, and linking.rs:180-182 `translated_context_group = hit.ctx_idx / 3`
(used by rev_compare_hsps_tbx 196-197, rev_compare_hsps_transl 216-217, fwd_compare_hsps_transl 233-234, and the frame-group split
636/640) therefore mis-groups every later query of the batch when the number of contexts before it is not a multiple of 3.
NCBI: 6 contexts per query always (BLAST_GetNumberOfContexts = NUM_FRAMES, blast_util.c:1373-1376, blast_def.h:88;
BlastQueryInfoNew last_context = nq*6-1, blast_query_info.c:77; SetupQueryInfo_OMF fills kNumContexts=6 contexts per query
even for empty frames, blast_setup_cxx.cpp:164,193-225,268, s_QueryInfo_SetContext :64-90 marks len-0 contexts invalid at :79,:85), so
hsp->context/3 (link_hsps.c:168-169, 246-247, 343-344, strand_factor grouping 513, 523-524) == 2*query_index + (frame<0).

### Intermediate values, batch [>a ACG, q2 (3000 nt)] (in/pq2_n3.fa), -evalue 10
LOSAT (instrumented, LOSAT_SHORTQ_DBG): contexts=8: idx0 a/+1 aa_len1 lambda .277865 K .081043 H .188476 eff_sp 2500 cutoff 11;
 idx1 a/-1 aa_len1 (ideal lambda .317606) cutoff 12; idx2..7 = q2 +1,+2,+3,-1,-2,-3 (aa_len 1000,999,999,1000,999,999).
 avg_query_length 500, min-lambda block = idx0 (.277865,.081043,.188476), cutoff_score_min 11,
 link cutoffs s1 small/big 22/37, s2 25/38 (q2 alone: avg 1000, min lambda .316192, cutoff_score_min 33, cutoffs 33/37, 33/37).
NCBI (source): 12 contexts, query_length [1,0,0,1,0,0,1000,999,999,1000,999,999], offsets via s_QueryInfo_SetContext
 [0,2,2,2,4,4,4,1005,2005,3005,4006,5006]; avg = (5006+999-1)/12 = 500 (blast_parameters.c:1023-1026; smallest-lambda block :1013, s_BlastFindSmallestLambda :92-112; cutoff_score_min :316-416, used :1069) == LOSAT.
 NCBI -outfmt 0 on ACG alone prints Lambda .278 K .0810 H .188, eff. search space 2500 == LOSAT idx0, so NCBI also has the 1-aa
 context valid and uses it for the smallest-lambda block and cutoff_score_min exactly as LOSAT.  => batch-wide values agree.
The two linked HSPs (q 870-956 frame +3, and q 847-873 frame +1; LOSAT_TRACE_HSP): LOSAT ctx_idx 4 (group 1) and 2 (group 0)
 -> never linked (E 0.57 / 3.6).  q2 alone: ctx 2 and 0 (both group 0) -> linked, 0.032.  NCBI contexts 8 and 6 -> both group 2.
Ablation: overriding avg (->1000), smallest-lambda block (exclude 1-aa contexts) and cutoff_score_min (->33) to the q2-alone
 values together does NOT change LOSAT's wrong output.  Replacing the group function fixes it.

### Batch-wide values verified against NCBI on sensitive fixtures (LOSAT-fixed == NCBI):
 [TGG,q2] -evalue 1000 and [ACG,q1] / [ACGT,q1] -evalue 1000/100000: NCBI output differs from "tail alone" (so NCBI does use the
 short query's contexts for min lambda / cutoff_score_min / avg), LOSAT-fixed reproduces it exactly, and LOSAT with
 min-lambda excluding 1-aa contexts, cutoff_score_min=33, or avg 100/2500 does NOT (avg sweep: match only for 500..1100, i.e. the
 12-context average 750, not q1-alone 1500).

### Proposed change (not applied): shortq/proposed_fix.diff  (one function)
--- a/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs
+++ b/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs
 fn translated_context_group(hit: &UngappedHit) -> usize {
-    hit.ctx_idx / 3
+    hit.q_idx as usize * 2 + usize::from(hit.q_frame < 0)
 }
 (follows link_hsps.c:168-169,246-247,343-344,513,523-524 with NCBI's fixed 6-contexts-per-query numbering). q_idx is the
 batch-local index at link time (offset by `start` only after search_query_batch returns, run_impl.rs run_in_pool).
 Identical to ctx_idx/3 whenever all earlier queries have 6 contexts (verified: no output change on q1/q2/q1+q2/q2+q1+q2/5-6 nt).
 Unit test test_final_link_sorts_preserve_payloads_addresses_and_all_ties (linking.rs ~3085) sets only ctx_idx; it should set
 q_idx/q_frame so its groups are exercised (passes either way). Other ctx_idx uses (purge_endpoints.rs, hsp_culling.rs,
 run_impl.rs:467-469) only compare/equate and the compact order is monotone in (query, frame) => unaffected.
 blastx indexes contexts q*6+i (query_setup.rs:290, results.rs:1402), so it is not affected.

### Scope of the divergence (169 NCBI/LOSAT-orig/LOSAT-fixed comparisons, w/final_matrix.log; fixed == NCBI in 169/169, orig differs in 105)
 Trigger: a later query exists in the same batch and the contexts of all earlier queries total != 0 (mod 3), where
 contexts(len) = 0 (<=2 nt), 2 (3 nt), 4 (4 nt), 6 (>=5 nt). All-N and partly-N queries count like any other (NNN->2, NNNN->4;
 N>=5 or frames with N -> 6, no problem). Examples: [3],[4],[N3],[N4],[1,3],[3,1],[3,3],[4,4],[q1,3],[q1,4],[3,q1,3] break;
 [1],[2],[5],[6],[3,4],[4,3],[3,3,3],[N3,N4],[q2,3] (short query LAST) do not.  The short query's own HSPs are not affected
 (1-aa frames cannot link).  Position: breaks every query AFTER the short one in the batch, not those before it.
 BATCH_SIZE: short query alone in its batch (BATCH_SIZE<=3 for 3 nt, <=4 for 4 nt) or last in its batch -> no effect; first-but-not-alone
 or in the middle -> breaks (e.g. BATCH_SIZE=4,5,3002..10002 for [ACG,q2]; 4504,4507,7503 for [q1,ACG,q2]).
