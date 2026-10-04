# S08+a：E2e の第 1 回の独立監査の残件

- ブランチ `feature/losat-web-gui-s08pa`（worktree `LOSAT-web-gui-s08pa`）、起点 `d96412265`。
- NCBI のソースは `598d8ae6`（写し `~/.cache/losat-web-gui-target/s08p/ncbi/c++`）、比較の NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin`）。
- 作業の写しと trace は `~/.cache/losat-web-gui-target/s08pa/`（`tn2/` など。`/mnt/c` の外）。

## TN-2：hard mask と少ない `-max_target_seqs` の和の統計の e-value

### 最初に値が食い違う箇所

再現：`audit/round1/tblastn_inputs/gen_q_lc.faa` × `gen_s.fna`、`-comp_based_stats 0 -lcase_masking -seg no -max_target_seqs 1 -outfmt 6`。`gl0 gs20` の最初の行が NCBI `9.33e-177`、LOSAT `9.00e-177`。

1. NCBI 自身が `-max_target_seqs` で変わる。`gl0 gs20` は mts 20 以上で `9.00e-177`、mts 12 以下で `9.33e-177`。LOSAT は常に `9.00e-177`。
2. TLOSAN Stage D の shim（`ncbi_d_call_trace.c`）で、`gl0 × gs20` の traceback の後の `BLAST_LinkHsps` の subject の長さが mts 1 で 489（その中の `Blast_HSPListGetEvalues` は 163）、mts 500 で 488（162）。context 0 の `eff_searchsp` 5529312・length adjustment 57、traceback の後の HSP の一覧（9 個）は同じ。
3. この長さは `Blast_TracebackFromHSPList` の `stat_length`（blast_traceback.c:294,425-433）：最後に翻訳した HSP の `translated_length`。`gs20`（1467 塩基）の全翻訳の長さは frame ±1 で 489、±2・±3 で 488。`gl0` は最初の traceback で fence に触れ（blast_traceback.c:516-520）、全翻訳でやり直す（blast_traceback.c:1660-1700、`kFullTranslation`）。
4. 追加の shim（`s08pa/tn2/tb_trace.c`：`Blast_TracebackFromHSPList` の入力の HSP の一覧と fence、`Blast_HSPGetTargetTranslation` の frame と長さ）で、traceback に入る HSP の**順**が違うことを確かめた。mts 500 は得点の順（1233, 74, 73, 70, 70, 68, 65, 59, 46。最後は frame 2 → 488）、mts 1 は e-value の順（1233, 74, 73, 70, 70, 59, 46, 68, 65。59 と 46 は連結した組で e-value 0.00244 が 68 の 0.00251 より小さい。最後は frame 1 → 489）。
5. 順を変えるのは予備の hit list：`Blast_HitListUpdate`（blast_hits.c:3266-3284）は hit list が一杯になると（`hsplist_count >= hsplist_max`、予備の大きさは `MIN(MAX(2*hitlist_size,10), hitlist_size+50)`、blast_hits.c:43-70）heap にする時に保存済みの全 HSP list を `Blast_HSPListSortByEvalue` で並べ、以後に来る list も e-value の順にする。traceback（blast_traceback.c:358-365）は `_DEBUG` の時だけ得点の順に並べ直し、release は並べ直さない。collector（hspfilter_collector.c:84-160）と stream の close（blast_hspstream.c:133-208）も並べ直さない。
6. LOSAT は予備の hit list（`stage_d_results.rs` の `KappaResultHitList::update`）で NCBI と同じく e-value の順にするが、traceback（`search_gapped.rs` の `full_translation_traceback_with_matrix_and_events_with_mask_mode_owned`）が入力を `sort_gapped_score_if_needed` で得点の順に並べ直していた。そのため最後に翻訳する HSP と `stat_length` が NCBI と違った。

### 直し方

traceback の入口の並べ直しを除いた（NCBI の blast_traceback.c:358-365 と blast_hits.c:3269-3284 を直上に引用）。Stage C の観察用の wrapper（`full_translation_traceback_with_matrix_and_events_with_mask_mode`、試験だけが使う）は、1 つの subject の予備の HSP を渡すので、今までどおり得点の順にしてから呼ぶ（blast_engine.c:539-552）。予備の hit list が heap にならない場合、LOSAT の linking の出力（`stage_d_linking.rs` の `score_compare` の並べ替え、link_hsps.c:1802-1803）と連結しない場合の出力は既に得点の順なので、変わらない。

### 確かめ

`gen_q_lc.faa`・`gen_q.faa` × `gen_s.fna`、`-comp_based_stats 0 -lcase_masking -seg no`、mts 1・2・5・10・11・12・20・100・500、BLOSUM62 と BLOSUM45（`-matrix BLOSUM45 -word_size 2`）：全て NCBI と一致（変更前は mts 1〜12 で 5 行が違った）。
