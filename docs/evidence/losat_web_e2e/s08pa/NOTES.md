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

## TN-4：巨大な `-xdrop_gap`・`-xdrop_gap_final` の時間と記憶

### 原因

- NCBI の 3 つの gapped DP（`ALIGN_EX` blast_gapalign.c:457-470、`Blast_SemiGappedAlign` 797-808、`s_RestrictedGappedAlign` 1040-1050）は、`num_extra_cells = x_dropoff / gap_extend + 3` が `dp_mem_alloc` を超えると `malloc`（帯を広げる時は `realloc`、923-931 など）で `dp_mem_alloc` 個の cell を取る。DP は cell に書いてから読み、`N`（subject の残りの長さ）より先の cell を使わない（最初の行 811-821、帯の延長と番兵 945-958 は `b_size <= N`）。そのため触れるのは最初の `MIN(dp_mem_alloc, N + 1)` 個で、残りの page は確保されない。traceback の行の状態（`s_GapGetState` 70-115）も、chunk を `malloc` して使った分だけ触れる。
- LOSAT（`blastp/gapalign.rs` の `gap_dp_reserve_initial`・`gap_dp_reserve_band`）は `Vec::resize` で `dp_mem_alloc` 個の全ての cell を書いていた。`-xdrop_gap 5e8` の raw の X-drop は約 1.3e9 で、1 回の resize が約 10 GB を書く。しかも TBLASTN の traceback は HSP ごとに新しい `GapAlignScratch` を作る（`search_gapped.rs`）ので、HSP ごとに書き直した。traceback の行の `Vec` も `num_extra_cells` の容量（最大 1.7 GB の仮想記憶）を予約していた。

### 直し方

`dp_mem_alloc` は NCBI の値のまま（再確保の時点と、LOSAT の範囲の検査 `b_size < dp_mem_alloc` に使う）にし、`dp_mem` の実際の長さを `MIN(dp_mem_alloc, N + 1)`（`gap_dp_ensure_cells`）にした。後の伸長の `N` が長ければ伸ばす。traceback の行の容量の予約も、行に入りうる `N + 2 - first_b_index` 個までにした（容量は性能のための予約で、出力に関わらない）。DP が読む cell と値は変わらない（cell は書いてから読む。`s_RestrictedGappedAlign` の最初の行の上限 `dp_mem.len() - 1` も、`len2` か以前と同じ `dp_mem_alloc - 1`）。BLASTP・TBLASTN・BLASTX が共有する経路で、出力は変わらない。

## TN-5：BLOSUM45 の組の大きな `-evalue` で同じ得点の HSP の frame

### 最初に値が食い違う箇所

再現：`e2e_protein_query.faa` × `e2e_tblastn_subject.fna`、`-matrix BLOSUM45 -word_size 2 -comp_based_stats 0 -evalue 1e4 -outfmt 6`。29900 行中 4 行で subject の座標だけが違う（NCBI `4712 4695`、LOSAT `4710 4693`）。

1. 違う行は、subject の frame f と f+2（1 と 3、−1 と −3）で protein の座標・得点・query の範囲が同じ 2 つの HSP の、どちらが残るか。`4712 4695` は frame −1 の `s.offset` 96、`4710 4693` は frame −3 の同じ 96（`LvMJNV_160001_165000` は 5000 塩基）。
2. traceback に入る HSP の一覧（`s08pa/tn2/tb_trace.c` の `Blast_TracebackFromHSPList` の入力、BDT62620.1 × subject 2 の 662 個）は同じ集合で、5 組だけ順が逆（NCBI は frame −1 が先、LOSAT は −3 が先）。後の方は、containment の検査と common endpoint の purge で落ちる。
3. 予備の段の linking の入力（Stage D の shim の `BLAST_LinkHsps` の `link_before`、5928 個）で既に 12 か所（6 組）の順が逆。NCBI の `s_BlastUnevenGapLinkHSPs`（link_hsps.c:1613-1757）は別の配列で連結し、`hsp_array` の順を変えない。最後の `Blast_HSPListSortByScore`（link_hsps.c:1802-1803）は安定。
4. 順を作るのは frame ごとの HSP の追加：`Blast_HSPListAppend` → `s_BlastHSPListsCombineByScore`（blast_hits.c:2749-2768）は新しい frame の HSP を後ろに足して `Blast_HSPListSortByScore`（blast_hits.c:1374-1383、`qsort` と `ScoreCompareHSPs`）。`ScoreCompareHSPs`（blast_hits.c:1330-1356）は frame を比べないので、f と f+2 の同じ座標の HSP は同順位になる。固定した NCBI BLAST+ 2.17.0 は glibc 2.39 の上で動き、その `qsort` は安定な merge sort（同順位は追加の順、つまり frame 1, 2, 3, −1, −2, −3 の順に残る）。
5. LOSAT の `search_gapped.rs` の `sort_preliminary_by_score`（frame の追加の後の並べ替え）は `sort_unstable_by` で、同順位の順を保たなかった。traceback の `purge_traceback_common_endpoints`（blast_hits.c:2478-2486,2504 の `qsort`。比較は frame を見ず、削除は隣の HSP が同じ frame の時だけ）、chunk の併合、初めの hit の並べ替え（blast_extend.c:306-310）も `sort_unstable_by` だった。

### 直し方

TBLASTN だけが使う、NCBI の `qsort` に当たる並べ替え（`search_gapped.rs` の 7 か所と `search_init.rs` の 1 か所）を安定な `sort_by` にした（NCBI の `qsort` の行を直上に引用）。BLASTP（`blastp/hsp.rs`、`blast_engine.rs`）と BLASTX と共有する composition の窓の並べ替え（`redo_alignment.rs`）は変えていない（下の「残り」）。

### fixture と確かめ

- fixture `e2e.tblastn.b45_frame_ties`（outfmt 6）と `e2e.tblastn.b45_frame_ties_fmt0`（outfmt 0）：BDT62620.1 × `LvMJNV_160001_165000`（`e2e_frame_tie_query.faa`・`e2e_frame_tie_subject.fna`）、`-matrix BLOSUM45 -word_size 2 -comp_based_stats 0 -evalue 1000`。変更前は 3 行（outfmt 0 は 6 行）違う。
- 監査の再現：`-evalue` 5000・1e4・1e5・1e10 で全て一致（19033〜69454 行）。
