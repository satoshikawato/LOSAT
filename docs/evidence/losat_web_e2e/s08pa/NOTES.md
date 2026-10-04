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

## RP-4：60000 残基の query の TBLASTN の bit score

### 最初に値が食い違う箇所

再現：`audit/round1/rp4/` の `bigq.faa`（60000 残基）× `bigs4.fna`（sbig4 だけ）、既定の option、`-outfmt 6`。最初の HSP の bit score が NCBI 10383、LOSAT 10404（5 subject の `bigs.fna` では NCBI 10404、LOSAT 10412）。

1. TLOSAN Stage D の composition の shim（`ncbi_kappa_composition_trace.c`）で、最初の `Blast_AdjustScores` の subject の長さ（窓の長さ）が NCBI 11000（sbig4 の frame の全長）、LOSAT 5649。呼び出しの数も NCBI 330、LOSAT 143。窓は kappa に入る HSP から作る（redo_alignment.c:683-731、`kWindowBorder` 200）ので、HSP の一覧が既に違う。
2. 原因は予備の段の query の分割：NCBI は gapped の検索の query の batch を `SplitQuery_CalculateNumChunks`（split_query_aux_priv.cpp:73-138）で分ける。chunk の大きさは tblastn 20000、blastp 10000（local_blast.cpp:74-76,91-94）、重なり 100（split_query_aux_priv.cpp:53-60）。`SplitQuery_ShouldSplit` は blastp・tblastn で真（73-97）。60000 残基は 3 つの chunk になり、chunk ごとの予備の検索の HSP を合わせて traceback と composition の調整に進む（prelim_stage.cpp:232-233）。
3. 確かめ：NCBI を `CHUNK_SIZE=100000`（分割しない）で動かすと、LOSAT と完全に一致する（`bigs4.fna` 8 行、`bigs.fna` 7 行）。BLASTP も同じ：`bigq.faa` × `e2e_many_subject.faa`、`-evalue 1000` で、NCBI の分割あり（136 行）と分割なし（126 行）が 11 行違い、LOSAT は分割なしと一致（TBLASTN は `e2e_many_subject.fna` で 78 行中 6 行、LOSAT は分割なしの 72 行と一致）。
4. LOSAT は BLASTN の分割（S07++、`blastn/query_split.rs`）を移植したが、BLASTP と TBLASTN には無い。`check_query_split_environment` の注釈（「blastp と tblastn は protein の query を分割しない」）は誤りだった。
5. 認証の範囲：TLOSAN Stage G と E2 の fixture の protein の query の batch（20000 残基ごと）は最大 25,086 残基（`TrcuMJNV.faa`）で、TBLASTN の閾値 39,800 に届かない。BLASTP の batch（10000 残基ごと）は 16,729 残基未満で、閾値 19,800 に届かない。RP-4 は認証の範囲の外。

### 直し方（query の分割の移植）

最初は分割される batch を明示的に拒否したが、NCBI の分割が結果を変えない入力（BLASTP の 30000 残基の query など、以前は NCBI と一致していたもの）まで拒否するので、保守者の指示（2026-10-04）で移植に替えた。

TBLASTN（`tblastn/stage_d_pipeline.rs` の `split_preliminary_hitlists`、`stage_d_results.rs` の `merge_query_chunk`、`common/protein_query_split.rs`）：

1. chunk の作成（split_query_cxx.cpp:145-171,196-247,415-418）：`split_protein_batch`。chunk の範囲は BLASTN の移植（`blastn/query_split.rs`）と同じで、protein の query は 1 query 1 context（frame 0）。chunk の部分の offset は plus の側の式（split_query_cxx.cpp:630-643、`GetStartingChunk` split_query_aux_priv.cpp:266-284）。`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE` は NCBI と同じく読む（`app::protein_query_split_sizes`）。
2. chunk の設定（`SplitQuery_CreateChunkData` split_query_aux_priv.cpp:185-210）：chunk の部分配列を query の集合として、全体と同じ設定の関数（`query_set_setup`、非分割の経路からくくり出した）で、SEG・Karlin の値・lookup・cutoff を chunk について計算する。lower-case の mask は全体の query の mask を部分に制限したもの（`RestrictToSeqInt` の癖で 1 残基長い、`protein_chunk_part`）。有効な探索空間は、全体の batch で計算した値を option として与える（`SplitQuery_SetEffectiveSearchSpace` split_query_aux_priv.cpp:150-183）。chunk の `BLAST_CalcEffLengths` はその値を chunk の context の番号で引き（`s_GetEffectiveSearchSpaceForContext` blast_setup.c:676-697。値が 1 つなら全ての context が `[0]`）、length adjustment は chunk の部分の長さで計算する（`local_subject_effective_lengths_with_search_spaces`）。複数の query の batch では、後ろの chunk の context 0 が前の query の探索空間を使う（NCBI の決まった結果として再現）。全ての context が Karlin の計算に失敗する chunk は、NCBI が例外を無視するので何も見つけない（prelim_stage.cpp:270-283）。
3. chunk の予備の段：非分割と同じ subject の loop（`run_preliminary_subjects`、くくり出した）で、chunk の linking と予備の e-value の reap、chunk の予備の hit list（大きさは全体と同じ）。
4. 合わせ方（`BlastHSPStreamMerge` blast_hspstream.c:399-534）：HSP を全体の query に移し（query の offset に chunk の offset を足す）、`Blast_HitListMerge`（blast_hits.c:2132-2217。chunk の hit list の大きさの新しい hit list に OID の順で `Blast_HitListUpdate`、同じ OID は offset が正なら query 側の `Blast_HSPListsMerge` blast_hits.c:2857-3035、そうでなければ `Blast_HSPListAppend`）。最後に全ての HSP list を得点の順に並べる（blast_hspstream.c:527-534）。`s_BlastMergeTwoHSPs` は subject の chunk の併合の移植（`merge_two_chunk_hsps`）を使い、gap を許す（`GetGappedMode`、split_query_cxx.cpp:864-865）。
5. traceback と kappa は今までどおり全体の query と全体の統計で行う（prelim_stage.cpp:286-296）。
6. NCBI が落ちる組は明示的な拒否：chunk をもう一度分けるほど重なりが大きい組（NCBI は null の参照で止まる）、batch を分ける負の `CHUNK_SIZE` の組（BLASTN の DW-16 と同じ扱い）。

BLASTP は同じ移植を続けて行う（下）。

### 確かめ（TBLASTN）

- RP-4 の入力：`bigq.faa` × `bigs4.fna`・`bigs.fna`・`e2e_many_subject.fna`、既定・`-comp_based_stats 0`・`-evalue 1000` の 9 組が全て NCBI と一致（最大 802 行）。
- 小さい `CHUNK_SIZE`（1000、500、300/重なり 50、2000/重なり 300、700/重なり 0、1500/重なり 10）：`e2e_protein_query.faa`・`e2e_many_query.faa` × `e2e_tblastn_subject.fna`・`e2e_many_subject.fna`、`gen_q_lc.faa` × `gen_s.fna`（`-lcase_masking`）、cbs 0/2、sum stats の有無、`-evalue 1000`、少ない `-max_target_seqs`、`-soft_masking true`、BLOSUM45、outfmt 0/6/7 の 53 組で標準出力が全て一致（`-num_threads 4` の 3 組は NCBI のスレッドの警告だけが違う、承認済みの例外 1）。
- `CHUNK_SIZE=300 OVERLAP_CHUNK_SIZE=250` と `CHUNK_SIZE=-5 OVERLAP_CHUNK_SIZE=-10`：NCBI は CCoreException（null の参照）、LOSAT は明示的な拒否。
- fixture：`e2e.tblastn.query_split`（outfmt 6）、`query_split_fmt0`、`query_split_two`（19000 と 21000 残基の 2 query、`-comp_based_stats 0 -evalue 1000`）。変更前の実行ファイルは 3 件とも違う。
- TBLASTN の単体試験 110 件、分割の単体試験 3 件、TBLASTN の fixture 20 件。

## `-out -version`（BP-8 の残り）

### NCBI の振る舞い

`CNcbiApplicationAPI::AppMain` の argv の前処理（ncbiapp.cpp:926-1001）は、`--`（`s_ArgDelimiter`、ncbiargs.cpp:89）より前の全ての語を、他の option の値の位置でも調べる：`-version`・`-version-full`・`-version-full-xml`・`-version-full-json` は version を出して終了 0、`-dryrun` は argv から取り除く（`-out -dryrun` は `-out` の値が無くなり USAGE）、`-logfile`・`-conffile` は次の語を値として取る。LOSAT は `-out -version` の `-version` を `-out` の値として読み、`-version` という名のファイルを作っていた。

### 直し方（明示的な拒否）

`cli.rs` の `try_parse_from` が、program の名の後の語を NCBI と同じく `--` まで調べ（`ncbi_preparsed_toolkit_word`）、最初に見つかった上の語（と `-conffile=…`）を、option の位置と同じ文言で拒否する（`unknown_option_error` に既存の分岐をくくり出して共有。`-version` は「the NCBI BLAST+ option -version is not supported by LOSAT's <PROGRAM>」、他は「the NCBI C++ Toolkit option …」、終了コード 2）。BLASTN・BLASTP・TBLASTN・TBLASTX（BLASTX は DW-10 で変えない）。アダプタの `validate` も同じ `try_parse_from` で読むので、同じ拒否になる。決定 D9（toolkit の option と `-version` は拒否）の範囲で、`-version` を出す NCBI と同じにはしない。

試験 `tests/cli_v2.rs` の `toolkit_words_in_an_option_value_are_rejected_as_options`（4 program × 5 つの組）。`-out o.txt` は今までどおり。
