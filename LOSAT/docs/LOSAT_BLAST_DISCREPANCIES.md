## LOSAT BLAST Discrepancies (NCBI Parity Audit)

### 目的

`TBLASTX_NCBI_PARITY_DEVLOG.md` の **直近5セッション**で行われた修正について、
**NCBI BLAST 実装（ground truth）**: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/`
と LOSAT 実装を突き合わせて検証し、**非準拠（出力差分につながる/仕様解釈が違う/移植が未完了/誤植や誤読）** を網羅的に列挙する。

### 運用ルール

- **NCBI 実装が唯一の正解**。ここに書くのは「NCBI と一致していない点」のみ。
- 推測は避け、可能な限り
  - **NCBI 側の該当ファイル/関数/行**
  - **LOSAT 側の該当ファイル/関数**
  - **差分の根拠（条件分岐/定数/データ構造差による論理差）**
  を明記する。
- 不確実なものは **Needs verification** として別枠に隔離し、確証が取れた時点で確定セクションへ移す。

### 対象セッション（devlogの直近5セッション）

devlog 末尾から（ファイル順）:

- `2026-01-02 追記: sum_stats_linking NCBI完全準拠修正`
- `2026-01-01 セッション: 侵入型リンクリスト実装と NCBI パリティ修正`
- `2026-01-02 追記: best 選択バグの修正と can_skip 最適化の安全化`
- `2026-01-02 セッション終了: sum_stats_linking の根本的なバグ発見`
- `2026-01-02 追記: 重大バグ発見と修正 - extension で masked sequence を使用していた`

---

## Discrepancies（確定）

### 1) ~~devlogの主張がNCBI実装と矛盾: 「extensionはunmasked(sequence_nomask)を使う」~~ **[修正完了: 2026-01-02 Plan D]**

- **devlog**: `2026-01-02 追記: 重大バグ発見と修正 - extension で masked sequence を使用していた` にて
  - 「LOSAT は extension 時に SEG masked された `aa_seq` を使用していた。NCBI は `sequence_nomask` (unmasked) を使う」
  - という主張＋ `aa_seq_nomask` を使う修正案が書かれている。
- **NCBI ground truth**:
  - マスク設定で `query_blk->sequence_start_nomask` / `query_blk->sequence_nomask` を作るのは事実だが、これは **非マスク配列の保持**であり（例: identity 計算や出力用途）、
    extension で参照される配列は **`query->sequence`（マスク済み）**。
  - `aa_ungapped.c` の `s_BlastAaExtendTwoHit()` は `Uint1 *q = query->sequence;` を参照し、`matrix[q[q_right_off + i]]...` を評価する（`query->sequence_nomask` は参照しない）。
- **NCBI refs**:
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1350-1406`
    - `BlastSetUp_MaskQuery()` が `sequence_start_nomask` を作りつつ、マスクは `query_blk->sequence`（working buffer）に適用している。
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158`
    - `Uint1 *q = query->sequence;` を用いて ungapped extension を実施。
- **LOSAT現状（確認済み 2026-01-02）**:
  - `LOSAT/src/algorithm/tblastx/utils.rs` の normal/neighbor-map の two-hit extension は `aa_seq`（masked working buffer）を渡している。
  - `aa_seq_nomask` は identity 計算（出力用）にのみ使用されている。
  - これは **NCBI と完全に整合**する。
- **修正内容（Plan D）**:
  - devlog の誤記を「⚠️ CORRECTION」セクションとして訂正（`TBLASTX_NCBI_PARITY_DEVLOG.md` 2026-01-02 追記部）
  - NCBI 根拠コード（`aa_ungapped.c`, `blast_filter.c`）を verbatim 引用
  - 現在の LOSAT 実装（extension=masked, nomask=identity）が正しいことを確認

### 2) ~~`sum_stats_linking` の link cutoffs 計算が NCBI `CalculateLinkHSPCutoffs()` と一致していない~~ **[修正完了: 2026-01-02]**

- **対象**:
  - LOSAT: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs:calculate_link_hsp_cutoffs_ncbi()` (新規追加)
  - NCBI: `blast_parameters.c:CalculateLinkHSPCutoffs()`
- **修正内容（Plan A + Plan B 実施）**:
  - **query_length の算出方法**: ✅ `compute_avg_query_length_ncbi()` を追加し NCBI 式を再現
  - **`cutoff_small_gap` の下限**: ✅ `cutoff_score_min.max(cutoff_small_raw)` で NCBI と同一
  - **`sbp->scale_factor` の適用**: ✅ `scale_factor` パラメータを追加、cutoffs に乗算
  - **`gap_prob` の扱い**: ✅ `LinkHspCutoffs` で `gap_prob` を返却し linking に伝播、small search space で 0 に設定
  - **`y_variable` の計算分岐**: ✅ `db_length > subject_length` 分岐を実装（-subject モードでは else 分岐）
  - **expected length の丸め**: ✅ `blast_nint()` 関数を追加し NCBI の `BLAST_Nint` を再現
- **NCBI refs**:
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082`
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:437-441`
- **NCBI 完全一致確認（2026-01-02 検証）**:
  - `blast_nint()`: NCBI `BLAST_Nint` と同一ロジック（`x + 0.5` for positive, `x - 0.5` for negative → trunc）
  - `compute_avg_query_length_ncbi()`: NCBI `s_QueryInfo_SetContext` のオフセット構築と平均化を再現
  - `calculate_link_hsp_cutoffs_ncbi()`: NCBI C コードをコメントとして引用、verbatim port
  - `gap_prob` 伝播: NCBI の `link_hsp_params->gap_prob = 0` を `LinkHspCutoffs.gap_prob` で再現
  - subject 単位 linking: neighbor-map モードで hits を `s_idx` 単位に分割して適用

### 3) ~~`Blast_HSPListPurgeHSPsWithCommonEndpoints` 相当処理が NCBI と同等になっていない（適用単位/キー）~~ **[修正完了: 2026-01-02]**

- **対象**:
  - LOSAT: `LOSAT/src/algorithm/tblastx/utils.rs:purge_hsps_with_common_endpoints()`
  - NCBI: `blast_hits.c:Blast_HSPListPurgeHSPsWithCommonEndpoints()` + `s_QueryOffsetCompareHSPs()` / `s_QueryEndCompareHSPs()`
- **修正内容（Plan C 実施）**:
  - **比較キーを NCBI 完全一致に置換**:
    - Phase 1 sort: `ctx_idx ASC → q_aa_start ASC → s_aa_start ASC → raw_score DESC → q_aa_end DESC → s_aa_end DESC`
    - Phase 1 duplicate key: `(ctx_idx, q_aa_start, s_aa_start, s_frame)`
    - Phase 2 sort: `ctx_idx ASC → q_aa_end ASC → s_aa_end ASC → raw_score DESC → q_aa_start DESC → s_aa_start DESC`
    - Phase 2 duplicate key: `(ctx_idx, q_aa_end, s_aa_end, s_frame)`
    - NCBI C コードをコメントとして引用（`blast_hits.c` lines 2267-2387, 2482-2513）
  - **適用単位を BlastHSPList 相当に修正**:
    - neighbor-map mode: `all_ungapped` への混在 purge を撤去
    - normal mode: もともと subject 単位で閉じているため変更なし
  - **NCBI parity 発見**: ungapped tblastx では NCBI は purge を**呼ばない**
    - `blast_engine.c` line 545 の purge は `if (score_options->gapped_calculation)` 内のみ
    - tblastx は ungapped のため、この purge は実行されない
    - LOSAT も同様に purge を呼ばないことで NCBI parity を達成
  - **テスト追加**:
    - `test_purge_ncbi_parity`: NCBI `testCheckHSPCommonEndpoints` 相当
    - `test_purge_different_s_frame_not_duplicate`: s_frame は重複判定にのみ使用
    - `test_purge_different_context_not_duplicate`: context は重複判定に使用
    - `test_mixed_subject_purge_is_wrong`: 混在 purge 禁止の回帰テスト
- **NCBI refs**:
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2534`
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2400`（比較関数）
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:545`（gapped-only 呼び出し）

### 4) ~~devlogの「`--max-hsps-per-subject` を追加」記述に対して、実装が未接続（dead option）~~ **[修正完了: 2026-01-02 Plan E]**

- **devlog**: `2026-01-02 追記: sum_stats_linking NCBI完全準拠修正` にて
  - `--max-hsps-per-subject` を追加した旨の記述があった。
- **問題**:
  - `LOSAT/src/algorithm/tblastx/args.rs` に定義はあったが、実装に接続されていなかった（dead option）。
- **修正内容（Plan E）**:
  - `args.rs` から `--max-hsps-per-subject` オプションを削除
  - `TBLASTX_NCBI_PARITY_DEVLOG.md` の該当セクションを「削除済み (Plan E)」と注記
  - `NEXT_SESSION_PROMPT_2026-01-02_SUM_STATS_LINKING_PERF.md` から関連タスク/参考実装を除去
  - NCBI parity 経路を汚さないため、非パリティ最適化は将来必要時に再実装すること

### 5) devlogの「`ordering_method` を保持する」主張と現コードが一致していない（HSPチェーン復元の欠落の可能性）

- **devlog**: `2026-01-01 セッション: 侵入型リンクリスト実装と NCBI パリティ修正` にて
  - `HspLink` に `ordering_method` を追加し、チェーン抽出時に `ordering` を保存する旨が書かれている。
- **NCBI ground truth**:
  - `link_hsps.c` は `LinkHSPStruct.ordering_method` を保持し（line 973 付近）、最終的に HSP を `next/prev` で「チェーン順」に接続する処理で参照する（line 1027 付近）。
- **LOSAT現状（確定）**:
  - `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` の `HspLink` に `ordering_method` は存在しない。
  - `next/prev` での「NCBI同等のチェーン復元」も実装されていない（現状は `Vec<UngappedHit>` に E-value を付与して返すのみ）。
- **影響（Needs verification）**:
  - LOSAT の出力仕様が「NCBI の HSP チェーン出力順/連結」を要求する場合、**出力順や HSP セット構造が一致しない**要因になりうる。

### 6) ~~normal mode の ungapped cutoff (`cutoff_score`) 計算が NCBI `word_params->cutoffs[curr_context].cutoff_score` と一致していない~~ **[修正完了: 2026-01-02]**

- **対象**:
  - LOSAT: `LOSAT/src/algorithm/tblastx/utils.rs`（normal mode の scan/extension 前に `cutoff_scores` を自前計算）
  - NCBI: `aa_ungapped.c:s_BlastAaWordFinder_TwoHit()` 内の `cutoffs = word_params->cutoffs + curr_context;` と `if (score >= cutoffs->cutoff_score)`  
    + `blast_parameters.c:BlastInitialWordParametersUpdate()`（`cutoff_score` を context ごとに更新）
- **修正内容（Plan 0 実施）**:
  - `LOSAT/src/algorithm/tblastx/ncbi_cutoffs.rs` を新規追加し、NCBI C コードを verbatim port:
    - `gap_trigger_raw_score()`: ungapped params (kbp_std) を使用、trunc 丸め
    - `compute_eff_searchsp_subject_mode_tblastx()`: gapped params (kbp_gap) を使用、length adjustment 計算
    - `cutoff_score_from_evalue()`: gapped params、ceil 丸め、kSmallFloat=1e-297 clamp
    - `compute_tblastx_cutoff_score()`: 上記を組み合わせて `MIN(gap_trigger, cutoff_score_max)` を計算
  - `LOSAT/src/stats/tables.rs` に `lookup_protein_params_ungapped()` を追加
  - `LOSAT/src/algorithm/tblastx/utils.rs` の normal mode / neighbor-map mode 両方で `compute_tblastx_cutoff_score()` を呼び出すように変更
- **NCBI refs**:
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:340-345`（gap_trigger は kbp_std を使用）
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:860-861`（cutoff_score_max は kbp_gap を使用）
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:368-374`（final cutoff = MIN(gap_trigger, cutoff_score_max)）
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4040-4063`（BlastKarlinEtoS_simple）
- **テスト結果**:
  - `test_gap_trigger_raw_score`: BLOSUM62 ungapped, bits=22 → **41** (trunc, not ceil)
  - `test_cutoff_score_cap_effective`: 小さい subject で cap が効くことを確認
  - `test_cutoff_score_gap_trigger_wins`: 大きい subject で gap_trigger が勝つことを確認
- **NCBI 完全一致確認（2026-01-02 再検証）**:
  - **gap_trigger 計算**: NCBI `(Int4)((gap_trigger * LN2 + logK) / Lambda)` と bit-perfect 一致
    - 丸め: `(Int4)` = trunc（LOSAT の `as i32` と同一）
    - パラメータ: `kbp_std` (ungapped) を使用
  - **cutoff_score_max 計算**: NCBI `BlastKarlinEtoS_simple` と bit-perfect 一致
    - 式: `ceil(log(K * searchsp / E) / Lambda)`
    - E clamp: `MAX(E, 1e-297)` = `e.max(K_SMALL_FLOAT)`
    - パラメータ: `kbp_gap` (gapped) を使用
  - **BLOSUM62 パラメータ**: NCBI `blast_stat.c:blosum62_values` と完全一致
    - ungapped: λ=0.3176, K=0.134, H=0.4012, α=0.7916, β=-3.2
    - gapped (11,1): λ=0.267, K=0.041, H=0.14, α=1.9, β=-30.0
  - **GAP_TRIGGER_BIT_SCORE**: NCBI `BLAST_GAP_TRIGGER_PROT = 22.0` と一致
  - **eff_searchsp 計算**: NCBI `db_length/3` + `BLAST_ComputeLengthAdjustment` と同一ロジック
  - **最終 cutoff_score**: NCBI `MIN(gap_trigger * scale_factor, cutoff_score_max)` と一致（scale_factor=1.0）
  - **結論**: 全計算式・パラメータ・丸め方法が NCBI C コードと bit-perfect に一致することを確認済み

---

## Needs verification（未確定）

- **gapped vs ungapped linking params**:
  - LOSAT `sum_stats_linking` は `GAP_PROB_UNGAPPED=0.5` / `GAP_DECAY_RATE_UNGAPPED=0.5` を固定で使用している。
  - NCBI は `BlastLinkHSPParametersNew(program, gapped_calculation)` により
    - ungapped: `BLAST_GAP_PROB=0.5`, `BLAST_GAP_DECAY_RATE=0.5`
    - gapped: `BLAST_GAP_PROB_GAPPED=1.0`, `BLAST_GAP_DECAY_RATE_GAPPED=0.1`
    を切り替える。
  - LOSAT が「NCBI default（gapped）出力」を目標にするなら、ここは一致していない可能性。
- **per-context Karlin params**:
  - NCBI `link_hsps.c` は `kbp[H->hsp->context]` を参照して `Lambda/logK` を context ごとに使用する。
  - LOSAT は `params` を単一で使用している。composition-based stats 等で差が出るか要確認。
- **出力順（determinism）**:
  - LOSAT は `apply_sum_stats_even_gap_linking()` を group 単位で並列処理し、最終出力は主に bit score で sort している。
  - NCBI の outfmt 出力順（query→subject→HSP の安定順）と一致するか要検証（同点 tie-break など）。

---

## 修正計画（本ドキュメントに基づく網羅計画）

### 方針

- **NCBI 実装が ground truth**。ロジック差分は「NCBI の該当 C コードを Rust の直上に引用」し、verbatim port で潰す。
- “性能のための拡張” は **デフォルトOFF** に隔離し、パリティ経路を汚さない（Strict Protocol）。

### 優先順位（NCBIパリティ実現のための順序）

1. **Plan 0**（Discrepancy #6）: **ungapped cutoff (`cutoff_score`) を NCBI の `word_params->cutoffs[curr_context].cutoff_score` に一致**
2. **Plan A + Plan B**（Discrepancy #2）: `sum_stats_linking` の link cutoffs（`CalculateLinkHSPCutoffs`）/`gap_prob` 伝播を完全一致
3. **Plan C**（Discrepancy #3）: endpoint purge の「適用単位/キー/タイミング」を NCBI と一致（または parity 経路から除外）
4. **Plan F（必要なら）**（Discrepancy #5 / Needs verification）: chain 復元・出力順が parity 条件なら実装
5. **Plan D / Plan E**（Discrepancy #1/#4）: 誤記・dead option の整理（将来の誤移植防止/再現性）

### Plan 0: ungapped cutoff (`cutoff_score`) を NCBI に完全一致させる（最優先） **[実施完了: 2026-01-02]**

- **対象ファイル**: `LOSAT/src/algorithm/tblastx/utils.rs`
- **対象箇所**: normal mode / neighbor-map mode の「score >= cutoff なら HSP を保存」の cutoff 生成部
- **実施内容**:
  - `LOSAT/src/algorithm/tblastx/ncbi_cutoffs.rs` を新規追加
  - NCBI `BlastInitialWordParametersUpdate()` の cutoff 計算ロジックを Rust へ verbatim port
  - ungapped params (kbp_std) と gapped params (kbp_gap) を明確に分離
  - `gap_trigger` は trunc 丸め（NCBI の `(Int4)` キャストを再現）
  - `cutoff_score_max` は ceil 丸め + kSmallFloat clamp
  - normal mode / neighbor-map mode 両方で `compute_tblastx_cutoff_score()` を使用
- **検証結果**:
  - 5 つの unit test すべて pass（gap_trigger=41, cap 動作確認済み）

### Plan A: `sum_stats_linking` の link cutoffs 計算を NCBI `CalculateLinkHSPCutoffs()` に完全一致させる **[実施完了: 2026-01-02]**

- **対象ファイル**: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`
- **対象関数**: `calculate_link_hsp_cutoffs_ncbi()` (新規追加)
- **実施内容**:
  - NCBI `blast_parameters.c:CalculateLinkHSPCutoffs()` (lines 998-1082) を **verbatim port**:
    - `query_length` の **平均化式**: `compute_avg_query_length_ncbi()` を新規追加
      - NCBI `s_QueryInfo_SetContext` のオフセット計算を再現
      - 式: `(last_offset + last_length - 1) / num_contexts`
    - `expected_length = BLAST_Nint(log(K*q*s)/H)`: `blast_nint()` 関数を追加
    - `db_length > subject_length` 分岐での `y_variable` 計算: 実装済み（-subject モードでは db_length=0 なので else 分岐）
    - small search space 時の `gap_prob=0`, `cutoff_small_gap=0`: `LinkHspCutoffs` 構造体で返却
    - `cutoff_small_gap = MAX(word_params->cutoff_score_min, ...)`: `cutoff_score_min` パラメータを追加
    - `cutoff_* *= sbp->scale_factor`: `scale_factor` パラメータを追加（現在は 1.0）
  - `LinkingParams` 構造体を追加し、linking に必要なパラメータをまとめて渡す
  - `utils.rs` で各 subject ごとに `cutoff_score_min` を計算（全 context の最小値）
  - neighbor-map モードでは subject 単位で linking を適用するよう修正
- **NCBI refs**:
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082`
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:437-441` (BLAST_Nint)
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp` (s_QueryInfo_SetContext)

### Plan B: `gap_prob` の取り扱いを NCBI の状態遷移どおりにする **[実施完了: 2026-01-02]**

- **対象**: `sum_stats_linking.rs`（prob 補正）
- **実施内容**:
  - `calculate_link_hsp_cutoffs_ncbi()` から `LinkHspCutoffs` 構造体で `gap_prob` を返却
  - small search space (`search_sp <= 8*window^2`) の場合、`gap_prob = 0` を設定
  - `apply_sum_stats_even_gap_linking()` で cutoffs を計算し、`gap_prob` を `link_hsp_group_ncbi()` に伝播
  - `link_hsp_group_ncbi()` 内の prob 補正で、伝播された `gap_prob` を使用:
    - `prob[0] /= gap_prob` (gap_prob=0 の場合は INT4_MAX)
    - `prob[1] /= (1-gap_prob)` (1-gap_prob=0 の場合は INT4_MAX)

### Plan C: endpoint purge を NCBI と同等にする（適用単位/キー/適用タイミング） **[実施完了: 2026-01-02]**

- **対象ファイル**: `LOSAT/src/algorithm/tblastx/utils.rs`
- **対象**: `purge_hsps_with_common_endpoints()`
- **実施内容**:
  - **comparator を NCBI 完全一致に置換** (`s_QueryOffsetCompareHSPs` / `s_QueryEndCompareHSPs`):
    - Phase 1 sort: `ctx_idx ASC → q_aa_start ASC → s_aa_start ASC → raw_score DESC → q_aa_end DESC → s_aa_end DESC`
    - Phase 2 sort: `ctx_idx ASC → q_aa_end ASC → s_aa_end ASC → raw_score DESC → q_aa_start DESC → s_aa_start DESC`
    - **注意**: `s_frame` はソートキーに含まれない（重複判定にのみ使用）
  - **重複判定キーを NCBI 一致**:
    - Phase 1: `(ctx_idx, q_aa_start, s_aa_start, s_frame)`
    - Phase 2: `(ctx_idx, q_aa_end, s_aa_end, s_frame)`
  - **混在 purge を撤去**: neighbor-map mode の `all_ungapped` への purge を削除
  - **NCBI parity 発見**: ungapped tblastx では NCBI は purge を**呼ばない**ため、LOSAT も呼ばない
  - **テスト追加**: NCBI `testCheckHSPCommonEndpoints` 相当 + 混在禁止回帰テスト（4件 pass）

### Plan D: devlog の誤りを修正（将来の誤移植を防ぐ） **[実施完了: 2026-01-02]**

- **対象**: `TBLASTX_NCBI_PARITY_DEVLOG.md`
- **実施内容**:
  - 「NCBI は extension に `sequence_nomask` を使う」という誤記を「⚠️ CORRECTION」セクションで訂正
  - NCBI 根拠コード（`aa_ungapped.c:1089-1094`, `blast_filter.c:1350-1406`）を verbatim 引用
  - 正しい NCBI 挙動: extension は `query->sequence`（masked）、`sequence_nomask` は identity 等の保持用
  - 現在の LOSAT 実装が NCBI と整合していることを確認

### Plan E: "非パリティ最適化" の整理（dead option を解消） **[実施完了: 2026-01-02]**

- **対象**: `LOSAT/src/algorithm/tblastx/args.rs` の `max_hsps_per_subject`
- **実施内容**:
  - `--max-hsps-per-subject` オプションを `args.rs` から削除
  - `TBLASTX_NCBI_PARITY_DEVLOG.md` に「削除済み (Plan E)」と注記
  - `NEXT_SESSION_PROMPT_2026-01-02_SUM_STATS_LINKING_PERF.md` から関連タスク/参考実装を除去
  - NCBI parity 経路を汚さないことを優先

### Plan F: HSP チェーン復元/出力順のパリティ要否を確定し、必要なら実装

- **対象**: `sum_stats_linking.rs`（`ordering_method`/chain order）と最終出力ソート
- **やること**:
  - NCBI の `link_hsps.c` が最終的に `next/prev` で HSP を連結する出力順を、LOSAT の outfmt 仕様で要求するか決める。
  - 要求するなら:
    - `ordering_method` を保持し、NCBI 同等に chain を “next で辿れる形” に復元するか、
    - もしくは outfmt の出力順ソートを NCBI と一致させる（tie-break 含む）。


