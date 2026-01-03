# TBLASTX NCBI Parity Status Report

**作成日時**: 2026-01-03  
**更新日時**: 2026-01-03  
**現象**: LOSATがNCBI BLAST+より多くのヒットを出力する  
**目標**: 出力を1ビットの狂いもなく一致させる

---

## 1. 修正完了済み (Completed Fixes)

### 1.1 DUST フィルタリングの削除
- **状態**: ✅ 完了
- **内容**: TBLASTXはNCBIでSEGのみ使用。LOSATからDUST関連コードを削除済み。
- **ファイル**: `args.rs`, `utils.rs`, `lookup.rs`

### 1.2 diag_offset の更新式
- **状態**: ✅ 完了
- **内容**: `diag_offset += s_aa_len + window` に修正済み。オーバーフロー処理 (`INT4_MAX/4` リセット) も実装済み。
- **ファイル**: `utils.rs`

### 1.3 Sum-Statistics Linking の主要ロジック
- **状態**: ✅ 完了
- **確認済み項目**:
  - グルーピングキー `(q_idx, s_idx, q_strand, s_strand)` ✅
  - ソート順序 (`s_RevCompareHSPsTbx` と同等 - reverse query position) ✅
  - `lh_helper` 配列構造の再現 ✅
  - `next_larger` skip-list 最適化 ✅
  - `linked_to` カウンター管理 ✅
  - `changed` フラグ管理 ✅
  - `gap_prob` / `(1 - gap_prob)` 適用 (num > 1 の場合) ✅
  - `cutoff_small_gap` / `cutoff_big_gap` 計算 (NCBI `CalculateLinkHSPCutoffs` ポート) ✅
  - E-valueのチェイン全体への適用 ✅
- **ファイル**: `sum_stats_linking.rs`

---

## 2. 修正が必要と判明している点 (Known Required Fixes)

### 2.1 ⚠️ X-drop 計算: 固定値 vs 動的計算
- **状態**: ⚠️ 要修正 (CRITICAL)
- **問題**: 
  - **LOSAT**: `constants.rs` で `X_DROP_UNGAPPED: i32 = 16` として固定値を使用
  - **NCBI**: `blast_parameters.c:219-221` で **Lambda を使って動的に計算**:
    ```c
    p->cutoffs[context].x_dropoff_init =
        (Int4)(sbp->scale_factor *
               ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    ```
  - NCBI は `BLAST_UNGAPPED_X_DROPOFF_PROT = 7` (bits) を入力とし、Lambda で raw score に変換
  - BLOSUM62 (Lambda ≈ 0.3176): `ceil(7 * 0.693 / 0.3176) = ceil(15.27) = 16`
  - **現状LOSATはたまたま16で一致しているが、異なる scoring matrix や異なる Lambda では不一致となる**
- **NCBIコード場所**: 
  - `blast_parameters.c:219-221` (初期化)
  - `blast_parameters.c:380-383` (per-subject 更新)
- **修正方針**: 
  - `utils.rs` の extension 呼び出し時に `ceil(X_DROP_BITS * LN2 / lambda)` を計算して渡す
  - または初期化時に ungapped_params.lambda を使って計算

### 2.2 ⚠️ Per-Subject Cutoff Score 更新
- **状態**: ⚠️ 要修正
- **問題**: NCBIは `BlastInitialWordParametersUpdate` でサブジェクト長に基づき cutoff を再計算
- **NCBIコード**: `blast_parameters.c:279-419`
  ```c
  Int2 BlastInitialWordParametersUpdate(EBlastProgramType program_number, 
     const BlastHitSavingParameters* hit_params, 
     const BlastScoreBlk* sbp, 
     BlastQueryInfo* query_info, Uint4 subj_length,
     BlastInitialWordParameters* parameters)
  {
     // ...
     BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp, searchsp, TRUE, gap_decay_rate);
     // ...
     new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
     curr_cutoffs->cutoff_score = new_cutoff;
     curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;
  }
  ```
- **LOSATの現状**: `ncbi_cutoffs.rs:compute_tblastx_cutoff_score` で計算しているが、per-context 更新ロジックが NCBI と完全一致するか未確認
- **修正方針**: NCBI と同様に各サブジェクト処理前に cutoff を再計算

### 2.3 ⚠️ X-dropoff の Per-Context 適用
- **状態**: ⚠️ 要修正
- **問題**: NCBIは context ごとに `cutoffs->x_dropoff` を持ち、extension 時にそれを参照
- **NCBIコード**: `aa_ungapped.c:579`
  ```c
  cutoffs = word_params->cutoffs + curr_context;
  score = s_BlastAaExtendTwoHit(matrix, subject, query,
                                last_hit + wordsize,
                                subject_offset, query_offset,
                                cutoffs->x_dropoff, ...);  // <-- context-specific
  ```
- **LOSATの現状**: `utils.rs:833` で固定の `dropoff` を渡している
  ```rust
  let dropoff = X_DROP_UNGAPPED;  // line 419
  // ... (line 833 で extend_hit_two_hit に dropoff を渡す)
  ```
- **修正方針**: context ごとの x_dropoff を計算し、extension 呼び出し時に正しい値を使用

---

## 3. "Might Need Adjustments" レベルの相違点

### 3.1 🔶 Sentinel バイト値の違い
- **状態**: 🔶 要検証 (低優先度)
- **問題**: 
  - **NCBI**: `NULLB = 0` を sentinel として使用 (`blast_encoding.c:120`)
  - **LOSAT**: `SENTINEL_BYTE = 255` を使用 (`constants.rs:98`)
- **影響分析**: 
  - 両者とも BLOSUM62 マトリックスで対応する残基がないため、`defscore = -4` が返る
  - LOSAT `extension.rs:24-29` では sentinel を明示的にチェックして `SENTINEL_PENALTY = -4` を返す
  - **機能的には同等のはず**だが、NCBI の extension コードでは `matrix[q[i]][s[i]]` で直接参照
- **調査状況**: 
  - NCBI `aa_ungapped.c:847` のマトリックス参照を確認
  - NCBI は NULLB (0) がマトリックス範囲外なら負のスコアを返す想定
- **結論**: 両者で X-drop による終了タイミングが同じなら問題なし。X-drop が正しく計算されていれば影響なし。

### 3.2 🔶 Frame Base 計算の Sentinel 考慮
- **状態**: 🔶 部分修正済み・要確認
- **変更内容**: `lookup.rs` で `base += frame.aa_seq.len().saturating_sub(1)` に修正
- **問題**: NCBI は frame 間に sentinel 1バイトを挿入し、次のフレームの開始は `+length+1`
- **調査状況**: 
  - NCBI `blast_setup.c` / `blast_util.c` の `BLAST_GetTranslation` を確認したが、LOSAT との完全比較は未完了
  - LOSAT は frame ごとに独立した aa_seq を持つため、NCBI の concatenated buffer とは構造が異なる
- **影響**: sum_stats_linking での abs_coords 計算に影響する可能性

### 3.3 🔶 HSP ソート順序の細部
- **状態**: 🔶 確認済み (問題なし)
- **NCBIコード** (`link_hsps.c:331-375`):
  ```c
  static int s_RevCompareHSPsTbx(const void *v1, const void *v2) {
     // 1. context/(NUM_FRAMES/2) で比較 (query strand+index)
     // 2. SIGN(subject.frame) で比較 (subject strand)
     // 3. query.offset descending
     // 4. query.end descending
     // 5. subject.offset ascending
     // 6. subject.end ascending
  }
  ```
- **LOSATコード** (`sum_stats_linking.rs:517-524`):
  ```rust
  group_hits.sort_by(|a, b| {
      bqo.cmp(&aqo)           // query offset descending
          .then(bqe.cmp(&aqe)) // query end descending
          .then(bso.cmp(&aso)) // subject offset ascending? (NO! should be aso.cmp(&bso))
          .then(bse.cmp(&ase)) // subject end ascending? (NO! should be ase.cmp(&bse))
  });
  ```
- **問題発見**: **subject の比較が逆順になっている!**
  - NCBI: `if (h1->subject.offset < h2->subject.offset) return 1;` → ascending order
  - LOSAT: `bso.cmp(&aso)` → descending order
  - **これは修正が必要**

### 3.4 🔶 E-value 計算の丸め処理
- **状態**: 🔶 確認済み (問題なし)
- **問題**: E-value から cutoff score への変換で NCBI は ceiling を使用
- **NCBIコード** (`blast_stat.c:4049-4063`):
  ```c
  S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ));
  ```
- **LOSATコード** (`ncbi_cutoffs.rs:152`):
  ```rust
  let score = ((gapped_params.k * searchsp / e).ln() / gapped_params.lambda).ceil();
  ```
- **結論**: ✅ 一致している

### 3.5 🔶 Gap Trigger スコアの計算
- **状態**: 🔶 確認済み (問題なし)
- **NCBIコード** (`blast_parameters.c:343-344`):
  ```c
  gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda);
  ```
- **LOSATコード** (`ncbi_cutoffs.rs:48`):
  ```rust
  let raw = (bit_trigger * NCBIMATH_LN2 + ungapped_params.k.ln()) / ungapped_params.lambda;
  raw as i32
  ```
- **結論**: ✅ 一致している (truncation = `as i32`)

### 3.6 🔶 Extension 終了条件
- **状態**: 🔶 確認済み (問題なし)
- **NCBIコード** (`aa_ungapped.c:859`):
  ```c
  if (score <= 0 || (maxscore - score) >= dropoff)
      break;
  ```
- **LOSATコード** (`extension.rs:149`):
  ```rust
  if right_score <= 0 || (max_score_total - right_score) >= x_drop {
      break;
  }
  ```
- **結論**: ✅ 一致している

### 3.7 🔶 Sum-Statistics の effective length 計算
- **状態**: 🔶 要確認
- **問題**: NCBI `link_hsps.c:560-571` で length_adjustment を context から取得
  ```c
  query_context = hp_start->next->hsp->context;
  length_adjustment = query_info->contexts[query_context].length_adjustment;
  query_length = query_info->contexts[query_context].query_length;
  query_length = MAX(query_length - length_adjustment, 1);
  subject_length = MAX(subject_length - length_adjustment, 1);
  ```
- **LOSATコード** (`sum_stats_linking.rs:532`):
  ```rust
  let search_space = SearchSpace::with_length_adjustment(query_len_aa, subject_len_aa, params);
  ```
- **調査状況**: `SearchSpace::with_length_adjustment` が NCBI と同じ計算をしているか確認が必要

---

## 4. 調査未着手の領域

### 4.1 ❓ Two-hit Window の詳細
- **状態**: ❓ 未調査
- **概要**: 2ヒット法の window / threshold 処理が NCBI と完全一致するか
- **関連NCBIコード**: `aa_ungapped.c:380-398`
  ```c
  diff = subject_offset - last_hit;
  if (diff >= window_size) {
      diag_array[diag_coord].last_hit = subject_offset + diag_offset;
      continue;
  }
  if (diff < wordsize) {
      continue;
  }
  ```
- **LOSAT**: `utils.rs` で同様のロジックを実装しているが line-by-line 比較は未実施

### 4.2 ❓ Lookup Table 構築の詳細
- **状態**: ❓ 未調査
- **概要**: Lookup table のワードサイズ、threshold 処理が NCBI と完全一致するか
- **関連NCBIコード**: `aa_lookup.c`

### 4.3 ❓ Masked Region の Extension 時処理
- **状態**: ❓ 未調査
- **概要**: SEG でマスクされた領域の extension 時の処理が NCBI と一致するか
- **関連NCBIコード**: `blast_seg.c`, `blast_filter.c`
- **LOSAT**: 
  - `utils.rs:481-492` でマスクされた残基を `X (21)` に置換
  - Extension 時にスコアが低くなり自然に終了する想定

### 4.4 ❓ HSP の重複排除 (Culling)
- **状態**: ❓ 未調査
- **概要**: HSP 間の重複排除ロジックが NCBI と一致するか
- **関連NCBIコード**: `link_hsps.c` の culling 関連関数

---

## 5. 発見した明確なバグ

### 5.1 🐛 HSP ソートの subject 比較順序が逆
- **ファイル**: `sum_stats_linking.rs:517-524`
- **問題**: 
  ```rust
  // 現在のコード (間違い)
  .then(bso.cmp(&aso)) // subject offset: descending
  .then(bse.cmp(&ase)) // subject end: descending
  
  // NCBI の正しい順序
  .then(aso.cmp(&bso)) // subject offset: ascending
  .then(ase.cmp(&bse)) // subject end: ascending
  ```
- **影響**: リンキングの順序が変わり、結果が異なる可能性

---

## 6. 推定される根本原因

LOSATがNCBIより多くのヒットを出力する原因として、以下が推定される:

1. **X-drop / Cutoff の不整合**: 
   - X-drop が固定値で、context-specific な値を使っていない
   - Per-subject cutoff 更新が NCBI と異なる可能性

2. **HSP ソート順序のバグ**: 
   - Subject offset/end の比較順序が逆
   - リンキング結果が変わり、E-value 計算に影響

3. **Sum-Statistics の Length Adjustment**: 
   - effective length の計算が NCBI と異なる可能性
   - Search space が異なれば E-value も異なる

---

## 7. 優先度順の修正作業リスト

| 優先度 | ID | 内容 | 推定工数 | ファイル |
|--------|-----|------|----------|----------|
| **1** | 5.1 | HSP ソート順序 subject ascending 修正 | 小 | `sum_stats_linking.rs` |
| **2** | 2.1 | X-drop 動的計算 (ceil(7*LN2/Lambda)) | 小 | `utils.rs`, `constants.rs` |
| **3** | 2.3 | X-drop の per-context 適用 | 中 | `utils.rs` |
| **4** | 2.2 | Per-subject cutoff 更新ロジック確認 | 中 | `utils.rs`, `ncbi_cutoffs.rs` |
| **5** | 3.7 | Sum-stats effective length 計算確認 | 小 | `sum_stats_linking.rs` |
| **6** | 3.2 | Frame base 計算の再検証 | 小 | `lookup.rs` |
| **7** | 4.1-4.4 | 未調査領域の調査 | 大 | 各種 |

---

## 8. 次のアクション

1. **🐛 HSP ソート順序を修正** (`sum_stats_linking.rs:517-524`)
2. **⚠️ X-drop を動的計算に変更** (`utils.rs`)
3. **⚠️ Per-context x_dropoff を extension に渡す** (`utils.rs`)
4. **差分確認テストを実行し、残存差異を特定**
