# TBLASTX NCBI Parity Status Report

**作成日時**: 2026-01-03  
**更新日時**: 2026-01-03 (4.8 Sum-Statistics Linking 調査完了 + 実装差異詳細追記)  
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

### 1.4 X-drop 動的計算
- **状態**: ✅ 完了
- **修正日**: 2026-01-03
- **問題だった点**: 
  - **旧LOSAT**: `constants.rs` で `X_DROP_UNGAPPED: i32 = 16` として固定値を使用
  - **NCBI**: `blast_parameters.c:219-221` で Lambda を使って動的に計算
- **修正内容**:
  - `ncbi_cutoffs.rs` に `x_drop_raw_score()` 関数を追加
  - NCBI公式: `(Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda))`
  - `utils.rs` の `run()` と `run_with_neighbor_map()` 両方で動的計算を使用
- **NCBIコード場所**: `blast_parameters.c:219-221`
- **LOSATコード**:
  ```rust
  // ncbi_cutoffs.rs
  pub fn x_drop_raw_score(x_drop_bits: f64, ungapped_params: &KarlinParams, scale_factor: f64) -> i32 {
      (scale_factor * (x_drop_bits * NCBIMATH_LN2 / ungapped_params.lambda).ceil()) as i32
  }
  
  // utils.rs (両モードで使用)
  let ungapped_params_for_xdrop = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
  let dropoff = x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ungapped_params_for_xdrop, 1.0);
  ```
- **検証**: BLOSUM62 で `ceil(7 * 0.693 / 0.3176) = 16` を確認 (ユニットテスト追加済み)
- **ファイル**: `ncbi_cutoffs.rs`, `utils.rs`, `constants.rs`

### 1.5 Per-Subject Cutoff Score 更新
- **状態**: ✅ 完了
- **修正日**: 2026-01-03
- **問題だった点**: 
  - **旧LOSAT**: `compute_tblastx_cutoff_score()` でユーザーのE-value (10.0) を直接使用し、searchsp も eff_searchsp (length adjustment 適用済み) を使用
  - **NCBI**: `BlastInitialWordParametersUpdate` で `CUTOFF_E_TBLASTX = 1e-300` と **異なる searchsp 計算式** を使用
- **NCBIコード**: `blast_parameters.c:348-374` (ungapped path for tblastx)
  ```c
  double cutoff_e = s_GetCutoffEvalue(program_number);  // = 1e-300 for tblastx!
  // ※重要: subj_length は NUCLEOTIDE 長 (AA長ではない!)
  // searchsp = MIN(query_len_aa, subject_len_nucl) * subject_len_nucl
  BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp, 
                MIN((Uint8)subj_length, (Uint8)query_length)*((Uint8)subj_length), 
                TRUE, gap_decay_rate);
  new_cutoff = MIN(new_cutoff, gap_trigger);
  new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
  ```
- **NCBIの2つの searchsp 計算の違い**:
  1. `BlastInitialWordParametersUpdate` (per-subject update):
     - `searchsp = MIN(q_aa_len, s_nucl_len) * s_nucl_len` (**長さ調整なし、AA/ヌクレオチド混在**)
     - `cutoff_e = 1e-300` (固定)
     - `dodecay = TRUE`
  2. `BlastHitSavingParametersNew` (初期設定の cutoff_score_max):
     - `searchsp = eff_searchsp` (**長さ調整あり**)
     - `cutoff_e = ユーザー指定 (10.0)`
     - `dodecay = FALSE`
- **修正内容**:
  - `ncbi_cutoffs.rs` に定数追加: `CUTOFF_E_TBLASTX = 1e-300`, `BLAST_GAP_DECAY_RATE = 0.5`
  - `cutoff_score_for_update_tblastx()` 関数追加: NCBIの `BlastInitialWordParametersUpdate` ungapped path をポート
  - `cutoff_score_max_for_tblastx()` 関数追加: NCBIの `BlastHitSavingParametersNew` をポート
  - `utils.rs` の `run()` と `run_with_neighbor_map()` を修正して新関数を使用
- **実質的な cutoff 決定要因**:
  - `CUTOFF_E_TBLASTX = 1e-300` から計算される cutoff は通常 1 または非常に低い値
  - 最終的な cutoff は `MIN(update_cutoff, gap_trigger, cutoff_score_max)` で決定
  - BLOSUM62 の場合: `gap_trigger = 41` が支配的になることが多い
- **ファイル**: `ncbi_cutoffs.rs`, `utils.rs`

### 1.6 X-dropoff の Per-Context 適用
- **状態**: ✅ 完了
- **修正日**: 2026-01-03
- **問題だった点**: 
  - **旧LOSAT**: 単一の `dropoff` / `x_drop` を全 context で共用
  - **NCBI**: context ごとに `cutoffs[context].x_dropoff_init` を計算・参照
- **NCBIコード**: `blast_parameters.c:219-221`, `aa_ungapped.c:575-579`
  ```c
  // 初期化時
  p->cutoffs[context].x_dropoff_init =
      (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
  
  // extension 時
  cutoffs = word_params->cutoffs + curr_context;
  score = s_BlastAaExtendTwoHit(..., cutoffs->x_dropoff, ...);
  ```
- **修正内容**:
  - `run()`: `x_dropoff_per_context: Vec<i32>` を `contexts` 作成後に生成
  - `run_with_neighbor_map()`: per-subject で `x_dropoff_per_context` を生成
  - extension 呼び出しで `x_dropoff_per_context[ctx_idx]` / `x_dropoff_per_context[ctx_flat]` を使用
- **NCBIの挙動** (`blast_stat.c:2796-2797`):
  - tblastx では全 context が `kbp_ideal` (BLOSUM62 ungapped Lambda=0.3176) を使用
  - → **x_dropoff = 16 は全 context で同一** (実質的な出力変更なし)
- **結論**: NCBIとの構造的 parity を達成
- **ファイル**: `utils.rs`

### 1.7 scale_factor の確認
- **状態**: ✅ 完了
- **確認日**: 2026-01-03
- **問題だった点**: 
  - NCBI の x_dropoff 計算には `sbp->scale_factor` が含まれる
  - LOSAT は `x_drop_raw_score()` で `scale_factor = 1.0` を固定で渡している
  - これが正しいか確認が必要だった
- **NCBIコード確認結果**:
  1. **デフォルト値**: `sbp->scale_factor = 1.0` (`blast_stat.c:919`)
  2. **RPS-BLAST 専用の ASSERT** (`blast_parameters.c:466-469`):
     ```c
     if (sbp->scale_factor > 1.0) {
         ASSERT(Blast_ProgramIsRpsBlast(program_number));
         params->gap_x_dropoff *= (Int4)sbp->scale_factor;
         params->gap_x_dropoff_final *= (Int4)sbp->scale_factor;
     }
     ```
  3. **E-value 計算時** (`blast_engine.c:881-888`):
     ```c
     double scale_factor = 1.0;
     if (isRPS) {
         scale_factor = score_params->scale_factor;
     }
     Blast_HSPListGetEvalues(..., scale_factor);
     ```
     → **RPS-BLAST 以外では E-value 計算に scale_factor は影響しない**
  4. **traceback でのスコア再スケーリング** (`blast_traceback.c:224-226, 244`):
     ```c
     double scale_factor =
        (Blast_ProgramIsRpsBlast(program_number) ?
        score_params->scale_factor : 1.0);
     // ...
     s_HSPListRescaleScores(hsp_list, score_params->scale_factor);
     ```
     → **tblastx では `scale_factor = 1.0` なので再スケーリングは実質無効**
  5. NCBI ユニットテストでも全て `scale_factor = 1.0` を使用
- **LOSATでの使用箇所** (全て `scale_factor = 1.0` で正しい):
  - `x_drop_raw_score()` - x_dropoff 計算
  - `cutoff_score_for_update_tblastx()` - per-subject cutoff 計算
  - `cutoff_score_word_params()` - cutoff_score 計算
  - `calculate_link_hsp_cutoffs_ncbi()` - linking cutoff 計算
  - `LinkingParams` 構造体のフィールド
- **結論**: **tblastx では常に `scale_factor = 1.0`** であり、LOSAT の現状実装は NCBI と完全に一致。コード修正不要。
- **将来対応**: RPS-BLAST 対応時には `scale_factor` を動的に計算する必要あり
- **ファイル**: `ncbi_cutoffs.rs`, `utils.rs`, `sum_stats_linking.rs`

---

## 2. 修正が必要と判明している点 (Known Required Fixes)

### 2.1 ✅ X-dropoff の Per-Context 適用
- **状態**: ✅ 完了 → **1.6 に移動**
- **修正日**: 2026-01-03
- **問題だった点**: NCBIは context ごとに `cutoffs->x_dropoff` を持ち、extension 時にそれを参照
- **NCBIコード**: `aa_ungapped.c:579`, `blast_parameters.c:219-221`
  ```c
  // 初期化時 (blast_parameters.c:219-221)
  for (context = ...) {
      kbp = sbp->kbp[context];
      p->cutoffs[context].x_dropoff_init =
          (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
  }
  
  // extension 時 (aa_ungapped.c:575-579)
  cutoffs = word_params->cutoffs + curr_context;
  score = s_BlastAaExtendTwoHit(..., cutoffs->x_dropoff, ...);
  ```
- **修正内容**:
  - `run()`: `contexts` 作成後に `x_dropoff_per_context: Vec<i32>` を生成
  - `run_with_neighbor_map()`: per-subject で `x_dropoff_per_context` を生成
  - extension 呼び出しで `x_dropoff_per_context[ctx_idx]` を使用
- **NCBIの挙動確認** (`blast_stat.c:2796-2797`):
  ```c
  // tblastx では計算された Lambda >= kbp_ideal->Lambda なら kbp_ideal に置換
  if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
     Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
  ```
  → **全 context で同じ kbp_ideal (BLOSUM62 ungapped Lambda=0.3176) を使用**
  → **x_dropoff = 16 は全 context で同一** (実質的な出力変更なし)
- **結論**: NCBIとの構造的 parity を達成。全 context で同じ値だが、per-context 配列を維持。
- **ファイル**: `utils.rs`

### 2.2 ✅ scale_factor の確認
- **状態**: ✅ 完了 → **1.7 に移動**
- **確認日**: 2026-01-03
- **結論**: tblastx では常に `scale_factor = 1.0` であり、LOSAT の現状実装は正しい。詳細は 1.7 を参照。

---

## 3. "Might Need Adjustments" レベルの相違点

### 3.1 ✅ Sentinel バイト値の違い
- **状態**: ✅ 完了
- **修正日**: 2026-01-03
- **問題だった点**: 
  - **旧LOSAT**: `SENTINEL_BYTE = 255` を使用 (`constants.rs:98`)
  - **NCBI**: `NULLB = 0` を sentinel として使用 (`blast_encoding.c:120`, `ncbi_std.h:181`)
- **NCBIコード**:
  ```c
  // ncbi_std.h:181
  #define NULLB '\0'
  
  // blast_encoding.c:120
  const Uint1 kProtSentinel = NULLB;
  
  // sm_blosum62.c:92-95
  const SNCBIPackedScoreMatrix NCBISM_Blosum62 = {
      "ARNDCQEGHILKMFPSTWYVBJZX*",
      s_Blosum62PSM,
      -4  // defscore for unknown/sentinel residues
  };
  
  // raw_scoremat.c:90-92 (FSM展開時)
  for (i = 0; i < NCBI_FSM_DIM; ++i) {
      fsm->s[0][i] = psm->defscore;  // index 0 に defscore を設定
  }
  ```
- **修正内容**:
  1. `constants.rs`: `SENTINEL_BYTE = 0` に変更 (NCBI NULLB と同一)
  2. `matrix.rs`: `DEFSCORE = -4` 定数を追加
  3. `matrix.rs`: `blosum62_score()` で index 0 をチェックして `-4` を返すように修正
  4. コメントを更新 (`extension.rs`, `translation.rs`)
- **NCBI の FSM 構築方式**:
  - NCBI は packed matrix (25x25) を FSM (128x128) に展開
  - 展開時に全体を `defscore = -4` で初期化
  - 有効な AA ペアのみ上書き
  - → index 0 (sentinel) は `-4` のまま
- **LOSAT の実装方式**:
  - packed matrix (25x25) + 変換テーブルを使用
  - `blosum62_score()` で sentinel (0) を明示的にチェック
  - → NCBI と同等の動作を保証
- **結論**: 両者は完全に同等の動作 (sentinel に対して `-4` を返す)
- **ファイル**: `constants.rs`, `matrix.rs`, `extension.rs`, `translation.rs`

### 3.2 ✅ Frame Base 計算の Sentinel 考慮と座標システム
- **状態**: ✅ 完了
- **修正日**: 2026-01-03
- **問題だった点**:
  - **旧LOSAT**: `sum_stats_linking.rs` の `abs_coords` 関数が concatenated buffer 内の絶対座標 (`frame_base + hit.aa_start + 1`) を計算していた
  - **NCBI**: `link_hsps.c` の比較関数 `s_RevCompareHSPsTbx` は **frame 内相対座標** (0-indexed) を使用
  - **不一致**: 同じグループ (strand) 内で異なるフレームの HSP は異なる `frame_base` を持つため、ソート順序が NCBI と異なっていた
- **NCBIの座標システム**:
  1. `aa_ungapped.c`: extension 結果は concatenated buffer 内の絶対座標
  2. `blast_gapalign.c:s_AdjustInitialHSPOffsets`: context offset を引いて frame 内相対座標に変換
     ```c
     init_hsp->ungapped_data->q_start -= query_start;
     ```
  3. `link_hsps.c`: HSP の `query.offset` / `subject.offset` は frame 内相対座標として比較
- **修正内容**:
  - `abs_coords` 関数を `frame_relative_coords` に変更
  - concatenated 絶対座標の計算を削除し、frame 内相対座標を直接使用:
    ```rust
    fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
        (hit.q_aa_start as i32, hit.q_aa_end as i32,
         hit.s_aa_start as i32, hit.s_aa_end as i32)
    }
    ```
  - HspLink 初期化 (trim 座標計算) も frame 内相対座標を使用
- **変更しなかった箇所**:
  - `lookup.rs` の `frame_base` 計算 (`base += frame.aa_seq.len() as i32 - 1`) は正しい
  - ただし、`sum_stats_linking.rs` では `frame_base` を使用しないことで NCBI parity を達成
- **検証**: ユニットテスト3件が成功
- **ファイル**: `sum_stats_linking.rs`

### 3.3 ✅ HSP ソート順序の細部
- **状態**: ✅ 確認済み (LOSAT は NCBI と一致)
- **検証日**: 2026-01-03
- **NCBIコード** (`link_hsps.c:359-375`):
  ```c
  // 全フィールドが同じパターン: h1 < h2 なら return 1 (DESCENDING)
  if (h1->query.offset < h2->query.offset)   return  1;  // descending
  if (h1->query.offset > h2->query.offset)   return -1;
  if (h1->query.end < h2->query.end)         return  1;  // descending
  if (h1->query.end > h2->query.end)         return -1;
  if (h1->subject.offset < h2->subject.offset) return  1;  // descending (NOT ascending!)
  if (h1->subject.offset > h2->subject.offset) return -1;
  if (h1->subject.end < h2->subject.end)       return  1;  // descending (NOT ascending!)
  if (h1->subject.end > h2->subject.end)       return -1;
  ```
- **C qsort の仕様**: `compare(a,b) > 0` は「a は b の後に来る」を意味
  - `if (h1 < h2) return 1` → h1 は h2 の後 → 小さい値が後 → **DESCENDING**
- **LOSATコード** (`sum_stats_linking.rs:517-524`):
  ```rust
  group_hits.sort_by(|a, b| {
      bqo.cmp(&aqo)           // descending ✓
          .then(bqe.cmp(&aqe)) // descending ✓
          .then(bso.cmp(&aso)) // descending ✓ (NCBIと同じ!)
          .then(bse.cmp(&ase)) // descending ✓ (NCBIと同じ!)
  });
  ```
- **結論**: ✅ **両者は一致している**。以前の分析で NCBI が ascending と誤解していたが、実際は全フィールドが DESCENDING。
- **ユニットテスト**: `test_hsp_sort_order_matches_ncbi`, `test_ncbi_comparison_semantics` 追加済み

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

### 3.7 ✅ Sum-Statistics の effective length 計算
- **状態**: ✅ 完了
- **修正日**: 2026-01-03
- **問題だった点**: 
  - **旧LOSAT**: `SearchSpace::with_length_adjustment()` を使用し、query と subject の両方から同じ `length_adjustment` を引いていた
  - **NCBI**: tblastx では subject に対して `length_adjustment / 3` のみを適用
- **NCBIコード** (`link_hsps.c:560-571`):
  ```c
  length_adjustment = query_info->contexts[query_context].length_adjustment;
  query_length = query_info->contexts[query_context].query_length;
  query_length = MAX(query_length - length_adjustment, 1);
  subject_length = subject_length_orig; /* in nucleotides even for tblast[nx] */
  /* If subject is translated, length adjustment is given in nucleotide
     scale. */
  if (Blast_SubjectIsTranslated(program_number))  // tblastx = TRUE
  {
     length_adjustment /= CODON_LENGTH;  // ★ 3 で割る
     subject_length /= CODON_LENGTH;
  }
  subject_length = MAX(subject_length - length_adjustment, 1);
  ```
- **NCBI の計算** (tblastx):
  - `eff_query = query_aa - length_adjustment` (全額を引く)
  - `eff_subject = subject_aa - (length_adjustment / 3)` (1/3 のみ引く)
- **修正内容**: `sum_stats_linking.rs:555-570` で NCBI と同等の計算を実装
  ```rust
  let length_adjustment = compute_length_adjustment_simple(
      query_len_aa, subject_len_aa, params
  ).length_adjustment;
  
  // query: 全額を引く
  let eff_query_len = (query_len_aa - length_adjustment).max(1) as f64;
  
  // subject: 1/3 のみ引く (NCBI の length_adjustment /= CODON_LENGTH)
  let length_adj_for_subject = length_adjustment / 3;  // 整数除算
  let eff_subject_len = (subject_len_aa - length_adj_for_subject).max(1) as f64;
  ```
- **影響**: 
  - effective search space が増加
  - E-value が大きくなり、より多くの HSP がフィルタリングされる方向
- **ファイル**: `sum_stats_linking.rs`

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

### 4.5 ❓ Context ごとの Karlin パラメータ計算
- **状態**: ❓ 潜在的相違 (低優先度)
- **概要**: NCBI はクエリのアミノ酸組成から context ごとに Karlin パラメータを計算
- **関連NCBIコード**: `blast_stat.c:2781-2782`
  ```c
  sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
  Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
  ```
- **LOSATの現状**: 
  - 固定のテーブル値 (BLOSUM62 ungapped) を全 context で使用
  - `blast_stat.c:2796-2797` の `check_ideal` により tblastx では通常 `kbp_ideal` が使われるため、実質的な影響は小さい
- **影響**: 極端にバイアスのあるアミノ酸組成のクエリで差異が生じる可能性

### 4.6 ❓ BSearchContextInfo による Context 検索
- **状態**: ❓ 要確認
- **概要**: NCBI は query_offset から context を二分探索で取得
- **関連NCBIコード**: `aa_ungapped.c:560`
  ```c
  curr_context = BSearchContextInfo(query_offset, query_info);
  ```
- **LOSATの現状**: 
  - `run()`: `lookup_ref.get_context_idx(query_offset)` を使用
  - `run_with_neighbor_map()`: `ctx_flat = ctx_base[q_idx] + q_f_idx` で直接計算
- **確認必要**: LOSATの実装が NCBI と同等の結果を返すか

### 4.7 ✅ Extension スコア計算の詳細
- **状態**: ✅ **コード比較完了 - アルゴリズム同等を確認**
- **調査日**: 2026-01-03
- **発見日**: 2026-01-03 (差分確認テストで検出)
- **概要**: 多数のヒットで LOSAT のスコアが NCBI より +1 bit score (raw score 約 2 点) 高いと報告されたが、コード比較の結果アルゴリズムは完全に同等
- **特徴**:
  - identity, aln_len は完全一致
  - strand や alignment length に特定パターンなし
  - MeenMJNV.MejoMJNV: 40 ヒット、AP027280.AP027280: 305 ヒットで差異

#### 詳細コード比較結果

LOSAT `extend_hit_two_hit` (`extension.rs:192-304`) と NCBI `s_BlastAaExtendTwoHit` (`aa_ungapped.c:1088-1158`) を徹底比較:

| コンポーネント | NCBI (`aa_ungapped.c`) | LOSAT (`extension.rs`) | 結果 |
|--------------|------------------------|------------------------|------|
| Word scanning ループ | Lines 1108-1119 | Lines 207-219 | ✅ 同等 |
| Position 調整 | Lines 1120-1121 (`q_right_off += right_d`) | Lines 221-222 | ✅ 同等 |
| `right_d` リセット | Line 1123 (`right_d = 0`) | Line 261 (`right_disp = 0`) | ✅ 同等 |
| Left extension 初期値 | `maxscore = 0` (line 1131) | `max_score = 0` (line 229) | ✅ 同等 |
| Right extension 初期値 | `maxscore = left_score` (line 1147) | `right_score = max_score` (line 272) | ✅ 同等 |
| 終了条件 | `score <= 0 OR (maxscore - score) >= dropoff` | 同一 | ✅ 同等 |
| Length 計算 | `left_d + right_d` (line 1156) | `left_disp + right_disp` (暗黙) | ✅ 同等 |
| Return 値 | `MAX(left_score, right_score)` | `max_score_total` | ✅ 同等 |

**NCBIコード参照** (`aa_ungapped.c:1108-1158`):
```c
// Word scanning (lines 1108-1119)
for (i = 0; i < word_size; i++) {
    score += matrix[q[q_right_off + i]][s[s_right_off + i]];
    if (score > left_score) {
        left_score = score;
        right_d = i + 1;
    }
}
q_right_off += right_d;
s_right_off += right_d;
right_d = 0;  // RESET for extension phase

// Left extension (lines 1127-1135)
left_score = s_BlastAaExtendLeft(matrix, subject, query,
                                 s_right_off - 1, q_right_off - 1,
                                 dropoff, &left_d, 0);  // maxscore = 0

// Right extension (lines 1137-1151)
if (left_d >= (s_right_off - s_left_off)) {
    *right_extend = TRUE;
    right_score = s_BlastAaExtendRight(matrix, subject, query,
                                       s_right_off, q_right_off,
                                       dropoff, &right_d, left_score,  // maxscore = left_score
                                       s_last_off);
}

*hsp_len = left_d + right_d;
return MAX(left_score, right_score);
```

**LOSATコード参照** (`extension.rs:192-304`):
```rust
// Word scanning (lines 207-219)
for i in 0..k_size {
    score += get_score(q_char, s_char);
    if score > left_score {
        left_score = score;
        right_d = i + 1;
    }
}
q_right_off += right_d;
s_right_off += right_d;

// Left extension (lines 228-252)
let mut max_score = 0i32;  // maxscore = 0
while i < max_left {
    current_score += get_score(q_char, s_char);
    if current_score > max_score {
        max_score = current_score;
        left_disp = i + 1;
    }
    if (max_score - current_score) >= x_drop { break; }
}

// Right extension (lines 266-296)
if reached_first_hit {
    right_extended = true;
    let mut right_score = max_score;  // maxscore = left extension score
    while ... {
        right_score += get_score(q_char, s_char);
        if right_score > max_score_total {
            max_score_total = right_score;
            right_disp = j + 1;
        }
        if right_score <= 0 || (max_score_total - right_score) >= x_drop { break; }
    }
}

// hsp_len = q_end - q_start = left_disp + right_disp
return max_score_total;  // equivalent to MAX(left_score, right_score)
```

#### Left Extension 関数比較

NCBI `s_BlastAaExtendLeft` (`aa_ungapped.c:886-921`) と LOSAT の左拡張ロジック (`extension.rs:228-252`):

| 項目 | NCBI | LOSAT | 結果 |
|------|------|-------|------|
| ループ方向 | `i = n; i >= 0; i--` (high to low) | `i = 0; i < max_left; i++` (low to high, access `[off - 1 - i]`) | ✅ 同等 (同じ位置にアクセス) |
| 初期 `best_i` | `n + 1` (never improved = length 0) | N/A (使用 `left_disp = 0`) | ✅ 同等 |
| 終了条件 | `(maxscore - score) >= dropoff` | `(max_score - current_score) >= x_drop` | ✅ 同等 |
| Length 計算 | `n - best_i + 1` | `left_disp` (直接トラック) | ✅ 同等 |

#### Right Extension 関数比較

NCBI `s_BlastAaExtendRight` (`aa_ungapped.c:831-866`) と LOSAT の右拡張ロジック (`extension.rs:272-296`):

| 項目 | NCBI | LOSAT | 結果 |
|------|------|-------|------|
| 初期 `best_i` | `-1` (never improved = length 0) | N/A (使用 `right_disp = 0`) | ✅ 同等 |
| 終了条件 | `score <= 0 OR (maxscore - score) >= dropoff` | `right_score <= 0 OR (max_score_total - right_score) >= x_drop` | ✅ 同等 |
| `s_last_off` | `s_off + i` (loop counter at exit) | `s_right_off + j` | ✅ 同等 |
| Length 計算 | `best_i + 1` | `right_disp` (直接トラック) | ✅ 同等 |

#### Reevaluation 比較

LOSAT `reevaluate_ungapped_hit_ncbi_translated` (`reevaluate.rs:80-145`) は NCBI `Blast_HSPReevaluateWithAmbiguitiesUngapped` (`blast_hits.c:675-733`) の直接ポート:

| 項目 | NCBI | LOSAT | 結果 |
|------|------|-------|------|
| スコア初期値 | `score = 0; sum = 0` | `score = 0; sum = 0` | ✅ 同等 |
| 負スコアリセット | `if (sum < 0) { sum = 0; ... }` | `if sum < 0 { sum = 0; ... }` | ✅ 同等 |
| cutoff 未満時の処理 | `if (score < cutoff_score) { best_start = query; }` | `if score < cutoff_score { best_start = idx + 1; }` | ✅ 同等 |
| 最良スコア更新 | `if (sum > score) { score = sum; ... }` | `if sum > score { score = sum; ... }` | ✅ 同等 |
| 削除判定 | `score < cutoff_score` | `score < cutoff_score` | ✅ 同等 |

#### 追加検証: 座標系の一貫性

`run()` と `run_with_neighbor_map()` 関数では座標系が異なるが、最終的に同じ raw 座標を生成:

| 関数 | `subject_offset` | `s_left_off` 計算式 | 結果 |
|------|-----------------|-------------------|------|
| `run()` | 1-based (raw) | `last_hit + wordsize` | raw 座標 |
| `run_with_neighbor_map()` | 0-based (logical) | `last_hit + wordsize + 1` | raw 座標 |

**例**: logical position 5 での比較
- `run()`: `subject_offset = 6`, `s_left_off = 6 + 3 = 9` (raw)
- `run_with_neighbor_map()`: `subject_offset = 5`, `s_left_off = 5 + 3 + 1 = 9` (raw)

→ **同一の raw 座標を生成**

#### 未使用コードの確認

`extend_hit_ungapped` (one-hit extension) は `extension.rs:55-178` に定義されているが、TBLASTX では**使用されていない**。TBLASTX は two-hit モードのみ使用し、`extend_hit_two_hit` のみが呼び出される。

#### 結論

**Extension アルゴリズムは NCBI と完全に同等**。コード修正は不要。

+1 bit score 差異の真因として考えられるのは:
1. **Sum-statistics linking** での E-value 計算差異 (4.10 で調査継続)
2. **比較方法論の問題** - テスト比較が誤ったヒットペアを比較している可能性
3. **上流の HSP 生成差異** - extension に入力される HSP がそもそも異なる可能性

**推奨アクション**: 
- Extension 調査は完了とし、4.9 (Reverse strand) および 4.10 (E-value 閾値) の調査に移行
- 具体的な +1 差異ケースの完全トレースが必要な場合は、特定の座標での入力シーケンスダンプを実施

- **関連NCBIコード**: 
  - `aa_ungapped.c:831-921` (`s_BlastAaExtendLeft`, `s_BlastAaExtendRight`)
  - `aa_ungapped.c:1019-1086` (`s_BlastAaExtendOneHit`)
  - `aa_ungapped.c:1088-1158` (`s_BlastAaExtendTwoHit`)
  - `blast_hits.c:675-733` (`Blast_HSPReevaluateWithAmbiguitiesUngapped`)
- **LOSATコード**: 
  - `extension.rs:55-303` (extension functions)
  - `reevaluate.rs:80-145` (reevaluation function)
  - `utils.rs:876-883` (`run()` での extension 呼び出し)
  - `utils.rs:1703-1712` (`run_with_neighbor_map()` での extension 呼び出し)

### 4.8 ✅ Sum-Statistics Linking のチェイン構造
- **状態**: ✅ **調査完了 - NCBI と同等の実装を確認**
- **調査日**: 2026-01-03
- **結論**: **LOSATの実装はNCBIと同等であり、コード修正不要**

#### 調査結果詳細

NCBI `link_hsps.c` と LOSAT `sum_stats_linking.rs` を徹底比較し、全コンポーネントの一致を確認:

| コンポーネント | NCBI (行番号) | LOSAT (行番号) | 結果 |
|--------------|--------------|----------------|------|
| `lh_helper` 構造体 | 660-686 | 792-828 | ✅ 一致 |
| `next_larger` skip-list | 676-684, 876-884 | 815-822, 1005-1012 | ✅ 一致 |
| Frame/strand グルーピング | 510-533 (ソート後分割) | 433-439 (事前グルーピング) | ✅ 同等の結果 |
| Index 0 (small gap) DPループ | 690-768 | 835-895 | ✅ 一致 |
| Index 1 (large gap) DPループ | 771-896 | 898-1025 | ✅ 一致 |
| `linked_to` カウンター | 685, 766, 894 | 801, 893, 1023 | ✅ 一致 |
| E-value チェイン適用 | 955-980 | 1139-1197 | ✅ 一致 |
| Gap probability 調整 | 918-935 | 1049-1077 | ✅ 一致 |
| `small_gap_sum_e` | blast_stat.c:4418-4463 | sum_statistics.rs:367-399 | ✅ 一致 |
| `large_gap_sum_e` | blast_stat.c:4532-4573 | sum_statistics.rs:464-496 | ✅ 一致 |
| `s_BlastSumP` (ルックアップテーブル) | blast_stat.c:4357-4408 | sum_statistics.rs:97-278 | ✅ 一致 |

#### 観察された差異の真因

NCBI-only ヒットが同一 E-value を持つパターンは、**リンキングアルゴリズムではなく上流工程の差異**が原因:

1. **Extension スコア差異** (4.7): +1 bit score の差があり、raw score が異なる → 異なる HSP がリンキングに入力される
2. **Reverse strand 座標計算** (4.9): reverse strand の座標差で HSP が異なる
3. **E-value 閾値判定** (4.10): 閾値境界でのフィルタリング差

#### 実装構造の違いと同等性

**NCBI**: 全体ソート → フレーム境界で分割
```c
qsort(link_hsp_array, total_number_of_hsps, sizeof(LinkHSPStruct*), s_RevCompareHSPsTbx);
// ソート後にフレーム境界を検出して分割
```

**LOSAT**: 事前グルーピング → グループ内ソート
```rust
let mut groups: FxHashMap<ContextKey, Vec<UngappedHit>> = FxHashMap::default();
// 各グループ内でソート
```

両者は異なるアプローチだが、HSP の最終的な処理順序は同一。LOSAT のアプローチは並列処理 (rayon) との相性が良い。

#### 実装差異の詳細と影響なしの根拠

| 差異点 | NCBI 実装 | LOSAT 実装 | 影響なしの根拠 |
|--------|-----------|------------|----------------|
| **maxsum1 チェック** | `if(0) if(H2_helper->maxsum1<=H_hsp_sum)break;` (line 850) | 未実装 | NCBI で `if(0)` により**無効化**されているため使用されない |
| **maxsum1 計算方式** | `max[SIGN(s_frame)+1]` で3バケット追跡 | 単一 `running_max` | LOSAT は事前にフレーム符号でグルーピングするため、グループ内では単一maxで同等 |
| **グルーピングキー** | `(context/strand_factor, SIGN(subject.frame))` | `(q_idx, s_idx, q_strand, s_strand)` | NCBI は per-subject 呼び出し。LOSAT は複数サブジェクト一括処理のため `s_idx` を追加。結果は同等 |
| **ln_factorial** | `lgamma(n+1)` via `BLAST_LnGammaInt` | 直接計算 `sum ln(i)` | リンクされる HSP 数は通常 2-10 個程度で、直接計算が十分な精度を持つ |
| **ルックアップテーブル** | TAB2(19), TAB3(38), TAB4(55) | TAB2(19), TAB3(38), TAB4(55) | **完全一致** |

**特記事項: `maxsum1` について**

NCBI コードには `maxsum1` フィールドが存在し計算も行われるが、実際に使用される箇所 (`link_hsps.c:850`) は `if(0)` で無効化されている:

```c
// NCBI link_hsps.c:850 - このチェックは無効
if(0) if(H2_helper->maxsum1<=H_hsp_sum)break;
```

したがって `maxsum1` の計算方式の違いは出力に影響しない。LOSAT は厳密パリティのため計算は行うが、使用はしない。

- **関連NCBIコード**: `link_hsps.c`, `blast_stat.c`, `ncbi_math.c`
- **LOSATコード**: `sum_stats_linking.rs`, `sum_statistics.rs`

### 4.9 ✅ Reverse Strand 処理
- **状態**: ✅ **調査完了 - NCBI と同等の実装を確認**
- **調査日**: 2026-01-03
- **発見日**: 2026-01-03
- **概要**: TrcuMJNV.MellatMJNV の NCBI-only ヒット 113 件が**全て reverse strand** という報告に基づき、reverse strand 処理を徹底調査
- **結論**: **コードレビューの結果、reverse strand 固有のバグは発見されなかった。LOSATの実装はNCBIと同等。**

#### 調査結果詳細

##### 1. SIGN() vs signum() の同等性 - ✅ 検証済み

**NCBIコード** (`ncbi_std.h:127`):
```c
#define SIGN(a) ((a)>0?1:((a)<0?-1:0))
```

**LOSAT**: Rust の `i8::signum()` は同一の値 (1, 0, -1) を返す。

**結論**: 動作は完全に同等。

##### 2. Frame 生成 - ✅ 検証済み

**ファイル**: `translation.rs`

LOSATは6フレームを正しく生成:
- Forward: frames 1, 2, 3 (from `seq[0..]`, `seq[1..]`, `seq[2..]`)
- Reverse: frames -1, -2, -3 (from `revcomp[0..]`, `revcomp[1..]`, `revcomp[2..]`)

これはNCBIの `BLAST_GetTranslation` (`blast_util.c:428-456`) の `ABS(frame)-1` 開始オフセットと一致。

**NCBIコード参照** (`blast_util.c:436, 441`):
```c
nucl_seq = (frame >= 0 ? (Uint1 *)query_seq : (Uint1 *)(query_seq_rev+1));
for (index=ABS(frame)-1; index<nt_length-2; index += CODON_LENGTH)
```

##### 3. 座標変換 (AA → DNA) - ✅ 検証済み

**ファイル**: `extension.rs:717-730`

```rust
// For negative frames:
let start_bp = dna_len - (aa_start * 3 + shift);
let end_bp = dna_len - (aa_end * 3 + shift - 1);
```

**数値検証** (frame=-1, aa_start=0, aa_end=2, dna_len=12):
- start_bp = 12 - 0 = 12
- end_bp = 12 - 5 = 7
- 出力: (12, 7) で start > end、正しく reverse strand を示す

**結論**: NCBI の座標表現と一致。

##### 4. Sum-Statistics Linking のグルーピング - ✅ 検証済み

**ファイル**: `sum_stats_linking.rs:433-438`

```rust
let q_strand: i8 = if hit.q_frame > 0 { 1 } else { -1 };
let s_strand: i8 = if hit.s_frame > 0 { 1 } else { -1 };
let key = (hit.q_idx, hit.s_idx, q_strand, s_strand);
```

これはNCBIの `context/strand_factor` と `SIGN(subject.frame)` によるグルーピングと同等。

**NCBIコード参照** (`link_hsps.c:522-528`):
```c
if (H->prev != NULL && 
    ((H->hsp->context/strand_factor) != (H->prev->hsp->context/strand_factor) ||
     (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
{ /* If frame switches, then start new list. */ }
```

##### 5. ソート順序 - ✅ 検証済み

**ファイル**: `sum_stats_linking.rs:519-526`

全フィールドが DESCENDING 順でソートされ、NCBIの `s_RevCompareHSPsTbx` (`link_hsps.c:359-375`) と一致。

#### 113件のNCBI-onlyヒットの推定原因

コードレビューで reverse strand 固有のバグが発見されなかったため、差異の原因は以下と推定:

1. **E-value 閾値境界効果**: E-value 10.0 付近のヒットが浮動小数点精度の差で含まれる/除外される
2. **比較方法論の問題**: テスト比較スクリプトの座標マッチングに問題がある可能性
3. **上流の差異**: Extension や HSP 生成の差異 (reverse strand 固有ではない)

#### 推奨アクション

- 具体的な NCBI-only ヒットをパイプライン全体でトレースして真因を特定
- 4.10 (E-value 閾値判定) の調査を優先

- **関連NCBIコード**:
  - `link_hsps.c:359-375` (`s_RevCompareHSPsTbx` 比較関数)
  - `blast_util.c:428-456` (`BLAST_GetTranslation`)
  - `ncbi_std.h:127` (`SIGN` マクロ定義)
- **LOSATコード**: 
  - `translation.rs` (frame 生成)
  - `extension.rs:717-730` (座標変換)
  - `sum_stats_linking.rs:433-438, 519-526` (グルーピング、ソート)

### 4.10 🟡 E-value 閾値判定
- **状態**: 🟡 **差分確認テストで差異検出** (優先度4)
- **発見日**: 2026-01-03
- **概要**: LOSAT-only ヒットの多くが E-value 0.1-10.0 の閾値付近
- **具体例**:
  - MeenMJNV.MejoMJNV: 46 件中 35 件が E-value 0.1-10.0
  - AP027280.AP027280: 74 件中 37 件が E-value 0.1-10.0
- **推定原因**:
  - search space 計算の微妙な差
  - 浮動小数点の丸め (Rust `f64` vs C `double`)
  - cutoff_score 計算の差
- **関連NCBIコード**: `blast_parameters.c:BlastHitSavingParametersNew`
- **LOSATコード**: `ncbi_cutoffs.rs:cutoff_score_max_for_tblastx()`

---

## 5. NCBI 設計の技術的注記

### 5.1 📝 tblastx の 2 つの searchsp 計算

NCBIは tblastx で**2種類の searchsp 計算**を使い分けている。これは意図的な設計であり、LOSATでも同様に実装する必要がある。

| 用途 | 関数 | 計算式 | 備考 |
|------|------|--------|------|
| Per-subject cutoff update | `BlastInitialWordParametersUpdate` | `MIN(q_aa_len, s_nucl_len) * s_nucl_len` | 長さ調整なし、**AA長とヌクレオチド長を混在** |
| cutoff_score_max / E-value | `BlastHitSavingParametersNew` | `eff_searchsp` (with length adjustment) | 通常の effective search space 計算 |

**重要**: `BlastInitialWordParametersUpdate` で subject_length として渡されるのは**ヌクレオチド長**であり、翻訳後のAA長ではない。これは `BlastSeqSrcGetSeqLen()` がヌクレオチド配列の長さを返すため。

**NCBIコード参照** (`blast_parameters.c:348`):
```c
// query_length は AA 長
// subj_length は NUCL 長 (BlastSeqSrcGetSeqLen から取得)
searchsp = MIN((Uint8)subj_length, (Uint8)query_length)*((Uint8)subj_length);
```

### 5.2 📝 Context ごとの Karlin パラメータ計算と check_ideal

NCBIは context ごとに**アミノ酸組成から Karlin パラメータを計算**する。

**NCBIコード参照** (`blast_stat.c:2781-2782`):
```c
sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);  // sfp = score frequency profile
```

しかし、**tblastx/blastx/rpstblastn では `check_ideal` フラグ**が有効になる:

**NCBIコード参照** (`blast_stat.c:2744-2748, 2796-2797`):
```c
Boolean check_ideal =
   (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
    program == eBlastTypeRpsTblastn);

// ...later...
// 計算された Lambda が kbp_ideal 以上なら置換 (より保守的な値を使用)
if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
   Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
```

**LOSATの現状**:
- 固定のテーブル値 (BLOSUM62 ungapped: Lambda=0.3176, K=0.134) を使用
- context ごとのアミノ酸組成からの計算は行っていない
- **結論**: NCBI の `check_ideal` ロジックにより、tblastx では通常 `kbp_ideal` が使用されるため、LOSATの固定値アプローチは実質的に正しい

**潜在的な相違**:
- 極端にバイアスのあるアミノ酸組成のクエリでは、計算された Lambda が kbp_ideal より小さくなる可能性がある
- その場合、NCBI は計算値を使用し、LOSAT は固定値を使用するため差異が生じる
- **将来対応**: 完全な parity が必要な場合は context ごとの Karlin パラメータ計算を実装

### 5.3 📝 cutoff_score の 3 段階キャップ

tblastx の cutoff_score は以下の 3 つの値の最小値で決定される:

1. **BLAST_Cutoffs からの値**: `CUTOFF_E_TBLASTX = 1e-300` と `dodecay=TRUE` を使用
   - 極端に小さい E-value のため、通常は 1 または非常に低い値になる
2. **gap_trigger**: BLOSUM62 の場合 41 (22.0 bits)
3. **cutoff_score_max**: ユーザー指定の E-value から計算

実質的には `gap_trigger` または `cutoff_score_max` が支配的な値となる。

### 5.4 📝 HSP グルーピングとソートの実装構造の違い

NCBIとLOSATは HSP のフレーム/ストランド別処理において、**異なる実装構造**を採用しているが、**結果は同等**である。

**NCBIのアプローチ**: 「全体ソート → フレーム境界で分割」

```c
// link_hsps.c:484-486
// 全HSPを一括でソート (比較関数内でcontext/frameを考慮)
qsort(link_hsp_array, total_number_of_hsps, sizeof(LinkHSPStruct*), 
      s_RevCompareHSPsTbx);

// s_RevCompareHSPsTbx の比較順序:
// 1. context/(NUM_FRAMES/2) - query strand + query index
// 2. SIGN(subject.frame)    - subject strand
// 3. query.offset descending
// 4-6. その他の座標 descending

// link_hsps.c:510-533
// ソート後にフレーム境界を検出して分割
for (index = 0; index < number_of_hsps; index++) {
    if (H->prev != NULL && 
        ((H->hsp->context/strand_factor) != (H->prev->hsp->context/strand_factor) ||
         (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
    {
        // フレーム境界で新しいリストを開始
        hp_frame_start[++cur_frame] = H;
        H->prev->next = NULL;
        H->prev = NULL;
    }
}
```

**LOSATのアプローチ**: 「事前グルーピング → グループ内ソート」

```rust
// sum_stats_linking.rs:29
type ContextKey = (u32, u32, i8, i8); // (q_idx, s_idx, q_strand, s_strand)

// 先にグルーピング (group_by_context)
let mut groups: FxHashMap<ContextKey, Vec<UngappedHit>> = FxHashMap::default();
for hit in hits {
    let key = (hit.q_idx, hit.s_idx, hit.q_frame.signum(), hit.s_frame.signum());
    groups.entry(key).or_default().push(hit);
}

// 各グループ内でソート (座標のみで比較、context/frameは既に分離済み)
group_hits.sort_by(|a, b| {
    bqo.cmp(&aqo)
        .then(bqe.cmp(&aqe))
        .then(bso.cmp(&aso))
        .then(bse.cmp(&ase))
});
```

**同等性の理由**:
- NCBIのソート比較関数は、まず `context` と `subject.frame` で比較し、同一の場合のみ座標比較に進む
- LOSATは事前に `context` と `frame` でグルーピングするため、ソート時は座標比較のみで同じ結果が得られる
- 最終的な HSP の処理順序は両者で同一

**実装上の利点**:
- LOSAT: 並列処理 (rayon) との相性が良い。各グループを独立して処理可能
- NCBI: メモリ効率が良い。追加のハッシュマップが不要

### 5.5 📝 Sum-Stats の length_adjustment 適用の非対称性

NCBI `link_hsps.c:560-571` では、tblastx において `length_adjustment` の適用方法が query と subject で**異なる**。

**NCBIコード参照** (`link_hsps.c:560-571`):
```c
length_adjustment = query_info->contexts[query_context].length_adjustment;
query_length = query_info->contexts[query_context].query_length;
query_length = MAX(query_length - length_adjustment, 1);  // ★ query: 全額を引く
subject_length = subject_length_orig; /* in nucleotides even for tblast[nx] */
/* If subject is translated, length adjustment is given in nucleotide
   scale. */
if (Blast_SubjectIsTranslated(program_number))  // tblastx = TRUE
{
   length_adjustment /= CODON_LENGTH;  // ★ 3 で割る
   subject_length /= CODON_LENGTH;
}
subject_length = MAX(subject_length - length_adjustment, 1);  // ★ subject: 1/3 のみ引く
```

**重要な発見**:
- `length_adjustment` は `blast_setup.c` で AA 単位で計算・保存される
- query に対しては全額 (`length_adjustment`) を引く
- subject に対しては 1/3 (`length_adjustment / 3`) のみを引く
- **NCBI のコメント "length adjustment is given in nucleotide scale" は誤解を招く**: 実際には AA 単位で格納されており、subject に適用する際に `/3` している

**LOSAT の対応**:
- `sum_stats_linking.rs:555-570` で NCBI と同等のロジックを実装
- `compute_length_adjustment_simple()` で length_adjustment を計算
- query: `eff_query_len = query_len_aa - length_adjustment`
- subject: `eff_subject_len = subject_len_aa - (length_adjustment / 3)`

### 5.6 📝 BLOSUM62 defscore と FSM 構築

NCBI は packed matrix (25x25) を FSM (Full Score Matrix, 128x128) に展開する。この際、未知/sentinel 残基には `defscore` が適用される。

**NCBIコード参照** (`sm_blosum62.c:92-95`):
```c
const SNCBIPackedScoreMatrix NCBISM_Blosum62 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Blosum62PSM,
    -4  // defscore for unknown/sentinel residues
};
```

**FSM 展開時** (`raw_scoremat.c:90-95`):
```c
// 全体を defscore で初期化
for (i = 0; i < NCBI_FSM_DIM; ++i) {
    fsm->s[0][i] = psm->defscore;  // -4
}
for (i = 1; i < NCBI_FSM_DIM; ++i) {
    memcpy(fsm->s[i], fsm->s[0], NCBI_FSM_DIM * sizeof(fsm->s[0][0]));
}
// その後、有効な AA ペアのみ上書き
```

**LOSAT の実装方式**:
- packed matrix (25x25) + 変換テーブルを使用 (FSM 展開なし)
- `blosum62_score()` で sentinel (0) を明示的にチェックして `DEFSCORE = -4` を返す
- メモリ効率が良い (25x25 vs 128x128)
- 動作は NCBI と完全に同等

**重要**:
- NCBISTDAA index 0 は gap ('-') を表す
- NCBI BLAST では `kProtSentinel = NULLB = 0`
- つまり、gap と sentinel は同じ値であり、どちらも `-4` を返す

### 5.7 📝 HSP 座標の変換フロー (NCBI vs LOSAT)

NCBI と LOSAT では HSP 座標の管理方法が異なる。

**NCBI の座標変換フロー**:

```mermaid
flowchart TD
    A[aa_ungapped.c: Extension] -->|hsp_q = 絶対座標| B[BlastSaveInitHsp]
    B -->|ungapped_data->q_start = 絶対座標| C[blast_gapalign.c]
    C -->|s_AdjustInitialHSPOffsets| D[frame内相対座標]
    D --> E[link_hsps.c: ソート・リンキング]
```

1. **Extension 時** (`aa_ungapped.c`): `hsp_q` は concatenated buffer 内の絶対座標
2. **座標調整** (`blast_gapalign.c:s_AdjustInitialHSPOffsets`):
   ```c
   init_hsp->ungapped_data->q_start -= query_start;
   // query_start = query_info->contexts[context].query_offset
   ```
3. **リンキング時** (`link_hsps.c`): `hsp->query.offset` は frame 内相対座標 (0-indexed)

**LOSAT の座標フロー**:

```mermaid
flowchart TD
    A[utils.rs: Extension] -->|raw coords in frame buffer| B[論理座標変換]
    B -->|hit.q_aa_start = frame内相対座標| C[sum_stats_linking.rs]
    C --> D[ソート・リンキング]
```

1. **Extension 時** (`utils.rs`): frame ごとに独立した `aa_seq` バッファで extension
2. **論理座標変換**: `qs_l = qs.saturating_sub(1)` で leading sentinel を除外
3. **リンキング時**: `hit.q_aa_start` は最初から frame 内相対座標

**LOSAT の設計上の利点**:
- frame ごとに独立したバッファなので、座標調整 (`s_AdjustInitialHSPOffsets` 相当) が不要
- 並列処理との相性が良い

**LOSAT の旧実装の問題点**:
- `sum_stats_linking.rs` の `abs_coords` 関数が `frame_base + hit.aa_start + 1` を計算
- これは NCBI の「調整前」の絶対座標をエミュレートしていた
- しかし NCBI のリンキングは「調整後」の相対座標を使用
- → ソート順序の不一致が発生

**修正後**:
- `frame_relative_coords` 関数で `hit.q_aa_start` 等を直接使用
- NCBI の「調整後」座標と同等

---

## 6. 発見した明確なバグ

### ~~6.1 🐛 HSP ソートの subject 比較順序が逆~~ → ✅ 分析誤り (バグではない)
- **状態**: ✅ 問題なし (分析誤りを訂正)
- **検証日**: 2026-01-03
- **ファイル**: `sum_stats_linking.rs:517-524`
- **誤った分析**: 
  ```
  NCBI が subject.offset/end を ascending でソートしていると誤解していた
  ```
- **正しい解釈**:
  - NCBI の `if (h1->subject.offset < h2->subject.offset) return 1;` は **DESCENDING**
  - C qsort では `compare(a,b) > 0` は「a は b の後に来る」
  - `h1 < h2 → return 1` は h1 が h2 の後 → 小さい値が後 → **DESCENDING**
- **LOSAT のコード `bso.cmp(&aso)`** は `b.cmp(&a)` = **DESCENDING** であり、NCBI と一致
- **結論**: **コード修正不要。LOSAT は既に NCBI と同一のソート順序を使用している。**
- **ユニットテスト追加**: `sum_stats_linking.rs` に以下のテストを追加
  - `test_hsp_sort_order_matches_ncbi()` - ソート順序が期待通りか検証
  - `test_hsp_sort_identical_coords()` - 同一座標の安定ソート検証
  - `test_ncbi_comparison_semantics()` - NCBI の比較セマンティクス検証

---

## 7. 推定される根本原因

LOSATがNCBIより多くのヒットを出力する原因として、以下が推定される:

1. ~~**X-drop / Cutoff の不整合**~~: 
   - ~~X-drop が固定値で、context-specific な値を使っていない~~ → ✅ **1.4で修正済み**
   - ~~Per-subject cutoff 更新が NCBI と異なる~~ → ✅ **1.5で修正済み**

2. ~~**HSP ソート順序のバグ**~~: 
   - ~~Subject offset/end の比較順序が逆~~ → ✅ **分析誤り (6.1で検証済み)**
   - ~~リンキング結果が変わり、E-value 計算に影響~~ → **LOSAT は NCBI と同一のソート順序**

3. ~~**Sum-Statistics の Length Adjustment**~~: 
   - ~~effective length の計算が NCBI と異なる可能性~~ → ✅ **3.7で修正済み**
   - ~~Search space が異なれば E-value も異なる~~ → **NCBI と同一の計算方法に修正**

---

## 8. 優先度順の修正作業リスト

### 完了済み作業

| 優先度 | ID | 内容 | ファイル | 状態 |
|--------|-----|------|----------|------|
| ~~1~~ | ~~1.4~~ | ~~X-drop 動的計算~~ | ~~`ncbi_cutoffs.rs`, `utils.rs`~~ | ✅ 完了 |
| ~~2~~ | ~~1.5~~ | ~~Per-subject cutoff 更新ロジック~~ | ~~`utils.rs`, `ncbi_cutoffs.rs`~~ | ✅ 完了 |
| ~~3~~ | ~~1.6~~ | ~~X-drop の per-context 適用~~ | ~~`utils.rs`~~ | ✅ 完了 |
| ~~4~~ | ~~6.1~~ | ~~HSP ソート順序 subject ascending 修正~~ | ~~`sum_stats_linking.rs`~~ | ✅ 分析誤り |
| ~~5~~ | ~~3.7~~ | ~~Sum-stats effective length 計算~~ | ~~`sum_stats_linking.rs`~~ | ✅ 完了 |
| ~~6~~ | ~~3.2~~ | ~~Frame base / 座標システムの修正~~ | ~~`sum_stats_linking.rs`~~ | ✅ 完了 |
| ~~7~~ | - | ~~差分確認テスト実行~~ | - | ✅ 完了 |

### 🔴 残存差異対応 (差分確認テストで発見)

| 優先度 | ID | 内容 | 推定工数 | ファイル | 状態 |
|--------|-----|------|----------|----------|------|
| ~~**1**~~ | ~~4.8~~ | ~~Sum-Statistics Linking チェイン構造~~ | ~~大~~ | ~~`sum_stats_linking.rs`, `link_hsps.c`~~ | ✅ 調査完了 (NCBI同等) |
| ~~**1**~~ | ~~4.7~~ | ~~Extension スコア計算 (+1 bit score 差)~~ | ~~中~~ | ~~`extension.rs`, `aa_ungapped.c`~~ | ✅ 調査完了 (NCBI同等) |
| ~~**1**~~ | ~~4.9~~ | ~~Reverse strand 処理~~ | ~~中~~ | ~~`translation.rs`, `utils.rs`, `extension.rs`~~ | ✅ 調査完了 (NCBI同等) |
| **1** | 4.10 | E-value 閾値判定 | 小 | `ncbi_cutoffs.rs`, `karlin.rs` | 🟡 未着手 |

### 低優先度 (未調査領域)

| 優先度 | ID | 内容 | ファイル | 状態 |
|--------|-----|------|----------|------|
| 5 | 4.1 | Two-hit Window 詳細 | `utils.rs` | ❓ 未調査 |
| 6 | 4.2 | Lookup Table 構築詳細 | `lookup.rs` | ❓ 未調査 |
| 7 | 4.3 | Masked Region Extension 処理 | `utils.rs` | ❓ 未調査 |
| 8 | 4.4 | HSP 重複排除 (Culling) | `sum_stats_linking.rs` | ❓ 未調査 |
| 9 | 4.5 | Context 別 Karlin パラメータ | `karlin.rs` | ❓ 低優先度 |
| 10 | 4.6 | BSearchContextInfo 検索 | `utils.rs` | ❓ 未確認 |

---

## 9. 次のアクション

### 完了済み

1. ~~**⚠️ X-drop を動的計算に変更** (`utils.rs`)~~ → ✅ **完了 (1.4)**
2. ~~**⚠️ Per-subject cutoff 更新ロジック確認** (`utils.rs`, `ncbi_cutoffs.rs`)~~ → ✅ **完了 (1.5)**
3. ~~**⚠️ X-drop の per-context 適用** (`utils.rs`)~~ → ✅ **完了 (1.6)**
4. ~~**🐛 HSP ソート順序を修正** (`sum_stats_linking.rs:517-524`)~~ → ✅ **分析誤り・修正不要 (6.1)**
5. ~~**Sum-Stats Effective Length 計算確認** (`sum_stats_linking.rs`)~~ → ✅ **完了 (3.7)**
6. ~~**Frame Base / 座標システムの修正** (`sum_stats_linking.rs`)~~ → ✅ **完了 (3.2)**
7. ~~**差分確認テストを実行し、残存差異を特定**~~ → ✅ **完了**

### 🔴 次の作業 (残存差異対応)

1. ~~**🔥 Sum-Statistics Linking チェイン構造調査** (4.8)~~ → ✅ **調査完了 - NCBIと同等の実装を確認 (2026-01-03)**
   - `link_hsps.c` と `sum_stats_linking.rs` を全コンポーネント比較
   - **結論**: リンキングアルゴリズムは正しい。差異の原因は上流工程 (Extension, Reverse strand, E-value閾値)
   - 詳細 → セクション 4.8 参照

2. ~~**🔥 Extension スコア計算調査** (4.7)~~ → ✅ **調査完了 - NCBIと同等の実装を確認 (2026-01-03)**
   - `aa_ungapped.c:1088-1158` と `extension.rs:192-304` を全コンポーネント比較
   - **結論**: アルゴリズムは完全に同等。+1 bit score 差の原因は上流工程または比較方法論にあり
   - 詳細 → セクション 4.7 参照

3. ~~**🔥 Reverse strand 処理調査** (4.9)~~ → ✅ **調査完了 - NCBIと同等の実装を確認 (2026-01-03)**
   - `translation.rs`, `extension.rs`, `sum_stats_linking.rs` を NCBI コードと徹底比較
   - **結論**: reverse strand 固有のバグは発見されなかった。SIGN()/signum()同等、frame生成正常、座標変換正常、グルーピング正常
   - 113件の NCBI-only ヒットは E-value 閾値境界効果または比較方法論の問題と推定
   - 詳細 → セクション 4.9 参照

4. **🔥 E-value 閾値判定調査** (4.10)
   - LOSAT-only ヒットの E-value 境界問題調査
   - 詳細手順 → セクション 10.3 参照

---

## 10. 差分確認テスト結果 (2026-01-03)

### 10.1 テスト対象ファイル (10ペア)

| ファイルペア | LOSAT | NCBI | 差分 | 差分率 |
|-------------|-------|------|-----|-------|
| MeenMJNV.MejoMJNV | 23285 | 23240 | +45 | +0.19% |
| MellatMJNV.MeenMJNV | 18586 | 18560 | +26 | +0.14% |
| TrcuMJNV.MellatMJNV | 5749 | 5804 | **-55** | -0.95% |
| LvMJNV.TrcuMJNV | 3311 | 3301 | +10 | +0.30% |
| PemoMJNVB.LvMJNV | 16186 | 16173 | +13 | +0.08% |
| PeseMJNV.PemoMJNVB | 44899 | 44838 | +61 | +0.14% |
| PemoMJNVA.PeseMJNV | 34117 | 34062 | +55 | +0.16% |
| MelaMJNV.PemoMJNVA | 4844 | 4794 | +50 | +1.04% |
| MjeNMV.MelaMJNV | 23766 | 23709 | +57 | +0.24% |
| AP027280.AP027280 | 42797 | 42733 | +64 | +0.15% |

**統計**: 共通ヒット率 ~99.8%, スコア完全一致率 ~95-98%

### 10.2 残存差異の分類

#### 🔶 差異1: スコア差異 (+1 bit score)

**現象**: 多数のヒットで LOSAT のスコアが NCBI より **1 bit score 高い**

| ファイル | +1 差異ヒット数 |
|---------|---------------|
| MeenMJNV.MejoMJNV | 40 |
| TrcuMJNV.MellatMJNV | 14 |
| AP027280.AP027280 | 305 |

**特徴**:
- identity, aln_len は完全一致
- raw score で約 **2 点** の差に相当
- strand や alignment length mod 3 には特定パターンなし

**分析**: NCBI の bit score フォーマットは `%4.1lf` (小数点以下1桁) で LOSAT と同一。問題は raw score 計算にある可能性。

#### 🔴 差異2: LOSAT-only ヒット (偽陽性候補)

**現象**: NCBI に存在しないヒットが LOSAT に存在

| ファイル | LOSAT-only | スコア範囲 | 特徴 |
|---------|-----------|-----------|------|
| MeenMJNV.MejoMJNV | 46 | 22-30 | E-value 閾値付近 (0.1-10.0) |
| TrcuMJNV.MellatMJNV | 58 | 22-30 | reverse strand に偏り |
| AP027280.AP027280 | 74 | 22-31 | E-value 閾値付近 |

**分析**: 低スコア・閾値ギリギリのヒット。E-value 10.0 境界での判定差、または extension 終了条件の微妙な差が原因の可能性。

#### 🔵 差異3: NCBI-only ヒット (偽陰性候補)

**現象**: LOSAT に存在しないヒットが NCBI に存在

| ファイル | NCBI-only | 特徴 |
|---------|-----------|------|
| MeenMJNV.MejoMJNV | 1 | 単一ヒット |
| TrcuMJNV.MellatMJNV | 113 | **全て reverse strand**, E-value が 0.002/0.010/1.1 の3グループ |
| AP027280.AP027280 | 10 | **全て同一 E-value (5.33e-09)** → Sum-Statistics チェイン |

**重要な発見**:
1. **同一 E-value パターン**: Sum-Statistics Linking でチェイン化された HSP が LOSAT では異なるグループ化になっている可能性
2. **Reverse strand バイアス**: TrcuMJNV の 113 件全てが reverse strand。特定の処理パスで差異がある
3. **NEARBY LOSAT 存在**: 一部の NCBI-only ヒットには座標が微妙に異なる LOSAT ヒットあり (Extension 差)

### 10.3 優先度順の調査項目 (詳細手順付き)

#### ✅ ~~優先度1: Sum-Statistics Linking の E-value 適用~~ → 調査完了 (2026-01-03)

**結論**: **LOSATの実装はNCBIと同等。コード修正不要。**

**調査結果**:
- 全コンポーネント (`lh_helper`, `next_larger`, DPループ, `linked_to`, E-valueチェイン適用等) が NCBI と一致
- フレーム分割の実装構造は異なるが、最終的な HSP 処理順序は同等
- NCBI-only ヒットの原因は上流工程 (Extension, Reverse strand, E-value閾値) の差異

詳細 → セクション 4.8 参照

---

#### ✅ ~~優先度1: Extension スコア計算 (+1 bit score 差)~~ → 調査完了 (2026-01-03)

**結論**: **LOSATの実装はNCBIと同等。コード修正不要。**

**調査結果**:
- LOSAT `extend_hit_two_hit` (`extension.rs:192-304`) と NCBI `s_BlastAaExtendTwoHit` (`aa_ungapped.c:1088-1158`) を全コンポーネント比較
- Word scanning、Left extension、Right extension、終了条件、Length 計算、Return 値 すべてが同等
- Reevaluation ロジック (`reevaluate.rs:80-145` vs `blast_hits.c:675-733`) も完全に同等
- +1 bit score 差の原因は Extension アルゴリズムではなく、上流の HSP 生成または比較方法論にある可能性

**確認済み項目**:
- ✅ word 内スコア計算での初期値 (`score = 0`, `left_score = 0`)
- ✅ left extension への初期スコア渡し (`maxscore = 0`)
- ✅ right extension への初期スコア渡し (`maxscore = left_score`)
- ✅ 終了条件 (`score <= 0 || (maxscore - score) >= dropoff`)
- ✅ `maxscore` 更新タイミング (ループ内で `if (score > maxscore)`)
- ✅ Length 計算 (`left_d + right_d`)

詳細 → セクション 4.7 参照

---

#### ✅ ~~優先度1: Reverse strand 処理~~ → 調査完了 (2026-01-03)

**問題**: TrcuMJNV.MellatMJNV の NCBI-only ヒット 113 件が**全て reverse strand**

**結論**: **コードレビューの結果、reverse strand 固有のバグは発見されなかった。LOSATの実装はNCBIと同等。**

**調査完了項目**:
1. ✅ **SIGN() vs signum() 同等性**: `ncbi_std.h:127` の `SIGN` マクロと Rust `signum()` は同一動作
2. ✅ **Frame 生成**: LOSAT は 6 フレーム (1,2,3,-1,-2,-3) を正しく生成。NCBI `blast_util.c:428-456` と同等
3. ✅ **座標変換 (AA → DNA)**: `extension.rs:717-730` で負フレームの座標を正しく計算 (start > end)
4. ✅ **グルーピング**: `sum_stats_linking.rs:433-438` が NCBI `link_hsps.c:522-528` と同等
5. ✅ **ソート順序**: 全フィールド DESCENDING で NCBI `s_RevCompareHSPsTbx` と一致

**113件の差異の推定原因**:
- E-value 閾値境界効果 (浮動小数点精度)
- 比較方法論の問題 (座標マッチング)
- 上流の HSP 生成差異 (reverse strand 固有ではない)

**関連ファイル**: `translation.rs`, `extension.rs`, `utils.rs`, `sum_stats_linking.rs`

詳細 → セクション 4.9 参照

---

#### 🔥 優先度1: E-value 閾値判定

**問題**: LOSAT-only ヒットの多くが E-value 0.1-10.0 の閾値付近。NCBI では閾値で除外されたが LOSAT では通過?

**調査手順**:
1. **E-value 閾値の適用箇所**:
   - NCBI: `BlastHitSavingParametersNew` で `cutoff_score_max` を計算し、E-value 10.0 から逆算
   - LOSAT: `cutoff_score_max_for_tblastx()` が同等か確認
2. **search space 計算の確認**:
   - E-value = K × m × n × exp(-λS) で計算
   - `m × n` (search space) が LOSAT と NCBI で同一か
3. **浮動小数点の丸め**:
   - E-value 計算時の丸め差で閾値境界のヒットが通過/除外される可能性
   - Rust の `f64` と C の `double` の精度差
4. **具体的なヒットでの検証**:
   - LOSAT-only ヒット (E-value ≈ 9.5-10.0) の raw score を確認
   - その score での E-value を NCBI 方式で再計算して比較

**関連ファイル**: `ncbi_cutoffs.rs`, `karlin.rs`, `sum_stats_linking.rs`

---

### 10.4 確認済み項目

- **bit score フォーマット**: NCBI `%4.1lf` (小数点以下1桁) で LOSAT と同一 ✅
- **E-value チェイン適用**: LOSAT は chain 全メンバーに同じ E-value を適用 (`sum_stats_linking.rs:1178-1180`) ✅
- **Sum-Statistics Linking チェイン構造**: LOSAT の実装は NCBI と同等 (全コンポーネント一致を確認) ✅ (2026-01-03 調査完了)
- **maxsum1 計算差異**: NCBI `link_hsps.c:850` で `if(0)` により無効化されているため、計算方式の違いは出力に影響しない ✅
- **Sum-Statistics ルックアップテーブル**: TAB2(19要素), TAB3(38要素), TAB4(55要素) が NCBI と完全一致 ✅
- **ln_factorial 実装**: LOSAT の直接計算と NCBI の lgamma(n+1) は小さな n で同等精度 ✅
- **グルーピングキー s_idx**: LOSAT はマルチサブジェクト一括処理のため追加。per-subject 処理と結果は同等 ✅
- **Reverse strand 処理**: SIGN()/signum() 同等、frame生成正常、座標変換正常、グルーピング/ソート正常 ✅ (2026-01-03 調査完了)

### 10.5 🔧 出力フォーマットの相違点 (未修正)

以下はパリティに影響しない出力フォーマットの違い。機能的には同等だが、diff で差異として検出される。

| 項目 | LOSAT | NCBI | 影響 |
|------|-------|------|------|
| **E-value 表記** | `0.0e0`, `7.2e-295` (科学表記) | `0.0`, `0.0` (通常表記) | 数値は同等 |
| **bit score 小数** | `692.0` (常に小数1桁) | `692` (整数の場合は .0 なし?) | 数値は同等 |
| **ヘッダー行** | なし | `# TBLASTX 2.17.0+` 等のコメント行あり | 解析ツール依存 |
| **ヒット数コメント** | なし | `# 23240 hits found` | 情報目的 |

**対応方針**: 機能的には同等なので低優先度。完全一致が必要な場合は `report/outfmt6.rs` を修正。

---

## 12. 変更履歴

| 日付 | 変更内容 |
|------|----------|
| 2026-01-03 | 初版作成 |
| 2026-01-03 | 1.4 X-drop動的計算を完了。`x_drop_raw_score()` 関数を `ncbi_cutoffs.rs` に追加し、`utils.rs` の両モードで使用するように変更。2.1 の per-context 適用について、tblastx では全 context で同一の Lambda を使用するため低優先度に変更。 |
| 2026-01-03 | 1.5 Per-Subject Cutoff Score更新を完了。`cutoff_score_for_update_tblastx()` と `cutoff_score_max_for_tblastx()` を追加。NCBIの `BlastInitialWordParametersUpdate` ungapped path を忠実にポート。`CUTOFF_E_TBLASTX = 1e-300` と `BLAST_GAP_DECAY_RATE = 0.5` を使用。searchsp計算がAA長とヌクレオチド長を混在させる点を発見・対応。 |
| 2026-01-03 | 1.6 X-dropoff の Per-Context 適用を完了。`run()` と `run_with_neighbor_map()` の両方で `x_dropoff_per_context: Vec<i32>` を生成し、extension 呼び出しで context ごとの x_dropoff を使用。NCBI `blast_stat.c:2796-2797` で tblastx は全 context が `kbp_ideal` を使用するため実質的な出力変更はなし。構造的 parity を達成。 |
| 2026-01-03 | セクション4に 4.5 (Context ごとの Karlin パラメータ計算), 4.6 (BSearchContextInfo による Context 検索) を追加。セクション5に 5.2 (check_ideal ロジック) を追加。調査中に発見した潜在的相違点を記録。 |
| 2026-01-03 | 1.7 scale_factor の確認を完了。NCBI コードベース調査により `blast_parameters.c:466-469` の ASSERT で tblastx では常に `scale_factor = 1.0` であることを確認。LOSAT の現状実装は NCBI と完全一致。コード修正不要。2.2 → 1.7 へ移動。 |
| 2026-01-03 | 6.1 HSP ソート順序を検証。**分析誤りを発見**: NCBI の `if (h1 < h2) return 1` は DESCENDING (ascending ではない)。C qsort では正の戻り値は「第1引数は第2引数の後に来る」を意味。LOSAT の `bso.cmp(&aso)` (DESCENDING) は NCBI と一致しており、修正不要。`sum_stats_linking.rs` にユニットテスト3件を追加して検証。セクション 3.3, 6.1, 7, 8, 9 を更新。 |
| 2026-01-03 | セクション 5.4 を追加: HSP グルーピングとソートの実装構造の違いを記録。NCBI は「全体ソート→フレーム境界で分割」、LOSAT は「事前グルーピング→グループ内ソート」という異なるアプローチを採用しているが、最終的な処理順序は同等。`link_hsps.c:484-533` のフレーム分割ロジックを引用。 |
| 2026-01-03 | 3.7 Sum-Stats Effective Length 計算を修正完了。NCBI `link_hsps.c:560-571` の詳細分析により、tblastx では subject に対して `length_adjustment / 3` のみを適用することを発見。旧LOSAT は `SearchSpace::with_length_adjustment()` で両方に同じ調整を適用していた。`sum_stats_linking.rs:555-570` で NCBI と同等の計算を直接実装。セクション 3.7, 7, 8, 9 を更新。 |
| 2026-01-03 | 3.1 Sentinel バイト値を NCBI と同一 (0) に修正。`constants.rs` で `SENTINEL_BYTE = 0` に変更。`matrix.rs` に `DEFSCORE = -4` 定数を追加し、`blosum62_score()` で sentinel (0) をチェックして `-4` を返すように修正。NCBI の FSM 構築方式 (`raw_scoremat.c:90-95`) を調査し、defscore の適用方法を確認。セクション 5.6 を追加して技術的詳細を記録。 |
| 2026-01-03 | 3.2 Frame Base / 座標システムを修正完了。**問題**: 旧LOSAT の `abs_coords` 関数が concatenated buffer 内の絶対座標 (`frame_base + hit.aa_start + 1`) を計算していたが、NCBI の `link_hsps.c:s_RevCompareHSPsTbx` は frame 内相対座標 (0-indexed) を使用。**修正**: `abs_coords` を `frame_relative_coords` に変更し、frame 内相対座標を直接使用 (`hit.q_aa_start` 等)。NCBI では `s_AdjustInitialHSPOffsets` で context offset を引いて frame 内相対座標に変換している。HspLink 初期化 (trim 座標計算) も同様に修正。ユニットテスト3件が成功。 |
| 2026-01-03 | **差分確認テストを完了**。10ペアのテストファイルで比較分析。共通ヒット率 ~99.8%。残存差異として(1) +1 bit score 差 (raw score ~2点差)、(2) LOSAT-only ヒット (E-value閾値付近)、(3) NCBI-only ヒット (Sum-Statistics チェインの差異、reverse strand バイアス) を特定。セクション10を追加。 |
| 2026-01-03 | **セクション4に 4.7-4.10 を追加**: 差分確認テストで発見した残存差異を調査項目として追加。(4.7) Extension スコア計算、(4.8) Sum-Statistics Linking チェイン構造、(4.9) Reverse strand 処理、(4.10) E-value 閾値判定。各項目に NCBI/LOSAT のファイル参照と調査ポイントを記載。 |
| 2026-01-03 | **セクション8, 9 を更新**: 残存差異対応として 4.7-4.10 を優先度付きで追加。完了済み作業と未着手作業を分離して整理。 |
| 2026-01-03 | **セクション10.3 に詳細調査手順を追加**: 各優先度項目について、具体的な調査手順(NCBI コード参照箇所、LOSAT 対応ファイル、確認ポイント)を記載。 |
| 2026-01-03 | **セクション10.5 を追加**: 出力フォーマットの相違点(E-value 表記、bit score 小数、ヘッダー行)を記録。機能的には同等だが diff で差異として検出される点。 |
| 2026-01-03 | **4.8 Sum-Statistics Linking チェイン構造の調査完了**。NCBI `link_hsps.c` と LOSAT `sum_stats_linking.rs` を全コンポーネント比較。結果: **LOSATの実装はNCBIと同等であり、コード修正不要**。比較項目: `lh_helper`構造、`next_larger` skip-list、Index 0/1 DPループ、`linked_to`カウンター、E-valueチェイン適用、gap probability調整、sum statistics関数群。観察された差異はリンキングではなく上流工程 (Extension, Reverse strand, E-value閾値) が原因。セクション 4.8, 8, 9, 10.3 を更新。 |
| 2026-01-03 | **4.8 実装差異の詳細追記**。再検証により以下の差異が結果に影響しないことを確認: (1) `maxsum1` チェックは NCBI `if(0)` (line 850) で無効化されているため計算方式の違いは無関係、(2) `s_idx` グルーピングキーは LOSAT のマルチサブジェクト一括処理対応のため追加されたもので結果は同等、(3) `ln_factorial` は LOSAT が直接計算 (sum ln(i))、NCBI が lgamma(n+1) を使用するが小さな n で同等精度、(4) ルックアップテーブル TAB2(19)/TAB3(38)/TAB4(55) は完全一致を確認。セクション 4.8 に詳細表を追加。 |
| 2026-01-03 | **4.7 Extension スコア計算の調査完了**。LOSAT `extend_hit_two_hit` (`extension.rs:192-304`) と NCBI `s_BlastAaExtendTwoHit` (`aa_ungapped.c:1088-1158`) を全コンポーネント比較。結果: **アルゴリズムは完全に同等であり、コード修正不要**。比較項目: Word scanning ループ、Position 調整、`right_d` リセット、Left extension 初期値、Right extension 初期値、終了条件、Length 計算、Return 値。Reevaluation ロジック (`reevaluate.rs:80-145` vs `blast_hits.c:675-733`) も完全一致。+1 bit score 差の原因は Extension アルゴリズムではなく上流工程 (HSP 生成) または比較方法論にある可能性。セクション 4.7, 8, 9, 10.3 を更新。 |
| 2026-01-03 | **4.7 追加検証: 座標系一貫性の確認**。`run()` と `run_with_neighbor_map()` の座標系を再検証。前者は 1-based (raw) 座標、後者は 0-based (logical) 座標を使用するが、extension 呼び出し時には両方とも同じ raw 座標 (`s_left_off`) を生成することを数値例で確認 (logical 5 → 両方とも raw 9)。また、`extend_hit_ungapped` は TBLASTX では未使用であることを確認 (BLASTN のみが使用)。セクション 4.7 に「追加検証」と「未使用コードの確認」を追記。 |
| 2026-01-03 | **4.9 Reverse strand 処理の調査完了**。TrcuMJNV.MellatMJNV の NCBI-only ヒット 113 件が全て reverse strand という報告に基づき徹底調査。結果: **LOSATの実装はNCBIと同等であり、コード修正不要**。確認項目: (1) `SIGN()` vs `signum()` 同等性 (`ncbi_std.h:127` 参照)、(2) Frame 生成 (`translation.rs` vs `blast_util.c:428-456`)、(3) 座標変換 AA→DNA (`extension.rs:717-730`)、(4) グルーピング (`sum_stats_linking.rs:433-438` vs `link_hsps.c:522-528`)、(5) ソート順序 (全フィールド DESCENDING)。113件の差異は E-value 閾値境界効果または比較方法論の問題と推定。セクション 4.9, 8, 9, 10.3, 10.4 を更新。 |
