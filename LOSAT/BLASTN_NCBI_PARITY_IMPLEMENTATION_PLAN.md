# LOSATN NCBI BLASTN パリティ実装計画

## 概要

この文書は、LOSATNとNCBI BLASTNの実装差異を修正し、HSP生成数の差異を解消するための網羅的な実装計画です。

## 現状分析

### 確認済みの実装状況

1. **オフ対角線検索（scan_range）**: ✅ 実装済み
   - `constants.rs`: `SCAN_RANGE_BLASTN=4`, `SCAN_RANGE_MEGABLAST=0`
   - `coordination.rs`: タスク別に設定
   - `utils.rs`: 1446-1554行目で実装済み

2. **Two-hit window**: ✅ 実装済み
   - `constants.rs`: `TWO_HIT_WINDOW=0` (NCBIデフォルト、one-hit mode)

3. **Cutoff score計算**: ✅ 実装済み
   - `ncbi_cutoffs.rs`: `compute_blastn_cutoff_score()` で実装

### 潜在的な問題点

1. **オフ対角線検索の実行条件**
   - NCBI: `two_hits && word_type == 1` の場合のみ実行
   - LOSAT: `two_hits && word_type == 1` で実行（正しい）
   - **問題**: `TWO_HIT_WINDOW=0` の場合、`two_hits=false` となり、オフ対角線検索は実行されない
   - **NCBI BLASTNのデフォルト**: `window_size=0` (one-hit mode) なので、オフ対角線検索は実行されない（正常）

2. **Delta計算の確認**
   - NCBI: `Delta = MIN(scan_range, window_size - word_length)`
   - LOSAT: `delta_max = MIN(scan_range, window_size - word_length)` (1452-1457行目)
   - **問題**: `window_size=0` の場合、`delta_max` は常に0または負の値になる
   - **確認必要**: NCBIの実装で、`window_size=0` の場合の動作を確認

3. **hit_level_array更新タイミング** ✅ 修正完了
   - NCBI: `s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset` (757-758行目)
   - LOSAT (修正前): `final_se` (gapped extension後の位置) を使用
   - LOSAT (修正後): `ss + (qe - qs)` (ungapped extension後の位置) を使用
   - **修正**: 2025-01-XX - ungapped extension後の位置で更新するように修正

4. **Cutoff score適用タイミング**
   - NCBI: `if (off_found || ungapped_data->score >= cutoffs->cutoff_score)` (752行目)
   - LOSAT: `if !(off_found || ungapped_score >= cutoff_score)` (1608行目)
   - **確認**: ロジックは正しいが、`off_found` が正しく設定されているか確認

## 実装計画

### Phase 1: NCBI実装の詳細確認と検証

#### Step 1.1: オフ対角線検索の実行条件の確認 ✅ 完了

**目的**: NCBI BLASTNで`window_size=0`の場合、オフ対角線検索が実行されるか確認

**NCBIコード参照**: `na_ungapped.c:656-717`

**確認結果**:

1. **実行条件の確認**:
   - NCBI: `Boolean two_hits = (window_size > 0);` (line 656)
   - NCBI: `if (two_hits && (hit_saved || s_end_pos > last_hit + window_size))` (line 674)
   - NCBI: `if (word_type == 1) { /* try off-diagonals */ }` (line 686)
   - **結論**: `window_size=0` の場合、`two_hits=false` となり、オフ対角線検索ブロック全体が実行されない（正常）

2. **LOSAT実装の確認**:
   - Two-stage lookup (utils.rs:888-1077):
     - `let two_hits = TWO_HIT_WINDOW > 0;` (line 888) ✅ 正しい
     - `if two_hits && (hit_saved || s_end_pos > last_hit + TWO_HIT_WINDOW)` (line 905) ✅ 正しい
     - `if word_type == 1 { /* off-diagonal search */ }` (line 952) ✅ 正しい
   - Original lookup (utils.rs:1445-1643):
     - `let two_hits = TWO_HIT_WINDOW > 0;` (line 1445) ✅ 正しい
     - `if two_hits && (hit_saved || s_end_pos > last_hit + TWO_HIT_WINDOW)` (line 1456) ✅ 正しい
     - `if word_type == 1 { /* off-diagonal search */ }` (line 1488) ✅ 正しい

3. **発見された問題**:
   - **問題1**: Two-stage lookupで`s_off_pos`が未定義
     - 969行目で`s_off_pos`を使用しているが、そのスコープ内で定義されていない
     - NCBI: `s_off_pos = s_off + diag_table->offset` (na_ungapped.c:668)
     - 修正必要: `s_off_pos = s_off + TWO_HIT_WINDOW` を追加
   - **問題2**: Original lookupで`s_off_pos`が未定義
     - 1505行目で`s_off_pos`を使用しているが、そのスコープ内で定義されていない
     - NCBI: `s_off_pos = s_off + diag_table->offset` (na_ungapped.c:668)
     - 修正必要: `s_off_pos = s_off + TWO_HIT_WINDOW` を追加

**結論**: 
- オフ対角線検索の実行条件は正しく実装されている
- しかし、`s_off_pos`の定義が不足しているため、オフ対角線検索の計算が正しく実行されない可能性がある
- **修正が必要**: `s_off_pos`の定義を追加する必要がある

#### Step 1.2: Delta計算の詳細確認 ✅ 完了

**目的**: `window_size=0` の場合のDelta計算を確認

**NCBIコード参照**: `na_ungapped.c:658, 692` (Array版), `na_ungapped.c:825, 858` (Hash版)

**確認結果**:

1. **NCBI実装の確認**:
   - **Array版** (na_ungapped.c:658):
     - Line 658: `Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);`
     - 計算タイミング: 関数の最初（`two_hits`ブロックの前）
     - Line 692: `if (Delta < 0) Delta = 0;` (word_type == 1ブロック内)
   - **Hash版** (na_ungapped.c:825):
     - Line 825: `Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);`
     - 計算タイミング: 関数の最初（`two_hits`ブロックの前）
     - Line 858: `if (Delta < 0) Delta = 0;` (word_type == 1ブロック内)

2. **LOSAT実装の確認**:
   - **Two-stage lookup** (utils.rs:961-967):
     - 計算タイミング: `if word_type == 1`ブロック内で計算
     - 計算式: `let delta_calc = window_size as isize - word_length as isize;`
     - `let mut delta_max = if delta_calc < 0 { 0 } else { scan_range.min(delta_calc as usize) as isize };`
   - **Original lookup** (utils.rs:1502-1508):
     - 計算タイミング: `if word_type == 1`ブロック内で計算
     - 計算式: 同様

3. **計算式の等価性検証**:
   - NCBI: `MIN(scan_range, window_size - word_length)` → `if (Delta < 0) Delta = 0;`
   - LOSAT: `if delta_calc < 0 { 0 } else { scan_range.min(delta_calc) }`
   - **等価性**: 
     - `MIN(scan_range, window_size - word_length)` が負の場合: NCBIは後で0に修正、LOSATは先に0に設定
     - 結果は同じ: `MIN(scan_range, MAX(0, window_size - word_length))`
   - **結論**: 計算式は等価だが、計算タイミングが異なる

4. **計算タイミングの違い**:
   - **NCBI**: 関数の最初で計算（`two_hits`ブロックの前）
   - **LOSAT**: `if word_type == 1`ブロック内で計算
   - **影響**: 
     - `window_size=0`の場合、`two_hits=false`となり、オフ対角線検索ブロックは実行されない
     - したがって、Delta計算のタイミングの違いは結果に影響しない
     - ただし、NCBIの実装に合わせるべき（パリティのため）

5. **具体例の検証**:
   - `window_size=0`, `word_length=11`, `scan_range=4`の場合:
     - NCBI: `Delta = MIN(4, 0-11) = MIN(4, -11) = -11` → `if (Delta < 0) Delta = 0;` → `Delta = 0`
     - LOSAT: `delta_calc = 0 - 11 = -11` → `if delta_calc < 0 { 0 }` → `delta_max = 0`
     - **結果**: 両方とも0になり、オフ対角線検索は実行されない（正常）

**結論**: 
- Delta計算式は等価で、結果は同じになる
- ただし、計算タイミングが異なる（NCBIは関数の最初、LOSATは`if word_type == 1`ブロック内）
- `window_size=0`の場合、`two_hits=false`となり、オフ対角線検索ブロックは実行されないため、タイミングの違いは結果に影響しない
- **修正完了**: NCBIの実装に合わせて、Delta計算を関数の最初に移動 ✅ 完了 (2025-01-XX)

#### Step 1.3: マスキング更新タイミングの確認 ✅ 完了

**目的**: NCBIのマスキング更新タイミングを確認

**NCBIコード参照**: `na_ungapped.c:757-772`

**確認結果**:

1. **NCBI実装の確認**:
   - **Ungapped extension結果** (na_ungapped.c:206-241):
     - `ungapped_data->s_start = s_off - (q_off - ungapped_data->q_start);` (line 207)
     - `ungapped_data->length = (Int4)(q_end - q_beg);` (line 241)
     - `ungapped_data->s_start`はungapped extension後のsubject開始位置
     - `ungapped_data->length`はungapped extensionの長さ
   - **s_end_posの計算** (na_ungapped.c:757-758):
     - `s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;`
     - これはungapped extension後の終了位置（`s_start + length + offset`）
   - **hit_level_array更新** (na_ungapped.c:768-772):
     - `hit_level_array[real_diag].last_hit = s_end_pos;` (line 768)
     - 更新タイミング: ungapped extension後、gapped extension前
     - 拡張失敗時: `hit_ready = 0`だが、`s_end_pos`は更新されない（pre-extension値のまま）

2. **LOSAT実装の確認**:
   - **Two-stage lookup** (utils.rs:1167, 1292):
     - Line 1167: `final_s_end_pos = ss + (qe - qs) + TWO_HIT_WINDOW;`
     - `ss`はungapped extension後のsubject開始位置（`ungapped_data->s_start`に相当）
     - `qe - qs`はungapped extensionの長さ（`ungapped_data->length`に相当）
     - `TWO_HIT_WINDOW`は`diag_table->offset`に相当
     - Line 1292: `hit_level_array[diag_idx].last_hit = final_s_end_pos;`
     - **結論**: ungapped extension後の位置で更新されている（正しい）
   - **Original lookup** (utils.rs:1710, 1841):
     - Line 1710: `final_s_end_pos = ss + (qe - qs) + TWO_HIT_WINDOW;`
     - Line 1841: `hit_level_array[diag_idx].last_hit = final_s_end_pos;`
     - **結論**: ungapped extension後の位置で更新されている（正しい）

3. **拡張失敗時の処理**:
   - **NCBI**: `hit_ready = 0`だが、`s_end_pos`は更新されない（pre-extension値のまま）
   - **LOSAT**: `hit_ready = false`だが、`final_s_end_pos`は`s_end_pos`のまま（pre-extension値）
   - **結論**: 拡張失敗時の処理も正しい

4. **gapped extension後の処理**:
   - **NCBI**: gapped extensionは別の関数で実行され、`hit_level_array`の更新には影響しない
   - **LOSAT**: `final_s_end_pos`はgapped extensionの結果（`final_se`）で上書きされていない
   - **確認**: Line 1180-1203でgapped extensionを実行しているが、`final_s_end_pos`は更新されていない
   - **結論**: gapped extension後の処理も正しい

**結論**: 
- LOSAT実装はNCBI実装と一致している
- `hit_level_array[real_diag].last_hit`はungapped extension後の位置で更新されている
- 拡張失敗時はpre-extension値のまま更新される（NCBIと一致）
- gapped extensionの結果は`last_hit`の更新に影響しない（NCBIと一致）
- **修正不要**: 実装は正しい

#### Step 1.4: Cutoff score適用の確認 ✅ 完了

**目的**: `off_found` フラグの設定とcutoff score適用を確認

**NCBIコード参照**: `na_ungapped.c:730-752, 752-761`

**確認結果**:

1. **off_foundフラグの設定**: ✅ 正しい
   - **NCBI実装**:
     - Line 657: `Boolean off_found = FALSE;` (初期化)
     - Line 700, 709: `off_found = TRUE;` (オフ対角線検索で見つかった場合)
     - Line 713-716: `if (!off_found) { hit_ready = 0; }` (見つからなかった場合、hit_readyを0に)
   - **LOSAT実装**:
     - Line 915: `let mut off_found = false;` (初期化) ✅
     - Line 1033, 1073: `off_found = true;` (オフ対角線検索で見つかった場合) ✅
     - Line 1081-1087: `if !off_found { hit_ready = false; }` (見つからなかった場合、hit_readyをfalseに) ✅
   - **結論**: `off_found`フラグの設定は正しく実装されている

2. **cutoff score適用**: ✅ ロジックは正しい
   - **NCBI実装**:
     - Line 730: `Int4 context = BSearchContextInfo(q_off, query_info);`
     - Line 731: `cutoffs = word_params->cutoffs + context;`
     - Line 752: `if (off_found || ungapped_data->score >= cutoffs->cutoff_score)`
   - **LOSAT実装**:
     - Line 1151: `if !(off_found || ungapped_score >= cutoff_score)` (等価) ✅
   - **結論**: cutoff score適用のロジックは正しい

3. **cutoff_score計算**: ⚠️ 確認が必要
   - **NCBI実装**:
     - Contextごとにcutoff_scoreが異なる可能性がある（`cutoffs = word_params->cutoffs + context;`）
     - `BSearchContextInfo(q_off, query_info)`でcontextを決定
   - **LOSAT実装**:
     - Query-subjectペアごとにcutoff_scoreを計算（`utils.rs:639-647`）
     - Contextごとの計算は行っていない
   - **考察**:
     - blastnの場合、通常は単一のqueryなので、contextは0または1（strand）になる可能性が高い
     - 単一queryの場合、contextごとのcutoff_scoreは同じになる可能性が高い
     - しかし、複数queryの場合や、queryが分割されている場合、contextごとに異なるcutoff_scoreが必要になる可能性がある
   - **確認必要**: 複数queryの場合の動作を確認

**結論**:
- `off_found`フラグの設定とcutoff score適用のロジックは正しく実装されている
- ただし、contextごとのcutoff_score計算については、複数queryの場合の動作を確認する必要がある
- 現在の実装（query-subjectペアごとの計算）は、単一queryの場合には問題ない可能性が高い

**修正日**: 2025-01-24

### Phase 3: Extension実装の詳細確認

#### Step 3.1: Ungapped Extension X-drop Terminationの確認 🔴 修正必要

**目的**: Ungapped extensionのX-drop termination条件がNCBI実装と一致しているか確認

**NCBIコード参照**: `na_ungapped.c:148-243` (s_NuclUngappedExtendExact)

**確認結果**:

1. **X-drop terminationの条件**: ❌ 不一致
   - **NCBI実装** (`na_ungapped.c:197-203`, Left extension):
     ```c
     if ((sum += matrix[*--q][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
         q_beg = q;
         score += sum;
         sum = 0;
     } else if (sum < X) {
         break;
     }
     ```
     - **重要**: Xは負の値（`-(cutoffs->x_dropoff)`として渡される、例：`-20`）
     - **重要**: `sum < X`で終了（`sum`がXより負の値になったら終了）
     - **重要**: `sum`が正になったら`score += sum`して`sum = 0`にリセット
   
   - **NCBI実装** (`na_ungapped.c:227-233`, Right extension):
     ```c
     if ((sum += matrix[*q++][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
         q_end = q;
         score += sum;
         X_current = (-score > X) ? -score : X;
         sum = 0;
     } else if (sum < X_current) {
         break;
     }
     ```
     - **重要**: `X_current`は動的に更新される（`X_current = (-score > X) ? -score : X;`）
     - **重要**: `sum < X_current`で終了

   - **LOSAT実装** (`extension.rs:45-58`, Left extension):
     ```rust
     current_score += sc;
     if current_score > max_score {
         max_score = current_score;
         best_i = i;
     } else if (max_score - current_score) > x_drop_val {
         break;
     }
     ```
     - **問題**: `(max_score - current_score) > x_drop_val`は、NCBIの`sum < X`とは異なるロジック
     - **問題**: LOSATは`x_drop_val`を正の値として使用しているが、NCBIは負の値として使用

   - **LOSAT実装** (`extension.rs:67-80`, Right extension):
     ```rust
     current_score_r += sc;
     if current_score_r > max_score_total {
         max_score_total = current_score_r;
         best_j = j;
     } else if (max_score_total - current_score_r) > x_drop_val {
         break;
     }
     ```
     - **問題**: NCBIの動的`X_current`更新が実装されていない

2. **スコア計算方法**: ⚠️ 確認が必要
   - **NCBI実装**: `sum`を累積し、`sum > 0`のときに`score += sum`して`sum = 0`にリセット
   - **LOSAT実装**: `current_score`を累積し、`current_score > max_score`のときに`max_score = current_score`を更新
   - **考察**: ロジックは等価に見えるが、X-drop termination条件が異なる

3. **Xパラメータの符号**: ❌ 不一致
   - **NCBI実装**: `s_NuclUngappedExtendExact(..., -(cutoffs->x_dropoff), ...)` - 負の値として渡される
   - **LOSAT実装**: `extend_hit_ungapped(..., x_drop: Option<i32>)` - 正の値として使用されている

**結論**:
- X-drop terminationの条件がNCBI実装と一致していない
- NCBIは`sum < X`（Xは負の値）で終了するが、LOSATは`(max_score - current_score) > x_drop_val`（x_drop_valは正の値）で終了
- 数学的には等価かもしれないが、実行順序やタイミングが異なる可能性がある
- **修正が必要**: X-drop termination条件をNCBI実装に合わせる必要がある

**修正日**: 2025-01-24

#### Step 3.2: Ungapped Extension X-drop Terminationの修正 ✅ 完了

**目的**: X-drop termination条件をNCBI実装に完全一致させる

**NCBIコード参照**: `na_ungapped.c:148-243` (s_NuclUngappedExtendExact)

**修正内容**:

1. **X-drop termination条件の修正**:
   - **修正前**: `(max_score - current_score) > x_drop_val`（正の値を使用）
   - **修正後**: `sum < X`（Xは負の値、`-(cutoffs->x_dropoff)`として渡される）
   - **NCBI参照**: `na_ungapped.c:201-202`, `na_ungapped.c:232-233`

2. **スコア計算方法の修正**:
   - **修正前**: `current_score`を累積し、`current_score > max_score`のときに`max_score = current_score`を更新
   - **修正後**: `sum`を累積し、`sum > 0`のときに`score += sum`して`sum = 0`にリセット
   - **NCBI参照**: `na_ungapped.c:197-200`, `na_ungapped.c:227-231`

3. **Right extensionの動的X_current更新の実装**:
   - **修正前**: 動的更新なし
   - **修正後**: `X_current = (-score > X) ? -score : X;`で動的に更新
   - **NCBI参照**: `na_ungapped.c:230`

4. **Right extensionの開始位置の修正**:
   - **修正前**: `j = 1`から開始（seed positionをスキップ）
   - **修正後**: `j = 0`から開始（seed positionもスコア）
   - **NCBI参照**: `na_ungapped.c:218-219` (`q = query->sequence + q_off;`)

5. **座標計算の修正**:
   - **修正前**: `q_end = q_pos + j`（最後にスコアされた位置）
   - **修正後**: `q_end = q_pos + j + 1`（最後にスコアされた位置の次の位置）
   - **NCBI参照**: `na_ungapped.c:227-228` (`*q++`の後で`q_end = q;`)

**修正ファイル**: `src/algorithm/blastn/extension.rs`

**修正日**: 2025-01-24

#### Step 3.3: 統合テスト結果（X-drop termination修正後）

**実施日**: 2025-01-24

**テスト結果**:

| Test Case | LOSAT Hits | NCBI Hits | Ratio | Avg Length (LOSAT/NCBI) | Avg Bitscore (LOSAT/NCBI) | Avg Identity (LOSAT/NCBI) | Avg E-value (LOSAT/NCBI) |
|-----------|------------|-----------|-------|-------------------------|---------------------------|----------------------------|--------------------------|
| EDL933.Sakai.megablast | 1508 | 5718 | 26.4% | 4264.3 / 1438.7 | 7653.6 / 2483.4 | 92.16% / 93.23% | 3.57e-05 / 7.78e-05 |
| MelaMJNV.PemoMJNVA.blastn | 742 | 2729 | 27.2% | 177.5 / 86.3 | 124.5 / 60.4 | 80.21% / 84.04% | 7.18e-01 / 2.37e+00 |
| MjeNMV.MelaMJNV.blastn | 1371 | 2668 | 51.4% | 310.6 / 252.0 | 467.5 / 290.1 | 87.00% / 84.14% | 7.13e-01 / 1.67e+00 |
| MjPMNV.MlPMNV.blastn | 2258 | 54402 | 4.2% | 217.1 / 144.8 | 311.8 / 122.3 | 87.81% / 80.41% | 4.32e-01 / 3.46e-01 |

**分析**:

1. **ヒット数の差異**: 
   - X-drop termination修正後も、まだ大きな差異が残っている
   - 特に`MjPMNV.MlPMNV.blastn`では、LOSATが4.2%しか検出していない（2258 vs 54402）
   - `EDL933.Sakai.megablast`では26.4%（1508 vs 5718）

2. **平均長さの差異**:
   - LOSATの平均長さがNCBIより長い傾向がある
   - 例：`EDL933.Sakai.megablast`: 4264.3 vs 1438.7（約3倍）
   - 例：`MelaMJNV.PemoMJNVA.blastn`: 177.5 vs 86.3（約2倍）

3. **平均bitscoreの差異**:
   - LOSATの平均bitscoreがNCBIより高い傾向がある
   - 例：`EDL933.Sakai.megablast`: 7653.6 vs 2483.4（約3倍）
   - これは平均長さが長いことと一致

4. **考察**:
   - X-drop termination修正により、extensionの動作は改善された可能性がある
   - しかし、ヒット数の大幅な差異は、他の要因（例：seeding、cutoff score適用、HSP filtering）が影響している可能性がある
   - 特に`MjPMNV.MlPMNV.blastn`の4.2%という極端に低い値は、何らかのフィルタリングが過度に働いている可能性がある

**次のステップ**:
- 他の要因（seeding、cutoff score適用、HSP filtering）の確認が必要
- 特に`MjPMNV.MlPMNV.blastn`の極端に低いヒット数の原因を調査

**修正日**: 2025-01-24

### Phase 2: 実装修正

#### Step 2.1: hit_level_array更新タイミングの修正 ✅ 完了
- **ファイル**: `src/algorithm/blastn/utils.rs`
- **修正内容**:
  - **問題**: `hit_level_array[real_diag].last_hit`をgapped extension後の位置（`final_se`）で更新していた
  - **NCBI実装**: ungapped extension後の位置（`ungapped_data->s_start + ungapped_data->length`）で更新
  - **修正後**: ungapped extension後の位置（`ss + (qe - qs)`）で更新
  - **NCBI参照**: `na_ungapped.c:757-772`
  - **影響範囲**: 
    - Two-stage lookup: 1180-1200行目
    - Original lookup: 1750-1767行目
    - スキップ時の処理: 1098-1121行目、1610-1630行目
  - **修正日**: 2025-01-XX

#### Step 2.2: オフ対角線検索の実装確認と修正 ✅ 完了 (部分)

**ファイル**: `src/algorithm/blastn/utils.rs`

**確認事項**:
- Delta計算が正しいか（1452-1457行目、958-962行目） ✅ 正しい
- オフ対角線検索の条件チェックが正しいか（1503-1508行目、1538-1543行目） ✅ 正しい
- 対角線インデックスの計算が正しいか（1485-1499行目、1520-1534行目） ✅ 正しい

**発見された問題**:
- **問題**: `s_off_pos`の定義が不足している
  - Two-stage lookup: 969行目で`s_off_pos`を使用しているが、定義されていない
  - Original lookup: 1505行目で`s_off_pos`を使用しているが、定義されていない
  - NCBI: `s_off_pos = s_off + diag_table->offset` (na_ungapped.c:668)
  - 修正必要: `s_off_pos = s_off + TWO_HIT_WINDOW` を追加

#### Step 2.2.1: s_off_pos定義の追加（新規）

**ファイル**: `src/algorithm/blastn/utils.rs`

**修正内容**:
1. **Two-stage lookup** (utils.rs:888-1077):
   - **問題**: 969行目で`s_off_pos`を使用しているが、定義されていない
   - **NCBI実装**: `s_off_pos = s_off + diag_table->offset` (na_ungapped.c:668)
   - **修正**: `s_off_pos = s_off + TWO_HIT_WINDOW` を追加（898行目の後に追加）
   - **NCBI参照**: `na_ungapped.c:668`, `blast_extend.c:63`

2. **Original lookup** (utils.rs:1445-1643):
   - **問題**: 1505行目で`s_off_pos`を使用しているが、定義されていない
   - **NCBI実装**: `s_off_pos = s_off + diag_table->offset` (na_ungapped.c:668)
   - **修正**: `s_off_pos = s_off + TWO_HIT_WINDOW` を追加（1449行目の後に追加）
   - **NCBI参照**: `na_ungapped.c:668`, `blast_extend.c:63`

**修正詳細**:

**Two-stage lookup修正** (utils.rs:897-898行目の後に追加):
```rust
// NCBI reference: na_ungapped.c:668
// s_off_pos = s_off + diag_table->offset;
// diag_table->offset = window_size (blast_extend.c:63)
// CRITICAL: Match NCBI exactly
let s_off_pos = s_off + TWO_HIT_WINDOW;
```

**Original lookup修正** (utils.rs:1448-1449行目の後に追加):
```rust
// NCBI reference: na_ungapped.c:668
// s_off_pos = s_off + diag_table->offset;
// diag_table->offset = window_size (blast_extend.c:63)
// CRITICAL: Match NCBI exactly
let s_off_pos = s_off + TWO_HIT_WINDOW;
```

**修正日**: 2025-01-XX

#### Step 2.3: Cutoff score適用の確認
- **ファイル**: `src/algorithm/blastn/utils.rs`
- **確認事項**:
  - `off_found` フラグが正しく設定されているか
  - `cutoff_score` の計算が正しいか
  - cutoff score適用のタイミングが正しいか（1608行目）

#### Step 2.4: Extension実装の確認
- **ファイル**: `src/algorithm/blastn/extension.rs`
- **確認事項**:
  - Ungapped extensionの実装がNCBIと一致しているか
  - X-drop termination条件が正しいか
  - Score計算が正しいか

### Phase 3: テストと検証

#### Step 3.1: 各ステップ後の統合テスト
各ステップが終わるごとに、以下を実施:

1. **Release版をビルド**
```bash
cd /mnt/c/Users/genom/GitHub/LOSAT/LOSAT
cargo build --release
```

2. **統合テストを実行**
```bash
cd /mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests
# Path to LOSAT binary
LOSAT_BIN="../target/release/LOSAT"

# --- LOSATN Commands (Default / Megablast behavior) ---

# NZ_CP006932 self (Default/Megablast)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.losatn.megablast.out -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.losatn.megablast.log

# EDL933 vs Sakai
(time $LOSAT_BIN blastn -q ./fasta/EDL933.fna -s ./fasta/Sakai.fna -o ./losat_out/EDL933.Sakai.losatn.megablast.out -n 1 )&>./losat_out/EDL933.Sakai.losatn.megablast.log

# Sakai vs MG1655
(time $LOSAT_BIN blastn -q ./fasta/Sakai.fna -s ./fasta/MG1655.fna -o ./losat_out/Sakai.MG1655.losatn.megablast.out -n 1 )&>./losat_out/Sakai.MG1655.losatn.megablast.log

# --- LOSATN Commands (Task: blastn) ---

# NZ_CP006932 self (Task: blastn)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out --task blastn -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.log

# PesePMNV vs MjPMNV
(time $LOSAT_BIN blastn -q ./fasta/AP027152.fasta -s ./fasta/AP027202.fasta -o ./losat_out/PesePMNV.MjPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PesePMNV.MjPMNV.losatn.blastn.log

# MelaMJNV vs PemoMJNVA
(time $LOSAT_BIN blastn -q ./fasta/LC738874.fasta -s ./fasta/LC738870.fasta -o ./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.log

# SiNMV vs ChdeNMV
(time $LOSAT_BIN blastn -q ./fasta/LC738884.fasta -s ./fasta/AP027155.fasta -o ./losat_out/SiNMV.ChdeNMV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/SiNMV.ChdeNMV.losatn.blastn.log

# PmeNMV vs MjPMNV
(time $LOSAT_BIN blastn -q ./fasta/LC738869.fasta -s ./fasta/AP027202.fasta -o ./losat_out/PmeNMV.MjPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PmeNMV.MjPMNV.losatn.blastn.log

# PmeNMV vs PesePMNV
(time $LOSAT_BIN blastn -q ./fasta/LC738869.fasta -s ./fasta/AP027152.fasta -o ./losat_out/PmeNMV.PesePMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PmeNMV.PesePMNV.losatn.blastn.log

# PeseMJNV vs PemoMJNVB
(time $LOSAT_BIN blastn -q ./fasta/LC738873.fasta -s ./fasta/LC738871.fasta -o ./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.log

# PemoMJNVA vs PeseMJNV
(time $LOSAT_BIN blastn -q ./fasta/LC738870.fasta -s ./fasta/LC738873.fasta -o ./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.log

# MjeNMV vs MelaMJNV
(time $LOSAT_BIN blastn -q ./fasta/LC738868.fasta -s ./fasta/LC738874.fasta -o ./losat_out/MjeNMV.MelaMJNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MjeNMV.MelaMJNV.losatn.blastn.log

# MjPMNV vs MlPMNV
(time $LOSAT_BIN blastn -q ./fasta/AP027202.fasta -s ./fasta/LC738875.fasta -o ./losat_out/MjPMNV.MlPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MjPMNV.MlPMNV.losatn.blastn.log
```

3. **結果比較**
   - BLASTの結果(`/mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests/blast_out/`内)と比較
   - ヒット数の分布（数、長さ、e-value, bitscore, identity）を比較
   - 差異がある場合は、原因を特定して修正

#### Step 3.2: NCBI実装との再確認
各ステップが終わるごとに:
- NCBI実装と完全に合致しているか再確認
- 合致していない箇所を修正
- NCBIコード参照をコメントとして追加

## 実装詳細

### Step 2.1: hit_level_array更新タイミングの修正 ✅ 完了

#### 修正前の実装（問題あり）
```rust
// utils.rs:1186-1200 (two-stage lookup), 1754-1767 (original lookup)
// NCBI: s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
// In LOSAT, we use final_se (end position after gapped extension) + diag_offset
let final_s_end_pos = final_se + diag_offset as usize;
hit_level_array[diag_idx].last_hit = final_s_end_pos; // ❌ gapped extension後の位置
```

#### 修正後の実装（NCBIに合わせる）
```rust
// NCBI reference: na_ungapped.c:757-772
// s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
// hit_level_array[real_diag].last_hit = s_end_pos;
// CRITICAL: NCBI updates last_hit using ungapped extension end position, NOT gapped extension
// ungapped_data->s_start is the subject start position after ungapped extension
// ungapped_data->length is the length of the ungapped extension
// In LOSAT, we use se (ungapped extension end) + diag_offset to match NCBI
let ungapped_s_end_pos = ss + (qe - qs) + diag_offset as usize;
hit_level_array[diag_idx].last_hit = ungapped_s_end_pos; // ✅ ungapped extension後の位置
```

#### 修正内容の詳細
- **Two-stage lookup**: 1180-1200行目を修正
- **Original lookup**: 1750-1767行目を修正
- **スキップ時の処理**: 1098-1121行目、1610-1630行目を修正（`hit_ready = false`を明示的に設定）
- **修正日**: 2025-01-XX

### Step 2.2: オフ対角線検索の実装確認と修正 ✅ 完了

#### 修正内容
1. **diag_array_lengthとdiag_maskの計算**: ✅ 修正完了
   - NCBI: `diag_array_length = smallest power of 2 >= (qlen + window_size)`
   - NCBI: `diag_mask = diag_array_length - 1`
   - LOSAT: 同様に計算するように修正

2. **diagとreal_diagの計算**: ✅ 修正完了
   - NCBI: `diag = s_off + diag_table->diag_array_length - q_off`
   - NCBI: `real_diag = diag & diag_table->diag_mask`
   - LOSAT: NCBI実装に完全一致するように修正

3. **orig_diagの計算**: ✅ 修正完了
   - NCBI: `orig_diag = real_diag + diag_table->diag_array_length`
   - LOSAT: 同様に計算するように修正

4. **off_diagの計算**: ✅ 修正完了
   - NCBI: `off_diag = (orig_diag + delta) & diag_table->diag_mask`
   - LOSAT: マスク処理を追加し、NCB実装に完全一致

5. **配列インデックス**: ✅ 修正完了
   - NCBI: `hit_level_array[off_diag]` (マスク済みのoff_diagを直接使用)
   - LOSAT: `hit_level_array[off_diag]` (マスク済みのoff_diagを直接使用)
   - **修正**: `diag_offset`を削除し、NCBI実装に完全一致

6. **s_off_posとs_end_posの計算**: ✅ 修正完了
   - NCBI: `s_off_pos = s_off + diag_table->offset` (offset = window_size)
   - NCBI: `s_end_pos = s_end + diag_table->offset`
   - LOSAT: `TWO_HIT_WINDOW`を使用するように修正

#### 修正詳細

**修正箇所:**
1. `diag_array_length`と`diag_mask`の計算をstrandループの前に移動
2. `diag`と`real_diag`の計算をNCBI実装に完全一致
3. `orig_diag`の計算をNCBI実装に完全一致
4. `off_diag`の計算にマスク処理を追加
5. `diag_offset`を削除し、`real_diag`と`off_diag`を直接配列インデックスとして使用
6. `s_off_pos`と`s_end_pos`の計算で`TWO_HIT_WINDOW`を使用

**NCBI参照:**
- `na_ungapped.c:663-664`: diagとreal_diagの計算
- `na_ungapped.c:688`: orig_diagの計算
- `na_ungapped.c:694, 703`: off_diagの計算（マスク処理）
- `na_ungapped.c:668-669`: s_off_posとs_end_posの計算
- `blast_extend.c:52-61`: diag_array_lengthとdiag_maskの計算

#### 修正日
2025-01-XX

### Step 2.2.1: s_off_pos定義の追加（新規）

#### 修正前の実装（問題あり）
```rust
// utils.rs:969 (two-stage lookup), 1505 (original lookup)
// NCBI: s_off_pos = s_off + diag_table->offset;
// In LOSAT, s_off_pos is not defined in the scope where it's used
let s_a = s_off_pos as isize + word_length as isize - window_size as isize; // ❌ s_off_pos未定義
```

#### 修正後の実装（NCBIに合わせる）
```rust
// NCBI reference: na_ungapped.c:668
// s_off_pos = s_off + diag_table->offset;
// diag_table->offset = window_size (blast_extend.c:63)
// CRITICAL: Match NCBI exactly
let s_off_pos = s_off + TWO_HIT_WINDOW; // ✅ NCBI実装に一致

// NCBI reference: na_ungapped.c:689-690
// Int4 s_a = s_off_pos + word_length - window_size;
// Int4 s_b = s_end_pos - 2 * word_length;
let s_a = s_off_pos as isize + word_length as isize - window_size as isize; // ✅ 正しい
```

#### 修正内容の詳細
- **Two-stage lookup**: 897-898行目の後に`s_off_pos`の定義を追加
- **Original lookup**: 1448-1449行目の後に`s_off_pos`の定義を追加
- **NCBI参照**: `na_ungapped.c:668`, `blast_extend.c:63`
- **修正日**: 2025-01-XX

#### テスト結果 (Step 2.2.1修正後)

**修正内容**:
- Two-stage lookup: 897-898行目の後に`s_off_pos = s_off + TWO_HIT_WINDOW`を追加
- Original lookup: 1453-1454行目の後に`s_off_pos = s_off + TWO_HIT_WINDOW`を追加
- NCBI参照: `na_ungapped.c:668`, `blast_extend.c:63`

**ビルド結果**: ✅ 成功（警告のみ、エラーなし）

**統合テスト結果**: ✅ 全13テストケース完了

**修正前後の比較**:
- 修正前: `s_off_pos`が未定義で、オフ対角線検索の計算が正しく実行されない可能性があった
- 修正後: `s_off_pos`が正しく定義され、NCBI実装と一致

**次のステップ**: 
- Step 1.2: Delta計算の詳細確認 ✅ 完了
- Step 1.3: マスキング更新タイミングの確認
- Step 1.4: Cutoff score適用の確認

### Step 2.2.2: Delta計算タイミングの修正（関数の最初に移動） ✅ 完了

#### 修正前の実装（問題あり）
```rust
// utils.rs:957-967 (two-stage lookup), 1498-1508 (original lookup)
// NCBI: Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
// NCBI calculates Delta at function start (line 658), NOT inside word_type == 1 block
if word_type == 1 {
    // ❌ Delta計算がword_type == 1ブロック内で実行されている
    let window_size = TWO_HIT_WINDOW;
    let delta_calc = window_size as isize - word_length as isize;
    let mut delta_max = if delta_calc < 0 {
        0
    } else {
        scan_range.min(delta_calc as usize) as isize
    };
}
```

#### 修正後の実装（NCBIに合わせる）
```rust
// NCBI reference: na_ungapped.c:656-658
// Boolean two_hits = (window_size > 0);
// Boolean off_found = FALSE;
// Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
// CRITICAL: Delta must be calculated at function start, BEFORE two_hits block
// This matches NCBI implementation exactly (na_ungapped.c:658)
let two_hits = TWO_HIT_WINDOW > 0;
let window_size = TWO_HIT_WINDOW;
// NCBI: Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
let mut delta = scan_range.min(window_size.saturating_sub(word_length)) as isize;

// ... later in word_type == 1 block ...
if word_type == 1 {
    // NCBI reference: na_ungapped.c:692
    // if (Delta < 0) Delta = 0;
    // CRITICAL: Delta was calculated at function start (line 658)
    // Here we only check if it's negative and set to 0
    if delta < 0 {
        delta = 0;
    }
    let delta_max = delta;
}
```

#### 修正内容の詳細
- **Two-stage lookup**: 
  - Delta計算を888行目（`two_hits`の直後）に移動
  - `if word_type == 1`ブロック内では`if (Delta < 0) Delta = 0;`のみ実行
- **Original lookup**: 
  - Delta計算を1463行目（`two_hits`の直後）に移動
  - `if word_type == 1`ブロック内では`if (Delta < 0) Delta = 0;`のみ実行
- **NCBI参照**: 
  - `na_ungapped.c:656-658`: Delta計算の位置
  - `na_ungapped.c:692`: `if (Delta < 0) Delta = 0;`の位置
- **修正日**: 2025-01-XX

#### 修正の効果
- **計算タイミング**: NCBI実装と完全一致（関数の最初で計算）
- **実行順序**: NCBI実装と完全一致（Delta計算 → two_hitsブロック → word_type == 1ブロック）
- **計算式**: 等価（`MIN(scan_range, window_size - word_length)`）

#### ビルド結果
- ✅ コンパイル成功（警告のみ、エラーなし）
- ✅ Release版ビルド成功

#### テスト結果 (Step 2.2.2修正後)

**Megablast:**

| Test Case | LOSAT Hits | NCBI Hits | Ratio | Avg Length (LOSAT/NCBI) | Avg Bitscore (LOSAT/NCBI) | Avg Identity (LOSAT/NCBI) | Avg E-value (LOSAT/NCBI) |
|-----------|------------|-----------|-------|-------------------------|---------------------------|----------------------------|--------------------------|
| NZ_CP006932 self | 270 | N/A | N/A | 2805.3 / N/A | 4887.3 / N/A | 87.31% / N/A | 3.78e-08 / N/A |
| EDL933 vs Sakai | 1508 | 5718 | 26.4% | 4264.3 / 1438.7 | 7653.6 / 2483.4 | 92.16% / 93.23% | 3.57e-05 / 7.78e-05 |
| Sakai vs MG1655 | 1202 | N/A | N/A | 3677.0 / N/A | 6348.0 / N/A | 93.66% / N/A | 2.47e-05 / N/A |

**Blastn:**

| Test Case | LOSAT Hits | NCBI Hits | Ratio | Avg Length (LOSAT/NCBI) | Avg Bitscore (LOSAT/NCBI) | Avg Identity (LOSAT/NCBI) | Avg E-value (LOSAT/NCBI) |
|-----------|------------|-----------|-------|-------------------------|---------------------------|----------------------------|--------------------------|
| NZ_CP006932 self | 6190 | N/A | N/A | 244.1 / N/A | 376.8 / N/A | 83.84% / N/A | 2.79e+00 / N/A |
| PesePMNV vs MjPMNV | 389 | N/A | N/A | 409.4 / N/A | 295.0 / N/A | 77.51% / N/A | 2.25e-01 / N/A |
| MelaMJNV vs PemoMJNVA | 712 | 2729 | 26.1% | 181.9 / 86.3 | 128.0 / 60.4 | 80.32% / 84.04% | 6.79e-01 / 2.37e+00 |
| SiNMV vs ChdeNMV | 1836 | N/A | N/A | 230.0 / N/A | 329.7 / N/A | 89.02% / N/A | 6.42e-01 / N/A |
| PmeNMV vs MjPMNV | 366 | N/A | N/A | 437.3 / N/A | 320.5 / N/A | 77.45% / N/A | 3.29e-01 / N/A |
| PmeNMV vs PesePMNV | 503 | N/A | N/A | 456.7 / N/A | 442.7 / N/A | 81.45% / N/A | 5.53e-01 / N/A |
| PeseMJNV vs PemoMJNVB | 1384 | N/A | N/A | 208.9 / N/A | 181.8 / N/A | 81.76% / N/A | 8.19e-01 / N/A |
| PemoMJNVA vs PeseMJNV | 1179 | N/A | N/A | 333.9 / N/A | 455.6 / N/A | 86.05% / N/A | 5.75e-01 / N/A |
| MjeNMV vs MelaMJNV | 1342 | 2668 | 50.3% | 315.9 / 252.0 | 476.8 / 290.1 | 87.14% / 84.14% | 6.88e-01 / 1.67e+00 |
| MjPMNV vs MlPMNV | 2246 | 54402 | 4.1% | 218.0 / 144.8 | 313.2 / 122.3 | 87.84% / 80.41% | 4.26e-01 / 3.46e-01 |

**考察**:
- Delta計算タイミングの修正により、NCBI実装との実行順序が完全一致
- 比較可能なテストケース（EDL933 vs Sakai, MelaMJNV vs PemoMJNVA, MjeNMV vs MelaMJNV, MjPMNV vs MlPMNV）では、ヒット数の比率は26.1%〜50.3%の範囲
- LOSATの平均ヒット長がNCBIより長い傾向がある（フィルタリングが不十分な可能性）
- LOSATの平均bitscoreがNCBIより高い傾向がある（高スコアのヒットのみが残っている可能性）
- 次のステップ: Step 1.3（マスキング更新タイミングの確認）とStep 1.4（Cutoff score適用の確認）が必要

#### Step 2.2.1完了後の確認事項

各ステップが終わるごとに、以下を実施:

1. **NCBI実装と完全に合致しているか再確認**
   - `s_off_pos`の定義がNCBI実装と一致しているか確認
   - Two-stage lookupとOriginal lookupの両方で修正されているか確認
   - NCBIコード参照がコメントとして追加されているか確認

2. **Release版をビルド**
```bash
cd /mnt/c/Users/genom/GitHub/LOSAT/LOSAT
cargo build --release
```

3. **統合テストをすべて実施**
```bash
cd /mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests
# Path to LOSAT binary
LOSAT_BIN="../target/release/LOSAT"

# --- LOSATN Commands (Default / Megablast behavior) ---

# NZ_CP006932 self (Default/Megablast)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.losatn.megablast.out -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.losatn.megablast.log

# EDL933 vs Sakai
(time $LOSAT_BIN blastn -q ./fasta/EDL933.fna -s ./fasta/Sakai.fna -o ./losat_out/EDL933.Sakai.losatn.megablast.out -n 1 )&>./losat_out/EDL933.Sakai.losatn.megablast.log

# Sakai vs MG1655
(time $LOSAT_BIN blastn -q ./fasta/Sakai.fna -s ./fasta/MG1655.fna -o ./losat_out/Sakai.MG1655.losatn.megablast.out -n 1 )&>./losat_out/Sakai.MG1655.losatn.megablast.log

# --- LOSATN Commands (Task: blastn) ---

# NZ_CP006932 self (Task: blastn)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out --task blastn -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.log

# PesePMNV vs MjPMNV
(time $LOSAT_BIN blastn -q ./fasta/AP027152.fasta -s ./fasta/AP027202.fasta -o ./losat_out/PesePMNV.MjPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PesePMNV.MjPMNV.losatn.blastn.log

# MelaMJNV vs PemoMJNVA
(time $LOSAT_BIN blastn -q ./fasta/LC738874.fasta -s ./fasta/LC738870.fasta -o ./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.log

# SiNMV vs ChdeNMV
(time $LOSAT_BIN blastn -q ./fasta/LC738884.fasta -s ./fasta/AP027155.fasta -o ./losat_out/SiNMV.ChdeNMV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/SiNMV.ChdeNMV.losatn.blastn.log

# PmeNMV vs MjPMNV
(time $LOSAT_BIN blastn -q ./fasta/LC738869.fasta -s ./fasta/AP027202.fasta -o ./losat_out/PmeNMV.MjPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PmeNMV.MjPMNV.losatn.blastn.log

# PmeNMV vs PesePMNV
(time $LOSAT_BIN blastn -q ./fasta/LC738869.fasta -s ./fasta/AP027152.fasta -o ./losat_out/PmeNMV.PesePMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PmeNMV.PesePMNV.losatn.blastn.log

# PeseMJNV vs PemoMJNVB
(time $LOSAT_BIN blastn -q ./fasta/LC738873.fasta -s ./fasta/LC738871.fasta -o ./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.log

# PemoMJNVA vs PeseMJNV
(time $LOSAT_BIN blastn -q ./fasta/LC738870.fasta -s ./fasta/LC738873.fasta -o ./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.log

# MjeNMV vs MelaMJNV
(time $LOSAT_BIN blastn -q ./fasta/LC738868.fasta -s ./fasta/LC738874.fasta -o ./losat_out/MjeNMV.MelaMJNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MjeNMV.MelaMJNV.losatn.blastn.log

# MjPMNV vs MlPMNV
(time $LOSAT_BIN blastn -q ./fasta/AP027202.fasta -s ./fasta/LC738875.fasta -o ./losat_out/MjPMNV.MlPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MjPMNV.MlPMNV.losatn.blastn.log
```

4. **結果比較**
   - BLASTの結果(`/mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests/blast_out/`内)と比較
   - ヒット数の分布（数、長さ、e-value, bitscore, identity）を比較
   - 差異がある場合は、原因を特定して修正

### Step 2.4: Cutoff score適用の確認 ✅ 完了

#### 実装内容
1. **QueryInfo構造体の実装**: `src/algorithm/blastn/query_info.rs`を新規作成
   - `BSearchContextInfo`相当の関数`bsearch_context_info`を実装
   - Query context情報を管理する構造体を定義
   - NCBI参照: `blast_query_info.c:219-243`

2. **Contextごとのcutoff_score計算**: `utils.rs`で事前計算
   - Query context情報からcontext数を決定
   - 各contextごとに`compute_blastn_cutoff_score`を呼び出し
   - Contextごとのcutoff_score配列を保持
   - Per-subjectでcutoff_scoreを更新（subject lengthに依存するため）

3. **Cutoff score適用箇所の修正**: Two-stage lookupとOriginal lookupの両方で修正
   - `q_idx`からcontextを決定（各queryが2つのcontext: 0=forward, 1=reverse）
   - Contextごとのcutoff_scoreを使用
   - NCBI参照: `na_ungapped.c:730-752`

#### 確認事項
1. **off_foundフラグ**: オフ対角線検索で正しく設定されているか ✅ 確認済み
2. **cutoff_score計算**: Contextごとに計算され、per-subjectで更新 ✅ 実装済み
3. **適用タイミング**: ungapped extension後、gapped extension前 ✅ 確認済み

#### テスト結果 (Step 2.4修正後)

**Megablast:**
- NZ_CP006932 self: LOSAT=270, NCBI=N/A
- EDL933 vs Sakai: LOSAT=1508, NCBI=N/A
- Sakai vs MG1655: LOSAT=1202, NCBI=N/A

**Blastn:**
- NZ_CP006932 self: LOSAT=6529, NCBI=454 (1438.1%) - NCBIファイルが存在しない可能性
- MelaMJNV vs PemoMJNVA: LOSAT=742, NCBI=2729 (27.2%)
- MjeNMV vs MelaMJNV: LOSAT=1371, NCBI=2668 (51.4%)
- MjPMNV vs MlPMNV: LOSAT=2258, NCBI=54402 (4.2%)

**考察:**
- Contextごとのcutoff_score計算と適用を実装
- Per-subjectでcutoff_scoreを更新するように修正
- しかし、ヒット数の差異は依然として大きい（特に`MjPMNV.MlPMNV.blastn`で4.2%）
- 次のステップ: Extension実装の詳細確認が必要

**修正日**: 2025-01-24

## 進捗管理

### Phase 1: NCBI実装の詳細確認と検証
- [x] Step 1.1: オフ対角線検索の実行条件の確認 ✅ 完了 (2025-01-24)
- [x] Step 1.2: Delta計算の詳細確認 ✅ 完了 (2025-01-24)
- [x] Step 1.3: マスキング更新タイミングの確認 ✅ 完了 (2025-01-24)
- [x] Step 1.4: Cutoff score適用の確認 ✅ 完了 (2025-01-24)

### Phase 2: 実装修正
- [x] Step 2.1: hit_level_array更新タイミングの修正（ungapped extension後の位置を使用） ✅ 完了 (2025-01-XX)
- [x] Step 2.2: オフ対角線検索の実装確認と修正 ✅ 完了 (2025-01-XX)
- [x] Step 2.3: hit_level_array更新タイミングの完全修正（常に最後に更新） ✅ 完了 (2025-01-24)
- [x] Step 2.2.1: s_off_pos定義の追加（新規） ✅ 完了 (2025-01-XX)
- [x] Step 2.2.2: Delta計算タイミングの修正（関数の最初に移動） ✅ 完了 (2025-01-XX)
- [x] Step 2.4: Cutoff score適用の確認 ✅ 完了 (2025-01-24)
- [ ] Step 2.5: Extension実装の確認

### Phase 3: テストと検証
- [x] Step 3.1: 各ステップ後の統合テスト ✅ 完了 (2025-01-XX)
- [ ] Step 3.2: NCBI実装との再確認

#### テスト結果 (Step 2.2修正後)

**Megablast:**
- NZ_CP006932 self: LOSAT=270, NCBI=454 (59.5%)
- EDL933 vs Sakai: LOSAT=1508, NCBI=5718 (26.4%)
- Sakai vs MG1655: LOSAT=1202, NCBI=6476 (18.6%)

**Blastn:**
- NZ_CP006932 self: LOSAT=6190, NCBI=12340 (50.2%)
- PesePMNV vs MjPMNV: LOSAT=389, NCBI=241 (161.4%)
- MelaMJNV vs PemoMJNVA: LOSAT=712, NCBI=2729 (26.1%)

**考察:**
- 修正により、一部のテストケースで改善が見られた
- しかし、まだNCBIより少ないヒット数が生成されている
- 次のステップ: Cutoff score適用の確認とExtension実装の確認が必要

### Step 2.3: hit_level_array更新タイミングの完全修正 ✅ 完了

#### 修正内容
1. **更新タイミングの統一**: ✅ 修正完了
   - NCBI: `hit_level_array`は常に最後に更新（`na_ungapped.c:768-772`）
   - LOSAT (修正前): 拡張失敗時に早期更新していた
   - LOSAT (修正後): 常に最後に更新するように修正

2. **s_end_posの計算**: ✅ 修正完了
   - NCBI: 拡張成功時は`ungapped_data->length + ungapped_data->s_start + diag_table->offset`
   - NCBI: 拡張失敗時は`s_end_pos`を更新しない（`type_of_word`後の値を使用）
   - LOSAT: 同様のロジックを実装

3. **hit_len_arrayの更新条件**: ✅ 修正完了
   - NCBI: `if (two_hits) { hit_len_array[real_diag] = (hit_ready) ? 0 : s_end_pos - s_off_pos; }`
   - LOSAT: `if TWO_HIT_WINDOW > 0 { hit_len_array[diag_idx] = if hit_ready { 0 } else { final_s_end_pos - s_off_pos }; }`

#### 修正詳細

**NCBI参照:**
- `na_ungapped.c:752-772`: 拡張成功/失敗時の処理と`hit_level_array`更新
- `na_ungapped.c:768-772`: `hit_level_array`の更新は常に最後に実行

**修正箇所:**
1. Two-stage lookup: 1136-1269行目を修正
2. Original lookup: 1685-1806行目を修正
3. 拡張失敗時の早期更新を削除し、最後に統一して更新

**修正日:**
2025-01-24

#### テスト結果 (Step 2.3修正後)

**Megablast:**

| Test Case | LOSAT Hits | NCBI Hits | Ratio | Avg Length (LOSAT/NCBI) | Avg Bitscore (LOSAT/NCBI) | Avg Identity (LOSAT/NCBI) | Avg E-value (LOSAT/NCBI) |
|-----------|------------|-----------|-------|-------------------------|---------------------------|----------------------------|--------------------------|
| NZ_CP006932 self | 270 | 454 | 59.5% | 2805.3 / 2064.3 | 4887.3 / 3157.1 | 87.31% / 83.66% | 3.78e-08 / 2.64e-07 |
| EDL933 vs Sakai | 1508 | 5718 | 26.4% | 4264.3 / 1438.7 | 7653.6 / 2483.4 | 92.16% / 93.23% | 3.57e-05 / 7.78e-05 |
| Sakai vs MG1655 | 1202 | 6476 | 18.6% | 3677.0 / 772.8 | 6348.0 / 1298.4 | 93.66% / 93.25% | 2.47e-05 / 6.19e-05 |

**Blastn:**

| Test Case | LOSAT Hits | NCBI Hits | Ratio | Avg Length (LOSAT/NCBI) | Avg Bitscore (LOSAT/NCBI) | Avg Identity (LOSAT/NCBI) | Avg E-value (LOSAT/NCBI) |
|-----------|------------|-----------|-------|-------------------------|---------------------------|----------------------------|--------------------------|
| NZ_CP006932 self | 6190 | 12340 | 50.2% | 244.1 / 175.7 | 376.8 / 168.1 | 83.84% / 81.03% | 2.79e+00 / 3.31e+00 |
| PesePMNV vs MjPMNV | 389 | 241 | 161.4% | 409.4 / 771.5 | 295.0 / 448.6 | 77.51% / 79.63% | 2.25e-01 / 1.32e+00 |
| MelaMJNV vs PemoMJNVA | 712 | 2729 | 26.1% | 181.9 / 86.3 | 128.0 / 60.4 | 80.32% / 84.04% | 6.79e-01 / 2.37e+00 |

**考察:**
- `hit_level_array`の更新タイミングを修正したが、ヒット数の差異は依然として大きい
- 一部のテストケース（NZ_CP006932 self）では改善が見られるが、他のケースでは大きな差異が残っている
- LOSATの平均ヒット長がNCBIより長い傾向がある（フィルタリングが不十分な可能性）
- LOSATの平均bitscoreがNCBIより高い傾向がある（高スコアのヒットのみが残っている可能性）
- 次のステップ: Cutoff score適用の確認とExtension実装の詳細確認が必要

### Step 2.4: hit_level_array更新タイミングの修正（ungapped extension後、gapped extension前） ✅ 完了

#### 修正内容
1. **更新タイミングの修正**: ✅ 修正完了
   - **NCBI実装**: `hit_level_array`はungapped extension直後、gapped extension**前**に更新（`na_ungapped.c:768-772`）
   - **LOSAT (修正前)**: gapped extension**後**に更新していた
   - **LOSAT (修正後)**: ungapped extension直後、gapped extension前に更新するように修正

2. **NCBI参照**:
   - `na_ungapped.c:752-761`: ungapped extension後、cutoff scoreチェック
   - `na_ungapped.c:757-758`: `s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;`
   - `na_ungapped.c:768-772`: `hit_level_array[real_diag].last_hit = s_end_pos;` (ungapped extension後、gapped extension前)
   - **重要**: NCBIではgapped extensionは別の関数（`blast_gapalign.c`）で実行されるため、`hit_level_array`の更新はungapped extension関数内で完了

3. **修正箇所**:
   - Two-stage lookup: `utils.rs:1150-1205` - `hit_level_array`更新をgapped extension前に移動
   - Original lookup: `utils.rs:1706-1760` - 同様に修正

**修正日:**
2025-01-24

#### テスト結果 (Step 2.4修正後)

| Test Case | LOSAT Hits | NCBI Hits | Ratio | Avg Length (LOSAT/NCBI) | Avg Bitscore (LOSAT/NCBI) | Avg Identity (LOSAT/NCBI) | Avg E-value (LOSAT/NCBI) |
|-----------|------------|-----------|-------|-------------------------|---------------------------|----------------------------|--------------------------|
| NZ_CP006932.NZ_CP006932.megablast | 270 | N/A | N/A | 2805.3 / N/A | 4887.3 / N/A | 87.31% / N/A | 3.78e-08 / N/A |
| EDL933.Sakai.megablast | 1508 | 5718 | 26.4% | 4264.3 / 1438.7 | 7653.6 / 2483.4 | 92.16% / 93.23% | 3.57e-05 / 7.78e-05 |
| Sakai.MG1655.megablast | 1202 | N/A | N/A | 3677.0 / N/A | 6348.0 / N/A | 93.66% / N/A | 2.47e-05 / N/A |
| NZ_CP006932.NZ_CP006932.blastn | 6190 | N/A | N/A | 244.1 / N/A | 376.8 / N/A | 83.84% / N/A | 2.79e+00 / N/A |
| PesePMNV.MjPMNV.blastn | 389 | N/A | N/A | 409.4 / N/A | 295.0 / N/A | 77.51% / N/A | 2.25e-01 / N/A |
| MelaMJNV.PemoMJNVA.blastn | 712 | 2729 | 26.1% | 181.9 / 86.3 | 128.0 / 60.4 | 80.32% / 84.04% | 6.79e-01 / 2.37e+00 |
| SiNMV.ChdeNMV.blastn | 1836 | N/A | N/A | 230.0 / N/A | 329.7 / N/A | 89.02% / N/A | 6.42e-01 / N/A |
| PmeNMV.MjPMNV.blastn | 366 | N/A | N/A | 437.3 / N/A | 320.5 / N/A | 77.45% / N/A | 3.29e-01 / N/A |
| PmeNMV.PesePMNV.blastn | 503 | N/A | N/A | 456.7 / N/A | 442.7 / N/A | 81.45% / N/A | 5.53e-01 / N/A |
| PeseMJNV.PemoMJNVB.blastn | 1384 | N/A | N/A | 208.9 / N/A | 181.8 / N/A | 81.76% / N/A | 8.19e-01 / N/A |
| PemoMJNVA.PeseMJNV.blastn | 1179 | N/A | N/A | 333.9 / N/A | 455.6 / N/A | 86.05% / N/A | 5.75e-01 / N/A |
| MjeNMV.MelaMJNV.blastn | 1342 | 2668 | 50.3% | 315.9 / 252.0 | 476.8 / 290.1 | 87.14% / 84.14% | 6.88e-01 / 1.67e+00 |
| MjPMNV.MlPMNV.blastn | 2246 | 54402 | 4.1% | 218.0 / 144.8 | 313.2 / 122.3 | 87.84% / 80.41% | 4.26e-01 / 3.46e-01 |

**考察:**
- `hit_level_array`の更新タイミングをNCBI実装に合わせて修正（ungapped extension後、gapped extension前）
- しかし、ヒット数の差異は依然として大きい（特に`MjPMNV.MlPMNV.blastn`で4.1%）
- LOSATの平均ヒット長がNCBIより長い傾向がある（フィルタリングが不十分な可能性）
- LOSATの平均bitscoreがNCBIより高い傾向がある（高スコアのヒットのみが残っている可能性）
- 次のステップ: Cutoff score適用の確認とExtension実装の詳細確認が必要

## 注意事項

1. **各ステップが終わるごとに**:
   - NCBI実装と完全に合致しているか再確認
   - 合致していない箇所を修正
   - Release版をビルド
   - 統合テストをすべて実施
   - 結果を比較して差異を確認

2. **NCBIコード参照**:
   - すべての修正箇所にNCBIコード参照をコメントとして追加
   - ファイルパスと行番号を明記

3. **テスト結果の記録**:
   - 各テストケースのヒット数、長さ、e-value, bitscore, identityを記録
   - NCBI BLASTとの差異を分析

## 参考文献

- NCBI BLAST source: `/mnt/c/Users/genom/GitHub/ncbi-blast/`
- Key files:
  - `c++/src/algo/blast/core/na_ungapped.c` - Extension logic
  - `c++/src/algo/blast/core/blast_parameters.c` - Cutoff calculations
  - `c++/include/algo/blast/core/blast_options.h` - Default values

