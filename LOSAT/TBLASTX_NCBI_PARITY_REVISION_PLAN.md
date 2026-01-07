# TBLASTX NCBI実装完全一致改修計画

**作成日**: 2026-01-06  
**目的**: LOSAT TBLASTX実装をNCBI BLAST+実装に完全に一致させるための網羅的な改修計画  
**NCBI codebase**: `/mnt/c/Users/genom/GitHub/ncbi-blast/`  
**原則**: アルゴリズムの実行順序、タイミング、コンテキストをNCBI実装と完全に一致させる

---

## 1. 改修項目の分類

### 1.1 優先度分類

- **🔴 最優先**: 出力結果に直接影響する不一致
- **🟡 高優先度**: アルゴリズムの実行順序・タイミングの不一致
- **🟢 中優先度**: 構造的な差異（出力には影響しないが、NCBIとの構造的一致のため）
- **⚪ 低優先度**: 最適化・リファクタリング（出力に影響しない）

---

## 2. 最優先改修項目

### 2.1 リンキング処理のフレームグループ並列化

**状態**: 🔴 未実装（並列化可能だが逐次処理）

**NCBI実装**: `link_hsps.c:553`
```c
for (frame_index=0; frame_index<num_query_frames; frame_index++)
{
    // Process each frame group independently
    // No shared state between groups
}
```

**LOSAT実装**: `sum_stats_linking.rs:682-699`
```rust
// Step 3: Process each frame group SEQUENTIALLY (NCBI does not parallelize)
let mut results: Vec<UngappedHit> = Vec::new();
for group_hits in frame_groups {
    let processed = link_hsp_group_ncbi(...);
    results.extend(processed);
}
```

**改修内容**:
- `for` ループを `par_iter()` に変更
- NCBIコード引用をコメントとして追加
- 出力順序は最終ソートで保証されるため問題なし

**詳細**: `TBLASTX_PARALLELIZATION_PLAN.md` を参照

---

### 2.2 HSP生成数の過剰問題（長い配列）

**状態**: 🔴 未解決

**現象**: 600kb+の配列でLOSATが338,859 HSPsを生成 vs NCBIの30,000-45,000（約7倍）

**根本原因候補**:
1. **Seeding/Extension段階での過剰生成**
   - Cutoff score計算の不一致
   - X-drop termination条件の不一致
   - Two-hit window処理の不一致

2. **Reevaluation段階でのフィルタリング不足**
   - NCBI: `Blast_HSPReevaluateWithAmbiguitiesUngapped` で厳密にフィルタ
   - LOSAT: 同等の実装だが、パラメータや条件が異なる可能性

**調査が必要な項目**:
- [ ] Cutoff score計算の詳細比較（`ncbi_cutoffs.rs` vs `blast_parameters.c`）
- [ ] Extension X-drop条件の詳細比較（`extension.rs` vs `aa_ungapped.c`）
- [ ] Two-hit window処理の詳細比較（`utils.rs` vs `blast_extend.c`）
- [ ] Reevaluation条件の詳細比較（`reevaluate.rs` vs `blast_hits.c`）

**参考**: `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md`

---

## 3. 高優先度改修項目

### 3.1 処理フローの順序とタイミング

#### 3.1.1 Reevaluationのタイミング

**NCBI実装**: `blast_engine.c:1492-1497`
```c
// Step 1: Extension → Save (absolute coordinates)
BlastSaveInitHsp(ungapped_hsps, ...);

// Step 2: Convert to context-relative coordinates
BLAST_GetUngappedHSPList(...);

// Step 3: Batch reevaluation
Blast_HSPListReevaluateUngapped(hsp_list, ...);
```

**LOSAT実装**: `utils.rs:2070-2093`
```rust
// Step 1: Extension → Save (absolute coordinates)
let ungapped_hits = get_ungapped_hsp_list(...);

// Step 2: Batch reevaluation
let ungapped_hits = reevaluate_ungapped_hsp_list(...);
```

**状態**: ✅ 実装済み（順序は一致）

**確認事項**:
- [ ] 座標変換のタイミングがNCBIと一致しているか
- [ ] Reevaluationの条件がNCBIと完全に一致しているか

**参考**: `TBLASTX_SEEDING_EXTENSION_DIFFERENCES.md`

---

#### 3.1.2 座標変換の実装

**NCBI実装**: `blast_gapalign.c:4719-4775`
- 絶対座標（concatenated buffer）→ context-relative座標への変換
- `s_AdjustInitialHSPOffsets` 関数を使用

**LOSAT実装**: `utils.rs:620-682`
- `get_ungapped_hsp_list()` で座標変換を実装

**状態**: ✅ 実装済み

**確認事項**:
- [ ] 座標変換の計算式がNCBIと完全に一致しているか
- [ ] Frame base offsetの計算がNCBIと一致しているか

---

### 3.2 アルゴリズムパラメータの一致

#### 3.2.1 Cutoff Score計算

**状態**: ✅ 実装済み（`ncbi_cutoffs.rs`）

**確認済み項目**:
- ✅ Per-subject cutoff更新: `cutoff_score_for_update_tblastx()`
- ✅ Cutoff score max: `cutoff_score_max_for_tblastx()`
- ✅ Gap trigger計算: `gap_trigger_raw_score()`
- ✅ X-drop計算: `x_drop_raw_score()`

**確認事項**:
- [ ] 計算式がNCBIと完全に一致しているか（浮動小数点精度含む）
- [ ] エッジケース（非常に長い配列、極端なE-value）での動作が一致しているか

---

#### 3.2.2 Sum-Statistics Linkingパラメータ

**状態**: ✅ 実装済み（`sum_stats_linking.rs`）

**確認済み項目**:
- ✅ Gap size: `GAP_SIZE = 40`
- ✅ Overlap size: `OVERLAP_SIZE = 9`
- ✅ Window size: `WINDOW_SIZE = 50`
- ✅ Trim size: `TRIM_SIZE = 5`
- ✅ Gap decay rate: `BLAST_GAP_DECAY_RATE = 0.5`

**確認事項**:
- [ ] すべての定数値がNCBIと一致しているか
- [ ] 計算式がNCBIと完全に一致しているか

---

## 4. 中優先度改修項目

### 4.1 構造的な差異（出力には影響しない）

#### 4.1.1 フレームグループの処理順序

**NCBI実装**: フレームグループを逐次処理

**LOSAT実装**: フレームグループを並列処理可能（改修予定）

**状態**: 🟡 改修予定（`TBLASTX_PARALLELIZATION_PLAN.md`参照）

**影響**: 出力順序が変わる可能性があるが、最終的にNCB順序でソートされるため問題なし

---

#### 4.1.2 Subject処理の並列化

**NCBI実装**: Subjectを逐次処理

**LOSAT実装**: Subjectを並列処理（`utils.rs:1612`）

**状態**: ✅ 実装済み

**影響**: 出力順序が変わる可能性があるが、最終的にNCB順序でソートされるため問題なし

**確認事項**:
- [ ] 並列化が出力結果に影響していないか
- [ ] 最終ソートが正しく実装されているか

---

### 4.2 コード構造の一致

#### 4.2.1 関数名と構造の対応

**状態**: ✅ 概ね一致

**確認事項**:
- [ ] 主要関数の名前がNCBIと対応しているか
- [ ] 関数の分割がNCBIと類似しているか

---

## 5. 低優先度改修項目

### 5.1 デバッグ出力の整理

**状態**: ⚪ 多数のデバッグ出力が残っている

**改修内容**:
- デバッグ出力を環境変数で制御可能にする
- 本番ビルドではデバッグ出力を無効化

---

### 5.2 コメントの充実

**状態**: ⚪ NCBIコード引用が一部不足

**改修内容**:
- すべての移植関数にNCBIコード引用を追加
- ファイルパスと行番号を明記

---

## 6. NCBI実装との詳細比較

### 6.1 処理フローの比較

#### NCBI BLAST+ TBLASTX処理フロー

```
1. Query読み込み・翻訳
   ↓
2. Lookup table構築
   ↓
3. Subject読み込み・翻訳（逐次）
   ↓
4. Per-subject処理（逐次）:
   a. Cutoff score計算
   b. Seeding (two-hit window)
   c. Extension
   d. HSP保存（絶対座標）
   ↓
5. 座標変換（絶対→context-relative）
   ↓
6. Reevaluation（一括）
   ↓
7. Sum-statistics linking（フレームグループ逐次）
   ↓
8. HSP culling（オプション）
   ↓
9. 出力（NCB順序でソート）
```

#### LOSAT TBLASTX処理フロー

```
1. Query読み込み・翻訳
   ↓
2. Lookup table構築
   ↓
3. Subject読み込み・翻訳（並列化可能）
   ↓
4. Per-subject処理（並列化）:
   a. Cutoff score計算
   b. Seeding (two-hit window)
   c. Extension
   d. HSP保存（絶対座標）
   ↓
5. 座標変換（絶対→context-relative）
   ↓
6. Reevaluation（一括）
   ↓
7. Sum-statistics linking（フレームグループ逐次→並列化予定）
   ↓
8. HSP culling（オプション）
   ↓
9. 出力（NCB順序でソート）
```

**差異**:
- Subject処理: NCBIは逐次、LOSATは並列化
- フレームグループ処理: NCBIは逐次、LOSATは並列化予定
- **出力結果への影響**: なし（最終ソートで保証）

---

### 6.2 アルゴリズムの詳細比較

#### 6.2.1 Two-Hit Window処理

**NCBI実装**: `blast_extend.c:103-200`
- `last_hit = -window` で初期化
- `diff >= window` でtwo-hit判定

**LOSAT実装**: `utils.rs:51-200`
- `last_hit = 0` で初期化
- `diff >= window` でtwo-hit判定

**状態**: ✅ 動作は同等（初期化値の違いは出力に影響しない）

**確認事項**:
- [ ] 動作が完全に一致しているか

---

#### 6.2.2 Extension処理

**NCBI実装**: `aa_ungapped.c:575-591`
```c
score = s_BlastAaExtendTwoHit(matrix, subject, query,
                              last_hit + wordsize,
                              subject_offset, query_offset,
                              cutoffs->x_dropoff, 
                              &hsp_q, &hsp_s,
                              &hsp_len, use_pssm,
                              wordsize, &right_extend,
                              &s_last_off);

if (score >= cutoffs->cutoff_score)
    BlastSaveInitHsp(...);
```

**LOSAT実装**: `extension.rs` + `utils.rs`
- `extend_hit_two_hit()` 関数で実装
- Cutoff check後に保存

**状態**: ✅ 実装済み

**確認事項**:
- [ ] X-drop termination条件が完全に一致しているか
- [ ] Score計算が完全に一致しているか

---

#### 6.2.3 Sum-Statistics Linking処理

**NCBI実装**: `link_hsps.c:553-1013`
- フレームグループを逐次処理
- `lh_helper` 配列を使用
- `changed` フラグで最適化

**LOSAT実装**: `sum_stats_linking.rs:719-1945`
- フレームグループを逐次処理（並列化予定）
- `lh_helpers` 配列を使用
- `changed` フラグで最適化

**状態**: ✅ 実装済み（並列化は改修予定）

**確認事項**:
- [ ] すべての最適化が実装されているか
- [ ] チェーン形成ロジックが完全に一致しているか

---

## 7. 実装チェックリスト

### 7.1 並列化関連

- [ ] リンキング処理のフレームグループ並列化（`sum_stats_linking.rs:682-699`）
- [ ] NCBIコード引用の追加
- [ ] 並列化前後で出力が一致することを確認
- [ ] 性能評価を実施

### 7.2 HSP生成数問題の調査

- [ ] Cutoff score計算の詳細比較
- [ ] Extension X-drop条件の詳細比較
- [ ] Two-hit window処理の詳細比較
- [ ] Reevaluation条件の詳細比較
- [ ] 根本原因の特定と修正

### 7.3 アルゴリズムの検証

- [ ] すべての計算式がNCBIと一致しているか
- [ ] エッジケースでの動作が一致しているか
- [ ] 浮動小数点精度がNCBIと一致しているか

### 7.4 コード品質

- [ ] すべての移植関数にNCBIコード引用を追加
- [ ] デバッグ出力の整理
- [ ] コメントの充実

---

## 8. 検証計画

### 8.1 機能検証

1. **短い配列（300kb）**
   - 入力: AP027280
   - 期待: NCBIと一致（42,733 hits）
   - 現在: LOSAT 42,797 vs NCBI 42,733（差: +64, 0.15%）

2. **長い配列（600kb+）**
   - 入力: 長い配列
   - 期待: NCBIと一致（14,871 hits）
   - 現在: LOSAT 29,766 vs NCBI 14,871（約2倍）

3. **複数クエリ**
   - 入力: 複数のクエリ配列
   - 期待: 各クエリの結果がNCBIと一致

### 8.2 性能検証

1. **並列化効果**
   - 並列化前後の処理時間を比較
   - スレッド数による性能変化を測定

2. **メモリ使用量**
   - 並列化によるメモリ使用量の変化を確認

---

## 9. 参考資料

### NCBI実装ファイル

- `c++/src/algo/blast/core/link_hsps.c` - Sum-statistics linking
- `c++/src/algo/blast/core/aa_ungapped.c` - Extension logic
- `c++/src/algo/blast/core/blast_parameters.c` - Cutoff calculations
- `c++/src/algo/blast/core/blast_extend.c` - Two-hit window
- `c++/src/algo/blast/core/blast_hits.c` - Reevaluation
- `c++/src/algo/blast/core/blast_engine.c` - Main processing flow

### LOSAT実装ファイル

- `src/algorithm/tblastx/sum_stats_linking.rs` - Sum-statistics linking
- `src/algorithm/tblastx/extension.rs` - Extension logic
- `src/algorithm/tblastx/ncbi_cutoffs.rs` - Cutoff calculations
- `src/algorithm/tblastx/utils.rs` - Main processing flow

### 関連ドキュメント

- `TBLASTX_PARALLELIZATION_PLAN.md` - 並列化計画
- `TBLASTX_NCBI_PARITY_STATUS.md` - パリティ状態レポート
- `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md` - 長い配列でのHSP過剰生成問題
- `TBLASTX_SEEDING_EXTENSION_DIFFERENCES.md` - Seeding/Extension段階での違い

---

## 10. 結論

LOSAT TBLASTX実装をNCBI BLAST+実装に完全に一致させるためには、以下の改修が必要です：

1. **最優先**: リンキング処理のフレームグループ並列化（出力結果に影響なし）
2. **最優先**: HSP生成数の過剰問題の調査と修正（長い配列での不一致の根本原因）
3. **高優先度**: アルゴリズムパラメータの詳細検証
4. **中優先度**: 構造的な差異の確認（出力には影響しない）
5. **低優先度**: コード品質の向上

並列化は比較的簡単に実装可能で、性能向上が期待できます。HSP生成数の過剰問題は、詳細な調査が必要ですが、根本原因を特定できれば修正可能です。

すべての改修は、**NCBI実装を唯一の正解（GROUND TRUTH）として**、アルゴリズムの実行順序、タイミング、コンテキストを完全に一致させることを目標とします。

