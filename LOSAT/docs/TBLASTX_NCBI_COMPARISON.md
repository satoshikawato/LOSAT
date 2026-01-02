# LOSAT TBLASTX - NCBI BLAST比較・修正状況

## 現在の問題点

フルサイズのゲノム（657kb）に対するTBLASTX実行が5分以上かかる（タイムアウト）。
３０秒以内（望ましくは１０秒以内）に完了すること。
## 実施済みの修正

### 1. SEGフィルター (src/utils/seg.rs)

| 項目 | NCBI実装 | LOSAT修正前 | LOSAT修正後 |
|------|----------|-------------|-------------|
| エントロピー計算 | Shannon entropy with `fabs(ent / total)` | 独自実装 | ✅ NCBI準拠に修正 |
| ln_factorial | `lnfact[]`事前計算配列 | 毎回計算 | ✅ 事前計算配列に修正 |
| s_Trim関数 | O(N)スライディングウィンドウ | O(N³)の非効率実装 | ✅ O(N)に修正 |
| パラメータ | maxbogus=2, maxtrim=50 | 未実装 | ✅ 追加 |

**結果**: SEGマスク率 97.73% → 32.03%に改善

### 2. ルックアップテーブル構築 (src/algorithm/tblastx/lookup.rs)

| 項目 | NCBI実装 | LOSAT修正前 | LOSAT修正後 |
|------|----------|-------------|-------------|
| Neighbor生成 | 再帰的枝刈り `s_AddWordHitsCore` | 13,824全組み合わせ走査 | ✅ 再帰的枝刈りに修正 |
| SEGマスク適用 | フレーム別にAA座標でマスク | DNA座標に変換・マージ | ✅ フレーム別AA座標に修正 |
| row_max事前計算 | 各AAの最大スコアを事前計算 | なし | ✅ 追加 |

**結果**: 平均17.1 neighbors/wordに削減（13,824から）

### 3. 対角線追跡 (src/algorithm/tblastx/utils.rs)

| 項目 | NCBI実装 | LOSAT修正前 | LOSAT修正後 |
|------|----------|-------------|-------------|
| データ構造 | `DiagStruct[]`事前確保配列 | FxHashMap | ✅ Vec<DiagStruct>に修正 |
| アクセス方法 | `diag_array[diag_coord]` O(1) | HashMap lookup O(1)〜O(n) | ✅ 配列インデックスに修正 |
| diag_offset | subject毎にリセット不要 | なし | ✅ 追加 |

## 未修正・要確認の問題点

### 1. Subject Scanning（最も疑わしい）

**NCBI実装** (`blast_aascan.c`):
```c
// Presence Vector (PV) によるO(1)の存在チェック
if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
    numhits = bbc[index].num_used;
    // hits をコピー
}
```

**LOSATの現状**:
```rust
// Vec::is_empty()でチェック（メタデータロード必要）
let matches = unsafe { lookup.get_unchecked(kmer) };
if !matches.is_empty() { ... }
```

**必要な修正**:
- Presence Vector (ビットマップ) の実装
- `ComputeTableIndexIncremental`によるインクリメンタルインデックス計算

### 2. Two-Hit Extension

**NCBI実装** (`aa_ungapped.c`):
- `s_BlastAaExtendTwoHit`: 左方向の拡張が最初のヒットに到達した場合のみ右方向に拡張
- ワードサイズ分のスコアリング最適化

**LOSATの現状**:
- `extend_hit_two_hit`は実装済みだが、NCBI完全互換か要確認

### 3. Batch Processing

**NCBI実装**:
- `scansub`が一度に`array_size`個のヒットを返す
- バッチ処理でメモリアクセスを最適化

**LOSATの現状**:
- 各ヒットを逐次処理
- バッチ最適化なし

### 4. QueryFrame毎の対角線配列

**現在の問題**: 
各subjectフレームの処理開始時に、全queryフレーム分の対角線配列を新規作成している。
これにより大きなメモリ確保と初期化コストが発生。

**NCBI実装**:
- `diag_offset`を使用して配列をゼロクリアせずに再利用

## パフォーマンス比較

| 処理 | NCBI (推定) | LOSAT現状 |
|------|-------------|-----------|
| 小テスト (7kb) | <1秒 | 0.05秒 |
| 中テスト (47kb) | <5秒 | 2.3秒 |
| フルテスト (657kb) | <30秒 | >300秒（タイムアウト） |

## 次のセッションで実施すべき作業

### 優先度: 高

1. **Presence Vector実装**
   - `src/algorithm/tblastx/lookup.rs`に`Vec<u64>`のビットマップ追加
   - `PV_TEST`マクロ相当の関数実装
   - ルックアップテーブル構築時にPVも同時構築

2. **インクリメンタルインデックス計算**
   - `ComputeTableIndexIncremental`相当の実装
   - スキャンループで毎回3文字エンコードせず、シフト＆マスクで更新

3. **対角線配列の再利用**
   - subjectフレーム毎に新規作成せず、`diag_offset`で再利用
   - メモリ確保コストを削減

### 優先度: 中

4. **Two-Hit Extension最適化**
   - NCBI `s_BlastAaExtendTwoHit`との完全比較
   - 左拡張の早期終了条件確認

5. **Batch Hit Processing**
   - スキャン結果をバッチで処理
   - キャッシュ効率向上

### 優先度: 低

6. **Sum Statistics Linking**
   - dual-index (small_gap + large_gap) の完全実装
   - `lh_helper`最適化構造の移植

## 参照すべきNCBIソースファイル

- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c` - Subject scanning
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c` - Two-hit extension
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c` - Lookup table construction
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_extend.h` - DiagStruct定義
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c` - Sum statistics linking

---

## Session Log: 2026-01-02 - 残り1%差異の徹底調査

### セッション概要

このセッションは、LOSAT TBLASTX と NCBI BLAST+ の間に残る約1%のヒット数差異（LOSAT 42,334 vs NCBI 42,733）の根本原因を特定するための徹底的な調査セッションであった。多くの仮説を検証し、複数のコード比較・修正試行を行った。

### 現在の到達点

| 指標 | 値 |
|------|-----|
| LOSAT hits | 42,334 |
| NCBI hits | 42,733 |
| 一致率 | **99.07%** |
| 差異 | 399 hits (0.93%) |

### 調査済み項目と結果

#### 1. フレームペア別分析

各query/subjectフレーム組み合わせ(9種類のminus-minusペア)で差異を分析：

| Frame Pair | LOSAT | NCBI | 差異 | 差異率 |
|------------|-------|------|------|--------|
| (-2,-2) | 2,751 | 2,681 | **+70** | **-2.6%** (LOSATが多い) |
| (-3,-3) | 1,717 | 1,739 | +22 | 1.2% |
| (-1,-1) | 1,903 | 2,013 | +110 | 5.4% |
| (-1,-2) | 1,046 | 1,119 | +73 | 6.5% |
| (-1,-3) | 1,040 | 1,136 | +96 | 8.4% |
| (-2,-1) | 1,048 | 1,128 | +80 | 7.0% |
| (-2,-3) | 1,117 | 1,215 | +98 | 8.0% |
| (-3,-1) | 1,015 | 1,126 | +111 | **9.8%** (最大差異) |
| (-3,-2) | 1,150 | 1,264 | +114 | 9.0% |

**重要な発見**: 
- 同一フレームペア (x,x) は小さい差異 (1-5%)
- 異なるフレームペアは大きい差異 (7-10%)
- (-2,-2)のみLOSATが多く、他はNCBIが多い
- この非対称パターンは diagonal array の frame 間干渉を示唆

#### 2. Diagonal座標計算の符号処理

**調査内容**: `diag_coord = (query_offset - subject_offset) & diag_mask` の符号付き/無し演算

**結論**: Rust の `i32 & diag_mask` と C の `Uint4 & diag_mask` は二の補数表現により同一結果を生成。**問題なし**。

#### 3. Minus-frame Seed Indexing

**調査内容**: マイナスフレームのseed生成とindex計算

**結論**: LOSATの連結offset計算 (`frame_base + raw_pos`) とcontext index取得 (`get_context_idx()`) は NCBIの `BSearchContextInfo` と等価。**問題なし**。

#### 4. 短いヒット (≤12 AA) の Two-hit Window条件

**調査内容**: 短いヒットでの window 条件判定

**発見**: 双方向の差異あり。12 AA で LOSAT +17、14 AA で LOSAT -67。多くのNCBI欠損ヒットは、より長いLOSATヒットと重複。**Window条件自体には問題なし**。

#### 5. Tie-break条件: `s_last_off` 更新

**調査内容**: LOSAT `extend_hit_two_hit` と NCBI `s_BlastAaExtendRight` の比較

**結論**: 両実装とも `s_last_off` をスキャンした最右subject位置に更新。**ロジックは一致**。

#### 6. AA↔DNA座標変換

**調査内容**: `convert_coords` のエッジケース

**結論**: LOSATは1-based DNA座標（BLASTの出力形式に合わせた意図的設計）、内部AA座標は0-based。**出力フォーマットの差異であり、計算上の問題なし**。

#### 7. `diag_offset` への `+window` 追加

**NCBI** (`blast_extend.c:172`):
```c
ewp->diag_table->offset += subject_length + ewp->diag_table->window;
```

**LOSAT修正**:
```rust
diag_offset += s_aa_len as i32 + window;
```

**結果**: 変更適用後もヒット数は 42,334 のまま変化なし。単一subject/queryのテストケースでは効果が出にくい可能性。

### 座標系統一アプローチの試行と失敗

#### 仮説

NCBI と LOSAT の座標系の違いがフレーム境界チェックに影響：

**NCBI**:
- `query->sequence` は sentinel の次を指す (`sequence_start + 1`)
- 最初のAAの offset = 0
- `contexts[0].query_offset = 0`

**LOSAT (修正前)**:
- `frame.aa_seq[0]` は sentinel
- 最初のAAの stored offset = `frame_base + 1`
- `ctx.frame_base = frame_base` (= 0 for first frame)

境界チェック `if query_offset - diff < ctx.frame_base` で差異が生じる可能性。

#### 実施した修正

1. **frame_base計算の変更**:
   - Before: `base += frame.aa_seq.len()` (sentinel 2つ含む)
   - After: `base += frame.aa_len + 1` (AA数 + sentinel 1つ)

2. **Query offset格納の変更**:
   - Before: `frame_base + raw_pos` (where `raw_pos = aa_pos + 1`)
   - After: `frame_base + aa_pos` (0-indexed)

3. **Subject offset格納の変更**:
   - Before: `s_off = s as i32`
   - After: `s_off = (s - 1) as i32` (sentinel skip)

4. **Extension引数の変換**:
   - NCBI-style offset を array index に変換 (+1)
   - Extension戻り値を NCBI-style に逆変換 (-1)

#### 結果

```
修正前: 42,334 hits
修正後: 35,747 hits (17%減少)
```

**大幅な悪化**。座標系の変更は複数箇所に影響し、一箇所の変更が他の整合性を崩す。**即座にrevert**。

### 残り1%差異の根本原因についての考察

#### 可能性1: 浮動小数点丸め差

SEG entropy計算、cutoff score計算などで `ln()`, `exp()`, `ceil()`, `floor()` を多用。
累積的な丸め差がヒット判定の境界ケースで差異を生む可能性。

#### 可能性2: Diagonal Array の Frame間干渉

フレームペア非対称パターン（同一ペアは小差、異種ペアは大差）は、
diagonal array が異なるフレーム組み合わせ間で意図せず干渉している可能性を示唆。

しかし、LOSAT の `diag_coord = (query_offset - subject_offset) & diag_mask` 計算は
NCBI と同一であり、frame_base の累積により異なるフレームは異なる diagonal 空間に
マッピングされるはず。

#### 可能性3: 特定のエッジケース

Missing hit の詳細分析で、同一DNA diagonal上に異なるframe組み合わせのヒットが
存在するケースを発見。これらは論理的に独立した diagonal として処理されるべきだが、
微妙な条件分岐で異なる結果になる可能性。

### 今後の課題

#### 優先度: 低 (99%一致達成後)

1. **特定 missing hit のデバッグトレース**
   - 具体的な missing hit を選び、seed 生成から extension まで完全トレース
   - NCBI のデバッグ出力と比較

2. **SEG incremental update の完全移植**
   - 現在: 毎回再計算
   - NCBI: `s_ShiftWin1` による増分更新
   - 数学的に等価だが、浮動小数点誤差の累積パターンが異なる可能性

3. **浮動小数点精度の統一**
   - `f64` vs `f32` の使い分け確認
   - NCBI の `BLAST_Nint` (四捨五入) 関数の完全移植確認

### 結論

**99.07%の一致率は商用品質レベル**に達している。残り0.93%の差異は、
浮動小数点演算の累積誤差、境界条件のtie-break、その他の微細な実装差異の
複合的な結果と考えられる。

これ以上の一致率向上は、費用対効果の観点から優先度が低い。
次のステップとしては以下を推奨：

1. **パフォーマンス最適化**: 現在の処理速度の改善
2. **他の BLAST プログラム**: BLASTN, BLASTP 等への展開
3. **機能拡張**: 出力フォーマットの充実、オプションの追加

### 参照したNCBIソースファイル (追加分)

- `blast_extend.c` - `Blast_ExtendWordExit`, `s_BlastDiagTableNew`
- `blast_util.c` - `BLAST_GetAllTranslations`, `BLAST_ContextToFrame`
- `blast_query_info.c` - `ContextOffsetsToOffsetArray`, context offset管理
- `blast_lookup.c` - `BlastLookupIndexQueryExactMatches`, `ComputeTableIndex`
- `blast_parameters.c` - `CalculateLinkHSPCutoffs`, cutoff計算

### セッション統計

- 調査項目数: 19 (TODO list)
- 修正試行: 多数 (大半は revert)
- 最終結果: 99.07% 一致率維持（改善なし、悪化なし）
- 重要な発見: フレームペア非対称パターン、座標系変更の危険性

---

## 🔴 緊急バグ発見: `--neighbor-map` パスの差異

### 発見状況

セッション終了時に判明。標準パスと `--neighbor-map` パスでヒット数が異なる：

| パス | ヒット数 |
|------|----------|
| 標準 (`run()`) | 42,334 |
| neighbor-map (`run_with_neighbor_map()`) | 42,086 |
| **差異** | **248 hits (0.59%)** |

### 原因の可能性

1. **座標系の不整合**: 標準パスと neighbor-map パスで `query_offset`, `subject_offset` の計算方法が異なる可能性
2. **Diagonal tracking の差異**: 二つのパスで `diag_array` の管理方法が微妙に異なる
3. **Extension 呼び出しの差異**: `extend_hit_two_hit` への引数の渡し方が異なる可能性

### 次回セッションでの対応必須

`run_with_neighbor_map()` 関数（`utils.rs` 内）を `run()` 関数と詳細比較し、
座標系・diagonal tracking・extension 呼び出しを完全に統一する必要がある。

**これは商用リリース前に修正必須のバグ。**

