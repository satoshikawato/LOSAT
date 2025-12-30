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

