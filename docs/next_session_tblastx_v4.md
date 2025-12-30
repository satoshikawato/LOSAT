# LOSAT TBLASTX 改善セッション引継ぎ（v4）

## 今回のセッションでの変更

### 1. MIN_UNGAPPED_SCORE の変更
- **変更前**: 22
- **変更後**: 14
- **理由**: より多くの低スコアヒットを捕捉するため

### 2. X_DROP_UNGAPPED の変更
- **変更前**: 11
- **変更後**: 7 (NCBIデフォルト)
- **理由**: NCBI BLASTとの互換性向上

### 3. sum_stats_linking.rs の改善
- small_gap_sum_e と large_gap_sum_e の両方を使用
- 最も良い（低い）E-valueを選択
- window_size = 50 (NCBI default: gap_size + overlap_size + 1)
- 効率改善（処理時間短縮）

## 現在の状態

### パフォーマンス
| 比較タイプ | LOSAT | NCBI | 座標共通率 |
|-----------|-------|------|----------|
| Self (NZ_CP006932) | 100,091 | 62,053 | **7.6%** |
| Cross (AP027078 vs AP027131) | 6,029 | 30,175 | **1.9%** |

### 根本的な問題

**座標の共通率が非常に低い**（2-8%）
- NCBIとLOSATが異なる座標でヒットを生成
- これはsum-statisticsリンキングの問題ではない
- extensionロジックまたはマスキングロジックの違いが原因

### 発見された違い

1. **ヒットパターンの違い**:
   - NCBI: 同じ領域で複数の短いヒット（16-24 aa）を出力
   - LOSAT: 同じ領域で1つの長いヒット（33-115 aa）を出力

2. **E-value分布の違い**:
   - NCBI: 9,789ヒットがE-value=0.0（sum-statisticsリンキングの効果）
   - LOSAT: 2-4ヒットがE-value=0.0

3. **マスキングの影響**:
   - LOSATは1つのヒットを出力後、その領域をマスクして追加ヒットを抑制
   - NCBIは同じ領域で複数のヒットを出力

---

## 次のステップ（優先度順）

### 優先度1: マスキングロジックの調査
NCBIは同じ領域で複数のヒットを出力しているが、LOSATは1つの長いヒットのみ。
- `utils.rs` のマスク更新ロジックを見直す
- `aa_ungapped.c` のdiag_array更新ロジックと比較

調査ポイント:
```
// LOSAT (utils.rs)
mask.insert(mask_key, se_ungapped.max(s_last_off.saturating_sub(k_size - 1)));

// NCBI (aa_ungapped.c)
diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
```

### 優先度2: Extensionロジックの比較
NCBIのextension結果（短いヒット）とLOSATの結果（長いヒット）が異なる理由を調査。
- `extend_hit_ungapped` と `extend_hit_two_hit` の実装を `aa_ungapped.c` と比較
- X-drop終了条件の詳細確認

### 優先度3: フレーム処理の確認
TBLASTXでは36フレームの組み合わせがある。フレーム処理が正しいか確認。

---

## テストコマンド

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT

# ビルド
cargo build --release

# Cross比較テスト
./target/release/losat tblastx \
  --query tests/fasta/AP027078.fasta \
  --subject tests/fasta/AP027131.fasta \
  --query-gencode 4 --db-gencode 4 \
  --num-threads 1 --out /tmp/cross_test.out

# Self比較テスト
./target/release/losat tblastx \
  --query tests/fasta/NZ_CP006932.fasta \
  --subject tests/fasta/NZ_CP006932.fasta \
  --query-gencode 4 --db-gencode 4 \
  --num-threads 1 --out /tmp/self_test.out

# 座標比較
awk -F'\t' 'NR>0 {print $7"-"$8" "$9"-"$10}' /tmp/cross_test.out | sort -u > /tmp/losat_coords.txt
grep -v "^#" tests/blast_out/AP027078.AP027131.tblastx.n1.out | awk -F'\t' '{print $7"-"$8" "$9"-"$10}' | sort -u > /tmp/ncbi_coords.txt
comm -12 /tmp/losat_coords.txt /tmp/ncbi_coords.txt | wc -l  # Common
```

---

## 重要なファイル

### メインソースファイル
- `LOSAT/src/algorithm/tblastx/utils.rs` - メイン処理（マスキングロジック含む）
- `LOSAT/src/algorithm/tblastx/extension.rs` - ungapped extension
- `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` - sum-statisticsリンキング
- `LOSAT/src/algorithm/tblastx/constants.rs` - 定数定義

### NCBI参照ファイル
- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c` - extension処理
- `ncbi-blast/c++/src/algo/blast/core/link_hsps.c` - HSPリンキング

---

## 注意事項

1. **座標共通率が改善目標**: ヒット数よりも座標の一致率が重要
2. **マスキングが鍵**: 同じ領域で複数ヒットを出力するかどうかが大きな違い
3. **git状態**: `MIN_UNGAPPED_SCORE`, `X_DROP_UNGAPPED`, `sum_stats_linking.rs` が変更済み

