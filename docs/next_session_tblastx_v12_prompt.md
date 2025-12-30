# TBLASTXセッション v12 - cutoff_score問題の修正

## 発見された根本原因

**LOSATの`cutoff_score`が高すぎて、NCBIが出力する低bit scoreのヒットを全て落としている。**

### 証拠

| 項目 | NCBI | LOSAT |
|------|------|-------|
| 出力ヒット数 | 62,053 | 9,077 |
| 最小bit score | 22.1 | 32.2 |
| bit score 22〜32のヒット数 | 約31,000件 | 0件 |

NCBIは22.1〜32.2のbit scoreで31,000件以上を出力しているが、LOSATは全てフィルタリングしている。

### 差分の内訳
- NCBI 62,053件 - LOSAT 9,077件 = 約53,000件の差
- うち約31,000件はbit score 22〜32の範囲（cutoff_scoreで落とされた）

## 今までやったこと（v11セッション）

### 1. sum_stats_linking最適化（完了・効果あり）
- `changed`フラグによる差分更新
- `linked_to`カウンター
- `path_changed`フラグ
- `use_current_max`最適化
- `next_larger`による効率的なスキップ
- **結果**: 処理時間 2分超 → 11秒に短縮

### 2. cutoff_scoreにgap_decay_rate適用（実装済み・効果不十分）
- `raw_score_from_evalue_with_decay`関数を追加
- `dodecay=true`、`gap_decay_rate=0.5`を適用
- **結果**: ヒット数変わらず

### 3. cutoff_big_gapチェック追加（実装済み）
- `calculate_cutoff_big_gap`関数を追加
- linking内部でのスコアチェックに使用
- **結果**: パフォーマンス最適化のみ、ヒット数変わらず

### 4. E-value計算の修正（実装済み）
- `small_gap_sum_e`と`large_gap_sum_e`の両方を計算
- `gap_prob`による調整を追加
- より小さいE-valueを選択
- **結果**: ヒット数変わらず

### 5. index=0（small gap）の実装（試行→リバート）
- NCBIのindex=0ループを実装
- ユーザーがリバート

## 問題箇所の特定

`LOSAT/src/algorithm/tblastx/utils.rs` の817行付近：

```rust
if ungapped_score < cutoff_score {
    // HSPが落とされる
}
```

ここで`cutoff_score`が高すぎて、bit score 22〜32のヒットが全て落とされている。

## 次のセッションでやるべきこと

### 1. cutoff_score計算のデバッグ
現在の`cutoff_score`の値を出力し、NCBIと比較する。

```rust
// utils.rs 715行付近
let cutoff_score = *cutoff_score_cache.entry(context_key).or_insert_with(|| {
    // この計算が正しいか確認
});
```

### 2. NCBIのcutoff_score計算を確認
NCBIの`blast_parameters.c`の`BlastHitSavingParametersNew`を詳細に確認：
- `cutoff_score`の計算式
- `gap_trigger`との関係
- `cutoff_score_max`との関係

### 3. 考えられる初歩的な原因

1. **gap_triggerの値が間違っている**
   - NCBIの`BLAST_GAP_TRIGGER_PROT = 22.0` (bit score)
   - LOSATで正しく変換されているか？

2. **scale_factorが適用されていない/間違っている**
   - TBLASTXではscale_factor = 3.0
   - cutoff_score計算に正しく適用されているか？

3. **cutoff_score_maxがcutoff_scoreを上書きしている**
   - 787行: `cutoff = cutoff.min(cutoff_score_max)`
   - cutoff_score_maxの計算が間違っている可能性

4. **検索空間の計算が違う**
   - NCBIは`MIN(q_len, s_len) * s_len`を使用
   - LOSATで同じ計算をしているか？

5. **そもそもcutoff_scoreの適用タイミングが違う**
   - NCBIはどの段階でcutoff_scoreを適用しているか？

## 関連ファイル

- `LOSAT/src/algorithm/tblastx/utils.rs`: cutoff_score計算と適用（715-817行）
- `LOSAT/src/algorithm/tblastx/constants.rs`: 定数定義
- `LOSAT/src/stats/karlin.rs`: raw_score_from_evalue関数
- `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c`: NCBI実装

## 確認コマンド

```bash
# bit score分布の比較
awk -F'\t' '{print $12}' tests/ncbi_out/*.out | sort -n | uniq -c | head -30
awk -F'\t' '{print $12}' tests/losat_out/*.out | sort -n | uniq -c | head -30

# cutoff_score計算のデバッグ出力を有効化
# utils.rs 809行のeprintln!を有効化してビルド・実行
```

## 目標
- 出力ヒット数: 62,053件（NCBI同等）
- 最小bit score: 22.1（NCBI同等）
- 処理時間: 10秒以下（維持）


