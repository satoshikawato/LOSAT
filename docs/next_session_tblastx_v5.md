# LOSAT TBLASTX 改善セッション引継ぎ（v5）

## 今回のセッションでの変更

### 1. 定数の変更
- `MIN_UNGAPPED_SCORE`: 22 → 14
- `X_DROP_UNGAPPED`: 11 → 7 (NCBIデフォルト)

### 2. sum_stats_linking.rs の改善
- `small_gap_sum_e` と `large_gap_sum_e` の両方を使用
- 最も良い（低い）E-valueを選択
- 効率的なper-frame linkingを維持

### 3. 試した改善策（効果なし）
- **Cross-frame linking**: 処理が遅すぎて実用的でない
- **マスキング緩和**: ヒット数は増えるが座標の一致率は改善せず

## 現在の状態

### パフォーマンス
| 比較タイプ | LOSAT | NCBI | 処理時間 |
|-----------|-------|------|---------|
| Cross (AP027078 vs AP027131) | 5,964 | 30,175 | 39秒 |
| Self (NZ_CP006932) | ~100,000 | 62,053 | - |

### 座標の共通率
**約2%**（575/30,175）- マスキング緩和でも改善しなかった

## 根本的な問題（重要）

### 座標の不一致が本質的な問題

マスキングやsum-statisticsリンキングの調整では**座標の一致率が改善しない**ことが判明。

**原因**: extensionロジックの根本的な違い
- NCBIとLOSATは同じシード位置から異なる座標のHSPを生成
- これはsum-statisticsやマスキングの問題ではなく、**ungapped extension**の実装の違い

### マスキング緩和実験の結果

| MASK_OVERLAP_ALLOW | ヒット数 | ユニーク座標 | 共通座標 |
|-------------------|---------|------------|---------|
| 0 (元の値) | 5,964 | ~5,800 | 575 |
| 20 | 14,920 | 6,167 | 575 |
| 50 | 21,385 | 6,160 | 575 |

→ ヒット数は増えるが**座標の一致率は変わらない**

## 次のセッションでの優先タスク

### 優先度1: Extensionロジックの詳細比較（最重要）

NCBIの`s_BlastAaExtendTwoHit`とLOSATの`extend_hit_two_hit`を詳細比較：

**NCBI** (`aa_ungapped.c`):
```c
score = s_BlastAaExtendTwoHit(matrix, subject, query,
                              last_hit + wordsize,  // ← 開始位置
                              subject_offset, query_offset,
                              cutoffs->x_dropoff, 
                              &hsp_q, &hsp_s, &hsp_len, ...);
```

**確認すべき点**:
1. X-dropの適用方法
2. 最大スコア位置の決定方法
3. 左右extension のマージ方法
4. アライメント境界の決定

### 優先度2: シード処理の比較

NCBIのseed処理（`scansub`関数）とLOSATの比較：
- k-merの生成方法
- マッチング条件
- 対角線の計算方法

## 技術的詳細

### Sum-statisticsリンキングの限界

現在の実装では、ほとんどのチェーンが2つのHSPしか含まない：
- 15,336チェーン: 2 HSPs
- 543チェーン: 3 HSPs
- 以降は急減

これは、LOSATが長いHSPを1つ生成するため、リンク可能なHSPが少ないことが原因。

### xsum計算の問題

xsum=20程度では：
```
pair_search_space.ln() ≈ 24.5 (4.53e10)
adjusted_xsum = 20 - 24.5 - ... = 負の値
→ E-value = 2.15e9 (INT4_MAX)
```

より多くのHSPをリンクしてxsumを増やす必要があるが、そもそもリンク可能なHSPが不足。

## 変更されたファイル

| ファイル | 変更内容 |
|---------|---------|
| `constants.rs` | MIN_UNGAPPED_SCORE=14, X_DROP_UNGAPPED=7 |
| `sum_stats_linking.rs` | small/large gap両方使用、効率化 |
| `utils.rs` | マスキング緩和を試したが元に戻した |

## 結論

**座標の一致率を改善するには、extensionロジックの詳細な比較と修正が必要**。
sum-statisticsやマスキングの調整は効果がない。
