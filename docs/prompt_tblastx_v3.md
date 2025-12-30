# LOSAT TBLASTXセッションプロンプト（v3）

## 背景

LOSATのTBLASTXアルゴリズムをNCBI BLASTと整合させる作業を行っています。

## 現在の状態

- **Self比較**（NZ_CP006932 vs 自身）: NCBI比96.7%の整合性 ✅
- **Cross比較**（異なる配列間）: NCBI比11-27%の整合性 ❌

## 発見された問題

### 1. 座標の不一致（最重要）
NCBIとLOSATで生成されるヒットの座標が全く異なる（共通率約2%）。
extensionロジックまたはフレーム計算に根本的な違いがある可能性。

### 2. E-value計算の差異
NCBIはsum-statisticsリンキングでほぼ全ヒットをE-value=0.0に改善しているが、LOSATは高いE-valueのまま。

### 3. ヒット分布の問題
グラフで確認：BLAST+は短いアライメント（20-100 aa）・低identity（20-60%）に多くのヒットを持つが、LOSATは不足。

## 今回のタスク

1. **extensionロジックの調査**: NCBIと座標が一致しない原因を特定
2. **sum-statisticsリンキングの改善**: E-value計算をNCBI互換に
3. **パラメータチューニング**: 必要に応じて閾値を調整

## 参照ファイル

詳細な状況は以下を参照：
- `docs/next_session_tblastx_v3.md`

主要ソースファイル：
- `LOSAT/src/algorithm/tblastx/utils.rs`
- `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`
- `LOSAT/src/algorithm/tblastx/extension.rs`

NCBI参照結果：
- `LOSAT/tests/blast_out/*.tblastx.n1.out`

## 目標

- Cross比較でもNCBI比80%以上の整合性を達成
- 実行速度は現状維持（約8秒）
- 全ゲノムヒットは出力しない（現状の正しい挙動を維持）


