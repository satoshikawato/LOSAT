# LOSAT TBLASTX 改善セッション引継ぎ（v3）

## 現状サマリー

### パフォーマンス
- **実行時間**: 約8秒（NZ_CP006932 self比較、1スレッド）
- **ビルド状態**: 正常

### NCBI BLASTとの整合性

| 比較タイプ | LOSAT | NCBI | 一致率 |
|-----------|-------|------|--------|
| Self (NZ_CP006932) | 59,984 | 62,053 | **96.7%** ✅ |
| Cross (AP027131 vs AP027133) | 4,035 | 14,871 | **27.1%** ❌ |
| Cross (AP027132 vs NZ_CP006932) | 55,586 | 62,575 | **88.8%** |
| Cross (AP027078 vs AP027131) | 3,380 | 30,175 | **11.2%** ❌ |

### グラフで確認された問題
1. **Accumulated Length vs Alignment Length**: BLAST+のピーク（20-100 aa）にLOSATが追いついていない
2. **Accumulated Length vs Identity**: BLAST+は40-60% identityに大きなピーク、LOSATは不足
3. **短いアライメント・低identityのヒットが不足**

---

## 発見された根本的な問題

### 問題1: 座標の不一致（共通率約2%）
```
AP027078 vs AP027131:
- NCBI unique coords: 30,175
- LOSAT unique coords: 3,380
- Common: 736 (約2%)
```
→ extensionロジックが根本的に異なる

### 問題2: E-value計算の差異
同じ座標のヒットでも：
| 座標 | NCBI E-value | LOSAT E-value |
|------|-------------|---------------|
| 101408-101563 | **0.0** | 4.8e-60 |
| 101642-102004 | **0.0** | 4.9e-77 |
| 102711-102764 | **0.0** | 6.6e0 |

→ NCBIはsum-statisticsリンキングでほぼ全ヒットをE-value=0.0に改善

### 問題3: bit_score分布の差異
```
NCBI (AP027078 vs AP027131): 最小bit_score 22.1、22-26範囲に多数
LOSAT: 最小bit_score 13.0、13-20範囲に少数
```
→ LOSATは低スコアヒットを生成しているが、NCBIの22-26範囲が不足

---

## これまでに試したこと

### 1. f593910スタイルのE-value計算
- alignment-length-based search space（min_space=60M）
- **結果**: ヒット数が85,091に増加（多すぎ）

### 2. diagonal計算の修正
- `diag = s_pos - q_pos`（f593910スタイル）
- mask_keyの計算をf593910スタイルに変更
- **結果**: ヒット数に大きな変化なし

### 3. mask更新タイミングの変更
- extension直後にmask更新（スコアチェック前）
- **結果**: 効果限定的

### 4. sum-statisticsリンキングの実装
- diagonal-basedのHSPリンキング
- large_gap_sum_e関数でE-value計算
- **結果**: Self比較は96.7%一致、Cross比較は不十分

---

## 重要なファイル

### メインソースファイル
- `LOSAT/src/algorithm/tblastx/utils.rs` - メイン処理（914行）
- `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` - リンキング実装（171行）
- `LOSAT/src/algorithm/tblastx/extension.rs` - extension処理
- `LOSAT/src/algorithm/tblastx/chaining.rs` - HSPチェーニング
- `LOSAT/src/algorithm/tblastx/constants.rs` - 定数定義

### 参照用NCBIファイル
- `LOSAT/tests/blast_out/*.tblastx.n1.out` - NCBI BLAST結果

---

## 今後の方針

### 優先度1: extensionロジックの調査と修正
座標の共通率が2%しかないのは根本的な問題。以下を調査：
1. **X-drop値の確認**: NCBIデフォルトとの比較
2. **seed処理**: threshold、word_size、two-hit windowの設定
3. **座標変換**: アミノ酸→ヌクレオチド変換の検証
4. **フレーム計算**: 6フレーム翻訳の検証

### 優先度2: sum-statisticsリンキングの改善
NCBIがE-value=0.0を実現している方法を調査：
1. **link_hsps.cの詳細分析**: NCBIの実装を正確に理解
2. **リンキング条件の見直し**: window_size、gap条件
3. **E-value計算式の検証**: large_gap_sum_eの実装確認

### 優先度3: パラメータチューニング
- `MIN_UNGAPPED_SCORE` (現在22)
- `seed_score`閾値 (現在11, 30)
- `TWO_HIT_WINDOW` (現在40)
- `X_DROP_UNGAPPED` (現在11)

---

## テストコマンド

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT

# ビルド
cargo build --release

# Self比較テスト
./target/release/losat tblastx \
  --query tests/fasta/NZ_CP006932.fasta \
  --subject tests/fasta/NZ_CP006932.fasta \
  --query-gencode 4 --db-gencode 4 \
  --num-threads 1 --out /tmp/test.out

# Cross比較テスト
./target/release/losat tblastx \
  --query tests/fasta/AP027078.fasta \
  --subject tests/fasta/AP027131.fasta \
  --query-gencode 4 --db-gencode 4 \
  --num-threads 1 --out /tmp/cross_test.out

# 結果比較
wc -l /tmp/test.out
grep -v "^#" tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out | wc -l
```

---

## 参考リンク

### NCBI BLASTソースコード参照箇所
- `ncbi-blast/c++/src/algo/blast/core/link_hsps.c` - sum-statisticsリンキング
- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c` - ungapped extension
- `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c` - パラメータ設定

---

## 注意事項

1. **全ゲノムヒットは不正確** - f593910は全ゲノムヒット（657,099 bp）を出力していたが、これは誤り。現在のビルドは正しくこれを除外している。

2. **git状態**: 多くのファイルが変更されている状態。必要に応じて`git stash`でクリーンにする。

3. **診断モード**: `LOSAT_TBLASTX_DIAGNOSTICS=1`環境変数で診断出力が有効になる。

