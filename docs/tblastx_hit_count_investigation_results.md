# TBLASTX ヒット数不足の調査結果

## 調査日
2024年12月30日

## 問題
- **NCBI BLAST+**: 62,053個のヒット
- **LOSAT**: 9,789個のヒット（約1/6）
- **差**: 約6倍

## 調査結果

### 1. シードスコア閾値の影響
- **実装**: `--min-seed-score`引数を追加（デフォルト: 11）
- **結果**: シードスコア範囲は13-31（平均14.17）のため、閾値11ではフィルタリングされていない
- **結論**: シードスコア閾値は問題ではない

### 2. Two-hit要件の影響
- **Two-hit (window=40)**: 6,516個のヒット
- **One-hit (window=0)**: 8,783個のヒット（+35%）
- **Two-hit (window=16)**: 6,380個のヒット
- **診断出力**: 67,200,016個のシードがtwo-hit要件でフィルタリング
- **結論**: Two-hit要件が主な要因の一つだが、それだけでは説明できない

### 3. SEGマスク率
- **LOSAT**: 47.50% (アミノ酸), 92.04% (DNA座標)
- **NCBI segmasker**: マスクなし（出力: "0 - 657100"）
- **結論**: SEGマスクが過剰にマスクしている可能性が高い

## 診断出力の詳細

### デフォルト設定（window=40, SEG有効）
```
K-mer matches found:        76,447,033
Seeds filtered (low score): 0
Seed score statistics:
  Count:                      76,351,066
  Range:                      13 - 31
  Average:                    14.17
Seeds suppressed (mask):    95,967
Seeds filtered (two-hit):   67,200,016  ← 主な要因
Seeds passed to extension:  448,034
Final hits:                 6,516
```

### One-hit mode（window=0, SEG有効）
```
K-mer matches found:        76,447,033
Seeds filtered (low score): 0
Seeds suppressed (mask):   8,836,736
Seeds filtered (two-hit):   0
Seeds passed to extension: 67,610,297
Final hits:                 8,783
```

## 次のステップ

### 優先度1: SEGマスク率の検証
- NCBI BLAST+の実際のSEGマスク率を確認
- LOSATのSEG実装が正しいか検証
- SEGパラメータ（locut, hicut, window）の調整

### 優先度2: Two-hit要件の詳細分析
- NCBI BLAST+のtwo-hit実装と比較
- Window sizeの影響をさらに調査
- 隣接ワード閾値（threshold=13）の影響

### 優先度3: その他の要因
- Cutoff scoreの影響
- E-value計算の違い
- 拡張アルゴリズムの違い

## 関連ファイル
- 診断ログ: `LOSAT/tests/losat_out/NZ_CP006932.*.log`
- 出力ファイル: `LOSAT/tests/losat_out/NZ_CP006932.*.out`
- 実装: `LOSAT/src/algorithm/tblastx/utils.rs`
- 引数: `LOSAT/src/algorithm/tblastx/args.rs`


