# TBLASTX SEGフィルター実装完了 - 次のステップ

## 今回のセッションで達成したこと

### 1. SEGフィルターの実装
- **場所**: `LOSAT/src/utils/seg.rs`
- **機能**: NCBI BLAST互換のSEGアルゴリズム（低複雑度領域フィルター）
- **パラメータ**: window=12, locut=2.2, hicut=2.5（NCBI BLASTデフォルト）

### 2. TBLASTXへの統合
- **場所**: `LOSAT/src/algorithm/tblastx/utils.rs`, `args.rs`
- **動作**:
  1. 翻訳後の各フレームのアミノ酸配列にSEGを適用
  2. マスク領域のアミノ酸を'X'（コード23）に置換
  3. DNA座標に変換してlookupテーブルでもフィルタリング
- **効果**: 'X'は他のアミノ酸とのスコアが低いため、拡張が自然にマスク領域で停止

### 3. 結果
**NZ_CP006932 self-comparison:**
- **最大アラインメント長**: 113,828 aa → **2,253 aa** に改善！
- **目標の2,500 aa以下を達成**
- runaway extensionの問題は解決

## 残存する問題

### ヒット数の不足
| ツール | ヒット数 |
|--------|----------|
| NCBI BLAST+ | 62,053 |
| LOSAT | 9,789 |
| **差** | **約6倍** |

### グラフから見える傾向
1. **低identity（<80%）のヒットが大幅に不足**
2. **短いヒット（<30 aa）が不足**
3. 分布の形状はBLAST+と類似しているが、量が少ない

## 考えられる原因

### 1. SEGフィルターが過剰にマスク
- 現在: **41.33%のアミノ酸がマスク**（543,175 / 1,314,198 aa）
- DNA座標でのマスク: **87.87%**
- NCBI BLASTの実際のマスク率と比較が必要

### 2. Seed検出閾値が厳しすぎる
- 現在のseed_score最小値: 11（f593910で追加）
- NCBI BLASTはもっと低いスコアのシードも許容？

### 3. 隣接ワード閾値（neighboring words threshold）
- 現在: threshold=13
- NCBI BLASTと同じだが、実装が異なる可能性

### 4. Two-hit window size
- 現在: 40（デフォルト）
- NCBI BLAST互換モード: 16

## 実装済みの改善（最新セッション）

### 1. シードスコア閾値の調整可能化
- **実装**: `LOSAT/src/algorithm/tblastx/args.rs`に`--min-seed-score`引数を追加
- **デフォルト**: 11（従来通り）
- **使用方法**: `--min-seed-score 0`で無効化、`--min-seed-score 9`で緩和など
- **効果**: シードフィルタリングの影響を調査可能

### 2. 診断出力の強化
- **実装**: `LOSAT/src/algorithm/common/diagnostics.rs`にシードスコア統計を追加
- **追加情報**:
  - シードスコアの最小値・最大値・平均値
  - シードスコア分布の追跡
- **使用方法**: `LOSAT_DIAGNOSTICS=1`環境変数を設定して実行
- **効果**: シードフィルタリングの影響を詳細に分析可能

### 3. SEGマスク率比較ツール
- **実装**: `docs/compare_seg_mask.sh`スクリプトを作成
- **機能**: NCBI segmaskerとLOSATのマスク率を比較
- **使用方法**: `./docs/compare_seg_mask.sh <fasta_file>`

## 次のセッションでの調査・対策

### 優先度1: SEGマスク率の検証
```bash
# SEGマスク率比較スクリプトを使用
./docs/compare_seg_mask.sh tests/data/NZ_CP006932.fasta

# または手動でNCBI segmaskerを実行
segmasker -in tests/data/NZ_CP006932.fasta -outfmt interval -out seg_mask.intervals
```

### 優先度2: シードスコア閾値の調整テスト
```bash
# デフォルト（seed_score >= 11）
cd LOSAT
time ./target/release/losat tblastx \
  -q ../tests/data/NZ_CP006932.fasta \
  -s ../tests/data/NZ_CP006932.fasta \
  -o tests/losat_out/default.out \
  --seg 2>&1 | tee tests/losat_out/default.log

# 閾値を緩和（seed_score >= 9）
time ./target/release/losat tblastx \
  -q ../tests/data/NZ_CP006932.fasta \
  -s ../tests/data/NZ_CP006932.fasta \
  -o tests/losat_out/min_seed_9.out \
  --seg --min-seed-score 9 2>&1 | tee tests/losat_out/min_seed_9.log

# 閾値を無効化（seed_score >= 0、すべてのシードを許可）
time ./target/release/losat tblastx \
  -q ../tests/data/NZ_CP006932.fasta \
  -s ../tests/data/NZ_CP006932.fasta \
  -o tests/losat_out/min_seed_0.out \
  --seg --min-seed-score 0 2>&1 | tee tests/losat_out/min_seed_0.log

# 診断出力を有効にして詳細を確認
LOSAT_DIAGNOSTICS=1 time ./target/release/losat tblastx \
  -q ../tests/data/NZ_CP006932.fasta \
  -s ../tests/data/NZ_CP006932.fasta \
  -o tests/losat_out/diagnostic.out \
  --seg --min-seed-score 0 2>&1 | tee tests/losat_out/diagnostic.log
```

### 優先度3: 診断出力の分析
- `LOSAT_DIAGNOSTICS=1`で実行して以下を確認:
  - シードスコア分布（最小・最大・平均）
  - フィルタリングされたシード数（low score, mask, two-hit）
  - 拡張数とヒット数の関係

### 優先度4: パラメータ感度分析
- SEG locut/hicut を調整
- window sizeを変更

## 関連ファイル

### LOSAT実装
- `LOSAT/src/utils/seg.rs` - SEGアルゴリズム
- `LOSAT/src/algorithm/tblastx/utils.rs` - TBLASTX統合
- `LOSAT/src/algorithm/tblastx/args.rs` - コマンドライン引数

### テスト結果
- `tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out`
- `tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out`

### NCBI BLAST参照
- `ncbi-blast/c++/src/algo/blast/core/blast_seg.c`
- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c`

## 比較プロット
プロットは`tests/`ディレクトリに保存されている。
- TBLASTXの全体トレンド
- NZ_CP006932 self-comparison詳細比較

