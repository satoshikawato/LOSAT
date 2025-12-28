# TBLASTX NCBI BLAST互換性改善 - セッション引継ぎ文書

## 概要
LOSATのTBLASTX実装がNCBI BLASTと比較して大幅にヒット数が少ない問題の調査と修正を行った。

## 現在の状況

### 数値比較（自己比較テスト: NZ_CP006932.fasta）
| 指標 | LOSAT | NCBI BLAST | 比率 |
|------|-------|------------|------|
| 総ヒット数 | 1,595 | 62,053 | 2.6% |
| 最小ビットスコア | 32.2 | 22.1 | - |
| 20-29ビットスコア範囲のヒット | 0 | 27,985 | 0% |
| 30-39ビットスコア範囲のヒット | 65 | 10,045 | 0.6% |
| off-diagonalヒット | 293 | 49,358 | 0.6% |
| Two-hitフィルタ通過率 | 3% | 不明 | - |

### 診断結果（LOSAT two-hitモード）
```
K-mer matches found:        1,685,593
Seeds filtered (two-hit):   1,637,755 (97%)
Seeds passed to extension:  2,370
Filtered (cutoff_score):    667 (score range: 11-41)
E-value passed:             1,595
```

## 実装済み修正

### 1. Karlin-Altschulパラメータ（重要）
**ファイル**: `src/algorithm/tblastx/utils.rs` (line 106-115)
- TBLASTXは非ギャップアライメントを使用するため、非ギャップパラメータに変更
- λ=0.3176, K=0.134 (以前: λ=0.267, K=0.041)
- cutoff_score: 42 (以前: 46)

### 2. Two-hitフラグメカニズム
**ファイル**: `src/algorithm/tblastx/utils.rs`
- NCBI BLASTの`diag_array[].flag`に相当する`diag_extended`ハッシュマップを実装
- 拡張後、次のヒットは新しいシーケンスとして扱う
- 参照: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:519-606`

### 3. 近傍ワード生成
**ファイル**: `src/algorithm/tblastx/lookup.rs`
- `generate_neighborhood_words()`関数を実装
- BLOSUM62スコア閾値11でシード生成
- 参照: NCBI BLASTのneighborhood word concept

### 4. K-merエンコーディング修正（重要）
**ファイル**: `src/algorithm/tblastx/utils.rs`
- 修正前: `c1 * 625 + c2 * 25 + c3` (base 25)
- 修正後: `c1 * 676 + c2 * 26 + c3` (base 26)
- `lookup.rs`との整合性を確保

### 5. One-hitモードサポート
**ファイル**: `src/algorithm/tblastx/args.rs`, `utils.rs`
- `--window-size 0`でone-hitモードを有効化
- NCBI BLASTの`-window_size 0`に相当

### 6. マスクチェック修正
**ファイル**: `src/algorithm/tblastx/utils.rs`
- 修正前: `s_pos <= last_end`
- 修正後: `s_pos < last_end`
- NCBI BLASTの厳密な`<`比較に合わせた

## 根本的な問題（未解決）

### 問題の本質
Two-hitフィルタが97%のシードを除外している。これにより：
1. 拡張回数が少ない（2,370回 vs 推定100k+）
2. off-diagonalヒットが見つからない（293 vs 49,358）
3. 低スコアヒット（22-29ビットスコア）が0個

### 仮説
1. **近傍ワード生成の違い**: NCBI BLASTはより多くの近傍ワードを生成している可能性
2. **ルックアップテーブル構造の違い**: エントリ数や分布が異なる可能性
3. **K-merスキャン順序の違い**: 同じk-merでもスキャン順序が結果に影響する可能性

## 重要なファイルと場所

### LOSAT
- `LOSAT/src/algorithm/tblastx/utils.rs` - メインパイプライン
- `LOSAT/src/algorithm/tblastx/extension.rs` - 拡張アルゴリズム
- `LOSAT/src/algorithm/tblastx/lookup.rs` - ルックアップテーブル
- `LOSAT/src/algorithm/tblastx/constants.rs` - 定数定義
- `LOSAT/src/stats/tables.rs` - Karlin-Altschulパラメータ

### NCBI BLAST（参照用）
- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c` - Two-hit拡張
- `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c` - cutoff計算
- `ncbi-blast/c++/include/algo/blast/core/blast_options.h` - 定数

## テストコマンド

### LOSAT
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release

# Two-hitモード（デフォルト）
./target/release/losat tblastx --query tests/fasta/NZ_CP006932.fasta --subject tests/fasta/NZ_CP006932.fasta --query-gencode 4 --db-gencode 4

# One-hitモード
./target/release/losat tblastx --query tests/fasta/NZ_CP006932.fasta --subject tests/fasta/NZ_CP006932.fasta --query-gencode 4 --db-gencode 4 --window-size 0

# 診断付き
LOSAT_DIAGNOSTICS=1 ./target/release/losat tblastx ...
```

### NCBI BLAST
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT

# Two-hitモード（デフォルト）
tblastx -query tests/fasta/NZ_CP006932.fasta -subject tests/fasta/NZ_CP006932.fasta -query_gencode 4 -db_gencode 4 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# 出力ファイルは既存
tests/ncbi_out/NZ_CP006932.NZ_CP006932.tblastx.out
```

## 次のセッションで試すべきこと

### 優先度高
1. **ルックアップテーブル診断の追加**
   - エントリ数、バケットごとのヒット数分布を出力
   - NCBI BLASTと比較

2. **特定ヒットのトレース**
   - NCBI BLASTの低スコアヒット（bit score 24.9など）を選択
   - そのヒットの座標でLOSATがどう処理しているか追跡

3. **対角線ごとのk-mer一致数比較**
   - 特定のoff-diagonalでk-mer一致がどれだけあるか確認

### 優先度中
4. **NCBI BLASTのneighborhood word実装詳細確認**
   - `lookup_wrap.c`や関連ファイルを調査

5. **高スコアシード単独拡張の調整**
   - 現在seed_score >= 22で単独拡張可能
   - この閾値の妥当性を検証

## 注意事項

1. **パフォーマンス**: NCBI BLAST one-hitモードは非常に遅い（20分以上）
2. **テストファイル**: `NZ_CP006932.fasta`は660kb細菌ゲノム
3. **遺伝暗号**: 両方とも`--query-gencode 4 --db-gencode 4`を使用

