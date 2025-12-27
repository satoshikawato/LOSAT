# 次のセッションに向けたプロンプト

## 背景

LOSATのTBLASTX実装をNCBI BLASTと比較し、ヒット数が少ない問題（8.6% - 27.6%）を調査しています。これまでの調査で、主要な実装はNCBI BLASTと一致していることが確認されましたが、ヒット数が少ない問題は解決されていません。

## コードベースの場所

### LOSAT
- パス: `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/`
- 主要ファイル:
  - `src/algorithm/tblastx/utils.rs` - メインの実行ロジック
  - `src/algorithm/tblastx/extension.rs` - Extension実装
  - `src/algorithm/tblastx/constants.rs` - 定数定義
  - `src/stats/search_space.rs` - Effective search space計算
  - `src/stats/karlin.rs` - E-value計算
  - `TBLASTX_IMPLEMENTATION_REVIEW.md` - 調査結果の詳細

### NCBI BLAST
- パス: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/`
- 主要ファイル:
  - `c++/src/algo/blast/core/aa_ungapped.c` - Two-hit extension実装
  - `c++/src/algo/blast/core/blast_setup.c` - Effective search space計算
  - `c++/src/algo/blast/core/blast_parameters.c` - パラメータ設定
  - `c++/include/algo/blast/core/blast_options.h` - 定数定義

## テストケース

以下の5つのテストケースを使用しています：

```bash
# テストケース1: NZ_CP006932 self
time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/NZ_CP006932.fasta -s ./LOSAT/tests/fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1

# テストケース2: AP027132 vs NZ_CP006932
time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/AP027132.fasta -s ./LOSAT/tests/fasta/NZ_CP006932.fasta -o ./losat_out/AP027132.NZ_CP006932.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1

# テストケース3: AP027078 vs AP027131
time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/AP027078.fasta -s ./LOSAT/tests/fasta/AP027131.fasta -o ./losat_out/AP027078.AP027131.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1

# テストケース4: AP027131 vs AP027133
time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/AP027131.fasta -s ./LOSAT/tests/fasta/AP027133.fasta -o ./losat_out/AP027131.AP027133.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1

# テストケース5: AP027133 vs AP027132
time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/AP027133.fasta -s ./LOSAT/tests/fasta/AP027132.fasta -o ./losat_out/AP027133.AP027132.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1
```

**テスト結果**:
- LOSAT hits: 2,609 - 15,121
- NCBI hits: 14,876 - 62,580
- 比率: 8.6% - 27.6%

## 完了した調査項目

1. ✅ **MIN_UNGAPPED_SCORE削除** - NCBI BLASTは固定閾値を使用していない
2. ✅ **Two-hit extension確認** - 実装は正しい
3. ✅ **Diagonal suppression mask更新タイミング修正** - E-valueチェック後に更新
4. ✅ **Effective search space計算確認** - 計算式は一致
5. ✅ **Two-hit要件確認** - 実装は一致（71.5%のseedがフィルタリングされるのは想定内）

## 現在の状況

**診断情報（NZ_CP006932 self）**:
- K-mer matches: 59,256,542
- Seeds filtered (low score): 16,742,533
- Seeds filtered (two-hit): 42,418,658 (71.5%)
- Seeds passed to extension: 95,351
- E-value passed: 13,850
- E-value failed: 81,501
- **E-value通過率: 14.5%**

**問題点**:
- 実装はNCBI BLASTと一致しているが、ヒット数が大幅に少ない
- E-value通過率が低い（14.5%）
- 原因は不明（E-value threshold、seed score threshold、その他の要因の可能性）

## 次の調査項目

### 1. E-value thresholdの確認

**調査内容**:
- NCBI BLASTのデフォルトE-value thresholdを確認
- LOSATのデフォルトE-value thresholdと比較
- 実際のNCBI BLASTの出力からE-value分布を確認

**参照ファイル**:
- `ncbi-blast/c++/include/algo/blast/core/blast_parameters.h` - `CUTOFF_E_TBLASTX`
- `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c` - E-value threshold設定

**確認方法**:
```bash
# NCBI BLASTのデフォルトE-valueを確認
grep -r "CUTOFF_E_TBLASTX\|evalue.*tblastx" /mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/

# LOSATのデフォルトE-valueを確認
grep -r "evalue\|E_VALUE" /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/src/algorithm/tblastx/
```

### 2. Seed score thresholdの確認

**調査内容**:
- 現在13を使用（NCBI BLASTと一致）
- 低スコアのseedがフィルタリングされている可能性
- 診断情報: 16,742,533 seeds filtered (low score)

**参照ファイル**:
- `ncbi-blast/c++/include/algo/blast/core/blast_options.h:115` - `BLAST_WORD_THRESHOLD_TBLASTX`
- `LOSAT/src/algorithm/tblastx/utils.rs:279` - seed score threshold

### 3. Extensionの実装の再確認

**調査内容**:
- X-dropoffの適用方法
- 最大extension長の制限
- 高アイデンティティ領域での長いヒットの生成

**参照ファイル**:
- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1088-1158` - `s_BlastAaExtendTwoHit`
- `LOSAT/src/algorithm/tblastx/extension.rs` - `extend_hit_two_hit`

## 診断情報の有効化

診断情報を有効にするには：

```bash
export LOSAT_DIAGNOSTICS=1
export LOSAT_BIN=./LOSAT/target/release/losat

# 診断情報付きで実行
$LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/NZ_CP006932.fasta -s ./LOSAT/tests/fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.diag.out --query-gencode 4 --db-gencode 4 -n 1 2>&1 | grep -E "Effective search space|Seeds|E-value|two.*hit|low.*score"
```

## ビルド方法

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release
```

## 重要な発見

1. **Effective search spaceの計算は正しい**
   - 計算式: `effective_space = (query_len - length_adj) * (db_len - length_adj)`
   - 各contextごとに正しく計算されている

2. **Two-hit要件の実装は正しい**
   - 71.5%のseedがフィルタリングされるのは想定内（NCBI BLASTの動作と一致）

3. **E-value計算式は正しい**
   - `E = effective_space * 2^(-bit_score)`

4. **問題は実装の違いではなく、設定の違いである可能性が高い**
   - E-value threshold
   - Seed score threshold
   - その他のパラメータ

## 推奨される次のステップ

1. **NCBI BLASTの実際の出力を確認**
   - E-value分布
   - Bit score分布
   - ヒット長分布

2. **LOSATとNCBI BLASTの設定を比較**
   - デフォルトE-value threshold
   - Seed score threshold
   - X-dropoff値

3. **実際のNCBI BLASTの出力から逆算**
   - Effective search spaceの値
   - E-value thresholdの値

## 参考資料

- `TBLASTX_IMPLEMENTATION_REVIEW.md` - 詳細な調査結果
- NCBI BLASTソースコード - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/`
- LOSATソースコード - `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/`

