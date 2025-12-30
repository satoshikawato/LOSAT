# テストヘルパーとユーティリティ

## 概要

このディレクトリには、ユニットテストで使用するヘルパー関数とユーティリティが含まれています。

## 重要な発見: 実装は既に移植済み

**✅ LOSATは既にNCBI BLASTの実効長計算ロジックを移植しています**

- **ファイル**: `src/stats/length_adjustment.rs`
- **関数**: `compute_length_adjustment_ncbi()`
- **参照**: NCBI BLAST `blast_stat.c` の `BLAST_ComputeLengthAdjustment`

これは、単純に「許容誤差を設ける」のではなく、**完全一致を目指す**ことが可能であることを意味します。

## ファイル構成

### `ncbi_reference.rs`
NCBI BLAST参照データと比較ユーティリティ

- `NcbiEvalueTestCase`: テストケース構造体
- `verify_effective_search_space()`: 実効探索空間の検証（0.1%許容誤差）
- `compare_evalue_with_ncbi()`: E-value比較
- `compare_bit_score_with_ncbi()`: Bit score比較

### `ncbi_reference_data.md`
NCBI BLAST参照データの準備方法（詳細ガイド）

### `TESTING_STRATEGY.md`
実装上の落とし穴と対策（4つの主要な問題点）

### `ncbi_blast_verification.md`
NCBI BLAST実装との完全一致検証方法

### `IMPLEMENTATION_STATUS.md`
実装移植状況の詳細

### `extract_ncbi_cases_enhanced.py`
NCBI BLAST出力からテストケースを抽出するスクリプト

## テスト戦略の更新

### 以前の理解（誤り）
- 実効長計算は複雑で、完全一致は不可能
- 1%の許容誤差を設ける必要がある

### 現在の理解（正しい）
- ✅ LOSATは既にNCBI BLASTの実装を移植済み
- ✅ 完全一致または0.1%以下の差が期待される
- ✅ より大きな差（>1%）は実装のバグを示す可能性がある

## 使用方法

### 1. NCBI BLAST参照データの準備

```bash
# BLASTNの場合
blastn -query test.fasta -subject db.fasta \
       -dust no \
       -reward 1 -penalty -2 \
       -gapopen 0 -gapextend 0 \
       -outfmt 7 > ncbi_output.txt

# テストケースを抽出
python extract_ncbi_cases_enhanced.py \
    ncbi_output.txt \
    <q_len> <db_len> <db_num_seqs> \
    blastn
```

### 2. テストケースの追加

生成されたRustコードを`ncbi_reference.rs`の関数に追加：

```rust
pub fn get_ncbi_blastn_evalue_cases() -> Vec<NcbiEvalueTestCase> {
    vec![
        // 生成されたテストケースをここに追加
    ]
}
```

### 3. テストの実行

```bash
# NCBI比較テストを実行
cargo test --test unit_tests evalue_ncbi

# 全てのユニットテストを実行
cargo test --test unit_tests
```

## 期待される結果

### 実効探索空間
- **許容誤差**: 0.1%以下（浮動小数点精度の違いのみ）
- **バグの可能性**: 1%以上の差は実装のバグを示す

### 長さ調整値
- **期待**: 完全一致
- **許容**: ±1単位（浮動小数点精度の違い）

### E-value
- **許容誤差**: 10-20%（実効探索空間の差が反映される）

## 注意事項

1. **フィルタリング**: 必ず `-dust no` (BLASTN) または `-seg no` (TBLASTX) を指定
2. **組成補正**: TBLASTXでは `-comp_based_stats 0` を指定
3. **パラメータ明示**: 全てのパラメータをコマンドライン引数で明示
4. **バージョン固定**: 同じNCBI BLASTバージョンを使用


