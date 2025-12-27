# NCBI BLAST実装検証ガイド

## 概要

このディレクトリには、LOSATの実装がNCBI BLASTと完全一致することを検証するためのツールが含まれています。

**重要な発見**: LOSATは既にNCBI BLASTの`BLAST_ComputeLengthAdjustment`関数を直接移植しています。
この検証プロセスは、実装が正しく動作し、NCBI BLASTと完全一致することを確認します。

## 検証プロセス

### ステップ1: NCBI BLASTで参照データを生成

```bash
cd tests/unit/helpers
./generate_ncbi_reference.sh \
    ../../fasta/NZ_CP006932.fasta \
    ../../fasta/NZ_CP006932.fasta \
    test_case_1 \
    blastn
```

このスクリプトは：
- NCBI BLASTを実行（適切なパラメータで）
- 実効探索空間と長さ調整値を含む出力を生成
- ログファイルを保存

### ステップ2: 自動検証スクリプトを実行

```bash
./verify_ncbi_implementation.sh \
    ../../fasta/NZ_CP006932.fasta \
    ../../fasta/NZ_CP006932.fasta \
    test_case_1 \
    blastn
```

このスクリプトは：
1. NCBI BLASTを実行
2. 実効探索空間と長さ調整値を抽出
3. LOSATの計算結果と比較
4. テストケースを生成

### ステップ3: テストケースを実行

生成されたテストケースを実行：

```bash
cd ../../..  # LOSAT root
cargo test --test unit_tests test_case_1_effective_search_space -- --nocapture
```

## 期待される結果

### 実効探索空間
- **許容誤差**: 0.1%以下（浮動小数点精度の違いのみ）
- **バグの可能性**: 1%以上の差は実装のバグを示す可能性がある

### 長さ調整値
- **期待**: 完全一致
- **許容**: ±1単位（浮動小数点精度の違い）

## 手動検証

### NCBI BLAST出力から情報を抽出

```bash
# 実効探索空間を抽出
grep "Effective search space" ncbi_reference_data/test_case_1.out

# 長さ調整値を抽出
grep "Length adjustment" ncbi_reference_data/test_case_1.out

# クエリ長を抽出
grep "^# Query:" ncbi_reference_data/test_case_1.out

# データベース情報を抽出
grep "^# Database:" ncbi_reference_data/test_case_1.out
```

### LOSATで計算

```rust
use LOSAT::stats::search_space::SearchSpace;
use LOSAT::stats::tables::KarlinParams;

let params = KarlinParams {
    lambda: 1.28,
    k: 0.46,
    h: 0.85,
    alpha: 1.5,
    beta: -2.0,
};

let search_space = SearchSpace::for_database_search(
    query_len,
    db_len,
    db_num_seqs,
    &params,
    true,
);

println!("Effective space: {}", search_space.effective_space);
println!("Length adjustment: {}", search_space.length_adjustment);
```

## トラブルシューティング

### NCBI BLASTが見つからない

```bash
# PATHに追加
export PATH=/path/to/ncbi-blast/bin:$PATH

# または、スクリプトを修正して直接パスを指定
```

### 実効探索空間が抽出できない

NCBI BLASTのバージョンによっては、出力形式が異なる場合があります。
`outfmt 7`を使用していることを確認してください。

### テストが失敗する

1. **Karlin-Altschulパラメータを確認**: NCBI BLAST出力から正しい値を抽出
2. **許容誤差を確認**: 0.1%以下であることを確認
3. **実装を確認**: `src/stats/length_adjustment.rs`の実装を確認

## 次のステップ

1. ✅ 実装を確認（既に移植済み）
2. ✅ 検証ツールを作成
3. ⏳ 実際のデータで検証
4. ⏳ テストケースを統合

