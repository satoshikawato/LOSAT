# NCBI BLAST実装移植の詳細

## 実装状況

### ✅ 実効長計算は既に移植済み

LOSATの`src/stats/length_adjustment.rs`は、NCBI BLASTの`BLAST_ComputeLengthAdjustment`関数を**直接移植**しています。

**移植元**:
- ファイル: `ncbi-blast/c++/src/algo/blast/core/blast_stat.c`
- 関数: `BLAST_ComputeLengthAdjustment(double K, double logK, double alpha_d_lambda, double beta, Int4 query_length, Int8 db_length, Int4 db_num_seqs, Int4 * length_adjustment)`

**移植先**:
- ファイル: `LOSAT/src/stats/length_adjustment.rs`
- 関数: `compute_length_adjustment_ncbi(query_length: i64, db_length: i64, db_num_seqs: i64, params: &KarlinParams) -> LengthAdjustmentResult`

## 実装の一致確認

### アルゴリズムの一致

両方とも同じbisection-likeイテレーションアルゴリズムを使用：

1. **初期値計算**: 二次方程式で`ell_max`を計算
2. **イテレーション**: `ell_bar = alpha_d_lambda * (logK + log(ss)) + beta`を計算
3. **収束判定**: `ell_bar - ell_min <= 1.0`
4. **最終値決定**: `ell_min`を直接キャスト（またはceil値に更新）

### コードの対応関係

| NCBI BLAST (C) | LOSAT (Rust) |
|----------------|--------------|
| `(Int4) ell_min` | `ell_min as i64` |
| `ceil(ell_min)` | `ell_min.ceil()` |
| `MAX(m, n)` | `m.max(n)` |
| `sqrt(...)` | `discriminant.sqrt()` |
| `log(...)` | `k.ln()` |

### 型の違い

- **NCBI BLAST**: `Int4` (32-bit signed integer)
- **LOSAT**: `i64` (64-bit signed integer)

**影響**: 値の範囲は同じ（実用的な範囲内では）ため、動作は同じです。

## 検証方法

### 1. 単体テスト

```bash
cargo test stats::length_adjustment
```

既存のテストが全て通ることを確認。

### 2. NCBI BLASTとの直接比較

実際のNCBI BLAST出力と比較：

```bash
# NCBI BLASTを実行
blastn -query test.fasta -subject db.fasta \
       -dust no -reward 1 -penalty -2 \
       -gapopen 0 -gapextend 0 \
       -outfmt 7 > ncbi_output.txt

# ヘッダーから実効探索空間を抽出
grep "Effective search space" ncbi_output.txt
```

### 3. 実装の詳細比較

NCBI BLASTのコードとLOSATのコードを1行ずつ比較：

```bash
# NCBI BLASTの実装を確認
grep -A 100 "BLAST_ComputeLengthAdjustment" \
    /mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c

# LOSATの実装を確認
cat LOSAT/src/stats/length_adjustment.rs
```

## 期待される結果

### 完全一致が期待される理由

1. **直接移植**: NCBI BLASTのコードを1行ずつ移植
2. **同じアルゴリズム**: 同じbisection-likeイテレーション
3. **同じ数式**: 同じ計算式を使用

### 許容誤差

- **実効探索空間**: 0.1%以下（浮動小数点精度の違いのみ）
- **長さ調整値**: 完全一致または±1単位
- **E-value**: 10-20%（実効探索空間の差が反映される）

### バグの可能性

1%以上の差が見られる場合は、以下の可能性があります：

1. **実装のバグ**: 移植時のミス
2. **数値計算の違い**: 浮動小数点の扱いの違い
3. **型変換の違い**: CとRustの型変換の違い

## 次のステップ

1. ✅ 実装を確認（既に移植済み） - **完了**
2. ✅ 細かい違いを修正（直接キャストに変更） - **完了**
3. ⏳ NCBI BLASTの出力から実効探索空間を抽出
4. ⏳ 実際のデータで比較テストを実行
5. ⏳ 完全一致を確認（または0.1%以下の差であることを確認）

## まとめ

**重要な発見**: LOSATは既にNCBI BLASTの実効長計算ロジックを移植しています。

**テストの役割**:
- 実装が正しく動作することを確認
- NCBI BLASTとの完全一致を検証
- 回帰を防ぐ

**許容誤差の理由**:
- 浮動小数点精度の違い（C vs Rust）
- 型の違い（Int4 vs i64、ただし値の範囲は同じ）

**完全一致を目指す**:
- 実装は既に移植済み
- 細かい違いを修正済み
- 実際のデータで検証が必要

