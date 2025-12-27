# NCBI BLAST実装との完全一致検証

## 実装状況

### ✅ 実効長計算は既に移植済み

LOSATの`src/stats/length_adjustment.rs`は、NCBI BLASTの`BLAST_ComputeLengthAdjustment`関数をRustに移植しています。

**参照**: 
- NCBI BLAST: `c++/src/algo/blast/core/blast_stat.c` の `BLAST_ComputeLengthAdjustment`
- LOSAT: `src/stats/length_adjustment.rs` の `compute_length_adjustment_ncbi`

## 実装比較

### NCBI BLAST (C) - 実際のコード

```c
if(converged) {
    *length_adjustment = (Int4) ell_min;  // 直接キャスト
    ell = ceil(ell_min);
    if( ell <= ell_max ) {
      ss = (m - ell) * (n - N * ell);
      if(alpha_d_lambda * (logK + log(ss)) + beta >= ell) {
        *length_adjustment = (Int4) ell;  // ceil値に更新
      }
    }
} else {
    *length_adjustment = (Int4) ell_min;  // 収束しなくても設定
}
```

### LOSAT (Rust) - 現在の実装

```rust
// Match NCBI BLAST's exact behavior: (Int4) ell_min (direct cast)
let mut length_adjustment = ell_min as i64;  // 直接キャスト（修正済み）

if converged {
    let ell_ceil = ell_min.ceil();
    if ell_ceil <= ell_max {
        let ss = (m - ell_ceil) * (n - n_seqs * ell_ceil);
        if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell_ceil {
            length_adjustment = ell_ceil as i64;  // ceil値に更新
        }
    }
}
// 収束しなかった場合もell_minが既に設定されている
```

## 実装の一致確認

### ✅ 一致している点

1. **イテレーションロジック**: 同じbisection-likeアルゴリズム
2. **収束判定**: `ell_bar - ell_min <= 1.0`
3. **初期値計算**: 同じ二次方程式を使用
4. **ceil値の検証**: 同じロジック

### ⚠️ 確認が必要な点

1. **浮動小数点精度**: CとRustで微妙に異なる可能性
2. **型変換**: `Int4` (32-bit) vs `i64` (64-bit) - 値の範囲は同じだが、型が異なる

## 検証方法

### ステップ1: NCBI BLASTで実効探索空間を取得

```bash
blastn -query test.fasta -subject db.fasta \
       -dust no \
       -reward 1 -penalty -2 \
       -gapopen 0 -gapextend 0 \
       -outfmt "7" \
       -out output.txt

# ヘッダーから実効探索空間を抽出
grep "Effective search space" output.txt
```

### ステップ2: LOSATで同じ計算を実行

```rust
use LOSAT::stats::search_space::SearchSpace;
use LOSAT::stats::tables::KarlinParams;

let params = KarlinParams { ... };
let search_space = SearchSpace::for_database_search(
    q_len, db_len, db_num_seqs, &params, true
);
println!("LOSAT: effective_space={}, length_adj={}",
         search_space.effective_space,
         search_space.length_adjustment);
```

### ステップ3: 比較

- **期待**: 完全一致または0.1%以下の差（浮動小数点精度の違いのみ）
- **許容誤差**: 0.1%以下
- **バグの可能性**: 1%以上の差が見られる場合は実装のバグを調査

## テスト戦略

### ユニットテスト

1. **実装の正確性**: アルゴリズムが正しく動作することを確認
2. **エッジケース**: 短い配列、長い配列、多数の配列など
3. **収束性**: 20イテレーション以内に収束することを確認

### NCBI BLASTとの比較テスト

1. **実効探索空間の比較**: NCBI BLASTの出力と比較（0.1%許容誤差）
2. **長さ調整値の比較**: 完全一致または±1単位の差
3. **様々なパラメータセット**: 複数のスコアリングパラメータで検証

## 次のステップ

1. ✅ 実装を修正（直接キャストに変更） - **完了**
2. ⏳ NCBI BLASTの出力から実効探索空間を抽出
3. ⏳ 実際のデータで比較テストを実行
4. ⏳ 完全一致を確認（または0.1%以下の差であることを確認）

## まとめ

**重要な発見**: LOSATは既にNCBI BLASTの実効長計算ロジックを移植しています。

**テストの役割**:
- 実装が正しく動作することを確認
- NCBI BLASTとの完全一致を検証（0.1%以下の許容誤差）
- 回帰を防ぐ

**許容誤差の理由**:
- 浮動小数点精度の違い（C vs Rust）
- 型の違い（Int4 vs i64、ただし値の範囲は同じ）

**バグの可能性**:
- 1%以上の差が見られる場合は、実装のバグを調査する必要がある
