# NCBI BLAST実装移植状況

## 実効長計算（Length Adjustment）

### ✅ 実装済み

LOSATは既にNCBI BLASTの`BLAST_ComputeLengthAdjustment`関数をRustに移植しています：

- **ファイル**: `src/stats/length_adjustment.rs`
- **関数**: `compute_length_adjustment_ncbi()`
- **参照**: NCBI BLAST `blast_stat.c` の `BLAST_ComputeLengthAdjustment`

### 実装比較

#### NCBI BLAST (C)
```c
*length_adjustment = (Int4) ell_min;  // 直接キャスト
```

#### LOSAT (Rust) - 修正前
```rust
let mut length_adjustment = ell_min.floor() as i64;  // floor()を使用
```

#### LOSAT (Rust) - 修正後
```rust
let mut length_adjustment = ell_min as i64;  // 直接キャスト（NCBI BLASTと同じ）
```

**注意**: 正の値では `f64 as i64` と `floor() as i64` は同じ結果になりますが、
完全一致を目指すため、NCBI BLASTと同じ直接キャストを使用します。

### 検証方法

1. **単体テスト**: `cargo test stats::length_adjustment`
2. **NCBI BLASTとの比較**: 実際のNCBI BLAST出力と比較
3. **許容誤差**: 0.1%以下（浮動小数点精度の違いのみ）

### 次のステップ

1. ✅ 実装を修正（直接キャストに変更） - **完了**
2. ⏳ NCBI BLASTの出力から実効探索空間を抽出してテスト
3. ⏳ 実装が完全一致することを確認

## 実効探索空間（Effective Search Space）

### ✅ 実装済み

LOSATは実効探索空間の計算も実装済みです：

- **ファイル**: `src/stats/search_space.rs`
- **関数**: `SearchSpace::for_database_search()`
- **参照**: NCBI BLASTの実効探索空間計算ロジック

### 計算式

```
effective_query_len = query_len - length_adjustment
effective_db_len = db_len - (db_num_seqs * length_adjustment)
effective_space = effective_query_len * effective_db_len
```

これはNCBI BLASTの計算式と一致しています。

## テスト戦略の更新

### 以前のアプローチ（許容誤差1%）
- 実装が複雑なため、1%の許容誤差を設ける
- 統計的に意味のある範囲内で一致することを確認

### 現在のアプローチ（実装移植済み）
- ✅ 実装は既に移植済み
- ✅ 完全一致または非常に近い値が期待される
- ✅ 許容誤差は0.1%以下（浮動小数点精度の違いのみ）
- ⚠️ より大きな差（>1%）が見られる場合は、実装のバグを調査

## まとめ

**良いニュース**: LOSATは既にNCBI BLASTの実効長計算ロジックを移植しています。

**次のステップ**:
1. 実装の細かい違いを修正（直接キャストに変更） - **完了**
2. NCBI BLASTの出力と実際に比較して検証
3. 完全一致を確認（または0.1%以下の差であることを確認）

**テストの役割**:
- 実装が正しく動作することを確認
- NCBI BLASTとの完全一致を検証
- 回帰を防ぐ


