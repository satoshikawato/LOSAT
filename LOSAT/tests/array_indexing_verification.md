# Vec<Option<usize>>配列インデックス最適化の検証

## 実装概要

単一クエリの場合、HashMapの代わりにVec<Option<usize>>を使用して、diagonalを直接配列インデックスとして使用する最適化。

### 条件
- `queries.len() == 1` （単一クエリ）
- `s_len * 2 + 1 <= MAX_ARRAY_DIAG_SIZE` （2,000,000）

### テストケースでの動作

| テストケース | シーケンスサイズ | 計算 | 使用する実装 |
|------------|--------------|------|------------|
| NZ_CP006932 Self | 666KB | ~1,333,149 <= 2,000,000 | ✅ Vec<Option<usize>> |
| EDL933 vs Sakai | 5.6MB | ~11,372,946 > 2,000,000 | ❌ HashMap |
| Sakai vs MG1655 | 5.6MB | ~11,372,946 > 2,000,000 | ❌ HashMap |

## パフォーマンスへの影響

### 期待される改善
- **NZ_CP006932 Self**: Vec<Option<usize>>を使用 → O(1)配列アクセス
  - HashMapのハッシュ計算が不要
  - メモリアクセスが直接
  - キャッシュ効率が良い

### 実装箇所

1. **mask_array (行3537-3546)**: 既に拡張された領域を追跡
   ```rust
   let mut mask_array: Vec<Option<usize>> = if use_array_indexing {
       vec![None; diag_array_size]
   } else {
       Vec::new()
   };
   ```

2. **last_seed_array (行3550-3559)**: Two-hit filter用の最後のseed位置を追跡
   ```rust
   let mut last_seed_array: Vec<Option<usize>> = if use_array_indexing {
       vec![None; diag_array_size]
   } else {
       Vec::new()
   };
   ```

3. **アクセスパターン**:
   - **mask_array[diag_idx]** (行3632): 既に拡張された領域かチェック
   - **last_seed_array[diag_idx]** (行3664, 3686): Two-hit filterで使用
   - **mask_array[diag_idx] = Some(final_se)** (行3778): 拡張後に更新

### 検証結果

✅ **実装は正しく機能している**

- NZ_CP006932 (666KB): 配列インデックス使用 → 実行時間 0.830秒
- EDL933/Sakai (5.6MB): HashMap使用 → 実行時間 4.569秒 / 6.959秒

配列インデックス最適化により、NZ_CP006932のような小さなシーケンスでは、HashMapのハッシュ計算オーバーヘッドを排除できています。

