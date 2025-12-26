# 最適化の検証まとめ

## 実装完了した最適化

### ✅ 1. Vec<Option<usize>>配列インデックス最適化（PR #18）

**実装内容:**
- 単一クエリ時、HashMapの代わりにVec<Option<usize>>を使用
- diagonalを直接配列インデックスとして使用（O(1)アクセス）
- 大きな配列（>500KB subject）ではHashMapにフォールバック

**条件:**
- `queries.len() == 1` （単一クエリ）
- `s_len * 2 + 1 <= MAX_ARRAY_DIAG_SIZE` (2,000,000)

**検証結果:**
| テストケース | シーケンスサイズ | 使用する実装 | 実行時間 |
|------------|--------------|------------|---------|
| NZ_CP006932 Self | 666KB | ✅ Vec<Option<usize>> | 0.830秒 |
| EDL933 vs Sakai | 5.6MB | ❌ HashMap | 4.569秒 |
| Sakai vs MG1655 | 5.6MB | ❌ HashMap | 6.959秒 |

**効果:**
- NZ_CP006932のような小さなシーケンスでは、HashMapのハッシュ計算オーバーヘッドを排除
- 直接メモリアクセスによりキャッシュ効率が向上

### ❌ 2. PackedSequence実装（リバート済み）

**問題点:**
- 初期化コストがO(n)で大きい
- get_base()のオーバーヘッド（除算、剰余、ビットシフト）
- is_ambiguous()のバイナリサーチコスト
- 元のENCODE_LUT実装の方が高速

**パフォーマンス比較:**
| テストケース | ベースライン | PackedSequence実装 | 差分 |
|------------|-----------|------------------|------|
| NZ_CP006932 Self | 0.830秒 | 0.915秒 | +10.2% |
| EDL933 vs Sakai | 4.569秒 | 5.307秒 | +16.2% |
| Sakai vs MG1655 | 6.959秒 | 7.011秒 | +0.7% |

**結論:**
- PackedSequence実装はリバートし、元のENCODE_LUT + rolling k-mer scanner実装を維持

## 修正したエラー

1. **Hit構造体のフィールド名修正**: 
   - `aln_len` → `length`
   - `mismatches` → `mismatch`
   - `gap_opens` → `gapopen`
   - `evalue` → `e_value`

2. **型の不一致修正**:
   - `as u32`キャストを削除し、`usize`に統一

## パフォーマンステスト結果

### ベースライン（リバート後 - 現在の実装）

**BLASTN (Megablast):**
- NZ_CP006932 Self: **0.830秒** ✅（目標: 1秒未満）
- EDL933 vs Sakai: 4.569秒
- Sakai vs MG1655: 6.959秒

**BLASTN (Task: blastn):**
- 10テストケースすべて完了

**TBLASTX:**
- 5テストケースすべて完了

## 結論

1. **Vec<Option<usize>>配列インデックス最適化**: ✅ 正しく機能している
2. **PackedSequence実装**: ❌ パフォーマンス低下のためリバート済み
3. **元の実装（ENCODE_LUT + rolling k-mer scanner）**: ✅ 最適化済みで高速

現在の実装は最適化されており、特にNZ_CP006932のような小さなシーケンスでは配列インデックス最適化により高速化されています。

