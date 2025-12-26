# NCBI BLAST vs LOSAT: 重要な発見

## 根本的な違い

### NCBI BLASTの実装（megablast、word_size=28）

1. **`lut_word_length = 8`** (lookup tableのインデックス化)
2. **`word_length = 28`** (extensionをトリガーする完全マッチ)
3. **`backbone_size = 1 << 16 = 65,536`** (8-letter word = 2^16)
4. **直接配列インデックス**: `thick_backbone[index]` でO(1)アクセス
5. **Presence Vector (PV)**: 事前フィルタリングで高速化
6. **`scan_step = 28 - 8 + 1 = 21`**: 大きなストライドでスキャン

### LOSATの実装

1. **`word_size = 28`** を直接使用してlookup
2. **HashMap使用**: `FxHashMap<u64, Vec<(u32, u32)>>`
3. **4^28 = 約72,000兆エントリ**: 直接配列は非現実的
4. **ハッシュ計算**: 各lookupでハッシュ計算が必要
5. **`scan_step = 4`**: LOSATの実装（word_size >= 16の場合）

## パフォーマンスへの影響

### NCBI BLASTのルックアップ
```c
// 2バイト（8塩基）を読み取り、直接インデックス化
Int4 index = s[0] << 8 | s[1];  // O(1) - ビット演算のみ
if (PV_TEST(pv, index, PV_ARRAY_BTS)) {  // O(1) - bit vector test
    num_hits = lookup->thick_backbone[index].num_used;  // O(1) - 直接配列アクセス
}
```

### LOSATのルックアップ
```rust
let current_kmer = /* 28-merを計算 */;  // O(1) - rolling k-mer
hash_lookup.get(&current_kmer)  // O(1) 平均だが、ハッシュ計算 + メモリアクセス
```

## 重要な違い

| 項目 | NCBI BLAST | LOSAT |
|------|-----------|-------|
| Lookup word length | 8 | 28 |
| Lookup table size | 65,536 (直接配列) | ~5.2M (HashMap) |
| アクセス方法 | `thick_backbone[index]` | `hash_lookup.get(&kmer)` |
| ハッシュ計算 | 不要 | 必要 |
| キャッシュ効率 | 高い（配列） | 低い（ハッシュテーブル） |
| scan_step | 21 (word_length - lut_word_length + 1) | 4 |

## NCBI BLASTのmegablast実装詳細

### lookup tableの構造（word_size=28）

```c
// blast_nalookup.c 行164-165
if (lookup_options->word_size >= 16) {
    *lut_width = 8;  // lookup tableは8-letter wordを使用
}

// 行392-397
lookup->word_length = 28;      // extensionをトリガーする完全マッチ
lookup->lut_word_length = 8;   // lookup tableのインデックス化
lookup->backbone_size = 1 << 16;  // 2^16 = 65,536（直接配列）
lookup->scan_step = 28 - 8 + 1 = 21;  // 大きなストライド
```

### スキャンループ（blast_nascan.c）

```c
// 2バイト（8塩基）を読み取り、直接インデックス化
for (; s <= s_end; s++) {
    Int4 index = s[0] << 8 | s[1];  // O(1) - ビット演算のみ
    
    // PV (Presence Vector) で事前フィルタリング
    if (PV_TEST(pv, index, PV_ARRAY_BTS)) {  // O(1) - bit vector test
        // 直接配列アクセス
        num_hits = lookup->thick_backbone[index].num_used;  // O(1)
        // hitsを取得
        s_BlastLookupRetrieve(lookup, index, ...);
    }
}
```

## LOSATとの比較

### LOSATの現在の実装

- `MAX_DIRECT_LOOKUP_WORD_SIZE = 13`
- `word_size = 28` → HashMap使用
- `scan_step = 4` (word_size >= 16の場合)
- 各lookupでハッシュ計算が必要

### 重要な違い

| 項目 | NCBI BLAST | LOSAT |
|------|-----------|-------|
| Lookup word | 8 | 28 |
| Table size | 65,536 (直接配列) | ~5.2M (HashMap) |
| scan_step | 21 | 4 |
| アクセス | `thick_backbone[index]` | `hash_lookup.get(&kmer)` |
| ハッシュ計算 | 不要 | 必要 |

## 解決策

### 推奨: 2段階lookup tableの実装

1. **短いlookup word (lut_word_length = 8-12)**で直接配列インデックス
   - `backbone_size = 1 << (2 * lut_word_length)` = 65,536 または 4,194,304
   - `thick_backbone[index]`でO(1)アクセス

2. **長いword (word_length = 28)**でextensionをトリガー
   - lookup table内に28-merの情報も保存
   - または、後で28-merマッチをチェック

3. **Presence Vector (PV)** を実装
   - 事前フィルタリングで不要なlookupを削減

4. **scan_stepを21に設定** (word_size=28, lut_word_length=8の場合)
   - NCBI BLASTと同じストライドを使用

## 期待される改善

- **Lookup速度**: HashMap → 直接配列で4-10倍改善
- **scan_step**: 4 → 21で約5倍のlookup削減
- **合計**: 約20-50倍の改善が期待される

**目標**: 4.6秒 → 1秒未満（実現可能）

## 実装の優先順位

### Phase 1: 2段階lookup tableの実装（最重要）

**変更点:**
- `lut_word_length = 8`（または12）で直接配列lookup tableを作成
- `word_length = 28`でextensionをトリガー
- `thick_backbone`構造を実装

**参考実装:**
- `ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c` (行392-411)
- `ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h` (行131-156)

### Phase 2: Presence Vector (PV) の実装

- ビットベクターで事前フィルタリング
- `PV_TEST`マクロの実装

### Phase 3: scan_stepの最適化

- `scan_step = word_length - lut_word_length + 1 = 21`に変更（word_size=28の場合）

### Phase 4: 2-bit packed sequenceの再検討

- 初期化コストを最小化した実装
- または、スキャンループ内でのみ使用

