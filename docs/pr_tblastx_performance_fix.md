# TBLASTX パフォーマンス回復 (f593910レベルへ)

## コミット Description

```
fix(tblastx): restore f593910-level performance (30s → 7s)

Major performance regression fixes:
- MAX_HITS_PER_KMER: 1000 → 200 (reduce k-mer matches by 5.5x)
- Simplify two-hit logic to match f593910 approach
- Add seed_score early filtering (< 11 skip, < 30 require two-hit)
- X_DROP_UNGAPPED: restore to 11 (was changed to 7)
- cutoff_score_max: use MIN_UNGAPPED_SCORE (22)

Performance on NZ_CP006932 self-comparison:
- Before: 30+ seconds, 128,994 hits
- After: 6.8 seconds, 45,599 hits
- f593910: 6.0 seconds, 60,145 hits
- NCBI: 62,059 hits
```

---

## 概要

f593910以降に追加された変更により発生していた重大なパフォーマンス劣化を修正しました。

## 問題

| 指標 | 修正前 | 修正後 | f593910 | NCBI |
|------|--------|--------|---------|------|
| 実行時間 | 30秒+ | **6.8秒** | 6.0秒 | - |
| ヒット数 | 128,994 | 45,599 | 60,145 | 62,059 |

テストケース: `NZ_CP006932` self-comparison (657 kbp)

## 主な変更点

### 1. `lookup.rs`
- `MAX_HITS_PER_KMER`: 1000 → 200 に修正
- k-merマッチが2.7億→5千万に削減（5.5倍改善）

### 2. `constants.rs`
- `X_DROP_UNGAPPED`: 7 → 11 に復元

### 3. `utils.rs`
- two-hitロジックをf593910のシンプルな方式に戻し
- `seed_score < 11` の早期フィルタを追加
- `seed_score < 30 && !two_hit` のフィルタを追加
- extension分岐を `if let Some(prev_s_pos) = two_hit_info` 方式に変更
- `cutoff_score_max = MIN_UNGAPPED_SCORE` (22) を使用

## 診断比較

```
f593910:
  K-mer matches: 49,863,028
  Extensions: 1,892,458
  Final hits: 60,145

修正後:
  K-mer matches: 49,863,028
  Extensions: 1,968,086
  Final hits: 45,599
```

## 残課題

- ヒット数がf593910/NCBIより少ない（45,599 vs 60,145/62,059）
  - 原因候補: extension.rsに追加されたstop codonチェックの影響
  - cutoff_scoreフィルタでf593910より多く弾かれている（582,681 vs 466,503）

## テスト

```bash
cd LOSAT && cargo build --release
time ./target/release/losat tblastx \
  --query ./tests/fasta/NZ_CP006932.fasta \
  --subject ./tests/fasta/NZ_CP006932.fasta \
  --query-gencode 4 --db-gencode 4 --num-threads 1 \
  --out /tmp/test.out
```


