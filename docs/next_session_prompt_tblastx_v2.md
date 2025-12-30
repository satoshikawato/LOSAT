# LOSAT TBLASTX 次セッション引継ぎ

## 現状

f593910レベルのパフォーマンスを回復しました（30秒+ → 6.8秒）。
しかし、ヒット数がf593910/NCBIより少ない問題が残っています。

### ベンチマーク結果 (NZ_CP006932 self-comparison)

| 指標 | 現在 | f593910 | NCBI |
|------|------|---------|------|
| 実行時間 | 6.8秒 | 6.0秒 | - |
| ヒット数 | 45,599 | 60,145 | 62,059 |

---

## 次の課題

### 1. ヒット数差異の調査（優先）

**現在 vs f593910:**
- Extensions: ほぼ同じ（1,968,086 vs 1,892,458）
- cutoff_scoreフィルタ: 多く弾かれている（582,681 vs 466,503）

**原因候補:**
1. `extension.rs`のstop codonチェック（f593910: 4箇所、現在: 9箇所）
2. extension内部ロジックの差異

**調査コマンド:**
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT
git diff f593910 HEAD -- LOSAT/src/algorithm/tblastx/extension.rs | less
```

### 2. sum-statistics再実装

現在の`sum_stats_linking.rs`は無効化されています（passthrough）。
NCBIスタイルのeven-gap linkingを正しく実装すると、E-value計算が改善される可能性があります。

**参照ファイル:**
- NCBI: `ncbi-blast/c++/src/algo/blast/core/link_hsps.c`
- LOSAT: `LOSAT/src/stats/sum_statistics.rs`（関数は実装済み、リンキングが未実装）

### 3. neighboring words

現在は完全一致のみ（f593910方式）。
NCBIはthreshold=13でneighboring wordsを使用。

**以前の実装は遅すぎた（構築時に全neighbor展開）。**

代替案:
- スキャン時にneighborテーブルを使用（以前試したが3分かかった）
- 事前計算neighbor展開を高速化

---

## 重要ファイル

| ファイル | 説明 |
|----------|------|
| `LOSAT/src/algorithm/tblastx/utils.rs` | メインロジック |
| `LOSAT/src/algorithm/tblastx/extension.rs` | ungapped extension |
| `LOSAT/src/algorithm/tblastx/lookup.rs` | k-mer lookup |
| `LOSAT/src/algorithm/tblastx/constants.rs` | パラメータ |
| `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` | リンキング（現在無効） |
| `LOSAT/src/stats/sum_statistics.rs` | sum-statistics関数 |

---

## テストコマンド

### 高速テスト
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT

time ./target/release/losat tblastx \
  --query ./tests/fasta/NZ_CP006932.fasta \
  --subject ./tests/fasta/NZ_CP006932.fasta \
  --query-gencode 4 --db-gencode 4 --num-threads 1 \
  --out /tmp/losat.out
```

### 診断付き
```bash
LOSAT_DIAGNOSTICS=1 ./target/release/losat tblastx \
  --query ./tests/fasta/NZ_CP006932.fasta \
  --subject ./tests/fasta/NZ_CP006932.fasta \
  --query-gencode 4 --db-gencode 4 --num-threads 1 \
  --out /tmp/losat.out 2>&1 | tail -40
```

### f593910でテスト
```bash
git stash && git checkout f593910
cd LOSAT && cargo build --release
time ./target/release/losat tblastx \
  --query ./tests/fasta/NZ_CP006932.fasta \
  --subject ./tests/fasta/NZ_CP006932.fasta \
  --query-gencode 4 --db-gencode 4 --num-threads 1 \
  --out /tmp/f593910.out
git checkout tblastx_fix && git stash pop
```

---

## NCBIリファレンス

### コードベース

```
C:\Users\kawato\Documents\GitHub\ncbi-blast
```

**重要ファイル:**
- `c++/src/algo/blast/core/link_hsps.c` - HSPリンキング
- `c++/src/algo/blast/core/blast_stat.c` - 統計計算、length adjustment
- `c++/src/algo/blast/core/aa_ungapped.c` - ungapped extension
- `c++/src/algo/blast/core/lookup_wrap.c` - lookup table

### テストコマンド

```bash
# NCBI tblastx
tblastx -query tests/fasta/NZ_CP006932.fasta \
  -subject tests/fasta/NZ_CP006932.fasta \
  -query_gencode 4 -db_gencode 4 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out \
  -num_threads 1
```

---

## 変更履歴

### 今回の修正 (パフォーマンス回復)

1. **`MAX_HITS_PER_KMER`**: 1000 → 200
2. **two-hitロジック**: f593910のシンプルな方式に戻し
3. **seed_score早期フィルタ**: `< 11` と `< 30 && !two_hit` を追加
4. **`X_DROP_UNGAPPED`**: 7 → 11 に復元
5. **`cutoff_score_max`**: MIN_UNGAPPED_SCORE (22) を使用

