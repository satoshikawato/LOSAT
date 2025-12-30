# 次のセッション用プロンプト（tblastx / sum-stats）

以下を次セッションのプロンプトとしてコピーして使ってください。

---

## プロンプト

LOSAT の `tblastx` が NCBI BLAST+ と比べてヒット数（特に 20–29 bit score 帯）が少ない問題の調査/修正を継続してください。**いま残っている主原因は sum statistics（even-gap linking）と E-value の付与方法が NCBI と一致していない点**です。NCBI のコードを一次情報として、同等実装に寄せてください。

### 現状（ベースライン）

#### NCBI（one-hit: `-window_size 0`）
- 総ヒット数: 80,339
- 最小 bit score: 22.1
- 20–29 bit score: 44,563
- 30–39 bit score: 11,529

#### LOSAT（one-hit: `--window-size 0`、現状）
- 総ヒット数: 58,101
- 最小 bit score: 22.1
- 20–29 bit score: 7,808
- 30–39 bit score: 15,053

#### 追跡対象の低スコア例（NCBIで残る）
- 座標: 450063–449974 vs 1135–1224
- bit score: 24.9
- frame: qframe=-3, sframe=1
- NCBI（`-sum_stats T`）で evalue ≈ 0.007
- LOSAT は `--evalue 1e9` にすると出るが evalue ≈ 3.1e3 で、`--evalue 10` だと落ちる

### 既に直した乖離（実装済み）
- threshold: NCBI の `BLAST_WORD_THRESHOLD_TBLASTX=13` に合わせた
- lookup: 頻出 k-mer バケット丸ごと削除（`MAX_HITS_PER_KMER` 超で `clear()`）を撤廃（NCBI は overflow で保持）
- ungapped extension: stop codon を境界に止める（one-hit側も修正）
- E-value フィルタ: 拡張直後ではなく後段で行う（cutoff_score を通ったHSPを保持→後段で刈る）

### 次の最優先タスク（ここが本丸）
1. **NCBIの even-gap linking をRustへ移植**（`link_hsps.c`）
   - tblastx用の並び（例: `s_RevCompareHSPsTbx` 相当）とフレーム分離
   - `overlap_size=9`, `gap_size=40` の window/trim の扱い
   - small-gap / large-gap の2パス（`cutoff_small_gap`, `cutoff_big_gap`）
2. **sum-statistics の式を移植**（`blast_stat.c`）
   - `BLAST_SmallGapSumE`, `BLAST_LargeGapSumE`
   - `s_BlastSumP`（r<=4のテーブル補間 + r>4のRomberg積分）
   - Romberg積分: `ncbi_math.c:BLAST_RombergIntegrate` を移植
3. set-level evalue を各HSPに付与して `--evalue` を適用（NCBIと同じ段階）
4. 追跡対象HSPが `--evalue 10` でも残ることを確認し、20–29帯が NCBI に近づくか再集計

### 重要ファイル
- LOSAT: `LOSAT/src/algorithm/tblastx/utils.rs`, `lookup.rs`, `extension.rs`
- NCBI参照:
  - `ncbi-blast/c++/src/algo/blast/core/link_hsps.c`
  - `ncbi-blast/c++/src/algo/blast/core/blast_stat.c`
  - `ncbi-blast/c++/src/algo/blast/core/ncbi_math.c`
  - `ncbi-blast/c++/src/algo/blast/api/tblastx_options.cpp`
  - `ncbi-blast/c++/include/algo/blast/core/blast_options.h`
- 引継ぎ文書:
  - `docs/session_handoff_tblastx_sumstats.md`

### テストコマンド
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release

# LOSAT one-hit
LOSAT_DIAGNOSTICS=1 ./target/release/losat tblastx \
  --query tests/fasta/NZ_CP006932.fasta --subject tests/fasta/NZ_CP006932.fasta \
  --query-gencode 4 --db-gencode 4 --window-size 0 \
  --out /tmp/losat_tblastx_windowsize0.out
rg "450063\t449974\t1135\t1224" /tmp/losat_tblastx_windowsize0.out

# NCBI (frames付き)
tblastx -query tests/fasta/NZ_CP006932.fasta -subject tests/fasta/NZ_CP006932.fasta \
  -query_gencode 4 -db_gencode 4 -window_size 0 -sum_stats T \
  -outfmt "6 qstart qend sstart send qframe sframe evalue bitscore" \
| rg "450063\t449974\t1135\t1224"
```

### 目標
NCBIの `tblastx` の **even-gap linking + sum-stats E-value** を再現し、
低スコア帯の取りこぼしとヒット数乖離を解消すること。

---

