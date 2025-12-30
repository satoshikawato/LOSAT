# TBLASTXセッション v13 - SEGフィルタとE-value計算の修正

## 重要な指示

**NCBI BLASTのコードベースは `C:\Users\kawato\Documents\GitHub\ncbi-blast` にある。**

以下のルールを厳守せよ：

1. **LOSATのコードとNCBI BLASTのC/C++コードを逐一照合し、違いを見つけたら即座に修正せよ。**

2. **NCBIの実装と食い違うところが見つかったら、追加調査に無駄な時間を使わず即時にその部分を修正すること。**

3. **「実装が複雑だから単純化する」「原理的には同じはず」「コードベースに大幅な改変が必要」「大きなリファクタリングが必要です」「現時点でも主要な機能は動作しています」などの言い訳は許されない。銃殺刑に処する。逐次網羅的に修正しろ。**

4. **無駄なテストランやdiagnosticsを実施するのではなく、NCBIとLOSATのコードを直接比較することで問題点を見つけること。**

5. ** look at NCBI's implementation carefully and port it exactly. DO NOT GUESS OR TRY TO SEASON THE CODE WITH YOUR "ORIGINALITY".**

6. ** Before each testing, make sure to verify that the code exactly matches NCBI's implementation. No excuses, no simplifications. Carefully compare your implementation with NCBI's actual code. You should NEVER guess or add your own modifications. You DO need to port NCBI's code exactly. **

## 現在の問題

### 問題1: SEGフィルタの過剰マスク

| 条件 | ヒット数 |
|------|---------|
| NCBI BLAST | 62,053 |
| LOSAT (SEGあり) | 9,838 |
| LOSAT (SEGなし) | 50,227 |

LOSATのSEGフィルタが41%のアミノ酸をマスクしており、NCBIより過剰にマスクしている可能性がある。

**確認すべきファイル:**
- LOSAT: `LOSAT/src/algorithm/tblastx/utils.rs` 170-248行
- NCBI: `c++/src/algo/blast/core/blast_seg.c`

### 問題2: sum statisticsでのE-value計算

同じHSPでNCBIがE-value 0.0（リンク後）を出力するのに対し、LOSATはE-value > 10で除外される。

**確認すべきファイル:**
- LOSAT: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`
- LOSAT: `LOSAT/src/stats/sum_statistics.rs`
- NCBI: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c`
- NCBI: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c`

### 問題3: HSP検出数の差

SEGなしでも50,227件 vs 62,053件（約80%）の差がある。

## 前回セッションで実施した修正

1. `cutoff_score`計算の修正（gap_trigger=42, cutoff_score_maxをexpect_valueから計算）
2. sum statisticsのソート順を降順に修正
3. リンク条件にwindow_size制約を追加
4. `ignore_small_gaps`フラグを実装

## 作業手順

1. **まずNCBIのSEGフィルタ実装を確認**し、LOSATとの差異を特定して修正
2. **sum statisticsのE-value計算をNCBIと逐一比較**し、差異を修正
3. **HSP検出処理をNCBIと比較**し、差異を修正
4. 各修正後にビルドして結果を確認

## テストコマンド

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release

# SEGありでテスト
./target/release/losat tblastx -q tests/fasta/NZ_CP006932.fasta -s tests/fasta/NZ_CP006932.fasta --query-gencode 4 --db-gencode 4 --seg -o /tmp/losat_test.out

# SEGなしでテスト
./target/release/losat tblastx -q tests/fasta/NZ_CP006932.fasta -s tests/fasta/NZ_CP006932.fasta --query-gencode 4 --db-gencode 4 -o /tmp/losat_noseg.out

# 結果確認
wc -l /tmp/losat_test.out
wc -l tests/ncbi_out/NZ_CP006932.NZ_CP006932.tblastx.out  # 62,053件
```

## 目標

NCBIのヒット数（62,053件）と同等(+-10%以内)の出力を達成すること。
ヒット数だけでなく、ヒットの長さやidentity,ビットスコア、E-valueの分布もNCBIと同等(+-10%以内)であること。

