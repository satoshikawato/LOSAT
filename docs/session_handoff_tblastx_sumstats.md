# TBLASTX NCBI互換性改善 - セッション引継ぎ（2025-12-28）

## 目的（変わらず）
LOSATの `tblastx` を **NCBI BLAST+（tblastx）に可能な限り忠実**に合わせ、ヒット数（特に低スコア帯）と座標/統計（E-value/bit score）を揃える。

---

## 背景（ベースライン）
### NCBI（one-hit: `-window_size 0`）
`tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.windowsize0.out` の集計:
- 総ヒット数: **80,339**
- 最小 bit score: **22.1**
- 20–29 bit score: **44,563**
- 30–39 bit score: **11,529**

### LOSAT（このセッション開始時点のone-hit）
（`--window-size 0` 実行）
- 総ヒット数: **1,238**
- 最小 bit score: **32.2**
- 20–29 bit score: **0**

---

## このセッションで確定した「NCBI一次情報」（重要）
このセッションでは **LOSATのコメントは信用せず**、`ncbi-blast/` の該当ソースを直接読んで挙動を固定した。

### 1) tblastxのデフォルト閾値（threshold）は 13
- NCBI: `BLAST_WORD_THRESHOLD_TBLASTX 13`（`ncbi-blast/c++/include/algo/blast/core/blast_options.h`）
- `tblastx` の pairwise 出力にも `Neighboring words threshold: 13` が出るのを確認

### 2) tblastxはデフォルトで sum statistics を有効化
- NCBI: `CTBlastxOptionsHandle::SetHitSavingOptionsDefaults()` が `SetSumStatisticsMode()` を呼ぶ  
  （`ncbi-blast/c++/src/algo/blast/api/tblastx_options.cpp`）
- `-sum_stats T` がデフォルト相当であることをローカル実行で確認

### 3) lookup（近傍ワード）側で「頻出k-merバケットを丸ごと捨てる」挙動はNCBIに無い
- NCBIのAA lookupは `blast_aalookup.c` で **可変長ヒットリストを overflow で保持**する（`AA_HITS_PER_CELL=3` など）
- LOSATがやっていた `MAX_HITS_PER_KMER 超で bucket.clear()` は **明確な乖離**で、感度を大きく落としていた

### 4) ungapped extension は stop codon を境界に止まる
- NCBI実装（`aa_ungapped.c`）の挙動に合わせ、LOSAT one-hit extension も stop で止める必要がある

### 5) two-hit の diag flag/last_hit 更新は NCBI 由来の規則がある
- NCBI: `aa_ungapped.c` で `diag_array[diag].flag` と `last_hit` を用い、
  `right_extend` の場合は `last_hit = s_last_off - (wordsize-1) + diag_offset` に更新する

---

## このセッションでLOSATに入れた修正（実装済み）
### A) lookup診断の追加
- `LOSAT/src/algorithm/tblastx/lookup.rs`
  - `LOSAT_DIAGNOSTICS=1` 時に、lookup 生成の統計（総エントリ、非空バケット、上位バケット等）を出力

### B) threshold の取り扱い修正（NCBI=13へ）
- `LOSAT/src/algorithm/tblastx/lookup.rs`: `build_direct_lookup()` のデフォルト閾値を **13** に変更
- `LOSAT/src/algorithm/tblastx/utils.rs`: lookup 構築で **`args.threshold`** を反映（CLIの `--threshold` が効く）

### C) 「頻出バケット丸ごと削除」を撤廃（NCBIのoverflow保持に寄せる）
- `LOSAT/src/algorithm/tblastx/lookup.rs`
  - `MAX_HITS_PER_KMER` 超で `bucket.clear()` していた処理を撤廃（保持する）

### D) one-hit ungapped extension の stop-codon 停止
- `LOSAT/src/algorithm/tblastx/extension.rs`
  - `extend_hit_ungapped()` の seed/左右拡張で stop codon に当たったら拡張を止める
  - これにより異常に長い拡張（自己比較で巨大bit score）が解消（max bitscore が正常化）

### E) E-value フィルタを「拡張直後」から「後段」へ移動
- `LOSAT/src/algorithm/tblastx/utils.rs`
  - NCBIに合わせて、拡張段階では `cutoff_score` のみで足切りして HSP を保持し、
    `--evalue` のフィルタは後段（writer側）で行うように変更

### F) 暫定の sum_stats 風 “set-level evalue” 付与（まだ近似）
- `LOSAT/src/algorithm/tblastx/utils.rs`
  - writer側で HSP を（query/subject/frame）でグループ化し、近傍を“雑に”連結して
    set-level evalue を各HSPへ付与してから `--evalue` で刈る暫定実装を追加
  - **ただし NCBI の even-gap linking / BLAST_SmallGapSumE/BLAST_LargeGapSumE を未移植**

---

## 改善した結果（確認できたこと）
### one-hit（`--window-size 0`）の大幅改善
stop-codon停止 + bucket clear撤廃 + 閾値13 + 後段E-value化で:
- LOSAT one-hit: **58,101 hits**
- 最小 bit score: **22.1**
- 20–29 bit score: **7,808**
- 30–39 bit score: **15,053**

→ 「低スコア帯が0」という致命的状況は解消し、NCBIの挙動に明確に近づいた。

---

## うまくいかなかったこと（残課題）と原因
### 1) 例の低スコアHSP（bit score 24.9）が `--evalue 10` で残らない
NCBIで観測される例:
- 座標: **450063–449974 vs 1135–1224**
- bit score: **24.9**
- frame: **qframe=-3, sframe=1**
- NCBI（sum_stats=T）では evalue ≈ **0.007** で残る

LOSATでは:
- `--evalue 1e9` にすると **そのHSP自体は出力される**（= seed/拡張はできている）
- しかし evalue が **約 3.1e3** と大きく、`--evalue 10` で落ちる

**原因（ほぼ確定）**:
- NCBIの本物の挙動は `link_hsps.c`（even-gap linking）+ `blast_stat.c`（`BLAST_SmallGapSumE/BLAST_LargeGapSumE`）に依存
- 現状LOSATはここを未移植で、暫定の近似（num>1の式/DP/trim/window/2パスの欠落）になっている

### 2) effective search space / evalue の計算経路がNCBIと一致していない可能性
LOSATは単純に `SearchSpace::with_length_adjustment(q_aa_len, s_aa_len)` を使っているが、
NCBIは `blast_setup.c` で `query_info->contexts[].eff_searchsp` を生成し、linkingでそれを使う。

この差も evalue のズレに寄与している可能性がある。

---

## 次のセッションでやるべきこと（優先度順）
### 最優先: NCBIの even-gap linking と sum-statistics 式の移植
1. `ncbi-blast/c++/src/algo/blast/core/link_hsps.c` の **even-gap linking**（tblastx対応）をRustに移植  
   - `s_RevCompareHSPsTbx` 相当の並び・フレーム分離
   - `overlap_size=9`, `gap_size=40` を使った `window_size` と trim の扱い
   - small-gap / large-gap 2パス（`cutoff_small_gap`, `cutoff_big_gap`）のロジック
2. `ncbi-blast/c++/src/algo/blast/core/blast_stat.c` の以下をRustに移植  
   - `BLAST_SmallGapSumE`
   - `BLAST_LargeGapSumE`
   - `BLAST_KarlinPtoE`（= `-log1p(-p)`）
   - `s_BlastSumP`（r<=4のテーブル補間 + r>4のRomberg積分）
   - Romberg積分は `ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:BLAST_RombergIntegrate` を移植
3. その set-level evalue を各HSPに割当てて `--evalue` で刈る（NCBIと同じ段階）

### 検証ゴール
- LOSAT `--window-size 0` で、ターゲットHSPが `--evalue 10` でも残ること
- 20–29 bit score 帯が **44,563** に近づくこと（NCBI one-hit）

---

## 再現コマンド（次セッション用）
### NCBI（one-hit、frames確認）
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
tblastx -query tests/fasta/NZ_CP006932.fasta -subject tests/fasta/NZ_CP006932.fasta \
  -query_gencode 4 -db_gencode 4 -window_size 0 -sum_stats T \
  -outfmt "6 qstart qend sstart send qframe sframe evalue bitscore" \
| rg "450063\t449974\t1135\t1224"
```

### LOSAT（one-hit、診断あり）
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release
LOSAT_DIAGNOSTICS=1 ./target/release/losat tblastx \
  --query tests/fasta/NZ_CP006932.fasta --subject tests/fasta/NZ_CP006932.fasta \
  --query-gencode 4 --db-gencode 4 --window-size 0 \
  --out /tmp/losat_tblastx_windowsize0.out
rg "450063\t449974\t1135\t1224" /tmp/losat_tblastx_windowsize0.out
```

### LOSAT（HSPは生成されているかの確認）
```bash
./target/release/losat tblastx ... --window-size 0 --evalue 1e9 --out /tmp/losat_tblastx_e1e9.out
rg "450063\t449974\t1135\t1224" /tmp/losat_tblastx_e1e9.out
```


