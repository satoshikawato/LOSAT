# NEXT SESSION PROMPT — 2025-12-31
## TBLASTX: “LOSATのヒットが微妙に少ない”問題をNCBIパリティで詰める

### ゴール
- 直近セッションで修正した「過剰に長い高一致HSP」の再発を防ぎつつ、
  **BLAST+（NCBI）と比べて LOSAT のヒットが少ない**差分を、根拠ベースで潰す。

### 重要な前提（すでに確定した事実）
- sum-statistics linking は **座標を結合しない**（E-valueを付与するだけ）。
- SEGに関しては NCBI と同様に、
  - workingクエリ配列を **Xで上書き**（seed/extensionに効かせる）
  - identity/mismatch は **nomask** で計算
  を LOSAT に実装済み。

### まずやること（再現・定量）
1. **差分が本当に「ヒット数」なのかを定量**する
   - 比較対象:
     - `LOSAT/tests/blast_out/*.tblastx.n1.out`（BLAST+）
     - `LOSAT/tests/losat_out/*.tlosatx.n1.out`（LOSAT）
   - 例: NZ_CP006932 self
     - BLAST+ 行数: outfmt7のコメント行除外後（`comment='#'`で読む）
     - LOSAT 行数: 同じ列定義で読む
2. “どの帯域で足りないか”を特定
   - identity bins（例: 20–30, 30–40, …）
   - length bins（log bin）
   - ここで **不足が集中するゾーン**を1つ決める（例: identity 30–60, length 30–150 aa など）

### 次にやること（代表的な欠損HSPを1つ選ぶ）
3. BLAST+にあってLOSATにないHSPを1つ選び、**近傍のLOSAT出力を探す**
   - 近傍とは:
     - 同じ(sseqid,qseqid)、同じ(フレーム組み合わせが推定できるなら)、
     - qstart/qend, sstart/send が近いもの
   - 「完全欠損」なのか、「あるが短い/スコアが低い/閾値で落ちた」のかを分類する

### “落ちた場所”を特定する（LOSAT側のパイプライン）
4. LOSATのtblastxパイプライン（現状）を段階別にチェックする
   - **lookup生成**: `LOSAT/src/algorithm/tblastx/lookup.rs` `build_ncbi_lookup()`
   - **subject走査**: `LOSAT/src/algorithm/tblastx/utils.rs` `s_blast_aa_scan_subject()` 呼び出し周辺
   - **two-hit ungapped extension**: `LOSAT/src/algorithm/tblastx/extension.rs` `extend_hit_two_hit()`
   - **cutoff判定**: `utils.rs` の `score >= cutoff` の条件
   - **sum-stats linking**: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`
   - **最終evalueフィルタ**: `h.e_value > evalue_threshold`（CLI `--evalue`）

### NCBIパリティ観点で疑うべき“具体項目”（仮説 → 要検証）
以下は**推測で直さない**。必ず「NCBI側コード/出力」と突き合わせて検証してから修正する。

#### A) スコア行列（特に stop `*`）がNCBIと一致しているか
- NCBI `sm_blosum62.c` では `*-* = +1` が含まれる。
- LOSAT 側 `LOSAT/src/utils/matrix.rs` に **NCBIと異なる改変が残っていないか**確認する。
- もし改変が残っているなら、まず「NCBIの値に戻したときに overlong tail が再発するか」を検証する。

#### B) lookupアルファベット / seedでの“非標準残基”扱い
- LOSAT は lookup を `0..=23`（stop除外）に絞っている（コメントあり）。
- NCBI `blast_aalookup.c` は `BLASTAA_SIZE` を使う（lookupのalphabetサイズ・圧縮alphabet含む）。
- 欠損HSPが「stopを含むseed近傍」由来で落ちている可能性があるため、
  - NCBIが実際に seed生成で `*` をどう扱うか
  - LOSATの除外がヒット数に効いているか
  を確認する。

#### C) high-frequency word suppression
- LOSAT: `MAX_HITS_PER_KMER` を超えるwordは捨てる（lookup構築時）。
- NCBIも類似の抑制を持つが、閾値や適用タイミングが異なる可能性がある。
- “ヒット数が少ない”が、特定のリピート/低複雑領域で顕著ならここを疑う。

#### D) cutoff（two-hit extension後）と丸め/閾値の差
- LOSATは evalue→raw score換算や gap decay を入れている。
- NCBI `blast_parameters.c` / `aa_ungapped.c` の該当部と式・rounding（floor/ceil/+1）を照合する。

### 実験手順（次セッションの作業順）
1. まず “不足帯域” を定量して、代表欠損HSPを1つ選ぶ
2. そのHSPが LOSAT のどの段階で落ちたかを特定（seed? extension? cutoff?）
3. NCBI側コード（該当関数）をピンポイントで参照し、差分があれば修正
4. `LOSAT/tests/run_comparison.sh` → `plot_comparison.py` で再確認
5. **overlong tail が再発していない**ことを必ず確認（最優先の安全柵）

### 便利コマンド（再現）
```bash
cd LOSAT/LOSAT
cargo build --release
cd tests
bash run_comparison.sh
python3 plot_comparison.py
python3 plot_overall_trend.py
```

### 関連資料
- 直近まとめ: `LOSAT/docs/SESSION_SUMMARY_2025-12-31_TBLASTX_QUERY_MASKING.md`

