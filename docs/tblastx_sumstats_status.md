### 目的（再確認）
- **LOSAT `tblastx` の hit 分布（特に 20–29 bit）と set-level E-value を、NCBI BLAST+ と極力一致**させる
- 普段の検証は **`--window-size 40`（BLAST+ protein系のデフォルト）**で回し、`--window-size 0` は必要時のみ

---

### これまでにやったこと（実装済み）
- **sum-statistics（NCBI `blast_stat.c` / `ncbi_math.c` 相当）移植**
  - `BLAST_RombergIntegrate`, `s_BlastSumPCalc`, `s_BlastSumP`（テーブル補間+積分）
  - `BLAST_SmallGapSumE`, `BLAST_LargeGapSumE`, `BLAST_KarlinPtoE`, `BLAST_LnFactorial`
  - 実装箇所: `LOSAT/src/stats/sum_statistics.rs`

- **even-gap linking（NCBI `link_hsps.c` 相当）をRust側に導入**
  - small-gap / large-gap の2パス、`gap_size=40`, `overlap_size=9` でDPリンク
  - 実装箇所: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`

- **NCBIの `/4` トリム heuristic を採用**
  - `trim = MIN(len/4, trim_size)` に修正（短いHSPを切り過ぎない）

- **tblastx の grouping を NCBI に寄せる方向で修正**
  - NCBI `s_RevCompareHSPsTbx` 前提（`context/3` + `SIGN(subject.frame)`）に合わせて、
    現状 Rust 側でも **(query strand sign, subject strand sign)** で group するよう変更

- **set-level evalue を各HSPへ付与し、`--evalue` を同段階で適用**
  - 生成HSP → linking → **set-level evalue を上書き** → `--evalue` で刈る
  - 実装箇所: `LOSAT/src/algorithm/tblastx/utils.rs`

- **デバッグ用オプション追加**
  - `--only-qframe/--only-sframe`
  - `--only-qstrand/--only-sstrand`（フレーム束の検証を短縮）
  - 実装箇所: `LOSAT/src/algorithm/tblastx/args.rs`, `LOSAT/src/algorithm/tblastx/utils.rs`

- **`--window-size` デフォルトを 40 に維持**
  - `LOSAT/src/algorithm/tblastx/args.rs` で `default_value_t = 40`

---

### 現状の課題（未解決）
- **追跡HSP（24.9 bits, 450063–449974 vs 1135–1224）が NCBI みたいに evalue が下がらない**
  - LOSAT は出るが **evalue ≈ 3.1e3 のまま**で、`--evalue 10` で落ちる（= linking がまだ効いていない）

- **link_hsps.c の“速度のための最適化”が未移植**
  - NCBI は `LinkHelpStruct`（`next_larger`, `changed`, 事前配列化）＋「前回ベストを再利用」＋早期break を使って **再計算を減らす**
  - 今のRust DPは **再計算が重くなりやすい**（特に `--window-size 0` で顕著）

---

### 今後の方針（優先順位つき）
- **最優先: NCBI `link_hsps.c` の最適化をそのまま移植**
  - `lh_helper` 相当の配列化（`sum/num/xsum`, `q_off_trim/s_off_trim`, `next_larger`, `maxsum1`）
  - `changed` / `use_current_max` / `path_changed` のロジックを入れて **DPの再計算を最小化**
  - これでまず **`--window-size 40` のテストが ~1分以内**に収まる状態を作る

- **次: tblastx用の並べ替えを NCBI と完全一致**
  - `s_RevCompareHSPsTbx` 相当（`context/3`、`SIGN(subject.frame)`、tie-break順）をRustで厳密化

- **次: 長さ・search space の与え方を NCBI に合わせ切る**
  - group（=frame束）内での `query_context` 選び方、`length_adjustment`、`eff_searchsp` の扱いを NCBI と同じにする

- **検証**
  - 普段は **`--window-size 40`** で回して分布比較
  - `--window-size 0` は **フレーム/ストランド制限付きで必要時のみ**
  - 目標: 追跡HSPが `--evalue 10` でも残り、20–29 bit 帯が NCBI の 44,563 に近づく


