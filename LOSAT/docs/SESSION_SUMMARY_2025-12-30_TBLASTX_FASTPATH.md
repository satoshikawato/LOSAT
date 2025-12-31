## LOSAT TBLASTX: セッションまとめ（2025-12-30）

### 目的（このセッションのゴール）
- **657kb self-comparison がタイムアウトする性能劣化を解消**し、NCBI BLAST に近いホットパスへ移植する
- `blast_aascan.c` / `aa_ungapped.c` の scan→two-hit→ungapped extension の流れを **Rust に忠実移植**する

---

### 1) このセッションで実施したこと（実装変更）

#### **TBLASTX ルックアップ/スキャン/延長のホットパスを NCBI 方式へ**
- **Lookup index を NCBI の bit-shift + mask 方式に変更**
  - 参照: `blast_lookup.h` の `ComputeTableIndex` / `ComputeTableIndexIncremental`
  - 変更先: `src/algorithm/tblastx/lookup.rs`
  - **backbone/pv は “backbone_size（max index + 1）” を確保**する方式へ
  - 現状は stop-codon（24）を lookup alphabet から除外しているため **alphabet_size=24（0..=23）**

- **Subject scanning を NCBI `s_BlastAaScanSubject` 構造へ寄せて再実装**
  - 参照: `ncbi-blast/.../blast_aascan.c:48-131`
  - 変更先: `src/algorithm/tblastx/utils.rs` の `s_blast_aa_scan_subject`
  - `PV_TEST` + incremental index を中心に実装（stop など invalid residue は rolling index を reset）

- **Two-hit extension を NCBI `s_BlastAaExtendTwoHit` 相当に統一**
  - 参照: `ncbi-blast/.../aa_ungapped.c:1089-1158`
  - `utils.rs` 側の簡略実装を撤去し、`extension.rs` の `extend_hit_two_hit()` を使用

- **diag 配列/offset_pairs をスレッド単位で再利用**
  - `rayon::par_iter().for_each_with` → `for_each_init`
  - `diag_array` / `offset_pairs` を worker thread で一度だけ確保して使い回し

#### **inner loop の重い処理（String clone / identity / bit/evalue）を排除**
- scan/extension の内側では `UngappedHit`（数値ID + 座標 + raw_score）を積むだけに変更
- `sum_stats_linking.rs` の grouping を **(q_idx, s_idx, strand)** に変更し、`String::clone()` を除去
- sum-stats linking 後に **最終出力する Hit のみに**:
  - identity/mismatch 計算
  - bit_score 計算
  - `query_id/subject_id` の `String` 生成

#### **関連ファイル（重要）**
- `src/algorithm/tblastx/utils.rs`
- `src/algorithm/tblastx/lookup.rs`
- `src/algorithm/tblastx/extension.rs`
- `src/algorithm/tblastx/translation.rs`
- `src/algorithm/tblastx/sum_stats_linking.rs`
- `src/algorithm/tblastx/chaining.rs`（`UngappedHit` 追加）

---

### 2) 実験・再現手順（比較/プロット）

#### **出力ディレクトリ**
- **LOSAT output**:
  - Windows: `C:\\Users\\kawato\\Documents\\GitHub\\LOSAT\\LOSAT\\tests\\losat_out`
  - WSL: `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests/losat_out`
- **BLAST output**:
  - Windows: `C:\\Users\\kawato\\Documents\\GitHub\\LOSAT\\LOSAT\\tests\\blast_out`
  - WSL: `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests/blast_out`

#### **比較実行スクリプト**
- `tests/run_comparison.sh`
  - TBLASTX（LOSAT）: `--query-gencode 4 --db-gencode 4 --seg -n 20`
  - BLAST+ パートはスクリプト内で comment-out されている（既存の `blast_out` を比較に使用）

#### **プロット生成**
- `python3 tests/plot_comparison.py` → `tests/plots/compare_*_TBLASTX.png` など
- `python3 tests/plot_overall_trend.py` → `tests/plots/overall_trend_comparison.png`
- `python3 tests/plot_execution_time.py` → `tests/plots/execution_time_comparison_all.png`

---

### 3) 現状の結果（観察された差分）

#### **性能**
- `tests/plots/execution_time_comparison_all.png` より、TBLASTX は **BLAST+ より大幅に短時間**で完走している（少なくとも今回の比較タスクでは）。

#### **品質（優先課題）**
- **高 identity の比較で、BLAST+ に比べて “長すぎる hit（length が極端に大きい）” が生成される**
  - 例:
    - `tests/plots/compare_AP027132_vs_NZ_CP006932_TBLASTX.png`
    - `tests/plots/compare_NZ_CP006932_Self_TBLASTX.png`
- **低 identity ヒットが若干少ない**
  - 例:
    - `tests/plots/compare_AP027078_vs_AP027131_TBLASTX.png`
    - `tests/plots/compare_AP027131_vs_AP027133_TBLASTX.png`

---

### 4) 次の優先課題と方針

#### **最優先: 高identityで“長すぎる hit”の原因究明と是正**
- **仮説**: NCBI は TBLASTX の翻訳配列を **stop codon（ORF境界）で分割（seq_ranges/fence）**しており、HSP が stop を跨いで伸びない。\n
  - LOSAT は stop codon を **“通常のAA（24）”として保持**しているため、X-drop をすり抜けて stop を跨ぐ超長尺 ungapped HSP が生成され得る。\n
- **方針**:
  - NCBI 側の翻訳/seq_ranges の実装を確認し、**stop を “境界（sentinel）” として扱う**（または seq_ranges を生成して scan/extension を分割）\n
  - まずは `NZ_CP006932_Self` の LOSAT 出力から最長HSPを1つ抽出し、\n
    - その HSP 範囲内に stop が多く含まれるか\n
    - BLAST+ の対応領域ではどのように分割されているか\n
    を確認して仮説を確定する\n

#### **次点: 低identityヒットの不足**
- 高identity問題を潰した後に、\n
  - seed/neighbor 生成の閾値\n
  - masking（SEG/DUST）適用範囲\n
  - cutoffs（E-value cutoff）\n
  の差分を点検する。

---

### 5) 備考
- `cargo test` は TBLASTX 変更部分のビルドは通るが、現状 `dust` / `packed_nucleotide` の既存テストが3件失敗している（今回の TBLASTX 変更とは別系統の可能性が高い）。\n
  ただし次セッションの主目的は **TBLASTX の品質差分修正**なので、必要になった時点で切り分ける。


