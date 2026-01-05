---
name: sum_stats_linking speedup
overview: Optimize `sum_stats_linking` hot loops to reduce allocation and bounds-check overhead while keeping NCBI control flow, validating each step with the user-provided AP027280 self-run benchmark.
todos:
  - id: baseline-ap027280
    content: 現状 `../target/release/LOSAT` で AP027280 self をユーザー指定コマンドで実行し、log の real/user/sys と hit 数を記録する
    status: pending
  - id: remove-hot-alloc-chainmembers
    content: "`sum_stats_linking.rs` のチェーン除去ループでデバッグ用 `chain_members` 収集を通常実行から排除し、alloc/push をゼロ化する"
    status: pending
    dependencies:
      - baseline-ap027280
  - id: stream-frame-groups
    content: "`apply_sum_stats_even_gap_linking` の frame_groups 二段階処理をストリーミング化して不要な Vec<Vec<_>> を廃止する"
    status: pending
    dependencies:
      - remove-hot-alloc-chainmembers
  - id: reuse-scratch-buffers
    content: "`hsp_links` / `lh_helpers` を scratch 化して capacity 再利用し、グループ毎の大きな初期化/再確保を減らす"
    status: pending
    dependencies:
      - stream-frame-groups
  - id: tighten-inner-loops
    content: INDEX0/INDEX1 内側ループの重複ロードと bounds check を段階的に削減（必要なら unsafe + debug_assert）
    status: pending
    dependencies:
      - reuse-scratch-buffers
---

# sum_stats_linking 高速化（NCBI準拠のまま Rust オーバーヘッド削減）

## ゴール
- `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` の **`apply_sum_stats_even_gap_linking` / `link_hsp_group_ncbi`** を中心に、NCBI `link_hsps.c` の制御フローを維持したまま、Rust 側の **不要な alloc / env 参照 / bounds check** を削って実行時間を落とします。
- 各ステップ後に、ユーザー指定のベンチコマンドで **必ず計測**します。

## ベンチ手順（各ステップ共通）
- `cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests`
- `LOSAT_BIN="../target/release/LOSAT"` を前提に以下を実行:

```bash
(time $LOSAT_BIN tblastx -q ./fasta/AP027280.fasta -s ./fasta/AP027280.fasta -o ./losat_out/AP027280.AP027280.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 --seg )&>./losat_out/AP027280.AP027280.tlosatx.n1.log
```

## 進め方（変更は `sum_stats_linking` に集中）
### Step 0: ベースライン確定
- 現状の `real/user/sys` を `./losat_out/AP027280.AP027280.tlosatx.n1.log` から控える。
- 参考として NCBI 出力 `tests/blast_out/AP027280.AP027280.tblastx.n1.out` の hit 数（ヘッダ行除外）と、LOSAT 出力の hit 数を控える（大きくズレていないことの確認）。

### Step 1: ホットループ内の **不要 alloc をゼロ**（出力非依存）
- `link_hsp_group_ncbi` の「チェーン抽出」部分で、デバッグ用途の `chain_members: Vec<_>` を **デバッグ有効時のみ**収集するように変更（通常実行では完全に skip）。
  - 現状はデバッグ無効でも `Vec` を毎チェーン確保・push しており、HSP 数が多いケースで致命的なオーバーヘッドになり得る。
  - NCBI 参考: `link_hsps.c` のチェーン除去部（`best[ordering_method]` から `linked_to=-1000` にしつつ unlink する箇所）に合わせ、Rust 側は観測用データ構築を外に逃がす。

### Step 2: グループ処理の alloc 削減（出力非依存）
- `apply_sum_stats_even_gap_linking` の `frame_groups: Vec<Vec<UngappedHit>>` を作る二段階処理をやめ、
  - グループ境界検出と同時に `link_hsp_group_ncbi(...)` を呼ぶ **ストリーミング処理**に変更して、グループ `Vec` を溜めない。

### Step 3: per-group / per-pass の再確保を抑制（出力非依存）
- `hsp_links` と `lh_helpers` をグループ毎に `vec![..; n]` で初期化している箇所を、
  - `reserve_exact` + `set_len`（または `resize`）で **capacity を使い回す scratch** に寄せ、
  - 「毎回ゼロ初期化が必要な領域」だけを必要範囲で上書きする。
- NCBI 側が `lh_helper` を再利用する意図（コメント済み）に合わせる。

### Step 4: インナーループの bounds check/重複ロード削減（出力同等）
- INDEX0/INDEX1 の内側ループで `hsp_links[i]` / `lh_helpers[idx]` を何度も index している部分を、
  - 事前に参照/値をローカルへ引き出す、
  - さらに必要なら `get_unchecked` / raw pointer で **境界チェックを除去**（不変条件は `debug_assert!` で担保）
  - という順で段階的に詰める。

### Step 5（必要なら）: `sum_statistics` の計算コスト見直し
- まだ `sum_stats_linking` が支配的なら、`small_gap_sum_e` / `large_gap_sum_e` の内部（`ln_factorial_int` 等）を NCBI 実装（`blast_stat.c` / `ncbi_math.c`）と突き合わせ、
  - **同じ数値経路**を保ったまま（or 近づけつつ）高速化する。

## 主要変更ファイル
- [`LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`](LOSAT/src/algorithm/tblastx/sum_stats_linking.rs)
- （必要なら）[`LOSAT/src/stats/sum_statistics.rs`](LOSAT/src/stats/sum_statistics.rs)
- 参照（NCBI）: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c`

## 完了条件
- 各ステップで上記ベンチを回し、`real` が単調改善（または少なくとも悪化しない）こと。
- 生成される `AP027280.AP027280.tlosatx.n1.out` が極端に崩れない（hit 数や分布が NCBI 側 `tests/blast_out/AP027280.AP027280.tblastx.n1.out` から大きく逸脱しない）こと。