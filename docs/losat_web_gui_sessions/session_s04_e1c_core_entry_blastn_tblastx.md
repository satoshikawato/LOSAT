# Session S04 — E1c：核の入口（BLASTN と TBLASTX）

## INSTRUCTION PROMPT

LOSAT の段階 E1c を実行する。S02・S03 と同じく、**出力を変えない**まとめ直しである。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S04 の行である。この段階では outfmt 0 は作らない（BLASTN は 6/7、TBLASTX は 6 のまま）。outfmt 0 は S06〜S08 で移植する。

1. 基準：S02 の最初に取った BLASTN と TBLASTX の基準（`docs/evidence/losat_web_e1a/baseline/` と Gate A の `LOSAT/tests/platform_native_v010_canonical.tsv`）を使う。性能の基準として、各 program の 1 つの fixture で、ネイティブと command-WASI の時間を 1 回の暖機と 3 回の計測で取る。
2. BLASTN：メモリ上の入口は wasm32 のときだけコンパイルされる（`LOSAT/src/algorithm/blastn/blast_engine/run.rs:4349`）。`run_local` を通すようにし、CLI の `run` と v1 の `run_web_pair` をその薄い層にする。
3. TBLASTX：出力を書く箇所が少なくとも 4 つあり、そのうち 1 つは writer thread（`LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1352`）、1 つは `wasm-threads` のときだけコンパイルされる（同 `:3061-3064`）。書く箇所を 1 つにまとめ、`run_local` の出力の配布につなぐ。NCBI の順序に集約した後で書く、という今の順序は変えない。
4. 両 program とも、`hits` は S07・S08 で `PairwiseHit` を作るまで空にしておき、観測者は outfmt 6/7 の行についてだけつなぐ。`docs/web/verification_cells.tsv` の両 program の行を埋める。
5. 試験：第 1 項の基準との一致（Gate A のハッシュを含む）。v1 の serial と threaded の reactor を既存の方法（`LOSAT/tests/build_wasi_artifacts.py`）で作り直し、v1 の reactor の検査（`check_wasi_reactor.js`、`check_wasi_api_limits.js`）が通ること。TBLASTX の `wasm-threads` 専用の出力箇所を通る threaded の実行で、serial と同じ出力になること。性能の中央値が基準の +5% 以内。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
6. 独立監査を受ける。

完了条件は計画 §7 の S04 の行による。記録は `docs/evidence/losat_web_e1c/README.md`。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S05 — アダプタと ABI v2](session_s05_e1d_adapter_abi_v2.md)。
