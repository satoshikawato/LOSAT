# Session S03 — E1b：核の入口（TBLASTN）

## INSTRUCTION PROMPT

LOSAT の段階 E1b を実行する。S02 と同じく、**出力を変えない**まとめ直しである。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S03 の行である。S02 のゲート記録（`docs/evidence/losat_web_e1a/README.md`）にある `run_local` と `ReportOutputs` の形、基準のファイルに合わせる。BLASTX には手を付けない（計画 DW-10）。

1. 基準：S02 の最初に取った TBLASTN の基準（TLOSAN 計画の Stage G の凍結バイト。再現手順は `docs/release/v0.2.0.md` の「Reproduction」と `docs/evidence/tlosan_stage_g/`）を使う。性能の基準として、TBLASTN の 1 つの fixture で、ネイティブと command-WASI の時間を 1 回の暖機と 3 回の計測で取る。
2. TBLASTN は、今はファイルを読む `TblastnArgs::run(self)`（`LOSAT/src/algorithm/tblastn/args.rs:382`）しか入口が無く、警告を stderr に直接書き、エラーのときに途中までの出力を残さないように、出力の全体をバッファしてから書く（同 `:525-566`）。CLI ではこの性質を保つ（Web では失敗した実行の出力を破棄するので、同じ結果になる）。`run_local` を通すようにし、CLI の `run` をその薄い層にする。20,000 残基の query バッチの扱い（`args.rs:468-472`）と、形式の分岐（`LOSAT/src/algorithm/tblastn/stage_e_report.rs` の `render`）の順序は変えない。警告は `diagnostics` に書く。
3. `hits` と観測者を S02 と同じ方法でつなぐ。TBLASTN の `PairwiseHit` は今 outfmt 0 のときだけ作られる（`stage_e_report.rs`）。最終の結果から 1 回だけ作り、`hits` と outfmt 0 の formatter の両方に渡す。`docs/web/verification_cells.tsv` の TBLASTN の行を埋める。
4. 試験：Stage G の凍結バイトとの一致（27 の遺伝暗号の升目を含む）。v1 の reactor の検査。3 形式の同時出力と単一形式の CLI との一致。観測者の範囲の一致。警告が 1 回だけ出ること。性能の中央値が基準の +5% 以内。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
5. 独立監査を受ける。

記録は `docs/evidence/losat_web_e1b/README.md`。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S04 — 核の入口：BLASTN と TBLASTX](session_s04_e1c_core_entry_blastn_tblastx.md)。
