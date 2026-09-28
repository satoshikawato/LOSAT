# Session S07 — E2a-2：BLASTN outfmt 0（実装とゲート）

## INSTRUCTION PROMPT

LOSAT の段階 E2a-2 を実行する。S06 で固定した権威と fixture に従って、BLASTN（LOSATN）に outfmt 0 を移植する。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S07 の行である。入力は `docs/evidence/losat_web_e2a/`（S06 の成果物）である。

1. S06 の経路の表の順に移植する。Rust の移植箇所の直上に、NCBI のファイル・行と C/C++ の断片を書く。`LOSAT/src/report/pairwise.rs` の既存の部品を使い、同じ働きの writer を二重に作らない。
2. BLASTN の最終の HSP 一覧（並べ終えたもの）から `PairwiseHit` を作る。整列文字列は、traceback の編集操作（`gap_info`）から作る。この一覧を `run_local` の `hits` と、outfmt 0 の formatter の両方に渡す（計画 §4.4）。観測者を outfmt 0 の節にもつなぐ。
3. `parse_blastn_output_format` が 0 を受け付けるようにする。未対応の組合せ（S06 の一覧）は明示的に拒否する。
4. 試験：S06 の fixture すべてで、NCBI とバイト一致すること。最初に違ったところは、query / subject の ID、鎖、座標、score、E 値、HSP の順序、表示の文字列に分けて記録し、該当するものをまとめて直す。既存の BLASTN の outfmt 6/7 の回帰ゲートと Gate A のハッシュが変わらないこと。`docs/web/verification_cells.tsv` の BLASTN の outfmt 0 の升目を埋め、BLASTN の全升目（0/6/7 × スレッド 1/2/4）の V-ABI が通ること。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
5. 独立監査を受ける。

完了条件は計画 §7 の S07 の行による。記録は `docs/evidence/losat_web_e2a/README.md`。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S08 — TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md)。BLASTN で作った部品のうち TBLASTX に使えるものを S08 の指示書に書き足す。
