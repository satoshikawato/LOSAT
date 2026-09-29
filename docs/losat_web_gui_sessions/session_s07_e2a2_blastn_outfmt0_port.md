# Session S07 — E2a-2：BLASTN outfmt 0（実装とゲート）

## INSTRUCTION PROMPT

LOSAT の段階 E2a-2 を実行する。S06 で固定した権威と fixture に従って、BLASTN（LOSATN）に outfmt 0 を移植する。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S07 の行である。入力は `docs/evidence/losat_web_e2a/`（S06 の成果物）である。

1. S06 の経路の表の順に移植する。Rust の移植箇所の直上に、NCBI のファイル・行と C/C++ の断片を書く。`LOSAT/src/report/pairwise.rs` の既存の部品を使い、同じ働きの writer を二重に作らない。
2. BLASTN の最終の HSP 一覧（並べ終えたもの）から `PairwiseHit` を作る。整列文字列は、traceback の編集操作（`gap_info`）から作る。この一覧を `run_local` の `hits` と、outfmt 0 の formatter の両方に渡す（計画 §4.4）。観測者を outfmt 0 の節にもつなぐ。
3. `parse_blastn_output_format` が 0 を受け付けるようにする。未対応の組合せ（S06 の一覧）は明示的に拒否する。
4. 試験：S06 の fixture すべてで、NCBI とバイト一致すること。最初に違ったところは、query / subject の ID、鎖、座標、score、E 値、HSP の順序、表示の文字列に分けて記録し、該当するものをまとめて直す。既存の BLASTN の outfmt 6/7 の回帰ゲートと Gate A のハッシュが変わらないこと。`docs/web/verification_cells.tsv` の BLASTN の outfmt 0 の升目を埋め、BLASTN の全升目（0/6/7 × スレッド 1/2/4）の V-ABI が通ること。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
5. 独立監査を受ける。

完了条件は計画 §7 の S07 の行による。記録は `docs/evidence/losat_web_e2a/README.md`。

### S05 の独立監査の残り（アダプタ。最初に行う）

S07 は reactor を作り直して full の V-ABI を再実行するので、S05 の独立監査の指摘のうち reactor のバイトを変えるものをここで直す（記録は `docs/evidence/losat_web_e1d/README.md` の「独立監査」）。

1. 計画 TD-11：`web/adapter/tools/build_reactors.py` で、checkout のパスと `CARGO_HOME` を `--remap-path-prefix` で置き換え、`tools/check_build_identity.py` はこの置き換えだけを rustflags の差として認める。別の場所に展開した同じ commit から同じバイトの reactor ができることを確かめ、`docs/web/abi_v2.md` §2 を直す。
2. `validate`：解析に渡す binary の名前を CLI と同じ `LOSAT` にする。`-num_threads` を `-out`・`-outfmt` と同じく拒否する（host が付ける）。`tests/v_abi.js` の `checkSurface` に、必須の引数が無い場合と `-num_threads` の場合を足す。
3. `register`：記録の照合を関数に分け、ID・長さ・記録の数が食い違う場合に失敗することの単体試験を足す（S05 の指示書 4）。
4. `RangeRecorder`：対応の無い `hsp_end` と重複した番号を、release のビルドでもエラーにする（今は `debug_assert` だけ）。
5. `tools/run_v_abi_parallel.py`：部分を実行する前に、その部分の古い結果を消す。`summary.json` に 3 つの成果物の SHA-256 を記録する。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S08 — TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md)。BLASTN で作った部品のうち TBLASTX に使えるものを S08 の指示書に書き足す。
