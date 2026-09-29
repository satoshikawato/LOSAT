# Session S08 — E2b：TBLASTX outfmt 0 と 7

## INSTRUCTION PROMPT

LOSAT の段階 E2b を実行する。TBLASTX（TLOSATX）に outfmt 0 と 7 を移植する。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S08 の行である。1 セッションで終わらなければ、S06・S07 と同じく「権威と fixture」（S08a）と「実装とゲート」（S08b）に分け、README の表を直す。

現状：TBLASTX は outfmt 6 だけを受け付ける（`LOSAT/src/blastinput/value_parsers.rs:313-320`）。最終の結果は平らな `Vec<Hit>` で、subject frame・整列文字列を持たない（`LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs` と `LOSAT/src/common.rs`）。承認済みの例外がある：局所 `-subject` の検索では、既定以外の `-db_gencode` を subject の翻訳・検索・表示に使う（ルートの `AGENTS.md` の「Mandatory Compliance Requirements」の 3）。

1. NCBI の経路を記録する：`c++/src/app/blast/tblastx_app.cpp` から `CBlastFormat` と `c++/src/objtools/align_format/`（`showalign.cpp`、`showdefline.cpp`、表形式の注釈は `tabular.cpp`）へ。翻訳検索に特有の表示（`Frame = +1/-2` のような両側の frame、ギャップの無いアラインメント、positives、翻訳した配列の行、核酸座標）と、outfmt 7 の query ごとの見出しと末尾を、ファイルと行の番号付きで表にする。
2. fixture を選ぶ（`LOSAT/tests/tblastx_v010_parity_manifest.tsv` の入力を基に、6 frame、複数 HSP、複数の query と subject、ヒット無し、遺伝暗号 1 と既定以外のコード）。NCBI BLAST+ 2.17.0 の outfmt 0 と 7 を oracle として固定する。既定以外の subject の遺伝暗号の場合は、承認済みの例外のとおりに差を分類する（ほかの差は認めない）。
3. 最終の HSP 一覧から `PairwiseHit` を作るのに必要な情報（subject frame、翻訳した配列の区間）を、NCBI と同じ時点で保持する。計算の順序は変えない。`PairwiseHit` を `hits` と outfmt 0 の formatter の両方に渡し、観測者を outfmt 0 と 7 につなぐ。
4. outfmt 0 と 7 を移植し、`tblastx_outfmt` が受け付けるようにする。Rust の移植箇所の直上に NCBI の参照を書く。S07 の部品を使う。
5. 試験：fixture すべてで NCBI とバイト一致（例外の分類を除く）。既存の outfmt 6 の回帰ゲート（`LOSAT/tests/audit_tblastx_v010.py` と Gate A のハッシュ）が変わらない。`docs/web/verification_cells.tsv` の TBLASTX の outfmt 0/7 の升目を埋め、TBLASTX の全升目（0/6/7 × スレッド 1/2/4）の V-ABI が通ること。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
6. アダプタ：TBLASTX の形式に 0 と 7 を足す（`web/adapter/src/run.rs` の `Program::formats`、`tests/v_abi.js` と `tools/v_abi_cases.py` の `FORMATS`）。これで全 program が既定の `-outfmt 0` を受け付けるので、`run.rs` の `parse` が挿入している `-outfmt 6` をやめ、`validate` の文言を CLI と同じにする（S05 の独立監査の指摘。`docs/web/abi_v2.md` §4 を直す）。
7. 独立監査を受ける。

完了条件は計画 §7 の S08 の行による。記録は `docs/evidence/losat_web_e2b/README.md`。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S09 — ブラウザでの実行基盤](session_s09_w1_browser_runtime.md)。
