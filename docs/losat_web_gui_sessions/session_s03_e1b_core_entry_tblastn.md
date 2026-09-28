# Session S03 — E1b：核の入口（TBLASTN）

## INSTRUCTION PROMPT

LOSAT の段階 E1b を実行する。S02 と同じく、**出力を変えない**まとめ直しである。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S03 の行である。S02 のゲート記録（`docs/evidence/losat_web_e1a/README.md`）にある `run_local` と `ReportOutputs` の形、基準のファイルに合わせる。BLASTX には手を付けない（計画 DW-10）。

1. 基準：S02 で取った基準 `docs/evidence/losat_web_e1a/baseline/hashes.tsv` を使う（236 件。TBLASTN の 162 件は Stage G の凍結値とも一致している）。基準は取り直さない。変更の後は、`docs/evidence/losat_web_e1a/capture_outputs.py run` を全 program で実行し（共有のコードを変えたときの影響も見るため）、`compare` で基準と比べる（差 0 が条件。`Sakai.MG1655.megablast` の Gate A との不一致は変更前からある既知の差で、基準と同じならよい）。変更前の成果物（ネイティブと serial / threaded の command-WASI）は、S02 の後のコミットの木から作る。性能は、変更前と変更後の成果物を `docs/evidence/losat_web_e1a/measure_perf.py run --before … --after …` で測る（1 回の暖機と 3 回の計測を、変更前と変更後を 1 回ずつ交互に行う。S02 で、時間を分けて測ると機械の状態の変化で同じ成果物の中央値が 20% 以上動くことが分かったため）。TBLASTN の case は `tblastn` と `tblastn-fmt0`（Stage G の benchmark の fixture）。
2. TBLASTN は、今はファイルを読む `TblastnArgs::run(self)`（`LOSAT/src/algorithm/tblastn/args.rs:382`）しか入口が無く、警告を stderr に直接書き、エラーのときに途中までの出力を残さないように、出力の全体をバッファしてから書く（同 `:525-566`）。CLI ではこの性質を保つ（Web では失敗した実行の出力を破棄するので、同じ結果になる）。`run_local` を通すようにし、CLI の `run` をその薄い層にする。20,000 残基の query バッチの扱い（`args.rs:468-472`）と、形式の分岐（`LOSAT/src/algorithm/tblastn/stage_e_report.rs` の `render`）の順序は変えない。警告は `diagnostics` に書く。
3. `hits` と観測者を S02 と同じ方法でつなぐ（S02 で確定した型は計画 §4.2、BLASTP の実装は `LOSAT/src/algorithm/blastp/blast_engine.rs` の `run_local` と出力の分岐、試験の書き方は `LOSAT/tests/run_local_blastp.rs`）。観測者があるときだけ、HSP の区切りで formatter のバッファを流す。TBLASTN の `PairwiseHit` は今 outfmt 0 のときだけ作られる（`stage_e_report.rs`）。最終の結果から 1 回だけ作り、`hits` と outfmt 0 の formatter の両方に渡す。`docs/web/verification_cells.tsv` の TBLASTN の行を埋める。
4. 試験：Stage G の凍結バイトとの一致（27 の遺伝暗号の升目を含む）。v1 の reactor の検査。3 形式の同時出力と単一形式の CLI との一致。観測者の範囲の一致。警告が 1 回だけ出ること。性能の中央値が基準の +5% 以内。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
5. 独立監査を受ける。

記録は `docs/evidence/losat_web_e1b/README.md`。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S04 — 核の入口：BLASTN と TBLASTX](session_s04_e1c_core_entry_blastn_tblastx.md)。
