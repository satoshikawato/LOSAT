# Session S11 — E2d：`-query_loc` / `-subject_loc`

## INSTRUCTION PROMPT

LOSAT の段階 E2d を実行する。これはエンジン（`LOSAT/`）の NCBI パリティ作業である。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S11 の行である。

目的：設計書 §2.1 の「領域は座標入力とプレビュー選択」を、NCBI の範囲指定の意味どおりに、BLASTN・BLASTP・TBLASTN・TBLASTX で使えるようにする（BLASTX は SX。計画 DW-10、DW-11）。範囲は、その役割のレコードが 1 つのときだけ画面で指定できる（計画 DW-9）が、CLI と同じく、エンジンは複数のレコードでも NCBI の意味どおりに扱う。UI で配列を切り出して検索することはしない。

1. 現状を確かめる：S08+（E2e）の後、4 つの program とも `-query_loc` / `-subject_loc` を `LOSAT/src/cli.rs` の `is_unported_blastn_arg`・`is_unported_blastp_arg`・`is_unported_tblastn_arg`・`is_unported_tblastx_arg`（2026-10-05 の時点で 215・283・493・551 行）が明示的に拒否する（文言は「… is not supported by LOSAT's <PROGRAM>」、ABI v2 の `validate` も同じ）。BLASTX の `is_unported_blastx_arg` は SX の範囲で、変えない。
2. NCBI の経路を追う：`c++/src/algo/blast/blastinput/blast_args.cpp` の `kArgQueryLocation`（定義は 1946 行付近、解析は 1996〜1997 行）と `kArgSubjectLocation`（2373 行付近）から、範囲が配列・統計・座標の表示にどの時点で効くかを program ごとに記録する。範囲は入力元が読むすべてのレコードに適用される（`c++/src/algo/blast/blastinput/blast_fasta_input.cpp:433-459`。逆向きの範囲はエラー、長さを超える終点は黙って切り詰め）。
3. 比較用の fixture を決める（4 program、範囲の端、範囲の内外にまたがるヒット、minus 鎖、翻訳の frame、範囲を全レコードに適用する複数レコードの入力、範囲の始点より短いレコード（NCBI ではエラー、`blast_fasta_input.cpp:446-452`））。NCBI BLAST+ の出力を oracle として固定し、SHA-256 を記録する。
4. 4 program に移植する。移植した箇所の直上に、NCBI のファイル・行と断片を書く。NCBI に無い組合せは明示的に拒否する。`run_local` と ABI v2 の `validate` から使えるようにする。
5. すべての fixture で outfmt 0/6/7 が NCBI とバイト一致すること、既存の認証済みの出力に退行が無いことを確かめ、`docs/web/verification_cells.tsv` に升目を足す。独立監査を受ける。

完了条件は計画 §7 の S11 の行による。記録は `docs/evidence/losat_web_e2d/README.md`。1 セッションで終わらなければ S11b に分ける。

## E2e（S08+・S08+a・S08+b）からの引き継ぎ（2026-10-05）

- E2e の保守者の判断（D11〜D15、web ABI v1 の BLASTP の誤りの文言）は 2026-10-05 に推奨の案で決まり、`PD-LOSAT-NCBI-DEFECTS` 版 1.3・`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 版 1.4・`AGENTS.md` に反映した（計画 DW-19）。E2e から残る完了条件は無い。範囲の指定で同じ種類の NCBI の不具合に当たったら、同じ方針（`PD-LOSAT-NCBI-DEFECTS` の規則）で扱う。
- 範囲と query の分割：BLASTP・TBLASTN は NCBI と同じ chunk（blastp 10000・tblastn 20000 残基、重なり 100、`split_query_aux_priv.cpp`）で query を分ける（S08+a、判断 D13）。範囲を指定した query がいつ分割されるか（範囲の長さか、レコードの長さか）を NCBI の経路で確かめ、境界をまたぐ範囲を fixture に入れる。TBLASTX は NCBI が 10002 塩基で分けるが出力は変わらず、LOSAT は分けない（E2e の残件 R2c-1、記憶だけ）。
- option の値の検査は 3 つの program とも NCBI の app の層（`blastinput/app.rs` の `check_options`、NCBI の順）にある。`-query_loc` / `-subject_loc` の文法の検査（`CArgAllow`、逆向きの範囲の誤り）はこの順に入れる。
- ゲートの script の出発点：[`docs/evidence/losat_web_e2e/gates/s08pb_gates.sh`](../evidence/losat_web_e2e/gates/s08pb_gates.sh)（`GATE` で build の directory を分ける）と `s08pb_gate_a.sh`。fixture は 151 件（E2e の終わり）。V-PERF の変更前は E2e の最後の成果物（`~/.cache/losat-web-gui-target/s08pb2-gate-*`、native `6f070575…`）。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S12 — 検索画面](session_s12_w3_search_ui.md)。範囲指定の意味（端の扱い、単位）を S12 の指示書に書き足す。
