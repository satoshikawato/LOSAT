# Session S11 — E2d：`-query_loc` / `-subject_loc`

## INSTRUCTION PROMPT

LOSAT の段階 E2d を実行する。これはエンジン（`LOSAT/`）の NCBI パリティ作業である。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S11 の行である。

目的：設計書 §2.1 の「領域は座標入力とプレビュー選択」を、NCBI の範囲指定の意味どおりに、BLASTN・BLASTP・TBLASTN・TBLASTX で使えるようにする（BLASTX は SX。計画 DW-10、DW-11）。範囲は、その役割のレコードが 1 つのときだけ画面で指定できる（計画 DW-9）が、CLI と同じく、エンジンは複数のレコードでも NCBI の意味どおりに扱う。UI で配列を切り出して検索することはしない。

1. 現状を確かめる：BLASTN・BLASTP・TBLASTX は `-query_loc` / `-subject_loc` を定義しておらず、CLI が未知の引数として拒否する。TBLASTN は `-query_loc` を未移植の名前として拒否し（`LOSAT/src/cli.rs:187-193`）、`-subject_loc` も拒否する（`LOSAT/src/algorithm/tblastn/args.rs:364`）。
2. NCBI の経路を追う：`c++/src/algo/blast/blastinput/blast_args.cpp` の `kArgQueryLocation`（定義は 1946 行付近、解析は 1996〜1997 行）と `kArgSubjectLocation`（2373 行付近）から、範囲が配列・統計・座標の表示にどの時点で効くかを program ごとに記録する。範囲は入力元が読むすべてのレコードに適用される（`c++/src/algo/blast/blastinput/blast_fasta_input.cpp:433-459`。逆向きの範囲はエラー、長さを超える終点は黙って切り詰め）。
3. 比較用の fixture を決める（4 program、範囲の端、範囲の内外にまたがるヒット、minus 鎖、翻訳の frame、範囲を全レコードに適用する複数レコードの入力、範囲の始点より短いレコード（NCBI ではエラー、`blast_fasta_input.cpp:446-452`））。NCBI BLAST+ の出力を oracle として固定し、SHA-256 を記録する。
4. 4 program に移植する。移植した箇所の直上に、NCBI のファイル・行と断片を書く。NCBI に無い組合せは明示的に拒否する。`run_local` と ABI v2 の `validate` から使えるようにする。
5. すべての fixture で outfmt 0/6/7 が NCBI とバイト一致すること、既存の認証済みの出力に退行が無いことを確かめ、`docs/web/verification_cells.tsv` に升目を足す。独立監査を受ける。

完了条件は計画 §7 の S11 の行による。記録は `docs/evidence/losat_web_e2d/README.md`。1 セッションで終わらなければ S11b に分ける。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S12 — 検索画面](session_s12_w3_search_ui.md)。範囲指定の意味（端の扱い、単位）を S12 の指示書に書き足す。
