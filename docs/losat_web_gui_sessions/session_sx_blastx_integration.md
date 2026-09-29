# Session SX — BLASTX の統合（条件付き）

## INSTRUCTION PROMPT

LOSAT の BLASTX（LOSATX）を LOSAT Web の構成に組み入れる。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の SX の行である。背景は計画の DW-10（BLASTX は LOSATX 計画の認証の後に扱う）と DW-11（BLASTX も領域の指定の対象にする）である。

入口の条件：LOSATX 計画（`docs/losatx_blastx_v0.2.0_plan.md`）の v0.2.0 の認証（Stage G）が合格し、その記録が `main` に入っていること。入っていなければ、何も変えずに止め、README の表の状態を「条件待ち」のままにする。保守者が、BLASTX の範囲の拡大（`-query_loc` / `-subject_loc`）を LOSATX 計画の範囲の記録に書いていることも確かめる。

1. `git merge origin/main` で取り込み、衝突の解消を独立したコミットにする。取り込んだ後、S02〜S08 のゲートのうち、取り込みで変わったコードに関わるものを再実行する。
2. 基準：取り込んだ後のこの worktree で、LOSATX 計画の比較ゲートの出力の SHA-256 と、BLASTX の 1 つの fixture の性能（1 回の暖機と 3 回の計測、ネイティブと command-WASI）を取る。
3. BLASTX を `run_local` に通す。BLASTX は 10,002 文字の query バッチごとに検索と整形を繰り返す（`LOSAT/src/algorithm/blastx/web.rs`、`native.rs`、`query_setup.rs`）。この順序（prolog → バッチごとの検索と整形 → epilog）を保つ。formatter の中で書いている query ごとの警告（`LOSAT/src/algorithm/blastx/report.rs`、`LOSAT/src/report/pairwise.rs`）を `diagnostics` に移し、形式の数だけ重ならないようにする。outfmt 0 の既定のアラインメント数（`LOSAT/src/algorithm/blastx/args.rs` の `num_alignments`）は形式ごとに解決する。CLI の `run` と v1 の BLASTX の経路を、`run_local` の薄い層にする。
4. `hits` と観測者をつなぎ、`docs/web/verification_cells.tsv` の BLASTX の行を埋める。ABI v2 の `register` と `scan`（NCBI 型の解析器）を BLASTX で使えるようにし、V-ABI と V-BR に BLASTX の升目を加える。
5. `-query_loc` / `-subject_loc` を BLASTX に移植する（S11 と同じ手順。NCBI の範囲指定の意味は `c++/src/algo/blast/blastinput/blast_args.cpp:1946-1997`、`:2373`、`c++/src/algo/blast/blastinput/blast_fasta_input.cpp:433-459`）。`LOSAT/src/cli.rs` の BLASTX の範囲外の一覧から 2 つを外す。比較用の fixture を固定し、NCBI とバイト一致させる。
6. 既定以外のオプション（TD-13）：S08+ と同じ手順で、BLASTX の検索のオプションの既定以外の値を NCBI と比べ、NCBI と同じにするか、明示的に拒否する（`docs/losat_web_gui_sessions/session_s08p_e2e_protein_options.md`）。直した値と拒否する値を S12 の指示書（または、S12 が終わっていれば検索画面）に反映する。
7. 試験：LOSATX 計画の比較ゲートが第 2 項の基準と一致すること（移植した範囲指定を除く）。範囲指定の fixture で NCBI とバイト一致。BLASTX の全升目（0/6/7 × スレッド 1/2/4）の V-ABI と V-BR。v1 の reactor の検査。性能の中央値が基準の +5% 以内。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
8. 独立監査を受ける。

記録は `docs/evidence/losat_web_sx/README.md`。1 セッションで終わらなければ、まとめ直し（SXa）と範囲指定（SXb）に分ける。

## 終了・引き継ぎ

README の規則 8 に従う。SX を終えた後は、中断していた順番のセッションに戻る。
