# Session S02 — E1a：核の入口（共通部と BLASTP）

## INSTRUCTION PROMPT

LOSAT の段階 E1a を実行する。これはエンジン（`LOSAT/`）の変更で、**出力を 1 バイトも変えない**まとめ直しである。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S02 の行である。設計は計画 §4.1〜§4.5 と TD-2〜TD-5 にある。`PD-LOSAT-WEB-APP-BOUNDARY`（Accepted）が、ここで作る配管を許可している。

制約：BLASTX の入口と、BLASTX だけが使う formatter の経路は変えない（計画 DW-10、§4.2）。BLASTX と共有する関数（`LOSAT/src/report/pairwise.rs`、`LOSAT/src/report/outfmt6.rs` の一部）に観測者が要るときは、関数の本体と引数を変えず、呼出し側か新しい関数で受ける。

1. **全 program の基準を 1 回だけ取る**（計画 §6.3）。凍結バイトがある升目はそれを使う（Gate A の `LOSAT/tests/platform_native_v010_canonical.tsv`、`docs/evidence/tlosan_stage_g/`）。無い升目は、変更の前にこの worktree で出力を取り、SHA-256 を `docs/evidence/losat_web_e1a/baseline/` に残す。対象は 5 program の既存の回帰ゲート（BLASTP：`LOSAT/tests/blastp_v010_parity_manifest.tsv` と `LOSAT/tests/audit_blastp_v010.py`、`LOSAT/tests/run_blastp_*.sh`。BLASTN：`LOSAT/tests/blastn_parity_manifest.tsv` と `LOSAT/tests/compare_blastn_parity.py`。TBLASTX：`LOSAT/tests/tblastx_v010_parity_manifest.tsv` と `LOSAT/tests/audit_tblastx_v010.py`。TBLASTN：TLOSAN 計画の Stage G。BLASTX：LOSATX 計画の比較）と、v1 の serial / threaded reactor の検査（`LOSAT/tests/check_wasm_threading.py` が呼ぶ `check_wasi_reactor.js`、`check_wasi_api_limits.js`）である。NCBI の oracle が要るゲートで、その実行ファイルがこの環境に無いものは、凍結した期待値との比較に限り、そのことを記録する。S03・S04 もこの基準を使う。
2. **性能の基準**（計画 §6.2 の V-PERF）：BLASTP の 1 つの fixture で、ネイティブと serial / threaded の command-WASI の時間を、1 回の暖機と 3 回の計測で取る（中央値と範囲）。
3. NCBI の構成を読む：`c++/src/app/blast/blastp_app.cpp`（と `c++/src/app/blast/blastx_app.cpp:197-297`）で、subject 集合を一度用意し、query バッチごとに `CLocalBlast` で検索して `CBlastFormat` で整形する流れを記録する。`c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978` から、出力形式で検索オプション（hitlist の大きさ）が変わる条件を記録し、LOSAT の BLASTP ではそれが起きない（`-num_descriptions` / `-num_alignments` を受け付けない）ことを確かめる。
4. `LOSAT/src/api/local_blast.rs` に核の入口 `run_local` を作る（計画 §4.2 の形）。どのターゲットでもコンパイルされるようにする。`ReportOutputs` は、形式ごとの writer、`diagnostics`（CLI が stderr に書く警告）、`hits`（並べ終えた最終一覧の `PairwiseHit` と HSP の ID）、観測者（HSP の行・節の開始と終了）を持つ。形式ごとに解決した検索オプションが食い違えば、明示的なエラーで止める（TD-4）。
5. BLASTP を `run_local` に通す。BLASTP にはすでにレコードを受け取る経路（`LOSAT/src/algorithm/blastp/blast_engine.rs` の `run_resolved_with_records`）があり、出力は最後の `match outfmt`（同ファイルの 6139 行付近）で 1 回だけ分かれている。CLI の `run` と v1 の `run_web_pair` / `run_web_pair_records` は、入力を変換して `run_local` を呼ぶだけの層にする。要求された形式ごとに、同じ最終結果から writer を呼ぶ。formatter には観測者の呼出しを足すが、書くバイトは変えない。どの変更にも、それが中継する NCBI の箇所の参照コメントを付ける。
6. `docs/web/verification_cells.tsv` を作る：program × 形式 × プロファイル × 実行経路の升目ごとに、期待値の出どころ（認証済みの凍結バイトか、同じ commit のネイティブ CLI の単一形式の出力か）とそのファイルを書く（TD-5）。BLASTP の行を埋め、ほかの program は行だけ用意する。
7. 試験：
   - 第 1 項の基準と、変更の後の出力の SHA-256 がすべて一致すること（BLASTP だけでなく、変更した共有のコードを使う全 program）。
   - v1 の serial / threaded reactor を既存の方法（`LOSAT/tests/build_wasi_artifacts.py`）で作り直し、v1 の reactor の検査が通ること。
   - `run_local` に 3 形式を同時に要求した出力が、それぞれ単一形式の CLI の実行と一致すること。
   - 観測者が報告した HSP の数と範囲が、outfmt 6 の行数・outfmt 0 の節と合うこと（`-max_target_seqs` の境界の fixture を含める）。
   - 警告が 1 回だけ出ること。
   - 第 2 項と同じ条件の時間の中央値が、基準の中央値の +5% 以内であること。
   - `cargo fmt --check`、`cargo clippy --all-targets --all-features -- -D warnings`、`cargo test --all-features`、`.github/workflows/web.yml` の `engine-web-api` が通ること。
8. 独立監査（README の「レビュー」）を受ける。

記録は `docs/evidence/losat_web_e1a/README.md`（ゲート記録、`evidence.sha256`、`run-<UTC>/`）。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S03 — 核の入口：TBLASTN](session_s03_e1b_core_entry_tblastn.md)。`run_local` と `ReportOutputs` の確定した形、観測者の呼び方、基準のファイルの場所、BLASTP で分かった注意点を、S03・S04・SX の指示書に書き足す。
