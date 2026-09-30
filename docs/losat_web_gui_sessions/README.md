# LOSAT Web GUI セッション用 INSTRUCTION PROMPTS

[総合計画書](../losat_web_gui_plan.md)から、セッションごとの実行指示を分けたものである。各ファイルの「INSTRUCTION PROMPT」は、1 つの作業セッションにそのまま渡す指示文である。セッションは表の順に 1 つずつ実行する。各段階の完了条件は計画書 §7 の表が正本で、指示書はそれに細部を足すだけである（条件を緩めない）。

## 全セッションに共通する規則

各指示書は、この節を前提にしている。

1. **ブランチと作業場所**：エンジン側（`LOSAT/`・`web/adapter/` を変えるセッション）のブランチは `feature/losat-web-gui`（上流は `origin/feature/losat-web-gui`）、作業ディレクトリは git worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui`。アプリ側（`web/app/` だけを変えるセッション）のブランチは `feature/losat-web-gui-app`（上流は `origin/feature/losat-web-gui-app`）、作業ディレクトリは git worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui-app`。2 本までを並行し（計画 DW-7）、アプリ側はセッションの終わりに `feature/losat-web-gui` へ merge する。アプリ側は `LOSAT/`・`web/adapter/`・計画・この README の表を変えず、変えるべきことはゲート記録に書く（エンジン側が merge のときに反映する）。これ以外の新しい clone や worktree（一時的なものを含む）は作らない。この worktree が無い環境に限り、既存の LOSAT のクローンで一度だけ `git fetch origin && git worktree add /mnt/c/Users/genom/GitHub/LOSAT-web-gui feature/losat-web-gui` を実行し（クローンも無ければ、一度だけ `git clone https://github.com/satoshikawato/LOSAT.git` を実行し）、以後それを使い続ける。
2. **開始時の確認**：worktree で `git branch --show-current` が `feature/losat-web-gui` であり、`git status` に未コミットの変更が無いことを確かめ、`git pull --ff-only` を行う（取り込むのは `origin/feature/losat-web-gui` だけ）。`main` の変更が必要なときだけ `git merge origin/main` を行い、衝突の解消を独立したコミットにする。
3. **最初に読むもの**：ルートの `AGENTS.md`、`web/AGENTS.md`、総合計画書、`docs/product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md`、前のセッションのゲート記録（`docs/evidence/losat_web_<stage>/README.md`）。要求の詳細は `docs/web/losat_web_design_v0.1.md`（設計書）と `docs/web/requirements_trace.tsv` にある。
4. **エンジン（`LOSAT/`）を変更するセッション**：ルートの `AGENTS.md` の必須規則と `.agents/skills/verify-ncbi-parity-and-speed/SKILL.md` に従う。NCBI のソースは `/mnt/c/Users/genom/GitHub/ncbi-blast/`（固定 commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`）で読み、NCBI の実行ファイルは比較のためだけに使う。Rust は `cargo +1.92.0` を使い、ビルドの出力先は worktree の外（`--target-dir /home/kawato/.cache/losat-web-gui-target/<用途>`）にする。`rustc` や `cargo` を直接呼ぶ既存のスクリプト（例：`LOSAT/tests/build_wasi_artifacts.py`）は、環境変数 `RUSTUP_TOOLCHAIN=1.92.0` を付けて実行する。出力を変えない変更の前後の比較には `docs/evidence/losat_web_e1a/capture_outputs.py`（全 program の出力の SHA-256）と `docs/evidence/losat_web_e1a/measure_perf.py`（性能の非退行。変更前と変更後を 1 回ずつ交互に測る）を使う。`cargo test --all-features` は診断用の feature を含むので、CI と同じく環境変数 `LOSAT_BLASTX_WORKER_LOG` に書き込めるファイルのパスを設定する。回帰の基準には、凍結バイト（例：`LOSAT/tests/platform_native_v010_canonical.tsv`）か、変更の前にこの worktree で取った出力を使う。
5. **アプリ（`web/`）を変更するセッション**：`web/AGENTS.md` に従う。コミットの前に `cd web/app && npm ci && npm run check && npm run e2e` を通す。エンジン側が V-PERF を測る間（ファイル `/home/kawato/.cache/losat-web-gui-target/vperf.lock` がある間）は、アプリ側の試験やビルドを始めず、実行中のものは終わるのを待つ（計画 DW-7）。
6. **証拠**：`docs/evidence/losat_web_<stage>/` に、ゲート記録の `README.md`、`evidence.sha256`、再現スクリプトを置く。実行の記録は `run-<UTC 時刻>/` に置き、作った後は書き換えない。
7. **保守者の判断**：PD の承認、公開、Cloudflare へのデプロイ、GA4 の設定、外部サービスへの送信は、セッションの中で行わない。必要になったら、ゲート記録に「保守者の判断待ち」として書き、最終回答でも示す。
8. **終了と引き継ぎ**：エンジンの変更とアプリの変更は別のコミットにする。コミットして `git push` を行い、ゲート記録に結果・コミット SHA・残件を書く。次のセッションの指示書を、実測や結果に合わせて更新し、この README の表の状態を直して、同じブランチにコミットしてプッシュする。最終回答には、コミット SHA、プッシュの結果、検証の結果、残件、次のセッションの指示書のパスと全文を示す。完了条件が残っているときは完了と書かず、その解消を次の指示書の最初の作業にする。

## レビュー

レビュー役の定義は保守者の Codex 設定（`~/.codex/agents/*.toml`）にある。使えない環境では、別のエージェントか人に、次の観点を読み取り専用で確かめてもらう。

- **計画レビュー**（`plan_critic`）：計画を現在のリポジトリと規約に照らし、古い前提、責任の境界の誤り、足りない依存、やり残しの削除、弱い完了条件、試験の汚染、破壊的な手順、満たされない要求を探す。見つけたものごとに、計画の節とリポジトリの根拠、影響、最小の修正を示し、重大な順に並べ、最後に「実行できる／修正すれば実行できる／判断が必要」のどれかで結ぶ。
- **独立監査**（`ncbi_parity_auditor`）：パリティや性能の主張を、コードを変えずに確かめる。fixture、オプション、LOSAT の commit、ターゲット、NCBI BLAST+ の版、適用した例外を確かめ、NCBI と LOSAT の呼出し経路（呼ぶ時点、エンコード、座標系、精度、並べ替え、刈り込み、整形、ネイティブと Wasm の一致、並列の集約の順序）を追う。ヒット数だけに基づく一致の主張や、最速の 1 回だけに基づく性能の主張は認めない。結論は「supported / unsupported / inconclusive」。
- **画面レビュー**（`visual_regression_reviewer`）：画面の記録を、同じ表示サイズで前の記録と比べる。はみ出し、重なり、配置、凡例、表記、狭い画面の配置、操作の対象を確かめ、見つけたものを重大な順に並べ、合否と、未解決の項目に必要な撮り直しやコードの証拠を示す。

## セッション一覧

| 順序 | セッション | 段階 | 主な完了条件 | 状態 |
|---|---|---|---|---|
| S01 | [契約と骨格](session_s01_w0_contract_skeleton.md) | W0 | FakeEngine の E2E、`crossOriginIsolated`、TBLASTX v1 の fail-fast | 完了（2026-09-29、[ゲート記録](../evidence/losat_web_w0/README.md)） |
| S02 | [核の入口：共通部と BLASTP](session_s02_e1a_core_entry_blastp.md) | E1a | 全 program の基準、変更したコードを使う全 program のゲート、v1 の検査、V-PERF | 完了（2026-09-29、[ゲート記録](../evidence/losat_web_e1a/README.md)） |
| S03 | [核の入口：TBLASTN](session_s03_e1b_core_entry_tblastn.md) | E1b | TLOSAN の Stage G のゲートが変わらない | 完了（2026-09-29、[ゲート記録](../evidence/losat_web_e1b/README.md)。V-PERF は保守者の確認待ち） |
| S04 | [核の入口：BLASTN と TBLASTX](session_s04_e1c_core_entry_blastn_tblastx.md) | E1c | 既存ゲートと Gate A が変わらない、v1 の検査 | 完了（2026-09-29、[ゲート記録](../evidence/losat_web_e1c/README.md)。V-PERF と CLI の 2 つの差は保守者の確認待ち） |
| S05 | [アダプタと ABI v2](session_s05_e1d_adapter_abi_v2.md) | E1d | 4 program の V-ABI、ビルドの同一性の検査 | 完了（2026-09-29、[ゲート記録](../evidence/losat_web_e1d/README.md)） |
| S06 | [BLASTN outfmt 0：権威と fixture](session_s06_e2a1_blastn_outfmt0_authority.md) | E2a-1 | NCBI の経路の表、固定した fixture | 完了（2026-09-29、[ゲート記録](../evidence/losat_web_e2a/README.md)） |
| S07 | [BLASTN outfmt 0：実装とゲート](session_s07_e2a2_blastn_outfmt0_port.md) | E2a-2 | NCBI とバイト一致、6/7 に退行なし | 完了（2026-09-29、[ゲート記録](../evidence/losat_web_e2a/README.md)） |
| S07+ | [BLASTN の得点のオプションと入力の読み方](session_s07p_e2c_blastn_scoring_options.md) | E2c | 既定以外の得点と入力の読み方が NCBI と同じか、明示的に拒否される（S06・S07 で追加。S12 の前に終える） | 未着手 |
| S07++ | [BLASTN の query の batch](session_s07pp_e2f_blastn_query_batches.md) | E2f | 複数の query が NCBI と同じ batch で検索される（S07+ で追加、TD-14。S12 の前に終える） | 未着手 |
| S08 | [TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md) | E2b | NCBI とバイト一致（承認済みの例外を除く） | 未着手 |
| S08+ | [BLASTP・TBLASTN・TBLASTX の既定以外のオプション](session_s08p_e2e_protein_options.md) | E2e | 既定以外のオプションが NCBI と同じか、明示的に拒否される（S07+ で追加、TD-13。S12 の前に終える） | 未着手 |
| S09 | [ブラウザでの実行基盤](session_s09_w1_browser_runtime.md) | W1 | V-BR、前処理の割合の実測（アプリ側。S10 の後、入口の条件は S08） | 未着手 |
| S09+ | 前処理キャッシュ（条件付き） | R2 | DW-8 の条件を満たした program だけ。S09 の結果で行を足す | 未定 |
| S10 | [データ層](session_s10_w2_data_layer.md) | W2 | BlockStore の契約試験、回収・保護・容量不足（アプリ側。S09 より先に行う） | 進行中（アプリ側、2026-09-30） |
| S11 | [領域の指定](session_s11_e2d_query_subject_loc.md) | E2d | `-query_loc` / `-subject_loc` が NCBI とバイト一致 | 未着手 |
| S12 | [検索画面](session_s12_w3_search_ui.md) | W3 | 研究作業と境界条件の E2E | 未着手 |
| S13 | [結果画面](session_s13_w4_results_ui.md) | W4 | 5 program の E2E、HSP と行・節の対応 | 未着手 |
| S14 | [抽出と候補](session_s14_w5_extraction_candidates.md) | W5 | 原配列との一致 | 未着手 |
| S15 | [出力と再現性](session_s15_w6_export_session.md) | W6 | 再計算しない再読込、明示的なつなぎ直し | 未着手 |
| S16 | [配信の仕上げ](session_s16_w7_delivery.md) | W7 | V-OFF・V-PRIV、プレビューでの隔離 | 未着手 |
| SX | [BLASTX の統合](session_sx_blastx_integration.md)（条件付き） | SX | LOSATX の v0.2.0 の認証が `main` に入った後の最初の区切りで実施。S17 の前に必ず終える | 条件待ち |
| S17 | [公開判定](session_s17_g_release_decision.md) | G | 初期の要求に未達が無い | 未着手 |

1 つのセッションで段階を終えられないときは、続きを同じ段階の次のセッション（例：S07b）とし、この表に行を足す。条件付きのセッション（S09+、SX）は、条件が満たされた後の最初のセッションの区切りで実行する。
