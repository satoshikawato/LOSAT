> Successor: [2026-09-13 threading remediation](../wasm_threading_remediation_plan_20260913.md) replaces the earlier pool-size and reactor lifecycle assumptions. Evidence below is historical.

# 進捗台帳 — LOSATN / LOSATP Wasm性能改善

本ファイルが、このパッケージに沿う開発の進捗の正本である。現時点では計画作成のみで、下記sessionは実行していない。S00が実在pathと初期baseを記入する。各担当は自分の状態・成果path・次の推奨だけを更新し、無関係な状態を成功へ変更しない。

## 実行環境

| 項目 | 値 |
|---|---|
| PLAN_DIR | 未設定。S00が本パッケージのrootを記録 |
| REPO_ROOT | 未設定。S00が実Git rootを記録 |
| WORK_DIR | 未設定。S00が永続evidence/report領域を記録 |
| コード読解SHA | `7db9bb0060e4e057f9f50807bf9edc2362f20133` |
| 初期固定base | 未設定 |
| 直近採用base / 統合candidate | 未設定 |
| gbdraw checkout / SHA | 未設定 |
| NCBIソース版・golden所在 | 未設定 |
| 主runtime・browser | 未設定 |
| 測定予算・採用基準の記録 | S00 report作成後に参照 |
| 独立レビュー担当 | 未設定 |

## セッション

| ID | 内容 | 状態 | 実report / candidate / 判断 |
|---|---|---|---|
| [S00](INSTRUCTION_PROMPTS/S00_baseline_contracts.md) | base・互換性契約・fixture・実行予算の固定 | NOT_STARTED | 未実施 |
| [S01](INSTRUCTION_PROMPTS/S01_profiling_selection.md) | 最小限のWasm計測と改善候補の選択 | NOT_STARTED | 未実施 |
| [S02](INSTRUCTION_PROMPTS/S02_n_traceback_copy.md) | LOSATN greedy tracebackの往復コピー削減 | NOT_STARTED | 未実施 |
| [S03](INSTRUCTION_PROMPTS/S03_n_query_pack.md) | LOSATN query再packの削減 | NOT_STARTED | 未実施 |
| [S04](INSTRUCTION_PROMPTS/S04_n_match_kernel.md) | LOSATN一致延伸カーネルの局所最適化 | NOT_STARTED | 未実施 |
| [S05](INSTRUCTION_PROMPTS/S05_n_existing_parallelism.md) | LOSATN既存並列仕事の効率化 | NOT_STARTED | 未実施 |
| [S06](INSTRUCTION_PROMPTS/S06_p_single_thread_kernel.md) | LOSATP単一スレッドの上位カーネル改善 | NOT_STARTED | 未実施 |
| [S07](INSTRUCTION_PROMPTS/S07_p_diagonal_reset.md) | LOSATP parallel preliminaryの初期化削減 | NOT_STARTED | 未実施 |
| [S08](INSTRUCTION_PROMPTS/S08_p_kappa_scratch.md) | LOSATP Kappa match-redoのscratch再利用 | NOT_STARTED | 未実施 |
| [S09](INSTRUCTION_PROMPTS/S09_p_bounded_redo.md) | LOSATP Kappa先行redoの有界化 | NOT_STARTED | 未実施 |
| [S10](INSTRUCTION_PROMPTS/S10_browser_runtime.md) | ブラウザ実測と限定的なruntime統合 | NOT_STARTED | 未実施 |
| [S11](INSTRUCTION_PROMPTS/S11_qualification_review.md) | 統合候補の独立検証と最終結論 | NOT_STARTED | 未実施 |
| [S12](INSTRUCTION_PROMPTS/S12_n_partition_design.md) | LOSATN単一仕事を分割する設計と隔離検証 | NOT_STARTED | 未実施 |

## 次の一手

最初はS00を実行する。S00/S01後の次候補は実profileを根拠に選ぶ。S02〜S09を番号順に全部実装することは要求しない。S12は所有者が選択した場合だけ開始する。

## 状態表記

`NOT_STARTED / READY / IN_PROGRESS / COMPLETE / ACCEPTED / REJECTED / SKIPPED / BLOCKED / INCONCLUSIVE`を使用する。意味はMASTER_PLAN §12に従う。独立レビューの`REVIEW_PENDING`とplatformの`NOT_RUN`はreport内で区別する。

## 運用メモ

WIPは原則一つ。計画を別場所へ移した場合はpathを更新する。初期配布版を再解凍して実行済みSTATUSを上書きしない。レポート本文やraw timingをこの台帳へ重複転記しない。
