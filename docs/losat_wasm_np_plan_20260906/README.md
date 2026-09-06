# LOSATN / LOSATP Wasm性能改善 — 実行パッケージ

版1.0 / 2026-09-06

**LOSATNはコピー・反復演算を削減し、LOSATPは重い計算と並列化の余分な仕事を削減する。** 結果の一致を保ちながらNativeとの性能差を詰めるための、総合計画とsession別の実行指示である。

このZIPには計画書・指示書・未記入テンプレートだけを含む。LOSAT/gbdrawのソース、NCBIソース、golden output、実行バイナリ、新しいベンチマーク結果は含まない。コード変更や性能検証はまだ実行していない。

## 最初の使い方

1. [MASTER_PLAN.md](MASTER_PLAN.md)を読み、対象範囲、既存互換性契約、採用基準、session依存を確認する。
2. 開発環境でこのフォルダを参照できる状態にし、[S00の指示書](INSTRUCTION_PROMPTS/S00_baseline_contracts.md)を新しい作業sessionへ渡す。repositoryの場所などはS00が実在環境から確認して[STATUS.md](STATUS.md)へ記録する。
3. S00とS01の後は、STATUSが指す必要なsessionの指示書を一つずつ渡す。新担当へは該当指示書に加え、このパッケージとSTATUSに登録したreport/evidenceへのアクセスを与える。過去の会話を別途説明する必要はない。

単独の指示書にも対象、目的、開始条件、制約、読取資料、手順、完了条件を記載している。ただし共通の正式基準はMASTER_PLAN、実行状態はSTATUS、根拠は実sourceとreportで管理する。指示書だけを切り出して必要資料が利用できない場合、担当者は不足物を記録し、未検証を合格にしない。

## ファイル構成

| ファイル | 用途 |
|---|---|
| [MASTER_PLAN.md](MASTER_PLAN.md) | 総合計画、アーキテクチャ、意味の不変条件、実験設計、採否・停止基準 |
| [SOURCES.md](SOURCES.md) | commit固定のコード参照と一次資料、確認範囲 |
| [STATUS.md](STATUS.md) | 実行path・base・session状態・次の作業を管理する唯一の台帳 |
| [templates/SESSION_REPORT.md](templates/SESSION_REPORT.md) | 各sessionの証拠・判断・引き継ぎ用の共通テンプレート |
| `INSTRUCTION_PROMPTS/` | 下表の13本の実行指示書 |

## セッション一覧

| ID | 指示書 | 実施条件 |
|---|---|---|
| S00 | [base・互換性契約・fixture・実行予算の固定](INSTRUCTION_PROMPTS/S00_baseline_contracts.md) | 必須／原則read-only |
| S01 | [最小限のWasm計測と改善候補の選択](INSTRUCTION_PROMPTS/S01_profiling_selection.md) | 必須／計測の欠落だけ実装可 |
| S02 | [LOSATN greedy tracebackの往復コピー削減](INSTRUCTION_PROMPTS/S02_n_traceback_copy.md) | 条件付き／局所実装 |
| S03 | [LOSATN query再packの削減](INSTRUCTION_PROMPTS/S03_n_query_pack.md) | 条件付き／局所実装 |
| S04 | [LOSATN一致延伸カーネルの局所最適化](INSTRUCTION_PROMPTS/S04_n_match_kernel.md) | 条件付き／局所実装 |
| S05 | [LOSATN既存並列仕事の効率化](INSTRUCTION_PROMPTS/S05_n_existing_parallelism.md) | 条件付き／スケジューリング限定 |
| S06 | [LOSATP単一スレッドの上位カーネル改善](INSTRUCTION_PROMPTS/S06_p_single_thread_kernel.md) | 条件付き／一カーネル限定 |
| S07 | [LOSATP parallel preliminaryの初期化削減](INSTRUCTION_PROMPTS/S07_p_diagonal_reset.md) | 条件付き／状態管理の局所変更 |
| S08 | [LOSATP Kappa match-redoのscratch再利用](INSTRUCTION_PROMPTS/S08_p_kappa_scratch.md) | 条件付き／run内の寿命に限定 |
| S09 | [LOSATP Kappa先行redoの有界化](INSTRUCTION_PROMPTS/S09_p_bounded_redo.md) | 条件付き／順序保持のスケジューリング変更 |
| S10 | [ブラウザ実測と限定的なruntime統合](INSTRUCTION_PROMPTS/S10_browser_runtime.md) | 必須の評価／host変更は条件付き |
| S11 | [統合候補の独立検証と最終結論](INSTRUCTION_PROMPTS/S11_qualification_review.md) | 必須／read-onlyレビューを主とする |
| S12 | [LOSATN単一仕事を分割する設計と隔離検証](INSTRUCTION_PROMPTS/S12_n_partition_design.md) | 拡張／設計・非production試験のみ |

S00・S01・S10・S11は基準固定、原因把握、ブラウザ評価、最終評価を担う。S02〜S09の8本は条件付きの実験で、効果や適用性が確認できない候補はSKIPPED/REJECTEDにする。S12は単一仕事のN検索に新しい分割が必要な場合の拡張設計であり、productionの既定経路を変更する許可ではない。

## 何を変えず、何を変えるか

検索感度・候補集合・scoring・統計・順序・既存errorを保ち、データの持ち方、専有scratchの寿命、必要な初期化、既存独立仕事の割当を改善する。NCBIは意味仕様とテストoracleであってruntime依存にはしない。既存のplatform別契約と固定raw outputを守り、比較器の緩和で結果差を隠さない。詳細と根拠は総合計画§5〜§8に記載している。

SOLID/KISS/DRY/YAGNIは、traitやframeworkの増設ではなく、責任の分離、同じ契約での置換、局所的な変更、一実験一判断、単一の進捗台帳、不要候補の棄却へ適用する。全面rewriteや新しい汎用schedulerを前提にしない。

## 実行成果の保存場所

S00が選ぶ永続`WORK_DIR`配下に`reports/`と`evidence/`を置く。各reportは`reports/Sxx.md`。この配布パッケージのSTATUSへ実在pathを記録する。`/tmp`の消失や別sessionのファイル不可視を前提に引き継がない。作業領域が別環境へ移る場合は、実在artifactとhashを確認してpathを更新する。

remoteへのpush、merge、release公開、配信設定の変更は、このパッケージだけでは許可しない。最終評価が完了しても、公開・採用の操作は所有者の判断と権限に従う。
