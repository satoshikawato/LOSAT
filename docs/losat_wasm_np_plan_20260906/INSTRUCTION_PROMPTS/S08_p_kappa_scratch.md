# S08 — LOSATP Kappa match-redoのscratch再利用

種別: 条件付き／run内の寿命に限定  
前提: S01で単一query match-redo経路の生成費を確認

## このセッションの目的

LOSATPの単一query match-redoにあるworkspace/scratch生成を、意味状態を混ぜずに再利用する。pool寿命をブラウザ全体へ拡大せず、まず同一検索run内に閉じる。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

`blastp/blast_engine.rs`の`use_parallel_kappa_match_redo`、`postprocess_preliminary_hits()`とworkspace構築を読む。読解SHAではmatch closure内でcomposition workspace、gap scratch、subject range cache等を作る。multi-query経路にも同じ問題があると一般化しない。[R10](../SOURCES.md#r10)

`map_init`/`for_each_init`だけではworkerにつき一度の生成を保証しない。初期化回数を実測する。[R14](../SOURCES.md#r14)

## 実施手順

1. 対象branchへの到達とmatch数、workspace生成回数、allocation bytes、初期化時間を確認する。準備費と実redo計算を分ける。
2. 各fieldをcapacity等の物理資源、query依存、subject/range依存、adjusted matrix依存、結果所有に分類する。reset時に必要なNCBI stateを表にする。
3. まず明示的な小batchまたは既存parallel job局所のscratch所有を選ぶ。結果は元のmatch indexに対応付ける。query共通不変情報を安全に共有できる場合も、既存の更新有無を確認する。
4. stale subject rangeやadjusted scoreが次matchへ漏れないresetを実装する。共有Mutex workspaceやcross-run/global cacheを導入しない。
5. S09のredo仕事数削減と同時に変更せず、仕事数を同じにした比較で再利用だけの効果を測る。

## 必須のfocused検証

一match、多数match、大小subject混在、同じqueryで異なるsubject、異なるqueryを連続run、空/失敗経路、adjusted matrixが変わるケースを確認する。single-query/multi-queryの既存raw outputとエラー結果を保持する。

通常releaseでWasm serial、threaded n1/n2/n4/n8、Nativeを比較し、scratch生成数の減少とtotal search時間の両方を示す。初回allocationとhigh-water markを分け、反復に比例したメモリ増加がないことを確認する。

## 採用判定

生成費が下がり、各状態のreset根拠が明確で、S00の非回帰基準を満たす場合ACCEPTED。小caseへの固定費増加が支配的、stateの寿命を説明できない、効果が仕事数変更と混ざった場合はREJECTED/INCONCLUSIVE。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
