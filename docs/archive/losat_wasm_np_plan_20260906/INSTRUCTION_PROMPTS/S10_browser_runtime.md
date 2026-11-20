# S10 — ブラウザ実測と限定的なruntime統合

種別: 必須の評価／host変更は条件付き  
前提: S00/S01、測定可能なgbdraw checkoutと既存host

## このセッションの目的

N/Pの改善を実際のブラウザ利用で検証し、Wasm計算とworker/転送/初期化の費用を分ける。Node上の成功をブラウザ全般の成功へ読み替えず、host側変更が必要な場合も観測した一つの費用に絞る。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

gbdrawはLOSATをブラウザから呼ぶ利用側アプリである。対象checkoutのSHA、実際にロードするserial/threaded artifact、該当service/workerを固定する。静的読解の参照には`gbdraw/web/js/services/losat.js`と`workers/losat-threaded-worker.js`があるが、現在の経路を優先する。[R13](../SOURCES.md#r13)

LOSATの新規公開API、永続shared executor、cross-run配列cache、配信基盤全面変更は本sessionの既定scopeに含めない。gbdrawの一般描画・ラベル・Undo等は触らない。配信環境の変更やremote公開は別許可が必要である。

## 実施手順

1. 主browserの版、OS、capability、`crossOriginIsolated`、Wasm feature、実際のworker spawnを確認する。Node/WSLとWindows browserの違いなど、完全に揃わない条件は明示する。[R18](../SOURCES.md#r18)
2. 配布artifact hashとbuild sourceを対応付ける。cacheが古いartifactを返していないか、task/threadsが実際に渡るか、serial fallbackの理由が何かを確認する。通常ユーザー経路を特殊設定へ置換しない。
3. worker起動、compile/instantiate、input/prepare、search、output/転送、teardownを測る。初回、同一worker再使用、同一instance再使用を明確に区別する。`_start`を同じinstanceで再度呼んで安全と仮定しない。
4. N-one-jobとN-MT-many、P-smallとP-heavyを通す。1ペアlatencyと複数ペアthroughputを別系列にする。外側pair並列数と内側thread予算の積、実worker数、memoryを記録する。
5. host費が主因なら既存module cache、入力転送、不要worker準備などから小さい変更を一つ選ぶ。既にあるmodule共有を新設しない。Rust poolの現在thread参加やregistry寿命を確認せず永続化しない。[R15](../SOURCES.md#r15)
6. 既存APIの再入性、結果所有、環境設定が初回で固定される状態を確認する。再利用できない部分は明示し、新APIが必要なら別決定へ渡す。hostだけでRustの意味stateを無理に初期化しない。

## 必須のfocused検証

連続run、入力変更、thread予算変更、通常の失敗処理、既存取消後の再実行、worker終了、memory高止まり/増え続ける状態を確認する。未実装のcancelを追加する必要はない。共有host変更は既存TBLASTXの小さい回帰fixtureも通す。

secure/cross-origin-isolated条件を満たす通常運用で測る。ブラウザ保護を無効化するflagsでthreaded成功を作らない。サポートする別engineの互換性を確認し、利用できないものはNOT_RUNとする。特殊JITモードの数値は診断のみ。

## 完了判定

host変更なしでも、実browserのresult/latency/memoryと制限が記録できればCOMPLETE。限定変更を採る場合はS00基準に従いACCEPTED/REJECTEDを併記する。browser環境がなければBLOCKEDとし、S11にNode-onlyの範囲を正確に渡す。未測定browserを含む「Native並み」を宣言しない。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
