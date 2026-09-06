# S07 — LOSATP parallel preliminaryの初期化削減

種別: 条件付き／状態管理の局所変更  
前提: S01でdiagonal fillまたはscratch初期化が有意

## このセッションの目的

LOSATPのparallel preliminaryで、subjectごとのdiagonal配列初期化や過剰なscratch生成を削減する。独立subjectの状態分離とdiagonal offsetの意味を保持する。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

`BlastpSubjectScratch`、`prepare_independent_subject()`、`finish_subject()`、`blastp_subject_diag_offsets()`、prelimの`for_each_init`を読む。読解SHAのparallel処理はsubjectごとにdiag配列をfillする。既存scratch自体は再利用されているため、新しいpoolの導入を既定案にしない。[R10](../SOURCES.md#r10)

`for_each_init`の初期化回数をOS worker数と同じだと仮定しない。固定Cargo.lockに対応するRayonの契約を読む。[R14](../SOURCES.md#r14)

## 実施手順

1. P-many-shortとdense/long-subject対照で、scratch生成回数、fill bytes、実際に触れたdiag要素の密度、stage時間を記録する。
2. NCBIのsubject遷移、diagonal clear、offset増分/overflow処理と、現在のserial/parallelでの独立化の方法を対応付ける。
3. 明示batchで生成回数を減らす案、touched-index reset案、generation案のいずれか最も小さいものを一つ選ぶ。全方式を同時に実装しない。
4. generation方式は全read/writeを対象に有効性確認を行い、世代overflowのfull resetを定義する。touched-index方式は初期登録と重複登録の費用を含める。stateを隠すグローバルcacheは作らない。
5. query/subjectごとの初期値とdiag offsetを維持し、既存serialが持つ効率を悪化させない。

## 必須のfocused検証

短subjectを多数順に実行、subject順の入力変更、denseアクセス、空/短配列、同じscratchへの長→短→長、forced generation-wrapの小さなtestを行う。generationを採らない場合、そのtestは不要と明記する。全thread数の出力順を確認する。

Wasm n1/n2/n4/n8とNativeのtime/memoryを比較する。fillが減っても各diagアクセスのbranch増加で悪化した場合は採用しない。cache localityと追加index bufferを含むpeak memoryも測る。

## 採用判定

対象ケースで初期化費用とsearch時間が改善し、dense/small guard・Nativeの基準を満たす場合ACCEPTED。全状態readを把握できない、serialとの意味差が解消できない場合REJECTED。対象負担が微小ならSKIPPED。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
