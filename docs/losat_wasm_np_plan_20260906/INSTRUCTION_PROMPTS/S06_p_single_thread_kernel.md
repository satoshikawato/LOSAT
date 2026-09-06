# S06 — LOSATP単一スレッドの上位カーネル改善

種別: 条件付き／一カーネル限定  
前提: S01の信頼できるWasm段階別profile

## このセッションの目的

LOSATPのgapped alignment、traceback、Kappa関連処理から、単一スレッドで実際に時間を使う一処理を選んで改善する。既存の特殊化やscratch設計を再実装しない。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

`blastp/blast_engine.rs`、`gapalign.rs`、profileで指された直近callerを読む。読解SHAではBLOSUM62専用経路とadjusted matrixのrow取得が存在する。[R10](../SOURCES.md#r10) [R11](../SOURCES.md#r11)

S01がWasm計測を得ていないなら、Native profileだけでWasmの主因を確定しない。S08のscratch再利用やS09の先行redoが主因なら、重複作業を避けて本sessionをSKIPPEDにする。

## 実施手順

1. 主対象caseとguardを指定し、時間割合/absolute時間、calls、処理されたcellやcopy bytes等の説明指標を一つ以上固定する。
2. 該当NCBI DP/traceback/Kappa処理についてmatrix、gap penalty、X-drop、restricted/exact retry、tie-break、edit script、統計の適用順を確認する。
3. DP行の局所メモリ配置、不要なclone/copy、生存領域の初期化、score lookupの残存dispatch等から最小の変更を一つ選ぶ。profileで未証明の大規模SIMD DPや別alignmentアルゴリズムへの置換は行わない。
4. Kappaのadjusted matrixと通常BLOSUM62を同じ契約で扱い、precisionと丸めを保持する。CBS/SEGを省略しない。速くするために候補HSP数を変えない。
5. 小パッチを通常releaseでA/B測定する。計測専用変更が速度に影響していないか確認する。

## 必須のfocused検証

同点のgap選択、長短alignment、restrictedからexact retryへの切替、標準/adjusted matrixが実際に通る既存ケース、異常入力の既存error、query/subjectを変えた反復を確認する。対象taskの既存raw output契約を通す。

Native n1、Wasm serial n1、threaded artifact n1の差を別々に示す。共通カーネル変更はP-MTもguardとして測り、各workerのscratch容量やpeak memoryが急増していないか確認する。

## 採用判定

S00の速度/非回帰条件を満たし、変更が既存意味と同じことを説明できればACCEPTED。kernelだけ速く検索全体が変わらないなら、必要性を再評価してREJECTED/INCONCLUSIVEにする。速度が出るまで別の仮説を同じsessionへ継ぎ足さない。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
