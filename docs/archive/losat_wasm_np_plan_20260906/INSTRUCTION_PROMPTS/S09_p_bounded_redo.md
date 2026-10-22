# S09 — LOSATP Kappa先行redoの有界化

種別: 条件付き／順序保持のスケジューリング変更  
前提: S01で安全に回避可能な先行計算が実在

## このセッションの目的

LOSATPのKappa match-redoで、後から不要と判定される仕事を無制限に先行計算しないようにする。既存のearly-termination、heap採用順序、結果・エラー契約を変えずに、元の順序の有界windowを検討する。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

`blastp/blast_engine.rs`のlocal_matches並列計算、`blast_compo_early_termination()`、heap replayと対応するNCBI Kappa処理を読む。read-onlyで依存関係を確認してから実装する。[R10](../SOURCES.md#r10)

対象が単一query match-redoであることを証明する。最終的にheapへ残らない結果が全て不要だったとは扱わない。redo自体が採用判定に必要な場合がある。S08が採用済みならそれを直近baseとし、scratch再利用を再実装しない。

## 実施手順

1. matchごとに元のindex、prelim情報、redo実行、early-skipped-before-use、heap判定、retainedを区別して小さな診断を行う。実行結果の意味を変えない。
2. どの判定がredo前に可能で、どれが既存heap state/redo結果を必要とするかをNCBIと現在の制御フローから確認する。先行して計算した関数の共有副作用・乱数・診断・error伝播も確認する。
3. 元のindex順の小windowを計算し、元の順序で既存heap判定を適用してから次windowを発行する候補を一つ作る。初期windowはthread数に比例した少数候補で比較し、fixture別のmagic numberを増やさない。
4. 既存が`continue`する条件を単調性の証明なしに`break`へ変えない。未計算matchを除外する規則は既存の安全な条件だけを用いる。全件必要なら元の仕事を全て実行する。
5. worker完了順での採用、共有heapの並列更新、検索感度の低下、経験則的なredo打切りはしない。一般化した非同期schedulerは作らない。

## 必須のfocused検証

early条件がほぼ効かないケース、強く効くケース、prelim同点、heap上限が作用する既存サポート条件、複数windowに跨る結果、threadsより少ないmatch、branch境界を試す。既存サポート内の条件だけを使う。

順序とerror propagationをunitレベルでも検証する。計算されないjobに潜在するerrorを、新経路だけ無条件に隠していないか確認する。現在の契約が未確定なら、その変更はBLOCKEDにする。window境界のbarrier、CPU作業量、未使用redo数、peak result memory、wall timeを記録する。

## 採用判定

安全に回避可能なredoが減り、実latency/memoryが改善し、全件必要なguardとNativeが基準内ならACCEPTED。CPU使用率が下がってもlatencyが改善する結果は許容する。仕事量だけ減って負荷分散が悪化する、または意味を証明できない場合REJECTED。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
