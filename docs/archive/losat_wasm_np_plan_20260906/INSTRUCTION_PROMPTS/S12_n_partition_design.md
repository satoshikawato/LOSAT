# S12 — LOSATN単一仕事を分割する設計と隔離検証

種別: 拡張／設計・非production試験のみ  
前提: S05等で単一仕事が障害と確認され、所有者が本sessionを選択

## このセッションの目的

LOSATNの単一query×単一subjectなど、既存subject/chunk並列化で独立仕事が足りないケースを扱う。検索の意味を維持した新規分割が成立するかを設計と最小の隔離試験で判定する。本sessionはproductionの分割条件を変更しない。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と境界

初期scopeのN-copy/N-pack/既存MT改善後にも、指定された実ワークロードで単一ペアlatencyの問題が残ることを確認する。複数pair throughputを別目標として記録する。所有者がこのsessionを起動していない場合は、S05から自動で着手しない。

`blastn/blast_engine/run.rs`、diagonal/lookup、chunk merge、traceback/purge、NCBIの対応call pathを読む。[R09](../SOURCES.md#r09) 新しいquery/strand/range分割を互いに同じものとして扱わない。

## 設計する内容

1. 実際の状態依存を表にする。query context、diagonal table、two-hit履歴、seed順序、extension終端、mask boundary、統計母集団、HSP採用、同点順序、subject ID/座標を含める。
2. queryごとの独立化、seed列挙だけの並列化とordered処理、既存意味単位を保った別分割などを比較する。query 1本にquery-record分割を提示して解決としない。
3. 各候補について必要な境界state、halo/overlapが十分である根拠、重複除去、安定した結果index、ordered reduction、最大memory、残存serial部分を示す。固定長overlapを置けば全てのX-drop/tracebackが安全とは仮定しない。
4. 出力byte一致に加え、中間候補集合と採用順の一致を確認する方法を定義する。単に全件sortして最終行数が同じという証明にしない。
5. 必要な場合だけ、既存test領域に隔離した最小の比較prototypeを作る。production entrypoint、既定CLI、既定chunk長、出力契約を変更しない。source対応を説明できない案は実装しない。

## 判定と成果

GOは、指定ケースに有効な独立単位があり、境界state/統計/順序を保存でき、計算重複とreduction費を含めても利益が見込める設計が存在することを意味する。性能達成やproduction採用の意味ではない。NO_GOは具体的な依存・費用・証拠不足を示す。

S12 report内に設計決定、代替案、必要テスト、危険境界、prototype結果、production変更がないことを記録する。GOの場合のみ、所有者承認用の次実装sessionのscopeと完了条件を同reportに提示する。新たな汎用executorや多数の補助文書を作らない。実装を進める判断が出た後は、その差分についてS11相当の独立検証が必要である。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
