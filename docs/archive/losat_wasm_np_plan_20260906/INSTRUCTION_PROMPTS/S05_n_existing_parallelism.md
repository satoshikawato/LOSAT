# S05 — LOSATN既存並列仕事の効率化

種別: 条件付き／スケジューリング限定  
前提: S00/S01、実際に複数subjectまたは既存chunkがある入力

## このセッションの目的

LOSATNの既存subject/chunk単位を維持して、Wasmの並列化判定、負荷分散、起動費用を改善する。新しい配列分割を導入せず、実際に独立した仕事がある場合だけ並列化する。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

`blastn/blast_engine/run.rs`の`blastn_wasi_parallel_decision()`、chunk見積もり、pool構築、既存reducerを読む。最大chunk長5,000,000塩基とWasm閾値の意味を確認する。[R09](../SOURCES.md#r09)

入力が単一仕事なら、本sessionの目的は「増やせない理由の確認」とguardの維持で完了してよい。thresholdを下げただけで新しい独立仕事が生まれると扱わない。

## 実施手順

1. 指定threads、実効threads、実spawn、subject/chunk数、worker別の仕事数・時間、pool起動時間を対象ケースで確定する。Rust診断とrunnerの観測を突き合わせる。
2. N-MT-many、N-one-job、小さい複数subject、長さ/負担が偏るsubjectの対照を用意する。単にCPU使用率が高いことを成功指標にしない。
3. 並列仕事は十分だがgateが過剰、仕事は分散したが細かすぎる、長い一仕事が支配、起動費が支配のどれかに分類する。原因に対応する変更を一つだけ選ぶ。
4. 既存work unitをbatchにまとめる、単純なwork量判定を修正する、結果indexを保った割当を変える範囲で比較する。query/subjectを物理分割せず、chunk境界・overlap・統計母集団を変更しない。
5. ordered reduction、common-endpoint purge、traceback、最終hitlistの既存順序を保つ。完了順に結果をappendする設計へ変更しない。

## 必須のfocused検証

単一仕事でworkerを無駄に作らないこと、多数仕事で実際に並列計算すること、chunk境界を跨ぐhitの一致、負strandと同点HSP順序を確認する。同じthreaded artifactでn1/n2/n4/n8、serial artifact、Nativeを比較する。要求threadsが大きい場合でも登録CPU/memory予算を超えない。

## 採用判定と引き継ぎ

対象でlatencyが改善し、小さいguard・Nativeが基準内、順序が一致する場合ACCEPTED。仕事不足が原因ならSKIPPEDとしてS12の設計候補へ渡す。単一ペアの課題を複数ペアthroughputの数値で解決済みにしない。新規executorや範囲partitionの実装は本sessionに含めない。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
