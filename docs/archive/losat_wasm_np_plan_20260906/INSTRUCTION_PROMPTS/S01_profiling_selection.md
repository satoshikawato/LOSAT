# S01 — 最小限のWasm計測と改善候補の選択

種別: 必須／計測の欠落だけ実装可  
前提: S00の環境・fixture・契約記録

## このセッションの目的

N/Pの遅さを起動費、計算、copy/初期化、並列仕事の偏りへ分解し、最初のN改善とP改善をそれぞれ一つ選ぶ。profileの結果が出る前に複数の最適化を実装しない。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 読む対象と根拠

N: `blastn/blast_engine/run.rs`、`alignment/greedy.rs`、`extension.rs`、`sequence_compare.rs`。P: `blastp/blast_engine.rs`と`gapalign.rs`。host: 既存WASI/browser runner。Pの`blastp_timing_env_enabled()`は読解SHAでwasm32に対しfalseを返す。既存環境変数を渡すだけで時間が得られると仮定しない。[R06〜R11](../SOURCES.md)

## 実施手順

1. S00で選んだcommandからproduction callerを辿る。似た名称の未使用helper、test-only経路、古い実装を最適化対象から除外する。各ケースのtask/feature/branchを記録する。
2. 既存のstage timerとdiagnosticsを再利用する。不足がある場合だけ、対応targetで安全な時計を段階・batch境界へ追加する。`wasm32-wasip1`と`unknown-unknown`で時計の供給が異なる可能性を確認し、後者のclockを当然利用可能とは扱わない。計測無効時にinner loopへ時計・log・重いatomicを追加しない。
3. Pの既存時間区分（準備、lookup、scan/ungapped、prelim gapped、Kappa、traceback/output）を使用する。Nは実際に必要なprepare/search/traceback程度の少数区分に留める。入れ子time、worker時間合計、wall timeを混ぜない。
4. N-copyのallocation/copy bytes、N-packの反復回数、N-MTの実仕事数を対象ケースで確認する。Pはdiag fill量、scratch生成数、redo計算数、実際にearly条件で不使用となった数、heap判定数を確認する。全allocatorの置換や全関数へのtimer挿入はしない。
5. threaded artifact n1/n2/n4/n8とserial artifactを分ける。requested/effective threads、実spawn、stageごとのuse_parallel/serial_reason、work item数を記録する。runnerへ環境変数が本当に渡るかも確認する。
6. profiling有効版は原因調査用、性能採用値は計測を外した通常releaseで取る。同じ結果が出ることを確認する。native/browserの未測定値を推定しない。
7. 各候補の到達、負担割合またはabsolute時間、変更の危険性、見込み上限を報告する。時間割合pと局所speedup sから`1/((1-p)+p/s)`を理論的な上限の目安として使えるが、実測改善として記録しない。

## 選択するもの

S02/S03/S04からNの第一候補を最大一つ選ぶ。Nの並列仕事が問題ならS05も候補になる。PはS06/S07/S08/S09から第一候補を最大一つ選ぶ。候補が弱ければSKIPPEDにする。起動/転送が主因ならS10を先行する。性能と無関係な大規模アーキテクチャ変更を改善案にしない。

## 検証と完了

計測のon/offでraw outputと終了状態が同じで、選択理由が数値とcall pathに結び付けばCOMPLETE。時計が供給できず信頼できるsamplingもない場合は対象profileをBLOCKEDとする。計測基盤の拡張だけでsessionを繰り返さず、確認できた最小の実験へ進める資料を残す。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
