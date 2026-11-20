# S02 — LOSATN greedy tracebackの往復コピー削減

種別: 条件付き／局所実装  
前提: S01でcopy経路の到達と有意な負担を確認

## このセッションの目的

LOSATNのgreedy tracebackが既存poolのデータを一時Vecへ複製し、後でpoolへ書き戻す費用を減らす。検索・tracebackの意味、行の更新順、保持されるセル値を変えない。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

`LOSAT/src/algorithm/blastn/alignment/greedy.rs`の`GreedyNonAffineMem::alloc_traceback_row()`、`NonAffineGreedyRow`、`persist()`、そのproduction callerを読む。読解SHAでは`to_vec()`と`copy_from_slice()`の往復があるが、実際に使われることをS01 evidenceで確認する。[R06](../SOURCES.md#r06)

影響範囲はこのrow/pool表現と直近のcaller・テストに限定する。別のgreedy変種やgapped DPを同時に全面改修しない。到達しない、費用が微小、既に改善済みならSKIPPED。

## 実施手順

1. NCBIの該当greedy memory確保・refresh・traceback読出しを確認し、Rustのrow所有者、base row、pool row、old cellのread-before-writeを対応付ける。初期化が未証明なら古い値をゼロ化して隠さず、差異/既存問題として別記する。
2. 現在のcopy bytesとallocation回数をcase別に固定する。結果採用とは別に、どのcopyを消すかを明示する。
3. `start/len/origin`等の小さいrange表現でpool内の値を直接扱う案を最初に評価する。growするVecへの生pointerを持ち続けず、短いborrow、必要時のindex access、事前reserveを使う。安全な借用を成立させるために汎用arenaやunsafe共有を作らない。
4. 一つの最小差分で実装する。rowアクセスの追加branchやindex変換がcopy削減を相殺しないか測る。既存capacityとhigh-water markの扱いを保持する。

## 必須のfocused検証

通常tracebackに加え、poolを拡大する長いalignment、短いalignmentの反復、長→短→長、base/pool storage切替、空に近い境界、同点・indel、負strandの該当ケースを試験する。同じscratchを再使用したときのedit scriptとraw outputを確認する。必要な境界fixtureは既存を優先する。

NativeとWasm両方の出力・search時間・allocation/copy bytes・peak memoryを比較する。kernel単体だけ速くても検索全体が改善しなければ速度目的の採用にはしない。棄却する場合は本sessionの差分のみ外す。

## 採用判定

S00の性能/非回帰基準に合格し、copy削減が確認でき、表現と所有権が旧版より複雑になりすぎていない場合にACCEPTED。borrow成立のために探索順やtraceback結果を変える必要があるならREJECTED。判断不能ならINCONCLUSIVE。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
