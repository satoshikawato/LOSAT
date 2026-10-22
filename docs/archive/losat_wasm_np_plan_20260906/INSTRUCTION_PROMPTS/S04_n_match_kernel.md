# S04 — LOSATN一致延伸カーネルの局所最適化

種別: 条件付き／局所実装  
前提: S01で実際のcallerとhotspotを確認

## このセッションの目的

LOSATNの完全一致検査または最初の不一致までの延伸を一つ選び、Wasm向けに高速化する。SIMDを使うこと自体ではなく、同じ結果を少ない時間で得ることを目的とする。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

`blastn/sequence_compare.rs`と呼び出すgreedy/word処理を読む。読解時にNative専用dispatchやscalar mismatch loopがあることは、当該helperがproductionで重い証明ではない。[R08](../SOURCES.md#r08)

対象helperを一つに限定する。実際には別のinline比較が支配的なら、そのcall pathを報告して対象を明記し直す。未使用helperの大量SIMD化、PやTBLASTXへの展開はしない。

## 実施手順

1. 元の比較predicateを確定する。単なるbyte equalityなのか、通常塩基のみの継続条件や曖昧文字停止条件を含むかをNCBIとcallerから確認する。
2. 既存release Wasmの生成コードを、利用可能なdisassembler/profilerで確認する。scalar sourceでも自動ベクトル化され得ることを考慮する。ツールがなければ断定を避け、実時間で小さいA/Bを行う。
3. 基準scalarと、scalar unrollまたは明示SIMDの少数候補を比較する。`v128`で16-byteを比較する場合、最初の不一致位置、reverse時の位置変換、全部一致のmask処理を厳密に保つ。
4. tailは範囲内loadと安全なscalar処理を優先する。28-byteのsliceから32-byteを読む方法をコピーしない。maskで不要laneを無視しても範囲外loadの安全性は回復しない。
5. target分岐をhelper入口へ閉じる。S03の前計算設計や新しいdispatch frameworkを同時に変更しない。

## 必須のfocused検証

長さ0/1/15/16/17/28/31/32/33、各laneの不一致、全一致、最初/最後の不一致、非整列開始、forward/reverse、allocation末端、曖昧文字を試す。最適化helperと契約を満たすscalar基準の返値を比較し、検索raw outputも確認する。

短い一致列と長い一致列の両方で測定する。microbenchmarkと検索全体を分け、deploymentと同じSIMD flagsで測る。S00の通常runtimeを特殊JIT設定に置換しない。

## 採用判定

kernelと対象検索全体の改善、全境界の一致、Native/他targetの非回帰が成立したときACCEPTED。SIMD版の方が遅い場合は残さない。速度差が微小ならINCONCLUSIVEまたはSKIPPEDとして、他の実ボトルネックへ判断を返す。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
