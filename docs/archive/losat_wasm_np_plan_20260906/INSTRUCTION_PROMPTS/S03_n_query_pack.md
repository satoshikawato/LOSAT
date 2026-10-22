# S03 — LOSATN query再packの削減

種別: 条件付き／局所実装  
前提: S01で既存4塩基extensionの反復費を確認

## このセッションの目的

LOSATNの既存4塩基単位ungapped extensionと同じ演算を保ち、queryを繰り返しpackする費用を減らす。前計算の時間と追加メモリまで含めて利益を判断する。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## 開始条件と変更範囲

`blastn/extension.rs`の`extend_hit_ungapped_approx_ncbi()`と`extend_hit_ungapped_exact_ncbi()`、query準備の直近caller、encodingを読む。既存approximate経路はNCBI互換処理の名称であり、新しい近似heuristicを追加する許可ではない。[R07](../SOURCES.md#r07)

既存のpacked queryまたは同等表現が再利用できるかを先に確認する。S02のscratch変更とは独立に評価する。extension回数が少なく回収不能、または異なる経路しか使われないならSKIPPED。

## 実施手順

1. queryのBLASTNA値、mask、sentinel、left/right offset、phase、既存pack式の型とbit演算結果を確認する。曖昧文字をACGT用2-bit値へ勝手に変換しない。
2. 代表fixtureでpack回数、query長、同じqueryの再利用回数を確認する。前計算を導入する場合のbyte数と初期化時間を見積もり、measurementに含める。
3. 既存表現の再利用を優先し、不足ならrun内に閉じる最小のphase別/offset別表現を一案だけ比較する。全入力に巨大cacheを作らず、単純な採用条件と既存fallbackで足りるか判断する。cross-run cacheや新CLI optionは導入しない。
4. 既存のscore_table参照順、X-drop評価、exact再計算への切替条件・開始位置を保つ。ブロックをまとめるために最初の停止位置を飛び越えない。

## 必須のfocused検証

query長0〜8程度の境界、phase 0/1/2/3、左右extension、sequence末端、mask内外、曖昧文字、負strand、exact切替閾値の直前/一致/直後を試験する。既存pack演算と前計算値は、サポートするencoding値の組合せを系統的に比較する。出力の順序を並べ替えて比較しない。

少数extensionと多数extensionの両方をtiming対象にする。全体prepare+search時間、準備済みの場合の時間、peak memoryを分ける。既存Native経路も悪化させない。

## 採用判定

前計算込みでS00の基準を満たす場合のみACCEPTED。計測外へ処理を動かした場合、追加メモリが登録予算を超える場合、曖昧文字で意味を保てない場合はREJECTED。適応条件の調整は少数候補に留め、fixture別の特例を入れない。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
