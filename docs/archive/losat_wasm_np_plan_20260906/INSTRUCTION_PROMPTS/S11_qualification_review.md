# S11 — 統合候補の独立検証と最終結論

種別: 必須／read-onlyレビューを主とする  
前提: 候補の採否・artifact・S10の実施状態が明確

## このセッションの目的

採用済みN/P変更だけを含む統合候補について、互換性、安全性、性能、アーキテクチャとワークフローの簡潔さを独立に確認する。未達や未測定も含む条件付きの結論を作り、リリース操作は行わない。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## レビュー独立性と入力

可能なら実装担当と別のreviewerが実行する。利用可能な`ncbi_parity_auditor`または別read-only sessionを使用する。同じ担当の再確認だけの場合はSELF_REVIEWと記録し、独立レビューが完了したとは書かない。[R02](../SOURCES.md#r02) [R03](../SOURCES.md#r03)

STATUS、S00の初期base/契約/fixture、各採否report、全採用patch/commit、統合candidate、実build flagsとartifact hashを読む。証拠の実在を確認し、reportの文章だけで合格としない。

## 実施手順

1. 棄却/未確認の差分がcandidateへ紛れていないかを確認する。baseと統合candidateを固定し、未commit変更がある場合はpatch/treeのidentityを確定する。
2. 変更されたNCBI意味と参照コメントを照合する。Gate A/B、raw golden、登録済み例外の範囲を確認し、expectedの更新や比較器緩和がないことを確認する。
3. owned scratch、reset、pool range、SIMD tail、generation overflow、Kappa ordered replay、error契約をread-onlyで検査する。新しい汎用型・依存・設定が効果に必要だったかを確認する。
4. focused regressionの通過後に、対象N/Pの既存certification範囲と必要な共通部回帰を実行する。未知のfingerprintや原本不足はBLOCKEDでありPASSではない。共有箇所を変更した場合のみTBLASTX回帰を含める。
5. 計測診断を外した通常releaseで、初期baseと統合candidateを同一条件で測る。各sessionの改善率を掛け合わせず最終値を実測する。Native/serial Wasm/threaded n1/n2/n4/n8とbrowserの実施範囲を明記する。
6. performance sampleとoutput hashを対応付ける。失敗runを黙って除外しない。latency/throughput/memory、cold/warm、platform差を別々に評価する。
7. 原則への実質的な準拠を確認する。準備/計算/reduction/hostの責任分離、最小の共通化、変更の局所性、WIP、単一進捗台帳、不要実験の棄却が成立しているかを見る。trait数や文書量を評価指標にしない。

## 最終成果物

`WORK_DIR/reports/S11.md`を最終reportとする。結論、検証範囲表、採用/棄却一覧、固定identity、初期base対比の数値、Native比、actual work分散、メモリ、安全性/互換性証拠、browser未実施一覧、再実行command、rollback、残余課題を含める。

「全体が速くなった」「Wasm/Native差が縮まった」「単一ペアMTが改善した」「複数ペアthroughputが改善した」を別々に判定する。N単一仕事のままならその事実を明記する。目標達成は測ったcase/targetの範囲でのみ宣言する。

## 不合格・完了の扱い

新しい不具合を見つけた場合、このレビューsessionで無関係な大改修をしない。blocking項目を該当実装sessionへ返し、修正版を再検証する。環境/権威不足と候補の実不合格を区別する。独立レビュー未了はREVIEW_PENDING、未測定はNOT_RUNとする。S11の完了はremote merge/release権限を与えない。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
