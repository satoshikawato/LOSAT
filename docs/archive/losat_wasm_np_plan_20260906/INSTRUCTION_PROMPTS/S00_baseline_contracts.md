# S00 — base・互換性契約・fixture・実行予算の固定

種別: 必須／原則read-only  
前提: 作業checkoutと参照資料へのアクセス

## このセッションの目的

LOSATN/LOSATPのWasm性能改善を再現可能にするため、対象コード、互換性の権威、少数の代表fixture、測定条件を固定する。ここでは高速化実装を行わない。性能改善の合格を宣言するのではなく、後続sessionが実際に走らせられる出発点を作る。


## 共通の着手条件

対象はRust実装のNCBI BLAST互換検索LOSATである。LOSATNは`blastn`、LOSATPは`blastp`を指す。WasmをNativeに近づけるための性能改善であり、検索感度や結果を変更する作業ではない。コード読解基準は`7db9bb0060e4e057f9f50807bf9edc2362f20133`。作業時のHEADが異なる場合は、差分を調べて現在の実装を権威として扱う。

まず対象repositoryの`AGENTS.md`、対象階層の局所指示、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`とその`references/`を読む。このパッケージの[総合計画](../MASTER_PLAN.md)、[進捗台帳](../STATUS.md)、本sessionに必要な既存reportを読む。一般基準を別文書へ複製せず総合計画を参照する。旧計画は実装済み項目を含み得るため、source/call pathを確認する。

`PLAN_DIR`は本パッケージのroot、`REPO_ROOT`は対象Git checkout、`WORK_DIR`はS00が登録した永続作業領域である。S00以外ではSTATUSに記録した実在pathを使う。前提が不足する場合は推測で埋めず、実施可能な調査と不足資料をreportに保存する。該当session以外へ自動で作業範囲を広げない。

## 実装・検証の共通制約

NCBIソースの意味、既存のraw output契約、scoring/統計/順序/エラーを維持する。現行AGENTSのGate A/Bと登録済み例外を確認し、platform-local NCBI出力でgoldenを更新しない。NCBIは検証oracleのみで、runtime/build/FFI/fallback依存は禁止。変更箇所のNCBI参照コメントは実在sourceから作り、出力差を比較器の緩和で隠さない。

一つの根本仮説と一つの局所差分に絞る。既存関数・構造体・sliceを優先し、汎用executor/arena/cache/trait体系を作らない。測定だけのセッションで不要なruntime変更をしない。手元の未commit変更・旧証拠を保護し、remote push/merge/releaseは行わない。

このsessionに明記したfocused試験と、候補完成後の該当互換性・性能試験を実行する。候補による結果差を確認したら性能試験の拡大を止め、同じ意味違反の解消に集中する。計測は同条件のbase/candidate、通常flags、同一sinkで行う。parallelの伸びは同一threaded artifactのn1/n2/n4/n8で測る。起動・準備を計測外へ移しただけの改善は採用しない。未実行をPASSと書かない。


## このセッション固有の初期化

S00だけはSTATUSのpathが未設定でよい。現在のcheckoutを`git rev-parse --show-toplevel`で確認し、planの実在rootを特定する。証拠はrepositoryの既存出力や一時的な`/tmp`を上書きせず、永続`WORK_DIR`へ保存する。適切な領域があればそこで`reports/`と`evidence/`を作る。STATUSに実在absolute path、移設時の相対関係、branch、HEAD、dirty状態を記録する。標準的な場所が使えない場合は、その理由をBLOCKEDとして示す。

## 読む対象

`Cargo.toml`、`Cargo.lock`、`.cargo/config.toml`、Native/Wasm buildとrunner、既存certification/manifest、`benchmarks/v0.1.0/README.md`と`plot_data.json`を読む。ソースのrootはrepository内の`LOSAT/`である。N/P各engineのentrypointと対象APIの対応条件を確認する。[R02〜R05・R09〜R11](../SOURCES.md)

## 実施手順

1. 作業branchを勝手にresetしない。コード読解SHAと実行HEADの差を確認し、初期baseを固定する。未commit差分があるなら所有範囲を確認してhash付きで記録し、他人の作業をcandidateへ混ぜない。toolchain、target、CPU/OS/WSL、Node/ブラウザ、NCBI版・source識別子を記録する。
2. 互換性の権威をケース別に固定する。raw goldenの実在path/hashと登録契約を読む。Native Gate A/Bを区別し、unknown platform fingerprintを許容しない。生のgoldenが利用できなければ、その検証はBLOCKEDとし、base/candidate比較だけで認定しない。
3. 現在サポートされるNのblastn/megablast、Pのoption・outfmt範囲を調べる。最初は既存Pese/Mj、Sakai/MG1655、P pairwise等の実在fixtureを探索の起点とするが、名前だけからpathや特徴を確定しない。
4. 総合計画§9.1の役割を満たす8〜12件以内を目安に選ぶ。実際の経路到達はS01で確定する。N-one-job、N-MT-many、P-one-query/P-multi-query、小さいguardと重い候補を区別する。境界unit testは既存を優先する。
5. 実行する既存scriptのargv、出力path、clean処理、環境変数の受け渡しを読む。古い出力を破壊するscriptをそのまま実行しない。実在するbuild target、lib/bin衝突、配布artifactを確認する。架空のrunner名や未対応flagをcommandに書かない。
6. 小さいN/Pでfresh buildとsmoke parityを実行する。既存不一致があれば最初の差を記録し、無関係なoptimizationを開始しない。初期wall timeを観察してtimeoutと最大反復数を決める。このsessionで全プログラム・全platformの巨大計測はしない。
7. 各測定系列のcold/warm定義、output sink、thread count、compile flags、Native対照、memory指標、採用/非回帰基準、反復予算を固定する。候補を見た後に基準を緩めない。

## 出力する最小の契約

S00 report内に、環境表、fixture表、実在command、出力/golden参照、実行予算をまとめる。既存manifestが同じ情報を持つなら参照と不足分だけを保存する。collectorを新設せず、必要最小限のscript変更案がある場合はS01へ渡す。

fixture表の各行はcase ID、input hash、program/task、完全argv、contract、expected hash/所在、役割、未確認経路を持つ。commandでは実在pathを使用し、利用者のPC固有pathを文書から無条件コピーしない。

## 完了・停止の判定

N/Pを少なくとも一つずつ再実行でき、authorityと実行予算が確認できればCOMPLETE。browserや一部大型fixtureの不足は個別にBLOCKED/NOT_RUNを記録できる。基礎となる検索を実行できない、raw契約が不明、checkoutが混在してbaseを定義できない場合は対象範囲をBLOCKEDにする。足りない資料を新しいgoldenや推測値で代替しない。


## セッション終了時の提出

[SESSION_REPORTテンプレート](../templates/SESSION_REPORT.md)を使い、`WORK_DIR/reports/SESSION_ID.md`へ保存する（`SESSION_ID`はこのファイルのSxx）。base/candidateのSHAまたはpatch hash、artifact hash、実行command、生証拠path、結果一致、timingとmemory、採用/棄却理由、残課題を記録する。変更なしの判断にも根拠を残す。

STATUSのこのsession行と次の推奨sessionだけを更新する。runtime候補はACCEPTED/REJECTED、適用なしはSKIPPED、資料不足はBLOCKED、測定不能はINCONCLUSIVE、非実装作業の完了はCOMPLETEを使う。棄却した差分を統合baseへ残さない。次の担当者には具体的なpathと再実行commandを渡し、過去のチャットや「いつもの設定」を参照させない。
