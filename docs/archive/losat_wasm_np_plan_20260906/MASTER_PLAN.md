# LOSATN / LOSATP Wasm性能改善・総合計画書

版: 1.0  
作成日: 2026-09-06（日本時間）  
対象: `satoshikawato/LOSAT` の `blastn`（LOSATN）と `blastp`（LOSATP）  
コード読解の基準: `7db9bb0060e4e057f9f50807bf9edc2362f20133`  
文書の状態: **実行前の計画。改善実装、現環境でのベンチマーク、合格認定は未実施。**

## 1. 目的と対象読者

LOSATは、NCBI BLAST互換の配列相同性検索をRustで実装するソフトウェアである。本計画のLOSATNは`blastn`サブコマンド、LOSATPは`blastp`サブコマンドを指す。独立した実行ファイルがそれぞれ存在するという意味ではない。対象は、既にサポートされている検索条件でのNative実行とWasm実行である。[R01](SOURCES.md#r01) [R02](SOURCES.md#r02)

目標は、結果と互換性を変えず、Wasmの単一スレッド待ち時間とマルチスレッド実行効率を改善することである。LOSATNでは不要なコピーと反復演算を優先して削減する。LOSATPではgapped alignmentとKappa redoの計算・状態管理を調べ、並列化に伴う初期化と不要な先行計算を減らす。

読者は、過去の議論を知らない開発者、コード実行エージェント、レビュー担当者である。着手に必要な文書はこのパッケージ、対象リポジトリの現行ルール、実際の検証資料だけであり、チャット履歴や個人の記憶は権威として使わない。

### 用語と比較対象

本書のNative対照は、原則として**同じLOSATコードを通常のCPU向けにビルドした実行物**であり、NCBI BLAST+そのものではない。NCBIは意味仕様・互換性検証の対照である。WasmはWebAssembly実行物、WASIはその実行物へシステム機能を供給するインターフェースであり、ブラウザで動かす場合は対応するhost実装を確認する。[R02](SOURCES.md#r02) [R19](SOURCES.md#r19)

| 用語 | 本書での意味 |
|---|---|
| fixture | 入力・全検索条件・期待結果を固定した検証ケース |
| HSP | 検索で得られる高スコアの局所alignment候補 |
| Kappa redo | タンパク質検索のcomposition-based処理に伴う再alignment・後処理 |
| scratch | 計算途中で再使用する作業用メモリ。保持容量と意味状態は別に管理する |
| heap | 結果の採用・順位管理に使うデータ構造。メモリ領域一般の意味ではない |
| golden / raw output | 登録契約が指定する比較用原本／整形し直していない出力byte列 |
| artifact | 実際に測定・配布するNative/Wasm実行物。sourceと別にhashを記録する |
| n1 / n8 | 指定スレッド数。実際のworker数や有効な仕事数と同じとは限らない |
| latency / throughput | 一つの検索が完了するまでの時間／複数検索を処理する量 |

用語はNCBIと既存LOSATの実装上の役割に沿って使用する。Kappa・heap等の具体的な状態と順序は実コードで確認する。[R10](SOURCES.md#r10)

## 2. 成果と成功の定義

### 2.1 三つの成果を別々に測る

| 成果 | 指標 | 誤った代替指標 |
|---|---|---|
| LOSATそのものの高速化 | 同じtarget・runtime・fixtureの変更前後時間 | Wasm/Native比だけを見る |
| WasmとNativeの差の縮小 | 同一候補コード・同条件のWasm時間÷Native時間 | 古い遅いNativeと新Wasmの比較 |
| マルチスレッド効率向上 | 同一threaded artifactのn1/n2/n4/n8 | serial専用artifactとthreaded n8だけの比較 |

共通Rustコードを改善するとNativeも速くなり得る。双方が速くなったが比率は変わらない結果を「差が縮小した」と記録しない。独立した複数ペアの処理量改善と、1ペアのlatency改善も分ける。

### 2.2 開発上の目安と採用ゲート

「Native並み」の最初の到達目安は、計算主体の代表ケースでwarm searchのWasm/Native比をn1で1.20以下、同一nの並列実行で1.30以下とする。**この数値は製品目標案であって予測・保証・既達成値ではない。** Native側の並列効果がない入力では比率だけでなく絶対時間を報告する。目標未達でも証拠のある個別改善は採用できるが、達成を宣言しない。

S00で候補を見る前にケースと採用基準を固定する。初期案は、主対象のWasm search中央値を5%以上改善、guardケースとNativeに`max(3%, 10 ms/検索)`を超える再現性のある悪化なし、peak memoryの10%超増加なしとする。メモリ増加を許容する設計は、変更前に別途上限と理由を登録する。軽量ケースでは絶対差を併記し、測定ノイズが閾値より大きければINCONCLUSIVEとする。メモリ削減のみを目的に変更を採用する場合は、速度向上とは別分類にする。

warm searchは、同一runtime内で必要な初期化条件を揃えた、検索呼び出し開始から結果バッファ完成までの時間とする。検索に必要な前計算を計測外へ移しただけの改善を認めない。準備済み入力を使う別APIは別の利用形態として測る。

## 3. 読解で確認した事実と、未検証の仮説

2026-09-06に公開mainのSHAを確認した。以下はそのSHAまたは明示した保存済み測定の情報である。実行開始時には手元のHEADと照合する。[R01](SOURCES.md#r01)

| 項目 | 確認した事実 | 性能について未確認の点 |
|---|---|---|
| リリースビルド | LTO、`codegen-units=1`、`opt-level=3`、Wasm `+simd128`設定が存在 | 実際の配布artifactに適用されたflagsと効果 |
| N greedy | `alloc_traceback_row()`にpool領域の`to_vec()`、rowの`persist()`にコピー戻しがある | 対象fixtureでの到達頻度、allocation/copy時間 |
| N ungapped | 既存4塩基extensionでqueryを内側ループでpackする | 前計算を含めた収支と曖昧塩基の扱い |
| N比較helper | 読解したdispatchにx86/AArch64経路、`find_first_mismatch_ex()`にscalar loopがある | production caller、LLVMの自動ベクトル化、Wasm専用実装の優位性 |
| N並列化 | subject/chunk単位。最大chunk長5,000,000塩基。Wasm用判定がある | runtime別の仕事分布と適切な閾値 |
| P計測 | `blastp_timing_env_enabled()`はwasm32でfalse | 計測を安全に有効化したときの各段階の負担 |
| P pool | lazyなper-run共通poolとWasm `use_current_thread()`が既にある | 複数runの寿命、実際のworker起動費用 |
| P DP | BLOSUM62の特殊化、adjusted matrixの行参照が既にある | 残存するDP/tracebackホットループ |
| P prelim | parallel subject処理前にdiagonal配列をfillする | 初期化割合、sparse resetの損益 |
| P Kappa | 単一queryのmatch並列経路にmatchごとのworkspace/scratch作成と全件先行redoがある | 経路の適用範囲、不要計算量、batch化の負荷分散損失 |

コード根拠: [R04](SOURCES.md#r04)、[R06](SOURCES.md#r06)〜[R11](SOURCES.md#r11)。**表の右列を測る前に、左列だけを根拠に改修を決めない。** 古い性能分析や計画書は探索の手掛かりとし、実装済みの項目を再実装しない。

### 3.1 保存済み数値は参考資料である

既存snapshotはコード`af3e2ea837afdb8a00cf19920f68be4f0bf3bfb5`、同一物理マシンのWSL2 Ubuntu 24.04.3 / Intel Core i9-14900HXで各5回のwall timeを保存している。以下はその中央値（秒）で、新規計測ではない。[R05](SOURCES.md#r05)

| ケース | Native n1 | Wasm serial | Native n8 | Wasm n8指定 | Wasm n8の分類 |
|---|---:|---:|---:|---:|---|
| PesePMNV.MjPMNV.task_blastn | 0.302596 | 0.623383 | 0.350770 | 0.498079 | 実効直列 |
| Sakai.MG1655.megablast | 0.953918 | 1.543504 | 0.991877 | 1.334621 | 実効直列 |
| pairwise_default_serial（blastp） | 0.654475 | 1.090408 | 0.244928 | 0.701468 | 実効並列 |

このsnapshotには異なる互換性契約のケースが含まれる。Sakaiの`SOURCE_UNDETERMINED_ACCEPTED`を無条件のNCBI完全一致ケースに読み替えない。短い3ケースだけから重いN/P全体のボトルネックを断定しない。snapshotファイル、旧ログ、旧図は上書きしない。

## 4. スコープと製品上の決定

### 4.1 今回実行する範囲

既存のLOSATN `blastn` / `megablast`、既存のLOSATP対応条件、Native・serial WASI・threaded WASI、およびgbdrawの既存ブラウザ呼び出しを対象とする。実際に対応しているoption・outfmtの範囲はS00でmanifestとコードから固定する。エンジンの名称だけから対応範囲を増やさない。

検索結果、感度、masking、スコア、統計、HSP選択・順序・出力、エラー契約を維持する。scratch寿命、メモリ配置、条件を保つカーネル差し替え、既存独立仕事のスケジューリングは変更候補とする。

### 4.2 今回自動的には実施しない範囲

全面rewrite、新しいcrate階層、汎用executor、全検索を包むtrait体系、GPU化、言語変更、未対応task追加、NCBI以外の探索heuristic導入、ブラウザAPI全面変更は対象外。依存更新、PGO、allocator総入替、コンパイラ変更をアルゴリズム最適化と同時に行わない。

TBLASTXを新しい最適化対象にしない。共通コードやworkerを変更した場合の回帰確認だけを行う。`tblastx-wasm-scalar`はN/P用の切り替えではない。

Nの新規query/range/seed分割は、既存subject/chunk並列化とは別の拡張判断とする。S12は設計と小さな隔離実験までを扱い、既定のproduction経路を変更しない。単一ペアのlatency目標が残る場合は「未達」と記録し、複数ペアthroughputで代用しない。

## 5. 互換性の権威と変更禁止条件

### 5.1 権威の順序

実行環境の上位指示、リポジトリの現行`AGENTS.md`と適用される局所指示、登録済み製品決定・certification・期待出力、この計画の順に確認する。衝突は勝手に解釈せず、対象ケースをBLOCKEDにして原文と不足物を記録する。この計画は既存の例外を増やす権限を与えない。[R02](SOURCES.md#r02) [R03](SOURCES.md#r03)

NCBIのC/C++ソースは意味仕様の権威である。参照先の版、ファイル、関数、行範囲、呼び出し順を確認し、変更箇所にはリポジトリ規約どおり参照コメントを付す。行番号やローカルパスを記憶から作らない。Rust固有の所有権・スケジューリング変更は、対応するNCBI上の意味・寿命との関係を記述する。無関係なNCBI断片を貼って根拠にしない。

### 5.2 プラットフォーム差を誤って「修正」しない

現行`AGENTS.md`の`PD-NCBI-PLATFORM-VARIANCE`では、Native Gate Aは固定PR5のraw bytesに対するLOSATの一致、Gate Bは指定された公式NCBI検索の登録済みplatform fingerprintとの一致である。実行環境のNCBI出力でLOSATのgoldenを置換しない。未知のfingerprintは合格ではない。[R02](SOURCES.md#r02)

N/Pの各ケースには、その実際の登録契約を付ける。TBLASTXのlocal-subject非標準遺伝暗号例外はN/Pへ転用しない。raw goldenがない場合、候補とbaseの一致は回帰証拠にはなるが、既存認定契約への合格の代わりにはならない。normalized alignmentデータからraw textのgoldenを捏造しない。

### 5.3 絶対にしない変更

NCBI実行ファイル・ライブラリ・FFI・subprocessをruntime/build/fallbackに導入しない。テストoracleとしての外部実行だけを許可する。CBS/SEG/DUSTを切る、seed数を減らす、X-dropを変える、top-kを縮める、float精度やtie-breakを変えることを高速化として採用しない。並列結果を完了順でheapへ投入しない。結果差を隠すsort・正規化・許容誤差を比較器に追加しない。

候補に起因する最初の差異が出たら、性能測定の拡大を止める。同じ意味違反を生む関連箇所の修正とfocused検証に集中し、無関係な最適化や試験を重ねない。既存不具合の修正が必要なら、性能差分と分離した記録・パッチにする。

## 6. 実装アーキテクチャ

### 6.1 責任の境界

新しい物理ディレクトリ構造を強制するのではなく、既存実装内で次の責任を明確にする。

```text
host/runtime adapter
  起動・入出力・worker/instance寿命・能力検出
          ↓ 値としての入力/検索条件/実行上限
NまたはPの既存run coordinator
  準備・段階の順序・仕事割当・既存reductionの呼び出し
          ↓ 不変入力と専有scratch
検索カーネル
  既存と同じ演算・停止条件・局所結果
          ↓ 安定した元の識別子付き結果
既存の順序依存reducer / formatter
  採用・間引き・統計・順序・出力
```

JSのWorker、WASI environment、ベンチマーク計測器をDP/extensionの内側へ直接持ち込まない。runtime能力や予算は必要な値だけを渡す。標準の関数・構造体・sliceで十分ならtrait objectやDI frameworkを作らない。

### 6.2 データの寿命

| データ | 原則的な寿命 | 注意点 |
|---|---|---|
| 生配列と安定ID | 入力または検索run | 元のindexと順序を維持 |
| packed/encoded query・subject | まず同一run | masking・strand・option依存を明示 |
| DP/diagonal/traceback scratch | batchまたは検証済みworker寿命 | capacity再利用と意味状態のresetを分ける |
| composition workspace | 再初期化可能な局所処理単位 | adjusted matrix等の別query/subjectへの漏出を防ぐ |
| 採用heapと結果列 | 既存のquery/subject範囲 | 順序依存stateをworker間で競合更新しない |
| コンパイル済みWasm module | hostの既存cache範囲 | instance再利用・pool再利用とは別 |

共有可変scratchをMutexで包んで全workerに配るのではなく、計算単位ごとの専有を優先する。プールを使っていても、内部でcopyやclearをしているかを測る。capacityのhigh-water markが残ることと、検索回数に比例するleakを区別する。

Rayonの`for_each_init`/`map_init`はOS workerごとに厳密に一度だけ初期化する保証として使わない。job分割に伴う初期化回数を数える。厳密な寿命が必要なら小さい明示batchで所有権を固定し、まずその費用を測る。[R14](SOURCES.md#r14)

### 6.3 安全性と意味の置換可能性

sliceの範囲外SIMD loadをtail maskで正当化しない。Vecの再allocationをまたいでpointer/referenceを保持しない。整数幅、overflow、sentinel、曖昧塩基、左右方向、最後の一致位置、X-dropの最初の成立位置を維持する。安全なRustで成立する局所変更を優先する。

scalar/optimized経路は同じ契約を満たすことが前提であり、別の「近似結果」を返す実装ではない。内部診断用の一時的な比較経路は許すが、採用時に不要なfeature/環境変数/重複実装を整理する。恒久的なscalar fallbackは、対応targetとテストoracleとして明確な役割がある場合に限る。

### 6.4 SOLID / KISS / DRY / YAGNIの適用

| 原則 | 実装への適用 | 開発ワークフローへの適用 |
|---|---|---|
| SRP | prepare、compute、ordered reduction、host寿命を混ぜない | 一つの実験、一つの意思決定、一つのreport |
| OCP | 現在必要なtarget分岐を局所化し、既存契約で置換 | 新しい証拠で候補順を変えられる条件付き計画 |
| LSP | scalar/SIMD、serial/parallelで結果・エラー契約を保持 | 合格基準を担当者やセッションごとに変えない |
| ISP | カーネルには必要なslice・設定・scratchだけを渡す | セッションには対象の資料と必要な証拠だけを渡す |
| DIP | 検索意味がJS、WASI、タイマーに依存しない | 記憶ではなく登録契約と再実行可能な証拠に依存 |
| KISS | 既存関数・構造体・明示batchを優先 | WIP=1、同一ファイルへの並行編集を避ける |
| DRY | 同じ意味・同じ寿命の重複だけを共通化 | 共通基準は本書、進捗はSTATUS、結果はsession report |
| YAGNI | 汎用pool/executor/cacheを先回りで新設しない | 適用条件が成立しないセッションを実行しない |

NとPのscratchは意味が異なるため、型が似ているだけで統合しない。原則への準拠をコード量やtrait数で評価しない。削除したcopy量、責任境界、状態の所有者、実測効果、レビュー可能性で評価する。

## 7. LOSATNの実験系列

### 7.1 N-copy: greedy tracebackの往復コピー

S02は`GreedyNonAffineMem`、`NonAffineGreedyRow`等のproduction callerを特定し、pool→Vec→poolのallocation回数とbyte数を測る。到達しないhelperは変更しない。[R06](SOURCES.md#r06)

変更候補は、必要なrowをpool内の安定したrange/indexで表現して直接更新すること。借用関係が複雑化するなら、十分な容量確保後の短いborrow、計算順を保った局所API、またはコピー削減の小さい代案を選ぶ。unsafeな汎用arenaを作らない。古いセル値の保持が結果に影響する可能性は、ゼロ化で隠さず初期化証明と既存契約で確認する。

### 7.2 N-pack: queryの4塩基表現再利用

S03は`extend_hit_ungapped_approx_ncbi()`の既存4塩基単位の演算を対象にする。候補は任意offset/phaseの既存pack式と同一の値をrun内で前計算すること。既に存在するpackedデータを使えるかを先に調べる。[R07](SOURCES.md#r07)

曖昧塩基を通常のACGTへ丸めない。左右extension、offset mod 4、mask/sentinel、exact再計算へ移る閾値と位置を保存する。前計算・追加メモリ・cache missも含めて測る。短いquery、少数extensionでは元経路が適切なら、単純な損益条件で選ぶ。外部設定の追加を前提にしない。

### 7.3 N-match: 完全一致延伸の局所最適化

S04は実際に呼ばれる一致検査/延伸を一つだけ対象にする。Native専用dispatchがあるという理由だけでWasmが必ず遅いとは判断しない。生成Wasmの自動ベクトル化を確認し、scalar/unroll/明示SIMDから少数の比較候補を選ぶ。[R08](SOURCES.md#r08)

16-byte境界、短いtail、reverse方向、buffer末端、曖昧文字の停止条件を重点検証する。広いDP SIMD化やNativeの32-byte overreadの模倣はしない。

### 7.4 N-MT: 既存の独立仕事を適切に使う

S05はsubject/chunkの実仕事数、Wasm閾値、各workerの負担、pool起動費用を確認する。閾値変更と新しい検索分割を同じ差分にしない。[R09](SOURCES.md#r09)

多数subjectまたは実際に複数chunkを持つ入力では、既存の意味上の単位を保ったbatch化・割当調整を比較する。query 1本×subject 1本で通常chunk条件により仕事が1個なら、n8を強制しても独立仕事は増えない。小仕事をserialにする正当な判断を残す。

同一run内の新規query分割・strand分割・座標分割はS12で別途扱う。単純なFASTA分割と出力連結で統計や境界が同じになると仮定しない。

## 8. LOSATPの実験系列

### 8.1 P-ST: gapped/traceback/Kappaの上位一処理

S06はS01のWasm段階別profileを根拠に一つのボトルネックを選ぶ。既存のBLOSUM62特殊化とadjusted matrixの行参照を確認し、未実装として作り直さない。[R10](SOURCES.md#r10) [R11](SOURCES.md#r11)

候補はDPの生存領域、traceback buffer、必要なscore lookup、clone/copy、局所的なscratch初期化の削減。汎用Smith–Watermanへの置換、CBS/SEG省略、数値精度変更はしない。S08/S09の所有範囲と重複した場合、その作業へ振り分けてS06をSKIPPEDとする。

### 8.2 P-reset: parallel preliminaryのdiagonal初期化

S07は`prepare_independent_subject()`のfill負担とscratch生成回数を測る。query集合が大きくsubjectが多数短いケースと、長くアクセスがdenseなケースを対照にする。[R10](SOURCES.md#r10)

候補は明示batchでのscratch再利用、変更済みindexだけのreset、またはgeneration方式のいずれか一つ。generation方式なら全read/writeの有効性判定、世代overflow時のfull reset、初期sentinelとの等価性を確認する。denseケースで追加branchが損になる場合は採用しない。serial経路で既に使われるdiag offsetの意味を崩さない。

### 8.3 P-scratch: Kappa match経路のworkspace再利用

S08は単一query match-redo経路で実際に発生するworkspace作成を対象とする。配列capacity等の再利用可能な部分と、query/subject/adjusted matrix依存の意味状態を分離する。まずrun内のbatchで寿命を閉じ、cross-run cacheへ広げない。[R10](SOURCES.md#r10)

multi-query経路も同じ問題と決めつけない。`map_init`だけでworker数と同数のallocationになると仮定しない。cold-small/warm-many、query/subjectを切り替える連続実行で漏出を検証する。

### 8.4 P-redo: 順序を保つ有界先行計算

S09はmatch-redoが全件を先行計算し、その後にearly-termination/heap判定をreplayする経路を対象とする。`computed`、`early-skipped-before-use`、`heap-rejected`、`retained`を区別する。**最終heapに残らない結果を全て不要計算と数えない。判定のために計算が必要だった結果もある。** [R10](SOURCES.md#r10)

変更候補は、元のlocal-match順の有界windowを計算し、同じ順で既存判定を適用してから次のwindowを発行する方式である。既存`continue`条件を根拠なく`break`へ変更しない。将来のmatchを計算前に安全に除外できる条件、heap更新、統計state、エラー伝播をNCBIと現行契約から証明する。

速いworker順の採用、未計算候補の経験則的cutoff、全queryの共有heapは導入しない。全件必要なケースでは追加barrierが悪化し得るため、early-terminationが効くケースと効かないケースを両方測る。必要なら単純な既存経路との選択に留める。複雑な非同期schedulerを作らない。

## 9. 実験と計測の設計

### 9.1 S00で固定するfixtureの役割

既存データを優先し、存在しないファイル名や入力サイズを作らない。初期の代表suiteは8〜12件程度を上限の目安とし、1ケースに複数の役割を持たせてよい。適用経路の存在を実測で確認してから固定する。

| 役割 | 必要な特徴 |
|---|---|
| N-small | 小さい既存blastnペア、起動費の対照 |
| N-greedy | megablast、greedy/tracebackコピー経路に到達 |
| N-many-extensions | pack費用を観察できるextension数 |
| N-MT-many | 複数subjectまたは本当に複数の既存chunk |
| N-one-job | query 1本×subject 1本、並列化しない判断の対照 |
| P-small | 小さい既存pairwise、serial gateの対照 |
| P-ST-heavy | gapped/Kappaの時間が十分な既存タンパク質集合 |
| P-many-short | 大きなquery集合と短い多数subject、reset費の候補 |
| P-one-query | match-redoへ実際に入る入力 |
| P-multi-query | query-redo、負荷分散とcache漏出の対照 |
| 境界・互換性 | 曖昧塩基、同点、重複HSP、短いtail、負strandなど |

新しいsynthetic入力は最小限の境界テストとして、生成方法・seed・hashを保存する。成果を良く見せるために本番から乖離したsynthetic workloadだけを主対象にしない。既存のサポート範囲にないoptionをfixture作成のために実装しない。

### 9.2 実行モード

Native n1/n2/n4/n8、serial Wasm n1、**同一threaded Wasm artifact**のn1/n2/n4/n8を区別する。ハードウェア予算がn8未満なら能力上限までとし、未測定を明記する。Rust既定releaseを主なNative対照とし、CPU特化ビルドを使うなら別系列としてflagsを示す。Nativeを意図的に遅くして比率を良くしない。

各candidateは直前に採用されたbaseと比較する。最終S11は初期固定baseとも比較する。中間候補の改善率を掛け算して最終値を作らない。棄却された候補をbaseへ混ぜない。

### 9.3 時間境界

`process/worker start → compile → instantiate → input/prepare → search → output transfer → teardown`を区別する。重なって進むphaseの時間を機械的に合計しない。parallel stageのwall time、worker workの合計、process CPU timeは異なる量である。

cold command E2E、ブラウザ初回E2E、同一instance warm search、入力準備再利用ありのwarm searchを別系列にする。V8のtieringがあるため、別processでのwarmupを同一instanceのwarmupと呼ばない。検索の途中で最適化版へ切り替わると仮定しない。[R17](SOURCES.md#r17)

subsecondの結果は倍率だけでなくms差を示す。通常ブラウザ設定を製品判定の対象とし、特殊JIT flags、debug runtime、profiler接続時の速度は診断値とする。

### 9.4 試行数・ノイズ・予算

診断はfocusedケースから始める。採用候補は少なくとも5回、通常7回のtimed samplesを採り、base/candidateを交互または固定seedで順序無作為化する。warmupはtimed samplesと分離して全件記録する。測定中にビルドや別ベンチマークを並列実行しない。温度・電源・CPU affinity・P/E-core・OS/WSL・ブラウザtab状態を揃える。

中央値とmin/max、必要ならIQRを示す。5〜7回でp95を強い結論に使わない。ノイズが大きければ追加の独立blockを一回測る。それでも判断不能ならINCONCLUSIVEとして候補を既定にしない。重い全suiteを候補ごとに無制限に再実行しない。S00でfixtureごとのtimeout・最大反復回数・同時実行数・メモリ上限を、baseを観察して記録する。

### 9.5 最小の証拠構造

既存runner/collectorを確認し、必要ならその一つだけを小さく拡張する。汎用ベンチマークframework、dashboard、複数の似たcollectorを新設しない。既存snapshotのrendererへ計測処理を混ぜない。[R05](SOURCES.md#r05)

推奨する作業領域は`WORK_DIR`配下の`reports/`と`evidence/`である。planルートの`STATUS.md`だけを進捗の正本とし、各reportには生証拠の実在パスを記す。再利用可能な測定commandはS00 report内に一度定義し、他reportから識別子で参照する。

fixture記録には、case ID、実在入力path/hash、program/task、全argv、サポート契約、expected raw hash/所在、入力特徴、適用経路を含める。timing記録には、source SHAまたはSHA+patch hash、artifact hash、target/features/flags、runner/runtime/browser版、thread指定/実効値、cold/warm、測定範囲、反復番号、wall/CPU/memory、出力hash、終了状態を含める。元データが得られない項目は数値を補わず`not available`と理由を残す。

## 10. ブラウザ統合の扱い

S10では実際のgbdraw側SHAとpackaged Wasm artifact hashを固定し、変更したartifactが呼ばれていることを確認する。コマンドWasmと直接API、pair単位parallelと検索内部parallelを混同しない。既存hostコードではmodule共有とthread worker準備が存在するため、一から同機能を作らない。[R13](SOURCES.md#r13)

共有memoryのthreaded経路ではsecure context、`crossOriginIsolated`、worker内の能力、実際のspawnを確認する。配信headerが意図どおりでない場合に、ブラウザ保護機構を無効化して合格としない。[R18](SOURCES.md#r18) [R19](SOURCES.md#r19)

まず既存runnerで入力転送、初期化、計算、結果転送を測る。host費用が支配的なら、観測した一つの費用だけを改善する。module cache、worker再利用、instance再利用、Rust pool再利用は寿命が異なる。`use_current_thread()`のregistryに関する注意を、固定Cargo.lockの実装でも確認する。同じ生存threadにlocal poolを繰り返し作ることを安全な再利用と見なさない。[R15](SOURCES.md#r15)

cross-runの準備済み配列cacheや新handle APIは自動導入しない。効果と寿命・無効化条件が必要な場合だけ別決定にする。共有状態へ複数トップレベル検索を同時投入してよいとは仮定しない。まず既存APIの再入性・結果所有権を読む。

全体CPU予算とメモリ予算を一つのhost側所有者が管理し、外側pair数×内側thread数が無制限に増えないようにする。既存取消・失敗処理がある場合は、取消後のrun、threads変更、worker終了、状態漏出も確認する。未実装の取消APIをこの計画で追加しない。

主評価browserは実際の配布対象からS00で選ぶ。少なくとも利用対象のChromium系で速度を確認し、他にサポートを主張するengineでは互換性を検証する。実行できないbrowser/OSはNOT_RUNとし、Nodeだけの測定からブラウザ全般の成功を宣言しない。

## 11. セッション構成と依存関係

S00/S01で計測条件と対象経路を固定した後、S02〜S09は証拠で必要とされたものだけ実行する。番号順の総当たりはしない。通常はNから一つ、Pから一つの主改善を採用して再profileし、残余ボトルネックがある場合だけ次を選ぶ。S10を早めに行うとhostが主因だと分かった場合、カーネル改善よりhost側の小変更を先行できる。

| ID | セッション | 前提 | 性格 |
|---|---|---|---|
| S00 | base・契約・fixture・実行予算の固定 | リポジトリと参照資料へのアクセス | 必須 |
| S01 | 最小限のWasm計測と実験選択 | S00の実行可能部分 | 必須 |
| S02 | N greedyの往復コピー削減 | S01で到達・負担を確認 | 条件付き |
| S03 | N query再pack削減 | S01で前計算の候補を確認 | 条件付き |
| S04 | N一致延伸カーネル | S01でproduction hotspotを確認 | 条件付き |
| S05 | N既存並列仕事の効率化 | 実際に複数仕事のfixture | 条件付き |
| S06 | P単一スレッドの上位カーネル | S01のgapped/Kappa profile | 条件付き |
| S07 | P preliminary resetの削減 | fill/初期化の負担あり | 条件付き |
| S08 | P Kappa scratch再利用 | 対象match経路に到達 | 条件付き |
| S09 | P Kappa先行計算の有界化 | 安全に避けられるredoあり | 条件付き |
| S10 | browser/hostの測定と限定統合 | S00/S01、測定可能なgbdraw環境 | 必須の評価。変更は条件付き |
| S11 | 統合候補の独立検証と結論 | 採否済み候補、S10の状態 | 必須 |
| S12 | N新規分割の設計・隔離試験 | 単一仕事の課題が残り、所有者が起動 | 拡張。production変更は別承認 |

通常の流れは`S00 → S01 → 選択したS02〜S09 → S10 → S11`。S12を実施しても、その設計完了だけで新しい並列実装を認定しない。必要なら実装セッションを所有者承認のもと追加し、再度S11相当の検証を行う。

同一repositoryでの同時編集は既定で一つとする。特にS06〜S09は同じP engineに触れる可能性があるため直列化する。read-onlyレビューは別担当に分けてよいが、共通実行機でtiming測定を競合させない。

## 12. セッションの契約・進捗・停止

### 12.1 着手と完了

各sessionは該当する`INSTRUCTION_PROMPTS/Sxx_*.md`を実行指示とする。reportは[共通テンプレート](templates/SESSION_REPORT.md)を使用し、既存証拠形式がある場合は内容をそこへ対応付ける。sessionごとの新しいテンプレートを増やさない。

着手時に`REPO_ROOT`、`PLAN_DIR`、`WORK_DIR`、base、対象case、変更上限、最重要不変条件、採用基準を確認する。手元の未commit変更は所有者のものとして保護する。勝手なreset/clean、旧証拠削除、認定golden更新はしない。remote push、PR merge、release公開、配信設定変更は本計画だけでは許可しない。

原則、一つの根本仮説、一つの局所差分、base/candidate対照、採否、引き継ぎで完了する。測定のみでNO_GOが確定した場合も有効な完了である。どの変更も採用されない場合、その理由と未達目標を提出する。

### 12.2 状態の定義

`NOT_STARTED`は未着手。`READY`は前提確認済み。`IN_PROGRESS`は実施中。`ACCEPTED`は所定の局所ゲート合格。`REJECTED`は試したが基準不達で既定差分から除外。`SKIPPED`は適用条件なし。`BLOCKED`は必要資料・権限・環境が不足。`INCONCLUSIVE`は証拠が不足し採否保留。`COMPLETE`はS00/S01等の非実装成果の完了に使う。

局所ACCEPTEDと製品全体のリリース認定を同一視しない。S11のread-only独立レビューが未実施なら`REVIEW_PENDING`を最終reportに残す。利用可能なら`ncbi_parity_auditor`を使用し、なければ別の担当者/セッションによるread-onlyレビューを要求する。自分の再確認を独立監査と呼ばない。[R03](SOURCES.md#r03)

### 12.3 見込み違い・環境不足への対応

旧候補が既に改善済み、対象関数に到達しない、計測で微小、単一仕事、必要なNCBIソース/goldenがない、browserで実行不能の場合、理由を具体的に記録する。実行可能なread-only調査と資料整理は続けられるが、未検証を合格にしない。

同じセッションで次々に別仮説へ広げない。最初の小変更とその合理的な修正を試し、根拠が崩れたらREJECTED/SKIPPEDにする。改善が小さくノイズの範囲なら比較の見せ方を変えずINCONCLUSIVEにする。測定系を作ること自体を目的にしない。

## 13. 最終検証と提出物

最終提出は、変更・非変更範囲、採用/棄却一覧、固定baseと統合候補の実在hash、Native/serial Wasm/threaded Wasm/browserの結果、raw output一致、thread指定と実効仕事、latency/throughput/memoryを分けた表、再実行command、未実施範囲、戻し方を含める。

試験はfocused unit/boundary → 対象taskの既存互換性suite → 変更共通部の回帰 → 最終性能の順に広げる。unsafe・共有memory・状態再利用変更は、利用可能な安全性検査と敵対的境界テストを加える。Miri等で対応しないtarget/intrinsicを検査済みと偽らない。

最終S11では、初期base、各採用差分、統合candidateを対応付ける。差分を一つずつ取り外せる単位にし、rollbackは採用commitまたは保存patchの明示的なrevertで行う。無関係な利用者変更を巻き戻す`reset --hard`を手順にしない。

**完了とは、全部の候補を実装することではなく、必要な改善だけを、結果・安全性・性能の証拠とともに残すことである。** Native相当目標の達成状況は、検証済み条件ごとに正確に報告する。
