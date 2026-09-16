# Wasmに残る計算・メモリアクセス・起動費用の修正実装計画

作成日: 2026-09-15  
状態: **第3回のI4を本体へ反映済み。主要2入力で通常Nodeのn8逆転を解消し、全性能・非回帰条件、反映後の確認、独立確認を完了。**

実行結果: [採用判断と最終結果](evidence/wasm_remaining_work_20260915/run-01/ADOPTION.md)、[再現用の証拠](evidence/wasm_remaining_work_20260915/run-01/README.md)。
2026-09-16のユーザー指示により、`Wasm / Native <= 1.20`は採用条件・完了条件から外し、報告の中心にも置かない。現状から速くなることを重視し、出力一致と必要な非回帰確認で採否を判断する。以下の1.20倍への言及は当初案の記録であり、実行時の判定基準ではない。

### 第3回の実装結果（2026-09-16）

I4を`linking.rs`へ反映した。通常Nodeのn8検索本体は、Mjeで22.224→10.833秒（51.25%短縮）、AP selfで37.739→17.591秒（53.39%短縮）。両入力とも候補n8が候補n1より速い。Native/serial、module再利用、同一reactor反復、実ブラウザ反復は、元の時間・メモリ非回帰基準をすべて通過した。ブラウザのMje n1反復は2.0〜2.6%増であり、その値も結果へ残した。

[第3回の全結果・測定条件・残る制限](evidence/wasm_remaining_work_20260915/run-03/I4-RESULTS.md)。以下の候補経緯・当初計画は記録として保存する。

### 第2回の追加候補（2026-09-15）

ユーザーから本体の改善を継続する指示を受け、D0の生成コードと通常設定／TurboFan固定設定の差を起点に検証する。第1回の証拠・棄却判断は保存する。

- C1: 独立した先頭frame groupを完了してから残りを並列処理する。スレッドpoolとgroup内の状態は維持する。初回だけ速く反復時に悪化する可能性があるため、主要入力のmodule再利用を早期の採否条件にする。
- C2: NCBIの既存のchain選択反復を単位に、Wasmの関数入口へ戻る。反復本体、group状態、浮動小数点演算順を維持する。C1とは別の候補とし、生成コード・反復境界の全状態・Native/serial/threaded・再利用を検証する。
- 両案とも実検索だけを実行する。人工的な検索、warm状態のcache、暗黙のthread上限は導入しない。
- 保存先: `evidence/wasm_remaining_work_20260915/run-02/`。1.20のNative比を必須条件にせず、現状からの速度改善と非回帰で採否を決める。

### 第3回: n8逆転を優先する追加検証（2026-09-16）

- 環境再開時に、第2回の未反映ソースと生データを置いた`/tmp`が失われた。第2回の進捗記録を採用根拠には使わない。詳細は[復旧記録](evidence/wasm_remaining_work_20260915/run-02/RECOVERY.md)。
- I1は、既存large-gap predecessor走査を実際のWasm関数に切り出し、Nativeではインライン展開を指定する。過去C-X3との差は、現在のb0早期continueの維持、呼出し元の初期状態の全引継ぎ、Nativeインライン指定。検索候補・順序・数値演算・スケジューラは変更しない。
- 通常設定の主要2入力について、n8／n1の固定3組を先に測定する。有望な候補だけを状態全体・Native/serial・module/reactor再利用・実ブラウザへ進める。これは仮説の選別であり、先行する速度測定だけでは採用しない。
- ソース、実行ファイル、生データは[第3回の保存先](evidence/wasm_remaining_work_20260915/run-03/README.md)へ逐次保存する。ビルドキャッシュだけを`/tmp`に置く。
- I1は通常Nodeのn8でMjeが22.630→12.360秒、APが40.763→21.125秒となり、両入力で候補のn1より速くなった。ただしmodule再利用の最初の固定セッションで9.681→10.311秒（+6.51%）となり、非回帰基準を超えたため採用しない。
- I2bでは既存threaded-WASI判定の内側だけで借用sliceによる後方走査を行う。NCBIの候補・訪問順・演算順は維持し、生成コード上のアドレス再計算の削減を確認した。802ケース×3target×trace有無と、NCBI raw/thread gate 12条件が一致。最初の速度条件はI1が失敗したMje n8のmodule再利用とし、通過後に残りの採用検証へ進む。[診断と採否の状態](evidence/wasm_remaining_work_20260915/run-03/DIAGNOSIS.md)。
- I2bはMje n8のmodule再利用を通過したが、serial比較で数値上の悪化が出た。最終sampleには外部負荷も重なり、採否には使わない。I3では関数分割をthreaded-WASIだけに限定し、serialは明示的にインライン化する。生成コードと802ケースのserial状態一致を確認済みで、性能評価は未完。[I3の条件](evidence/wasm_remaining_work_20260915/run-03/I3-POLICY.md)。
- I3はserial固定比較で13.541→16.847秒（+24.4%）となり不採用。I4ではNative/serialの走査を元の呼出元へ戻し、共通macroで選択処理を一か所に保った。802ケース×3target×trace有無とserial raw 6条件は一致。実際の計測用Wasmも変更前の命令構造に戻ったが、定数値・memory offsetには差があり、速度の同一性は未判定。[I4の条件](evidence/wasm_remaining_work_20260915/run-03/I4-POLICY.md)。

## 1. 目的と起点

NCBI BLASTと同じ検索候補、状態更新、順序、浮動小数点演算、出力を維持し、Wasmの実行時間を短縮する。最初にTBLASTXの並列実行で遅くなる原因を調べ、その結果から変更を選ぶ。

起点は採用済みコミット`8f23f774`と、実装開始時の作業ツリー。計画作成時にも多数の既存変更があるため、コミット番号だけでbaselineを識別しない。既存変更を保存し、ビルド入力のhashを固定する。

前回の[採用判断](evidence/wasm_four_program_20260915/run-01/ADOPTION.md)を維持する。前回の測定を再開したり、六つの時間基準違反を未解決の採用判断として扱ったりしない。本計画は新しいbaselineから行う別の評価である。

### 既存の観測

以下は前回採用版の起動込み時間。Node 26.8.2、3組の中央値であり、現在の作業ツリーの再測定値ではない。

| 入力・処理 | Native n1 | Native n8 | threaded Wasm n1 | threaded Wasm n8 |
|---|---:|---:|---:|---:|
| MjeNMV/MelaMJNV・TBLASTX | 14.441秒 | 8.645秒 | 22.016秒 | 35.640秒 |
| AP027280 self・TBLASTX | 26.452秒 | 16.063秒 | 39.361秒 | 63.585秒 |
| AP027132/NZ_CP006932・BLASTP | 45.826秒 | 8.997秒 | 45.466秒 | 10.780秒 |

出典: [前回の全条件・Native比](evidence/wasm_four_program_20260915/run-01/integrated-three-pair-results.md)。

TBLASTXのn8逆転は確認済みの観測だが、原因は未確定。前回の検索本体・module再利用・TurboFan固定条件の最終計測は未実施である。通常ブラウザの速度へも一般化しない。

## 2. 対象と到達目標

- 主対象: TBLASTXの検索本体とthreaded Wasm。副対象: megablastの一致走査、Node/WASIホストの起動費用。
- BLASTNとBLASTPは共通ホスト・共有コードの非回帰対象。今回、根拠のない新しいDP最適化を追加しない。
- 主実行環境は実装開始時に固定するNode/WASIの通常設定。ブラウザは実際の利用側で別途確認する。Wasmtimeや`wasm32-unknown-unknown`の速度は、実測しない限り対象外と明記する。
- 「Native並み」の目標案は、事前指定した主入力について、同じ版・同じ実スレッド数・同じ検索本体区間の`Wasm / Native <= 1.20`。n1とn8を別々に判定する。達成予測ではない。
- NativeとWasmの絶対時間、起動込み時間、反復時の時間も併記する。Nativeの悪化による比率改善を成果に数えない。

## 3. 実施順序

| 段階 | 内容 | 実装へ進む条件 | 成果物 |
|---|---|---|---|
| B0 | 現行ソース・入力・実行環境と出力の固定 | 必須 | baseline manifest、raw比較結果 |
| D0 | TBLASTX n1/n2/n4/n8の仕事量・時間・生成コードの診断 | 必須 | 原因別の費用表、変更候補の選定理由 |
| X1 | large-gap走査で常に読む項目を小さな連続配列へ移す | D0で対象ループが主要費用 | 単独差分、状態比較、速度・メモリ比較 |
| X2 | 並列リンク用scratchの過大確保を減らす | D0で確保・初期化・memory growthの費用を確認 | 単独差分、使用量・容量・寿命の記録 |
| M1 | megablastの非圧縮一致走査をまとめて比較する | 一致長分布とprofileが有効性を支持 | 単独差分、境界比較、全体時間 |
| H1 | ホストの検査用moduleを再利用する | 検査・compile要求の重複が実時間に影響 | ホストのみの差分、ABI・起動時間確認 |
| V0 | 採用候補を統合し、共有経路と実利用側を確認 | 各候補の個別判定後 | 最終比較表、採用・見送り・未解決の一覧 |

D0で別の主因が見つかった場合は、その関数・命令・待ち区間を示す局所的な計画に更新する。X1/X2の実装を目的化しない。原因不明のままスケジューラ全体を書き換えない。

## 4. 実装全体の制約

1. NCBI C/C++が検索挙動の唯一の根拠。実装変更には実在するNCBIソースのpath・行番号・snippetを付ける。Rustの所有権、メモリ配置、Wasmホスト固有の変更は、その区別もコメントに残す。
2. NCBI実行ファイルは比較oracleに限る。runtime、build、fallbackに導入しない。
3. 感度、word size、候補pruning、X-drop、HSP linkingの実施回数を高速化のために変えない。`linked_set && !start_of_chain`は従来どおり出力時に除く。
4. `xsum = (h_xsum + score * lambda) - logK`の評価順を維持する。fast-math、精度低下、演算の再結合は使わない。
5. 指定n本は呼び出し元1本＋子n−1本。黙ってnを下げない。plain `wasm32-wasip1`はserialのまま。並列結果は既存のNCBI順へ戻す。
6. Workerの検索終了待ち、失敗伝播、キャンセル後の回復、APIの所有権とmemory budgetを維持する。永続Worker poolや検索間の可変状態cacheを混ぜない。
7. 採用済みのN1/N2/N3/M1/P1/P2/X2を新規成果として数えない。旧計画の同名IDとの混同を避け、本計画の証拠は専用ディレクトリへ保存する。

## 5. B0: 再現条件と正しさを固定する

### ビルドと記録

- baseline/candidateのsource、Cargo.lock、build.rs、Cargo設定、runner、fixtureをhash付きで保存する。必要なuntracked入力も含める。既存の出力ファイルを上書きしない。
- Native、serial command、threaded commandを別targetディレクトリへrelease buildする。reactorは反復/APIの検証対象となる候補で追加する。
- rustc、LLVM、Node/V8、NCBI BLAST+の版、CPU、affinity、メモリ制約、実際のコマンド、環境変数を記録する。
- Wasmとrunnerは両版とも同じ種類のファイルシステムに置く。コード変更がなければrunnerの実パスも同じにする。過去のDrvFS/ext4差を再混入させない。

### 最初のfixture

| 目的 | 入力 |
|---|---|
| TBLASTX主入力 | MjeNMV/MelaMJNV、AP027280 self |
| TBLASTX対照 | MelaMJNV/PemoMJNVA、LC738874/LC738875のevalue 10/100/10000 |
| megablast主入力 | EDL933/Sakai、NZ_CP006932 self、明示的な`-task megablast` |
| BLASTN対照 | LC738874/LC738870、AP027202/LC738875、word7の不均等複数subject回帰 |
| BLASTP対照 | AP027078/AP027131、AP027132/NZ_CP006932、既存の単一query・32 matches |
| 境界 | 有効queryのno-hit、短配列、負方向座標、複数query/subject、thread数未満の仕事数 |

`comparison_cases.tsv`の表示名とFASTA名は必ず実行manifestに解決する。特にBLASTNの表示名は実ファイル名と異なる。genetic codeと全argvも明示する。

fresh oracleとraw bytesを比較し、差があれば最初の不一致フィールドとNCBI ownerを特定する。既知の差がある対象経路では性能候補の試験を進めず、必要な正しさの修正を別差分で完了する。

PR5固定出力のGate Aと、登録済みNCBI platform fingerprintのGate Bは分離する。新しいplatform出力でexpectedを置き換えない。TBLASTX local-subject非既定`db_gencode`の例外はその挙動だけに適用し、AP027131/AP027133 code4回帰では対応するNCBI DB oracleを使う。

統計的に無効なTBLASTX queryのexit-status差、reactor反復memory、TBLASTX reuse時間の既存課題は別に明示する。本計画が直していない課題を合格扱いしない。通常の有効queryによる限定的な性能評価と、全入力の互換認証を区別する。

## 6. D0: TBLASTXのn8逆転を切り分ける

対象: `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`、`sum_stats_linking/linking.rs`、`LOSAT/src/utils/threading.rs`、既存のWASI計測ホスト。

### 6.1 仕事量

検索ごとに以下を集計する。初回linkingと再評価後linkingは別の記録にする。

- subject/chunk/frame group別のHSP数、group数、groupサイズ分布。
- large-gap helper訪問数、sumで拒否した数、座標で拒否した数、predecessor更新数。
- pass数、既存の変更なしskipの回数、chain取り出し数。
- scratch生成・拡張回数、初期化bytes、コピーbytes、capacityと実際の使用量。
- 実際の処理slot数、各slotの仕事量、開始・終了時刻、Worker起動・終了待ち。

並列のslotはWorker identityと区別する。n1/n8でchunk overlapなどの正当な中間処理量の差があり得るため、単純なcounter不一致を欠陥と断定しない。共通のmerge後境界でも候補・順序を照合し、追加処理が何に必要か追う。

### 6.2 時間と生成コード

- 起動込み、検査/変換/compile、instantiate、Worker起動、入力準備、検索各段階、出力、Worker終了待ちを分ける。
- 時計がWASI上で実際に有効であることを確認する。`LOSAT_TIMING=1`だけで計装が動いたと判断しない。
- countersは診断buildのslot内に保持し、段階終了時に回収する。内側のループで時計や共有atomicを更新しない。時間の採否は診断なしbuildで判定する。
- 並列groupの合計CPU/経過時間をwall時間の内訳として足さない。重なりとクリティカルパスを明示する。
- emitted Wasmと実際に動くV8機械語を確認し、余分なload、index計算、bounds check、compile tierを特定する。Rustの見た目だけで除去可能な命令数を主張しない。
- 通常Nodeを主条件とする。補助診断として両版に同じ`--no-liftoff --no-wasm-tier-up`を与え、検索時間とcompile費用を分ける。flagの対応は固定したNodeで確認する。

V8には、実行中のWasm関数を後から完成した最適化コードへ切り替えない仕様がある。ただしLOSATの逆転原因であることは未証明。[V8公式のコンパイル説明](https://v8.dev/docs/wasm-compilation-pipeline)。DevToolsによるtier変更も記録し、profile時の挙動を通常実行と同一視しない。

### 6.3 実装候補の決定

| D0の結果 | 次の行動 |
|---|---|
| 同じhelper走査が支配し、load/データ配置の費用が大きい | X1を実装 |
| 過大なscratch確保・初期化や拡張が支配する | X2を実装 |
| n8で余分な検索仕事が増える | 増えた処理のNCBI上の必要性を追跡。不要な重複だけを別候補として設計 |
| compile tierが主因 | 通常設定で遅い実関数と生成コードを示して局所案を作る。過去に不採用のループ抽出をそのまま再採用しない |
| 起動/終了待ちが短い検索を支配する | H1を先行評価。スレッド数を黙って縮小しない |
| どれも確定しない | 測定事実と未解決点を残し、TBLASTXの変更を見送る |

## 7. X1: large-gap用データを小さくする

NCBI owner: `c++/src/algo/blast/core/link_hsps.c:100–110,659–688,827–895`。LOSAT owner: `sum_stats_linking/linking.rs`の`LhHelper`、helper構築、INDEX1走査、`next_larger`更新。

1. 常に参照する`sum[1]`と`next_larger`を、一つの小さな連続配列へ移す案を最初の候補とする。座標・HSP index・small-gap用状態は同じhelper indexで参照する。
2. 移した項目を元の構造体にも残す二重管理はしない。構築、sentinel、再構築、両ループからの更新箇所を一括で修正する。
3. 既存のsum不成立時のjumpと、座標読み出しを遅らせる挙動を維持する。訪問順、`<=`/`>=`、同点時の選択、`changed`、`linked_to`を変えない。
4. 配列サイズと使用量を各targetで測る。単なるfield削除や型名変更を高速化とみなさない。全体の計算量・メモリ量のオーダーを維持する。

検証は、旧配置と新配置のhelper訪問列、各更新後のsum/next_larger、predecessor、chainとE-value bitsを比較する。両linking段階、small/large-gap、同点、cutoff境界、空/単一group、sentinel付近、frame切替、反復pass、n1/n2/n4/n8、trace有効時を含める。

前回の約40億訪問・37.9% sum拒否は前回baselineの診断値であり、X1の見積もりにはD0の新しい値を使う。構築・更新費用を含む全体時間が改善しなければ見送る。

## 8. X2: 並列scratchの確保量を実作業に合わせる

現行の並列linkingはthread-local poolの初回確保に検索全体の`total_hits`を使う。各Workerが実際に処理するgroupとの比をD0で調べる。大きなcapacityだけでは、物理メモリ消費や遅延の原因を証明したことにならない。

NCBI owner: `link_hsps.c:452–454,553–558`。NCBIのhelper領域の確保とframe group処理を参照し、並列化のためのRust側の所有者数は別に説明する。

条件が成立した場合のみ、各poolが実際に受け取るgroup長と必要なsentinelに応じて確保・拡張する。必要な初期化と残存状態を維持し、groupごとのshrinkやglobal cacheを追加しない。X1と同時に変更せず、それぞれの効果を測る。

検証: 小→大→小group、Worker間の異なる配分、再利用、エラー、n1/n8、reactorの連続検索。確保回数、使用量、capacity、post-call live allocation、linear memory、RSSを区別する。既存のreactor memory課題を「pool削減で解決」と推定しない。

## 9. M1: megablastの非圧縮一致走査をまとめる

対象: `LOSAT/src/algorithm/blastn/alignment/greedy.rs:2175`付近の`find_first_mismatch_greedy`。NCBI owner: `c++/src/algo/blast/core/greedy_align.c:313–375`の`s_FindFirstMismatch`。

### 最初の実装

1. D0とは別の診断で、`rem == 4`の呼出数、forward/reverse別の一致長分布と処理時間を取得する。圧縮経路への効果を合算しない。
2. 非圧縮経路だけを対象に、Wasmの16-byte SIMD比較を第一候補とする。各laneで「query < 4かつquery == subject」を判定し、最初に条件を満たさない位置を得る。
3. 両配列の有効範囲に収まるblockだけ読み、端数は既存scalar処理で扱う。逆方向ではlaneと走査順の対応を明示する。終了後のfence検出を維持する。
4. SIMD条件は既存の`target_arch`/`simd128`を使う。新しい利用者向けオプションは作らない。現行scalar経路を正しさの対照として残す。

検証: 長さ0/1/15/16/17/31/32/33、各laneの不一致、曖昧塩基、queryの非ACGT値、fence、非整列開始位置、forward/reverse、score-only/traceback、非アフィン/アフィン、scratch再利用。返す一致長と`fence_hit`に加え、最終score・座標・edit script・raw bytesを比較する。圧縮経路も非回帰確認する。

短い一致でsetup費用が勝つ場合は見送る。thresholdやblock幅を多数試し、最速だけを選ぶ探索はしない。別の整数比較案に進む場合は、新しい根拠と条件を先に記録する。

## 10. H1: 検査と実行のmodule所有者を一つにする

対象: `LOSAT/tests/wasi_artifact.js`、`run_losat_wasi.js`、`wasi_thread_host.js`。検索kernelと別に変更・評価する。

### H1a: serial

現行は`inspectArtifact`内で`new WebAssembly.Module(bytes)`を作り、呼出し側が同じbytesを再度`WebAssembly.compile`へ渡す。まず実費用を測る。二つのAPI呼出しがあることから、機械語生成も二重だと決めつけない。

検査済みmoduleを実行側へ渡す内部APIに整理し、compile済みmoduleを再利用する。ABI判定の所有者は一つとし、検査CLIのJSON形式、artifact identity、command/reactor判別、imports/exports/memory条件、例外の伝播を維持する。

### H1b: threadedは条件付き

threadedでは元bytesの検査後に`guardSharedMemory`で別bytesを生成し、そのmoduleを実行する。H1aと同じ方法で元moduleを実行することはできない。

まずraw検査・guard変換・guard後compileを別計測する。無駄が実測できた場合だけ、元bytesの妥当性確認とABI解析を維持しながら実行用moduleの生成を一回にまとめる局所案を作る。既存parserを共用できず大きな新規parserが必要になるなら、本ラウンドでは見送る。

共有メモリguardを性能のために外さない。現在のNodeでの必要性の再検証は別のruntime互換作業であり、この最適化に混ぜない。

NCBI参照はアプリの終了・エラー伝播境界（`c++/src/app/blast/blastn_app.cpp:172–176`等）、Worker開始/終了境界（`c++/src/algo/blast/api/prelim_stage.cpp:145–188`）。module検査と共有メモリ対策はWasmホスト固有であり、NCBIに同じAPIがあるとは説明しない。

検証: 正常serial/threaded command/reactor、取り違え、壊れたWasm、欠けたexports、不正memory limits、guard後のABI、実際のshared-memory拡張、元来不正な範囲のtrap、Worker失敗、反復、終了コード。既存の`test_wasi_runners.js`、`test_wasi_shared_memory.js`等を変更範囲に合わせて使う。

この変更の速度改善はNode/WASIホストについて報告する。実際のブラウザ利用側に同じ処理があることを確認しない限り、ブラウザ改善とは呼ばない。

## 11. 測定量と採用条件

### 測定を絞る

- 最初の詳細profileはTBLASTX主入力2件のn1/n2/n4/n8に限定する。候補の速度測定は対象主入力のn1/n8を中心とし、n2/n4は順序・thread契約の確認に使う。
- 各速度条件はwarmup 1回＋交互のA/B 3組を基本案とする。全sampleと中央値を残す。3組では不確定なら「不確定」とし、自動で測定を延長しない。
- 検索kernel変更では両版に同じrunner、flag、入出力方式を使う。ホスト変更では同じWasm artifactと対称な配置で比較する。
- 時間計測中にビルド、他ベンチマーク、ブラウザQA、圧縮処理を並行実行しない。単調時計を使用する。
- 最初から全候補×全runtime×全fixtureの総当たりをしない。個別で棄却した候補の統合計測は行わない。
- 反復は独立したAB/BAの2セッション、各warmup 1回＋測定3回を初期案とする。これで無制限反復のmemory安定性を主張しない。memory修正の主張には既存16回window等の対応する検証を別途行う。

### 時間区間

| 区間 | 含める費用 |
|---|---|
| 起動込み | process開始から検索・出力・Worker終了・process終了まで |
| 検索本体 | Native/Wasmで同じRust入口から出口。前処理、scratch確保、linking等を恣意的に外さない |
| module再利用 | compile済みmoduleから新instance/memoryを作り、検索・結果取得・Worker終了まで |
| reactor反復 | 同一instanceへの入力コピー、実行、結果取得、入力解放、Worker終了まで。setup時間も別途報告 |

### 候補の判定

1. 対象範囲のraw output、終了状態、必要な内部状態・浮動小数点bits・順序が一致すること。
2. 事前指定した主条件で、対象区間の中央値5%以上の短縮を目標とする。局所counterの改善だけで採用しない。
3. 対照条件の時間悪化は`max(5%, 50 ms)`、peak RSS増加は`max(10%, 16 MiB)`以内を基準案とする。API/linear-memory予算は別に維持する。
4. kernel改善と起動改善を別に判定する。compiler flag条件だけの改善を通常設定の改善に数えない。
5. raw差、trap、thread契約違反、候補由来timeoutは不採用。時間・メモリの未測定を合格にしない。既存失敗を新しい例外にしない。
6. 個別効果を加算しない。統合版をbaselineと比較して採否を決める。各候補は採用・見送り・不確定のいずれかを明記する。

長いcode4 oracleは正しさの確認として先に一度取得し、出力hash・完全なargv・binary hashを固定する。速度sampleごとにoracleを再実行しない。timeoutは実行前に固定し、過去に600秒超を要した事実を考慮する。超過はTIMEOUTとして残し、自動延長しない。

## 12. V0: 最終確認と引き渡し

1. 候補ごとのfocused parity sweepを完了し、その後に該当program/taskの比較suiteを実行する。最後に統合版のRust test、clippy、format、影響するruntime harnessを実行する。
2. X1/X2ではTBLASTXの短い閾値sweepと長いcode4を維持する。M1ではBLASTN共有greedy経路を含める。H1では四処理とcommand/reactorの境界を含める。
3. ブラウザ向けartifactを変更する候補では、実利用側へ一時overlayして通常設定で検証する。`browser-offline-qa`に沿い、外部通信遮断、isolated threaded/nonisolated serial、n8予算、反復、キャンセル・不正入力後の回復、pagehide後の終了を確認する。配布物の更新とは分ける。
4. ブラウザ速度は主入力を別測定する。計測できなければ未測定とし、Node結果で代用しない。
5. ordering/pruning/浮動小数点に関係する変更、parity解決やrelease向け速度の主張は、AGENTS.mdに従って`ncbi_parity_auditor`の独立read-only確認を受ける。
6. production Rust、ホスト、検証コード、文書、生成artifactを分けてdiff reviewする。変更前後のsource/build/artifact hash、raw output、全sample、除外・失敗理由、残る課題を保存する。

新しい証拠の保存先案: `docs/evidence/wasm_remaining_work_20260915/run-01/`。前回の証拠と採用判断は変更しない。

最終報告は「何秒短縮したか」「どの無駄が減ったか」「Native比」「確認したruntime/入力」「未達・未解決」を一つの表にする。目標未達でも結果を明記し、未根拠の大規模改変へ拡大しない。

## 13. 本提案での推奨着手点

**B0 → D0 → 原因に対応するX1/X2 → M1 → H1 → V0**。

TBLASTXのn8逆転原因を特定することが最初の完了条件。そこから、同じ検索を実行するために必要なデータ移動・準備・命令を減らす。Native同等の達成は、実測した範囲だけで判断する。
