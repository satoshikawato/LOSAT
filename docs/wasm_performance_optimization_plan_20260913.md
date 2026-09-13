# Wasm性能改善計画：TBLASTX → BLASTN → BLASTP

作成日: 2026-09-13

状態: 計測基盤と TBLASTX 限定 TurboFan 比較実行設定を実装。ユーザーのメモリ増許容を受けた [追加報告](evidence/wasm_performance_20260913/run-20260913-02/REPORT.md)に再測定・採用と追加 Rust 3候補の棄却を記録。[初回報告](evidence/wasm_performance_20260913/run-20260913-01/REPORT.md)は凍結したまま保持。P4は既存 megablast 出力差により PARTIAL。

対象: LOSAT Native、WASI serial command、WASI threaded command。再利用時の検証には既存reactorも使用する。

## 1. 目的と実行順序

NCBIとの出力一致を維持しながら、次の順に改善する。

1. **TBLASTXのコンパイル方式を検証する。** V8の初期コンパイルと最適化への移行がn8の低速化にどれだけ寄与するかを確定する。
2. **BLASTNの並列DP効率を調整する。** 既存の先行traceback計算の同期・割当・メモリ費用を減らす。megablastを非回帰対象に含める。
3. **BLASTPの上位ボトルネックを最適化する。** Wasmのプロファイルから対象を選び、単一スレッドの処理効率とn8の実行時間を改善する。

最初に、失敗時間と古い出力が混在しない測定条件を整える。実装候補は一度に一つとし、各段階の採否を決めてから次へ進む。実行順序は `P0 → P1 → P2 → P3 → P4` とする。

本書は今回依頼された作業の計画である。旧計画全体の実行、リリース、公開、別リポジトリの変更は含めない。

## 2. 出発点となる観測と未確定事項

2026-09-13の診断は、HEAD `c0e123d771f3a58c7529aed4d3f684a1a81f4cb5` に未コミット変更がある作業ツリーからrelease buildして実施した。HEADだけではソースを識別できないため、実行時には作業ツリーの内容も固定する。

使用環境はRust/Cargo 1.92.0、Node 24.21.0、V8 `13.6.233.17-node.53`、NCBI BLAST+ 2.17.0。診断時の成果物SHA-256は次のとおり。

| 成果物 | SHA-256 |
|---|---|
| Native command | `834a99ed11f6a4415df6f119d11736a16e50515d1f0f851fda9a3f687688a35f` |
| Serial WASI command | `f687af514acb81a62219237c25425412d24c8d5a527d11a988e41dc114c53940` |
| Threaded WASI command | `13b88d349ec3449330b11428bb15e9a688448c4e144e756b6f038b4b9412a437` |

### 2.1. TBLASTXの反復測定

MelaMJNV/PemoMJNVA、遺伝暗号1、outfmt 6。各条件の診断実行後に、診断を無効にして新規プロセスで3回測定した。実行順を巡回させ、起動・コンパイル・worker生成・検索・出力・終了待ちを含めた。

| 条件 | n1中央値 | n8中央値 |
|---|---:|---:|
| Native | 4.188秒 | 2.767秒 |
| Wasm threaded、通常設定 | 5.878秒 | 7.531秒 |
| Wasm threaded、TurboFan固定 | 5.174秒 | 4.678秒 |

Serial WASI n1は6.180秒。TurboFan固定はNodeに `--no-liftoff --no-wasm-tier-up` を指定したもので、Wasmファイルと検索引数は同じである。この表の全出力は、その診断で実行したNCBI oracleとraw bytesが一致した。

V8は実行中のWasm関数を最適化コードへ途中置換しないと説明している。長い初回呼び出しが複数workerで重なることは原因候補だが、責任関数と実際のコンパイル段階の対応は未測定である。[V8 compilation pipeline](https://v8.dev/docs/wasm-compilation-pipeline)

### 2.2. 並列仕事と測定上の問題

| 対象 | 確認済みの事実 | まだ断定できないこと |
|---|---|---|
| TBLASTX MelaMJNV/PemoMJNVA | poolは8 worker、subjectは1本、linkingは4組 | 低速化の全量をlinkingやJITだけで説明できるか |
| BLASTN LC738874/LC738870 | poolは8 worker、subject/chunkは各1、最大16候補の先行DPバッチを選択 | 同期、負荷の偏り、後で不要になるDPの寄与率 |
| megablast EDL933/Sakai | subjectは1本、検索chunkは2個 | chunkごとの実仕事量と検索本体の最大短縮余地 |
| BLASTP WSSV/PajaWSV | subject探索86仕事、Kappa query redo177仕事を並列選択 | 大きなAP027xxx入力のNative/Wasm差を支配する関数 |

上記のn8診断のうちBLASTN、megablast、TBLASTXでは、workerのspawn要求からreadyまでの待ち時間合計が約0.7～1.0秒だった。これはworker稼働率ではない。BLASTPの小入力ではTurboFan固定による速度改善は観測されなかった。

添付の全体比較はAvCLPV/PsCLPVの途中でtimeoutし、その後の図は既存ログを読んでいた。BLASTNのMjPMNV/MlPMNVでNative n1の0.139秒は `-outfmt` 重複エラーであり、検索時間ではない。更新日が異なるログや、その時点で生成されていない出力を新しい比較へ混ぜてはならない。

診断原本は `/tmp/losat-wasm-diagnosis-20260913-x2i1a9bw/` の `REPORT.md`、`metadata.json`、`runs.jsonl`、`summary.tsv`、`checks.json` にある。これは一時保管先である。P0で原本を利用可能な永続領域へ保存し、消失していれば現在のfixtureを再実行する。本書の観測値を新しいbaselineの代わりにしない。

## 3. 共通制約

- [AGENTS.md](../AGENTS.md) と [verify-ncbi-parity-and-speed](../.agents/skills/verify-ncbi-parity-and-speed/SKILL.md) を適用する。実装前に対応するNCBIソースを読み、変更箇所へファイル名・行番号・関連snippetをコメントとして付ける。
- NCBIは挙動の唯一の権威であり、実行ファイルは検証oracleに限る。LOSATのruntime/build/FFI/fallbackから呼ばない。
- 候補HSP、seed、二-hit状態、X-drop、scoring、統計、丸め、tie-break、pruning、出力順を維持する。CBS/SEGの省略、近似計算、新しい候補削減で速くしない。
- local-subject TBLASTXの非既定 `db_gencode` は既存の承認済み例外を適用する。それ以外の差を例外へ含めない。
- Nativeの凍結PR5出力に対するGate Aと、登録済みNCBI platform fingerprintに対するGate Bを維持する。platform-local NCBIの出力でgoldenを作り直さない。
- Plain `wasm32-wasip1` はserial。並列は `wasm32-wasip1-threads` と `wasm-threads` を使い、同一threaded artifactのn1/n2/n4/n8を比較する。
- 明示されたn8を黙ってn1へ縮小しない。現在のsearch-scoped pool、worker終了、エラー伝播を維持する。段階ごとの仕事の選択とpoolの実数を区別する。
- shared-memory guardを推測で削除しない。永続worker pool、検索を跨ぐscratch cache、独自の配列分割、新しい公開APIは今回の初期実装範囲に含めない。
- 未コミット変更を保存し、production、計測・テスト、文書、生成物の差分を分ける。既存の証拠ファイルを上書きしない。

## 4. 測定と採用の共通規約

### 4.1. 比較条件

| 軸 | 条件 |
|---|---|
| Native | release n1/n2/n4/n8 |
| Wasm serial | serial artifact n1 |
| Wasm threaded | 同一artifact n1/n2/n4/n8 |
| コンパイル方式 | 通常設定、TurboFan固定。変更するのはこの軸だけ |
| 実行寿命 | 新規プロセス、compiled module再利用＋新規command instance、同一reactor instanceでの反復 |
| baseline | P0で固定した元baselineと、直前までの採用済み状態の両方を保持 |

全fixtureで全条件の直積を実行しない。各段階の主対象で候補を絞り、選んだ設定を他のfixtureと必要なtargetへ広げる。P2/P3のA/B比較ではNode/V8設定を固定し、JIT設定変更の効果をRust変更の効果へ混ぜない。Rust共通処理の候補は通常設定でも比較する。

### 4.2. 保存する証拠

- run ID、開始・終了時刻、task、全argv、入力hash、ソースsnapshot/patch hash、artifact hash、Cargo.lock、features、build flags、Node/V8・NCBI版、runner/guardのhash。
- CPU/OS、実行可能CPU数・制限、ファイルシステム、関連環境変数、output sink。独立した性能計測を同時に走らせず、測定中はbuildもしない。
- 各回のexit status、stdout/stderr、生出力hash、wall、CPU user/system、取得可能ならpeak RSSとlinear memory。取得不能な値は理由付きで未測定とする。
- Node準備、artifact検査、guard処理、compile、instance作成、worker起動、検索呼び出し、出力、終了待ち。重複する区間を単純加算しない。
- 要求threads、pool実数、spawn/ready/exited、段階の仕事数。worker別の仕事時間を測る場合は専用診断を使い、CPU/wallから稼働率を推定しない。

compiled-module再利用は最適化済みコードの再利用まで保証するとは限らない。初回検索、warmup、以降の検索を分け、moduleとinstanceのどちらを共有したかを明記する。APIコピーやI/Oを含む時間をengine単独時間と呼ばない。RSSのprocess lifetime high-waterを検索ごとのpeakと誤記しない。

### 4.3. 反復数・timeout・判定

1. 比較対象の出力一致と正常終了を先に確認する。既知の差やtrapがある対象では、無関係な性能試験へ広げず原因を切り分ける。
2. 探索は短中入力を中心に各条件3回。採用候補はwarmup 1回＋計測5回以上を、baseline/candidateの順序を交互にして実施する。cold測定では毎回新規プロセスとし、warmupを同一instanceの再利用と混同しない。
3. 長時間fixtureはまず対象設定の正しさを1回確認する。timeoutは短中入力300秒、長入力は初回3,600秒を出発点とし、観測後はbaselineの所要時間に対して十分な共通上限を実行前に固定する。timeoutを上限値の成功時間へ変換しない。
4. 暫定採用基準は、主対象の同じ測定境界のwall中央値を10%以上短縮すること。全サンプル、範囲とばらつきを示す。5回で判断不能なら追加5回まで行い、それでも不明なら `INCONCLUSIVE` とする。
5. 小入力、Native、他taskの非回帰対象で、`max(5%, 50 ms)` を超えるwall悪化を認めない。メモリ悪化の基準は `max(10%, 16 MiB)` とし、予算超過や反復時の継続増加は不採用。これらはNCBIの仕様ではなく、本計画の工学的な採否基準であり、P0で候補測定前に確定する。
6. raw bytes不一致、実thread契約の違反、crash、候補に起因するtimeoutは不採用。小さい差や単発の最速値だけで採用しない。無効・未実行の行も別表に残し、成功例だけで全体を改善済みとしない。

再利用経路の改善は、その経路の効果として報告する。準備費用を計測外へ移しただけでcold性能が改善したとは扱わない。n8が常にn1より速いことやNative完全同速は、達成を保証する目標にしない。

## 5. P0 — 測定条件を固定する

### 作業

1. 作業ツリーと既存変更を確認し、入力・ソース・runnerを含むbaseline snapshotを保存する。別target directoryで候補を作り、baseline artifactを差し替えない。
2. 結果をrun IDごとの新規ディレクトリへ出す。古い `.out` と新しい `.log` が組み合わさらないようにし、成功status・対応するargv/hashが揃った結果だけを比較表へ入れる。
3. [run_comparison.sh](../LOSAT/tests/run_comparison.sh) と [comparison_data.py](../LOSAT/tests/comparison_data.py) の必要箇所だけを修正する。statusがない旧ログは新しい性能比較では拒否し、過去資料として残す。失敗・timeout・同名の古い出力を使う検証を追加する。
4. 計測は既存の [wasm_performance.py](../LOSAT/tests/wasm_performance.py)、[benchmark_wasm_threading.py](../LOSAT/tests/benchmark_wasm_threading.py)、[benchmark_wasi_reuse.js](../LOSAT/tests/benchmark_wasi_reuse.js) を再利用する。新しい汎用benchmark基盤を作らず、必要なfixture選択・Node引数・run manifestを既存の担当箇所へ追加する。
5. Node引数はargv配列で記録する。`NODE_BIN` にオプション込み文字列を入れない。Node起動時に指定するコンパイル設定と、LOSATの検索引数を分ける。
6. 次の旧失敗を現在のartifactで再確認する。BLASTN `NZ_CP006932` selfのthreaded n8、BLASTP `AP027131.faa/NZ_CP006932.faa` のthreaded n8のout-of-boundsを、単なる速度欠測と扱わない。再現すれば正しさの修正を性能差分から分離し、影響する段階の採用を止める。

標準buildの入口は次のとおり。実行時はbaseline/candidateを識別できるtarget directoryへ分離する。

```bash
cd LOSAT
cargo build --release --bin LOSAT
cargo build --release --bin LOSAT --target wasm32-wasip1 --no-default-features --target-dir target/serial-command
cargo build --release --bin LOSAT --target wasm32-wasip1-threads --features wasm-threads --target-dir target/threaded-command
```

### 完了条件

新規runだけから、入力・版・設定・status・出力一致を説明できる比較表を再生成できる。旧ログを意図的に置いても成功結果へ混入しない。実行fixture、採否基準、timeout、証拠の永続保存先が固定されている。

## 6. P1 — TBLASTXのコンパイル方式を検証する

### 主対象と対照

| 用途 | query / subject |
|---|---|
| 最初の反復比較 | `MelaMJNV.fasta / PemoMJNVA.fasta` |
| 中程度の再現確認 | `MjeNMV.fasta / MelaMJNV.fasta`、`AP027280.fasta` self |
| n8低速化の強い対照 | `PemoMJNVB.fasta / LvMJNV.fasta` |
| 既存の長配列回帰 | `AP027131.fasta / AP027133.fasta`、query/db gencode 4 |
| 他taskへの副作用 | BLASTP `WSSV.faa / PajaWSV.faa`、megablast `EDL933.fna / Sakai.fna` |

### 作業

1. 同じartifact・入力・threadsで通常設定と `--no-liftoff --no-wasm-tier-up` を比較する。対象Nodeの `--v8-options` を保存し、実際のflag受理とworkerへの適用を確認する。
2. 主対象でn1/n2/n4/n8とserial n1を測る。初回cold結果とcompiled-module再利用時の初回・反復結果を分ける。既存hostはworkerへcompiled moduleを渡しているので、それ自体を未実装の改善案としない。
3. 別の診断実行でコンパイルeventとworkerの関数サンプルを取得し、scan、ungapped extension、linkingのどこがLiftoff/TurboFanで実行されたかを追う。必要なら既存段階timerを使うが、hot loopへの頻繁な時計呼び出しは避ける。
4. profilerがコンパイル段階や実行時間を変える可能性を確認する。診断中の時間を採用用のrelease timingへ混ぜない。main isolateだけのサンプルでworkerの支配処理を決めない。
5. コンパイル方式の効果と、linkingが4仕事しかないこと、負荷の偏り、buffer容量・初期化の費用を分けて説明する。必要なら「同じ最適化方式でのn1対n8」を追加し、並列化自体の費用を測る。
6. 短中入力で絞った設定を長入力と非回帰対象へ適用する。Node起動前の設定なので、別プログラムにも影響する全体の既定値変更を先に行わない。

### 判断と変更範囲

- 効果が再現できるNode条件は、再現可能な実行設定として文書化する。通常設定と明示的な実験設定の双方の結果を残す。
- 通常設定での残差が長い初回関数に帰属すると確認できた場合に限り、既存計算関数の呼び出し単位を調整する局所候補を検討する。検索区間、diagonal状態、HSP統合の境界は変えない。コンパイラによる再inline後も効果を確認する。
- ブラウザでNode flagsが使えるとは扱わない。WASI hostの検証からブラウザでの改善を主張しない。既存module再利用の効果は確認するが、ブラウザアプリへの統合は後続の別作業とする。
- runtimeの互換性を変えるguard削除や永続pool導入へ広げない。

主な担当箇所は [wasi_thread_host.js](../LOSAT/tests/wasi_thread_host.js)、既存benchmark runner、必要と判明した場合の [TBLASTX run_impl.rs](../LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs) と [linking.rs](../LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs)。

### 完了条件

同じartifactでのコンパイル方式・threads・実行寿命の比較を保存し、原因の確認済み範囲と残差を説明できる。出力一致を満たした設定・局所候補だけを採用する。効果が特定条件に限られる場合も範囲を確定してP2へ進み、効果が出るまで無関係な変更を続けない。

## 7. P2 — BLASTNの並列DP効率を調整する

### 主対象と対照

- 主対象: `LC738874.fasta / LC738870.fasta` と `AP027202.fasta / LC738875.fasta`、いずれも `-task blastn`。
- self・高密度・メモリ境界: `NZ_CP006932.fasta` self、`LC738873.fasta / LC738871.fasta`。
- 小入力: `AP027152.fasta / AP027202.fasta`。
- megablast非回帰: `NZ_CP006932.fasta` self、`EDL933.fna / Sakai.fna`、`Sakai.fna / MG1655.fna`。
- 仕事量の対照: 既存の複数query/subject・長さが偏るfixture。単一ペアのlatencyと複数ペアのthroughputを別評価する。

### 計測する境界

現行の `SPECULATIVE_TRACEBACK_BATCH_SIZE` と、scratch slotへ `skip(slot).step_by(num_threads)` で仕事を割り当てる経路を出発点とする。

| 指標 | 確認したいこと |
|---|---|
| batch数・batchごとのjob数 | 小さな並列起動を何回繰り返すか |
| DP計算数・結果を消費した数・順序付き判定で不要になった数 | 先行計算の無駄がどれだけあるか |
| worker/slot別の仕事時間と最後の完了待ち | 静的な割当が負荷の偏りを生むか |
| scratch生成・reserve/realloc・copy bytes・最大容量 | 同期以外の費用とメモリ増加 |
| 順序付きcontains/materialize/add、sort/purgeの時間 | 並列化できない部分がどれだけ残るか |

診断counterはworkerローカルに集め、段階の終わりで集約する。全workerの累積CPU時間とwallのcritical pathを区別する。通常実行へ常時の細粒度atomicやログを加えない。

### 候補を選ぶ順序

1. batch内jobが0/1のときの不要な並列iterator、結果用の一時vector、再確保の費用を測り、根拠がある箇所だけ除く。poolの要求数は維持する。
2. 既存scratchの容量・結果slotを再利用し、所有権と初期化範囲を明示する。検索終了後まで容量を無制限に保持しない。
3. 負荷の偏りが主因なら、同じjob集合をscratch付きの既存pool内で割り当て直す。結果は元のprelim indexに戻す。jobごとのscratch生成へ退行させない。
4. 同期回数が主因なら、固定batch候補を小さな範囲で比較する。例として8/16/32/64を上限とした探索から始め、batch増加による不要DPとpeak memoryも計測する。fixture名に依存した閾値や、結果を見てから候補集合を変える仕組みを作らない。

各候補を別々に評価し、不要な候補を重ねない。NCBIの順序付きcontainment判定、start offset調整、traceback、identity判定、interval-tree挿入は元の順序を守る。先行計算はscratch内の純粋な計算に限り、共有treeやhitlistを更新しない。採否・統計・最終HSPを先行確定しない。

### megablastの扱い

DP経路を通らないmegablastで、DP調整の効果を主張しない。2 chunkや1 chunkしかない入力では仕事量不足を記録する。NCBIの `MAX_DBSEQ_LEN`、overlap、mask境界を性能だけの理由で変更しない。さらに細かいscan分割が必要と判明した場合は、その状態継承の設計課題を残し、本段階へ混ぜない。

主な担当箇所は [BLASTN run.rs](../LOSAT/src/algorithm/blastn/blast_engine/run.rs) と [alignment/gapped.rs](../LOSAT/src/algorithm/blastn/alignment/gapped.rs)。greedy、purge、interval treeそのものは、プロファイルとNCBI参照で必要性が判明しない限り変更しない。

### 完了条件

主対象の採用用A/B測定と非回帰条件を満たし、元のprelim index順で同じHSPを採用することを説明できる。n1/n2/n4/n8、複数query/subject、負strand、同点・包含HSP、短長alignment、反復・エラーを検証する。trapが残る経路を成功扱いしない。

## 8. P3 — BLASTPの上位ボトルネックを最適化する

### 主対象と対照

- 大入力主対象: `AP027078.faa / AP027131.faa`、`AP027132.faa / NZ_CP006932.faa`。
- 正常終了の確認を兼ねる対象: `AP027131.faa / NZ_CP006932.faa`。
- 小入力: `WSSV.faa / PajaWSV.faa`、`SicyWSV.faa / CoBV.faa`。
- 経路対照: 既存のsingle-query/many-subjectとmany-query fixture、標準行列とcomposition-adjusted matrix、restricted→exact retryを通るfixture。

### 作業

1. 現行Native n1、serial Wasm n1、threaded Wasm n1/n8を再測定する。通常コンパイル設定をbaselineとし、P1の設定は別条件として扱う。
2. Wasm自身のプロファイルを得る。現行の `blastp_timing_env_enabled()` はwasm32でfalseになるため、`LOSAT_TIMING=1` を設定しただけで段階時間が取れたと扱わない。
3. V8 samplingでworkerの上位処理を調べ、必要な場合だけWASIで使用可能な時計による段階計測を加える。`wasm32-unknown-unknown` へ同じ時計が使えるとは仮定しない。NativeのプロファイルだけでWasmの主因を決めない。
4. query/subject準備、scan/ungapped、preliminary DP、Kappa redo、traceback、merge/formatを比較し、呼び出し数・処理cell数・copy bytes・allocation数などで説明する。
5. 主対象のwallへ十分な寄与がある一処理を選ぶ。処理時間の構成と改善可能範囲から、共通の10%短縮基準へ届く見込みを確認してから実装する。

### プロファイルで選ぶ局所候補

| 支配的だった費用 | 最初に検討する変更 | 維持する条件 |
|---|---|---|
| DP/traceback内のメモリアクセス | 既存row/scratch配置、繰り返すcopy、残存dispatchの削減 | 同じcell、X-drop、gap tie-break、edit script |
| subjectごとの初期化 | diagonal/scratchの必要範囲と再利用を見直す | subject間の状態混入なし、NCBIと同じ初期状態 |
| Kappa redoの準備・確保 | 既存worker内のworkspace再利用、不要cloneの削減 | adjusted matrix、CBS/SEG、丸め、retry順 |
| 結果保持・merge | 中間vectorとcopyを削減し、元indexでreplay | pruning前倒しなし、完了順でのhitlist更新なし |

既に存在するBLOSUM62特殊化、SIMD設定、search-scoped pool、subject/Kappa並列化を作り直さない。`for_each_init` / `map_init` の初期化回数をworker数と同じと仮定せず、実数を計測する。大規模SIMD DP、別alignment算法、汎用arena/cacheは初期候補にしない。

主な担当箇所は [BLASTP blast_engine.rs](../LOSAT/src/algorithm/blastp/blast_engine.rs)、[gapalign.rs](../LOSAT/src/algorithm/blastp/gapalign.rs)、選ばれた処理の直接のcaller/helper。Kappa match-redoとquery-redoを混同せず、変更した分岐を通るfixtureを用意する。

### 完了条件

選んだ処理の費用削減と検索全体の時間短縮の両方を示す。NativeとWasmのn1/n8を比較し、小入力・他の行列・retry・出力形式の非回帰を満たす。候補HSP数を減らした結果を高速化として採用しない。効果不足なら候補を外し、未確定の原因を明記する。

## 9. P4 — 統合検証と引き渡し

1. 採用した変更だけを直列に統合し、元baseline対最終候補、各段階のbaseline対候補を区別して示す。全targetを最終ソースから再buildする。
2. 変更したNCBI移植関数には境界を含むunit testを追加し、対応NCBI unit testがあれば参照する。batch境界の15/16/17件など、分岐を実際に通す検証を選ぶ。
3. 各段階のparity sweep後、該当taskのcomparison suiteと共有runtimeに必要な回帰を行う。outfmt 6のほか既存0/7/custom契約を適用し、header除去・並べ替え・数値正規化をせず比較する。
4. TBLASTXのLC738874/LC738875、E-value 10/100/10000と、AP027131/AP027133のgencode 4回帰を保持する。NCBI local `-subject` の非既定db gencodeをそのまま失敗oracleにせず、既存のDB oracle等の検証契約を使う。
5. threading/hostを変更した場合は、既存の `check_wasm_threading.py`、`check_wasm_threading_regressions.py`、host testsを使用し、反復1→2→4→8→2→1、engine切替、終了、部分spawn失敗、trap、次回呼び出しへの影響を確認する。
6. Rustを変更した場合は、該当focused検証後に `cargo test --all-features`、`cargo clippy --all-features --all-targets -- -D warnings`、`cargo fmt --check` と両WASI targetのrelease buildを行う。host/比較器変更には該当JS/Python検証を適用する。
7. リリース向けの性能・parity主張、ordering/pruning/浮動小数点へ関わる変更の採用前には、`ncbi_parity_auditor` へ独立したread-only監査を依頼する。未測定・延期・失敗を合格へ変更しない。

[既存remediation結果](wasm_threading_remediation_results_20260913.md)の凍結回帰は99/145成功、46件未完了という記録である。本計画のfocused合格でその残りを完了扱いしない。広い性能主張に必要な範囲と実際の実行範囲を照合し、正式リリース認証は既存のGate A/Bと認証lineageに従う。

### 最終成果物

- 変更の目的、NCBI owner、実装差分、採用・棄却理由を記載した段階ごとのMarkdown報告。
- run manifest、全raw timing、出力hash/byte diff、診断profile、thread契約、memory結果。
- 同じ実行条件を比較した表。coldと再利用時を分け、失敗・未実行の一覧を付ける。
- 使用可能なNode設定、適用範囲、残る制約、再実行コマンド。
- production、test/計測、文書、生成物を分けたdiff確認と、英語commit title・短いsummaryの引き渡し。push/tag/publishは行わない。

保存先はP0で `docs/evidence/wasm_performance_20260913/<run-id>/` 等の永続領域へ固定する。巨大なraw出力は保管場所とchecksumを記録したarchiveへまとめ、計画書へ貼り込まない。失敗試行も保持する。

## 10. 参照実装と既存資料

NCBIソースのrootは `/mnt/c/Users/genom/GitHub/ncbi-blast/`。以下は計画作成時の対応箇所であり、実装時に現在の行番号とcall pathを再確認する。

| 対象 | NCBI参照 | LOSAT参照 |
|---|---|---|
| pool・thread寿命 | `c++/src/algo/blast/api/prelim_stage.cpp:145–188` | `LOSAT/src/utils/threading.rs` |
| subject chunk・WordFinder | `c++/src/algo/blast/core/blast_engine.c:246–268,478–498` | TBLASTX `run_impl.rs`、BLASTN `run.rs` |
| translated frame順 | `c++/src/algo/blast/core/blast_engine.c:804–845` | TBLASTX `run_impl.rs` |
| linking group | `c++/src/algo/blast/core/link_hsps.c:497–558` | TBLASTX `sum_stats_linking/linking.rs` |
| containment・DP・replay | `c++/src/algo/blast/core/blast_traceback.c:403–405,509–512,598–601` | BLASTN `run.rs`、`alignment/gapped.rs` |
| BLASTP subject探索 | `c++/src/algo/blast/core/blast_engine.c:1409–1475` | BLASTP `blast_engine.rs` |
| Kappa redo・thread-local状態 | `c++/src/algo/blast/core/blast_kappa.c:3429–3449,3493–3503` | BLASTP `blast_engine.rs` |

参考資料:

- [V8 compilation pipeline](https://v8.dev/docs/wasm-compilation-pipeline)：2026-09-13確認。Liftoff/TurboFan、実行中の置換、実験flags、計測時のコンパイル段階の説明。
- [V8 sampling profiler](https://v8.dev/docs/profile)：2026-09-13確認。使用runtimeでworkerを含む実際の取得範囲を確認して使う。
- [Wasm threading remediation結果](wasm_threading_remediation_results_20260913.md)：既存pool、reactor、guard、残る認証範囲。
- [旧N/P計画](losat_wasm_np_plan_20260906/MASTER_PLAN.md) と [STATUS](losat_wasm_np_plan_20260906/STATUS.md)：局所最適化の参考。旧thread閾値やsingle-jobでworkerを作らない指示は、現在の明示thread契約より優先しない。旧STATUSを本書の進捗で上書きしない。
- [BLASTN traceback/prune旧計画](blastn_traceback_prune_performance_plan.md)：過去のprofile・hit数は現在のbaselineに流用しない。

## 11. 進捗と判断の記録

| 段階 | 状態 | 成果物・判断 |
|---|---|---|
| P0 測定固定 | COMPLETE | 新規 run・status/argv/hash の検証、baseline 再build、旧 trap の非再現。完全な source snapshot は `snapshot-v2.json`。 |
| P1 TBLASTXコンパイル方式 | ACCEPTED（限定） | run-02で TBLASTX の cold Node/WASI 比較実行へ TurboFan を採用。主対象 n8 は6.56→4.44秒、RSS +27.95 MiB。旧2対照は約2.8倍、+47.33/+36.85 MiB。ユーザー承認後のメモリ条件で再評価。汎用採用・局所 Rust 候補の棄却記録は保持。 |
| P2 BLASTN並列DP | COMPLETE | 5入力の raw parity と DP 診断完了。Wasm batch32 は大入力 8.03%短縮・もう一方18.72%悪化で REJECTED。元の batch16 を維持。 |
| P3 BLASTP上位処理 | COMPLETE | 5入力の raw parity と Wasm profile 完了。初回DP行書込みに加え、run-02でquery間SEG cache再利用、DP関数境界、adjusted matrix行特殊化を個別実装・計測。最大5.45%短縮に留まり全て REJECTED。 |
| P4 統合検証 | PARTIAL | 検索Rust・Wasm成果物はbaselineと同一。run-02は43 process PASS、実比較12件raw一致。既存閾値12件と前段回帰を保持。既存megablast2入力のraw差と正式認証の未完了分は未合格のまま。 |

run-02 の方針変更はユーザー回答「TBLASTX限定でメモリ増を許容する」に基づく。
TBLASTX cold Node/WASI のピーク RSS 増分だけを `max(基準の10%, 48 MiB)`
へ変更し、出力・スレッド・wall の基準は維持する。これは計測上の採用基準であり、
実行時のメモリ上限ではない。比較スクリプトでのみ既定適用し、直接 Node 起動には
明示フラグが必要。旧対照の時間は元の CLOCK_MONOTONIC 観測値を再評価したもので、
新測定ではない。ブラウザや再利用時へ同じ倍率を外挿しない。

非実装の検証完了は `COMPLETE`、候補採用は `ACCEPTED`、棄却は `REJECTED`、対象外の根拠を確認した場合は `SKIPPED`、測定が不確定なら `INCONCLUSIVE` と記録する。前段階の非採用は次段階の独立した調査を妨げないが、未解決の共有runtime障害や出力差が影響する候補は採用しない。

各段階の報告では、共通の[証拠テンプレート](../.agents/skills/verify-ncbi-parity-and-speed/references/evidence.md)を使い、次に使うbaselineのソース・artifact・runtime設定と、残る課題を明記する。計画の保存を、実装や検証の完了として記録しない。
