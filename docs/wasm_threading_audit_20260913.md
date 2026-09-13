# LOSAT Wasm threading・性能・再実行監査

監査日: 2026-09-13。**修正の実装は行わず、問題の特定と報告を成果物とする。** セッション前半に作成した修正候補は撤回し、対象ソースを着手前の保存内容へバイト単位で復元した。元から存在したユーザーの変更は保持している。本報告の測定対象は修正候補ではなく、復元した現行実装である。

## 結論

スレッド数の無視だけではない。次の問題を確認した。

1. スレッド対応 Wasm が入力サイズ・処理件数によって指定を実質的に無効化する。同じ入力で指定を 2 から 4 に増やすと並列実行が消えるケースも再現した。
2. Wasm のローカル Rayon プールが呼び出し元スレッドの登録状態を残し、次の並列検索でプールを作成できない。
3. BLASTN direct API の並列経路がメモリへの出力を選択せず、ファイル出力を試みて失敗する。
4. TBLASTX は複数 subject を入力しても subject 1 件分の統計を使い、E-value と採用ヒット数が NCBI と異なる。この問題は native にもあり、Wasm threading 起因ではない。
5. threaded library artifact は通常の直接API初期化後、最初の複数スレッド検索で停止する。command/libraryで起動契約が違ううえ、ビルド出力名も衝突する。
6. BLASTN の診断ログと性能集計が実際のワーカー起動を正しく表していない。

一方、**しきい値を取り除くだけで高速化するという結論は支持されない**。今回測定した小さい検索では、実際にワーカーを起動する方が遅かった。要求したスレッド数を尊重すること、起動コストを下げること、検索仕事を効率よく配分することは、別々の検証課題である。

## 監査範囲と判定方法

- BLASTN / megablast、BLASTP、TBLASTX の CLI → thread count → pool → 並列処理 → reduction → 出力を確認。
- native、serial command-Wasm、threaded command-Wasm、direct API の境界を区別。
- 現行ソース、NCBI 2.17.0 の C/C++、依存している Rayon / num_cpus の実ソースを参照。
- 実行ごとの終了コード、出力の SHA-256、実際の `worker start tid=`、診断の requested/effective/pool 数を照合。
- 生物学的出力は通常の NCBI の raw outfmt 6 と比較。行をソートして合格にすることはしていない。
- 同一インスタンスでの再実行、ワーカー起動条件、統計の query/subject 軸を分けて再現。
- 独立した `ncbi_parity_auditor` にソース経路と実測証拠を再確認してもらった。最終報告の事実確認と、40件の測定サンプルの正常終了・oracle一致・集計値の再計算も完了。

以下の「実測確認」は記載した fixture・artifact・ランタイムでの確認。「ソース確認」は分岐または依存実装を確認したもの。「要追加測定」は性能への影響量を確定していないものを指す。P1 は検索結果・実行可否・要求契約への直接影響、P2 は診断・制限・回帰防止、P3 は診断専用経路または性能改善候補の優先度であり、外部公開の severity 定義ではない。

NCBI の API は設定済み N スレッドを生成するが、NCBI CLI 自体には CPU 数への削減と local `-subject` 時の無視が**警告付き**で存在する。したがって今回の「指定数を維持、できなければ明示的エラー」はユーザーの LOSAT 実行契約として扱う。NCBI CLI が常に要求数を使うという主張はしない。生物学的出力の唯一の基準は引き続き NCBI である。

根拠: `ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:86–88,145–180`、`c++/src/algo/blast/blastinput/blast_args.cpp:3152–3187,3205–3236`。

## 確認した問題

| ID | 優先度 | 問題 | 確度・対象 |
|---|---|---|---|
| F01 | P1 | 複数 subject の TBLASTX 統計が 1 subject のまま | 実測・ソース確認。native / serial Wasm / threaded Wasm |
| F02 | P1 | BLASTN direct API の並列出力先が誤っている | 実測・ソース確認。threaded Wasm |
| F03 | P1 | 同一スレッドで次の Rayon pool を構築できない | 実測・依存ソース確認。3 エンジン |
| F04 | P1 | 入力しきい値でスレッド指定を無効化、指定数増加で直列化 | 実測・ソース確認。3 エンジン |
| F05 | P1 | BLASTP が host の CPU 推定値へ黙って制限 | 実測・ソース確認。特に独自 host |
| F06 | P2 | serial Wasm でも複数スレッド要求を成功扱い | 実測・ソース確認。3 エンジン |
| F07 | P1 | TBLASTX の subject 判定が既存 linking の並列化まで停止 | ソース確認。1 subject の場合も影響し得る |
| F08 | P2 | BLASTN の pool/worker 診断と性能記録が実態と不一致 | 実測・ソース確認。DP 経路 |
| F09 | P2 | Rayon の Wasm 上限 255 でさらに黙った縮小 | 依存ソース確認。大量 worker の実起動は未実施 |
| F10 | P3 | 診断用 PARALLEL_SCAN_CHUNKS が scan を並列化しない | ソース確認。既定で無効の診断経路 |
| F11 | P2 | direct API がエラー原因の chain を落とす | 実測・ソース確認 |
| F12 | P2 | スレッド要求・再実行を継続検証する gate が不足 | CI・テストソース確認 |
| F13 | P1 | threaded library の direct API 起動契約が成立していない | 停止を実測、artifact初期化経路も確認。停止命令は未特定 |
| F14 | P2 | command/library が同名 Wasm を上書きし得る | 実際の --bin / --lib build で確認 |

### 用語: Rust/Wasm と `libc` の関係

ここでいう `libc` はC標準ライブラリ、`CRT` は起動用の実行基盤を指す。Rustで書いたプログラムでも、対象OS・実行環境向けの標準ライブラリが内部でこれらを利用する。RustのWASIターゲットはWasm向けの`wasi-libc`と起動用objectを同梱し、今回のtoolchainにも`libc.a`、`crt1-command.o`、`crt1-reactor.o`がある。[Rust公式のWASI説明](https://doc.rust-lang.org/rustc/platform-support/wasm32-wasip1.html#requirements)、[threads targetの同梱runtime説明](https://doc.rust-lang.org/rustc/platform-support/wasm32-wasip1-threads.html#no-interop-with-c-required)。

F13が扱うのは「RustのWASI実行基盤で、スレッドを使う前の初期化が実行されるか」という境界。LOSATの検索処理をCで実装したという意味でも、NCBI BLASTのC/C++を実行部品として利用しているという意味でもない。`TLS`は各スレッド固有の状態を置く仕組みで、初期化にはその状態を参照するthread pointerの準備も含む。

### F01: TBLASTX の複数 subject 統計が現行 NCBI と異なる

`run_impl.rs:1602–1606` は現在処理中の subject の長さで有効長を計算する。`ncbi_cutoffs.rs:143–164` は `db_length = subject_len_nucl / 3`、`db_num_seqs = 1` を固定する。

現行 NCBI の通常 CLI は `InitializeSubject` → `CLocalDbAdapter(..., true)` → subject 全体の長さと件数 → length adjustment / effective search space という経路である。各 subject ごとの再計算は legacy の `db_length == 0` の場合だけ。

対応ソース:

- LOSAT: `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1602–1606`
- LOSAT: `LOSAT/src/algorithm/tblastx/ncbi_cutoffs.rs:143–164`
- NCBI: `c++/src/app/blast/tblastx_app.cpp:118–119`
- NCBI: `c++/src/app/blast/blast_app_util.cpp:204–210`
- NCBI: `c++/src/algo/blast/api/seqsrc_multiseq.cpp:175–180`
- NCBI: `c++/src/algo/blast/core/blast_setup.c:821–847,866–883`
- NCBI: `c++/src/algo/blast/core/blast_engine.c:1434–1443`

`small_test.fasta` の先頭 900 bp を異なる ID で複製した fixture の結果:

| query 件数 | subject 件数 | NCBI 行数 | LOSAT 行数 | native / serial Wasm / threaded Wasm |
|---:|---:|---:|---:|---|
| 1 | 3 | 99 | 111 | LOSAT 各 target の出力は互いに完全一致 |
| 3 | 1 | 111 | 111 | NCBI とも完全一致 |
| 3 | 3 | 297 | 333 | LOSAT 各 target の出力は互いに完全一致 |

3×3 の場合、E-value 列を除けば NCBI の全 297 行が LOSAT に存在し、LOSAT だけに 36 行ある。先頭の対応行の E-value は NCBI `7.59e-165`、LOSAT `2.59e-165`。対応行の identity、長さ、座標、bit score は同じ。全列を含む通常 NCBI との完全一致行は 0。4 本指定・しきい値 0 で実際に 3 workers を起動しても、LOSAT の出力は直列経路と同じだった。

診断目的だけで NCBI に `BL2SEQ_LEGACY=1` を設定すると、LOSAT と全 333 行が multiset 一致する。ただし subject 順序が違うため raw bytes は一致しない。これは legacy の subject 単位統計を使っていることの裏付けであり、**legacy を新たな期待値へ置き換えてよいという意味ではない**。

影響は過小な E-value と追加ヒット。後段の仕事も増え得るが、「実行時間が 36/297 だけ増える」といった速度換算はしていない。default genetic code なので承認済み `db_gencode` 例外は適用されない。新しい一般的な chaining 不具合として分類せず、まず複数 subject の統計所有箇所を対象とする。

### F02: BLASTN direct API は並列経路で Vec への出力を失う

`run_web_pair` は出力 Vec を用意し、それを返す。通常の最終出力は `in_memory` があれば `Writer(input.output)` を使う。しかし threaded Wasm の並列 subject 経路では `Path(&args.out)` を渡し、先に `return Ok(())` する。

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs:4585–4607`
- 同 `10590–10621`: `BlastnOutputTarget::Path(&args.out)` と早期 return
- 同 `10946–10949`: 到達すべき `Writer` の選択

F03で説明する、実行基盤を初期化したcommand moduleの診断用直接呼び出しで、n1 は成功して 549 bytes、並列条件を満たした最初の n2 は `Is a directory (os error 31)`、0 bytes となった。Web API が空の `PathBuf` を渡すことと対応する。閾値が小さな入力を直列へ戻すため、通常の小規模テストではこの障害が隠れる。

これは出力の受け渡しの不具合であり、速度低下として扱わない。F03 の次回プール作成失敗より先に起きる。

### F03: `use_current_thread()` の登録状態が残り、次の検索が失敗する

3 エンジンの Wasm pool builder は `use_current_thread()` を使用する。

- BLASTN: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:282–306`
- TBLASTX: `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:138–161`
- BLASTP: `LOSAT/src/algorithm/blastp/blast_engine.rs:1074–1099`

使用中の `rayon-core 1.13.0/src/lib.rs:534–546` は local pool の registry が leak すると明記。`src/registry.rs:295–311` は登録済みのスレッドで新しい current-thread pool を作ると `CurrentThreadAlreadyInPool` を返す。

以下はcommand moduleを`wasi.start(--help)`で起動し、その後direct exportsを呼ぶ**診断用の切り分け**。C系の起動処理を経由させることで、F13のlibrary初期化問題を迂回し、その先のengineの状態を確認した。productionでこの呼び方を推奨するものでも、reactor型direct APIの動作保証でもない。最初に初期化を省いたcommand probeは停止したため、この表の根拠から除外している。

その診断用command moduleの同一インスタンスでの結果:

| プログラム | n1 | 最初の n2 | 次の n4 | もう一度 n2 | 最後の n1 |
|---|---|---|---|---|---|
| BLASTP | 成功 | 成功 | pool 構築失敗 | pool 構築失敗 | 成功 |
| TBLASTX | 成功 | 成功 | pool 構築失敗 | pool 構築失敗 | 成功 |
| BLASTN | 成功 | F02 の出力エラー | pool 構築失敗 | pool 構築失敗 | 成功 |

各呼び出しの 100 ms 後、host の live worker 数は 0 だった。**動き続けるワーカーの leak は今回確認していない。確認したのは Rayon の呼び出し元登録状態の残留。** メモリサイズの増加だけで無制限メモリ leak とも判断しない。Wasm linear memory は通常縮まないためである。

毎回別プロセス・別インスタンスを作る command 実行ではこの問題は隠れる。通常の parallel iterator は `pool.install` 配下にあり、メインスレッドが仕事を拾わず常に遊ぶという別の不具合は確認できなかった。

### F04: 明示したスレッド数より独自の入力しきい値を優先する

対象:

- BLASTN `run.rs:718–755,4859–4898`
- TBLASTX `run_impl.rs:235–273,1190–1229`
- BLASTP `blast_engine.rs:1196–1227,1312–1339`

| プログラム | 主な Wasm 独自条件 |
|---|---|
| BLASTN | jobs > 1、jobs ≥ 2、subject 合計長 ≥ 262144 × N |
| TBLASTX | jobs > 1、jobs ≥ 2、subject 合計長 ≥ 262144 × N |
| BLASTP | jobs > 1、jobs ≥ 2 × N、subject 合計長 ≥ 4096 × N |

1 subject / 900 bp の megablast と TBLASTX、1 subject / 87 aa の BLASTP では n4 でも実際の起動 worker は 0。3 件の小さな subject でも size 条件により直列になる。通常の stderr には警告せず、`LOSAT_WASI_THREADS_DEBUG=1` で初めて理由が見える。

**BLASTN の DP は例外的に別の pool を作るため、「小さい BLASTN は必ず全処理が直列」とは言えない。F08 参照。**

指定数に比例するしきい値により、同じ入力で並列性が逆転する。LC738874 / LC738875 の 2 subject、合計 646708 bp、900 bp query を使う megablast では:

- n2: 下限 524288 bp、`parallel=true`、1 worker 起動。
- n4: 下限 1048576 bp、`parallel=false`、0 workers。

件数やサイズは起動費用を節約する判断材料にはなる。しかしその判断を明示的な N 指定より優先し、無警告で直列化するのは今回のユーザー要求に反する。

### F05: BLASTP だけが host の thread cap へ縮小する

`blastp_available_thread_cap` は `LOSAT_WASI_THREAD_CAP` を読み、なければ `num_cpus::get()` を使う。`blastp_effective_num_threads` は `requested.min(cpu_count)` を返す。

- `LOSAT/src/algorithm/blastp/blast_engine.rs:1005–1051`
- `LOSAT/tests/run_losat_wasi_threads.js:223–228`
- `num_cpus 1.17.0/src/lib.rs:435–454`: WASI は CPU discovery の実装対象外で、fallback は **1**。

実測では `LOSAT_WASI_THREAD_CAP=1`、要求 n4 が `effective_threads=1` に削減された。独自 host がこの LOSAT 固有環境変数を供給しない場合も 1 になる構造。BLASTN/TBLASTX の正の要求数には同じ cap 処理がなく、プログラム間で契約が違う。

runner の cap は worker spawn 関数全体の強制上限ではない。BLASTP の入力値として使われるため、名前から host 共通の制約だと受け取ると誤る。

native BLASTP にも `min(cpu_count)` があり、「native は全経路で常に指定数そのまま」という前提も正確ではない。NCBI は CPU 制限時に警告するが、この LOSAT 経路は通常の警告を出さない。

### F06: serial Wasm が複数スレッド指定を正常終了で無視する

BLASTN `run.rs:4638–4660` と TBLASTX `run_impl.rs:889–911` は serial Wasm の thread count を 1 にする。BLASTP は cfg で parallel 経路を除外し、serial の subject preparation は `_num_threads` を使用しない (`blast_engine.rs:718–737`)。

serial `wasm32-wasip1 --no-default-features` artifact に n4 を指定したところ、3 プログラムとも終了コード 0、n1 と同じ出力だった。n1 と出力が同じこと自体は正しいが、複数スレッドを使えないことを通知していない。

serial target で真の threads を実行できないという制約は維持すべきもの。将来の対応は未対応要求の明示的拒否を含めて設計する。なお `parallel` feature を除いた native build にも同種の cfg 経路があるが、今回その native variant の実行測定はしていない。

### F07: TBLASTX は subject 数を見て、別段階の linking まで直列化する

Wasm 判定の仕事数は `subject_records.max(subject_chunk_work_items)`。しかし同じ `use_parallel` が `apply_sum_stats_even_gap_linking_with_parallel` に渡される。linking には既存の複数 frame group を処理する `into_par_iter` がある。

- `run_impl.rs:243,2737–2745`
- `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs:633–695`

1 subject の scan をさらに分割できるかどうかと、複数 frame group の linking を並列にできるかどうかは別である。native は要求数で許可する一方、Wasm は subject 条件で linking の既存並列性も止める。これは利用可能な仕事の見積もりを別段階へ流用する問題。

並列 linking で何倍速くなるかは未測定。候補の pruning や HSP 構築を変更してよいという意味でもない。

### F08: BLASTN の診断が pool 0 と報告しても、実際には workers がいる

`speculative_traceback` は `requested_parallel && config.use_dp` などから決まり、pool は `use_parallel || speculative_traceback` で作る。診断は `decision.parallel` だけから pool 数と expected workers を作る。

- `run.rs:4756–4759,4923–4927`
- 同 `796–815`: 診断値の算出

900 bp、1 subject、`-task blastn -num_threads 4` の実測:

```text
rayon_pool_threads=0 expected_spawned_workers=0 parallel=false
実際の worker start tid= は 3 件
```

`LOSAT/tests/wasm_performance.py:854–868` は `parallel=false` を `effective_compute_threads=1` と記録するため、この誤表示を性能証拠へ持ち込む。`spawned_helpers` は別途数えるが、両者の矛盾を要求契約違反として失敗させない。

さらに、`effective_compute_threads` は現状では stage のログ値の最大を取ったもので、実際に同時に検索を処理した数の測定ではない。Node/V8 自体にも threads があり、process CPU / wall の比から Rayon の利用率を逆算するのも不適切。

### F09: Rayon の Wasm 最大 thread 数 255 を事前検証していない

`rayon-core 1.13.0/src/registry.rs:245–246` は要求数を `max_num_threads()` と比較して縮小する。`src/sleep/counters.rs:60–77` の 32-bit 設定では最大 255。

LOSAT の builder は要求値を直接渡し、構築後の `current_num_threads()` と要求数の一致も検査しない。従って engine のしきい値・BLASTP cap を通過しても、256 以上の要求を正確に反映できない。事前 gate が serial にした場合や先に thread 起動が失敗した場合には、この縮小まで到達しない。

大量 worker 起動による資源消費を避け、255 worker の実測は行っていない。これは依存実装に基づく確定した上限であり、OS でそれだけの worker を起動できることの保証ではない。

### F10: PARALLEL_SCAN_CHUNKS の scan は逐次処理

`LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS` を有効にすると subject の parallel traversal を通らなくなる。一方、scan interiors は通常の `into_iter`、chunks も逐次ループで処理する。

- `run_impl.rs:1784–1788,2290–2327,3165,3252`

名前と診断値だけで並列 scan の効果を測ったと判断してはいけない。linking は別途並列になり得るので実行全体が直列という断定もしない。既定で無効の診断専用経路であり、通常検索の減速原因として混ぜない。

### F11: direct API がエラー chain を表示しない

`LOSAT/src/web_api.rs:622,634,645,705` の `err.to_string()` は anyhow の外側の context だけを文字列化する。F03 は API 利用者には `failed to build ... thread pool` としか見えず、Rayon の登録済み状態なのか OS/host 起動失敗なのかを区別できない。

通常 CLI の anyhow 表示と同等の原因情報が API から得られないため、thread 指定・host 制約・再実行問題の診断を難しくする。エラー発生自体を隠して成功扱いしているわけではない。

### F12: 継続的な threading 検証の空白

- `.github/workflows/ci.yml:29–36` の Wasm job は serial `wasm32-wasip1` の build だけ。
- 同ファイルの `cargo test --all-features` は native target であり、`target_arch="wasm32"` の分岐は実行しない。
- 現在の release readiness workflow も serial Wasm を中心に検証する。これを threaded の証明として使えない。
- `test_wasi_runners.js` は小さな人工 Wasm を使って host の exit / trap / worker failure を検証する。実際の LOSAT の要求数・pool 選択・direct API の再実行を検証するものではない。
- `wasm_performance.py:798–873` は出力 hash / oracle / memory を gate にするが、要求数と実際の worker 数の一致を合否条件にしない。
- warm 測定は serial 限定 (`wasm_performance.py:924–928`)。`run_wasm_performance_warm.js:7–11` にも新規インスタンスを毎回作ると明記されており、resident/direct API の再実行証拠ではない。

既存の出力一致テストは必要だが、それだけでは F03/F04/F06/F08 を検出できない。

### F13: threaded library artifact の最初の並列呼び出しが停止する

復元した現行ソースを次の条件でbuildしたartifactを別名で保存して確認した。

```bash
cargo build --release --lib --target wasm32-wasip1-threads --features wasm-threads
```

このlibraryを標準のNode `wasi.initialize(instance)`で準備し、同じdirect APIを呼ぶと、3プログラムともn1は成功するが、続く最初のn2は`thread-spawn`を1回も呼ばず15秒でtimeoutとなった。F03の「次回のpool作成エラー」に到達する前の問題である。

| バイナリ上の確認 | command | library |
|---|---|---|
| `_start` export | あり | なし |
| `_initialize` export | なし | なし |
| Wasm start section | `__wasm_init_memory` | `__wasm_init_memory` |
| main thread初期化 | `_start`から`__wasi_init_tp`を呼ぶ | 対応呼び出し経路を確認できない |
| API呼び出しのwrapper | 通常のAPI export | `.command_export`、API後に`__wasm_call_dtors` |

独立監査でfunction bodyも確認した。command `_start`は`__wasi_init_tp` → constructors → main → destructors → proc_exitを実行する。library側の`__wasm_init_memory`は他の初期化関数を呼ばない。libraryでは`losat_web_run_pair`だけでなく`wasi_thread_start`のcommand export wrapperにもdestructor呼び出しがある。

したがって「libraryには初期化が一切ない」は誤りで、memory初期化はある。一方、通常のcommand/reactorが実施するmain threadの実行基盤初期化を、この呼び出し方では通らない。実toolchainの`crt1-command.o`と`crt1-reactor.o`はともに`__wasi_init_tp`を参照する。Node18の`wasi.initialize`は`_initialize`があれば呼ぶが、今回のlibraryにはその入口がない。

**初期化契約の不足は有力な原因。ただし実際に停止した命令・lockの特定は未完了であり、thread pointer初期化だけを直せば解消すると断定しない。** API終了のたびに呼ばれるdestructorと反復利用の契約も併せて調べる必要がある。

証拠: `direct-lib/*.stdout` / `*.stderr` / `*.command.json`、保存したlibrary hashは後述。F02/F03は初期化したcommand artifactで切り分けた追加障害であり、library版の最初の障害を解消した後も検証が必要になる。

### F14: command と library の出力名が衝突する

`LOSAT/Cargo.toml:1–2,40–41`の同名package/binと`cdylib`、`LOSAT/.cargo/config.toml:15–26`のbuild aliasにはartifact名の分離がない。この環境では`--bin LOSAT`と`--lib`がともに次のパスを生成した。

```text
target/wasm32-wasip1-threads/release/LOSAT.wasm
```

`--lib` buildの後には同じパスから`_start`が消え、direct APIのwrapperも変化した。bin/libを一度にbuildするとCargoがoutput filename collisionを警告する。threaded aliasは`--bin`も`--lib`も指定していない。

ファイル名だけを検査するpackagingや後続コマンドでは、初期化契約の違うartifactを取り違える危険がある。現在のcommand runnerは`_start`を確認してエラーにするため、この境界でのsilent fallbackは確認していない。実際のrelease配布物が誤っているという断定もしない。

今後は名前、export集合、thread imports、memoryのshared属性をartifact契約として結び付ける必要がある。今回の検証後は、作業ツリーのbuild出力パスに元のcommand artifactを戻している。

## 性能の測定結果と未確定事項

### 5 回測定: 指定尊重と高速化は別の課題

修正前の同じ threaded command artifact を使った。各設定 1 回 warmup、5 回測定、最後に 1 回 diagnostics。各反復内で n1 → n2 → n4 → n4/しきい値0 を交互に実行。wall は Node 起動から command 終了までの cold process 時間であり、compiled module の常駐実行時間ではない。計測期間に他の LOSAT build/benchmark は実行していない。

megablast は LC738874/LC738875 の 2 subject・計646708 bp、query は LC738874 の先頭900 bp。BLASTP は PajaWSV の最長 protein の先頭1600 aa を6件複製したsubject・計9600 aa、queryは同じproteinの先頭300 aa。全測定を通常 NCBI 2.17.0 の raw outfmt 6 と比較し、完全一致したものだけを集計した。

「しきい値0」は既存の診断用環境変数を0にして既存の並列経路を選んだもので、コードの改変や新しい分割アルゴリズムは使っていない。

| プログラム | 指定・条件 | 実際の起動 workers | wall 中央値 [最小–最大] 秒 | 最大 RSS MiB |
|---|---|---:|---|---:|
| megablast | n1 | 0 | 0.218 [0.217–0.221] | 117.9 |
| megablast | n2 | 1 | 0.281 [0.276–0.295] | 126.7 |
| megablast | n4 | 0 | 0.228 [0.218–0.233] | 117.7 |
| megablast | n4-thresholds0 | 3 | 0.407 [0.399–0.495] | 147.9 |
| BLASTP | n1 | 0 | 0.243 [0.228–0.275] | 128.8 |
| BLASTP | n2 | 1 | 0.304 [0.293–0.365] | 136.9 |
| BLASTP | n4 | 0 | 0.250 [0.221–0.254] | 125.3 |
| BLASTP | n4-thresholds0 | 3 | 0.422 [0.407–0.512] | 153.3 |

これらの小さな検索では、実際に n4 を使う方が n1 より遅い。しきい値による serial 化には起動費用の節約という背景があるが、ユーザーの明示指定を無視する契約問題は残る。ここから大きな検索の速度、native/Wasm 全般の倍率、thread count の最適値は結論できない。

TBLASTX の同じ2-subject測定は、初回 warmup の段階で F01 による NCBI 出力不一致を検出したため停止した。その時間を速度比較表へ混ぜていない。

### 起動内訳: 順番に立ち上がる Node workers

`run_losat_wasi_threads.js:232–278` は各 `new Worker` の起動完了を `Atomics.wait` で待ち、戻ってから次の Rust thread spawn を受ける。起動費用が順番に積み上がる構造。

外部の診断用 host コピーで時刻だけを追加した1回の n4 実行では:

| 区間 | 時間 |
|---|---:|
| shared-memory guard の検査・書換え | 86.4 ms |
| その後の WebAssembly.compile | 22.0 ms |
| worker 1 の生成→起動確認 | 62.0 ms |
| worker 2 の生成→起動確認 | 58.8 ms |
| worker 3 の生成→起動確認 | 61.1 ms |

この実行も出力は NCBI と完全一致。区間値は1回の診断であり、5回の performance samples と同一のサンプルではない。V8 JIT/コンパイル活動などもあるため、全CPU費用を engine の計算コストへ帰属させていない。

### shared-memory guard は除去可能な無駄と断定できない

`run_losat_wasi_threads.js:71–74` は各 command で guard と compile を行う。`wasi_shared_memory.js:128–138,200–216` は bulk `memory.fill` / `memory.copy` をチェック付き関数呼び出しへ置き換える。今回の module では125 fill sites、1202 copy sites。元の bulk 操作自体は保持されており、バイトごとのコピーに退化させているわけではない。

cold startup への実測寄与はある。一方、この guard は現在の host で shared memory が成長する境界を扱う既存対応であり、単純除去を提案しない。将来の検討対象は、同一 bytecode の再解析回避・同一 module の安全な再利用・ランタイム別に必要性を実証すること。guard内部の動的検査が検索全体を何%遅くするかは未測定。

### 追加測定が必要な性能候補

- **仕事の粒度と偏り**: pool N と同時に処理可能な仕事 N は異なる。1 subject/1 chunk の scan、少数の Kappa matches、順序依存の heap/pruning を要求数だけで分割できるわけではない。worker別の実仕事時間は今回未計測。
- **全結果の一時保持**: BLASTN は `run.rs:10535–10587`、TBLASTX は `run_impl.rs:3189–3227,3338–3355` で subject ごとの結果を集めてから順序復元・縮約する。大量HSPでのpeak memory候補。nativeにもバッチ保持があり、Wasm固有の欠陥とはまだ断定できない。メモリ不足の再現・Big Oの厳密比較は未実施。
- **scratch と一時bufferの増加**: thread数増加に伴うallocation/RSS増加を広いfixtureで調べる必要がある。上表のRSS差だけでは長い検索の上限は分からない。
- **resident実行の費用**: cold Node processと常駐hostを同じグラフで比較しない。F03など再実行の正しさを確立してから、pool/モジュール再利用の費用を測る。

## 将来の修正順序と必要な受け入れ条件

本節は提案であり、実装・承認済み計画ではない。

1. **F01/F02: 出力の正しさを先に確立。** TBLASTXは全subject統計の所有者を一本化し、BLASTNは出力先の選択を並列/直列で共有する。複数query・複数subject・outfmt・default/non-default gencodeを分けてoracle比較する。
2. **F03とdirect artifact起動: ライフサイクルの契約を決める。** ローカルpoolに`use_current_thread()`を使う設計を再検討。候補はN専用計算workerと待機caller、または寿命とサイズ変更を明示したresident pool。単純に呼び出し元登録を残す実装は再利用できない。通常のRayon `Drop`は停止要求であり同期Joinの保証ではない点も確認する。
3. **F04/F05/F06/F09: 要求数を使うか明示的に拒否。** serial非対応、Rayon上限、OS/host起動失敗を区別。各engineで隠れたcapや独自しきい値が契約を変えないようにする。nativeとの違いも説明可能にする。
4. **F07/F08/F10/F11: 段階別の仕事量と診断を一致させる。** 設定数・実pool数・起動数・稼働数を別の値として扱う。
5. **F12: 回帰gateへ接続。** native / serial Wasm / threaded Wasmのn1/n2/n4、単一/複数subject、同じインスタンスで1→2→4→2→1、serialにn2、上限超過、host起動失敗、raw出力とworker数を検証する。
6. **その後に費用を削減。** worker順次起動、module準備、allocation、仕事の偏りを個別に測定。NCBIと異なるseed/pruning/統計を導入して速さを得ることは認めない。

`use_current_thread()`を外す場合、現状の「N−1 spawned workers + callerがworker0」という診断・hostテストも見直す必要がある。すべてのparallel iteratorが正しいlocal pool配下にあることも確認対象。新しいschedulerを入れるだけではF01/F02/F08は解決しない。

## 再現用fixtureと代表コマンド

以下は repository root から実行する例。出力を一時ディレクトリへ置き、既存の比較出力を上書きしない。

```bash
cd LOSAT
cargo build --release --bin LOSAT
cargo build --release --bin LOSAT --target wasm32-wasip1 --no-default-features
cargo build --release --bin LOSAT --target wasm32-wasip1-threads --features wasm-threads
cd ..
```

```python
from pathlib import Path
out = Path('/tmp/losat-thread-audit-repro')
out.mkdir(exist_ok=True)
for kind, source, length in [
    ('nuc', 'small_test.fasta', 900),
    ('aa', 'PajaWSV.faa', 300),
]:
    raw = (Path('LOSAT/tests/fasta') / source).read_text()
    seq = ''.join(raw.split('>', 2)[1].splitlines()[1:])[:length]
    for suffix, count in [('single', 1), ('multi', 3)]:
        name = f'{kind}_{suffix}'
        (out / f'{name}.fasta').write_text(
            ''.join(f'>{name}_{i}\n{seq}\n' for i in range(count))
        )
```

PajaWSV の最初の protein は87 aaなので、上の単一protein fixtureも87 aa。配列長を300 aaと誤記しない。

```bash
# F04: poolが0になる例。出力は成功する。
LOSAT_WASI_THREADS_DEBUG=1 node LOSAT/tests/run_losat_wasi_threads.js \
  LOSAT/target/wasm32-wasip1-threads/release/LOSAT.wasm \
  blastn -task megablast \
  -query /tmp/losat-thread-audit-repro/nuc_single.fasta \
  -subject /tmp/losat-thread-audit-repro/nuc_single.fasta \
  -num_threads 4 -outfmt 6 -out /tmp/losat-thread-audit-repro/megablast.out

# F08: 同じcommandを -task blastn にすると、ログpool0でもworkers3を確認できる。
# F06: runnerをrun_losat_wasi.js、artifactをwasm32-wasip1へ替えてn4実行。
# F05: BLASTPにLOSAT_WASI_THREAD_CAP=1を与えてrequested4/effective1を確認。
```

F01は`tblastx`で`nuc_single`/`nuc_multi`を組み合わせ、native・両Wasm・NCBI2.17.0に同じ`-query/-subject/-outfmt 6/-num_threads 1`を渡す。通常NCBIの比較では`BL2SEQ_LEGACY`を設定しない。

直接APIのextra argsは空白区切りではなく**NUL区切り** (`web_api.rs:143–145`)。再実行probeでは同じinstanceの`losat_web_alloc`で文字列を置き、`losat_web_run_pair`を1→2→4→2→1で呼び、result/errorのptr/lenを読む。小さいfixtureで並列経路に確実に入るため、関連する既存`LOSAT_*_WASI_MIN_*`を0にした。これらの診断条件を通常利用時のdefaultと混同しない。

## 環境・証拠・制約

- HEAD: `78c7556da6fa18b44447c3a3fa73774a8ee5956f`。作業ツリーには着手前から変更があり、clean HEAD certificationという主張ではない。
- Rust: `1.92.0 (ded5c06cf 2025-12-08)`、LLVM21.1.3。
- Node: `v18.19.1`。Linux x86_64 / WSL、Intel Core i9-14900HX、OSから32 logical CPUs。
- NCBI oracle: `/home/kawato/micromamba/pkgs/blast-2.17.0-h66d330f_0/bin/{blastn,blastp,tblastx}`、2.17.0+。
- NCBI source root: `/mnt/c/Users/genom/GitHub/ncbi-blast/`。
- シリアル/並列Wasmのflag: `LOSAT/.cargo/config.toml` の既存`+simd128`。最適化条件や生物学的処理を今回変更していない。
- 広いrelease fixture全体、全genetic code、全outfmt、ブラウザ各社、Wasmtime、WASI SDK各版、長時間のresident運用をcertifyしたものではない。

実行証拠の一時保存先: `/tmp/losat-wasm-thread-request/`。一時領域なので、長期保存には別途アーカイブが必要。本Markdownに主要な条件・結果・根拠を記載している。

| パス | 内容 |
|---|---|
| `audit-metadata.json` | commit、toolchain、source/artifact/runner hashes |
| `audit/commands.json` | 41 commandの引数・終了コード・hash・worker数・診断 |
| `audit/*.out`, `audit/*.stderr` | 通常のthread request/threshold比較、および明示したlegacy診断のraw証拠 |
| `parity-axes/results.json` | TBLASTXのquery/subject軸18実行 |
| `parity-axes/*.out` | NCBI/native/serial/threadedのraw出力 |
| `direct/*.command.json`, `*.stdout`, `*.stderr` | 初期化済みcommand moduleでの診断用direct再実行probe |
| `direct-lib/*.command.json`, `*.stdout`, `*.stderr` | 標準wasi.initialize後のlibrary direct API停止probe |
| `measurement/samples.json`, `summary.json` | megablast 28実行と、TBLASTXの不合格warmup1実行 |
| `measurement-blastp/samples.json`, `summary.json` | BLASTP 28実行 |
| `profile-host.stderr` | 起動区間の追加診断、raw出力も保存 |
| `audit_commands.py`, `parity_axes.py`, `run_direct.py`, `direct_probe.js`, `measure.py`, `measure_blastp.py` | 一時診断スクリプト。製品コードへ追加していない |

主要SHA-256:

```text
native command:
384dadffd6d38c64c6c40492a7d795d4f2b4214481f6a928931035819665dc03
threaded command Wasm:
b49d437373ff3f6b28b8bb7d5a38f7832d60c06574ea85d35e98a7d956edb05d
threaded library Wasm (F13):
cacc08c67643df66a089437c28bb8c819580e25da19177acb0a734e547b7c7f0
serial command Wasm:
272805273d0b8108fdb077933786239e18948e61fda26c7d532cc3c69b7737f8
TBLASTX 3x3 通常NCBI output:
d4e8a675fc637a2117386cd858e28675404f45c2509496b7434fe8b80de43c1d
TBLASTX 3x3 LOSAT output (serial/parallel同一):
99e5c76e5d0097455f0ec2b3e9e41c7c264a31a0089c1344a1421c0905231b87
```

測定結果から修正後の高速化・完全parityを主張しない。未解決の問題を残したまま「threaded Wasmはnative同様に指定を反映する」と認定することもできない。
