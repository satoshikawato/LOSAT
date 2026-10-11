# Wasm のスレッド版の共有メモリの上限を外す（1 GiB → 65,535 ページ）：記録

- 対象：[総合計画書](../../losat_web_gui_plan.md) の DW-28（保守者の決定、2026-10-12：「これは意味わかんないね。取っ払おうぜそんな制限」）。この決定で TD-7、G12、§2.3、§5.5、§8 と [`docs/web/abi_v2.md`](../../web/abi_v2.md) §3・§6 を改訂した。4 GiB を超える 64 ビットの Wasm（Memory64）は保守者の依頼で調査した
- ブランチ：エンジン側 `feature/wasm-mem-max`（`feature/losat-web-gui` の `141615ec` から。変更は `a6b49721`・`95c1a7f7`・`fe2ca26b`・`dcf9d591`、検証した tree は `dcf9d591`）。アプリ側 `feature/wasm-mem-max-app`（`feature/losat-web-gui-app` の `81dbeec9` から。変更は `7172aaa0`。記録は同じ名前のフォルダの `app/README.md`）
- 作業：2026-10-12、2 セッション。環境は WSL2（x64、32 スレッド、31 GiB）、Rust 1.92.0、Node v26.8.2
- 状態：実装と検証は終わった。PR を 2 本開いて保守者の merge を待つ（アプリ側を先に。下の「残件」）。実機（iOS、Android など）での確認は残件

## 結論

1. 1 GiB は技術的な制限ではなく、Rust の `wasm32-wasip1-threads` ターゲットの既定の link 引数（`--max-memory=1073741824`）だった。`wasm-ld` は最後の `--max-memory` を使うので、rustflags に `-C link-arg=--max-memory=4294901760` を足して **65,535 ページ（4 GiB − 64 KiB）** を宣言する。
2. 変更の前後のモジュールの違いは、threaded の 3 つとも import 節の 3 バイト（共有メモリの最大値 16,384 → 65,535）だけで、code 節と data 節は同じ。serial のモジュールは 1 バイトも変わらない。出力は変更の前後と native ですべて同じで、速さも変わらない（A/B の 5 つの組で B/A の wall の中央値は 0.980〜1.031）。
3. 途中で、最大値とは別の不具合を見つけて直した：2 GiB 以上のポインタは JavaScript に負の数で届く。`node:wasi` はそれに `EINVAL` を返し（nodejs/node#62671）、アプリの ABI の binding と WASI の shim は例外を出した。これは最大値を持たない serial のモジュールにも前からあった。host はポインタと長さを `>>> 0` で符号なしにする（ABI v2 §6 に書いた）。
4. 宣言した最大値を確保できない端末がある。JavaScriptCore は `RangeError` を出し、guard region の無い V8（Android、アドレス空間の制限のある環境）は誤りを出さずに小さい予約に切り替え、最悪は初期ページだけになって伸びない。host は 65,535 → 49,152 → 32,768 → 24,576 → 16,384 ページの順に試し、各候補を「初期サイズより 1 ページ少なく作って 1 ページ伸ばす」で確かめる。16,384 ページ（これまでの最大値）も取れないときの動作は変更前と同じ（アプリは serial のモジュールで検索する。Node の host は失敗する）。
5. Memory64 は当面見送る。Chrome 133・Firefox 134・Node 24 は出しているが、Safari の安定版（27.0）は出しておらず、Rust には wasm64 でスレッドや WASI を使えるターゲットが無い（調査：[`memory64_study.md`](memory64_study.md)）。4 GiB を超える規模は、wasm32 のまま、データベースを linear memory の外に置いて分けて読む（計画 §2.3 のトリガーに書いた）。

## 変更

| コミット | 内容 |
|---|---|
| `a6b49721`（エンジン側） | `LOSAT/.cargo/config.toml` とその写し `web/adapter/.cargo/config.toml` の `wasm32-wasip1-threads` の rustflags に `-C link-arg=--max-memory=4294901760`。`LOSAT/tests/build_wasi_artifacts.py` と `web/adapter/tools/build_reactors.py`：rustflags の `--max-memory` とモジュールが宣言した最大値が違えば失敗する（引数を落とした build が黙って 1 GiB になるのを防ぐ）。Node の host：`node:wasi` を使う道具すべてで WASI の i32 の引数を u32 にして渡す（`wasi_artifact.js` の `unsignedWasiImports`。`wasi_thread_host.js`、`run_losat_wasi.js`、`check_wasi_api_limits.js`、`benchmark_wasi_reuse.js`、`run_wasm_performance_warm.js`）、export のポインタで memory を読む道具は `>>> 0`（`check_wasi_reactor.js`、`web/adapter/tests/v_abi.js`）。共有メモリの段階的な引き下げ。`LOSAT/tests/wasm_performance.py` は `LOSAT_WASM_MEMORY_MAXIMUM_PAGES=16384` の固定をやめた（書いた時点では宣言した最大値と同じ値だった）。試験：`test_wasi_runners.js`（2 GiB 超のポインタでの WASI の呼び出し、serial と threaded。引数の符号なし化）、`test_wasi_shared_memory.js`（引き下げ、2 GiB の境と 65,535 ページの最上部での fill と copy） |
| `95c1a7f7`（エンジン側） | アプリの試験の 1 行（`web/app/tests/unit/engine-runtime.test.ts`）：threaded の reactor の最大値を固定の 16,384 でなく、`build_reactors.py` が記録した値と比べる。アプリ側の `7172aaa0` と同じ変更（エンジン側の CI の app の job が新しい reactor で通るため） |
| `fe2ca26b`（エンジン側） | 計画書：DW-28 を足し、TD-7、G12、§2.3（YAGNI の行とトリガーの表）、§5.5、§8 を改訂。`docs/web/abi_v2.md` §3（最大値と host の引き下げ）と §6（ポインタと長さは符号なし） |
| `dcf9d591`（エンジン側） | `createSharedMemory` を `wasi_thread_host.js` から `wasi_artifact.js` に移した。`origin/main`（`perf/wasm-host` が `wasi_thread_host.js` の export の行を変えた）との merge が衝突しなくなる |
| `7172aaa0`（アプリ側） | `web/app/src/infra/reactor/shared-memory.ts`（新規。引き下げと伸長の確かめ）、`engine-worker/runtime.ts`（それを使う）、`reactor/abi.ts` と `reactor/wasi.ts`（ポインタと長さを `>>> 0`）、単体試験 `engine-memory.test.ts`、E2E `tests/e2e/engine-memory.spec.ts` と harness の worker。旧い reactor（16,384 ページ）でも新しい reactor でも動く |

## 検証

A は変更前（`141615ec`、最大 16,384 ページ）の build（`out-base`）、B は 2 つの `config.toml` の rustflags だけを変えた build（`out-new`。`141615ec` から `dcf9d591` までに Rust のソースの変更は無い）。native は同じ tree の release の実行ファイル。host の JavaScript と試験は `dcf9d591`。

| 項目 | 結果 | 記録 |
|---|---|---|
| A と B のモジュールの比較（節ごと） | threaded の command・reactor・Web の reactor は import 節の 3 バイト（最大値）だけが違う。command は name 節のモジュール名の文字列も違う。code・data 節は同じ。serial の Web の reactor は同じ bytes（`5f8faec5…`） | [`results/sections_base_vs_new.txt`](results/sections_base_vs_new.txt)、[`scripts/wasm_sections.py`](scripts/wasm_sections.py) |
| 最上部までの確保（`memprobe`：ブロックを確保できなくなるまで確保して印を確かめ、最上部で fill・copy、スレッド 4 本、時計、乱数、ファイルの読み書き） | 15/15 が通った（宣言 16,384・65,535・65,536 ページ × ブロック 256・64・8・4・1 MiB）。65,535 ページまで伸び、使えた最上位の番地は `0xfffeffb0`。65,536 ページと宣言しても使えた最上位の番地は同じ（最後のページは確保に使われない） | [`results/probe3_summary.txt`](results/probe3_summary.txt)、[`scripts/memprobe/`](scripts/memprobe/)、[`scripts/run_probe_matrix.sh`](scripts/run_probe_matrix.sh) |
| 引き下げ（アドレス空間を `ulimit -v` で絞った Node） | 1.8 GB の制限：B は 16,384 ページまで下がり、A と同じ 4,108 ページまで使えた（確かめの無い素の `WebAssembly.Memory` では、同じ制限で最大値 65,536 ページの memory は初期の 23 ページから伸びなかった。`wasm32_4g_study.md` の付録）。3 GB の制限：B は 65,535 ページのまま 16,398 ページまで使えた。制限なし：65,535 ページ | [`results/step_down_ulimit.txt`](results/step_down_ulimit.txt) |
| Node の host の単体試験（`test_wasi_runners.js`、`test_wasi_shared_memory.js`） | 24/24。変更前（`141615ec`）の host のファイルに新しい試験だけを置くと 4 件が失敗する：2 GiB 超のポインタでの WASI の呼び出し（serial と threaded）は `28 !== 0`（`EINVAL`）、残りの 2 件は新しい関数が無い | [`results/verify_summary.txt`](results/verify_summary.txt)、[`results/host_test_old_hosts.txt`](results/host_test_old_hosts.txt) |
| Python の単体試験（`test_wasm_performance.py`、`test_simple_comparison.py`） | 22 件と 22 件、失敗 0 | 同上 |
| TD-6 の同一性の検査（`check_build_identity.py`） | `identical: true`、失敗 0（2 つの `config.toml` の rustflags が同じ） | 同上 |
| V-ABI quick | 60 実行が native の CLI と一致、凍結ハッシュ 16/16（serial 16 行、threads 46 行） | 同上 |
| V-ABI quick、ヒープを 2.4 GiB 超に置く（`LOSAT_WASM_MEMORY_INITIAL_PAGES=40000`） | 60 実行が native の CLI と一致、凍結ハッシュ 16/16。threads の 46 行は 40,000 ページ（2.44 GiB）から始まる共有メモリで、確保はすべてその上に置かれる。serial の 16 行は memory を自分で定義するので上の行と同じ条件 | 同上 |
| v1 の WASI の行列（`check_wasm_threading.py`） | 通過：command と oracle の記録 306 件、reactor の lifecycle のゲート通過、形式の失敗 0 | 同上 |
| 高いメモリでの出力の一致（p03、d04、blastp、tblastn、blastn の 1・8 スレッド。native、A、B、ヒープ 2.4 GiB 超の B） | 10/10 で 4 つの出力の SHA-256 が同じ | [`results/high_memory.tsv`](results/high_memory.tsv)、[`scripts/ab_time.py`](scripts/ab_time.py) |
| A/B の時間（threaded の command、暖機 1 回の後に A と B を交互に） | p03（TBLASTX）の 1・8 スレッド、blastp の 1・8 スレッド、d04（TBLASTX）の 8 スレッド（d04 は 3 回、ほかは 5 回）。回ごとの B/A の中央値は wall が 1.007・1.022・1.002・0.980・1.031、CPU が 1.008・1.011・1.008・0.984・0.989。全実行で終了コード 0、出力は組ごとに 1 つの SHA-256。ほかのセッションも動いていた（load 2.5〜15.5）ので、回ごとの比は 0.76〜1.32 に散った。code 節と data 節が同じなので、差は計測の揺れと読む | [`results/ab_timing_summary.txt`](results/ab_timing_summary.txt)、[`results/ab_timing.tsv`](results/ab_timing.tsv) |
| アプリ（`npm run check`、FakeEngine の E2E、新旧の reactor の E2E、3 ブラウザでの最大値と 2 GiB 超の検索） | `npm run check` は通過（単体 916 件）。E2E は新しい reactor で 327 件通過・失敗 0、旧い reactor で 324 件通過・失敗 0。Chromium・Firefox・WebKit とも宣言どおり 65,535 ページの共有メモリを確保し、2.06 GiB から始めた共有メモリでの 4 プログラムの検索（4 スレッド）が通常の memory と同じ bytes。FakeEngine の E2E は 296 件通過、1 件失敗：WebKit の `offline.spec.ts` の試験 6 が boot の画面のボタンを待って時間切れ（エンジンのコードを通らない試験）。同じ試験は新旧の reactor の E2E の WebKit で通り、単独の 5 回の繰り返しも 5/5 で通った。再現しないので原因は特定していない | アプリ側のブランチの `app/README.md` |
| wasm64 の試作 build（`wasm64-unknown-unknown`、serial、ソースは変えない） | build できない。Rust 1.97.0 に `rust-src` を足し、`RUSTC_BOOTSTRAP=1` と `-Z build-std=std,panic_abort` で試した（固定の 1.92.0 は触っていない）。依存の `getrandom` 0.2.17 が `wasm*-unknown-unknown` を既定では受け付けず（`js` の feature を求める `compile_error!`）、そこで止まった。std の build と LOSAT 自身のコードがこのターゲットで通るかは、この試作では分からない。調査の経路 A の見積もり（WASI の置き換えが要る）と合う | [`results/wasm64_summary.txt`](results/wasm64_summary.txt)、[`scripts/wasm64_probe.sh`](scripts/wasm64_probe.sh) |
| GitHub の CI（Web） | エンジン側：run 38155449735（`fe2ca26b`）と 38155831282（`dcf9d591`）は 4 つの job（engine-web-api、adapter、app、browser-engine）とも成功。アプリ側：run 38155451498（`7172aaa0`）は 3 つの job が成功、browser-engine が 1 件の失敗（323 件通過。WebKit の `results.spec.ts:512` が dot plot の canvas の `data-titles` を 1 回だけ読んで null）。同じ失敗は、このブランチを作る前の `feature/losat-web-gui-app` の run 38109093750・38109446588・38114111384・38118880176 にもある既存の試験の競合で、この変更は結果の画面を触っていない。手元の WebKit では新旧の reactor とも通った（下の「残件」） | — |

計画のゲート（skill `losat-gates`）との対応：変更は `web/adapter` の設定と host、`LOSAT` の link 引数と試験の道具で、Rust のソースは変えていない。standard の段（reactor の build、同一性、V-ABI quick、v1 の WASI の行列）を行った。V-ABI full と capture は、code 節と data 節が変更の前後で同じなので行っていない。

## 調査

- [`wasm32_4g_study.md`](wasm32_4g_study.md)：1 GiB の出どころ、65,535 と 65,536 の違い、エンジンごとの予約（V8、SpiderMonkey、JavaScriptCore）、2 GiB 以上での JavaScript の落とし穴、Node。印（V：一次情報を読んだか手元で走らせた、S：二次情報、I：推測）付き。
- [`memory64_study.md`](memory64_study.md)：Memory64 の標準の状態、エンジンの対応と上限、速さ、Rust の道具、移植で作り直すもの、3 つの経路と見積もり。結論は「当面見送り」で、経路 C（wasm32 のまま、データベースを linear memory の外に置く）を勧める。見直す条件は、Safari の安定版が Memory64 を出し、かつ Rust に wasm64 でスレッドを使える経路（`wasm64-unknown-emscripten`、または wasm-bindgen の threads の memory64 対応）ができること。

## Owner-delegated の決定

保守者の常設の指示（2026-09-29、2026-10-07 に再掲）により、推奨の案を採った。保守者の決定は DW-28（上限を外す）と Memory64 の調査の依頼の 2 つ。

| 番号 | 決定 | 理由 |
|---|---|---|
| 3 | 宣言する最大値は 65,536 でなく **65,535 ページ**（`--max-memory=4294901760`） | 65,536 ページではバイト数が 2^32 になり、wasi-libc の `sbrk(0)` と dlmalloc の区間の終わりが 0 に回る（rust-lang/rust#144199 がこの回避を挙げる。Emscripten も 1 ページ手前で止める）。`memprobe` で 65,536 ページ目は確保に使われなかった。SpiderMonkey は 65,536 ページ未満のときだけ 32 ビットの境界検査を保つ。失うのは 64 KiB |
| 4 | link 引数は `build.rs` でなく、2 つの `.cargo/config.toml` の rustflags に置く | TD-6 の同一性の検査がこの 2 つを比べている。成果物の JSON が `target_rustflags` と `memory.maximum` を記録する（改訂前の §2.3 が求めた記録） |
| 5 | ブランチは `feature/wasm-mem-max` と `feature/wasm-mem-max-app`、専用の worktree。セッションは merge しない | `.worktrees/web-gui` に SFe の未コミットの変更がある。計画の決定の行（DW-28、TD-7）を変えるので保守者の merge を待つ |
| 6 | host は宣言した最大値を確保できないときに段階的に下げ、各候補を 1 ページ伸ばして確かめる。環境変数 `LOSAT_WASM_MEMORY_MAXIMUM_PAGES` を与えたときは下げない | JavaScriptCore は `RangeError`、guard region の無い V8 は誤り無しで小さい予約になる（上の表の `ulimit` の行） |
| 7 | `node:wasi` を使う道具すべてで WASI の引数を符号なしにし、export のポインタは `>>> 0` | `node:wasi` は負の数に `EINVAL` を返す（nodejs/node#62671 は未解決、修正の PR #66636 は未 merge）。serial の runner の不具合も直る |
| 8 | `wasm_performance.py` は最大値を 16,384 に固定しない | 計測は出荷する上限で行う |
| 9 | 計画書に DW-28 の行を足し、TD-7・G12・§2.3・§5.5・§8 から引く | 保守者の決定の記録 |
| 10 | アプリの instance の作り直しの閾値は 512 MiB のまま（「最大値の半分」でなく固定の値にする） | linear memory は縮まないので、大きな検索の後のメモリをブラウザに返すための値である |
| 11 | 引き下げは 16,384 ページ（1 GiB）で止める。それも取れないときは変更前と同じ動作。`createSharedMemory` は `wasi_artifact.js` に置く | 1 GiB 未満の threaded を許すかどうかは製品の判断なので保守者に残す。serial のモジュールのメモリは前もって予約されず、小さい共有メモリの最大値より先まで伸びうる |

## 残件

- **実機で未確認**：iPhone・iPad（Safari と、WebKit を使う iOS の Chrome）、Android の Chrome（64 ビットと 32 ビット）と Firefox、macOS の Safari、Windows での commit の量。どの段で確保できるかは端末で違いうる。確認したのは Linux の Chromium・Firefox・WebKit（Playwright）と Node だけ。
- **1 GiB 未満の threaded を許すかどうか**（製品の判断）：今は 1 GiB を確保できない端末では、変更前と同じく serial のモジュールで検索する。
- **メモリ不足の知らせに大きさを出す案**：検索の途中で `memory.grow` が失敗したときに、確保できた最大値とその時点の大きさを知らせに入れる（今の知らせは変更前と同じ）。
- **取り込みの順**：アプリ側の PR を先に merge する。エンジン側だけを先に入れると、旧いアプリは引き下げをしないので、4 GiB を予約できない端末で threaded が始まらず serial になる。SFe の未コミットの変更は `docs/web/abi_v2.md` の 86 行目を触るので、SFe を先にコミットしてからこのブランチを取り込む。
- **この変更の外の、アプリの E2E の揺れ 2 件**（アプリの track に渡す）：WebKit の `results.spec.ts` の 784〜785 行目は dot plot の `data-titles` を待たずに 1 回だけ読む（CI で失敗した。このブランチの前からある）。WebKit の `offline.spec.ts` の試験 6 は FakeEngine の E2E で 1 回だけ時間切れになり、再現しない。
- **Memory64**：上の条件が揃ったら見直す。それまでは大きなデータベースを linear memory の外に置く設計にする（`dbsearch-20261012` の作業に知らせた）。

## このフォルダ

| ファイル | 内容 |
|---|---|
| `wasm32_4g_study.md`、`memory64_study.md` | 調査の報告（英語、出典付き） |
| `scripts/` | build（`build_variant.sh`、`build_probe.sh`）、`memprobe/`、`run_probe_matrix.sh`、`wasm_sections.py`、`verify_engine.sh`、`ab_time.py`、`run_ab_timing.sh`、`wasm64_probe.sh`。パスはこの機械のもの |
| `results/` | 上の表の要約（小さいファイルだけ。実行ごとの出力は `$BUILD_ROOT/wasm-mem/` と作業フォルダに置き、commit していない） |
| `evidence.sha256` | このフォルダのファイルの SHA-256 |
