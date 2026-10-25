# dev への持ち込み（線 NP 系列 WW0、2026-10-25）：記録

- 対象：#137（エンジン側、merge commit `7b8af58e`）と #136（アプリ側、`a821d123`）。2026-10-12 に `feature/losat-web-gui` / `feature/losat-web-gui-app` にだけ merge され、dev にも main にもなかった（`7b8af58e` は origin/dev の祖先でなく、`max-memory=4294901760` は origin/dev に無かった）。線 NP の DECISIONS 23 h で dev に持ち込む項目 WW0 にした（LOSAT Web の WindowMasker の counts、人で 691 MB、が 1 GiB に収まらないため）
- ブランチ：`feature/np-ww0`（origin/dev `fb8b10ee` から）。10 commit を `cherry-pick -x` で載せ（元の commit は各 commit の本文の `cherry picked from`）、dev の app に合わせた試験の直しを 1 つ足した
- 作業：2026-10-25、1 セッション（losat-baselines-04）。環境は WSL2（x64、32 スレッド、31 GiB）、Rust 1.92.0、Node v26.8.2、Playwright の Chromium・Firefox・WebKit。作業の記録は `/home/kawato/losat-baselines/ncbi-parity-20261019/WW/`（STATE.md、DECISIONS.md 1〜5）

## dev との衝突と直し

| file | 何が起きたか | 直し |
|---|---|---|
| `LOSAT/tests/wasi_thread_host.js` | dev が後から `blockingStdioWrite`（stdout と stderr の EAGAIN を待つ）を足していた | 符号なしにした WASI の import の `fd_write` を `blockingStdioWrite` で包む（両方を残す）。export は dev の形（`createSharedMemory` は `dcf9d591` で `wasi_artifact.js` に移る） |
| `docs/losat_web_gui_plan.md` | app branch の DW-26・DW-27 の行が dev に無い | DW-28 の行だけを足した |
| `web/app/tests/e2e/harness/engine.ts` | dev に無い `memoryDownloader` の import | その import を落とした（使っていない） |
| `web/app/tests/unit/engine-memory.test.ts`、`tests/e2e/engine-memory.spec.ts`、`tests/e2e/harness/high-memory-worker.ts`（文字の衝突は無し、typecheck が失敗） | dev の app には self-test の表（`src/examples/self-test.json`、`src/domain/examples.ts`）も `fetchReactor` も無い | worker は dev の `loadReactor` で threaded の module を読む。E2E は V-BR の 4 cell（`compact.blastn`、`width.blastp`、`width.tblastn`、`tblastx.multi`）で、2 GiB より上から始めた memory の出力を、通常の memory、同じ commit の native CLI、manifest にある NCBI の凍結 hash と比べる。unit は serial reactor の heap を 2 GiB より上に置いた instance と通常の instance を比べる（`b73127fa`、試験だけ） |

## 検証（head `636f31ae`、手元。ws には Node が無い）

手順は [`ww0_gate.sh`](ww0_gate.sh) と [`ww0_webkit.sh`](ww0_webkit.sh)（パスはこの機械のもの）。各段の 1 行は [`gate_summary.txt`](gate_summary.txt)。

| 項目 | 結果 |
|---|---|
| dev の reactor（`fb8b10ee`、同じ target directory で build） | threads の `memory.maximum` は 16,384 ページ、rustflags は `+simd128` だけ |
| WW0 の reactor | threads の `memory.maximum` は 65,535 ページ、rustflags に `link-arg=--max-memory=4294901760` |
| 節ごとの比較（dev と WW0） | serial は全節が同じ。threads は import 節の 3 バイト（最大値 16,384 → 65,535）だけが違い、code・data 節は同じ（[`sections_dev_vs_ww0.txt`](sections_dev_vs_ww0.txt)） |
| build の同一性（`check_build_identity.py`）、adapter の fmt と `cargo test` | 通過、14 件通過 |
| host の試験（`test_wasi_runners.js`、`test_wasi_shared_memory.js`、`test_wasm_performance.py`） | 40/40、OK |
| V-ABI quick | 68 実行が native の CLI と一致、凍結ハッシュ 16/16 |
| V-ABI quick、threaded の共有メモリを 40,000 ページ（2.44 GiB）から始める | 68 実行が native の CLI と一致、凍結ハッシュ 16/16 |
| v1 の WASI の行列（`build_wasi_artifacts.py`、`check_wasm_threading.py`） | command と oracle の記録 419 件、reactor の lifecycle のゲート通過、形式の失敗 0 |
| アプリ（`npm run check`） | 単体 406 件通過、40 件 skip |
| vitest（`engine-runtime.test.ts`、`engine-memory.test.ts`、reactor 付き） | 73/73 |
| E2E、エンジンの build | Chromium・Firefox 94 件通過、WebKit 43 件通過、失敗 0 |
| E2E、FakeEngine の build | Chromium・Firefox 64 件通過、WebKit 28 件通過、失敗 0 |
| 3 ブラウザでの共有メモリと 2 GiB 超の検索（`engine-memory.spec.ts`、4 スレッド） | 3 ブラウザとも 65,535 ページの共有メモリを確保。2,214,592,512 バイトから始めた memory で、最初の高い番地は 2,214,592,528。4 program の全 stream が通常の memory と同じで、outfmt 0・6・7 は native の CLI と同じ（[`engine-memory-chromium.json`](engine-memory-chromium.json)、[`-firefox`](engine-memory-firefox.json)、[`-webkit`](engine-memory-webkit.json)） |

- WebKit は、この機械に無い共有ライブラリ（`libxslt.so.1` ほか、Ubuntu 24.04 の package 16 個）を `apt-get download` で取り、`$BUILD_ROOT/np-ww0/webkit-deps/` に展開して、bundle の MiniBrowser を `LD_LIBRARY_PATH` 付きで起動する launcher（`LOSAT_WEB_WEBKIT_EXECUTABLE`）で動かした。system には何も入れていない。最初の gate の WebKit の失敗（engine 44 件、fake 31 件）はすべて起動の失敗（bundle の `MiniBrowser` の script が `LD_LIBRARY_PATH` を置き換える）で、launcher を直して WebKit だけを走らせ直した。
- V-ABI full、capture、エンジンの full gate は走らせていない。Rust のソースは変わらず、変わるのは threaded の link 引数だけで、code・data 節が dev と同じだから（#137 の記録と同じ判断）。native の build は対象の rustflags の外で変わらない。
- 時間の A/B は取っていない。code 節が同じで、#137 の A/B（5 つの組で B/A の wall の中央値 0.980〜1.031）と同じ条件である。

## 残件

元の記録（[`../README.md`](../README.md) と [`../app/README.md`](../app/README.md)）の残件はそのまま：iOS・Android・macOS の Safari・Windows の実機、1 GiB 未満の threaded を許すか（製品の判断）、メモリ不足の知らせに大きさを出す案。
