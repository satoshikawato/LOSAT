# Wasm のスレッド版の共有メモリの上限を外す：アプリ側の記録

- 対象：計画 DW-28（保守者の決定、2026-10-12）と改訂した TD-7。エンジン側の記録は、エンジン側のブランチ `feature/wasm-mem-max` の `docs/evidence/wasm_memory_max_20261012/README.md`（結論、調査、Owner-delegated の決定 3〜11、残件）
- ブランチ：`feature/wasm-mem-max-app`（`feature/losat-web-gui-app` の `81dbeec9` から。変更は `7172aaa0`）
- 状態：実装と検証は終わった。PR は保守者の merge を待つ。**エンジン側の PR より先に merge する**

## 変更（`7172aaa0`）

- `src/infra/reactor/shared-memory.ts`（新規）：threaded のモジュールが宣言した最大値で共有メモリを作り、確保できないときは 49,152・32,768・24,576・16,384 ページの順に下げる。各候補は初期サイズより 1 ページ少なく作って 1 ページ伸ばし、両方が通った最初の候補を使う（作れても伸びない memory を使わない）。16,384 ページ（1 GiB）も取れないときは変更前と同じで、threaded は始まらず、検索は serial のモジュールで走る。`engine-worker/runtime.ts` がこれを使う。
- `src/infra/reactor/abi.ts`、`wasi.ts`：ポインタと長さを `>>> 0` で符号なしにする。2 GiB 以上のポインタは JavaScript に負の数で届き、変更前は `RangeError: Start offset -2080374736` になった（最大値を持たない serial のモジュールにも前からあった不具合）。
- 作り直しの閾値は 512 MiB のまま（決定 10）。`policy.ts`・`artifact.ts`・`block-store.ts` は注釈だけ。
- 試験：単体 `tests/unit/engine-memory.test.ts`（引き下げの順、`RangeError` を出すエンジンと伸びない memory を作るエンジン、1 GiB 未満での失敗、WASI の shim と ABI の binding の高い番地）、E2E `tests/e2e/engine-memory.spec.ts` と harness の worker（`high-memory-worker.ts`）。`engine-runtime.test.ts` は最大値を、固定の 16,384 でなく build が記録した値と比べる。

旧い reactor（16,384 ページを宣言）でも新しい reactor（65,535 ページ）でも動く。

## 検証（`7172aaa0`、2026-10-12、Node v26.8.2、Playwright の Chromium・Firefox・WebKit）

reactor と native は、エンジン側の `out-new`（65,535 ページ）と `out-base`（16,384 ページ、`141615ec`）。手順は [`verify_app.sh`](verify_app.sh)（パスはこの機械のもの）。

| 項目 | 結果 |
|---|---|
| `npm run check`（typecheck、lint、単体、build） | 通過。単体は 916 件通過、69 件 skip（FakeEngine の build。エンジンの要る試験は skip）。エンジンの build を与えた単体は、新しい reactor で 983 件、旧い reactor で 983 件通過（セッション 1） |
| 新しい reactor の E2E | 327 件通過、失敗 0（skip 189） |
| 旧い reactor の E2E | 324 件通過、失敗 0（skip 192。`engine-memory.spec.ts` の 3 件は、宣言が 2 GiB 以下なので skip） |
| FakeEngine の E2E | 296 件通過、1 件失敗（skip 219）。失敗は WebKit の `offline.spec.ts:417`（試験 6：固定した release をサイトがもう配らないときの boot の画面）で、452 行目の「Try again」のボタンを待って 180 秒で時間切れ。この試験はエンジンのコードを通らず、この変更は boot と Service Worker を触っていない。同じ試験は、続けて走らせた新旧の reactor の E2E の WebKit で通り（3.3 秒、3.2 秒）、FakeEngine の build で単独に 5 回繰り返すと 5/5 で通った（16.4 秒）。再現しないので原因は特定していない |
| 3 ブラウザでの共有メモリ（`engine-memory.spec.ts`、新しい reactor、4 スレッド） | Chromium・Firefox・WebKit とも、宣言どおり 65,535 ページの最大値で共有メモリを確保した。2.06 GiB（2,214,592,512 バイト、33,792 ページ）から始めた memory で、2 GiB より下の空きを使い切った後の最初の確保の番地は 2,214,592,528（2 GiB 以上）。その状態での 4 プログラム（blastn、blastp、tblastn、tblastx）の検索は、すべての stream（outfmt 0・6・7 と HSP の記録）が通常の memory での実行と同じ bytes で、outfmt 0・6・7 は自己試験の表の SHA-256 とも同じ。3 ブラウザの間でも同じ。記録：[`engine-memory-chromium.json`](engine-memory-chromium.json)、[`engine-memory-firefox.json`](engine-memory-firefox.json)、[`engine-memory-webkit.json`](engine-memory-webkit.json) |
| GitHub の CI（Web、run 38155451498） | engine-web-api・adapter・app の job は成功。browser-engine は 323 件通過、1 件失敗：WebKit の `results.spec.ts:512`（BLASTN）が、dot plot の canvas の `data-titles` を 1 回だけ読んで null だった。この属性は描画の関数が後から付ける（`src/ui/DotPlot.vue` の `dataset['titles']`）ので、試験の読み方の競合。同じ失敗は、このブランチの前の `feature/losat-web-gui-app` の run 38109093750・38109446588・38114111384・38118880176 にもある。この変更は結果の画面を触っていない |

## 残件

- iPhone・iPad、Android、macOS の Safari、Windows の実機では確かめていない（どの段で確保できるかは端末で違いうる）。
- 1 GiB を確保できない端末で、1 GiB 未満の threaded を許すかどうかは製品の判断（今は変更前と同じく serial）。
- 検索の途中でメモリが足りなくなったときの知らせに、確保できた最大値とその時点の大きさを出す案は入れていない。
- FakeEngine の E2E で 1 回だけ時間切れになった WebKit の `offline.spec.ts` の試験 6（上の表）は、再現せず原因が分かっていない。アプリの track でまた出たら調べる。
- 既存の試験の競合（この変更の外。アプリの track に渡す）：`results.spec.ts` の 784〜785 行目は `data-titles` を待たずに 1 回だけ読む。`expect.poll` か `toHaveAttribute` で待てば直る。
