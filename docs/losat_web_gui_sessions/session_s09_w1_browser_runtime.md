# Session S09 — W1：ブラウザでの実行基盤

## INSTRUCTION PROMPT

LOSAT Web の段階 W1 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S09 の行である。設計は計画 §3、§4.6〜§4.7、§5.5 にある。入口の条件は、S05〜S08 で BLASTN・BLASTP・TBLASTN・TBLASTX の升目が V-ABI を通っていることである（`docs/web/verification_cells.tsv`）。BLASTX は SX で加える（計画 DW-10）。

目的：FakeEngine を本物の Wasm エンジンに置き換え、ブラウザの中で 4 program（BLASTX は SX の後）を serial / threaded の両方で動かす。

このセッションはアプリ側（README の規則 1）で、S10（データ層）の後に行う（計画 DW-7）。S10 が port と試験用の実装で作った部分（Data worker の `scan_*` によるレコード表と `register` の照合、Engine worker から Data worker への結果のチャンクの MessagePort）を本物の reactor と Engine worker につなぎ、S10 の契約試験を本物の実装で再実行する。何をどこにつなぐかは、下の「S10 で作った port と、S09 が本物につなぐもの」にある。

S05 で確定した ABI（`docs/web/abi_v2.md`。version 2）と、それを Node で確かめる方法（V-ABI）の要点：

- 成果物は `web/adapter/tools/build_reactors.py` が作る `losat-web-serial.wasm` と `losat-web-threads.wasm`（どちらも WASI の reactor で、instance 化の後に `_initialize` を 1 回呼ぶ）。ビルドの同一性は `web/adapter/tools/check_build_identity.py` が確かめる。
- threaded の module は `env.memory`（共有、最大 16384 ページ）と `wasi.thread-spawn` と `losat_host.emit` を import する。スレッド用の worker の instance にも `losat_host.emit` を渡す必要があるが、呼ばれることはない（エンジンは、export を呼んだスレッドで整形する）。呼ばれたら例外にしてよい。Node では `LOSAT/tests/wasi_thread_host.js` の `createThreadHost(path, "threaded-reactor", [], { imports: { losat_host: { emit } } })` がこれを行う。
- `emit(stream, ptr, len)` のバイトは呼出しの間だけ有効なので写す。stream 0/6/7 は各形式の出力（1 MiB ごと）、1 は HSP レコード（JSON Lines、`out6`・`out0`・`out0_subject` の範囲つき）、2 は `describe`・`register`・`scan_end` の JSON、3 は警告。
- V-ABI の実行：`node web/adapter/tests/v_abi.js --native <LOSAT> --serial <…> --threads <…> --cases <cases.json>`（case は `web/adapter/tools/v_abi_cases.py --suite quick|full`）。ABI の結合の書き方（引数の確保と解放、エラーの読み方、ハンドルの扱い）は `v_abi.js` の `Reactor` がそのまま参考になる。

1. `web/app/src/infra/engine-worker/` に Engine worker を作る。WASI shim（gbdraw と同じ `@bjorn3/browser_wasi_shim`。版は固定する）、ABI v2 の結合、`losat_host.emit` の受け口を置く。`EngineGateway`（`src/ports/engine.ts`。S10 で `run(request, output, onPhase)` になり、`output` は Data worker の MessagePort）を実装し、`src/composition.ts` で FakeEngine と入れ替える。FakeEngine は単体試験のために残す。
2. ThreadHost を作る。`wasi.thread-spawn` を受けて、スレッド用の worker で module を共有メモリ付きで instance 化し、`wasi_thread_start` を呼ぶ。スレッド用の worker は再利用する。参考にする実装は、Node の `LOSAT/tests/wasi_thread_host.js` と、gbdraw の `gbdraw/web/js/workers/` と `gbdraw/web/js/services/losat.js`（<https://github.com/satoshikawato/gbdraw>、commit `538e9ec5`）である。ABI v2 に合わせて書き、コードを写さない。
3. 機能を実際に確かめて経路を選ぶ（`crossOriginIsolated`、`SharedArrayBuffer`、module が宣言する最大値での共有メモリの確保。最大値は認証済みの threaded ビルドと同じ 1 GiB で、host は大きくしない。計画 TD-7）。threaded にできないときは、argv の `-num_threads` を 1 にして serial のモジュールで実行し、理由を `RuntimeInfo.fallbackReason` に入れる（計画 §4.7）。Auto のスレッド数は、gbdraw の既定値（FASTA が 500,000 文字未満なら serial）から始めて測り直す。
4. 取消は、Engine worker とすべてのスレッド用 worker を終了して行う。終了の後に `runtimeGeneration` を増やし、古い世代のメッセージを捨てる。次の実行では instance を作り直し、Subject を原バイトから登録し直す（R1）。
5. instance を作り直す条件（実行回数の上限、またはメモリの高水位）を実測して決める。reactor を繰り返し呼ぶとメモリが増える既知の問題（`docs/wasm_reactor_memory_followup_20260914.md`）の推移を RunRecord に記録する。
6. Chromium で、共有メモリが別の isolate で伸びた後の境界の問題（Node では `LOSAT/tests/wasi_shared_memory.js` が回避している）が起きるかを試し、要るなら同じ変換を使う。
7. V-BR：実アプリの EngineGateway を通して、`verification_cells.tsv` の升目の一部（4 program × 対応する形式 × スレッド 1/2/4）を Chromium・Firefox・WebKit で実行し、期待値の SHA-256 と比べる。取消 → 次の実行の成功、instance の作り直しの後の一致も試す。Firefox と WebKit の Playwright のブラウザは CI で入れる。
8. DW-8 の判断材料：program ごとに、warm 実行の時間のうち前処理（エンコード、翻訳、lookup の構築のうち subject だけで決まる部分）が占める割合を、1 回の暖機と 3 回の計測で測る。20% 以上の program があれば、README の表に R2 の行（program ごと）を足し、その指示書を書く。取消の後の再準備の時間も測り、協調取消の閾値の案を保守者に示す。

9. S10 の port を本物の実装につなぎ、S10 の契約試験を本物の実装で再実行する（次の節）。

完了条件は計画 §7 の S09 の行による。記録は `docs/evidence/losat_web_w1/README.md`。S10 の契約試験（`record-scanner`・`engine-input`・`run-output`）が本物の reactor と Engine worker で通ることも、このセッションで確かめる（S10 のゲート記録で、計画 §7 の S09 の行に足すよう依頼した）。

## 前回の中断（2026-10-01）

最初の S09 の作業は使用量の上限で止まった。`/mnt/c/Users/genom/GitHub/LOSAT-web-gui-app` の作業ツリーに、未コミットで試験していない変更が残っている（`web/app/` の 25 ファイル：Engine worker・reactor の結合・harness など。ブランチは `origin/feature/losat-web-gui` を取り込み済みで、push していない）。最初の作業として、この変更を読み、使うか捨てるかを決める（`git stash` で退避してもよい）。入口の条件の判断（S08 の TBLASTX の outfmt 0/7 の升目を除いて、V-ABI を通る升目で始める）は前回と同じ推奨案とする。使用量を抑えるため、機械的な作業は Agent の `model: "sonnet"` に回す（[S07++b](session_s07ppb_e2f_close.md) の「進め方」）。

## S10 で作った port と、S09 が本物につなぐもの

S10（W2、[ゲート記録](../evidence/losat_web_w2/README.md)）は、本物の reactor に依存する部分を port と試験用の実装で作り、契約試験を書いた。どれも `web/app/` の中にある。

| port（interface） | S10 の試験用の実装 | 契約試験 | S09 がつなぐ本物の実装 |
|---|---|---|---|
| `RecordScanner`（`src/ports/scan.ts`）：ABI v2 の `scan_begin` / `scan_chunk` / `scan_end` | `FakeScanner`（`src/infra/fake/fake-fasta.ts`） | `tests/contract/record-scanner.contract.ts`（23 case） | Data worker の中の serial の adapter reactor |
| `EngineGateway.run(request, output, onPhase)`（`src/ports/engine.ts`）の `EngineInput.records` と `InputMismatchError`：`register` の照合 | `FakeEngine`（`src/infra/fake/fake-engine.ts`） | `tests/contract/engine-input.contract.ts`（4 case） | Engine worker の中の `register` の直後 |
| run の出力の経路（`src/ports/run-output.ts`、`RunOutputWriter` / `RunOutputReceiver`）：Engine worker から Data worker へ、worker どうしの MessagePort | FakeEngine（main thread）と、E2E の engine double worker（`tests/e2e/harness/engine-double-worker.ts`） | `tests/contract/run-output.contract.ts`（11 case） | Engine worker の `emit` の受け口 |

S10 で本物として動いているもの（S09 は変えずに使う）：Data worker（`src/infra/data-worker/`、`startDataWorker()`）、`DataService`（`src/infra/data/data-service.ts`）、BlockStore の OPFS / Memory の実装と契約試験（`tests/contract/block-store.contract.ts`）、Web Locks による所有と回収（`src/infra/data/session.ts`）、容量不足の扱い、使用量の表示。

S09 の作業：

1. **索引の走査**：Data worker の中で serial の reactor（`losat-web-serial.wasm`）を instance 化し、`RecordScanner` を実装する。`scan(parser, chunks)` は `losat_web2_scan_begin(parser)`、chunk ごとの `losat_web2_scan_chunk`、`losat_web2_scan_end` を呼び、stream 2 の *scan* JSON（ABI v2 §9）をそのまま `ScanResponse` として返す。失敗は `losat_web2_last_error_*` の文で reject する（S10 の Data worker は、その文を `Query FASTA: …` の形で利用者に出す）。`src/infra/data-worker/data-worker.ts` の `new FakeScanner()` を置き換える。Data worker の instance は Engine worker と別に持つ（計画 §3.1）。WASI shim と ABI の結合（確保・解放、誤りの読み方、stream 2 の受け取り）は、Engine worker と同じモジュールを使う（DRY）。
2. **`register` の照合**：Engine worker の中で `losat_web2_register` の *register* JSON の `records`（`id`、`length`）を、`request.query.records` と `request.subject.records` と `recordMismatch`（`src/domain/dataset.ts`）で比べる。違えば `InputMismatchError(role, detail)` で reject し、`losat_web2_run` を呼ばない（出力を 1 つも送らない）。main thread と Engine worker の間でも、誤りの `name` と `message` を保つ（`src/infra/data-worker/rpc.ts` と同じ規則）。契約試験はこの `name` を確かめる。
3. **出力の経路**：`output`（Data worker が `openRun` で作った MessagePort）を Engine worker へ transfer する。Engine worker は `RunOutputWriter`（`src/infra/run-output/writer.ts`）を作り、`emit(stream, ptr, len)` の stream 0・6・7・1・3 を `write(stream, view)` にそのまま渡す（`write` が写して transfer するので、呼出しの間だけ有効な wasm メモリの view を渡してよい。stream 2 は渡さない）。`losat_web2_run` が 0 を返した後、main thread に成功を返す**前に** `end()` を呼ぶ。失敗と取消では `end()` を呼ばない（coordinator が `discardRun` を呼び、Data worker が捨てる。S10 で実装済み）。HSP レコード（stream 1）は JSON Lines のまま保存され、`readHits` が読むときに解く。
4. **契約試験の再実行**：(a) `record-scanner.contract.ts` を本物の scanner で通す（Node なら `node:wasi` で serial の reactor を読み込む env、ブラウザなら `tests/e2e/harness/` に scanner の env を足す）。期待値の権威は adapter の scan（`web/adapter/src/scan.rs`、`web/adapter/tests/scan_properties.rs`。計画 TD-8）で、食い違ったらその根拠で期待値か FakeScanner を直す。kind 1 を拒否する case は SX で置き換える。(b) `engine-input.contract.ts` を本物の `EngineGateway` と、本物の argv と FASTA（2 レコード以上）で通す。(c) `run-output.contract.ts` を、writer を本物の Engine worker に置いた env で通す（E2E の `EngineDouble` の代わり）。V-BR のバイト比較は、Data worker に commit された出力（`readOutput`）で行う。
5. **取消と世代**：Engine worker を終了したら、coordinator は今のとおり `discardRun` を呼ぶ。port の相手が消えても Data worker は `end` を待ち続けない（discard が閉じる）。
6. **Subject の保持（R1）**：S10 の `RunSnapshot` は、役割ごとに入力のバイト、SHA-256、`revisionIds`、`records` を持つ。登録済みのハンドルを再利用するときは、`revisionIds` と SHA-256 を鍵にする。
7. **ブラウザ**：Firefox と WebKit を Playwright の project に足す。`tests/e2e/storage.spec.ts` と `tests/e2e/contracts.spec.ts` は、Chromium の CDP（クォータ、凍結）、ディスク上の Chromium のプロファイル（`tests/e2e/support/profile.ts`）と `/proc`（renderer の強制終了）を使うので Chromium だけで動く。Firefox と WebKit では、BlockStore と run output の契約のうち容量不足以外の case と、通常の E2E を通す（その env を harness に足す）。OPFS が使えないブラウザでは Memory に切り替わり、`Temporary storage` の欄に理由が出ることを確かめる。
8. **Memory の上限**：`MemoryBlockStore` は既定で上限が無い（`setCapacity` は試験で使う）。S09 のメモリの実測から、OPFS が使えないときの上限を決めるか、無制限のままにするかを決めて記録する。
9. **空白だけの BLASTN の入力**：adapter の `register` は BLASTN の空白だけの入力を「レコードが無い」として受け付け、CLI は NCBI の警告（query）か誤り（subject）を出す（`docs/web/abi_v2.md` §4）が、`scan` は `Expected > at record start.` で失敗するので、S10 のアプリはキューに入れる前に拒否する。推奨案：adapter の `scan` の扱いに合わせて決め、アプリで BLASTN の規則を作り直さない（エンジン側の変更が要るなら、エンジン側のセッションに頼む）。

S10 で分かった注意：

- Playwright の既定の context（incognito）の OPFS はメモリの上にあり、CDP の `Storage.overrideQuotaForOrigin` を無視する。容量不足の試験はディスク上のプロファイルで行う。CDP の上書きは、その DevTools の session を detach すると消える。`navigator.storage.estimate()` の usage は、Chromium が書き込みで確かめる usage より遅れるので、容量不足を起こすときはクォータを 0 にする。
- この環境（WSL2）では CDP の `Page.crash` が renderer を止めない（`crash` イベントが来ず、ロックも残る）。強制終了は renderer のプロセスに SIGKILL を送って起こす。
- 同期アクセスハンドルを `close` した直後の `removeEntry` は、ファイルがまだ使われているとして `NoModificationAllowedError` になることがある。OPFS の実装は、同じディレクトリの削除を順に行い、短い間隔で試し直す。

## 終了・引き継ぎ

README の規則 8 に従う（アプリ側：`feature/losat-web-gui` への merge と、計画・README の表の更新は、エンジン側が行う）。アプリ側の次は、R2 の行を足した場合はその最初のもの、そうでなければ [S12 — 検索画面](session_s12_w3_search_ui.md)（入口の条件は計画 §7 のとおり。S10 は完了）。実測値（メモリの上限、Auto の閾値、作り直しの条件、前処理の割合）を、後の指示書と計画 §10 に反映する（計画は、ゲート記録を通してエンジン側に依頼する）。
