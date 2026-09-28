# Session S09 — W1：ブラウザでの実行基盤

## INSTRUCTION PROMPT

LOSAT Web の段階 W1 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S09 の行である。設計は計画 §3、§4.6〜§4.7、§5.5 にある。入口の条件は、S05〜S08 で BLASTN・BLASTP・TBLASTN・TBLASTX の升目が V-ABI を通っていることである（`docs/web/verification_cells.tsv`）。BLASTX は SX で加える（計画 DW-10）。

目的：FakeEngine を本物の Wasm エンジンに置き換え、ブラウザの中で 4 program（BLASTX は SX の後）を serial / threaded の両方で動かす。

S05 で確定した ABI（`docs/web/abi_v2.md`。version 2）と、それを Node で確かめる方法（V-ABI）の要点：

- 成果物は `web/adapter/tools/build_reactors.py` が作る `losat-web-serial.wasm` と `losat-web-threads.wasm`（どちらも WASI の reactor で、instance 化の後に `_initialize` を 1 回呼ぶ）。ビルドの同一性は `web/adapter/tools/check_build_identity.py` が確かめる。
- threaded の module は `env.memory`（共有、最大 16384 ページ）と `wasi.thread-spawn` と `losat_host.emit` を import する。スレッド用の worker の instance にも `losat_host.emit` を渡す必要があるが、呼ばれることはない（エンジンは、export を呼んだスレッドで整形する）。呼ばれたら例外にしてよい。Node では `LOSAT/tests/wasi_thread_host.js` の `createThreadHost(path, "threaded-reactor", [], { imports: { losat_host: { emit } } })` がこれを行う。
- `emit(stream, ptr, len)` のバイトは呼出しの間だけ有効なので写す。stream 0/6/7 は各形式の出力（1 MiB ごと）、1 は HSP レコード（JSON Lines、`out6`・`out0`・`out0_subject` の範囲つき）、2 は `describe`・`register`・`scan_end` の JSON、3 は警告。
- V-ABI の実行：`node web/adapter/tests/v_abi.js --native <LOSAT> --serial <…> --threads <…> --cases <cases.json>`（case は `web/adapter/tools/v_abi_cases.py --suite quick|full`）。ABI の結合の書き方（引数の確保と解放、エラーの読み方、ハンドルの扱い）は `v_abi.js` の `Reactor` がそのまま参考になる。

1. `web/app/src/infra/engine-worker/` に Engine worker を作る。WASI shim（gbdraw と同じ `@bjorn3/browser_wasi_shim`。版は固定する）、ABI v2 の結合、`losat_host.emit` の受け口を置く。`EngineGateway`（`src/ports/engine.ts`）を実装し、`src/composition.ts` で FakeEngine と入れ替える。FakeEngine は単体試験のために残す。
2. ThreadHost を作る。`wasi.thread-spawn` を受けて、スレッド用の worker で module を共有メモリ付きで instance 化し、`wasi_thread_start` を呼ぶ。スレッド用の worker は再利用する。参考にする実装は、Node の `LOSAT/tests/wasi_thread_host.js` と、gbdraw の `gbdraw/web/js/workers/` と `gbdraw/web/js/services/losat.js`（<https://github.com/satoshikawato/gbdraw>、commit `538e9ec5`）である。ABI v2 に合わせて書き、コードを写さない。
3. 機能を実際に確かめて経路を選ぶ（`crossOriginIsolated`、`SharedArrayBuffer`、module が宣言する最大値での共有メモリの確保。最大値は認証済みの threaded ビルドと同じ 1 GiB で、host は大きくしない。計画 TD-7）。threaded にできないときは、argv の `-num_threads` を 1 にして serial のモジュールで実行し、理由を `RuntimeInfo.fallbackReason` に入れる（計画 §4.7）。Auto のスレッド数は、gbdraw の既定値（FASTA が 500,000 文字未満なら serial）から始めて測り直す。
4. 取消は、Engine worker とすべてのスレッド用 worker を終了して行う。終了の後に `runtimeGeneration` を増やし、古い世代のメッセージを捨てる。次の実行では instance を作り直し、Subject を原バイトから登録し直す（R1）。
5. instance を作り直す条件（実行回数の上限、またはメモリの高水位）を実測して決める。reactor を繰り返し呼ぶとメモリが増える既知の問題（`docs/wasm_reactor_memory_followup_20260914.md`）の推移を RunRecord に記録する。
6. Chromium で、共有メモリが別の isolate で伸びた後の境界の問題（Node では `LOSAT/tests/wasi_shared_memory.js` が回避している）が起きるかを試し、要るなら同じ変換を使う。
7. V-BR：実アプリの EngineGateway を通して、`verification_cells.tsv` の升目の一部（4 program × 対応する形式 × スレッド 1/2/4）を Chromium・Firefox・WebKit で実行し、期待値の SHA-256 と比べる。取消 → 次の実行の成功、instance の作り直しの後の一致も試す。Firefox と WebKit の Playwright のブラウザは CI で入れる。
8. DW-8 の判断材料：program ごとに、warm 実行の時間のうち前処理（エンコード、翻訳、lookup の構築のうち subject だけで決まる部分）が占める割合を、1 回の暖機と 3 回の計測で測る。20% 以上の program があれば、README の表に R2 の行（program ごと）を足し、その指示書を書く。取消の後の再準備の時間も測り、協調取消の閾値の案を保守者に示す。

完了条件は計画 §7 の S09 の行による。記録は `docs/evidence/losat_web_w1/README.md`。

## 終了・引き継ぎ

README の規則 8 に従う。次は R2 の行を足した場合はその最初のもの、そうでなければ [S10 — データ層](session_s10_w2_data_layer.md)。実測値（メモリの上限、Auto の閾値、作り直しの条件、前処理の割合）を、後の指示書と計画 §10 に反映する。
