# LOSAT Web W1（Session S09）ゲート記録

- 段階：W1 ブラウザでの実行基盤（[総合計画書](../../losat_web_gui_plan.md) §7 の S09、[指示書](../../losat_web_gui_sessions/session_s09_w1_browser_runtime.md)）
- ブランチ：`feature/losat-web-gui-app`（アプリ側。worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui-app`）。エンジンは `origin/feature/losat-web-gui` の `dfa65cfd8`（エンジンの最後のコミットは `24fcfe41b`）を早送りで取り込み、ゲートの前に文書だけの 2 コミット（`27cfdefd1` まで）を merge した（`bb6cdaad5`。衝突なし）
- 実行記録：[`run-20261003T053627Z/`](run-20261003T053627Z/)（commit `3669ad28f` の木。ゲートの記録）と [`run-20261003T043424Z/`](run-20261003T043424Z/)（commit `bb6cdaad5` の木。この実行の計測で TBLASTN のスレッドの不具合を見つけた。下の「計測で見つけて直したもの」）。作成後は書き換えない。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)、再現は [`run_gate.sh`](run_gate.sh)
- 判定：**完了条件を満たした**。2 回目のゲートの実行（`run_gate.sh`）のすべての段階が通った：V-ABI quick（52 の検索がネイティブの CLI と一致）、`npm run check`、reactor 付きの単体試験（184 件）、FakeEngine と本物のエンジンの両方のビルドでの E2E（3 ブラウザ。26 件と 56 件）、エンジンの E2E の 2 回の繰り返し（60 件）、共有メモリの試験、計測
- 保守者の判断待ち：2 件（下の「保守者の判断待ち」。TBLASTN の R2 の進め方、協調取消の閾値）

## コミット

| コミット | 内容 |
|---|---|
| `0b65619bf` | アプリ：W1 の実行基盤（Engine worker、ThreadHost、ABI v2 の結合、Data worker の本物の走査、V-BR、計測、CI） |
| `bb6cdaad5` | `origin/feature/losat-web-gui` の merge（文書だけ） |
| `3669ad28f` | アプリ：一つの検索でスレッドのプールを次々に作る場合の予備の thread worker。E2E のポートを既存のサーバーと共有しない（ゲートの記録の木） |
| `a1157dd02` | アプリ：`policy.ts` のコメントだけ（ゲートの計測の値に合わせた）。動作は `3669ad28f` と同じ |
| このゲート記録を含むコミット | 文書：このゲート記録、`run_gate.sh`、実行記録、共有メモリの試験、TBLASTN の profile、S09+（R2）の指示書、S12 の指示書への申し送り |

エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表、`docs/web/` は変えていない。依存ライブラリは `@bjorn3/browser_wasi_shim` 0.4.2（gbdraw と同じ WASI shim、版を固定）を足した。

## 完了条件と結果

| 完了条件（計画 §7 の S09） | 結果 | 証拠 |
|---|---|---|
| V-BR（BLASTX を除く 4 program、3 ブラウザ、n=1/2/4） | 通過。9 つの検索 × 1/2/4 スレッド = 27 の検索を、Chromium 149・Firefox 151・WebKit 26.5 のそれぞれで、アプリ（Coordinator、Data worker、`WasmEngine`、Engine worker）を通して実行し、commit された outfmt 0・6・7 と診断がすべてネイティブの CLI と一致した。NCBI の凍結バイトとの比較 21（認証済み 15、TBLASTX 0/7 の記録だけ 6）はすべて一致。1 スレッドは serial、2・4 は threaded の経路。HSP レコードの数と範囲も outfmt 6 の行・outfmt 0 と一致 | `run-…/npm-e2e-engine.log`、`e2e-engine-repeat.log`、`records/v-br-*.json` |
| 取消の後の実行が成功する | 通過。実行中の TBLASTX（ゲノムどうし）を取消すと、次の検索は新しい世代（`runtimeGeneration` 2）で、1 と 4 スレッドとも成功し、出力はネイティブの CLI と NCBI の凍結バイトに一致（3 ブラウザ）。作り直し（毎回）の後も一致し、Subject は毎回登録し直される | `records/cancel-*.json`、`measure/measure-cancel-*.json` |
| メモリの推移と前処理の割合を記録する | 記録した（下の「実測」） | `measure/` |
| S10 の契約試験（`record-scanner`・`engine-input`・`run-output`）が本物の reactor と Engine worker で通る | 通過。`record-scanner` 23 case を本物の serial の reactor でブラウザ（3 ブラウザ）と Node、`engine-input` 4 case を本物の `EngineGateway` で n=1 と n=4（3 ブラウザ）、`run-output` の容量不足以外の 10 case を書き手を本物の Engine worker に置いて（3 ブラウザ）。容量不足の case は、S10 のとおり Chromium の `contracts.spec.ts` で通る | `unit-cases.log`、`npm-e2e-engine.log` |

## 前回の中断の変更

最初の S09（2026-10-01、使用量の上限で中断）の未コミットの変更（`web/app/` の 25 ファイル）を読み、設計が計画 §3.1・§4.6・§4.7 に沿っていたので**使った**。そのときの記録では、Chromium で V-BR が通り、Firefox と WebKit では threaded の検索が serial に切り替わっていた。このセッションで、原因（下の判断 2）を直し、途中だった `serveHarness` の引数を整え、デバッグ用の 2 ファイル（`.try.mjs`、`zz-debug.spec.ts`）を消した。エンジンは `dfa65cfd8` を早送りで取り込み（アプリ側のブランチには独自のコミットが無かったので merge コミットは無い）、reactor とネイティブの CLI をその木から作り直した。

## 作業との対応

指示書の作業 1〜9：

1. Engine worker（`src/infra/engine-worker/`）、WASI shim（`src/infra/reactor/wasi.ts`）、ABI v2 の結合（`src/infra/reactor/abi.ts`）、`losat_host.emit` の受け口。`EngineGateway` を `WasmEngine`（`gateway.ts`）で実装し、`src/composition.ts` で入れ替えた。FakeEngine は、reactor を入れないビルドと単体試験のために残した。
2. ThreadHost（`thread-host.ts`、`thread-worker.ts`）。thread worker は前もって module を共有メモリで instance 化して待ち、スレッドが終わると instance を作り直して再利用する。N スレッドの検索に N−1 個を 2 組用意する（「計測で見つけて直したもの」）。gbdraw と Node の host（`LOSAT/tests/wasi_thread_host.js`）を参考に、ABI v2 に合わせて書いた。
3. 機能の確認と serial への切り替え（`runtime.ts` の `threadedSupport`、module が宣言する最大値（16384 ページ）での共有メモリの確保）。Auto のスレッド数は実測で決めた（判断 5）。
4. 取消は Engine worker（とその thread worker）の終了。`runtimeGeneration` を増やし、古い世代のメッセージを捨てる。次の実行は instance を作り直し、Subject を原バイトから登録し直す。
5. 作り直しの条件を実測で決めた（判断 6）。RunRecord に `memory`（検索の前後の linear memory、instance の実行回数）を記録する。
6. Chromium で共有メモリの境界の問題を試し、起きたので同じ変換を使った（判断 3）。
7. V-BR を 3 ブラウザで実行（下の「実測」）。取消 → 次の実行、作り直しの後の一致も試した。CI で Firefox と WebKit を入れる（`.github/workflows/web.yml`）。
8. DW-8 の前処理の割合と、取消の後の再準備の時間を測った。TBLASTN が 20% 以上なので、S09+（R2）の指示書を書いた（README の表の行は、規則 1 によりエンジン側が足す）。
9. S10 の port を本物につなぎ、契約試験を本物で再実行した（次の表）。

S10 の port と本物の実装（指示書の「S10 で作った port と、S09 が本物につなぐもの」の 1〜9）：

| 作業 | 結果 |
|---|---|
| 1 索引の走査 | Data worker が自分の serial の reactor（Engine worker とは別の instance）で `RecordScanner`（`src/infra/reactor/scanner.ts`）を実装し、FakeScanner を置き換えた。WASI と ABI の結合は Engine worker と同じモジュール。`describe` と `validate` も Data worker の serial の reactor が答える（検索の実行中でも答えられる） |
| 2 `register` の照合 | Engine worker が `register` の直後に `recordMismatch` で比べ、違えば `InputMismatchError(role, detail)` で、`losat_web2_run` を呼ばずに失敗する。誤りの `name`・`role`・文は main thread で作り直す |
| 3 出力の経路 | `output` を Engine worker に transfer し、`RunOutputWriter` に stream 0・6・7・1・3 を渡す（stream 2 は渡さない）。`end()` は `losat_web2_run` が 0 を返し、query のハンドルを解放した後、成功を返す前に 1 回だけ呼ぶ。失敗と取消では呼ばない |
| 4 契約試験の再実行 | (a) `record-scanner`（23 case）：ブラウザの dedicated worker と Node（`tests/unit/engine-runtime.test.ts`）で本物の serial の reactor。(b) `engine-input`（4 case）：本物の `EngineGateway` で n=1 と n=4、2 レコードずつの BLASTN。(c) `run-output`（容量不足を除く 10 case）：書き手を本物の Engine worker（試験用のビルドだけのフック）に置いた。V-BR のバイト比較は、Data worker に commit された出力（`exportOutput` / `readOutput`）で行う |
| 5 取消と世代 | Engine worker の終了の後、coordinator は `discardRun` を呼び、Data worker は `end` を待たずに捨てる（取消の E2E が止まらずに通る） |
| 6 R1 | 鍵は program、SHA-256、`revisionIds`（判断 9）。保持した Subject を使った検索は `RunRecord.subjectRetained` が `true` |
| 7 ブラウザ | Firefox と WebKit を Playwright の project に足した。`storage.spec.ts` と `contracts.spec.ts` は Chromium だけ。`browsers.spec.ts` が Firefox と WebKit で BlockStore（OPFS・Memory）と run の出力の契約（容量不足を除く）を通し、アプリの保存先を確かめる。Playwright の WebKit には OPFS が無く、アプリは Memory に切り替わり、`Temporary storage` の欄に理由（「this browser has no Origin Private File System」）を出す |
| 8 Memory の上限 | 512 MB にした（判断 7） |
| 9 空白だけの BLASTN の入力 | `scan` の扱いに合わせる（判断 8） |

## 作ったもの

```text
Main thread ── Coordinator ── DataGateway ── rpc ──► Data worker ─ DataService / BlockStore / Web Lock
                    │                                     └ serial reactor：scan（RecordScanner）、describe、validate
                    │ EngineGateway = WasmEngine（世代、取消、作り直し、Auto）
                    ▼
              Engine worker（世代ごと）─ EngineRuntime：warm な instance（serial | threaded）、保持した Subject
                    ├ ThreadHost ─► thread worker × 2(n−1)（前もって instance 化、再利用。次のプールのための予備の組）
                    └ RunOutputWriter ──MessagePort（openRun）──► Data worker
```

| 場所 | 内容 |
|---|---|
| `build/reactors.ts` | `LOSAT_WEB_REACTORS` の 2 つの module を `artifacts.json` の SHA-256 で確かめてビルドに入れ、`virtual:losat-engine` で URL と SHA-256 を渡す。threaded の module には `LOSAT/tests/wasi_shared_memory.js` の guard をかける（判断 3） |
| `src/infra/reactor/` | `abi.ts`（ABI v2 の結合：引数の確保と解放、誤りの文、stream の受け口、1 MiB ごとに届く stream 2 の結合）、`instance.ts`（取得と SHA-256 の確認、instance 化）、`artifact.ts`（module の種類と memory の上限をバイトから読む）、`wasi.ts`、`scanner.ts`、`control.ts`、`assets.ts` |
| `src/infra/engine-worker/` | `gateway.ts`（main thread の `WasmEngine`）、`engine-worker.ts`、`runtime.ts`、`thread-host.ts`、`thread-worker.ts`、`protocol.ts`、`policy.ts`（Auto と作り直しの値） |
| `src/infra/data/` | Memory の結果の上限と、その容量不足の文 |
| `tests/e2e/` | `engine.spec.ts`（契約、V-BR、R1、取消、作り直し、serial への切り替え）、`browsers.spec.ts`、`measure.spec.ts`（`LOSAT_WEB_MEASURE` のときだけ）、`support/harness-server.ts`（判断 2）、`support/native.ts`（V-BR の期待値） |
| `tests/unit/` | `engine-gateway.test.ts`（偽の worker で世代・取消・作り直し・障害）、`thread-host.test.ts`（起動の手順と時間切れ）、`engine-runtime.test.ts`（Auto と作り直しの規則、module の検査、stream 2 の結合、reactor があれば Node での走査の契約と 1 MiB を超える応答） |
| `.github/workflows/web.yml` | app の job で 3 ブラウザを入れる。新しい `browser-engine` の job が adapter の job の reactor とネイティブの CLI を受け取り、エンジン入りのビルドで E2E 全体と Node の reactor の試験を行う |

## 判断（推奨案で進めたもの）

1. **前回の変更を使う**（上の「前回の中断の変更」）。
2. **E2E のハーネスは本物の HTTP サーバーで配る**（`tests/e2e/support/harness-server.ts`）。Playwright の route による配信では、Firefox で worker が作る worker（thread worker）の要求が route を通らず `vite preview` の HTML が返り、WebKit では route で返した worker が cross-origin isolated にならない（`SharedArrayBuffer` が無い）。Node の HTTP サーバーで同じ応答ヘッダー（`public/_headers`）を全ての応答に付けると、3 ブラウザとも worker の中で `crossOriginIsolated`、`SharedArrayBuffer`、入れ子の module worker、最大 1 GiB の共有メモリ、compile 済み module と共有メモリの postMessage が使えた。アプリの欠陥ではなく試験の配り方の問題だった。
3. **threaded の module に共有メモリの guard をかける**。下の「共有メモリの境界」のとおり Chromium 149 で起きたので、Node の host と同じ関数（`LOSAT/tests/wasi_shared_memory.js` の `guardSharedMemory`）を、ビルドの時に threaded の module にかける（ブラウザで Node と同じバイトを動かす。検索のコードは変えない）。配る module の SHA-256 は guard の後のもので、RunRecord の `engineBuild` は元の artifact の SHA-256 と「(shared-memory guard)」を示す（例：`losat-web-threads.wasm sha256:a83376643e459cb1 (shared-memory guard)`）。serial の module は共有メモリを使わないので、そのまま配る。
4. **V-BR の期待値**（計画 TD-5、§6.1）：同じ commit のネイティブの CLI を、アプリが実行したのと同じ argv で形式ごとに 1 スレッドで実行した出力（V-ABI と同じ方法。`tests/e2e/support/native.ts`）。`LOSAT/tests/outfmt0_manifest.tsv` の NCBI の凍結バイトがあり、その升目が認証済みのものは、それとも比べる。TBLASTX の outfmt 0/7 は `verification_cells.tsv` で `pending`（S08b が確かめる）なので、凍結バイトとの比較は記録だけにした（すべて一致した）。S08b が `checked` にしたら、`native.ts` の `UNCERTIFIED` から外す。
5. **Auto のスレッド数**（`policy.ts`）：入力（query と subject の FASTA）が 20,000 バイト未満なら serial、それ以上なら論理プロセッサの半分（最大 4）。gbdraw の既定値（500,000 文字未満なら serial）から測り直した（下の「Auto」）。4 スレッドは、仕事が多くのレコードに分かれる検索（BLASTP の proteome どうし）で 2 万残基から約 3 倍速く、1 本の genome どうしの検索では絶対値でわずかしか変わらない（最大で 69 ms 遅い）。gbdraw の値では、20 万残基の BLASTP が serial のまま（13.4 s。4 スレッドなら 4.4 s）になる。上限 4 は V-ABI と V-BR が確かめるスレッド数（1/2/4）に合わせた。手動の指定はそのまま使う。
6. **instance の作り直し**（`policy.ts`）：検索の後の linear memory が 512 MiB（threaded の最大値の半分）以上なら作り直す。実行回数では作り直さない。下の「メモリ」のとおり、linear memory は最初の 1〜3 回の検索で水準に達し、その後は 20 回の繰り返しと、全 program を順に 30 回でも増えなかった。最大の検索（大腸菌の 2 つの genome）で約 240 MiB（serial 237.4、4 スレッド 242.6）。
7. **Memory の上限**：OPFS が使えないときの結果の保存（Memory の BlockStore）を、1 つのタブで 512 MB までにした（`MEMORY_RESULTS_CAPACITY_BYTES`）。その場合メモリが結果の唯一の写しなので、使い切ってタブが落ちるとすべての結果を失う。上限を超えた実行は、OPFS が一杯のときと同じく理由付きで失敗し（「…this browser keeps results in memory, at most 512 MB in a tab. Earlier results are kept.」）、前の結果は残る。実測の出力の大きさ（大腸菌の genome どうしで outfmt 0 が 20.7 MiB、6/7 が 0.5 MiB）なら、約 20 回分が入る。エンジンのメモリ（threaded は最大 1 GiB）は別にかかる。
8. **空白だけの BLASTN の入力**：アダプタの `scan` は「Expected > at record start.」で失敗し、`register` はレコードが無いとして受け付ける（CLI は query なら NCBI の警告「Query is Empty!」、subject なら誤り）。推奨案のとおり `scan` に合わせ、キューに入れる前に `Query FASTA: Expected > at record start.` で拒否する（E2E：`smoke.spec.ts`）。アプリで BLASTN の規則は作り直さない。CLI と同じにするなら、エンジン側で `scan` を変える（下の「エンジン側への依頼」）。
9. **R1 の鍵**は指示書のとおり program、SHA-256、`revisionIds`。今の `enqueue` は入力ごとに新しい revision を作るので、アプリの画面からは再利用が起きない（独立レビューの指摘 M1）。S12 で `DatasetStore` から同じ revision を続く検索に渡すと効く（S12 の指示書に書いた）。レビューは「program と SHA-256 だけを鍵にする」案を勧めた（`register` の結果は program とバイトだけで決まり、レコード表は毎回比べるので、正しさは変わらない）。S12 の後も再利用が起きにくければ、その案に変える。
10. **E2E は既存のサーバーを使わない**：エンジン側も同じ機械の別の worktree で同じ試験を動かす（DW-7）。Playwright がポート 4173 にあった別のサーバー（同じ時間に動いていたエンジン側の試験のものと推定。cross-origin isolated でなく、FakeEngine の帯も無かった）を使い回して試験したことがあった（FakeEngine の E2E が 18 件失敗。`3669ad28f` の前）。既存のサーバーは使わず（ポートが使われていれば失敗する）、`LOSAT_WEB_E2E_PORT` でポートを変えられるようにした（既定は 4173。このセッションのゲートは 4273）。
11. **E2E はどちらのエンジンでも通る**：`vite preview` のビルドは `LOSAT_WEB_REACTORS` を読むので、FakeEngine の出力を前提にした試験を、ビルドにエンジンがあるかで切り替えた（`BUILD_HAS_ENGINE`、`OUTFMT7_MARK`）。
12. **計測のための環境変数は持たない**：エンジンの `LOSAT_TIMING` は、BLASTP では wasm32 で無効（`LOSAT/src/algorithm/blastp/blast_engine.rs:410-418`）、TBLASTN には無いので、reactor に環境変数を渡す仕組みは作らなかった（ABI v2 §3 のとおり空）。

## 実測

計測は `tests/e2e/measure.spec.ts`（ゲートの実行の `measure/`）。warm な検索（instance と登録済みの Subject を再利用）の検索の段階（`losat_web2_run` の時間）を、1 回の暖機と 3 回の計測の中央値で示す。同じ機械で、エンジン側（S08b）の V-ABI full などが並行して動いていた時間がある（32 論理プロセッサ）ので、絶対値には揺れがあり、比べるのは同じ実行の中の比である（1 回目のゲートの実行と、試しの計測でも同じ傾向だった）。

### ブラウザの機能

Playwright 1.61.1 の Chromium 149.0.7827.55、Firefox 151.0、WebKit 26.5（WPE。この機械ではシステムライブラリを root なしで展開した起動 script を `LOSAT_WEB_WEBKIT_EXECUTABLE` で使った。CI は `npx playwright install --with-deps` で入れる）。3 つとも、COOP / COEP を返すサーバーで、worker の中で `crossOriginIsolated`、`SharedArrayBuffer`、入れ子の module worker、最大 16384 ページの共有メモリ、compile 済みの module と共有メモリの postMessage が使えた。OPFS は Chromium と Firefox で使え、Playwright の WebKit の context には無い（アプリは Memory にし、理由を示す）。Firefox は、この機械で同じ検索が Chromium と WebKit の 4〜6 倍遅い（例：BLASTP の V-PERF の fixture で 3.7 s と 0.59 s）。

### 共有メモリの境界（作業 6）

[`shared_memory_growth_probe.mjs`](shared_memory_growth_probe.mjs)（`LOSAT/tests/test_wasi_shared_memory.js` の fixture）：一つの worker が Wasm の中で待つ間に別の worker が共有メモリを伸ばし、伸びたページへの `memory.fill` / `memory.copy` を、各 20 回、guard なしと guard ありで行った（`run-…/shared-memory-growth-probe.log`）。

| ブラウザ | guard なし | guard あり |
|---|---|---|
| Chromium 149 | `memory access out of bounds`：fill 1、copy（書き込み先）1 / 60 回（試しの実行でも 2 / 60、1 回目のゲートで 1 / 60） | 0 / 60 |
| Firefox 151 | 0 / 60 | 0 / 60 |
| WebKit 26.5 | 0 / 60 | 0 / 60 |

### DW-8：subject だけで決まる前処理の割合（作業 8）

V-PERF の program ごとの代表の fixture（`docs/evidence/losat_web_e1a/measure_perf.py`）で、同じ Subject（保持して再利用）に、本当の query と、その先頭の数残基だけの query（qmin：BLASTP と TBLASTN は 30 残基、BLASTN は 60 塩基、TBLASTX は 90 塩基）を当てた。subject だけで決まる仕事（エンコード、翻訳など）は qmin の検索にもすべて含まれ、それに subject の走査と実行の固定費が加わるので、qmin と本当の query の検索の段階の比は、その割合の**上限**である。エンジンの `LOSAT_TIMING` は reactor では使えない（判断 12）。

| program（fixture） | Chromium n=1 / n=4 | Firefox n=1 / n=4 | WebKit n=1 / n=4 | 判定 |
|---|---|---|---|---|
| BLASTP（SicyWSV.faa / PajaWSV.faa、`-max_hsps 1`） | 0.9% / 1.4% | 1.0% / 1.2% | 1.1% / 1.6% | 20% 未満 |
| BLASTN megablast（AP027152 / LC738884） | 6.5% / 7.5% | 5.7% / 5.7% | 7.3% / 7.9% | 20% 未満 |
| TBLASTX（LC738884 / LC741431） | 1.3% / 1.3% | 1.3% / 1.4% | 1.5% / 1.6% | 20% 未満 |
| TBLASTN（Stage G の 1 本の protein / AvCLPV） | 59.9% / 58.3% | 58.8% / 56.3% | 56.5% / 56.4% | 上限が 20% を超える → 下の profile |

TBLASTN の上限は判定に使えないので、ネイティブ（同じ commit）の命令数を callgrind で役割に分けた（[`tblastn_profile/findings.md`](tblastn_profile/findings.md)）。subject だけの前処理（6 フレームの翻訳 `generate_frames` と、subject の ncbi2na への変換）は全体の **40.9%**（下限 36.1%、上限は約 50%）。subject の長さを半分にした実行と、wall time の比べ方も同じ値を示した。**TBLASTN は DW-8 の条件（20% 以上）を満たす**。その命令の 76% は codon ごとの `GeneticCode::get`（codon あたり約 159 命令）で、表引きにすれば約 16% になる見込み（推定）。→ 保守者の判断待ち 1、[S09+（R2）の指示書](../../losat_web_gui_sessions/session_s09p_r2_tblastn_subject_cache.md)。

### メモリ（作業 5）

Chromium、同じ検索を 20 回（入力は毎回索引し直し、Subject は毎回登録し直す。今のアプリと同じ）と、全 program を順に 6 周（30 回）。数字は検索の後の linear memory（MiB）。

| 検索 | スレッド | 1〜3 回目 | 18〜20 回目 | 最大 | outfmt 0 / 6 / 7（MiB） |
|---|---:|---|---|---:|---|
| BLASTP SicyWSV / PajaWSV | 1 | 22.6, 22.6, 22.6 | 22.6, 22.6, 22.6 | 22.6 | 0.4 / 0.0 / 0.1 |
| TBLASTN AvCLPV protein / AvCLPV | 1 | 12.7, 13.1, 13.1 | 13.1, 13.1, 13.1 | 13.1 | 0.0 / 0.0 / 0.0 |
| BLASTN megablast AP027152 / LC738884 | 1 | 73.6, 73.6, 73.7 | 73.7, 73.7, 73.7 | 73.7 | 0.0 / 0.0 / 0.0 |
| TBLASTX LC738884 / LC741431 | 1 | 67.2, 67.2, 67.2 | 67.2, 67.2, 67.2 | 67.2 | 0.8 / 0.2 / 0.2 |
| BLASTN megablast Sakai / MG1655 | 1 | 237.4, 237.4, 237.4 | 237.4, 237.4, 237.4 | 237.4 | 20.7 / 0.5 / 0.5 |
| 全 program を順に | 1 | 22.6, 22.6, 73.9 | 237.2, 237.2, 237.2 | 237.2 | |
| BLASTP SicyWSV / PajaWSV | 4 | 25.6, 25.7, 25.7 | 25.7, 25.7, 25.7 | 25.7 | 0.4 / 0.0 / 0.1 |
| TBLASTN AvCLPV protein / AvCLPV | 4 | 15.9, 15.9, 15.9 | 16.0, 16.0, 16.0 | 16.0 | 0.0 / 0.0 / 0.0 |
| BLASTN megablast AP027152 / LC738884 | 4 | 76.6, 76.6, 76.6 | 76.8, 76.8, 76.8 | 76.8 | 0.0 / 0.0 / 0.0 |
| TBLASTX LC738884 / LC741431 | 4 | 70.2, 71.4, 71.4 | 71.5, 71.5, 71.5 | 71.5 | 0.8 / 0.2 / 0.2 |
| BLASTN megablast Sakai / MG1655 | 4 | 240.3, 240.3, 242.6 | 242.6, 242.6, 242.6 | 242.6 | 20.7 / 0.5 / 0.5 |
| 全 program を順に | 4 | 25.6, 25.6, 78.6 | 242.2, 242.2, 242.2 | 242.2 | |

linear memory は最初の 1〜3 回で水準に達し、その後は増えない（`docs/wasm_reactor_memory_followup_20260914.md` の Node での観察と違い、この範囲では増え続ける系列は無かった）。水準は、その instance で行った最大の検索で決まる。→ 判断 6。

### Auto（作業 3）

Chromium、検索の段階の中央値（ms）。入力は query と subject に半分ずつ（BLASTN：Sakai / MG1655 の先頭、TBLASTX：LC738884 / LC741431 の先頭、BLASTP：NZ_CP006932 の proteome の先頭のレコードどうし、TBLASTN：その proteome の先頭 1/8 と genome の先頭 7/8）。

| 残基（合計） | BLASTN 1 / 2 / 4 | TBLASTX 1 / 2 / 4 | BLASTP 1 / 2 / 4 | TBLASTN 1 / 2 / 4 |
|---:|---|---|---|---|
| 20,000 | 3 / 3 / 3 | 18 / 18 / 17 | 298 / 190 / 105 | 58 / 58 / 57 |
| 60,000 | 6 / 5 / 6 | 41 / 41 / 40 | 1,926 / 1,097 / 629 | 205 / 204 / 203 |
| 200,000 | 15 / 15 / 15 | 199 / 266 / 268 | 13,431 / 7,842 / 4,354 | 1,138 / 1,151 / 1,151 |
| 500,000 | 37 / 38 / 37 | 1,094 / 1,107 / 1,059 | 39,106 / 22,601 / 12,265 | 3,359 / 3,530 / 3,519 |

LOSAT の並列の単位はレコード（の組）なので、1 本の genome どうし（BLASTN、TBLASTX）や、1 本の genome に対する TBLASTN ではスレッドがほとんど効かず、proteome どうしの BLASTP では 4 スレッドで約 3 倍速い。threaded の損は絶対値で小さい（最大 69 ms）。→ 判断 5。

### 取消と再準備（作業 4、8）

実行中の TBLASTX（LC738884 / LC741431）を、検索の段階に入って 500 ms 後に取消し、同じ Subject の BLASTN を新しい runtime で行った。準備の時間（要求から検索の段階まで：worker と instance の起動、query と Subject の登録）、ms。

| ブラウザ | Subject | スレッド | warm（Subject を保持） | 取消の後（新しい runtime、Subject を登録し直す） |
|---|---|---:|---:|---:|
| Chromium | LC738884（365 kb） | 1 / 4 | 3 / 3 | 11 / 27 |
| Chromium | MG1655（4.6 Mb、query Sakai 5.5 Mb） | 1 / 4 | 45 / 47 | 93 / 106 |
| Firefox | LC738884 | 1 / 4 | 16 / 17 | 42 / 58 |
| Firefox | MG1655 | 1 / 4 | 306 / 317 | 576 / 606 |
| WebKit | LC738884 | 1 / 4 | 3 / 11 | 13 / 65 |
| WebKit | MG1655 | 1 / 4 | 45 / 58 | 87 / 139 |

取消そのものは、Engine worker の終了ですぐに効く（取消から取消の状態まで 0〜1 ms）。取消の後の再準備は最大 0.61 s（Firefox、大腸菌の genome）で、warm より最大 0.29 s 長い。→ 保守者の判断待ち 2。

## 独立レビュー

コミット前の変更を、別のエージェントが読み取り専用でレビューした（判定「ready after fixes」）。指摘と対応（すべて `0b65619bf` に含む）：

| 重さ | 指摘 | 対応 |
|---|---|---|
| 高 | stream 2 の応答（*scan*・*register*）も 1 MiB ごとに届くのに、最後の塊だけを読んでいた。約 3,000 レコード以上の走査、約 20,000 レコード以上の登録が失敗し、読めなかった登録のハンドルが残る | 呼出しごとに塊を集めてつなぐ。応答を読めない登録はハンドルを解放する。単体試験（塊に分けた応答）と、本物の reactor で 25,000 レコードの走査と登録（数 MiB の応答）を足した |
| 中 | R1 はアプリの画面からは起きない。作り直しの試験の `subjectRetained` の確認は何も示していない | 判断 9 に記録。作り直しの試験を、同じ revision を使う経路にし、作り直さないときは保持され（`[false, true, false, true]`）、作り直すと毎回登録し直す（`[false, false, false, false]`）ことを確かめる |
| 中〜低 | Data worker の module の取得が 1 回失敗すると、タブを開き直すまで走査・`describe`・`validate` が失敗する | 失敗した取得を覚えず、次の利用で取り直す |
| 低〜中 | 時間切れで諦めた起動の後に thread worker がスレッドを始めると、解放された引数で動き、slot も健全に戻る | 起動を compare-and-swap の手順にした（待機 → 開始 / 放棄。放棄の後は始めない）。時間切れは障害として runtime を終える。単体試験を足した |
| 低 | 一時的な起動の失敗でも、その世代のすべての検索が serial になる | 機能が無い（隔離されていない、共有メモリを確保できない）ときだけ世代を通して serial にし、thread worker の準備の失敗は次の検索で作り直して試す |
| 低 | serial と threaded の切り替えのたびに instance を捨てる | 時間だけの損失（下の「取消と再準備」の作り直しの時間と同じ程度）。記録して、そのままにした |
| 細 | query のハンドルの解放が `end()` の後 | `end()` の前に解放する |
| 整理 | 使っていないコード、古いコメント、`findReactors()` の重複、guard の読み込み | 消した。`findReactors()` は設定で 1 回、guard はエンジン入りのビルドのときだけ読む |
| 試験 | HSP レコードの照合が、ヒットが無いと素通りする | レコードの数 = outfmt 6 の行の数、`out0`・`out0_subject` の範囲が outfmt 0 の中にあることも確かめる |

レビューで正しいとされた点：世代の選り分けと遅れたメッセージ、compile 中と実行中の取消、世代ごとの障害の通知、共有メモリの view の取り直しと写し、確保と解放、`end` と `done` の順、stream 2 を渡さないこと、照合が出力の前であること、走査の誤りの文、`claim` と `spawn` が必ず時間切れで終わること、web/AGENTS.md の規則（BLAST の値を計算しない、層、依存の固定、FakeEngine の帯）、CSP と COOP / COEP、ABI v2（`-num_threads`、`env.memory` の上限、版、thread の `emit`、thread の instance で `_initialize` を呼ばない）、V-BR の比較の方法。

## 計測で見つけて直したもの

TBLASTN は query を 20,000 残基ごとの batch に分け、batch ごとにスレッドのプールを作る（`LOSAT/src/algorithm/tblastn/args.rs:587-589`、`LOSAT/src/utils/threading.rs` の `with_search_pool`）。Auto の計測（62,500 残基の query、4 つの batch）で、新しい runtime の最初の 2 スレッドの検索が「failed to build tblastn pool with 2 threads; … Resource temporarily unavailable (os error 6)」で失敗した（試しの計測で 2 回、最初のゲートの実行 `run-20261003T043424Z/measure/chromium.log` で 1 回）。

- 原因（記録を足して再現した）：前のプールのスレッドは、次のプールの `thread-spawn` が戻るまで JavaScript に戻らない（エンジンのスレッドの library が `thread-spawn` の間に持つ lock を、終わるスレッドが必要とする）。ThreadHost が thread worker を N−1 個しか用意していないと、次のプールの `thread-spawn` は前のスレッドの thread worker を待ち、互いに待って準備の時間切れ（30 s）の後に −1 を返していた。最初の修正（実行中の thread worker も待つ）は、この待ち合いのため効かなかった。
- 修正（`3669ad28f`）：N スレッドの検索に N−1 個の thread worker を 2 組用意し、次のプールは前のプールが使わなかった組で始める。`thread-spawn` は、準備し直している thread worker だけを待ち、スレッドが実行中のものは待たない（残りがそれだけなら、すぐに −1）。Node の host（`LOSAT/tests/wasi_thread_host.js`）は、`thread-spawn` のたびに新しい Node の worker を作るので、この問題は無い。
- 試験：単体試験（予備の組、実行中を待たないこと、準備し直しを待つこと）、E2E（このような検索を、新しい runtime の最初の検索として 2 と 4 スレッドで、3 ブラウザ。出力はネイティブの CLI と一致）、V-BR の case `tblastn.batches`（64,069 残基）、計測の Auto の TBLASTN（2 回目のゲートの実行で失敗なし）。E2E は時間の具合で再現しないことがある（修正の前のコードでも 30 s 遅れて通ることがあった）ので、規則は単体試験で固定した。

## 保守者の判断待ち

1. **TBLASTN の R2 の進め方**（DW-8）：TBLASTN の subject だけの前処理は 20% 以上（上の「DW-8」）なので、S09+（R2）の[指示書](../../losat_web_gui_sessions/session_s09p_r2_tblastn_subject_cache.md)を書いた。その大部分（命令数の 76%）は codon ごとの `GeneticCode::get` の遅さによる。**推奨案**：R2 の最初の作業として翻訳を表引きにし（出力は変えない。TBLASTX・BLASTX も速くなる）、測り直して 20% 未満になればキャッシュは入れない（KISS。DW-8 の条件どおり）。20% 以上のままならキャッシュを移植する。別の案：表引きをせずにキャッシュだけを入れる。指示書は推奨案の順にした。
2. **協調取消の閾値**（計画 §2.3、§10）：**推奨案**：協調取消（NCBI の `TInterruptFnPtr` の移植、計画 §2.3）を入れる条件を「公開する対応規模の最大の入力で、取消の後の再準備（instance の作り直しと Subject の登録し直し）が 1 秒を超えたとき」とする。実測は最大 0.61 s（Firefox、4.6 Mb の Subject と 5.5 Mb の query）、Chromium と WebKit では 0.14 s 以下で、Subject の長さにほぼ比例するので、数 Mb の Subject では条件を満たさない。S17 で対応規模を公開するときに、その最大の入力で測り直す。別の案：閾値を 0.5 s にする（今の Firefox の大腸菌の genome で条件を満たし、協調取消を S17 の前に入れることになる）。

## 合流のときにエンジン側（コーディネーター）が行うこと

1. README の表：S09 を「完了（2026-10-03、ゲート記録）」にし、S09+ の行を「[R2：TBLASTN の subject の前処理](session_s09p_r2_tblastn_subject_cache.md)（条件付き。DW-8 の条件を TBLASTN が満たした。保守者の判断 1 の後）」にする。
2. 計画：§7 の S09 の行を完了にする。§3.2 のコード配置に `infra/engine-worker/`、`infra/reactor/`、`build/reactors.ts` を足す。§5.5 の Auto と作り直しの値（判断 5・6）。§2.3 の Memory の上限の行（判断 7：512 MB）と協調取消の行（実測値）。§10 の「instance を作り直す条件」を決定済みに、「§2.3 の閾値（取消の後の再準備）」を判断待ち 2 に。§8 の Chromium の共有メモリの境界のリスクの行に「起きる（Chromium 149）。guard で対処」。S09+ の行に TBLASTN を書く。
3. `docs/web/verification_cells.tsv`：ブラウザの升目（実行経路「browser, application EngineGateway, threads 1/2/4」、期待値「the native CLI rows; NCBI frozen bytes where certified」、根拠 `web/app/tests/e2e/engine.spec.ts` と `docs/evidence/losat_web_w1/run-20261003T053627Z/records/`）を、BLASTP・TBLASTN・BLASTN の 0/6/7 と TBLASTX の 6 は `checked`、TBLASTX の 0/7 は S08b の升目に合わせて足す（S08b が同じファイルを変えているので、アプリ側では変えなかった）。
4. `docs/web/abi_v2.md`：§5 の「each output stream in chunks of 1 MiB」に stream 2 も含むことを明記する（判断の元になった指摘）。判断 8 の空白だけの BLASTN の入力：`scan` と `register`・CLI の違いをどうするか（`scan` に「レコードが無い」を足すか、今のまま拒否するか）を決める。
5. S08b が TBLASTX の outfmt 0/7 を `checked` にしたら、アプリ側の次のセッションで `web/app/tests/e2e/support/native.ts` の `UNCERTIFIED` から外す（または merge の時に外す）。

## S12 への申し送り

[S12 の指示書](../../losat_web_gui_sessions/session_s12_w3_search_ui.md)の「S09（W1）から引き継ぐこと」に書いた（本物のエンジンのビルド、R1 が効く条件、Auto、診断情報の項目、空白だけの入力、TBLASTX の `-max_target_seqs`、ブラウザごとの試験）。

## 残件と注意

- ThreadHost の障害の経路（時間切れ、thread worker の失敗）は、偽の worker の単体試験だけで確かめた（本物のブラウザで起こす方法が無い）。スレッドのプールを次々に作る検索の不具合は、時間の具合で E2E では再現しないことがあるので、規則は単体試験で固定した（「計測で見つけて直したもの」）。
- serial と threaded の切り替え（Auto が 20,000 バイトの前後の入力を交互に受けたとき）では、毎回 instance と保持した Subject を捨てる（レビューの指摘 L3。時間だけの損失で、上の「取消と再準備」と同じ程度）。
- R1 は、S12 で同じ dataset revision を続く検索に渡すまで、アプリの画面からは効かない（判断 9）。
- Firefox の Wasm の遅さ（Chromium の 4〜6 倍の時間）は、この機械と Playwright の Firefox 151 での値。対応環境の表（S17）と V-MOB（実機）で確かめる。
- WebKit は Playwright の WPE の MiniBrowser で試した。Safari（macOS / iOS）での確認は V-MOB と S17。
- 1 回目のゲートの実行（`run-20261003T043424Z/`）は、計測の段階で TBLASTN の不具合により止まり（`set -e`）、Firefox と WebKit の計測を行っていない。その前に、`run_gate.sh` がエンジンの変数を FakeEngine の段階にも渡していた誤りに気付いて、別の実行を V-ABI と check の後で止めて捨てた（記録は残していない。script を直した）。
- 2 回目のゲートの後に、`policy.ts` のコメントだけを直した（`a1157dd02`）。その前に `npm run check` と両方のビルドの E2E を通した。
- V-PERF の lock（`vperf.lock`）は、このセッションの間に一度も現れなかった。
