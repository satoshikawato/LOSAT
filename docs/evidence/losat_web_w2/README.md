# LOSAT Web W2（Session S10）ゲート記録

- 段階：W2 データ層（[総合計画書](../../losat_web_gui_plan.md) §7 の S10、[指示書](../../losat_web_gui_sessions/session_s10_w2_data_layer.md)）
- ブランチ：`feature/losat-web-gui-app`（アプリ側。worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui-app`、基点は `feature/losat-web-gui` の `9df105107`。計画 DW-7）
- 実行記録：[`run-20260930T081839Z/`](run-20260930T081839Z/)（commit `d28d65ac7` の木。作成後は書き換えない）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)、再現は [`run_gate.sh`](run_gate.sh)
- 判定：**完了条件を満たした**。本物の reactor に依存する 2 つの部分（ABI v2 の `scan_*` によるレコード表と `register` の照合、Engine worker から Data worker への結果のチャンク）は、コーディネーターの指定どおり port と試験用の実装で作り、契約試験を書いた。本物への接続と契約試験の再実行は S09 で行う（[S09 の指示書](../../losat_web_gui_sessions/session_s09_w1_browser_runtime.md)に書き足した）

## コミット

| コミット | 内容 |
|---|---|
| `62fa9ec3c` | アプリ：W2 のデータ層（Data worker、OPFS / Memory の BlockStore、データセット、run の staging、保存領域の所有と回収、使用量の表示、契約試験、E2E） |
| `d28d65ac7` | アプリ：独立レビューの指摘への対応（下の「独立レビュー」） |
| このゲート記録を含むコミット | 文書：このゲート記録、`run_gate.sh`、実行記録、S09 の指示書 |

エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表は変えていない。依存ライブラリも足していない。

## 完了条件と結果

| 完了条件（計画 §7 の S10） | 結果 | 証拠 |
|---|---|---|
| 両方の実装（OPFS / Memory）が契約試験を通る | 通過。`tests/contract/block-store.contract.ts` の 12 case（書き込み、確定（seal）、破棄、読み出しの範囲、確定前の読み出しの拒否、48 MiB の大きなブロック、容量不足など）を、Memory は Vitest（Node）と Chromium の dedicated worker で、OPFS は Chromium の dedicated worker（ディスク上のプロファイル）で実行し、すべて通った | `run-20260930T081839Z/unit-cases.log`、`npm-e2e.log`（`pass [opfs]` 12 行、`pass [memory]` 12 行） |
| 強制終了の後の回収 | 通過。タブの renderer のプロセス（Data worker を含む）に SIGKILL を送ると、そのタブのロックは外れ、`tmp/<session>/` は残る。次に開いたタブがそれだけを消し（`data-removed-sessions="1"`、画面に「Removed the temporary data left by 1 closed tab.」）、自分の検索は通る。通常に閉じたタブも同じ | `npm-e2e.log`（`the next tab removes the data of a tab that is killed / is closed`） |
| 2 つのタブの保護 | 通過。順に開いた 2 つのタブは互いのデータを消さず、1 つを閉じた後に開いた 3 つ目のタブは閉じたタブのデータだけを消す。同時に開いた 2 つのタブは、放置された 1 つのセッションをちょうど 1 回だけ消し、互いのデータを残す。停止（CDP で frozen）したタブのデータとロックは、別のタブの起動の後も残り、再開したタブは結果を読める | `npm-e2e.log`（`two open tabs …`、`two tabs that start at the same time …`、`a stopped (frozen) tab …`） |
| 容量不足の試験 | 通過。CDP でクォータを 0 にすると、次の実行は `failed` になり「Not enough temporary storage for the results of this run: the browser refused to store more data. Earlier results are kept.」を示す。前の実行の結果は読めて、書き出せる（バイトはそのまま）。クォータを戻すと次の実行は完了する。BlockStore と run の出力の契約の容量不足の case も、OPFS と Memory で通った | `npm-e2e.log`（`a run that runs out of storage …`）、`unit-cases.log` |

指示書の作業 1〜7 との対応：1 Data worker と `DataGateway`（Memory の実装を Data worker の中へ移した）、2 OPFS の実装と配置 `tmp/<session-token>/runs/<run-token>/{out0,out6,out7,hits,diagnostics}`（名前は token だけ。E2E で確かめた）と機能を試す経路の選択、3 `File.slice` での読み出し・`RecordScanner` port によるレコード表・レコードごとの SHA-256・DatasetRevision・`register` の照合（port と試験用の実装、契約試験）、4 staging と破棄、worker どうしの MessagePort（port と試験用の実装、契約試験。E2E では別の worker から本物の Data worker へ送った）、5 Web Locks と回収、6 容量不足と使用量の表示、7 Playwright の試験。すべて行った。

## 作ったもの

### 構成（計画 §3.1、§3.3）

```text
Main thread ── Coordinator（application） ── DataGateway（port）── rpc.ts ──► Data worker（src/infra/data-worker）
                    │                                                         ├ DataService：source（File の参照）、レコード表、登録簿
                    │ EngineGateway（port）                                   ├ BlockStore：OPFS | Memory
                    ▼                                                         └ Web Lock：losat-web:tmp:<session-token>
              FakeEngine（S09 で Engine worker に替わる）──MessagePort（openRun が作る）──► run の出力のチャンク
```

| 場所 | 内容 |
|---|---|
| `src/domain/dataset.ts` | `IndexedRecord`（ABI v2 の *scan* の記録そのまま）、`DatasetRecord`（＋SHA-256）、`DatasetRevision`（不変。除外を変えると新しい revision）、`RecordKey`、`recordMismatch`（`register` の照合） |
| `src/ports/data.ts` | `DataGateway` = `DatasetStore`（`addSource`・`indexSource`・`reviseDataset`・`buildRunInput`）＋ `RunStore`（`openRun` は MessagePort を返す、`commitRun`・`discardRun`・読み出し・`deleteRun`）＋ `storageInfo` |
| `src/ports/scan.ts` | `RecordScanner`：ABI v2 の `scan_begin/chunk/end` |
| `src/ports/run-output.ts` | run の出力の経路の型（ABI v2 の stream 0/6/7/1/3、`chunk` と `end`） |
| `src/ports/engine.ts` | `EngineRunRequest` の入力は `{ bytes, records }`、`InputMismatchError`、`run(request, output: MessagePort, onPhase)` |
| `src/infra/data/` | `block-store.ts`（契約）、`memory-block-store.ts`、`opfs-block-store.ts`（同期アクセスハンドル）、`data-service.ts`、`session.ts`（ロック、回収、起動） |
| `src/infra/data-worker/` | Data worker（`data-worker.ts`）、main thread 側（`gateway.ts`）、RPC（`rpc.ts`）、方法の一覧（`methods.ts`） |
| `src/infra/run-output/` | `RunOutputWriter`（エンジン側）、`RunOutputReceiver`（Data worker 側） |
| `src/infra/fake/fake-fasta.ts` | `FakeScanner`（ABI v2 の scan の試験用の実装）。FakeEngine の `register` の照合もこれで読む |
| `src/ui/StorageStatus.vue` | 使用量の表示（操作を妨げない）：保存先、このタブの使用量、ブラウザの見積り、OPFS が使えない理由、回収の状態 |
| `tests/contract/` | 契約試験（BlockStore 12、RecordScanner 23、run の出力 11、エンジンの入力 4 の case）。test runner に依存しないので、Vitest、ブラウザの worker、S09 の本物の実装で同じ case を実行できる |
| `tests/e2e/` | `storage.spec.ts`（7 件）、`contracts.spec.ts`（3 件）。`harness/` は契約試験をブラウザの worker で動かす頁で、E2E の中で Vite がメモリ上にビルドし、Playwright の route で配る（アプリのビルドには入らない）。`support/profile.ts` はディスク上の Chromium のプロファイル |

### 判断（推奨案で進めたもの）

1. **本物の reactor に依存する部分**：`RecordScanner` と `EngineInput.records` / `InputMismatchError` と run の出力の経路を port にし、FakeScanner、FakeEngine、E2E の engine double を試験用の実装にした（コーディネーターの指定）。
2. **`DataGateway` の平らな API**：S01 の `RunStaging`（方法を持つオブジェクト）は worker の境界を越えられないので、`openRun`（MessagePort を返す）・`commitRun`・`discardRun` にした。
3. **run の出力の順序**：出力の経路と RPC は別の MessagePort で、両者の間の順序は保証されない。そこでエンジン側が最後に `end`（送った chunk とバイトの数）を送り、`commitRun` は `end` が届いて数が合うまで待つ。合わなければ失敗にして捨てる。
4. **確定はファイルを動かさない**（計画 §2.3）：登録簿が確定かどうかを持つ。run のブロックは `openRun` で 5 つ作る（追記を同期にするため）。
5. **HSP レコード**：stream 1 の JSON Lines をそのままブロックに保存し、`readHits` が読むときに解く。数は空行でない行の数（BLAST の値は計算しない）。
6. **貼り付けた入力**：`query.fa` / `subject.fa` という名前の File にして、ファイルの入力と同じ経路で扱う（設計書 §6.1）。
7. **run の入力のバイト**：除外の無い revision は元のファイルをそのまま、除外があれば含むレコードの `[header_offset, end_offset)` を元の順に並べる。複数のファイルでは、改行で終わらないファイルの後に改行を 1 つ補う（計画 §5.3）。レコードごとの SHA-256 は `[header_offset, end_offset)` のバイトの値。
8. **所有**：作業ごとの Web Lock だけを所有の印にし、設計書 §7.1 の例の `owner.meta` は作らない。ロックを持てないタブ（Web Locks が無い、または要求が失敗した）は OPFS を使わず Memory にする（他のタブがデータを放置と見なさないため）。ロックの名前 `losat-web:tmp:` と `tmp/` は、旧版のタブを守るため版をまたいで変えない。
9. **回収**：起動時に、ロックが空いている `tmp/*` だけを、そのロックを持ったまま消す。OPFS が開けないとき（放置されたデータでクォータが一杯のときなど）は、先に回収してからもう一度開く。
10. **容量不足**：計画 §5.6 のとおり、その実行を理由付きで失敗にし、途中で Memory に切り替えない。失敗した run のブロックはすぐに消して容量を返す。
11. **レコード表は Data worker のメモリに置く**（設計書 §7.1 の `datasets/<revision>/index.blocks` は作らない）。後から入れる条件：S12 の実測で、レコード表がメモリの負担になったとき。
12. **Memory の BlockStore の上限**：既定は無制限（`setCapacity` は試験で使う）。S09 のメモリの実測で決める。
13. **キューに入れるときに索引を作る**：今の画面は貼り付けだけなので、`enqueue` が入力ごとに source と revision を作る。ファイルの選択、レコード一覧、除外は S12 で `DatasetStore` を使って作る。
14. **BLASTX**：ABI v2 は SX まで BLASTX を拒否する。FakeEngine の `validate` も同じ文（`blastx is not available in LOSAT Web ABI v2 yet`）で拒否するようにした（S01 の FakeEngine は受け付けていた）。画面の program の一覧には残っている（表示の扱いは S12）。
15. **E2E の方法**：容量不足と強制終了は、ディスク上の Chromium のプロファイルで行う（下の「実測で分かったこと」）。強制終了は renderer のプロセスへの SIGKILL（Linux の `/proc` で探す）。Playwright の worker は 2 つまで（DW-7）。

## 試験

- 単体試験（Vitest）：9 ファイル 128 件（S01 の 37 件から増えた）。うち契約試験は、BlockStore（Memory）12、RecordScanner（FakeScanner）23、run の出力（DataService と Memory、同じスレッド）11、エンジンの入力（FakeEngine）4。ほかに DataService 14、セッションの所有 9、RPC 7、データセット 5、coordinator 14。
- E2E（Playwright、Chromium 149、headless shell）：12 件。smoke 2、`storage.spec.ts` 7、`contracts.spec.ts` 3（ブラウザの worker で契約の 35 case：OPFS 12、Memory 12、run の出力 11。run の出力の case は、書き手を別の worker に置き、受け手は本物の Data worker（OPFS）にした）。
- 安定性：`storage.spec.ts` と `contracts.spec.ts` を 3 回ずつ繰り返し、30 件すべてが通った（`e2e-repeat.log`）。開発中にも全 E2E を 3 回と 2 回繰り返して通した。

## 実測で分かったこと

1. Playwright の既定の context（incognito）の OPFS はメモリの上にあり、CDP の `Storage.overrideQuotaForOrigin` を無視する（書き込みは通り続ける）。ディスク上のプロファイル（`launchPersistentContext`）では、同期アクセスハンドルの `write` が `QuotaExceededError`（「No space available for this operation」）になる。
2. CDP のクォータの上書きは、その DevTools の session を detach すると消える。上書きの間は session を保つ。
3. `navigator.storage.estimate()` の usage は、Chromium が書き込みで確かめる usage より遅れる（また、上書きしたクォータを返さない）。容量不足はクォータを 0 にして起こす。
4. この環境（WSL2）では CDP の `Page.crash` が renderer を止めない（`crash` イベントが来ず、ロックも外れない）。renderer への SIGKILL では、`crash` イベントが来てロックが外れる。
5. 同期アクセスハンドルを `close` した直後の `removeEntry` は `NoModificationAllowedError` になることがある。OPFS の実装は、同じディレクトリの削除を順に行い、25 ms ごとに 1 秒まで試し直す。これが無いと、容量不足の後の片付けの誤りが容量不足の理由を隠していた（E2E で見つけて直した）。
6. 同期アクセスハンドルで書いた後の `SyncAccessHandle.write` の戻り値は要求したバイト数と一致した。短い書き込みは容量不足として扱う。

## 独立レビュー

コミット `62fa9ec3c` を、別のエージェントが読み取り専用でレビューした。指摘と対応：

| 重さ | 指摘 | 対応（`d28d65ac7`） |
|---|---|---|
| 中 | 放置されたデータで OPFS が一杯だと、起動の確認が失敗して回収も行われず、以後のタブはずっと Memory になる | 回収を先に行ってから OPFS をもう一度開く。理由の文は「the browser storage of this site is full」にした |
| 中 | BLASTX が FakeScanner の分かりにくい文で拒否される | ABI v2 と同じ文で `validate` が拒否する（判断 14） |
| 中 | ロックの要求が失敗するとデータ層全体が動かない。ロックが無くても OPFS を使っていた | ロックを持てないときは Memory にする（判断 8） |
| 低 | 同じ run の `commitRun` を 2 回呼ぶと run が消えうる | 1 つの commit を共有する |
| 低 | 空の chunk と空行で HSP の数が `readHits` と食い違いうる | 空の chunk を無視し、空行でない行を数える |
| 低 | FakeScanner が見出しの先頭の U+FEFF を落とす。UTF-8 を入力全体で先に確かめるので、誤りの順序が bio と違う | U+FEFF を保ち、bio が読む行だけを行ごとに確かめる。契約試験に 3 case を足した |
| 低 | 除外で空になったファイルの後にも改行を足す | 空の部分の前には足さない |
| 低 | worker の `error` イベントで RPC の client が二度と使えなくなる | 最初の応答の前（起動の失敗）だけを扱う |
| 低 | ディレクトリの作成の容量不足が `StorageFullError` にならない | 変換する |
| 試験 | 遅れて届く chunk を確かめていない。否定の待ち時間が短い。同時に起動するタブを試していない。case の名前の誤り | 使用量が元に戻ることを確かめる。待ち時間を 300 ms にした。同時起動の E2E を足した。名前を直した |

問題が無いとされた点：ロックを取ってからディレクトリを作る順序と回収、確定した run が破棄・取消で消えないこと、出力の経路の順序と数の照合、OPFS のハンドルの扱い、run の入力のバイト、`register` の照合が出力の前であること、走査の契約の期待値が `web/adapter/src/scan.rs` と一致すること、RPC の transfer、層の規則、BLAST の値を計算していないこと。

## 合流のときにエンジン側（コーディネーター）が行うこと

1. 計画の状態の行と §7 の S10 の行を「完了（2026-09-30）」にし、このゲート記録へのリンクを付ける。README の表の S10 を「完了（2026-09-30、ゲート記録）」にする。
2. 計画 §3.2 のコード配置：`infra/memory/` は無くなり、`infra/data/`（BlockStore の OPFS / Memory、DataService、セッション）、`infra/data-worker/`、`infra/run-output/` ができた。`tests/contract/`、`tests/e2e/harness/`、`tests/e2e/support/` も足した。
3. 計画 §7 の S09 の行の完了条件に「S10 の契約試験（`record-scanner`・`engine-input`・`run-output`）が本物の reactor と Engine worker で通る」を足す（S09 の指示書には書いた）。
4. 計画 §5.6 に、判断 4・8・9 を短く足す（確定は登録簿、所有はロックだけでロックの無いタブは Memory、OPFS が開けないときは先に回収）。§2.3 に判断 11・12 の後から入れる条件を足す。
5. S12 の指示書に：入力は `DatasetStore`（`addSource`・`indexSource`・`reviseDataset`・`buildRunInput`）を使い、ファイルを選んだ時点で索引を作る（今は `enqueue` が作る。貼り付けは Data worker の File と snapshot のバイトの 2 つの写しになる）。レコード一覧は `DatasetRevision.records`、除外は `reviseDataset`、Combined は `buildRunInput([...])`。BLASTX の表示（SX までの扱い）。多数のレコードの索引とハッシュの時間の実測。
6. S16 の指示書に：旧版のタブのデータを新しい版が消さないこと（`tmp/` と `losat-web:tmp:` を変えないこと）の E2E（設計書 §7.2）。
7. エンジン側への申し送り（S09 で決める）：BLASTN の入力が空白だけのとき、adapter の `register` は「レコードが無い」として受け付け、CLI は NCBI の警告（query）か誤り（subject）を出す（ABI v2 §4）が、`scan` は `Expected > at record start.` で失敗するので、アプリはキューに入れる前に拒否する。CLI と同じにするなら、`scan` にこの規則を持たせるか、Data worker が BLASTN のときに扱いを合わせる必要がある。推奨案：S09 で adapter の `scan` の扱いに合わせて決め、アプリで BLASTN の規則を作り直さない。

## 保守者の判断待ち

このセッションで新しく生じたものは無い。

## 残件と注意

- S09：`RecordScanner` を serial の reactor で、`register` の照合と出力の経路を Engine worker で本物にし、S10 の契約試験を本物で再実行する（S09 の指示書の「S10 で作った port と、S09 が本物につなぐもの」）。
- Firefox と WebKit では未試験（S09 でブラウザを足す）。`storage.spec.ts` と `contracts.spec.ts` は CDP、ディスク上の Chromium のプロファイル、`/proc` を使うので Chromium（Linux）だけで動く。OPFS が使えないブラウザでの Memory への切り替えは単体試験で確かめた。
- Data worker そのものが落ちた（メモリ不足など）ときの検出は無い。`messageerror` も扱っていない（S09 の worker の監視で扱う）。
- 旧版のタブの保護の E2E は S16（上の 6）。
- Memory の上限、レコード表の置き場所（判断 11・12）は、実測の後に決める。
- 開発中の実験は worktree の外（scratchpad）で行い、一時的な `vite preview` は終了させた。
