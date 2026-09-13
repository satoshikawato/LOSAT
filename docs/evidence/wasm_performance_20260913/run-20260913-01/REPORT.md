# Wasm性能改善計画の実装・測定結果

測定基盤の修正を実装し、P1〜P3の候補を実測して採否を確定した。TurboFan固定と3つのRust局所候補は採用基準を満たさず棄却した。最終の検索実装と4成果物は元baselineと同一であり、製品の既定性能が改善したとは主張しない。P4は既存のmegablast出力差2件により **PARTIAL**。

## TurboFanでどれだけ速くなり、メモリが増えるか

Node24.21.0 / V8 13.6.233.17-node.53、同じthreaded WASI成果物、8threads、通常設定と `--no-liftoff --no-wasm-tier-up` を比較した。warmup1回＋計測5回、baseline/candidate交互、毎回新規プロセス。時間は起動・コンパイル・worker生成・検索・出力・終了待ちを含む中央値。メモリは計測5回の最大process peak RSSであり、Wasm linear memoryそのものではない。表の比較対象はすべて現在のNCBI oracleとraw bytes一致。

| 入力 | 通常 → TurboFan（秒） | 速度比 | peak RSS（MiB） | RSS増加 |
|---|---:|---:|---:|---:|
| MelaMJNV.PemoMJNVA.tlosatx | 6.81 → 3.99 | 1.71× | 293.9 → 309.3 | +15.5 MiB (+5.3%) |
| MjeNMV.MelaMJNV.tlosatx | 39.77 → 14.36 | 2.77× | 298.3 → 345.7 | +47.3 MiB (+15.9%) |
| AP027280.AP027280.tlosatx | 70.53 → 25.57 | 2.76× | 309.9 → 346.8 | +36.9 MiB (+11.9%) |
| WSSV.PajaWSV.losatp | 1.22 → 2.47 | 0.49× | 281.1 → 384.9 | +103.8 MiB (+36.9%) |
| EDL933.Sakai.losatn.megablast | 2.28 → 3.14 | 0.73× | 454.6 → 454.8 | +0.1 MiB (+0.0%) |

TBLASTXの3入力では **1.71〜2.77倍、約15〜47 MiB増**。主入力Mela/PemoAは計画のメモリ基準内だが、Mje/Melaの+15.9%、AP027280 selfの+11.9%は各入力の `max(10%,16 MiB)` を超える。汎用採用では小さなBLASTPとmegablastの時間悪化も基準を超えるため、全体の既定設定にもTBLASTX共通の既定設定にも採用しない。これは利用不能という意味ではなく、今回の工学的な採用上限による判断である。

全サンプル・範囲は `p1-adoption/summary.json` と `p1-controls/summary.json`、表の計算値は `turbofan-n8-comparison.json`。冷起動以外の探索は `p1-reuse-summary.json` に分離した。長入力の単発正しさ診断を5回中央値へ混ぜていない。ブラウザ性能への外挿もしない。

## 統合した変更

- 新規runディレクトリ、成功status、対応argv、出力/入力/runner/toolのhashと版を要求し、古い出力・失敗・timeout・片側だけの成功が比較表へ混入する経路を閉じた。
- Node引数を配列で記録し、既存benchmarkのfixture・実行寿命・反復・timeout選択と失敗記録を拡張した。再利用測定の検索境界と証拠処理、検索ごとの値とprocess lifetime RSSを区別した。
- Python watchdogの早期起床時にmonotonic期限を再確認する。realtime/boottime/GNU elapsedも診断として保存する。後退するrealtimeによる初回起床遅延の限界は `clock-audit.md` に記載した。
- 比較・計測の回帰テストと利用手順を更新した。NCBI実行はテストoracleに限り、検索実装やbuildへの依存は追加していない。

## 段階別判断

| 段階 | 結果 | 理由・証拠 |
|---|---|---|
| [P0](P0.md) | COMPLETE / 計測変更を統合 | Python36件、実Wasm再利用4シナリオPASS。旧trapは現在のbaselineで非再現 |
| [P1](P1.md) | COMPLETE / 設定・局所候補REJECTED | TurboFanはメモリ/他task基準超過。小関数呼出し候補は1.78%短縮のみ |
| [P2](P2.md) | COMPLETE / batch32 REJECTED | 大入力8.03%短縮、他方18.72%悪化。元のbatch16を維持 |
| [P3](P3.md) | COMPLETE / traceback行書込み候補REJECTED | 主入力4.39%悪化/2.46%短縮。診断はDPに処理集中を示すが、十分な短縮なし |
| [P4](P4.md) | PARTIAL |4成果物のbaseline同一性、TBLASTX閾値12件と実比較/描画5件PASS。megablast2入力はraw差が残る |

候補が探索で不採用となった後の拡大採用試験は未実行であり、合格扱いしていない。P2の物理copy/allocation、P3の正確な命令別費用やretry経路など、取得しなかった量は未測定。worker profileは全9isolateを含むが、未解決PCも多く、サンプル比をCPU占有率と呼ばない。

## 残る出力差

NZ_CP006932 self megablastは454行中1行、Sakai/MG1655 megablastは6,476行中5行でidentity/length/mismatch/gap表示に差がある。元baselineと最終成果物、Native/serial/threadedのLOSAT出力は各入力内で同一。今回の回帰ではないが、現在のraw gateでは不合格のまま保存した。

Sakaiは既存の狭いSOURCE_UNDETERMINED登録と同じ差。NZはその登録範囲外で、greedy traceback・gap reduction・endpoint選択のどこで初めて分かれるかは未確定。独立監査の根拠と次の診断境界は [audit-p4.md](audit-p4.md)。新しい例外、独自tie-break、goldenの差替えは行っていない。正式GateA/B認証、旧99/145からの46件解消、全体release parityは主張しない。

## 再実行と証拠

通常の比較手順は `LOSAT/tests/README.md`。`NODE_BIN` には実行ファイルだけを渡す。実験設定の配列は次のとおり。

```bash
export NODE_ARGS_JSON='["--no-liftoff", "--no-wasm-tier-up"]'
```

冷起動の詳細A/Bは保存した `p1-adoption-command.json` / `p1-controls-command.json` と既存 `LOSAT/tests/benchmark_wasm_threading.py` を使い、必ず新しい `--output-dir` を指定する。各manifestにargv・入力hash・source/artifact・Node/V8・oracle・環境がある。元の絶対パスを移す場合は入力と成果物のhashを照合してからパスを差し替える。実験コードの復元は `experimental-source-deltas.json` が指定するbaseline-v2＋delta archiveを使う。診断用driverは `experiment-drivers.tar.gz`、全workerのraw profileは `p1-raw-profiles.tar.gz` と `p3-raw-profiles.tar.gz` に保存する。

`SHA256SUMS` はこのrunの最終ファイル一覧とhash。初期dirty patch archive末尾の33バイトの非patchデータも記録し、原本を維持した。失敗・timeout・欠落build.rsによる無効な初期候補・初期profile失敗も保持した。完全なbaselineはv2で、元archiveを上書きしていない。`current-task.patch` と `worktree-preservation.json` はユーザーの初期dirty変更と今回の差分を区別する。

## Commit handoff

Title: `Harden Wasm benchmark provenance and document optimization decisions`

Summary: Reject stale and incomplete comparison evidence, make Wasm timing and reuse boundaries explicit, and preserve reproducible optimization experiments. Keep rejected performance candidates out of production and report existing megablast parity failures without broadening exceptions.

ローカルcommit、push、tag、publishは実行していない。
