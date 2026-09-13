# Wasm optimization follow-up — TBLASTX compilation profile accepted

TBLASTX の cold Node/WASI 比較実行に `--no-liftoff --no-wasm-tier-up` を既定適用した。主対象 n8 は再測定で **6.56→4.44秒、32.31%短縮**、ピーク RSS は **+27.95 MiB**。ユーザー回答「TBLASTX限定でメモリ増を許容する」に基づく採用である。追加の Rust 3候補も独立に実装・計測し、10%短縮基準未達で取り込まなかった。

## Adopted behavior and scope

- `LOSAT/tests/comparison_data.py` が共通 `NODE_ARGS_JSON` と TBLASTX 専用 `NODE_TBLASTX_ARGS_JSON` を別々に解析する。後者の既定値は上記2フラグ。JSON 配列以外、非文字列、NUL は明示的に失敗する。
- `run_comparison.sh` / それを呼ぶ `run_wasm_comparison.sh` が検索ごとにプログラムの Node argv を選ぶ。native / NCBI の引数、BLASTN / BLASTP の既定設定は同じ。`NODE_TBLASTX_ARGS_JSON='[]'` で専用追加分を解除できる。共通フラグは残る。
- manifest の `node_argv_by_program`、NUL 区切り argv、各実行の `ordered_argv` を対応付ける。
- 直接 `node run_losat_wasi*.js ...` を呼ぶ場合は Node フラグを明示する必要がある。Rust / Wasm 自体や JS host の既定動作を書き換えたものではない。
- Node 24.21.0 / V8 `13.6.233.17-node.53` で確認した cold 起動の結果。ブラウザ、全入力、別runtime、再利用済み module/reactor に同じ速度を保証するものではない。

## Measurement and memory policy

全時間はプロセス起動・コンパイル込みの CLOCK_MONOTONIC 秒。各条件で warmup 1 回を除外し、5 回の中央値を比較した。各 repetition で条件順を反転する A/B・B/A 測定で、診断実行の時間は中央値に混ぜていない。RSS は各 measured process の GNU `%M` を byte 換算し、その5回の最大値を条件間で比較する。JS heap や Wasm linear memory 単独の値ではない。

従来の採用基準は主対象10%以上短縮、対照悪化 `max(5%,50ms)` 以下、メモリ増 `max(10%,16 MiB)` 以下。ユーザー承認により **TBLASTX cold Node/WASI の RSS 増分だけ `max(10%,48 MiB)`** へ変更した。約37〜47 MiBという承認対象をカバーする値であり、実行時のメモリ上限や将来の全入力の上限ではない。raw 出力・thread・wall gate は変更していない。[承認と方針](tblastx-policy.json)、[集計と依存資料hash](tblastx-evaluation.json)、[最終採用判定](tblastx-decision.json)。

| 入力 / target | 記録 | 通常 秒 | TurboFan 秒 | 倍率 | 通常 RSS MiB | TurboFan RSS MiB | RSS 増分 MiB |
|---|---|---:|---:|---:|---:|---:|---:|
| MelaMJNV.PemoMJNVA.tlosatx / serial n1 | 今回 | 5.937 | 5.041 | 1.178× | 154.762 | 140.762 | -14.000 |
| MelaMJNV.PemoMJNVA.tlosatx / threaded n1 | 今回 | 6.059 | 5.178 | 1.170× | 161.707 | 172.398 | +10.691 |
| MelaMJNV.PemoMJNVA.tlosatx / threaded n8 | 今回 | 6.563 | 4.443 | 1.477× | 288.012 | 315.961 | +27.949 |
| MjeNMV.MelaMJNV.tlosatx / threaded n8 | 旧記録再評価 | 39.768 | 14.364 | 2.769× | 298.324 | 345.652 | +47.328 |
| AP027280.AP027280.tlosatx / threaded n8 | 旧記録再評価 | 70.532 | 25.572 | 2.758× | 309.949 | 346.801 | +36.852 |

[今回43 process](tblastx-adoption/samples.json) は全件 PASS、除外0。oracle 1、diagnostic 6、cold warmup 6、cold measured 30。serial n1、threaded n1/n8 の全条件で出力バイトが同じ。6 diagnostic の pool 数・tid ごとの spawn→ready→exit(code 0) を検証した。通常の cold 行は軽量に測定しており、個々の worker event の直接証拠は diagnostic 行にある。

今回43件では monotonic / boottime が一致し、最大差は約5.16 µs。realtime 不一致7件は元の方針どおり診断値として保持した。旧2入力の測定には boottime の同境界照合がないため、照合 PASS とは書き換えていない。元の CLOCK_MONOTONIC 観測値として [旧 clock audit](../run-20260913-01/clock-audit.md) の制約を継承する。

旧2入力は fresh 測定ではなく [p1-controls](../run-20260913-01/p1-controls/summary.json) の完全な各5組を、新メモリ基準で再評価した。Wasm artifact、Node executable / V8、Node argv、Rust source、JS execution runners、入力hashは同一。計測器の版は異なり、旧時計の制約も残る。旧 run の全体は別入力 Pemo の TIMEOUT により PARTIAL のまま。全体の失敗・旧判定・raw は上書きしていない。Mje のメモリ増は48 MiBに対して余裕約0.672 MiBしかなく、幅広い余裕があるとの主張はしない。

## Three additional Rust candidates

通常 Node の BLASTP n8、大入力2組をそれぞれ3反復で探索した。候補は同じ baseline から別々に構築し、重ねていない。各候補18 process（oracle 2、diagnostic 4、cold 12）は raw / lifecycle PASS。3反復は探索の根拠であり、5反復の採用証拠とは区別する。

| 候補 | AP027078/AP027131 | AP027132/NZ_CP006932 | 判断 |
|---|---:|---:|---|
| [query 間の subject SEG cache 再利用](cache-decision.json) | 6.426→6.630秒（3.17%悪化） | 9.580→9.648秒（0.70%悪化） | REJECTED |
| [Wasm DP の関数境界を保持](dp-decision.json) | 6.400→6.583秒（2.86%悪化） | 9.229→9.501秒（2.95%悪化） | REJECTED |
| [adjusted matrix の row dispatch 特殊化](matrix-decision.json) | 6.607→6.620秒（0.19%悪化） | 9.779→9.246秒（5.45%短縮） | REJECTED |

Cache は query-parallel branch の folder 内で subject 固定の SEG 結果を再利用した。hit 数は0→18,908 / 32,838に増えたが、検索全体は遅くなった。`map_init` の初期化数を worker 数と同じとは扱っていない。query 依存の near-identical / bias 判定は毎回行う。汎用の新 cache や pruning 変更はない。

DP 境界候補は Wasm に限った `inline(never)` 属性2箇所で、DP の順序や tie-break を変更していない。機械語の最終関数境界や V8 内部の inlining を証明したとは扱わない。Matrix 候補は既存行列の分類を外側で特殊化し、cell の計算順・値を保ったが、最大5.45%短縮だった。全て非採用なので本作業の production Rust 差分は0。

各候補の source patch、release build argv/log/artifact/hash、全 timing/usage/output/worker 記録を `cache-*` / `dp-*` / `matrix-*` に保持する。独立した source audit と matrix 候補の参照コメントの行番号訂正は [audit-summary.md](audit-summary.md) と [matrix-decision.json](matrix-decision.json)。測定済み source を後から修正して同じ成果物の由来と称していない。候補を将来復活させる場合は参照訂正・必要な境界unit・5反復以上の採用検証・native/serial/retry/行列対照を新しい版として行う必要がある。

## Verification and parity boundaries

- Python focused/full harness gate: [39 tests / OK](final-python-tests.log)。新規3件はプログラム別設定、解除・不正値拒否、manifest/NUL argv一致を検証。
- [実比較12件](simple-verification.json): TBLASTX default / optout、BLASTP SicyWSV/CoBV、BLASTN AP027152/AP027202。それぞれ serial / threaded n8 / current NCBI の3出力がバイト一致。実際の Node argv の位置と manifest も一致。BLASTP の既存描画は3行で成功した。これは各1回の呼出し確認で、採用中央値には混ぜない。
- [shell syntax / Python compile](final-syntax-checks.json)。JS host と production Rust に今回の変更はないため、その再build・全Rust test・host全回帰の追加実行は対象外。前段の最終4artifactと同じsource/hashを再利用し、今回の実検索は同じartifactで実行した。
- NCBI source owner: `c++/src/algo/blast/blastinput/cmdline_flags.cpp:46–75` の query / subject / out / num_threads 引数。コンパイル設定は Node host の引数であり、NCBI の検索機能を追加していない。NCBI 2.17.0 executable はテスト oracle のみ。
- TBLASTX current fixture は標準 code 1 / outfmt 6 / E-value既定。main raw SHA-256: `bbef153094025db89dd4a69b2a3c4246b0b8853396fe740bc7e1bd654db105b6`（351,051 bytes）。他の全出力hashは各 summary / sample / result に保存。
- 初回 run の E-value 10/100/10000、gencode 4 DB oracle、長入力・再利用の検証を保持する。gencode 例外は既存の非既定 local-subject 契約だけであり、今回拡張していない。
- P4 は引き続き PARTIAL。既存 megablast の NZ_CP006932 self は未解決、Sakai/MG1655 の既存 SOURCE_UNDETERMINED 差も未合格のまま。native Gate A/B や凍結回帰99/145の残り46件を今回の focused PASS で完了扱いしていない。 [前段 P4](../run-20260913-01/P4.md)。

## Provenance, reproduction, and handoff

HEAD は `c0e123d771f3a58c7529aed4d3f684a1a81f4cb5`、開始時から多数の未コミット変更がある。Rust/Cargo 1.92.0、NCBI source HEAD `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`。HEAD 単独を測定ソースの識別には使わない。

`snapshot.json` / `start-source.tar.gz` は最初の207ファイルを保持。開始 archive に shell / README / 計画が漏れていたため、前段の凍結 `final-task-source.tar.gz` から補った `snapshot-v2.json` / `start-source-v2.tar.gz` を別途作った。補完3ファイルの出典と限界を記録し、元々開始 archive に存在したとラベルし直していない。shell の4つの既知置換を前段版へ適用すると実装版と一致し、独立監査済み。

`current-task.patch` は今回だけの5ファイル差分。production は0、比較設定2ファイル、test1ファイル、文書2ファイルに分けて確認する。`final-task-source.tar.gz` / `final-task-source.json` と `worktree-preservation.json` で前段のユーザー変更との境界を保持する。生成物はこの run directory のみ。旧 run の全 checksum を再照合した結果も保存する。

各 `*-command.json` が正確な呼出し条件。再実行時は新しい出力先を使う。`experiment-drivers.tar.gz` に隔離候補作成・build・測定・集計driverと、その検索前に失敗した初版2件を保存する。候補 source は前段 `baseline-source-v2.tar.gz` と各 `*-source.patch` から復元できる。絶対パスを移す場合は source / input / runtime / artifact hash を照合する。NCBI は測定器の oracle に限り、LOSAT runtime/buildへの依存にしない。

通常コンパイルの比較 baseline を取る場合は両設定を空にする:

```bash
NODE_ARGS_JSON='[]' NODE_TBLASTX_ARGS_JSON='[]' BENCHMARK_PROGRAMS=tblastx \
  BENCHMARK_CASE=MelaMJNV.PemoMJNVA.tlosatx bash LOSAT/tests/run_comparison.sh
```

TurboFan 既定の比較は `NODE_TBLASTX_ARGS_JSON` を設定しない。実行可能な Node/Wasm/NCBI のパスは [tests README](../../../../LOSAT/tests/README.md) の環境変数で指定する。直接 Node 呼出しの例も同READMEにある。詳細A/Bは `tblastx-adoption-command.json` の引数を使い、基準と候補を独立に明示する。

Commit title: `Enable a TBLASTX-only TurboFan profile for WASI comparisons`

Summary: Apply program-specific Node compilation flags, preserve exact argv and output evidence, and document the accepted TBLASTX memory tradeoff plus three rejected Rust experiments. No commit, push, tag, or publication performed.
