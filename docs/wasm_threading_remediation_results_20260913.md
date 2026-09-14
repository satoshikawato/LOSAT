# Wasm threading 修正の実装・検証記録

> 2026-09-14: N 本の専用子 worker + 待機 caller というスレッド数の契約は、ユーザー指定により[親を含む合計 N 本](wasm_total_threads_20260914.md)へ変更する。以下は9月13日時点の履歴。

2026-09-13。[実装計画](wasm_threading_remediation_plan_20260913.md)の F01〜F14 を実装した。ユーザーが追加承認した BLASTP の既存書式差も修正した。実装と完了済み検証を確認し、ユーザーの指示により長時間の固定回帰は後日に延期した（`DEFERRED_BY_USER`）。

## 対象と根拠

開始時 HEAD は `78c7556da6fa18b44447c3a3fa73774a8ee5956f`。既存の変更を含む作業ツリーを対象にし、開始時のソースコピー・SHA-256 と最終ソースの SHA-256 を別に記録した。この報告を含む実装変更をユーザーの指示により commit する。関連する既存の WASI memory guard と比較 runner/catalog は実装の前提として含め、無関係な既存変更は作業ツリーに保持する。package 検証専用の隔離コピーの commit は `6a849993dbe26c7af8da7421720db3fc2a26f4bf` である。

最終検証は Rust/Cargo **1.92.0**、LLVM **21.1.3**、Node **24.21.0** / V8 **13.6.233.17-node.53**、Linux x86_64 / Intel Core i9-14900HX で実施した。開始時の Node 18.19.1 による診断とは区別する。NCBI ソースは `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`、最終比較 oracle は公式 NCBI BLAST+ **2.17.0+** の Linux 配布バイナリを使用した。NCBI は比較・ソース確認専用である。

| oracle | SHA-256 |
|---|---|
| blastn | `33b64bc67d3149cee2459b2f7766b363323df632cf12c099546de00aea9698b5` |
| blastp | `5ce267c04e4988c265357bfbedc64e809545b6fcfae7ff6775266fabbee8ba0e` |
| tblastx | `583e5d60bbd444ac455d20e0956c5aa0aeef675da8daee8204d8f9376ddb8804` |

## 変更と指摘の対応

| 指摘 | 実装と根拠 | 最終検証記録の場所 |
|---|---|---|
| F01 | 全 subject の総塩基数・件数を一度集計。総長を合計後に `/3`。有効 context の統計・word/hit cutoff を共有し、subject ごとの link cutoff を scan 前に更新。NCBI `seqsrc_multiseq.cpp:175–180`、`blast_setup.c:716–847`、`blast_engine.c:1372–1455`、`blast_parameters.c:998–1082`。 | `matrix-repaired-formats/`、`audit-ncbi-values.json`、`lc-thresholds/` |
| F02 | BLASTN の memory/file writer 選択を共有し、native/WASI の並列 return 経路も同じ出力先へ接続。NCBI `blast_format.cpp:68–96,770–832`。 | matrix の CLI file/stdout、reactor の blastn DP/megablast |
| F03 | 3 engine 共通の検索単位 scoped pool。caller 登録を使わず、部分的な spawn 失敗も含め Rust thread を join。NCBI `prelim_stage.cpp:145–188`、rayon-core 1.13.0 `build_scoped`（rayon 1.12.0）。 | reactor 132 呼び出し、失敗後 recovery、worker exited |
| F04 | 小入力しきい値による明示 N の無効化を除去。 | 小入力・旧しきい値 fixture の n1/n2/n4 |
| F05 | BLASTP 固有の暗黙 cap を除去。明示 cap は全 engine で検証して超過を拒否。 | native/CLI/API の cap 検証 |
| F06 | serial Wasm と parallel 無効 native の N > 1 を検索前に明示拒否。 | matrix の unsupported、serial reactor |
| F07 | subject traversal と linking の gate を分離。単一 subject の frame-group linking も指定 pool 内で実行。NCBI `link_hsps.c:553–558,959–982`。 | tblastx-linking、long-single の stage 診断 |
| F08 | 実 pool 数、host の attempt/spawned/ready/exited、stage の仕事数を別々に記録。未測定の稼働率は null。集計器は矛盾を拒否。 | matrix 112 契約、benchmark 35 契約、Python 集計器テスト |
| F09 | Rayon 上限を事前検証し、構築した pool の実数も検証。黙った縮小を認めない。 | N=256 の WASI 拒否、実 pool 数 |
| F10 | 診断用逐次 scan を `LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS` に改名。linking の並列性と分離。 | unit tests、stage 診断、AGENTS/README |
| F11 | direct API の `{err:#}` による原因 chain と失敗時の結果消去。FASTA error を伝播。 | reactor の不正入力・部分 spawn 失敗 |
| F12 | 実 Wasm 用の reusable CI workflow、固定版 toolchain、通常 oracle・既存回帰・異常系 gate を接続。 | `.github/workflows/wasm-threading.yml`、ローカルの matrix/異常系は合格、固定回帰は一部延期 |
| F13 | WASI cdylib のみに実 Rust sysroot の reactor CRT と `_initialize` を指定。main を一度初期化し、child は thread entry へ入る。API wrapper の process destructor を除去。 | `final-crt-proof/proof.json`、最小 thread probe、正規 reactor の反復 |
| F14 | command/reactor と serial/threaded を別 target directory・別配布名へ分離。imports/exports/memory を検査して誤種別を拒否。 | `format-artifacts/`、逆順 build、`validation-package-smoke-final/` |

F01 の開始時 1×3 fixture は、LOSAT 111 行に対し oracle 99 行で、先頭 E-value は `2.59e-165` 対 `7.59e-165` だった。件数や出力後の補正ではなく、統計と cutoff の入力・呼び出し時点を修正した。短い subject は全体統計には含め、NCBI の `<11` 条件に従って検索を skip する。無効な query context の長さ・cutoff・Lambda 参照も NCBI に合わせた。

追加の BLASTP 書式修正は、outfmt 0 のタイトル、68-byte wrapping、空白・改行・列幅、outfmt 7 の query ごとの header/no-hit/footer、custom fields の accession 名を対象とする。NCBI `showdefline.cpp:668–959`、`showalign.cpp:340–359,2385–2471`、`ncbistr.cpp:5088–5340`、`tabular.cpp:415–439,845–859,1111–1217` を参照した。local FASTA ID の version 風の末尾は保持する。`stitle` は NCBI の local FASTA 経路と同じ `N/A` である。bitscore の小数表示の最小幅も `align_format_util.cpp:952,986–993` に合わせた。raw output の比較後に header 除去や正規化は行わない。

## 検証結果

| gate | 結果 |
|---|---|
| 新規 command/oracle matrix | **300 記録 PASS**。直接の candidate/oracle 比較 160 件はすべて raw bytes 一致。1×1/1×3/3×1/3×3、長さの偏り、`/3` の端数、無効 context、形式 0/6/7/custom を含む。 |
| 同一 threaded reactor | **132 呼び出し PASS**。96 成功・36 想定失敗。1→2→4→2→1、engine 切り替え、部分 spawn 失敗後の再検索を確認。 |
| serial/cap API | **48 呼び出し PASS**。拒否時に空の結果、原因、worker 0 を確認。 |
| 回復不能な host 障害 | 起動 timeout、起動 error、trap、異常終了の 4 ケースを subprocess で検証。すべて期待した異常終了となり、timeout を成功扱いしない。 |
| LC738874/LC738875 閾値 | **12 記録 PASS**。E-value 10/100/10000 の公式 oracle、native n1、serial n1、threaded n4。 |
| 凍結 PR5 回帰 | **99/145 ケース成功、46 件はユーザー指示で延期**。native 34/43・serial 32/41・threaded 33/43。retained Linux oracle 0/6・repeatability 0/12 は未完了。初回は 98 PASS・5 TIMEOUT・3 中断を保存。3,600 秒上限の再実行で serial p11 が成功し、残る実行中 7 件をユーザー指示で停止した。計 114 試行（99 成功・5 timeout・10 中断）を保存し、失敗を成功へ変更していない。固定回帰全体は未合格。 |
| Rust tests | **613 passed、3 既存 ignored、0 failed**。`cargo test --all-features`。 |
| Rust lint/format | `cargo clippy --all-features --all-targets -- -D warnings`、`cargo fmt --check` PASS。 |
| host JS | **17 tests PASS**。実 Worker、memory trap、guard の index/境界・共有 memory growth を含む。 |
| Python gate tools | **92 tests PASS**。隔離コピーで集約実行。期待される不合格・不足証拠も拒否できることを検証。 |
| pure Rust runtime boundary | **PASS**。最終ソース検査。build.rs は正規化 SHA を固定した Rust CRT 設定の観測対象であり、汎用 build script 許可ではない。NCBI delegation allowlist は空。 |
| 配布検証 | 両 build 順序で 4 Wasm の SHA が一致。source package 356 files の作成・展開後再ビルド・最終 Rust source bytes 一致。配布 archive 展開後の CLI 9 件と reactor 18 回が PASS。 |

繰り返し stress の最後の 12 回では Wasm linear memory が **7,864,320 bytes** で一定だった。これはこの fixture の観測であり、全入力のヒープ上限や leak 不在の一般的証明ではない。各検索の worker 終了と host event の排出も別に検証した。

独立した `ncbi_parity_auditor` は、production 差分、NCBI 統計の抽出 C 計算、出力の hash、worker events、性能集計を read-only で確認した。formatter の wrapping 計算量・backspace 条件・accession の指摘は修正して再検証した。最終監査でも実装上の blocker は見つからず、固定回帰の成功 99 件の raw hash・argv・worker 契約を独立確認した。未完了 46 件と正式なリリース認証は、この完了判断の対象外である。

## 性能と残る費用

最終計測の各条件は **warmup 1 回、timed 5 回**。実行順を交互化した。cold process は Node 起動、artifact 検査、guard、compile、worker、検索、出力、終了処理を含む。CPU/wall 比から worker 稼働率は推定していない。

`performance-final/` に 448 記録（診断 70、warmup/timed 378）、63 条件の集計を保存した。内訳は出力・実行契約を検証した candidate 35 条件と、出力一致を確認した baseline 28 条件である。baseline の暗黙の pool 縮小は現行の実行契約に合格したとは扱わず、以前の挙動と修正費用の記録として区別する。出力が誤っている baseline の複数 subject 7 条件は測定対象から除外した。壊れた baseline reactor に対する speedup も算出していない。

修正後 cold process の中央値（ms）。全条件の範囲・peak RSS・hash は [cold-performance.tsv](evidence/wasm_threading_20260913/cold-performance.tsv) に記載する。

| fixture | native n1 | n2 | n4 | serial n1 | threaded n1 | n2 | n4 |
|---|---:|---:|---:|---:|---:|---:|---:|
| small | 4.514 | 4.738 | 5.407 | 101.937 | 166.675 | 319.606 | 512.136 |
| multi-subject | 8.242 | 8.302 | 7.976 | 97.779 | 177.663 | 303.885 | 484.998 |
| old-threshold | 3.137 | 3.671 | 4.160 | 109.311 | 180.759 | 314.579 | 467.582 |
| long-single | 13.731 | 13.169 | 12.961 | 127.445 | 204.629 | 360.909 | 527.341 |
| dense-biased | 9.654 | 7.537 | 6.838 | 134.188 | 194.547 | 327.407 | 471.605 |

小入力では worker 起動費用が支配的であり、今回の修正は速度向上を示さない。旧しきい値 megablast の baseline threaded n1/n2/n4 の中央値は 118.605/149.863/118.160 ms、修正後は 180.759/314.579/467.582 ms だった。baseline n4 は暗黙に直列化されていたため、n4 の前後比を同じ並列処理の速度差として説明できない。dense-biased の native n4 は n1 より短かったが、5 回の局所的な観測に留める。

compiled module の再利用は command ごとに新しい instance を作る測定、reactor は同一 instance の測定として分けた。210 呼び出し（175 timed）の出力・worker 契約が一致した。API 用の入力文字列は測定開始前に各 fixture 一度だけ読み、計 6,189,606 bytes を保持した。API invocation 時間はコピー・I/O を含み、engine 単独時間と呼ばない。RSS は process lifetime の high-water であり、個々の検索だけに帰属させない。

再利用準備の費用（ms、各 process で一度）。検査時間は raw module compile を含む。

| mode | 検査 | guard | compile | process peak RSS (MiB) |
|---|---:|---:|---:|---:|
| serial-command-compiled-module | 7.861 | 0.000 | 2.898 | 187.43 |
| threaded-command-compiled-module | 5.796 | 56.907 | 3.975 | 395.14 |
| threaded-reactor-same-instance | 4.821 | 36.841 | 2.446 | 220.34 |

各条件の中央値・範囲・process RSS・worker 起動時間・invocation/終了待ちは [reuse-performance.tsv](evidence/wasm_threading_20260913/reuse-performance.tsv) に保存した。worker 起動時間は host の spawn_attempt→ready の壁時計差であり、検索 worker の稼働時間ではない。

`cost-profile/` は最終 BLASTP 書式修正前のソースを複製して計測コードを加えた診断専用 build である。別の source/artifact hash と instrumentation patch を保存した。24 出力は oracle と一致した。最終 production build の速度試料には混ぜていない。

BLASTN DP の threaded n4 では scratch constructor 18 回の集計 0.790 ms、DP band reserve/check 23,355 回の集計 2.996 ms、traceback row reuse/allocation 13,977 回の集計 13.687 ms、subject 順序復元・hitlist replay 1 回 0.171 ms を観測した。TBLASTX 3×3 threaded n4 の subject sort/flatten は 0.020 ms、選択境界での結果 vector capacity は最大 44,416 bytes だった。

これらは instrumented scope の累積時間と特定境界での vector capacity である。allocator 単独時間、完全な memory peak、worker 稼働率ではない。nested gap scripts、共有文字列、allocator metadata、stack、realloc 中の一時的重複は除外する。未計測の BLASTP 等の counter=0 は費用ゼロを意味しない。

全 subject の統計集計は一度の O(S)。context 統計は subject loop 外で共有し、NCBI の pruning を前倒ししていない。既存の subject-index sort と HSP 保持・replay を維持した。外側の subject-index sort と batch traversal/flatten は O(S log S + H)、外側 vector/HSP record の保持は O(S + H) である。この計算量には既存の hitlist/pruning 処理を含めない。BLASTN の hitlist replay は subject ごとの query slot 走査と query 内 HSP sort を行うため、replay 全体を H に対して線形とは扱わない。保持量にも可変長の nested scripts/strings は含めない。formatter wrapping は newline を再走査し続けない線形実装とした。常駐 pool、新しい candidate pruning、scan 分割、無制限 cache は追加していない。明示的な compiled-module 再利用 factory は呼び出し側が所有する単一 artifact の解析・compile のみを再利用する。

## 最終成果物と再現

| artifact | SHA-256 |
|---|---|
| native Linux command | `834a99ed11f6a4415df6f119d11736a16e50515d1f0f851fda9a3f687688a35f` |
| serial command | `f687af514acb81a62219237c25425412d24c8d5a527d11a988e41dc114c53940` |
| threaded command | `13b88d349ec3449330b11428bb15e9a688448c4e144e756b6f038b4b9412a437` |
| serial reactor | `76c4cc492ac274c716f07815e45b61744d731466b599b61533cd1d65767aac3d` |
| threaded reactor | `f853c69b6ea9070b3d86a9b9be3dd8e9c539c40475a4da84f37089dd4bb8975a` |

Cargo.lock SHA は `7a81afba8adc0fee7efe2868103dd468fe274c7e28629df4dd360725336b26e8`、build.rs のレビュー済み正規化 SHA は `7d70a92bfe68a76662ed485e0d98ea82e4294a84c959714fb62ad331ef808830`。target/features/link flags/CRT と runtime JS の SHA は各 artifact JSON に保存した。4 Wasm は build 順序を逆転しても同じ bytes となった。

再現コマンドは [tests README](../LOSAT/tests/README.md) にまとめた。使用する主要 runner は `build_wasi_artifacts.py`、`check_wasm_threading.py`、`check_wasm_threading_regressions.py`、`check_tblastx_thread_thresholds.py`、`benchmark_wasm_threading.py`。固定回帰の CI は `--jobs 3`、ローカルの継続 runner は独立した 6 subprocess を併走させた。各検索の要求スレッド数は同じであり、この回帰の時間を性能比較には使わない。`--timeout-seconds` は正の値を受け取り、既定値は 3,600 秒。最初の run の 900 秒 watchdog による失敗は、記録された wall_seconds と同一の測定時間とは扱わない。再実行では artifact・凍結 FASTA・Node・host JS の同一性を別途確認する。

証拠の作業場所は `/tmp/losat-wasm-remediation-20260913/`。証拠は [保管先](evidence/wasm_threading_20260913/README.md) にアーカイブ済み。[index](evidence/wasm_threading_20260913/index.json) に 7,598 members・49,960,799 bytes を記録し、全 member の hash を検証した。archive SHA-256 は `284296ec044db0d5ea4cf1be8f731a1b11cec12cf9e360bacfddfe8f5e380039`。`frozen-deferred.json` は成功 99 件の根拠と未完了 46 件の argv・環境・期待 hash を保持する。古い `/tmp/losat-wasm-thread-request/` は上書きせず hash 記録を残した。途中の `frozen-regression/` は最終 artifact へ切り替えるため 97 記録時点で中止した不完全な測定であり、最終結果には含めない。旧 run の p09 付近で native artifact が一時的に差し替わった履歴も保存し、最終 gate は新しい出力ディレクトリで最初から実行した。

## 認証と引き渡しの境界

凍結 PR5 の Gate A と登録済み native NCBI fingerprint の Gate B は変更していない。local-subject TBLASTX 非既定 `db_gencode` の承認済み例外も維持する。retained Linux oracle の一致は、この作業ツリーの回帰証拠であって macOS/Windows の hosted Gate B 認証ではない。

CI workflow は実装したが、このセッションから remote CI を実行していない。既存 RC helper は runtime/build 変更後の旧認証流用を `RC_HANDOFF_FAILED` として正しく拒否した。今回の archive はローカル検証用成果物である。正式リリースには、候補 SHA を固定した新しい認証 lineage と既存 hosted Gate B が必要となる。

英語 commit title: `Fix WASI reactor startup and exact search-scoped threading`。
