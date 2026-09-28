# LOSAT Web E1a（Session S02）ゲート記録

- 段階：E1a 核の入口：共通部と BLASTP（[総合計画書](../../losat_web_gui_plan.md) §7 の S02、[指示書](../../losat_web_gui_sessions/session_s02_e1a_core_entry_blastp.md)）
- ブランチ：`feature/losat-web-gui`。変更前は `bdeaa9bc6`。エンジンと試験の変更はコミット `afd57e7cb`（この記録の直前のコミット）
- 実行記録：[`run-20260928T185126Z/`](run-20260928T185126Z/)、基準 [`baseline/`](baseline/)、変更後 [`after/`](after/)。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)（`sha256sum --check evidence.sha256` をこのディレクトリで実行する）
- 判定：**完了条件を満たした**（独立監査の結論は下の「独立監査」の節）

## 変更の内容

| ファイル | 内容 |
|---|---|
| `LOSAT/src/api/local_blast.rs` | 出力の配布の共通の型：`OutputSink`（stdout / `-out` のファイル / 呼び出し側の writer）、`FormatOutput`、`HspIndex`、`FormatObserver`、`FormatProbe`、`ReportOutputs`（形式ごとの出力、`diagnostics`、`hits`、`observer`。どれも `Send`） |
| `LOSAT/src/algorithm/blastp/blast_engine.rs` | `run_local`（どのターゲットでもコンパイルされる BLASTP の入口）と `resolve_for_formats`（形式ごとに検索オプションを解決し、形式間で食い違えばエラー）。CLI の `run` と v1 の `run_web_pair_records` は、どちらも `run_local` を呼ぶ。最後の出力の分岐は、要求された形式ごとに回り、`PairwiseHit` は必要なときに 1 回だけ作る。`open_output_writer` と `run_internal_with_records` は不要になったので消した |
| `LOSAT/src/report/pairwise.rs` | `write_blastp_pairwise_report`（BLASTP だけが使う）に観測者を足した |
| `LOSAT/tests/run_local_blastp.rs`、`LOSAT/tests/run_local_support/mod.rs` | 試験 3 件と、各 program の `run_local` の試験が共有する道具（形式ごとの出力先、HSP レコードの受け取り、観測者の範囲の記録と検査） |

HSP の ID は、最終の HSP 一覧（`collect_hits_from_hit_lists` の順）の中での番号 `HspIndex` で、outfmt 0・6・7 のどれでも同じ HSP を指す。観測者が報告する範囲は、outfmt 6・7 では 1 行、outfmt 0 ではスコアの行とアラインメントである。outfmt 0 の範囲は、NCBI の `x_DisplayAlnvecInfo` 1 回分の出力から、各 subject の最初の HSP の前に書かれる subject の見出し（`x_DisplayAlnvecInfo` から呼ばれる `x_ShowAlnvecInfo` の中、`showalign.cpp:3613-3632`。呼出しは `3956-3975`）を除いたものに当たる。観測者があるときだけ、HSP の区切りで formatter のバッファを流す。CLI と v1 には観測者が無いので、バッファの動作は変わらない。

CLI の `run` は、オプションの解決と対応範囲の検査を入力の読み込みの前に行い（今までと同じ順序）、その後で `run_local` を呼ぶ。`run_local` は同じ解決と検査をもう一度行うが、結果は同じである。各形式の `-outfmt` は、今までと同じく、検索のスレッドプールの中で検索の前に解析する（CLI は、コマンド行を解析するときにすでに検査している）。

### v1 の振る舞い

v1（`run_web_pair_records`）の振る舞いは変わらない。二つの誤りを含む要求で先に報告される誤りも同じである。変更前と変更後の threaded reactor に同じ要求 14 件（誤りの組合せ 5 件、単独の誤り 4 件、成功 5 件）を送り、状態・結果の長さと SHA-256・エラーの文がすべて一致した（`v1_requests.js`、`run-…/v1-requests-{before,after}.jsonl`）。`check_wasi_reactor.js` の 344 件の出力とエラーの記録も一致した（`run-…/v1-reactor-records-compare.txt`）。serial reactor も、`check_wasi_api_limits.js` の 25 件の記録が変更前と一致した（`run-…/v1-serial-reactor-records-compare.txt`）。変更前の reactor は、`bdeaa9bc6` の木（`git archive`）から作り、その command モジュールのハッシュは以前に作ったものと一致した。

reactor の線形メモリの大きさ（`runs.json` の `memory_bytes`、172 回の読み取り）は、変更後に 169 回で大きく、最大で 3 ページ（192 KiB）増えた。小さくなったのは 2 回（最大 4 ページ）、同じだったのは 1 回である。最大値は 11,862,016 から 11,993,088 バイト（+1.1%）になった。出力には関係しない。

## 完了条件と結果

| 完了条件（計画 §7 の S02） | 結果 | 証拠 |
|---|---|---|
| 変更したコードを使う全 program の既存のゲートの出力が変わらない | 変更前と変更後で 236 件すべて、コマンド・終了コード・出力の SHA-256・stderr の SHA-256 が一致（差 0）。対象：BLASTN 14、BLASTP 28（manifest の 9 件 × outfmt 6/0/7 と、カスタム列 1 件）、TBLASTX 20、TBLASTN 162（Stage G の matrix）、BLASTX 12 | `baseline/hashes.tsv`、`after/hashes.tsv`、`run-…/capture-compare.txt`、スクリプト `capture_outputs.py` |
| Gate A などの凍結ハッシュが変わらない | 変更前も変更後も同じ：TBLASTN 162/162（Stage G）、TBLASTX 20/20、BLASTP 9/9、BLASTN 13/14（Gate A）。一致しない 1 件 `Sakai.MG1655.megablast` は、変更前の `main` でも同じハッシュで一致しない既知の差である（下の「既存の差」） | 同上 |
| 変更前の基準の出どころ | 変更前のネイティブ（SHA-256 `bda9b7dd…`）は、`bdeaa9bc6` の木から作り直しても同じハッシュになった。変更前の WASI の 4 つの成果物も、同じ木から作った | `run-…/binaries.txt` |
| v1 の serial / threaded reactor の検査が通る | 通過。変更後の 4 つの WASI 成果物で `LOSAT/tests/check_wasm_threading.py` を serial の互換の検査も含めて実行：409 件の command / oracle の記録、reactor の寿命の検査 PASS、形式の失敗 0。v1 の振る舞いが変わらないことは上の節 | `run-…/wasm-threading.log`、`wasm-threading-runs.json`、`wasm-threading-metadata.json` |
| wasm32 だけでコンパイルされる v1 の試験 | 通過（`web_api::tests` の 4 件） | `run-…/wasm32-web-api-tests.log` |
| V-NAT：1 回の検索で複数の形式を出しても、その形式だけの CLI の実行と一致する | 通過。outfmt 0・6・7・カスタム列を同時に要求した `run_local` の出力が、それぞれ CLI（`LOSAT blastp … -outfmt N -out …`）の出力と一致し、`diagnostics` は CLI の stderr と一致した。既定と `-max_target_seqs 1`（8 件 × 8 件の fixture で query と subject の組が 24 から 7 に減る）の両方で確かめた。観測者の範囲は、outfmt 6・7・カスタム列の 1 行ずつ、outfmt 0 の ` Score =` で始まる節に正確に一致し、同じ番号の HSP は、どの形式でも同じ座標と raw score を持つ。未対応の形式は、何も書かずに拒否される | `LOSAT/tests/run_local_blastp.rs`（3 件） |
| V-PERF の非退行（中央値が基準の +5% 以内） | 調査の後、退行は無いと判断した（下の「性能の計測」）。**保守者の確認を求める**：計画の手順の計測で +5% を超えた回がある（まとめて測った回の outfmt 0 の 4 件と、交互に測った 1 回目の serial WASI の outfmt 0 の ×1.080。`perf-check-1.txt` の FAIL の行）。BLASTP（SicyWSV.faa × PajaWSV.faa、`-max_hsps 1`）の outfmt 6 と outfmt 0、ネイティブ・serial WASI・threaded WASI（n=4）。出力のハッシュはすべて同じ | `run-…/perf-{1,2}.json`、`perf-check-{1,2}.txt`、`run-…/perf-investigation/`、スクリプト `measure_perf.py` |
| `cargo fmt --check`、`clippy -D warnings`、`cargo test --all-features` | 通過。clippy の終了コード 0。試験は合計 786 件が通過（ライブラリ 589、CLI と結合 176、`run_local` 3 ほか）。`--all-features` は診断用の `blastx-worker-probe` を含むので、CI（`.github/workflows/ci.yml:26-29`）と同じく `LOSAT_BLASTX_WORKER_LOG` を設定して実行した | `run-…/fmt.log`、`clippy.log`、`cargo-test.log`（先頭の行が実行したコマンド） |
| `docs/web/verification_cells.tsv` | 作成。BLASTP の行を埋め、ほかの program の行を用意した | [`../../web/verification_cells.tsv`](../../web/verification_cells.tsv) |

成果物のハッシュと道具の版は `run-…/binaries.txt` にある。`after/binary.json` と `wasm-threading-metadata.json` に書かれた成果物のパス（`native-after/`、`wasi-after-artifacts/`）は、その後の S03 のビルドで上書きされた。S02 の成果物は `s02-final/` に同じハッシュで残してある（`binaries.txt`）。`run-…` は `run-20260928T185126Z` である。

### 性能の計測

計画の手順（1 回の暖機と 3 回の計測。`AGENTS.md`）で測った結果と、超えた場合の調査である。

1. 初めに、変更前と変更後をそれぞれまとめて測り、それを 2 回繰り返した（変更前、変更後、変更前、変更後。`perf-investigation/perf-{before,after}-{1,2}.json`、`perf-compare-{1,2}.txt`、そのときのスクリプト `measure_perf_blocks.py`）。outfmt 6 は 6 件すべて +5% 以内だったが、outfmt 0 は 6 件中 4 件が超えた（×1.107、×1.099、×1.123、×1.260）。一方で、同じ成果物の中央値が回によって大きく動いた（変更前のネイティブの outfmt 6 は 0.536 秒と 0.443 秒、変更後の threaded の outfmt 0 は 0.498 秒と 0.635 秒）。outfmt 0 の経路で変えたのは、HSP の番号の受け渡しと、観測者が無いときに何もしない分岐だけで、計算の量は変わらない。
2. そこで、3 回の計測では結論が出ないと判断し（`AGENTS.md` の「demonstrably inconclusive」）、変更前と変更後を 1 回ずつ交互に測るように `measure_perf.py` を直した。暖機はそれぞれ 1 回、計測はそれぞれ 3 回のままである。これで 2 回測った（`perf-1.json`、`perf-2.json`、判定は `perf-check-{1,2}.txt`）。1 回目は 6 件中 5 件が +5% 以内で、serial WASI の outfmt 0 が ×1.080 だった（変更前 0.657・0.756・0.632 秒、変更後 0.709・0.923・0.621 秒で、同じ成果物の 3 回の中にも 49% の開きがある）。2 回目は 6 件すべてが ×0.967〜×1.030 だった。
3. 調査として、同じ交互の方法で計測を 12 回ずつにした（`perf-investigation/perf-investigation.json`、`perf-investigation-check.txt`）。6 件すべてが ×0.989〜×1.008 で、serial WASI の outfmt 0 は変更前 0.620〜0.642 秒、変更後 0.610〜0.638 秒（中央値 ×0.994）だった。

結論：性能の退行は無い。まとめて測った回の超過と、交互に測った 1 回目の超過は、機械の状態の変化と 1 回ごとのばらつき（約 0.5 秒の実行）によると判断した。この原因は、同じ成果物の中央値の変動と、交互に測ったときの一致から推し量ったもので、CPU 時間や負荷は記録していない。これからのセッションは、交互に測る `measure_perf.py` を使う。

## 既存の差（この変更とは無関係。記録だけ）

- `Sakai.MG1655.megablast` は、変更前の `main` でも Gate A の凍結ハッシュと一致しない（変更前後とも `ac817776…`）。LOSATX の Stage G の権威の記録（`docs/evidence/losatx_stage_g_authority_v3/run-20260928T121516Z/GATE_STATUS.md` の 3 行目）が、PR5 の Gate A の Sakai の段階を「historical HARD_FAIL」とし、新しい Sakai の期待値を別の版として登録している。
- LOSAT の BLASTP は、NCBI が hitlist の大きさが 5 未満のときに出す警告「Examining 5 or more matches is recommended」（`blast_args.cpp:2975-2976`）を出さない（BLASTX は出す）。`-max_target_seqs 1` の NCBI 2.17.0 の実行でこの警告が出ることを、独立監査が確かめた。試験は `diagnostics` を CLI の stderr と比べるので、この警告を移植したときも書き直さずに済む。
- LOSAT は `-out` のファイルを検索の後に作る。NCBI は引数を処理するときに開く（`blast_args.cpp:3478-3480`）。この変更はその時点を変えていない。
- outfmt 0 の既定の表示数（NCBI の 500 件の説明と 250 件のアラインメント、`format_flags.cpp:219-221`）を LOSAT の BLASTP が適用しているかは、まだ確かめていない。適用していない場合、v2 の `out0` の範囲を持つ HSP の集合が NCBI と違う。

## 独立監査

`ncbi_parity_auditor` の定義に従う読み取り専用の監査を受けた。

- 1 回目：出力のバイトが変わるような欠陥は無く、バイトの同一性は supported とした。そのうえで、主張全体は unsupported とした。理由は、変更したブロックの一部に NCBI の参照が無いこと（AGENTS.md の規則 4）、`OutputSink` のコメントが NCBI の `-out` を開く時点を誤って書いていたこと、CLI の `run` が `run_local` を通っていなかったこと（計画 TD-2、PD 決定項目 4.1）である。証拠の不足として、CLI と `run_local` の直接の比較、`-max_target_seqs` の境界と形式をまたぐ中身の確かめ、変更後の wasm32 の試験、ゲート記録を挙げた。
  - 対応：参照を足してコメントを直し、CLI を `run_local` に通し、試験を足した。
- 2 回目：主張（出力のバイトが変わらないこと、複数形式の出力が単一形式の CLI と一致すること）は supported とした。独立に作り直した成果物のハッシュの一致、CLI の 15 件の誤りと境界の場合の変更前後の一致、reactor の 344 件の記録の一致を確かめたうえでの結論である。残った指摘と対応は次のとおり。
  - NCBI の行番号の誤り（`blast_args.cpp:3476-3478` は `3478-3480`）：直した。
  - outfmt 0 の観測の単位の説明（`x_DisplayAlnvecInfo` 1 回分と書いていたが、subject の見出しを含まない）：コメント、`abi_v2.md`、計画 §4.5 を直した。見出しの範囲の扱いは S05 の指示書に書いた。
  - `run` の事前の検査と `needs_pairwise_hits` に NCBI の参照が無い：足した。
  - v1 で `-outfmt` をオプションの解決と対応範囲の検査より先に検査するようになり、二つの誤りを含む要求で先に報告される誤りが変わっていた：`run_local` での `-outfmt` の事前の検査をやめ、変更前と同じ順序に戻した（上の「v1 の振る舞い」）。
  - 証拠の不足（`evidence.sha256` が無い、成果物の出どころのコミットが無い、v1 の変更前との比較、fmt と clippy の記録、性能の計測を交互に行った記録、2 回目の監査の結論が空欄）：すべて足した。性能の JSON は開始と終了の時刻を持つ。
  - 追加の観察：試験が NCBI との既知の差（上の警告）を固定していた。`diagnostics` を CLI の stderr と比べる形に直した。性能の計測が outfmt 6 だけだった。outfmt 0 の計測（`blastp-fmt0`）を足した。足した計測で +5% を超えた場合の調査と、計測の方法の変更は、上の「性能の計測」に書いた。
- 3 回目（直した後の確認）：supported。独立に作り直したネイティブと threaded reactor のハッシュが記録と一致し、変更前・作り直した変更後・記録の変更後の 3 つの reactor で `v1_requests.js` の出力が一致した（2 回目の監査のときの reactor では、誤りを二つ含む最初の 4 件が違い、スクリプトが以前の欠陥を検出できることも確かめた）。CLI の 10 件の誤りの場合も変更前と一致した。`evidence.sha256` の検査、fmt、clippy、`run_local_blastp` の試験も通った。性能の扱いは `AGENTS.md` と V-PERF に沿うとした。残った指摘と対応は次のとおり。
  - `local_blast.rs` の `showalign.cpp:3613-3630` は、引用した `out << "\n";` が 3632 行目にあるので `3613-3632` が正しく、呼出し元の `3956-3975` も引くべき：S03 のコミットで直す（この記録と `abi_v2.md`、S05 の指示書は直した）。
  - serial reactor の変更前との比較が無い：足した（上の「v1 の振る舞い」）。
  - 記録にある成果物のパスが上書きされた：上に書いた。
  - V-PERF の表で FAIL の行を明示し、保守者の確認を求めるべき：表に書いた。
  - 追加の観察：`diagnostics` を LOSAT の CLI とだけでなく NCBI の stderr とも比べること、outfmt 0 の節の検査で座標とアラインメントの中身も見ること。これは次のセッション以降の課題とした（TBLASTN の試験は NCBI の保存した stdout と stderr と比べている）。

## 実行の方法（再現）

```bash
# 変更の前後で同じ 236 件を実行し、比べる（入力は凍結時と同じパス文字列で渡す）
python3 docs/evidence/losat_web_e1a/capture_outputs.py run --losat <LOSAT> --out <dir> --jobs 6
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare <before>/hashes.tsv <after>/hashes.tsv
# 性能の非退行（変更前と変更後を 1 回ずつ交互に）
python3 docs/evidence/losat_web_e1a/measure_perf.py run --before <LOSAT>,<serial.wasm>,<threaded.wasm> --after <LOSAT>,<serial.wasm>,<threaded.wasm> --cases blastp,blastp-fmt0 --out <file.json>
python3 docs/evidence/losat_web_e1a/measure_perf.py check <file.json>
# v1 の WASI の検査（NCBI BLAST+ 2.17.0 を比較用の oracle として使う）
RUSTUP_TOOLCHAIN=1.92.0 python3 LOSAT/tests/build_wasi_artifacts.py --target-dir <dir> --output-dir <artifacts> --include-serial
RUSTUP_TOOLCHAIN=1.92.0 python3 LOSAT/tests/check_wasm_threading.py --native ... --native-serial ... --serial ... --threaded ... --reactor ... --serial-reactor ... --oracle-dir <ncbi-bin> --output-dir <dir>
# v1 の変更前との比較（変更前の reactor は変更前の木から作る）
node docs/evidence/losat_web_e1a/v1_requests.js <reactor.wasm> <dir>/fixtures/aa3.fasta
node LOSAT/tests/check_wasi_reactor.js <reactor.wasm> <dir>/fixtures <records>
# wasm32 だけでコンパイルされる v1 の試験
CARGO_TARGET_WASM32_WASIP1_RUNNER="node web/tools/wasi-test-runner.mjs" cargo +1.92.0 test --lib --target wasm32-wasip1 --no-default-features -- web_api::tests
```

Gate A の入力は `certify_platform_native_v010.py` の `HISTORICAL_LEXICAL_ROOT`（`/tmp/losat-pr5-runtime-cert-5845d22/LOSAT`）を通して渡す。このセッションでは、LOSATX の再認証の作業が作った同じ名前のシンボリックリンク（このリポジトリの元のクローンの `docs/evidence/losatx_stage_g_recertification/run-20260927T152852Z/pr5_fixtures` への参照）を、WSL の再起動で消えた後に同じ参照先で作り直して使った。`capture_outputs.py` は、実行の前に、そこから読むすべての入力がこの checkout のファイルとバイト単位で同じであることを確かめる。

## 引き継ぎ（S03 以降へ）

- `run_local` の形と使い方は計画 §4.2 に書いた。program をまたぐ振り分けは S05 で足す。
- 基準は `baseline/hashes.tsv` を S03・S04 でも使う（取り直さない）。変更後の比較は `capture_outputs.py compare` で行う。
- 各 program の `run_local` の試験は `LOSAT/tests/run_local_support/mod.rs` を使う（`run_formats`、`assert_observer_ranges`）。
- BLASTX が使う共有の関数（`LOSAT/src/report/pairwise.rs` の `write_hsp_info`、`write_alignment` など）の本体と引数は変えていない。
- HSP レコードを求める（`hits`）とアラインメントを描画する。表形式だけを求める使い方では、描画できない HSP があると CLI の outfmt 6 より先に失敗しうる（S05 の指示書に書いた）。
- この環境では、`rustc` を直接呼ぶスクリプトに `RUSTUP_TOOLCHAIN=1.92.0` が要る（既定の stable は不完全）。`cargo test --all-features` には `LOSAT_BLASTX_WORKER_LOG` が要る。
