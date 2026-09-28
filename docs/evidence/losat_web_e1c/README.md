# LOSAT Web E1c（Session S04）ゲート記録

- 段階：E1c 核の入口：BLASTN と TBLASTX（[総合計画書](../../losat_web_gui_plan.md) §7 の S04、[指示書](../../losat_web_gui_sessions/session_s04_e1c_core_entry_blastn_tblastx.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S03 の後の `a0ebbb849`（エンジンは `00fdfca14`）。エンジンと試験の変更はコミット `458cbe110`（この記録の直前のコミット）
- 実行記録：[`run-20260928T203929Z/`](run-20260928T203929Z/)、変更後の出力のハッシュ [`after/`](after/)。基準は S02 の [`../losat_web_e1a/baseline/`](../losat_web_e1a/baseline/)（取り直していない）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)（`sha256sum --check evidence.sha256` をこのディレクトリで実行する）
- 判定：**完了条件を満たした**。ただし、CLI の振る舞いの差が 2 つある（下の「CLI の振る舞いの差」。出力のバイトは変わらない）。V-PERF とあわせて保守者の確認を求める

## 変更の内容

| ファイル | 内容 |
|---|---|
| `LOSAT/src/algorithm/blastn/blast_engine/run.rs` | `run_local`（どのターゲットでもコンパイルされる BLASTN の入口）。CLI の `run` は、入力を読む前の検査（スレッド数と `-limit_lookup` の検査。`check_blastn_lookup_options` にまとめた）を行い、query を読み、query が空なら以前と同じく subject を読まずに終わり、subject を読んで `run_local` を呼ぶ。v1 の `run_web_pair` も `run_local` を呼ぶ。`post_process_hits_and_write` は要求された形式ごとに書く。以前の `BlastnInMemoryRun`・`BlastnOutputTarget` は不要になったので消した。検索は query の記録の写しを持つ（`prepare_sequence_data` が所有するため） |
| `LOSAT/src/algorithm/blastn/hsp.rs` | 表形式の各行に観測者の印を付けた。ファイルを開く `write_output_blastn_hitlists` は使われなくなったので消した |
| `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs` | `run_local`。CLI の `run` は、スレッド数を検査し、query と subject を読み（以前と同じ順。その間に失敗しうる処理は無い）、`run_local` を呼ぶ。出力を書く箇所を 1 つ（`write_tblastx_outputs`）にまとめた。ネイティブの並列の検索の集約のスレッドは、書く代わりに最終の HSP を返し、`wasm-threads` だけの集約の経路と、順次の経路と同じく、書くのは 1 か所になった。NCBI の順序に並べてから書く、という順序は変えていない |
| `LOSAT/src/common.rs` | TBLASTX だけが使う `write_output_ncbi_order_evalue_hsp_order_to_writer` とその実装に観測者を足した。ファイルを開く `write_output_ncbi_order_evalue_hsp_order` は使われなくなったので消した |
| `LOSAT/src/api/local_blast.rs`、各 `mod.rs` | `run_local_blastn`、`run_local_tblastx` の再 export |
| `LOSAT/src/algorithm/tblastn/args.rs`、`stage_e_report.rs`、`blastp/blast_engine.rs`、`local_blast.rs` | S03 の独立監査の指摘：`SearchSettings` に NCBI の参照を足し、省略していた NCBI の引用に省略の印（`...`）と、区切り文字の行を足した（コメントだけ） |
| `LOSAT/tests/run_local_blastn.rs`、`run_local_tblastx.rs` | 試験 3 件ずつ |
| `LOSAT/tests/run_local_support/mod.rs` | 表形式だけの観測者の範囲の検査（`assert_tabular_ranges`）、一時の FASTA（`TempFasta`）。S03 の独立監査の指摘により、outfmt 0 の節の検査に、Query と Sbjct の座標と、節がその HSP の query の区切りの中にあることを足した（BLASTP と TBLASTN の試験が使う） |

両 program とも、`hits` は S07・S08 で `PairwiseHit` を作るまで呼ばない。観測者は outfmt 6/7 の行について知らせる。HSP の番号は、書いた行の順（最終の HSP 一覧の順）である。

## 完了条件と結果

| 完了条件（計画 §7 の S04） | 結果 | 証拠 |
|---|---|---|
| 既存のゲートと Gate A のハッシュが変わらない | 変更後の 236 件（全 program）が S02 の基準と一致（差 0）。凍結ハッシュ：BLASTN 13/14、TBLASTX 20/20、BLASTP 9/9、TBLASTN 162/162。一致しない 1 件は S02 からの既知の `Sakai.MG1655.megablast` | `after/hashes.tsv`、`run-20260928T203929Z/capture-compare.txt` |
| v1 の reactor の検査 | 通過。`check_wasm_threading.py`：409 件の command / oracle の記録、reactor の寿命の検査 PASS、形式の失敗 0。S03 の成果物との比較：threaded reactor の 344 件（BLASTN の 2 つの task と TBLASTX の、1・2・4 スレッドの実行を含む）と serial reactor の 25 件の記録、`v1_requests.js` の 14 件の応答が一致 | `run-20260928T203929Z/wasm-threading.log`、`wasm-threading-metadata.json`、`wasm-threading-runs.json`、`v1-reactor-records-compare.txt`、`v1-requests-{before,after}.jsonl` |
| TBLASTX の `wasm-threads` だけの出力の経路を通る threaded の実行が serial と同じ | 通過。2 つの query × 4 つの subject（`LC738884` と `LC741431` の一部）で、threaded command-WASI の 4 スレッド（`LOSAT_WASI_THREADS_DEBUG=1` の診断 `threaded-4.debug.log` が、4 つの subject を並列に処理したことを示す）が、serial command-WASI、ネイティブの 1 と 4 スレッド、S03 の成果物の出力と一致（164 行） | `run-20260928T203929Z/tblastx-threaded-direct/` |
| V-NAT | 通過。BLASTN：3 つの query × 4 つの subject で、outfmt 6 と 7 を同時に要求した `run_local` の出力が、それぞれの形式だけの `run_local` と CLI に一致し、`diagnostics` は CLI の stderr に一致した（既定、`-max_target_seqs 1`（行が減ることを確かめる）、`-task blastn`、`-num_threads 2`）。TBLASTX：outfmt 6 を 2 つ同時に要求し、どちらも CLI に一致（1 と 2 スレッド。2 スレッドはネイティブの集約のスレッドを通り、1 スレッドの結果とも一致）、`-max_target_seqs 1`。観測者の範囲は各行に正確に一致する。未対応の形式は何も書かずに拒否される | `LOSAT/tests/run_local_blastn.rs`、`run_local_tblastx.rs` |
| V-PERF の非退行 | 調査の後、退行は無いと判断した（下の「性能の計測」）。**保守者の確認を求める**：計画の手順の 2 回の計測で +5% を超えた件がある（1 回目の `blastn` の serial WASI ×1.065、2 回目の `blastn-large-fmt7` の serial WASI ×1.081 と threaded WASI ×1.063。`perf-check-{1,2}.txt` の FAIL の行）。出力のハッシュはすべて同じ | `run-20260928T203929Z/perf-*.json`、`perf-check-*.txt`、`run-20260928T203929Z/perf-investigation/`、`perf_cases.py` |
| `cargo fmt --check`、`clippy -D warnings`、`cargo test --all-features` | 通過。clippy は `--all-features`、`--no-default-features`、`wasm32-wasip1`、`wasm32-wasip1-threads --features wasm-threads` の 4 つの構成で通過（`wasm-threads` だけでコンパイルされる TBLASTX の経路を含む）。試験は合計 796 件（S03 の 790 件と、BLASTN と TBLASTX の `run_local` の 6 件） | `run-20260928T203929Z/fmt.log`、`clippy.log`、`cargo-test.log` |
| wasm32 だけでコンパイルされる v1 の試験 | 通過（4 件） | `run-20260928T203929Z/wasm32-web-api-tests.log` |
| `docs/web/verification_cells.tsv` の両 program の行 | 埋めた | [`../../web/verification_cells.tsv`](../../web/verification_cells.tsv) |

成果物のハッシュと道具の版は `run-20260928T203929Z/binaries.txt` にある。コミットの木から作り直した 2 つのネイティブと 4 つの WASI の成果物は、ゲートに使った成果物と同じハッシュになった。同じハッシュの複製を `s04-final/` に残した。

### 性能の計測

変更前は S03 の成果物（`s03-final/`）、変更後はこの記録の成果物で、変更前と変更後を 1 回ずつ交互に測る（`measure_perf.py`、S02 で決めた方法）。case は、指示書の `blastn`（AP027152 × LC738884 の megablast。ヒットが無く、ネイティブで約 0.05 秒）と `tblastx`（LC738884 × LC741431、約 0.8 秒）に、`perf_cases.py` で Gate A の EDL933 × Sakai の megablast（約 0.8 秒、query 5.6 Mb）の outfmt 6（`blastn-large`）と outfmt 7（`blastn-large-fmt7`）を足した 4 つである。計画の手順の計測を 2 回、その後に 12 回の計測を 1 回行うことを、計測を始める前に決めた。すべての標本を残してある。

1. 計画の手順（1 回の暖機と 3 回の計測。`perf-{1,2}.json`、`perf-check-{1,2}.txt`）：1 回目は 12 件中 11 件が +5% 以内で、`blastn` の serial WASI が ×1.065（0.107 秒と 0.114 秒の、7 ミリ秒の差）だった。2 回目は 12 件中 10 件が +5% 以内で、`blastn-large-fmt7` の serial WASI が ×1.081、threaded WASI が ×1.063 だった。2 回目のこの 2 件は、変更前と変更後の両方が 1 回目より大きく遅い（serial WASI の変更前は 1 回目の 1.065 秒から 1.675 秒に、変更後は 1.071 秒から 1.811 秒になった）ので、その間の機械の状態が変わっていたと判断した。
2. 12 回の計測（`perf-investigation/perf-12.json`、`perf-12-check.txt`）：12 件すべてが ×0.909〜×1.018 だった。この計測は、3 回の計測で結論が出ない場合に備えて事前に決めていた。`AGENTS.md` は結論が出ない場合に限って回数を足すとする。1 回目には、同じ成果物の標本が二峰に分かれる件（変更後のネイティブの `blastn`：0.054〜0.280 秒）があり、結論が出ない場合に当たる。

計測の範囲：TBLASTX の case は subject が 1 つで、ネイティブは 1 スレッドなので、変更した集約の経路（ネイティブの集約のスレッド、threaded WASI の複数 subject の経路）の時間は測っていない（その経路の出力の一致は `tblastx-threaded-direct/` と `run_local_tblastx.rs` で確かめた）。CPU 時間とメモリの最大値は記録していない（BLASTN の検索は query の記録の写しを持つので、メモリは query の大きさの分だけ増える）。

結論：BLASTN と TBLASTX の経路に性能の退行は無い。超えた件は、起動が大半を占める短い実行のばらつきと、2 回目の途中の機械の状態の変化によると判断した。BLASTN の検索が query の記録の写しを持つことによる差は、5.6 Mb の query の `blastn-large` で見えなかった（12 回の計測でネイティブ ×0.994）。CPU 時間は記録していない。

## CLI の振る舞いの差（独立監査が見つけた。出力のバイトは変わらない）

1. **書き込みの失敗を報告するようになった。** BLASTN と TBLASTX は、出力を書き終えた後に `writer.flush()?` を呼ぶ（`blastn/blast_engine/run.rs`、`tblastx/blast_engine/run_impl.rs` の書き出し）。以前の書き出し（`blastn/hsp.rs` と `common.rs` のファイルを開く関数）は `BufWriter` を捨てるときに書き込みの失敗を黙って無視していた。報告の全体が 8 KiB のバッファより小さいとき、例えば `-out /dev/full` や、読み手が閉じたパイプへの stdout で、S03 は終了コード 0 で何も言わず、S04 は `No space left on device (os error 28)` や `Broken pipe (os error 32)` で終了コード 1 になる。BLASTP（S02 から）と TBLASTN（以前から）はすでにこの振る舞いである。書いたはずの出力が失われるのを黙って成功とするより正しいので、この差を残す（推奨）。**保守者の確認を求める。**
2. **診断用の stderr の出る順序が変わった。** CLI が入力を検索のスレッドプールの前に読むようになったため、`LOSAT_WASI_THREADS_DEBUG=1` の `[losat-thread-pool]` の行が BLASTN の `-verbose` の `Reading query & subject...` の後に出るようになり、BLASTN の query が空のときと、入力の読み込みに失敗したときには出なくなった。TBLASTX の `[DEBUG SCAN_OFF]` も読み込みに失敗したときには出ない。TBLASTX の `LOSAT_TIMING` は `read_subjects: 0.000s` を出し、`total` にファイルの読み込みを含まない。どれも LOSAT の開発用の診断で、NCBI の出力ではない。TBLASTX の `run` のコメント（読み込みの間に失敗しうる処理は無い）は、この診断の出力に触れていなかったので、S05 のコミットで直す。

## 既存の差と観察（この変更とは無関係。記録だけ）

- NCBI の `tblastx_app.cpp:114-124` は、subject を用意してから query を読む。LOSAT の TBLASTX と BLASTN は query を先に読む（以前から）。query と subject の両方が読めないときに報告される誤りが NCBI と違いうる。この変更はその順序を変えていない。
- `measure_perf.py` の `blastn` の case（AP027152 × LC738884 の megablast）は、ヒットが無く、ネイティブで約 0.05 秒で、起動の時間だけを測る。そのため `perf_cases.py` で Gate A の EDL933 × Sakai を足した。
- NCBI 2.17.0 は、`-subject` のファイルが無いと、query が空でも引数の処理の時点で拒否する（"File is not accessible"、終了コード 1）。LOSAT の BLASTN は、query が空だと subject を読まずに終了コード 0 で終わる（以前から。独立監査が確かめた）。これは承認された例外ではないので、後の課題とする。
- BLASTN の `-verbose` の出力は `ReportOutputs::diagnostics` を通らず、stderr に直接書かれる（以前から）。v2 のアダプタでは stream 3 に出ない（S05 の記録で扱う）。
- この記録を作る途中で、v1 の WASI の検査を一度、別の clone（このリポジトリの元の clone）のスクリプトで実行してしまった（シェルの `cd` が片方の背景の処理にしか効いていなかった）。その実行は記録の `head` と runner のハッシュで見分けられ、捨てて、この worktree のスクリプトでやり直した（`wasm-threading-metadata.json` の `head` は `a0ebbb849`）。元の clone のファイルは変わっていない。

## 独立監査

`ncbi_parity_auditor` の定義に従う読み取り専用の監査を受けた（コミット `458cbe110` と、この記録）。

- 結論：出力のバイトが変わらないこと、TBLASTX の出力の一本化、DW-10（BLASTX の経路を変えていない）、観測者の行の順序、成果物の出どころは supported。CLI の振る舞いは、上の 2 つの差（書き込みの失敗の報告、診断の stderr）を除いて supported。V-PERF は inconclusive で、保守者の判断が要る。
- 監査が独立に確かめたこと：コミットから作り直した 6 つの成果物のハッシュの一致、BLASTN と TBLASTX の 34 件の出力の S02 の基準との一致（`after/hashes.tsv` は 236 件すべて一致）、clippy、`run_local` の試験、試験への 4 つの変異（BLASTN の番号を query ごとに振り直す、集約のスレッドが HSP を 1 つ落とす、TBLASTX の警告を形式ごとに書く、観測者の `begin` を飛ばす）をどれも試験が検出すること、NCBI BLAST+ 2.17.0 との一致（BLASTN の outfmt 6/7 を 2 query × 5 subject、TBLASTX の 2 query × 4 subject の 2 つの fixture をネイティブの 1・2・4 スレッド、serial WASI、threaded WASI の 4 スレッドで。threaded の診断は `parallel_selected=true work_items=4` で、`wasm-threads` だけの経路を通った）、S03 と S04 の CLI の 110 件の誤りと境界の場合（上の差を除いてすべて一致）。
- 残った指摘と対応：
  - 上の 2 つの差は、この記録に書いた。1 は残すことを推奨し、保守者の確認を求める。2 の TBLASTX のコメントは S05 で直す。
  - `tblastx-threaded-direct/` に、`wasm-threads` の経路を通った証拠の診断が無い：`threaded-4.debug.log` を足した（`stage=subjects work_items=4 parallel_selected=true`）。
  - 凍結ハッシュの行はすべて 1 スレッドで、ネイティブの集約のスレッドは `run_local_tblastx.rs` の 2 スレッドと `tblastx-threaded-direct/` だけが通る。
  - 性能の計測の範囲と、12 回の計測を事前に決めたこと：上の「性能の計測」に書いた。
  - `v1_requests.js` は BLASTP だけを送る。BLASTN と TBLASTX の v1 は reactor の 344 件の記録で確かめている。
  - BLASTN と TBLASTX の試験の fixture は警告を出さないので、`diagnostics` の比較は空の stderr どうしである（変異で、形式ごとの警告は検出される）。
  - NCBI の引用の細部：`tabular.cpp:1100-1108` の引用が行の注釈（`// Add tab ...`）を省略の印なしに省いている、`blast_formatter.cpp:429-465` の最後の `}` は 466-467 行：S05 のコミットで直す。
- 例外として扱ったもの：S02 から続く `Sakai.MG1655.megablast`、`-out` のファイルを作る時点、警告を形式の数によらず 1 回だけ書くこと（計画 §4.3）。

## 実行の方法（再現）

S02 の記録（[`../losat_web_e1a/README.md`](../losat_web_e1a/README.md)）の「実行の方法」と同じ。変更前の成果物は S03 の成果物（`00fdfca14` の木から作ったもの）。性能は次のとおり。

```bash
python3 docs/evidence/losat_web_e1c/perf_cases.py run --before <S03 native>,<S03 serial.wasm>,<S03 threaded.wasm> --after <native>,<serial.wasm>,<threaded.wasm> --cases blastn,tblastx,blastn-large,blastn-large-fmt7 [--repeat 12] --out <file.json>
python3 docs/evidence/losat_web_e1c/perf_cases.py check <file.json>
```

## 引き継ぎ（S05 以降へ）

- 4 つの program の入口は `api::local_blast::run_local_{blastp,tblastn,blastn,tblastx}`（計画 §4.2）。program をまたぐ振り分けは S05 で足す。BLASTX は SX。
- BLASTN と TBLASTX は `hits` を呼ばない（S07・S08 まで）。v2 のアダプタは、HSP レコードが無い program を扱えること。
- BLASTN の `-outfmt` は、今は 6 と 7 だけで、検索のオプションは形式によらない。S07 で outfmt 0 を足すときは、NCBI が outfmt 0 のときに hitlist の大きさを `-num_descriptions` / `-num_alignments` から決める（`blast_args.cpp:2894-2978`）ことに合わせ、形式ごとの解決（BLASTP の `resolve_for_formats`）が要るかを確かめる。
