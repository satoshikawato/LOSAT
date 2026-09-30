# LOSAT Web E1b（Session S03）ゲート記録

- 段階：E1b 核の入口：TBLASTN（[総合計画書](../../losat_web_gui_plan.md) §7 の S03、[指示書](../../losat_web_gui_sessions/session_s03_e1b_core_entry_tblastn.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S02 の後の `fd06d5115`（エンジンは `afd57e7cb`）。エンジンと試験の変更はコミット `00fdfca14`（この記録の直前のコミット）
- 実行記録：[`run-20260928T194657Z/`](run-20260928T194657Z/)、変更後の出力のハッシュ [`after/`](after/)。基準は S02 の [`../losat_web_e1a/baseline/`](../losat_web_e1a/baseline/)（取り直していない）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)（`sha256sum --check evidence.sha256` をこのディレクトリで実行する）
- 判定：**完了条件を満たした**。ただし V-PERF は、計画の手順の計測で +5% を超えた回があるので、保守者の確認を求める（下の「性能の計測」）

## 変更の内容

| ファイル | 内容 |
|---|---|
| `LOSAT/src/algorithm/tblastn/args.rs` | `run_local`（どのターゲットでもコンパイルされる TBLASTN の入口）。入力を使わないオプションの検査を `search_settings` にまとめ、CLI の `run` と `run_local` の両方がそれを呼ぶ。CLI の `run` は、オプションを検査し、ファイルを読み、`run_local` を呼ぶ薄い層になった。CLI は以前と同じく出力の全体をメモリに書いてから stdout / `-out` に書く（`OutputSink::Writer`）ので、エラーのときに途中までの出力を残さない。無効な query の警告は `diagnostics` に 1 回だけ書く |
| `LOSAT/src/algorithm/tblastn/stage_e_report.rs` | `render` が `ReportOutputs` を受け取り、並べ替えた最終の結果から、要求された形式を順に書く。`PairwiseHit` の組み立て（`pairwise_hits`）を outfmt 0 の整形（`write_pairwise`）から分け、outfmt 0 か `hits` が求めるときに 1 回だけ作る。表形式の各行に観測者の印を付けた |
| `LOSAT/src/report/pairwise.rs` | `write_tblastn_pairwise_report`（TBLASTN だけが使う）に観測者を足した。節は、スコアの行、アラインメント、`x_DisplayAlnvecInfo` の最後の空行（`showalign.cpp:3980`）である |
| `LOSAT/src/algorithm/tblastn/mod.rs`、`LOSAT/src/api/local_blast.rs` | `run_local_tblastn` の再 export。`local_blast.rs` の `showalign.cpp` の引用の範囲を直した（S02 の 3 回目の独立監査の指摘） |
| `LOSAT/tests/run_local_tblastn.rs` | 試験 4 件（S02 の `run_local_support` を使う） |

HSP の番号 `HspIndex` は、並べ替えた最終の結果を query、subject、HSP の順にたどった順番で、表形式の行の順と同じである。`PairwiseHit` もこの順に作る。20,000 残基の query バッチの扱いと、並べ替え（`sort_hsps_for_report`）は変えていない。

## 完了条件と結果

| 完了条件（計画 §7 の S03） | 結果 | 証拠 |
|---|---|---|
| TLOSAN 計画の Stage G のゲートが変わらない | 変更後の 236 件（全 program）が S02 の基準と一致（差 0）。TBLASTN の 162 件（27 の遺伝暗号 × outfmt 0/6/7 × スレッド 2 通り）は Stage G の凍結ハッシュとすべて一致。一致しないのは S02 と同じ既知の `Sakai.MG1655.megablast`（BLASTN）だけ | `after/hashes.tsv`、`run-…/capture-compare.txt` |
| v1 の reactor の検査 | 通過。`check_wasm_threading.py`：409 件の command / oracle の記録、reactor の寿命の検査 PASS、形式の失敗 0。v1 には TBLASTN の入口が無いが、共有のコードを変えたので、S02 の成果物と比べた：threaded reactor の 344 件と serial reactor の 25 件の記録、`v1_requests.js` の 14 件の要求への応答が、S02 と一致 | `run-…/wasm-threading.log`、`wasm-threading-metadata.json`、`wasm-threading-runs.json`、`v1-reactor-records-compare.txt`、`v1-requests-{before,after}.jsonl` |
| V-NAT：3 形式の同時出力と単一形式の CLI との一致 | 通過。4 つの query（うち 1 つは無効）と 6 つの subject（同じ subject に複数の HSP を持つものを含む）で、outfmt 0・6・7 を同時に要求した `run_local` の出力が、それぞれの形式だけの `run_local` と CLI の出力に一致し、`diagnostics` は CLI の stderr に一致した。既定と `-max_target_seqs 1`（HSP の数が減ることを試験で確かめる）の両方で確かめた。未対応の形式（`5`、`6 qseq sseq`、`0 qseqid`）は、何も書かずに拒否される | `LOSAT/tests/run_local_tblastn.rs` |
| NCBI のバイトとの一致（複数形式を同時に出す経路） | 通過。Stage G の保存した NCBI BLAST+ 2.17.0 の出力（`batch_boundary` の 2 件と `native_real/av_code1`）と、3 形式を同時に出した `run_local` の出力が一致した。入力のパスは NCBI が読んだパスを表示名として渡した | 同上 |
| 警告が 1 回だけ出ること | 通過。3 形式を同時に出しても、`diagnostics` は NCBI の保存した stderr（無効な query 1 件または 2 件の警告）と一致した | 同上 |
| 観測者の範囲の一致 | 通過。outfmt 6 の行は HSP の番号の順に出力を隙間なく覆い、座標が HSP レコードと一致する。outfmt 7 の行は outfmt 6 の行と同じ。outfmt 0 の節は ` Score =` で始まり、その HSP の raw score を持ち、空行で終わる | 同上 |
| V-PERF の非退行 | 調査の後、退行は無いと判断した（下の「性能の計測」）。**保守者の確認を求める**：計画の手順の 2 回の計測で、各 1 件が +5% を超えた（1 回目の serial WASI の outfmt 0 が ×1.102、2 回目の threaded WASI の outfmt 6 が ×1.059。`perf-check-{1,2}.txt` の FAIL の行）。出力のハッシュはすべて同じ | `run-…/perf-*.json`、`perf-check-*.txt`、`run-…/perf-investigation/` |
| `cargo fmt --check`、`clippy -D warnings`、`cargo test --all-features` | 通過。試験は合計 790 件（S02 の 786 件と `run_local_tblastn` の 4 件）。`initial_tabular_oracle_bytes`（Stage E の保存したバイト）は、3 形式を 1 回の `render` で出す形に直して通過 | `run-…/fmt.log`、`clippy.log`、`cargo-test.log` |
| wasm32 だけでコンパイルされる v1 の試験 | 通過（4 件） | `run-…/wasm32-web-api-tests.log` |
| `docs/web/verification_cells.tsv` の TBLASTN の行 | 埋めた | [`../../web/verification_cells.tsv`](../../web/verification_cells.tsv) |

成果物のハッシュと道具の版は `run-…/binaries.txt` にある（`run-…` は `run-20260928T194657Z`）。ゲートは、コミットの最後の編集（`local_blast.rs` のコメントだけ）の前の木から作った成果物で実行した。コミットの木から作り直した 2 つのネイティブと 4 つの WASI の成果物は同じハッシュになった。`after/binary.json` などに書かれた成果物のパスは次のセッションのビルドで上書きされうるので、同じハッシュの複製を `s03-final/` に残した。

### 性能の計測

変更前は S02 の成果物（`s02-final/`）、変更後はこの記録の成果物で、`measure_perf.py` が変更前と変更後を 1 回ずつ交互に測る（S02 で決めた方法）。すべての計測の JSON は、各回の開始と終了の時刻と、捨てていないすべての標本を持つ。

1. 計画の手順（1 回の暖機と 3 回の計測）を、Stage G の benchmark の fixture（最初の AvCLPV の蛋白質 1 本 × AvCLPV のゲノム、`tblastn` は outfmt 6、`tblastn-fmt0` は outfmt 0）で 2 回行った（`perf-{1,2}.json`、`perf-check-{1,2}.txt`）。1 回目は serial WASI の outfmt 0 が ×1.102（変更後の標本 0.136・0.271・0.127 秒のうち 1 つが飛び離れている）、2 回目は threaded WASI の outfmt 6 が ×1.059（変更前 0.288・0.290・0.290 秒、変更後 0.307・0.311・0.298 秒）で、ほかの 10 件は ×0.935〜×1.027 だった。この fixture は、ネイティブで約 0.045 秒、serial WASI で約 0.12 秒、threaded WASI で約 0.29 秒で、起動の時間が大半を占める。
2. 同じ fixture で、計測を 12 回ずつにした（`perf-investigation/perf-12.json`。この計測は、1. を始める前に、1. の直後に続けて行うと決めていた）。6 件すべてが ×0.972〜×1.022 だった。threaded WASI の outfmt 6 は ×1.022 である。
3. 1. と 2. の結果を見た後に、調査として、TBLASTN の検索と整形が時間の大半を占める大きな fixture を足した。Stage G の full real code-1 cross（AvCLPV の蛋白質 120 本 × PsCLPV のゲノム、`perf_investigation.py` の `tblastn-full` と `tblastn-full-fmt0`）で、ネイティブで約 2.3 秒、serial WASI で約 2.7 秒、threaded WASI で約 3.5 秒である。3 回の計測（`perf-investigation/perf-full-3.json`）は 6 件すべてが ×0.953〜×1.009、12 回の計測（`perf-full-12.json`）は 6 件すべてが ×0.994〜×1.003 だった。

結論：TBLASTN の経路に性能の退行は無い。小さな fixture での超過は、起動が大半を占める 0.05〜0.3 秒の実行のばらつきによると判断した。ただし、threaded WASI の小さな fixture の outfmt 6 は、3 回の計測（×1.027、×1.059、×1.022）のすべてで 1 を超えていて、ばらつきだけでは説明しきれていない。大きな fixture では同じ場合が ×0.983 と ×0.997 なので、TBLASTN の検索と整形の差ではなく、起動（モジュールの準備やスレッドの起動）の数ミリ秒の差の可能性がある。これは確かめていない。変更後の WASI のモジュールは S02 より少し小さい（threaded command は 3,470,293 から 3,468,250 バイト）。CPU 時間は記録していない。

## 独立監査

`ncbi_parity_auditor` の定義に従う読み取り専用の監査を受けた（コミット `00fdfca14` と、この記録）。

- 結論：出力のバイトが変わらないこと、CLI の振る舞いが変わらないこと、複数形式の出力が NCBI の保存したバイトと単一形式の CLI に一致することは supported。HSP の番号の主張は、コードを読んだ結果としては supported だが、試験が部分的にしか確かめていない。性能は、計画の手順では inconclusive で、保守者の判断が要る（扱い方は V-PERF に沿う）。出力を変える欠陥は無かった。
- 監査が独立に確かめたこと：コミットから作り直した成果物のハッシュの一致、BLASTP・TBLASTN・BLASTX の 202 件の出力の S02 の基準との一致、S02 と S03 と作り直したネイティブの 3 つでの CLI の 69 件の誤りと境界の場合の一致（stdout、stderr、終了コード、`-out` のファイルができるか。すべて無効な query、20,001 残基の無効な query、書けない `-out`、壊れた FASTA などを含む）、`search_settings` が 2 回呼んでも同じであること、HSP の順序、移したコードが以前と同じであること、NCBI の引用の行、BLASTX に触れていないこと。試験への変異（形式ごとに警告を書く、`-outfmt` の事前の検査を外す、表形式の番号を query ごとに振り直す、outfmt 0 の番号を逆にする）は、どれも試験が検出した。
- 残った指摘と対応：
  - `struct SearchSettings`（`args.rs`）に NCBI の参照が無い。`tabular.cpp:1100-1108` と `blast_formatter.cpp:429-465` の引用に省略の印が無い：S04 のコミットで直す。
  - outfmt 0 の節の試験は raw score だけを比べるので、同じ raw score の HSP の番号の入れ替えを検出できない：S04 で、節の中の Query と Sbjct の座標と、節がその HSP の query の区切りの中にあることを確かめる検査を `run_local_support` に足す（TBLASTN と BLASTP の試験が使う）。
  - `evidence.sha256` が無い、判定と監査の欄が空、serial のネイティブの出どころ：足した（serial のネイティブもコミットから作り直して同じハッシュ）。
  - CLI が失敗したときに出力を残さないことの自動の試験は無い（監査が手で確かめた）。
  - TBLASTN のネイティブと WASI の一致は、性能の計測の fixture の出力のハッシュだけで確かめている。V-ABI は S05 で行う。
- 例外として扱ったもの：`PD-TLOSAN-LOCAL-GENCODE-32`（遺伝暗号 32）、S02 から続く `Sakai.MG1655.megablast` の Gate A との不一致、`-out` のファイルを作る時点（S02 から記録）、警告を形式の数によらず 1 回だけ書くこと（計画 §4.3 の設計。NCBI は `CBlastFormat` ごとに出す）。

## 実行の方法（再現）

S02 の記録（[`../losat_web_e1a/README.md`](../losat_web_e1a/README.md)）の「実行の方法」と同じ。変更前の成果物は S02 の成果物（`afd57e7cb` の木から作ったもの）である。TBLASTN の性能は `--cases tblastn,tblastn-fmt0` で測る。調査の大きな fixture は `perf_investigation.py`（`measure_perf.py` に case を 2 つ足したもの）で測る。

```bash
python3 docs/evidence/losat_web_e1a/measure_perf.py run --before <S02 native>,<S02 serial.wasm>,<S02 threaded.wasm> --after <native>,<serial.wasm>,<threaded.wasm> --cases tblastn,tblastn-fmt0 --out <file.json>
python3 docs/evidence/losat_web_e1b/perf_investigation.py run --before ... --after ... --cases tblastn-full,tblastn-full-fmt0 [--repeat 12] --out <file.json>
```

## 引き継ぎ（S04 以降へ）

- TBLASTN の入口は `api::local_blast::run_local_tblastn(args, queries, subjects, outputs)`。表示名は `-query` / `-subject` の値から取る（ラベルの引数は持たない。計画 §4.2）。
- `docs/evidence/losat_web_e1a/` の `measure_perf.py` と `capture_outputs.py` は、S02 の `evidence.sha256` が覆うので変えない。case を足すときは、この記録の `perf_investigation.py` のように、読み込んで足す。
- Stage G の benchmark の fixture（`tblastn`、`tblastn-fmt0`）は、ネイティブで約 0.05 秒、threaded WASI で約 0.3 秒で、起動の時間が大半を占める。
