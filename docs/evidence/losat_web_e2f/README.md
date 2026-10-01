# LOSAT Web E2f（Session S07++ / S07++b）ゲート記録

- 段階：E2f BLASTN の query の batch と query の分割（[総合計画書](../../losat_web_gui_plan.md) §7 の S07++、指示書 [S07++](../../losat_web_gui_sessions/session_s07pp_e2f_blastn_query_batches.md)・[S07++b](../../losat_web_gui_sessions/session_s07ppb_e2f_close.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S07+ の記録の `e8de6861b`（エンジンは `bbc9c7809`、実行ファイルは S07+ の成果物 `s07p-final/LOSAT`、SHA-256 `7e5cf78a0db5…`）
- エンジン：`93c58a2e4`（batch）、`bbb869a9e`（query の分割）、`74223fee3`（予備の段階の hit list。独立監査の第 1 回の修正）。ゲートを実行した `7693a9c73` から HEAD まで、エンジン（`LOSAT/`）と `web/adapter/` は変わっていない（`git diff --stat 7693a9c73 HEAD -- LOSAT web/adapter` が空）ので、ゲートの実行は最後のエンジンの実行である
- アダプタ：S07++ で変えていない（`git log 93c58a2e4~1..HEAD -- web/adapter` が空。最後の変更は S07+ の `e4409ccfa`）
- 権威の記録：[`AUTHORITY.md`](AUTHORITY.md)（NCBI の batch の経路 §A、batch ごとに決まるもの §B、LOSAT の対応 §C、query の分割 §D、予備の段階の hit list §E）
- 実行記録：[`run-20260930T163727Z/`](run-20260930T163727Z/)（ゲート。`head.txt` は `7693a9c73`）、[`run-20261001T130537Z/`](run-20261001T130537Z/)（独立監査の第 1 回の再現、`head.txt` は `7a2160754`）、[`audit_round2/`](audit_round2/)（独立監査の第 2 回の比較）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)（リポジトリのルートで `sha256sum --check docs/evidence/losat_web_e2f/evidence.sha256` を実行する）
- 判定：**完了条件を満たした**。独立監査は 2 回目で supported（1 回目は unsupported）。V-PERF の最初の計測で 15 件すべてが閾値の内だった（下の「性能の計測」）

## 変更の内容

### エンジン

| 範囲 | 内容（NCBI の経路は `AUTHORITY.md`） |
|---|---|
| query の batch（`93c58a2e4`） | `CBatchSizeMixer` と `GetNextSeqBatch` を `blastinput/query_batch.rs` に移植した（`(Int4)` の変換は x86-64 の実行ファイルと同じ）。最初の batch は約 5000 残基で、次の batch の大きさは前の batch の `good_init_extends`（subject の chunk ごとに保存した初期の hit の数）から決める。`run.rs` の `run_in_pool` は subject を 1 度だけ読み符号化し、batch ごとに query の塊・lookup table・対角線の表・gapped の X-drop・予備の段階の `-subject_besthit` を作って検索し、全 query の結果を 1 度に書く（§A〜§C） |
| S07+ の拒否の置き換え | 最初の batch だけを推定した `first_query_batch`、検索されなかった batch を推定した `unsearched_queries`、表を超える gap で query が最初の batch に収まらない場合の拒否、Karlin-Altschul の表の失敗の文言の繰り返しの回数を、本物の batch で置き換えた（§C） |
| query の分割（`bbb869a9e`） | batch の query の総文字数が 2 × (chunk − 100) 以上のとき（`-task blastn` で 1,999,800、megablast で 9,999,800 文字以上）、NCBI は batch を重なり 100 の塊に分けて予備の段階を塊ごとに行い、HSP を batch の context に移して合わせる。`blastn/query_split.rs`（塊の数・範囲・部分・context・補正・mask の制限）と、`run.rs` の `search_query_chunks`・`merge_query_chunk`・`merge_prelim_hit_lists` で移植した。`Blast_HSPListsMerge` の帯の HSP の入れ替えは、subject の分割の枝にも入れた（§D） |
| 予備の段階の hit list（`74223fee3`） | NCBI は予備の段階で subject の HSP の一覧を query ごとに大きさ `prelim_hitlist_size` の hit list（既定で 550、`-max_target_seqs` 1〜5 で 10）に保存し、traceback は残った一覧だけを読む。LOSAT は subject ごとに予備の段階と traceback を続けて行い、大きさの制限を traceback の後の e-value に適用していた。collector（`collect_prelim_hit_lists`）と `Blast_HitListUpdate`・`Blast_HitListMerge` の移植に置き換え、traceback の後の hit list の大きさは `hitlist_size` にした。この差は S07++ の前からあり、分割しない batch にも出ていた（§E） |
| 明示的な拒否 | 1 つだけ残る：得点の表が無い得点で、最初の batch が無効な query だけのとき（NCBI は、その batch の結果と警告を書いた後に、次の batch で誤りを出す。§C）。LOSAT のそのほかの明示的な拒否は S07+ の `AUTHORITY.md` §G のまま |
| 試験 | エンジンの試験は 823 件（S07+）から 830 件になった（`cargo-test.log`） |

### アダプタと文書

- `web/adapter/`：S07++ で変えていない。
- `docs/evidence/losat_web_e2f/`：`AUTHORITY.md`、再現スクリプト `batch_sweep.py`・`check_inputs.py`・`split_check.py`（`check_inputs.py` と `batch_sweep.py` は S07+ の `check_inputs.py`・`slice_sweep.py` を読み込む）。
- 計画：TD-14（S07++、BLASTN の query の batch）、DW-12（NCBI の経路を棚卸しして未移植の部分を一括で transpile する方針。query の分割を拒否でなく移植で入れた根拠）。

## 完了条件と結果

| 完了条件（計画 §7 の S07++、指示書） | 結果 | 証拠 |
|---|---|---|
| 複数の query の sweep が NCBI とバイト一致（batch の境で結果が変わる入力） | 通過。`batch_sweep.py`（subject の窓 5〜40 kb と 2〜40 本の query、5 つのオプション × 150 件 = 750 件）は差 0。S07+ の実行ファイルでは `-task blastn` で 1 件が違った（NCBI 123 行、S07+ の LOSAT 124 行。既定のオプションは 150 件とも一致）。NCBI に `BATCH_SIZE` を与えた比較は行っていない（`batch_sweep.py` は NCBI の既定の batch と比べる） | `batch-sweep-batches_*.tsv`、`batch-sweep-batches_task_blastn-before.tsv` |
| query の分割が NCBI とバイト一致 | 通過。`split_check.py` 213 件が差 0。そのうち 15 件は、NCBI の出力が分割で変わる case（探索空間、mask の細部、境の HSP）で、LOSAT も同じに変わる。第 1 回の指摘の再現 7 件（`r1_`）は HEAD の実行ファイルで差 0 | `split-check.tsv`、`run-20261001T130537Z/split-check-r1-head.tsv` |
| 予備の段階の hit list が NCBI とバイト一致（独立監査） | 通過。`audit_round2/` の 325 件のうち 321 件がバイト一致、違い 0、LOSAT の明示的な拒否 4 件（下の「独立監査」） | `audit_round2/results.tsv` |
| S07+ の検査に退行なし | 通過。`check_inputs.py` 294 件がすべて期待どおり（S07+ は 292 件。S07+ で LOSAT が拒否していた 14 件が NCBI とバイト一致の `same` になり、分割の 2 件を足した）。得点の sweep は outfmt 0・6・7 のそれぞれで same 300 / same-error 580。S06 の `scoring_sweep.py` は 80 組合せのうち 46 が outfmt 6 のバイト一致、34 が NCBI と同じ終了コードの拒否。word size の sweep 256 件・切り出しの sweep 17 の設定は差 0。title の sweep 1023 件は 957 がバイト一致、66 が NCBI の落ちる定義行の明示的な拒否（S07+ と同じ数） | `check-inputs.tsv`、`scoring-sweep-fmt{0,6,7}.tsv`、`scoring-sweep-s06.tsv`、`word-size-sweep.tsv`、`slice-sweep-*.tsv`、`title-sweep.tsv` |
| S07 の fixture に退行なし | 通過。47 件すべてが `-num_threads` 1・2・4 で stdout と stderr まで NCBI BLAST+ 2.17.0 とバイト一致、outfmt 6 も一致（`precheck.tsv`）。NCBI の再実行は manifest のハッシュと一致 | `check-losat-n{1,2,4}.tsv`、`precheck.tsv`、`oracle-check.log` |
| 既存の BLASTN のゲートに退行なし | 通過。全 program の 236 件の出力・stderr・終了コードが S02 の基準と一致（差 0）。Gate A の BLASTN の行を含む | `capture-compare.txt`、`capture/hashes.tsv` |
| BLASTN の V-ABI | 通過。full：115 の検索を serial reactor（1 スレッド）と threaded reactor（1・2・4 スレッド）で実行した 460 件すべてで、全形式・HSP レコード・診断がネイティブの CLI と一致。凍結ハッシュは 676 件中 672 件が一致し、一致しない 4 件は S02 からの既知の `Sakai.MG1655.megablast` の outfmt 7。quick：52 件すべて一致 | `v-abi-full/summary.json`、`v-abi-quick/v-abi.log` |
| V-PERF の記録 | 通過（下の「性能の計測」） | `perf-1.json`、`perf-check-1.txt`、`split-timing.txt` |
| （規則 4）既存のゲート | 通過。v1 の WASI の検査（423 件の command / oracle の記録、形式の失敗 0）、`v1_requests.js` の 14 件の応答が S05 の記録と一致、wasm32 だけでコンパイルされる v1 の試験 5 件が通過。v1 の reactor の記録は、threaded の 172 件のうち 21 件の誤りの文言が S07+ で足した句だけ違い（S07+ と同じ 21 件）、serial の 12 件は一致 | `wasm-threading.log`、`v1-reactor-records-compare.txt`、`v1-requests-compare.txt`、`wasm32-web-api-tests.log` |
| `cargo fmt --check`、`clippy -D warnings`、`cargo test --all-features` | 通過。エンジンは 4 つの構成の clippy と 830 件の試験、アダプタは 3 つの構成の clippy と 7 件の試験 | `fmt.log`、`clippy.log`、`adapter-clippy.log`、`cargo-test.log`、`adapter-test.log` |
| ビルドの同一性（TD-6、TD-11） | 通過。両方の reactor で、リンクの引数・profile・共有の依存・target の rustflags が一致（失敗 0） | `reactors/build-identity.json` |
| 独立監査 | 第 2 回で supported | 下の「独立監査」 |

## 独立監査

読み取り専用の独立監査（役割 `ncbi_parity_auditor`、コードを変えない別のエージェント）を 2 回受けた。見つかったものの NCBI の経路は `AUTHORITY.md` §E にある。

| 回 | 結果と主な指摘 | 修正 | 記録 |
|---|---|---|---|
| 1（`b951afeb2` で実施） | unsupported。指摘 1（高）：`prelim_hitlist_size` の制限を traceback の前でなく後に適用していた。S07++ の前からの差で、分割しない batch にも出る（既定の設定で 560 の subject、`-max_target_seqs` 1 と 3、分割しない 300,000 文字の query で再現）。指摘 2（低）：`split_check.py` が megablast の分割を実行していなかった（2 つの query は NCBI では 2 つの batch になる）。注記：分割した batch の lookup table を NCBI は作らない（`blast_aux_priv.cpp:206-207`）ので注釈を直した | `74223fee3`、`7693a9c73` | §E、`run-20261001T130537Z/` |
| 2（`7693a9c73` のエンジン、2026-10-01） | supported。下の「第 2 回の比較」 | — | `audit_round2/` |

**第 1 回の再現（`run-20261001T130537Z/`）。** 再現の入力のスクリプトは一時ディレクトリにあって失われていたので、`split_check.py` の `r1_` 7 件にした（(a) 560 の subject と境をまたぐ 1 本、(b) `-max_target_seqs` 1・3 と対照の 10・11、(c) 分割した megablast、(d) 分割しない batch）。HEAD の実行ファイルは 7 件とも NCBI とバイト一致。`bbb869a9e` の実行ファイル（予備の hit list の移植の前、`s07pp-split1-native`）は 5 件（(a)、(b) の 1 と 3、(c)、(d)）で違い、対照の 10・11 は一致した。

**第 2 回の比較（`audit_round2/`）。**

- 比較 325 件（`results.tsv`）：321 件がバイト一致、違い 0、LOSAT が明示的な誤りで拒否したもの 4 件（`-xdrop_gap 10` と `-ungapped` が 2 件ずつ。NCBI は実行する）。ゲートと同じ実行ファイル（`a4ff5abb3a61…`）で実行した。`results-rerun.tsv` は `W_other_options` の 15 件を除く 310 件の再実行で、325 件の該当行と同じ。
- 入力が古い振る舞いと修正を分ける：`74223fee3` の親の実行ファイルでは、同じ 325 件のうち 95 件が NCBI と違った（`results.tsv` の `prefix_binary_matches_ncbi=NO` は 99 行で、95 件の違いと 4 件の拒否）。NCBI の hit subject が 560 本以上の入力が 166 件（550 を超えるものは 172 件）。
- 範囲（`results.tsv` の category と argv から数えた）：`-max_target_seqs` 1〜5 で hit subject が 10 を超える 154 件、6・10・11 の 22 件、100・550・600 の 13 件。`-subject_besthit` 27 件、`-max_hsps` 26 件、`-num_threads 4` 20 件、`-evalue` を与えた 17 件。outfmt 7 が 16 件、outfmt 0 が 14 件。megablast 54 件、`-task blastn` 271 件。query を分割する batch は category `G_split`・`S_split_gap`・`T_chunk_overflow` の計 67 件（EDL933[0:2M]、EDL933 と Sakai をつないだ 11 Mb の megablast、EDL933 全体の 5 つの塊）。複数の query の category `H_multi_query` 13 件。
- NCBI のソースと LOSAT を比べて一致と判断したもの：collector（`hspfilter_collector.c:83-161` と `run.rs` の `collect_prelim_hit_lists`）、heap（`blast_hits.c` の `s_EvalueComp`・`s_EvalueCompareHSPLists`・`s_Heapify`・`s_CreateHeap`・`Blast_HitListNew`・`Blast_HitListUpdate` と `hsp.rs`。同点は e-value（どちらも 1e-180 未満なら等しい）、最初の HSP の得点、subject の番号の順）、予備の e-value の key と合わせた HSP の e-value、`Blast_HitListMerge`、traceback の読む順（subject の番号の昇順）と最後の hit list の大きさ、オプションの組合せ。
- 読んだが実行していない NCBI の枝：組成補正（`ADAPTIVE_CBS`）の予備の hit list の大きさ、ungapped の枝（`-ungapped` は LOSAT が拒否する）、database の mode の複数スレッドの予備の順、最良の e-value の HSP が一覧の最初の HSP と違う場合（サポートする得点の blastn の CLI からは出ない）、`hsp_num_max` が `INT4_MAX` 未満の場合。LOSAT の他の明示的な拒否（`-qcov_hsp_perc`、`-culling_limit`、`-best_hit_overhang`、`-best_hit_score_edge`、`-min_raw_gapped_score`、`-xdrop_gap_final`、`-strand`、`-searchsp`、`-dbsize`、`-window_size`、`-off_diagonal_range`、`-xdrop_ungap`、`-no_greedy`）は LOSAT の側だけを確かめた。入力は EDL933 と Sakai だけである。

## 性能の計測

変更前は S07+ の成果物（`s07p-final/`）、変更後はこの記録の成果物で、1 回の暖機の後に 3 回ずつ交互に測った（`perf_cases.py`、`--repeat` は 3）。case は E2c と同じ：`blastn`、`blastn-large`・`blastn-large-fmt7`・`blastn-large-fmt0`（EDL933 × Sakai の megablast）、`blastn-many`（260 の query、`-task blastn`）。

- **V-PERF（`perf-check-1.txt`）：** 15 件すべてが +5% の内（×0.736〜×1.050、出力はすべて同じ）。最大は `blastn-many` の native（×1.050。比は 1.0496、変更前 [0.127, 0.131] 秒、変更後 [0.125, 0.138] 秒で範囲が重なる）、次が `blastn` の serial WASI（×1.047）。最小は `blastn-many` の serial WASI（×0.736）。`blastn-many`（batch ごとに subject を検索し直す case）は、serial WASI が 0.347 秒から 0.256 秒、threaded WASI が 0.361 秒から 0.358 秒（×0.994）。
- **再計測は行っていない：** 閾値を超えた case が無いので、AGENTS.md の手順（結論が出ないときだけ回数を増やす）の 5 回の再計測は要らなかった。
- **分割の時間（`split-timing.txt`）：** EDL933（5,528,445 文字、5 つの塊）× Sakai[1000000:1500000]、`-task blastn -outfmt 6`、1 スレッド、3 回交互。NCBI 1.56・1.64・1.58 秒（1.56〜1.64）、LOSAT 1.71・1.67・1.51 秒（1.51〜1.71）、出力は同じ（1919 行）。
- 独立監査の後にエンジンは変えていないので、このゲートの計測が最後のエンジンの計測である。計測の間、アプリ側の試験とビルドは止めた（計画 DW-7、`vperf.lock`）。

## 残件と扱い

- **単体試験の不足（重大度 低、独立監査の第 2 回の指摘）：** `74223fee3` は `#[test]` を足していない。予備の hit list の単体試験は `hsp.rs` の `test_hitlist_update_keeps_best` だけで、大きさ 1 の list と異なる e-value の場合に限る。subject の番号による同点の順、e-value・得点が等しい場合、1e-180 の規則、2 つ以上の list の heap、collector、`merge_prelim_hit_list`、最後の hit list の大きさの上限は固定されていない。回帰を守っているのは、NCBI の実行ファイルを要る比較（`split_check.py` の `r1_` と `audit_round2/`）である。S07++b では足さなかった。試験だけの変更でもソースの行が動いて実行ファイルが変わり、ゲート全体をやり直すことになるため。**扱い：** S07+++（E2g）で、次のエンジンの変更と同時に足す。
- **S07+ の実行ファイルの差 1 件：** `batch_sweep.py` の `-task blastn` で、S07+ の実行ファイルは 1 件が NCBI と違った（NCBI 123 行、LOSAT 124 行）。S07++ の前の振る舞い（batch の境の近似の ungapped 伸長）で、S07++ の実行ファイルでは一致する。対応は要らない。
- **残る明示的な拒否 1 件：** 得点の表が無い得点で、最初の batch が無効な query だけのとき（§C）。NCBI は、その batch の結果を書いた後に次の batch で誤りを出す。**扱い：** 拒否のまま。S07+++ の棚卸しの表で、拒否か承認済みの例外として載せる。
- **証拠の範囲：** 第 2 回の入力は EDL933 と Sakai だけである。`BATCH_SIZE` を与えた NCBI との比較は行っていない（`batch_sweep.py` と `split_check.py` は NCBI の既定の batch と、`split_check.py` は `CHUNK_SIZE=20000000` も比べる）。
- **S07+++（E2g）への引き継ぎ（DW-12）：** BLASTN の経路の棚卸しと、未移植・差のある移植の一括の transpile。確かめる候補は、ほかの `qsort` を Rust の安定な並べ替えにした箇所（オラクルの glibc 2.39 の `qsort` は安定な merge sort）と、`run.rs` の最初の hit の並べ替えが使う `sort_unstable_by`。
- **merge：** `main` への PR と merge は、この記録の commit の後に行う（保守者は merge を承認済み）。
- **アプリ側：** 中断した S09 の作業（`LOSAT-web-gui-app`）には触っていない。S07+++ の後に別のセッションで行う。

## 実行の方法（再現）

```bash
# NCBI との比較（NCBI BLAST+ 2.17.0 を比較のときに実行する）
python3 docs/evidence/losat_web_e2f/batch_sweep.py --bin-dir <NCBI> --losat <LOSAT> --work <dir> \
  --options="<options>" --cases 150 --jobs 8
python3 docs/evidence/losat_web_e2f/split_check.py --bin-dir <NCBI> --losat <LOSAT> --work <dir> \
  [--only r1_] --jobs 6
python3 docs/evidence/losat_web_e2f/check_inputs.py --bin-dir <NCBI> --losat <LOSAT> --work <dir>
# 第 2 回の比較（audit_round2/ から。NCBI と LOSAT のパスは driver.py の先頭にある）
python3 driver.py [name-prefix] [out.tsv]
# S07+ の sweep、fixture、回帰の出力、v1 の WASI の検査、V-ABI：S07+ の記録の「実行の方法」と同じ
python3 docs/evidence/losat_web_e2c/perf_cases.py run --before <S07+ の 3 つ> --after <S07++ の 3 つ> \
  --cases blastn,blastn-large,blastn-large-fmt7,blastn-large-fmt0,blastn-many --out <file.json>
```

`<NCBI>` は NCBI BLAST+ 2.17.0 の `bin`、`<LOSAT>` は `cargo +1.92.0 build --release --locked` の実行ファイル。ゲートの実行ファイルのハッシュは `run-20260930T163727Z/artifacts.sha256`。

## 引き継ぎ（S07+++ へ）

- 次は S07+++（E2g、[指示書](../../losat_web_gui_sessions/session_s07ppp_e2g_blastn_inventory.md)）。S07++ のコミット、`AUTHORITY.md` §A〜§E、上の「残件と扱い」を前提にする。
- 変更前の成果物は、この記録の成果物（ハッシュは `run-20260930T163727Z/artifacts.sha256`）。
