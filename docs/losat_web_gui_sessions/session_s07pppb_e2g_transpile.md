# Session S07+++b — E2g の続き：一括の transpile、試験、独立監査

## INSTRUCTION PROMPT

LOSAT の段階 E2g の続きを実行する。S07+++ は、CI の整備と、NCBI の BLASTN の経路の棚卸し（`docs/evidence/losat_web_e2g/INVENTORY.tsv`、1006 行）までを終えた。このセッションでは、棚卸しで決めた一括の transpile（と 2026-10-02 の保守者の決定で足した項目）を行い、試験・ゲート・V-PERF・独立監査を通して、E2g の完了条件を満たす。完了条件の正本は総合計画書 §7 の S07+++ の行で、[S07+++ の指示書](session_s07ppp_e2g_blastn_inventory.md)がその細部である（条件を緩めない）。

先に次を読む：
- [セッション README](README.md) の共通規則（特に規則 4）
- S07+++ の指示書（「進め方」と 4.「試験」）
- [E2g のゲート記録](../evidence/losat_web_e2g/README.md)
- `docs/evidence/losat_web_e2g/stage2/REVIEW.md`（項目ごとの確認と判断）

### 0. CI の残り（最初に行う。CI の設定と検査の道具だけを変え、エンジンとは別のコミットにする）

1. **毎晩の WASI の TIMEOUT を解く。**
   - 事実：`wasm-threading.yml` の `workflow_dispatch`（run `36882087513`）で、許可リストは働いた（Sakai の 4 つの label が allowed、閾値の 9 件が PASS）。しかし `threaded/tblastx/p11_avclpv_psclpv` が TIMEOUT（1 検索の期限 3600 秒）で失敗した。native は PASS。前回の赤い run `36745543144` でもこの label は PASS を出していない。最後に緑だった run `35201795582`（2026-09-17、`implement-v010-distribution`）では 1 時間以内に PASS した。
   - 速度の退行か runner の差かを確かめる：
     - 手元で threaded Wasm の p11 を、`check_wasm_threading_regressions.py` と同じコマンド（4 スレッド）で、HEAD の成果物と、残っている前の成果物（`~/.cache/losat-web-gui-target/s07p-final/` など）で 1 回ずつ測る。
     - run `35201795582` と `36882087513` のログから、threaded p11 の開始と終了の時刻を取る。
     - 測定は sonnet の agent に回してよい。
   - 退行なら：原因の範囲を記録し、S08（TBLASTX）の指示書の最初の作業にする。
   - runner の差なら：`wasm-threading.yml` の `--timeout-seconds` を根拠つきで上げる。
   - どちらでも、毎晩の検査が何で赤いかを記録に残す。
2. **緑の PR の実行時間（cache あり）を測る。** このセッションの最初の PR の run を記録する。
3. **`main` への PR を作る。** S07+++ の CI の変更（`db67e02cc` 以降）と、このセッションの成果を含める。merge は保守者が承認済みの流れ（S07++b と同じ）に従う。`nightly.yml` は `main` に入ってから `workflow_dispatch` できる。merge の後に 1 回 dispatch し、緑を確かめる。

### 1. 一括の transpile（設計と実装はメインの Opus。エンジンのソースは 1 つずつ変える）

各項目の確認の根拠は `stage2/REVIEW.md` にある。行番号は S07+++ の終わり（`2fabad16c`）の値。

| 項目 | NCBI | LOSAT | 移植すること | NCBI と比べる入力 |
|---|---|---|---|---|
| T1 初期の hit の順 | `blast_extend.c:274-310`（`score_compare_match`）、`na_ungapped.c:1690-1691` | `run.rs:2811` `score_compare_ungapped_hits`、`run.rs:9719` `sort_unstable_by` | `q_start` を連結した query の位置（`query_context_offset + qs`。2 つの作り手 `run.rs` ~8439・~9486 は context の中の位置を入れる）で比べ、安定な `sort_by` にする（glibc 2.39 の `qsort` は安定な merge sort） | X と revcomp(X) をつないだ query など、両鎖に得点・subject の始まり・長さの等しい hit が出る入力 |
| T2 gapped の開始点 | `blast_gapalign.c:3323-3389`（`BlastGetStartForGappedAlignmentNucl`） | `alignment/gapped.rs:493`（呼出し `run.rs:10596`） | `offset`・`q_start`・`s_start`・`q_len` を `Int4` の符号付きで計算する（`usize` の差は負で折り返す）。範囲の守りは残す | `-task blastn` の小さい `-word_size`（4〜7）で、DP の開始点の近くに mismatch がある入力（word size の sweep） |
| T3 予備の e-value の刈り込み | `blast_engine.c:640-672` | `run.rs` ~10458-10470 | `prelim_evalue <= threshold` を C と同じ `!(prelim_evalue > threshold)` にする（NaN は LOSAT が拒否するので出力は変わらない） | 単体試験だけ |
| T4 traceback の順 | `blast_traceback.c:358-365`（得点の並べ直しは `_DEBUG` だけ） | `run.rs` ~10530 `sort_prelim_hits_by_score` | traceback の前の並べ直しを除く。heap にしない list は chunk ごとの得点の並べ替え（`blast_engine.c:555`、LOSAT `run.rs` 2953・10161）で既に得点の順になっている。heap にした list（550 を超える subject）は `HitList` の e-value の順を保つ。subject のすべての query で 1 つの区間木を使う点（REVIEW の E2-note）も確かめる | hit のある subject が 550 を超え、表を超える gap（context ごとに gapped block が違う）の得点。`scoring_sweep` の case から選ぶ |
| T5 lookup の cell の順 | `blast_lookup.c:33-75`（`BlastLookupAddWordHit` は後ろに足す）、`BlastLookupIndexQueryExactMatches`、`blast_nalookup.c:197-300`（`s_BlastSmallNaLookupFinalize`）と `s_BlastNaLookupFinalize`。megablast の table は新しい順（`blast_nalookup.c:1090-1091`） | `lookup.rs`（小さい table を megablast の鎖で持つ。~640-870） | eSmallNaLookupTable と eNaLookupTable の cell で、query の位置を昇順に並べる。megablast の table は変えない | 既定の blastn（word 11、lut 8）で 4〜16 kb の query（対角線の hash を使う）。1 つの subject の位置で seed が複数出る反復と、複数の query |
| T6 対角線の配列 | `blast_parameters.c:225-231` | `run.rs:7239-7242`（`queries.len() == 1` と `MAX_ARRAY_DIAG_SIZE`） | query の塊が 8000 以下なら、query の数によらず配列を使う。NCBI の「最後の context の位置 + 長さ」と LOSAT の `query_concat_length` が同じものか、配列の経路が複数の query で連結した位置の対角線を使うかを確かめる | 2〜20 本の短い query（片鎖の計 4000 以下）。小さい query の `batch_sweep.py` |
| T7 batch と分割の環境変数 | `blast_input_aux.cpp:85-91`（`BATCH_SIZE`。0 は無いのと同じ：`blastn_app.cpp:261-300`）、`local_blast.cpp:54-62`（`CHUNK_SIZE`、`IsBlank`）、`split_query_aux_priv.cpp:51-61`（`OVERLAP_CHUNK_SIZE`）、`ncbistr.cpp:635-643,798-870`（`StringToInt`：符号と 10 進数、`int` の範囲、それ以外は `CStringException`） | `run.rs:5448-5469`（今の拒否）、`query_split.rs:23-39`（定数）、`run.rs` `run_in_pool`（5598 から） | NCBI と同じく読む。NCBI の経路が普通の batch と分割の計算になる値は移植する。整数でない値（S07+++ の range A の agent の実行では、NCBI は `CStringException` で終了コード 255）・負の値・NCBI の `size_t` の計算が折り返す値は、NCBI の文言と終了コードを再現できればそうし、できなければ明示的に拒否する。E2c §G の 2 つの拒否がなくなる | NCBI にも同じ変数を与えた `batch_sweep.py`・`split_check.py`（`BATCH_SIZE` 100・1000・5000・100000、`CHUNK_SIZE` 40000・300000、`OVERLAP_CHUNK_SIZE` 6・50・500）。S07++ の証拠の範囲の穴（`BATCH_SIZE` を与えた比較が無い）も埋まる |
| T8 得点の表の失敗（後の batch） | `setup_factory.cpp:170-185`、`blast_setup.c:73-107`、`blastn_app.cpp` の batch の繰り返しと例外の経路 | `run.rs:6296-6305`、`run_in_pool`、`post_process_hits_and_write`（`run.rs:4531`、全 query を 1 度に書く） | 最初の batch が無効な query だけで、後の batch（その batch の context がすべて有効）が表の失敗に当たったら、NCBI と同じく前の batch の報告を書いてから NCBI の誤り（終了コード 3）を出す。前の batch の報告では、outfmt 7 の終わりの行と outfmt 0 の epilog を書かない（NCBI の例外の経路で確かめる）。失敗する batch に無効な query がある分岐（NCBI が落ちる）は拒否のまま | 得点の表に無い得点で、最初の batch が無効な query だけ、後の batch が有効な query の入力（outfmt 0・6・7） |
| T9 batch ごとの警告 | `fasta.cpp:1617-1679`（title の警告は batch を読むとき）、無効な query の警告（batch の結果とともに） | `input.rs:180` `write_title_warnings`、`run.rs:5416` | batch ごとに、NCBI と同じ順で警告を書く（T8 の batch ごとの書き出しと合わせる） | 2 つ以上の batch で、title の警告と無効な query の両方がある入力（stderr を比べる） |
| T10 `-subject` が無い | `blast_args.cpp:2560-2562`（`CBlastDatabaseArgs::ExtractAlgorithmOptions`） | `args.rs:21`（clap が必須にする）、`cli.rs` | NCBI の文言と終了コード 1 を、NCBI と同じ時点（`-query`・`-out` を開く前）で出す。アダプタの `describe`（clap から作る）が正しいままか確かめる | `-subject` の無いコマンド（`-out` があるとき・無いとき） |
| T11 `CTOOLKIT_COMPATIBLE` | `showdefline.cpp:80`（`kBits` の静的な初期化で `getenv`） | `report/pairwise.rs:748,826` | 変数があれば「(bits)」。全 program の outfmt 0 の説明の一覧が共有する書き手で直す | 変数を与えた outfmt 0（BLASTN と、共有するほかの program） |
| T12 `PRE_FETCH_SEQS_LIMIT` | `blast_app_util.cpp:731-751`（`s_PreFetchSeqs`） | 無し | 整数でない値は、NCBI と同じ時点で NCBI の誤り（終了コード 255）を出す。整数は prefetch を切り替えるだけで、出力は変わらない | 変数を与えたコマンド |
| T13 `%#8.3g` | `align_format_util.cpp:578-620`（`PrintKAParameters`） | `report/pairwise.rs:957` `format_ncbi_ka_value` | C の `%#.3g` と同じにする：3 桁に丸めた後の指数で書き方を選び、指数の書き方は符号つき 2 桁以上、`#` で小数点と末尾の 0 を残す。全 program の outfmt 0 が共有する | 単体試験（Python の `'%#.3g' % x` は C と同じ）。0.99996・9.9996・999.5・1.23e-5 など |
| T14 outfmt 0 の書き込みの失敗 | `blast_app_util.hpp:242-254`（`CIOException::eFlush` は「BLAST failed to write output: <msg>」、`std::ios::failure` は「BLAST failed to write output」、どちらも終了コード 6 `BLAST_OUTPUT_ERROR`） | BLASTN の報告の書き出し（`run.rs` の `writer.flush()?` と `cli.rs` の `NativeError`） | outfmt 0 の書き込みの失敗で、NCBI の文言と終了コード 6 を出す（NCBI がどちらの例外になるかを oracle の `-out /dev/full` で確かめる）。outfmt 6/7 は承認済みの例外（`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 3）で、今の誤りの報告のまま | `-out /dev/full` と、読み手が閉じたパイプへの stdout（outfmt 0） |
| R1 `BL2SEQ_LEGACY` | `blast_app_util.cpp:206-211`、`local_blast.cpp:189,289`、`blast_format.cpp`（`m_IsBl2Seq && !m_IsDbScan`） | 無し | 明示的に拒否する（別の検索の mode と報告。値が空でも `getenv` は真）。今の環境変数の検査（`run.rs` ~5462）と同じ所で | 変数を与えたコマンド |
| R2 NCBI の診断と registry の変数、`.ncbirc` | `ncbiapp.cpp:1050-1135`、`ncbidiag.cpp`、registry の読み込み（`.ncbirc` を探す場所と順序） | 無し | blastn の stdout・stderr を変える環境変数（`DIAG_POST_LEVEL`、`DIAG_POST_PREFIX`、`NCBI_CONFIG__*` など）と、NCBI が読む `.ncbirc` の出力を変えるキー（`[BLAST] LONG_SEQID` など）を NCBI のソースから挙げ、oracle で確かめ、明示的に拒否する。出力に効かないキーだけの `.ncbirc`（`BLASTDB` など）は受け付ける（DW-13） | 変数ごと、キーごとの oracle の実行 |
| R3 UTF-8 でないファイル名 | `GetSubjectFile`（`-subject` の文字列をそのまま出す）、NCBI はファイル名を生のバイトで扱う | `report` の `Database:` の行、`cli.rs` の `inaccessible`（U+FFFD にする） | `-query`・`-subject`・`-out` のパスが UTF-8 でなければ、「not supported by LOSAT」で明示的に拒否する（DW-13）。拒否の時点は NCBI がそのファイルを扱う時点に合わせる | UTF-8 でないパスの入力 |
| V1 得点の幅 3000 | `blast_stat.c:2698-2734` | `scoring.rs:174`（`MAX_SCORE_RANGE`）、`scoring.rs:125`（`check_losat_limits`） | 制限を外した LOSAT を NCBI と比べる（reward は 32767 まで、penalty は −32768 まで、表の最大公約数の形を含む。出力と時間・メモリ）。一致すれば拒否を外し、しなければ見た理由を書いて残す | 得点の sweep の大きい値 |

**順序**
1. 小さい項目を 1 つずつ、それぞれ別のコミットで：T3、T13、T11、T2、T1、T6、T5、T4。
2. 次に T7、R1、R2、R3、T12、T10、T14。
3. 次に T8 と T9（batch ごとの書き出し）。
4. 最後に V1。

**各項目の後に行うこと**
- 直上に NCBI のファイル・行と逐語の断片を書き、`~/.cache/losat-web-gui-target/s07p-resume/verify_refs.py <変えた .rs>` を通す。
- `cargo fmt --check`、`clippy -D warnings`、`cargo test --all-features` を通す。
- `LOSAT/tests/ci_fast_regressions.py --losat <bin> --out <dir>`（capture の S02 の基準、凍結ハッシュ、outfmt 0 の fixture、BLASTN の fixture）を通す。
- その分岐の入力を NCBI と比べる（sonnet の agent。複数を並行してよい）。
- 移植した分岐の case を `LOSAT/tests/blastn_regression_fixtures.py` の `CASES` に足して `freeze` し（NCBI 2.17.0、`/home/kawato/micromamba/bin/blastn`）、PR の CI に守らせる。`freeze` は全件を作り直すので、既存の出力が変わらないことを確かめる。
- `build_inventory.py` の `RESULTS` に結果（例：`faithful after <commit>; <証拠>`）を書いて `INVENTORY.tsv` を作り直す。
- capture の case の出力が変わったら、それが NCBI との差を直した結果であることを NCBI との比較で示して記録する。基準は書き換えず、記録と対応させる。

注釈を足すと `inventory_refs.py` の固定した数（`EXPECTED`：54 ファイル、1829 の注釈、390 の関数）が変わり、最初の検査が失敗する。transpile の最後に新しい数を測り、`EXPECTED` と記録を同じコミットで直す。

### 2. 試験とゲート（S07+++ の指示書の 4.）

- **ゲートの script：** `~/.cache/losat-web-gui-target/s07p-resume/s07pp_gates2.sh` を写して E2g 用にする。
  - run は `docs/evidence/losat_web_e2g/run-<UTC>`。
  - 変更前は S07++ の成果物（`s07pg-*`。native の SHA-256 は `a4ff5abb3a61…80f5`。エンジンの検索のコードは `7693a9c73` から S07+++ の終わりまで変わっていない。`2fabad16c` は試験だけ）。
  - 実行の前に `date -u +%Y%m%dT%H%M%SZ > …/s07p-resume/ts2.txt`。
- **退行が無いことを確かめる検査：**
  - S07+ と S07++ の全検査（`check_inputs.py`、`scoring_sweep.py`、`word_size_sweep.py`、`slice_sweep.py`、`title_sweep.py`、`batch_sweep.py`、`split_check.py`）
  - Gate A と capture（236 件、S02 の基準）
  - outfmt 0 の fixture（1・2・4 スレッド）
  - v1 の WASI の検査、V-ABI（full と quick）
- **追加の比較：** T7 の後は `BATCH_SIZE`・`CHUNK_SIZE` を与えた比較も行う。
- **V-PERF：**
  - 変更前は S07++ の 3 つの成果物。`perf_cases.py run … --repeat 3` で、case は E2f と同じにする。
  - 閾値を超えた case は `--repeat 5` で測り直す。
  - V-PERF の段階はアプリ側の lock（`vperf_lock.sh`）を取る。

### 3. 独立監査

読み取り専用の監査（`ncbi_parity_auditor` の役割）を、`INVENTORY.tsv` を基準に行う。

- **確かめる観点：**
  - (a) 経路の網羅：経路にあるのに表に無い NCBI の関数は無いか。
  - (b) 移植の忠実さ：T・R・V の各項目。
  - (c) 拒否のまま残した理由の妥当さ。
- **入力：**
  - `docs/evidence/losat_web_e2f/audit_round2/{gen,cases,driver}.py` を使い回す。
  - 別の生物の配列（`LOSAT/tests/fasta` のウイルスなど）、IUPAC の文字の多い配列を含める。
  - T7 の後は `BATCH_SIZE` を与えた比較も含める。
- **進め方：** 観点ごとに sonnet の agent を並行させる（作業ディレクトリと出力のファイルを分ける）。
- 結論は supported / unsupported / inconclusive。unsupported なら直して、もう 1 回受ける。

### 4. 記録と引き継ぎ

README の規則 8 に従う。

- **ゲート記録（`docs/evidence/losat_web_e2g/README.md`）：**
  - 完了条件の表、transpile の結果、ゲート、V-PERF、独立監査
  - `evidence.sha256`
  - CI の残りの結果
- **計画と README：** 計画の状態と README の表を直す。
- **`main` への PR：** 0.3 のとおり。
- **次のセッション：** [S08 — TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md)。0.1 が TBLASTX の速度の退行なら、それを S08 の最初の作業にする。完了条件が残れば、続きを S07+++c にする。

## S07+++ からの引き継ぎ（2026-10-02 の実測）

### 済んだこと

- **CI（`db67e02cc`、`eb262bac2`）：**
  - PR の速い検査：`rust` 約 4〜5 分、fast output regressions 約 4.5 分（cache なし）。
  - 許可リスト：`LOSAT/tests/frozen_mismatch_allowlist.json` に、`Sakai.MG1655.megablast` の 1 項目。
  - `nightly.yml` を新設した。
  - 検出は PR #109 の run `36881993150` で確かめた（PR は閉じ、ブランチは消した）。
- **Sakai の差の中身（`4bb3be2a6`）：**
  - 凍結した `3f27c1f1…` は v0.1.0 の LOSAT の出力。
  - 今の LOSAT の `ac8177764f…` は NCBI 2.17.0 とバイト一致する。
  - 6482 行のうち 5 行が違い、違うのは得点の等しい greedy な traceback の経路で、query の分割とは関係しない。
- **棚卸しの第 1 段（`c6c3386cf`）：**
  - `inventory_refs.py`：54 ファイル、1829 の注釈、390 の関数。
  - 前回の 2180 と 669 は、数え方が失われていて再現できなかった。ファイルの数 54 は再現した。
- **第 2 段（`0ce717faa`〜`2dc7275b4`）と `INVENTORY.tsv`（`d1f96b309`）：** 内容はゲート記録の表。
- **予備の hit list の単体試験（`2fabad16c`）：** 試験だけ。

### 道具

- **`LOSAT/tests/ci_fast_regressions.py`**
  - PR の速い検査。手元では `--losat <bin> --out <dir> [--programs …] [--all-cases]`。
  - Gate A の字句のパス `/tmp/losat-pr5-runtime-cert-5845d22/LOSAT/tests/fasta` は自動で用意する。`capture_outputs.py` を直接使うときは、このパスが要る。
- **`LOSAT/tests/blastn_regression_fixtures.py`：** `generate`、`freeze --oracle`、`check --losat`。入力はコミットしてある。
- **`docs/evidence/losat_web_e2g/build_inventory.py`：** 規則と `RESULTS` から `INVENTORY.tsv` を作る。判断の無い divergent・unported の行があれば失敗する。

### 進め方（保守者の指示）

- **オーケストレータ：** Opus（メインのセッション）。設計を伴う作業だけを行う：移植の設計と実装、監査の指摘の判断、agent の成果の確認。機械的な作業は Agent ツールに `model: "sonnet"` で回す。
- **並行：** Sonnet の agent は、互いに独立な作業なら 4 つを超えて並行してよい（保守者、2026-10-02）。エンジンのソースの変更、ゲートの実行、V-PERF は 1 つずつ行う。agent ごとに出力のファイルと `--target-dir` を分ける。
- **使用量の上限：** agent には結果を途中でファイルに書かせる。上限に達したら並行数を 2 に減らす。
- **アプリ側：** S09（`LOSAT-web-gui-app` の未コミットの作業）には触らない。

### 保守者の決定（2026-10-02、計画 DW-13、`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`、AGENTS.md）

S07+++ の終わりに、判断待ちだった 8 件を保守者が決めた（どれも推奨の案）。

- **承認済みの例外（全 program）：**
  - 引数の構文の誤りと `-help`：LOSAT の文言と終了コード 2。引数の解析の後の誤りは NCBI と同じにする。
  - `-subject` での `-num_threads`：LOSAT は並列のまま、NCBI の 2 つのスレッドの警告を出さない。
  - outfmt 6/7 の書き込みの失敗：NCBI は abort、LOSAT は誤りを報告して 0 でない終了コード。
  - メモリ不足：LOSAT は abort。
  - `INVENTORY.tsv` の X-cli・X-threads・X-oom。
- **このセッションで移植：** outfmt 0 の書き込みの失敗（T14）。
- **このセッションで明示的に拒否：** UTF-8 でないファイル名（R3）、出力を変える `.ncbirc` のキー（R2）。
- **FASTA の読み方（TD-12）：** S17 の前の専用のセッション SF で移植する。このセッションでは拒否のまま。
- **E1a〜E1c の V-PERF と E1c の CLI の 2 つの差：** 記録どおり承認した（閉じた）。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S08 — TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md)。S08+ は同じ棚卸しの方式で行う。
