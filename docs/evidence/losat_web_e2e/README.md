# LOSAT Web E2e（Session S08+・S08+a・S08+b）ゲート記録

- 段階：E2e BLASTP・TBLASTN・TBLASTX の既定以外の option（[総合計画書](../../losat_web_gui_plan.md) §7 の S08+、計画 TD-13・DW-12・DW-15・DW-16・DW-17。指示書 [S08+](../../losat_web_gui_sessions/session_s08p_e2e_protein_options.md)、[S08+a](../../losat_web_gui_sessions/session_s08pa_e2e_audit_open_items.md)、[S08+b](../../losat_web_gui_sessions/session_s08pb_e2e_close.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S08+ の開始の `78c06fe61`（エンジンは SD の最後の `eeeea4fb2`。native の SHA-256 `487ac387…`、SD の最後の成果物 `~/.cache/losat-web-gui-target/sd-final-*`）。変更後は最後のエンジンのコミット `75cc8e565`（最後のゲートは記録と script を足した `8c905ccc5` で、native `4d6c036f…`）。S08+ は 2026-10-03〜04、S08+a は S08+ の第 1 回の監査の残件を別のブランチ `feature/losat-web-gui-s08pa` で並行して直し（2026-10-04）、S08+b が merge して閉じた（2026-10-04）
- 状態：（最後に書く）

## 完了条件（計画 §7 の S08+ の行）

| 条件 | 状態 | 根拠 |
|---|---|---|
| 各 program の sweep の全組合せが、NCBI と同じ拒否、outfmt 0/6/7 のバイト一致、明示的な拒否のどれか | 満たした（1 件は保守者の判断待ち） | 最後のゲートの sweep：BLASTP 1196・TBLASTN 1427・TBLASTX 641 で DIFF 0・timeout 0（変更前は DIFF 987・1073・194）。TBLASTN の既定以外の `-db_gencode` は承認済みの例外で、C++ API の oracle と一致。2 回の独立監査の約 13 万の比較でも、一致、同じ誤り、承認済みの例外、明示的な拒否だけ。保守者の判断待ちは D11〜D14（明示的な拒否）と D15（NCBI が落ちる検索で LOSAT は NCBI の valgrind の下の出力を出す。推奨は承認済みの例外。断られれば明示的な拒否にする）（下の「保守者に諮ること」） |
| 直した組合せの fixture で NCBI とバイト一致 | 満たした | outfmt 0/7 の manifest の E2e の 56 行（1・2・4 スレッドで差 0）、TBLASTX の `culling.*` 11 件。変更前の実行ファイルは足した行で違う（S08+ の [`audit/ROUND1.md`](audit/ROUND1.md) の確かめ、S08+a の `s08pa/checks/fixtures-base-e2e.log`、S08+b の `s08pb/r2a2_1/check-losat-old-n1.tsv`） |
| 各 program の既存のゲートと S07・S08 の fixture に退行なし | 満たした（Gate A は下の「Gate A」） | capture 236 件が S02 の基準と差 0（S08+ が NCBI とバイト一致にした BLASTP の 11 行の hash を更新）、CI の速い検査の全件 236 件で失敗 0（Gate A と TLOSAN Stage G の凍結ハッシュを含む）、TBLASTX 84・BLASTN 183 の回帰 fixture、題の sweep、BLASTN の入力、`CTOOLKIT_COMPATIBLE`、v1 の WASI の行列で形式の失敗 0 |
| V-PERF の非退行 | （V-PERF の後に書く） | 下の「V-PERF」 |
| 変えた program の V-ABI | 満たした | V-ABI full 800 の実行で失敗した部分 0（E2e の fixture の 50 の検索を含む）、凍結ハッシュ 1080/1084（違う 4 件は既知の Sakai）、V-ABI quick 60 の実行と凍結ハッシュ 16/16。`docs/web/verification_cells.tsv` に E2e の 4 行を `checked` で足した |
| 独立監査 | 満たした | 第 1 回（S08+、5 観点）の指摘を直し、第 2 回（S08+b、4 観点）で 4 観点とも supported（(a) は 3 回目） |

## NCBI の経路の記録

[`AUTHORITY.md`](AUTHORITY.md)。§A 3 つの program に共通の app の流れ（`CBlastAppArgs`、`NativeError` の「BLAST query/options error:」と終了コード、標準入力、`-out`）、§B 引数の読み方（`CArg_Integer`・`CArg_Double`・真偽値、`-outfmt` の書き方）、§C task の既定値（`blastp-fast`・`tblastn-fast` の threshold 20、word size による lookup の切り替え）、§D 引数の組の処理の順、§E `BLAST_ValidateOptions`、§F lookup の threshold と `(Int4)` の変換（x86-64 の `cvttsd2si` で `INT_MIN`）、§G hit list の大きさ（予備の `int` の回り込み）、§H 報告（500 の説明・250 の整列、蛋白の題、Method の文、epilog の `%g`）、§I 環境変数、§J 入力、§K LOSAT が対応する値と拒否する値（S12 の検索画面の入力）、§L TBLASTX の culling、§M 承認済みの例外と判断、§N 残り。ソースの引用は固定 commit（598d8ae6）のファイルと行で、コードの中の引用は `verify_refs.py`（このセッションで足した行の誤り 0、下の「ゲート」）で確かめた。

## 棚卸し（DW-12）

[`INVENTORY.tsv`](INVENTORY.tsv)（1063 行）。NCBI の 3 つの program の経路を 13 の範囲（PA・NA・XA：BLASTP・TBLASTN・TBLASTX の app と引数、IP・IN：蛋白と核酸の入力、PV・NV：option の検査、PS・NS：得点と統計、LK：lookup と word finder、XC：TBLASTX の culling、RP・RN：BLASTP と TBLASTN の報告）に分け、読み取り専用の agent（sonnet）が移植の前に行を作った（`inventory/{範囲}.tsv`、規則は `inventory/COMMON.md`、範囲ごとの指示と覚え書きは `inventory/{範囲}_brief.md`・`_notes.md`）。列 `status_before` は S08+ の前の LOSAT（`78c06fe61`）。

| | faithful | divergent | unported | rejected | exception | n/a |
|---|---|---|---|---|---|---|
| 移植の前（`status_before`） | 259 | 366 | 270 | 49 | 78 | 41 |

棚卸しから作業の束（WP1 CLI の入口と `NativeError`、WP2 引数の値と検査の順、WP3 入力、WP4 得点の表と行列、WP5 lookup と word finder、WP6 composition の mode、WP7 ungapped、WP8 TBLASTX の追加の option、WP9 報告、WP10 明示的な拒否・アダプタ・記録）を作り、移せる行を移し、残りを明示的に拒否した（§K）。行ごとの移植の後の状態は付けていない。移植の後の振る舞いは、下の sweep（3 つの program の全組合せ）と 2 回の独立監査（約 4 万の比較）で確かめた。

## 移植（S08+）

| コミット | 内容 |
|---|---|
| `b4bd7768e` | 最初の作業 1（DW-17）：TBLASTX・TBLASTN の outfmt 0 の句読点だけの subject の題を BLASTN と同じに書く（承認済みの例外 2）。[`title_sweep.py`](title_sweep.py) |
| `08396413b`・`391d8c383` | TBLASTX の `-culling_limit`（判断 D1、§L：`hspfilter_culling.c` の writer と pipe）。アダプタの `validate` |
| `4a473b228` | BLASTP・TBLASTN の最短の subject の下限（`BLAST_SEQSRC_MINLENGTH`）、TBLASTN の cbs 0 の最終の X-drop（`MAX`） |
| `68469680b` | BLASTP の outfmt 0 の Method の文を HSP ごとの行列の調整から |
| `8bd76f731`・`88430d947` | BLASTP の NCBI の app の層と option の意味（WP1・WP2）：`blastinput/app.rs`、`stats/protein_options.rs`・`protein_tables.rs`（[`gen_protein_tables.py`](gen_protein_tables.py) が NCBI の表から作る）、NCBI の順の `check_options`（アダプタの `validate` も使う）、`-outfmt` の field、`is_unported_blastp_arg` |
| `80e248e40`・`5c14ac1ee` | TBLASTN の app の層、option の意味、subject の読み方（最初の作業 2：TBLASTX と同じ `blastn/input.rs` の部品、判断 D2）。アダプタ |
| `7569e1f21` | TBLASTN：得点 0 の HSP、word size による lookup の切り替え、`Int4` の X-drop |
| `2271bedb0` | TBLASTX の NCBI の引数の文法と読んだ後の検査 |
| `1ddccbbec` | TBLASTX：NaN の e-value は全ての HSP を残す |
| `80e3188ae`・`6e2e0a15c`・`9ad9a3fbf`・`b2ef1e5cd` | 報告と入力：BLASTP・TBLASTN の outfmt 0 の既定の 500 の説明・250 の整列、BLASTP の蛋白の題（`x_InitDeflineTable`、50 文字の題の警告）、TBLASTN の Sbjct の B・Z・J の再翻訳、query の警告・無効な query・残基の無いレコード（O を X に読む警告を含む） |
| `c17856738`・`6ddc06327`・`2389859e9`・`a81abb584`・`2d26499ff` | NCBI BLAST+ 2.17.0 で凍結した fixture（下の「fixture」） |
| `5ac1dbd85` | 記録：`AUTHORITY.md`、棚卸し、変更前の sweep、ゲートの script |
| `f7ce8199f` | S12 に option の値、SX に BLASTX の項目を渡した |
| `b1f164529`・`d96412265` | TBLASTN の既定以外の `-db_gencode` の確かめ：NCBI の C++ API の oracle（[`gencode_api/`](gencode_api/)、[`gencode_api_check.py`](gencode_api_check.py)、[`gates/build_api_oracle.sh`](gates/build_api_oracle.sh)）。sweep の分類 |
| `a2b43ee12`・`3a3e4d1a8` | 第 1 回の独立監査の指摘の修正と記録（[`audit/ROUND1.md`](audit/ROUND1.md)、判断 D11・D12） |

## S08+a（第 1 回の監査の残件、[`s08pa/NOTES.md`](s08pa/NOTES.md)）

S08+ の第 1 回の監査で残った 5 件を、S08+ の閉じる作業と並行して別のブランチで直した（起点 `d96412265`、S08+b が `e1d4f98a0` で merge）。

| コミット | 内容 |
|---|---|
| `36572eace` | TN-2：traceback の入口で HSP を得点の順に並べ直さない（予備の hit list が heap になると NCBI は e-value の順のまま traceback に渡し、最後に翻訳する HSP が `stat_length` を決める。blast_hits.c:3266-3284、blast_traceback.c:358-365） |
| `20856865e` | TN-4：gapped DP は NCBI が触れる cell だけを書く（`dp_mem_alloc` は NCBI の値、実際の長さは `MIN(dp_mem_alloc, N + 1)`）。TBLASTN の `-xdrop_gap_final 1e8` が 26.9 秒・3.9 GB から 0.12 秒・8.4 MB（NCBI 0.23 秒・68 MB）。BLASTP・TBLASTN・BLASTX の共有の経路、出力は変わらない |
| `5c5132987` | TN-5：TBLASTN の NCBI の `qsort` に当たる並べ替えを安定に（`ScoreCompareHSPs` は frame を比べず、glibc 2.39 の `qsort` は安定な merge sort） |
| `7d77ee408` → `44ccde85c`・`61d3bf4f5` | RP-4：原因は NCBI の query の分割（split_query_aux_priv.cpp:73-138、chunk は tblastn 20000・blastp 10000、重なり 100）。最初の拒否の版を、保守者の指示で TBLASTN と BLASTP の移植に替えた（判断 D13） |
| `6fa06cb6a` | `-out -version`（BP-8 の残り）：値の位置の toolkit の語も明示的に拒否（判断 D14） |
| `68154c73e` | TN-1 の残り（監査の再現で見つかった）：窓の右端の traceback の `N`、`N <= 0` で DP を回さない、窓の後ろの番兵（blast_traceback.c:463-464,509-512、blast_gapalign.c:429-432,563-578） |
| `f3c8d978d`・`b9b2a6f06` | rustfmt と参照の行 |
| `c394f1381`・`c6fb13915` | 記録（[`s08pa/`](s08pa/)） |

確かめ（`68154c73e`、[`s08pa/checks/`](s08pa/checks/)）：fixture 147 件が 1・2・4 スレッドで差 0（変更前は足した 10 件だけが違う）、TBLASTN・BLASTP の sweep は DIFF 0・timeout 0、`cargo test --all-features` 909 件、監査の再現（5 つの agent、[`s08pa/audit_rerun/`](s08pa/audit_rerun/)）で説明の付かない差は TN-1 の残りだけ（直した）。native の性能は上限 ×1.05 の内（BLASTP は 2 回とも約 1〜2 % 遅い側）。事故：BLASTP の再現の agent が `pkill -f` と広い `kill` を使い、TBLASTX の再現の agent の 1 件と `bash -c` の 2 件を止めた（S08+ のゲートのプロセスは見当たらなかった）。

## S08+b（このセッション）

| コミット | 内容 |
|---|---|
| `76b61bf08` | 指示書と README の表 |
| `e1d4f98a0` | S08+a の merge（衝突なし。merge の後の native は S08+a の `LOSAT-final` と同じ `6f5e268d…`、fixture 147 件 × 1・2・4 スレッドで差 0、`cargo test` 909 件。[`s08pb/merge_check/`](s08pb/merge_check/)） |
| `dbc14613a` | 申し送り 2.1：アダプタの試験の BLASTP の登録を tab の無い定義行に（`store.rs`。S08+ の IN-10 から失敗していた） |
| `642b2669e` | 申し送り 2.2：`cli.rs` の 2 つの NCBI の参照の範囲 |
| `1f3a8fc18` | 申し送り 2.3：BLASTP の NCBI の `qsort` に当たる 6 つの並べ替え（`hsp.rs` の得点・e-value・purge の 2 つ、`blast_engine.rs` の chaining と得点）を安定に。`hsp.rs` の hit list の並べ替えは OID で決まるので変えない。変更の前後で 340 の比較が同じ（[`s08pb/stable_sorts/`](s08pb/stable_sorts/)）。BLASTX と共有の `redo_alignment.rs` は SX |
| `491292327` | ゲートの sweep の timeout 1200 秒と Gate A の script |
| `0e59e2f64` | 第 2 回の監査 R2A-1・R2A-2：BLASTP の compressed の lookup の短い subject と one-hit の負の長さで abort していた（[`audit/ROUND2.md`](audit/ROUND2.md)）。fixture 2 件 |
| `13493774f` | R2c-2：最後の `--` は何もしない（ncbiargs.cpp:2866-2872） |
| `efa445b19` | web ABI v1 の BLASTP の誤りの順と tabular の field（計画 TD-1、下の「ゲート」） |
| `141e49567` | R2b-1：`-evalue` が DBL_MAX 以上を拒否（判断 D12 を広げた） |
| `75cc8e565` | R2D-1：注釈の行 |
| `7f64c9cba`・`3373ae5fe`・`8c905ccc5` | 記録、最後のゲートの script（[`gates/s08pb_gates.sh`](gates/s08pb_gates.sh)、自分の build の directory、TBLASTX の sweep の並列 3）、監査 (a) の 2 回目の指示 |
| `3b84ce6c4` | 第 2 回の監査 (a) の 2 回目の R2A2-1：`-task blastp-fast` の chaining は ungapped の整列が 1 つの list にも行う（NCBI が飛ばすのは空の list だけ、blast_gapalign.c:3707-3708、drop の検査 3628-3636）。fixture 2 件 |
| `9e07a403e`・`875bfe114`・`70a1bb5a3`・`3bc4d3454` | 記録、ゲートの script の build の directory の接頭辞（`GATE`。前のゲートの Gate A が動いている間に次のゲートを組むため） |

## fixture

| 集まり | 件数 | 内容 |
|---|---|---|
| `LOSAT/tests/outfmt0_manifest.tsv` の `e2e.*`・`method.blastp`・`punct.tblastn`・`punct.tblastx` | 56（S08+ 42、S08+a 10、S08+b 4） | NCBI の outfmt 0 か 7 を凍結（`run_oracle.py`、stderr も固定した行がある）。BLASTP：`blastp-fast`、word 5、threshold の実数と `+inf`、`-seg` の小文字と窓、行列の名の大小、`-max_target_seqs` 3・280・回り込み、全て X、O の同一性と警告の順、無効な query、空の subject、300 の subject、蛋白の題、Method の文、query の分割、compressed の短い subject、one-hit の負の長さ、blastp-fast の 1 つの ungapped の整列（outfmt 0 と 7）。TBLASTN：BLOSUM45、同じ得点の frame、hard mask と少ない `-max_target_seqs`、fence と X-drop、最終と予備の X-drop、IUPAC の subject、真偽値の flag、O の警告、300 の subject、query の分割、窓の右端。TBLASTX：threshold の実数、句読点の題 |
| `LOSAT/tests/tblastx_regression_fixtures.py` の `culling.*` | 11 | TBLASTX の `-culling_limit`（§L） |

fixture の数は 151（SD の 95 に 56）。凍結と確かめは各ゲートの `check-losat-n{1,2,4}.tsv`・`oracle-check-gate.log`。

## sweep

[`option_sweep.py`](option_sweep.py) が、3 つの program の option の値と組（行列・gap、word size、threshold、window、`-comp_based_stats`、`-seg`、`-evalue` の書き方、`-max_target_seqs`、遺伝暗号、`-outfmt`、真偽値、X-drop、culling など）を [`make_inputs.py`](make_inputs.py) の入力で、NCBI BLAST+ 2.17.0 と outfmt 0・6・7 で比べる（stdout、stderr、終了コード）。TBLASTN の既定以外の `-db_gencode` は承認済みの例外で、C++ API の oracle で確かめる（[`gencode_api/check.tsv`](gencode_api/check.tsv)）。蛋白の題は [`protein_title_sweep.py`](protein_title_sweep.py)、句読点の題は [`title_sweep.py`](title_sweep.py)。

| sweep | 件数 | 変更前（[`sweeps/before-*.tsv`](sweeps/)、`78c06fe61`） | 最後のゲート |
|---|---|---|---|
| BLASTP | 1196 | DIFF 987、一致 140、引数の構文の誤り 63、LOSAT の拒否 6 | DIFF 0・timeout 0：一致 231、同じ誤り 269、LOSAT の拒否 633、引数の構文の誤り 63 |
| TBLASTN | 1427 | DIFF 1073、一致 192、引数の構文の誤り 111、拒否 6、`-db_gencode` の例外 45 | DIFF 0・timeout 0：一致 258、同じ誤り 269、拒否 720、引数の構文の誤り 102、`-db_gencode` の例外 78（C++ API の oracle で 0 failures） |
| TBLASTX | 641 | DIFF 194、一致 219、引数の構文の誤り 123、拒否 24、同じ誤り 3、例外 75、timeout 3 | DIFF 0・timeout 0：一致 270、同じ誤り 83、拒否 90、引数の構文の誤り 123、`-db_gencode` の例外 75（`-db` の oracle） |

## 決めたこと

[`AUTHORITY.md`](AUTHORITY.md) §M。推奨の案で進め、記録した（保守者の常の指示）。

- D1 TBLASTX の `-culling_limit` を移植。D2 TBLASTN の subject は `blastn/input.rs` の部品で読む（`CFastaReader` の全体の移植は SF）。D3 NCBI が落ちる組は明示的な拒否。D4 `double` から `Int4` への未定義の変換は `INT_MIN`。D5 unified P は拒否。D6 実数は `CArg_Double` の書き方。D7 未知の行列は NCBI の engine error。D8 NCBI が終わらない window は拒否。D9 `-remote`・`-db` の一族・toolkit の option は拒否。D10 パイプの空の query は拒否。
- D11・D12（第 1 回の監査）、D13・D14（S08+a）は保守者の確認を待つ（下の「保守者に諮ること」）。

## ゲート

最後のゲートは `70a1bb5a3`（エンジンは最後のエンジンのコミット `3b84ce6c4`。後のコミットは記録と script だけ）：run [`run-20261004T163746Z/`](run-20261004T163746Z/)（`head.txt`）。script は [`gates/s08pb_gates.sh`](gates/s08pb_gates.sh)（S08+ の [`gates/s08p_gates.sh`](gates/s08p_gates.sh) を写し、sweep の timeout を 1200 秒、TBLASTX の sweep の並列を 3 にした。build の directory は `GATE=s08pb2-gate` の `s08pb2-gate-*`）、Gate A は [`gates/s08pb_gate_a.sh`](gates/s08pb_gate_a.sh)。成果物のハッシュは `artifacts.sha256`（native `6f070575…`、reactor `cb0734e3…`・`2080b5b4…`。第 2 回の監査 (a) の 3 回目の実行ファイルと同じ）。変更前は SD の最後の成果物（native `487ac387…`）。

その前の 2 つのゲートは、どちらも後でエンジンを直したので記録に使わない（`~/.cache/losat-web-gui-target/s08pb/superseded-run-*`）：`491292327` の run（監査の第 2 回の 1 回目の指摘の前）と `8c905ccc5` の run（R2A2-1 の前。全ての検査が通り、Gate A は EXACT_TEXT 14・`-db_gencode` の HSP_SET_DIFF 6）。最後のゲートの結果は、R2A2-1 の 2 つの fixture のほかは `8c905ccc5` の run と同じ（sweep の行の分類、capture の 236 の hash、v1 の記録）。

| 検査 | 結果 |
|---|---|
| `cargo fmt --check`（LOSAT、アダプタ）、clippy `-D warnings` の 4 構成とアダプタの 3 構成 | すべて終了コード 0 |
| `cargo test --all-features`、アダプタ、wasm32 の web API の試験 | 911 件通過・失敗 0（無視 3）、アダプタ 5 件と 2 件、web API 6 件（`blastp_web_pair_keeps_v1_error_order_and_fields` を含む） |
| pure-Rust の境界、`ci_fast_regressions.py` の単体試験、`protein_tables` の検査、このセッションで足した行の NCBI の参照（`verify_added.py`、BASE `78c06fe61`） | 通過、誤り 0（`verify-refs-session-added.txt`） |
| outfmt 0 の fixture（`check_losat.py`、1・2・4 スレッド） | 151 件、差 0（3 つのスレッド数とも。承認済みの `-db_gencode` の例外 2 行 `tblastx.code4.0`・`code4.7` は分類どおり） |
| `run_oracle.py`（NCBI の凍結の確かめ） | 151 件、差 0 |
| TBLASTX の回帰 fixture（`culling.*` 11 件を含む）、BLASTN の回帰 fixture | 84 件差 0、183 件差 0 |
| `CTOOLKIT_COMPATIBLE`（`ctoolkit_compare.py`）、句読点の定義行（`punct-defline.tsv`） | 147 行 441 実行で差 0、予期しない 0 |
| 句読点だけの subject の題（[`title_sweep.py`](title_sweep.py)） | TBLASTX・TBLASTN とも 1023 の定義行で一致 957、例外 2 が 66、予期しない 0 |
| BLASTN の入力（E2g の `check_inputs.py`） | 300 件、予期しない 0 |
| CI の速い検査の全件（`ci_fast_regressions.py --all-cases`。Gate A の凍結ハッシュと TLOSAN Stage G を含む） | 236 件、失敗 0、許可した既知の不一致 1（`Sakai.MG1655.megablast`） |
| sweep（[`option_sweep.py`](option_sweep.py)、`sweeps/after-*.tsv`） | BLASTP 1196・TBLASTN 1427・TBLASTX 641 で DIFF 0・timeout 0（上の「sweep」）。TBLASTX の `-threshold +inf` の outfmt 0・6・7 も一致（下の注） |
| TBLASTN の `-db_gencode` の C++ API の oracle（[`gencode_api_check.py`](gencode_api_check.py)）、蛋白の題（[`protein_title_sweep.py`](protein_title_sweep.py)） | failures 0、一致 1014・LOSAT の拒否 69 |
| capture（`capture_outputs.py`） | 236 件が S02 の基準と差 0（S08+ が BLASTP の 11 行の hash を NCBI とバイト一致の出力に更新した、[`s02_baseline_updates.tsv`](s02_baseline_updates.tsv)）。SD の最後の成果物とは BLASTP の 11 件（12 の hash）が違い、11 件とも stdout が NCBI とバイト一致、stderr は `pairwise_limit1` の 3 件が NCBI と同じ警告になり、4 スレッドの 1 件は NCBI のスレッドの警告だけ（承認済みの例外 1、SD と同じ）（[`s08pb/capture_vs_ncbi_final.tsv`](s08pb/capture_vs_ncbi_final.tsv)） |
| v1 の WASI の行列（`check_wasm_threading.py`） | 433 の記録、reactor の lifecycle の gate 通過、形式の失敗 0。`v1-requests` は E1d と 14 件中 10 件が同じで、4 件は誤りの文言だけが違う（下の注） |
| V-ABI quick | 60 の実行が native の CLI と一致、凍結ハッシュ 16/16 |
| V-ABI full | 200 の検索（BLASTN 77、TBLASTN 49、TBLASTX 38、BLASTP 36。E2e の fixture の検索は BLASTP 27・TBLASTN 21・TBLASTX 2 の 50）× 4（serial n1、threads n1・n2・n4）= 800 の実行、失敗した部分 0。各実行の 0/6/7 の stream が native の CLI と一致。凍結ハッシュ 1084 件中 1080 件一致、違う 4 件は既知の `Sakai.MG1655.megablast` outfmt 7 |

**TBLASTX の `-threshold +inf`（第 2 回の監査 R2c-1）。** 最初のゲート（`491292327`、記録には使わない）では、TBLASTX の sweep のこの 2 行で LOSAT が並列の実行の中で記憶不足で止まった（signal 9）。単独では outfmt 0・6・7 とも NCBI とバイト一致（LOSAT 8.6 GB・約 225 秒、NCBI 5.9 GB・約 175 秒、`s08pb/thrinf/`）。最後のゲートは TBLASTX の sweep の並列を 3 にし、3 行とも一致した。記憶の差は残件（AUTHORITY §N）。

**web ABI v1 の BLASTP（計画 TD-1）。** S08+ が BLASTP を NCBI の app の層に移したため、v1 の BLASTP の誤りの順と、書かない tabular の field の扱いが E1d の記録から変わっていた。`efa445b19` で v1 の経路だけ、option を先に確かめ、書かない field を以前の文言で拒否するようにした。残る差は 4 件の誤りの文言で（[`s08pb/v1_requests_vs_e1d.txt`](s08pb/v1_requests_vs_e1d.txt)：未知の行列、BLOSUM45、`-ungapped`）、エンジンの今の option の検査の文言（NCBI の文言か LOSAT の明示的な拒否）である。状態（−1）、出力の無さ、検査の順は E1d と同じ。TD-1 は v1 の引数・形式・誤りの扱いを凍結し、エンジンの修正が v1 に及ぶことと、共有の誤りの文言が変わること（S07+ のスレッドの上限）を認めているので、この扱いにした（保守者に諮る、下）。

## Gate A

（ゲートの後に書く）

## V-PERF

（ゲートの後に書く）

## 独立監査

### 第 1 回（S08+、`f7ce8199f`）

[`audit/ROUND1.md`](audit/ROUND1.md)。5 つの観点（BLASTP の app と option、TBLASTN、TBLASTX、BLASTP・TBLASTN の報告、入力）を sonnet の監査役 5 人が読み取り専用で並行して、約 7000 の argv で NCBI と比べた（報告は [`audit/round1/`](audit/round1/)）。指摘 53 件（重複を除くと 46 件）を、直した（`a2b43ee12`）、明示的な拒否（D11・D12）、受け入れ（承認済みの例外、v1 の凍結、両方が終了コード 1 の失敗の順）、S08+b に残した 5 件（TN-2・TN-4・TN-5・RP-4・`-out -version`）に分けた。残した 5 件は S08+a が直した（上）。

### 第 2 回（S08+b）

[`audit/ROUND2.md`](audit/ROUND2.md)。第 1 回の 5 つの報告を (a) BLASTP、(b) TBLASTN、(c) TBLASTX、(d) 報告と入力 の 4 観点にまとめ、sonnet の監査役 4 人が読み取り専用で並行して、第 1 回の全ての再現の命令と harness を繰り返し、S08+・S08+a・S08+b の変更の周り（query の分割、巨大な X-drop、BLOSUM45 の同点、hard mask、窓の右端、toolkit の語、BLASTP の同順位）を新しく調べた（指示は [`audit/round2/brief/`](audit/round2/brief/)、報告は [`audit/round2/`](audit/round2/)）。比べた数は約 12 万。

| 観点 | 結論 | 指摘と対応 |
|---|---|---|
| (a) BLASTP | 1 回目 unsupported（`491292327`）→ 2 回目 unsupported（`75cc8e565`）→ **3 回目 supported**（`3b84ce6c4`） | R2A-1 compressed の lookup の短い subject と R2A-2 one-hit の負の長さの abort（`0e59e2f64`、fixture 2 件）、R2A2-1 blastp-fast の 1 つの ungapped の整列の chaining（`3b84ce6c4`、fixture 2 件）。R2A2-2 は NCBI が配列の外を読んで落ちる検索で、LOSAT は NCBI の `valgrind -q` の下の出力と一致（D15、保守者の判断待ち） |
| (b) TBLASTN | supported（1 回目） | R2b-1 DBL_MAX の `-evalue` で NCBI が落ちる（D12 を広げた、`141e49567`）、R2b-2 小さい `CHUNK_SIZE` の速さ（残件） |
| (c) TBLASTX | supported（1 回目） | R2c-1 完全な近傍の記憶（残件、出力は同じ）、R2c-2 最後の `--`（`13493774f`）、R2c-3 文言（受け入れ） |
| (d) 報告と入力 | supported（1 回目） | R2D-1 注釈の行（`75cc8e565`） |

## 保守者に諮ること

推奨の案で進め、記録した。判断が出たら、`PD-LOSAT-NCBI-DEFECTS`・`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` と `AGENTS.md` の承認済みの例外の記述を、承認の範囲で更新する。出なければ今のまま残件として次のエンジン側のセッションに渡す。

1. **D11**（第 1 回の監査 TX-1）：query の長さ＋window が 2^31 − 1 を超えると NCBI の `Int4` が回り込み hit が無い（決まった結果）。実用が無く再現の費用が高いので、D8 と合わせて和が 2^30 を超える値を BLASTP・TBLASTX で明示的に拒否した。推奨：このまま。
2. **D12**（BP-4・TN-6・R2b-1）：無限大と DBL_MAX 以上の `-evalue` で NCBI の blastp・tblastn は入力によって SIGSEGV。LOSAT はいつ落ちるかを事前に知れないので明示的に拒否した（1 つ下の有限の値は実行し NCBI と一致）。推奨：このまま。別案：近い有限の値の NCBI の出力と一致を確かめた LOSAT の結果を承認済みの例外にする。
3. **D13**（S08+a、RP-4）：BLASTP・TBLASTN は NCBI と同じ chunk で query を分ける。NCBI が落ちる分割の設定（chunk をもう一度分けるほど大きい重なり、batch を分ける負の `CHUNK_SIZE`）は明示的な拒否（BLASTN は前者を承認済みの例外 1 で検索している）。推奨：このまま。
4. **D14**（S08+a、BP-8 の残り）：option の値の位置の `-version` などの toolkit の語も、option の位置と同じ文言で拒否する（NCBI は version を出して終了 0）。推奨：このまま。
5. **D15**（第 2 回の監査 R2A2-2）：BLASTP の one-hit と `-threshold` 1〜3 の一部で NCBI は配列の外を読んで SIGSEGV。LOSAT は配列の外を NULLB として読み、NCBI を `valgrind -q` の下で走らせた出力と 23 件ともバイト一致。推奨：承認済みの例外（NCBI の不具合の方針の「確かめられる妥当な結果」）。別案：窓が配列の外を読む検索の明示的な拒否。
6. **web ABI v1 の BLASTP の誤りの文言**（計画 TD-1）：4 つの要求の誤りの文言がエンジンの今の文言（NCBI の文言か LOSAT の明示的な拒否）に変わった。状態、出力の無さ、検査の順は E1d と同じ。推奨：TD-1 の S07+ のスレッドの上限の文言と同じ扱いで受け入れる。
7. V-PERF：（V-PERF の後に書く）

S08b の V-PERF の `tblastx-multi`（E2b のゲート記録）は今も保守者の確認待ち。

## アプリ側（S09・S12）への注意

- 検索画面の option の値は [`AUTHORITY.md`](AUTHORITY.md) §K（S12 の指示書の S08+ の項、S08+b で更新）：BLASTP・TBLASTN の長い query は NCBI と同じに分けて検索する。`-version` などの toolkit の語は値の位置でも拒否する。`-evalue` は DBL_MAX 以上を拒否する。最後の `--` は何もしない。
- ABI v2 の `validate` は BLASTP・TBLASTN・TBLASTX の NCBI の option の検査と LOSAT の上限を CLI と同じ順と文言で返す（S08+）。
- BLASTP・TBLASTN の outfmt 0 は NCBI の既定の 500 の説明と 250 の整列を出す（S08+）。300 の subject の検索で報告は大きくなる。
- web ABI v1 は凍結のまま（TD-1）。BLASTP の 4 つの誤りの文言だけが E1d と違う（上の「ゲート」の注）。
- verification の升目：`docs/web/verification_cells.tsv` に E2e の 4 行（fixture、TBLASTX の culling、sweep、V-ABI）を `checked` で足した。

## 残件と引き継ぎ

- **`main` への PR**：（PR の後に書く）
- **残件**（[`AUTHORITY.md`](AUTHORITY.md) §N。出力に関わるものは無い）：R2c-1 TBLASTX の完全な近傍の記憶（NCBI の最大 2.3 倍。NCBI は tblastx の query を 10002 塩基の chunk に分ける）、R2b-2 小さい `CHUNK_SIZE` の TBLASTN の速さ、TBLASTN の 300 の subject の速さ（前から）、BLASTP の分割される batch の全体の検索の費用、R2c-3 文言。
- **SX**（[指示書](../../losat_web_gui_sessions/session_sx_blastx_integration.md)の「S08+a・S08+b からの引き継ぎ」）：BLASTX と共有の gapped DP の確保と compressed の走査の修正、`read_end_sentinel`、`redo_alignment.rs` の `sort_unstable_by`、blastx の query の分割（重なり 297）、toolkit の語、最後の `--`、one-hit の負の長さ。
- **S12**（[指示書](../../losat_web_gui_sessions/session_s12_w3_search_ui.md)）：S08+b の 2 つの点を足した（長い query の分割、値の位置の toolkit の語）。
- **次のエンジン側のセッション**：[S11 — 領域の指定](../../losat_web_gui_sessions/session_s11_e2d_query_subject_loc.md)。保守者の判断（上）が出たら、その反映を最初の作業にする。
