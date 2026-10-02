# LOSAT Web E2b（Session S08）ゲート記録

- 段階：E2b TBLASTX の outfmt 0 と 7（[総合計画書](../../losat_web_gui_plan.md) §7 の S08、指示書 [S08](../../losat_web_gui_sessions/session_s08_e2b_tblastx_outfmt0_7.md)、計画 DW-6・DW-10・DW-12・DW-15・TD-1）
- ブランチ：`feature/losat-web-gui`。変更前はセッションの開始の `2bcb86b1f`（エンジンは E2g の最後の `f4057718a` と同じ。native の SHA-256 `331fcba36447…44d8`）。変更後はこの記録の時点で `1117e8c17`（最後のエンジンのコミット。全ゲートは S08b で実行する）
- 状態：**未完了（2026-10-03、S08b に続く）**。移植・fixture・棚卸し・独立監査の第 1・2 回とその指摘への対応は済んだ。最後のコミットでの全ゲート（lint と試験は済み）、独立監査の第 3 回（確認）、Gate A、V-PERF、`main` への PR は、指示書 [S08b](../../losat_web_gui_sessions/session_s08b_e2b_final_gates.md) で行う（WSL の `/mnt/c` の I/O の誤りでゲートが止まり、保守者の指示でここで区切った）。下の「完了条件」の表

## 完了条件（計画 §7 の S08 の行）

| 条件 | 状態 | 根拠 |
|---|---|---|
| 固定した fixture で NCBI とバイト一致（承認済みの遺伝暗号の例外を除く） | 最後のコミットでは未確認（S08b） | `7fbbfad96` のゲート（部分、`~/.cache/losat-web-gui-target/s08/gate-7fbbfad96-partial/`）：outfmt 0 の fixture 74 件（1・2・4 スレッド）、TBLASTX の回帰 fixture 69 件、BLASTN の fixture 107 件がすべて一致。その後の 70 件目（`cutoff.floor_evalue_1e10`）を含む 70 件は `aa4b6f8c3` の手元の build で一致。下の「ゲート」 |
| 既存の 6 に退行なし | 最後のコミットでは未確認（S08b） | `968fa98c8` と `7fbbfad96` のゲート：Gate A 以外の検査（capture 236 件が S02 の基準と差 0、速い検査の全件 236 件で失敗 0、v1 の WASI の行列）。Gate A（`audit_tblastx_v010.py`）は S08b |
| TBLASTX の全升目の V-ABI | 最後のコミットでは未確認（S08b） | 監査 (d) の第 1・2 回（serial と threaded の reactor、0・6・7 × スレッド 1・2・4、凍結の NCBI のハッシュ 24/24・72/72）。ゲートの V-ABI full は S08b |
| 独立監査 | 第 1・2 回の指摘に対応済み、第 3 回（確認）は S08b | 下の「独立監査」。第 2 回：(b)・(d) は supported、(a)・(c) は直す前の実行ファイルで unsupported（指摘はすべて直したか記録した） |

## NCBI の経路の記録（指示書の 1.）

[`AUTHORITY.md`](AUTHORITY.md)（英語）。§A アプリの経路（`tblastx_app.cpp` の `CTblastxApp::Run` から、subject の読み込み、query の batch、`CLocalBlast`、`CBlastFormat`）、§B 報告が検索から受け取るもの（batch、hit list、ncbi2na の乱数の塩基、sum statistics の `num`、HSP の並べ替え、linking の組、`BLAST_LargeGapSumE`、SEG）、§C outfmt 0（`showdefline.cpp` の説明の一覧と `N` の列、`showalign.cpp` の両側の frame と翻訳の行、`blast_format.cpp` の query の末尾と epilog）、§D outfmt 7 と 6（`tabular.cpp`）、§E 表示の翻訳と遺伝暗号（承認済みの例外の分類はデータベースのオラクル）、§F このセッションで決めたこと、§G 棚卸し。ソースの引用はすべて固定 commit（598d8ae6）のファイルと行で、引用の文字列が行と一致することを機械で確かめた。

## fixture（指示書の 2.）

| 集まり | 件数 | 内容 |
|---|---|---|
| `LOSAT/tests/outfmt0_manifest.tsv` の `tblastx.*` | 26 | NCBI BLAST+ 2.17.0 の outfmt 0 と 7 を凍結（`run_oracle.py`）。6 frame、複数の HSP、複数の query と subject（4 query・3 batch、260 subject）、ヒット無し、曖昧な文字、`-max_target_seqs` 255、遺伝暗号 1 と既定以外（`-query_gencode 4`、`-db_gencode 4`）。`tblastx.code4.*`（既定以外の `-db_gencode`）は承認済みの例外で、期待値は subject を `makeblastdb -dbtype nucl` にした NCBI の `-db` の出力（データベースの行と Posted date を正規化）。ほかの差は認めない |
| `seg.blastp.6`（同じ manifest） | 1 | `blastp -seg yes`（SEG の左の再帰。下の「調査」） |
| `LOSAT/tests/fixtures/tblastx_regression/`（`LOSAT/tests/tblastx_regression_fixtures.py`） | 70 | BLASTN の fixture と同じ方式（case ごとの環境変数、`.merged` は stderr を stdout に合わせる）。`BATCH_SIZE` 700・22000・100000・−1・0・整数でない値、検索されない batch の警告、`-max_target_seqs` 1・2・5・250 と既定、曖昧な文字（2・4 スレッド）、`-threshold`・`-window_size`・`-seg no`、`CTOOLKIT_COMPATIBLE`、害の無い `.ncbirc`、空の入力、RNA（`U`）、3 nt と 4 nt の query、同点の HSP（4 件）、`BLAST_LargeGapSumE`、SEG の左の再帰、`/dev/full`。独立監査の後に 17 件：`-evalue 0`（3）、最初の空行の subject と空の query、NCBI が decode しない題と表示されない subject の題（2）、DEL（2）、`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`（5）、0 以下の `-seg` の値（3）、`BLAST_Cutoffs` の下限。入力は 1 MB ほどでコミットした |

outfmt 0 の fixture の NCBI による凍結と確かめの記録は [`run-20261002T084101Z/`](run-20261002T084101Z/)（`oracle-freeze.log`、`oracle-check.log`）。

CI の速い検査（`LOSAT/tests/ci_fast_regressions.py`）は、TBLASTX を選んだとき、outfmt 0 の fixture の TBLASTX の行と上の 70 件を確かめる。TBLASTX は BLASTN の hit list と入力の読み方を使うので、BLASTN の変更でも TBLASTX を選ぶ。`ci.yml` と `nightly.yml` は `*-regression-fixtures.tsv` を成果物にする。

区別の確認：変更前の実行ファイル（`331fcba36447…`）は最初の 53 件のうち 50 件で違い（outfmt 0/7 の case は引数の誤り、outfmt 6 の case は batch・hit list・曖昧な文字の差）、棚卸しの結果の修正の前の `0d533ba76` は、その後に足した 13 件すべてで違う。監査の後に足した 17 件は、修正の前の実行ファイル（`968fa98c8`、`cutoff.*` は `7fbbfad96`）で、効かないことを示す 3 件を除く 14 件すべてが違う。環境変数を使う 21 件は、変数を外すと結果が変わる（7 件は、変数が効かないことを示すための case。`env_discrimination.py`）。

## 移植（指示書の 3.〜6.）

| コミット | 内容 |
|---|---|
| `fd896a52c` | outfmt 0 と 7 を 1 つの検索から書く。NCBI の経路を一括で移した（DW-12）：10002 nt の query の batch（`BATCH_SIZE` は NCBI と同じに読む）と batch ごとの検索（linking の cutoff も batch ごと）、有効な query の無い batch は検索しない（警告、Karlin block −1、「hits found」の行が無い）。`-max_target_seqs` の hit list（`Blast_HitListUpdate`、あふれたときの e-value の順、同点は subject の順）、説明 500・整列 250 の既定。予備の検索は NCBI の ncbi2na の subject（曖昧な文字は乱数の塩基）、再評価は ncbi4na の翻訳。sum statistics の `hsp->num`（`Expect(n)`、説明の一覧の `N` の列）。表示の翻訳（B・Z・J・X、query と db の遺伝暗号）、表示の文字からの identities の数え直し、query の行の SEG の小文字。outfmt 0 の prolog を最初の batch の前に書く、書き込みの失敗（「BLAST failed to write output」、終了コード 6）、「Query is Empty!」と空の subject の誤り。最終の HSP 一覧から `PairwiseHit` を作り、`hits` と outfmt 0 の formatter の両方に渡し、観測者を outfmt 0 と 7 につないだ。TBLASTN も E2a §G.3 の説明の一覧と NCBI の核酸の title を使う（凍結した出力は変わらない）。アダプタ：TBLASTX の形式 0/6/7、`parse` は `-outfmt` を挿入しない、`validate` は CLI の文言（`docs/web/abi_v2.md` §4・§5・§8） |
| `0d533ba76` | NCBI を凍結した TBLASTX の回帰 fixture（40 件）と CI の速い検査への組み込み |
| `06953c68d` | 棚卸しの結果の修正：入力を BLASTN と同じに読む（`blastn/input.rs` を program の名前つきにした。`-subject` を先に読み、query を開き、`-out` を作り、それから query を読む。`U` を `T`。CFastaReader の title の警告。NCBI が違う読み方をする定義行・文字・中身の無い record・UTF-8 でないファイル名・パイプからの空の query は明示的な拒否）。sum statistics の linking を query と鎖ごとに（NCBI は query ごとに 6 context）。subject の塊ごとの init hit list を `score_compare_match` で並べる（`aa_ungapped.c:234-235`）。TBLASTX の `BLAST_LargeGapSumE` を NCBI の評価の順に（BLASTX は SX まで前の関数、DW-10）。`-culling_limit` 1 以上の拒否。epilog の threshold を C++ の stream の `%g` で |
| `a46dd63fd` | RNA と短い query の fixture（6 件）。`run_oracle.py` はデータベースのオラクルの Posted date を固定の文字列にしてハッシュを取る |
| `2c9f71aa7` | 最初の `BLAST_LinkHsps` の後の `Blast_HSPListSortByScore`（`link_hsps.c:1802-1803`） |
| `2049103b0` | 同点の HSP と `BLAST_LargeGapSumE` の fixture（5 件） |
| `aee365eb5` | SEG の左の再帰は先頭の区間だけを残す（`blast_seg.c:2086-2101`）。BLASTP・TBLASTN・TBLASTX。BLASTX は `keeping_all_left_segments` で SX まで前のまま（DW-10） |
| `968fa98c8` | SEG の fixture（TBLASTX 2 件、BLASTP 1 件） |
| `7fbbfad96` | 独立監査の第 1 回の指摘への対応（下の「独立監査」）：outfmt 0 の題の検査を表示される subject だけに、NCBI の `HtmlDecode` と同じ判定に（BLASTN も共有）。`-evalue 0`。`-seg` の分け方と 0 以下の値（BLASTP・TBLASTN・TBLASTX）。`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`・`BL2SEQ_LEGACY`。`BATCH_SIZE` を query の前に読む。DEL。help の文。アダプタの `validate`。`abi_v2.md`。`check_wasm_threading.py` |
| `32fd67a52` | 注釈だけ（NCBI の行の範囲） |
| `aa4b6f8c3` | `BLAST_Cutoffs` の下限 1、`PRE_FETCH_SEQS_LIMIT`（監査 (b)） |
| `6017f0ea2` | アダプタの `register` の検査の順、`abi_v2.md` の文言（監査 (d) 第 2 回） |
| `5fa53b3f0` | `bio` が読めない subject の後回しを NCBI が黙って読むものだけに（監査 (a)・(c) 第 2 回） |
| `1117e8c17` | 閉じた標準出力の見分け方を外す（監査 (c) 第 2 回の N2。保守者に諮る）。`-culling_limit` の help |

## 棚卸し（DW-12）

[`INVENTORY.tsv`](INVENTORY.tsv)（447 行、[`build_inventory.py`](build_inventory.py) が作る）。NCBI の tblastx の outfmt 0/6/7 の経路を 7 つの範囲（A アプリと orchestration、B 整列、C 説明の一覧、D outfmt 7/6、E 報告が検索から受け取るもの、F 文字と遺伝暗号、G hit list）に分け、読み取り専用の agent（sonnet）が移植の前に行を作り（`inventory/{A..G}.tsv`、規則は `inventory/COMMON.md`。列 `status` は S08 の前の LOSAT）、別の agent が移植の後の `0d533ba76` で行ごとの結果を確かめた（`inventory/result_{A..G}.tsv`、規則は `inventory/RESULT_COMMON.md`。NCBI との比較の実行は合わせて 3000 を超える）。GAP とされた 15 行は直すか拒否した（`build_inventory.py` の `RESOLUTION`、列 `s08_final`）。

| | faithful | reused | ported | rejected | exception | deferred | n/a |
|---|---|---|---|---|---|---|---|
| 移植の前（`status_before`） | 64 | 92（reusable） | — | 8 | — | — | 88 |
| 最後（`s08_final`） | 63 | 64 | 195 | 23 | 6 | 8 | 88 |

移植の前の残り：divergent 29、missing 61、needs-param 80（行の数は 447。result の agent が足した 2 行と、独立監査が足した範囲 H の 23 行は移植の前の状態が無い）。`deferred` の 8 行は S08+ の範囲（NCBI の実数・整数の引数と `-outfmt` の文字列の読み方、`-subject` の欠落、NCBI のファイルの誤りの文言、`-query -`、引数の解析の後の検査の文言など）と、保守者に諮る閉じた標準出力。

## 棚卸しの結果の後の調査（4 件）

結果の確かめと、その後の広い比較で、LOSAT の報告に出ていなかった差が 4 つ見つかった。どれも S08 より古く、outfmt 6 にも出る。原因を NCBI のソースと計装（qsort の shim、`LD_PRELOAD` の SEG の shim、LOSAT の段階の dump）で特定し、直した（[`investigations/`](investigations/)。各 `.md` が記録、`.diff` が提案した修正）。

| 調査 | 原因 | 修正と確認 |
|---|---|---|
| [`shortq`](investigations/shortq.md) | linking の組を `ctx_idx / 3` で作っていた。NCBI は query ごとに 6 context を必ず作る（長さ 0 の context を含む）が、LOSAT の context は詰めて番号を振るので、3 nt や 4 nt の query の後の query の frame の組がずれた | 組を `2 × query + (frame < 0)` に（`06953c68d`）。169 の入力で NCBI と一致 |
| [`linktie`](investigations/linktie.md) | subject の塊ごとの init hit list の `score_compare_match`（絶対の query の位置で同点を分ける）が無く、frame の違う同点の HSP が走査の順で linking に届いた | `sort_init_hsps_by_score_ncbi`（`06953c68d`） |
| [`nisland`](investigations/nisland.md) | 最初の `BLAST_LinkHsps` の後の得点の並べ替えが無く、再評価で同点になった HSP（subject の N の島の乱数の塩基、二重の hit の端の切り詰め）が鎖の順のままで、2 回目の linking で鎖を入れ替えた | `2c9f71aa7`。乱数の入力 17000 で後退なし |
| [`segmask`](investigations/segmask.md) | NCBI の `s_SegSeq` は左の再帰の区間の一覧を `leftsegs->next = *segs` でつなぐので、先頭の 1 つだけが残る（決まった結果）。LOSAT は全部を残し、一部の低複雑度の区間を多く mask した | `aee365eb5`。SEG を直接呼ぶ約 12 万の入力で NCBI と一致。BLASTP・TBLASTN・TBLASTX の差の集まりは 0。BLASTX は SX まで前のまま |

`BLAST_LargeGapSumE` の評価の順（`blast_stat.c:4560-4561`、1 ULP の差で HSP の順が変わった）も棚卸しの結果の行（E-X2）から直した。

## 決めたこと（推奨の案。AUTHORITY.md §F）

1. TBLASTN も E2a §G.3 の説明の一覧の規則と NCBI の核酸の title（`ncbi_nucleotide_title`）を使う。凍結した TBLASTN の出力は変わらず、棚卸しの範囲 C の 2 つの再現（p417、p716：最初の HSP が最良でない）が NCBI と一致するようになった。BLASTP は S08+ まで今の一覧（protein の title は棚卸ししていない）。
2. batch、hit list、ncbi2na の乱数の塩基、表示の文字の identities を S08 で移した（DW-12）。これらは一部の入力で outfmt 6 も変える（NCBI と一致するようになる）。Gate A と S02 の case の capture は変わらない（下の「ゲート」）。
3. 明示的な拒否（「… not supported by LOSAT's TBLASTX」）：`-window_size 0`（one-hit の word finder は未移植）、NCBI が違う読み方をする定義行・文字・record、outfmt 0 の subject の HTML の文字参照、整数でない `BATCH_SIZE`、0・6・7 以外と欄の指定のある `-outfmt`、`-culling_limit` 1 以上（LOSAT の culling は NCBI と違う HSP を残す）。
4. **保守者に諮る（下）：** NCBI の `x_CleanAndCompress` が文字列の終わりを越えて読む outfmt 0 の subject の title。
5. ABI v1 は凍結のまま（TD-1）：TBLASTX は outfmt 6 だけと前の文言（`tblastx_v1_outfmt`）、hit list の大きさ無し。
6. 例外の分類はデータベースのオラクルで行う（手で書いた差でなく）。
7. 調査の 4 件と `BLAST_LargeGapSumE` は S08 で直した。BLASTX は共有の SEG と sum statistics の前の振る舞いを SX まで保つ（DW-10）。

## 保守者に諮ること

2 つの問いをまとめて諮る（保守者の指示）。どちらも今の実装は推奨でない方（明示的な拒否、または例外のない差）のままで、判断の後に S08+ で合わせる。

**1. 句読点だけの outfmt 0 の題（TBLASTX と TBLASTN、`PD-LOSAT-NCBI-DEFECTS` 版 1.1）。** 例外 2 は BLASTN だけを対象にしている。NCBI tblastx と tblastn も、`x_CleanAndCompress`（`src/objmgr/util/create_defline.cpp:219-312`）が文字列の終わりを越えて読む題（例：`, ,`、`;~ ;`）で、blastn と同じく SIGSEGV で落ちる（[`punct_defline.py`](punct_defline.py)：tblastx と tblastn の outfmt 0 で NCBI が落ち、outfmt 6 と 7 は NCBI と LOSAT が一致）。今は TBLASTX と TBLASTN が outfmt 0 でその subject を明示的に拒否する（`report/defline.rs` の `ncbi_nucleotide_title_reads_past_end`）。BLASTN が例外の前にしていたのと同じである。TBLASTN は S08 の前は、その題を整えずに出していた。

- 推奨：例外 2 を TBLASTX と TBLASTN に広げる（BLASTN と同じ代わりの題の置き換え。妥当な結果を近い入力の NCBI の出力と一致で示せる。BLASTN の `title_sweep.py` の方式で確かめる）。実装は S08+ の最初の作業にする。
- 代わり：明示的な拒否のまま。

**2. 起動の時に閉じた標準出力（`>&-`、全 program、`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`）。** NCBI は最初の書き込みで失敗する（outfmt 0 は「BLAST failed to write output」と終了コード 6、outfmt 6/7 は abort で 134）。Rust の runtime は main の前に閉じた標準の記述子へ `/dev/null` を開くので、LOSAT はそれを、呼び出し側が開いた `/dev/null` と区別できない（main の前に確かめるには `extern "C"` の関数が要り、pure-Rust の境界の規則が許さない）。LOSAT は `/dev/null` に書いて終了コード 0 になる。報告はどちらでも捨てられ、違うのは終了コードだけ。

- 推奨：承認済みの例外にする（`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` に「起動の時に閉じた標準出力では、LOSAT は報告を捨てて成功する」を足す）。
- 代わり：Linux で `/proc/self/fdinfo/1` を見る近似（読み書きの `/dev/null` を閉じたものとみなす）。呼び出し側の `subprocess.DEVNULL` などを誤って失敗させる（S08 の監査 (c) 第 2 回の N2）。

## ゲート

最後のコミット（`1117e8c17`）の全ゲートは S08b で行う。script は `~/.cache/losat-web-gui-target/s08/s08_gates.sh`（写しは [`gates/`](gates/)。`s08_gate_a.sh`、`s08_perf.sh`、`verify_added.py` も）。

- [`run-20261002T162101Z/`](run-20261002T162101Z/)（`1117e8c17`）：lint と試験の段は済んだ（`cargo fmt --check`、clippy `-D warnings` の 4 構成と adapter の 3 構成、`cargo test --all-features` 875 件通過・失敗 0、adapter と wasm32 の web API の試験、pure-Rust の境界の検査、`ci_fast_regressions.py` の単体試験、このセッションで足した行の NCBI の参照の誤り 0）。「build wasi」の段で WSL の `/mnt/c` の I/O の誤りで止まった（`NOTE.txt`）。
- 途中で止めたゲート（記録には入れない。`~/.cache/losat-web-gui-target/s08/` の `gate-968fa98c8-partial/`、`gate-7fbbfad96-partial/`）。どちらも後のコミットで置き換わったので止めた。済んだ検査はすべて通った：
  - `968fa98c8`：lint と試験（872 件）、outfmt 0 の fixture 74 件（1・2・4 スレッド）、`run_oracle.py` と `precheck_hits.py`（差は承認済みの例外の `tblastx.code4.*` 2 件だけ）、TBLASTX の fixture 53 件、BLASTN の fixture 107 件、`CTOOLKIT_COMPATIBLE` の比較 216 実行、BLASTN の得点の sweep（outfmt 0/6/7 各 300 一致・580 同じ誤り）、題の sweep（1023 の定義行、一致 957・例外 2 が 66）、速い検査の全件 236 件（失敗 0、許可した既知の不一致 1）、capture 236 件が S02 の基準とも変更前とも差 0、v1 の WASI の行列（433 の記録、形式の失敗 0）、V-ABI quick 52 実行。
  - `7fbbfad96`：lint と試験（873 件）、上の fixture と sweep のすべて（TBLASTX の fixture 69 件）、HTML の題の sweep 3354 実行、句読点の題、BLASTN の `check_inputs.py` 300 件（予期しない 0）、閉じたパイプと閉じた標準出力。

## V-PERF

S08b で行う（`~/.cache/losat-web-gui-target/s08/s08_perf.sh`。変更前はセッションの開始の実行ファイル、case は `perf_cases.py` の TBLASTX 3 つと TBLASTN・BLASTP・BLASTN）。

## 独立監査（指示書の 7.）

`ncbi_parity_auditor` の役割を、観点ごとに sonnet の agent 4 つが読み取り専用で並行して行った（共通の指示 `~/.cache/losat-web-gui-target/s08-audit/COMMON.md`、観点ごとの指示 `ANGLE_{A,B,C,D}.md`、第 2 回 `ROUND2.md`。写しは [`audit/`](audit/)）。比べた全件の入力と出力は作業ディレクトリ `~/.cache/losat-web-gui-target/s08-audit/{a,b,c,d,r2a,r2b,r2c,r2d}/` にある（大きいので入れない）。基準は `INVENTORY.tsv`。

### 第 1 回（`968fa98c8`）

| 観点 | 結論 | 内容 |
|---|---|---|
| (a) 経路の網羅 | unsupported | NCBI を callgrind の下で 29 の場面で実行し（届いた記号は約 6500、関係するライブラリで 1511）、`getenv` を gdb で追い（140 の名前）、棚卸しと突き合わせた（`COVERAGE.tsv` 515 行）。差の比較は約 985 件。**指摘：** F1 `-evalue 0` を受け付ける（中）、F2 `-seg` を任意の空白で分ける（低）、F3 0 以下の `-seg` の値を NCBI は既定のままにし LOSAT は切り詰める（中、黙った誤り）、F4 `CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`・`BL2SEQ_LEGACY` を TBLASTX が黙って無視する（低）、F5 引数の解析の後の検査の文言と終了コード（低、例外 1 の外）、F6 help の文と `-num_threads` の上限（低）。棚卸しに行の無い領域：option の検証の層、`BlastSetUp_Filter` の本体、検索の核の 114 関数（差は無し）、配列の源。検索の核は差分の実行で差が無かった |
| (b) 移植の忠実さ | supported（低い未解決 2 件） | 34418 件（すべて `968fa98c8`）。差のある 2020 件をすべて分類した：承認済みの例外 164、明示的な拒否（理由が成り立つものと後回し）423、後回しの option と引数の差 88、BLASTP だけの差（範囲外、S08 の前から）703、黙った差 642（すべて極端な値：`-seg` の 0 以下の cut-off、極端な `-evalue`。通常の値と S08 で移した部分には無い）。差のあった 2020 件と、一致した 10039 件の無作為の標本、新しい 501 件を `7fbbfad96` で再実行：584 件が一致に変わり、一致していたものの後退は 0。項目ごと（batch、hit list、ncbi2na の乱数、sum statistics と linking、説明の一覧と整列の表示、遺伝暗号、SEG、epilog と outfmt 7、書き込みの失敗、入力の読み方）に supported。**指摘：** F-1 `BLAST_Cutoffs` の下限 1 が無く、短い query に大きな `-evalue` で cutoff が負になり連結が変わる（低、S08 の前から、`7fbbfad96` でも未解決）、F-2 `-seg` の 0 以下の値（中、`7fbbfad96` で解決）、F-3 `-evalue 0`（解決）、F-4 環境変数（解決。F-4b `PRE_FETCH_SEQS_LIMIT` は未解決）、F-5 `-seg` の分け方（解決）、F-6 help（解決）、F-7 拒否の前の余計な題の警告（情報） |
| (c) 拒否の理由 | unsupported（僅か） | 977 件。拒否の理由はすべて NCBI のソースとオラクルで成り立ち、LOSAT の拒否はすべて明示的（stdout は空、`BATCH_SIZE=0` は NCBI と同じ）。句読点の題の判定は TBLASTX 15624、TBLASTN 3124 の定義行で NCBI の SIGSEGV と完全に一致。`-db_gencode` の例外は 26 のコードでデータベースのオラクルと一致。**指摘：** F1 出力に出ない subject の題でも拒否する（中）、F2 HTML の判定が NCBI の上位集合（低）、F3 DEL を制御文字として拒否（低）、F4・F5 NCBI の出力が再現できる入力の拒否（低、入力の読み方は S08+）、F6 誤りの順（低）、F7 閉じた標準出力で終了コード 0（中）、F8 `-evalue 0`（中）、F9 `-seg` の文字列（低）、F10 `-num_threads` 65535 以上（低）、F11 標準入力とパイプ（低、既知）、F12 help（低）、F13 `DIAG_POST_LEVEL`（低、既知）、F14 TBLASTN の subject の定義行を黙って読み違える（中、S08 の範囲外で S08 の前から） |
| (d) ABI v2 と構造化結果 | supported（製品の振る舞い） | 133 の検索（serial n1、threaded n1・2・4、一部 n3・8、600 subject、250 の整列の切り詰め、1 MiB を越える outfmt 0、遺伝暗号 6 種）で、1 回の実行の 0・6・7 の stream が CLI のその形式の実行とバイト一致、stream 3 も一致。HSP の記録を outfmt 0 の本文と全項目で照合し、観測者の範囲も正確（8 つの陰性対照をすべて検出）。ABI v1 は凍結のとおり。**指摘：** F1 `check_wasm_threading.py` が TBLASTX の outfmt 0/7 の失敗を期待したまま（中、毎晩の WASI の検査が赤になる）、F2・F3・F6・F8 `abi_v2.md` の記述の不足（低）、F4 outfmt 0 を必ず書く ABI の実行は、題の拒否で outfmt 6/7 だけなら通る subject でも失敗する（低、正当な拒否と記述の不足）、F10 `-evalue 0`（= (c) F8） |

### 第 1 回の指摘への対応（`7fbbfad96`）

| 指摘 | 対応 |
|---|---|
| (c) F1・F2、(d) F4 | NCBI は最終の hit list の subject の題だけを作る（`blast_format.cpp:1540`）。TBLASTX は検索の後、outfmt 0 の prolog の後に（NCBI も prolog の最初の行を flush してから落ちる）、TBLASTN は報告の前に、ヒットのある subject だけを確かめる（`check_shown_subject_titles`）。HTML の判定は NCBI の `NStr::HtmlDecode` の走査と実体の表（ncbistr.cpp:4223-4590）を、`GenerateDefline` が decode する題（最後の `.,;~ ` を落とした後）に当てる（`ncbi_nucleotide_title_is_decoded`、BLASTN も共有）。[`html_titles.py`](html_titles.py)：1118 の定義行 × 3 program、3354 実行で NCBI とバイト一致 2499、NCBI が実際に decode して LOSAT が拒否 855、不要な拒否 0。fixture `title.kept_fmt0`（`R&D; x`、`a&foo;b`、`a&amp;`、NCBI が decode しない `&xi;`・`&X41;`）、`title.hitless_fmt0` |
| (c) F3 | DEL は NCBI と同じに保つ（題は最初の 0x20 未満の byte で終わる、fasta_reader_utils.cpp:215-225）。fixture `input.del_deflines_fmt{0,7}` |
| (c) F6 | `BATCH_SIZE` を query の前、後回しにした subject の検査の前に読む（tblastx_app.cpp:136-137）。`bio` が読めない subject（最初の空行。NCBI は黙って読む）は検索の始まりで拒否する。fixture `input.subject_blank_first_empty_query`、試験 `a_batch_size_that_is_not_an_integer_is_rejected`。subject の残基の拒否（o15・o16）は、NCBI が読むときに警告を出すので、読むときのままにした（NCBI の結果も stderr が違う） |
| (c) F7 | Rust の runtime は閉じた標準の記述子に `/dev/null` を読み書きで開く。`7fbbfad96` は Linux でそれを `/proc/self/fdinfo/1` で見分けて（shell の `> /dev/null` は書き込みだけ）書き込みを失敗させたが、第 2 回の監査 (c) の N2 のとおり、呼び出し側が読み書きで開いた `/dev/null`（Python の `subprocess.DEVNULL`、Node の `stdio: 'ignore'`）も失敗させたので外した（`1117e8c17`、試験 `a_standard_output_on_dev_null_is_written`）。main の前に確かめるには `extern "C"` の関数が要り、pure-Rust の境界の検査が許さない。閉じた標準出力では LOSAT は `/dev/null` に書いて成功する。保守者に諮る（下の「保守者に諮ること」の 2） |
| (c) F8、(a) F1、(d) F10 | NCBI の hit saving の検査（blast_options.c:1518-1523）を移植した。`Query is Empty!` の前、`-max_target_seqs` の警告の後。fixture `options.evalue0_*`、adapter の `validate` |
| (a) F2・F3、(c) F9 の一部 | `-seg` を 1 つの空白で 3 つに分け、窓は `int`、0 以下の窓・locut・hicut は NCBI の既定のまま（blast_filter.c:1147-1154）。BLASTP・TBLASTN・TBLASTX で共有。outfmt 6 は 9 の値 × 3 program で NCBI と一致（BLASTP の outfmt 0 は S08 の前からの `Method:` の文言の差だけ）。fixture `options.seg_*`。誤りの文言と終了コード 1 は S08+ |
| (a) F4 | `CHUNK_SIZE` が 3 で割れない（`size_t` として）と NCBI と同じ「BLAST engine error: Split query chunk size must be divisible by 3」（終了コード 3、prolog と最初の batch の題の警告の後）。NCBI が変換できない `CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE` と `BL2SEQ_LEGACY` は明示的な拒否。fixture `env.chunk*`、`env.overlap_negative_fmt6` |
| (a) F6、(c) F12 | help の文を直した |
| (d) F1 | `check_wasm_threading.py` は TBLASTX の outfmt 0/7 を NCBI のオラクルと比べる（S07 の BLASTN と同じ） |
| (d) F2・F3・F6・F8 | `abi_v2.md` を直した。`validate` は TBLASTX の option だけの検査（`check_options`）も行う |
| (a) F5、(c) F4・F5・F9 の残り・F10・F11・F13・F14 | S08+ に記録した（[S08+ の指示書](../../losat_web_gui_sessions/session_s08p_e2e_protein_options.md)の「S08 からの引き継ぎ」、`INVENTORY.tsv` の範囲 H と A24・E50）。(c) F14 は S08+ の最初の作業 |
| (a) 棚卸しの行 | 範囲 H の 14 行を足し、A24・E50・A69・D51 の最後の状態を直した（`build_inventory.py` の `AUDIT_FINAL`・`AUDIT_ROWS`） |
| (b) F-1・F-4b（`aa4b6f8c3`） | `BLAST_Cutoffs` は呼び出し側の 1（blast_parameters.c:321、916）と E からの得点の大きい方を返す（blast_stat.c:4097、4126-4129）。TBLASTX の 2 つの呼び出し（`ncbi_cutoffs.rs` の `blast_cutoffs_from_one`）に移した（BLASTX は SX まで自分の呼び出しのまま）。監査の再現 12 件が NCBI と一致（ゲートの実行ファイルは 6 件で違う）、fixture `cutoff.floor_evalue_1e10`。NCBI が変換できない `PRE_FETCH_SEQS_LIMIT` は BLASTN と同じに拒否する。試験 `environment_that_ncbi_cannot_convert_or_reports_otherwise_is_rejected` |
| (b) F-7 | 情報（拒否の前の題の警告。実行は失敗する）。S08+ の入力の読み方で扱う |

### 第 2 回（`7fbbfad96`）

| 観点 | 結論 | 内容 |
|---|---|---|
| (a) 経路の網羅 | unsupported（`7fbbfad96` について。N2・N3 は `aa4b6f8c3` で直した） | 第 1 回の 1014 行のうち 918 行を同じ引数で再実行（他は新しい行列 221 件で再実行）：第 1 回に一致した 789 行は今も一致（1 行は再実行の引用符の誤り）、26 行が一致に変わり、90 行は記録した種類の差。後退 0。第 1 回の F1〜F7 はすべて直したか記録した後回し。新しい差分の実行は 1600 の無作為の組（0 差）、SEG の 288 件、`CHUNK_SIZE` の 341 件、HTML の 324 件（不要な拒否も見逃しも 0）、`x_CleanAndCompress` の 400 件、接頭辞 136 件、DEL 18 件。`getenv` は LD_PRELOAD の shim で 148 の名前。**指摘：** N3 `BLAST_Cutoffs` の下限（中、= (b) F-1）、N2 `PRE_FETCH_SEQS_LIMIT`（低、= (b) F-4b）、N1 空の query と NCBI 自身が拒否する subject（BOM など）で LOSAT が終了コード 0（低）、N4 `>` だけの subject（低、後回し）。H13 の一括の「faithful」の行は安全でないので、統計の関数を関数ごとの行に分けるよう勧めた |
| (b) 移植の忠実さ | supported | 第 1 回の再実行（上の表）を `7fbbfad96` で行った。F-1・F-4b は `aa4b6f8c3` で直した |
| (c) 拒否の理由 | unsupported（僅か。N1・N2 は第 1 回の修正が入れた差で、`5fa53b3f0`・`1117e8c17` で直した） | 第 1 回の 959 件をすべて再実行：一致から差への後退 0、差から一致 42。第 1 回の 14 の指摘はすべて直したか記録した後回し。HTML の判定は NCBI の `HtmlDecode` と行ごとに同じ、実体の表 280 も同じで、14082 の題（TBLASTX と TBLASTN）で不要な拒否も見逃しも 0。句読点の題の fuzz は第 1 回と同じ（TBLASTX 272 = 272、TBLASTN 66 = 66）。`CHUNK_SIZE` の 40 件、変換できない値の 34 × 2、`-seg` の 135 件、`-evalue` の順の 21 の組、DEL と制御文字の 312 件。**指摘：** N2 呼び出し側が読み書きで開いた `/dev/null` を閉じた標準出力と取り違える（中）、N1 NCBI が読むときに拒否する subject（BOM など）と空の query で終了コード 0（低）、N3 `PRE_FETCH_SEQS_LIMIT`（低、`aa4b6f8c3` で直した）、N4 TBLASTN の `-max_target_seqs` 5 未満の警告が無い（低、範囲外、S08 の前から） |
| (d) ABI v2 と構造化結果 | supported | 第 1 回の組を新しい reactor で再実行：core 33、edge 74、big 25、LC738874 × LC738875、スレッド 3・8、`v_abi.js` の TBLASTX（凍結の NCBI のハッシュ 24/24・72/72）、quick、8 つの陰性対照、40 回の繰り返し。後退 0。第 1 回の F1・F2・F6・F10 は直り、F3・F4・F8 は記述どおり、F5 は正当な拒否、F7 は試験の不足（欠陥でない）。ABI を通した題の 84 件と 600 subject の順位の試験で、拒否は最終の hit list の subject だけ。**指摘：** N1 空の query の run より `validate` が厳しい（低）、N2 2 つの問題のある subject で `register` と CLI の文言の順が違う（低、S08 の前から）、N3 `abi_v2.md` の文言（低）、N4 `aa4b6f8c3` の後の成果物の作り直しと V-ABI（手順） |

### 第 2 回の指摘への対応

| 指摘 | 対応 |
|---|---|
| (a) N3・(b) F-1、(a) N2・(b) F-4b | `aa4b6f8c3`（上の「第 1 回の指摘への対応」の表の (b) の行）。監査 (a) 第 2 回の下限の 56 件と環境の 81 件を直した実行ファイルで再実行：56 件が一致、環境は一致 79・明示的な拒否 2 |
| (a) N1 | `5fa53b3f0`：`bio` が読めない subject を検索の始まりまで後回しにするのは、最初の定義行の前の行がすべて NCBI の飛ばす行（空白だけか、`!`・`#`・`;` で始まる注釈、fasta.cpp:375-384）か、定義行が UTF-8 でないときだけ。BOM、NCBI が止まる行、定義行の無いレコードは読むときに拒否する。空の query との組み合わせ 8 種で、NCBI が終了コード 0 の 5 種は一致、NCBI が終了コード 1 の 2 種は LOSAT も拒否（終了コード 1）、NCBI が定義行の無いレコードとして読む 1 種は LOSAT が拒否（明示的） |
| (a) N4 | S08+ に記録した（中身の無いレコードの subject 側） |
| (a) 棚卸しの行 | H13 を「`BLAST_Cutoffs` を除く」とし、`BLAST_Cutoffs`・`BlastKarlinEtoS_simple`・`BLAST_KarlinStoE_simple`・`BLAST_SmallGapSumE`/`BLAST_UnevenGapSumE`・`BLAST_ComputeLengthAdjustment`・`s_GetNextSubjectChunk`（5 Mnt を越える subject）・lookup の bone の種類・`s_PreFetchSeqs`・最初の定義行の前の行の 10 行を足した（範囲 H は 24 行） |
| (d) N1・N3 | `6017f0ea2`：`abi_v2.md` に、空の query の run より `validate` が厳しいこと、各入力を別々に確かめること、0x20 未満の byte、最終の hit list の subject、TBLASTN の文言を書いた |
| (d) N2 | `6017f0ea2`：`register` は CLI と同じ順（query は定義行が先、subject は残基の後）で確かめる |
| (d) N4 | 最後のゲート（下の「ゲート」）で、最後のコミットから成果物を作り直し、V-ABI を全升目で再実行した |
| (c) N1 | (a) N1 と同じ（`5fa53b3f0`） |
| (c) N2 | `1117e8c17`（上の (c) F7 の行） |
| (c) N4 | S08+ に記録した |

### 第 3 回

S08b で行う（指示 [`audit/ROUND3.md`](audit/ROUND3.md)。最後の実行ファイルで、第 2 回の未確認の指摘と後退を確かめる）。

## アプリ側（S09）への注意

- TBLASTX は ABI v2 で形式 0・6・7 を出し、stream 1 に HSP の記録を出す（翻訳した行、SEG で mask した query の文字は小文字、両側の frame。`docs/web/abi_v2.md` §8）。
- `describe` の TBLASTX の `-max_target_seqs` に `default` が無くなった（省略すると hit list は 500 で、outfmt 0 は説明 500・整列 250。値を与えると両方がその値。clap の既定値では表せないので、help の文に「(default: 500)」と書いた）。アプリの検索画面の既定値の表示はこの help か、未指定として扱う。
- `parse` は `-outfmt` を挿入しなくなった。argv の誤りは、未知の program を含めて CLI と同じ文言になる。
- TBLASTX の register は BLASTN と同じに入力を確かめる（NCBI が違う読み方をする record は `not supported by LOSAT's TBLASTX` の誤り、`U` は `T`）。

## 残件と引き継ぎ

- **S08b（次のセッション、[指示書](../../losat_web_gui_sessions/session_s08b_e2b_final_gates.md)）**：最後のコミットでの全ゲート、独立監査の第 3 回、Gate A、V-PERF、`docs/web/verification_cells.tsv` の TBLASTX の outfmt 0/7 の行（26・28）、この記録と `evidence.sha256` の仕上げ、`main` への PR（`4fc67f9ab` と `2bcb86b1f` を含む。merge しない）。S08 の完了条件はここで満たす。
- **保守者の判断（上の「保守者に諮ること」の 1 と 2）**：句読点だけの outfmt 0 の題（TBLASTX と TBLASTN に例外 2 を広げるか）、起動の時に閉じた標準出力（承認済みの例外にするか）。判断の後に S08+ で合わせる。
- **S08+**（[指示書](../../losat_web_gui_sessions/session_s08p_e2e_protein_options.md)の「S08 からの引き継ぎ」）：最初の作業は TBLASTN の subject の定義行（監査 (c) の F14、黙った読み違え）。ほかに、引数の解析の後の NCBI の検査の文言と終了コード 1（`-threshold 0`、壊れた `-seg`、負の `-evalue`）、TBLASTX の入力の読み方の安く移せる拒否、標準入力とパイプ、`-num_threads` 65535 以上、TBLASTX の `check_ncbi_application_settings`、BLASTN の HTML の題の検査を表示される subject に狭めること、BLASTP の差（説明の一覧、`Method:` の文言、SEG の小文字、全体が mask された query の Karlin の警告、`-evalue 1000` の 1 残基の HSP）、TBLASTN の `-max_target_seqs` 5 未満の警告、中身の無い subject のレコード。
- **SX**（[指示書](../../losat_web_gui_sessions/session_sx_blastx_integration.md)の「S08 からの引き継ぎ」）：BLASTX は共有の SEG・`BLAST_LargeGapSumE`・説明の一覧・`-seg` の値・`BLAST_Cutoffs` の呼び出しを、取り込むまで前のまま（DW-10）。
- **アプリ側（S09）**：上の「アプリ側（S09）への注意」。S09 の入口の条件（S08 の完了）は S08b で満たす。
