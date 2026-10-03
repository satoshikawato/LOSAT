# LOSAT Web E2b（Session S08）ゲート記録

- 段階：E2b TBLASTX の outfmt 0 と 7（[総合計画書](../../losat_web_gui_plan.md) §7 の S08、指示書 [S08](../../losat_web_gui_sessions/session_s08_e2b_tblastx_outfmt0_7.md)、計画 DW-6・DW-10・DW-12・DW-15・TD-1）
- ブランチ：`feature/losat-web-gui`。変更前はセッションの開始の `2bcb86b1f`（エンジンは E2g の最後の `f4057718a` と同じ。native の SHA-256 `331fcba36447…44d8`）。変更後は `24fcfe41b`（最後のエンジンのコミット。S08 の移植は `1117e8c17` まで、S08b の `bc521f450` は独立監査の第 3 回の指摘の修正、`24fcfe41b` は注釈だけ。native `2a46c2e4…`）。セッションは S08（2026-10-02〜03）と S08b（2026-10-03、指示書 [S08b](../../losat_web_gui_sessions/session_s08b_e2b_final_gates.md)）
- 状態：**完了（2026-10-03、S08b）**。最後のエンジンのコミットの全ゲート、Gate A、V-ABI full、独立監査の第 3 回（4 観点とも supported、その指摘の修正の確かめも supported）を通した。V-PERF は 1 つの case（`tblastx-multi`、NCBI の batch の費用）の判断を保守者に確かめる（下の「V-PERF」）。保守者の 2 つの判断（DW-17）を記録した

## 完了条件（計画 §7 の S08 の行）

| 条件 | 状態 | 根拠 |
|---|---|---|
| 固定した fixture で NCBI とバイト一致（承認済みの遺伝暗号の例外を除く） | 満たした | 最後のゲート（下の「ゲート」、`run-20261003T021918Z/`）：outfmt 0 の fixture 74 件（1・2・4 スレッド）差 0、NCBI の凍結の確かめ差 0、`precheck` の差は承認済みの `tblastx.code4.*` だけ、TBLASTX の回帰 fixture 73 件差 0、BLASTN の fixture 107 件差 0 |
| 既存の 6 に退行なし | 満たした | Gate A（`audit_tblastx_v010.py`）20 case が v0.1.0 と同じ分類（parity 14 が `EXACT_TEXT`、承認済みの 6 が `HSP_SET_DIFF`）、capture 236 件が S02 の基準とも変更前とも差 0、CI の速い検査の全件 236 件で失敗 0、v1 の WASI の行列で形式の失敗 0 |
| TBLASTX の全升目の V-ABI | 満たした | V-ABI full 524 実行（TBLASTX の 36 の検索を含む 131 の検索 × serial n1・threads n1・n2・n4、形式 0・6・7）で失敗した部分 0、quick 52 実行。`docs/web/verification_cells.tsv` の 26・28 行を `checked`（S08） |
| 独立監査 | 満たした | 第 3 回（`1117e8c17`）で 4 観点とも supported。その低い指摘（(a) N1-ii・N1-iii、(b) N-b1、(d) L1）は `bc521f450` で直し、4 観点の確かめ（`2a46c2e4…`）も supported。残りは明示的な拒否か記録した後回し（下の「独立監査」） |
| V-PERF（指示書 S08b の 3.） | 1 case を保守者に確かめる | 27 組のうち `tblastx-multi` の native ×1.204・serial-WASI ×1.106（`--repeat 5`）。NCBI の batch の構造の費用（下の「V-PERF」）。ほかは閾値以内 |

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
| `1117e8c17` | 閉じた標準出力の見分け方を外す（監査 (c) 第 2 回の N2。S08b で承認済みの例外 6）。`-culling_limit` の help |
| `bc521f450`（S08b） | 最初の定義行の前に NCBI が飛ばす行（空白、`!`・`#`・`;` の注釈）だけがある subject：そのような行だけのファイルは subject の無いファイル（`Empty CBlastQueryVector`、終了コード 3）、その後のレコードは読むときに残基を確かめ title の警告を書く（監査 (a)・(b)・(d) 第 3 回）。fixture 3 件と CLI の試験 |
| `24fcfe41b`（S08b） | 注釈だけ（`cli.rs` の閉じた標準出力は承認済みの例外 6） |

## 棚卸し（DW-12）

[`INVENTORY.tsv`](INVENTORY.tsv)（447 行、[`build_inventory.py`](build_inventory.py) が作る）。NCBI の tblastx の outfmt 0/6/7 の経路を 7 つの範囲（A アプリと orchestration、B 整列、C 説明の一覧、D outfmt 7/6、E 報告が検索から受け取るもの、F 文字と遺伝暗号、G hit list）に分け、読み取り専用の agent（sonnet）が移植の前に行を作り（`inventory/{A..G}.tsv`、規則は `inventory/COMMON.md`。列 `status` は S08 の前の LOSAT）、別の agent が移植の後の `0d533ba76` で行ごとの結果を確かめた（`inventory/result_{A..G}.tsv`、規則は `inventory/RESULT_COMMON.md`。NCBI との比較の実行は合わせて 3000 を超える）。GAP とされた 15 行は直すか拒否した（`build_inventory.py` の `RESOLUTION`、列 `s08_final`）。

| | faithful | reused | ported | rejected | exception | deferred | n/a |
|---|---|---|---|---|---|---|---|
| 移植の前（`status_before`） | 64 | 92（reusable） | — | 8 | — | — | 88 |
| 最後（`s08_final`） | 63 | 64 | 195 | 23 | 7 | 7 | 88 |

移植の前の残り：divergent 29、missing 61、needs-param 80（行の数は 447。result の agent が足した 2 行と、独立監査が足した範囲 H の 23 行は移植の前の状態が無い）。`deferred` の 7 行は S08+ の範囲（NCBI の実数・整数の引数と `-outfmt` の文字列の読み方、`-subject` の欠落、NCBI のファイルの誤りの文言、`-query -`、引数の解析の後の検査の文言など）。閉じた標準出力（H11）は S08b で承認済みの例外 6 になった（`exception` の 7 行目）。

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
4. **保守者の判断（下、S08b）：** NCBI の `x_CleanAndCompress` が文字列の終わりを越えて読む outfmt 0 の subject の title は、例外 2 を TBLASTX と TBLASTN に広げた（実装は S08+）。
5. ABI v1 は凍結のまま（TD-1）：TBLASTX は outfmt 6 だけと前の文言（`tblastx_v1_outfmt`）、hit list の大きさ無し。
6. 例外の分類はデータベースのオラクルで行う（手で書いた差でなく）。
7. 調査の 4 件と `BLAST_LargeGapSumE` は S08 で直した。BLASTX は共有の SEG と sum statistics の前の振る舞いを SX まで保つ（DW-10）。

## 保守者の判断（S08b、2026-10-03、計画 DW-17）

S08 が記録した 2 つの問いを S08b でまとめて諮り、保守者はどちらも推奨の案を選んだ。

**1. 句読点だけの outfmt 0 の題（TBLASTX と TBLASTN、`PD-LOSAT-NCBI-DEFECTS` 版 1.1）。** 例外 2 は BLASTN だけを対象にしている。NCBI tblastx と tblastn も、`x_CleanAndCompress`（`src/objmgr/util/create_defline.cpp:219-312`）が文字列の終わりを越えて読む題（例：`, ,`、`;~ ;`）で、blastn と同じく SIGSEGV で落ちる（[`punct_defline.py`](punct_defline.py)：tblastx と tblastn の outfmt 0 で NCBI が落ち、outfmt 6 と 7 は NCBI と LOSAT が一致）。今は TBLASTX と TBLASTN が outfmt 0 でその subject を明示的に拒否する（`report/defline.rs` の `ncbi_nucleotide_title_reads_past_end`）。BLASTN が例外の前にしていたのと同じである。TBLASTN は S08 の前は、その題を整えずに出していた。

- **判断：例外 2 を TBLASTX と TBLASTN に広げる**（`PD-LOSAT-NCBI-DEFECTS` 版 1.2、AGENTS.md）。BLASTN と同じ代わりの題の置き換えの移植と、BLASTN の `title_sweep.py` の方式での確かめは S08+ の最初の作業。それまでは明示的な拒否のまま（承認済みの例外の範囲の、より厳しい振る舞い）。

**2. 起動の時に閉じた標準出力（`>&-`、全 program、`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`）。** NCBI は最初の書き込みで失敗する（outfmt 0 は「BLAST failed to write output」と終了コード 6、outfmt 6/7 は abort で 134）。Rust の runtime は main の前に閉じた標準の記述子へ `/dev/null` を開くので、LOSAT はそれを、呼び出し側が開いた `/dev/null` と区別できない（main の前に確かめるには `extern "C"` の関数が要り、pure-Rust の境界の規則が許さない）。LOSAT は `/dev/null` に書いて終了コード 0 になる。報告はどちらでも捨てられ、違うのは終了コードだけ。

- **判断：承認済みの例外 6 にする**（`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 版 1.3、AGENTS.md）。実装の変更は無い（注釈だけ `24fcfe41b`）。証拠 [`closed_stdout/`](closed_stdout/)：BLASTN・TBLASTX・BLASTP・TBLASTN・BLASTX の outfmt 0/6/7 の 15 の組で、LOSAT は終了コード 0 と空の stderr、NCBI は outfmt 0 で 6、6/7 で 134。

## ゲート

最後のエンジンのコミット `24fcfe41b`（`bc521f450` の後の注釈だけのコミット）の全ゲート：run [`run-20261003T021918Z/`](run-20261003T021918Z/)（`head.txt`）。script は `~/.cache/losat-web-gui-target/s08/s08_gates.sh` と `s08_gate_a.sh`（写しは [`gates/`](gates/)）。S08b で、Gate A の字句のパス（`/tmp/losat-pr5-runtime-cert-5845d22/LOSAT/tests/fasta`、WSL の再起動で消える）を CI と同じ `ci_fast_regressions.py` の `stage_lexical_fixtures` で用意する段を script の始めに足した。成果物のハッシュは `artifacts.sha256`（native `2a46c2e4…`。独立監査の第 3 回の確かめの実行ファイルと同じ）。変更前はセッションの開始の実行ファイル（`~/.cache/losat-web-gui-target/s08/bin/LOSAT-base`、`331fcba36447…44d8`）。

| 検査 | 結果 |
|---|---|
| `cargo fmt --check`（LOSAT、adapter）、clippy `-D warnings` 4 構成と adapter の 3 構成 | すべて終了コード 0 |
| `cargo test --all-features`、adapter、wasm32 の web API の試験 | 876 件通過・失敗 0、adapter 7 件、web API 5 件 |
| pure-Rust の境界、`ci_fast_regressions.py` の単体試験、このセッションで足した行の NCBI の参照 | 通過、誤り 0（`verify-refs-session-added.txt`） |
| outfmt 0 の fixture（`check_losat.py`、1・2・4 スレッド） | 74 件、差 0（3 つのスレッド数とも） |
| `run_oracle.py`（NCBI の凍結の確かめ）、`precheck_hits.py` | 差 0。`precheck` の差は承認済みの `-db_gencode` の例外の `tblastx.code4.0`・`tblastx.code4.7` だけ |
| TBLASTX の回帰 fixture（`tblastx_regression_fixtures.py check`） | 73 件、差 0。区別の確認：変更前の実行ファイルは 68 件、`0d533ba76` は 30 件で違う。`1117e8c17` は S08b で足した 3 件で違う（上の「独立監査」） |
| 環境変数の区別（`env_discrimination.py`） | 21 件、予期しない 0 |
| BLASTN の fixture（`blastn_regression_fixtures.py check`） | 107 件、差 0 |
| `CTOOLKIT_COMPATIBLE`（`ctoolkit_compare.py`） | 216 実行、差 0 |
| 句読点だけの題（`punct_defline.py`） | 予期しない 0（tblastx・tblastn の outfmt 0 で NCBI が落ち、LOSAT は明示的に拒否。outfmt 6/7 は一致） |
| HTML の題（`html_titles.py`） | 1118 の定義行 × 3 = 3354 実行：一致 2499、NCBI が decode して LOSAT が拒否 855、不要な拒否と差 0 |
| BLASTN の入力（E2g の `check_inputs.py`） | 300 件、予期しない 0 |
| BLASTN の得点の sweep（E2c、outfmt 0/6/7） | 各 300 一致・580 同じ誤り |
| BLASTN の題の sweep（E2g の `title_sweep.py`） | 1023 の定義行：一致 957・例外 2 が 66、予期しない 0 |
| 閉じたパイプ（`closed-pipe.txt`） | outfmt 0 は終了コード 6「BLAST failed to write output」、6/7 は 1（承認済みの例外 3・5）。起動の時に閉じた標準出力は LOSAT が 0、NCBI が 6 と 134（承認済みの例外 6） |
| CI の速い検査の全件（`ci_fast_regressions.py --all-cases`） | 236 件、失敗 0、許可した既知の不一致 1（`Sakai.MG1655.megablast`） |
| capture（`capture_outputs.py`） | 236 件が S02 の基準とも変更前とも差 0 |
| v1 の WASI の行列（`check_wasm_threading.py`） | 433 の記録、reactor の lifecycle の gate 通過、形式の失敗 0。reactor の記録の S05 との差 21 は E2g と同じ内容（LOSAT の文言の差だけ）。`v1-requests` は E1d と一致 |
| V-ABI quick | 52 の実行が native の CLI と一致 |
| V-ABI full | 131 の検索（BLASTN 58、TBLASTX 36、TBLASTN 28、BLASTP 9）× 4（serial n1、threads n1・n2・n4）= 524 の実行、失敗した部分 0。各実行の 0/6/7 の stream が native の CLI と一致。凍結ハッシュ 776 件中 772 件一致、違う 4 件は既知の `Sakai.MG1655.megablast` outfmt 7。TBLASTX の 36 の検索は、凍結の NCBI のハッシュを outfmt 0 で 14、6 で 20、7 で 10 の検索が持つ |
| Gate A（`audit_tblastx_v010.py`、TBLASTX v0.1.0 の outfmt 6） | 20 case：parity の 14 件が `EXACT_TEXT`、承認済みの `-db_gencode` の 6 件が `HSP_SET_DIFF`（契約 PASS）。v0.1.0 の認証と同じ分類（`audit-tblastx-v010/`）。Gate A の出力のハッシュは capture（上）で S02 の基準と同じ |

reactor（`d2db178a…`・`8b179d93…`）は `bc521f450` で変わった（TBLASTX の library を含む）。`bc521f450` が変えたのは CLI の `run` だけで、ABI の確かめは上の V-ABI quick と full、独立監査 (d) の第 3 回の確かめ（`1117e8c17` の reactor と新しい native）で行った。

途中で止めた run（記録には入れない。`~/.cache/losat-web-gui-target/s08/stopped-runs/`）：`run-20261003T001302Z`（`1117e8c17`、lint と試験の後、V-ABI の case を作る段で Gate A の字句のパスが無く止まった）、`run-20261003T002556Z`（`1117e8c17`、v1 の WASI の行列までの全検査が通った後、独立監査の第 3 回の指摘を直すために V-ABI と Gate A を止めた）、S08 の `run-20261002T162101Z`（`1117e8c17`、WSL の I/O の誤り。この記録から外した）。S08 の部分の run は `~/.cache/losat-web-gui-target/s08/` の `gate-968fa98c8-partial/`・`gate-7fbbfad96-partial/`。

## V-PERF

`~/.cache/losat-web-gui-target/s08/s08_perf.sh` と `s08_perf_rerun.sh`（写しは [`gates/`](gates/)）、V-PERF の lock を取って（アプリ側の S09 は止まる）。変更前はセッションの開始の実行ファイル（native `331fcba36447…`、E2g の WASI の成果物）、変更後は上のゲートの成果物。case は [`perf_cases.py`](perf_cases.py) の 9 つ × native・serial-WASI・threaded-WASI（4 スレッド）、変更前と変更後を 1 回ずつ交互に（暖機 1 回）。先に、9 つの case の出力が変更前と変更後で同じことを native の 1・4 スレッドで確かめた（`1117e8c17` の実行ファイルで 18 組すべて同じ、`perf-precheck-1117e8c17.tsv`。S08 が出力を変えた入力は case に無い。最後の実行ファイルでは、V-PERF の計測そのものが 27 組すべてで出力の同じことを確かめた）。計算機は静かでなかった（同じ計算機でアプリ側の S09 の Playwright と別のセッションの試験。交互の計測で両方に同じだけ掛かる）。

| 計測 | 結果 |
|---|---|
| `perf-1`（`--repeat 3`、27 組） | 22 組が閾値（×1.05）以内、出力はすべて同じ。超えた 5 組：`tblastx` serial-WASI ×1.064、`tblastx-multi` native ×1.135・threaded-WASI ×1.079、`blastp-fmt0` serial-WASI ×1.180・threaded-WASI ×1.577（変更前の範囲 0.515〜1.304 秒） |
| `perf-2`（超えた 3 つの case を `--repeat 5`） | `tblastx` は 3 つとも以内（×1.009〜1.029）、`blastp-fmt0` も以内（×0.755〜1.007）。`tblastx-multi` は native ×1.204（0.100 → 0.120 秒）、serial-WASI ×1.106（0.206 → 0.228 秒）、threaded-WASI ×0.961 |
| outfmt ごとの費用（`perf-formats.txt`、変更後の native） | LC738884 × LC741431：outfmt 6 は 0.900 秒、0 は 0.944 秒、7 は 0.913 秒。LC738874 × LC738875：6 は 2.846 秒、0 は 2.925 秒、7 は 2.938 秒 |

**`tblastx-multi` の判断（推奨の案、保守者の確認待ち）：** S08 の移植した NCBI の batch のため。4 つの query（10000・300・10000・10000 nt）は NCBI の 10002 nt の batch で 2 回の検索になり（変更前は 1 回）、batch ごとに lookup を作り直し（`LOSAT_TIMING`：8 + 11 ms、変更前は 14 ms）、244 kb の subject の準備（6 frame の翻訳、ncbi2na の乱数の塩基）を繰り返す。走査と ungapped の拡張の合計は変わらない。NCBI も batch ごとに検索全体を行う（`tblastx_app.cpp` の batch の loop と `CLocalBlast`）。batch は他の入力で出力を NCBI と同じにするのに要る（linking の cutoff、検索されない batch、hit list）。1 つの batch の検索（全ゲノムの `tblastx`、260 subject の `tblastx-many`）と、TBLASTN・BLASTP・BLASTN に退行は無い。推奨：NCBI の batch の構造の費用として認め、batch の間で subject の準備を使い回す最適化（出力が同じなら許される）を後のセッションの項目にする（S08+ の指示書の「S08 からの引き継ぎ」）。

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

### 第 3 回（`1117e8c17`、S08b）

指示 [`audit/ROUND3.md`](audit/ROUND3.md)。対象は最後のゲートの成果物の写し（`~/.cache/losat-web-gui-target/s08-audit/r3-bin/`：native `55a7e574…`、reactor `d4a1bc43…`・`905cc994…`。ゲートの `artifacts.sha256` と同じ）。作業ディレクトリは `s08-audit/r3{a,b,c,d}/`。途中で Claude Code の再起動で止まり、同じ文脈で再開した。保守者の判断（上）の後は、閉じた標準出力を承認済みの例外 6、句読点だけの題を例外 2 の範囲の明示的な拒否として分類させた。

| 観点 | 結論 | 内容 |
|---|---|---|
| (a) 経路の網羅 | supported（低い残り 2 件、`bc521f450` で直した） | 第 2 回の 5039 実行のうち 4896（97%）を再実行：一致のまま 4289、記録した種類の差のまま 605、一致から差 2（定義行の無い配列だけの subject と空の query：`5fa53b3f0` が読むときに明示的に拒否するようにした、記録した後回しの範囲）。新しい実行は約 6500（`BLAST_Cutoffs` の下限の sweep 444 件はすべて一致）。第 2 回の N3 は直った、N2 は理由の成り立つ拒否、N1 は BOM・UTF-16・NUL で直った、N4 は記録した後回し、H13 は関数ごとの行に分けた（`s_BlastSumP`・`BLAST_GapDecayDivisor`・`BLAST_Powi` は H13 のまま、差分の実行で差なし）。**指摘：** N1-ii 最初の定義行の前に NCBI が飛ばす行があり、もっともらしくない配列の行がある subject と空の query（NCBI 終了コード 1、LOSAT は `Query is Empty!` で 0、低）、N1-iii 注釈の行だけの subject と空の query（NCBI は `Empty CBlastQueryVector` で 3、LOSAT は 0、低） |
| (b) 移植の忠実さ | supported（低い 1 件、`bc521f450` で直した。もう 1 件は記録した後回し） | 第 2 回の組から 9095 件（第 1 回の 34418 件の 26%）を再実行：一致から差 0、差から一致 14（F-1 の `Ev2.*.sq`）。F-1（`BLAST_Cutoffs` の下限）は直った（新しい 1182 件も一致、`7fbbfad96` の実行ファイルは全ゲノムの subject の 18 件すべてで違う）、F-4b（`PRE_FETCH_SEQS_LIMIT`）は理由の成り立つ拒否（28 の値で、NCBI が変換する値と LOSAT が受け付ける値が同じ）。`5fa53b3f0` の 494 件。**指摘：** N-b1 注釈の行だけの subject と空の query（= (a) N1-iii、低）、N-b2 最初の subject が中身の無い `>` のレコード（低、記録した後回し「中身の無いレコード」） |
| (c) 拒否の理由 | supported（低い 1 件と情報 2 件。どれも明示的か記録した後回し） | 第 2 回の 2156 行をすべて再実行：第 1 回の 959 件で一致から差 0、第 2 回から差に変わった 8 件はすべて下の L1 の種類、ほかの変化 36 件は改善か文言だけ。第 2 回の N1 は終了コードが直った（文言は後回し）、N2（呼び出し側の `/dev/null`）は直った、N3 は理由の成り立つ拒否、N4（TBLASTN の警告）は記録した後回し。閉じた標準出力は例外 6 の記述どおり（24 実行）。句読点だけの題：TBLASTX・TBLASTN とも長さ 5 までの 3124 の題で、NCBI が落ちる 66 と LOSAT の拒否が同じ、残り 3058 はバイト一致。`HtmlDecode` は 6826 の題で不要な拒否も見逃しも 0。拒否 1572 件はすべて明示的。**指摘：** L1 最初の定義行の前に NCBI が読む文字列（定義行の無い配列、` >s` など）がある subject と空の query で、NCBI は 0、LOSAT は読むときに拒否して 1（低、明示的、記録した後回し）、L2 `PRE_FETCH_SEQS_LIMIT=abc` と `CHUNK_SIZE=10` の組の誤りの順（情報）、L3 TBLASTN が `PRE_FETCH_SEQS_LIMIT=abc` を受け付ける（範囲外、情報） |
| (d) ABI v2 と構造化結果 | supported（低い 1 件、CLI だけ、`bc521f450` で直した） | 第 2 回の 630 実行と比べて 624 が同じ、6 は第 2 回の N2 の修正による文言か harness の差。`v_abi.js` の TBLASTX（serial 15、threads 45、凍結の NCBI のハッシュ 24/24・72/72）、quick 52、`BLAST_Cutoffs` の下限の新しい 32 case（CLI と NCBI 32/32、ABI と CLI 32/32・64/64）、8 つの陰性対照。第 2 回の N1・N3（`abi_v2.md`）と N2（`register` の順）は直った。**指摘：** L1 飛ばす行の後の subject に無効な残基があり query が空のとき、NCBI の残基の警告が無い（低、CLI だけ） |

### 第 3 回の指摘への対応

| 指摘 | 対応 |
|---|---|
| (a) N1-ii・N1-iii、(b) N-b1、(d) L1 | `bc521f450`：NCBI は query を見る前に subject を読み（tblastx_app.cpp:118-132）、飛ばす行（fasta.cpp:375-384）の後のレコードをほかのファイルと同じに読む。そのような行だけのファイルは subject の無いファイル（`Empty CBlastQueryVector`、終了コード 3、NCBI とバイト一致）、その後のレコードは読むときに残基を確かめ（拒否、終了コード 1）title の警告を書く。fixture `input.subject_comment_only_{empty_query,fmt0}`・`input.subject_comment_first_title_empty_query`（`1117e8c17` は 3 件とも違う）、試験 `subjects_after_skipped_lines_are_read_before_the_query_is_found_empty` |
| (b) N-b2、(c) L1・L2・L3 | S08+ に記録した（指示書の「S08 からの引き継ぎ」の入力の読み方と環境の項） |

### 第 3 回の確かめ（`bc521f450`、S08b）

`bc521f450` の実行ファイル（`s08-audit/r4-bin/LOSAT`、`2a46c2e4…`。最後のゲートの native と同じ、下の「ゲート」）で、4 つの観点が同じ文脈で確かめた（`s08-audit/r3{a,b,c,d}/round3b/`）。

| 観点 | 結論 | 内容 |
|---|---|---|
| (a) | supported | N1-iii は直った（注釈だけのファイルの CRLF・tab・CR だけ・末尾の改行なしを含めて NCBI とバイト一致）、N1-ii は直った（無作為の 400 ファイルで NCBI が失敗し LOSAT が 0 だった 18 件のうち 15 件が両方 1）。入力の読み方の組 2241 件：改善 2 件、一致から差 2 件（飛ばす行の後の配列の行の `;` や `!` の行：飛ばす行の無いファイルと同じに読むときに明示的に拒否、記録した後回し）。**残り：** N1-iv 空の定義行のレコード（`>` だけ）の後にもっともらしくない配列の行がある subject と空の query で、NCBI は 1、LOSAT は 0（`bio` の読み方による。`>` で始まるファイルでは S08 の前から同じ。低、S08+ の中身の無いレコードの項に記録した） |
| (b) | supported | N-b1 の 77 件はすべて NCBI とバイト一致（`r3-bin` は 0/77）。`5fa53b3f0` の 494 件で LOSAT が 0 で NCBI が失敗するものは 0（前は 5）。回帰の 3471 件：一致から差 5（すべて 1 つの入力：飛ばす行の後の、レコードの間の `;` の行。飛ばす行の無いファイルと同じ明示的な拒否、記録した後回し）、差から一致 90。新しい 450 件で、LOSAT が 0 で NCBI が失敗するもの・終了コードの違い・明示的でない失敗は 0 |
| (c) | supported | 第 3 回の 2613 行で LOSAT の stdout・stderr・終了コードの変化 0（拒否の 1962 行を含む）。新しい 456 件：注釈・空白だけの subject はどの query でも NCBI とバイト一致（57 行）、飛ばす行の後の title の警告は空の query で NCBI とバイト一致（58 行が差から一致）、無効な残基は読むときの明示的な拒否。**指摘：** L1b 飛ばす行の後の配列の行の中の空白と空の query（NCBI 0、LOSAT 1、明示的、記録した後回し「配列の行の中の空白」） |
| (d) | supported | L1 は直った（飛ばす行の後の無効な残基は、空の query でも読むときの明示的な拒否。ほかの subject ファイルと同じ）。飛ばす行の後の title の警告と注釈だけの subject は NCBI とバイト一致。`sub5fa` の 156 行で変わったのは意図した 1 行だけ。ABI と CLI：`abi_diff.js` の core・edge の 354 実行が第 3 回と同じ、`reg_vs_cli.js` で CLI が `register` に近づいた 2 行のほかは同じ。**指摘：** 注釈の行だけの subject を `register` は `bio` が読めない FASTA として拒否し、CLI は `Empty CBlastQueryVector`（低、記述の不足。どちらも失敗する。`abi_v2.md` の `register` の行に書いた） |

## アプリ側（S09）への注意

- TBLASTX は ABI v2 で形式 0・6・7 を出し、stream 1 に HSP の記録を出す（翻訳した行、SEG で mask した query の文字は小文字、両側の frame。`docs/web/abi_v2.md` §8）。
- `describe` の TBLASTX の `-max_target_seqs` に `default` が無くなった（省略すると hit list は 500 で、outfmt 0 は説明 500・整列 250。値を与えると両方がその値。clap の既定値では表せないので、help の文に「(default: 500)」と書いた）。アプリの検索画面の既定値の表示はこの help か、未指定として扱う。
- `parse` は `-outfmt` を挿入しなくなった。argv の誤りは、未知の program を含めて CLI と同じ文言になる。
- TBLASTX の register は BLASTN と同じに入力を確かめる（NCBI が違う読み方をする record は `not supported by LOSAT's TBLASTX` の誤り、`U` は `T`）。注釈の行（`!`・`#`・`;`）だけの subject は `register` が `bio` の読めない FASTA として拒否する（CLI は `Empty CBlastQueryVector`。`abi_v2.md`）。
- outfmt 0 の題が句読点だけの subject（NCBI が読み過ぎて落ちる）は、TBLASTX・TBLASTN の outfmt 0 で明示的に拒否する。保守者は例外 2 を広げた（DW-17）ので、S08+ で BLASTN と同じ代わりの題を出すようになる。

## 残件と引き継ぎ

- **`main` への PR**：S08 と S08b のコミット（`4fc67f9ab`・`2bcb86b1f` を含む）。merge は保守者。
- **保守者の確認**：V-PERF の `tblastx-multi`（上の「V-PERF」。推奨は NCBI の batch の費用として認め、batch の間で subject の準備を使い回す最適化を後の項目にする）。
- **SD（次のエンジン側のセッション、[指示書](../../losat_web_gui_sessions/session_sd_e2i_blastn_dc_megablast.md)、DW-18）**：保守者の依頼で、BLASTN の `-task dc-megablast` と `-task blastn-short` を NCBI と同じにする。
- **S08+**（[指示書](../../losat_web_gui_sessions/session_s08p_e2e_protein_options.md)の「S08 からの引き継ぎ」、SD の後）：最初の作業は、TBLASTX・TBLASTN の句読点だけの題の例外 2 の実装（DW-17）と TBLASTN の subject の定義行（監査 (c) の F14）。ほかに、第 3 回の監査の後回し（中身の無いレコードと `bio` の読み方（(a) N1-iv、(b) N-b2）、最初の定義行の前の NCBI が読む文字列と配列の行の中の空白・`;`（(c) L1・L1b、(b) N-r4-1）、垂直タブ、TBLASTN の `PRE_FETCH_SEQS_LIMIT`（(c) L3））、引数の解析の後の NCBI の検査の文言、TBLASTX の入力の読み方の安く移せる拒否、標準入力とパイプ、`-num_threads` 65535 以上、TBLASTX の `check_ncbi_application_settings`、BLASTN の HTML の題の検査を表示される subject に狭めること、BLASTP の差、TBLASTN の `-max_target_seqs` 5 未満の警告、batch の間の subject の準備の使い回し（V-PERF）。
- **SX**（[指示書](../../losat_web_gui_sessions/session_sx_blastx_integration.md)の「S08 からの引き継ぎ」）：BLASTX は共有の SEG・`BLAST_LargeGapSumE`・説明の一覧・`-seg` の値・`BLAST_Cutoffs` の呼び出しを、取り込むまで前のまま（DW-10）。入口の条件（LOSATX の v0.2.0 の認証が `main` に入ること）はまだ満たされていない。
- **アプリ側（S09）**：S08b と並行して 2026-10-03 に始めた（保守者の指示、[指示書](../../losat_web_gui_sessions/session_s09_w1_browser_runtime.md)の「S08b と並行して始める」）。TBLASTX の outfmt 0/7 の升目（26・28 行）はこの記録で `checked` になったので、取り込んで V-BR に加える。
