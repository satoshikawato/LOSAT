# LOSAT Web E2g（Session S07+++）ゲート記録

- 段階：E2g BLASTN の経路の棚卸しと一括の移植（[総合計画書](../../losat_web_gui_plan.md) §7 の S07+++、指示書 [S07+++](../../losat_web_gui_sessions/session_s07ppp_e2g_blastn_inventory.md)、計画 DW-12）
- ブランチ：`feature/losat-web-gui`。変更前は S07++ の記録（エンジンは `7693a9c73`、実行ファイルは `~/.cache/losat-web-gui-target/s07pg-native/release/LOSAT`、SHA-256 `a4ff5abb3a61…80f5`、ハッシュの一覧は `../losat_web_e2f/run-20260930T163727Z/artifacts.sha256`）
- 状態：**完了条件は未達（S07+++b に続く）**。S07+++ で済んだもの：CI の整備（指示書の 0.。毎晩の WASI の 1 件の TIMEOUT が残る）、棚卸しの第 1 段と第 2 段（`INVENTORY.tsv` 1006 行）、予備の hit list の単体試験（`2fabad16c`）。残るもの：一括の transpile（13 項目）、新しい明示的な拒否 2 つ、拒否の見直し 1 つ、試験・ゲート・V-PERF、独立監査（下の「残件と引き継ぎ」）

## CI の整備（指示書の 0.）

CI の設定とその検査の道具だけを変えた。エンジン（`LOSAT/src/`）は変えていない（コミット `db67e02cc`、`eb262bac2`）。

### 赤の原因

`wasm / integration` の「Frozen default and approved genetic-code regressions」の失敗（run [`36745543144`](https://github.com/satoshikawato/LOSAT/actions/runs/36745543144)、`7693a9c73` の push）は、ログ全文で `native/blastn/Sakai.MG1655.megablast` の出力が `ac8177764f35…`、凍結ハッシュが `3f27c1f1396b…` の不一致だった。`threaded/…`・`repeatability/…/2`・`…/3` の同じ case も PASS を出していない（最初の不一致の例外は、ほかの case が終わるまで表に出なかった）。これは S02 からの既知の差と同じ組である：

- S02 の基準（`../losat_web_e1a/baseline/hashes.tsv` の 12 行目）の LOSAT の出力が `ac8177764f35…`、Gate A の期待値が `3f27c1f1396b…`、`matches_expected` が false（`../losat_web_e1a/README.md` の「既存の差」）。
- LOSATX の Stage G の権威 v3（`../losatx_stage_g_authority_v3/run-20260928T121516Z/registry.json` の `sakai_new_gate_a` の 11・52・53、`approval` は `USER_APPROVED_EXACT_FOUR_PROPOSALS_2026-09-28`）が、新しい Sakai の期待値を `ac8177764f35…` として承認し、PR5 の Gate A の 11・52・53 を historical HARD_FAIL のまま残している（`GATE_STATUS.md`）。

差の中身（2026-10-02 に確かめた）：凍結した `3f27c1f1…` は v0.1.0 の LOSAT の出力で（`LOSAT-v0.1.0-candidate-20260829T051719Z` の実行ファイル、`ca61245` で再現）、NCBI と違っていた。今の LOSAT の `ac8177764f…` は、同じコマンドの NCBI BLAST+ 2.17.0 の出力とバイト一致し（`retained-linux-oracle` の実行）、登録済みの Linux の NCBI の指紋（`LOSAT/tests/ncbi_platform_variance_v010.json` の `retained_linux_raw_sha256`）とも同じである。v0.1.0 の出力との差は outfmt 7 の 6482 行のうち 5 行で、座標・e-value・bit score は同じで、gap の開始の数・長さ・mismatch・一致率だけが違う（megablast の greedy な traceback が、得点の等しい別の経路を選んでいた）。5 つの HSP は query の 48 kb・495 kb・4.66 Mb（2 つ）・4.94 Mb にあり、query の分割とは関係しない（Sakai は 5,498,578 文字で、megablast の分割の閾値 9,999,800 より短い）。v0.1.0 の後、S02 の前に `main` で直った。

よって回帰でもパリティの例外でもなく、古い凍結の期待値（NCBI と違っていた v0.1.0 の出力）との不一致で、新しい期待値が承認済みの別の版である。許可リストに載せた（凍結ハッシュは書き換えていない）。

### 変えたもの

| ファイル | 内容 |
|---|---|
| `LOSAT/tests/frozen_mismatch_allowlist.json`、`frozen_allowlist.py` | 既知の凍結ハッシュの不一致の許可リスト。項目は `Sakai.MG1655.megablast` の 1 つ（凍結 `3f27c1f1…`、許す出力 `ac817776…`、上の根拠）。許可するのは記録した出力だけで、別の出力は失敗。項目の case が凍結ハッシュに一致するようになったとき、または全件の実行でその case が実行されなかったときも失敗にして、リストの掃除を促す |
| `LOSAT/tests/check_wasm_threading_regressions.py` | 最初の不一致で止めず、全件を実行して `runs.json` と `summary.json` に集め、許可リストにない不一致・実行の失敗・スレッドの証拠の失敗・古い項目があれば最後に失敗する |
| `.github/workflows/ci.yml` | `push` は `main` だけ、PR は `pull_request`。`concurrency`（`${{ github.workflow }}-${{ github.ref }}`、`cancel-in-progress: true`）。`rust` のジョブに `Swatinem/rust-cache`、`cargo test` は `CARGO_PROFILE_TEST_OPT_LEVEL=1`（test の profile は dev の debug assertions と overflow checks を保つ）。`cargo build --release` は下の新しいジョブに移した。新しいジョブ「fast output regressions」（下） |
| `LOSAT/tests/ci_fast_regressions.py` | PR の速い検査（NCBI の実行ファイルは使わない）。変更したパスから program を選び（program のディレクトリはその program と、それを `crate::algorithm::` で読み込む program。それ以外のエンジン・fixture のパスは全 program。文書だけなら無し）、`../losat_web_e1a/capture_outputs.py` の case を S02 の基準（出力・正規化した stderr・終了コード・コマンド）と凍結ハッシュ（Gate A、TLOSAN Stage G。許可リストつき）で比べ、選んだ program の outfmt 0 の fixture（`../losat_web_e2a/check_losat.py`、1 と 4 スレッド）と、BLASTN を選んだときは下の BLASTN の fixture を確かめる。TBLASTX は 30 秒未満の 7 case だけ（ほかの 13 case は毎晩）。`--all-cases` で全件 |
| `LOSAT/tests/blastn_regression_fixtures.py`、`LOSAT/tests/fixtures/blastn_regression/` | NCBI BLAST+ 2.17.0 の BLASTN の出力を凍結した 29 件（下の「検出の確認」で足した）。入力は決定的に作ってコミットした（436 KB）。600 の subject で予備の hit list があふれる場合（既定、`-task blastn`、`-max_target_seqs` 1・2・3・5・6、`-subject_besthit`、`-max_hsps`、`-evalue`、outfmt 0/7、`-num_threads 4`、得点 3/−4 と 1/−1）、得点と e-value が等しい 600 の subject（subject の順の同点）、40 の query の batch、2 Mb の query の分割した batch（実行のときに `EDL933.fna` から切り出す）、IUPAC の曖昧な文字、小文字の mask。NCBI は `freeze` でだけ実行する。S07++ の実行ファイルは 29 件すべて一致 |
| `.github/workflows/nightly.yml` | 新設（`schedule` 毎日 18:17 UTC、`workflow_dispatch`）。全 case の出力（TBLASTX の長い case を含む）と、WASI の行列（`wasm-threading.yml`） |
| `.github/workflows/wasm-threading.yml` | `Swatinem/rust-cache`。TBLASTX のスレッドの閾値の検査を別の段階にし、凍結の段階が失敗しても実行する。CI の PR からは外した（毎晩、`workflow_dispatch`、release-readiness から呼ぶ） |
| `LOSAT/tests/test_ci_fast_regressions.py` | 許可リストの判定、基準と凍結の検査、パスからの選択、BLASTN の fixture の manifest の単体試験（`rust` のジョブで実行） |

### 実行時間（前と後）

| | 前（run `36745543144`、`7693a9c73` の push） | 後（PR [#109](https://github.com/satoshikawato/LOSAT/pull/109) の run [`36881993150`](https://github.com/satoshikawato/LOSAT/actions/runs/36881993150)、cache なし） |
|---|---|---|
| 1 つの変更で走る回数 | 2（`push` と `pull_request`） | 1（`pull_request`。`main` の push で 1） |
| `rust` | 14 分 38 秒（`cargo test` が 11 分 36 秒） | 3 分 43 秒（`cargo test` の lib の試験で失敗。通るときの残りは約 1 分） |
| `wasm / integration` | 1 時間 59 分 58 秒（失敗） | PR では走らない（毎晩） |
| fast output regressions | — | 4 分 22 秒（release のビルドと検査。失敗） |

手元（32 コア）で `cargo test --all-features` は、debug の profile の試験の実行が約 5 分 20 秒、`opt-level=1` ではビルドを含めて 1 分 55 秒（830 件すべて通過）。速い検査は 4 並列で 67 秒（223 case と fixture）。

### 検出の確認

1. まず指示書の例（`prelim_hitlist_size` を 1 つずらす：`get_prelim_hitlist_size` の `max(2h, 10)` を 9、`h + 50` を 49 にした実行ファイル）を手元で確かめた。capture の 216 case と outfmt 0 の fixture では**検出できなかった**（全 case が基準と一致）。予備の hit list の大きさを 1 つ変えても、出力が変わるのは、予備の段階の順位が 10 番目（または 550 番目）の subject が traceback の後に最終の上位に入る場合だけで、既存の case には subject が 11 以上ヒットするものが無い。そこで上の BLASTN の fixture を足したが、それでも 1 つずらした大きさは出力を変えなかった（29 件とも一致）。query の分割の重なりを 100 から 101 にしても同じだった。どちらも NCBI の設計上、出力に出にくい変更である。
2. 出力を変える変更として、予備の hit list の heap の同点の順（`s_EvalueCompareHSPLists` の subject の番号の比較、`hsp.rs`）を逆にしたブランチ `ci/e2g-detection-check` を push し、draft PR #109 を作った。run [`36881993150`](https://github.com/satoshikawato/LOSAT/actions/runs/36881993150) で、fast output regressions（[job](https://github.com/satoshikawato/LOSAT/actions/runs/36881993150/job/110435801031)）が outfmt 0 の fixture（1 と 4 スレッド）と BLASTN の fixture 14 件の差で失敗し、`rust`（[job](https://github.com/satoshikawato/LOSAT/actions/runs/36881993150/job/110435800623)）も単体試験 `word7_uses_subject_set_cutoffs` で失敗した。どちらも 5 分以内。PR は閉じ、ブランチは消した（変更は捨てた）。

### 毎晩の検査の確認

- 全 case（236）：手元で S07++ の実行ファイルの全 case を S02 の基準と比べ、差 0（`capture_outputs.py compare`）。`nightly.yml` は `main` に入るまで `workflow_dispatch` できない。
- WASI の行列：`wasm-threading.yml` を `feature/losat-web-gui` で `workflow_dispatch` した（run [`36882087513`](https://github.com/satoshikawato/LOSAT/actions/runs/36882087513)、`c6c3386cf`、2 時間）。凍結の 104 件のうち、Sakai の 4 つの label（native、threaded、repeatability の 2 つ）が許可リストどおり「allowed known mismatch」、TBLASTX の閾値の検査 9 件が PASS（凍結の段階の失敗の後も、別の段階として実行された）。**1 件が失敗した：`threaded/tblastx/p11_avclpv_psclpv` が TIMEOUT**（1 検索の期限 3600 秒。native は PASS）。前回の赤い run `36745543144` でもこの label は PASS を出していない（最初の不一致の例外に隠れていた）。最後に `wasm / integration` が緑だった run [`35201795582`](https://github.com/satoshikawato/LOSAT/actions/runs/35201795582)（2026-09-17、`implement-v010-distribution`）では、この label は 1 時間以内に PASS し、104 件すべてが通っていた。threaded Wasm の TBLASTX の速度の退行か、runner の差かは未確認（BLASTN の範囲の外。S07+++b の最初の作業）。手元では、Sakai の 2 case に絞った実行で、許可リストの 4 つの label が「ALLOWED」、ほかが PASS、終了コード 0 だった。

## 棚卸しの第 1 段

`inventory_refs.py`（このディレクトリ）が、BLASTN から届く Rust のファイルの「NCBI reference」の注釈を NCBI の関数に対応させる。

- 届くファイル：`LOSAT/src/algorithm/blastn/` の 26 ファイルから、`use` と `crate::`・`super::`・`self::` のパス（子のモジュールで始まるパスを含む）をたどり、ほかの program のディレクトリ（`algorithm/{blastp,blastx,tblastn,tblastx}`）には入らない。**54 ファイル**で、前回の数（54）を再現した。
- 注釈：「NCBI reference」を含む行。**1829**。前回の数 2180 は再現できなかった（前回の数え方は失われた。行の数え方を変えた 10 通り以上を試したが、2180 になるものは無かった。例：`file:line` の出現 2008、NCBI の file を引く注釈の行 2072）。この版の数を固定し、script の最初の検査にした（`EXPECTED`）。
- NCBI の関数：注釈の行（無ければ続く 3 行のコメント）の `<file>.<ext>:<行>` を、固定 commit の NCBI のファイルで、行の範囲の中の最初の関数の本体に対応させた。参照 1820、NCBI のファイル 111、**関数 390**（前回は 669。数え方が違う）、関数でない対象（struct、`#define`、file scope）74、解決できないファイル 0。
- 出力：`stage1_refs.tsv`（注釈ごと）、`stage1_functions.tsv`（関数ごとの参照の数と LOSAT の位置）、`stage1_summary.json`。

## 棚卸しの第 2 段

指示書の 7 つの範囲（A：app と引数、B：API、C：統計と DUST、D：lookup・scan・ungapped、E1：gapped、E2：traceback と hit の保存、F：整形）ごとに、読み取り専用の agent（sonnet）が `stage2/<範囲>.tsv` に行を追記しながら分類した（共通の指示は `stage2/PROMPT_COMMON.md`。並行は 4 つまで、2 回に分けた）。agent は NCBI のソースと LOSAT を読み、A・C・F の agent は小さな入力で NCBI BLAST+ 2.17.0 も実行して文言・終了コードを確かめた。オーケストレータ（Opus）が、divergent・unported・rejected の行をすべて NCBI のソースと LOSAT で確かめ、判断を `stage2/REVIEW.md` に書いた。

| 範囲 | 行 | faithful | n/a | rejected | divergent | unported |
|---|---|---|---|---|---|---|
| A | 236 | 118 | 62 | 42 | 6 | 8 |
| B | 175 | 127 | 41 | 5 | 1 | 1 |
| C | 174 | 130 | 36 | 8 | 0 | 0 |
| D | 88 | 64 | 17 | 0 | 7 | 0 |
| E1 | 68 | 52 | 11 | 1 | 3 | 1 |
| E2 | 119 | 85 | 33 | 0 | 1 | 0 |
| F | 146 | 98 | 29 | 10 | 1 | 8 |
| 計 | 1006 | 674 | 229 | 66 | 19 | 18 |

`INVENTORY.tsv`（`build_inventory.py` が 7 つの表と判断の規則から作る）は、行ごとに id、NCBI のファイル・行・関数・分岐、LOSAT の位置、状態、影響、`e2g_action`（faithful と n/a 以外の行の扱い）、`e2g_result`（移植の結果。S07+++b で埋める）を持つ。影響が high の差は無かった。

| 扱い | 行 | 内容（`build_inventory.py` の `ACTIONS`） |
|---|---|---|
| T1〜T13（transpile） | 25 | T1 初期の hit の並べ替え（連結した query の位置、安定な並べ替え）、T2 `BlastGetStartForGappedAlignmentNucl` の `Int4`、T3 予備の e-value の刈り込みの比較、T4 traceback の順（heap の e-value の順を保つ）、T5 小さい・標準の lookup table の cell の query の位置の昇順、T6 query の塊が 8000 以下なら対角線の配列、T7 `BATCH_SIZE`・`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`、T8 最初の batch が無効な query だけのときの得点の表の失敗、T9 batch ごとの警告、T10 `-subject` が無いときの誤り、T11 `CTOOLKIT_COMPATIBLE`、T12 `PRE_FETCH_SEQS_LIMIT` の誤り、T13 Lambda・K・H の `%#8.3g` |
| R1、R2（新しい明示的な拒否） | 11 | `BL2SEQ_LEGACY`（別の検索の mode と報告）、NCBI の診断と registry の環境変数（`DIAG_POST_LEVEL` など）と、出力を変える `.ncbirc` の設定（2026-10-02 の決定） |
| V1（見直し） | 1 | reward − penalty が 3000 を超える得点の拒否（NCBI と比べて一致すれば外す） |
| K-*（拒否のまま） | 60 | FASTA の読み方（S17 の前のセッション SF で移植、DW-13）21、NCBI の C++ の層か範囲外の機能 20、別の mode 15、NCBI が落ちる 2、NCBI の 32 ビットの溢れ 2 |
| X-*（承認済みの例外、2026-10-02） | 5 | `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`：引数の構文の誤りと `-help`（X-cli 2）、`-subject` での `-num_threads` の 2 つの警告（X-threads 2）、メモリ不足（X-oom 1） |
| T14、R3（2026-10-02 の決定で足した） | 2 | T14 outfmt 0 の書き込みの失敗（`BLAST failed to write output`、終了コード 6。outfmt 6/7 は例外）、R3 UTF-8 でないファイル名の明示的な拒否 |

確かめたことの要点（詳細は `stage2/REVIEW.md`）：NCBI の traceback は `_DEBUG` のときだけ得点で並べ直す（オラクルは release）。NCBI は subject の chunk ごとに得点で並べる（blast_engine.c:555）ので、LOSAT の並べ直しが効くのは heap にした（550 を超える subject の）list だけ。小さい lookup table の cell は NCBI では query の位置の昇順、megablast の table は新しい順で、LOSAT は小さい table を megablast の鎖で持つ。NCBI は query の塊が 8000 以下なら query の数によらず対角線の配列を使う（blast_parameters.c:225-231）。`BL2SEQ_LEGACY` は `CLocalBlast` の dbscan mode を切る（local_blast.cpp:189,289）。`.ncbirc` の `[BLAST] LONG_SEQID` も出力を変えるが、LOSAT は `.ncbirc` を読まない（記録だけ）。

## 予備の hit list の単体試験

S07++ の独立監査の第 2 回の指摘（重大度 低）を `2fabad16c` で入れた（試験だけ。検索のコードは変えていない）：subject の番号による同点（等しい list は番号の大きい方が残り、先に来る）、1e-180 の規則、最初の HSP の得点、多くの list からの heap の選択、heap にするときの list の e-value の並べ替え、空の list の除去と最後の大きさの上限、collector（query ごとの `prelim_hitlist_size`）、`merge_prelim_hit_list`（`Blast_HitListMerge`）。`cargo fmt --check`、`clippy --all-targets --all-features -D warnings`、`cargo test --all-features`（lib 626 件、ほか全件）、`verify_refs.py` が通過。

## 残件と引き継ぎ（S07+++b へ）

- **毎晩の WASI の TIMEOUT：** `threaded/tblastx/p11_avclpv_psclpv`（上の「毎晩の検査の確認」）。S07+++b の最初に、速度の退行か runner の差かを確かめ、毎晩の検査を緑にする。
- **一括の transpile：** T1〜T13、R1、R2、V1（上の表と `stage2/REVIEW.md`）。エンジンは変えていない（`2fabad16c` は試験だけ）。
- **試験・ゲート・V-PERF・独立監査：** 指示書の 4.。移植した分岐ごとの NCBI との比較と、BLASTN の fixture（`blastn_regression_fixtures.py`）への case の追加。
- **CI の残り：** 緑の PR の実行時間（cache あり）を測る。`main` への PR（この記録の後）。
- **保守者の決定（2026-10-02、DW-13、`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`、AGENTS.md）：** 判断待ちだった 8 件を、保守者が推奨の案で決めた。
  - 承認済みの例外：引数の構文の誤りと `-help`、`-subject` での `-num_threads`（並列のまま、2 つの警告を出さない）、outfmt 6/7 の書き込みの失敗、メモリ不足。
  - S07+++b で移植：outfmt 0 の書き込みの失敗（T14）。S07+++b で明示的に拒否：UTF-8 でないファイル名（R3）、出力を変える `.ncbirc` の設定（R2 に含める）。
  - FASTA の読み方：S17 の前の専用のセッション SF で移植（それまで拒否のまま）。
  - E1a〜E1c の V-PERF と E1c の CLI の 2 つの差：記録どおり承認。
- 次の指示書：[`session_s07pppb_e2g_transpile.md`](../../losat_web_gui_sessions/session_s07pppb_e2g_transpile.md)。
