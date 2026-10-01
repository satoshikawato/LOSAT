# LOSAT Web E2g（Session S07+++）ゲート記録

- 段階：E2g BLASTN の経路の棚卸しと一括の移植（[総合計画書](../../losat_web_gui_plan.md) §7 の S07+++、指示書 [S07+++](../../losat_web_gui_sessions/session_s07ppp_e2g_blastn_inventory.md)、計画 DW-12）
- ブランチ：`feature/losat-web-gui`。変更前は S07++ の記録（エンジンは `7693a9c73`、実行ファイルは `~/.cache/losat-web-gui-target/s07pg-native/release/LOSAT`、SHA-256 `a4ff5abb3a61…80f5`、ハッシュの一覧は `../losat_web_e2f/run-20260930T163727Z/artifacts.sha256`）
- 状態：**進行中**。CI の整備（指示書の 0.）は済んだ。棚卸しの第 1 段は済み、第 2 段の分類を実行している

## CI の整備（指示書の 0.）

CI の設定とその検査の道具だけを変えた。エンジン（`LOSAT/src/`）は変えていない（コミット `db67e02cc`、`eb262bac2`）。

### 赤の原因

`wasm / integration` の「Frozen default and approved genetic-code regressions」の失敗（run [`36745543144`](https://github.com/satoshikawato/LOSAT/actions/runs/36745543144)、`7693a9c73` の push）は、ログ全文で `native/blastn/Sakai.MG1655.megablast` の出力が `ac8177764f35…`、凍結ハッシュが `3f27c1f1396b…` の不一致だった。`threaded/…`・`repeatability/…/2`・`…/3` の同じ case も PASS を出していない（最初の不一致の例外は、ほかの case が終わるまで表に出なかった）。これは S02 からの既知の差と同じ組である：

- S02 の基準（`../losat_web_e1a/baseline/hashes.tsv` の 12 行目）の LOSAT の出力が `ac8177764f35…`、Gate A の期待値が `3f27c1f1396b…`、`matches_expected` が false（`../losat_web_e1a/README.md` の「既存の差」）。
- LOSATX の Stage G の権威 v3（`../losatx_stage_g_authority_v3/run-20260928T121516Z/registry.json` の `sakai_new_gate_a` の 11・52・53、`approval` は `USER_APPROVED_EXACT_FOUR_PROPOSALS_2026-09-28`）が、新しい Sakai の期待値を `ac8177764f35…` として承認し、PR5 の Gate A の 11・52・53 を historical HARD_FAIL のまま残している（`GATE_STATUS.md`）。

よって回帰ではなく、承認済みの差として許可リストに載せた（凍結ハッシュは書き換えていない）。

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
- WASI の行列：`wasm-threading.yml` を `feature/losat-web-gui` で `workflow_dispatch` した（run [`36882087513`](https://github.com/satoshikawato/LOSAT/actions/runs/36882087513)）。結果は下の「残件」に書く。手元では、Sakai の 2 case に絞った実行で、許可リストの 4 つの label（native、threaded、repeatability の 2 つ）が「ALLOWED」、ほかが PASS、終了コード 0 だった。

## 棚卸しの第 1 段

`inventory_refs.py`（このディレクトリ）が、BLASTN から届く Rust のファイルの「NCBI reference」の注釈を NCBI の関数に対応させる。

- 届くファイル：`LOSAT/src/algorithm/blastn/` の 26 ファイルから、`use` と `crate::`・`super::`・`self::` のパス（子のモジュールで始まるパスを含む）をたどり、ほかの program のディレクトリ（`algorithm/{blastp,blastx,tblastn,tblastx}`）には入らない。**54 ファイル**で、前回の数（54）を再現した。
- 注釈：「NCBI reference」を含む行。**1829**。前回の数 2180 は再現できなかった（前回の数え方は失われた。行の数え方を変えた 10 通り以上を試したが、2180 になるものは無かった。例：`file:line` の出現 2008、NCBI の file を引く注釈の行 2072）。この版の数を固定し、script の最初の検査にした（`EXPECTED`）。
- NCBI の関数：注釈の行（無ければ続く 3 行のコメント）の `<file>.<ext>:<行>` を、固定 commit の NCBI のファイルで、行の範囲の中の最初の関数の本体に対応させた。参照 1820、NCBI のファイル 111、**関数 390**（前回は 669。数え方が違う）、関数でない対象（struct、`#define`、file scope）74、解決できないファイル 0。
- 出力：`stage1_refs.tsv`（注釈ごと）、`stage1_functions.tsv`（関数ごとの参照の数と LOSAT の位置）、`stage1_summary.json`。

## 棚卸しの第 2 段

指示書の 7 つの範囲（A：app と引数、B：API、C：統計と DUST、D：lookup・scan・ungapped、E1：gapped、E2：traceback と hit の保存、F：整形）ごとに、読み取り専用の agent（sonnet）が `stage2/<範囲>.tsv` に行を追記しながら分類する。共通の指示は `stage2/PROMPT_COMMON.md`。

（実行中）

## 残件

- WASI の行列の `workflow_dispatch`（run `36882087513`）の結果。
- 緑の PR の実行時間（cache あり）は、この記録の後の `main` への PR で測る。
