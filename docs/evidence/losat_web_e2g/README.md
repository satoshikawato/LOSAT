# LOSAT Web E2g（Session S07+++・S07+++b）ゲート記録

- 段階：E2g BLASTN の経路の棚卸しと一括の移植（[総合計画書](../../losat_web_gui_plan.md) §7 の S07+++、指示書 [S07+++](../../losat_web_gui_sessions/session_s07ppp_e2g_blastn_inventory.md) と [S07+++b](../../losat_web_gui_sessions/session_s07pppb_e2g_transpile.md)、計画 DW-12・DW-13・DW-14・DW-15）
- ブランチ：`feature/losat-web-gui`。変更前は S07++ の記録（エンジンは `7693a9c73`、実行ファイルは `~/.cache/losat-web-gui-target/s07pg-native/release/LOSAT`、SHA-256 `a4ff5abb3a61…80f5`、ハッシュの一覧は `../losat_web_e2f/run-20260930T163727Z/artifacts.sha256`）。変更後は S07+++b の最後のエンジン `f4057718a`（native の SHA-256 `331fcba36447…44d8`、[`run-20261002T042123Z/artifacts.sha256`](run-20261002T042123Z/artifacts.sha256)）
- 状態：**完了（2026-10-02、S07+++b）**。S07+++ で CI の整備と棚卸し、S07+++b で一括の transpile、保守者の判断（NCBI の不具合の扱い）、試験・ゲート・V-PERF、独立監査（第 4 回で supported）。下の「完了条件」の表

## 完了条件（計画 §7 の S07+++ の行、S07+++b の指示書）

| 条件 | 結果 | 根拠 |
|---|---|---|
| 棚卸しの表が経路を網羅する | 満たす | `INVENTORY.tsv` 1015 行（監査で 9 行を足した）。監査 (a) が NCBI を callgrind の下で 235 回実行し、届いた 4410 の定義と突き合わせた。出力に効く漏れ F1 は例外 1 に |
| 未移植と差のある移植が残らない（残すものは明示的な拒否か承認済みの例外） | 満たす | 状態が faithful・n/a でない 104 行：移植して faithful 25、移植（`-evalue`）1、承認済みの例外 7、明示的な拒否 71（理由は監査 (c) と第 3 回で確かめた）。下の「一括の transpile」「保守者の判断」 |
| S07+ と S07++ の全検査、既存の BLASTN のゲート、fixture に退行なし | 満たす | 下の「ゲート」（`run-20261002T042123Z`、capture 236 件が S02 の基準と差 0、fixture 107 件、sweep 3710 件ほか） |
| V-PERF の非退行 | 満たす | 下の「V-PERF」 |
| 棚卸しの表を基準にした独立監査 | 満たす | 下の「独立監査」（第 1〜3 回の指摘を直し、第 4 回で supported） |
| 毎晩の WASI の TIMEOUT の解消（S07+++b の 0.） | 満たす（原因と変更）。`main` での dispatch の結果は下 | runner の差。期限を 7200 秒（`beea0bddf`） |

# S07+++（CI の整備と棚卸し）

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

## S07+++ の残件と引き継ぎ（S07+++b へ。この記録の時点）

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

# S07+++b（一括の transpile、保守者の判断、試験、独立監査）

## CI の残り（S07+++b、指示書の 0.）

### 毎晩の WASI の TIMEOUT

sonnet の agent が CI の 3 つの run の記録（`gh run view --log` と、成果物 `wasi-threading-evidence` の `frozen-regressions/runs.json`）を集め、手元で threaded Wasm の p11 を測った（`~/.cache/losat-web-gui-target/e2g-wasi-timeout/RESULTS.md`）。

- 最後に緑だった run `35201795582`（2026-09-17、`6bfb1b09b`）でも、`threaded/tblastx/p11_avclpv_psclpv` は 3358 秒（CPU 6617 秒）で、期限 3600 秒の 7% 手前だった。赤い 2 つの run（`36745543144`、`36882087513`。同じ wasm と native）では、どちらも 3600 秒で TIMEOUT。
- 赤い run では、ほぼすべての case（native も threaded も）が緑の run より 1.1〜1.9 倍遅い（同じ実行ファイルの 2 つの赤い run の間でも native p11 が 1413 秒と 2060 秒）。`--jobs 3`（4 vCPU の runner で 3 つの検索が並ぶ）は 3 つの run とも同じで、threaded p11 の間は常にほかの 2 つの検索が走っていた。
- 手元の 4 コアで、CI と同じ argv・環境（`RAYON_NUM_THREADS=1`、`-num_threads 4`）で p11 を測ると、HEAD 相当の threaded Wasm が 2164 秒（CPU 6327 秒）、緑の run の成果物が 2809 秒（CPU 8000 秒）。TBLASTX の速度の退行は無い（`ci_wasi_timeout/`）。

結論：runner の差（共有の 4 vCPU の上で 3 つの検索が並び、threaded p11 だけ CPU を多く使う）。退行ではないので S08 の最初の作業にはしない。`wasm-threading.yml` の検索ごとの期限を 3600 秒から 7200 秒にした（`beea0bddf`、CI だけのコミット）。

### 緑の PR の実行時間（cache あり）

PR [#110](https://github.com/satoshikawato/LOSAT/pull/110)。最初の run `36942594965`（`baf180fbb`）は `rust` のジョブが pure-Rust の境界の検査で失敗した（T14 で足した SIGPIPE の `extern "C"` の import。`48a9ee0c2` で外し、閉じたパイプは保守者の判断で承認済みの例外 5、DW-14）。直した後の run [`36943546622`](https://github.com/satoshikawato/LOSAT/actions/runs/36943546622)（`48a9ee0c2`、cache あり）が緑：全体 4 分 44 秒、fast output regressions 3 分 11 秒、`rust` 4 分 39 秒（S07+++ の cache なしの run は 4〜5 分と 4.5 分）。20 分の目標の内。

### `main` への PR

PR #110（S07+++ の CI の変更 `db67e02cc` 以降と、このセッションの成果）。最後の push（`a9f604ddc`）の CI は緑（run `36964410624`：`rust` 2 分 42 秒、fast output regressions 2 分 58 秒）。merge と毎晩の検査は下の「`main` への merge と毎晩の検査」。

## 一括の transpile（S07+++b、指示書の 1.）

設計と実装はオーケストレータ（Opus）が行い、エンジンのソースは 1 項目ずつ変えた。各項目の後に、`verify_refs.py`、`cargo fmt --check`、`clippy --all-targets --all-features -D warnings`、`cargo test --all-features`（`opt-level=1`）、`ci_fast_regressions.py`（capture の S02 の基準、凍結ハッシュ、outfmt 0 の fixture、BLASTN の fixture）を通した（記録は `~/.cache/losat-web-gui-target/e2g-items/<項目>/`）。途中から CI の `rust` のジョブと同じ `check_pure_rust_runtime_boundary.py` も通す。NCBI との比較は、項目ごとに `LOSAT/tests/blastn_regression_fixtures.py` の `CASES` に足して NCBI 2.17.0 で `freeze` し（既存の凍結出力は変わらないことを毎回確かめた）、変更前の実行ファイルでも `check` して、その case がその分岐を区別するかを記録した。fixture は 29 件から 84 件になった（case ごとの環境変数の列 `env` を足した）。

| 項目 | コミット | 内容 | NCBI との比較 |
|---|---|---|---|
| T3 | `35e3e6d7f` | 予備と最後の e-value の刈り込みを `!(evalue > cutoff)` に（blast_engine.c:662-664、blast_hits.c:1996-1998） | 単体試験（NaN だけが違う。NaN の `-evalue` は拒否） その後 `30713884f` で NaN・無限大の `-evalue` を移植したので、NaN の分岐にも届く（下の「保守者の判断」#10）。 |
| T13 | `a99527f20` | Lambda・K・H を C の `%#8.3g` どおりに（align_format_util.cpp:584-603）。全 program の outfmt 0 が共有 | 単体試験（C の文字列 32 値：0.99996、9.9996、999.5、1.23e-5、9.999e-5、99.95 など） |
| T11 | `e4b5c4a63` | `CTOOLKIT_COMPATIBLE`（値が空でも）で説明の一覧の見出しを「(bits)」に（showdefline.cpp:80） | fixture `ctoolkit.*` 2 件（変更前は差）。全 program の outfmt 0 の 55 件 × 3 環境：51 件一致、TBLASTX の 4 件は LOSAT に outfmt 0 が無い（S08） |
| T2 | `b5940fba9`、fixture `de42604ff` | `BlastGetStartForGappedAlignmentNucl` を `Int4` で（blast_gapalign.c:3327-3389） | overflow を検査する build で 2762 の BLASTN のコマンドを走らせ（`overflow_hunt/`）、見つかった overflow はこの 1 か所だけ（`-task blastn -word_size 4` の 188 件）。188 件とも NCBI・変更前・変更後が一致（release の折り返しが同じ開始点を選んでいた）。word_size_sweep 256 と slice_sweep 7 種 670 の計 926 の組合せで差 0。fixture `t2.word4` |
| T1 | `019ba8ea4` | 初期の hit の並べ替えを、連結した query の位置と安定な並べ替えに（blast_extend.c:274-310、na_ungapped.c:206） | fixture `pal.*`（X + revcomp(X) の 40 query）。変更前も一致（両鎖の同点は後の並べ替えで隠れる） |
| T6 | `0b3b851c7` | query の塊が 8000 以下なら query の数によらず対角線の配列（blast_parameters.c:225-231、blast_engine.c:1002-1003） | `query->length` と LOSAT の `query_concat_length` が同じこと、両方の scan の経路が連結した位置で対角線を作ることを確かめた。fixture `sq.*`（12 query、塊 5155）。変更前も一致 |
| T5 | `435f97afd` | 小さい・標準の lookup table の cell を、索引を作った順（query の位置の昇順）に（blast_lookup.c:74-76、blast_nalookup.c:289-293・522-526）。megablast の table は新しい順のまま（`TaskConfig::mb_lookup`） | 単体試験。fixture `rep.*`・`rep2.*`（塊が 8000 を超え、小さい table、繰り返し）。変更前も一致 |
| T4 | `7eb018f70` | traceback の前の得点の並べ直しを除き、保存した順で追う（blast_traceback.c:358-365、blast_hits.c:3272-3284） | BLASTN では 1 つの list の e-value の順は得点の順と同じ（両鎖は同じ Karlin block と探索空間、1e-180 の規則は得点に戻る）。subject のすべての query で 1 つの区間木は、query ごとの木と同じ（連結した位置は重ならない）。fixture `prelim.gaps10_*`、`ties.gaps10_1_1`。変更前も一致 |
| T7 | `8ccf079ba` | `BATCH_SIZE`・`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE` を NCBI と同じく読む（`StringToInt`、`size_t`・`TSeqPos`・`Int4` の計算）。E2c §G の 2 つの拒否がなくなる。拒否のまま：整数でない値（NCBI の `CStringException` の文言はオラクルの build のソースのパスを含む。終了コード 255）、`BATCH_SIZE` なしの `CHUNK_SIZE=1000`（最初の batch が空になる。「BLAST engine error: Empty CBlastQueryVector」、終了コード 3）、負の `CHUNK_SIZE` と負の `OVERLAP_CHUNK_SIZE` の組 | fixture `env.*` 20 件（`CHUNK_SIZE=2000`・`500` は既定と違う出力になり、一致）。batch_sweep 22 回 × 60 case と split_check 4 回を変数つきで：差 0 その後：`5dac71f72` で `CHUNK_SIZE=1000` を NCBI の誤りの再現に（監査 (c)）、`993310891` で塊をもう一度分ける設定を拒否（監査 (a)）→ `d846e9bbe` で承認済みの例外 1、`ad5fa9c85` で負の組の拒否を batch を分けるときだけに（監査 (b) D2）。 |
| R1 | `7b63980b2` | `BL2SEQ_LEGACY` を明示的に拒否（値が空でも。blast_app_util.cpp:204-210） | 空の query では NCBI と同じく「Query is Empty!」で終わる |
| T12 | `650b02771` | `PRE_FETCH_SEQS_LIMIT` の整数は受け付け（出力は変わらない）、整数でない値（空を含む）は明示的に拒否（blast_app_util.cpp:732-737） | 0・5・2147483647 で NCBI と一致。fixture `env.prefetch*` |
| T10 | `1a0fd98c1` | `-subject` が無いとき NCBI の誤り（「BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified」、終了コード 1）を、`-query`・`-out` を開く前に（blast_args.cpp:2558-2562）。アダプタの `describe` は `-subject` を出さないので変わらない | fixture `nosubject.*`（変更前は clap の終了コード 2） |
| R2 | `219c2c49e` | NCBI の診断と registry の設定を明示的に拒否：`DIAG_*`、`NCBI_CONFIG_*`（出力に効かない項目を除く）、`ABORT_ON_THROW`、stack trace と `LOG_*` の parameter、Boolean でない `BLAST_USAGE_REPORT`。`blastn.ini` と `.ncbirc` を NCBI の探す場所と順序で探し（`NCBI_CONFIG_PATH`、`.`、`$HOME`、`$NCBI`、`/etc`、実行ファイルの場所。最初の 1 つ）、出力に効かない項目（`[BLAST] BLASTDB` など）だけなら受け付ける（DW-13）。CLI の blastn の入口で検査する（アダプタには NCBI の application の層が無い） | 調べは sonnet の agent（`~/.cache/losat-web-gui-target/e2g-r2r3/FINDINGS.md`、オラクルの実行約 300）。fixture `ncbirc.harmless_fmt0`（出力に効かない `.ncbirc` のある `HOME` で NCBI の出力は変わらず、LOSAT も一致） |
| R3 | `6182cef73` | UTF-8 でない `-subject`・`-query`・`-out` を、NCBI がそのファイルを開く順（subject、query、out）で明示的に拒否 | NCBI は名前のバイトをそのまま `Database:` の行と誤りの文言に書く（オラクル） |
| T14 | `5cd9cc3cf`、`48a9ee0c2` | outfmt 0 の書き込みの失敗で「BLAST failed to write output」、終了コード 6（blast_format.cpp:118-119、blast_app_util.hpp:252-255）。outfmt 6/7 は承認済みの例外のまま | `-out /dev/full` で NCBI と一致（fixture `write.devfull_fmt0`）。閉じたパイプ：下の「保守者の判断待ち」 閉じたパイプは保守者の判断で承認済みの例外 5（DW-14）。SIGPIPE の既定への戻しは pure-Rust の境界の検査が拒否するので外した（`48a9ee0c2`）。その後 `e37099f44`・`f4057718a` で、書き込みの失敗が警告より前に（NCBI の prolog の最初の flush と同じ位置で）止まる。 |
| T8・T9 | `baf180fbb` | batch ごとに、読むときに title の警告、報告とともに無効な query の警告（blastn_app.cpp:277-318、blast_format.cpp:1450-1452）。得点の表が無い得点で、無効な query だけの batch の後の batch が表の誤りに当たったら、前の batch の報告を書き（epilog なし：outfmt 0 の末尾、outfmt 7 の最後の行、最後の query の後の空行を書かない）、NCBI の誤り（終了コード 3）を出す。失敗する batch に無効な query がある分岐（NCBI が落ちる）は拒否のまま | fixture `kaerror.later_batch.fmt{0,6,7}`、`warnings.batches.fmt{0,6}`、`warnings.batch1000.fmt7`（6 件とも変更前は差）。E2c の `scoring_error.first_batch_invalid` は「same-error」になった その後 `e37099f44` で、警告を query の報告の間に NCBI の順で書く（stdout と stderr を 1 つにしたとき。監査 (b) D3）。 |
| V1 | `82e8c7593` | reward − penalty が 3000 を超える拒否を外し、reward 32767・penalty −32768 だけを拒否（`BlastScoreBlkMaxScoreSet` は `BLAST_SCORE_MAX`・`MIN` を得点の範囲から外し、`BlastScoreFreqCalc` は reward 32767 を配列の外に数える。blast_stat.c:1506-1520、2175-2179） | 制限を外した実行ファイルで、NCBI が実行するすべての case が一致（最大公約数で表に当たる組、24000/−30000 まで、表の gap・0/0・表を超える gap、両 task、outfmt 0/6）。LOSAT は最大 0.02 秒・8 MB（NCBI 77 MB）。差は megablast の gap が 32767 を超える既存の拒否と、境界の 2 値だけ。fixture `scores.*` 4 件（変更前は拒否） |

棚卸しの表：`build_inventory.py` の `RESULTS` に各項目の結果を書き、`INVENTORY.tsv` を作り直した（`c318d038c`、`81eb47efd`）。`inventory_refs.py` の固定した数は、最後に 54 ファイル、1881 の注釈、398 の関数（変更前は 1829 と 390。監査の修正と保守者の判断の後の数、`a9f604ddc`）。棚卸しの表は 1015 行（監査の後に 9 行を足した）。

## 保守者の判断：NCBI の不具合に当たる挙動（DW-15、`PD-LOSAT-NCBI-DEFECTS` 版 1.0）

監査 (a) の F1（query の塊を NCBI がもう一度分ける設定）について、保守者は「妥当な結果が出るならば承認済みの例外」とした。その後、ほかにも NCBI 側の不具合に近い挙動がないかまとめて判断を求められたので、棚卸しの表と監査の記録から 10 項目を出し、2026-10-02 に次のとおり決まった（#8・#9 は一度「意図どおりの結果を出す例外」とされ、同日に「再現のまま」に改められた）。

| # | 入力 | NCBI | 決定 | 結果 |
|---|---|---|---|---|
| F1 | 塊をもう一度分ける `CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE` | CCoreException（null pointer）、終了コード 3 | 承認済みの例外 1 | `d846e9bbe`。LOSAT は塊を 1 回ずつ検索 |
| 1 | 句読点だけで終わりに空白と区切りが続く subject の定義行（outfmt 0） | 題の整形が文字列の外を読み、hit があると SIGSEGV | 承認済みの例外 2 | `c72452236`。整形を文字列の終わりで止める（`, ,` は `, `） |
| 2 | K-A の表に無い得点系で、batch の先頭に無効な query | gapped の block なしで進み SIGSEGV | 明示的な拒否のまま | 変更なし |
| 3 | `-max_target_seqs` が 2^31−51 を超える | 予備の hit list の大きさが負になり落ちる | 明示的な拒否のまま | 変更なし |
| 4 | reward 32767・penalty −32768 | 頻度の配列の外に数える／全 query が無効 | 明示的な拒否のまま | 変更なし（reward 32768 以上が 16 bit で 0 以下に回る場合も同じ扱い。確認待ち） |
| 5 | subject の合計 2^31 文字以上 | 32 bit の合計が回る | 明示的な拒否のまま | 変更なし |
| 6 | megablast の gap cost 32767 超 | greedy の距離の 32 bit の計算が回りうる | 明示的な拒否のまま | 変更なし |
| 7 | 負の `CHUNK_SIZE` が負の `OVERLAP_CHUNK_SIZE` より大きい | `size_t` で回る | 明示的な拒否のまま | `ad5fa9c85`。監査で、batch を分けないとき（例 −1 と −2147483648）は NCBI が普通に検索すると分かったので、拒否を batch を分ける場合に絞った（分けると NCBI は CCoreException か、隙間のある塊で検索して hit を落とす） |
| 8 | `BATCH_SIZE` なしの `CHUNK_SIZE=1000` | batch の大きさが 0、prolog の後「Empty CBlastQueryVector」、終了コード 3 | 再現のまま | 変更なし |
| 9 | `-max_target_seqs` 2^30〜2^31−51 | 予備の hit list の大きさが回って 10 | 再現のまま | 変更なし |
| 10 | `-evalue` の `+inf`・`-nan`・`+nan(1)`・`1e999` | 受け付け、最大の e-value と同じく検索 | NCBI どおり移植 | `30713884f` |

検証：

- 例外 1（`docs/evidence/losat_web_e2g/resplit/`）：30〜60 kb の 3 query × blastn・megablast × 塊と重なりの 7 組 = 42 構成。NCBI はすべて終了コード 3。LOSAT の出力は、35 構成で NCBI の分けない検索と一致し、84 の実行（outfmt 6・0）のうち 74 で NCBI のもう一度分けない最大の重なりの出力と一致。残りの 10 は塊の境の HSP の差だけで、NCBI 自身の分割でも同じ種類の HSP が落ちる。fixture `env.resplit_*`（新しい列 `oracle_env`：NCBI はその最大の重なりで凍結）。
- 例外 2（`punct_defline/`）：hit のある subject に 8 つの句読点の定義行、megablast・blastn、`-max_target_seqs 3` の有無の 4 実行。NCBI はすべて SIGSEGV。LOSAT の報告は、定義行を同じ長さの代わりの文字列にした NCBI の報告で代わりの文字列を LOSAT の題に置き換えたものと一致。hit のない句読点の subject（NCBI が動く）は fixture `punct.nohit_*` で NCBI と一致（変更前は拒否）。
- 移植 #10：fixture `evalue.*` 4 件、走査 `e2g-sweep-AC5` 440 件（`+inf`・`-nan`・`+nan(1)`・`1e999`、megablast・blastn・word 7・`-subject_besthit`・IUPAC）で差 0。`BLAST_Cutoffs` は INT_MIN／`e > 0` の偽で cutoff 1 のまま（単体試験）。web API も CLI と同じ読み方（裸の `inf` は NCBI と同じく引数の誤り）。

## 監査 (b) の指摘への対応

| 指摘 | 内容 | 対応 |
|---|---|---|
| D1 | 大きい重なりで NCBI が落ち、LOSAT が検索 | F1 と同じ。承認済みの例外 1 |
| D2 | 負の組の拒否が広すぎる | `ad5fa9c85`（上の #7） |
| D3 | stdout と stderr を 1 つにしたとき（`2>&1`）、LOSAT は警告をすべて先に出す。NCBI は batch・query ごとに挟む | `e37099f44`。NCBI は警告を cerr に出し、cerr は cout に tie されているので、警告の前に書いた分の stdout が出る。LOSAT は query ごとの警告を保持し、報告の writer が各 query の前で報告を flush してから書く（`QueryWarnings`）。outfmt 0 の `"\n\n"` の前置きは警告の後 |
| D4 | outfmt 0 の書き込みの失敗で、LOSAT は警告を先に出す | `e37099f44`。NCBI の最初の flush（版の行の後）で失敗し、query を読む前に止まる。LOSAT は警告の前の flush で止まる |
| D5 | 閉じたパイプ | 承認済みの例外 5（DW-14） |

D3・D4 の検証（`order/`）：19 の入力 × outfmt 0/6/7 × 1 つの流れ・別々・stdout を `/dev/full`・`-out /dev/full` × `BATCH_SIZE` なし・1・200 の 456 実行がすべて NCBI と一致（変更前の実行ファイルは 229 が差）。NCBI の flush の位置は gdb で write の system call を数えて確かめた（`flushmap.py`）。fixture `*.merged`（新しい形：stderr を stdout に合わせる）5 件と `write.devfull_warnings_fmt0`。

## ゲート（指示書の 2.）

run [`run-20261002T042123Z/`](run-20261002T042123Z/)（`head.txt` = `f4057718a`、最後のエンジン）。script は `~/.cache/losat-web-gui-target/s07p-resume/e2g_gates.sh`（`s07pp_gates2.sh` を E2g 用にしたもの）、記録は `gates-e2g-4.log`。変更前は S07++ の成果物（`s07pg-*`、native `a4ff5abb3a61…80f5`）。成果物のハッシュは `artifacts.sha256`（native `331fcba36447…44d8`、WASI の 4 つ、reactor の 2 つ）。

ゲートは 4 回走らせた：第 1 回（`run-20261002T002347Z`）と第 2 回（`run-20261002T021016Z`）は独立監査の指摘の修正のため、第 3 回（`run-20261002T031510Z`）は第 3 回の監査の D4 の残りの修正のために止め、部分の run は `~/.cache/losat-web-gui-target/e2g-aborted-runs/` に移した（記録には入れない）。

| 検査 | 結果 |
|---|---|
| `cargo fmt --check`（LOSAT、adapter）、clippy `-D warnings` 4 構成（all-targets all-features、no-default-features、wasm32-wasip1、wasm32-wasip1-threads）、adapter の clippy 3 構成 | すべて終了コード 0 |
| `cargo test --all-features` | 854 件通過、失敗 0。adapter の試験、wasm32 の web API の試験も通過 |
| Gate A と capture（236 件） | S02 の基準（`../losat_web_e1a/baseline/hashes.tsv`）と 236 件とも差 0（`capture-compare.txt`）。凍結ハッシュとの不一致は既知の `Sakai.MG1655.megablast` の 1 件だけ（許可リスト） |
| CI の速い検査の全件（`ci_fast_regressions.py --all-cases`） | 236 件、失敗 0、許可した既知の不一致 1 |
| BLASTN の fixture（`blastn_regression_fixtures.py check`） | 107 件、差 0 |
| outfmt 0 の fixture（1・2・4 スレッド、`check_losat.py`）と `precheck_hits.py` | 47 件、差 0（3 つのスレッド数とも） |
| `run_oracle.py`、`scoring_sweep.py`（S06） | 一致 |
| `scoring_sweep.py`（E2c、outfmt 0/6/7） | 各 300 件一致、580 件は同じ誤り |
| `word_size_sweep.py` | 256 件、差 0 |
| `slice_sweep.py` 17 種と `batch_sweep.py` 11 種（`BATCH_SIZE` 1000・100000、`CHUNK_SIZE` 40000、`CHUNK_SIZE` 2000 + `OVERLAP_CHUNK_SIZE` 50 を含む、変更前の実行ファイルの 2 種も） | 3710 件、差 0 |
| `split_check.py` | 220 件、差 0（NCBI の分割で出力が変わる 19 件も一致） |
| `check_inputs.py`（E2g の期待） | ゲートの run では E2g の期待のまま 6 件が「予期しない」になった：保守者の判断で `-evalue` の無限大・NaN（3 件）と句読点の題（3 件）の扱いが変わったため（`check-inputs.tsv`）。期待を直して同じ実行ファイルで走らせ直し、300 件とも期待どおり（`check-inputs-e2g.tsv`。NCBI が落ちる題の 2 件は例外 2） |
| `title_sweep.py` | ゲートの run は S07+ の期待（落ちる題は拒否）のままで 66 件が「予期しない」（`title-sweep.tsv`）。E2g の `title_sweep.py`（例外 2：NCBI が落ち、LOSAT の報告が代わりの題の NCBI の報告の置き換えと一致）で同じ実行ファイルを走らせ直し、1023 の定義行が一致 957・例外 2 が 66、予期しない 0（`title-sweep-e2g.tsv`）。script は以後これを使う |
| v1 の WASI の検査（`check_wasm_threading.py`） | 423 の記録、reactor の lifecycle の gate 通過、形式の失敗 0。reactor の記録の S05 との差 21 は、E2f と同じ LOSAT の文言の差だけ（`v1-reactor-records-compare.txt`）。`v1-requests` は E1d と一致 |
| V-ABI quick | 52 の実行が native の CLI と一致 |
| V-ABI full | 115 の検索 × 4（serial n1、threads n1・n2・n4）= 460 の実行、失敗した部分 0。凍結ハッシュ 676 件中 672 件一致、違う 4 件は既知の `Sakai.MG1655.megablast` outfmt 7 |

## V-PERF（指示書の 2.）

アプリ側の lock（`vperf_lock.sh`）を取り、静かな計算機で測った。変更前は S07++ の 3 つの成果物（`s07pg-native`、`s07pg-wasi-artifacts` の serial と threaded の command）、変更後は run の成果物。case は E2f と同じ 5 つ（`blastn`、`blastn-large`、`blastn-large-fmt7`、`blastn-large-fmt0`、`blastn-many`）× 3 つの実行の形。

| case | native | serial WASI | threaded WASI |
|---|---|---|---|
| blastn（`--repeat 3`） | ×1.013 | ×0.996 | ×2.601（FAIL、変更後の範囲 0.291〜1.124 秒） |
| blastn（`--repeat 5` で測り直し） | ×1.006 | ×1.048 | ×0.983 |
| blastn-large | ×1.000 | ×0.979 | ×0.855 |
| blastn-large-fmt7 | ×1.029 | ×1.005 | ×1.017 |
| blastn-large-fmt0 | ×0.994 | ×1.021 | ×0.991 |
| blastn-many | ×0.754 | ×1.023 | ×0.978 |

中央値の比（変更後 ÷ 変更前）。出力はすべて同じ。閾値を超えたのは `--repeat 3` の threaded の `blastn` だけで（変更後の最小 0.291 秒は変更前の中央値 0.301 秒より速い）、`--repeat 5` で測り直すと ×0.983（`perf-1.*`、`perf-2.*`、`perf-check-*.txt`）。非退行。

分割の時間（`split-timing.txt`、EDL933 × Sakai の 50 万文字、5 つの塊、`-task blastn`、1 スレッド、交互に 3 回）：NCBI 1.76・1.35・1.94 秒、LOSAT 1.36・1.46・1.45 秒、出力は同じ（1919 行）。

## 独立監査（指示書の 3.）

記録：[`audit/`](audit/)（共通と観点・回ごとの指示、第 2〜4 回の NOTES、監査 (c) の FINDINGS）。比較の全件と差の出力は作業ディレクトリ `~/.cache/losat-web-gui-target/e2g-audit/{a,b,c,r2,r3a,r3b,r4}/` にある（大きいので入れない）。

`ncbi_parity_auditor` の役割を、観点ごとに sonnet の agent 3 つが読み取り専用で並行して行った（共通の指示 `~/.cache/losat-web-gui-target/e2g-audit/COMMON.md`、観点ごとの指示 `ANGLE_{A,B,C}.md`。対象は `82e8c7593` の実行ファイル）。基準は `INVENTORY.tsv`。

### 第 1 回

| 観点 | 結論 | 内容 |
|---|---|---|
| (a) 経路の網羅 | unsupported（中程度 1 件） | NCBI を callgrind の下で 235 回実行して、届いた関数を記号の表から集め（main は strip されているので `blastn_app.cpp` は読んだ）、4410 の NCBI の定義を `INVENTORY.tsv` の行の範囲と突き合わせた。届いて行の無い関数 447 のうち、出力に効くものは 1 つだけで、ほかは記憶の管理、accessor、option の getter・setter、行のある関数の下請けなど。NCBI と LOSAT の比較 339 件（272 件が一致、67 件は明示的な拒否か承認済みの例外）、`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`・`BATCH_SIZE` の組の走査 228 件、すべての ASCII の文字を含む query 67 件、MG1655 × Sakai、EDL933 × MG1655（5 Mb を超える subject の分割）、10.9 Mb の複数の query で、LOSAT が受け付ける case はすべて一致。**中程度の差（F1）：** `OVERLAP_CHUNK_SIZE` が塊の大きさの約半分を超えると、NCBI は `CCoreException`（eNullPtr）で止まり（終了コード 3）、LOSAT は検索した。低い指摘：`blast_options.c` の検証関数、`GetDisplayIds`、`~CBlastFormat` の行が無い、古い参照 |
| (b) 移植の忠実さ | unsupported（重大度 低） | 項目ごとに NCBI のソースと読み比べ、3728 件を比べた（2790 件が一致、631 件が明示的な拒否（612 件は理由が成り立つ、19 件は広すぎる）、93 件が承認済みの例外、4 件は時間切れ）。T1〜T6・T8・T10〜T13・R1〜R3・V1 は試したすべてで忠実。指摘：D1（大きい重なりで NCBI が落ち LOSAT が検索。= F1）、D2（負の組の拒否が広すぎる）、D3（stdout と stderr を 1 つにしたときの警告の順）、D4（outfmt 0 の書き込みの失敗の前に警告を出す）、D5（閉じたパイプ。承認済みの例外 5） |
| (c) 拒否の理由 | unsupported（重大度 低、誤った出力の経路なし） | 拒否・例外の 86 行と、拒否の文言を出すコードを突き合わせ、手で作った FASTA 180、変異させた FASTA 6000、変異させた defline 30,800、52 の未対応のオプション × 3 つの書き方、60 を超える環境変数・得点・形式で、LOSAT が終了コード 0 で NCBI と違う出力を出す case は無かった。理由の一部が誤りか過大だった：HTML の文字参照は NCBI が落ちるのでなく decode する（K-crash でなく K-layer）、負の `CHUNK_SIZE` と `OVERLAP_CHUNK_SIZE` は塊が重なりより大きいときだけ NCBI が落ちる、`CHUNK_SIZE=1000` の NCBI の誤りは再現できる、R2 の一部の変数は通常の実行で効かない（族ごとの拒否）、penalty −32768 は NCBI が落ちるのでなく全 query が無効になる、greedy の gap の上限 32767 は控えめ（NCBI は約 600000 から失敗）、BOM の文言、`-num_threads` 65535 以上の失敗（全 program の CLI。S08+） |

### 第 1 回の指摘への対応

- `5dac71f72`：`BATCH_SIZE` なしの `CHUNK_SIZE=1000` は拒否をやめ、NCBI と同じく outfmt 0 の prolog の後に「BLAST engine error: Empty CBlastQueryVector」（終了コード 3）。負の組の拒否は、塊が重なりより大きいときだけに狭めた。R2 の文言を、拒否の実際の理由に直した。
- `993310891`：F1 の原因を NCBI のソースで特定した。`SplitQuery_CreateChunkData`（split_query_aux_priv.cpp:190-201）は塊ごとに新しい query の splitter で setup し、その塊がもう一度分かれないことを debug の build でだけ `_ASSERT` で確かめる。release の build で分かれると、`BlastSetupPreliminarySearchEx` は塊の lookup table を作らず（blast_aux_priv.cpp:206-208）、後で空の参照をたどる。LOSAT はこの条件（`chunk_would_be_split`：塊の長さ ÷（塊の大きさ − 重なり）が 2 以上）で明示的に拒否する。塊の大きさ 1500・3000・8000・20000 で、境界の両側とも NCBI と一致（NCBI が動く側は出力も一致、落ちる側は LOSAT が拒否）。
- `19a2f8397` ほか：棚卸しの表に 6 行を足し（1012 行）、分類と結果を直した。

### 第 2 回

第 2 回（`AC2` の実行ファイル、3418 件）：第 1 回の 3 つの修正（`5dac71f72`、`993310891`、`19a2f8397`）は、ソースとオラクルで確かめられた（塊をもう一度分ける条件は、size_t の回り込みを含む独自のモデルと 1152 件中 1152 件で一致）。結論は unsupported（低）：D4 と同じ。負の組の拒否の文言が、NCBI が動く場合を含むことの指摘。

### 保守者の判断と第 2 回までの指摘への対応

上の「保守者の判断」と「監査 (b) の指摘への対応」の表。

### 第 3 回（最後の実行ファイル `AC7` = `ad5fa9c85`）

| 部分 | 結論 | 内容 |
|---|---|---|
| A（監査 (b) の 3728 件の再実行と D3・D4 の拡大） | unsupported（低 1 件） | 3728 件：2973 件が一致、ほかはすべて承認済みの例外、理由の成り立つ明示的な拒否、保守者が残した拒否。第 1 回で一致し今違うもの（後退）は 0。第 1 回の D1・D3・D4・D2 の行は、例外 1 か一致に変わった。D3・D4 を広げた 1336 件（batch 1〜1000、title の警告、全 N、得点の表のない得点、`/dev/full`、`-out /dev/full`）は 1 つを除いて一致。**残り（D4 の残り）：** gapped の K-A の値のない得点で outfmt 0 の書き込みが失敗すると、LOSAT は保持した title の警告を出してから「BLAST failed to write output」（39 件）。また無効な query が先頭の batch では LOSAT は拒否（終了コード 1）、NCBI は prolog の書き込みで終了コード 6（40 件） |
| B（保守者の判断と第 2 回の T7 の 1823 件） | supported | 例外 1：343 件（NCBI の stack はすべて塊の null pointer）。65 件を NCBI の分けない最大の重なりと比べ 61 件がバイト一致、残りは塊の境の HSP か、query の多い batch の弱い HSP（NCBI 自身の重なりの間の差と同程度）。負の組は batch を分けるときだけ拒否、分けないときはバイト一致。例外 2：手で作った 54 組 432 実行、`,` `;` `~` 空白（`.`）の定義行をすべて試す 12392 実行で、NCBI が落ちる 272 の定義行は代わりの文字列の置き換えと一致、ほかはバイト一致。`-evalue` の 520 件。低い指摘：G1（塊に query の無い分割の拒否が PD に無い。理由は成り立つ）、G2（`-evalue 0x1p3`・`+0x10`・`1e` は NCBI が受け付け LOSAT が拒否。棚卸し A-207、既存） |

D4 の残りへの対応：`f4057718a`。NCBI の `PrintProlog`（blastn_app.cpp:258）は mixer と最初の batch の前で、書き込みの失敗はその最初の flush で起きる。LOSAT も outfmt 0 の prolog を batch の前に書いて flush する（file の sink は、報告を書くときに作るので従来どおり）。第 3 回 A の 1336 件を再実行して、残る差は承認済みの例外 3 と保守者が残した K-A の拒否だけ。fixture `write.devfull_*kaerror*`。G1 は PD に記録、G2 は既存の拒否として残す（下の残件）。

### 第 4 回（`AC8` = `f4057718a`）

supported。CLI 3958 件と web の経路 966 件（`run_local_blastn` を adapter と同じに呼び、outfmt 0/6/7 と診断を NCBI と比べ、observer の HSP の範囲 9156 も確認）。3280 件が一致、ほかは承認済みの例外 1（引数の構文）100、例外 2 が 12、例外 3 が 187、棚卸しに記録した明示的な拒否 379。未分類 0、後退 0。prolog を前に出したことで、Query is Empty!、`-max_target_seqs` の警告、`CHUNK_SIZE=1000`、K-A の誤り、title の警告、`-subject` の欠落・空、`-out` の file・`/dev/full` のどれも NCBI と一致。

## `main` への merge と毎晩の検査

この記録の commit の後に PR [#110](https://github.com/satoshikawato/LOSAT/pull/110) を S07++b と同じ流れ（merge commit、保守者は承認済み）で `main` に merge し、`nightly.yml` を `main` で 1 回 `workflow_dispatch` する。結果はこの節に追記する。

## 残件と引き継ぎ（S08 へ）

- **保守者の確認（2026-10-02、DW-16、推奨の案）：** DW-15 の後に分かった 3 つの細部。
  1. 負の `CHUNK_SIZE` が負の `OVERLAP_CHUNK_SIZE` より大きく、batch を分けるがもう一度は分けない場合（NCBI は隙間のある塊で検索して hit を落とす、決まった結果）：明示的な拒否のまま（`PD-LOSAT-NCBI-DEFECTS` 版 1.1）。
  2. reward 32768 以上（16 bit で 0 以下に回り、NCBI は全 query を無効にして hit 無し）：明示的な拒否のまま（同 版 1.1）。
  3. outfmt 6/7 の閉じたパイプで、LOSAT の書き込みが読み手が閉じる前に済むと LOSAT は終了コード 0、NCBI は後の flush で SIGPIPE（141）：承認済みの例外 3・5 の時間に依る部分（`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 版 1.2、AGENTS.md）。
- **`-evalue` の 16 進と `1e`：** NCBI は `0x1p3`・`+0x10`・`1e` を読む（`strtod`／`StringToDoublePosix`）。LOSAT は明示的に拒否のまま（棚卸し A-207、K-layer。アプリの出す書き方ではない）。BLASTP などの e-value の書き方は TD-13（S08+）。
- **整数でない環境変数の値：** NCBI の `CStringException` の文言はオラクルの build のソースのパスを含むので再現せず、明示的に拒否する（計画 TD-15、`check_inputs.py` の `e2g.t7.batch_size_text`）。
- **R2 の族ごとの拒否：** `DIAG_*` などは値によって NCBI の出力が変わらない場合も拒否する（監査 (b)：472 件中 125 の場面、第 2 回 74 件）。文言は族の理由を言う。誤った出力の経路は無い。
- **時間切れ（両方とも遅い）：** `-reward 32766 -penalty -32767` の megablast（NCBI と LOSAT とも 3300 秒を超える）、重なりが塊の大きさに近い 60 kb の query（1 万〜6 万の塊、第 3 回 B の 7 件）。
- **S08：** [`session_s08_e2b_tblastx_outfmt0_7.md`](../../losat_web_gui_sessions/session_s08_e2b_tblastx_outfmt0_7.md)（「S07+++b からの引き継ぎ」の節）。毎晩の WASI の TIMEOUT は runner の差だったので、S08 の最初の作業にはしない。
