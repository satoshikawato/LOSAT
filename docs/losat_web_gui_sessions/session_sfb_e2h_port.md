# Session SFb — E2h の続き：読み込み器の接続と移植、fixture、ゲート

## INSTRUCTION PROMPT

LOSAT の段階 E2h を SF から続ける。指示書の正本は [SF の指示書](session_sf_e2h_blastn_fasta_reader.md)（背景、NCBI の経路、手順 1〜9、「保守者の判断」1〜5、終了・引き継ぎ）で、この指示書は SF で終わったことと残りの順を足すだけである。完了条件は総合計画書 §7 の SF の行である。先に [セッション README](README.md) の共通規則を読む。エンジン側のセッションで、アプリ側（S13 の続き）と並行する（計画 DW-7）。

セッションは Linux の clone `/home/kawato/losat-work`（`$WORK_ROOT`）で起動する。そうすると `.claude/settings.local.json` の env（`WORK_ROOT`・`BUILD_ROOT`・`NCBI_SRC`・`NCBI_BIN`・`ORACLE_JOBS`）、hook、agent（`losat-implementer`・`losat-inventory`・`losat-mechanic`・`losat-reviewer`・`losat-test-runner`）が読み込まれる。SF は別の folder で起動したので、これらが無く、env を手で設定した（`/home/kawato/losat-baselines/sf-e2h-20261008/env.sh`）。場所は `CLAUDE.local.md`、手順は skill `losat-campaign`・`losat-worktree`・`losat-gates`・`losat-oracle-runs` に従う。

SF の結果（[ゲート記録](../evidence/losat_web_e2h/README.md)）：棚卸し [`INVENTORY.tsv`](../evidence/losat_web_e2h/INVENTORY.tsv)（222 行）、権威の記録 [`AUTHORITY.md`](../evidence/losat_web_e2h/AUTHORITY.md)、移植の計画 [`PORT_PLAN.md`](../evidence/losat_web_e2h/PORT_PLAN.md)（S0〜S10）、共有の読み込み器 `LOSAT/src/blastinput/fasta_reader/`（S0 まで。program にはまだ接続していない）、NCBI 2.17.0 を固定した fixture `LOSAT/tests/fasta_input_fixtures.py`（3932 行。変更前の LOSAT では same 1131、differs 13、拒否 2764、pending 24）。推奨の案で進めた判断はゲート記録の「推奨の案で進めた判断」1〜10 にある（ABI v1 は `from_bio` で変えない、Seq-id の判定、`register` の data loader など）。これに従う。

0. **準備**
   - 規則 2 の確認（worktree `$WORK_ROOT/.worktrees/web-gui`、`feature/losat-web-gui`、SF の最後のコミットから）。
   - タスクフォルダ `/home/kawato/losat-baselines/sfb-e2h-<yyyymmdd>/` を skill `losat-campaign` の手順で作り、セッションを bind する。SF の `STATE.md`・`DECISIONS.md`（`/home/kawato/losat-baselines/sf-e2h-20261008/`）を読み、`DECISIONS.md` は新しい folder に写して追記する。
   - アプリ側の S13 が終わり、`feature/losat-web-gui-app` が push されていれば `feature/losat-web-gui` へ merge する（README 規則 1。衝突の解消は別のコミット）。アプリ側のゲート記録に書かれた README の表と計画の変更も反映する。S13 が続いていれば、区切りのよいところでまた確かめる。
1. **移植 S1〜S10**（`PORT_PLAN.md` §2 の順。`losat-implementer` に 1 段階ずつ任せ、1 体ずつ）
   - 核酸の入力を先に終える：S1 共有の部品（`read_subjects`、`stream_is_empty`、`DATA_LOADERS` の実効値、`InputRecord`、`from_bio`）、S2 報告の byte の層（`report/defline.rs`：GenerateDefline、HtmlDecode、GuessEncoding）、S3 BLASTN の読み込みと報告、S4 BLASTN の中身の無いレコードと batch の時点の誤り、S5 TBLASTX、S6 TBLASTN の subject。
   - 次に蛋白：S7 BLASTP、S8 TBLASTN の query と後片付け（`blastn/input.rs` の拒否の関数と `fasta_input` の別名を消す。`BLASTINPUT_GEN_DELTA_SEQ` の確かめ直し）。
   - 最後にアダプタ：S9 `scan` の種類 1（核酸の旗）と 2（蛋白の旗）と性質試験、S10 `register`・`run`・`docs/web/abi_v2.md`。
   - 各段階の検証は `PORT_PLAN.md` §2 の Q と F に、変えた program の `fasta_input_fixtures.py check --rows <program> --jobs 3` を足す。その段階で移植した種類の行は `same` になり、ほかの行は変わらない（変わる行は理由を記録）。`PORT_PLAN.md` の「受け付ける入力で変わるのは `Subject_` で始まる subject の ID だけ」を確かめる。
   - 移植した箇所の直上に NCBI のファイル・行と断片を書く。BLASTX（`algorithm/blastx/`）は変えない。
   - 読み込み器の速度（S0 の後、bio の 1.2〜1.6 倍）は、接続の後の V-PERF で確かめる。
2. **fixture の仕上げ**（指示書の手順 5 の残り）
   - Seq-id の 24 行に LOSAT の明示的な拒否の文言（`not supported by LOSAT's <PROGRAM>` を含む）と終了コードを固定する。NCBI では流さない（ネットワークに出る）。
   - TBLASTX と BLASTP の `-lcase_masking` の 34 行：E2e で承認済みの例外か明示的な拒否として記録したかを `docs/evidence/losat_web_e2e/` で確かめる。記録が無ければ移植する（mask は読み込み器の小文字から作れる。推奨）。記録があれば、明示的な拒否として fixture の期待を直す。
   - `LOSAT/tests/ci_fast_regressions.py` に program ごとの `fasta_input_fixtures.py check` を足し、path の絞り込みに `LOSAT/tests/fasta_input_fixtures*` と `LOSAT/tests/fixtures/fasta_input/` を足す。全行が `same` か明示的な拒否になってから足す。
   - 変更前の `check` の結果（`$BUILD_ROOT/sf-e2h/fixtures/check-before.tsv`）と移植の後の結果をゲート記録に並べる。
3. **以降**は SF の指示書の手順 6〜9 をそのまま行う：sweep、ゲートの全工程、V-PERF、独立監査（4 観点とも supported まで）。
4. **ゲートの規模**（SF の指示書 7.）
   - 読み込み器は全 program が通るので、S11 で省いた工程（TBLASTX の option の sweep、capture、Gate A の 20 組、V-ABI full の全部）も行う。
   - 工程は並行させず、順に流す。各工程は `flock "$BUILD_ROOT/oracle.lock"` の下で、`--jobs` を明示する（V-ABI full と capture は 4、sweep は 3、TBLASTN・TBLASTX の option の sweep は 2）。
   - capture の変更前の基準は `~/.cache/losat-web-gui-target/sf/capture-before/hashes.tsv`、変更前の成果物は `~/.cache/losat-web-gui-target/sf/bin/` である。
   - 待つ間は watcher 1 本で、30 分以内の間隔で見る。V-PERF の lock は `~/.cache/losat-web-gui-target/s07p-resume/vperf_lock.sh` を使う（`s07p-resume/` の `verify_refs.py` と `vperf_lock.sh` は git に無いので消さない）。
   - SF の作業 folder `$BUILD_ROOT/sf-e2h/`（棚卸しの scratch、fixture の NCBI の出力 71 MB、S0 の速度の入力）は、ゲート記録が参照するので、この段階が `main` に入るまで消さない。
5. **NCBI の不具合**：`AUTHORITY.md` §K の 12 件は、各行の推奨の規則で移植する。保守者の答えがこのセッションの中で来たら従う。来なければ、最終回答で一括の問いとしてまた示す。
6. **区切り**
   - 核酸の入力（S1〜S6 と、その fixture の行が `same`）が終わった時点、またはコンテキストが 300k に近づいた時点で、残りを SFc に分けてよい（README の表に行を足す）。
   - 区切る前に、`STATE.md`・ゲート記録・次の指示書を更新して push する。
   - アプリ側の worktree（`.worktrees/web-gui-app`）とそのビルドには触れない。

最終回答は README 規則 8 のとおり（コミット SHA、push の結果、検証の結果、残件、次の指示書のパスと全文）。推奨案で進めた判断は `DECISIONS.md` とゲート記録に書き、最終回答に一覧で示す。終わったら、SF の指示書の「終了・引き継ぎ」（計画、`abi_v2.md`、`verification_cells.tsv`、README の表、S17 と SX の指示書、アプリ側への引き継ぎ）を行う。

## 進め方

- 移植は `losat-implementer` に 1 段階ずつ任せる。メインは、設計の判断と結果の確認だけをする。NCBI と LOSAT のコードを自分で大量に読まない。長い実行は `losat-test-runner` に任せ、ログはファイルに置く。
- 同じ worktree でコミットする agent は 1 体だけにする。ほかの agent は file を書くだけにし、コミットはメインが path を指定して行う（`git commit -- <path>`）。
- 検索を流す agent は同時に 2 体まで。始める前に `uptime` を見て、load average が 12 を超えていれば待つ（SF の間、アプリ側と別の project の Playwright の E2E で一時 66 になった）。
- 新しいビルドは `$BUILD_ROOT/<用途>` に置き、既存の名前（`native`、`test`、`wasi`、`adapter`、`gate-sf`）を使い回す。
