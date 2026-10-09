# Session SFd — E2h の続き：最後のコミットでのゲート、Gate A、V-PERF、修正の再監査、終了

## INSTRUCTION PROMPT

LOSAT の段階 E2h を SFc から続ける。指示書の正本は [SF の指示書](session_sf_e2h_blastn_fasta_reader.md)（背景、NCBI の経路、手順 1〜9、「保守者の判断」1〜5、終了・引き継ぎ）、[SFb の指示書](session_sfb_e2h_port.md)（ゲートの規模、進め方）と [SFc の指示書](session_sfc_e2h_protein_gates.md)（進め方）である。この指示書は、SFc で終わったことと残りの順を足すだけである。完了条件は総合計画書 §7 の SF の行である。先に [セッション README](README.md) の共通規則を読む。エンジン側のセッションで、アプリ側（S13b）と並行する（計画 DW-7）。

セッションは Linux の clone `/home/kawato/losat-work`（`$WORK_ROOT`）で起動する（`.claude/settings.local.json` の env、hook、agent が読み込まれる）。SF〜SFc は別の folder で起動したので、env を手で設定した（`/home/kawato/losat-baselines/sfc-e2h-20261009/env.sh`。agent には `agent_preamble.md` を渡した）。

SFc の結果（[ゲート記録](../evidence/losat_web_e2h/README.md) の「SFc で行ったこと」）：移植 S7〜S10 で 4 つの program の全入力とアダプタの `scan`（種類 1・2）・`register`・`run` が NCBI の `CFastaReader` の移植で読む。fixture 3932 行（監査の修正で 3993 行）は全行が NCBI と一致（same）か明示的な拒否（Seq-id 24 行、TBLASTX・BLASTP の `-lcase_masking` 34 行）。2 つの sweep の記録の無い拒否・差・時間切れは 0。CI の速い検査が fixture を確かめる。ゲートの全工程を監査の修正の前の `9dfe7efd` で試走し（run `20261009T145447Z`、記録は `/home/kawato/losat-baselines/sfc-e2h-20261009/gate-trial-record/`）、失敗は環境と script の不具合だけだった（直した）。独立監査 4 観点（`9dfe7efd`）の指摘 B-1・A-3・A-1 を直し、D-1 は文書、D-2 は v1 の変化が NCBI の byte であることを確かめて記録した（`321aa9ca`、`0be1afec`、`34dc8f99`、178 の再現が NCBI とバイト一致）。保守者の指示（2026-10-09「とりあえず main に適宜マージして」）で、検証の済んだ区切りごとに `main` へ PR を出して merge している（#116〜#118、SFc の最後の PR は下の最終回答を見る）。推奨の案で進めた判断はゲート記録の 16〜23 と `/home/kawato/losat-baselines/sfc-e2h-20261009/DECISIONS.md`。これに従う。

0. **準備**
   - 規則 2 の確認（worktree `$WORK_ROOT/.worktrees/web-gui`、`feature/losat-web-gui`、SFc の最後のコミットから）。
   - タスクフォルダ `/home/kawato/losat-baselines/sfd-e2h-<yyyymmdd>/` を skill `losat-campaign` の手順で作り、セッションを bind する。SFc の `STATE.md`・`DECISIONS.md`・`handoff-S7〜S10.md`・`handoff-fix.md` を読み、`DECISIONS.md` は新しい folder に写して追記する。
   - アプリ側（S13b）が終わり push されていれば `feature/losat-web-gui` へ merge する（README 規則 1。S13 は SFc で merge 済み、`39a563f1`）。
1. **修正の再監査**：SFc の監査の指示（`/home/kawato/losat-baselines/sfc-e2h-20261009/audit-brief/`、報告は `$BUILD_ROOT/sfc-e2h/audit/{A,B,C,D}.md`）を、監査するコミットを最後のコミットに変えて使う。(a)・(b)・(c)・(d) の各観点で、修正した点（B-1、A-3/C-1、A-1、D-1 の文書、D-2 の一覧）と、修正が触れた周り（`ncbi_environment.rs` の層の順、BLASTP の tabular の書き出し、gap の行の修飾子）を確かめ、4 観点とも supported まで続ける。検索を流す agent は同時に 2 体まで。
2. **ゲートの全工程**（SFb の指示書の 4. のまま、最後のコミットで、`FRESH=1`）：`docs/evidence/losat_web_e2h/gates/sf_gates.sh`。試走から直したこと：TBLASTN の Stage G は `/mnt` の tmpfs の中の bind（既定）、`oracle-check` は `TMPDIR=/tmp`、`capture` は既知の不一致 `Sakai.MG1655.megablast` を許可、`collect` の下位の directory。試走の時間：build 7 分、fast-all 23 分、option の sweep 40 分、capture 26 分、V-ABI full 3 時間 21 分（TBLASTX が 3 時間 17 分）、全体で約 6 時間。各工程は `flock "$BUILD_ROOT/oracle.lock"` の下で順に。
   - check_inputs の 1 行（outfmt 0 の `Database:` の行が worktree の path を含む）は `$WT` から流すと一致する（`handoff-fix.md`）。
3. **Gate A**（`sf_gate_a.sh`、20 組、単独、3〜5.6 時間）。`/tmp/losat-pr5-runtime-cert-*` を置き直す（WSL の再起動で消える）。
4. **V-PERF**（`sf_perf.sh`、変更前 `~/.cache/losat-web-gui-target/sf/bin/` と交互）：読み込みの重い case（10 万レコードの query、5 Mb の 1 行と 80 桁の subject、BLASTP の多数の query）の大きさは、ネイティブで 1 回測ってから決める。`-`（標準入力、1 byte ずつの補充）の case も足す。lock は `~/.cache/losat-web-gui-target/s07p-resume/vperf_lock.sh`（アプリ側の試験を止める）。閾値を超えたら静かな計算機で `--repeat 5`。その後 `FORCE=1 STAGES=collect sf_gates_resume.sh`。
5. **終了**：完了条件（計画 §7 の SF の行）を満たしたら、ゲート記録を「完了」にし、SF の指示書の「終了・引き継ぎ」（計画 §0.5 の TD-8・TD-12、§5.4、§7、§10 の FASTA の読み方の行、状態の行、`docs/web/verification_cells.tsv` の SF の升目、README の表、S17 の指示書の 1. と 6.、SX の指示書、アプリ側への引き継ぎ）を行う。アプリ側への引き継ぎにはゲート記録の「アプリ側への引き継ぎ（SFc で増えたもの）」を含める。`main` への PR を出し、CI が通れば merge する（保守者の指示 2026-10-09）。
6. **NCBI の不具合と保守者への問い**：`AUTHORITY.md` §K の 12 件、BLASTP の句読点だけの題（承認済みの例外 2 を広げるか）、ABI v1 の出力の変化（判断 12 と D-2、`/home/kawato/losat-baselines/sfc-e2h-20261009/logs/fix/v1-changes.md`）、S13 の判断 1〜3（アプリ側のゲート記録）。保守者の答えがこのセッションの中で来たら従う。来なければ、最終回答で一括の問いとしてまた示す。

区切り：コンテキストが 300k に近づいたら残りを次のセッションに分けてよい（README の表に行を足す）。区切る前に `STATE.md`・ゲート記録・次の指示書を更新して push する。アプリ側の worktree（`.worktrees/web-gui-app`）とそのビルドには触れない。最終回答は README 規則 8 のとおり。推奨案で進めた判断は `DECISIONS.md` とゲート記録に書き、最終回答に一覧で示す。

## 進め方（SF〜SFc で分かったこと）

- 長い実行は background で始めて 1 回だけ待つ（短い間隔で見に行かない）。ゲートの最中は `$WT` にコミットしない（`guard` が HEAD を確かめる）。どうしても branch を進めるときは一時的な worktree（`$WORK_ROOT/.worktrees/` の下）でコミットし、`git push origin HEAD:feature/losat-web-gui` で push する（SFc で 2 回、ゲートの worktree は動かさなかった）。
- PR の CI：`watch_prs.py` は check が 1 つ動いている間に green を報告したことがある。merge の前に必ず `gh pr view` で全 check の完了を確かめる（SFc の `scripts/wait_checks.py`）。`browser-engine`・`app` の firefox と webkit の E2E で 5 秒の待ちの timeout が繰り返し出た（毎回別の試験、`gh run rerun --failed` で通過）。
- TMPDIR を `$TASK_DIR/tmp` にすると、path を出力に含む凍結（`-db_gencode` の代わりの DB の `# Database:`）が変わる。その工程だけ `TMPDIR=/tmp` で流す。
- V-ABI の実行ファイルの名前は `LOSAT` でなければならない（clap の usage の行）。
- 検索を流す agent は同時に 2 体まで。始める前に `uptime` を見て、load average が 12 を超えていれば待つ。
- 新しいビルドは `$BUILD_ROOT/<用途>` に置き、既存の名前（`native`、`test`、`wasi`、`adapter`、`gate-sf`）を使い回す。`$BUILD_ROOT/sf-e2h/`・`sfb-e2h/`・`sfc-e2h/` は、ゲート記録が参照するので、この段階が `main` に入るまで消さない。
