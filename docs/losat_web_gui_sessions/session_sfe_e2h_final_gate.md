# Session SFe — E2h の続き：修正の再監査の最後の回、最後のコミットでのゲート、Gate A、V-PERF、終了

## INSTRUCTION PROMPT

LOSAT の段階 E2h を SFd から続ける。指示書の正本は [SF の指示書](session_sf_e2h_blastn_fasta_reader.md)（背景、NCBI の経路、手順、「保守者の判断」1〜5、終了・引き継ぎ）、[SFb の指示書](session_sfb_e2h_port.md)（ゲートの規模、進め方）、[SFc の指示書](session_sfc_e2h_protein_gates.md) と [SFd の指示書](session_sfd_e2h_final_gate.md)（進め方）である。この指示書は、SFd で終わったことと残りの順を足すだけである。完了条件は総合計画書 §7 の SF の行である。先に [セッション README](README.md) の共通規則を読む。エンジン側のセッションで、アプリ側（S14）と並行できる（計画 DW-7）。

セッションは Linux の clone `/home/kawato/losat-work`（`$WORK_ROOT`）で起動する（`.claude/settings.local.json` の env、hook、agent が読み込まれる）。SFd は `/home/kawato/losat-baselines` で起動したので、env を手で設定した（`/home/kawato/losat-baselines/sfd-e2h-20261010/env.sh`。agent には `agent_preamble.md` を渡した）。

SFd の結果（[ゲート記録](../evidence/losat_web_e2h/README.md) の「SFd で行ったこと」と判断 25〜29）：

- アプリ側の S13b を merge した（`dcbb149d`。README・計画 DW-25・DW-26 は `2f954c66`）。
- 最後のコミット `2f954c66` でゲートの全工程を始め（run `20261009T232829Z`）、修正の再監査の指摘でエンジンを変えることになったので `fast-all` で止めた（それまでの工程は全部終了コード 0。最後のコミットの run ではない）。
- 修正の再監査の第 1 回（`2f954c66`）：C（蛋白の入力、BLASTP の報告の順、拒否の理由）は supported（3,787 件）。B は FASTA の読み方の 666 件が全部一致し、A・B の指摘は `[BLAST] DATA_LOADERS` を決める registry の層だけ（R1〜R4）。D は ABI v1 の定義行の変化（D-1）。
- registry の層（`LOSAT/src/blastinput/ncbi_environment.rs`）を 4 回直した：`f5729e60`（R1〜R4）、`68cc3ba5`（第 2 回の監査 R の R-1〜R-9、R-8 は `AUTHORITY.md` §K-13）、`201651b8`（第 3 回の監査 S の S-1〜S-5、S-4 は §K-14）。再現は task folder の `scripts/registry_repros.py`（398 件）、`registry_repros3.py`（477 件）、`s_harness/`（193 件）で、どれも差 0 か記録した拒否。監査の回ごとに指摘は特殊な起動条件（重複した環境変数、`exec -a`、長い passwd の行）に移った。
- ABI v1 の層（`LOSAT/src/web_api/v1_bio.rs`）を 3 回直した：`844de3d8`（D-1）、`a5d3421b`（`>?` の行の後の番号、outfmt 0 の最後の byte の脱落）、`2e82c250`（第 2 回の監査 V の V-1：`bio` が読まない record の定義行を調べない。V-2：単独の CR の後で NCBI の record の数は行ごとには数えられないので、その後の変わった record を出るところで拒否する、範囲は判断 12 の (3)。V-3：役割の ID か題を出す形式だけ。hit の無い record は v1 の層から見えないので拒否のまま、判断 30）。決着しない 2 種類は保守者への問い（判断 28 の (4)）。
- V-PERF の読み込みの重い case の大きさをネイティブの計測で決め、標準入力の case を足した（`cf809a87`：50 Mb の subject、2000 の蛋白、`blastn-q100k-stdin-file`・`-stdin-pipe` はネイティブだけ）。
- 推奨の案で進めた判断は `/home/kawato/losat-baselines/sfd-e2h-20261010/DECISIONS.md` の「SFd」とゲート記録の 25〜29。これに従う。
- `main` への PR は SFd では出していない（最後の PR は #119、`dd6fadb8`）。

0. **準備**
   - 規則 2 の確認（worktree `$WORK_ROOT/.worktrees/web-gui`、`feature/losat-web-gui`、SFd の最後のコミットから）。
   - タスクフォルダ `/home/kawato/losat-baselines/sfe-e2h-<yyyymmdd>/` を skill `losat-campaign` の手順で作り、セッションを bind する。SFd の `STATE.md`・`DECISIONS.md`・`handoff-fix.md`・`handoff-fix2.md`・`handoff-fix3.md`・`handoff-v1fix.md` を読み、`DECISIONS.md` は新しい folder に写して追記する。
   - SFd の監査と修正の build（`$BUILD_ROOT/sfd-e2h/`：`audit*`、`fix*`、`v1fix`）はゲート記録が参照するので、この段階が `main` に入るまで消さない。
1. **修正の再監査の最後の回**（最後のコミットで、2 観点を並行。検索を流す agent は同時に 2 体まで）
   - (R) registry の層：`/home/kawato/losat-baselines/sfd-e2h-20261010/audit-brief/re3_registry_brief.md` を、監査するコミットを最後のコミットに、前の回を第 3 回の監査 S（`$BUILD_ROOT/sfd-e2h/audit3/S.md`）と第 4 回の修正（`handoff-fix3.md`、`scripts/s_harness/`）に変えて使う。
   - (V) ABI v1：`.../audit-brief/re2_v1_brief.md` を、前の回を第 2 回の監査 V（`$BUILD_ROOT/sfd-e2h/audit-v1/V.md`）と v1 の第 3 回の修正に変えて使う。
   - 監査の対象は source の写し（`git archive`）、ゲートと同じ build のバイナリと reactor に固定する。
   - 指摘があれば直して、その観点だけ再監査する。特殊な起動条件の指摘は、SFd と同じく、はっきりした安い規則は移植し、残りは範囲を記録した明示的な拒否（§J）か NCBI の不具合（§K、保守者への問い）にして閉じる。
   - B（核酸の読み方）・C（蛋白の読み方と BLASTP の報告）・D の v2 の部分は、第 1 回の後に変わったのが registry と v1 の層だけなので、再監査しない（理由をゲート記録に書く）。
2. **ゲートの全工程**（最後のコミットで、`FRESH=1`）：`docs/evidence/losat_web_e2h/gates/sf_gates.sh`。SFd の指示書の 2. のとおり（約 6 時間、background で始めて 1 回だけ待つ）。ゲートの最中は `$WT` にコミットしない。
3. **Gate A**（`sf_gate_a.sh`、20 組、単独、3〜5.6 時間）。
4. **V-PERF**（`sf_perf.sh`）：大きさは決めてある（`cf809a87`）。入力は script が作り直す（`$BUILD_ROOT/sfb-e2h/perf-inputs/` は SFd で消した）。静かな計算機で（script は load average 1 未満を 20 分まで待つ。他のセッションの重い処理が終わってから始める）。SFd の大きさを決める 1 回の計測（負荷の高い中、`/home/kawato/losat-baselines/sfd-e2h-20261010/logs/perf-sizing.log`）では 5 Mb の subject（0.28→0.55 秒）と標準入力の pipe（1.76→2.82 秒）に遅くなった兆しがあった（交互の 5 回では再現しなかった）。この 2 つは結果をよく見る。閾値を超えたら `--repeat 5`。その後 `FORCE=1 STAGES=collect sf_gates_resume.sh`。
5. **終了**：完了条件（計画 §7 の SF の行）を満たしたら、ゲート記録を「完了」にし、SF の指示書の「終了・引き継ぎ」（計画 §0.5 の TD-8・TD-12、§5.4、§7、§10 の FASTA の読み方の行、状態の行、`docs/web/verification_cells.tsv` の SF の升目、README の表、S17 の指示書の 1. と 6.、SX の指示書、アプリ側への引き継ぎ）を行う。アプリ側への引き継ぎは S14 の指示書に節として足す（ゲート記録の「アプリ側への引き継ぎ（SFc で増えたもの）」を含める）。`main` への PR を出し、CI が通れば merge する（保守者の指示 2026-10-09）。
6. **NCBI の不具合と保守者への問い**：`AUTHORITY.md` §K の 14 件（SFd で K-13 の parameter の変数の異常終了、K-14 の `=` の無い環境の項目の SIGSEGV を足した）、BLASTP の句読点だけの題（承認済みの例外 2 を広げるか）、ABI v1 の出力の変化と v1 の新しい拒否（判断 12・26・28、とくに 28 の (4)）、S13 の判断 1〜3、S13b の判断 1〜88（DW-25）。保守者の答えがこのセッションの中で来たら従う。来なければ、最終回答で一括の問いとしてまた示す。

区切り：コンテキストが 300k に近づいたら残りを次のセッションに分けてよい（README の表に行を足す）。区切る前に `STATE.md`・ゲート記録・次の指示書を更新して push する。アプリ側の worktree（`.worktrees/web-gui-app`）とそのビルドには触れない。最終回答は README 規則 8 のとおり。推奨案で進めた判断は `DECISIONS.md` とゲート記録に書き、最終回答に一覧で示す。

## 進め方（SF〜SFd で分かったこと）

- SFd の指示書の「進め方」はそのまま有効（長い実行は 1 回だけ待つ、ゲートの最中は `$WT` にコミットしない、`wait_checks.py` で全 check の完了を確かめる、TMPDIR、V-ABI の実行ファイルの名前は `LOSAT`、load average 12、`$BUILD_ROOT` の名前）。
- ゲートは最後のコミットでしか意味が無い。エンジンを変える見込みのある監査は、ゲートの前か並行で行う。止めるときはゲートのプロセス群に `kill -TERM -- -<PGID>`（`ps -eo pid,pgid,args` で `sf_gates.sh` の PGID を見る。`pkill -f` は使わない）。
- 監査の回ごとに、registry の層の指摘は細かい起動条件へ移った。NCBI の規則が source ではっきりしていて安いものだけ移植し、残りは範囲を記録した拒否か §K にして閉じる。実装役には、直した後に条件の行列を自分で NCBI と突き合わせさせると往復が減った。
- 2 体の実装役を並行させるときは、片方を別の worktree（`.worktrees/<topic>`）と別の `$BUILD_ROOT/<topic>` で動かし、push の前に `origin/feature/losat-web-gui` へ rebase させる（SFd で 2 回、衝突は文書だけ）。
- `execve` で `argv[0]` と envp を正確に渡す監査の harness は `/home/kawato/losat-baselines/sfd-e2h-20261010/scripts/s_harness/`。
- 同じ計算機の他のセッション（アプリ側 S14 など）が task folder を間違えて書き換えたことがある（SFd、2026-10-10）。`DECISIONS.md` と `STATE.md` は変更の後に読み直す。
