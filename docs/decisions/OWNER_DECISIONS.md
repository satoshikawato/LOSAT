# Owner decisions cited by the repository

The repository cites the Owner's decisions in two ways: by date ("Owner, 2026-10-17") and by an entry
of a work line's decision log ("drop-in survey `DECISIONS.md` 12", "v4-C DECISIONS 25"). The logs are
kept outside the repository. This file copies each cited decision: its date, the Owner's words as
written (in Japanese) with a summary in English, and where the repository cites it. A few entries
that the repository does not cite are kept because a cited entry depends on them ("related entry").

As of dev 8a5656f7f (2026-11-05).

- Kind: Owner (the Owner's own words, or the Owner's answer to a question), Owner, in part (an entry
  that records both), Owner-delegated (a working session decided it under the Owner's standing
  delegation of 2026-09-29, restated 2026-10-07; the entry has no words of the Owner), Record (an
  entry that records a fact, not a decision), Not found (the repository cites the Owner and a date,
  but no words of the Owner for it were found in the work records; the row gives what they say).
- An Owner-delegated decision is not the Owner's approval where a rule asks for one (Owner,
  2026-11-05); [`SWITCH_APPROVALS.md`](SWITCH_APPROVALS.md) applies this to the stop rule for
  conditional switches.
- Line and stage IDs are in [`LINE_IDS.md`](LINE_IDS.md).
- When the repository cites a new decision, add its row here in the same pull request.

## Decisions cited by date

| Date | Rule | Kind | Owner's words | Decision | Scope | Cited in | Source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2026-11-05 | Web builds outfmt 6/7 column specifications (WA2-OP-03) | Owner | 「WA2-OP-03: 列指定を Web で扱う」 | The Web also builds outfmt 6 and 7 with a field list; an outfmt 6 or 7 in Additional options with fields sets the output format. | LOSAT Web argv and output format | 3: web/adapter/src/run.rs:1532, web/app/tests/e2e/results.spec.ts:2049 | web-audit-20261104/w5/DECISIONS.md:19 |
| 2026-11-05 | An Owner-delegated decision is not the Owner's approval | Owner | (a) | The Owner chose option (a), 「Owner-delegated では Owner の承認に代えられない」: a decision that a working session made under the standing delegation does not replace the Owner's approval that the stop rule for conditional switches asks for. Approvals are kept in the Owner's words in `docs/decisions/`; a switch without one is approved by the Owner one by one or becomes a candidate for structural line SF, stage S3. | approvals that a rule asks of the Owner | 1: docs/decisions/SWITCH_APPROVALS.md:18 | swebok-improve-20261105/OWNER_QUESTIONS.md:49, swebok-improve-20261105/OWNER_QUESTIONS.md:53 |
| 2026-11-04 | No "verification/certified" wording; the mark is Supported or Not yet supported | Owner | 「verificationっていう表現なくそうよ。前も言ったけど。Readmeも含め。Supported かNot yet supportedだけでしょ。」 | User-facing text, README and CLI help do not say verification, verified or certified; the run mark is "Supported" or "Not yet supported". | Web UI, README, CLI help | 2: changelog.d/fix-supported-wording-cli-pd.md:1, changelog.d/fix-web-supported-wording.md:2 | web-audit-20261104/prompts/W1.md:7, web-sx-blastx-20261102/DECISIONS.md:38,43 |
| 2026-11-04 | A reversed range offers a swap | Owner | 「範囲の逆順は『ひっくり返しますか？』って聞いてあげれば？」 | When a start is after the stop, the form offers the swap beside the problem; outside the record it offers none. | LOSAT Web region picker | 2: web/app/src/ui/RegionPicker.vue:27, web/app/tests/e2e/search.spec.ts:325 | web-audit-20261104/w3/STATE.md:74 |
| 2026-11-04 | The Web app does not read .bz2 files | Owner | 「bz2やめとく？」「bz2ウェブはやめよう」 | The browser reads .gz and .zst inputs but not .bz2 because a pure-JS decoder is slow and memory-hungry for a rare format. | LOSAT Web compressed inputs | 1: docs/web/ncbi_ui_mapping.md:267 | gzinput-20261103/DECISIONS.md:46 |
| 2026-11-03 | Keep searching while you're away (wake lock) is on by default | Owner | 「OK. 規定ONにしとこうぜ」 | The wake lock is on by default; the checkbox shows only on handheld devices under the label "Keep searching while you're away". | LOSAT Web attention / wake lock | 2: web/app/src/application/attention.ts:5, web/app/tests/e2e/search.spec.ts:1018 | web-refresh-20261102/DECISIONS.md:208 |
| 2026-11-03 | Explicit fetch may send an identifier and file URL to NCBI, EBI and CNCB | Owner | 「研究データは送らない。明示に頼んだ取得（accession、Example）のときだけ、その識別子と公開ファイルの URL を取得先へ送る」 | Research data is never sent; only on an explicit fetch (accession or Example) the identifier and public file URL go to NCBI, EBI or CNCB. | LOSAT Web privacy promise, fetch | 2: docs/web/losat_web_design_v0.1.md:477, docs/web/ncbi_ui_mapping.md:32 | web-fetch-20261103/DECISIONS.md:3-5 |
| 2026-11-02 | Whole Web UI redesigned in NCBI's three layers | Owner | 「配列を入れる画面もデザインを更新して。てかWebUI全体」「box、背景白いとちょっと明るすぎる」「背景とかも含めてデザインしなおそう」 | The whole UI uses a soft blue ground, grey groups with visible edges and accent name tabs, white inputs, tables and outputs, and no card shadows. | LOSAT Web styling | 1: changelog.d/feature-web-params-readable.md:8 | web-refresh-20261102/DECISIONS.md:54 |
| 2026-11-02 | Runs that read no time do not wait in the timed queue | Owner | 「予約だけでほぼ使っていない枠」 | Runs that read no time (the Wasm trap-counting modes) leave the timed queue and run on another node under a shared marker. | benchmark scheduling | 1: docs/evidence/wasm_oob_20261101/scripts/wo_ab2_q2.sh:3 | wasm-oob-20261101/STATE.md:37, wasm-oob-20261101/DECISIONS.md:17 |
| 2026-11-02 | Web copy follows stop-describing-ui; NCBI sentences stay | Owner | 「OK基本これでいこう。NCBI 由来の説明文は残す」 | LOSAT Web strips its own explanatory copy and shows what can be done, but keeps sentences taken from NCBI. | LOSAT Web copy | 2: changelog.d/feature-web-copy-stop-describing.md:1, web/app/src/ui/SearchPanel.vue:9 | the Owner's chat with a working session (no task-folder record) |
| 2026-11-02 | Run LOSAT opens the run's results directly, like NCBI BLAST | Owner | 「あと検索開始したら結果に直接飛ぶようにして。NCBIBLASTみたいに。」 | After Run LOSAT the page goes straight to the results page of the run it started and shows its status until it completes. | LOSAT Web search flow | 1: web/app/src/ui/AppView.vue:32 | web-refresh-20261102/STATE.md:3-8 |
| 2026-11-02 | Opened FASTA file shows name, size and count; record list starts closed | Owner | 「このファイル読み込むとレコードがずらっと出てくるのいらない。NCBIは出ないじゃん」 | An opened FASTA file shows its name, size and record count; the record list and first lines stay closed until opened, as NCBI shows only the file. | LOSAT Web file card | 2: changelog.d/feature-web-mode-term.md:8, web/app/src/ui/SourceCard.vue:150 | web-refresh-20261102/DECISIONS.md:86 |
| 2026-11-02 | Algorithm parameters form reads like NCBI's | Owner | 「NCBIと比べるとやっぱ見づらい。なんで？次ノセッションで対応したい。プロンプト作って」 | Every parameter control shows the value the search uses, NCBI's selects and Mask checkboxes are used, and labels, groups and rows follow NCBI's readability. | LOSAT Web Algorithm parameters | 1: changelog.d/feature-web-params-readable.md:1 | web-refresh-20261102/STATE.md:34 |
| 2026-11-02 | Okabe-Ito colour-universal palette for accents and score colours | Owner | 「案 4 でいきましょう」 | The accent, status colours and the Graphic Summary's five score colours use the Okabe-Ito colour-universal palette; the Dot Plot keeps blast2dotplot.py's blue and orange. | LOSAT Web colours | 2: changelog.d/feature-web-gbdraw-look.md:5, web/app/src/domain/plot-scale.ts:144 | web-refresh-20261102/STATE.md:32 |
| 2026-11-02 | Say "mode" for BLASTN etc. and "Task Selection" | Owner | 「Modeにしてください」「Program Selection … これはTask selectionにすべき」 | LOSAT Web calls BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX modes and names the block that chooses the task "Task Selection"; display names and output files do not change. | LOSAT Web wording | 1: changelog.d/feature-web-mode-term.md:1 | web-refresh-20261102/DECISIONS.md:66 |
| 2026-11-02 | Look B of seven: dark blue header band, grey-blue mode tabs | Owner | 「灰色で見ずらい。ていうか全体的に見ずらいから…モック作って提案して」「あと4つ案を出して」「Bにします」 | LOSAT Web uses a dark blue header band, grey-blue mode tabs on an accent rule, white groups and darker secondary text; usable text meets 4.5:1. | LOSAT Web styling | 1: changelog.d/feature-web-mode-term.md:4 | web-refresh-20261102/DECISIONS.md:77 |
| 2026-11-02 | Web look follows gbdraw (Inter, rounded cards, blue accent) | Owner | 「色とかボタンを gbdraw の GUI っぽくしてくれないか」「フォントも Inter で」 | LOSAT Web takes colours, buttons and the Inter font from gbdraw's Web GUI while keeping NCBI's layout, order and words. | LOSAT Web styling | 1: changelog.d/feature-web-gbdraw-look.md:1 | the Owner's chat with a working session (no task-folder record) |
| 2026-11-02 | Dot Plot chooses its subject in its own tool band | Owner | 「Dotplotの対象配列選択の挙動が謎なんだけどこれわかりづらくない？」 | The Dot Plot has a Subject list in its tool band, in Descriptions order; a Descriptions row and the list are one selection. | LOSAT Web Dot Plot | 3: changelog.d/feature-web-dotplot-pick.md:2, web/app/src/ui/DotPlot.vue:12 | web-refresh-20261102/DECISIONS.md:92 |
| 2026-11-02 | Start BLASTX in LOSAT Web without waiting for the DW-10 certification | Owner | 「DW-10 の認証を待たずに始めてよい」 | BLASTX work in LOSAT Web starts without waiting for the v0.2.0 certification named as the entry condition of DW-10. | LOSAT Web, BLASTX | 1: docs/losat_web_gui_plan.md:88 | web-sx-blastx-20261102/DECISIONS.md:3 |
| 2026-11-01 | `LOSAT_X_BPDBWIDE` is not a switch that runs under a rule | Owner | 「規則で動く switch ではない」 | `LOSAT_X_BPDBWIDE` (#670) does not run under a rule, so the stop rule's deadline fallback is not needed for it. | blastp, `LOSAT_X_BPDBWIDE` (#670) | 1: docs/decisions/SWITCH_APPROVALS.md:42 | dropin-survey-20261029/DECISIONS.md:204 |
| 2026-10-31 | Merge #652 (`LOSAT_X_TBXFEWGROUPS`) too | Owner | 「これも」 | Continuing 「直してマージして」 (fix and merge), said of another pull request, the Owner asked for #652 to be fixed and merged as well. #652's body named `LOSAT_X_TBXFEWGROUPS` a deadline fallback and asked for the Owner's OK; the merging session read this answer as that OK. No deadline date was given. | tblastx, `LOSAT_X_TBXFEWGROUPS` (#652) | 1: docs/decisions/SWITCH_APPROVALS.md:31 | dropin-survey-20261029/TBX/STATE.md:70 |
| 2026-10-28 | blastn query-batch slots removed; time gain must cover the memory rise | Owner | 「1割速くするために1割メモリ増やすのは100歩譲って受け入れられても、1割速くするために2割メモリ増やすのはダメ。最悪の場合でも1割-1割ルール」 | The blastn query-batch slots (QBATCHPAR) and their tuning switches are removed because max RSS rose more than the time saved; the per-case rule is RSS rise percent at most time saved percent. | blastn query-batch slots, speed versus memory rule | 1: changelog.d/revert-qbatchpar-remove.md:4 | v4-W/DECISIONS.md:738 |
| 2026-10-27 | Speed-up PRs name the structural cause; no new rule-limited switch | Owner | 「提案に賛成」 | A speed-up PR names the structural cause it removes and the inputs where its form fails; no new rule-limited switch while three structural items are open; thresholds are folded before graduation. | speed-up pull requests, switches | 3: AGENTS.md:237, AGENTS.md:244 | speedup-plan-v4-20261017/DECISIONS.md:382, pr-audit-20261027/REPORT.md:139 |
| 2026-10-25 | Threaded Wasm n16 out-of-bounds stops the v0.2.0 GO until the cause is known | Owner | 「原因が出るまで GO を止める」 | The release GO waits for a cause of the threaded-Wasm memory access out of bounds at thread start; the Owner chose this over shipping threaded Wasm as experimental. | v0.2.0 release, threaded Wasm | 1: docs/evidence/wasm_oob_20261025/README.md:3 | release-v020-20261025/OWNER_QUESTIONS.md:133 |
| 2026-10-25 | macOS libm may differ in the last digit of a full-precision E value (PD-LOSAT-PLATFORM-LIBM) | Owner | (a) | For v0.2.0, the last digit of a full-precision E value (outfmt 11 and 8) may differ on macOS because of Apple libm; a bit-identical exp/log port is a later line. | macOS, outfmt 11 and 8 | 2: changelog.d/fix-macos-libm-last-digit.md:1, docs/release/v0.2.0_certification_plan.md:162 | release-v020-20261025/DECISIONS.md:42 |
| 2026-10-23 | LOSAT_X_DPNEON2 joins the default set | Owner-delegated | 「どっちがいいと思いますか？」 | The Owner asked which option is better for the remaining opt-in speed-up; the session adopted its recommendation, that DPNEON2 joins the default set. | aarch64 NEON X-drop kernels | 1: LOSAT/src/utils/xswitch.rs:110 | switch-tidy-20261022/DECISIONS.md:173-184 |
| 2026-10-22 | Wide integer instead of reproducing NCBI undefined behaviour (culling limit sum) | Owner | 「広い型になおしましょう」 | LOSAT adds the culling limit in i64 and does not reproduce NCBI's Int4 overflow, so every limit above 1 is three more. | HSP culling limit | 1: LOSAT/src/core/hspfilter_culling.rs:160 | ncbi-parity-20261019/DECISIONS.md:115 |
| 2026-10-22 | Unread evidence goes to GitHub release assets | Owner | 「推奨通り」 | Evidence that no code, CI, script or skill reads moves to a GitHub release asset; the repository keeps a manifest and the asset hash. | docs/evidence housekeeping | 1: AGENTS.md:475 | slim-20261022/DECISIONS.md:9-10 |
| 2026-10-17 | Timed runs need only their own idle CPUs | Owner | 「けっこう空いてる気がするけど」「ほかの子にも伝えて」 | A timed benchmark run needs only its own physical cores and their hyper-thread siblings to be idle; other cores of the shared benchmark machine may be busy. | benchmark scheduling | 3: benchmarks/db6-small/scripts/cpu_busy.py:2, benchmarks/db6-small/scripts/ws_db6_bench.sh:6 | perf-mb-20261014/STATE.md:796, dbsearch-20261012/DECISIONS.md:855 |
| 2026-10-17 | Speed-ups land default ON in their own PR | Owner | 「てか速くする改善、実装する時opt-inじゃなくてデフォルトでONにしてくれない？毎回二重にCI通すじゃん。みんなに伝えて。」 | An exact speed-up joins the strict set in the pull request that adds it, with LOSAT_X_<NAME>=0 as the off path; no opt-in PR followed by a default-set PR. | speed-up switches (LOSAT_X_*) | 3: .claude/skills/losat-ship-pr/SKILL.md:40, AGENTS.md:211 | the Owner's chat with a working session (no task-folder record) |
| 2026-10-17 | A speed-up that cannot be default ON is still adopted | Owner | 「改善策、基本デフォルトでONにしてとは言ったけどデフォルトでONにできなければ採用しないとは言ってないからね」 | Default ON is the normal form, not a condition of acceptance: a speed-up that loses on a measured condition lands opt-in with its win and loss conditions recorded. | speed-up switches (LOSAT_X_*) | 1: .claude/skills/losat-ship-pr/SKILL.md:51 | the Owner's chat with a working session (no task-folder record) |
| 2026-10-17 | dev is the integration branch | Owner | 「gbdrawみたいに、あるていどブランチ？ごとにまとめてからPRしようよ」「devって必要？」「ごめん必要だわ、続けて」 | dev is the integration branch as in gbdraw: work PRs target dev and merge there, and dev goes to main in one promotion pull request. | branch model, pull requests | 1: .claude/skills/losat-ship-pr/SKILL.md:13 | followups-20261014/DECISIONS.md:626 |
| 2026-10-16 | BPCLUT and BPLKPFLAT join the default set | Owner | 「メリットがあるなら入れてください」 | The two blastp-fast lookup-table switches that were added after the default changed join the strict default set when they show a merit. | blastp-fast lookup table switches | 1: docs/agents/diagnostics.md:649 | perf-mb-20261014/DECISIONS.md:299 |
| 2026-10-15 | Reading of NCBI-parity rules 8 and 9 for performance work | Not found | (not found) | Performance work may reorder or parallelise when every value and decision behind the output is unchanged (rule 8), and may do less work than NCBI if complexity does not get worse (rule 9). | performance paths, AGENTS.md rules 8 and 9 | 2: AGENTS.md:201, AGENTS.md:205 | speedup-plan-v4-20261017/PLAN_v4.md:48, v4-P0/pr/AGENTS_md_proposal.md |
| 2026-10-14 | Build for aarch64 and other common platforms (Linux arm64 target) | Owner | 「今後はaarch64とか、ほかにもユーザーが多そうなプラットフォーム向けに作らなきゃね」 | LOSAT is also built for aarch64 and other platforms with many users; Linux arm64 is a required v0.2.0 release target, built natively on an arm64 runner. | platform support, v0.2.0 targets | 2: docs/release/v0.2.0.md:23, docs/release/v0.2.0.md:35 | followups-20261014/DECISIONS.md:11 |
| 2026-10-13 | -db output equals NCBI -num_threads 1 at every thread count | Owner | 「DB-3: NCBI の 1 スレッド。LOSAT はどのスレッド数でも同じ出力」 | LOSAT output is the same at every thread count and for every -mt_mode value, and equals NCBI's -num_threads 1 output. | all programs, -db and -subject | 4: LOSAT/src/algorithm/blastn/args.rs:237, LOSAT/src/algorithm/blastp/args.rs:401 | dbsearch-20261012/DECISIONS.md:290-292 |
| 2026-10-10 | TBLASTX chain sums beyond Int4 not reproduced | Owner | 「Int4は気にせず続けてください」 | TBLASTX linking keeps chain sums in i64 and does not reproduce NCBI's Int4 wrap (C undefined behaviour) beyond INT4_MAX; real inputs stay far below it. | TBLASTX sum-statistics linking | 3: LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs:424, LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs:565 | perf-linkidx-20261010/DECISIONS.md:26 |
| 2026-10-10 | Results page stacked like NCBI classic, not behind tabs | Owner | 「Graphic Summaryの下にDescriptionが直接出てくるやつ。いちいちタブを開かないといけないのはめんどくさい」 (2026-10-10); 「ヒット分布が赤とか緑とか青のリボンで出るボックスとかトップヒットのリストとかをタブで分けるんじゃなくて上下に配置して（NCBIの以前の方式）」 (2026-11-02) | Graphic Summary, Descriptions and the selected subject's Alignments are stacked on one page like NCBI's classic results; Dot Plot, Run details and Outputs stay tabs. | LOSAT Web results page | 2: web/app/src/ui/ResultsPanel.vue:6, changelog.d/feature-web-classic-results.md:7 | web-refresh-20261102/STATE.md:3-8 |
| 2026-10-10 | Gate A skipped when the capture matches the baseline | Owner | 「GateAは省きましょう。規則も変えていいです。私が承諾します」 | Gate A (20 TBLASTX pairs) is not run when the capture of the same commit matches the S02 baseline and the Gate A frozen hashes for all 20 pairs. | TBLASTX full gate tier | 1: .claude/skills/losat-gates/SKILL.md:87 | perf-sa-20261010/DECISIONS.md:70, #124 |
| 2026-10-08 | Gate tiers table and oracle-run cap in the agent skills | Owner | 「Q1,Q2,Q3は推奨通り。Q4はモノによります。コミットされなかったファイルは、どこまで再利用可能か？」 | The Owner approved the recommended agent setup: gate tiers (skill losat-gates) as the default for choosing gates, and the skills including the machine-wide cap on oracle search processes. | agent skills (gates, oracle runs) | 2: .claude/skills/losat-gates/SKILL.md:11, .claude/skills/losat-oracle-runs/SKILL.md:22 | session-audit-20261008/DECISIONS.md:10 |
| 2026-10-06 | Program display names (DW-21) | Not found | (not found) | The Web UI shows the program names BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX, matching NCBI BLAST+ and the LOSAT CLI. | LOSAT Web display names | 1: docs/evidence/losat_web_w3/README.md:8 | the Owner's chat with a working session (no task-folder record) |
| 2026-10-03 | Split parallel sessions and give the full prompt | Owner | 「S09/S10/SXとかのなかで、現在並行して進められるものがあれば進めて構いません。というかそれができるならばセッションを分けましょう。プロンプトの全文を出してください。」 | A session that can run alongside the current one gets its own session; update its instruction file, push, and give the Owner the complete prompt. | campaign workflow (engine and app tracks) | 1: .claude/skills/losat-campaign/SKILL.md:57 | the Owner's chat with a working session (no task-folder record) |
| 2026-10-02 | Treatment of NCBI BLAST+ defects (PD-LOSAT-NCBI-DEFECTS) | Owner | 「妥当な結果が出るならば承認済みの例外」「やっぱ再現のままにして」 | Where NCBI crashes or wraps, LOSAT gives the valid result as an approved exception when it can be checked against a nearby configuration; deterministic NCBI bugs that still give an NCBI result stay byte-reproduced. | NCBI defect handling, BLASTN first | 1: docs/evidence/losat_web_e2g/check_inputs.py:8 | the Owner's chat with a working session (no task-folder record) |
| 2026-10-02 | CLI differences outside search results (DW-13) | Owner | (no words recorded) | Owner accepted the recommended option for the eight items left open by the BLASTN inventory: approved exceptions for argument syntax errors, -help, -subject thread warnings, outfmt 6/7 write failures, out of memory. | CLI non-search behaviour, all programs | 2: docs/evidence/losat_web_e2g/README.md:116, docs/losat_web_gui_sessions/session_s07pppb_e2g_transpile.md:147 | the Owner's chat with a working session (no task-folder record) |
| 2026-10-01 | V-PERF settling rerun uses 5 repetitions | Owner | 「10回は以後5回でいいです」 | The alternating before/after rerun that settles V-PERF cases over +5% uses --repeat 5, not 10. | perf gate (V-PERF) | 2: .claude/skills/losat-gates/SKILL.md:91, docs/evidence/losat_web_e2h/gates/sf_perf.sh:7 | the Owner's chat with a working session (no task-folder record) |

## Decisions cited by log entry

The heading names the log (the work line's folder name, as the repository cites it); Entry is the
entry's number in that log. 123 entries: 15 Owner, 3 Owner in part,
99 Owner-delegated, 6 records.

### `blastdb-ids-20261021`: line DBID

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 2026-10-21 | Owner-delegated |  | Plan answers: ID1 first, NCBI crash inputs refused in LOSAT's words, compact accession guide built from the 2.17.0 table, environment variables accepted only at defaults and refused otherwise (Q8). | blastdb-ids plan (Q1, Q4, Q6-Q9) | `LOSAT/src/algorithm/tblastn/database.rs:527`, `LOSAT/src/algorithm/blastx/database.rs:808` and 3 more |
| 5 | 2026-10-21 | Owner | 「すべて推奨通り」 | Q2 NCBI's gb/emb/dbj type swap is not reproduced; Q3 v5 -parse_seqids LMDB is a faithful port; Q5 default -blastdb_version becomes 5 once ID7a lands. Also approves entry 3. | blastdb-ids Q2, Q3, Q5 | `LOSAT/src/seqid/tests.rs:336`, `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:227` and 1 more |
| 9 | 2026-10-26 | Owner-delegated |  | ID2a: parsed-id DB support for blastp and blastx with LOSAT refusals where NCBI's formatter merges hits or cannot look up a Seq-align id, plus ten more details. | parsed Seq-id databases, ID2a | `LOSAT/src/algorithm/blastx/native.rs:1117`, `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:368` |
| 10 | 2026-10-27 | Owner-delegated |  | ID2b extends the ID2a parts and refusals to blastn, tblastn and tblastx; punctuation-only titles of nucleotide parsed subjects follow approved exception 2. | parsed-id databases, blastn/tblastn/tblastx | `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:253` |
| 13 | 2026-10-27 | Owner-delegated |  | Taxdb messages are written for outfmt 0/5/6/7/17 too; PRE_FETCH_SEQS_LIMIT is not imitated, and an empty taxdb.btd is refused (rc 1) instead of reproducing NCBI's segfault. | taxdb messages, item T | `LOSAT/tests/outfmt16_blastn.rs:529`, `LOSAT/tests/outfmt16_tblastn.rs:605` and 2 more |
| 19 | 2026-11-01 | Owner-delegated |  | Pre-search check of tblastx is bounded for parsed-id databases only, with a constant shared with blastn; larger databases skip the check. | tblastx report streaming bound | `LOSAT/src/algorithm/tblastx/report_stream.rs:309`, `LOSAT/src/algorithm/tblastx/database.rs:756` |
| 20 | 2026-11-01 | Owner-delegated |  | psiblast opens parsed-id databases like blastp, reusing the blastp reporting path with the same refusals; refusals name PSIBLAST. | psiblast parsed-id databases | `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:368` |
| 21 | 2026-11-01 | Owner-delegated |  | Subjects that NCBI's formatter resolves to a different OID (version-less accessions, gb vs emb) are refused for Bioseq-reading formats through a port of SeqidToOid. | parsed-id databases, accession to OID (N1) | `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:368` |

### `cli-compat-20261029`: line CLI

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | 2026-10-29 | Owner-delegated |  | -version prints NCBI's two-line form for the six programs and makeblastdb, read anywhere in argv; -version-full*, -dryrun and -conffile stay refused. | losat CLI -version | `LOSAT/src/blastdbcmd/mod.rs:315`, `LOSAT/tests/blastdbcmd_info.rs:341` |
| 5 | 2026-10-29 | Owner-delegated |  | -mt_mode is accepted (0 to 2) for blastn, blastp, blastx and tblastn, with the same output as -num_threads 1; NCBI's query-split output is not reproduced. | losat CLI -mt_mode | `LOSAT/src/algorithm/tblastn/args.rs:267`, `LOSAT/src/algorithm/blastx/args.rs:400` and 3 more |

### `cli-compat-20261029/C2`: line CLI, stage C2

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 2026-10-29 | Owner-delegated |  | With taxdb found and a name field asked for a subject with nonzero taxid, PR 1 and 2 refuse; PR 3 adds the XML2-style reader. Without taxdb, write N/A and NCBI's warning. | custom tabular taxonomy fields | `LOSAT/src/algorithm/tblastx/report_stream.rs:193` |
| 5 | 2026-10-29 | Owner-delegated |  | Runs with a custom field list (6, 7, 10, 20) do not stream by batch; the report is written once after the search so refusals come before the first byte. | report streaming with custom tabular | `LOSAT/src/algorithm/tblastx/report_stream.rs:89` |
| 6 | 2026-10-29 | Owner-delegated |  | Outfmt 20: write the taxdb warning before the search and the field-name line with the first report, so byte streams match NCBI. | outfmt 20 header and warning | `LOSAT/src/algorithm/common/tabular.rs:7` |
| 7 | 2026-10-29 | Owner-delegated |  | run_local (Web/ABI) has no taxdb search path and treats taxdb as not found; the CLI searches the same path as NCBI; refusal applies only to nonzero looked-up taxids. | custom tabular taxdb path | `LOSAT/src/algorithm/tblastn/args.rs:283`, `LOSAT/src/algorithm/blastx/args.rs:801` and 3 more |

### `dropin-survey-20261029`: the drop-in survey (items DI-n, line BDC)

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 2026-10-29 | Owner-delegated |  | Refuse blastdbcmd combinations whose NCBI output is wrong, with an "is not supported" sentence; four of them were replaced by decision 12, the rc 255 build-path case stays refused. | blastdbcmd (line BDC) | `LOSAT/src/blastdbcmd/taxonomy.rs:71` |
| 9 | 2026-10-29 | Owner-delegated |  | blastdbcmd refuses -outfmt %d and %b (ASN.1 text of the defline set or Bioseq); no downstream tool uses them. | blastdbcmd -outfmt | `LOSAT/src/blastdbcmd/format.rs:113`, `LOSAT/tests/blastdbcmd_retrieve.rs:803` |
| 12 | 2026-10-29 | Owner | 「これはあえて直してお出しするのは？TBLASTX みたいに」 | Four wrong NCBI blastdbcmd combinations become approved differences that output the correct result, built from NCBI's own correct neighbouring path, instead of being refused. | blastdbcmd, PD-LOSAT-BLASTDBCMD-DEVIATIONS | `LOSAT/src/seqdb/lookup.rs:270`, `LOSAT/src/blastdbcmd/fasta.rs:229` and 7 more |
| 31 | 2026-11-01 | Owner-delegated |  | BDC3a: taxid lookups share the DI-8 LMDB readers; -taxids 1 is answered as in DI-8; NCBI rc 255 exceptions that print the build path (non-numeric taxids and similar) are refused with LOSAT sentences. | blastdbcmd taxid options (BDC3a) | `LOSAT/tests/blastdbcmd_batch_line_ends.rs:260` |
| 33 | 2026-11-01 | Owner, in part | 「差 5・6 の承認: 承認します」 | BDC3b reads mask column files; differences 5 and 6 (protein FASTA mask applied to the first third only, stale mask in batch lines) are fixed and presented as approved differences; the Owner approved them. | blastdbcmd masks (BDC3b) | `docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md:15` |
| 34 | 2026-11-02 | Owner-delegated |  | BDC4 lifts the refusal of several names in -db and of quotes in -metadata; the difference 6 extension (mask carried to another sequence of the same Seq-id) was approved by the Owner on 2026-11-02. | blastdbcmd -db "a b" (BDC4) | `docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md:16`, `docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md:18` |
| 36 | 2026-11-02 | Owner | 「承認する」 | The Owner approves the difference 6 extension of the blastdbcmd approved differences (mask carried over to another sequence with the same accession when several names are in -db). | blastdbcmd approved differences, item 6 | related entry (not cited) |
| 44 | 2026-11-03 | Owner | 「(b) 検索前の検査をやめる」 | Drop the pre-search definition-line and ambiguity checks of report streaming; check only the reported subjects at report time, as NCBI does. | NCMEM1, report streaming (RPTSTREAM) | `LOSAT/src/algorithm/blastn/blast_engine/x_rpt_stream.rs:71`, `LOSAT/src/algorithm/blastn/blast_engine/x_rpt_stream.rs:550` and 9 more |

### `dropin-survey-20261029/DI7`: the drop-in survey, item DI-7

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 7 | 2026-10-30 | Owner-delegated |  | Refuse input whose reading changes with molecule-type inference (id accession type conflicts with the program and the letters differ between readings). | -parse_deflines molecule type | `LOSAT/src/blastinput/fasta_reader/believed.rs:15`, `LOSAT/src/blastinput/fasta_reader/believed.rs:427` and 2 more |
| 8 | 2026-10-30 | Owner-delegated |  | Refuse inputs where scope side effects change how another sequence is shown: duplicate query ids, a subject id equal to a query id, or gnl\|BL_ORD_ID\| queries. | -parse_deflines scope | `LOSAT/src/algorithm/common/database.rs:1305`, `LOSAT/src/blastinput/fasta_reader/believed.rs:15` and 2 more |
| 11 | 2026-10-31 | Owner-delegated |  | BLASTN with the same accession as unversioned and versioned queries is refused, since NCBI warns on stderr; other programs pass. | -parse_deflines versions | `LOSAT/src/blastinput/fasta_reader/believed.rs:238` |
| 12 | 2026-11-01 | Owner-delegated |  | -parse_deflines with -subject: refuse equal subject ids, equal query/subject ids and versioned/unversioned accession pairs, because NCBI puts both in one scope. | -parse_deflines -subject (PR3b) | `LOSAT/src/blastinput/fasta_reader/believed.rs:238`, `LOSAT/src/blastinput/fasta_reader/believed.rs:703` and 2 more |
| 15 | 2026-11-01 | Owner | (no words recorded) | A record present in both query and subject (same id, title and letters) passes with -parse_deflines; only combinations where NCBI's display breaks are refused, judged by sequence hash. | -parse_deflines same record | `LOSAT/src/blastinput/fasta_reader/believed.rs:261`, `LOSAT/src/blastinput/fasta_reader/believed.rs:703` and 4 more |
| 16 | 2026-11-01 | Owner | (no words recorded) | Move the id sets used for refusals (QUERY_KEYS, SUBJECT_KEYS) from process-wide statics into a per-search context passed to the reader, in the same PR as decision 15. | -parse_deflines id sets | `LOSAT/src/blastinput/fasta_reader/believed.rs:770` |
| 19 | 2026-11-01 | Owner-delegated |  | Scope of the shared-record case: BLASTX/TBLASTN stay refused (NCBI itself fails), BLASTN/BLASTP/TBLASTX outfmt 0 refuse only when a mask is present, other formats pass. | -parse_deflines same record, scope | `LOSAT/src/blastinput/fasta_reader/believed.rs:326`, `LOSAT/src/blastinput/fasta_reader/believed.rs:703` and 7 more |
| 21 | 2026-11-02 | Owner-delegated |  | After independent review: gb\|A.1\|N and gb\|\|N pairs are refused (NCBI scope finds one for the other), and gim\|n is refused because NCBI stops with rc 3. | -parse_deflines review fixes | related entry (not cited) |

### `dropin-survey-20261029/DI8`: the drop-in survey, item DI-8

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 2026-10-29 | Owner-delegated |  | When taxonomy4blast.sqlite3 is missing, port NCBI as is: one warning, search the whole database unrestricted, exit 0. | DI-8 taxid options | `LOSAT/src/blastinput/db_filter.rs:304` |
| 3 | 2026-10-29 | Owner-delegated |  | An accession with more than 500 OIDs does not reproduce NCBI's lookup error; LOSAT returns the first duplicate, as an approved difference. | DI-8b seqid list, PD-LOSAT-NCBI-DEFECTS | `LOSAT/src/seqdb/filter.rs:1102`, `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:205` |
| 6 | 2026-10-29 | Owner-delegated |  | Alias keys SEQIDLIST, TAXIDLIST, GILIST, OIDLIST and STATS_* stay refused in DI-8; only TITLE and DBLIST of distribution aliases are read. | DI-8 alias keys | `LOSAT/src/seqdb/filter.rs:322` |
| 7 | 2026-10-29 | Owner-delegated |  | A restricted database is presented inside SeqDb as ordinals over the included OIDs (OID bitmap plus rank table), not as holes, zero lengths or an OID-to-ordinal Vec. | DI-8 database restriction design | `LOSAT/src/seqdb/filter.rs:9` |
| 10 | 2026-10-29 | Owner-delegated |  | Replaces 9: -taxids 1 and other taxids on a parent cycle are searched by collecting descendants while skipping visited taxids, rather than refused. | DI-8a taxid expansion, PD-LOSAT-NCBI-DEFECTS | `LOSAT/src/seqdb/tax_sqlite.rs:817`, `LOSAT/src/seqdb/tax_sqlite.rs:833` and 11 more |
| 11 | 2026-10-29 | Owner | 「NCBI と同じ動作のまま（決定 2）」 | Owner confirmed decision 2: without the sqlite file, behave as NCBI, with one warning and no restriction. | DI-8 taxid options | `LOSAT/src/blastinput/db_filter.rs:304`, `LOSAT/tests/db_search_manifest.tsv:1589` and 4 more |
| 17 | 2026-10-30 | Owner-delegated |  | When a positive seqid list removes every definition line of a kept OID, LOSAT refuses with a sentence instead of silently dropping its hits as NCBI does. | DI-8b seqid list, PD-LOSAT-NCBI-DEFECTS | `LOSAT/src/seqdb/filter.rs:376`, `LOSAT/src/seqdb/filter.rs:2911` and 5 more |
| 18 | 2026-10-30 | Owner-delegated |  | A text list whose ids are all gi ids: on a v5 database NCBI warns "Empty seqid list" and searches the whole database, which LOSAT reproduces; on v4 nothing is excluded. | gi-only seqid lists (DI-8b) | `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:281`, `LOSAT/tests/db_search_manifest.tsv:1694` and 2 more |
| 20 | 2026-10-30 | Owner-delegated |  | The missing -seqidlist file case (rc 255, exception text with the build path) is not a manifest row; it is checked in a cli_v2 test by exit code and sentence parts. | -seqidlist missing file test | `LOSAT/tests/db_seqidlist.rs:1` |
| 21 | 2026-10-30 | Owner-delegated |  | With a user list, measure the minimum subject length as NCBI does (otherwise the constant 10), because it enters the Spouge cutoff of the preliminary stage. | DI-8b list statistics | `LOSAT/src/algorithm/blastx/subjects.rs:181` |
| 22 | 2026-10-30 | Owner-delegated |  | A negative list naming OIDs with only gi/ti definition lines is not refused up front; it is refused when a hit lands on such an OID, keeping hit-free searches byte-identical to NCBI. | negative seqid list (DI-8b) | `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:332` |
| 24 | 2026-10-31 | Owner-delegated |  | A database with 0 OIDs goes through the normal path with NCBI's statistics; the only refusal is -dbsize with -outfmt 11 and 0 OIDs. | empty subset databases (DI-8c) | `LOSAT/tests/db_search_manifest.tsv:1738`, `LOSAT/tests/db_search_manifest.tsv:1758` and 3 more |
| 29 | 2026-11-01 | Owner-delegated |  | DI-8f: accept several names in -db for all five programs, resolved as NCBI does, with reports kept per name. | DI-8f multiple -db names | `LOSAT/src/seqdb/tests.rs:1181` |
| 31 | 2026-11-01 | Owner-delegated |  | Mixing v5 databases with and without the accession table (no -parse_seqids) is refused before the search, because NCBI fails depending on the hits; the Owner is asked to confirm in the PR body. | mixed v5 accession tables (DI-8) | `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md:345` |

### `gzinput-20261103`: line GZ

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 2026-11-03 | Owner-delegated |  | A compressed input is recognised by a name ending in .gz (case-insensitive), as NCBI makeblastdb does, not by gzip magic; - is never decompressed. | gzip input recognition | `LOSAT/src/blastinput/compressed.rs:3`, `docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md:113` and 1 more |
| 2 | 2026-11-03 | Owner-delegated |  | Only gzip is read; .bz2 and .zst stay refused naming gzip as the supported format. | compressed inputs (superseded by 8) | `LOSAT/src/blastinput/compressed.rs:3` |
| 3 | 2026-11-03 | Owner-delegated |  | Search programs inflate a .gz query or subject into an unnamed temporary file first and read it as a plain file; a damaged gzip is refused before the search. | gzip input mechanism | `LOSAT/src/blastinput/input_files.rs:70`, `LOSAT/src/blastinput/input_files.rs:444` and 2 more |
| 4 | 2026-11-03 | Owner-delegated |  | Concatenated gzip members (bgzip) are read whole, as NCBI makeblastdb does. | gzip members | `LOSAT/src/blastinput/compressed.rs:3`, `LOSAT/src/blastinput/compressed.rs:286` and 2 more |
| 5 | 2026-11-03 | Owner-delegated |  | A .gz that is not gzip or is cut short stops before any output with an explicit sentence, exit 1; NCBI's exit 6 text is not reproduced. | compressed input errors, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES | `LOSAT/src/blastinput/compressed.rs:310`, `LOSAT/src/blastinput/input_files.rs:70` and 5 more |
| 8 | 2026-11-04 | Owner, in part | (no words recorded) | .bz2 and .zst join .gz for search inputs and makeblastdb (supersedes 2), using the pure-Rust bzip2 (libbz2-rs-sys) and ruzstd crates. | bz2 and zst inputs | `LOSAT/src/blastinput/compressed.rs:3`, `docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md:113` and 1 more |
| 9 | 2026-11-04 | Owner-delegated |  | Concatenated bzip2 streams and zstd frames are all read; NCBI reads only the first and drops the rest silently, which LOSAT does not reproduce. | bz2/zst multi-stream, PD-LOSAT-NCBI-DEFECTS | `LOSAT/src/blastinput/compressed.rs:3`, `LOSAT/src/blastinput/compressed.rs:286` and 6 more |
| 12 | 2026-11-04 | Owner | 「3もやろうよ」 | The inflated copy of a compressed -query / -subject is written beside the input file, falling back to TMPDIR, so a large query does not use RAM-backed /tmp. | gzip temporary copy location | `LOSAT/src/blastinput/input_files.rs:84`, `LOSAT/src/blastinput/input_files.rs:489` and 1 more |

### `makedb-v5-20261029`: `losat makeblastdb` version 5 databases

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 2026-10-29 | Owner-delegated |  | LMDB page size equals the OS page size as in NCBI (4096 on Linux x86-64, 16384 on macOS arm64), and 4096 on wasm32. | makeblastdb LMDB page size | `LOSAT/src/writedb/writer.rs:686`, `LOSAT/src/writedb/platform.rs:116` |
| 9 | 2026-10-29 | Owner-delegated |  | S1 ports mdb.c at the buffer level; paths survey A found unreachable are internal errors; NCBI's LMDB environment variables are read as makeblastdb does. | makeblastdb LMDB writer | `LOSAT/src/writedb/mdb/mod.rs:12` |
| 10 | 2026-10-29 | Owner-delegated |  | NCBI v5 abort cases (gi-only record before another; duplicate in the first record) are refused with rc 1 naming the record, not reproduced. | makeblastdb v5 aborts | `LOSAT/src/writedb/tests.rs:1541`, `LOSAT/src/writedb/writer.rs:1007` and 1 more |
| 11 | 2026-10-29 | Owner-delegated |  | Error paths match makeblastdb in exit status and final message, without byte parity for leftover files or exception texts containing NCBI build paths. | makeblastdb error paths | `LOSAT/src/writedb/tests.rs:1406`, `LOSAT/src/seqid/parse.rs:54` |
| 12 | 2026-10-29 | Owner-delegated |  | .ncbirc LONG_SEQID and MAX_SEQID_LENGTH and the makeblastdb scan-bioseq variable are unsupported: an explicit error when set. | makeblastdb registry settings | `LOSAT/src/writedb/makeblastdb.rs:475`, `LOSAT/src/writedb/makeblastdb.rs:1683` and 2 more |
| 20 | 2026-10-29 | Owner-delegated |  | With several -in files, an unreadable file is an error (exit 255) that removes what earlier files added; -in "" is refused. | makeblastdb -in files | `LOSAT/src/writedb/tests.rs:1235` |
| 28 | 2026-10-30 | Owner-delegated |  | A FASTA id token ends at a byte >= 0x80 as with x86-64 signed char, on every platform. | makeblastdb FASTA id parsing | `LOSAT/src/writedb/lmdb.rs:635` |
| 29 | 2026-10-30 | Owner-delegated |  | CSeqIdException reports name the official Linux release build tree like LOSAT's other NCBI exception reports; the test gate masks the tree. | makeblastdb exception text | `LOSAT/src/writedb/makeblastdb.rs:635`, `LOSAT/src/writedb/tests.rs:1385` |

### `ncbi-parity-20261019`: line NP

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 7 | 2026-10-19 | Owner | 「-max_hsps とか -max_target_seq? とかも絶対必要だと思う」 | -max_hsps for tblastn and tblastx joins item O1, and -num_descriptions / -num_alignments (refused by all five programs) becomes a new item O1b. | ncbi-parity O1/O1b | `docs/losat_web_gui_plan.md:108` |
| 10 | 2026-10-20 | Owner-delegated |  | O1b details: argument errors are clap errors (rc 2), hit list 0 gives NCBI's "No hits are being saved", the "Examining 5 or more" notice follows the hit list, and Web refuses -num_* for now. | -num_descriptions / -num_alignments | `docs/losat_web_gui_plan.md:108` |
| 17 | 2026-10-22 | Owner | 「広い型になおしましょう」 / 「よろしい」 | -culling_limit arithmetic uses i64 instead of reproducing NCBI's Int4 overflow, done in the O3b PR for blastn, blastx and tblastx. | -culling_limit, hspfilter_culling | `LOSAT/src/algorithm/blastx/culling.rs:19`, `LOSAT/src/algorithm/blastx/results.rs:1015` and 4 more |
| 19 | 2026-10-22 | Owner | 「これがまずダメだよな。ネイティブで使えるものはWebでも使えなきゃ。次のセッションで改めよう」 | Whatever the native CLI accepts must also work in the Web adapter; the Web refusal of -num_descriptions / -num_alignments is to be fixed, and other native-only options are to be listed. | Web / native option parity | `docs/losat_web_gui_plan.md:108` |
| 22 | 2026-10-25 | Owner | 「この方針（D3 は代案に変える）で実装に進めてよい」 | LOSAT Web accepts -num_descriptions / -num_alignments and splits the search only when outfmt 0 and 6/7 hit lists differ; -num_alignments 0 behaves as native; design O/webopt D1, D2, D4-D6 adopted. | Web -num_* (O-webopt) | `docs/losat_web_gui_plan.md:108` |
| 23 | 2026-10-25 | Owner-delegated |  | Plan after three surveys: unified P stays refused, BLASTMAT matrices are no item, blastp CBS 0 moves to SC4; item h brings the threaded-Wasm 65,535-page shared-memory commits (#136/#137) into dev as WW0. | ncbi-parity scoring/masking plan; item h = WW0 | `docs/evidence/wasm_memory_max_20261012/dev_ww0/README.md:3` |

### `ncbi-parity-20261019/M`: line NP, stage M (WindowMasker)

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 9 | 2026-10-20 | Owner-delegated |  | Refuse an alias whose volumes have differing algorithm sets or a missing mask-algorithm column, and a big-endian-only mask volume; NCBI reads uninitialised or wrong data there. | window-masker column reader (series M, M4) | `LOSAT/src/blastdbcmd/retrieve.rs:185` |
| 14 | 2026-10-25 | Owner-delegated |  | Where NCBI exception reports are copied, write the official Linux x64 2.17.0 tarball build root; conda builds differ only in that root. | NCBI exception text, NCBI_BUILD_ROOT | `LOSAT/src/utils/winmask/ncbi_report.rs:30`, `LOSAT/tests/blastn_window_masker.rs:3` and 1 more |

### `ncbi-parity-20261019/WW`: line NP, stage WW0 (the Wasm memory maximum on dev)

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 2026-10-25 | Owner-delegated |  | WW0 brings in the ten commits of #136 / #137 with cherry-pick -x, not rewritten; only the DW-28 plan row is added. | WW0 shared memory commits | `docs/evidence/wasm_memory_max_20261012/dev_ww0/README.md:5` |
| 2 | 2026-10-25 | Owner-delegated |  | In the wasi_thread_host.js conflict both dev's blockingStdioWrite and #137's unsignedWasiImports are kept, wrapping the unsigned fd_write. | WW0 conflict | `docs/evidence/wasm_memory_max_20261012/dev_ww0/README.md:5` |
| 3 | 2026-10-25 | Owner-delegated |  | The three test files of #136 are adapted to dev's app: high-memory worker uses loadReactor, E2E runs four V-BR cells at high and normal memory, unit compares serial reactor heaps. | WW0 tests | `docs/evidence/wasm_memory_max_20261012/dev_ww0/README.md:5` |
| 4 | 2026-10-25 | Owner-delegated |  | The WW0 gate runs locally without Node; the full V-ABI, capture and engine gates are skipped because only threaded link arguments change, checked by comparing the reactor sections. | WW0 gate | `docs/evidence/wasm_memory_max_20261012/dev_ww0/README.md:5` |
| 5 | 2026-10-25 | Owner-delegated |  | The local WebKit run uses 16 downloaded Ubuntu 24.04 packages unpacked into a build directory and a launcher script, installing nothing system-wide. | WW0 WebKit check | `docs/evidence/wasm_memory_max_20261012/dev_ww0/README.md:5` |

### `outfmt-20261019`: line OF

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 7 | 2026-10-19 | Owner-delegated |  | SAM @PG CL: carries LOSAT's actual command line instead of a rebuilt NCBI argv; fixture comparisons skip only that item and the difference goes into a product decision. | outfmt 17 (SAM) | `LOSAT/tests/outfmt17_blastn.rs:12`, `LOSAT/tests/db_search_fixtures.py:163` |
| 8 | 2026-10-19 | Owner-delegated |  | Reproduce the outfmt 13 master file name (lists a nonexistent x.json_1.json, no final newline) as defined NCBI behaviour. | outfmt 13 master file, PD-LOSAT-NCBI-DEFECTS | `LOSAT/src/report/blast_xml2/files.rs:27`, `LOSAT/tests/db_search_manifest.tsv:457` and 1 more |
| 10 | 2026-10-19 | Owner-delegated |  | LOSAT Web gets the new output formats after the CLI as a separate item (OF10), with a format list required of ABI v2 and renumbered streams; the v1 API stays frozen. | Web ABI v2, outfmt | `web/adapter/tools/v_abi_cases.py:116` |
| 13 | 2026-10-19 | Owner-delegated |  | Environment: ARCHIVE_FORMAT text and bin in OF8; registry values only at defaults (FastWriteDouble, LONG_SEQID) and refused otherwise. | outfmt environment and registry values | `LOSAT/src/report/serial/real.rs:21` |
| 19 | 2026-10-20 | Owner-delegated |  | OF4 design: eight recommended answers adopted, including the OF4a/OF4b split, refusal of non-ASCII titles, shared XML encoder, and TBLASTX matrix/gap costs written as constants. | outfmt 5 (XML) design | `LOSAT/src/algorithm/tblastx/blast_xml.rs:66`, `LOSAT/tests/outfmt5_blastx.rs:4` and 8 more |
| 20 | 2026-10-20 | Owner-delegated |  | BLASTN report streaming is off whenever outfmt 5 is requested, so every refusal comes before the first byte. | outfmt 5 and RPTSTREAM | `LOSAT/src/algorithm/blastn/blast_engine/run.rs:10043` |
| 21 | 2026-10-20 | Owner-delegated |  | BLASTX outfmt 5 pre-scans query titles before the search and holds the XML in memory only when a subject title would be refused, so every non-ASCII refusal comes before the first byte. | BLASTX outfmt 5 | `LOSAT/tests/outfmt5_blastx.rs:467`, `LOSAT/tests/outfmt5_blastx.rs:622` and 3 more |
| 22 | 2026-10-20 | Owner-delegated |  | BLASTX outfmt 5 with a query that is not a regular file is refused before the search, because the pre-scan reads the file twice. | BLASTX outfmt 5 pipe queries | `LOSAT/src/algorithm/blastx/native.rs:663`, `LOSAT/tests/outfmt5_blastx.rs:466` and 3 more |
| 23 | 2026-10-20 | Owner-delegated |  | Keep TBLASTX displayed_rows and TBLASTN aligned_sequences unchanged; XML gets its own _into copies so default paths stay identical. | outfmt 5 TBLASTX/TBLASTN row helpers | `LOSAT/src/algorithm/tblastn/stage_e_report.rs:1572` |
| 26 | 2026-10-21 | Owner-delegated |  | OF5 design: Q1-Q2 and Q4-Q11 adopted as recommended (shared model, non-ASCII refusal, 13/14 CLI-only until ABI v2); Q3 was superseded by entry 27. | outfmt 13-16 design (OF5) | `LOSAT/src/report/blast_xml2/framing.rs:128`, `web/adapter/tools/v_abi_cases.py:106` and 1 more |
| 28 | 2026-10-21 | Owner | 「入れてください。許可します。」 | LOSAT reads NCBI's taxdb and writes the taxonomy names NCBI writes (13-16 sciname, tabular name fields, 18), reproducing defined edge behaviour and refusing otherwise. | taxdb reader and names | `LOSAT/src/report/blast_xml/records.rs:127`, `LOSAT/tests/outfmt13_14_files.rs:962` and 4 more |
| 30 | 2026-10-21 | Owner-delegated |  | The taxdb reader lands with its first consumer (outfmt 13-16 sciname); formats printing no names do not open taxdb; records past the last taxdb.btd page, empty taxdb.btd and non-ASCII scinames are refused. | taxdb reader, outfmt 13-16 | `LOSAT/tests/outfmt16_blastn.rs:6`, `LOSAT/tests/db_search_manifest.tsv:320` and 14 more |
| 32 | 2026-10-22 | Owner-delegated |  | BLASTX 15/16 with taxdb: decode the DB taxids once before the first byte and hold output only when some taxid would refuse; otherwise write incrementally. | BLASTX outfmt 15/16 taxdb refusal | `LOSAT/src/algorithm/blastx/native.rs:694`, `LOSAT/src/algorithm/common/database.rs:3117` and 6 more |
| 34 | 2026-10-22 | Owner-delegated |  | Answer "can a taxdb lookup refuse?" from the smaller side (walk the taxdb index or decode DB taxids), bounding the pre-pass cost. | taxdb pre-pass cost | `LOSAT/src/algorithm/blastp/x_query_blocks.rs:158`, `LOSAT/src/algorithm/blastp/blast_engine/database.rs:1766` |
| 35 | 2026-10-22 | Owner-delegated |  | BLASTX 13/14 checks each side file where NCBI would open it; a side-file write error maps to "BLAST failed to write output", rc 6 in all programs. | outfmt 13/14 side files | `LOSAT/src/algorithm/blastx/blast_xml.rs:496`, `LOSAT/src/algorithm/blastx/blast_xml.rs:1067` and 8 more |
| 38 | 2026-10-22 | Owner-delegated |  | OF-BPQB2 review: plan arm admits only Single15/16, taxid pre-check shares the title range visit, and both PRs use the smaller-side taxdb check. | blastp query-block outfmt 15/16 | `LOSAT/src/algorithm/blastp/x_query_blocks.rs:158`, `LOSAT/src/algorithm/blastp/blast_engine/database.rs:1766` |
| 44 | 2026-10-23 | Owner-delegated |  | blastp outfmt 11 takes the query-block search from the start, joining the existing strict-set switch LOSAT_X_BPQBLOCKSA with a shared run-level ArchiveHandle. | blastp outfmt 11 archive | `LOSAT/src/algorithm/blastp/blast_engine.rs:6376` |
| 46 | 2026-10-23 | Owner-delegated |  | BLASTX strict-set flow is bounded by letters, not batches, for large BATCH_SIZE, keeping max RSS at or below NCBI; the query title prescan runs per batch; archive rows differing only in algorithm-options are dropped. | BLASTX BATCH_SIZE, outfmt 11 manifest | `LOSAT/tests/outfmt0_manifest.tsv:9`, `LOSAT/tests/outfmt0_manifest.tsv:10` |

### `perf-si-20261012`: performance stage S-I (the strict set as the default)

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 175 | 2026-10-17 | Owner | 「LOSAT_X_TBNSCANBF を既定の組に入れる。同じ分岐なしの走査を blastp / blastx / tblastx にも入れる」 | LOSAT_X_TBNSCANBF joins the strict set, and the scalar branch-free scan is added to blastp, blastx and tblastx, each opt-in first, measured against its own scan, then strict set. | branch-free scan (TBNSCANBF, BXSCANBF) | `docs/agents/diagnostics.md:7088` |

### `release-v020-20261025`: line RL (the v0.2.0 release)

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 52 | 2026-11-01 | Owner-delegated |  | gbdraw's intake of the threaded-Wasm memory.copy trap in Chromium becomes a v0.2.0 GO item; the fix is a build step that guards browser-bound threaded artifacts (guardSharedMemory), worked in a separate session. | v0.2.0 GO, threaded Wasm guard | `web/app/tests/e2e/guard-measure.spec.ts:2`, `web/app/tests/e2e/engine-guard.spec.ts:2` and 7 more |
| 54 | 2026-11-01 | Owner, in part | 「これを v0.2.0 の GO の妨げにする」 「その案でいいです」 「GO 項目を外す」 「x64 の自然な率の扱い: OK」 「Android: 受け入れる」 | The rare threaded-Wasm worker out-of-bounds trap blocks the v0.2.0 GO until a fix PR merges and repeated runs show zero traps; lifted 2026-11-02 after #703, x64 natural rate accepted, Android accepted. | v0.2.0 GO, wasm_oob | `docs/evidence/wasm_oob_20261101/README.md:3`, `docs/evidence/wasm_oob_20261101/README.md:267` and 1 more |

### `sf-e2h-20261008`: LOSAT Web stage E2h (the TBLASTX and BLASTP engine entries)

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 2026-10-08 ABI v1 | 2026-10-08 | Owner-delegated |  | ABI v1 keeps the bio reader and its checks and enters run_local through FastaRecord::from_bio, so v1 outputs and messages do not change (planner Q1). | Web ABI v1 input reading | related entry (not cited) |

### `v4-A`: plan v4 line A

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 19 | 2026-10-19 | Owner-delegated |  | x_unit_records keeps its own body and the sliced path calls a copy, because a call into a shared helper changed inlining of the default scan. | tblastx unit scan, x_records_for | `LOSAT/src/algorithm/tblastx/x_pair_par.rs:1613` |
| 53 | 2026-10-26 | Owner-delegated |  | Session 14 opens with its own PR for the MEM s19 n1 row, LOSAT_X_BXPOOLN1: a one-thread BLASTX search builds lookup tables without the pool; emptying the pool after the last batch was rejected. | LOSAT_X_BXPOOLN1 | `docs/agents/diagnostics.md:7282` |

### `v4-A14`: plan v4 item A14

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 14 | 2026-10-23 | Owner-delegated |  | The cross scan (Avx2Cross) becomes the default for fused tables on AVX2+POPCNT CPUs: LOSAT_X_BXCROSS joins the strict set; =0 or LOSAT_X_DEFAULT=0 restores the old form. | LOSAT_X_BXCROSS | `docs/agents/diagnostics.md:7006` |
| 26 | 2026-10-24 | Owner-delegated |  | Deferred candidate calls with an in-band marker (variant B) are a go after the probe counts, proceeding to a prototype. | LOSAT_X_BXCROSSDEFER probe | `docs/agents/diagnostics.md:7054` |
| 27 | 2026-10-24 | Owner-delegated |  | The prototype is a go and a default-on PR follows; form r3 (marker handled without calls in the batch path) is chosen. | LOSAT_X_BXCROSSDEFER PR | `docs/agents/diagnostics.md:7054` |

### `v4-C`: plan v4 line C

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 25 | 2026-10-19 | Record |  | Probe numbers of the blastp composition redo on WSSV x PajaWSV: erfc, composition reads and Karlin parameter costs and repeated inputs, the basis for step C4. | blastp redo C4 probe | `docs/agents/diagnostics.md:1008`, `docs/agents/diagnostics.md:1143` |
| 31 | 2026-10-19 | Record |  | Record of probe counts for tblastn/blastx redo (SEG trim dominates, repeated trim inputs 118,765 calls of 12,584 kinds). | C5 SEG trim memo evidence | `LOSAT/src/utils/x_seg_trim_memo.rs:33`, `docs/agents/diagnostics.md:950` |
| 32 | 2026-10-19 | Record |  | Record of the trim memo sizing: cap at 32,768 entries, doubled as needed, cleared at the cap. | LOSAT_X_SEGTRIMMEMO size | `LOSAT/src/utils/x_seg_trim_memo.rs:56`, `docs/agents/diagnostics.md:950` |
| 33 | 2026-10-19 | Owner-delegated |  | Go for step C5 = LOSAT_X_SEGTRIMMEMO: a per-thread memo of SEG trim results keyed by the whole byte string. | SEG trim memo (C5) | `docs/agents/diagnostics.md:950` |
| 34 | 2026-10-19 | Owner-delegated |  | After review the SEG trim memo is opt-in per masker (tblastn redo windows, blastx query) so other maskers keep no table and RSS stays at or below NCBI. | SEG trim memo scope | `docs/agents/diagnostics.md:950` |

### `v4-D12`: plan v4 item D12

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 59 | 2026-10-27 | Owner-delegated |  | bp_wssv -subject n1 stays outside the memory floor after #576, so the u16-hit block-table prototype (BPLK16) goes ahead; the floor is recounted as 12,365 KiB. | blastp memory, BPLK16 | `docs/agents/diagnostics.md:1997` |
| 60 | 2026-10-27 | Owner-delegated |  | The remaining ur50e3 n32 memory is live data (thread scratch, window block results, lookup tables), not freed arena space. | blastp memory census | `docs/agents/diagnostics.md:1997` |
| 61 | 2026-10-27 | Owner-delegated |  | LOSAT_X_BPLK16 wins and goes into the strict set in the same PR (perf/bp-lk16, default on, =0 restores the old path). | LOSAT_X_BPLK16 | `docs/agents/diagnostics.md:1997` |
| 62 | 2026-10-28 | Owner-delegated |  | No-go on shrinking init-HSP capacity after a large subject; the peak is the sum of per-thread maxima, so next borrow large arrays from a shared pool. | blastp memory (init HSP arrays) | `LOSAT/src/algorithm/blastp/blast_engine/x_ih_pool.rs:25`, `docs/agents/diagnostics.md:1904` |
| 63 | 2026-10-28 | Owner-delegated |  | The shared initial-HSP pool saves 15-22 MiB on ur50e3 n32 at no time cost, so LOSAT_X_BPIHPOOL proceeds to the strict set. | LOSAT_X_BPIHPOOL | `docs/agents/diagnostics.md:1904` |
| 66 | 2026-10-28 | Owner-delegated |  | Restricted-gap arrays grow to the reached band as NCBI's dp_mem does (about 18-20 MiB saved on ur50e3 n32), so LOSAT_X_BPGAPBAND proceeds to the strict set. | LOSAT_X_BPGAPBAND | `docs/agents/diagnostics.md:5928` |
| 67 | 2026-10-28 | Record |  | The board's "O0 cycles x1.057 after #601" splits into x1.005-1.01 from #601 and binary placement. | BPIHPOOL cost analysis | `docs/agents/diagnostics.md:5928` |
| 68 | 2026-10-29 | Owner-delegated |  | The first BPGAPBAND form is not made a PR (O0 n1 cycles x1.0019); the comparison moves into the existing column-copy branch to keep the row loop identical to dev. | BPGAPBAND implementation form | `docs/agents/diagnostics.md:5928` |
| 69 | 2026-10-29 | Owner-delegated |  | Row functions stay character-identical to dev; the band limit is passed through the col length and the arrays double when a row stops. | BPGAPBAND implementation form | `docs/agents/diagnostics.md:5928` |
| 70 | 2026-10-30 | Owner-delegated |  | The d78bb8a8 form passes the rule, so a docs commit and a draft PR (perf/bp-gap-band) follow. | BPGAPBAND PR | `docs/agents/diagnostics.md:5928` |
| 71 | 2026-10-30 | Owner-delegated |  | Census shows small thread residues of trees and per-query lists; first implement LOSAT_X_BPHITSFLAT (counting sort, no per-query lists), trees later. | blastp memory (prelim hit lists) | `LOSAT/src/algorithm/blastp/blast_engine/x_hits_flat.rs:24`, `docs/agents/diagnostics.md:1938` |
| 72 | 2026-10-30 | Owner-delegated |  | LOSAT_X_BPHITSFLAT works and goes to a full gate and draft PR in the strict set; the n1 instruction increase of x1.0003 is accepted. | LOSAT_X_BPHITSFLAT | `docs/agents/diagnostics.md:1938` |

### `v4-W`: plan v4 line W

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 33 | 2026-10-20 | Record |  | TBXHITBF wins on Wasm in all cases (scan x0.62-0.75), so neither a rule nor a PR is needed. | TBXHITBF on Wasm | `docs/agents/diagnostics.md:6287` |
| 36 | 2026-10-20 | Record |  | Same-process pair measurement of the wasm32 split of x_hits_bf: scan stage x0.92-0.95, no losing case, so a PR follows. | wasm32 x_hits_bf split | `docs/agents/diagnostics.md:6287` |

### `wasm-oob-20261101`: line WO

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 15 | 2026-11-02 | Owner-delegated |  | Clearing check for the wasm_oob GO item: x64 widened module dc-megablast and blastp n16 at N=40 each, and ARM64 blastp n16 at N=100; no natural-rate x64 run is made. | wasm_oob clearing runs | `docs/evidence/wasm_oob_20261101/README.md:267`, `docs/evidence/wasm_oob_20261101/scripts/wo_clear.sh:2` |

### `web-audit-20261104`: line WA

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 7 | 2026-11-04 | Owner-delegated |  | The web audit work is split into five sessions (W1-W5) by touched files; user-visible "verification/certified" wording is removed from screens, titles, notices and READMEs, internal names and past records stay unchanged. | LOSAT Web wording, session split | `docs/web/ncbi_ui_mapping.md:227` |

### `web-sx-blastx-20261102`: line SX

| Entry | Date | Kind | Owner's words | Decision | Scope | Cited in |
| --- | --- | --- | --- | --- | --- | --- |
| 7 | 2026-11-02 | Owner | 「SXd を始める」 「Certified profile っていう表現やめた方がいいと思う。Certify って別に公式じゃないじゃん。使えます、だけでいいじゃん」 「既存の不具合 2 つ：なおしましょう」 | Owner starts SXd (BLASTX -query_loc / -subject_loc), drops "certified/certification" from the verification badge, and asks that two existing defects (a screens spec wait, an operator-precedence bug in v_abi_cases.py quick()) be fixed. | Web SXd, wording, two bug fixes | `web/adapter/tools/v_abi_cases.py:245` |
| 10 | 2026-11-03 | Owner-delegated |  | SXd1 shape for BLASTX ranges: NCBI behaviour for empty ranges, UTF-8 option text as other programs, formats 0/6/7/10 only, outfmt 0 lowercase quirk reproduced, error order per NCBI. | BLASTX -query_loc/-subject_loc (SXd1) | `LOSAT/src/cli.rs:1191` |
