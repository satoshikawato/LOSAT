# Session S08b — E2b の仕上げ：最後のゲート、独立監査の第 3 回、V-PERF、`main` への PR

## INSTRUCTION PROMPT

LOSAT の段階 E2b（TBLASTX の outfmt 0 と 7）の仕上げを行う。S08 は移植・fixture・棚卸し・独立監査の第 1・2 回とその指摘への対応までを終え、最後のエンジンのコミット `1117e8c17` の全ゲートの途中（「build wasi」）で WSL の `/mnt/c` の I/O の誤りに当たり、保守者の指示でそこで区切った。このセッションで、最後のコミットの全ゲート、独立監査の第 3 回（確認）、Gate A、V-PERF を通し、記録を仕上げ、`main` への PR を作って、S08 の完了条件を満たす。完了条件の正本は総合計画書 §7 の S08 の行で、[S08 の指示書](session_s08_e2b_tblastx_outfmt0_7.md)がその細部である（条件を緩めない）。

先に次を読む：
- [セッション README](README.md) の共通規則（特に規則 4 と 8）
- [S08 の指示書](session_s08_e2b_tblastx_outfmt0_7.md)
- [E2b のゲート記録](../evidence/losat_web_e2b/README.md)（状態、完了条件、「ゲート」「独立監査」「保守者に諮ること」「残件と引き継ぎ」）
- [E2b の権威の記録](../evidence/losat_web_e2b/AUTHORITY.md) の §F（決めたこと）

作業場所：worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui`、ブランチ `feature/losat-web-gui`。`RUSTUP_TOOLCHAIN=1.92.0`、`--target-dir` は worktree の外（`~/.cache/losat-web-gui-target/...`）。NCBI のソース `/mnt/c/Users/genom/GitHub/ncbi-blast/c++`（598d8ae6）が唯一の正本、`/home/kawato/micromamba/bin` の NCBI BLAST+ 2.17.0 は比較だけに使う。アプリ側（S09、worktree `LOSAT-web-gui-app`）には触れない。`docs/evidence/losat_web_e2g/README.md` と `evidence.sha256` は S07+++b の記録なので直さない。コミットと push の前に毎回 `git pull --rebase origin feature/losat-web-gui` を行う。

### 0. 状態の確かめ

- `ls /mnt/c/Users/genom/GitHub/LOSAT-web-gui` が通ること（通らなければ止めて保守者に知らせる。`/mnt/c` の作り直しは同じ計算機の他のセッションに及ぶ）。
- `git log --oneline 2bcb86b1f..HEAD` と `git status --short`：エンジンの最後のコミットは `1117e8c17`、その後は記録だけのコミット。未コミットの変更が無いこと。
- エンジンを変えずにゲートを始める。ゲートの途中で直すことが見つかったら、ゲートを止め、直してコミットし、新しい run の directory でやり直す（ゲートとエンジンの変更は 1 つずつ）。

### 1. 最後のゲート

1. run の directory を新しくする：`echo "docs/evidence/losat_web_e2b/run-$(date -u +%Y%m%dT%H%M%SZ)" > ~/.cache/losat-web-gui-target/s08/rundir`。
2. `~/.cache/losat-web-gui-target/s08/s08_gates.sh > ~/.cache/losat-web-gui-target/s08/gates.log 2>&1` を background で実行する（写しは `docs/evidence/losat_web_e2b/gates/`。約 2 時間。進み具合は 30 分以内の間隔で見る）。各段の期待：
   - lint と試験：fmt・clippy（4 構成と adapter の 3 構成）・`cargo test --all-features`（S08 の最後は 875 件通過）・adapter・wasm32 の web API の試験・pure-Rust の境界・`verify-refs-session-added.txt` の誤り 0。
   - fixture：`check-losat-n{1,2,4}.tsv` 74 件差 0、`oracle-check-gate.log` 差 0、`precheck.tsv` の差は `tblastx.code4.0`・`tblastx.code4.7`（承認済みの `-db_gencode` の例外）だけ、`tblastx-fixtures.log` 70 件差 0、`tblastx-fixtures-before.log` と `-0d533ba76.log` は差がある（区別の確認）、`tblastx-env-discrimination.log` 予期しない 0、`blastn-fixtures.log` 107 件差 0、`ctoolkit-compare.tsv` 差 0、`punct-defline.tsv` 予期しない 0、`html-titles.tsv` 不要な拒否 0・差 0、`blastn-check-inputs.tsv` 予期しない 0、`blastn-scoring-sweep-fmt{0,6,7}.tsv` 300 一致・580 同じ誤り、`blastn-title-sweep.tsv` 予期しない 0、`closed-pipe.txt`（outfmt 0 は終了コード 6。閉じた標準出力の行は、保守者に諮っている差のとおり LOSAT が終了コード 0）、`fast-regressions-all.log` 失敗 0（許可した既知の不一致 1：Sakai.MG1655.megablast）。
   - capture：`capture-compare.txt`（S02 の基準）と `capture-compare-before.txt`（セッションの開始の実行ファイル）が 236 件差 0。
   - v1 の WASI の行列：`wasm-threading.log` の形式の失敗 0、`v1-reactor-records-compare.txt`、`v1-requests-compare.txt`。
   - V-ABI quick と full（`v-abi-quick/v-abi.log`、`v-abi-full/run.log`：失敗した部分 0、凍結ハッシュとの不一致は既知の Sakai の 4 件だけ）。
3. `~/.cache/losat-web-gui-target/s08/s08_gate_a.sh`（Gate A：`audit_tblastx_v010.py`、outfmt 6 の v0.1.0 の parity。変わらないこと）。
4. 期待と違う結果は、原因を NCBI のソースで確かめてから直すか記録する（推奨の案で進め、判断を記録する）。

### 2. 独立監査の第 3 回（確認）

最後のゲートが native と reactor を作った後（`gates.log` の「v-abi full」の行の後）に、観点ごとに sonnet の agent 4 つを読み取り専用で並行して動かす。指示は `~/.cache/losat-web-gui-target/s08-audit/COMMON.md`、`ANGLE_{A,B,C,D}.md`、`ROUND3.md`（写しは `docs/evidence/losat_web_e2b/audit/`）。作業ディレクトリは `~/.cache/losat-web-gui-target/s08-audit/r3{a,b,c,d}/`、前の回の作業は `{a,b,c,d}/` と `r2{a,c,d}/` にある。agent には結果を途中でファイルに書かせ（`.md` の Write が拒まれたら `FINDINGS.txt`）、最後の返答に全指摘を書かせる。4 つとも supported になるまで、指摘を直し（エンジンの変更は 1 つずつ、変えたら 1. の該当する検査をやり直す）、その観点を再び監査する。閉じた標準出力（保守者に諮っている差）は「pending the maintainer」として扱わせる。

### 3. V-PERF

監査が終わり、計算機が静かなときに、アプリ側の lock を取って `~/.cache/losat-web-gui-target/s08/s08_perf.sh` を実行する（変更前はセッションの開始の実行ファイル `~/.cache/losat-web-gui-target/s08/bin/LOSAT-base` と `e2g-gate-wasi-artifacts`、変更後は 1. の成果物。case は `docs/evidence/losat_web_e2b/perf_cases.py` の `tblastx,tblastx-multi,tblastx-many,tblastn,tblastn-fmt0,blastp,blastp-fmt0,blastn-large-fmt0,blastn-large`）。先に、比べる case の出力が変更前と変更後で同じかを確かめる（S08 は batch・hit list・曖昧な文字・SEG・cutoff で一部の出力を NCBI に合わせて変えた。違う case は、速度の比較から外す理由を記録するか、同じ出力の入力に替える）。閾値を超えた case は `--repeat 5` で測り直す。outfmt ごとの費用（`perf-formats.txt`）も記録する。

### 4. 検証の升目

`docs/web/verification_cells.tsv` の TBLASTX の outfmt 0/7 の行（26：native CLI、28：Wasm ABI v2 serial / threaded、スレッド 1・2・4）を、1. の結果（fixture の件数、run の directory の名前、V-ABI full）で埋め、状態を `checked`、`filled_by` を `S08` にする。

### 5. 記録

- ゲート記録 `docs/evidence/losat_web_e2b/README.md`：状態（完了と日付）、完了条件の表（各条件の結果と根拠）、「ゲート」（最後の run の表。E2g の記録の形）、「V-PERF」、「独立監査」の第 3 回、残件。途中で止めた run は記録に入れない（`~/.cache/losat-web-gui-target/s08/` に移す）。
- `docs/evidence/losat_web_e2b/evidence.sha256`（E2g と同じ形：記録の directory の全ファイルの SHA-256。run の directory を含む）。
- 総合計画書 `docs/losat_web_gui_plan.md` の冒頭の状態の行と §7 の S08 の行の状態、[セッション README](README.md) の表の S08 と S08b の行。
- S08+ の指示書の「S08 からの引き継ぎ」を、このセッションの実測に合わせて直す。

### 6. 保守者の判断

ゲート記録の「保守者に諮ること」の 2 つ（TBLASTX・TBLASTN の句読点だけの outfmt 0 の題に例外 2 を広げるか、起動の時に閉じた標準出力を承認済みの例外にするか）を、まとめて保守者に諮る（推奨はどちらも記録のとおり）。判断が出たら、`PD-LOSAT-NCBI-DEFECTS` または `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` と AGENTS.md の承認済みの例外の記述を保守者の承認の範囲で更新し、実装は S08+ の最初の作業にする（判断が出なければ、今の明示的な拒否と記録した差のまま、S08+ に残件として渡す）。

### 7. コミット、push、`main` への PR

エンジンの変更とアプリの変更は別のコミットにする（アプリは変えない）。`git pull --rebase origin feature/losat-web-gui` の後に push し、`main` への PR を作る（`4fc67f9ab`・`2bcb86b1f`・S08 と S08b のコミットを含む。本文の終わりは `🤖 Generated with [Claude Code](https://claude.com/claude-code)`）。merge はしない（保守者が行う）。PR の CI（`rust`、fast output regressions）が緑であることを確かめる。

完了条件は計画 §7 の S08 の行による。記録は `docs/evidence/losat_web_e2b/README.md`。

## S08 からの引き継ぎ（2026-10-03 の実測）

- **コミット**（`git log --oneline 2bcb86b1f..HEAD`）：移植 `fd896a52c`、fixture と CI `0d533ba76`、棚卸しの結果の修正 `06953c68d`・`a46dd63fd`・`2c9f71aa7`・`2049103b0`、SEG `aee365eb5`・`968fa98c8`、監査の対応 `7fbbfad96`・`32fd67a52`（注釈）・`aa4b6f8c3`・`6017f0ea2`・`5fa53b3f0`・`1117e8c17`、記録のコミット。
- **済んだゲート**：`1117e8c17` の lint と試験の段（`docs/evidence/losat_web_e2b/run-20261002T162101Z/`、875 件通過）。`968fa98c8` と `7fbbfad96` の部分のゲート（`~/.cache/losat-web-gui-target/s08/gate-968fa98c8-partial/`、`gate-7fbbfad96-partial/`）で、済んだ検査はすべて通った（ゲート記録の「ゲート」）。
- **独立監査**：第 2 回は (b)・(d) が supported、(a)・(c) は直す前の実行ファイル（`7fbbfad96`）で unsupported。その指摘はすべて `aa4b6f8c3`・`5fa53b3f0`・`1117e8c17` で直したか、S08+ に記録したか、保守者に諮った。手元の build での確かめ：監査 (a) 第 2 回の下限の 56 件と環境の 81 件、監査 (b) の再現 12 件、空の query と読めない subject の 8 種、`subprocess.DEVNULL` の実行。第 3 回で最後の実行ファイルについて確かめる。
- **成果物とスクリプト**：変更前の実行ファイル `~/.cache/losat-web-gui-target/s08/bin/LOSAT-base`（SHA-256 `331fcba36447…44d8`）と `LOSAT-0d533ba76`、ゲートの script `~/.cache/losat-web-gui-target/s08/{s08_gates.sh,s08_gate_a.sh,s08_perf.sh,verify_added.py}`、V-PERF の lock `~/.cache/losat-web-gui-target/s07p-resume/vperf_lock.sh`、監査の指示 `~/.cache/losat-web-gui-target/s08-audit/`。
- **進め方**：機械的な作業（比較の実行、sweep、記録の下書き、監査）は Agent ツールに `model: "sonnet"` で回し、読み取り専用の agent は 4 つを超えて並行してよい。エンジンの変更、ゲート、V-PERF は 1 つずつ。V-PERF の測り直しは `--repeat 5`。長い実行は 30 分以内の間隔で見る。保守者に判断を仰ぐことはまとめて諮り、ほかは推奨の案で進めて記録する。
- **注意**：監査の agent の `git status` が worktree の index の lock を取ることがある（`/mnt/c` は遅い）。commit が `index.lock` で失敗したら、数秒おいてやり直す。

## 終了・引き継ぎ

README の規則 8 に従う。次は [SD — BLASTN の dc-megablast](session_sd_e2i_blastn_dc_megablast.md)（保守者の依頼、DW-18）、その次が [S08+ — BLASTP・TBLASTN・TBLASTX の既定以外のオプション](session_s08p_e2e_protein_options.md)。
