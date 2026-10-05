# Session S08+b — E2e の仕上げ：S08+a の merge、最後のゲート、独立監査の第 2 回、V-PERF、`main` への PR

## INSTRUCTION PROMPT

LOSAT の段階 E2e（BLASTP・TBLASTN・TBLASTX の既定以外のオプション、TD-13）を仕上げる。S08+ は移植・fixture・sweep・第 1 回の独立監査とその修正（`d96412265` まで）を終え、最後のゲートと記録を進めていた。並行の S08+a（別の worktree・ブランチ `feature/losat-web-gui-s08pa`）は、第 1 回の監査の残件（TN-2・TN-4・TN-5・RP-4・`-out -version`）と、監査の再現で見つかった TN-1 の残りを NCBI と同じにした（ブランチの先頭 `c6fb13915`、エンジンの最後のコミット `68154c73e`）。このセッションでは、S08+a を本線に merge し、申し送りを片付ける。そのうえで merge 後の最後のコミットで、全ゲート・Gate A・V-ABI・v1 の WASI の行列・独立監査の第 2 回・V-PERF を通し、記録を仕上げて、S08+ の完了条件を満たす。完了条件の正本は総合計画書 §7 の S08+ の行で、[S08+ の指示書](session_s08p_e2e_protein_options.md)がその細部である（条件を緩めない）。

先に読む：
- [セッション README](README.md) の共通規則（特に規則 4 と 8）
- [S08+ の指示書](session_s08p_e2e_protein_options.md)、[S08+a の指示書](session_s08pa_e2e_audit_open_items.md)
- S08+a の記録 [`s08pa/NOTES.md`](../evidence/losat_web_e2e/s08pa/NOTES.md)（merge の前は `git show origin/feature/losat-web-gui-s08pa:docs/evidence/losat_web_e2e/s08pa/NOTES.md`）
- [第 1 回の独立監査](../evidence/losat_web_e2e/audit/ROUND1.md)、[E2e の権威の記録](../evidence/losat_web_e2e/AUTHORITY.md) の §M（判断 D1〜D12）と §N（残り）
- S08+ が書いていれば、E2e のゲート記録 `docs/evidence/losat_web_e2e/README.md` と、S08+ の作業の覚え書き `~/.cache/losat-web-gui-target/s08p/NOTES.md`

### 作業場所と約束

- **作業場所とブランチ**：worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui`、ブランチ `feature/losat-web-gui`。
  - コミットと push の前に、毎回 `git pull --rebase origin feature/losat-web-gui` を行う。merge のコミットを作った後は、merge を平らにしないよう `git pull --ff-only` にする。
  - アプリ側の worktree `LOSAT-web-gui-app`（S09）には触れない。
  - S08+a の worktree `LOSAT-web-gui-s08pa` は読むだけにする（merge の後に消すかは保守者に任せる）。
- **ツールと NCBI**：
  - `RUSTUP_TOOLCHAIN=1.92.0`。`--target-dir` は `~/.cache/losat-web-gui-target/s08pb-*`（ゲートの script が使う `s08p-gate-*` はそのまま使ってよい）。
  - NCBI のソースは固定 commit `598d8ae6`（`/mnt/c/Users/genom/GitHub/ncbi-blast/c++`）。同じ commit の写し `~/.cache/losat-web-gui-target/s08p/ncbi/c++` を読む方が、9p の誤りに当たらない。
  - 比較用の NCBI BLAST+ 2.17.0 は `/home/kawato/micromamba/bin`。`-remote` は使わない。
- **WSL の 9p**：S08+ の最後のゲートの実行（`run-20261004T034405Z`）は、`/mnt/c` の Input/output error で壊れた（下の「S08+ の状態」）。
  - agent（Sonnet）には `/mnt/c` を読ませない。必要な写しは `~/.cache/losat-web-gui-target/s08pb/` に置く。
  - ゲートの途中で I/O の誤りが出たら、その段をやり直す。繰り返すなら止めて保守者に知らせる（`/mnt/c` の作り直しは、同じ計算機の他のセッションに及ぶ）。
- **agent の約束**：
  - `pkill`・`killall`、名前や `-f` の条件での `kill` を使わせない。止めてよいのは、自分が起動した PID だけ（S08+a で、agent の `pkill -f` が他の agent の実行を止めた）。
  - 機械的な作業（比較の実行、sweep、記録の下書き、監査）は、Agent ツールに `model: "sonnet"` で回す。読み取り専用の agent は、4 つを超えて並行してよい。
  - エンジンの変更、ゲート、V-PERF は 1 つずつ行う。
  - 長い実行は 30 分以内の間隔で見る。
- **V-PERF の lock** `/home/kawato/.cache/losat-web-gui-target/vperf.lock`：
  - このセッションが V-PERF の間に置く（`s08p_perf.sh` が `~/.cache/losat-web-gui-target/s07p-resume/vperf_lock.sh` で取る）。
  - 他が置いている間は、ビルド・試験・比較を始めない。
- **大きなファイルをコミットしない**：
  - stage は path を指定して行う（`git add -A` を使わない）。
  - コミットの前に stage した blob の大きさを確かめる。1 MB を超えるものがあれば、止めて保守者に諮る。
  - sweep の作業ディレクトリ、NCBI と LOSAT の生の出力、実行ファイル、trace は `~/.cache/losat-web-gui-target/` に置く。記録に入れるのは要約とログだけで、大きさは既存の `sweeps/*.tsv` までにする。
- **`/tmp`**：要らないファイルは消してよい。ただし `/tmp/losat-pr5-runtime-cert-*`（Gate A の字句の入力）と使用中のものは消さない。

### 0. 状態の確かめ

1. `ls /mnt/c/Users/genom/GitHub/LOSAT-web-gui` が通ること。
2. S08+ のセッションが終わっていること。次の 3 つを確かめる。S08+ が作業中なら、始めずに保守者に知らせる。
   - 本線に未コミットの変更が無い（`git status --short`）。S08+ の壊れた run の directory `docs/evidence/losat_web_e2e/run-20261004T034405Z/` は、未追跡で残っている場合がある。
   - `ps aux | grep -E 's08p_gates|s08p_perf'` に、実行中のものが無い。
   - `vperf.lock` が無い。
3. S08+b の指示書の扱い：
   - S08+ が S08+b の指示書（`session_s08pb_*.md`）を書いていれば、それを正本として読む。この指示文の S08+a の部分（1.、2.、記録の該当箇所）を足して行う。
   - 書いていなければ、この指示文を `docs/losat_web_gui_sessions/session_s08pb_e2e_close.md` に置く。README の表には S08+b の行を足し、S08+a の行を完了にする。
4. S08+ の記録の状態を確かめ、S08+ が済ませていない段を、このセッションの作業に入れる（`git log --oneline d96412265..HEAD`、`docs/evidence/losat_web_e2e/README.md` の有無、S08+ が終えた段）。
5. 未追跡の `run-20261004T034405Z/` は記録に入れない。S08+ がその後に使っていないことを確かめてから、`~/.cache/losat-web-gui-target/s08p/aborted-run-20261004T034405Z/` に移す。

### 1. S08+a の merge

1. `git fetch origin` の後、`git merge --no-ff origin/feature/losat-web-gui-s08pa` を行う。S08+a は `d96412265` から分かれた。本線の側には、セッションの指示書 `16ac97583` と S08+ の記録がある。
   - 衝突は、S08+a が `LOSAT/tests/outfmt0_manifest.tsv` の末尾に足した 10 行（`e2e.tblastn.hardmask_*`・`b45_frame_ties*`・`query_split*`・`window_end`、`e2e.blastp.query_split*`）を残す形で解く。
   - 解き方は merge のコミットの本文に書く。
2. merge の後、S08+a の確かめの一部を、この worktree の build で繰り返す。
   - `check_losat.py` の全 fixture：147 件。S08+ が後で足した行があれば、その分を加える。
   - `cargo test --all-features`：S08+a の最後は 909 件が通過した。

### 2. S08+a からの申し送り

1. **アダプタの試験の失敗**：`web/adapter/src/store.rs:253` の `store::tests::blastn_register_rejects_records_that_ncbi_reads_differently` が、`d96412265` から失敗している。
   - 原因：S08+ の IN-10（`a2b43ee12`）で、BLASTP の `register` が定義行の tab を拒否するようになった。それでもこの試験は、`register("blastp", ..., b">q\tt\nACXGT\n")` の成功を前提にしている。
   - 直し方：BLASTP の登録を、tab の無い定義行（例 `b">q t\nACXGT\n"`）にする。試験だけの変更で、BLASTN の検査の意図は変えない。
2. **NCBI の参照の誤り 2 件**（`verify-refs-session-added.txt`）：
   - `LOSAT/src/cli.rs:15` の `cmdline_flags.cpp:107-143`：例えば `const string kArgQuery("query");` が範囲に無い。
   - `cli.rs:172`（merge の前の行。S08+a のくくり出しで `:209` から移った）の `tblastn_args.cpp:63-110`：`m_PsiBlastArgs.Reset(new CPsiBlastArgs(CPsiBlastArgs::eNucleotideDb));` が範囲に無い。
   - NCBI の `598d8ae6` の正しい行に直し、`verify_refs.py` と `verify_added.py`（BASE `78c06fe61`）で誤り 0 にする。
3. **TN-5 で変えなかった BLASTP の `sort_unstable_by`**（`blastp/hsp.rs`、`blast_engine.rs`。行は NOTES の「残り」）：
   - これらは NCBI の `qsort`（glibc 2.39 では安定な merge sort）に当たる並べ替えである。
   - 比べ方が全てを区別しないものについて、同順位の順が出力に届くかを NCBI のソースで確かめる。届きうるなら、安定な並べ替えにする。対象は BLASTP だけが使う箇所である。
   - BLASTX と共有する `core/composition_adjustment/redo_alignment.rs` は変えず、SX の指示書に書く（DW-10）。
4. **保守者に諮る項目に S08+a の 2 つを足す**（S08+a は推奨の案で進めた）：
   - (a) NCBI が落ちる query の分割の組は、BLASTP・TBLASTN で明示的な拒否のままにした。組は、chunk をもう一度分けるほどの重なりと、batch を分ける負の `CHUNK_SIZE` である。BLASTN では前者を承認済みの例外 1 で検索しているので、それを BLASTP・TBLASTN に広げるかを諮る。
   - (b) option の値の位置にある toolkit の語（`-out -version` など）も、明示的に拒否した（D9 の範囲）。

### 3. 最後のゲート

エンジンを変えずにゲートを始める。途中で直すことが見つかったら、ゲートを止めて直し、コミットしてから、新しい run の directory でやり直す（ゲートとエンジンの変更は 1 つずつ）。

1. run の directory を新しくする：`echo "docs/evidence/losat_web_e2e/run-$(date -u +%Y%m%dT%H%M%SZ)" > ~/.cache/losat-web-gui-target/s08p/rundir`。
2. `docs/evidence/losat_web_e2e/gates/s08p_gates.sh > ~/.cache/losat-web-gui-target/s08p/gates-final.log 2>&1` を background で実行する。
   - 実行の前に、script の option sweep の `--timeout 300` を `1200` にする。NCBI の blastp・tblastn の `-word_size 6 -threshold 1` と、tblastx の `-word_size 4 -threshold 1` は数分かかる。tblastn は 1 つで約 5 GB の記憶を使う。S08+ の run と S08+a の 1 回目の sweep では、この組の NCBI の側が timeout になった。BLASTP の組は、NCBI が `malloc(): corrupted top size` で落ち、LOSAT は拒否する。
   - 記憶が足りなければ、sweep の `--jobs` を下げる（32 core、31 GB）。

   期待する結果は次のとおり（S08+ の `d96412265` の run のうち壊れていない段と、S08+a の実測による）。
   - **lint と試験**：
     - fmt、clippy の 4 構成とアダプタの 3 構成が通る。
     - `cargo test --all-features` は 909 件以上が通過し、失敗 0。
     - アダプタの試験は 5 件が通過（2.1 の後）。
     - wasm32 の web API の試験と pure-Rust の境界が通る。`protein-tables-check.log` は終了 0。
     - `verify-refs-session-added.txt` の誤りは 0（2.2 の後）。
   - **fixture**：
     - `check-losat-n{1,2,4}.tsv` は 147 件で差 0。`oracle-check-gate.log` は差 0。
     - `tblastx-fixtures.log` は 84 件で差 0。`blastn-fixtures.log` は 183 件で差 0。
     - `ctoolkit-compare.tsv` は差 0。`punct-defline.tsv` は予期しない 0。
     - `title-sweep.tsv`：tblastx・tblastn とも 1023 の定義行で、一致 957、例外 2 が 66、予期しない 0。
     - `blastn-check-inputs.tsv` は 300 件で予期しない 0。
     - `fast-regressions-all.log` は 236 件で失敗 0。許可した既知の不一致は 1 件（Sakai.MG1655.megablast）。この検査は Gate A と TLOSAN Stage G の凍結ハッシュを含む。
   - **sweep**（`sweeps/after-*.tsv`）：DIFF 0、timeout 0。
     - BLASTP 1196 件：一致 231、LOSAT の拒否 633。
     - TBLASTN 1427 件：一致 258、拒否 720、承認済みの `-db_gencode` の例外 78。
     - TBLASTX 641 件：S08+ の run の DIFF 95 のうち 88 件は、ファイルを読めない誤り（`File is not accessible`）だった。残る 7 件は原因を確かめていない。内訳は、`-threshold +inf` で LOSAT が signal 9 で終わった 1 件と、`-db_gencode` 2・5・6・9 などの 6 件である。S08+ のその前の sweep では、TBLASTX は DIFF 0 だった。再び出たら、NCBI のソースで原因を確かめる。
     - `gencode-api-check.tsv` は終了 0。`protein-title-sweep.tsv` は一致 1014・拒否 69。
   - **capture**：`capture-compare.txt` と `capture-compare-before.txt` を見る。
     - 前者の基準は `docs/evidence/losat_web_e1a/baseline/hashes.tsv`（S08+ が BLASTP の 11 行を更新した）。後者の基準は SD の最後の成果物。
     - 差は、NCBI と同じ出力に変えた行だけであることを、NCBI の出力で確かめて記録する。
     - 基準の hash の列を更新するのは、NCBI とバイト一致の行だけ（S08+ と同じ）。
   - **v1 の WASI の行列**：`wasm-threading.log` の形式の失敗 0、`v1-requests-compare.txt` 一致。
   - **V-ABI quick と full**：失敗した部分 0。凍結ハッシュとの不一致は、既知の Sakai の 4 件だけ。S08+ の run の `tblastx-34` の失敗は、I/O の誤りの中で起きた。再び出たら原因を確かめる。
3. Gate A：`~/.cache/losat-web-gui-target/s08/s08_gate_a.sh` の形で、`LOSAT/tests/audit_tblastx_v010.py --losat-bin ~/.cache/losat-web-gui-target/s08p-gate-native/release/LOSAT --output-dir <scratch>` を実行する。出力は run の directory の `audit-tblastx-v010*` に置く。TBLASTX の outfmt 6 の v0.1.0 の parity が変わらないこと。
4. 期待と違う結果は、原因を NCBI のソースで確かめてから、直すか記録する（推奨の案で進め、判断を記録する）。

### 4. 独立監査の第 2 回

最後のゲートが native と reactor を作った後に、AUTHORITY §N の 4 観点で、sonnet の agent を読み取り専用で並行して動かす。

- **観点**：S08+ が分け方を決めていなければ、第 1 回の 5 つの報告を (a) BLASTP、(b) TBLASTN、(c) TBLASTX、(d) 報告と入力 の 4 つにまとめる。
- **指示と作業場所**：
  - 指示は `~/.cache/losat-web-gui-target/s08pb-audit/COMMON.md` と `ANGLE_*.md` に書く。形は S08 の `~/.cache/losat-web-gui-target/s08-audit/COMMON.md` に倣い、写しを `docs/evidence/losat_web_e2e/audit/` に置く。
  - 作業ディレクトリは `~/.cache/losat-web-gui-target/s08pb-audit/<観点>/`。比べる実行ファイルは、ゲートの native の写し。

agent には次を求める。
- 第 1 回の `audit/round1/*.md` の全ての再現の命令と harness を、最後の実行ファイルで繰り返す。S08+a が `/mnt/c` の外に直した harness は `~/.cache/losat-web-gui-target/s08pa/audit-rerun/<観点>/` に、その報告は `s08pa/audit_rerun/` にある。
- S08+a の変更の周りを新しく調べる：
  - BLASTP・TBLASTN の query の分割：9,800 残基・19,800 残基を超える query、複数の query の batch、環境変数 `CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`・`BATCH_SIZE`、lower-case の mask、`-comp_based_stats` の値、少ない `-max_target_seqs`。
  - 巨大な X-drop の時間と記憶。
  - BLOSUM45 で同じ得点の HSP の frame。
  - hard mask と少ない `-max_target_seqs`。
  - 窓の右端の traceback：`-comp_based_stats 0`、巨大な最終の X-drop、大きい `-evalue` の組。
  - option の値の位置の toolkit の語。
  - BLASTP の同順位の並べ替え（2.3）。
- 結果を途中でファイルに書き（`.md` の Write が拒まれたら `FINDINGS.txt`）、最後の返答に全指摘を書く。
- `pkill` などを使わない（上の約束）。

4 つとも supported になるまで、指摘を直してから、その観点を再び監査する。エンジンの変更は 1 つずつ行い、変えたら 3. の該当する検査をやり直す。保守者の判断待ちの項目（D11、D12、2.4 の 2 つ）は、「pending the maintainer」として扱わせる。記録は `audit/ROUND2.md` と `audit/round2/` に置く。

### 5. V-PERF

監査が終わり、計算機が静かなときに、`docs/evidence/losat_web_e2e/gates/s08p_perf.sh` を実行する。

- **比べるもの**：
  - 変更前は SD の最後の成果物（`~/.cache/losat-web-gui-target/sd-final-native/release/LOSAT`、`sd-final-wasi-artifacts`）。変更後は 3. の成果物。
  - case は既定の `blastp,blastp-fmt0,tblastn,tblastn-fmt0,tblastx,tblastx-multi,tblastx-many,blastn-large`。
- **出力が同じかを先に確かめる**：S08+ は、outfmt 0 の既定の 500 の説明・250 の整列や蛋白の題などで、BLASTP・TBLASTN の出力を NCBI に合わせて変えた。出力が違う case は、速度の比較から外す理由を記録するか、同じ出力の入力に替える。
- **閾値**：`measure_perf.py` の ×1.05 を超えた case は、`REPEAT=5` で測り直す。
- **参考に記録するもの**：
  - S08+a の native だけの測定（`s08pa/checks/perf-final*.json`）：BLASTP は 2 回とも約 1〜2% 遅い側に出たが、上限の内。
  - TBLASTN の `e2e_many_subject.fna`（300 subject）は、LOSAT が約 2〜4 秒、NCBI が 0.2〜0.4 秒。SD の成果物から同じなので、退行ではない（残件として記録する）。

### 6. 記録

- **E2e のゲート記録** `docs/evidence/losat_web_e2e/README.md`（E2b・E2i の記録の形）：
  - 状態（完了と日付）。
  - 完了条件の表：計画 §7 の S08+ の行の各条件について、結果と根拠。
  - 「ゲート」：最後の run の表。
  - 「V-PERF」、「独立監査」（第 1 回と第 2 回）、「S08+a」（NOTES の要約とコミット）、「保守者に諮ること」、残件。
  - 途中で止めた run は記録に入れない。
- **監査と権威の記録**：
  - `audit/ROUND1.md` の TN-2・TN-4・TN-5・RP-4・BP-8（`-out -version`）の行の対応を、「直した（S08+a、コミット、`s08pa/NOTES.md`）」に直す。TN-1 の残り（S08+a の再現で見つかり、`68154c73e` で直した）の行を足す。
  - `AUTHORITY.md` の §M に、S08+a の 2 つの判断（2.4）を足す。§N は、このセッションの後の残りにする。
- **`docs/evidence/losat_web_e2e/evidence.sha256`**：E2g と同じ形で、記録の directory の全ファイルを入れる（run の directory と `s08pa/` を含む）。
- **計画と README の表**：総合計画書 `docs/losat_web_gui_plan.md` の冒頭の状態の行と、§7 の S08+ の行の状態。[セッション README](README.md) の表の S08+・S08+a・S08+b の行。
- **[S12 の指示書](session_s12_w3_search_ui.md)**：S08+ が `f7ce8199f` で書いた option の値を、S08+a で変わった点に合わせる。
  - BLASTP・TBLASTN の長い query は、分割して NCBI と同じに検索する。
  - `-version` などの toolkit の語は、値の位置でも拒否する。
- **[SX の指示書](session_sx_blastx_integration.md)**：BLASTX の項目を足す。
  - BLASTX と共有する gapped DP の確保（TN-4。出力は変えていない）と、`align_ex_protein` の `read_end_sentinel`（BLASTX は偽）。
  - `redo_alignment.rs` の `sort_unstable_by`。
  - NCBI は blastx の 1 つの query の batch も分割する（`split_query_aux_priv.cpp` の `SplitQuery_ShouldSplit`、重なり 297）。これを LOSATX が移植しているか。

### 7. 保守者の判断

次をまとめて諮る（推奨は、どれも記録のとおり）。
- D11
- D12（推奨の別案の承認済みの例外を含む）
- 2.4 の 2 つ
- V-PERF で上限を超えた case があれば、その扱い

判断が出たら、`PD-LOSAT-NCBI-DEFECTS`・`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` と、AGENTS.md の承認済みの例外の記述を、保守者の承認の範囲で更新する。判断が出なければ、今の明示的な拒否のまま、残件として渡す。

### 8. コミット、push、`main` への PR

1. エンジンの変更（2.1〜2.3 と監査の修正）と記録は、別のコミットにする（アプリは変えない）。
2. `git pull --rebase origin feature/losat-web-gui`（merge の後は `--ff-only`）の後に push する。
3. `main` への PR：
   - SD が作った PR が開いていれば、S08+・S08+a・S08+b のコミットはそれに入る。本文に E2e の要約を足す。
   - 閉じていれば、新しく作る。本文の終わりは `🤖 Generated with [Claude Code](https://claude.com/claude-code)`。
   - merge はしない（保守者が行う）。
4. PR の CI（`rust`、fast output regressions）が緑であることを確かめる。

完了条件は計画 §7 の S08+ の行による。記録は `docs/evidence/losat_web_e2e/README.md`。

## S08+a からの引き継ぎ（2026-10-04 の実測）

- **ブランチ**：`origin/feature/losat-web-gui-s08pa` の先頭は `c6fb13915`（起点 `d96412265`）。
  - エンジンのコミット：
    - TN-2 `36572eace`、TN-4 `20856865e`、TN-5 `5c5132987`
    - RP-4 の拒否の版 `7d77ee408`（後の移植で置き換えた）
    - `-out -version` `6fa06cb6a`
    - rustfmt と参照の行 `f3c8d978d`・`b9b2a6f06`
    - RP-4 の query の分割の移植：TBLASTN `44ccde85c`、BLASTP `61d3bf4f5`
    - TN-1 の残り `68154c73e`
  - 記録のコミット：`c394f1381`・`c6fb13915`（`docs/evidence/losat_web_e2e/s08pa/` だけ、33 ファイル・約 840 KB）。
- **確かめ**（`68154c73e`、`s08pa/checks/`）：
  - fixture は 147 件が 1・2・4 スレッドで差 0（変更前の実行ファイルは、足した 10 件だけが違う）。
  - TBLASTN と BLASTP の sweep は DIFF 0・timeout 0（`--timeout 1200`）。
  - fmt、clippy の 4 構成とアダプタの 3 構成が通過。`cargo test --all-features` は 909 件通過。
  - アダプタの試験は 1 件失敗（2.1。起点から失敗）。
  - NCBI の参照の誤りは、S08+ の 2 件だけ（2.2）。
  - TN-4：TBLASTN の `-xdrop_gap_final 1e8` で、26.9 秒・3.9 GB が 0.12 秒・8.4 MB になった（`s08pa/tn4/`）。
  - 監査の再現（`6fa06cb6a` で 5 つの agent、`s08pa/audit_rerun/`）で、説明の付かない差は TN-1 の残りだけ。
- **やっていないこと**：全ゲート、V-ABI、v1 の WASI の行列、Gate A、正式の V-PERF、独立監査の第 2 回、capture と fast regressions。全体の監査の再現は、最後の実行ファイルではやり直していない（分割と TN-1 の残りで変わる組は、最後の実行ファイルで繰り返した）。
- **成果物**：
  - `~/.cache/losat-web-gui-target/s08pa/bin/LOSAT-base`（`d96412265`、SHA-256 `a3e54fcf5c8a2a333a555049f5ac842d35a2d66a969ef0cee09c88d9bbc86f84`）。
  - `LOSAT-final`（`68154c73e`、SHA-256 `6f5e268d0765ee649f131b1f888a8bb5c77656f9151d8774d7caea0ff6983389`）。
  - 作業の写しと trace は `~/.cache/losat-web-gui-target/s08pa/`（`tn2/tb_trace.c` など）、sweep の作業は `s08pa-sweep-*`。
- **残件**（NOTES の「残り」）：
  - TBLASTN の 300 subject の速さ（前から）。
  - BLASTP で分割される batch の query を、全体の検索でも一度検索する費用（約 9,800 残基を超える query だけに掛かる）。
  - TN-5 で変えなかった並べ替え（2.3）。
- **事故**：S08+a の BLASTP の再現の agent が、`pkill -f` と `kill` を広い条件で使った。止めたのは、TBLASTX の再現の agent の 1 件と `bash -c` の 2 件。S08+ のゲートのプロセスは見当たらなかったが、止まったものは特定できない。なお S08+ の最後のゲートの実行は、S08+a の開始より前に I/O の誤りで壊れていた。

## S08+ の状態（2026-10-04 に S08+a から見たもの）

- 本線の先頭は `16ac97583`（S08+a の指示書）。
- S08+ の最後のゲートの実行 `docs/evidence/losat_web_e2e/run-20261004T034405Z/`（`d96412265`、未追跡）は、記録には使えない。
  - lint と試験・build・fixture・sweep の段は動いた（`cargo test` 904 件通過、fixture 137 件差 0 など）。
  - `/mnt/c` の I/O の誤りで、`gencode-api-check`・capture・v1 の WASI の行列・V-ABI quick が失敗した。
  - TBLASTX の sweep の DIFF 95 のうち 88 件は、ファイルを読めない誤りだった（残る 7 件は 3.2）。
  - V-ABI full は `tblastx-34` が失敗した。
- S08+ の作業の覚え書き（`~/.cache/losat-web-gui-target/s08p/NOTES.md`）の最後は、「Remaining for this session: final gates, V-PERF, README, S09/S12/SX docs, audit, next instruction (S08+b if conditions remain)」。S12 と SX の指示書への引き継ぎは、`f7ce8199f` で済んでいる。これより後に S08+ が進めたものは、0. で確かめる。

## 終了・引き継ぎ

README の規則 8 に従う。エンジン側の次は [S11 — 領域の指定](session_s11_e2d_query_subject_loc.md)（アプリ側の S09 は並行して続いている）。完了条件が残れば、その解消を次の指示書（S08+c）の最初の作業にする。
