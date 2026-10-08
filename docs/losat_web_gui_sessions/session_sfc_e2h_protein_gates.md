# Session SFc — E2h の続き：蛋白の入力の移植（S7〜S10）、fixture の仕上げ、ゲート

## INSTRUCTION PROMPT

LOSAT の段階 E2h を SFb から続ける。指示書の正本は [SF の指示書](session_sf_e2h_blastn_fasta_reader.md)（背景、NCBI の経路、手順 1〜9、「保守者の判断」1〜5、終了・引き継ぎ）と [SFb の指示書](session_sfb_e2h_port.md)（手順 0〜6、ゲートの規模、進め方）である。この指示書は、SFb で終わったことと残りの順を足すだけである。完了条件は総合計画書 §7 の SF の行である。先に [セッション README](README.md) の共通規則を読む。エンジン側のセッションで、アプリ側（S13 の続き）と並行する（計画 DW-7）。

セッションは Linux の clone `/home/kawato/losat-work`（`$WORK_ROOT`）で起動する。そうすると `.claude/settings.local.json` の env、hook、agent（`losat-implementer`・`losat-mechanic`・`losat-reviewer`・`losat-test-runner`）が読み込まれる。SF と SFb は別の folder で起動したので、これらが無く、env を手で設定した（`/home/kawato/losat-baselines/sfb-e2h-20261008/env.sh`）。

SFb の結果（[ゲート記録](../evidence/losat_web_e2h/README.md) の「SFb で行ったこと」）：移植 S1〜S6 で BLASTN・TBLASTX・TBLASTN の subject が `CFastaReader` の移植で読む（`7e49a1be`、`e4de0af6`、`ebaeeaf5`、`6c5d4220`、`8f8e8e08`、`c7061c08`、`5fbdf4c1`）。fixture は BLASTN same 1969 と Seq-id 12、TBLASTX same 625 と `-lcase_masking` 16 と Seq-id 2、TBLASTN same 297 と query の側 191 と Seq-id 4、BLASTP は変更前のまま（same 243、rejects 567、Seq-id 6）。受け付ける入力で変わった出力は `Subject_` の sseqid だけ。sweep の script（`fasta_sweep.py`、`check_inputs.py`、NCBI 側は固定済み）は evidence にある。ゲートの script の下書きは `docs/evidence/losat_web_e2h/gates/`（未実行）。推奨の案で進めた判断はゲート記録の 11〜15 と `/home/kawato/losat-baselines/sfb-e2h-20261008/DECISIONS.md`。これに従う。

0. **準備**
   - 規則 2 の確認（worktree `$WORK_ROOT/.worktrees/web-gui`、`feature/losat-web-gui`、SFb の最後のコミットから）。
   - タスクフォルダ `/home/kawato/losat-baselines/sfc-e2h-<yyyymmdd>/` を skill `losat-campaign` の手順で作り、セッションを bind する。SFb の `STATE.md`・`DECISIONS.md`・`handoff-S1〜S6.md`（`/home/kawato/losat-baselines/sfb-e2h-20261008/`）を読み、`DECISIONS.md` は新しい folder に写して追記する。段階ごとの検査の script（`scripts/s5_standard.sh`、`s5_compare_fixtures.py`、`s5_replay/`）も同じ folder にある。
   - アプリ側の S13 が終わり、`feature/losat-web-gui-app` が push されていれば `feature/losat-web-gui` へ merge する（README 規則 1。衝突の解消は別のコミット）。アプリ側のゲート記録に書かれた README の表と計画の変更も反映する。S13 が続いていれば、区切りのよいところでまた確かめる。
1. **移植 S7〜S10**（`PORT_PLAN.md` §2 の順。`losat-implementer` に 1 段階ずつ任せ、1 体ずつ）
   - S7 BLASTP（両方の役割が蛋白）：`check_inputs.py` の `EXPECT_OVERRIDE_PORTED` に `blastp` を足す。v1 は今の検査を前の位置で行う（S5 の `V1Records` の形）。
   - S8 TBLASTN の query（蛋白）と後片付け：`blastn/input.rs` の拒否の関数、`InputRecord` の bio の実装、`fasta_input` の別名を消す。`grep` で `fasta::Record` が `web_api.rs`（v1）、アダプタの種類 0、`algorithm/blastx/` にだけ残ることを示す。`BLASTINPUT_GEN_DELTA_SEQ` を確かめ直す（RD の 6 実行）。
   - S9 `scan` の種類 1（核酸の旗）と 2（蛋白の旗）と性質試験 `web/adapter/tests/scan_ncbi_properties.rs`。
   - S10 `register`・`run`・`docs/web/abi_v2.md`（`register` の行と §8・§9、種類 1 は「BLASTX、SX」でなくなる）、`v_abi.js` の 239〜241 行、V-ABI full の case に SF の入力の場合を足す。ABI v2 の file の出力で、最初の batch の失敗の前に outfmt 0 の prolog が書かれるようになった（S5）。これを V-ABI で確かめる。
   - 各段階の検証は SFb と同じ：quick と standard の段階、`fasta_input_fixtures.py check --rows <program> --jobs 3`（その段階で移植した種類の行は `same` になり、ほかの行は変わらない。変わる行は理由を記録）、2 つの sweep、棚卸しの oracle の再生。移植した箇所の直上に NCBI のファイル・行と断片を書く。BLASTX（`algorithm/blastx/`）は変えない。
2. **fixture の仕上げ**（SF の指示書の手順 5 の残り。S8 の後）
   - Seq-id の 24 行に、LOSAT の明示的な拒否の文言（§G4。`not supported by LOSAT's <PROGRAM>` を含む）と終了コード（1）を固定する。NCBI では流さない。
   - TBLASTX と BLASTP の `-lcase_masking` の 34 行は、E2e の明示的な拒否（`AUTHORITY.md` §A、終了コード 2）として期待を直す（SFb の判断 11）。
   - `LOSAT/tests/ci_fast_regressions.py` に program ごとの `fasta_input_fixtures.py check` を足し、path の絞り込みに `LOSAT/tests/fasta_input_fixtures*` と `LOSAT/tests/fixtures/fasta_input/` を足す。全行が `same` か明示的な拒否になってから足す。
   - 変更前の `check` の結果（`$BUILD_ROOT/sf-e2h/fixtures/check-before.tsv`）と移植の後の結果をゲート記録に並べる（S6 の後の結果は `/home/kawato/losat-baselines/sfb-e2h-20261008/fixtures-after-S6/`）。
3. **以降**は SF の指示書の手順 6〜9 をそのまま行う：sweep（`fasta_sweep.py` の `LISTED_REJECTIONS` を最後の文言に合わせる。全 program で、記録の無い拒否・差・時間切れが 0）、ゲートの全工程、V-PERF、独立監査（4 観点とも supported まで）。
4. **ゲートの規模**（SFb の指示書の 4. のまま）
   - S11 で省いた工程（TBLASTX の option の sweep、capture、Gate A の 20 組、V-ABI full の全部）も行う。工程は並行させず順に流し、各工程を `flock "$BUILD_ROOT/oracle.lock"` の下で `--jobs` を明示して流す（V-ABI full と capture は 4、sweep は 3、TBLASTN・TBLASTX の option の sweep は 2）。
   - 下書き `docs/evidence/losat_web_e2h/gates/`（`sf_gates.sh`・`sf_gates_resume.sh`・`sf_gate_a.sh`・`sf_perf.sh`・`NOTES.md`）は未実行である。先に `STAGES=fast-all` だけで試す。TBLASTN の Stage G の case は `/mnt/c/Users/genom/GitHub/LOSAT/...` の path を名指すので、下書きは `STAGE_G_BIND=1`（`unshare -rm` で `$WT` をその path に重ねる）を入れた。動かなければ、case の path を `$WT` に読み替える方法を選び、`/mnt/c` は読まない。
   - Gate A は `/tmp/losat-pr5-runtime-cert-*` を読む。WSL の再起動で消えるので、ゲートの前に置き直す。
   - capture の変更前の基準は `~/.cache/losat-web-gui-target/sf/capture-before/hashes.tsv`、変更前の成果物は `~/.cache/losat-web-gui-target/sf/bin/` である。
   - V-PERF：読み込みの重い case（10 万レコードの query、5 Mb の 1 行と 80 桁の subject、BLASTP の多数の query）の入力の大きさは、ネイティブで 1 回測ってから決める。S3 から `-`（標準入力）を NCBI の `cin` と同じく 1 byte ずつの補充で読むので、標準入力の case も足す。lock は `~/.cache/losat-web-gui-target/s07p-resume/vperf_lock.sh`。
   - 待つ間は watcher 1 本で、30 分以内の間隔で見る。`$BUILD_ROOT/sf-e2h/` と `$BUILD_ROOT/sfb-e2h/`（fixture と sweep の NCBI の出力、S0 の速度の入力、段階ごとの成果物）は、ゲート記録が参照するので、この段階が `main` に入るまで消さない。
5. **NCBI の不具合**：`AUTHORITY.md` §K の 12 件は、各行の推奨の規則で移植する（SFb までに (1) RP-20 は明示的な拒否、(4) file と標準入力の尾の消失は両方を再現）。保守者の答えがこのセッションの中で来たら従う。来なければ、最終回答で一括の問いとしてまた示す。
6. **区切り**：コンテキストが 300k に近づいたら、残りを SFd に分けてよい（README の表に行を足す）。区切る前に、`STATE.md`・ゲート記録・次の指示書を更新して push する。アプリ側の worktree（`.worktrees/web-gui-app`）とそのビルドには触れない。

最終回答は README 規則 8 のとおり（コミット SHA、push の結果、検証の結果、残件、次の指示書のパスと全文）。推奨案で進めた判断は `DECISIONS.md` とゲート記録に書き、最終回答に一覧で示す。終わったら、SF の指示書の「終了・引き継ぎ」（計画、`abi_v2.md`、`verification_cells.tsv`、README の表、S17 と SX の指示書、アプリ側への引き継ぎ）を行う。

## 進め方（SFb で分かったこと）

- 移植は `losat-implementer` に 1 段階ずつ任せる。メインは設計の判断と結果の確認だけをする。implementer は自分で quick と standard の段階を流し、コミットして push する（同じ worktree でコミットする agent は 1 体だけ）。ほかの agent は file を書くだけにし、コミットはメインが path を指定して行う（`git commit -- <path>`）。
- implementer の中から別の agent を起こさせない（SFb の S5 では、見えない worker と 2 体目の implementer が重なった）。TBLASTX のような大きい段階は、コンテキストが約 250k を超えたら、確かめた部分をコミットして区切るよう指示する。
- implementer が NCBI で確かめられなかった端の場合は、メインが `unshare -rn` の下で少数流して決着させる（SFb で 4 件、すべて一致）。v1 の WASI の行列は NCBI を自分で流すので `unshare -rn` の下で流す。
- 検索を流す agent は同時に 2 体まで。始める前に `uptime` を見て、load average が 12 を超えていれば待つ（アプリ側の Playwright の E2E で 9〜13 になる）。
- 新しいビルドは `$BUILD_ROOT/<用途>` に置き、既存の名前（`native`、`test`、`wasi`、`adapter`、`gate-sf`）を使い回す。段階ごとの成果物は `$BUILD_ROOT/sfb-e2h/S<n>/` の形にする。
