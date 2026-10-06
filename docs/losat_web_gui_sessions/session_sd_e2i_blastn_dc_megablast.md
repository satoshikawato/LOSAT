# Session SD — E2i：BLASTN の dc-megablast と blastn-short

## INSTRUCTION PROMPT

LOSAT の段階 E2i を実行する。BLASTN の `-task dc-megablast`（discontiguous megablast）と `-task blastn-short` を NCBI BLAST+ と同じにする。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`.agents/skills/verify-ncbi-parity-and-speed/SKILL.md`）に従う。完了条件の正本は、総合計画書 §7 の SD の行である（計画 DW-18）。エンジン側のセッションで、S08b の後、S08+ の前に行う。

背景（2026-10-03 の調べ）：

- v0.1.0 の CLI の前（`adbcaf7fb` の前、2026-09-11 まで）の `-task` は任意の文字列を受け付け、`dc-megablast` は blastn の既定値（word size 11、reward 2、penalty −3、gap 5/2）を使い、megablast の lookup を強いるだけだった（`coordination.rs` の `discontig_template`）。discontiguous の template（`-template_type`・`-template_length`）、その lookup と走査は無く、NCBI の dc-megablast とは違う結果だった。
- v0.1.0 の CLI（`adbcaf7fb`）は `-task` を `megablast` と `blastn` に絞り、S07+（`c33731d08`）は NCBI のほかの task（`blastn-short`・`dc-megablast`・`rmblastn`）を明示的に拒否した（`LOSAT/src/blastinput/value_parsers.rs` の `blastn_task`）。`-template_type`・`-template_length` は LOSAT に無い NCBI の option として拒否している（`LOSAT/src/cli.rs`）。
- 残っている部品：`LOSAT/src/algorithm/blastn/blast_engine/run.rs` の `discontig_template`、`LOSAT/src/algorithm/blastn/coordination.rs` の lookup の選び方（`mb_template_length > 0` で megablast の lookup を強いる NCBI の分岐、`blast_nalookup.c` の `BlastChooseNaLookupTable` の移植）。どちらも今は届かない。
- 保守者の依頼（2026-10-03）：dc-megablast をできるようにする。続けて blastn-short も。
- blastn-short は NCBI では blastn の program に task の既定値を重ねたもの（`c++/src/algo/blast/api/blast_options_handle.cpp:343-360`：`SetMatchReward(1)`、`SetMismatchPenalty(-3)`、`SetEvalueThreshold(1000)`、`SetWordSize(7)`、`ClearFilterOptions()`）。LOSAT の blastn は word size 4〜28 と NCBI の得点の表を既に NCBI と同じにしている（E2c の `word_size_sweep.py`・`scoring_sweep.py`）ので、残りは task の既定値と、NCBI の CLI がそれと `-evalue`・`-dust`・`-lcase_masking` などの指定を組み合わせる順序（`blastn_args.cpp`、`blast_args.cpp` の `ExtractAlgorithmOptions`）。

NCBI の経路（固定 commit 598d8ae6 のソースで確かめ、表にする。blastn-short は下の「task と既定値」「引数」と、blastn の経路のうち既定値で変わる分岐）：

- task と既定値：`c++/src/algo/blast/api/blast_options_handle.cpp`（`CBlastOptionsFactory::CreateTask` の `dc-megablast`）、`c++/src/algo/blast/api/disc_nucl_options.cpp`（`CDiscNucleotideOptionsHandle`：`SetTemplateType(0)`、`SetTemplateLength(18)`、`SetWordSize(BLAST_WORDSIZE_NUCL)`、`SetWindowSize(BLAST_WINDOW_SIZE_DISC)`（40）、ungapped と gapped の X-drop、gap trigger、`eDynProgScoreOnly` と `eDynProgTbck`、得点は `CBlastNucleotideOptionsHandle::SetScoringOptionsDefaults`）、program `eDiscMegablast`。
- 引数：`c++/src/algo/blast/blastinput/blast_args.cpp:696-760`（`CDiscontiguousMegablastArgs`：`-template_type` は `coding`・`optimal`・`coding_and_optimal`、`-template_length` は 16・18・21、互いに依存）、`blastn_args.cpp`。
- 検査：`c++/src/algo/blast/core/blast_options.c:1242-1270`（`s_DiscWordOptionsValidate`：word size 11・12、template の長さ）と `:1400-1406` の呼び出し。
- lookup：`c++/src/algo/blast/core/blast_nalookup.c`（`GetDiscTemplateType`（606-642）、discontiguous の megablast の表の作り方、2 つの template）、`mb_lookup` の構造（`blast_nalookup.h`）。
- 走査：`c++/src/algo/blast/core/blast_nascan.c:2202-2626`（`s_MB_DiscWordScanSubject_1`、`_TwoTemplates_1`、`_11_18_1`、`_11_21_1`、選び方）。
- 拡張と以後：discontiguous の hit の ungapped 拡張（`na_ungapped.c`）、two-hit の窓、gapped の段階、traceback、報告。BLASTN の既存の移植（megablast と blastn）と違う分岐だけを、E2g の棚卸し（`docs/evidence/losat_web_e2g/INVENTORY.tsv`）の行に照らして洗い出す。
- 報告：outfmt 0 の prolog と epilog、outfmt 7 の見出しの program と task の表記、`Method:` の有無。

1. **棚卸し（計画 DW-12）**：NCBI の dc-megablast の経路の関数を、E2g の棚卸しと同じ方式（`docs/evidence/losat_web_e2g/` の `inventory_refs.py`・`build_inventory.py`・`stage2/`）で表にする（`docs/evidence/losat_web_e2i/INVENTORY.tsv`）。megablast・blastn と共有で E2g が faithful とした行は、dc-megablast で通る分岐が同じかだけを確かめる。読み取り専用の sonnet の agent に範囲ごとに分けてよい。
2. **一括の移植**：未移植と差のある移植を、NCBI の関数ごとに簡略化せずに transpile する。Rust の移植箇所の直上に NCBI のファイル・行と断片を書く。速度のための LOSAT の方式（詰めた配列の走査、並列化）は、出力が同じなら使ってよい。`blastn_task` で `dc-megablast` と `blastn-short` を受け付け、`-template_type`・`-template_length` を NCBI と同じ制約と依存で足す（`cli.rs` の拒否の一覧から外す）。`rmblastn`（行列の得点と masklevel）は明示的な拒否のまま。
3. **fixture**：NCBI BLAST+ 2.17.0 の出力を oracle として固定する（`LOSAT/tests/outfmt0_manifest.tsv` と `LOSAT/tests/blastn_regression_fixtures.py` の方式。outfmt 0/6/7、スレッド 1/2/4）。dc-megablast は template の種類（coding・optimal・両方）× 長さ（16・18・21）× word size（11・12）、blastn-short は 50 塩基より短い query（primer の長さ）と既定値の上書き（`-word_size`・`-reward`・`-penalty`・`-evalue`・`-dust`）、両方の鎖、複数の query の batch、曖昧な文字、`-lcase_masking`、`-dust`、遠い種どうしの組（`LOSAT/tests/fasta` の genome の窓）、ヒット無し。変更前の実行ファイルでも `check` して、fixture が分岐を区別することを記録する。
4. **sweep**：E2c の `scoring_sweep.py`・`word_size_sweep.py`・`check_inputs.py`、E2f の batch の sweep を `-task dc-megablast` と `-task blastn-short` に広げ、NCBI と同じ拒否、outfmt 0/6/7 のバイト一致、明示的な拒否のどれかになることを確かめる。
5. **アダプタ**：ABI v2 の `describe` に task と新しい option を出し、`validate` が CLI と同じ文言になることを確かめる（`docs/web/abi_v2.md`）。`docs/web/verification_cells.tsv` に BLASTN の dc-megablast と blastn-short の行（native CLI と Wasm ABI v2 serial / threaded、スレッド 1/2/4）を足して埋める。
6. **ゲート**：BLASTN の既存のゲート（Gate A、S02 の capture、E2c・E2f・E2g の sweep と fixture、CI の速い検査）に退行なし、v1 の WASI の行列、V-ABI quick と full、`cargo fmt --check`・clippy・`cargo test --all-features`。V-PERF は megablast と blastn の case が非退行であること（dc-megablast は新しい経路なので、NCBI との時間の比を記録する）。
7. **独立監査**：観点ごとに sonnet の agent（読み取り専用）で、経路の網羅、移植の忠実さ、拒否の理由、ABI を確かめる（S08 の `~/.cache/losat-web-gui-target/s08-audit/` の指示の形）。4 つとも supported になるまで続ける。

完了条件は計画 §7 の SD の行による。記録は `docs/evidence/losat_web_e2i/README.md`（ゲート記録、`AUTHORITY.md`、`INVENTORY.tsv`、`evidence.sha256`、run の directory）。

## 進め方

- 作業場所：worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui`、ブランチ `feature/losat-web-gui`。S08b の PR を作った後に始める。`RUSTUP_TOOLCHAIN=1.92.0`、`--target-dir` は `~/.cache/losat-web-gui-target/sd-*`。NCBI のソース `/mnt/c/Users/genom/GitHub/ncbi-blast/c++`（598d8ae6）が唯一の正本、`/home/kawato/micromamba/bin` の NCBI BLAST+ 2.17.0 は比較だけに使う。
- アプリ側（S09、worktree `LOSAT-web-gui-app`）と並行する（計画 DW-7）。V-PERF の間は `vperf.lock` を置く（`~/.cache/losat-web-gui-target/s07p-resume/vperf_lock.sh`）。
- 機械的な作業（NCBI との比較、sweep、fixture の凍結、記録の下書き、監査）は Agent の `model: "sonnet"` に回し、互いに独立なら 4 つを超えて並行してよい。エンジンの変更、ゲート、V-PERF は 1 つずつ。agent には結果を途中でファイルに書かせ、作業ディレクトリと `--target-dir` を分ける。V-PERF の測り直しは `--repeat 5`。長い実行は `setsid nohup … &` で始め、30 分以内の間隔で見る。
- NCBI の不具合に当たる挙動は `PD-LOSAT-NCBI-DEFECTS` の規則で扱い、まとめて保守者に諮る。ほかは推奨の案で進めて記録する。
- ゲートの script は S08 のもの（`~/.cache/losat-web-gui-target/s08/s08_gates.sh`。最初に Gate A の字句のパスを `ci_fast_regressions.py` の `stage_lexical_fixtures` で用意する）を写して使う。

## S08b からの引き継ぎ（2026-10-03 の実測）

- **状態**：S08 と S08b は完了した（[E2b のゲート記録](../evidence/losat_web_e2b/README.md)）。エンジンの最後のコミットは `24fcfe41b`（native `2a46c2e4…`）。`main` への PR は S08b が作った（merge は保守者）。SD は PR の merge を待たずに同じブランチで始めてよい（S07+++ の後と同じ）。V-PERF の `tblastx-multi`（NCBI の batch の費用）は保守者の確認待ち。
- **変更前の基準**：S08b のゲートの成果物（`~/.cache/losat-web-gui-target/s08-gate-native/release/LOSAT`、`s08-gate-native-serial`、`s08-gate-wasi-artifacts`、`s08-gate-reactors`。ハッシュは `docs/evidence/losat_web_e2b/run-20261003T021918Z/artifacts.sha256`）。SD を始める前に `~/.cache/losat-web-gui-target/sd/bin/` に写し、capture（`docs/evidence/losat_web_e1a/capture_outputs.py`）を取る。
- **ゲートの script**：`docs/evidence/losat_web_e2b/gates/`（`s08_gates.sh`、`s08_gate_a.sh`、`s08_gate_a_after_build.sh`、`s08_perf.sh`、`s08_perf_rerun.sh`、`perf_precheck.py`、`verify_added.py`）。ゲートは約 3 時間（lint と試験 10 分、fixture 40 分、capture 35 分、V-ABI full 2.5 時間）、Gate A は約 3.5 時間で、別に並行して走らせてよい（native のビルドの後）。
- **環境の注意**：WSL の再起動で `/tmp` が空になる（Gate A の字句のパスは script が用意する）。Claude Code が再起動すると background のコマンドが止まった扱いになるが、実際には動き続けることがある：`ps` で確かめ、重ねて起動しない。長い実行は `setsid nohup … &` で起動し、待つのは log の until-loop（2 時間の上限があるので、切れたら起動し直す）。`pkill -f <パターン>` は自分の shell の command line にも一致するので、PID で止める。
- **監査**：S08 の 4 観点の指示（`~/.cache/losat-web-gui-target/s08-audit/COMMON.md`、`ANGLE_{A,B,C,D}.md`、`ROUND{2,3}.md`、写しは `docs/evidence/losat_web_e2b/audit/`）を BLASTN の task に直して使う。agent には監査する実行ファイルの写しを渡す（ゲートが同じ path を作り直すので）。agent が止まったら SendMessage で同じ文脈のまま再開できる。
- **アプリ側**：S09 が worktree `LOSAT-web-gui-app` で並行している。その directory（`app-s09-*`）には触れない。V-PERF の間は lock を置く（S09 は lock の間、試験とビルドを止める）。

## 終了・引き継ぎ

README の規則 8 に従う。`main` への PR を作り（merge は保守者）、CI が緑であることを確かめる。次は [S08+ — BLASTP・TBLASTN・TBLASTX の既定以外のオプション](session_s08p_e2e_protein_options.md)。
