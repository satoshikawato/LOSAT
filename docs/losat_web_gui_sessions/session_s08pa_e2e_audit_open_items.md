# Session S08+a — E2e：第 1 回の独立監査の残件（S08+ の仕上げと並行）

## INSTRUCTION PROMPT

LOSAT の段階 E2e（BLASTP・TBLASTN・TBLASTX の既定以外のオプション）の、第 1 回の独立監査で見つかり S08+ で直さなかった指摘を、NCBI とバイト一致にするか明示的に拒否する。S08+ のセッションは同じ時間に、本線の worktree で最後のゲート・V-PERF・記録を進めている（2026-10-04 に保守者の指示で分けた）。このセッションは別の worktree と別のブランチで作業し、本線の worktree・ゲート・V-PERF・記録に触れない。全ゲート、独立監査の第 2 回、V-PERF、S08+ の完了は、このセッションの後の S08+b（S08+ が最後に書く指示書）で、このブランチを本線に merge してから行う。先に [セッション README](README.md) の共通規則（特に規則 4・5）を読む。

先に読む：
- [第 1 回の独立監査の記録](../evidence/losat_web_e2e/audit/ROUND1.md)（指摘と対応の表）と観点ごとの報告 [`audit/round1/`](../evidence/losat_web_e2e/audit/round1/)（再現の命令）。
- [E2e の権威の記録](../evidence/losat_web_e2e/AUTHORITY.md) の §F（lookup と `(Int4)`、fence の行）、§M（判断 D1〜D12）、§N（残り）。

### 作業場所と約束

- worktree：`git -C /mnt/c/Users/genom/GitHub/LOSAT-web-gui worktree add /mnt/c/Users/genom/GitHub/LOSAT-web-gui-s08pa -b feature/losat-web-gui-s08pa d96412265`（`d96412265` は S08+ の監査の修正と記録の後のコミット。既に worktree があれば使う）。push は `origin feature/losat-web-gui-s08pa`。
- 本線の worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui` とブランチ `feature/losat-web-gui` には書かない（S08+ のゲートが走っている）。アプリ側の worktree `LOSAT-web-gui-app` にも触れない。
- 記録は新しい `docs/evidence/losat_web_e2e/s08pa/`（`NOTES.md` と、確かめの結果）だけに書く。`AUTHORITY.md`、`README.md`、`evidence.sha256`、`audit/`、セッションの指示書と README の表は S08+ が仕上げるので直さない（S08+b が merge の後に合わせる）。fixture の行は `LOSAT/tests/outfmt0_manifest.tsv` の末尾に足してよい（S08+ は以後この表を変えない）。
- `RUSTUP_TOOLCHAIN=1.92.0`。`--target-dir` は `~/.cache/losat-web-gui-target/s08pa-*`（S08+ の `s08p-*` を使わない）。NCBI のソースの写しは `~/.cache/losat-web-gui-target/s08p/ncbi/c++`（598d8ae6、読むだけ）、比較用の NCBI BLAST+ 2.17.0 は `/home/kawato/micromamba/bin`。`-remote` は使わない。
- **V-PERF の lock**：S08+ はゲートの後に `/home/kawato/.cache/losat-web-gui-target/vperf.lock` を置いて性能を測る。このファイルがある間は、ビルド・`cargo test`・NCBI の比較の実行を始めず、実行中のものは止めて、lock が消えるのを待つ（ソースを読むことと編集は続けてよい）。
- `/tmp` の要らないファイルは消してよいが、`/tmp/losat-pr5-runtime-cert-*`（Gate A の字句の入力）と使用中のものは消さない。
- agent（Sonnet）は読み取り専用の調べ物と比較の実行に使い、`/mnt/c` を読ませない（WSL の 9p の不安定）。必要な写しは `~/.cache/losat-web-gui-target/s08pa/` に置く。

### 1. 残件（優先の順）

各指摘は、まず NCBI のソースと trace で最初に値が食い違う箇所を突き止めて `s08pa/NOTES.md` に書き（NCBI のファイル・行と値）、それから直す。直す箇所の直上に NCBI のファイル・行と断片を書く。直せない値は、LOSAT が対応しないことを示す文言（`... is not supported by LOSAT's <PROGRAM>`）で明示的に拒否する（アダプタの `validate` にも出るように）。NCBI の値は TLOSAN Stage D の trace の shim（`docs/evidence/tlosan_stage_d/ncbi_d_call_trace.c` を `gcc -shared -fPIC ... -ldl` で作り `LD_PRELOAD`。`Blast_HSPListGetEvalues`・`BLAST_LinkHsps`・`Blast_HSPListReapByEvalue` の引数と HSP の一覧を出す）と、TLOSAN Stage E の C++ API の oracle（`~/.cache/losat-web-gui-target/s08p/api-oracle/tblastn_stage_e_local_oracle`、作り方は `docs/evidence/losat_web_e2e/gates/build_api_oracle.sh`）で取れる。

1. **TN-2**（高）：TBLASTN の `-comp_based_stats 0 -lcase_masking -seg no`（hard mask）と少ない `-max_target_seqs`（1〜11 で再現、500 では一致）で、和の統計の e-value が 1〜4 % 違う。`-sum_stats false` では一致。再現：`docs/evidence/losat_web_e2e/audit/round1/tblastn_inputs/gen_q_lc.faa` × `gen_s.fna`、`-comp_based_stats 0 -lcase_masking -seg no -max_target_seqs 1 -outfmt 6`（例：`gl0 gs20` の最初の行 NCBI `9.33e-177`、LOSAT `9.00e-177`）。S08+ で取った NCBI の trace（`~/.cache/losat-web-gui-target/s08p/tn1/n2.trace`）：`gs20` の traceback の `BLAST_LinkHsps` は subject の長さ 489、その中の `Blast_HSPListGetEvalues` は 163（489/3）、context 0 の `eff_searchsp` 5529312、length adjustment 57。LOSAT の `stage_d_pipeline.rs` の mode 0 の `link_preliminary_hsps`（`post_link`）に渡る長さ、Spouge の `db_length`、`lengths` を同じ点で比べ、`-max_target_seqs` に依る値（予備の段で残る subject から作る値など）を探す。
2. **TN-4**（高）：巨大な `-xdrop_gap`・`-xdrop_gap_final`（1e7〜5e8 bit）で LOSAT は数十秒・数 GB（`-xdrop_gap 5e8` で 20 GB、`-xdrop_gap_final 5e8` は強制終了）、NCBI は 1 秒未満。出力は同じ。原因の見込み：`LOSAT/src/algorithm/blastp/gapalign.rs` の `gap_dp_reserve_initial`・`gap_dp_reserve_band` が `Vec::resize` で全 cell を書く（NCBI の `blast_gapalign.c:797-808,923-931` は `malloc` で、触れない page は確保されない）。traceback の状態（`GapAlignScratch` の行）も確かめる。読む前に書く cell だけを使う形（初期化しない確保）にするか、それが NCBI の読み方と同じにできないなら、閾値で明示的に拒否する。BLASTP・TBLASTN（と BLASTX）が共有する経路なので、BLASTX の出力を変えない（DW-10）。速度の比べ方は `docs/evidence/losat_web_e2b/perf_cases.py` の BLASTP・TBLASTN の case を変更前（`d96412265` の成果物）と交互に（V-PERF の正式の測定は S08+b）。
3. **TN-5**（高）：BLOSUM45 の組（`-matrix BLOSUM45 -word_size 2 -comp_based_stats 0`）で `-evalue` 5000 以上のとき、同じ得点の HSP の subject の座標・frame が違う（`e2e_protein_query.faa` × `e2e_tblastn_subject.fna`、`-evalue 1e4 -outfmt 6` で 29900 行中 4 行、例 NCBI `4712 4695`、LOSAT `4710 4693`）。2000 以下と BLOSUM62 は一致。同点の並べ替え・heap・格納の順（`blast_hits.c`、`blast_gapalign.c` の `s_BlastGapAlignStruct...` の得点の比べ方）を NCBI の trace と比べる。
4. **RP-4**（前からの差、既定の option）：60000 残基の query の 4998 残基の一致で、TBLASTN の bit score が NCBI 10404（raw 26998）、LOSAT 10412（raw 27020）。subject を 1 つ（`bigs4.fna`）にすると NCBI 10383、LOSAT 10404（NCBI は他の subject に依る）。`-comp_based_stats 0` では一致、BLASTP の同じ蛋白は一致。入力は `docs/evidence/losat_web_e2e/audit/round1/rp4/`（gzip）。composition の行列の調整（`blast_kappa.c`、`composition_adjustment/`）の窓・組成・尺度を trace で比べる。TLOSAN の認証の範囲（Stage G）にあった入力でないことを確かめ、範囲の外なら記録する。
5. **`-out -version`**（低）：NCBI は `-version` を toolkit の version の option として読み version を出す。LOSAT は `-out` の値にして `-version` という名のファイルを作る。NCBI の `ncbiapp.cpp` の argv の前処理（`s_ArgVersion` を argv のどこでも探すか）を確かめ、同じにするか明示的に拒否する。

### 2. 試験と確かめ（このセッションの範囲）

- 直した組合せは fixture にし、`docs/evidence/losat_web_e2a/run_oracle.py --bin-dir /home/kawato/micromamba/bin --out LOSAT/tests/fixtures/outfmt0 --freeze --only <id>` で NCBI の出力を固定する。`docs/evidence/losat_web_e2a/check_losat.py --losat <新> --threads 1/2/4` で全 fixture の差 0、変更前（`d96412265`）の実行ファイルは新しい fixture で違うこと。
- TBLASTN と BLASTP の sweep（`docs/evidence/losat_web_e2e/option_sweep.py --program tblastn|blastp ... --api ~/.cache/losat-web-gui-target/s08p/api-oracle/tblastn_stage_e_local_oracle`、作業ディレクトリは `~/.cache/losat-web-gui-target/s08pa-sweep-*`）で DIFF 0、timeout 0。監査の再現の命令（`audit/round1/*.md`）を全て繰り返す。
- `cargo fmt --check`、clippy の 4 構成（`--all-targets --all-features`、`--all-targets --no-default-features`、`--lib --target wasm32-wasip1 --no-default-features`、`--lib --target wasm32-wasip1-threads --features wasm-threads`、`-D warnings`）と adapter、`LOSAT_BLASTX_WORKER_LOG=... cargo test --all-features`。このセッションで足した行の NCBI の参照は `~/.cache/losat-web-gui-target/s08p/verify_refs.py` と `docs/evidence/losat_web_e2e/gates/verify_added.py`（BASE は `78c06fe61`）で誤り 0。
- 全ゲート、V-ABI、v1 の WASI の行列、V-PERF の正式の測定、独立監査は行わない（S08+b）。

### 3. 終わり方

エンジンの変更は指摘ごとにコミットし、`origin feature/losat-web-gui-s08pa` に push する。`s08pa/NOTES.md` に、指摘ごとの原因（NCBI のファイル・行、値）、直し方、fixture、確かめの結果、残ったもの（直せず拒否にしたものと理由）を書く。保守者の判断が要るもの（NCBI の不具合に当たる挙動：`PD-LOSAT-NCBI-DEFECTS` の扱い、推奨の案）は、推奨の案で進めて記録し、まとめて諮る項目として NOTES に書く。最終回答には、コミット SHA、push の結果、確かめの結果、残件を示す。S08+ の本線への merge は S08+b が行う（このセッションは本線のブランチを変えない）。
