# Session S08+ — E2e：BLASTP・TBLASTN・TBLASTX の既定以外のオプション

## INSTRUCTION PROMPT

LOSAT の段階 E2e を実行する。BLASTP・TBLASTN・TBLASTX の検索のオプションのうち、認証済みの既定の値以外を、NCBI と同じにするか、明示的に拒否するセッションである（TD-13）。S07+（E2c）が BLASTN に行ったことを、ほかの 3 つの program に行う。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S08+ の行である。

現状：認証済みの fixture は、既定の値だけを使う（例：`LOSAT/tests/blastp_parity_options.sh` は BLOSUM62、gap 11/1、word size 3、threshold 11、window 40、`-comp_based_stats 2`、`-seg no` を明示するだけ）。例外は、承認済みの subject の遺伝暗号（`AGENTS.md`）だけである。一方、CLI と ABI v2 の `describe` は次を受け付ける（`LOSAT/src/algorithm/*/args.rs`）：

- BLASTP：`-matrix`、`-gapopen`、`-gapextend`、`-threshold`、`-word_size`、`-window_size`、`-comp_based_stats`、`-seg`、`-ungapped`、`-use_sw_tback`、`-evalue`、`-max_hsps`、`-max_target_seqs`
- TBLASTN：上の得点のオプションに加えて `-db_gencode`、`-max_intron_length`、`-xdrop_gap`、`-xdrop_gap_final`、`-sum_stats`、`-lcase_masking`、`-soft_masking`
- TBLASTX：`-threshold`、`-word_size`、`-window_size`、`-seg`、`-query_gencode`、`-db_gencode`、`-culling_limit`、`-evalue`、`-max_target_seqs`

S07+ で BLASTN に見つかった種類の差（NCBI が拒否する値を実行する、task の既定値の上書きの誤り、表に無い組の黙った代用、NCBI の検査の順序、16 ビットなどの値の型、e-value の書き方）が、これらにもあるかは確かめていない。これらのオプションは、アプリの検索画面（S12）に出る。

1. **範囲を決める。** program ごとに、受け付けるオプションと、NCBI の引数の制約（`c++/src/algo/blast/blastinput/blast_args.cpp` など）、検査（`c++/src/algo/blast/core/blast_options.c` の `BLAST_ValidateOptions`）、task の既定値（`c++/src/algo/blast/api/blast_prot_options.cpp`、`blast_advprot_options.cpp` など）、表（`c++/src/algo/blast/core/blast_stat.c` の行列ごとの gap の表、`Blast_KarlinBlkGappedLoadFromTables`）の対応を記録する（`docs/evidence/losat_web_e2e/AUTHORITY.md`）。
2. **sweep を作る。** `docs/evidence/losat_web_e2c/scoring_sweep.py` の形で、program ごとに、行列 × gap の組（NCBI の表にあるもの、無いもの、境界）、threshold と word size、`-comp_based_stats` の値、`-seg` の値、遺伝暗号（承認済みの例外の扱いは `AGENTS.md`）、e-value の書き方を、outfmt 0/6/7 で NCBI と比べる。比べる前に、今の commit の結果を記録する。
3. **NCBI の検査を移植する。** NCBI が拒否する値を、NCBI と同じ順序、文言、終了コードで拒否する（S07+ の `LOSAT/src/algorithm/blastn/scoring.rs` と `LOSAT/src/cli.rs` の `NativeError` を参考にし、共有できる部品は共有する）。
4. **違いの原因を調べて直す。** NCBI が受け付けて結果が違う組合せは、最初に値が食い違う箇所（lambda、K、H、生の得点、X-drop、cutoff、`eff_searchsp`）をソースで追って記録してから直す。直す箇所の直上に NCBI のファイル・行と断片を書く。直した組合せは fixture にし、`docs/evidence/losat_web_e2a/run_oracle.py` で固定する。
5. 直せなかった値は、LOSAT が対応していないことを示す文言（`not supported by LOSAT's <PROGRAM>`）で明示的に拒否する。黙って NCBI と違う結果を出さない。拒否がアダプタの `validate` にも出ることを確かめる。
6. BLASTX は範囲に入れない（DW-10）。BLASTX の同じ確認は SX で行う（[SX の指示書](session_sx_blastx_integration.md) に書き足す）。
7. 試験：sweep の全組合せが、NCBI と同じ拒否、outfmt 0/6/7 のバイト一致、明示的な拒否のどれかになる。各 program の既存のゲート（v0.1.0・v0.2.0 の manifest、Gate A、TLOSAN の Stage G）と S07・S08 の fixture が変わらない。エンジンの変更の V-PERF の非退行。変えた program の全升目の V-ABI。`cargo fmt --check`・`clippy`・`cargo test --all-features`。独立監査。

記録は `docs/evidence/losat_web_e2e/`（`README.md`、`evidence.sha256`、変更前と変更後の sweep の結果）。変更前の成果物は S08 の成果物である。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S09 — ブラウザでの実行基盤](session_s09_w1_browser_runtime.md)。対応するオプションの値（直したものと拒否するもの）を [S12 の指示書](session_s12_w3_search_ui.md) に書き足す。
