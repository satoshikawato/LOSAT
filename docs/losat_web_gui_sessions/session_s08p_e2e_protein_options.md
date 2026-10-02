# Session S08+ — E2e：BLASTP・TBLASTN・TBLASTX の既定以外のオプション

## INSTRUCTION PROMPT

LOSAT の段階 E2e を実行する。BLASTP・TBLASTN・TBLASTX の検索のオプションのうち、認証済みの既定の値以外を、NCBI と同じにするか、明示的に拒否するセッションである（TD-13）。S07+（E2c）が BLASTN に行ったことを、ほかの 3 つの program に行う。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S08+ の行である。

現状：認証済みの fixture は、既定の値だけを使う（例：`LOSAT/tests/blastp_parity_options.sh` は BLOSUM62、gap 11/1、word size 3、threshold 11、window 40、`-comp_based_stats 2`、`-seg no` を明示するだけ）。例外は、承認済みの subject の遺伝暗号（`AGENTS.md`）だけである。一方、CLI と ABI v2 の `describe` は次を受け付ける（`LOSAT/src/algorithm/*/args.rs`）：

- BLASTP：`-matrix`、`-gapopen`、`-gapextend`、`-threshold`、`-word_size`、`-window_size`、`-comp_based_stats`、`-seg`、`-ungapped`、`-use_sw_tback`、`-evalue`、`-max_hsps`、`-max_target_seqs`
- TBLASTN：上の得点のオプションに加えて `-db_gencode`、`-max_intron_length`、`-xdrop_gap`、`-xdrop_gap_final`、`-sum_stats`、`-lcase_masking`、`-soft_masking`
- TBLASTX：`-threshold`、`-word_size`、`-window_size`、`-seg`、`-query_gencode`、`-db_gencode`、`-culling_limit`、`-evalue`、`-max_target_seqs`

S07+ で BLASTN に見つかった種類の差が、これらにもあるかは確かめていない（BLASTP は、BLOSUM62 と gap 11/1 以外の行列と gap、`-comp_based_stats 2` 以外を、既に明示的に拒否する。`LOSAT/src/algorithm/blastp/blast_engine.rs`）。S07+ の差は次のとおり（`docs/evidence/losat_web_e2c/AUTHORITY.md`）：

- NCBI が拒否する値を実行する、task の既定値の上書きの誤り、表に無い組の黙った代用、16 ビットなどの値の型。
- 引数の読み方：NCBI の整数の引数は `0x` の 16 進数も読み（`value_parsers.rs` の `ncbi_integer`）、実数の引数は最初の文字の規則と `strtod` で読む（`ncbi_double`）。BLASTN だけがこれらを使う。ほかの program の引数（`positive_usize` などの共有の parser）も同じにする。
- 順序：`-outfmt` の解析、subject のファイル（レコードが無ければ engine error）、query のファイル、`-out`、警告、検査、「Query is Empty!」、LOSAT の上限の順（BLASTN の `run`）。開けないファイルは NCBI の文言（`cli.rs` の `inaccessible`）。`-` は標準入力・標準出力、`-query` の既定値は `-`。
- `-outfmt` の NCBI の文言と終了コード（`blastn/hsp.rs` の `parse_blastn_output_format`）。
- NCBI の program に無いオプションは削る（`AGENTS.md` の規則 5。BLASTN では `-verbose` など 4 つ）。
- 制約のある整数の引数は、制約が `NStr::StringToDouble` で読み直すので `0x` だけを拒否する（`ncbi_constrained_integer`）。`-dust`・`-seg` などの文字列の引数は、NCBI の区切り方と誤り（`parse_dust_filtering`）。`-out` は 1 度だけ開き、名前は 256 バイト未満。予備の hit list の大きさなどの `Int4` の計算は NCBI と同じく折り返す（`get_prelim_hitlist_size`）。
- outfmt 0 の subject の title：NCBI は `CDeflineGenerator::GenerateDefline` で作る（BLASTN は `report/defline.rs` の `ncbi_nucleotide_title` で移植済み。S07+ の第 10 回の監査）。TBLASTN（核酸の subject）はこれを使い、BLASTP・BLASTX（蛋白の subject）は `x_CleanAndCompress` の蛋白の規則と `x_AdjustProteinTitleSuffix` を移植するか、明示的に拒否する。
- 引数の誤り（NCBI は USAGE と終了コード 1、LOSAT は clap の終了コード 2）の扱いを決める。
- S07+ の第 12 回の監査（`AUTHORITY.md` §O）：NCBI は得点 0 の HSP を Seq-align にしない（`blast_seqalign.cpp:672-674`、全 program）。subject の小文字の区間（`SetupSubjects_OMF`、全 program）は、前の区間の最後の文字から走査の区間を始める。どちらも、ほかの program に当てはまる場合（例えば TBLASTN の subject の曖昧な文字と得点 0、`-lcase_masking` の subject）を確かめ、NCBI と同じにするか、明示的に拒否する。
- 全 program に共通の CLI の誤りの経路：出力の書き込みの失敗（NCBI は「BLAST failed to write output」と終了コード 6、表形式は abort。BLASTX は既にそうしている）と、UTF-8 でないファイル名の表示（`cli.rs` の `inaccessible`）と引数の値（`-dust`・`-outfmt` など）を、BLASTN を含めて NCBI と同じにする（S07+ の第 5・6 回の監査、`docs/evidence/losat_web_e2c/AUTHORITY.md` §I の末尾）。

これらのオプションは、アプリの検索画面（S12）に出る。

0. **棚卸しの方式（計画 DW-12）。** 監査の指摘を 1 つずつ直すのでなく、S07+++ と同じく、program ごとに NCBI の実行経路に現れる関数を棚卸しし（`docs/evidence/losat_web_e2e/INVENTORY.tsv`）、未移植と差のある移植を先に一括で transpile してから、sweep と監査で確かめる。LOSAT が速度のために NCBI と違う実装にしている箇所は、出力が同じなら新しく移植する部分にも使ってよい。
1. **範囲を決める。** program ごとに、受け付けるオプションと、NCBI の引数の制約（`c++/src/algo/blast/blastinput/blast_args.cpp` など）、検査（`c++/src/algo/blast/core/blast_options.c` の `BLAST_ValidateOptions`）、task の既定値（`c++/src/algo/blast/api/blast_prot_options.cpp`、`blast_advprot_options.cpp` など）、表（`c++/src/algo/blast/core/blast_stat.c` の行列ごとの gap の表、`Blast_KarlinBlkGappedLoadFromTables`）の対応を記録する（`docs/evidence/losat_web_e2e/AUTHORITY.md`）。
2. **sweep を作る。** `docs/evidence/losat_web_e2c/scoring_sweep.py` の形で、program ごとに、行列 × gap の組（NCBI の表にあるもの、無いもの、境界）、threshold と word size、`-comp_based_stats` の値、`-seg` の値、遺伝暗号（承認済みの例外の扱いは `AGENTS.md`）、e-value の書き方を、outfmt 0/6/7 で NCBI と比べる。比べる前に、今の commit の結果を記録する。
3. **NCBI の検査を移植する。** NCBI が拒否する値を、NCBI と同じ順序、文言、終了コードで拒否する（S07+ の `LOSAT/src/algorithm/blastn/scoring.rs` と `LOSAT/src/cli.rs` の `NativeError` を参考にし、共有できる部品は共有する）。
4. **違いの原因を調べて直す。** NCBI が受け付けて結果が違う組合せは、最初に値が食い違う箇所（lambda、K、H、生の得点、X-drop、cutoff、`eff_searchsp`）をソースで追って記録してから直す。直す箇所の直上に NCBI のファイル・行と断片を書く。直した組合せは fixture にし、`docs/evidence/losat_web_e2a/run_oracle.py` で固定する。
5. 直せなかった値は、LOSAT が対応していないことを示す文言（`not supported by LOSAT's <PROGRAM>`）で明示的に拒否する。黙って NCBI と違う結果を出さない。拒否がアダプタの `validate` にも出ることを確かめる。
6. BLASTX は範囲に入れない（DW-10）。BLASTX の同じ確認は SX で行う（[SX の指示書](session_sx_blastx_integration.md) に書き足す）。
7. 試験：sweep の全組合せが、NCBI と同じ拒否、outfmt 0/6/7 のバイト一致、明示的な拒否のどれかになる。各 program の既存のゲート（v0.1.0・v0.2.0 の manifest、Gate A、TLOSAN の Stage G）と S07・S08 の fixture が変わらない。エンジンの変更の V-PERF の非退行。変えた program の全升目の V-ABI。`cargo fmt --check`・`clippy`・`cargo test --all-features`。独立監査。

記録は `docs/evidence/losat_web_e2e/`（`README.md`、`evidence.sha256`、変更前と変更後の sweep の結果）。変更前の成果物は S08 の成果物である。

## S08 からの引き継ぎ（2026-10-02・03 の実測）

S08（E2b、[ゲート記録](../evidence/losat_web_e2b/README.md)、[権威の記録](../evidence/losat_web_e2b/AUTHORITY.md)、棚卸し `docs/evidence/losat_web_e2b/INVENTORY.tsv`）が残したもの。どれも S08 で明示的に拒否しているか、S08 の前からの差で、S08 では直していない。S08 の独立監査の作業ディレクトリ（`~/.cache/losat-web-gui-target/s08-audit/{a,b,c,d}/`、比べた全件の入力と出力）に再現がある。

- **最初の作業：TBLASTN の subject の定義行を黙って読み違える**（S08 の監査 (c) の F14、中程度、S08 の前から）。TBLASTN は subject を `bio` だけで読み（`algorithm/tblastn/args.rs`）、tab や制御文字のある定義行の outfmt 0 の title（NCBI は最初の 0x20 未満の byte で切る）、空の定義行や空白で始まる定義行の outfmt 6 の `sseqid`（NCBI は `Subject_2` など）が NCBI と違う。TBLASTX と同じく `blastn/input.rs` の部品（`check_deflines_of` ほか、program の名前つき）で読むか、NCBI の読み方を移植する。

- **TBLASTX の入力の読み方は BLASTN と同じになった**（`LOSAT/src/algorithm/blastn/input.rs` の `open_input`・`read_records`・`parse_fasta`・`check_deflines_of`・`check_residues_of`・`check_utf8_file_name`。program の名前を受け取る）。BLASTP・TBLASTN は今も `bio` で読むだけなので、同じ部品を使って、定義行・残基・空のレコード・UTF-8 でないファイル名・`-out` を開く順序・`U`・CFastaReader の title の警告を NCBI と同じにするか拒否する（TBLASTN は核酸の subject、BLASTP は蛋白なので、蛋白の残基の規則は `fasta.cpp` から別に移植する）。
- **TBLASTX の引数**：`-threshold`（NCBI は実数。LOSAT は整数だけ。outfmt 0 の epilog は S08 で `%g` に直した）、`-evalue`（NCBI の `CArg_Double` は `nan` を USAGE と終了コード 1 で拒否し、`+inf` などを読む。LOSAT は有限の数だけ）、整数の引数の `0x`（`-max_target_seqs 0x10` を NCBI は読む）、`-outfmt` の `+6`・`06`・`7 std` などの文言、`-subject` が無いときの NCBI の誤り、`-query -`（標準入力）、`-word_size` 2・4（LOSAT は 3 だけを受け付ける）、`-window_size 0`（一つの hit の word finder。S08 で明示的に拒否）、`-num_descriptions`・`-num_alignments`・`-sorthits`・`-sorthsps`・`-sum_stats`・`-line_length`（LOSAT の TBLASTX に無い）。
- **TBLASTX の `-culling_limit`**：S08 で明示的に拒否した（LOSAT の culling は NCBI と違う HSP を残す。例：`-culling_limit 1` で NCBI の 18 行に対し LOSAT は 0 行。`docs/evidence/losat_web_e2b/inventory/result_G_notes.md`・`result_E_notes.md`）。移植するなら `hspfilter_culling.c` と、culling のときに hit list の切り詰めが traceback の段に移ること（範囲 G の行 32）。
- **`-seg` の値**：S08 の監査の後（`7fbbfad96`）、NCBI と同じに 1 つの空白で 3 つに分け、窓は `int`、0 以下の窓・locut・hicut は NCBI の既定（12、2.2、2.5）のままにする（`blast_filter.c:1147-1154`、`value_parsers.rs` の `SegSpec::params`。BLASTP・TBLASTN・TBLASTX 共有）。outfmt 6 は 3 つの program の 27 の組で NCBI と一致した（ゲート記録）。残り：誤りの文言と終了コード（NCBI は `BLAST query/options error: Invalid number of arguments to filtering option` または `Invalid input for filtering parameters`、終了コード 1。LOSAT は clap の誤りで 2。BLASTN の `-dust` は E2g で NCBI と同じにした）、`nan`・`inf` の locut・hicut（LOSAT は拒否）、NCBI の `StringToDouble` が読む書き方。
- **引数の解析の後の NCBI の検査の文言**（監査 (a) の F5）：`-threshold 0`（`Non-zero threshold required`）、`-word_size`、壊れた `-seg` は、NCBI が `BLAST query/options error: …` と `Please refer to the BLAST+ user manual.`、終了コード 1。LOSAT は引数の解析器で止め、clap の文言と終了コード 2（承認済みの例外 1 は構文の誤りだけを対象にする）。`-evalue 0` は S08 で NCBI と同じにした（`check_hit_saving_options`）。負の `-evalue`（NCBI は同じ誤り、LOSAT は解析器で 2）も直す。
- **TBLASTX の環境と設定**：`check_ncbi_application_settings`（`LOSAT/src/blastinput/ncbi_environment.rs`：`DIAG_*`、`NCBI_CONFIG_*`、出力を変える `.ncbirc`、NCBI が値で止まる `BLAST_USAGE_REPORT`・`DEBUG_CATCH_UNHANDLED_EXCEPTIONS`・`NCBI_STATIC_ARRAY_COPY_WARNING`・`NCBI_USAGE_REPORT_*`・`OBJMGR_BLOB_CACHE`・`OBJMGR_SCOPE_*` など）は BLASTN だけが呼ぶ（監査 (c) の F13：`DIAG_POST_LEVEL=Error` で NCBI は `-max_target_seqs 1` の警告を出さない、LOSAT は出す。監査 (a) 第 2 回：NCBI が `getenv` で読む 148 の名前の一覧は `~/.cache/losat-web-gui-target/s08-audit/r2a/envtrace/all_names.txt`）。`BL2SEQ_LEGACY`・`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`・`PRE_FETCH_SEQS_LIMIT` は S08 で NCBI と同じか拒否にした。
- **TBLASTX の入力の読み方の拒否のうち、安く移せるもの**（監査 (c) の F4・F5、監査 (a) の F7）：`>` の後の空白を飛ばす、ID は最初の 0x20 以下の byte まで、空の ID は `Query_N`（subject は `Subject_N`）、最初の定義行の前の空行と `;` の注釈、配列の行の中の空白と `;` 以降（NCBI は黙って飛ばす）、`-`・`*`・`X` などの無効な残基（NCBI は `FASTA-Reader: Ignoring invalid residues at position(s): On line N: M` の警告で落とす）、中身の無いレコード（NCBI は警告して続けるか、query だけなら終了コード 3。subject の `>` だけのレコードは、NCBI が `Subject sequence contains no data` と `The average subject length is too short`、prolog の後に終了コード 3。LOSAT は `Empty CBlastQueryVector`、prolog 無しで終了コード 3。監査 (a) 第 2 回の N4）、最初の定義行の前の残基の文字列（NCBI は定義行の無いレコードとして読む。`CheckDataLine` で止まる行は `CFastaReader: Near line N, there's a line that doesn't look like plausible data…`、終了コード 1。LOSAT は S08 で読むときに拒否する）。BLASTN と共有の `blastn/input.rs` なので、BLASTN の S07+ の記録（`docs/evidence/losat_web_e2c/AUTHORITY.md` §E・§M・§N）と合わせて直す。
- **標準入力**：`-query -`・`-subject -` は clap の誤り（文言が LOSAT を名指さない）。空のパイプの query は NCBI が空とみなさず、query の無い報告を書く（outfmt 0 は 615 バイト、outfmt 7 は `# BLAST processed 0 queries`）。LOSAT は明示的に拒否（監査 (c) の F11）。
- **`-num_threads` 65535 以上**：全 program の CLI が Rayon の上限で拒否する。NCBI は CPU の数に切り詰め、`-subject` では無視する。
- **BLASTN の outfmt 0 の題の HTML の拒否**：S08 で NCBI の `HtmlDecode` と同じ判定にした（BLASTN も共有）が、BLASTN は今もヒットの無い subject も含めて全 subject を確かめる（TBLASTX・TBLASTN は報告に出る subject だけ）。BLASTN の `search` の検査を、TBLASTX の `check_shown_subject_titles` のように最終の hit list の subject に狭める。
- **閉じた標準出力**：起動の時に閉じた標準出力（`>&-`）で、NCBI は最初の書き込みで失敗する（outfmt 0 は終了コード 6、6/7 は abort）。Rust の runtime は main の前にそこへ `/dev/null` を開くので、LOSAT は呼び出し側が読み書きで開いた `/dev/null`（Python の `subprocess.DEVNULL`）と区別できず、書いて成功する。S08 の `7fbbfad96` の見分け方は後者を誤って失敗させたので外した（`1117e8c17`）。保守者に承認済みの例外にするか諮った（S08 のゲート記録）。
- **BLASTP の S08 の監査 (b) で見つかった差**（S08 の前から）：`-seg yes` で query 全体が mask されると NCBI は Karlin の警告を stderr に出し、LOSAT の BLASTP は出さない（TBLASTX・TBLASTN は出す。703 件の stderr だけの差）。`blastp -evalue 1000` で LOSAT が 1 残基の HSP を 1 つ多く報告する（`~/.cache/losat-web-gui-target/s08-audit/b/` の `P2.blastp.3383.f6`）。
- **BLASTP・TBLASTN の `BLAST_Cutoffs` の下限**：TBLASTX は S08 で、呼び出し側の 1 と E からの得点の大きい方にした（`ncbi_cutoffs.rs` の `blast_cutoffs_from_one`）。TBLASTN（`stage_d_stats.rs`）と BLASTX は自分の呼び出しで `.max(1)` を付けている。BLASTP の呼び出しも確かめる。
- **拒否の前の余計な題の警告**（監査 (b) の F-7、情報）：最後に 20 以上の塩基と空白がある subject の定義行で、LOSAT は CFastaReader の題の警告を出してから拒否する（NCBI はその空白のため警告しない）。入力の読み方を移すときに直す。
- **BLASTP の説明の一覧**：BLASTP（と BLASTX）は今も最初の HSP と素の幅の `write_subject_summary_table_with_sum_n` を使う（`report/pairwise.rs`）。NCBI の規則は E2a §G.3 の `x_InitDeflineTable`（最も高い bit score とその e-value、最後の行の合計の幅）。BLASTN・TBLASTN・TBLASTX は `write_blastn_description_table` を使う。蛋白の subject の title（`x_CleanAndCompress` の蛋白の規則、`x_AdjustProteinTitleSuffix`）と合わせて移すか拒否する。
- **BLASTP・TBLASTN の epilog の threshold**：TBLASTX は S08 で C++ の stream の既定の書き方（`%g`、精度 6。1000000 は `1e+06`）にした。BLASTP・TBLASTN の `write_*_epilog` は今も整数のまま書く。
- **TBLASTN の表示の翻訳**：TBLASTN の Sbjct の行は検索の翻訳（曖昧な codon は `X`）から作る（`algorithm/tblastn/stage_e_report.rs`）。NCBI の表示は `CTrans_table` で B・Z・J を出す（TBLASTX は S08 で移植、`algorithm/tblastx/report.rs` の `display_residue`）。
- **BLASTP の同点**：BLASTP の `score_compare_match` の移植（`algorithm/blastp/blast_engine.rs`）は `sort_unstable_by` で、glibc の `qsort`（安定な merge sort）と違い、同点の順が変わりうる。TBLASTX は S08 で `Blast_InitHitListSortByScore` の安定な並べ替えを足して、同点の HSP の連結の差を直した（`algorithm/tblastx/blast_gapalign.rs` の `sort_init_hsps_by_score_ncbi`）。
- **BLASTP・TBLASTN の outfmt 0 の S08 で見つかった差**（S08 の SEG の調べ、`docs/evidence/losat_web_e2b/investigations/segmask.md`）：BLASTP の `-seg yes` の outfmt 0 は SEG で mask した query の残基を小文字にしない（NCBI はする。TBLASTN・TBLASTX・BLASTX はする）。BLASTP の outfmt 0 の `Method:` の名前が一部の alignment で違う（NCBI `Composition-based stats`、LOSAT `Compositional matrix adjust`）。TBLASTN の outfmt 0 の Lambda・K・H の行が、乱数の 3000 の入力のうち約 35 で最後の桁が違う。どれも S08 の前からの差。
- **SEG の左の再帰**：S08 で NCBI と同じにした（`utils/seg.rs`。BLASTP・TBLASTN・TBLASTX が共有）。fixture `seg.blastp.6`（`LOSAT/tests/outfmt0_manifest.tsv`）。BLASTP・TBLASTN の sweep でも確かめる。
- **outfmt 0 の punctuation だけの subject の title**（NCBI の `x_CleanAndCompress` が文字列の外を読み、tblastx・tblastn も落ちる）：TBLASTX と TBLASTN は S08 で明示的に拒否した。保守者の判断（BLASTN の例外 2 を TBLASTX・TBLASTN に広げるか）の結果に合わせる（S08 のゲート記録の「保守者への質問」）。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S09 — ブラウザでの実行基盤](session_s09_w1_browser_runtime.md)。対応するオプションの値（直したものと拒否するもの）を [S12 の指示書](session_s12_w3_search_ui.md) に書き足す。
