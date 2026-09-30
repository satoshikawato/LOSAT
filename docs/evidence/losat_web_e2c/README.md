# LOSAT Web E2c（Session S07+）ゲート記録

- 段階：E2c BLASTN の得点のオプションと入力の読み方（[総合計画書](../../losat_web_gui_plan.md) §7 の S07+、[指示書](../../losat_web_gui_sessions/session_s07p_e2c_blastn_scoring_options.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S07 の記録の `1bf597c1a`（エンジンは S07 の `4404a082c`）
- エンジン：`60ebd226b`（移植）、独立監査の各回の修正 `5377fc1bd`・`2fd0ee2c1`・`c1e933757`・`637426f30`・`7d6fee2a7`・`c33731d08`・`010bee067`・`6b176da78`・`d46906312`・`f8e4b0426`・`9adac5b1e`・`07e8c4fa8`・`4cab5c3fb`・`181c18ac1`・`bbc9c7809`、v1 の WASI の検査の文言の追従 `2b97f22f9`
- アダプタ：`5ea3f9767`、`0c34edebe`、`c27fbf93b`、`295cf0b40`、`0c7985789`、`6c8826068`、`722e6b95d`、`a9e17ca43`、`e4409ccfa`
- 権威の記録：[`AUTHORITY.md`](AUTHORITY.md)（NCBI の経路の表 §A〜§F、LOSAT が明示的に拒否するもの §G、独立監査で見つかったもの §H〜§R）
- 実行記録：[`run-20260930T022908Z/`](run-20260930T022908Z/)（`head.txt` がゲートを実行した commit）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)（リポジトリのルートで `sha256sum --check docs/evidence/losat_web_e2c/evidence.sha256` を実行する）
- 判定：**完了条件を満たした**。ゲートは `e8de6861b`（`head.txt`）で実行した。独立監査は 16 回目で supported。V-PERF の最初の計測で 3 件が +5% を超えたが、10 回の交互の計測ではすべて閾値の内だった（下の「性能の計測」）

## 変更の内容

### エンジン

| 範囲 | 内容（NCBI の経路は `AUTHORITY.md`） |
|---|---|
| 得点のオプション | `-reward`・`-penalty`・`-gapopen`・`-gapextend`・`-word_size` は省略できるキーで、与えた値は megablast の既定値と同じでも残る（S06 §D.5 の reward 1 の差の原因、§A）。reward と penalty は NCBI と同じく 16 ビットに折り返す。NCBI のオプションの検査を、NCBI の順・文言・終了コードで移植した（§B）。Karlin-Altschul の表（`s_GetNuclValuesArray`、最大公約数、線形の行、表を超える gap、`round_down`、表に無い組の文言）を移植した（§C） |
| query の context ごとの block | 鎖ごとの組成から ungapped と gapped の block を求め、cutoff・X-drop・gap trigger・探索空間・e-value・bit score に使う。無効な context は NCBI と同じく検索しない。subject ごとの HSP を e-value で並べ直す（§D） |
| 引数の読み方と順序 | 整数は NCBI の `CArg_Integer`（10 進数、`0x` の 16 進数）、実数は NCBI の最初の文字の規則。`-outfmt` の解析、subject を読む、query と `-out` を開く、`-dust`、警告、検査、「Query is Empty!」、LOSAT の上限の順を NCBI の引数の処理に合わせた。`-query` の既定値 `-`、`-subject -`、`-out -`、名前付きパイプ、ディレクトリ、`-out` の名前の長さ。`-dust` の読み方。NCBI blastn に無い LOSAT だけのオプション（`-verbose` など）を CLI から除いた（§B、§I） |
| 入力の読み方 | `U` を `T` として読む。`bio` と NCBI の `CFastaReader` で読み方が違う入力（定義行の先頭の空白・制御文字・非 ASCII・`>?`、IUPAC 以外の残基、配列の行の非 ASCII のバイト、残基の無いレコード、`bio` が読めない FASTA）を明示的に拒否する（計画 TD-12）。NCBI の「Title ends with at least 20 valid nucleotide characters」の警告（§E、§M、§N） |
| query の batch と表形式 | 最初の batch の後の短い無効な query の連なりを検索する。outfmt 7 は検索されなかった batch の query に件数の行を出さない。subject の総文字数が 2^31 以上、予備の hit list の大きさが負になる `-max_target_seqs` は明示的に拒否する（§F） |
| 独立監査で見つかった、以前からの不具合 | greedy な伸長のメモリ（§H）、1 つの query の対角線の表の大きさ（§J）、gap の縮小の subject の先頭の panic、予備の段階の `-subject_besthit`、bit score の幅（§K）、subject の曖昧な文字の `CRandom`（§L。TBLASTN と共有）、lookup table の見積もり（§L、§M）、outfmt 0 の subject の title（`CDeflineGenerator::GenerateDefline`、§N）、`-lcase_masking` の subject の走査の区間と得点 0 の HSP（§O）、1 文字の HSP の ` Strand=`（§P）、小さい lookup table の word の延長と確かめ（§Q）、小さい e-value の cutoff の `Int4` の変換（§R） |
| 明示的な拒否 | NCBI と同じにできないもの・NCBI が落ちるものは「not supported by LOSAT's BLASTN」か「which LOSAT does not reproduce」で失敗する（一覧は §G） |
| ABI v1（計画 TD-1） | 引数・形式・誤りの扱いは S07+ の前のまま（fail-fast の修正を除く）。v1 は `run_local` で CLI と同じエンジンを使うので、エンジンを NCBI に合わせる修正は v1 の検索結果にも及ぶ |
| fixture | `LOSAT/tests/outfmt0_manifest.tsv` に 6 件を足した（47 件）：reward 1 の得点、表を超える gap（megablast、IUPAC の query）、公約数を持つ得点、RNA の query、最初の batch の後の短い無効な query |

### アダプタと文書

- `web/adapter/src/store.rs`：BLASTN の `register` が CLI と同じ入力の検査（定義行、配列の行、残基、残基の無いレコード、`bio` が読めない FASTA）を行う。入力の handle を program に結びつける。
- `web/adapter/src/run.rs`：`validate` が BLASTN の `-dust` と得点を CLI と同じく検査する。
- `docs/web/abi_v2.md`：BLASTN の拒否、1 文字の HSP の鎖（§8）。`docs/web/verification_cells.tsv`：BLASTN の S07+ の升目（fixture、得点の sweep、入力の検査、word size・切り出し・title の sweep）。
- `LOSAT/tests/check_wasm_threading.py`：BLASTN の独自の形式の拒否の文言（「not supported by LOSAT」）を受け付ける（`2b97f22f9`。拒否は以前と同じく、worker を起動する前で、出力を書かない）。
- 計画：TD-12（CFastaReader）、TD-13（S08+、ほかの program の既定以外のオプション）、TD-14（S07++、BLASTN の query の batch）、TD-1 の補足（v1 の結果と誤りの文言）、DW-7 の改訂（エンジン側とアプリ側の 2 本、2026-09-30 の保守者の承認）。S12 の指示書に BLASTN の得点と入力の扱いを、S13 の指示書に 1 文字の HSP の鎖を書き足した。

## 完了条件と結果

| 完了条件（計画 §7 の S07+） | 結果 | 証拠 |
|---|---|---|
| S06 の `scoring_sweep.py` の全組合せが、NCBI と同じ拒否、outfmt 6 のバイト一致、明示的な拒否のどれかになる | 通過。80 の組合せのうち 46 が outfmt 6 のバイト一致、34 が NCBI と同じ終了コードの拒否（1 が 10、3 が 24。S06 の LOSAT はこの 34 を受け付け、17 で結果が違った） | `scoring-sweep-s06.tsv` |
| 既定以外の得点（S07+ の sweep）| 通過。22 の reward / penalty の組 × 10 の gap × 2 つの task × 2 つの query の 880 の組合せが、outfmt 0・6・7 のそれぞれで、300 がバイト一致、580 が NCBI と同じ stdout・stderr・終了コードの拒否（明示的な拒否 0）。S07 の成果物では outfmt 6 で 806 が違った | `scoring-sweep-fmt{0,6,7}.tsv`、`scoring-sweep-fmt6-before.tsv` |
| 直した組合せの fixture で NCBI とバイト一致 | 通過。47 件（S07+ で足した 6 件を含む）すべてが `-num_threads` 1・2・4 で stdout と stderr まで NCBI BLAST+ 2.17.0 とバイト一致。47 件の outfmt 6 も一致（`precheck.tsv`）。NCBI の再実行は manifest のハッシュと一致（`oracle-check.log`） | `check-losat-n{1,2,4}.tsv`、`precheck.tsv`、`oracle-check.log` |
| 入力の読み方（指示書 6・7） | 通過。`check_inputs.py` の 292 件すべてが期待どおり（same、same-error、LOSAT の明示的な拒否、引数の誤り、both-fail）。outfmt 0 の句読点の title の 1023 件は、957 がバイト一致、66 が NCBI の落ちる定義行の明示的な拒否（`title-sweep.tsv`） | `check-inputs.tsv`、`title-sweep.tsv` |
| 独立監査で見つかった以前からの差の sweep | 通過。word size の sweep 256 件（S07 の成果物で 22 件、第 6 回の修正の前で 16 件が違った）、切り出しの sweep 17 の設定（viral 6 × 120 組、ambiguity 5 × 100 または 40 組、lcase 3 × 120 組、iupac 3 × 300 組）がすべて差 0。修正の前の LOSAT での差は `*-before*.tsv`（例：lcase の `-task blastn -word_size 7` で 20 組、iupac の既定で 3 組）。各回の LOSAT の実行ファイルのハッシュは `round-binaries.sha256` | `word-size-sweep.tsv`、`word-size-sweep-before.tsv`、`slice-sweep-*.tsv` |
| 既存の BLASTN のゲートと S07 の fixture に退行なし | 通過。全 program の 236 件（BLASTN 14、BLASTP 28、TBLASTX 20、TBLASTN 162、BLASTX 12）の出力・stderr・終了コードが S02 の基準と一致（差 0）。Gate A の BLASTN の行を含む | `capture-compare.txt`、`capture/hashes.tsv`、`check-losat-n*.tsv` |
| 拒否がアダプタの `validate` にも出る、V-ABI | 通過。full：115 の検索を serial reactor（1 スレッド）と threaded reactor（1・2・4 スレッド）で実行した 460 件すべてで、全形式・HSP レコード・診断がネイティブの CLI と一致。凍結ハッシュは 676 件中 672 件が一致し、一致しない 4 件は S02 からの既知の `Sakai.MG1655.megablast` の outfmt 7。quick：52 件すべて一致。アダプタの `validate` と `register` の拒否はアダプタの試験（`adapter-test.log`）で確かめた | `v-abi-full/`、`v-abi-quick/` |
| （規則 4）エンジンの変更の既存のゲート | 通過。v1 の WASI の検査（423 件の command / oracle の記録、形式の失敗 0）、`v1_requests.js` の 14 件の応答が S05 の記録と一致、wasm32 だけでコンパイルされる v1 の試験 5 件が通過。v1 の reactor の記録は、threaded の 172 件のうち 21 件の誤りの文言が S07+ で足した句だけ違い（下の「残件と扱い」）、serial の 12 件は一致 | `wasm-threading.log`、`v1-reactor-records-compare.txt`、`v1-requests-compare.txt`、`wasm32-web-api-tests.log` |
| （規則 4）V-PERF の非退行 | 通過（下の「性能の計測」） | `perf-1.json`、`perf-check-1.txt`、`perf-2.json`、`perf-check-2.txt` |
| `cargo fmt --check`、`clippy -D warnings`、`cargo test --all-features` | 通過。エンジンは 4 つの構成の clippy と 823 件の試験、アダプタは 3 つの構成の clippy と 7 件の試験 | `fmt.log`、`clippy.log`、`adapter-clippy.log`、`cargo-test.log`、`adapter-test.log` |
| ビルドの同一性（TD-6、TD-11） | 通過。reactor には `/losat` と `/cargo` に置き換えたパスだけが残る | `reactors/build-identity.json` |
| 独立監査 | 16 回目で supported | 下の「独立監査」 |

## 独立監査

読み取り専用の独立監査（役割 `ncbi_parity_auditor`、コードを変えない別のエージェント）を 16 回受けた。1〜15 回目はどれも unsupported で、見つかったものを直すか明示的に拒否し、次の回で確かめた。16 回目は supported。見つかったものの NCBI の経路と直し方は `AUTHORITY.md` にある。

| 回 | 主な指摘（◎は S07+ の前からの、得点のオプション以外にも及ぶ差） | 修正 | 記録 |
|---|---|---|---|
| 1 | ◎ greedy な伸長のメモリと時間が gap のコストに比例する。無効な context の扱い。検査の順（word size 100 超、e-value 0、`-out`）。LOSAT の上限（reward 0、3000 を超える得点、32767 を超える greedy の gap）、`BATCH_SIZE`・`CHUNK_SIZE`、空の定義行、残基の無いレコード | `5377fc1bd`、`0c34edebe` | §G、§H |
| 2 | reward と penalty の 16 ビット、空の subject の engine error、e-value の読み方、ファイルを開く順、パイプ | `2fd0ee2c1`、`c27fbf93b` | §A、§B、§E |
| 3 | subject を読む時点（引数の処理の中）、入力を 1 度だけ開く、ディレクトリ、`-evalue` の最初の文字、`bio` が読めない FASTA の文言 | `c1e933757`、`295cf0b40` | §B、§E |
| 4 | 整数の 16 進数、`-outfmt` の解析の時点と文言、`-query` の既定値 `-`、`-out -`、LOSAT だけのオプション、上限の時点 | `637426f30`、`0c7985789` | §B、§I |
| 5 | `delim` の検査、`0x` だけの値、予備の hit list の 32 ビット、`-dust` の読み方、`-out` を 1 度だけ開く、名前の長さ、残基の無い subject の時点 | `7d6fee2a7`、`6c8826068` | §B、§F |
| 6 | ◎ 1 つの query の two-hit の対角線の表を一方の鎖から求めていた（既定のオプションでも HSP が欠けた）。NCBI の他の task と 53 のオプションの明示的な拒否 | `c33731d08`、`722e6b95d` | §G、§J |
| 7 | ◎ gap の縮小が subject の先頭を越えて panic（既定の megablast）。◎ 予備の段階の `-subject_besthit`。◎ bit score の幅。複数の query の batch（S07++ に移した） | `010bee067` | §K、TD-14 |
| 8 | ◎ subject の曖昧な文字を固定の対応で詰めていた（NCBI は `CRandom`、既定のオプションでも差）。lookup table の見積もり。ABI v1 の空の query | `6b176da78` | §L |
| 9 | `>?` の gap の行、「Title ends with at least 20 valid nucleotide characters」の警告、覆われた context の見積もり | `d46906312` | §M |
| 10 | 配列の行の Unicode の空白、◎ outfmt 0 の subject の title（`GenerateDefline`） | `f8e4b0426`、`a9e17ca43` | §N |
| 11 | NCBI が落ちる句読点だけの title | `9adac5b1e`、`e4409ccfa` | §N |
| 12 | ◎ `-lcase_masking` の subject の走査の区間（既定のオプションでも HSP が欠けた）。得点 0 の HSP（第 8 回の修正の後）。文言 3 件 | `07e8c4fa8` | §O |
| 13 | 1 文字の minus 鎖の HSP の ` Strand=`（第 8 回の修正の後に出る） | `4cab5c3fb` | §P |
| 14 | ◎ query の曖昧な文字：小さい lookup table の word の延長（圧縮した query）と word の確かめ（既定の megablast で、NCBI が捨てる HSP を出した）。ABI v2 の 1 文字の HSP の鎖、引用 | `181c18ac1` | §Q |
| 15 | ◎ 1e-297 以下の `-evalue` の大きな探索空間の cutoff（NCBI の x86-64 の `Int4` の変換） | `bbc9c7809` | §R |
| 16 | supported。主張の範囲に黙った差は無し。第 15 回の変換がオラクルの実行ファイルの `cvttsd2si` と同じことを逆アセンブルで確かめた。1 件（`-task blastn -evalue 1e10 -dust "5 10 2"` の 2 つのゲノム）は NCBI が 43 分で終わらず比べていない | — | — |

### 性能の計測

変更前は S07 の成果物（`s07-final/`）、変更後はこの記録の成果物で、1 回ずつ交互に測った（`perf_cases.py`、1 回の暖機と 3 回の計測）。case は、BLASTN の `blastn`（起動の時間が大半）、`blastn-large`・`blastn-large-fmt7`・`blastn-large-fmt0`（Gate A の EDL933 × Sakai の megablast）、`blastn-many`（260 の query、`-task blastn`）。

- **最初の計測（`perf-1`）：** 15 件のうち 12 件が +5% の内で、3 件が超えた（出力はすべて同じ）。`blastn` の serial WASI（×1.169、変更後の範囲 [0.139, 0.226]）、`blastn-large-fmt0` の threaded WASI（×1.271、範囲 [1.705, 2.292]）、`blastn-many` の threaded WASI（×1.070、範囲 [0.321, 0.389]）。どれも変更前と変更後の範囲が重なる。
- **切り分け（`perf-2`）：** AGENTS.md の手順（結論が出ないときだけ回数を増やす）に従い、この 3 case を 10 回ずつ交互に測った。9 件すべてが +5% の内（×0.623〜×1.038）で、出力は同じ。同じ実行ファイルでも 2 回の計測で絶対時間が大きく違い（例：`blastn-large-fmt0` の threaded WASI の変更前の中央値が 1.722 秒と 3.729 秒）、計算機の負荷の揺れが最初の超過の原因と判断した。
- 計測の間、アプリ側の試験とビルドは止めた（計画 DW-7、`vperf.lock`）。

## 残件と扱い

- **S07++（TD-14）**：複数の query を NCBI の batch に分けて検索する（最初の batch の後の batch の大きさ、batch の端の近似の ungapped 伸長、lookup table と対角線の表の大きさ）。S07+ の batch に依存する明示的な拒否を置き換える。NCBI の query の分割（`CQuerySplitter`、batch の query が 2 × (chunk − 100) 以上）も S07++ で調べる（S07+ の終わりの確認では、分割の境をまたぐ 2.1 Mb の query で NCBI とバイト一致した）。指示書に、S07+ の終わりに読んだ `CBatchSizeMixer` の式と `good_init_extends` を数える箇所を書いた。
- **S08+**：出力の書き込みの失敗、UTF-8 でないファイル名と引数の値の表示、`-num_threads` の NCBI の警告（計画 §10）、LOSAT が対応しない出力形式の NCBI 自身の誤り、引数の誤り（NCBI の USAGE と終了コード 1、LOSAT の clap と終了コード 2）、ほかの program の既定以外のオプションと outfmt 0 の title（TD-13）。
- **計画 TD-12（§10）**：NCBI の `CFastaReader` の移植（今は読み方が違う入力を明示的に拒否する）。
- **ABI v1（TD-1）**：v1 の reactor の記録は、スレッドの上限の誤りの文言の 21 件が、S07+ で足した句（`, which is not supported by LOSAT`、`; a thread count that the system cannot start is not supported by LOSAT`）だけ S05 と違う。誤りの扱い（検索の前に失敗し、スレッドを起動しない）は同じ。計画の TD-1 に書いた。
- **保守者への報告（判断は不要）**：◎ の差（§J、§K、§L、§N、§O）は、BLASTN の以前からの不具合で、既定のオプションでも結果が NCBI と違いうるものだった（1 つの query の逆位の反復、subject の曖昧な文字、tandem repeat の中で始まる subject、`-lcase_masking` の小文字の subject など）。BLASTN の回帰の fixture（Gate A と S02 の基準）はこれらの条件に達しないので、凍結ハッシュは変わらない（capture の差 0）。TBLASTN（TLOSAN v0.2.0 の認証済みのコード）は、`CRandom` を共有にしただけで出力は変わらない。bit score の幅（§K）は全 program が共有する関数で、既定の得点では出ない。
- **ABI v2 の 1 文字の HSP の鎖**：HSP レコードの座標は、1 文字の HSP の鎖を表せない（`docs/web/abi_v2.md` §8）。鎖の欄を足すかを S13 で決める（指示書に書いた）。
- **比べていない 1 件**：第 16 回の監査の、`AP027131 × AP027132 -task blastn -evalue 1e10 -dust "5 10 2" -max_target_seqs 1 -outfmt 0` は、NCBI が 43 分で終わらなかったので比べていない（LOSAT は終わる）。極端に大きい e-value と弱い DUST の組合せで、主張の範囲の中だが、比べられる時間で終わる入力ではない。
- **ゲートの実行**：修正のたびにゲートを実行し直し、後の修正で置き換えた実行の記録は残していない（この記録の実行だけを残す）。修正の前の LOSAT での差は、各回の実行ファイル（`round-binaries.sha256`）で同じ実行の中で測り直した（`*-before-round*-fix.tsv`）。
- **保守者の確認待ち（以前から）**：E1a〜E1c の V-PERF の判断、E1c の CLI の 2 つの振る舞いの差。

## 実行の方法（再現）

```bash
# 入力を作る（既存のファイルと同じバイトになる）
python3 docs/evidence/losat_web_e2c/make_inputs.py
# NCBI との比較（NCBI BLAST+ 2.17.0 を比較のときに実行する）
python3 docs/evidence/losat_web_e2c/scoring_sweep.py --bin-dir <NCBI> --losat <LOSAT> --outfmt {0,6,7} --jobs 8
python3 docs/evidence/losat_web_e2c/check_inputs.py --bin-dir <NCBI> --losat <LOSAT> --work <dir>
python3 docs/evidence/losat_web_e2c/word_size_sweep.py --bin-dir <NCBI> --losat <LOSAT> --jobs 8
python3 docs/evidence/losat_web_e2c/slice_sweep.py --bin-dir <NCBI> --losat <LOSAT> --work <dir> \
  [--pool viral|ambiguity|lcase|iupac] --options="<options>" [--cases N] --jobs 8
python3 docs/evidence/losat_web_e2c/title_sweep.py --bin-dir <NCBI> --losat <LOSAT> --work <dir> --jobs 8
# S06 の sweep、fixture、回帰の出力、v1 の WASI の検査、V-ABI、性能：S07 の記録の「実行の方法」と同じ
python3 docs/evidence/losat_web_e2a/scoring_sweep.py --bin-dir <NCBI> --losat <LOSAT>   # docs/evidence/losat_web_e2a から
python3 docs/evidence/losat_web_e2a/check_losat.py --losat <LOSAT> [--threads N]
python3 docs/evidence/losat_web_e2c/perf_cases.py run --before <S07 の 3 つ> --after <S07+ の 3 つ> \
  --cases blastn,blastn-large,blastn-large-fmt7,blastn-large-fmt0,blastn-many --out <file.json> [--repeat 10]
```

`<NCBI>` は NCBI BLAST+ 2.17.0 の `bin`、`<LOSAT>` は `cargo +1.92.0 build --release --locked` の実行ファイル。

## 引き継ぎ（S07++ へ）

- S07++ の指示書に、S07+ の終わりに NCBI のソースで読んだ batch の経路（`CBatchSizeMixer` の式と `Int4` の変換、`good_init_extends`、query の分割）を書き足した。変更前の成果物は、この記録の成果物（`/home/kawato/.cache/losat-web-gui-target/s07p-final/`、ハッシュは `run-20260930T022908Z/artifacts.sha256`）。
- S08+ の指示書には、S07+ で分かった NCBI の振る舞い（引数の読み方と順序、`report/defline.rs`、CLI の誤りの経路、第 12 回の得点 0 の HSP と subject の小文字の区間）を、ほかの program でも確かめることを書いてある。
- アプリ側（計画 DW-7）：S10 を別の worktree で並行して進めている。S07++ のゲートの V-PERF の間は `vperf.lock` でアプリ側の試験を止める（ゲートの script が取る）。
