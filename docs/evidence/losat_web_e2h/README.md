# E2h（SF・SFb・SFc・SFd）ゲート記録：FASTA の読み方（NCBI `CFastaReader` の移植）

状態：**途中**（2026-10-10、SFc を監査の修正の後で区切った）。完了条件（計画 §7 の SF の行）はまだ満たしていない。続きは [SFd の指示書](../../losat_web_gui_sessions/session_sfd_e2h_final_gate.md)。

## SF（2026-10-06 開始、10-06 に `/mnt/c` の障害で中断、10-08 に Linux の clone で再開）で行ったこと

| 手順 | 結果 | コミット |
|---|---|---|
| 0. 準備 | `origin/main`（#113〜#115）を merge（衝突なし）。変更前の成果物（S11 の 8 つ）は `~/.cache/losat-web-gui-target/sf/bin/`、capture の基準は `sf/capture-before/hashes.tsv`（236 件、10-06 に取得） | `00fb6c0e` |
| 取り残された作業の復元 | 10-06 に未コミットだった読み込み器の module と棚卸しの brief を、セッション監査の写し（`stranded-backup-20261008`、SHA-256 を照合）から戻した | `eba45f75` |
| 読み込み器の単体試験 | 初めて流した 9 件のうち 4 件が失敗していた。4 件とも試験の期待が NCBI と違い、読み込み器は NCBI と同じだった。期待を NCBI のソースと oracle の観察に合わせた | `04cc336d` |
| 2. 棚卸し | [`INVENTORY.tsv`](INVENTORY.tsv)（222 行）。LR 行の読み込みと流れ 30、RD レコード 49、BI blastinput の層 55、RP 報告 56、AD アダプタと ABI 32。`status_before`（`f3048ffde`）：rejected-now 74、faithful 52、divergent 46、n/a 26、unported 22、exception 2。範囲ごとの指示と覚え書きは `inventory/`、oracle の入力と出力は `$BUILD_ROOT/sf-e2h/inventory/scratch_*`（追跡しない） | `e1eed966` |
| 1. 権威の記録 | [`AUTHORITY.md`](AUTHORITY.md)（§A〜§L）。§K に NCBI の不具合に見える挙動 12 件（保守者への一括の問い、下記） | `2a903cd7` |
| 3. 移植の計画 | [`PORT_PLAN.md`](PORT_PLAN.md)：段階 S0〜S10、記録の型、危険、未決の問い（判断は下記） | この記録と同じコミット |
| 3. 移植 S0 | 読み込み器の下準備（program への接続は無く、出力は変わらない）。Seq-id の判定を `CSeq_id::Set` の解析に合わせた（`0123`・`ACGT ACGT`・`contig1` は FASTA として読む。広めの拒否は accession の guide の形だけ）。行の読み方の試験（BLASTX の C++ oracle 12180 行、LR の癖）、`from_bytes` と file の一致の性質試験、まとめて読む経路（R1） | `24f56b17`、`905597ba` |
| 5. fixture（NCBI 側） | `LOSAT/tests/fasta_input_fixtures.py` と program ごとの manifest（`fasta_input_fixtures/*.tsv`、528 KB）、入力 683（387 KB）。3932 行：NCBI 2.17.0 を固定した 3908 行と、Seq-id の拒否の 24 行（文言は移植で埋める）。3 回の凍結でハッシュが一致 | `352ad49f` |

### 検証

- `cargo fmt --check`、release build、clippy `--all-targets --all-features -D warnings`、`cargo test --all-features`（`fasta_reader`：20 通過、速度の試験 1 は ignored）。ログは `/home/kawato/losat-baselines/sf-e2h-20261008/logs/`。
- fixture の `check` を変更前の実行ファイル（S11）で流した：same 1131、differs 13（BLASTN の `>` だけの subject：NCBI は `Subject sequence contains no data` の警告と終了コード 3 の `The average subject length is too short`、LOSAT は `Empty CBlastQueryVector`。移植の S4 で直す）、LOSAT の拒否 2764、pending 24。fixture は変更前の拒否と移植を区別する。
- 読み込み器だけの速度（release、3 回の中央値、`bio::io::fasta` との比）：100 Mb の genome（80 桁）1.6 倍、1 行の genome 1.2〜1.3 倍、10 万 × 300 塩基 1.5 倍（まとめて読む経路の前は 4.6〜6.0 倍）。V-PERF は program に接続した後。

## SFb（2026-10-08〜09、エンジン側。アプリ側の S13 と並行）で行ったこと

移植は `PORT_PLAN.md` §2 の段階ごとに 1 体の implementer に任せ、各段階で quick と standard の段階（clippy 4 構成、`cargo test --all-features`、`check_losat` 1/2/4 スレッド、`ci_fast_regressions.py` 全 program、BLASTN 183・TBLASTX 84・範囲 35 の回帰 fixture、`web_api`・アダプタ・`run_local` を変えた段階は V-ABI quick・v1 の WASI の行列・`v1_requests.js`）と、`fasta_input_fixtures.py check`・sweep・棚卸しの oracle の再生を流した。NCBI の新しい実行は、implementer が確かめられなかった端の場合だけをメインが `unshare -rn` の下で少数流した（下の「oracle で決着させた端」）。

| 段階 | 内容 | コミット |
|---|---|---|
| S1 共有の部品 | `read_subjects`、`stream_is_empty`、`QueryRecords`/`QueryEnd`、`FastaRecord::cut`・`from_bio`、`DATA_LOADERS` の実効値（環境変数 > `<prog>.ini` > `.ncbirc`）、`InputRecord`。出力は変わらない | `7e49a1be` |
| 引用の行番号 | S0 の読み込み器のコメントの NCBI の引用 14 か所を固定 commit に合わせた（`verify_refs` 14 → 0） | `e4de0af6` |
| S2 報告の byte の層 | `report/defline.rs`（byte の GenerateDefline、HtmlDecode、GuessEncoding、Windows-1252、`x_CleanAndCompress`）、折り返し・`Query=`・`# Query:`・警告の切り詰めの byte の形。受け付ける入力の出力は変わらない | `ebaeeaf5` |
| sweep | [`fasta_sweep.py`](fasta_sweep.py)（2280 件）と [`check_inputs.py`](check_inputs.py)（1032 件）、NCBI 側を 2 回固定して一致（[`SWEEPS.md`](SWEEPS.md)、[変更前の数](sweeps-before-counts.md)） | `60cf5823` |
| S3 BLASTN の読み込みと報告 | BLASTN が `CFastaReader` の移植で読む。空の pipe（零 batch）、空白の pipe、`-` を NCBI の `cin` と同じに読む（§K-4 の尾の消失を含む）、Seq-id の行の明示的な拒否（§G4 の文言、終了コード 1）、RP-20 の明示的な拒否、`Subject_` の sseqid | `6c5d4220` |
| S4 BLASTN の中身の無いレコードと batch の時点の誤り | 空の query・subject・範囲、batch ごとの読み込みの誤り、`Empty CBlastQueryVector`、`The average subject length is too short` | `8f8e8e08` |
| S5 TBLASTX | 両方の役割を核酸の読み込み器で。v1 は今の検査を前の位置で行う（`V1Records`） | `c7061c08` |
| S6 TBLASTN の subject | subject を引数の処理の中で読む。全部空の subject の `Effective search space used`（blast_setup.c:716-732）も移植 | `5fbdf4c1` |

### 検証（S6 の後）

- 各段階の最後の commit で、fmt、clippy 4 構成、`cargo test --all-features`（S6：967 通過、失敗 0）、`check_losat` 204 × スレッド 1/2/4、回帰 fixture 183/84/35、`ci_fast` 223 件（失敗 0、許可した既知の不一致 1：`Sakai.MG1655.megablast`）、アダプタの試験、wasm32 の `web_api`、reactor とビルドの同一性、v1 の WASI の行列 433 件、`v1_requests.js`、V-ABI quick 60/60（凍結ハッシュ 16/16）。受け付ける入力で変わった出力は、意図した outfmt 6/7 の `Subject_` で始まる subject の sseqid だけ（RP v01・v09、NCBI と一致）。ログは `/home/kawato/losat-baselines/sfb-e2h-20261008/logs/S1〜S6/`、段階ごとの引き継ぎは同じ folder の `handoff-S1〜S6.md`。
- fixture（`fasta_input_fixtures.py check`、program ごとの行の結果。変更前は S11 の成果物）：

| program | 変更前（S11） | S6 の後 | 残りの行 |
|---|---|---|---|
| BLASTN | same 561、differs 13、rejects 1395、pending 12 | same 1969、pending 12 | Seq-id の行 12（明示的な拒否、期待の固定は SFc） |
| TBLASTX | same 176、rejects 465、pending 2 | same 625、rejects 16、pending 2 | `-lcase_masking` 16（E2e の明示的な拒否）、Seq-id 2 |
| TBLASTN | same 151、rejects 337、pending 4 | same 297、rejects 191、pending 4 | query の側（S8）、Seq-id 4 |
| BLASTP | same 243、rejects 567、pending 6 | 変わらない | S7 |

- sweep（S6 の後、移植した program と役割）：BLASTN・TBLASTX・TBLASTN の subject は、一致、同じ誤り、記録した明示的な拒否（Seq-id の行）、承認済みの例外（句読点の題で NCBI が落ちる）だけ。`check_inputs.py` の想定外は 0。
- oracle で決着させた端（NCBI 2.17.0、`unshare -rn`、記録は `/home/kawato/losat-baselines/sfb-e2h-20261008/probes/`）：registry の file の空の `DATA_LOADERS=` は項目なしと同じ（data loader は有効。SF の判断どおり、S1 の implementer の別の読みは誤り。環境変数の空の値は別：判断 5、SFc の監査 B-1）、split される query の前の空のレコード（`CHUNK_SIZE=20000`、59 行で一致）、TBLASTX の空の query と短い query の batch（outfmt 0/6/7 で一致）、TBLASTN の `-subject_loc` が record の終わりを越える 3 通り（一致）。

## SFc（2026-10-09〜10、エンジン側。アプリ側の S13 は SFc で merge 済み）で行ったこと

SFb と同じく、移植は段階ごとに 1 体の implementer に任せ、各段階で quick と standard の段階、`fasta_input_fixtures.py check`、2 つの sweep、棚卸しの oracle の再生を流した。保守者の指示（2026-10-09「とりあえず main に適宜マージして」）で、検証の済んだ区切りごとに `main` へ PR を出して merge した（計画 DW-20 の「今は PR を作らない」を置き換える）。

| 段階 | 内容 | コミット |
|---|---|---|
| S7 BLASTP | 両方の役割を蛋白の旗の読み込み器で読む。v1 は今の検査を前の位置で（`V1Records`）。中身の無い subject を走査しない（buffer の 1 語先を読んでいた不具合も直った、blast_engine.c:1429-1432） | `0c6721ea` |
| S8 TBLASTN の query と後片付け | query を batch ごとに蛋白の旗で読む。`blastn/input.rs` と `InputRecord` の bio の実装を消し、CLI の file の部品は `blastinput/input_files.rs` へ。`from_bio` は題の先頭の空白を落とす。`BLASTINPUT_GEN_DELTA_SEQ` は出力を変えない（RD の 6 実行、有無とも NCBI と 12/12 一致） | `34db06aa` |
| S9 `scan` の種類 1・2 | `web/adapter/src/scan/ncbi.rs`（chunk で与える読み込み器の行と record の移植）、性質試験 `scan_ncbi_properties.rs`。1430 の入力 × 2 種類 × 10 の chunk の大きさでエンジンの読み込み器と差 0 | `27206434` |
| S10 `register`・`run`・文書 | `register` は種類 1・2 の scan の後にエンジンの読み込み器で読む。`run` は `FastaRecord` を直接渡す。`abi_v2.md` §4・§8・§9。V-ABI の `fasta_input` の群（528 run、NCBI の凍結ハッシュ 796/796）。ABI v1 の `bio` の読み込みをエンジンから `web_api/v1_*.rs` へ（`fasta::Record` は v1 の層と種類 0 の scan にだけ残る） | `e57ef708`、`02277def` |
| fixture の仕上げ | Seq-id の 24 行と `-lcase_masking` の 34 行を `losat-rejection` の行として LOSAT の終了コード・文言・ハッシュで固定（`freeze-rejections`）。`ci_fast_regressions.py` が program ごとに fixture を確かめる（1 program 約 2 秒） | `32dfe90c` |
| ゲートの script | TBLASTN の Stage G の path を `/mnt` の tmpfs の中に作って worktree を bind（`/mnt/c` を読まない）、wasm32 の試験の絞り込み `web_api::` | `9dfe7efd` |
| アプリの試験（tests only） | エンジンの契約の試験 3 file を ABI v2 の新しい契約に合わせた（種類 1 の拒否 → 未知の種類 3 の拒否、`register` が NCBI と同じく読む入力、engine の E2E の拒否の入力を `>?10` に）。アプリの `src/` は変えない | `4596b82f` |
| 監査の修正 | B-1（空の環境変数 `DATA_LOADERS=` は項目あり、data loader は無効。env_reg.cpp:157-167）、A-3（BLASTP の outfmt 6/7 と custom の tabular は書いたときに flush する）、A-1（gap-type は最初の NUL までで引く）の修正と試験。D-1 の文書（`abi_v2.md` の `register` の行）と D-2 の一覧（判断 12）。ゲートの script の試走からの修正。監査の再現 178 件が NCBI とバイト一致（`unshare -rn`）。fixture は 61 行増えて 3993 行（same 3935、rejection 58） | `321aa9ca`（B-1・A-3・A-1 と試験）、`0be1afec`（D-1・D-2）、`34dc8f99`（ゲートの script） |

`main` への merge：#116（S11・W1・W3・SF・SFb・S7、`03947fb9`）、#117（S8、`e31236bf`）、#118（S9・S10・fixture・ゲートの script・アプリの試験、`feea1cea`）。CI の firefox と webkit の E2E に 5 秒の待ちの timeout が 1 件ずつ 4 回出た（毎回別の試験、再実行で通過、手元の E2E 92/92）。

### fixture（`fasta_input_fixtures.py`、3932 行。監査の修正の後は 3993 行）

| program | 変更前（S11） same / rejects / differs / pending | S6 の後 | 最後（02277def） same / rejection |
|---|---|---|---|
| BLASTN | 561 / 1395 / 13 / 12 | 1969 / 0 / 0 / 12 | 1969 / 12 |
| TBLASTX | 176 / 465 / 0 / 2 | 625 / 16 / 0 / 2 | 625 / 18 |
| TBLASTN | 151 / 337 / 0 / 4 | 297 / 191 / 0 / 4 | 488 / 4 |
| BLASTP | 243 / 567 / 0 / 6 | 243 / 567 / 0 / 6 | 792 / 24 |
| 計 | 1131 / 2764 / 13 / 24 | 3134 / 774 / 0 / 24 | 3874 / 58 |

`rejection` は Seq-id の 24 行（終了コード 1、§G4 の文言）と TBLASTX・BLASTP の `-lcase_masking` の 34 行（終了コード 2）。sweep：`fasta_sweep.py` 2280 件は same 2134・same-error 68・明示的な拒否 70・承認済みの例外 8、記録の無い拒否・差・時間切れ 0（変更前は記録の無い拒否 1513）。`check_inputs.py` 1032 件は same 841・same-error 149・明示的な拒否 23・承認済みの例外 19、想定外 0（変更前は 485）。

### ゲートの試走（run `20261009T145447Z`、HEAD `9dfe7efd`、監査の修正の前）

run `20261009T145447Z`（HEAD `9dfe7efd`）は 23 工程を流した。build 7 分、fast-all 236 件 失敗 0（既知の許可 1、Stage G を含む）、lint・clippy・tests・pychecks・quick-fixtures・regression-fixtures・sf-fixtures・sf-sweeps・input-sweeps・v1-wasi・vabi-quick・option の sweep（BLASTP 10 分、TBLASTN 10 分、TBLASTX 20 分）・V-ABI full（BLASTN 3 分、BLASTP 1 分、TBLASTN 1 分、TBLASTX 3 時間 17 分）は終了コード 0。失敗は 3 つで、どれも LOSAT の差ではない：`oracle-check` 2 回（TMPDIR。`-db_gencode` の代わりの DB の報告の `# Database:` の行が `$TMPDIR` を含み、凍結は `/tmp`。`/tmp` に戻すとハッシュは manifest と一致）、`capture`（236 件、S02 の基準とも S11 の build とも差 0。終了コード 1 は既知の許可 `Sakai.MG1655.megablast` だけ）。`collect` は `sweeps/` の下位の directory を作らずに止まった。script を直した（`34dc8f99`）。一部の run の記録は `/home/kawato/losat-baselines/sfc-e2h-20261009/gate-trial-record/`。監査の修正の後のコミットでのゲートは SFd。

### 独立監査（4 観点、監査したコミット `9dfe7efd`、sonnet、読み取り専用）

| 観点 | 比較の数 | 結論 | 見つけたこと |
|---|---|---|---|
| (a) 棚卸しに沿った経路の網羅 | 約 13,200 run、Seq-id の分類 約 8,000 | unsupported | 222 行：ported 183、rejected 18、exception 19、divergent 2、missing 0。NCBI の引用 363 件の機械的な照合は通過。A-3（中）BLASTP の outfmt 6/7 で、後の batch の読み込みの警告が前の batch の行より先に出る（`2>&1` の順）。A-1（低）`[gap-type=…\0…]` の NUL 以降を NCBI は無視する。A-2（低）guide の表に無い prefix の accession の形の行（`NB798107`）を拒否する（記録済みの広めの拒否、下の判断） |
| (b) 核酸の入力 | 約 6,370 | unsupported（狭い） | B-1（中〜低）空の環境変数 `NCBI_CONFIG__BLAST__DATA_LOADERS=` は NCBI では項目あり（data loader は無効、env_reg.cpp:150-160）。LOSAT は項目なしとして Seq-id に見える行を拒否していた |
| (c) 蛋白の入力と残す拒否 | 約 7,200 と最初の行 1,300 | supported（A-3 を除く） | A-3 を確かめ広げた（全ての警告の種類、threads 1/2、BATCH_SIZE 60/100/既定）。残す拒否の文言・終了コード・時点・理由・範囲は記録どおり |
| (d) アダプタ | 約 190,000 | unsupported | `scan` の種類 1・2、`register`、`q_idx`/`s_idx`、`run` と CLI、種類 0、v1（変更前の reactor と比較）は一致。D-1（低）両方の入力に不備があるとき、CLI は subject の誤りを先に出し、Web は `register` を呼んだ順。D-2（中〜低）ABI v1 の出力が判断 12 に書いた 2 種類より広く変わる（BLASTP の空・空白だけ・先頭が空白の定義行の outfmt 0/6/7、TBLASTX の outfmt 6）。共有の報告の層が NCBI の移植になったため |

修正は `321aa9ca`（B-1・A-3・A-1）と `0be1afec`（D-1 は文書、D-2 は NCBI との一致の確認と記録）。修正の再監査と、最後のコミットでのゲートの全工程・Gate A・V-PERF は SFd で行う。

### アプリ側への引き継ぎ（SFc で増えたもの）

- `register` に残る拒否（`>?` の行、Seq-id の行、読み込みの誤り）の文言は行の番号を名指し、`query record N` を名指さない。アプリの `RECORD_IN_MESSAGE`（`draft.ts:114`）で拒否されたレコードの印と除外のボタンが出ない。行の番号から自分の索引でレコードを引くか、`register` の誤りの JSON にレコードの番号を足す（エンジン側の変更、次のエンジンのセッションで決める）。
- 定義行の無い最初のレコードの offset は 0/0（`data-service.ts:365` の `header_offset < sequence_offset` の検査）。種類 1・2 の一様な配置の `eol` は 1・2 以外の値も取る。注釈だけの入力はレコード 0。`FastaParserKind = 0 | 1`（`dataset.ts`）に 2 を足す。
- 両方の入力に不備があるとき CLI と同じ最初の誤りを出すには、subject を先に `register` する（D-1）。
- CI の firefox・webkit の E2E の 5 秒の待ちの timeout（`smoke.spec.ts:10`、`:37`、`search.spec.ts:391`）が繰り返し出る。
- ABI v1 の host には、空・空白だけ・先頭が空白の定義行で `unknown` の代わりに `Query_1`・`Subject_1`・`unnamed`・最初の語が返る（D-2）。subject を `records: []` で登録したときは、run の `Empty CBlastQueryVector` の誤りとして扱う。
- `screens.spec.ts` の検索画面の 2 つの試験（`LOSAT_WEB_SCREENS` を付けたときだけ流れる、CI には無い）は、エンジンの build で `query-source-0-exclude-refused` を待って止まる（10 分の timeout）。`register` の拒否の文言が行を名指すため（上の 1 つ目）。結果の画面の試験は通る。
- S13 の `results.spec.ts` の「失敗した run」は、TBLASTX の題の HTML の文字参照（E2h で NCBI と同じく復号されて完了する）から、中身の無い subject だけの TBLASTX（NCBI の `The average subject length is too short`）に変えた（tests only）。

## SFd（2026-10-10〜、エンジン側）で行ったこと

修正の再監査（4 観点、監査したコミット `2f954c66`）の指摘を直した。

| 修正 | 内容 | コミット |
|---|---|---|
| 再監査の修正（registry の層、`ncbi_environment.rs`） | R1（A-1 = B-3）：`<prog>.ini` の空の `DATA_LOADERS=` は項目（loader 無効）。R2（A-2）：`.ncbirc` の空の値は、`BLAST_USAGE_REPORT` が偽か `<prog>.ini` があるとき項目（アプリが file を自分で読む。既定は usage report の cache の写しで落ちる）。R3（A-3 = B-4）：空の `NCBI_CONFIG_PATH=` は探索 path を空にし、registry の file を読まない。R4（A-4 + B-2）：registry の file を `IRWRegistry::x_Read` のとおりに byte で読む（C locale の白空白、名前の文字、全ての値の引用符と `ParseEscapes`、継続行、UTF-8 の BOM、section より前の項目は捨てる）。NCBI が構文の誤りとする file と UTF-16 の file は明示的に拒否（`AUTHORITY.md` §J-8）。usage report だけが読む `.ncbirc` の構文の誤りと Boolean でない `BLAST_USAGE_REPORT` も拒否。`AUTHORITY.md` §G1・§J・§L1・§L2-1 と判断 5 の理由を直した。再現 398 件（再監査 A・B の環境と registry の全件と新しい件、program の directory の 4 件を含む。NCBI は `unshare -rn`）：修正の前は差 145、後は差 0（same 198、loader が有効で Seq-id の行を §J-1 で拒否 106、構文の誤りの拒否 85、その他の記録済みの拒否 9：`-conffile`、`[NCBI] .Inherits`、`[NCBI] DATA_LOADERS`、Boolean でない `BLAST_USAGE_REPORT` 4（NCBI は abort）、UTF-16 2）。fixture に 38 行（`*.seqid0_registry_*`、`*.registry_cfgpath_empty.o6`、registry の file は `LOSAT/tests/fixtures/fasta_input/registry/`。修正の前は 38 行とも拒否） | `f5729e60` |
| 再監査の修正（ABI v1、第 1 回の D-1・D-2） | D-1：v1 の BLASTP（`run_pair` と FASTA の handle）と TBLASTX は、`bio` が ID を空にし `from_bio` の record が NCBI の読み込み器と違う定義行（先頭か白空白だけの中の非 ASCII の白空白、行の途中の単独の CR、行頭の白空白の後の制御文字、末尾の非 ASCII の白空白）を検索の後で明示的に拒否する（`web_api/v1_bio.rs` の検査、第 1 回の名前は `check_empty_id_deflines_of`、今は `check_bio_deflines_of`、判断 12・26）。単体の試験は NCBI の読み込み器（`fasta_reader`）との全数の照合（12 文字からなる 4 文字までの定義行、蛋白と核酸）を含む。監査 D の格子を監査したコミットと比べた：D-1 の 285 run は全て拒否、残りは同じ（例外は判断 12 の BLASTP の outfmt 0 の query の 9 run と、fuzz（seed 7）で D-1 と同じ種類の 135 run。どちらも NCBI と違う出力だった）。D-2：判断 12 の白空白を C の `isspace`（空白・tab・VT・FF）と書き直した。`web_api` は wasm32 だけなので native の実行ファイルは変わらない（sha256 `54bab76e…0a4f`）。第 1 回で未解決とした 2 種類は第 2 回で直した | `844de3d8` |
| 再監査の第 2 回の修正（registry の層、R-1〜R-9、`ncbi_environment.rs`） | R-1：program の directory は探索 path を作る時点で決まる（usage report が `.ncbirc` を読まないとき、つまり `BLAST_USAGE_REPORT` が偽か `NCBI_DONT_USE_NCBIRC`・`NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` があるときは、`argv[0]` からの名前の directory が link を解いた directory の前に入る。symlink から起動したとき）。R-2：`HOME` が無いときは passwd の home（passwd に利用者の項目も無いときは、探索が home に届けば拒否）。R-3：`NCBI_CONFIG_PATH` は byte で分ける。R-4：`<prog>.ini` は BOM を 2 回検査する（`EF`・`FE`・`FF` + `BB BF`。ちょうど 2 byte の `EF|FE|FF BB` は空）。R-5：1 byte の `EF`・`FE`・`FF` の file は空で、NCBI の `Error: (110.4) Error reading the registry after line 1: ` を同じ回数出す。R-6：開けない file は registry なし（探索はそこで止まる）。R-7：空の `NCBI_CONFIG_OVERRIDES=` は無いのと同じ。R-8：`AUTHORITY.md` §K-13（NCBI の parameter の環境変数の解析できない値で NCBI は異常終了、LOSAT は普通に走る）。R-9：`AUTHORITY.md` §G1・§J-8・§K・§L1・§L2-1、判断 5・25・27、`registry_search_path` の説明を直した。再現 477 件（第 2 回の監査の 117 件、新しい 226 件、9 次元の条件の行列 134 件。NCBI は `unshare -rn`、passwd は `unshare -rnm`、開けない file は両方 `unshare -cn`）：差 0（same 276、§J-1 の Seq-id の拒否 133、構文の誤りの拒否 45、§K-13 の記録 13、home の拒否 4、UTF-16 の拒否 3、その他の記録済みの拒否 3：`.ncbirc` の Boolean でない `BLAST_USAGE_REPORT`（項目の許可表）2、`-version`）、model の判定 121/121。第 1 回の再現 398 件は判定が全て同じ | `68cc3ba5` |
| 再監査の修正（ABI v1、第 2 回） | 第 1 回の残り 2 種類を v1 の拒否にした（判断 12 の (2)・(3)、判断 28）：`>?` の行の後の題の無い record の番号が出るところ（TBLASTX、BLASTP の query の ID の field のある tabular）と、BLASTP の outfmt 0 の subject で `bio` の読み方が NCBI と違い題が非 ASCII で終わるもの（原因は共有の層の符号付き `char` の `x_CleanAndCompress`、CLI も使うので案 (b)）。題の無いまま残る record は名前の出るところだけ拒否する（判断 26 の 9 run は S11 の byte に戻った）。単体の試験は `>?` の規則と NCBI の読み込み器の局所 ID の全数の照合（9 種の定義行の 4 record まで）を含む。監査 D の格子：D-1 の 285 run だけが拒否、他は監査したコミットと同じ（fuzz seed 7 は 138 が拒否、どれも NCBI と違った）。差分の生成器 24000 run：S11 と同じ 10131、拒否 8429（生成器の形式では S11 と同じだった run は 0。hit の無い record と、形式が出さない役割の拒否は第 3 回の監査 V の V-3、判断 30）、NCBI と一致 781、決着しない 2 種類 4100・559（判断 28 の (4)）。native の実行ファイルは変わらない | `a5d3421b` |
| 再監査の第 3 回の修正（registry の層、S-1〜S-5、`ncbi_environment.rs`） | S-1：同じ名前の環境変数は NCBI の読み手ごとに最初（`getenv`）か最後（`CNcbiEnvironment`：`NCBI_CONFIG__…`、program を探す `PATH`）の項目を使い、値で決まる拒否もその値だけを見る。S-5：相対の `argv[0]` の末尾の `/` を `CFile` の検査で落とす。S-2：`<prog>.ini` の名前は NCBI の program が自分の名前で起動された場合のもの（範囲を §G1・§J-8 に記録）、空の `argv[0]` は拒否。S-3：home の拒否の文言と §J-8 を本当の条件に（`std::env::home_dir` が項目を返さない：項目が無いか `getpwuid_r` の buffer に収まらない）。S-4：§K-14（`=` の無い環境の項目で NCBI は SIGSEGV、LOSAT は普通に走る）。監査 S の harness（`execve`）で S の 161 件と近くの 32 件：S-1・S-5 の全件が same か §J-1 の Seq-id の拒否、S-2 は空の `argv[0]` の拒否 9 件（`argv[0]` が `foo` の 1 件は範囲の外）、S-3 は home の拒否 2 件、S-4 は §K-14 の 6 件、そのほかは記録済みの拒否。第 1 回の 398 件・第 2 回の 477 件の判定は変わらない（home の拒否の 4 件は文言だけ） | `201651b8` |
| 再監査の修正（ABI v1、第 3 回、監査 V の V-1〜V-4） | V-1：`bio` が止まる record より後の定義行を調べない。V-2：定義行の途中の単独の CR より後の record で byte が変わるもの（ID が空のもの、BLASTP の outfmt 0 の非 ASCII で終わる題）を拒否する（NCBI の行の読み手の行の終わりの型に依るので、行ごとの模型では決まらない。判断 12 の (3)）。V-3：その役割の ID か題を形式が出すときだけ拒否する（hit の無い record は拒否する）。V-4：名前と注釈。単体の試験は `bio` の停止と、番号の規則と NCBI の読み込み器の局所 ID の照合（12 種の定義行の 4 record まで。CR の無い入力は一致、CR のある入力は NCBI の番号が違えば必ず拒否）を含む。監査 V の格子 9457 run：S11 と同じ 3811、拒否 3994（S11 と同じ byte だった 74 はどれも hit の無い record）、NCBI と一致 1488、record ごとに NCBI の値 156、NCBI が拒む 7、既存の句読点の拒否 1。監査 D の格子は第 2 回と同じ。差分の生成器 24000 run：S11 と同じ 10131、拒否 8745（S11 と同じだった run は 0）、NCBI と一致 777、決着しない 2 種類 3992・355（判断 28 の (4)）。native の実行ファイルは変わらない | `2e82c250` |

## 推奨の案で進めた判断（保守者に委ねられた判断、2026-09-29 の常設の指示、10-07 に再掲）

1. 移植の順は `PORT_PLAN.md` の S0〜S10。`AUTHORITY.md` を S1 の前に書いた。
2. **ABI v1**：v1 は `bio` と今の検査のまま、`FastaRecord::from_bio` で `run_local` に入る（v1 の出力と文言は変わらない）。保守者の判断 4 の目的（TD-1 の凍結）を保つ。判断 4 の前提「v1 が受け付ける入力は bio と NCBI で読み方が同じ」は TBLASTX と BLASTP の v1 で成り立たない（棚卸し AD-24・AD-25）。
3. **Seq-id の行**：`CSeq_id` の解析（BI-15〜19・55）を移した。広めに拒否するのは accession の guide の形（文字と数字の数の 26 の形、先頭の byte が `A-Z`・`_`・`?` のときだけ）。guide の表（1458 の規則）は移さない。`ZZ123456` などは NCBI では FASTA だが LOSAT は拒否する（判断 2 の範囲）。
4. Web の `register` は data loader を有効とし（CLI の既定）、読み込みの誤りと注釈だけの query を NCBI の文言で早く拒否する。
5. 空の `DATA_LOADERS=`（空白だけ・`""` を含む）は、その項目を持つ層では項目があり、両方の data loader を無効にする：環境変数 `NCBI_CONFIG__BLAST__DATA_LOADERS=`（`env_reg.cpp:157-167`。SFc の監査 B-1）、`<prog>.ini`（`CCompoundRegistry::FindByContents` が `fCountCleared` を付けて尋ねる、`ncbireg.cpp:1235-1246`・`:984-991`。SFd の再監査 A-1・B-3）、`.ncbirc` は `BLAST_USAGE_REPORT` が偽の Boolean のときか `<prog>.ini` があるとき（開けない `<prog>.ini` は無いのと同じ、空や 1 byte の `<prog>.ini` は在る。第 2 回の再監査 R-5・R-6）（アプリが file を自分で読む。既定ではアプリは usage report が cache に読んだ registry を `Write` と `Read` で写し、`Write` が空の値を落とす、`metareg.cpp:152-171`・`ncbireg.cpp:226-234`。SFd の再監査 A-2。`AUTHORITY.md` §G1）。`blastdb` も `genbank` も含まない空でない値は両方を無効にする。以前は「registry の file の空の値は項目なし（`ncbireg.cpp:984-991`）」としていた（その行は `fCountCleared` のある場合を読み落としていた。SFd で訂正）。
6. 読み込み器の API：Seq-id の拒否の後は、拒否した 1 行の次から読める（局所 ID の番号は進まない）。program は拒否で止まる。
7. `from_bytes` は同じ byte の通常の file と同じに読む（NCBI の stream の補充の大きさによる尾の消失を含む）。
8. まとめて読む経路は `FastaStream::bulk` の旗の後ろに置く（program では常に有効。試験では 1 byte ずつの参照の経路と比べる）。
9. fixture の manifest は情報を落とさず小さくした（全長の SHA-256、行の定義は script の生成器、program ごとに 4 file）。CI の速い検査への接続は移植の後（SFb）。TBLASTX と BLASTP の `-lcase_masking` の 34 行は残し、移植か明示的な拒否かを SFb で決める。
10. `PORT_PLAN.md` の Q4〜Q6・Q8・Q10・Q11 は推奨のとおり。BI の覚え書き §10 の推奨（query は batch ごとに読む、subject にも同じ Seq-id の規則、空のレコードと空の batch を 4 program で移植、空の pipe の零 batch、`-parse_deflines` は拒否のまま）も同じ。

11. （SFb）TBLASTX と BLASTP の `-lcase_masking` の fixture の 34 行は、E2e の `AUTHORITY.md` §A が LOSAT に無い option として拒否を記録しているので、明示的な拒否として期待を直す（移植しない）。
12. （SFb）ABI v1 は今の検査を前と同じ位置で行い、`from_bio` で新しい型に入る。共有の報告の層が NCBI の移植になったので、v1 の BLASTP が受け付け CLI が拒否する 2 種類の入力（非 ASCII の byte で終わる subject の題、行頭に空白のある query の定義行）の出力は NCBI の byte に変わる。計画 TD-1（エンジンを NCBI に合わせる修正は v1 の検索結果にも及ぶ）の範囲として受け入れた。v1 の引数・形式・受け付ける入力・誤りの扱いは変わらない。SFc の監査 D-2 で、同じ原因で変わる入力が他にもあると分かった：BLASTP の空・C の `isspace` の白空白（空白・tab・VT・FF。CR だけの `>\r\r` も）だけの定義行（query は outfmt 6・7 の行が `unknown` から `Query_1`、subject は outfmt 0・6・7 が `unnamed protein product`・`unnamed`）、行頭に C の白空白（VT・FF を含む。その後が非 ASCII の白空白でもよい）のある定義行の outfmt 6・7 の行（query・subject、`unknown` から最初の語）と subject の outfmt 0、TBLASTX の outfmt 6 の同じ定義行（`Query_1`・`Subject_1`・最初の語）。どれも NCBI BLAST+ 2.17.0 の CLI の出力と byte で一致する（44 実行、`unshare -rn`。一覧は task folder の `logs/fix/v1-changes.md`。SFd の再監査 D で FF だけ・行頭の VT/FF の定義行、空白の後の非 ASCII の白空白、`>\r\r` も一致）。BLASTN の v1 と TBLASTX の v1 の outfmt 0・7 は前と同じ拒否。TBLASTX の v1 は平均の subject の長さの誤りに届かない（中身の無い subject を先に拒否する）。SFd の再監査 D-1（第 1 回）とその後の検証（第 2 回、第 3 回の監査 V）で、共有の報告の層が NCBI の移植になったため v1 の出力が S11 の byte とも NCBI の byte とも違う種類が分かった。v1 の BLASTP（`run_pair` と FASTA の handle）と TBLASTX は次の定義行を明示的に拒否する（`web_api/v1_bio.rs` の `check_bio_deflines_of`。拒否は v1 の他の検査と検索の後で、他の誤りとその順は変わらない。query の入力に record があるときだけ。調べるのは `bio` が読む record だけで、`bio` が止まる record（Unicode の白空白を除くと空の定義行で、次の `>` の行か終わりまでに残基の行が無いもの、bio-1.6.0 `src/io/fasta.rs:1034`）より後は調べない）。どの規則も、その役割について規則が変える物を形式が出すときだけ働く：query は ID（`std`・`qseqid`・`qacc`・`qaccver` のある tabular）と題（outfmt 0・7 の query の行）、subject は ID（outfmt 0 と、`std`・`sseqid`・`sacc`・`saccver` のある tabular）と題（outfmt 0 と `stitle`）、TBLASTX は両方の ID。検索が報告しない record（hit の無い query と subject）も形式だけで判定して拒否する（その record は出力に出ないので、S11 と同じ byte の run も拒否になる）。(1) `bio` が ID を空にする定義行（`char::is_whitespace` の文字で始まるか、それだけの定義行）のうち、`from_bio` の record が NCBI の読み込み器の record と違うもの：先頭か白空白だけの中の非 ASCII の白空白 U+0085・U+00A0・U+1680・U+2000〜200A・U+2028・U+2029・U+202F・U+205F・U+3000 を NCBI は題に残す、行の途中の単独の CR で NCBI の行が終わる、行頭の白空白の後の制御文字で NCBI の題が終わる、行頭の白空白の後の題の末尾の非 ASCII の白空白を NCBI は残す（S11 は `unknown`、拒否の前は最初の語・`Query_N`・`Subject_N`・`unnamed`）。題のある record は ID か題が出るとき、題の無いまま残る record は ID（名前）が出るときだけ（BLASTP の outfmt 0・7 の query の行は空の題だけなので S11 の byte のまま）。文言は `<role> record <n> has a defline that starts with white space and <理由>; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLASTP`（または `TBLASTX`）。(2) `>?` の行（`>?_?` を含み、最初の行を除く。NCBI は前の record の gap として読み、record にしない）の後の題の無い record：NCBI の番号が `bio` より小さい（TBLASTX の outfmt 6 で S11 `unknown`、拒否の前 `Subject_3`、NCBI `Subject_2`）。局所 ID の出るところだけ：TBLASTX の両方の役割と、query の ID の field のある BLASTP の tabular の query（BLASTP の題の無い subject は `unnamed`）。文言は `… has a defline without a title after a line that starts with '>?', which NCBI BLAST+ reads as a gap, not as a record (NCBI BLAST+ numbers this record <k>); …`。(3) 定義行の途中の単独の CR（その行の最後の byte でない CR）より後の record で byte が変わるもの（ID が空のもの、BLASTP の outfmt 0 の subject で題が非 ASCII で終わるもの）を、それが出るところで拒否する。NCBI の行の読み手は、CR の後の残りを押し戻すときに読んだ行の終わりを落とし（line_reader.cpp:249-258）、最初の行が CR で終わる file では次の CR まで読むので、後の行を別の record、定義行の一部、残基として読むことがある（query `>s1 a\r>\n<P>\n>\n<P>` の 3 番目の record を S11 は `unknown`、拒否の前は `Query_2`、NCBI は `Query_3` と書く。`>s1\r>?5\n…\n>s1 a\r>\n…` では NCBI は `>s1 a>` を 1 行に読む）。行ごとの模型では NCBI の読み方が決まらないので、NCBI の番号や読み方が `bio` と同じになる入力も拒否する。文言は `<role> record <n> comes after a defline with a carriage return before the end of its line, after which NCBI BLAST+ reads the lines otherwise (as other records or as residues); …`。(4) BLASTP の outfmt 0 の subject で、ID のある定義行の `from_bio` の record が NCBI と違い（tab などの制御文字や非 ASCII の白空白で NCBI の題が違う、行の途中の CR、`>?` の gap の行、NCBI が落とす `?_` の前置き）、題（末尾の `.,;~ ` を除く）が非 ASCII の文字で終わるもの。原因は共有の層の `x_CleanAndCompress` の移植（create_defline.cpp:304-306 の `curr > 0` は符号付きの `char` で、非 ASCII の最後の byte を落とす。S11 の層は落とさなかった）で、CLI も同じ層を使う（判断 28）。例 `>s1\tx é`：S11 `s1 x é`、拒否の前 `s1 x \xc3`、NCBI `s1`。文言は `subject record <n> has a defline that NCBI BLAST+ reads with another title (it <理由>) and that ends with a non-ASCII character, whose last byte NCBI BLAST+'s outfmt 0 drops; this is not supported by LOSAT's BLASTP outfmt 0 (…)`。検証：監査 D の格子は D-1 の 285 run だけが拒否に変わり、他は監査したコミット `2f954c66` と同じ（題 2460、Unicode の白空白 756 のうち 495、CR 50 のうち 26、誤りの行列 1000、fuzz seed 11 4500。fuzz seed 7 3600 は 3462 が同じで 138 が拒否、どれも S11 から変わり NCBI と違った run）。第 3 回の監査 V の格子 9457 run（S11、`2f954c66`、修正後。NCBI は `unshare -rn`）：S11 と同じ 3811、拒否 3994（そのうち S11 と同じ byte だった 74 はどれも hit の無い record：subject 71、query 3）、NCBI と一致 1488、変わった部分が record ごとに NCBI の値 156、NCBI が入力を拒む 7、既存の句読点の拒否 1（`2f954c66` から同じ。`bio` の題が句読点だけになる (1) の定義行）。tab・制御文字・非 ASCII（UTF-8 の途中の切り詰め、幅 60・65・68 の近くの長い題）・番号のある/無い `>?` の行・空の ID・CRLF に偏らせた生成器の差分（6 seed × 4000 run、S11 の reactor・監査したコミット・修正後、NCBI 4173 実行）：24000 run のうち S11 と同じ 10131、拒否 8745（S11 と同じだった run は 0、標本 300 のうち 1 は監査したコミットの出力が NCBI と同じ：(3) の広めの拒否）、NCBI と一致 777。残りの 2 種類は run 全体では S11 の byte でも NCBI の byte でもない（保守者への問い、判断 28 の (4)）：変わった部分は NCBI の値で、残りは S11 の凍結の byte の run 3992（同じ入力に、`bio` と NCBI で読み方の違う別の record、例えば `>?` の gap の行や tab のある題、がある）と、NCBI が `Near line N, there's a line that doesn't look like plausible data` で失敗する入力の run 355（v1 は S11 から受け付けている）。
13. （SFb）Seq-id の行の拒否の文言は `AUTHORITY.md` §G4、終了コード 1。RP-20 は tabular の形式でその subject に HSP があるときだけ、出力の前に拒否する。`-` は NCBI の `cin`（`in_avail` 0）と同じに読む。
14. （SFb）中身の無い query は NCBI の配置のまま batch に残す。満ちた batch の直後の eEOF は空の次の batch として `Empty CBlastQueryVector` で止まる。平均の subject の長さの検査は有効な query のある batch だけ。
15. （SFb）各段階の細部（API の形、検査の順、共有する型の置き場所）は `DECISIONS.md`（`/home/kawato/losat-baselines/sfb-e2h-20261008/`）に段階ごとに記録した。

16. （S7）BLASTP の CLI の順：`Query is Empty!` → BATCH_SIZE → 中身の無い subject の警告 → XInclude → v1 の検査 → LOSAT の制限 → 環境 → prolog → batch。検索は NCBI が終える batch の query をまとめて 1 回。outfmt 0 の句読点だけの subject の題は、承認済みの例外 2 が BLASTN・TBLASTX・TBLASTN だけなので、BLASTP では明示的な拒否のまま（保守者への問い）。ABI v1 の BLASTP の題は NCBI の定義行の読み方と同じく先頭の空白を落とす。
17. （S8）BLASTP の batch の読み方を TBLASTN と共有。文字の無い TBLASTN の query は NCBI と同じく無効として "Sequence contains no data"。LOSAT 独自の空の `-query_loc` 区間の拒否は、NCBI に無いので消した。`from_bio` は全 program で題の先頭の空白を落とす（v1 の TBLASTX の outfmt 6 の qseqid が NCBI と同じになる）。
18. （S9）Seq-id の最初の行は scan でも拒否（エンジンの読み込み器で判定）。read-forward の規則で追えない行の連結は、一様でないレコードだけ明示的な拒否。定義行の無い最初のレコードの offset は 0/0。注釈だけの入力はレコード 0。
19. （S10）`register` は種類 1・2 の scan を先に流す（Web だけの `>?` の拒否）。読み込みの誤りは CLI の文言（`BLAST query error: …`）、Seq-id と `>?` は LOSAT Web の文言。注釈だけの query は `BLAST engine error: Empty CBlastQueryVector`、中身の無い subject は登録できる。`q_idx` はアダプタが `-query_loc` から写す。v1 は `V1Checks` の hook で検査の順を保つ。V-ABI は実行ファイルの名前が `LOSAT` であること。
20. （fixture）拒否のハッシュは同じ manifest に固定。`rejection` は終了コードと 3 つのハッシュが一致するときだけ。
21. （ゲート）TBLASTN の Stage G は `unshare -rm` の中で `/mnt` に tmpfs を置き、path を作って worktree を bind する（`/mnt/c` を読まない）。
22. （アプリの試験）PR #118 の CI の `browser-engine` で失敗したアプリの試験は、エンジンの契約の変更に合わせて別のコミット（tests only）で直した。ゲートの HEAD を動かさないよう、一時的な worktree から push した。
23. （監査 A-2）accession の guide の表（`accguide2.inc`、1.95 MB、25,779 行）は移さず、記録済みの広めの拒否のままにする（約 1 MB を超える file の規則。保守者の判断 2 の「忠実に移せない部分は広めに拒否」）。
24. （監査の修正）B-1 は環境変数の層だけ（file の空の値は項目なしのまま。この file の部分は SFd の再監査で誤りと分かり、25 で直した）、D-1 は文書だけ（ABI の呼び方を変えない）、D-2 は v1 を変えない（変わった 44 run は NCBI の byte、`/home/kawato/losat-baselines/sfc-e2h-20261009/logs/fix/v1-changes.md`）。

25. （SFd の再監査の修正）R1〜R3 は NCBI の規則をそのまま移した（R2 の条件「アプリが `.ncbirc` を自分で読むか」は、環境変数 `BLAST_USAGE_REPORT` の偽の Boolean と `<prog>.ini` の有無で全部決まり、拒否に回す部分は無い）。R4 は registry の file を `IRWRegistry::x_Read` のとおりに byte で読み、NCBI が構文の誤りとする file は拒否する（`.ncbirc` の文言は NCBI の build の source の path と行を含み、`<prog>.ini` の文言は診断の層の例外の報告（`NCBI_REPORT_EXCEPTION_X`、この build では 2 回）で、どちらも LOSAT は移さない。終了コードは他の registry の拒否と同じ 1）。UTF-16 の file も拒否（`ReadIntoUtf8` を移さない）、UTF-8 の BOM と section より前の項目は NCBI と同じく扱う（BOM の検査は `IRWRegistry::Read` の `GetTextEncodingForm` で、`.ncbirc` は 1 回、`<prog>.ini` は 2 回。この修正は両方 1 回としていた。第 2 回の再監査 R-4 で直した、27）。`<prog>.ini` の `[NCBI] DONT_USE_NCBIRC` でアプリが読まない `.ncbirc` も usage report が読むので、構文の誤りと `BLAST_USAGE_REPORT` だけ調べる。fixture の行は registry の file を `LOSAT/tests/fixtures/fasta_input/registry/` に置き、相対の `NCBI_CONFIG_PATH`（R3 は空の値と `HOME`）で指す。構文の誤りの file は fixture にしない（NCBI の stderr が build の path を含む。証拠は task folder の再現）。
26. （SFd の再監査の修正、ABI v1）D-1 は v1 の BLASTP・TBLASTX で明示的に拒否する（推奨の案。NCBI の byte を出す案は `bio` の分割の横に生の定義行が要り、v1 の変更が大きいので採らない）。拒否する範囲は、`bio` が ID を空にする定義行のうち `from_bio` の record（題、CR のときは残基も）が NCBI の読み込み器の record と違うもの全部（D の 2 種類に加え、行頭の白空白の後の制御文字と末尾の非 ASCII の白空白）。第 3 回から、その役割の ID か題を形式が出すときだけ拒否する（題の無いまま残る record は ID が出るときだけ。S11 と同じ byte だった BLASTP の outfmt 0 の query の 9 run は第 2 回で S11 の byte に戻った。判断 28・30）。検索が報告しない record（hit の無い query と subject）は出力から分からないので、形式だけで判定して拒否する。拒否は v1 の検索の後に置き、v1 の他の誤りとその順を変えない（BLASTN の v1 の検査は検索の前だが、その位置では他の誤りより先に出る入力がある）。query の入力に record が無いときは調べない（NCBI の `Query is Empty!`、BLASTN の v1 と同じ）。
27. （SFd の再監査の第 2 回の修正、registry の層、R-1〜R-9）推奨の案のとおり、NCBI の規則がはっきりしているものは移した：R-1 program の directory は探索 path を作る時点で決まる（`<prog>.ini` の名前は program 名のまま：29。usage report が `.ncbirc` を読まないとき、つまり環境変数 `BLAST_USAGE_REPORT` が偽か `NCBI_DONT_USE_NCBIRC`・`NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` があるときは、`FindProgramExecutablePath` の名前（`argv[0]`：絶対、作業 directory から、`PATH` の base 名）の directory が先、link を解いた directory が次。それ以外は Linux で解いた directory だけ。macOS・Windows の分は source からで、oracle では確かめていない）。R-2 `HOME` が無いときは passwd の利用者の項目の home（Rust の `std::env::home_dir` は `HOME` が無いときだけ呼ぶ。`HOME=` の空の値を `getpwuid_r` に回すため）。passwd に項目が無いときの login 名の `getpwnam` は、`libc` への依存を足さないと移せないので、探索がその home に届くときだけ拒否する（§J-8。項目が `getpwuid_r` の buffer に収まらないときも同じ拒否になる：29）。R-3 `NCBI_CONFIG_PATH` は byte で分ける。R-4 `<prog>.ini` は BOM を 2 回検査（`EF`・`FE`・`FF` + `BB BF` を UTF-8 の BOM とする、ちょうど 2 byte の `EF|FE|FF BB` の `<prog>.ini` は 2 回目の検査で stream が失敗したまま空になることも移した）。UTF-8 の BOM の後の UTF-16 は UTF-16 の拒否のまま。R-5 1 byte の file の `Error: (110.4) Error reading the registry after line 1: ` は `.ncbirc` も含め再現する（文言が固定で path を含まず、回数は usage report とアプリのどちらが読むかで決まり、usage report・`BLAST_USAGE_REPORT` 偽・`<prog>.ini` あり・`[NCBI] DONT_USE_NCBIRC` の組み合わせを oracle で確かめた）。R-6 開けない file は registry なしで探索も止まる（開けたが読めない file は拒否）。R-7 空の `NCBI_CONFIG_OVERRIDES=` は無いのと同じ。R-8 は FASTA の読み方と関係しない NCBI の parameter の解析の異常終了なので `AUTHORITY.md` §K-13 に記録し、LOSAT は普通に走る（有効な値は出力も R2 の条件も変えないことを oracle で確かめた）。R-9 は文書（`AUTHORITY.md` §G1・§J-8・§K・§L1・§L2-1、判断 5・25、`registry_search_path` の説明）。再現は task folder の `scripts/registry_repros3.py`（第 2 回の監査の 117 件、新しい 226 件、生成した 9 次元の条件の行列 134 件、差 0。NCBI は `unshare -rn`、passwd は `unshare -rnm`、開けない file は両方 `unshare -cn`）。
28. （SFd の再監査の修正、ABI v1 の第 2 回）(1) `>?` の行の後の題の無い record は、番号の出るところだけ拒否する（調整役の推奨のとおり最も狭く）。(2) outfmt 0 の最後の byte の脱落は共有の報告の層（CLI も使う NCBI の `x_CleanAndCompress` の移植）に原因があるので、案 (b)：v1 の狭い拒否（BLASTP の outfmt 0 の subject で、ID のある定義行の record が NCBI と違い、題が非 ASCII で終わるもの）。案 (a) は採れない：S11 の byte には共有の層に v1 用の動きが要り、NCBI の byte には v1 が自分では読まない題を渡すことになり、CR の分割では NCBI の record にならず、出力の前の HtmlDecode の検査も別の題で行われて他の v1 の誤りが変わる。(3) 題の無いまま残る record は名前の出るところだけ拒否する（判断 26 の 9 run を S11 の byte に戻す）。(4) 差分の検証で決着しない 2 種類（変わった部分は NCBI の値で残りは S11 の凍結の byte の run、NCBI が失敗する入力の run）は、判断 2（v1 は `bio` で読む入力を凍結して受け付ける）を変えずには決着しない。BLASTN の v1 のように `bio` と NCBI で読み方の違う定義行を全部拒否すれば決着するが、S11 と同じ byte の run も拒否になる。推奨は今のまま（変わる byte は record ごとに NCBI の byte、残りは S11 の凍結の byte）で、保守者に問う。
29. （SFd の再監査の第 3 回の修正、registry の層、S-1〜S-5）推奨の案のとおり：S-1 同じ名前の環境変数が複数あるとき、NCBI の `getenv` の読み手（`BLAST_USAGE_REPORT`・`NCBI_CONFIG_OVERRIDES`・`NCBI_CONFIG_PATH`・`NCBI_DONT_USE_LOCAL_CONFIG`・`NCBI_DONT_USE_NCBIRC`・`NCBI`・`HOME`）は最初、アプリの `CNcbiEnvironment`（`NCBI_CONFIG__…` の registry の層と `FindProgramExecutablePath` の `PATH`）は最後の項目を使うので、LOSAT も読み手ごとにそうし、値で決まる拒否は NCBI が使う値だけを見る（LOSAT の他の環境変数 `BATCH_SIZE`・`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`・`BL2SEQ_LEGACY`・`PRE_FETCH_SEQS_LIMIT`・`OLD_FSC`・`ADAPTIVE_CBS` は NCBI も `getenv` で読み、LOSAT の `std::env::var_os` も最初の項目なので変えない）。S-5 相対の `argv[0]` の末尾の `/` は NCBI の `CFile` の検査（`CDirEntry::Reset`）のとおり落とす。S-2 `<prog>.ini` の名前は NCBI の program が自分の名前で起動された場合のものとし（LOSAT の `argv[0]` は program の directory だけ）、別の名前で起動した NCBI は範囲の外として `AUTHORITY.md` §G1・§J-8 に記録した。LOSAT の `argv[0]` が空のときは写す起動が無いので明示的に拒否する（NCBI は `ncbi` と名乗り `ncbi.ini` を読む）。S-3 `std::env::home_dir` は `getpwuid_r` を 1 回だけ呼び `ERANGE` で広げ直さず、std に別の方法も無いので、拒否の文言と §J-8 を本当の条件（home を決められない：項目が無いか buffer に収まらない）に直した（`libc` への依存は足さない）。S-4 `=` の無い環境の項目で NCBI が SIGSEGV で落ちるのは NCBI の不具合として §K-14 に記録し、LOSAT は普通に走る。再現は task folder の `scripts/s_harness/`（第 3 回の監査 S の harness の写し：`execve` で `argv[0]` と envp を正確に渡す）で、S の全件（`/proc/self/mem` の 3 件を除く 161 件）と近くの 32 件。
30. （SFd の再監査の修正、ABI v1 の第 3 回、監査 V の V-1〜V-4）(1) V-1：`bio` が止まる record（bio-1.6.0 `src/io/fasta.rs:1034` の `Records::next` は、ID も説明も残基も無い record で終わる）より後の定義行は調べない（v1 はその record を読まない）。(2) V-2：定義行の途中の単独の CR の後で NCBI が何を record として読むかは、NCBI の行の読み手の行の終わりの型と押し戻し（`x_AdvanceEOLUnknown`、`x_AdvanceEOLSimple` が落とす行の終わり、line_reader.cpp:219-258）に依り、行ごとの模型では正確に決まらない（反例 `>s1\r>?5\n…\n>s1 a\r>\n…` では NCBI は `>s1 a>` を 1 行に読む）。調整役の推奨の後半（正確に書けなければその狭い種類を拒否）に従い、そのような CR より後の record で byte が変わるもの（ID が空のもの、BLASTP の outfmt 0 の非 ASCII で終わる題）を、出るところで拒否する（判断 12 の (3)）。NCBI と同じ読み方になる入力も拒否するが、その範囲を文書に書いた。`>?` の行の番号の規則は正確なまま。(3) V-3：調整役の推奨（狭める）に従い、その役割の ID か題を形式が出すときだけ拒否する。hit の無い record は v1 の層からは分からないので、形式だけで判定して拒否する（文書に書いた）。(4) V-4：第 1 回の行の関数名（今は `check_bio_deflines_of`）と、`v1_blastp.rs`・`v1_tblastx.rs` の注釈を直した。

## 保守者への一括の問い（`PD-LOSAT-NCBI-DEFECTS`）

`AUTHORITY.md` §K の 14 件。移植は各行の推奨の規則で進める（再現：2〜7・9・11 の切り詰め、拒否：1・8・11 の巨大な gap・12、Linux の意味：10、LOSAT は普通に走る：13・14）。保守者の答えで変わるのは主に次の 3 つ：(1) RP-20 tabular の `Subject_` の語と復号できない題で NCBI が終了コード 255 で落ちる → 明示的な拒否、(4) 改行の無い最後の行の尾が届き方（file と pipe）で失われたり残ったりする → file の挙動を再現、(3) 行の途中の単独の CR が行の終わりを失わせる → 再現（BLASTX の移植と同じ）。

- BLASTP の outfmt 0 の句読点だけの subject の題（NCBI は SIGSEGV）：承認済みの例外 2 を BLASTP にも広げるか（推奨：広げる。今は明示的な拒否）。
- ABI v1 の出力の変化（判断 12 と D-2）：v1 の BLASTP の空・空白だけ・先頭が空白の定義行と非 ASCII の byte で終わる題（outfmt 0/6/7）、v1 の TBLASTX の outfmt 6 の同じ定義行。いずれも NCBI の byte に変わる（TD-1 の「NCBI に合わせる修正は v1 の検索結果にも及ぶ」の範囲として受け入れた）。NCBI の byte にならない定義行（先頭が非 ASCII の白空白、行の途中の単独の CR とその後の record、`>?` の行の後の題の無い record の番号、BLASTP の outfmt 0 の非 ASCII で終わる題で `bio` の読み方が違うもの、SFd の再監査 D-1 と第 2・3 回）は、形式がそれを出すときに拒否に変えた（判断 12・26・28・30。hit の無い record も拒否する）。同じ入力に `bio` と NCBI で読み方の違う別の record があるか、NCBI が入力を拒む run は、変わった部分だけが NCBI の値になる（判断 28 の (4)）。この扱いでよいか。
- `AUTHORITY.md` §K の 14 件（前の節）。13（第 2 回の再監査 R-8）は NCBI の parameter の環境変数の解析できない値で NCBI が異常終了するもの、14（第 3 回の再監査 S-4）は `=` の無い環境の項目で NCBI が SIGSEGV で落ちるもので、どちらも LOSAT は普通に走る。記録に留めるか拒否に加えるか。
- 保守者の指示（2026-10-09「とりあえず main に適宜マージして」）で、検証の済んだ区切りごとに `main` へ PR を出した（#116 `03947fb9`、#117 `e31236bf`、#118 `feea1cea`。計画 DW-20 の「今は PR を作らない」を置き換える、計画 DW-24）。

## 残件（SFd の最初の作業）

- 修正の再監査（4 観点とも supported まで）。
- 最後のコミットでのゲートの全工程（`FRESH=1`）。
- Gate A、V-PERF。
- 終了・引き継ぎ（SF の指示書）、`main` への PR。
- アプリ側の S13 は SFc で merge 済み（`39a563f1`）。
