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

## 推奨の案で進めた判断（保守者に委ねられた判断、2026-09-29 の常設の指示、10-07 に再掲）

1. 移植の順は `PORT_PLAN.md` の S0〜S10。`AUTHORITY.md` を S1 の前に書いた。
2. **ABI v1**：v1 は `bio` と今の検査のまま、`FastaRecord::from_bio` で `run_local` に入る（v1 の出力と文言は変わらない）。保守者の判断 4 の目的（TD-1 の凍結）を保つ。判断 4 の前提「v1 が受け付ける入力は bio と NCBI で読み方が同じ」は TBLASTX と BLASTP の v1 で成り立たない（棚卸し AD-24・AD-25）。
3. **Seq-id の行**：`CSeq_id` の解析（BI-15〜19・55）を移した。広めに拒否するのは accession の guide の形（文字と数字の数の 26 の形、先頭の byte が `A-Z`・`_`・`?` のときだけ）。guide の表（1458 の規則）は移さない。`ZZ123456` などは NCBI では FASTA だが LOSAT は拒否する（判断 2 の範囲）。
4. Web の `register` は data loader を有効とし（CLI の既定）、読み込みの誤りと注釈だけの query を NCBI の文言で早く拒否する。
5. registry の file（`<prog>.ini`・`.ncbirc`）の空の `DATA_LOADERS=` は項目なし（data loader は有効、`ncbireg.cpp:984-991`）。`blastdb` も `genbank` も含まない空でない値は両方を無効にする。環境変数の空の `NCBI_CONFIG__BLAST__DATA_LOADERS=` は項目があり、両方を無効にする（`env_reg.cpp:157-167`。SFc の監査 B-1 で訂正。以前は file の規則を環境変数にも当てはめていた）。
6. 読み込み器の API：Seq-id の拒否の後は、拒否した 1 行の次から読める（局所 ID の番号は進まない）。program は拒否で止まる。
7. `from_bytes` は同じ byte の通常の file と同じに読む（NCBI の stream の補充の大きさによる尾の消失を含む）。
8. まとめて読む経路は `FastaStream::bulk` の旗の後ろに置く（program では常に有効。試験では 1 byte ずつの参照の経路と比べる）。
9. fixture の manifest は情報を落とさず小さくした（全長の SHA-256、行の定義は script の生成器、program ごとに 4 file）。CI の速い検査への接続は移植の後（SFb）。TBLASTX と BLASTP の `-lcase_masking` の 34 行は残し、移植か明示的な拒否かを SFb で決める。
10. `PORT_PLAN.md` の Q4〜Q6・Q8・Q10・Q11 は推奨のとおり。BI の覚え書き §10 の推奨（query は batch ごとに読む、subject にも同じ Seq-id の規則、空のレコードと空の batch を 4 program で移植、空の pipe の零 batch、`-parse_deflines` は拒否のまま）も同じ。

11. （SFb）TBLASTX と BLASTP の `-lcase_masking` の fixture の 34 行は、E2e の `AUTHORITY.md` §A が LOSAT に無い option として拒否を記録しているので、明示的な拒否として期待を直す（移植しない）。
12. （SFb）ABI v1 は今の検査を前と同じ位置で行い、`from_bio` で新しい型に入る。共有の報告の層が NCBI の移植になったので、v1 の BLASTP が受け付け CLI が拒否する 2 種類の入力（非 ASCII の byte で終わる subject の題、行頭に空白のある query の定義行）の出力は NCBI の byte に変わる。計画 TD-1（エンジンを NCBI に合わせる修正は v1 の検索結果にも及ぶ）の範囲として受け入れた。v1 の引数・形式・受け付ける入力・誤りの扱いは変わらない。SFc の監査 D-2 で、同じ原因で変わる入力が他にもあると分かった：BLASTP の空・空白だけ（空白・tab・VT）の定義行（query は outfmt 6・7 の行が `unknown` から `Query_1`、subject は outfmt 0・6・7 が `unnamed protein product`・`unnamed`）、行頭に空白のある定義行の outfmt 6・7 の行（query・subject、`unknown` から最初の語）と subject の outfmt 0、TBLASTX の outfmt 6 の空・空白だけ・行頭に空白のある定義行（`Query_1`・`Subject_1`・最初の語）。どれも NCBI BLAST+ 2.17.0 の CLI の出力と byte で一致する（44 実行、`unshare -rn`。一覧は task folder の `logs/fix/v1-changes.md`）。BLASTN の v1 と TBLASTX の v1 の outfmt 0・7 は前と同じ拒否。TBLASTX の v1 は平均の subject の長さの誤りに届かない（中身の無い subject を先に拒否する）。
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
24. （監査の修正）B-1 は環境変数の層だけ（file の空の値は項目なしのまま）、D-1 は文書だけ（ABI の呼び方を変えない）、D-2 は v1 を変えない（変わった 44 run は NCBI の byte、`/home/kawato/losat-baselines/sfc-e2h-20261009/logs/fix/v1-changes.md`）。

## 保守者への一括の問い（`PD-LOSAT-NCBI-DEFECTS`）

`AUTHORITY.md` §K の 12 件。移植は各行の推奨の規則で進める（再現：2〜7・9・11 の切り詰め、拒否：1・8・11 の巨大な gap・12、Linux の意味：10）。保守者の答えで変わるのは主に次の 3 つ：(1) RP-20 tabular の `Subject_` の語と復号できない題で NCBI が終了コード 255 で落ちる → 明示的な拒否、(4) 改行の無い最後の行の尾が届き方（file と pipe）で失われたり残ったりする → file の挙動を再現、(3) 行の途中の単独の CR が行の終わりを失わせる → 再現（BLASTX の移植と同じ）。

- BLASTP の outfmt 0 の句読点だけの subject の題（NCBI は SIGSEGV）：承認済みの例外 2 を BLASTP にも広げるか（推奨：広げる。今は明示的な拒否）。
- ABI v1 の出力の変化（判断 12 と D-2）：v1 の BLASTP の空・空白だけ・先頭が空白の定義行と非 ASCII の byte で終わる題（outfmt 0/6/7）、v1 の TBLASTX の outfmt 6 の同じ定義行。いずれも NCBI の byte に変わる（TD-1 の「NCBI に合わせる修正は v1 の検索結果にも及ぶ」の範囲として受け入れた）。この扱いでよいか。
- `AUTHORITY.md` §K の 12 件（前の節）。
- 保守者の指示（2026-10-09「とりあえず main に適宜マージして」）で、検証の済んだ区切りごとに `main` へ PR を出した（#116 `03947fb9`、#117 `e31236bf`、#118 `feea1cea`。計画 DW-20 の「今は PR を作らない」を置き換える、計画 DW-24）。

## 残件（SFd の最初の作業）

- 修正の再監査（4 観点とも supported まで）。
- 最後のコミットでのゲートの全工程（`FRESH=1`）。
- Gate A、V-PERF。
- 終了・引き継ぎ（SF の指示書）、`main` への PR。
- アプリ側の S13 は SFc で merge 済み（`39a563f1`）。
