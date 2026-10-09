# E2h（SF・SFb・SFc）ゲート記録：FASTA の読み方（NCBI `CFastaReader` の移植）

状態：**途中**（2026-10-09、SFb を核酸の入力の移植の後で区切った）。完了条件（計画 §7 の SF の行）はまだ満たしていない。続きは [SFc の指示書](../../losat_web_gui_sessions/session_sfc_e2h_protein_gates.md)。

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
12. （SFb）ABI v1 は今の検査を前と同じ位置で行い、`from_bio` で新しい型に入る。共有の報告の層が NCBI の移植になったので、v1 の BLASTP が受け付け CLI が拒否する 2 種類の入力（非 ASCII の byte で終わる subject の題、行頭に空白のある query の定義行）の出力は NCBI の byte に変わる。計画 TD-1（エンジンを NCBI に合わせる修正は v1 の検索結果にも及ぶ）の範囲として受け入れた。v1 の引数・形式・受け付ける入力・誤りの扱いは変わらない。TBLASTX の v1 は平均の subject の長さの誤りに届かない（中身の無い subject を先に拒否する）。
13. （SFb）Seq-id の行の拒否の文言は `AUTHORITY.md` §G4、終了コード 1。RP-20 は tabular の形式でその subject に HSP があるときだけ、出力の前に拒否する。`-` は NCBI の `cin`（`in_avail` 0）と同じに読む。
14. （SFb）中身の無い query は NCBI の配置のまま batch に残す。満ちた batch の直後の eEOF は空の次の batch として `Empty CBlastQueryVector` で止まる。平均の subject の長さの検査は有効な query のある batch だけ。
15. （SFb）各段階の細部（API の形、検査の順、共有する型の置き場所）は `DECISIONS.md`（`/home/kawato/losat-baselines/sfb-e2h-20261008/`）に段階ごとに記録した。

## 保守者への一括の問い（`PD-LOSAT-NCBI-DEFECTS`）

`AUTHORITY.md` §K の 12 件。移植は各行の推奨の規則で進める（再現：2〜7・9・11 の切り詰め、拒否：1・8・11 の巨大な gap・12、Linux の意味：10）。保守者の答えで変わるのは主に次の 3 つ：(1) RP-20 tabular の `Subject_` の語と復号できない題で NCBI が終了コード 255 で落ちる → 明示的な拒否、(4) 改行の無い最後の行の尾が届き方（file と pipe）で失われたり残ったりする → file の挙動を再現、(3) 行の途中の単独の CR が行の終わりを失わせる → 再現（BLASTX の移植と同じ）。

## 残件（SFc の最初の作業）

- 移植 S7（BLASTP）、S8（TBLASTN の query と後片付け、`BLASTINPUT_GEN_DELTA_SEQ` の確かめ直し）、S9（`scan` の種類 1・2 と性質試験）、S10（`register`・`run`・`abi_v2.md`）。
- fixture の仕上げ：Seq-id の 24 行の期待（LOSAT の明示的な拒否の文言と終了コード）、`-lcase_masking` の 34 行の期待、CI の速い検査への接続、変更前と移植の後の表。
- 指示書の手順 6〜9（sweep、ゲートの全工程、V-PERF、独立監査）。ゲートの script の下書きは [`gates/`](gates/)（SFb で書いたが未実行）。
- アプリ側（`feature/losat-web-gui-app`、S13 の途中、`e200e60b`）の merge は S13 の終了後。
