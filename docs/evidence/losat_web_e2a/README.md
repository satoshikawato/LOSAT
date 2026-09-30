# LOSAT Web E2a-1（Session S06）ゲート記録

- 段階：E2a-1 BLASTN outfmt 0：権威と fixture（[総合計画書](../../losat_web_gui_plan.md) §7 の S06、[指示書](../../losat_web_gui_sessions/session_s06_e2a1_blastn_outfmt0_authority.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S05 の記録の後の `ccdc471e0`。LOSAT のコードは変えていない（エンジンは S05 の `edf825d36` のまま）。この記録と fixture は同じコミットに入る
- 権威の記録：[`AUTHORITY.md`](AUTHORITY.md)（英語。NCBI の経路の表、LOSAT の部品の対応、見つかった不具合、S07 のゲート）
- 実行記録：[`run-20260929T000046Z/`](run-20260929T000046Z/)。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)（リポジトリのルートで `sha256sum --check docs/evidence/losat_web_e2a/evidence.sha256` を実行する）
- 判定：**完了条件を満たした**。ただし、fixture を確かめる途中で、outfmt 0 とは別の既存の不具合を 5 つ見つけた（下の「見つかった既存の不具合」）。どれも推奨案で扱いを決め（S07 と、新しい S07+）、TLOSAN v0.2.0 の認証済みのコードに関わる 1 つを保守者に報告する

## 成果物

| ファイル | 内容 |
|---|---|
| [`AUTHORITY.md`](AUTHORITY.md) | NCBI の `blastn -outfmt 0` の呼出し経路（`blastn_app.cpp` → `CBlastFormat` → `showdefline.cpp` / `showalign.cpp`）をファイルと行の番号と逐語の断片で表にしたもの（§A）、BLASTP の報告との違い（§B）、`LOSAT/src/report/pairwise.rs` の部品の使い道（§C）、見つかった不具合と扱い（§D）、oracle の出力の例（§E）、S07 のゲートと拒否を続けるオプションと移植の順（§F） |
| [`LOSAT/tests/outfmt0_manifest.tsv`](../../../LOSAT/tests/outfmt0_manifest.tsv) | fixture の manifest。34 件（BLASTN 32、BLASTP 1、TBLASTN 1）。列は ID、program、task、query、subject、追加の引数、stdout の SHA-256 と大きさ、stderr の SHA-256、確かめる内容。すべて `LOSAT/` から実行する。`program` の列があるので、S08 の TBLASTX の fixture も同じ manifest に足せる |
| [`LOSAT/tests/fixtures/outfmt0/`](../../../LOSAT/tests/fixtures/outfmt0/) | 固定した NCBI の出力（`<ID>.out`、stderr があるものは `<ID>.err`）。合計 1.8 MB |
| [`LOSAT/tests/fasta/outfmt0/`](../../../LOSAT/tests/fasta/outfmt0/) | 作った入力 19 ファイル（`make_inputs.py` が固定の乱数の種からバイト単位で同じものを作る） |
| [`make_inputs.py`](make_inputs.py) | 入力を作る |
| [`run_oracle.py`](run_oracle.py) | manifest の全件を NCBI で実行する。既定は manifest のハッシュとの照合、`--freeze` は manifest のハッシュの列を書く（S06 だけ）。NCBI の出力を変える環境変数と `.ncbirc` があると実行しない。`search_argv` は S07 が LOSAT に同じ引数を渡すのに使う |
| [`precheck_hits.py`](precheck_hits.py) | 各 fixture を outfmt 6 で LOSAT と NCBI の両方で実行して比べる（outfmt 0 で一致するには、同じ HSP が見つかっている必要がある） |
| [`scoring_sweep.py`](scoring_sweep.py) | 既定以外の reward / penalty / gap の組合せを outfmt 6 で比べる（§D.5） |

## 完了条件と結果

| 完了条件（計画 §7 の S06） | 結果 | 証拠 |
|---|---|---|
| 経路の対応表 | 作成した。`AUTHORITY.md` §A の表の NCBI の断片は、引用した行の範囲にそのまま含まれることを機械的に確かめた（出力の文字列やオプションの名前を除く）。引用した NCBI のファイルは、固定 commit `598d8ae6` と改行以外に差が無い | `AUTHORITY.md` §A |
| 固定した fixture と SHA-256 | 34 件を NCBI BLAST+ 2.17.0 で固定した。2 回目の実行で全件が同じバイトになった。最長の実行は 0.43 秒 | `run-20260929T000046Z/oracle-freeze.log`、`oracle-check.log`、manifest |
| （指示書 5）S07 のゲートと、明示的に拒否するオプションの一覧 | `AUTHORITY.md` §F に書いた。一覧のオプションはすべて LOSAT の解析器がすでに拒否することを `edf825d36` で確かめた | `AUTHORITY.md` §F |

### fixture が確かめる範囲

指示書 3 の各項目と fixture の対応（詳細は manifest の `covers` の列）：

| 項目 | fixture |
|---|---|
| plus / minus 鎖 | `strand.*`、`multi.*`、`LC738874_LC738870.megablast`、`AP027152_AP027202.blastn` |
| 1 つの subject に複数の HSP | `strand.*`、`multi.*`（`msD` は plus と minus を 1 つずつ）、2 つのゲノムの組 |
| 複数の query と subject | `compact.*`（2 × 2）、`multi.*`（3 × 6。ヒットの無い query を挟む）、`many.*`（1 × 260） |
| ヒット無し | `compact.nohit.*`、`multi.*` の `mq2`、`edge.*` |
| `-lcase_masking`、`-dust` | `mask.*`（DUST、`-dust no`、query と subject の小文字、その組合せ）、`edge.alllower_lcase.blastn` |
| ギャップ | `strand.*`（行頭のギャップ）、`LC738874_LC738870.megablast` |
| 長い defline | `longdef.megablast`、`strand.*`、`LC738874_LC738870.megablast` |
| `-max_target_seqs` の境界 | `multi.mts1/mts3`（5 未満の警告）、`many.blastn`（説明 260・アラインメント 250）、`many.mts255`、`many.mts500`（省略と 500 の区別） |
| task megablast / blastn | ほぼすべての入力を両方の task で |
| そのほか | 座標の桁の境界（`width.*`）、無効な query（`edge.allN.blastn`）、HSP の選別のオプション（`multi.maxhsps1/perc97/besthit`）、既定以外の得点と 0 のギャップの末尾の式（`multi.r2p3g00.megablast`） |

入力は `LOSAT/tests/blastn_parity_manifest.tsv` の既存の入力（compact と 2 つのゲノムの組）をまず使い、既存の入力に無い場合（鎖、マスク、境界など）だけ小さな入力を作った。

除外した候補：`LC738873_LC738871.blastn`（7.2 MB、11,668 HSP）、`LC738874_LC738870.blastn`（1.3 MB、2,729 HSP）、`LC738874_LC738870.nodust.megablast`（1.1 MB、3,454 HSP）。どれも同じ書式の分岐を小さな fixture が確かめているので、リポジトリの大きさを優先した。outfmt 6 では後の 2 つも LOSAT と NCBI が一致した。

## 見つかった既存の不具合（outfmt 0 の writer の外。詳細と NCBI の根拠は `AUTHORITY.md` §D）

扱いは、セッションの規則（保守者に尋ねず推奨案で進め、記録する）に従って決めた。

| | 内容 | 扱い |
|---|---|---|
| D.1 | BLASTP と TBLASTN の outfmt 0 の座標の桁数が、NCBI（0 始まりの最大値）と違い 1 始まりの最大値から求められている。最大の座標がちょうど 10 の累乗のとき、空白が 1 つ多い（`width.blastp`、`width.tblastn`）。BLASTX は正しい | S07 で、NCBI の規則の 1 つの関数にまとめ、BLASTP・TBLASTN・BLASTN（S08 で TBLASTX）から使う。BLASTX は変えない（DW-10）。変わる凍結ハッシュはすべて列挙し、それぞれ NCBI の出力に一致することを示す |
| D.2 | BLASTN の `-task blastn -word_size 7`（または 8）で、ギャップ付き伸長の開始点を探す処理が subject の先頭を越えて panic する（`gapped.rs:559`）。NCBI は subject の両端の番兵のバイトで止まる | S07 の最初の作業（エンジンの変更、回帰試験は `multi.ws7_e1000.blastn`）。同じ修正を試しに入れたビルドは、この入力の word size 4〜8 と 34 件の fixture すべてで outfmt 6 が NCBI とバイト一致した（`probe-gapped-start-bound.diff`、`precheck-probe.tsv`。worktree のファイルは元に戻した） |
| D.3 | BLASTN は NCBI の警告（`-max_target_seqs` が 5 未満、無効な query）を stderr に出さない（どの形式でも） | S07 で足し、stderr を記録した 3 件で stderr もゲートにする |
| D.4 | BLASTN の `-max_target_seqs` は既定値 500 を持つので、省略と 500 を区別できない（outfmt 0 の表示数に必要） | S07 |
| D.5 | 認証済みの既定の得点（blastn 2/−3・5/2、megablast 1/−2・0/0）以外で：NCBI が拒否する組合せ 34 を LOSAT は受け付けて実行し、NCBI が受け付ける組合せのうち 17 は結果が違う（違うものはすべて reward 1 で gap が 0 でない） | 新しいセッション **S07+**（段階 E2c、S07 の直後）で、NCBI の検査（同じ拒否と文言）を移植し、違いの原因を調べる。終わっても一致しない組合せは明示的に拒否し、アプリが認証されていない BLASTN の得点を出さないようにする。S07 の fixture は認証済みの範囲に留める |

**保守者への報告（判断は不要）**：D.1 の TBLASTN の桁数は、TLOSAN v0.2.0 で認証した TBLASTN の outfmt 0 と同じコードにある。認証の fixture はこの境界に達しないので、認証の主張の範囲は変わらない。このブランチで S07 に直す。

## 実行の方法（再現）

```bash
# 入力を作る（既存のファイルと同じバイトになる）
python3 docs/evidence/losat_web_e2a/make_inputs.py
# NCBI の出力を manifest と照合する（固定するときだけ --freeze と --out LOSAT/tests/fixtures/outfmt0）
python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir /home/kawato/micromamba/bin --out <dir>
# LOSAT と NCBI の outfmt 6 の比較、得点の組合せの比較
python3 docs/evidence/losat_web_e2a/precheck_hits.py --bin-dir /home/kawato/micromamba/bin --losat <LOSAT>
python3 docs/evidence/losat_web_e2a/scoring_sweep.py --bin-dir /home/kawato/micromamba/bin --losat <LOSAT>
```

`<LOSAT>` は `edf825d36` から作ったネイティブの release の実行ファイル（`cargo +1.92.0 build --release --locked`）。

## 引き継ぎ（S07 へ）

- S07 の指示書に、`AUTHORITY.md` §F の移植の順（D.2 の修正と D.1 の関数から始める）とゲートを書き足した。
- 新しいセッション S07+ の指示書を作り、計画 §7 と README の表に行を足した。
- S08 は、TBLASTX の fixture をこの manifest（`program` の列）に足し、`run_oracle.py` で固定する。

---

# LOSAT Web E2a-2（Session S07）ゲート記録

- 段階：E2a-2 BLASTN outfmt 0：移植とゲート（[総合計画書](../../losat_web_gui_plan.md) §7 の S07、[指示書](../../losat_web_gui_sessions/session_s07_e2a2_blastn_outfmt0_port.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S06 の `63ddf7380`
- エンジン：`032a4b96c`（移植）、`1219c8990`（ABI v1 の固定）、`11b123e2d`（rustfmt）、`0388d1064`（独立監査の指摘）、`4404a082c`（V-PERF）
- アダプタ：`e228a0d8d`、`cb6876d24`、`61efc6a47`
- 実行記録：[`run-20260929T012339Z/`](run-20260929T012339Z/)（`head.txt` がゲートを実行した commit）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)
- 判定：**完了条件を満たした**。独立監査は 2 回目で supported。監査が見つけた、以前からの BLASTN の差のうち outfmt 0 の移植の外のもの（outfmt 7 の無効な query の行、FASTA の読み方、警告の時点）は S07+ に移した（下の「残件と扱い」）

## 変更の内容

### エンジン

| ファイル | 内容 |
|---|---|
| `LOSAT/src/report/pairwise.rs` | `write_blastn_pairwise_report`（BLASTN の報告の driver）、`write_blastn_description_table`（NCBI の `x_InitDeflineTable` の規則、`AUTHORITY.md` §G.3）、核酸の query の末尾（Gumbel の列なし、無効な query）、epilog（`blastn matrix`、gap extension の式）、megablast の文献。`write_hsp_info` の BLASTN の行（表示の行から数えた Identities と Gaps、` Strand=`）。`write_alignment_with_sequences` を核酸（`\|` の中線、下る subject の座標、空の行の規則）に広げ、座標の桁数を NCBI の 0 始まりの規則の `coordinate_width` にした（BLASTP と TBLASTN も、計画 TD-9）。BLASTX が呼ぶ関数は変えていない（DW-10） |
| `LOSAT/src/algorithm/blastn/pairwise.rs` | 最終の HSP 一覧から `PairwiseHit` を作る（編集操作から表示の行を作り、minus 鎖の HSP は query の plus 鎖で見せ、マスクを小文字にする）。query の組成から ungapped Karlin block を求める（`query_ungapped_karlin`、§G.2） |
| `LOSAT/src/algorithm/blastn/blast_engine/run.rs` | outfmt 0 の受け付けと報告の組み立て（`BlastnReportInputs`、`blastn_pairwise_report`、NCBI の query の batch の規則 `unsearched_queries`、§G.1）、`run_local` の `hits`、観測者、NCBI の 2 つの警告（`-max_target_seqs` が 5 未満、無効な query）。ABI v1 の BLASTN は outfmt 0 を以前と同じ文言で拒否する（TD-1） |
| `LOSAT/src/algorithm/blastn/alignment/gapped.rs` | ギャップ付き伸長の開始点の探索が subject の端で止まる（S06 §D.2 の panic の修正） |
| `LOSAT/src/algorithm/blastn/lookup.rs` | `reverse_complement` が小文字を大文字と同じ塩基として扱う（以前からの不具合、§G.6） |
| `LOSAT/src/algorithm/blastn/args.rs`、`hsp.rs` | `-max_target_seqs` は省略時に値を持たない（outfmt 0 の表示数、S06 §D.4）。outfmt 0 の受け付け |
| `LOSAT/src/blastinput/query_batch.rs`、`LOSAT/src/report/query_warnings.rs` | query の batch と、無効な query・少ない一致数の警告を共有にした（TBLASTN もこれを使う） |
| 試験 | `LOSAT/tests/run_local_blastn.rs`（0/6/7 と観測者、無効な query の batch、小文字の query）、`cli_v2.rs`、`unit/blastn/args.rs`、`gapped.rs`・`pairwise.rs`・`lookup.rs`・`query_batch.rs`・`query_warnings.rs` の単体試験、`web_api.rs`（ABI v1）、`LOSAT/tests/check_wasm_threading.py`（BLASTN の outfmt 0 を NCBI と比べる） |
| fixture | `LOSAT/tests/outfmt0_manifest.tsv` に 7 件を足した（41 件）：無効な query を挟む batch、IUPAC の記号の query（2 task）、最初の batch を埋める全 N の query（2 task）、minus 鎖の小文字の query（2 task） |

### アダプタと文書

- `web/adapter/src/run.rs`：BLASTN は 0/6/7 を書き、HSP レコードを出す。解析の binary の名前は CLI と同じ `LOSAT`。`RangeRecorder` は対応の無い・重複・終わらない HSP の出来事で run を失敗させる。
- `web/adapter/src/store.rs`：`register` の照合を関数に分け、残基の数も比べ、食い違いの単体試験を足した。
- `web/adapter/tools/build_reactors.py`：checkout のパスと `CARGO_HOME` を置き換える（計画 TD-11）。同一性の記録に checkout・`CARGO_HOME`・置き換えを残す。
- `web/adapter/tests/v_abi.js`、`tools/v_abi_cases.py`、`tools/run_v_abi_parallel.py`：BLASTN の outfmt 0 と HSP レコード、full に outfmt 0 の fixture（NCBI の凍結ハッシュ）、古い部分の結果を消す、成果物のハッシュ。
- `docs/web/abi_v2.md`：BLASTN の形式と HSP レコード、パスの置き換え、host の `-num_threads` の扱い。
- `AUTHORITY.md` §G：S07 で分かった NCBI の振る舞い（無効な query と batch、組成に依存する ungapped Karlin block、説明の一覧表の規則、ABI v1、S08 への部品、minus 鎖の小文字）。

## 完了条件と結果

| 完了条件（計画 §7 の S07） | 結果 | 証拠 |
|---|---|---|
| 固定した fixture で NCBI とバイト一致（stderr を記録したものは stderr も） | 通過。41 件（S06 の 34 件と S07 で足した 7 件）すべてが、`-num_threads` 1・2・4 で stdout と stderr まで NCBI BLAST+ 2.17.0 とバイト一致。S07 で足した 7 件の NCBI の出力は、2 回の実行で同じだった。41 件の outfmt 6 も NCBI と一致した（`precheck.tsv`） | `check-losat-n{1,2,4}.tsv`、`precheck.tsv`、`oracle-check.log`、`oracle-freeze-*.log` |
| 既存の 6/7 に退行なし | 通過。全 program の 236 件（BLASTN 14、BLASTP 28、TBLASTX 20、TBLASTN 162、BLASTX 12）の出力・stderr・終了コードが S02 の基準と一致（差 0）。凍結ハッシュの不一致は、S02 からの既知の `Sakai.MG1655.megablast` の 1 件だけ | `capture-compare.txt`、`capture/hashes.tsv` |
| TD-9 で変わる BLASTP・TBLASTN の凍結出力は NCBI と一致する | 変わる凍結出力は無かった（上の 236 件で差 0）。境界の fixture `width.blastp` と `width.tblastn` は NCBI とバイト一致 | 同上、`check-losat-n1.tsv` |
| BLASTN の全升目（0/6/7 × スレッド 1/2/4）の V-ABI | 通過。full：109 の検索（回帰の 68 と outfmt 0 の fixture の 41）を serial reactor（1 スレッド）と threaded reactor（1・2・4 スレッド）で実行した 436 件すべてで、全形式・HSP レコード・診断がネイティブの CLI と一致。BLASTN は 208 件で、HSP レコード（範囲を含む）も確かめた。凍結ハッシュは 652 件中 648 件が一致し、41 件の fixture の outfmt 0 は 164 件すべてが NCBI の凍結ハッシュと一致。一致しない 4 件は既知の `Sakai.MG1655.megablast` の outfmt 7。quick：52 件すべて一致 | `v-abi-full/summary.json`、`v-abi-full/v-abi-results.json`、`v-abi-quick/` |
| 独立監査 | 2 回目は supported（下の「独立監査」） | — |
| （規則 4）エンジンの変更の既存のゲート | v1 の WASI の検査が通過（423 件の command / oracle の記録、形式の失敗 0。BLASTN の outfmt 0 も NCBI と比べるようにした）。v1 の reactor の記録（threaded 172 件、serial 12 件）と `v1_requests.js` の 14 件の応答が S05 と一致。wasm32 だけでコンパイルされる v1 の試験 5 件が通過（ABI v1 の BLASTN が outfmt 0 を拒否し続けることを含む） | `wasm-threading.log`、`wasm-threading-metadata.json`、`v1-reactor-records-compare.txt`、`v1-requests-compare.txt`、`wasm32-web-api-tests.log` |
| （規則 4）V-PERF の非退行 | 通過（下の「性能の計測」） | `perf-1.json`、`perf-check-1.txt` |
| `cargo fmt --check`、`clippy -D warnings`、`cargo test --all-features` | 通過。エンジンは 4 つの構成の clippy と 805 件の試験、アダプタは 3 つの構成の clippy と 6 件の試験 | `fmt.log`、`clippy.log`、`adapter-clippy.log`、`cargo-test.log`、`adapter-test.log` |
| ビルドの同一性（TD-6、TD-11） | 通過。reactor には `/losat` と `/cargo` に置き換えたパスだけが残る | `reactors/build-identity.json`、`reactors/losat-web-*.json` |

ゲートは `61efc6a47`（`head.txt`）で実行した。LOSAT/ と web/adapter/ に commit との差は無い。成果物のハッシュは `artifacts.sha256` にあり、同じハッシュの複製を `/home/kawato/.cache/losat-web-gui-target/s07-final/` に残した（S07+ の変更前の成果物）。

### 性能の計測

変更前は S05 の成果物（`s05-final/`）、変更後はこの記録の成果物で、1 回ずつ交互に測った（`perf_cases.py`、1 回の暖機と 3 回の計測）。case は、BLASTN の `blastn`（起動の時間が大半）、`blastn-large`・`blastn-large-fmt7`（Gate A の EDL933 × Sakai）、座標の桁数を変えた BLASTP と TBLASTN の outfmt 0（`blastp-fmt0`、`tblastn-fmt0`）の 5 つ。BLASTN の outfmt 0 は変更前に無いので、比べる対象が無い。

15 件すべてが ×0.972〜×1.050 で、出力は同じ。この前のゲートの実行では、`blastn` が serial WASI で ×1.12〜1.14、threaded WASI で ×1.06 だった。

- **切り分け：** エンジンの時間（`LOSAT_TIMING` の `total`）は変わらなかった。Node の中の各段階（読み込み・compile・instantiate・実行）にも差は無かった。差は Wasm のコード生成の量（遅延 compile、eager compile、TurboFan のみのいずれでも差が出た）にあった。
- **原因：** outfmt 0 の関数が共有の後処理に inline され、どの BLASTN の run でも compile されていたこと。
- **対処：** `#[inline(never)]` でこれらを分けた（`4404a082c`）。10 回の計測で ×0.990（serial）と ×1.004（threaded）、この記録の計測で ×1.000 と ×0.981 になった。
- **大きい入力：** Sakai × MG1655 の megablast の outfmt 0 は、監査の指摘（M2）の前の 46.5 秒から 0.87 秒になった（outfmt 6 は 0.74 秒）。

## 独立監査

読み取り専用の独立監査（役割 `ncbi_parity_auditor`、コードを変えない別のエージェント）を 2 回受けた。

**1 回目（`11b123e2d` まで）：unsupported。** fixture の 37 件は一致していたが、次の指摘があった。

| | 指摘 | 扱い |
|---|---|---|
| M1（重大、S07 で入った） | 無効な query の footer を決める query の batch を、固定の大きさ（blastn 100,000、megablast 5,000,000 残基）と誤った NCBI の引用で決めていた。NCBI の blastn は `CBatchSizeMixer` で batch を決める | NCBI の最初の batch を正確に移植し、後の batch が決まらない場合は outfmt 0 を書く前に明示的に失敗させた（`0388d1064`、`AUTHORITY.md` §G.1）。fixture `edge.batch_allN.{blastn,megablast}` |
| M2（重大、S07 で入った） | 表示の行を作るときに HSP ごとに query と subject の全体を大文字にし、逆相補にしていた（O(HSP × 長さ)）。Sakai × MG1655 の megablast の outfmt 0 が 46.5 秒（outfmt 6 は 2.2 秒）。アダプタは毎回 outfmt 0 と HSP レコードを求めるので、Web の BLASTN のすべての run に効く | HSP の区間だけを読み、マスクは配列ごとに 1 度まとめて二分探索する（`0388d1064`）。0.87 秒で、NCBI とバイト一致 |
| P1（重大、以前から） | 小文字の query が minus 鎖で壊れる（`reverse_complement` が小文字を gap にしていた）。`-lcase_masking` の無い soft-mask の query の minus 鎖の HSP が、どの形式でも NCBI と違っていた | 大文字にしてから逆相補にする（`0388d1064`、§G.6）。fixture `mask.lcase_minus.{blastn,megablast}` と `run_local` の試験 |
| m1 | 一部の NCBI の引用が逐語でない、範囲の外 | 直した |
| m2 | V-ABI が BLASTN の HSP レコードを確かめていない | 確かめるようにした（`cb6876d24`、`61efc6a47`） |
| m3 | S05 の残り（`-num_threads`、`RangeRecorder`） | `RangeRecorder` を直した。`-num_threads` は host が付けるので拒否せず、host が実行する argv を検証する（`abi_v2.md` §7）と決め、繰り返しの試験を足した |
| m4 | `-outfmt` の help が古い | 直した |
| m5・m6（以前から） | 警告の時点と空の query、FASTA の読み方（tab、UTF-8 の折り返し、`U`、`X`） | S07+ で扱う（指示書に書いた） |

**2 回目（`0388d1064`・`cb6876d24`）：supported**（ゲートの結果待ち）。条件は、§G.1 の「表形式は影響を受けない」を直すことと、N1 を記録するか直すこと。M1・M2・P1 が直ったこと、41 件の一致、退行の無いこと（監査の 44 の入力のうち 42 が一致し、残り 2 つは N2 の拒否）を確かめた。

| | 指摘 | 扱い |
|---|---|---|
| N1（以前から） | outfmt 7 で、検索されなかった batch の query に LOSAT は `# 0 hits found` を出すが、NCBI は出さない（`tabular.cpp:1277-1282`） | §G.1 を直し、S07+ で outfmt 0 と同じ規則を outfmt 7 に使う |
| N2 | 最初の batch の後の 100 残基未満の無効な query の連なりで後ろに有効な query があるものは必ず検索されるが、今の fail-fast はこれも拒否する | S07+ |
| N3 | subject の総文字数が 2^31 以上では、NCBI の `int` の変換を再現しない | S07+ |

監査の後、ゲートの V-PERF で、起動の時間だけを測る `blastn` の case が WASI で遅くなった（下の「性能の計測」）。原因を直した `4404a082c` は、出力を変えない（`#[inline(never)]` だけ）。

## 残件と扱い

- **S07+ に移したもの**：無効な query と batch の N1〜N3、FASTA の読み方と警告の時点（m5・m6）、query の組成に依存する ungapped Karlin block をエンジンの gap trigger と X-drop にも使うか（§G.2）、既定以外の得点（S06 §D.5）。S07+ の指示書に書いた。
- **S08 に移したもの**：BLASTP・TBLASTN・BLASTX の説明の一覧表は、今も各 subject の最初の HSP と普通の幅の規則を使う（§G.3）。TBLASTX の一覧表で NCBI の規則を使い、BLASTP と TBLASTN を切り替えるかを決める。`validate` が挿入する `-outfmt 6` をやめ、文言を CLI と同じにする。
- **計画 TD-11**：reactor に機械のパスは残らなくなったが、バイトは checkout のパスに依存する（Cargo の metadata hash）。同一性の記録に checkout のパスを残し、再現は同じパスで行う。
- **保守者への報告（判断は不要）**：P1 は、以前からの BLASTN の不具合で、`-lcase_masking` の無い小文字を含む query の minus 鎖の HSP の結果（全形式）が NCBI と違っていた。回帰の fixture に小文字の query は無いので、Gate A と Stage G のハッシュは変わらない。S06 §D.1 の TBLASTN の座標の桁数（TLOSAN v0.2.0 の認証済みのコード）は、このブランチで直した。

## 実行の方法（再現）

```bash
# fixture（NCBI の凍結出力）との比較
python3 docs/evidence/losat_web_e2a/check_losat.py --losat <LOSAT> [--threads N]
# 回帰の出力（S02 の基準と比べる）、v1 の WASI の検査、V-ABI、性能：S02・S05 の記録の「実行の方法」と同じ
python3 docs/evidence/losat_web_e1a/capture_outputs.py run --losat <LOSAT> --out <dir> --jobs 8
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare docs/evidence/losat_web_e1a/baseline/hashes.tsv <dir>/hashes.tsv
python3 web/adapter/tools/v_abi_cases.py --suite full --out <cases.json>
python3 web/adapter/tools/run_v_abi_parallel.py --native <LOSAT> --serial <reactors>/losat-web-serial.wasm \
  --threads <reactors>/losat-web-threads.wasm --cases <cases.json> --out <dir> --jobs 12
python3 docs/evidence/losat_web_e1c/perf_cases.py run --before <S05 の 3 つ> --after <S07 の 3 つ> \
  --cases blastn,blastn-large,blastn-large-fmt7,blastp-fmt0,tblastn-fmt0 --out <file.json>
```

## 引き継ぎ（S07+・S08 へ）

- S07+ の指示書に、独立監査の N1〜N3、m5・m6、§G.2 を書き足した。S07+ は S12（検索画面）の前に終える（計画 TD-10）。
- S08 の指示書に、S07 の部品（`coordinate_width`、`write_sequence_row`、§G.3 の一覧表の規則、`query_warnings.rs`、`query_batch.rs`、`run_oracle.py`・`check_losat.py`・manifest）と、`validate` の `-outfmt 6` の挿入をやめることを書いた。
- `docs/web/verification_cells.tsv` の BLASTN の outfmt 0 の升目（ネイティブと V-ABI）を埋めた。
