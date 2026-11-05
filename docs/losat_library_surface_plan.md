# LOSAT ライブラリ化計画（公開面の設計メモ）

状態：**LIB-1 は `dev` に入った（#845）。段階 1（公開面を固める）は `dev` に入った（#846）。段階 2（結果と実行制御、#849）と段階 3（`pylosat`、#850）は積んだ PR が feature branch に merge されたので、1 本の PR で `dev` へ入れる（2026-11-05）。** 設計メモは LIB-0（2026-10-31、保守者「OK」）。v0.2.0 自体は crates.io に出さない（release の決定 5、2026-10-25）。LIB-1 で決めたことは §9、順番の変更と段階 1 で決めたことは §10、段階 2 は §11、段階 3 は §12。

| 項目 | 内容 |
|---|---|
| 目的 | LOSAT を Rust crate（crates.io）、npm、PyPI から使えるようにするときの**公開面**（public API surface）を先に決める |
| 対象 | `LOSAT/`（エンジン crate）の公開面、`web/adapter/` との分担、配布の 3 チャネルの形と順番 |
| 対象外 | wasm-bindgen、ネイティブ C ABI、`psiblast` と `makedb` の API 化（§8） |
| 基点 | `origin/dev` `4b2588e01` で確認し、`61fe55d2d` で関係するファイルに差分が無いことを確かめた（2026-10-31） |
| 規則 | ルートの [`AGENTS.md`](../AGENTS.md)（エンジンの変更は NCBI の権威と parity gate）、[`web/AGENTS.md`](../web/AGENTS.md)、[`PD-LOSAT-WEB-APP-BOUNDARY`](product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md)、[`PD-PURE-RUST-RUNTIME-AUTHORITY`](product_decisions/PD-PURE-RUST-RUNTIME-AUTHORITY.md) |
| 記録 | `losat-baselines/library-20261031/`（DECISIONS、STATE、次のセッションの指示） |

---

## 0. 結論

1. 三つのチャネルの前提は同じで、**公開面を `losat::api` の 1 モジュールに絞る**ことである（`cli` も隠す、§10.2）。今の lib は内部モジュールが全部 `pub` で（§1.2）、このまま publish すると docs.rs が内部実装の目録になり、互換の約束ができない。
2. 順番は **公開の前に Python で API を確かめる**：公開面を固める → 結果と実行制御 → `pylosat` で検証 → 公開（§10.1。2026-11-05 に変更。初版の順番は crates.io → npm → PyPI）。
   - crates.io は、Cargo.toml のメタデータが既に揃っていて、`cargo package --list` が通る（§1.1）。公開面を絞る LIB-1 の後なら、残りは名前とメタデータだけである。
   - npm は **wasm-bindgen を使わず**、既にある 2 つの WASI リアクタと ABI v2 のホスト側 TypeScript をそのまま包む（§4.2）。他のチャネルと独立に進められる。
   - PyPI は **pyo3 の binding `pylosat`**（§4.3、§9.1、§10.1）。
3. 公開面の形は LOSAT Web の ABI v2（[`docs/web/abi_v2.md`](web/abi_v2.md)）を Rust の型にしたものにする（§3）。検索条件は CLI と同じ argv 1 本、入力はメモリ上の FASTA、出力は outfmt ごとのバイト列 + HSP レコード + 診断。ABI v2 がこの形で 4 program を既に動かしているので、設計の新しさはない。
4. エンジンの動作は変えない。blastn / blastp / tblastn / tblastx の façade（§3）は既存の入口（各 program の `run_local`、`ReportOutputs`、`FormatObserver`）を包むだけである。blastx には `run_local` が無く（§1.2）、ABI v2 も blastx を断っているので、blastx は LIB-1b で入口を足す（§5、§7 Q6）。CLI の出力バイトは既存の gate で不変を確かめる（§6）。

## 1. 現状（dev `4b2588e01` で確認）

### 1.1 crate

- `LOSAT/Cargo.toml`: package 名 `LOSAT`、version `0.2.0`、`rust-version = "1.92"`、`[lib] crate-type = ["rlib", "cdylib"]`、`[[bin]] name = "losat"`、`readme = "../README.md"`、license / keywords / categories / exclude あり。feature は `parallel`（既定）、`wasm-threads`、`xstats`、`diagnostics`、`blastx-worker-probe`、`tblastx-wasm-scalar`。
- `cd LOSAT && cargo package --list --allow-dirty --no-verify` は通る（497 ファイル、`README.md` と `build.rs` を含む）。`exclude` の `"/matrix"` は、そのディレクトリがもう無いので死んでいる。非 test の `include_str!` に package の外を読むものは無い。
- `build.rs` は `wasm32-wasip1` / `wasm32-wasip1-threads` 以外では何もしない（cfg の宣言だけ）ので、docs.rs（Linux、既定 feature）でそのまま build できる。
- README §Installation（114 行目）は「crates.io、Homebrew、PyPI、Docker Hub、npm には出していない。Bioconda だけ」と明記している。この計画の各段が終わるごとに、この文を直す。
- 名前 `losat` は crates.io、PyPI、npm のどれも未取得である（2026-10-31 確認、HTTP 404）。

### 1.2 公開面

- `LOSAT/src/lib.rs` は `#![allow(warnings, clippy::all)]` で、19 のモジュール（`algorithm`、`core`、`seqdb`、`writedb`、`blastdbcmd`、`cli` など。`web_api` は wasm32 のみ）を全部 `pub mod` にしている。非公開のモジュールは無い。`pub` 項目 約 2,600、`pub(crate)` を含む行 約 1,600（`src/` 約 49 万行）。
- `api` モジュールは `mod.rs`（12 行）と `local_blast.rs`（353 行）だけで、中身は出力の束ね方（`OutputSink`、`FormatOutput`、`ReportOutputs`（`formats`、`diagnostics`、`hits`、`observer`）、`FormatObserver`、`HspIndex`）。
- 検索の入口は program ごとの `run_local(args, query_records, subject_records, outputs)` で、blastn（`blast_engine/run.rs`）、blastp（`blast_engine.rs`。query / subject のラベルが 2 つ余分）、tblastn（`tblastn_app.rs`。`-db` を断る: Web のホストにはファイルが無い）、tblastx（`blast_engine/run_impl.rs`）の 4 つにある。**blastx には無い**: 入口は `native::run`（CLI）と ABI v1 用の `web::run_pair`（`pub(crate)`）だけで、ABI v2 も blastx を断っている（`abi_v2.md` の冒頭、adapter `run.rs`）。
- `args` は clap 派生の `BlastnArgs` 等。`records` は **LOSAT 自前の** `blastinput::fasta_reader::FastaRecord`（`local_id`、`title: Vec<u8>`、`sequence`、`warnings`、`source`）で、reader は `FastaInputSource::from_bytes(bytes, ReaderConfig::query / subject(...))`。CLI と adapter（`store.rs`）が同じ reader を使う。rust-bio（`bio`）は 2026-11-05 に依存から外した（計画の DW-48）。
- 結果は返り値ではなく、`ReportOutputs` の writer（outfmt ごと）と `FormatObserver`（HSP の位置）と `hits` callback（`PairwiseHit`）に流れる。1 回の `run_local` の outfmt は同じ検索条件でなければならず（blastp は違うと断る）、adapter は hit list の大きさが違う outfmt を `search_groups`（`run.rs`）で別の検索に分けている。
- `cli` モジュール: `Cli` / `Commands`（blastx、blastn、blastp、tblastx、tblastn、psiblast、makedb の 7 つ。`blastdbcmd` は `Commands` の外で `main.rs` が扱う）、`try_parse_from`、`render_message`、`ncbi_version_request`、`exit_on_native_error`（`main.rs` 以外で `process::exit` を呼ぶ唯一の場所）、`report_standard_output`、`ReportStream`、`NativeError`、`inaccessible`、`command_line_bytes`。

### 1.3 利用者

| 利用者 | 入口 | 備考 |
|---|---|---|
| CLI `losat` | `LOSAT/src/main.rs`（323 行） | `cli::try_parse_from` → 各 program |
| LOSAT Web（`web/app`） | ABI v2（`web/adapter`、export 13 個、`losat_web2_*`） | WASI preview1 リアクタ 2 本（`wasm32-wasip1`、`wasm32-wasip1-threads`）。認定ビルドと同一性を検査（plan TD-6、`check_build_identity.py`） |
| gbdraw Web | ABI v1（`LOSAT/src/web_api.rs`、`losat_web_run_pair` など export 13 個） | 凍結（fail-fast の修正だけ、plan TD-1） |
| gbdraw ネイティブ | subprocess（解決順: 明示の path、conda、managed、PATH の `losat`、NCBI BLAST+。`--losat_bin` あり） | `gbdraw/comparisons/losat_runtime.py`。生の出力をキャッシュする |

### 1.4 ライブラリとして邪魔になる振る舞い

| 振る舞い | 場所 | 扱い（§3.4） |
|---|---|---|
| 進捗バー（indicatif 0.18、必須 dependency）を stderr に描く | blastn `run.rs` の 1 か所だけ（`subjects_parallel` = 並列 かつ subject 2 本以上 のとき無条件に `ProgressBar::new`）。tblastx は `ProgressBar::hidden()` | indicatif は stderr が端末でなければ描かない（`draw_target.rs` の `term()`、確認済み）。それでもライブラリの既定では作らない: 端末を持つ呼び出し側（例: Python の REPL）でも描かないため |
| stderr への診断 | `[TIMING]`（`LOSAT_TIMING`）、`print_summary`（`LOSAT_DIAGNOSTICS`）、`[losat-thread-pool]` / `[losat-thread-activity]`（`LOSAT_WASI_THREADS_DEBUG`）、tblastx の統計（`LOSAT_X_STATS`）、`BLEMIR_DEBUG`、`LOSAT_DEBUG_BLASTN`、`-verbose` | 全部 switch で gate されていて既定の経路では出ない（標本 15 か所。網羅は §6 の stderr 検査で）。`ignored_notices` の出力は `main.rs` だけ。そのまま |
| `process::exit` | `cli::exit_on_native_error` | `main.rs` 専用。façade は呼ばない（Web adapter も呼んでいない） |
| rayon の pool | `utils/threading.rs` `with_search_pool` | 検索ごとに local pool を作る（`build_global` は無い）。`-num_threads` 2 以上のときは rayon の worker の中からの呼び出しを断る（`current_thread_index().is_none()`）。1 のときは worker の中でも動く。API の文書に書く |
| 大域の状態 | `OnceLock` / `LazyLock` / `Mutex` / atomic（127 ファイル。`lazy_static` は無い）。`utils/xtiming.rs` の `START` と `SEARCH_SPANS` は process 全体で共有。`web_api.rs` の FASTA store は wasm32 のみ | 既定 build（`diagnostics` feature なし）で 2 検索を同時に走らせて出力が変わらないことを LIB-1 で検証し、文書に書く（§6）。診断を有効にした build では計時が混ざりうる |
| 環境変数 | `LOSAT_X_*`（性能 switch と診断） | 読む。診断系は `diagnostics` feature の build だけが読む（`utils/xdiag.rs`）。ライブラリの利用者には「診断用、既定で off」と文書に書く |

## 2. 方針

- **P1 検索条件は argv 1 本。** CLI と同じ clap の定義で parse する（`cli::try_parse_from`）。NCBI のオプション互換の実装を 1 か所に保つ。ABI v2 §7 と同じ。
- **P2 入力はメモリ上の FASTA か、`-db` の path。** メモリ入力は LOSAT 自前の `FastaRecord` を公開して受け、CLI / adapter と同じ reader（`FastaInputSource` と `ReaderConfig`）で読む。`-db` はファイルシステムのあるホストだけ（ネイティブ）。
- **P3 出力は、要求した outfmt ごとのバイト列 + HSP レコード + 診断。** HSP レコードの項目は ABI v2 §8 と同じ。バイト列は CLI の出力と一致する（§6）。
- **P4 公開面は `api` だけ（§10.2）。** 他のモジュールは `#[doc(hidden)]` を付け、semver の保証外と README に書く。`pub(crate)` 化は別の線（dead-code 警告が今は `allow(warnings)` で隠れていて、一斉に出る）。
- **P5 ライブラリは stdout / stderr に触れず、`process::exit` しない。** 出力は呼び出し側の writer だけに書く。
- **P6 新しい BLAST の動作を入れない。** 4 program の façade は既存の入口を包むだけ。blastx の入口（LIB-1b）は `native::run` を `run_local` の形に分けるもので、動作は変えない。`AGENTS.md` の規則（NCBI の権威、参照コメント、parity gate）はそのまま適用する。参照は NCBI の API 層（`c++/src/algo/blast/api/`、`local_blast`、`blast_results`）で、`api/mod.rs` の既存の参照と同じ。
- **P7 版。** crate の version は CLI の version と同じ。0.x の間は `api` の破壊的変更で minor を上げる。ABI v2 §10 と同じ考え方。

## 3. 公開面（案）

最終の signature は LIB-1 の PR の review で決める。ここでは形と境界を決める。実装した形は §9.2 と §10.2。

### 3.1 `losat::api`

```rust
pub mod api {
    /// 既存の reader をそのまま公開する（CLI と adapter が使うもの）。
    pub use crate::blastinput::fasta_reader::{FastaRecord, FastaInputSource, ReaderConfig};
    pub enum Role { Query, Subject }
    /// adapter の register と同じ ReaderConfig::query / subject で bytes を読む。
    pub fn read_fasta(program: Program, role: Role, bytes: &[u8]) -> Result<Vec<FastaRecord>, Error>;

    /// API で使える program。blastx は LIB-1b から（§5）。psiblast と makedb は CLI のみ（§8）。
    pub enum Program { Blastn, Blastp, Tblastn, Tblastx, Blastx }
    /// program の option と outfmt の一覧（adapter の describe.rs を移す。clap の定義から生成）。
    pub fn describe(program: Program) -> Description;

    /// argv を CLI と同じ parser で検証した結果。失敗の文言は cli::render_message と同じ。
    pub fn validate(program: Program, argv: &[&str]) -> Result<Search, Error>;
    pub struct Search { /* 検証済みの Commands。-out / -outfmt / -query / -subject は §3.3 */ }

    pub enum Subject<'a> { Records(&'a [FastaRecord]), Db(&'a std::path::Path) }
    pub struct Inputs<'a> { pub query: &'a [FastaRecord], pub subject: Subject<'a> }

    /// 低水準: 1 つの検索条件で、既存の ReportOutputs に書く（writer を呼び出し側が持つ。ストリーム）。
    pub fn run(search: Search, inputs: Inputs<'_>, outputs: &mut ReportOutputs<'_>) -> Result<(), Error>;
    /// 高水準: 要求した outfmt のバイト列と HSP レコードと診断を集めて返す。
    /// hit list の大きさが違う outfmt は、adapter の search_groups と同じく別の検索に分ける。
    pub fn run_collect(search: Search, inputs: Inputs<'_>, formats: &[&str]) -> Result<Outcome, Error>;
    pub struct Outcome { pub formats: Vec<(String, Vec<u8>)>, pub hsps: Vec<HspRecord>, pub diagnostics: Vec<u8> }

    /// ABI v2 §8 の項目をそのまま。
    pub struct HspRecord {
        pub index: usize, pub q_idx: usize, pub s_idx: usize, pub rank: usize,
        pub raw_score: i32, pub bit_score: f64, pub e_value: f64,
        pub q_start: usize, pub q_end: usize, pub s_start: usize, pub s_end: usize,
        pub query_frame: Option<i32>, pub subject_frame: Option<i32>, pub subject_length: Option<usize>,
        pub query_aligned: Option<String>, pub subject_aligned: Option<String>,
        pub out6: Option<(u64, u64)>, pub out0: Option<(u64, u64)>, pub out0_subject: Option<(u64, u64)>,
    }

    pub enum Error { Argv(String), Input(String), Search(String) }   // 文言は CLI と同じ（ABI v2 §4 の規則）
    // 既存のまま公開: OutputSink, FormatOutput, ReportOutputs, FormatObserver, HspIndex
}
```

- `run_collect` は、`search_groups` の分け方で `run` を呼び、`Vec<u8>` の writer と observer で集める。`search_groups` は program 非依存なので `api` に移す（LIB-1 でモデルにし、LIB-2 で adapter がそれを使う）。
- スレッド数は argv の `-num_threads`（CLI と同じ検証）。2 以上を rayon の worker の中から呼ぶと `Error::Search`（§1.4）。
- 同時に 2 検索を走らせてよいかは LIB-1 で検証して文書に書く（§1.4、§6）。
- 中断（cancel）はバッチの境界で扱う（段階 2、§10.1）。ブラウザではホストが worker を止める（§8）。

### 3.2 `losat::cli`

~~今のまま公開する。~~ `#[doc(hidden)]` にした（§10.2）：`Commands` は全オプションの引数構造体を含み、オプションを足すたびに破壊的変更になる。argv を取るホストには `validate` がある。Web adapter、`main.rs`、tests は `LOSAT::cli` をそのまま使える。

### 3.3 argv の扱い

ABI v2 §7 を基にするが、入力の渡し方が違う。

- `-outfmt` と `-out` は argv で受け取らず、`run` / `run_collect` の引数で指定する。argv にあれば `Error::Argv`（ABI v2 と同じ規則。adapter の文言を流用。保守者委任、§7 Q5）。
- ABI v2 は `-query <name>` / `-subject <name>` に登録した入力の名前を要求するが、API は `Inputs` でレコードを直接渡す。`Inputs` にレコードを渡すときに argv に `-query` / `-subject` があれば `Error::Argv`。`Subject::Db` のときの `-db` は argv に書く（path は呼び出し側のファイルシステム）。LIB-2 で adapter は、名前 → `Inputs` の対応を自分の store で持つ。
- `-num_threads`、`-task`、`-evalue` などの検索条件は argv のまま。既定値は CLI と同じ。

### 3.4 隠すもの

- `lib.rs` の `api` 以外の `pub mod` に `#[doc(hidden)]` を付ける。項目は今までどおり使えるが、docs.rs に出ず、semver の保証外（README と crate の doc に書く）。`api` が再 export する `FastaRecord` 等は `blastinput` の中にあるので、再 export の側だけを文書化する。
- 進捗バーは `ReportOutputs` に `progress: Option<...>` を足して CLI（`main.rs`）だけが渡すか、indicatif を optional dependency にして CLI の feature だけで有効にする（LIB-1 で決める）。
- `anyhow::Error` は `api::Error` に包む。`clap` は argv が API なので必須 dependency のまま。

## 4. 配布の 3 チャネル

### 4.1 crates.io（LIB-3）

| 項目 | 案 |
|---|---|
| crate 名 | **`losat`（小文字）に改名**（保守者委任、§7 Q1）。bioconda / PyPI / npm と同じ綴りになる。crates.io は大文字小文字を区別せずに名前を予約するので、`LOSAT` で publish しても `losat` は使えないし逆も同じ。`[[bin]] name = "losat"` は今のまま。`LOSAT::` の path を使うファイルは 49（`src` 1、`tests` 43、`web/adapter` 5）で、機械的な置換（Sonnet の仕事） |
| メタデータ | `rust-version` あり。`readme = "../README.md"` は `cargo package` が `README.md` として入れることを確認済み。死んだ `exclude` の `"/matrix"` を消す。`[package.metadata.docs.rs]`（既定 feature、Linux）を足す |
| feature | 今の 6 つはそのまま。進捗バーは §3.4 |
| gate | `cargo publish --dry-run --locked`（`release-readiness.yml` に足す）、`cargo doc --no-deps` が `api` だけを示す、tarball の大きさ（`cargo package` の出力） |
| publish | 保守者が行う（外向きの操作）。最初の版は cut 後の dev の checkpoint から |
| 版 | CLI と同じ（`0.2.x` の patch か `0.3.0`。どちらかは保守者、§7 Q4） |

### 4.2 npm（LIB-4）

**wasm-bindgen は使わない。** 理由: (a) ターゲットが `wasm32-unknown-unknown` になり、`wasm32-wasip1-threads` のスレッド版が無くなる（wasm-bindgen のスレッドは `wasm-bindgen-rayon` と build-std の組み合わせで、認定ビルドと別物になる）、(b) ABI v1、v2 に続く第 3 の ABI になる、(c) plan TD-6 の「認定ビルドとの同一性」が成り立たない。

代わりに、**既にあるものを包む**:

| 内容 | 出どころ |
|---|---|
| `losat-web-serial.wasm`、`losat-web-threads.wasm` | `web/adapter/tools/build_reactors.py`（同一性の記録つき） |
| ABI v2 のホスト側 | `web/app/src/infra/reactor/abi.ts`（`losat_web2_*` の呼び出し）、`infra/reactor/wasi.ts`（WASI shim。npm の `@bjorn3/browser_wasi_shim` に依存）、`infra/engine-worker/*`（worker、`thread-host.ts`）、`ports/engine.ts`（型と `EngineGateway`）。アプリが使う層と同じものを package にし、重複させない（切り出し方は LIB-4 で決める。`web/AGENTS.md` の規則の下） |
| 型 | `HspRecord`、`ProgramDescription`、`ValidationResult`（`engine.ts` に既にある） |
| README | cross-origin isolation（COOP / COEP）が threads 版に必須、メモリの候補（ABI v2 §3）、`abi_version` |

- package 名は **`losat`（スコープなし）**（保守者委任、§7 Q2）。版は engine の version、`losat_web2_abi_version` を同梱の表に書く。
- gbdraw Web が ABI v1 から v2 へ移る道になる。移るかどうかは gbdraw の製品判断で、この計画の対象外。

### 4.3 PyPI（LIB-5）

保守者の決定（§9.1、§10.1）：**pyo3 の binding `pylosat`** を `api` の上に薄く作る。初版の案（maturin の bin wheel だけ）は採らない。

- 置き場は `python/`（独自の Cargo root）で、`web/` と同じ境界規則を当てる（[`PD-LOSAT-PYTHON-BINDING-BOUNDARY`](product_decisions/PD-LOSAT-PYTHON-BINDING-BOUNDARY.md)、[`python/AGENTS.md`](../python/AGENTS.md)）。その PR は保守者が merge する。NCBI にない要素（Arrow など）はそこに置き、`LOSAT/` には入れない。
- `panic = "unwind"` はその crate の profile で設定する。LOSAT の既定の `abort` のままだと、panic が Python の process ごと落とす。
- 段階 3 は path 依存のまま動かし、Python 版の gbdraw で確かめる。関門は gbdraw が subprocess なしで動くこと。
- 配布は wheel のマトリクス（manylinux x86_64 / aarch64、macOS arm64 / x86_64、Windows x64）。release のバイナリと同じ build（`release-readiness.yml` / release workflow）から作る。詳細は LIB-5。

## 5. 実装の段（cut 後、PR はすべて `dev` へ）

順番は §10.1 の 4 段階。LIB-1 と段階 1、2 の後に `pylosat`（LIB-5 の先行版）で API を確かめ、LIB-3、LIB-5、LIB-4 を公開する。LIB-1c は要望に応じて前後させ、LIB-2 は段階 3 の後で、公開の前提ではない。

| 段 | 内容 | 触る場所 | gate |
|---|---|---|---|
| LIB-0 | この計画（このメモ） | `docs/` | — |
| LIB-1 | façade（4 program）: `api::{read_fasta, Program, describe, validate, Search, Inputs, run, run_collect, HspRecord, Error}` と `search_groups` の移設、`#[doc(hidden)]`、進捗バーの切り離し、`examples/library_search.rs` | `LOSAT/src/lib.rs`、`LOSAT/src/api/`、`LOSAT/src/main.rs`（進捗バー） | 既存の parity gate 全段（出力バイト不変）、adapter の V-ABI quick、`cargo doc --no-deps`、例の出力 = CLI の出力、同時 2 検索（§6） |
| LIB-1b | ~~blastx の入口~~ 不要（Web の線の SX で入った、§9） | — | — |
| LIB-1c | `-db` の入口: 各 program の `run_local` が `-db` を断る（§9.2）。ネイティブのホストに database の subject を開く入口を足す（`Subject::Db`、§3.1）。動作は CLI の `-db` と同じ | `LOSAT/src/algorithm/*/`、`LOSAT/src/api/` | `-db` の parity gate（既存の database fixtures） |
| LIB-2 | adapter を `api` の上に載せ替え: program 非依存の部分（`run.rs` の parse / validate / `search_groups` / HSP の収集、`describe.rs`）を `LOSAT/src/api` へ移し、adapter は WASI の export、store（名前 → `Inputs`）、scan、emit だけにする | `web/adapter/`、`LOSAT/src/api/` | V-ABI full（定義は `docs/losat_web_gui_plan.md` と `web/adapter/tools/v_abi_cases.py`）、同一性の検査（TD-6） |
| LIB-3 | crate 名、メタデータ、`cargo publish --dry-run` の gate、README の文、changelog | `LOSAT/Cargo.toml`、49 ファイルの path、`.github/workflows/release-readiness.yml`、`README.md` | `cargo publish --dry-run --locked`、CI 全緑 |
| LIB-4 | npm package | `web/`（切り出し方は LIB-4 で決める）、`.github/workflows/web.yml` | V-ABI quick を package の build 物で、`npm pack` の内容の確認 |
| LIB-5 | PyPI の `pylosat`（pyo3、§4.3）。先行版は段階 3 | `LOSAT/` の外の crate、release workflow | gbdraw が subprocess なしで動く（段階 3）、wheel で fixture 1 件 |

LIB-1、LIB-1b、LIB-2 は `LOSAT/src` を触るので `AGENTS.md` の規則の対象（NCBI の参照コメント、parity gate）だが、動作の変更は無い。LIB-1 が設計の判断を含む唯一の段で、残りは機械的。

## 6. 受け入れ条件

- CLI の出力バイトが不変（`fast output regressions` と skill `losat-gates` の full tier）。
- Web の V-ABI（LIB-1 は quick、LIB-2 は full）。
- `cargo doc --no-deps` が `api`（と `#[doc(hidden)]` でない再 export）だけを示す。
- `examples/library_search.rs`: 2 つのメモリ上のレコードで blastp を走らせ、`run_collect` の outfmt 6 と outfmt 0 のバイト列が、同じ fixture の CLI の出力と一致する（`tests/` の 1 ケースにする）。stderr に何も出ない（subprocess で確かめる）。
- 同時 2 検索: 既定 build で 2 スレッドから `run_collect` を同時に呼び、1 検索ずつの出力と一致する（`tests/` の 1 ケース）。
- `cargo package --list` と `cargo publish --dry-run --locked` が通る。

## 7. 保守者への問い（委任した推奨と、違う選択肢）

常設指示（2026-09-29、2026-10-07 再確認）により、推奨を採って進め、ここに記録する。違う選択をするときは保守者が言う。

| # | 問い | 委任した推奨 | 他の選択肢 |
|---|---|---|---|
| Q1 | crate 名 | `losat`（小文字）に改名 | `LOSAT` のまま（code の変更なし、綴りが他のチャネルと違う） |
| Q2 | npm の名前 | スコープなしの `losat` | `@satoshikawato/losat` など（org / user が要る） |
| Q3 | Python | ~~maturin の bin wheel だけ~~ → **pyo3 の binding `pylosat` を作る（保守者 2026-11-05、§9.1、§10.1）** | maturin の bin wheel だけ |
| Q4 | 時期と版 | cut（11-01）後に LIB-1 を dev で始め、`0.3.0` で出す | `0.2.x` の patch に入れる |
| Q5 | argv の `-outfmt` / `-out` | ABI v2 と同じく断る | 受け取って既定にする |
| Q6 | blastx | LIB-1 は 4 program、blastx は LIB-1b で入口を足してから | LIB-1 に含める（LIB-1 が blastx の engine の変更を抱える） |

## 8. 対象外と決めないこと

- wasm-bindgen（§4.2）。pyo3 は §9.1 で対象に入った（LIB-5、段階 3）。
- ネイティブの C ABI（`cdylib` は wasm のために既にあるが、C のヘッダと FFI は作らない）。
- `psiblast` と `makedb` の API 化（CLI のみ。`blastdbcmd` も同じ）。
- 検索の途中での中断。バッチの境界で止めることは対象に入った（段階 2、§10.1）：Python の Ctrl-C は同じ process の中で止めるほかない。ブラウザではホストが worker を止める。バッチの中で止めるには NCBI の `SetInterruptCallback` に当たるものをエンジンに入れる（まだやらない）。
- 内部モジュールの `pub(crate)` 化（別の線。`#[doc(hidden)]` で足りる間はやらない）。

## 9. LIB-1 の決定（2026-11-05、PR `feature/library-api`）

基点は `origin/dev` `701657e69`。dev では blastx にも `run_local` が入っていて（`algorithm/blastx/local.rs`、Web の線の SX）、ABI v2 も blastx を動かす。LIB-1b（§5）はもう要らず、LIB-1 は 5 program を包む。

### 9.1 保守者の決定

| # | 問い | 決定 |
|---|---|---|
| Q3 | Python | **pyo3 の binding を作る**（2026-11-05。§4.3 の推奨を変更）。`api` の上に薄く、LIB-5 で。保守者の呼び名は `pylosat`；PyPI の名前は LIB-5 で確定する |
| Q1 | crate の形 | 計画どおり（別の façade crate ではなく、エンジン crate の `api` に足す。改名は LIB-3） |
| — | 段 | LIB-1 を 1 PR で。npm（LIB-4）と PyPI（LIB-5）は後の PR |

### 9.2 公開面（実装したもの；§3 の案との違い）

`run`、`spans`、`search_groups`、`Inputs` の field、`cli` の扱いは段階 1 で変えた（§10.2）。

`LOSAT/src/api/`：`program.rs`（`Program`）、`input.rs`（`Role`、`read_fasta`、`Inputs`）、`describe.rs`（`describe`、`Description`、`Parameter`）、`search.rs`（`validate`、`Search`、`API_OWNED`）、`run.rs`（`run`、`run_collect`、`search_groups`、`Outcome`、`Report`、`HspRecord`）、`error.rs`（`Error`）。`local_blast.rs` は不変。

| 項目 | §3 の案 | 実装 | 理由 |
|---|---|---|---|
| 入力の名前 | argv の `-query` / `-subject` は `Error::Argv`、名前は無し | `Inputs { query, subject, query_name, subject_name }`（既定 `query` / `subject`、`Inputs::named`） | 報告は `-subject` の名前を印字する（outfmt 7 の `# Database: User specified sequence set (Input: NAME)`、outfmt 0 の `Database:`）。CLI と同じバイトを出すには名前が要る。`run` は argv の前に `-query NAME -subject NAME` を置いて CLI の parser に渡す |
| `-db` | `Subject::Db(path)` | 無し。argv の `-db` は `Error::Argv`（`API_OWNED`） | 5 program の `run_local` はどれも `-db` を断る（`-db (a database search) is not supported in the browser yet`）。`-db` の入口はエンジンの変更（別の段、LIB-1c として §5 に足す） |
| `-parse_deflines` | — | `Error::Argv`（`validate`） | adapter と同じ理由：レコードはオプションを知る前に読む（`read_fasta`） |
| `validate` | parser + option の検査 | parser と API の規則だけ。範囲・scoring・hit list の検査は `run` / `run_collect` が CLI の文言で返す | 検査は format と入力に依る（adapter は group ごとに検査している）。重複を作らない |
| `Outcome.formats: Vec<(String, Vec<u8>)>` | | `Outcome.reports: Vec<Report { outfmt, bytes }>`、要求した順 | 名前付き |
| `HspRecord.out6 / out0 / out0_subject` | | `spans: Vec<Option<Range<usize>>>` と `subject_headings: Vec<Option<Range<usize>>>`、`reports` と同じ順 | API は任意の `-outfmt`（`6 qseqid ...`、`7`、`10`、…）を取るので、format 番号で固定しない。npm の JSON はこれを ABI v2 の形に写せる |
| `run_collect` の失敗 | | どの検索の失敗も `Error::Search`（fail fast） | adapter の「hit list 0 で断られた検索を許す」規則は Web の 3 format 固定の都合。API では要求した format が CLI でも失敗する |
| HSP レコードの出どころ | | 最初の検索（tabular 群）。他の群の report には span を付けない | adapter の `link_format0` は Web の UI 用。API は report ごとの正しい span だけ返す |
| `describe` | adapter の JSON を移す | 型付き（`Description`、`Parameter`）。adapter の `defaults`（661 行の解決）は LIB-2 で | `defaults` は Web のフォーム用 |
| `search_groups` | adapter から移す | `api::search_groups(command, numbers) -> Vec<Vec<usize>>`（format 番号の列の位置で群を返す、program 非依存） | adapter は LIB-2 でこれを使う |
| 進捗バー（§3.4） | `ReportOutputs` に足すか feature | `BlastnArgs::progress_bar`（`#[arg(skip)]`）。CLI の `run` が `true` にし、`run_local` は触らない | `ReportOutputs` の構築箇所が 20 を超え（blastp の block、tblastn の stage E、adapter、ABI v1）、field 追加は engine の広い編集になる。feature は `parallel` と抱き合わせになる |
| エラーの型 | `Error { Argv, Input, Search }` | 同じ。`#[non_exhaustive]`。`Error::message()` | |

受け入れ条件（§6）の実装：`LOSAT/examples/library_search.rs`（`<program> <query.fa> <subject.fa> [options...]`、outfmt 6 を標準出力へ）、`LOSAT/tests/library_api.rs`（`common` suite：CLI とのバイト一致、同時 2 検索、例の stderr が空、`-db` の拒否）。

### 9.3 この PR で検証できたこと／残ること

- この PR の作業機は Rust 1.97.0（1.92.0 と wasm32 target を取得できない網）、NCBI BLAST+ 無し。quick tier（fmt、clippy、`cargo test --lib api::`、`cargo test --test common library_api::`、`cargo doc --no-deps`）はここで、`fast output regressions` と `rust` は CI で、V-ABI quick は `web.yml` で。
- 残り（保守者の機）：full tier の parity gate（出力バイト不変の確認。変更は `api/` の追加、`lib.rs` の `#[doc(hidden)]`、`BlastnArgs::progress_bar` の 3 行、`web_api.rs` の field 1 行）。

## 10. 段階 1：公開面を固める（2026-11-05、PR `feature/library-api-values`）

基点は `origin/dev` `4fa76ca`。もとの調査は保守者の文書「相同性検索ライブラリのAPI調査とLOSATへの提言」（2026-11-05）。

### 10.1 保守者の決定と順番

| # | 問い | 決定 |
|---|---|---|
| 1 | 中断 | 対象に入れる。バッチの境界で止める（§8） |
| 2 | binding の置き場 | `LOSAT/` の外で、`web/` と同じ境界規則（§4.3）。PD は段階 3 で書き、その PR は保守者が merge する |
| 3 | PyPI の名前 | `pylosat` |
| 4 | 最初の公開の `-db` | 含めない。`Inputs` は足せる形にする（§10.2）。LIB-1c は要望に応じて |
| 5 | 順番 | 公開の前に `pylosat` で確かめる（下の表）。LIB-2 は段階 3 の後で、公開の前提ではない |
| 6 | 段階 3 の検証先 | Python 版の gbdraw |

| 段階 | 内容 | 関門（次の段に進む条件） |
|---|---|---|
| 1 公開面を固める（この PR） | §10.2 | 型付きの値と `-outfmt` の列の一致テスト |
| 2 結果と実行制御 | `queries()` と `subjects()`、`run_batches`（逐次返却、進捗、境界での中断）、`Parameter` の `value_kind`、`-parse_deflines` | 分割実行と一括実行の全行一致 |
| 3 `pylosat` で検証（LIB-5 の先行版） | PyO3 で薄く包む（path 依存）、Arrow、GIL の解放、Ctrl-C、pickle。gbdraw とノートブックで API の不足を洗い出す | gbdraw が subprocess なしで動く |
| 4 公開（LIB-3、LIB-5、LIB-4） | crate 名 `losat`、crates.io と PyPI。npm は ABI v2 を包むだけなので独立に進められる | 公開面の差分検査（cargo-semver-checks、cargo-public-api）を CI に足す |

### 10.2 公開面の変更（§9.2 からの差分）

| 項目 | 実装 | 理由 |
|---|---|---|
| `HspRecord` の数 | `length`、`num_ident`、`num_positives`、`mismatch`（`i64`）、`gap_open`、`gaps`、`query_length`、`subject_length`（`usize`） | 数えるために報告を解析しなくてよくする。値は `-outfmt` の列（`length nident positive mismatch gapopen gaps qlen slen`）と同じ |
| `num_positives` | `PairwiseHit::positives`、無ければ `num_ident` | BLASTN の hit には無く、列は一致数を印字する（行列が無い） |
| `query_length` | record の長さ、`-query_loc` があれば範囲の長さ | `Hit::query_length` は program で違う（TBLASTX は 0）。`qlen` は `m_QueryRange.GetLength()`（`tabular.cpp:880-887`）で、範囲が record を越えてもそのまま。`slen` は `-subject_loc` でも record の長さ |
| 鎖 | `query_strand`、`subject_strand`：`Option<Strand>`、タンパク質は `None` | `sstrand` の規則（`tabular.cpp:1061-1067`）。1 文字の行は plus |
| 枠 | `query_frame`、`subject_frame`：翻訳する側だけ `Some` | 文書どおりにする（BLASTP の hit の `Some(1)` は写さない）。列の `qframe` / `sframe` は翻訳しない側にも数を印字する（BLASTN 1 / ±1、BLASTP 1 / 1、BLASTX の `sframe` 0、TBLASTN の `qframe` 0）が、それは報告の形で、鎖は `Strand` が持つ |
| 区間 | `HspRecord::query_range()`、`subject_range()`：plus 鎖の上の 0 始まり半開 | `Bio.Align`、PAF、BED の形。NCBI の `BlastSeg` の `offset` と `end` も同じ数え方 |
| 連結と組成補正 | `sum_n`、`comp_adjust_method` | E 値とスコアの読み方を変える値。`PairwiseHit` にある |
| 報告の中の位置 | `spans`、`subject_headings` を `Report` へ。`hsps` と同じ長さで、`HspRecord::index` で引く | HSP ごとの `Vec` 2 本をやめる。範囲は報告の中の位置なので報告に置く |
| 検索統計 | `Outcome::query_stats: Vec<Option<QueryStats>>`（`eff_space`、`hsp_len`、`lambda`、`kappa`、`entropy`）、`Inputs` の query の順。無効な query は `None` | NCBI の `CBlastAncillaryData`。E 値は検索全体の値なので、結果と一緒に返す |
| 統計の渡し方 | `FormatObserver::query_stats`（既定は何もしない）。5 program の XML の値を作る部分を関数にして、observer があるときだけ呼ぶ | `ReportOutputs` に field を足さない（§9.2 の進捗バーと同じ理由）。observer が無い CLI の経路は通らない |
| `run` | `(outfmt, writer)` の列と `diagnostics` を取る。hit list の大きさが違う format は `run_collect` と同じく別の検索に分ける | `ReportOutputs` は `PairwiseHit` を公開する（`HitsSink`）。`ReportOutputs` を取る `run` を使う呼び出し側は無かった |
| 隠すもの | `api::local_blast`（`ReportOutputs` 系と `FormatObserver`）と `cli` を `#[doc(hidden)]`。`api` の子モジュールは private（公開の path は `api::X` だけ）。`search_groups` は `pub(crate)`（adapter は LIB-2 まで自分のものを使う）。`FastaInputSource`、`ReaderConfig` は `api` から外す | 隠しモジュールの型が公開面に出ないようにする。後から公開するのは破壊的変更にならないが、逆はなる |
| `FastaRecord` | `#[non_exhaustive]`。`new`、`seq`、`title` と `local_id`、`title`、`sequence`、`warnings` だけが公開面。ほかは `#[doc(hidden)]` | engine 用のメソッドが `report` の型を連れてくる |
| `Inputs` | field を `pub(crate)` に。`-db` は別の constructor にする | 最初の公開の後に `-db` を足せる形にする（決定 4） |
| 足せる形 | `Report`、`Program`、`FastaRecord` に `#[non_exhaustive]`、`API_OWNED` は `&[&str]` | 公開の後の追加を破壊的変更にしない |
| format なしの呼び出し | `Error::Argv`（`no output format was requested`） | engine が断る（文言も同じ）。検索の前に揃えて断る |

engine の変更は 2 つで、どちらも CLI の出力に触れない：統計の受け渡し（上の表）と、BLASTX の field list つき tabular の行に observer の印を足したこと（`algorithm/blastx/tabular.rs`。`run_collect` の `spans` がその行に無かった。CLI では `probe` が `None`）。

残ること：

- BLASTX の `run_local` は outfmt 0、6、7 だけで、ほかは `not supported in the browser` と断る（段階 2 で広げた、§11.1）。
- `HspRecord` は `pident` を持たない。`100 * num_ident / length` で、印字は報告の列を読む。
- `Outcome` の型付きの結果だけが欲しい呼び出しでも、報告を 1 つ要求し整列行を描く（文書の「計算段階」。後回し）。

### 10.3 検証

- 関門のテストは `LOSAT/tests/library_api.rs`：5 program の 13 検索（`-query_loc`、`-subject_loc`、minus 鎖を含む）の全 HSP で、`length nident positive mismatch gapopen gaps qlen slen qstart qend sstart send qframe sframe sstrand score bitscore evalue` の列と `HspRecord` が一致し、`query_range` と `subject_range` が整列行の文字と一致する。統計は CLI の outfmt 5 の `Statistics` と query ごとに一致する（バッチに分かれる入力、BLASTP の block の経路、有効な context が無い query、`-query_loc` が飛ばす query を含む）。
- CLI の出力バイトは `fast output regressions`（CI）。full tier は保守者の機で（§9.3 と同じ）。

## 11. 段階 2：結果と実行制御（2026-11-05、PR `feature/library-api-batches`）

基点は段階 1 の `feature/library-api-values`。

### 11.1 公開面に足したもの

| 項目 | 実装 | 理由 |
|---|---|---|
| ビュー | `Outcome::queries()`、`Outcome::query(q_idx)`、`Batch::queries()` が `QueryHits`（`q_idx`、`stats`、`hsps`）を返し、`QueryHits::subjects()`、`QueryHits::subject(s_idx)` が `SubjectHits`（`q_idx`、`s_idx`、`hsps`）を返す。`hsps` の範囲でコピーしない | NCBI の `CSearchResultSet::operator[]` と `GetResults(qi, si)`。hit list は query の順で、query の中で subject ごとに隣接する（全ケースをテストで固定） |
| 逐次返却 | `run_batches(search, inputs) -> Batches`。`Batches` は `Iterator<Item = Result<Batch, Error>>` で、`next` が 1 バッチを検索する。`Batch` は `query_indices`（inputs の添字の範囲）、`hsps`、`query_stats`、`diagnostics` | 結果を 1 バッチ分ずつ持てばよい。NCBI の `for (; !input.End(); ...)` そのもの |
| 進捗と中断 | `Batch::query_indices` と `Batches::query_count()` が進捗で、次を求めなければ止まる（iterator を捨てる）。callback は無い | Python の generator にそのまま写る。Ctrl-C は `next` の間に効く（§10.1 の決定 1）。バッチの中は止まらない |
| 型付きのみ | `run_batches` の対象は `HspRecord` と `QueryStats` だけ。報告は無い | outfmt 0 と 7 は run 全体の頭と尾を持つ |
| オプションの型 | `Parameter::value_kind: ValueKind`（`Flag`、`Boolean`、`Integer`、`Real`、`Path`、`Text`）。parser が作る型から決める | NCBI の `CArgDescriptions::EType`。Python の引数と型スタブをここから作る |
| `-parse_deflines` | `read_inputs(search, query, subject)` が `(query, subject)` の records を返す。subject を先に、id の scope を共有して読む。`validate` は受け付ける（BLASTX を除く）。`read_fasta` の records で走らせると `Error::Argv` | scope は 1 組の入力のものなので、`Search` には持たせない。NCBI も subject を先に読む |
| BLASTX の format | `run_local` が outfmt 5、8、9、10、12、15、16、20 も書く（CLI の writer）。断るのは 11 と 13、14 だけで、文言から `in the browser` を除いた | 「各報告は CLI のバイトと同じ」の約束を BLASTX でも守る。LOSAT Web は adapter が 0、6、7 に絞る |

### 11.2 文書の見込みとの違い

文書は「`run_batches` はエンジンを変えず、バッチごとに `run_local` を呼ぶだけで作れる」としたが、3 点で足りなかった。

1. **run の警告がバッチごとに繰り返される。** subject の reader の警告、`Subject sequence contains no data`、オプションの警告（`-max_target_seqs 1` など）は `run_local` の頭で毎回書かれる。CLI は 1 回。`FormatObserver::queries_begin`（既定は何もしない）で run の分の終わりを教え、2 番目以降のバッチではそこまでを捨てる。
2. **BLASTN のバッチの大きさは固定でない。** `CBatchSizeMixer` が直前のバッチの extension 数から次の大きさを決める（`blastn_app.cpp:261-300`）。バッチが違えば lookup table の選び方も変わりうるので、CLI と同じ境界にする。`QueryBatchSizes`（`blast_engine/run.rs`）が同じ計算をし、`FormatObserver::batch_extensions` で各バッチの数を受け取る。入れ子の `run_local` が自分の mixer で切り直さないよう、`BlastnArgs::one_query_batch`（`#[arg(skip)]`、CLI は設定しない）で渡した records を 1 バッチにする。

3. **警告の `Query_N` は入力全体の番号。** BLASTN、BLASTP、TBLASTN、TBLASTX の警告は query の位置から番号を作るので、入れ子の `run_local` ではバッチの中の番号になる。各 program の引数に `first_query_ordinal`（`#[arg(skip)]`、CLI では 0。NCBI の `SetLocalIdCounterInitValue` に当たる）を足した。BLASTX の警告は record の id を使うので要らない。

どれも observer があるとき、または `run_batches` からのときだけ値が変わり、CLI の経路の出力は変わらない。

### 11.3 残ること

- バッチごとに `run_local` を 1 回呼ぶので、subject の準備もバッチごとにやり直す。文書の測定では 20 分割で 0.75 s が 0.87 s。`-db` の後で測り直す（再利用状態の項）。
- BLASTX の `-parse_deflines`：`run_local` が records の id を持ち越さない（`blastx_record`）。
- `Batches` は `Search` と records を借りる。`pylosat` の generator は所有する形を要するかもしれない（段階 3 で決める）。→ API は変えなかった（§12.1）。

### 11.4 検証

`LOSAT/tests/library_api.rs`：

- `run_batches` の `hsps`、`query_stats`、`diagnostics` を連結すると `run_collect` の値と全項目で一致する。5 program の 13 検索と、バッチに分かれるコピー入力（BLASTN は `-task blastn`、BLASTP は 2 ブロック分）で。
- run の警告は 5 program とも 1 回だけ（`-max_target_seqs 1` と残基のない subject）。
- BLASTN のバッチ数は、全 query の検索が engine から報告する数と一致する（`blastn`、`megablast`）。
- 全 HSP が `queries()` と `subjects()` のちょうど 1 か所に、hit list の順で現れる。
- `-parse_deflines` の報告（outfmt 6、0、7）が CLI とバイト一致する（BLASTN、BLASTP、TBLASTN、TBLASTX）。
- outfmt 5、8、9、10、12、15、16、20 の報告が 5 program で CLI とバイト一致する。

## 12. 段階 3：`pylosat` で確かめる（2026-11-05、PR `feature/pylosat`）

基点は段階 2 の `feature/library-api-batches`。使い方は [`python/README.md`](../python/README.md)。

### 12.1 形

| 項目 | 実装 | 理由 |
|---|---|---|
| 二層 | Rust は `run`、`Batches`、`HspTable`、`Hsp`、`describe`、`validate` と例外（`python/src/`）。ほかは Python（`pylosat/_api.py`） | dataclass なら等価と pickle を書かずに済み、公開の形を型つきの Python のまま読める |
| 引数 | `Search("blastn", evalue=1e-5)` は `describe` の表で `-evalue 1e-05` にして parser に渡す。`Search.from_argv` は語のまま | 検索条件の定義は argv 1 か所のまま。値の種類が違えば `TypeError`。`repr(float)` は同じ値に戻る |
| 型スタブ | `pylosat/_options.pyi` を `tools/generate_option_stubs.py` が同じ表から作る。`Search.__init__` は program ごとの overload に `Unpack[TypedDict]` | 表と stub のずれはテストが見る |
| 入力 | `Sequences.from_fasta`、`from_text`、`from_pairs`。bytes か path を持ち、検索の開始時に `read_inputs` で読む | reader は program と `-parse_deflines` で変わるので、作るときには読めない。`from_pairs` は FASTA の bytes にする |
| 表 | `HspTable` の中身は Arrow の `RecordBatch` 1 個。`Hsp` はその 1 行、slice と `queries()` はコピーしない | 保持を 1 か所にして、表、行、pickle を同じものから作る。`HspRecord` は `#[non_exhaustive]` で、外から作り直せない |
| 列名 | `-outfmt` の名前で、既定の報告の 12 列が先頭（`qaccver saccver pident ...`）。`Hsp` は `HspRecord` の名前 | 文書の方針どおり。文書の例の `qseqid` ではなく `qaccver`：既定の列がそれで、`-parse_deflines` のとき値が違う |
| Arrow | `arrow-array` の `ffi` で `__arrow_c_stream__` を直に作る。pyo3-arrow は使わない | 使うのは capsule 1 個で、pyo3-arrow は numpy、chrono-tz、arrow-cast を連れてくる。要求された schema は無視する（interface は best effort）。文字列は `LargeUtf8`（`qseq` の合計が 2 GiB を越えても切れない） |
| pickle | 表は Arrow IPC の bytes、`Hsp` は 1 行の表、ほかは dataclass | process pool で結果を返せる。`arrow-ipc` は `zstd` なし |
| `Batches` の所有 | binding が `self_cell` で `Search` と records を一緒に持つ。API は借りる形のまま | Rust の呼び出し側に `Arc` を強いる理由が無い。`unsafe` も書かずに済む |
| Ctrl-C | `Batches.__next__` と、報告なしの `run` の各バッチの前に `check_signals` | 決定 1。報告なしの `run` は `run_batches` で動かす（§11.4 で全項目一致）。報告ありは `run_collect` 1 回で、戻るまで効かない |
| GIL | 入力の読み込み、検索、表の作成、IPC を `Python::detach` の中で | 検索中も他の Python スレッドが動く |
| panic | `python/Cargo.toml` の release profile は LOSAT の写しで `panic = "unwind"` だけ違う | §4.3 |
| wheel | `abi3-py310` | wheel が Python の版ごとに要らない |

### 12.2 gbdraw とノートブックで見つかった API の不足

足したもの（`LOSAT/src/api`、規則は `LOSAT/` のまま）：

| 不足 | 足したもの | 根拠 |
|---|---|---|
| 行に id が無い。表が `qaccver` を持てない | `HspRecord::query_id`、`subject_id`（`Arc<str>`、record ごとに 1 回作る） | engine の `shown_id`、`subject_tabular_id`。列との一致は `check_columns` と `-parse_deflines` のテスト |
| `pident` は binding が計算できない | `HspRecord::percent_identity()`。式は `report/outfmt6.rs` の 1 か所 | `x_PrintPercentIdentical`（`tabular.hpp:437-442`） |
| ファイルは bytes にしてからで、`.gz` は呼び出し側が開けない | `read_input_file(role, path)`。CLI の開き方と文言 | `open_file_input`、`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` の例外 8 |
| エンジンの版が取れない（gbdraw の runtime 記録） | `api::VERSION` | `CBlastVersion`（`version.hpp`） |

足さなかったもの：

- 報告つきの検索は止まらない。`run_collect` は engine の 1 回の呼び出しで、バッチの境界を外に見せない。gbdraw は outfmt 6 の text が要るのでこの形になる。
- 1 query は 1 バッチなので、ゲノム 2 本の比較には境界が無い。gbdraw の検索はほぼこれで、subprocess なら kill できたものが、process の中では終わるまで待つことになる。NCBI の `SetInterruptCallback` に当たるものがエンジンに入れば解ける（文書の「後回し」のまま）。gbdraw が移るかどうかの判断材料になる。
- BLASTN の `-task` は `Parameter::choices` に出ない。parser が `String` で受け、engine が後で断るため（BLASTP は parser の list）。gbdraw の `losat_task_values` は `--help` を読んでいて、binding では検索の開始時のエラーになる。
- 警告は `diagnostics` の text のままで、深刻度も query の番号も無い（文書の「後回し」）。Python の `warnings` には流さない。

### 12.3 関門

- gbdraw が subprocess なしで動く：`docs/evidence/pylosat_gbdraw/`。gbdraw `033a3aa` に小さな patch を当て（gbdraw の repository は変えない）、`subprocess.Popen` を断った状態で linear と circular の 7 シナリオ（`losatn`、`tlosatx`、`losatp`）を実行し、同じ commit の `losat` を使う実行と SVG と raw の bytes が一致した。gbdraw 自身の T-CLI-07 の図の test は通る。T-CLI-09 の表は同じ commit の `losat` でも一致しない：表を作った v0.1.0 の遺伝暗号表は、コード 2、3、5、13、21 の `ATA` とコード 31 の `TGA` が `gc.prt` と違っていた（今の 27 表は `gc.prt` と同じ）。ヒトの配列の 3107 番の `N` をまたぐ行も v0.1.0 と違う。
- ノートブックの流れ（進捗つきの `batches`、`polars` と `pandas` と `pyarrow` への受け渡し、`queries()`、`spans`、pickle）を `python/tests` と手元の script で確認した。
- `python/tests`：5 program で、報告（outfmt 6、7、0、10、列指定）のバイトと `diagnostics` が `losat` と一致し、表の値が列と一致する。Ctrl-C、GIL、2 スレッド同時、process pool、pickle、Arrow の 3 consumer。
- `LOSAT/tests/library_api.rs`：`qaccver`、`saccver`、`pident` の列と全 HSP で一致（§10.3 の 13 検索と `-parse_deflines`）。
- `.github/workflows/python.yml`：`python/` と `LOSAT/src/api` の変更で走る。必須の check にはしていない。

### 12.4 残ること

- 公開（段階 4）：wheel のマトリクス、release workflow、`cargo-semver-checks` と `cargo-public-api`。
- gbdraw の製品判断：`pylosat` を runtime に入れるか、その順番と記録の形。
- `panic = "unwind"` の経路はテストしていない（公開の API から panic を起こす手段が無い）。release profile の build は通り、`python/tests` と gbdraw の比較はその build でも通る。
