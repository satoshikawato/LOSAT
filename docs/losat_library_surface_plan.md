# LOSAT ライブラリ化計画（公開面の設計メモ）

状態：**設計メモ（LIB-0）。実装は v0.2.0 の cut の後に始める（保守者 2026-10-31 「OK」）。** v0.2.0 自体は crates.io に出さない（release の決定 5、2026-10-25）。作成 2026-10-31。

| 項目 | 内容 |
|---|---|
| 目的 | LOSAT を Rust crate（crates.io）、npm、PyPI から使えるようにするときの**公開面**（public API surface）を先に決める |
| 対象 | `LOSAT/`（エンジン crate）の公開面、`web/adapter/` との分担、配布の 3 チャネルの形と順番 |
| 対象外 | wasm-bindgen、pyo3、ネイティブ C ABI、`psiblast` と `makedb` の API 化、検索の中断（§8） |
| 基点 | `origin/dev` `4b2588e01` で確認し、`61fe55d2d` で関係するファイルに差分が無いことを確かめた（2026-10-31） |
| 規則 | ルートの [`AGENTS.md`](../AGENTS.md)（エンジンの変更は NCBI の権威と parity gate）、[`web/AGENTS.md`](../web/AGENTS.md)、[`PD-LOSAT-WEB-APP-BOUNDARY`](product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md)、[`PD-PURE-RUST-RUNTIME-AUTHORITY`](product_decisions/PD-PURE-RUST-RUNTIME-AUTHORITY.md) |
| 記録 | `losat-baselines/library-20261031/`（DECISIONS、STATE、次のセッションの指示） |

---

## 0. 結論

1. 三つのチャネルの前提は同じで、**公開面を `losat::api` と `losat::cli` の 2 モジュールに絞る**ことである。今の lib は内部モジュールが全部 `pub` で（§1.2）、このまま publish すると docs.rs が内部実装の目録になり、互換の約束ができない。
2. 順番は **crates.io → npm → PyPI**。
   - crates.io は、Cargo.toml のメタデータが既に揃っていて、`cargo package --list` が通る（§1.1）。公開面を絞る LIB-1 の後なら、残りは名前とメタデータだけである。
   - npm は **wasm-bindgen を使わず**、既にある 2 つの WASI リアクタと ABI v2 のホスト側 TypeScript をそのまま包む（§4.2）。
   - PyPI は **maturin の bin wheel**（CLI を入れるだけ）から始め、pyo3 の binding は作らない（§4.3）。
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
- `args` は clap 派生の `BlastnArgs` 等。`records` は **LOSAT 自前の** `blastinput::fasta_reader::FastaRecord`（`local_id`、`title: Vec<u8>`、`sequence`、`warnings`、`source`）で、reader は `FastaInputSource::from_bytes(bytes, ReaderConfig::query / subject(...))`。CLI と adapter（`store.rs`）が同じ reader を使う。rust-bio（`bio`）は ABI v1（`v1_bio.rs`）と `bio::alphabets::dna`（translation、redo_alignment）に残るだけで、公開面の型には出ていない。
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
- **P4 公開面は `api` と `cli` だけ。** 他のモジュールは `#[doc(hidden)]` を付け、semver の保証外と README に書く。`pub(crate)` 化は別の線（dead-code 警告が今は `allow(warnings)` で隠れていて、一斉に出る）。
- **P5 ライブラリは stdout / stderr に触れず、`process::exit` しない。** 出力は呼び出し側の writer だけに書く。
- **P6 新しい BLAST の動作を入れない。** 4 program の façade は既存の入口を包むだけ。blastx の入口（LIB-1b）は `native::run` を `run_local` の形に分けるもので、動作は変えない。`AGENTS.md` の規則（NCBI の権威、参照コメント、parity gate）はそのまま適用する。参照は NCBI の API 層（`c++/src/algo/blast/api/`、`local_blast`、`blast_results`）で、`api/mod.rs` の既存の参照と同じ。
- **P7 版。** crate の version は CLI の version と同じ。0.x の間は `api` の破壊的変更で minor を上げる。ABI v2 §10 と同じ考え方。

## 3. 公開面（案）

最終の signature は LIB-1 の PR の review で決める。ここでは形と境界を決める。

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
- 中断（cancel）は API に入れない。LOSAT Web と同じく、ホストが worker / process を止める（§8）。

### 3.2 `losat::cli`

今のまま公開する（Web adapter が既に使う）。安定とみなすのは `Cli`、`Commands`、`try_parse_from`、`render_message`、`ncbi_version_request`。`main.rs` 用の `exit_on_native_error`、`report_standard_output`、`ReportStream`、`NativeError`、`inaccessible`、`command_line_bytes` は `#[doc(hidden)]`。

### 3.3 argv の扱い

ABI v2 §7 を基にするが、入力の渡し方が違う。

- `-outfmt` と `-out` は argv で受け取らず、`run` / `run_collect` の引数で指定する。argv にあれば `Error::Argv`（ABI v2 と同じ規則。adapter の文言を流用。保守者委任、§7 Q5）。
- ABI v2 は `-query <name>` / `-subject <name>` に登録した入力の名前を要求するが、API は `Inputs` でレコードを直接渡す。`Inputs` にレコードを渡すときに argv に `-query` / `-subject` があれば `Error::Argv`。`Subject::Db` のときの `-db` は argv に書く（path は呼び出し側のファイルシステム）。LIB-2 で adapter は、名前 → `Inputs` の対応を自分の store で持つ。
- `-num_threads`、`-task`、`-evalue` などの検索条件は argv のまま。既定値は CLI と同じ。

### 3.4 隠すもの

- `lib.rs` の `api` と `cli` 以外の `pub mod` に `#[doc(hidden)]` を付ける。項目は今までどおり使えるが、docs.rs に出ず、semver の保証外（README と crate の doc に書く）。`api` が再 export する `FastaRecord` 等は `blastinput` の中にあるので、再 export の側だけを文書化する。
- 進捗バーは `ReportOutputs` に `progress: Option<...>` を足して CLI（`main.rs`）だけが渡すか、indicatif を optional dependency にして CLI の feature だけで有効にする（LIB-1 で決める）。
- `anyhow::Error` は `api::Error` に包む。`clap` は argv が API なので必須 dependency のまま。

## 4. 配布の 3 チャネル

### 4.1 crates.io（LIB-3）

| 項目 | 案 |
|---|---|
| crate 名 | **`losat`（小文字）に改名**（保守者委任、§7 Q1）。bioconda / PyPI / npm と同じ綴りになる。crates.io は大文字小文字を区別せずに名前を予約するので、`LOSAT` で publish しても `losat` は使えないし逆も同じ。`[[bin]] name = "losat"` は今のまま。`LOSAT::` の path を使うファイルは 49（`src` 1、`tests` 43、`web/adapter` 5）で、機械的な置換（Sonnet の仕事） |
| メタデータ | `rust-version` あり。`readme = "../README.md"` は `cargo package` が `README.md` として入れることを確認済み。死んだ `exclude` の `"/matrix"` を消す。`[package.metadata.docs.rs]`（既定 feature、Linux）を足す |
| feature | 今の 6 つはそのまま。進捗バーは §3.4 |
| gate | `cargo publish --dry-run --locked`（`release-readiness.yml` に足す）、`cargo doc --no-deps` が `api` と `cli` だけを示す、tarball の大きさ（`cargo package` の出力） |
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

- 需要は実質 gbdraw で、gbdraw は subprocess + CLI の NCBI 互換で既に動き、生の出力をキャッシュしている。pyo3 の利点（process 起動なし、メモリ内 I/O、構造化した結果）はこの用途では小さく、代わりに wheel のマトリクス（manylinux x86_64 / aarch64、macOS arm64 / x86_64、Windows x64）と bioconda と二重の配布が増える。
- **maturin の bin wheel**（`[tool.maturin] bindings = "bin"`。`pip install losat` で CLI `losat` が入るだけ、API 設計なし）から始める（保守者委任、§7 Q3）。gbdraw の解決順（PATH の `losat`）はそのまま満たす。
- pyo3 の binding は作らない。必要になるのは、メモリ内の結果を Python で直接使う利用者が現れたときで、そのときも §3 の `api` の上に薄く作る（Rust の形を 2 言語分凍結することになるので、`api` が安定してから）。
- wheel は release のバイナリと同じ build（`release-readiness.yml` / release workflow）から作る。詳細は LIB-5。

## 5. 実装の段（cut 後、PR はすべて `dev` へ）

| 段 | 内容 | 触る場所 | gate |
|---|---|---|---|
| LIB-0 | この計画（このメモ） | `docs/` | — |
| LIB-1 | façade（4 program）: `api::{read_fasta, Program, describe, validate, Search, Inputs, run, run_collect, HspRecord, Error}` と `search_groups` の移設、`#[doc(hidden)]`、進捗バーの切り離し、`examples/library_search.rs` | `LOSAT/src/lib.rs`、`LOSAT/src/api/`、`LOSAT/src/main.rs`（進捗バー） | 既存の parity gate 全段（出力バイト不変）、adapter の V-ABI quick、`cargo doc --no-deps`、例の出力 = CLI の出力、同時 2 検索（§6） |
| LIB-1b | blastx の入口: `native::run` を `run_local(args, query, subject, outputs)` の形に分ける（動作は不変）。ABI v2 の「BLASTX joins in session SX」と同じ作業なので、Web の線と共有する | `LOSAT/src/algorithm/blastx/` | blastx の parity gate 全段、ABI v2 の blastx（SX の V-ABI） |
| LIB-2 | adapter を `api` の上に載せ替え: program 非依存の部分（`run.rs` の parse / validate / `search_groups` / HSP の収集、`describe.rs`）を `LOSAT/src/api` へ移し、adapter は WASI の export、store（名前 → `Inputs`）、scan、emit だけにする | `web/adapter/`、`LOSAT/src/api/` | V-ABI full（定義は `docs/losat_web_gui_plan.md` と `web/adapter/tools/v_abi_cases.py`）、同一性の検査（TD-6） |
| LIB-3 | crate 名、メタデータ、`cargo publish --dry-run` の gate、README の文、changelog | `LOSAT/Cargo.toml`、49 ファイルの path、`.github/workflows/release-readiness.yml`、`README.md` | `cargo publish --dry-run --locked`、CI 全緑 |
| LIB-4 | npm package | `web/`（切り出し方は LIB-4 で決める）、`.github/workflows/web.yml` | V-ABI quick を package の build 物で、`npm pack` の内容の確認 |
| LIB-5 | PyPI の bin wheel | `packaging/`、release workflow | wheel から `losat -version` と fixture 1 件 |

LIB-1、LIB-1b、LIB-2 は `LOSAT/src` を触るので `AGENTS.md` の規則の対象（NCBI の参照コメント、parity gate）だが、動作の変更は無い。LIB-1 が設計の判断を含む唯一の段で、残りは機械的。

## 6. 受け入れ条件

- CLI の出力バイトが不変（`fast output regressions` と skill `losat-gates` の full tier）。
- Web の V-ABI（LIB-1 は quick、LIB-2 は full）。
- `cargo doc --no-deps` が `api` と `cli`（と `#[doc(hidden)]` でない再 export）だけを示す。
- `examples/library_search.rs`: 2 つのメモリ上のレコードで blastp を走らせ、`run_collect` の outfmt 6 と outfmt 0 のバイト列が、同じ fixture の CLI の出力と一致する（`tests/` の 1 ケースにする）。stderr に何も出ない（subprocess で確かめる）。
- 同時 2 検索: 既定 build で 2 スレッドから `run_collect` を同時に呼び、1 検索ずつの出力と一致する（`tests/` の 1 ケース）。
- `cargo package --list` と `cargo publish --dry-run --locked` が通る。

## 7. 保守者への問い（委任した推奨と、違う選択肢）

常設指示（2026-09-29、2026-10-07 再確認）により、推奨を採って進め、ここに記録する。違う選択をするときは保守者が言う。

| # | 問い | 委任した推奨 | 他の選択肢 |
|---|---|---|---|
| Q1 | crate 名 | `losat`（小文字）に改名 | `LOSAT` のまま（code の変更なし、綴りが他のチャネルと違う） |
| Q2 | npm の名前 | スコープなしの `losat` | `@satoshikawato/losat` など（org / user が要る） |
| Q3 | Python | maturin の bin wheel だけ | pyo3 の binding を今作る |
| Q4 | 時期と版 | cut（11-01）後に LIB-1 を dev で始め、`0.3.0` で出す | `0.2.x` の patch に入れる |
| Q5 | argv の `-outfmt` / `-out` | ABI v2 と同じく断る | 受け取って既定にする |
| Q6 | blastx | LIB-1 は 4 program、blastx は LIB-1b で入口を足してから | LIB-1 に含める（LIB-1 が blastx の engine の変更を抱える） |

## 8. 対象外と決めないこと

- wasm-bindgen（§4.2）、pyo3（§4.3）。
- ネイティブの C ABI（`cdylib` は wasm のために既にあるが、C のヘッダと FFI は作らない）。
- `psiblast` と `makedb` の API 化（CLI のみ。`blastdbcmd` も同じ）。
- 検索の中断（ホストが止める）。部分結果のストリーム（`run` の writer は既にストリームだが、`run_collect` は集めて返す）。
- 内部モジュールの `pub(crate)` 化（別の線。`#[doc(hidden)]` で足りる間はやらない）。
