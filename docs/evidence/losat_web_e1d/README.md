# LOSAT Web E1d（Session S05）ゲート記録

- 段階：E1d アダプタと ABI v2（[総合計画書](../../losat_web_gui_plan.md) §7 の S05、[指示書](../../losat_web_gui_sessions/session_s05_e1d_adapter_abi_v2.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S04 の後の `4482616a6`（エンジンは `458cbe110`）。エンジンの変更はコミット `edf825d36`、アダプタ・ABI・アプリの型・CI はコミット `20d26bdc3`、V-ABI の並列実行の道具（`tools/run_v_abi_parallel.py`、`v_abi.js` の `--engines`）と文書の補足はその次のコミット
- 実行記録：[`run-20260928T221453Z/`](run-20260928T221453Z/)。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)
- 判定：**完了条件を満たした**。独立監査は supported。監査の指摘のうち reactor のバイトを変えるもの（パスの置き換え、`validate` の文言、`register` の照合の試験など）は、reactor を作り直す S07 の最初の作業にした（下の「独立監査」）

## 変更の内容

### アダプタ（`web/adapter/`、新しい crate）

| ファイル | 内容 |
|---|---|
| `src/lib.rs` | ABI v2 の export（`docs/web/abi_v2.md` §4）。ABI v1 の `losat_web_*` はエンジンからそのまま link され、状態を共有しない（計画 TD-1） |
| `src/run.rs` | argv の解析（`validate`。エンジンの CLI の解析器をそのまま使い、`-out` と `-outfmt` を拒否する）と `run`：1 回の `run_local` で、program が対応するすべての形式、HSP レコード、警告を出す。観測者の出来事から `out6`・`out0`・`out0_subject` の範囲を記録する |
| `src/store.rs` | 登録した入力（`register`、`release`）。program の解析器（`bio::io::fasta`）で読み、同じバイトの索引の走査と ID と長さを照合し、食い違えば止める（計画 TD-8） |
| `src/scan.rs` | 索引の走査（`scan_*`）。`bio::io::fasta` 1.6.0 の規則を、任意の大きさの chunk で再現する |
| `src/describe.rs` | `describe`。エンジンの clap の定義から生成し、遺伝暗号の候補はエンジンの解析器に尋ねる |
| `src/emit.rs`、`src/json.rs` | host への出力（`losat_host.emit`、1 MiB ごと）、JSON |
| `tests/scan_properties.rs` | 走査の性質試験（下） |
| `tests/v_abi.js` | V-ABI（下） |
| `tools/check_build_identity.py` | ビルドの同一性（計画 TD-6） |
| `tools/build_reactors.py` | 2 つの reactor を作り、同一性の記録を残す |
| `tools/v_abi_cases.py` | V-ABI の検索の一覧（quick：CI、full：この記録） |
| `tools/run_v_abi_parallel.py` | full の V-ABI を、engine と program ごと（TBLASTX は検索ごと）の部分に分けて並列に実行し、結果をまとめる |

reactor の起動（`crt1-reactor.o`、`--entry=_initialize`）は、`LOSAT/build.rs` の `rustc-cdylib-link-arg` が依存先の cdylib にも渡ることで、エンジンの reactor と同じになる。指示書は `build = "../../LOSAT/build.rs"` での共有を挙げていたが、それを加えると `_initialize` が二重に定義されて link が失敗したので、アダプタは build script を持たない。

### エンジン（コミット `edf825d36`。出力のバイトは変えない）

| ファイル | 内容 |
|---|---|
| `LOSAT/src/api/local_blast.rs` | `FormatObserver` に `subject_begin` / `subject_end`（既定は何もしない）、`FormatProbe` に同じもの |
| `LOSAT/src/report/pairwise.rs` | BLASTP と TBLASTN の outfmt 0 で、subject の見出し（NCBI の `x_ShowAlnvecInfo` が subject の最初の HSP の前に書く defline と `Length=` の行、`showalign.cpp:3613-3632`）の範囲を、最初の HSP の番号で知らせる。BLASTX と共有する `write_subject_header` の本体と引数は変えず、呼出し側で知らせる（DW-10） |
| `LOSAT/tests/wasi_thread_host.js` | WASI 以外の host の関数を渡す `imports`。主の instance には渡した関数を、スレッド用の instance には呼ばれたら例外にする関数を渡す |
| `LOSAT/tests/run_local_support/mod.rs` | subject の見出しの検査（query の subject ごとに 1 つ、最初の HSP の番号を持ち、その HSP の節の始まりで終わる） |
| 各所のコメント | S04 の独立監査の指摘（TBLASTX の `run` のコメント、NCBI の引用の細部） |

### 文書とアプリ

- `docs/web/abi_v2.md`：実装に合わせて確定した（状態は version 2）。`out0_subject` を足した。`emit` は export を呼んだスレッドで呼ばれること、reactor の起動の出どころ、`scan` の各項目の意味、エラーの文言の出どころを書いた。
- `web/app/src/ports/engine.ts`：HSP レコードに `out0_subject`、`ProgramDescription` に遺伝暗号の候補。FakeEngine も合わせた。
- `.github/workflows/web.yml`：`adapter` の job（同一性の検査、アダプタの試験、2 つの reactor、quick の V-ABI）。
- S09 の指示書に、確定した ABI と V-ABI の実行方法を書き足した。

## 完了条件と結果

| 完了条件（計画 §7 の S05） | 結果 | 証拠 |
|---|---|---|
| V-ABI：BLASTP・TBLASTN・BLASTN・TBLASTX の、その時点で対応する全形式 × スレッド 1/2/4 が期待値と一致 | 通過。full：Gate A・Stage G の行列・BLASTP の manifest の回帰の case をまとめた 68 の検索と TBLASTX の直接の検索を、serial reactor（1 スレッド）と threaded reactor（1・2・4 スレッド）で実行した 272 件すべてで、program が対応する全形式（BLASTP・TBLASTN は 0/6/7、BLASTN は 6/7、TBLASTX は 6）と、BLASTP・TBLASTN の HSP レコード（範囲 `out6`・`out0`・`out0_subject` を含む）と診断が、同じ commit のネイティブの CLI と一致した。凍結ハッシュは 488 件中 484 件が一致し、一致しない 4 件は S02 からの既知の `Sakai.MG1655.megablast` の outfmt 7（S02 の基準も同じく一致しない）。quick（CI と同じ 13 の検索）：52 件すべてが一致 | `run-20260928T221453Z/v-abi-full/`、`run-20260928T221453Z/v-abi-quick/` |
| ビルドの同一性の検査が通る | 通過。`[profile.release]`、Wasm の rustflags、2 つの `Cargo.lock` に共通の 144 の依存の版が一致し、アダプタは build script を持たない。threaded の reactor の共有メモリの最大値は、認証済みの threaded reactor と同じ 16384 ページ（1 GiB、計画 TD-7） | `run-20260928T221453Z/reactors/artifacts.json`、`losat-web-{serial,threads}.json` |
| `scan` の性質試験 | 通過。生成した 2,000 の入力（エラーになる入力 134、行長が揃った記録 2,068、チェックポイントの記録 2,143、65,536 × 2 を超える記録 479、CRLF を含む入力 1,090、非 ASCII を含む入力 1,495）と 17 の境界の場合で、走査の結果が `bio::io::fasta` の結果（記録、ID、長さ、エラー）と、すべての残基の位置を記録する参照の走査（位置、行の配置、残基の内訳）に一致し、1 バイト・1〜7 バイト・1〜4096 バイトの chunk で走査しても同じ | `run-20260928T221453Z/adapter-test.log` |
| `register` での照合 | `register` が走査と解析器の結果を照合し、食い違えば止める（`src/store.rs`）。V-ABI が `register` と chunk ごとの `scan` の一致を reactor で確かめる | `tests/v_abi.js` |
| エンジンの変更の既存のゲート | 変更後の 236 件が S02 の基準と一致（差 0）。v1 の WASI の検査が通過（409 件、形式の失敗 0。変更した `wasi_thread_host.js` で実行）。v1 の reactor の記録（threaded 344 件、serial 25 件）と `v1_requests.js` の応答が S04 と一致 | `run-20260928T221453Z/capture-compare.txt`、`wasm-threading.log`、`v1-reactor-records-compare.txt` |
| V-PERF（エンジンの変更） | 通過。変更した outfmt 0 の経路（`blastp-fmt0`：SicyWSV.faa × PajaWSV.faa、`tblastn-fmt0`：TLOSAN の Stage G の benchmark）を、変更前（S04 の成果物）と変更後で 1 回ずつ交互に測った（1 回の暖機と 3 回の計測、`measure_perf.py`）。6 件すべてが ×0.983〜×1.042 で、出力は同じ。独立監査の指摘（m6）で測った。変更は観測者があるときだけ flush と通知を行い（`pairwise.rs` の `if let Some(probe)`）、観測者の無い CLI では subject ごとに最初の HSP の番号を 1 回引くだけである | `run-20260928T221453Z/perf-1.json`、`perf-check-1.txt` |
| `cargo fmt --check`、`clippy -D warnings`、`cargo test` | 通過。エンジンは 4 つの構成の clippy と 796 件の試験、アダプタは 3 つの構成（ネイティブ、`wasm32-wasip1`、`wasm32-wasip1-threads --features threads`）の clippy と 4 件の試験 | `run-20260928T221453Z/fmt.log`、`clippy.log`、`adapter-clippy.log`、`cargo-test.log`、`adapter-test.log` |
| `web/app/src/ports/engine.ts` の型が ABI と一致 | `out0_subject` と遺伝暗号の候補を足した。`npm run check` が通過 | `web/app` |

成果物のハッシュと道具の版は `run-20260928T221453Z/binaries.txt` にある。コミットの木から作り直した 2 つのネイティブ、4 つの WASI の成果物と 2 つのアダプタの reactor は、ゲートに使った成果物と同じハッシュになった。アダプタの reactor のバイトは、ビルドしたディレクトリの絶対パスに依存する（LOSAT を path 依存として使うので、Cargo が記号の hash にその場所を含める）。セッションの規則の worktree のパスでビルドすれば同じになる（`docs/web/abi_v2.md` §2）。同じハッシュの複製を `s05-final/` に残した。

## 独立監査

読み取り専用の独立監査（役割 `ncbi_parity_auditor`、コードを変えない別のエージェント）を、full の V-ABI の最後の部分が走っている間に受けた。結論は **supported**（S05 は完了条件を満たす。V-ABI は最後の部分の結果待ち。その部分はその後に通過した）。

監査が確かめたこと：エンジンの変更の前後で 236 件の出力・stderr・終了コードが同じ（BLASTX の 12 件を含む）。flush と通知はすべて観測者があるときだけで、BLASTX と共有する `write_subject_header` の本体と引数は変わらない。subject の見出しの範囲が NCBI の `showalign.cpp:3613-3632` と一致する（reactor で BLASTP 7,093、TBLASTN 91 の見出しを確かめた）。reactor の export・import・形式・`-out`/`-outfmt` の拒否・alloc の対応・`emit` の呼ばれるスレッド・`register` が食い違えば止まること。`scan` が `bio::io::fasta` 1.6.0 と一致すること（監査が独自に 200,040 の入力で差分の試験を行い、差 0。先頭の空白、空白だけの行、NEL、U+2028、BOM、NUL、単独の CR、行の途中の `>`、不正な UTF-8 を含む）。V-ABI の比較の網羅と、並列の部分の結果がすべてまとめられること。同一性の検査の再実行。NCBI の参照と、NCBI の実行ファイルを実行時に使っていないこと。

指摘と扱い：

| | 指摘 | 扱い |
|---|---|---|
| M1（重大） | reactor が作った場所に依存する理由の説明が不十分。名前の区間だけでなく、コードとデータの区間も違う。panic の位置に checkout の絶対パス（約 75 か所）と `CARGO_HOME` のパス（32〜44 か所）が入っている。公開すれば、ビルドした機械のパスが漏れる | `docs/web/abi_v2.md` §2 と `binaries.txt` の説明を直した。`build_reactors.py` の同一性の記録に checkout と `CARGO_HOME` を足した。パスの置き換え（`--remap-path-prefix`）を計画の TD-11 とし、reactor を作り直す S07 の最初の作業にした。S16 の V-PRIV で公開物にパスが無いことを確かめる |
| m1 | `validate` のエラーの文言が CLI と一部違う（binary の名前が `losat`、使い方の行に `-outfmt` が出る、未知の program の文言） | `abi_v2.md` §4 に違いを書いた。binary の名前は S07、`-outfmt 6` の挿入をやめるのは全 program が outfmt 0 を受け付ける S08（それぞれの指示書に書いた） |
| m2 | `register` の照合が食い違ったときに止まることの試験が無い（指示書 4） | 照合の処理を関数に分けて単体試験を足すのは reactor のバイトを変えるので、S07 の最初の作業にした。止まることは監査がコードで確かめた |
| m3 | BLASTN と TBLASTX が HSP レコードを出さないことが ABI の文書に無い | `abi_v2.md` §4・§5・§8 に書いた |
| m4 | アプリの port に stream 3（診断）の受け口が無い | `RunSink.diagnostics`、`DataGateway.readDiagnostics` と Memory の実装、単体試験を足した（`npm run check`、`npm run e2e` 通過） |
| m5 | 同一性の検査が、`build` の鍵の無い `build.rs`、親のディレクトリと `CARGO_HOME` の Cargo の設定、環境変数を見ていない | 検査に足した。再実行して通過（`reactors/build-identity-recheck.json`） |
| m6 | エンジンの変更の V-PERF を測っていない（セッションの規則 4） | 測った（上の表と「性能の計測」） |
| 注記 | 走査の鍵の空白（`abi_v2.md` の記述）、`-num_threads` を `validate` が通すこと、WASI の import の一覧、`RangeRecorder` の release での検査、並列の道具の古い結果と成果物のハッシュ、TBLASTX の wasm-threads だけの経路の再記録、S05 と S06 の変更の分け方、`verification_cells.tsv` の升目 | 文書で済むものは直した（`abi_v2.md` §3・§9）。`-num_threads`、`RangeRecorder`、並列の道具は S07 の最初の作業にした（この記録の full の実行は空の `parts/` から始めたので古い結果は無い）。コミットは S05 と S06 で分けた。升目は埋めた |

**保守者への報告（判断は不要）**：認証済みの LOSAT の Wasm の成果物（v1 の command と reactor）にも、ビルドした機械の `CARGO_HOME` のパス（`/home/<user>/.cargo/registry/...`）が入っている（`losat-serial-reactor.wasm` で 32 か所。以前から）。アプリが配信するのは v2 のアダプタの reactor だけで、それは TD-11 で扱う。

## 実行の方法（再現）

```bash
# ビルドの同一性、アダプタの試験、2 つの reactor（同一性の記録つき）
python3 web/adapter/tools/check_build_identity.py
(cd web/adapter && cargo +1.92.0 test --locked -- --nocapture)
RUSTUP_TOOLCHAIN=1.92.0 python3 web/adapter/tools/build_reactors.py --target-dir <dir> --output-dir <reactors>
# V-ABI（full は Gate A の入力の綴りを使うので、docs/evidence/losat_web_e1a/README.md の lexical root が要る）
python3 web/adapter/tools/v_abi_cases.py --suite full --out <cases.json>
python3 web/adapter/tools/run_v_abi_parallel.py --native <LOSAT> --serial <reactors>/losat-web-serial.wasm \
  --threads <reactors>/losat-web-threads.wasm --cases <cases.json> --out <dir> --jobs 12
# quick（CI と同じ）は 1 つの v_abi.js で
node web/adapter/tests/v_abi.js --native <LOSAT> --serial <reactors>/losat-web-serial.wasm \
  --threads <reactors>/losat-web-threads.wasm --cases <cases.json> --out <dir>
```

エンジンの変更のゲート（出力の捕捉、v1 の WASI の検査、v1 の比較、性能）は S02 の記録の「実行の方法」と同じ。変更前の成果物は S04 の成果物（`s04-final/`）。

## 引き継ぎ（S06 以降へ）

- S07 の最初の作業：独立監査の指摘のうち reactor のバイトを変えるもの（計画 TD-11 のパスの置き換え、`validate` の binary の名前と `-num_threads` の拒否、`register` の照合の試験、`RangeRecorder` の release での検査、並列の道具の古い結果と成果物のハッシュ）。S07 の指示書に書いた。

- ABI は version 2 で確定した（`docs/web/abi_v2.md`）。S09 の指示書に、ブラウザの host が知るべき点（threaded の import、`emit` の扱い、V-ABI の実行方法）を書き足した。
- BLASTN と TBLASTX は、S07・S08 で `PairwiseHit` を作るまで HSP レコード（stream 1）を出さない。outfmt 0 を足すときは、観測者の `subject_begin` / `subject_end` も呼ぶ（BLASTP と TBLASTN の `pairwise.rs` の書き方）。
- BLASTX は SX で ABI v2 に加える（`Program::parse`、`describe` の形式、parser kind 1 の `scan`）。
- `v_abi.js` の case は `tools/v_abi_cases.py` が作る。新しい形式を足したら `FORMATS`（両方のファイル）を直す。
