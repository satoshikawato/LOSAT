# Wasm threading 監査指摘の実装計画

> 2026-09-14: N 本の専用子 worker + 待機 caller というスレッド数の契約は、ユーザー指定により[親を含む合計 N 本](wasm_total_threads_20260914.md)へ変更する。以下は9月13日時点の履歴。

作成日: 2026-09-13。状態: **実装済み・長時間回帰はユーザー指示で後日に延期**。

対象は [2026-09-13 の監査](wasm_threading_audit_20260913.md) の F01〜F14。検索結果の正しさ、command/library の起動、反復実行、スレッド指定、診断の順に契約を成立させ、その後に性能を測定する。

本書作成時の HEAD は `78c7556da6fa18b44447c3a3fa73774a8ee5956f`。作業ツリーには既存の変更が多数ある。本書は修正済み・認証済みという主張ではなく、監査後のソース確認を反映した実装計画である。実装開始時には HEAD に加えて対象ソース・依存・artifact の hash を取り直す。

追加承認（2026-09-13）: ユーザーの指示により、検証で見つかった BLASTP outfmt 0/7/custom fields の既存書式差も修正対象に含めた。元の sink 修正に限定する制約を、この追加修正の禁止としては扱わない。

終了指示（2026-09-13）: 実装と完了済み検証を確認して commit し、長時間回帰は後日実行する。固定回帰は 99/145 ケース成功、46 件未完了として保持する。[実装・検証記録](wasm_threading_remediation_results_20260913.md)を参照。

## 1. 完了時に満たす契約

1. TBLASTX の複数 subject 検索は、NCBI の通常検索と同じ全 subject 統計を使う。統計を参照する cutoff、linking、最終 E-value まで同じ入力を渡す。
2. CLI と direct API は、直列・並列のどちらでも指定された出力先へ同じ検索結果を返す。
3. command と direct API 用 reactor は、ビルド出力・配布名・初期化方法を区別できる。別 artifact の上書きや取り違えを検出する。
4. 同一 reactor instance で `1→2→4→2→1` の検索、プログラム切り替え、回復可能な失敗後の再検索が成立する。
5. 有効な検索の明示的な `N > 1` 指定は、N 本の計算プールを構築して実行するか、理由を示して失敗する。入力サイズや CPU 推定値による無警告の縮小を行わない。
6. プール数、host worker 数、段階別の仕事量を区別して記録し、要求契約の検証に使える。
7. 修正経路を実際の Wasm artifact で CI 検証する。性能結果は出力・実行契約を満たす測定だけから作る。

「N 本」は検索に利用する計算プールのサイズであり、全段階で N 本が同時稼働する保証ではない。仕事が 1 件しかない段階を、見かけの稼働率のために新しい方法で分割しない。

## 2. 維持する境界と旧計画との関係

- [AGENTS.md](../AGENTS.md) を最上位のリポジトリ規則とし、実装時は `verify-ncbi-parity-and-speed` に従う。NCBI はソースと比較 oracle に限定し、runtime・build・fallback に組み込まない。
- コード変更には、対応する NCBI のファイル・行番号・snippet を直上コメントとして付ける。WASI/Rayon 固有の実装理由は、NCBI が定義する検索・順序・thread lifecycle と区別し、依存ソースの根拠も併記する。無関係な CPU 制限の snippet を TLS 初期化の根拠にしない。
- 浮動小数点精度、演算順序、候補集合、HSP 構築、pruning、reduction、出力順序を維持する。出力のソート・正規化・許容誤差で不一致を合格にしない。
- TBLASTX local subject の非既定 `db_gencode` は承認済みの挙動を維持する。例外を統計や順序の差へ広げない。
- [PD-NCBI-PLATFORM-VARIANCE](product_decisions/PD-NCBI-PLATFORM-VARIANCE.md) の Gate A（凍結 PR 5 LOSAT bytes）と Gate B（登録済み native NCBI fingerprint）を維持する。修正で既存 Gate A に差が出た場合、原因・影響 fixture を記録し、既存認証は失敗のまま扱う。本計画は期待値や authority の自動更新を許可しない。
- `wasm32-wasip1` は serial のままとする。threaded command/reactor は `wasm32-wasip1-threads` と `wasm-threads` を使い、feature 名だけで実行能力があるとみなさない。
- 同一 instance への API 呼び出しは直列とする。複数 caller による同時検索、新しい検索オプション、汎用 scheduler、新しい scan 分割方式は本修正に含めない。
- [旧 current-thread 計画](blastp_wasm_current_thread_pool_plan.md) の caller 登録と N−1 worker 前提、および [旧改善計画](wasm_multithreading_improvement_plan.md) の小入力直列化・BLASTP 固有 cap は、本計画の修正方針へ置き換える対象とする。旧文書は履歴として残し、実装時に冒頭へ後継計画の案内を付ける。
- 明示指定の維持は監査に記載された LOSAT の実行契約である。NCBI CLI は CPU 制限や local subject の指定無視を警告付きで行うため、「NCBI CLI が常に要求 N 本を使う」と説明しない。

## 3. 指摘と実装段階の対応

| 指摘 | 主な修正 | 段階 | 完了を示す証拠 |
|---|---|---|---|
| F01 | 全 subject 統計とその利用箇所の修正 | S01 | 複数 subject の通常 NCBI raw 出力一致、統計・cutoff の比較 |
| F02 | BLASTN の出力先選択を共有 | S01 | direct API 並列出力の成功・bytes 一致 |
| F03 | current-thread 登録を使わず、検索単位で pool 終了待ち | S03 | 同一 instance 反復、失敗後の再実行、worker 回収 |
| F04 | 入力しきい値による要求数無効化を解消 | S03 | 小入力・従来の逆転 fixture で n2/n4 の pool 一致 |
| F05 | BLASTP の暗黙 cap を除去 | S03 | host cap 未設定、明示上限、native の要求数検証 |
| F06 | serial build の N > 1 を明示拒否 | S03 | CLI 非ゼロ終了、API error、worker 起動なし |
| F07 | subject traversal と linking の実行条件を分離 | S04 | 1 subject・複数 frame group で linking が指定 pool を使う |
| F08 | 実 pool・起動確認に基づく診断と集計 | S03/S04 | BLASTN DP の診断・host 記録一致、集計の矛盾検出 |
| F09 | Rayon 上限の事前検証と pool 構築後検証 | S03 | 上限超過の事前失敗、黙った縮小なし |
| F10 | 診断 scan の名称・表示を実処理に合わせる | S04 | 逐次 scan と linking の診断を区別 |
| F11 | API error chain を保持 | S01 | 外側 context と原因の両方を取得 |
| F12 | 実 artifact による CI・回帰 gate | S05 | target 別・反復・失敗注入・出力の実行結果 |
| F13 | reactor の初期化と終了処理を成立させる | S02 | 正規初期化後の最初の並列検索成功、起動経路の証拠 |
| F14 | command/library の build・packaging 分離 | S02 | build 順序を変えても別 artifact を保持、誤種別を拒否 |

推奨実装順序は `S00 → S01 → S02 → S03 → S04 → S05 → S06`。S01 の統計修正と I/O 修正は独立した変更単位に分けられる。S02 への着手に F02 の正式な API 合格は要求せず、I/O 修正を実装した上で、S02 の runtime 初期化確立時に両方を検証する。これにより S01/S02 間の検証待ちを循環させない。S03 の Rayon lifecycle は正規 reactor で判定する。

各段階に必要な回帰テストを同時に追加する。S05 はテスト作成の開始点ではなく、既存・追加テストを CI と認証へ接続する段階である。

## 4. 実装段階

### S00: 現行差分・再現条件・ソース対応を固定する

**作業**

1. working tree の状態、HEAD、対象ソースと runner の hash、Cargo.lock、Rust/LLVM、Node/V8、NCBI の版・実体を記録する。既存変更は復元・整形・上書きしない。
2. `/tmp/losat-wasm-thread-request/` の旧証拠が存在すれば hash を照合して保存する。存在しなければ新しい一時ディレクトリへ再現し、旧測定と同一結果だったとは記載しない。
3. native、serial command、threaded command、threaded library の baseline を別 target directory へビルドする。baseline library は既知の不完全な起動契約を持つ検証対象として識別する。
4. F01 の query/subject 軸、小入力の thread 指定、F02/F03/F13 の停止・失敗を現在の source/artifact で再確認する。F13 の正規 library probe と、監査の command 初期化後の診断 probe を別分類にする。
5. F01 の統計と cutoff の全呼び出し、3 engine の全 pool builder・parallel iterator・reducer・出力経路を列挙する。変更対象の NCBI 所有関数と呼び出しタイミングを対応させる。

**baseline build の例（リポジトリルートから）**

```bash
LOSAT_PLAN_BASELINE_DIR="$(mktemp -d /tmp/losat-wasm-plan-baseline.XXXXXX)"
(
  set -e
  cd LOSAT
  cargo build --locked --release --bin LOSAT \
    --target-dir "$LOSAT_PLAN_BASELINE_DIR/native-command"
  cargo build --locked --release --bin LOSAT \
    --target wasm32-wasip1 --no-default-features \
    --target-dir "$LOSAT_PLAN_BASELINE_DIR/serial-command"
  cargo build --locked --release --bin LOSAT \
    --target wasm32-wasip1-threads --features wasm-threads \
    --target-dir "$LOSAT_PLAN_BASELINE_DIR/threaded-command"
  cargo build --locked --release --lib \
    --target wasm32-wasip1-threads --features wasm-threads \
    --target-dir "$LOSAT_PLAN_BASELINE_DIR/threaded-library"
)
```

crate 内で Cargo を呼び、`LOSAT/.cargo/config.toml` の既存 target flags を適用する。実行環境の `RUSTFLAGS` 等も記録し、baseline と candidate の SIMD/LTO 条件を変えない。library の reactor 化をこの baseline コマンドだけで達成したと扱わない。

**受け入れ条件**

- 各 artifact を hash と実際の exports/imports/memory で識別できる。
- F01 は件数差だけでなく最初の E-value 差を記録する。既知の不一致を baseline の正常期待値に昇格させない。
- 再現用の入力・argv・環境が保存され、旧出力フォルダを上書きしていない。
- この段階は対象を絞った比較診断に限定し、既知の差が残る状態で無関係な広域テストや benchmark を走らせない。

### S01: TBLASTX 統計と出力・エラー境界を修正する

**対象ファイル**

- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`、必要な import を持つ `blast_engine/mod.rs`
- `LOSAT/src/algorithm/tblastx/ncbi_cutoffs.rs`
- `LOSAT/src/algorithm/tblastx/sum_stats_linking/{params,cutoffs,linking}.rs`
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/web_api.rs`
- 上記の既存 unit tests、`LOSAT/tests/unit/tblastx/`、対象比較 fixture

**F01 の作業**

1. subject 読み込み時に統計用の総塩基数・配列数を確定する。worker ごとに再集計しない。NCBI の入力検証・subject range が長さへ及ぼす影響も対応させる。
2. `compute_eff_lengths_subject_mode_tblastx` の single-subject 固定を解消し、DB 総長と配列数を受け取る一つの計算関数へ整理する。総塩基数を合計してから 3 で除算する。配列ごとに切り捨てた値を足し合わせない。
3. query context ごとに `length_adjustment` と `eff_searchsp` を算出し、subject loop の外で保持する。NCBI の context 有効性、整数型、精度、演算順序、既存オプションの上書き規則を保つ。
4. 初期 cutoff、hit-saving cutoff、再評価、sum statistics、最終 E-value への伝播をまとめて修正する。`compute_eff_searchsp_subject_mode_tblastx` と `compute_tblastx_cutoff_score` などの wrapper も調査し、single-subject 固定の別所有者を残さない。
5. `linking.rs` の `calculate_link_hsp_cutoffs_ncbi(..., 0, ...)` を全 subject の `db_length` に対応させる。`LinkingParams` など既存の受け渡しを必要最小限で拡張する。NCBI が subject ごとに更新する値は subject ごとに維持する。
6. `CalculateLinkHSPCutoffs` の計算タイミングと結果の共有を NCBI に合わせる。全 DB 長と個別 subject 長を区別し、引数を一括置換しない。`db_num_seqs = 1 for -subject` 等の誤った説明も修正する。
7. 出力後の E-value 補正、追加ヒットの後付け削除、`BL2SEQ_LEGACY` を使った期待値変更は実装しない。

**F02/F11 の作業**

1. BLASTN の出力先選択を共通化し、並列の早期 return 経路でも `in_memory` があれば既存 `Writer` を使う。reducer と formatter の処理順序は維持する。
2. native を含む同種の出力先選択を調査し、別の経路に同じ欠落があれば同じ修正単位で対応する。direct API のために一時ファイルへ迂回しない。
3. direct API の engine error 変換を共通化し、`format!("{err:#}")` 相当で原因 chain を保持する。失敗時に前回の result bytes を返さないことも既存 ABI のまま確認する。

**受け入れ条件**

- 同じ入力・オプションの native/serial/threaded で出力が一致し、default gencode は通常 NCBI oracle の raw bytes に一致する。
- 1×1、1×3、3×1、3×3 に加え、異なる subject 長、総長の `/3` の端数、短い翻訳 context、異なる E-value 閾値を含む。
- 統計関数の unit test は NCBI 計算を参照し、移植後の実装をそのまま期待値生成に使わない。DB 長が linking cutoff の分岐へ影響する fixture も含める。
- default/non-default gencode は第 5 節の別 oracle 契約で確認する。既存の長配列・閾値回帰を維持する。
- F02 の最終検証は S02 の正規 reactor で行い、n1 と最初の n2 の API output が一致する。F11 は context と下位原因を取得できる。
- 変更対象の parity sweep を完了してから focused tests を実行し、その後に対象プログラムの比較へ広げる。

### S02: command/reactor の成果物・初期化契約を修正する

**対象ファイル**

- `LOSAT/Cargo.toml`、`LOSAT/.cargo/config.toml`
- `LOSAT/tests/run_losat_wasi.js`、`LOSAT/tests/run_losat_wasi_threads.js`
- `LOSAT/tests/wasi_shared_memory.js`（共通処理が必要な場合のみ）
- `LOSAT/tests/prepare_release_candidate_v010.py`、関連 tests・build 呼び出し元
- `README.md`、`RELEASE.md`、`LOSAT/tests/README.md`
- direct API 検証用 host と artifact 検査（新規 test ファイル。製品 API を増やさない）

**F14 の作業**

1. command は `--bin LOSAT`、reactor は `--lib` を明示し、target directory を分離する。crate 名の全面変更より、この build 境界の分離を優先する。
2. `build-web-threaded` の役割を明確にし、library 用 alias と出力を分ける。serial の command/library にも同じ衝突が起きないようにする。
3. 配布名は種別を表す別名にする。artifact metadata へ target、features、Rust/LLVM、link flags、hash、exports、imports、shared memory 属性を記録する。
4. runner と packaging は種別を検査し、同名だったことや `_start` がないことだけで reactor と推定しない。既存 workflow・比較 script の artifact 参照先を更新する。

**F13 の作業**

1. 実 toolchain の CRT と linker 出力を確認し、reactor の `_initialize` から main thread の実行基盤・constructors へ到達する正式な経路を成立させる。具体的な linker 設定はこの調査結果で決め、推測のフラグを計画上の確定解にしない。
2. API export および `wasi_thread_start` の wrapper を逆アセンブルし、API 呼び出しごと・child entry ごとに process 全体の destructor が実行される構造を解消する。通常の thread 終了に必要な処理は維持する。
3. host は main instance の初期化を一度だけ行う。child instance の thread entry と main の `_initialize` を混同せず、child 作成時に共有状態を再初期化しない。
4. まず同じ target/link 設定の最小 Rust thread spawn/join probe で起動基盤を検証し、その後に実 LOSAT の最初の並列呼び出しを確認する。停止が残る場合は命令・lock・呼び出し位置を特定する。
5. command を `--help` で起動してから exports を呼ぶ監査用手法を製品の direct API 起動方式にしない。TLS の手動初期化だけで完了としない。

**artifact の検査項目**

| artifact | 起動入口 | thread 契約 | 用途 |
|---|---|---|---|
| serial command | `_start` | thread spawn 不要 | CLI n1 |
| threaded command | `_start`、child 用 `wasi_thread_start` | shared memory、thread-spawn import | CLI n1/nN |
| threaded reactor | `_initialize`、API exports、child 用 `wasi_thread_start` | shared memory、thread-spawn import | 正規初期化後の反復 API |

serial direct library を配布する場合も、専用の出力・初期化契約を検査する。browser の `wasm32-unknown-unknown` を、この WASI reactor と同じ契約で認証しない。

**受け入れ条件**

- command→library、library→command の両 build 順序で先に生成した artifact の hash が変わらない。
- 誤った artifact を渡すと検索前に明示的に失敗する。正しい artifact は標準の起動手順で動作する。
- reactor で最初の thread spawn/join が成功し、3 engine の最初の並列検索が timeout せず完了する。F02 の output gate もここで実行する。
- この段階では旧しきい値が残るため、確実に並列経路へ入る複数 subject fixture または監査と同じ診断用しきい値 0 を使い、worker 起動を確認する。S03 以降の最終 gate は通常設定で実行する。
- F13 の原因が未特定・起動未成立なら S02 は未完了とする。一時的な unsupported error は安全措置であり、F13 解決の代替ではない。

### S03: スレッド数とプールの寿命を共通化する

**対象ファイル**

- 3 engine の pool builder、thread count 解決、parallel gate
- `LOSAT/src/utils/` の小さな共通 thread 実行 helper と module 登録（既存の所有箇所で足りなければ新規）
- CLI/direct API から上記 helper への入口
- threaded host、要求数・再実行 tests、旧 thread 計画の案内

**設計**

- n1 は caller で直列実行し、追加計算 worker を作らない。
- nN は N 本の専用計算 worker を持つ local pool を一検索につき一つ作る。caller はその完了を待ち、計算 worker 数には含めない。
- `use_current_thread()` と global pool への暗黙依存を除去する。BLASTP の preparation、search、Kappa redo も一つの pool を共有する。
- `ThreadPoolBuilder::build_scoped` を第一候補とし、pool 破棄後の Rust worker join を検索スコープ終了までに完了させる。使用 target で成立しない場合は原因を特定し、所有する JoinHandle と明示 join による同じ寿命契約を実装する。join なしの Drop やリークへ迂回しない。
- 共通 helper の役割は要求検証、pool の構築・実数確認・終了管理・実状態の記録までとする。検索候補の生成や algorithm の判断を移さない。

**要求数の規則**

| 条件 | 処理 |
|---|---|
| 未指定 | 既存既定値 1 を維持 |
| 0・不正な数値 | 既存の引数検証で拒否。新しい auto モードを追加しない |
| serial Wasm / parallel 無効 native で N > 1 | 検索・worker 起動前に unsupported error |
| 対応 target で N > 1 | N をそのまま pool builder へ渡す |
| N が `rayon::max_num_threads()` を超える | worker を起動せず、要求値と上限を示して失敗 |
| pool 構築後の実数が N と異なる | pool を回収し、明示的に失敗 |
| OS/host が途中の worker 起動を拒否 | 起動済み worker を回収し、原因付きで失敗 |

現行 rayon-core 1.13.0 の wasm32 上限は 255 だが、検証には公開 API の値を使い、独立した定数を増やさない。255 本の実起動を通常 CI の合格条件にしない。

BLASTP の `requested.min(cpu_count)` は native/Wasm とも要求検証へ統合する。`num_cpus` や Node の CPU 推定値を暗黙 cap にしない。`LOSAT_WASI_THREAD_CAP` は host が明示した資源上限として扱い、設定がある場合は全 engine で同じく「超過なら拒否」とする。runner による CPU 推定値の自動注入は除去する。未設定は上限 1 を意味しない。

明示 cap の 0・不正値を無視して既定値へ戻さず、設定エラーとして報告する。API と CLI の両方が通る共通入口で検証し、途中の engine 呼び出しで別の値へ解決し直さない。

既存の Wasm MIN_JOBS/MIN_BASES 等は、明示 N を無効化する入口から除去する。撤去する診断変数は docs/tests と併せて整理し、新しい暗黙 auto 経路を作らない。

**lifecycle の受け入れ条件**

- 正常な nN 検索は `pool_threads=N`、新規起動・ready を確認した計算 worker は N 本。旧 N−1 期待値を runner/test/集計で同時に更新する。
- すべての parallel iterator が指定 local pool の配下にある。nested linking/DP/Kappa が別 pool を作らない。
- 同一 instance の `1→2→4→2→1` と、BLASTN→TBLASTX→BLASTP の切り替えが成功し、各出力が独立した n1 の同条件出力に一致する。
- 1 本目と途中の thread-spawn 失敗を test host で注入し、API が原因を返した後に有効な再検索が成功する。
- Rust 側は検索スコープ終了時に join 済み。host 側の Worker exit event は JavaScript のイベントループへ戻った後に確認し、固定の sleep ではなく終了通知と期限を使って残存 worker が 0 になることを検証する。
- guest trap・abort は回復可能な spawn error と区別する。破損した instance の再利用を要求せず、host の明示失敗・回収を検証する。
- worker 異常終了の検出を、Wasm 内で待機している main thread の JavaScript callback だけに依存させない。既存の即時失敗処理と外側の timeout を維持し、起動途中・計算途中の異常を実行テストで覆う。
- 反復で worker 数・未回収 allocation が増え続けないことを確認する。縮まない linear memory の high-water mark だけで leak と判定しない。

### S04: 段階別の並列処理と診断を修正する

**対象ファイル**

- BLASTN/TBLASTX/BLASTP の並列段階と診断
- `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
- `LOSAT/tests/run_losat_wasi_threads.js`
- `LOSAT/tests/wasm_performance.py`、`test_wasm_performance.py`、関連 summary 処理
- scan 診断変数の参照元と説明

**作業**

1. プールが利用可能かという条件と、subject traversal・DP・linking・Kappa の仕事が複数あるかという条件を分ける。pool の有無を subject 数から再推定しない。
2. TBLASTX 1 subject でも複数 frame group の linking を既存の並列経路へ投入できるようにする。subject traversal が直列でも linking を `pool.install` の外へ出さない。
3. 完了順で hitlist/heap を更新せず、subject、query、frame group の NCBI 順序へ戻す既存 reducer を維持する。
4. 診断は構築済み pool と実行した stage の状態から生成する。BLASTN DP だけが並列となる場合も正しく表現する。
5. F10 の変数は逐次 scan 分割を示す名称へ改め、参照元・説明を更新する。旧名を受け付ける必要がある場合は警告付きで同じ一つの診断経路へ接続する。scan 本体を新しく並列化しない。
6. performance 集計から `parallel=false → effective_compute_threads=1` の推定を除去し、診断不足や Rust/host の矛盾を合格扱いしない。旧形式の証拠を新しい意味へ自動変換しない。

**記録する値**

| 値 | 意味・取得元 |
|---|---|
| requested threads | 検証前の利用者指定 |
| pool threads | 構築済み pool の実サイズ。n1 の追加 pool は 0 |
| spawned / ready / exited workers | host の実イベント。spawn 試行と成功を別に数える |
| caller participates | 新設計では nN の計算 pool に参加しない |
| stage / work items / parallel selected | 実行段階、その独立した仕事数、選択経路 |
| measured activity | 必要な profiling 実行で測った worker 別仕事数・実仕事時間。未計測なら null |

**受け入れ条件**

- 1 subject の `blastn` DP、megablast、TBLASTX linking、BLASTP preparation/redo で、段階の違いが記録に現れる。
- 小入力や不均等な仕事量でも、要求 N と pool N の契約を満たす。全 worker の同時稼働は通常 gate にしない。
- 診断は通常無効・stderr のみ。診断の有無で検索出力が変わらない。
- host 起動数を意図的に欠落・改変した test record が集計で失敗し、誤った性能結果を作らない。

### S05: 回帰 gate と CI・配布検証へ接続する

**作業**

1. `test_wasi_runners.js` / `test_wasi_shared_memory.js` の host 検証を維持する。その上に、実 LOSAT の要求数・正規 reactor・反復を検証する専用 integration runner を追加する。
2. `.github/workflows/ci.yml` に serial command、threaded command、threaded reactor の build と実行を追加する。native の `cargo test --all-features` を Wasm 実行の代用にしない。
3. Rust と Node は検証した版を CI 設定と evidence に明示する。監査の Node 18 再現と、配布でサポートする host の合格を区別する。
4. release readiness の build 呼び出し・artifact 選択にも種別検査を接続する。threaded artifact の対応を主張する場合、その正規起動と実行 gate を必須にする。serial だけの既存結果を threaded 認証に流用しない。
5. 各修正の focused verification 後、対象プログラムの既存比較・CLI tests・pure Rust runtime boundary・format/clippy/test を実行する。shared core の変更があるため native も対象とする。
6. 既存 canonical/exception manifest の失敗は記録し、自動再生成しない。新規 fixture の期待値は固定 NCBI oracle と source 対応から作る。
7. README と旧計画への後継案内を最終実装に合わせる。production、test、documentation、generated evidence の差分を分けてレビューする。

**受け入れ条件**

- 第 5 節の必須行がすべて実行され、skip/timeout/未知 artifact を成功に含めない。
- 正常系は raw output と実行契約の両方に合格する。異常系は期待する失敗・原因・回収に合格する。
- リリースに関わる parity・順序・性能の結論を確定する前に、`ncbi_parity_auditor` の独立した read-only 確認を受ける。

### S06: 修正後の費用を測り、必要な最適化を選定する

**必須の測定**

1. baseline と candidate の release build、同一入力・argv・出力先で比較する。既知の F01 不一致や F13 停止を含む baseline を、正常な検索との速度倍率に混ぜない。
2. native/threaded Wasm の n1/n2/n4 と serial Wasm の n1 を比較する。serial Wasm の n2/n4 は拒否テストであり性能試料にしない。cold process、compiled module 再利用、同一 instance の反復は別集計にする。
3. 各条件 1 回以上の warmup と 5 回以上の計測を行い、実行順を交互化する。中央値・範囲・peak RSS・出力 hash・worker 契約を保存する。詳細診断は timed samples と分ける。
4. 小入力、既存のしきい値逆転入力、長い単一 subject、複数 subject、大量 HSP/偏った仕事量を含める。CPU/wall 比だけで計算 worker の利用率を推定しない。
5. module guard/compile、worker 起動、検索段階、scratch allocation、順序復元、一時結果保持を個別に測る。NCBI と異なる Big O や不必要な全量保持を導入していないか確認する。

**測定後に限定して検討する改善**

- 同一 bytecode の guard 解析・compile の再利用。cache key に artifact hash、guard 版、設定を含める。
- 独立した worker 初期化の費用削減。ただし起動失敗の伝播、共有 memory 初期化、thread entry の ready 契約を保つ。
- worker scratch の適切な再利用と一時結果の保持量削減。順序依存の pruning を前倒ししない。
- 常駐 pool は検索単位 pool の正しさを基準として、寿命・サイズ変更・終了管理を持つ別変更として評価する。単に `use_current_thread()` を戻さない。

shared-memory guard の単純除去、新しい seed/候補 pruning、scan 分割、無制限 pool cache は採用しない。最適化候補ごとに source 根拠・出力一致・実測改善が得られなければ採用しない。F01〜F14 の修正完了に、未測定の高速化率や全入力での speedup を条件として付けない。

## 5. 必須検証マトリクス

全組合せを無目的に増やさず、まず次の各行で触る分岐を覆い、その後に既存比較 suite へ広げる。

| 検証軸 | 必須ケース | 判定 |
|---|---|---|
| プログラム | blastn DP、megablast、tblastx、blastp | それぞれの pool と出力を検証 |
| target | native parallel、native no-default-features、serial command、threaded command、threaded reactor | 対応/非対応と artifact 種別を明示 |
| 正常な要求数 | n1/n2/n4。serial は n1 | raw output、一検索の pool サイズ、worker lifecycle |
| 拒否する要求 | serial n2/n4、0、不正値、Rayon 上限+1、明示 host cap 超過 | 検索開始前の error、出力・起動なし |
| F01 の入力 | 1×1、1×3、3×1、3×3、異なる subject 長、端数を持つ総長 | NCBI の有効長・cutoff・E-value・最終 raw bytes |
| スケジューリング | 小さい単一/複数 subject、従来の n2 並列/n4 直列 fixture、1 subject の linking | 要求を縮小せず、実行段階の診断が正しい |
| 出力先 | CLI file、既存の stdout 経路、API memory | sink 選択の一致、失敗時の stale result なし |
| 出力形式 | 対象 engine が対応する outfmt 0/6/7 と既存 custom fields | 既存の形式別 oracle 契約で bytes 比較 |
| 反復 | 同一 instance の 1→2→4→2→1、同じ N の連続実行、3 engine 切り替え | 各回の結果、worker 回収、再登録 error なし |
| 回復可能な失敗 | 最初/途中の spawn 拒否、出力先エラー後の正しい要求 | 原因 chain、終了待ち、次回検索成功 |
| 回復不能な失敗 | guest trap、worker 異常終了、起動 timeout | 成功扱いせず、host が残存 worker を回収 |
| artifact | 両 build 順序、誤種別、欠落 export/import、非 shared memory | 正しい種別だけを受理 |
| 診断/集計 | BLASTN DP のみ並列、欠落イベント、誤った起動数 | 推定で埋めず、矛盾を fail |
| 既存回帰 | AP027131/AP027133 gencode 4、LC738874/LC738875 の E-value 10/100/10000、対象既存 manifest | 承認済み境界を維持 |

### oracle と期待値の扱い

- 通常の default-gencode fixture は、監査と同じ通常 NCBI BLAST+ 2.17.0 の controlled raw output を基準にする。`BL2SEQ_LEGACY` は通常比較から外す。
- 非既定 subject gencode は既存の承認済み fixture に加え、対応する非既定コードの translation/search/reporting が維持されることを検証する。NCBI local `-subject` の遺伝暗号差だけを不合格理由にしない。
- 非既定コードの新規 oracle が必要な場合は、同じ配列を使う NCBI database 検索または対応 NCBI source の検証経路を先に特徴付ける。database と local subject の統計・masking・識別子等の条件を確認し、単に `-db` へ替えた出力を無条件に期待値にしない。NCBI database の作成・検索は tests のみに置く。
- 形式 0/7 のラベル・パス・header が比較に影響する場合、既存の形式別 fixture と固定入力名を使う。比較後の header 除去やソートを新しい合格条件にしない。
- API と CLI に別の既存ラベル契約がある場合は各契約の期待 bytes を固定し、それぞれで n1/nN を比較する。sink 修正のために無関係なラベルを変更しない。
- 既存 native 認証は Gate A/B をそのまま実行する。新規 fixture の通常 NCBI 比較と、既存認証の platform fingerprint 検査を混ぜない。

## 6. 実装時に参照するソース

行番号は計画作成時の入口。実装直前に使用 source snapshot の版・hash と行番号を確認してコメントへ反映する。

| 項目 | NCBI / 依存ソース | 対応する根拠 |
|---|---|---|
| subject 全体の集計 | `c++/src/app/blast/blast_app_util.cpp:204–210`、`c++/src/algo/blast/api/seqsrc_multiseq.cpp:175–180` | 通常 CLI の dbscan mode と全 subject 長 |
| 有効長 | `c++/src/algo/blast/core/blast_setup.c:716–740,770–847,852–884` | DB 長・件数、context、翻訳長換算、search space |
| subject ごとの更新条件 | `c++/src/algo/blast/core/blast_engine.c:1407–1455` | `db_length == 0` の場合だけ再計算、link cutoff の入力 |
| word/hit cutoff | `c++/src/algo/blast/core/blast_parameters.c:340–374,942–976` | cutoff 入力と cap、更新条件 |
| linking cutoff | `c++/src/algo/blast/core/blast_parameters.c:998–1082` | DB 総長と個別 subject 長を別々に利用 |
| linking 順序 | `c++/src/algo/blast/core/link_hsps.c:553–558,959–982` | frame group と chain の更新 |
| NCBI thread lifecycle | `c++/src/algo/blast/api/prelim_stage.cpp:145–188` | N threads 作成・起動・Join |
| CLI thread 指定 | `c++/src/algo/blast/blastinput/blast_args.cpp:3152–3187,3205–3236` | 制約、CPU 削減、local subject の警告 |
| output stream | `c++/src/algo/blast/format/blast_format.cpp:68–96,770–832` | caller が渡す ostream と formatter 選択 |
| current-thread の残留 | `rayon-core 1.13.0/src/lib.rs:534–546`、`src/registry.rs:295–311` | caller 登録の残留と再登録失敗 |
| Rayon 上限 | `rayon-core 1.13.0/src/registry.rs:245–246`、`src/sleep/counters.rs:60–77` | builder 内の縮小、wasm32 上限 |
| scoped pool | `rayon-core 1.13.0/src/lib.rs:322–343`、`src/thread_pool/mod.rs` の Drop | scoped spawn/join と Drop の役割 |
| WASI 初期化 | 使用 toolchain の `crt1-command.o`、`crt1-reactor.o`、生成 Wasm の起動関数・wrapper | F13 の実際の link 結果。未確認の停止原因を断定しない |

NCBI source root は `/mnt/c/Users/genom/GitHub/ncbi-blast/`。Rayon は現在の Cargo.lock と実際の registry source を使う。runtime の仕様確認には [Rayon ThreadPoolBuilder](https://docs.rs/rayon/latest/rayon/struct.ThreadPoolBuilder.html)、[Rust WASI threads target](https://doc.rust-lang.org/rustc/platform-support/wasm32-wasip1-threads.html)、[Node WASI](https://nodejs.org/api/wasi.html) も参照し、オンライン最新版を固定 toolchain の実装と同一視しない。

## 7. 証拠・引き渡し・完了確認

各段階で [スキルの evidence template](../.agents/skills/verify-ncbi-parity-and-speed/references/evidence.md) に沿った記録を残す。最低限、入力 hash、順序付き argv、環境、source/artifact/runner hash、NCBI 版、終了コード、raw output と hash、stderr、差分、worker events、実行時間・RSS を保存する。

baseline、candidate、oracle、diagnostics、timed samples は別に識別する。生の出力と異常系の証拠も保存し、要約 JSON だけを証拠にしない。`/tmp` の証拠は必要なものを保管先へアーカイブし、アーカイブ hash と再現コマンドを記録する。保管先が失われた場合、再現可能な記述と実際に保管された証拠を区別する。

実装は、統計修正、I/O、artifact/CRT、pool/要求数、段階診断、CI の変更単位に分ける。必要な tests と説明は各変更に含める。共有作業ツリーの既存変更を混ぜて commit しない。publish・push・tag・deploy は本計画の成果物に含めない。

- [x] S00: 現行 artifact と再現条件を固定した。
- [x] S01: 全 subject 統計、link cutoff、出力先、error chain を修正した。
- [x] S02: command/reactor の衝突と正規初期化を解決した。
- [x] S03: 要求数・一検索一 pool・終了待ち・失敗回復を成立させた。
- [x] S04: 段階別の利用可否と実状態の診断を一致させた。
- [ ] S05: CI gate 実装・ローカル matrix・独立監査は完了。固定 regression の残り 46 件はユーザー指示で延期。remote CI は未実行。
- [x] S06: 正しい出力を対象とする性能比較と、残る費用・制限を記録した。
- [x] F01〜F14 の各行に実行証拠を結び付け、未解決・未対応を完了扱いしていない。
- [x] 凍結 canonical/authority、承認済み gencode 例外、既存ユーザー変更を維持した。

F13 の停止、raw output の差、反復失敗、worker 数の矛盾、未知の native fingerprint が残る場合は、その項目と依存する認証を未完了とする。修正候補が速いことや、serial の結果だけが合っていることを理由に完了へ変更しない。
