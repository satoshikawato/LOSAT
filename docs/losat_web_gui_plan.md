# LOSAT Web GUI 総合実装計画

状態：**S01（W0）から S07++（E2f）までは完了条件を満たした（[W0](evidence/losat_web_w0/README.md)、[E1a](evidence/losat_web_e1a/README.md)、[E1b](evidence/losat_web_e1b/README.md)、[E1c](evidence/losat_web_e1c/README.md)、[E1d](evidence/losat_web_e1d/README.md)、[E2a-1・E2a-2](evidence/losat_web_e2a/README.md)、[E2c](evidence/losat_web_e2c/README.md)、[E2f](evidence/losat_web_e2f/README.md) のゲート記録。アプリ側の S10（W2）も完了した（[W2](evidence/losat_web_w2/README.md)）。E1a〜E1c の V-PERF の判断と、E1c の CLI の 2 つの振る舞いの差は、2026-10-02 に保守者が記録どおり承認した（DW-13））。`PD-LOSAT-WEB-APP-BOUNDARY` は 2026-09-29 に承認された。S07+++ と S07+++b（E2g）も完了条件を満たした（[E2g](evidence/losat_web_e2g/README.md)。棚卸し `INVENTORY.tsv` 1015 行の一括の transpile、NCBI の不具合の扱いの保守者の判断 DW-15、独立監査は第 4 回で supported）。S08・S08b（E2b：TBLASTX outfmt 0/7）も完了条件を満たした（[E2b](evidence/losat_web_e2b/README.md)。独立監査は第 3 回で supported、保守者の判断 DW-17。V-PERF の 1 case は保守者の確認待ち）。SD（E2i：BLASTN の dc-megablast と blastn-short、DW-18）も完了条件を満たした（[E2i](evidence/losat_web_e2i/README.md)。棚卸し `INVENTORY.tsv` 285 行の一括の移植、独立監査は第 1 回で 4 観点とも supported）。S08+・S08+a・S08+b（E2e：BLASTP・TBLASTN・TBLASTX の既定以外の option、TD-13）も完了条件を満たした（[E2e](evidence/losat_web_e2e/README.md)。棚卸し `INVENTORY.tsv` 1063 行、sweep 3264 組で差 0、独立監査は第 2 回で 4 観点とも supported、保守者の判断 DW-19）。S11（E2d：`-query_loc` / `-subject_loc`、DW-9〜DW-11）も完了条件を満たした（[E2d](evidence/losat_web_e2d/README.md)。棚卸し `INVENTORY.tsv` 248 行、範囲の fixture 53 と回帰 fixture 35 で NCBI とバイト一致、範囲の sweep 1860 件で差 0、独立監査は 4 観点と再監査 2 回、保守者の判断 DW-20）。エンジン側の残りは SF（E2h、FASTA の読み方。アプリ側の S13 と並行できる）と SX（条件待ち）。アプリ側の S09（W1：ブラウザでの実行基盤、V-BR）と S12（W3：検索画面、研究作業と境界条件の E2E）も完了条件を満たした（[W1](evidence/losat_web_w1/README.md)、[W3](evidence/losat_web_w3/README.md)。program の表示名は保守者の判断 DW-21）。W1 が諮った 2 件（TBLASTN の R2 の進め方、協調取消の閾値）は推奨の案で決まった（DW-22）。次はアプリ側の S13（W4）。** 作成 2026-09-28、改訂 2026-10-06。

| 項目 | 内容 |
|---|---|
| 作業ブランチ | **`feature/losat-web-gui`**。エンジン側のセッションはこのブランチを使い、アプリ側のセッションは `feature/losat-web-gui-app` で行ってセッションの終わりにこのブランチへ merge する（DW-7、2026-09-30 の改訂）。基点は `origin/main` の `0627c88f5`（2026-09-28 時点の最新）。リポジトリに `dev` ブランチは無いため、`main` を基点にした |
| 作業ディレクトリ | git worktree **`/mnt/c/Users/genom/GitHub/LOSAT-web-gui`**（アプリ側は `/mnt/c/Users/genom/GitHub/LOSAT-web-gui-app`）。ブランチと作業場所の規則（新しい clone や worktree を作らない、など）は、[セッション README](losat_web_gui_sessions/README.md) の規則 1〜2 にある |
| 要求の原典 | [`docs/web/losat_web_design_v0.1.md`](web/losat_web_design_v0.1.md)（以下「設計書」）と、それを要求 ID に整理した [`docs/web/requirements_trace.tsv`](web/requirements_trace.tsv) |
| 境界の規則 | [`PD-LOSAT-WEB-APP-BOUNDARY`](product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md)、[`web/AGENTS.md`](../web/AGENTS.md) |
| セッション指示書 | [`docs/losat_web_gui_sessions/README.md`](losat_web_gui_sessions/README.md) |

---

## 0. この計画の読み方

### 0.1 背景

LOSAT は、NCBI BLAST+ を純 Rust で再実装したものである。宣言した fixture では、NCBI と出力バイトが一致することを目標にしている。Rust crate は `LOSAT/` にあり、作業規約はルートの `AGENTS.md` が定める。要点は次の 3 つ。

- 動作の権威は NCBI の C/C++ ソースだけである。
- エンジンのコードを変更するときは、該当する NCBI ソースの位置と断片をコメントに書く。
- NCBI に無い機能はエンジンに入れない。

この計画は、LOSAT をブラウザの中で動かす、インストール不要の Web アプリ（以下「LOSAT Web」）を作るためのものである。何を作るかは設計書が定める。設計書に出てくる Q001–Q060・R01–R30 は、設計書を作るときに使った質問票の回答番号で、回答そのものはリポジトリに無い。この計画ではそれらを使わず、設計書の節番号（「設計書 §2.1」）と、要求トレース表の ID（`REQ-01` など）で参照する。設計書の中で定義されている T01（設計書 §0）と D01〜D10（設計書 §17）は、そのまま使う。

### 0.2 用語

| 用語 | 意味 |
|---|---|
| outfmt 0 / 6 / 7 | BLAST+ の出力形式。0 はペアワイズのアラインメント表示、6 は表形式、7 は注釈行付きの表形式 |
| LOSATN / LOSATP / LOSATX / TLOSAN / TLOSATX | LOSAT の中の各 program の実装名。順に BLASTN、BLASTP、BLASTX、TBLASTN、TBLASTX。CLI はどれも `LOSAT <program>` |
| LOSATX 計画・TLOSAN 計画 | BLASTX と TBLASTN の v0.2.0 の実装と認証の計画（`docs/losatx_blastx_v0.2.0_plan.md`、`docs/tlosan_tblastn_v0.2.0_plan.md`） |
| 認証済みプロファイル | NCBI との一致を記録と独立監査で確かめた、program・task・オプション・実行経路の組。記録は `docs/release/`、`docs/evidence/`、`LOSAT/tests/*manifest*.tsv` にある |
| 凍結バイト | 認証済みプロファイルで NCBI と一致すると確かめ、SHA-256 で固定した LOSAT の出力。例：`LOSAT/tests/platform_native_v010_canonical.tsv`（`PD-NCBI-PLATFORM-VARIANCE` の Gate A）、`docs/evidence/tlosan_stage_g/` |
| command / reactor | WASI の Wasm モジュールの種類。command は `_start` で一度だけ実行する。reactor は `_initialize` の後、export された関数を何度でも呼べる |
| COOP / COEP / CORP | 共有メモリ（`SharedArrayBuffer`）を使うのに必要な HTTP ヘッダー。揃うと `crossOriginIsolated` が真になる |
| OPFS | ブラウザが origin ごとに持つ、非公開のファイル領域 |
| PD | 製品決定の文書（`docs/product_decisions/`） |
| gbdraw | 同じ保守者のゲノム描画 Web アプリ（<https://github.com/satoshikawato/gbdraw>、この計画の作成時に参照したのは commit `538e9ec5`）。LOSAT の Web ABI v1 を使っている。ブラウザでの LOSAT の実行は `gbdraw/web/js/services/losat.js` と `gbdraw/web/js/workers/` にあり、threaded を使う入力の大きさの既定値は `losat.js` の `DEFAULT_THREADED_MIN_FASTA_CHARS = 500000` |
| Web ABI v1 / v2 | v1 は `LOSAT/src/web_api.rs` の `losat_web_*` で、gbdraw が使う。v2 はこの計画で新設する `losat_web2_*`（[`docs/web/abi_v2.md`](web/abi_v2.md)） |
| 計画レビュー・独立監査・画面レビュー | それぞれ `plan_critic`、`ncbi_parity_auditor`、`visual_regression_reviewer` という名前のレビュー役。定義は保守者の Codex 設定（`~/.codex/agents/`）にあり、その要点はセッション README の「レビュー」の節に書き写してある |

### 0.3 解釈の規則

1. 設計書で「確定している要求」とされたものは、§0.4 の決定で変更したものを除いて削らない。
2. 設計書で「提案」とされたもの（内部API、Worker分割、保存形式など）は、KISS と YAGNI に従って最小構成で始める。後回しにしたものには、導入する条件（トリガー）を付ける（§2.3）。
3. 検索の意味と出力バイトの権威は、NCBI BLAST だけである。Web アプリはこれを変えない。
4. 各段階の完了条件は §7 の表が正本である。セッションの指示書は、それに細部を足すことはできるが、条件を緩めることはできない。

### 0.4 保守者の決定

| ID | 日付 | 決定 | この計画への影響 |
|---|---|---|---|
| DW-1 | 2026-09-28 | 同一リポジトリで、規約の適用範囲を分ける | `LOSAT/` は現行の `AGENTS.md` のまま。`web/` には `web/AGENTS.md` を置き、`PD-LOSAT-WEB-APP-BOUNDARY` で、検索の意味と出力バイトを変えないアプリ機能を認める |
| DW-2 | 2026-09-28 | Vue 3 + Vite + TypeScript | gbdraw（Vue 3）と知識を共有する。アプリの状態は Vue に依存しない TS で書き、Vue は表示だけを担う |
| DW-3 | 2026-09-28 | Cloudflare で配信する | `_headers` で COOP / COEP / CORP を設定する（gbdraw と同じ）。GA4 用の計測文書が必要になれば、別サブドメイン（別 origin）に置く |
| DW-4 | 2026-09-28 | パラメーターの初期値は BLAST+ CLI の既定値 | 画面の項目名と並びは NCBI Web BLAST に寄せ、値は CLI の既定値にする |
| DW-5 | 2026-09-28 | query ごとの途中経過は作らない | 設計書 §2.1「途中結果」行の「完了queryの結果は仮表示」、§3.1 の「完了query数」、§9.3 の仮表示を取り下げる（`REQ-09`）。結果は実行が正常に完了した後にまとめて表示する。取消・失敗した実行の結果を破棄する規則は維持する（`REQ-10`） |
| DW-6 | 2026-09-28 | outfmt 0 の不足分はこのブランチで整備する | BLASTN の outfmt 0、TBLASTX の outfmt 0 と 7 を NCBI から移植する（§7 の S06〜S08） |
| DW-7 | 2026-09-28、改訂 2026-09-30 | エンジン側とアプリ側の 2 本まで | セッションは、エンジン側（`LOSAT/`・`web/adapter/` を変えるもの）とアプリ側（`web/app/` だけを変えるもの）の 2 本までを並行して、それぞれの中では表の順に 1 つずつ実行する。アプリ側は別の worktree とブランチで行い、セッションの終わりにエンジン側のブランチへ merge する。V-PERF を測る間はアプリ側の試験を止める（2026-09-30、保守者の承認。当初は単一ブランチ・単一 worktree で 1 つずつ） |
| DW-8 | 2026-09-29 | 前処理の再利用は、実測で必要な program だけ | 設計書 §2.1「反復検索」の再利用のうち、初期に入れるのは入力の参照・解析・warm なエンジンまで（`REQ-07`）。エンコードや翻訳などの前処理のキャッシュは、ブラウザでの実測で前処理が warm 実行時間の 20% 以上を占めた program にだけ入れる（§4.6 の R2） |
| DW-9 | 2026-09-29 | 領域の指定は、その役割のレコードが 1 つのときだけ | NCBI の `-query_loc` / `-subject_loc` は読み込むすべてのレコードに同じ範囲を適用する（§4.9）。レコードごとに別の範囲を指定するには検索を分ける必要があり、統計が変わるので行わない（`REQ-06`） |
| DW-10 | 2026-09-29 | TBLASTN は今、BLASTX は認証の後 | 核の入口へのまとめ直し（§4.2）は、BLASTP・TBLASTN・BLASTN・TBLASTX を先に行う。TBLASTN は TLOSAN 計画の認証のゲートをこのブランチで再実行する。BLASTX は、LOSATX 計画の v0.2.0 の認証が `main` に入った後に、`main` を取り込んでから扱う（§7 の SX）。進行中の LOSATX の作業とぶつからないようにするためである |
| DW-11 | 2026-09-29 | BLASTX も領域の指定の対象にする | LOSATX 計画は v0.2.0 の範囲外として `-query_loc` / `-subject_loc` を拒否している（`LOSAT/src/cli.rs:238-246`）。v0.2.0 の認証の後に、このブランチで BLASTX に移植する（SX）。v0.2.0 の後の範囲の拡大として扱う |
| DW-13 | 2026-10-02 | CLI の検索以外の差の扱いと FASTA の読み方の移植の時期 | S07+++ の棚卸しが保守者の判断に残した項目を決めた（`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`、版 1.0、AGENTS.md）。承認済みの例外：引数の構文の誤りと `-help`（LOSAT の文言と終了コード 2）、`-subject` での `-num_threads`（LOSAT は並列のまま、NCBI の 2 つのスレッドの警告を出さない）、outfmt 6/7 の書き込みの失敗（NCBI は abort、LOSAT は誤りを報告）、メモリ不足（LOSAT は abort）。移植：outfmt 0 の書き込みの失敗（`BLAST failed to write output`、終了コード 6）。明示的な拒否：UTF-8 でないファイル名、出力を変える `.ncbirc` の設定。FASTA の読み方（TD-12）は S17 の前の専用のセッション（SF）で移植する。E1a〜E1c の V-PERF と E1c の CLI の 2 つの差は記録どおり承認 |
| DW-14 | 2026-10-02 | outfmt 0 で閉じたパイプへの書き込み | NCBI は SIGPIPE の既定の動作で終わる（終了コード 141、文言なし）。SIGPIPE を既定に戻すには native の `signal` が要り、pure-Rust の境界の検査が拒否するので、LOSAT は outfmt 0 のほかの書き込みの失敗と同じく「BLAST failed to write output」、終了コード 6 を出す。`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 版 1.1 の例外 5、AGENTS.md（S07+++b、E2g T14） |
| DW-15 | 2026-10-02 | NCBI の不具合に当たる挙動の扱い | NCBI BLAST+ が落ちる・debug でだけ防ぐ前提が崩れる・buffer の外を読む入力は、近い入力の NCBI の出力と一致を確かめられる妥当な結果を LOSAT が出せるとき承認済みの例外にする（`PD-LOSAT-NCBI-DEFECTS` 版 1.0、AGENTS.md）。BLASTN：query の塊を NCBI がもう一度分ける CHUNK_SIZE / OVERLAP_CHUNK_SIZE（NCBI は CCoreException、LOSAT は塊を 1 回ずつ検索）、句読点だけの定義行の outfmt 0 の題（NCBI は文字列の外を読み落ちる、LOSAT は文字列の終わりで止める）。NCBI が決まった結果を出すもの（CHUNK_SIZE=1000 の空の batch、-max_target_seqs 2^30〜2^31−51 の 10 への回り込み）は再現する。妥当な結果を確かめられない NCBI の失敗（K-A の表に無い得点系と先頭の無効な query、-max_target_seqs が 2^31−51 より大きい、reward 32767 / penalty −32768、subject の合計 2^31 文字以上、megablast の gap cost 32767 超、batch を分ける負の CHUNK_SIZE の組）は明示的な拒否のまま。NCBI が受け付ける -evalue の +inf / -nan / +nan(1) / 1e999 は移植する（S07+++b、E2g） |
| DW-16 | 2026-10-02 | DW-15 の後に分かった 3 つの細部 | 監査で分かった細部を、保守者が推奨の案で決めた：負の CHUNK_SIZE の組で batch を分けるがもう一度は分けない場合（NCBI は隙間のある塊で検索して hit を落とす、決まった結果）と、16 bit で 0 以下に回る reward 32768 以上（NCBI は全 query を無効にする）は、実用の無い設定なので明示的な拒否のまま（`PD-LOSAT-NCBI-DEFECTS` 版 1.1）。outfmt 6/7 で LOSAT の書き込みが読み手が閉じる前に済んだときの終了コード 0（NCBI は後の flush で SIGPIPE）は、承認済みの例外 3・5 の時間に依る部分（`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 版 1.2、AGENTS.md）（S07+++b） |
| DW-17 | 2026-10-03 | S08 が諮った 2 つの NCBI との差 | 保守者が推奨の案で決めた。(1) outfmt 0 の句読点だけの subject の題（NCBI の `x_CleanAndCompress` が文字列の外を読み、tblastx・tblastn も SIGSEGV で落ちる）：BLASTN の承認済みの例外 2 を TBLASTX と TBLASTN に広げる（`PD-LOSAT-NCBI-DEFECTS` 版 1.2、AGENTS.md）。代わりの題の移植と、BLASTN と同じ `title_sweep.py` の方式での確かめは S08+ の最初の作業で、それまでは明示的な拒否のまま。(2) 起動の時に閉じた標準出力（`>&-`、全 program）：NCBI は最初の書き込みで失敗する（outfmt 0 は終了コード 6、6/7 は abort）が、Rust の runtime が main の前に開く `/dev/null` を LOSAT は呼び出し側の `/dev/null` と区別できないので、報告を捨てて検索の終わりの終了コードで終わる。承認済みの例外 6（`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 版 1.3、AGENTS.md、証拠 `docs/evidence/losat_web_e2b/closed_stdout/`）（S08b） |
| DW-18 | 2026-10-03 | BLASTN の dc-megablast | 保守者の依頼で、BLASTN の `-task dc-megablast` を NCBI と同じにする。v0.1.0 の CLI の前の LOSAT は `dc-megablast` を受け付けたが、blastn の既定値で megablast の lookup を使うだけで、discontiguous の template（`-template_type`・`-template_length`、`blast_nalookup.c` の表、`blast_nascan.c` の走査）は無かった。v0.1.0 の CLI が task を 2 つに絞り、S07+ が明示的に拒否した。エンジン側のセッション SD（段階 E2i）として、S08b の後、S08+ の前に、DW-12 の棚卸しと一括の移植で行う。続けて保守者の依頼で `blastn-short`（NCBI では blastn に reward 1・penalty −3・e-value 1000・word size 7・filter なしを重ねた task、`blast_options_handle.cpp:343-360`）も SD に入れた。`rmblastn` は明示的な拒否のまま（S08b）。SD で完了した（2026-10-04、[E2i のゲート記録](evidence/losat_web_e2i/README.md)。`-template_type`・`-template_length` はどの task にも NCBI と同じに効く。task が内部で決める `-window_size` などは E2c からの明示的な拒否のままで、S08+ の指示書に引き継いだ） |
| DW-19 | 2026-10-05 | S08+・S08+a・S08+b（E2e）が諮った NCBI との差 | 保守者が全て推奨の案で決めた（S08+b）。(1) D11（query の長さ＋window が 2^31 − 1 を超えると NCBI の `Int4` が回り込みヒット無し、決まった結果）は、D8（NCBI が終わらない (2^30, 2^31 − 1]）と合わせて和が 2^30 を超える値の明示的な拒否のまま（BLASTP・TBLASTX）。(2) D12：無限大と DBL_MAX 以上の `-evalue`（NCBI の blastp・tblastn は入力によって SIGSEGV）は明示的な拒否のまま。(3) D13：NCBI が落ちる query の分割の設定は BLASTP・TBLASTN で明示的な拒否のまま（BLASTN の承認済みの例外 1 を広げない）。(4) D14：option の値の位置の NCBI C++ Toolkit の語も明示的な拒否（`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 版 1.4）。(5) D15：BLASTP の one-hit の gapped の始点の窓が配列の外を読む検索（NCBI は配列の外を読み、落ちることがある）は、LOSAT が配列の外を番兵として読む結果（NCBI の valgrind の下の出力とバイト一致）を承認済みの例外 3 にする（`PD-LOSAT-NCBI-DEFECTS` 版 1.3、AGENTS.md）。(6) web ABI v1 の BLASTP の 4 つの要求の誤りの文言がエンジンの今の文言に変わったこと（状態・出力の無さ・検査の順は E1d と同じ）は、TD-1 の S07+ のスレッドの上限の文言と同じ扱いで受け入れる（[E2e のゲート記録](evidence/losat_web_e2e/README.md)） |
| DW-20 | 2026-10-06 | S11（E2d）が諮ったこと | 保守者が全て推奨の案で決めた。(1) R1：`NStr::StringToInt` が読めない範囲の部分（NCBI は build のパスを含む文言で終了コード 255）は明示的な拒否（TD-15 と同じ扱い）。(2) R2：文字の無い区間（範囲の始まりがレコードの長さ + 1）は明示的な拒否（BLASTP の subject は NCBI と同じ）。(3) A-1：BLASTN・BLASTP・TBLASTN・TBLASTX の option の値の UTF-8 でないバイトは、引数の解析の後の明示的な拒否（`-help` と parser の誤りが先）。(4) O-1：ABI v2 の `validate` はレコードを持たないので、引数とレコードの両方に誤りがあると、CLI と違い引数の誤りを先に返す（`docs/web/abi_v2.md`）。(5) 最後のゲートで省いた工程（TBLASTX の option の sweep、capture、Gate A の後半 10 組、V-ABI full の TBLASTX の残り）を受け入れる。(6) V-PERF は非退行として受け入れる（`tblastx` threaded-WASI の 2 つの時間は有意でない）。`main` への PR は今は作らない（[E2d のゲート記録](evidence/losat_web_e2d/README.md)） |
| DW-21 | 2026-10-06 | program の表示名 | BLASTN 系（BLASTN・BLASTP・BLASTX・TBLASTN・TBLASTX）。NCBI BLAST+ と LOSAT の CLI の program 名に合わせる。§10 の項目を閉じる（S12 の最初に保守者が決めた。[W3 のゲート記録](evidence/losat_web_w3/README.md)） |
| DW-22 | 2026-10-06 | W1 が諮った 2 つ（S12 の後に保守者が推奨の案で決めた） | (1) TBLASTN の R2（S09+）：最初に翻訳を表引きにし（`GeneticCode::get`、出力は変えない）、測り直して subject だけの前処理が 20% 未満ならキャッシュは入れない（DW-8 の条件どおり）。20% 以上のままならキャッシュを移植する（[S09+ の指示書](losat_web_gui_sessions/session_s09p_r2_tblastn_subject_cache.md)の作業の順）。(2) 協調取消（§2.3）を入れる閾値は 1 秒：公開する対応規模の最大の入力で、取消の後の再準備（instance の作り直しと Subject の登録し直し）が 1 秒を超えたとき。S09 の実測は最大 0.61 s。S17 で対応規模を公開するときに、その最大の入力で測り直す（[W1 のゲート記録](evidence/losat_web_w1/README.md)の「保守者の判断待ち」） |
| DW-23 | 2026-10-06 | SF（E2h、FASTA の読み方）の範囲と細部 | 保守者が全て推奨の案で決めた（[SF の指示書](losat_web_gui_sessions/session_sf_e2h_blastn_fasta_reader.md)の「保守者の判断」）。(1) 範囲：BLASTN・TBLASTX・TBLASTN・BLASTP の全入力を一度に NCBI の読み込みの経路で読む（DW-13 の BLASTN の TD-12 に、E2b・E2e の同じ種類の拒否を加える。BLASTX は SX）。(2) NCBI が Seq-id として data loader で取り寄せる最初の行は明示的な拒否（LOSAT はネットワークも BLAST DB も使わない）。配列だけの行は NCBI と同じく配列として読む。(3) `>?` の gap の行は、CLI と `run_local` は移植し、ABI v2 の `register` と `scan` だけが明示的に拒否する（gap の残基に元の byte が無く、索引から切り出せない）。(4) ABI v1 の入力の扱いは変えない（TD-1 の凍結）。(5) UTF-8 でない定義行は移植する（題を byte のまま持ち、報告も byte のまま。ABI の JSON だけ置き換え文字） |
| DW-12 | 2026-09-30 | NCBI の経路の棚卸しと一括の transpile | エンジンの NCBI との一致は、独立監査の指摘を 1 つずつ直すのでなく、アプリが出すオプションの範囲で NCBI の実行経路に現れる関数を program ごとに棚卸しし（忠実な移植・差のある移植・未移植・明示的な拒否）、未移植と差のある移植を NCBI の関数ごとに簡略化せずに transpile してから、棚卸しの表を基準に監査する。NCBI の C++ の層（object manager、ASN.1、`CFastaReader` など）は経路にある関数だけを移植する。LOSAT が速度のために NCBI と違う実装にしている箇所（詰めた配列の走査、並列化など）は、出力が同じなら、新しく移植する部分にも同じ方式を使ってよい（S07+ の 16 回の監査の後の、保守者の指示） |
| DW-12 | 2026-09-29 | `PD-LOSAT-WEB-APP-BOUNDARY` を承認する | 状態を Accepted（版 1.0）にした。S02 の入口の条件を満たす |

### 0.5 この計画で下した技術判断

保守者の決定ではなく、この計画の判断である。根拠が変われば、この表を改訂して変更できる。

| ID | 判断 | 理由 |
|---|---|---|
| TD-1 | Web ABI v2 は `losat_web2_*` という別の名前にし、v1 と同じ Wasm モジュールに共存させる。v1 は、fail-fast の欠陥の修正を除いて凍結する（凍結するのは v1 の引数・形式・誤りの扱いで、v1 は `run_local` で CLI と同じエンジンを使うので、エンジンを NCBI に合わせる修正は v1 の検索結果にも及ぶ。S07+ の第 9 回の監査で確認。全 program が共有するスレッドの上限の誤りの文言は、S07+ で末尾に LOSAT が対応しないことを示す句を足したので、v1 でも同じ句が付く。誤りの扱い（検索の前に失敗し、スレッドを起動しない）は変わらない）。v1 の廃止は、gbdraw が v2 へ移る時点で別に決める | v1 を feature で外すと、既存の serial ビルド（`--no-default-features`）から v1 の export が消え、gbdraw が壊れる（`LOSAT/.cargo/config.toml` の `build-web-serial`、`LOSAT/tests/build_wasi_artifacts.py` の serial-reactor） |
| TD-2 | 各 program の検索と出力は、1 本の核の入口 `run_local` だけを通す。CLI の `run`、v1 の `run_web_pair*`、v2 は、それぞれの入力の形を変換して `run_local` を呼ぶ薄い層にする | 入口が 3 本に分かれると、互いにずれていく（`CONTRIBUTING.md` は並行した実装経路を拒否の理由に挙げている） |
| TD-3 | HSP と outfmt 6 の行・outfmt 0 の節との対応は、formatter が「どの HSP を書き始め、書き終えたか」を観測者に知らせることで取る。出力の順番や、テキストの行頭で推測しない | formatter ごとに並べ方と表示する範囲が違う。例：BLASTX の outfmt 0 は、既定では先頭 250 subject のアラインメントだけを示す（`LOSAT/src/algorithm/blastx/args.rs:918-922`） |
| TD-4 | 1 回の検索から複数の形式を出す。形式ごとに解決した検索オプションが食い違う場合は、`run_local` が明示的なエラーで止める（形式ごとに別に検索する仕組みは作らない） | NCBI では hitlist の大きさが出力形式によって決まる場合がある（`blast_args.cpp:2894-2978`）。ただしそれは `-num_descriptions` / `-num_alignments` を使う場合で、LOSAT はどの program もこれを受け付けない（TBLASTN と BLASTX は明示的に拒否し（`cli.rs:192-193`、`:245-246`）、ほかの program は定義していない）ので、今は食い違いが起きない。起きたときの扱いは §2.3 のトリガーにした |
| TD-5 | 検証の期待値は、表の升目（program × 形式 × プロファイル × 実行経路）ごとに出どころを決める。認証済みの凍結バイトがあればそれを使い、無ければ同じ commit のネイティブ CLI を単一形式で実行した出力を使う。後者は「native-equivalent, not NCBI-certified」と表示する | 既存の manifest は 1 件につき 1 形式しか固定していない（例：`blastp_v010_parity_manifest.tsv` と `tblastx_v010_parity_manifest.tsv` はすべて outfmt 6） |
| TD-6 | アダプタは別の crate（`web/adapter/`）にし、reactor の起動は `LOSAT/build.rs` のものを使う（S05 で確かめた方法：その build script の `rustc-cdylib-link-arg` は依存先の cdylib にも渡るので、アダプタは build script を持たない。パスで共有すると `_initialize` が二重に定義される）。`[profile.release]` と Wasm の rustflags は写しを持ち、ビルドの同一性の検査（`web/adapter/tools/check_build_identity.py`：依存の版、profile、rustflags、アダプタに build script が無いこと）で、ずれを機械的に検出する | 規約の範囲を分けたまま（DW-1）、認証済みの command-WASI と同じ条件でビルドする。別の Cargo root は、LOSAT の profile、`LOSAT/.cargo/config.toml` の `+simd128`、`Cargo.lock` を引き継がない |
| TD-7 | スレッド版の共有メモリの最大値は、認証済みの threaded ビルドと同じ（16384 ページ、1 GiB）にする。host が 2 GiB や 4 GiB を試すことはしない | 取り込むメモリの最大値は、モジュールが宣言した最大値を超えられない。最大値を上げるには link の引数を変える必要があり、TD-6 の同一性とぶつかる。上げる条件は §2.3 に置いた |
| TD-8 | 索引用の FASTA の走査（`scan`）は、アプリの抽出のためだけの例外として置く。検索に渡す入力は、各 program の解析器が読む。`scan` は解析器の種類（`bio::io::fasta` 型か、BLASTX の NCBI 型か）を受け取り、その解析器とだけ性質試験で照合する。食い違いの最終的な判定は、`register` の時点のエンジンの解析結果との照合で行う | BLASTX の解析器は NCBI の CFastaReader に従い、`;` などで始まる行を読み飛ばす（`LOSAT/src/algorithm/blastx/input.rs`）。`bio` はそうしない。1 つの走査で両方と一致させることはできない |
| TD-9 | 既存の不具合のうち、BLASTP と TBLASTN の outfmt 0 の座標の桁数（NCBI は 0 始まりの最大値、LOSAT は 1 始まりの最大値から求める）は、S07 で NCBI の規則の 1 つの関数にまとめて直す。BLASTX の関数は変えない（DW-10） | BLASTN の outfmt 0 も同じ規則を使うので、共有の部品として直すのが最も小さい。境界（最大の座標が 10 の累乗）でだけ出力が変わり、変わる凍結出力は NCBI と一致することを示せる（S06 の `docs/evidence/losat_web_e2a/AUTHORITY.md` §D.1） |
| TD-10 | BLASTN の得点のオプション（既定以外の reward / penalty / gap）は、S07+ で NCBI と同じにするか、明示的に拒否する。アプリが認証されていない BLASTN の得点を出さないよう、S12 の前に終える | S06 で、NCBI が拒否する 34 の組合せを LOSAT が実行し、NCBI が受け付ける 17 の組合せで結果が違うことが分かった（`AUTHORITY.md` §D.5）。outfmt 0 の移植（S07）とは原因が別なので、段階を分ける |
| TD-11 | reactor は、ビルドした checkout のパスと `CARGO_HOME` を `--remap-path-prefix` で固定の名前に置き換えてビルドする。TD-6 の rustflags の比較は、この置き換えだけを差として認める。S07 で入れ（reactor を作り直し、full の V-ABI を再実行するため）、S16 の V-PRIV で公開物にローカルのパスが無いことを確かめる | アダプタは LOSAT を path 依存として使うので、エンジンのコードの panic の位置に checkout の絶対パスが入り、登録簿の crate の位置には `CARGO_HOME` が入る（S05 の独立監査、`docs/web/abi_v2.md` §2）。公開物にビルドした機械のパスを残さない。置き換えはコードの働きを変えない。S07 で、置き換えの後の reactor にローカルのパスが無いことを確かめた。バイトは checkout のパスにはなお依存する（path 依存の crate の Cargo の metadata hash が記号の名前に入る）ので、同一性の記録に checkout のパスを残し、再現は同じパスで行う |
| TD-12 | BLASTN の入力は、LOSAT の `bio` の読み方と NCBI の `CFastaReader` の読み方が同じものだけを受け付ける。違いは `U`（NCBI は `T` として読む）だけを移植し、ほかの違い（定義行の先頭の空白・制御文字・非 ASCII、IUPAC の文字以外の残基）は明示的に拒否する。v2 のアダプタは `register` で同じ検査をする。NCBI の読み込み器の移植は未決事項にした（§10） | NCBI の読み込み器の移植には、アダプタの索引の走査（TD-8）に新しい解析器の種類が要り、S07+ の範囲を超える。拒否する入力は研究で使う FASTA では稀で、黙って違う結果を出すより安全である（S07+、`docs/evidence/losat_web_e2c/AUTHORITY.md` §E） |
| TD-13 | BLASTP・TBLASTN・TBLASTX の既定以外のオプション（行列、gap、threshold、word size、window、`-comp_based_stats`、`-seg`、遺伝暗号、e-value の書き方など）は、S08+ で NCBI と同じにするか、明示的に拒否する。BLASTX は SX で同じことを行う。S12 の前に終える | 認証済みの fixture は既定の値だけを使う（`LOSAT/tests/blastp_parity_options.sh` など。例外は承認済みの遺伝暗号だけ）。S07+ で BLASTN の既定以外の値に多くの差が見つかった（`docs/evidence/losat_web_e2c/AUTHORITY.md`）が、ほかの program は比べていない。これらのオプションはアプリの検索画面に出るので、S07+ と同じ扱いにする（2026-09-29、推奨案で進める保守者の指示による）|
| TD-14 | BLASTN は、複数の query を NCBI と同じ query の batch（`CBatchSizeMixer`、`good_init_extends` による後の batch の大きさ）に分けて検索する。S07++ で移植し、S07+ の batch に依存する明示的な拒否をなくす。S12 の前に終える | S07+ の第 7 回の独立監査で、NCBI の近似の ungapped 伸長が query の塊（batch）の端で止まるため、batch の境の query の端で、既定のオプションでも LOSAT の HSP が増えることが分かった（`docs/evidence/losat_web_e2c/AUTHORITY.md` §K）。以前からの差で、S07+ の範囲（得点のオプションと入力の読み方）の外の、エンジンの構造の変更なので、段階を分ける（2026-09-30、推奨案で進める保守者の指示による）|
| TD-15 | BLASTN の `BATCH_SIZE`・`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`・`PRE_FETCH_SEQS_LIMIT` の整数でない値は、明示的に拒否する（NCBI の `NStr::StringToInt` の失敗と同じ値を、NCBI の誤りの文言を再現せずに止める） | NCBI の誤りは `CStringException::what()` で、オラクルの build のソースのパスと行の番号（`/opt/conda/conda-bld/…/ncbistr.cpp`）を含み、NCBI の公式の build とも違う build に固有の文字列なので、バイト一致の対象にできない（S07+++b の技術判断。終了コードも NCBI の 255 と違う 1） |

---

## 1. 現状とのギャップ

`0627c88f5`（`origin/main`）時点のコードを調べた結果である。行番号はこの commit のものである。

| # | 設計書の要求 | 現状（根拠） | 必要な作業 |
|---|---|---|---|
| G1 | 5 program（`REQ-01`） | 5 program すべてが `main` にある。BLASTX は PR #106 で統合されたが、BLASTX v0.2.0 の受入判定は HARD_FAIL、リリースは HOLD のままである（merge commit のメッセージと `docs/evidence/losatx_stage_g_authority_v3/`） | BLASTX の認証は LOSATX 計画で完了させる。このブランチでの BLASTX の扱いは、その後に行う（DW-10、SX） |
| G2 | 全 program で outfmt 0（`REQ-01`） | BLASTN は 6/7 だけ（`LOSAT/src/algorithm/blastn/hsp.rs:30-43`）。TBLASTX は 6 だけ（`LOSAT/src/blastinput/value_parsers.rs:313-320`）。BLASTP・TBLASTN・BLASTX は 0/6/7 | BLASTN の outfmt 0、TBLASTX の outfmt 0/7 を移植する（DW-6） |
| G3 | 1 回の検索から複数の出力（設計書 §10.1） | どの入口も、1 回の実行で 1 形式だけを書く（`LOSAT/src/web_api.rs:614-760`）。出力の書き先も program ごとにばらばらである。TBLASTX は少なくとも 4 か所から書き、そのうち 1 つは writer thread（`LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1352`）、1 つは `wasm-threads` のときだけコンパイルされる（同 `:3061-3064`） | 出力を 1 か所に集め、1 回の結果を複数の formatter に配る（§4.3） |
| G4 | 構造化結果（設計書 §10） | program をまたぐ共通型は `common::Hit`（`LOSAT/src/common.rs:78-183`）だけで、subject frame・整列文字列・ID を持たない。`report::PairwiseHit`（`LOSAT/src/report/pairwise.rs:76-102`）は outfmt 0 の入力で、BLASTP・TBLASTN・BLASTX が outfmt 0 のときだけ作る。BLASTN と TBLASTX は作らない | 全 program が、並べ終えた最終の HSP 一覧から `PairwiseHit` を作る（§4.4） |
| G5 | 入口の構成 | BLASTP にはレコードを受け取る共通の経路があり、CLI と v1 はそこを通る（`LOSAT/src/algorithm/blastp/blast_engine.rs` の `run_resolved_with_records`）。ほかの program のメモリ上の入口は wasm32 のときだけコンパイルされる（例：`blastn/blast_engine/run.rs:4349`、`tblastx/blast_engine/run_impl.rs:572`）。TBLASTN にはファイルを読む `run(self)` しか無い（`LOSAT/src/algorithm/tblastn/args.rs:382`）。BLASTX は 10,002 文字ごとの query バッチで検索と整形を繰り返す（`LOSAT/src/algorithm/blastx/web.rs:100-147`、`query_setup.rs:81`） | ターゲットに依存しない 1 本の入口 `run_local` にまとめる（TD-2、§4.2） |
| G6 | 警告 | 警告はプロセスの stderr に直接書かれる。BLASTX は formatter の中で query ごとの警告を書く | 警告の書き先を 1 つにし、形式を増やしても 1 回だけ出す（§4.3） |
| G7 | Subject の保持（`REQ-07`） | 呼出しをまたいで保持するのは、解析済み FASTA のハンドルだけで、対応するのは BLASTP と BLASTX（`web_api.rs:695-760`）。過去の Wasm 性能計画には「永続 Worker pool や検索間の可変状態 cache を混ぜない」とある（`docs/wasm_remaining_work_implementation_plan_20260915.md:86`） | 出力バイトの同一性を条件に、PD で許可した（§4.6） |
| G8 | 取消と進捗（`REQ-08`） | エンジンにも ABI にも無い。gbdraw は Worker を終了して取消している | 初期は Worker の終了で取消す（§2.3）。進捗は段階と経過時間だけ（DW-5） |
| G9 | 領域の指定（`REQ-06`） | BLASTN・BLASTP・TBLASTX は `-query_loc` / `-subject_loc` を定義しておらず、CLI が未知の引数として拒否する。TBLASTN は `-query_loc` を未移植の名前として拒否し（`LOSAT/src/cli.rs:187-193`）、`-subject_loc` も拒否する（`LOSAT/src/algorithm/tblastn/args.rs:364`）。BLASTX は v0.2.0 の範囲外として両方を拒否する（`cli.rs:238-246`） | 4 program へ移植する（S11）。BLASTX は SX で移植する（DW-11） |
| G10 | 引数（設計書 §11.2） | BLASTX の Web 入口は、CLI と同じ clap パーサに argv を渡す（`LOSAT/src/algorithm/blastx/web.rs` の `parse_args`）。BLASTN・BLASTP・TBLASTX の v1 は、手書きの小さな引数パーサと、それ自身の既定値を持つ（`web_api.rs:227-610`）。TBLASTN は Web に公開されていない。TBLASTX の v1 は outfmt を検証せず常に 6 を出していた（S01 で修正済み） | v2 は argv と clap で統一する。v1 のパーサは凍結する（TD-1） |
| G11 | ブラウザでの並列実行（設計書 §9.1） | rayon のプールは検索ごとに作られ、N-1 本を thread-spawn する（`LOSAT/src/utils/threading.rs:93`）。serial のビルドは `-num_threads` が 1 を超えると拒否する（同 `:64-77`）。ブラウザ用の host は gbdraw にしか無い | スレッド版 reactor 用のブラウザ host を新設する（S09）。serial に切り替えるときは argv の `-num_threads` を 1 にする（§4.7） |
| G12 | メモリ（設計書 §9.2） | スレッド版の共有メモリは最大 16384 ページ（1 GiB）。reactor を繰り返し呼ぶとメモリが増える問題が未解決（`docs/wasm_reactor_memory_followup_20260914.md`） | 最大値は据え置き（TD-7）、実測する。instance を作り直す方針を持つ（§5.5） |
| G13 | 検証の表示（`REQ-18`） | v0.2.0 のリリース判定では browser/reactor ABI は認証の範囲外（`docs/release/v0.2.0.md`）。既存の manifest は 1 件につき 1 形式だけを固定している | 升目ごとに期待値の出どころを決める（TD-5、§6） |
| G14 | フロントエンド | リポジトリにフロントエンドもブラウザの CI も無かった | S01 で `web/app` と `.github/workflows/web.yml` を新設した |

**変えてはいけない制約**：BLASTN・BLASTP・TBLASTX・TBLASTN は、subject を外側に回すループで検索し、全検索が終わってから query 順に出力する。BLASTX は query バッチごとに検索と整形を繰り返す。TBLASTX は全 query の平均長を使い（`run_impl.rs:1003-1005`）、TBLASTN は 20,000 残基の query バッチ単位で統計を取る（`tblastn/args.rs:468-472`）。したがって、**query や subject を分割して独立に検索し、結果をつなぐことはしない**（結果が変わる）。

---

## 2. SOLID / KISS / DRY / YAGNI の適用

### 2.1 アーキテクチャでの適用

| 原則 | 決定 |
|---|---|
| **S** 単一責任 | 変更の理由ごとに 3 つに分ける。`LOSAT/`：NCBI パリティ（変わる理由は NCBI との差だけ）。`web/adapter/`：ABI の契約（ポインタの受け渡し、ハンドル、直列化、FASTA 索引）。`web/app/`：利用者の体験（画面、ジョブ、保存、書き出し） |
| **O** 開放閉鎖 | program の違いは `ProgramDescriptor` の表に書く。program を足すときは、descriptor 1 行と試験 manifest 1 つを足すだけにし、画面・Worker・書き出しのコードは変えない |
| **L** 置換可能性 | `BlockStore`（OPFS / Memory）と `EngineRuntime`（threaded / serial）は、同じ契約試験を通る実装どうしとして入れ替えられる。「同じ入力なら同じバイト」を契約試験で保証する |
| **I** インターフェース分離 | UI は wasm のポインタも Worker のメッセージも知らない。アプリ層が使うのは `EngineGateway`・`DataGateway`・`Downloader` の小さな interface だけ（計測を入れるときに `Analytics` を足す） |
| **D** 依存性逆転 | アプリ層は上の interface にだけ依存し、実装は起動時に組み立てる。試験では FakeEngine などを差し込む。これにより、エンジン側が出来る前からアプリ側を進められる |
| **KISS** | 物理構成は UI・Data worker・Engine worker（＋スレッド用 worker）の 3 種類だけ。検索条件の表現は LOSAT CLI の argv ひとつ。取消は Worker の終了だけ。進捗は段階と経過時間だけ |
| **DRY** | ① program ごとに核の入口は `run_local` の 1 本（TD-2）。② argv ひとつで、実行・入力検証・CLI コマンドの表示・セッションをまかなう。③ 既定値・選択肢・help は clap の定義から取る。④ 検索に渡す FASTA は各 program のエンジンの解析器が読む（索引用の走査は TD-8 の例外）。⑤ NCBI 形式の数値整形と詳細アラインメントは、既存の formatter の出力をそのまま使い、TS で作り直さない。⑥ 応答ヘッダーは `web/app/public/_headers` の 1 か所で定義する |
| **YAGNI** | 結果アーカイブ形式、協調取消、部分読み出し、前処理キャッシュ、形式ごとの別の検索、共有メモリの拡大、複数の計算プール、query ごとの途中経過は作らない。そのうち前の 6 つは §2.3 のトリガーが立ってから入れる |

**`web/` の中の線引き**：BLAST が定義する値（スコア、E 値、identity、被覆率、NCBI 形式の数値やアラインメントの整形）を計算・整形するコードは `web/` に置かない。必要なものはエンジンから受け取る。こうすると、NCBI 参照コメントの要らないアプリのコードと、NCBI が権威となるコードとが混ざらない。

### 2.2 ワークフローでの適用

| 原則 | 決定 |
|---|---|
| S | 1 セッションで扱うのは 1 段階で、成果は 1 つのゲート記録。エンジンのパリティ作業とアプリの作業を同じコミットに混ぜない。エンジンの中でも、入口のまとめ直し（出力を変えない変更）と移植（出力を変える変更）を別のセッションにする |
| O | 新しい program や実行経路は、manifest と検証の升目の表に行を足すだけで検証の対象に入る |
| KISS | ゲートの種類は「バイトの同一性」と「E2E 試験の合格」の 2 つに揃え、エンジンのまとめ直しには性能の非退行を足す。ブランチと worktree は、エンジン側とアプリ側の 1 つずつ（DW-7） |
| DRY | 共通の手順（ブランチの確認、コミット、引き継ぎ、レビュー）はセッション README に一度だけ書き、計画書と各指示書はそれを参照する。期待値は既存の凍結バイトを参照し、コピーしない |
| YAGNI | 各セッションの指示書には、目的・入口の条件・作業・完了条件だけを書く。実測に依存する細部は、前のセッションが引き継ぎのときに書き足す |

### 2.3 簡素化と後から入れる条件（トリガー）

| 設計書の提案、または起こり得る要求 | 初期の実装 | 後から入れる条件 |
|---|---|---|
| 共有メモリ上の協調取消フラグ（設計書 §9.4） | Engine worker とスレッド用 worker をすべて終了する。Data worker の保持物は残す | 取消の後の再準備（instance の生成と Subject の再登録）の実測時間が、保守者の決める閾値を超えたとき。そのときは NCBI の `TInterruptFnPtr` を移植する（§4.9）。S09 の実測：取消そのものは 0〜1 ms、取消の後の再準備は最大 0.61 s（Firefox、4.6 Mb の Subject と 5.5 Mb の query）、Chromium と WebKit は 0.14 s 以下で、Subject の長さにほぼ比例する。閾値は 1 秒（DW-22）：公開する対応規模の最大の入力で、取消の後の再準備が 1 秒を超えたとき（S17 で測り直す） |
| SequenceProvider、部分読み出し、巨大な単一レコードへの範囲供給（設計書 §6.3、§8.2） | 選択したレコードの原バイトを、wasm メモリへまとめて渡す。扱える規模の上限は実測して公開する | 目標とするデータ規模が、実測したメモリ上限を超えるとき |
| 共有メモリの最大値の拡大 | 1 GiB のまま（TD-7） | 目標とする入力の実測で 1 GiB が足りないとき。そのときは `--max-memory` だけを link の引数の差として認め、成果物の JSON に記録する |
| PreparedSubject（エンコード・翻訳のキャッシュ）（設計書 §5、§8.1） | 解析済みのレコードと warm な instance を保持する（§4.6 の R1） | DW-8：ブラウザでの実測で、前処理が warm 実行時間の 20% 以上を占めた program |
| 形式ごとに解決した検索オプションが違う場合の、形式ごとの別の検索 | `run_local` が明示的なエラーで止める（TD-4） | いずれかの program が `-num_descriptions` / `-num_alignments` を移植したとき |
| 構造化バッチ、結果アーカイブからの再整形（設計書 §10.1） | 実行の完了時に、その program が対応する形式と構造化レコードを同時に作り、保存する | 生成していない形式（例：別のカスタム列）を後から出す要求が確定したとき |
| 一時領域の staging から committed への移動（設計書 §7.1） | 移動しない。Data worker の登録簿で「確定」かどうかを区別する | 無し |
| 残存領域の回収（設計書 §7.2） | Web Locks の保持状態だけで判断する。時刻では判断しない。Web Locks が使えない環境では回収しない | 無し |
| レコード表の保存（設計書 §7.1 の `datasets/<revision>/index.blocks`） | Data worker のメモリに置く（S10）。S12 の実測：10 万レコード（34 MB）の索引は 3〜10 秒、表は 1 つの写しで約 51 MiB（V8 での見積り）。負担にならないと判断した | S17 で公開する対応規模が 10 万レコードを大きく超えるとき、または V-MOB の実機で負担が見えたとき |
| Memory の BlockStore の上限 | 1 つのタブで 512 MB（S09。`MEMORY_RESULTS_CAPACITY_BYTES`）。超えた実行は理由付きで失敗し、前の結果は残る | 実機（V-MOB）で 512 MB が合わないとき |
| 進捗率と残り時間（設計書 §3.1） | 出さない。段階と経過時間だけを出す | エンジンが測定できる進捗を持ったとき |
| 図と一覧の集約表示（設計書 §11.1） | 仮想化した一覧と Canvas での描画だけ | 描画件数の実測値が、保守者の決める閾値を超えたとき |
| 複数の計算プール、独立した直列 Wasm の並列実行（設計書 §9.1） | 作らない | 作らない（単一 query の高速化には効かないため） |
| 別 origin の GA4 計測文書（設計書 §13.2） | S16 まで作らない。CSP は `'self'` だけにしておく | S16 で方式を比べて決める（§5.10） |

---

## 3. 全体アーキテクチャ

### 3.1 構成要素と責任

```text
Main thread ─ Vue UI（表示だけ）
            └ Application coordinator（Vue に依存しない TS）
                 Draft / Queue / Runs / 結果の索引 / ViewState / Candidates
                 │ MessagePort                         │ MessagePort
                 ▼                                     ▼
Data worker（作業中はずっと動く）          Engine worker（取消のときに終了し、作り直す）
  ├ SourceStore：File の参照、レコード表     ├ adapter wasm（threaded / serial の reactor）
  │   └ serial の adapter wasm で索引を作る  ├ WASI shim と ThreadHost ─► スレッド用 worker（再利用する）
  ├ BlockStore：OPFS | Memory               └ 結果のチャンクを Data worker へ直接送る
  ├ Run の登録簿（仮 / 確定）、GC              （worker どうしの MessagePort）
  └ Web Lock（作業ごとの名前空間）
Service Worker：release を固定した資産キャッシュだけ
```

- 検索の実行中でも次のジョブのためにファイルを開けるように、FASTA の索引は Data worker で作る。そのため、Data worker は serial の adapter instance を別に持つ。
- Engine worker は、作業状態の唯一のコピーを持たない（設計書 D01）。終了させても、失うのは warm な instance とそこに登録した Subject だけで、どちらも Data worker から作り直せる。

### 3.2 コード配置

```text
web/
  AGENTS.md                  # web/ 以下の規約
  adapter/                   # S05 で新設。Rust cdylib。ABI v2。wasm32-wasip1 と wasm32-wasip1-threads の reactor
    Cargo.toml               # LOSAT = { path = "../../LOSAT", default-features = false }、features は parallel / wasm-threads を明示
  app/                       # Vite + Vue 3 + TypeScript（S01 で新設）
    public/_headers          # 応答ヘッダーの唯一の定義（COOP / COEP / CORP / CSP）
    build/headers.ts         # _headers を読み、Vite と試験に渡す
    src/
      domain/                # program、argv、出力形式、run の型（Vue 非依存、何にも依存しない）
      application/           # coordinator（キュー、状態機械、グループ）、draft（検索画面の下書き、S12）、
                             # attention（Wake Lock、離れる前の注意、復帰時の確認、S12）、Store（Vue 非依存）
      ports/                 # EngineGateway / DataGateway / Downloader / InputChecker / PagePort の interface
      infra/                 # fake/（FakeEngine、FakeScanner、describe.json）、data/（BlockStore の OPFS / Memory、DataService、セッション）、
                             # data-worker/、run-output/（Engine worker から Data worker への出力の経路）、browser/（ダウンロード、SHA-256、ページ）、
                             # engine-worker/（Engine worker、ThreadHost、Auto と作り直しの規則、S09）、reactor/（ABI v2 の結合、入力の検査、S09・S12）
                             # 後で sw/ を足す
    build/reactors.ts        # LOSAT_WEB_REACTORS の reactor をビルドに入れる（無ければ FakeEngine）、threaded の module の共有メモリの guard（S09）
      ui/                    # Vue コンポーネント
      composition.ts         # 実装を選んで組み立てる唯一の場所
    tests/  unit/ contract/ e2e/（harness/、support/）
  tools/                     # Node で wasm32-wasip1 の試験を動かす runner など
docs/web/                    # 設計書、要求トレース表、ABI v2 の契約
docs/product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md
docs/evidence/losat_web_<stage>/
```

使う Wasm のターゲットは、認証済みの command-WASI と同じ系統の `wasm32-wasip1` と `wasm32-wasip1-threads` だけにする。`wasm32-unknown-unknown` と wasm-bindgen は使わない。

### 3.3 依存の向き

`ui → application → ports ← infra` とし、`domain` は他のどこにも依存しない。`ui` は `application` を通して操作し、表示のために `domain` と `ports` の型と定数だけを import できる。`src/main.ts` と `src/composition.ts` だけが全層を import できる。`infra/engine-worker → web/adapter → LOSAT` の向きにする。これらは `web/app/eslint.config.js` の `no-restricted-imports` で検査し、`web/app/tests/unit/layers.test.ts` がその規則自体を試験する。

---

## 4. エンジン側の契約

### 4.1 方針

エンジンへの変更は、Web に必要な「核の入口」「出力の配布」「観測者」「警告の書き先」と、DW-6・S11・SX の移植に限る。検索の本体と、認証済みの formatter のバイト出力は変えない。変更はすべて、ルートの `AGENTS.md` の規約（NCBI 参照コメント、パリティのゲート、独立監査）と `.agents/skills/verify-ncbi-parity-and-speed/SKILL.md` に従う。これらの LOSAT 側の配管は、`PD-LOSAT-WEB-APP-BOUNDARY` が名前を挙げて許可している。

### 4.2 核の入口 `run_local`（TD-2）

`LOSAT/src/api/local_blast.rs`（今は再 export だけ）に、NCBI の局所検索（subject 集合を一度用意し、query バッチごとに `CLocalBlast` を実行して整形する、§4.9）に相当する入口を置く。

S02 で次の形に確定した（`LOSAT/src/api/local_blast.rs`）。

```rust
pub enum OutputSink<'a> {
    Stdout,                                  // CLI（-out なし）。バッファする
    File(&'a Path),                          // CLI（-out）。整形を始めるときに作る
    Writer(&'a mut (dyn Write + Send)),      // 呼び出し側の writer。ここではバッファしない
}
pub struct FormatOutput<'a> { pub outfmt: &'a str, pub sink: OutputSink<'a> }
pub type HspIndex = usize;                   // 最終の HSP 一覧の中での番号。全形式で共通の ID
pub trait FormatObserver {                   // 行・節の開始と終了。書くバイトは変えない
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex);
    fn hsp_end(&mut self, format: usize, hsp: HspIndex);
}
pub struct ReportOutputs<'a> {               // すべて Send（整形はスレッドプールの中で走る）
    pub formats: Vec<FormatOutput<'a>>,                                // CLI は 1 つ
    pub diagnostics: &'a mut (dyn Write + Send),                      // CLI の stderr に当たる警告
    pub hits: Option<&'a mut (dyn FnMut(&[PairwiseHit]) + Send)>,     // 最終の HSP 一覧
    pub observer: Option<&'a mut (dyn FormatObserver + Send)>,
}

// program ごとの入口（S02 で BLASTP、S03 で TBLASTN、S04 で BLASTN と TBLASTX。BLASTX は SX）
pub fn run_local(args: BlastpArgs, queries: &[fasta::Record], subjects: &[fasta::Record],
                 query_label: &str, subject_label: &str, outputs: &mut ReportOutputs<'_>) -> Result<()>;
pub fn run_local(args: TblastnArgs, queries: &[fasta::Record], subjects: &[fasta::Record],
                 outputs: &mut ReportOutputs<'_>) -> Result<()>;   // api::local_blast::run_local_tblastn
pub fn run_local(args: BlastnArgs, queries: &[fasta::Record], subjects: &[fasta::Record],
                 outputs: &mut ReportOutputs<'_>) -> Result<()>;   // api::local_blast::run_local_blastn
pub fn run_local(args: TblastxArgs, queries: &[fasta::Record], subjects: &[fasta::Record],
                 outputs: &mut ReportOutputs<'_>) -> Result<()>;   // api::local_blast::run_local_tblastx
```

BLASTN と TBLASTX は、S07・S08 で `PairwiseHit` を作るまで `hits` を呼ばない。観測者は outfmt 6/7 の行について知らせる。

BLASTP の `query_label` / `subject_label` は v1 のためにある（空なら `-query` / `-subject` の値を使う）。ほかの program の入口は、表示名を `-query` / `-subject` の値からだけ取り、ラベルの引数を持たない（BLASTN と TBLASTX の v1 はラベルを渡していなかった）。TBLASTN の CLI は、エラーのときに途中までの出力を残さないために、`OutputSink::Writer` でメモリに書き、成功したときだけ stdout / `-out` に書く（以前と同じ）。

program をまたいで振り分ける関数は、それを呼ぶアダプタを作る S05 で足す（使う側の無い振り分けを先に作らない）。

- CLI の `run` は「ファイルを読んで `run_local` を呼び、stdout / `-out` に書く」だけになる。v1 の `run_web_pair*` も `run_local` を呼ぶ。こうして program ごとの核の経路を 1 本にする。
- `run_local` はターゲットに依存しない。そのため `cargo test` で、CLI と同じ経路の出力を fixture で確かめられる。`wasm-threads` のときだけコンパイルされる TBLASTX の出力箇所（G3）は、ネイティブでは通らないので、threaded の V-ABI で確かめる。
- `-query` / `-subject` の値は表示名にだけ使う（outfmt 7 の `# Database:` 行などに出る）。CLI で比較するときも同じ名前を渡す。
- まとめ直しは出力を変えない変更である。完了条件は、変更の前に取った全 program の基準（凍結バイトを優先する）との一致、v1 の reactor の検査（`LOSAT/tests/check_wasi_reactor.js`、`check_wasi_api_limits.js`）、性能の非退行である（§6.2 の V-PERF）。共有の formatter（`LOSAT/src/report/`）を変えたときは、それを使うすべての program のゲートを実行する。
- BLASTX の入口と、BLASTX だけが使う formatter の経路は、SX まで変えない（DW-10）。BLASTX と共有する関数を変える必要があるときは、関数の本体と引数を変えずに、呼出し側か新しい関数で観測者を受ける。

### 4.3 1 回の検索から複数の出力（TD-4）

- 形式ごとに検索オプションを解決する。解決した検索オプションがすべての形式で同じなら、1 回の検索の結果を、形式ごとの formatter に順に渡す。違う形式があれば、`run_local` は明示的なエラーで止める（§2.3）。
- 形式ごとにしか効かない表示のオプションは、その formatter にだけ渡す。
- 警告は `diagnostics` に 1 回だけ書く。formatter の中で書かれている警告（BLASTX）は、SX で、形式の数だけ重ならないように整理する。
- 根拠：NCBI は 1 つの結果（`CSearchResultSet`）を、再検索せずに `CBlastFormat` で整形できる（`blast_formatter`）。複数の形式は、同じ結果に対する複数の `CBlastFormat` に当たる（§4.9）。
- ゲート：どの形式のバイトも、その形式だけを指定した CLI の実行のバイトと一致すること。

### 4.4 構造化結果：`PairwiseHit`

新しい結果の型は作らない。全 program が、並べ終えた最終の HSP 一覧から `PairwiseHit` を作り（BLASTN と TBLASTX は S07・S08 で作るようにする）、それを構造化結果にする。各 HSP の ID は、その実行の最終の HSP 一覧の中での番号（`HspIndex`）で、どの出力形式でも同じ HSP を指す（S02 で確定）。ABI の HSP レコードは、この番号と、そこから求めた `q_idx`・`rank` を持つ。アダプタが直列化する内容は、[`docs/web/abi_v2.md`](web/abi_v2.md) の「HSP record」にある。

**表に出す値は、outfmt 6 の該当行を分割して使う。** NCBI 形式の数値整形を TS で作り直さないためである。並べ替えには原値を使う。

### 4.5 formatter の観測者（TD-3）

formatter は、HSP の行（outfmt 6/7）や節（outfmt 0 のスコアの行とアラインメント。subject の見出しは含まない）を書き始めるときと書き終えるときに、その HSP の ID を観測者に知らせる。書くバイトは変えない。アダプタは、そのときの書き込み位置から、HSP ごとのバイト範囲を記録する。outfmt 0 に現れない HSP（例：BLASTX の既定で 251 番目以降の subject）は、範囲を持たない。詳細画面には、選んだ HSP の subject 見出しと節を outfmt 0 の原文のまま表示する（見出しの範囲を得る方法は S05 で決める。S02・S03 の観測者は見出しを節に含めない）。midline（`|` や `+`）やマスクの表示規則を TS で作り直すことはしない。

### 4.6 Subject の保持

- **R1（初期）**：アダプタのハンドル表に、登録した Subject を保持する（v1 の `FASTA_STORE` と同じ考え方を、全 program 共通にする）。コンパイル済みの `WebAssembly.Module` と warm な instance も保持する。instance を失ったら、Data worker が原バイトから登録し直す。
- **R2（DW-8 のとき）**：program ごとの前処理キャッシュ（エンコード済みの subject など）。NCBI では、subject 集合から `CLocalDbAdapter` を一度作り、query バッチごとに `CLocalBlast` を作り直す構成がこれに当たる（§4.9）。query に依存する構造（例：全 query から作る lookup）は再利用しない。
- **ゲート**：同じ Subject に対して query や条件を変えて連続で実行したとき、どの出力も、毎回新しく実行した CLI の出力と一致すること。

### 4.7 取消・進捗・serial への切り替え

- 取消：初期は Engine worker とスレッド用 worker を終了するだけ（§2.3）。
- 進捗：host から段階（準備・検索・整理）を観測し、経過時間と一緒に出す。query ごとの途中経過は出さない（DW-5）。
- serial への切り替え：threaded を使えないとき、host は argv の `-num_threads` を 1 にして serial のモジュールで実行し、理由を RunRecord に書く。RunSnapshot には利用者が要求した値を残す（スレッド数は出力を変えない）。

### 4.8 ABI v2

契約は [`docs/web/abi_v2.md`](web/abi_v2.md) にある（S01 で下書き、S05 で確定）。`describe` は program ごとに対応する出力形式を返し、`run` はその形式だけを出す。アプリ側から見た型は `web/app/src/ports/engine.ts` で、HSP レコードの項目を変えるときは同じコミットで両方を変える。成果物は `losat-web-serial.wasm` と `losat-web-threads.wasm` の 2 つで、ビルドの条件は TD-6、共有メモリの最大値は TD-7 のとおりにする。

### 4.9 NCBI 側の拠り所（固定 commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`、`/mnt/c/Users/genom/GitHub/ncbi-blast/`）

| 用途 | NCBI の箇所 |
|---|---|
| subject 集合を一度用意し、query バッチごとに検索して整形する | `c++/src/app/blast/blastx_app.cpp:197-297`（`CLocalDbAdapter` を一度作り、`PrintProlog`、バッチごとの `CLocalBlast(...).Run()` と `PrintOneResultSet`、最後に `PrintEpilog`） |
| 再検索せずに結果を整形する | `c++/src/app/blast/blast_formatter.cpp:429-465`、`c++/src/algo/blast/format/blast_format.cpp:1411` |
| hitlist の大きさと出力形式 | `c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978` |
| 範囲の指定 | `c++/src/algo/blast/blastinput/blast_args.cpp:1946-1997`（`-query_loc`）、`:2373`（`-subject_loc`）。範囲は入力元が読むすべてのレコードに適用される（`c++/src/algo/blast/blastinput/blast_fasta_input.cpp:433-459`） |
| 中断コールバック（トリガーが立ったときだけ） | `c++/include/algo/blast/core/blast_def.h:335-354`、`c++/src/algo/blast/core/blast_engine.c:571` |
| Subject 一覧の集約列（採用した列だけ） | `c++/src/objtools/align_format/showdefline.cpp:89`、`:1215` |

実装するセッションは、これを起点に呼出し経路を読み、Rust の直上に引用を置く。ここに挙げた行は起点であり、それだけで移植の根拠になるという意味ではない。

---

## 5. アプリ側の設計

### 5.1 レイヤー

- `domain`：純粋な型と関数。座標（1 始まりで両端を含む、strand、単位）の変換はここ 1 か所に集める。
- `application`：coordinator と状態機械。副作用は `ports` を通してだけ起こす。
- `infra`：Worker、wasm、OPFS、Service Worker。
- `ui`：Vue。application の状態を表示し、操作を application に渡すだけ。UI・ヘルプ・エラーは英語（`REQ-02`）。

### 5.2 ドメインモデル（設計書 §5 を簡約したもの）

| オブジェクト | 内容 | 設計書からの変更 |
|---|---|---|
| SourceRef | File の参照、または貼り付けたテキスト | そのまま |
| DatasetRevision | source、レコード表、選択・除外、レコードごとの SHA-256。領域は、その役割のレコードが 1 つのときだけ持つ（DW-9） | 不変。編集すると新しい revision になる |
| RunSnapshot | program、正規化した argv、入力の名前・バイト・SHA-256、要求スレッド数、グループ ID | argv を中心に据えた（S01 で実装） |
| RunRecord | 実際の経路（threaded / serial）、実スレッド数、切り替えの理由、エンジンの build、段階の時刻、エラー | そのまま（S01 で実装） |
| ResultSet | 形式ごとのテキストのブロック参照、HSP レコード、query ごとの状態（ヒット無し、上限到達の可能性） | そのまま |
| ViewState、Candidate | 設計書のとおり | そのまま |
| PreparedSubject | アプリのオブジェクトとしては持たない。Engine worker の中のハンドルだけ | R1 に合わせて縮小 |
| RunState | 状態は §5.5 のとおり | 途中経過を持たない（DW-5） |

Combined は、複数ファイルのレコードをつないだ 1 回の実行にする。Separate は、ファイルごとの RunSnapshot を同じグループ ID で積む。「グループを取消」は、実行中のものを取消し、同じグループで待機中のものを除く操作とする（原子性は求めない）。

### 5.3 検索条件の表現は CLI の argv ひとつ

- フォームは、`ProgramDescriptor`（program ごとに、表示する引数・節・ラベルを並べた表）と、`describe` の JSON から組み立てる。既定値・選択肢・help・対応する出力形式はエンジンから取る（DW-4）。
- 「Add to queue」を押したとき、`validate(argv)` が通った場合にだけ RunSnapshot を確定する。エラーは CLI と同じ文で出す（S01 で実装）。
- この argv を、実行、CLI コマンドの表示、NCBI との比較用コマンド、セッションで共通に使う。`-query`・`-subject`・`-out`・`-outfmt`・`-num_threads` はアプリが管理し、利用者はパラメーターとして入力できない（S01 で実装）。
- 表示する CLI コマンドは、出力形式ごとに作り、必ず `-outfmt N` を付ける。いくつかの program は CLI の既定の形式を受け付けないためである（S01 で実装）。
- `-query` / `-subject` の名前の規則：ファイルが 1 つならそのファイル名。貼り付けなら `query.fa` / `subject.fa`。複数ファイルをつないだ場合は `combined_query.fa` / `combined_subject.fa`（S12）。つなぐとき、末尾に改行の無いファイルには改行を 1 つ補う。
- レコードを除外した実行やファイルをつないだ実行では、「この実行に使った入力 FASTA」（エンジンに渡したバイトそのもの）を書き出せるようにする。CLI での再現には、このファイルを使う。
- 領域は、その役割のレコードが 1 つのときだけ指定でき、`-query_loc` / `-subject_loc` として argv に入る（DW-9）。

### 5.4 入力と索引

- 元の File は読むだけで、コピーしない。Data worker は `File.slice` で必要な範囲だけを読む。
- 検索に渡すのは、選んだレコードの原バイトをそのまま並べたものである。検索の入力は、その program のエンジンの解析器が読む（`bio::io::fasta`、BLASTX は `LOSAT/src/algorithm/blastx/input.rs`）。
- 抽出のための索引は、アダプタの `scan` が作る（レコードの位置、行の配置、残基の内訳。TD-8）。行長が揃っていれば計算式で、揃っていなければ 64 Ki 残基ごとのチェックポイントで位置を引く。`scan` は解析器の種類を受け取り、その解析器の結果との一致を性質試験で確かめる。
- `register` の時点で、エンジンが解析したレコードの ID と長さを `scan` の表と照合し、食い違えば止める。索引の誤りが、抽出の誤りとして表に出ないようにするためである。
- 不正なレコードかどうかは、エンジン自身の入力変換で判定する（アダプタが、レコードごとのエラーを返す）。TS で検証規則を作り直さない。配列種別の警告（例：「タンパク質に見える」）だけはアプリの推定で出し、program は変更しない。

### 5.5 実行・キュー・状態機械

- 状態は `Queued → Preparing → Running → Finalizing → Completed | Cancelled | Failed`。実行は常に 1 つで、待ち行列は FIFO（S01 で実装）。
- 取消と確定の順序：実行は `Finalizing` に入るまで取消せる。`Finalizing` の後の取消は受け付けず、その実行は完了する（S01 で実装）。
- 段階のイベントは、実行中で取消されていない実行のものだけを受け付ける。S09 で Worker を作り直すようになったら、イベントに `runtimeGeneration` を付け、古い世代のものを捨てる。
- Auto のスレッド数（S09 で決めた）：入力（query と subject の FASTA）が 20,000 バイト未満なら serial、それ以上なら論理プロセッサの半分（最大 4）。手動の指定はそのまま使う（`src/infra/engine-worker/policy.ts`）。
- instance の作り直し（S09 で決めた）：検索の後の linear memory が 512 MiB（threaded の最大値の半分）以上なら作り直す。実行回数では作り直さない（G12 の対策。linear memory は最初の 1〜3 回の検索で水準に達し、増え続けなかった）。
- 設計書 §9.4 の Wake Lock の選択、終了前の注意、スリープやバックグラウンドからの復帰時の状態確認（S12、`REQ-22`）：Wake Lock は利用者が選び、Run が実行中か待っている間、ページが見えている間だけ持つ。Run が実行中か待っている間は離れる前に注意する。Run があるまま 1 秒以上隠れたページが戻ると、隠れていた時間、Run ごとの前と今の状態、Data worker の応答を示す。タブを閉じても検索が続くとは言わない。

### 5.6 保存層と寿命

- `BlockStore` は、ブロック単位で追記して読むだけの契約にする。実装は OPFS（Data worker の同期アクセスハンドル）と、OPFS が使えない環境向けの Memory の 2 つ。両方に同じ契約試験をかける。
- 配置は `tmp/<session-token>/runs/<run-token>/…`。ファイル名に入力名を入れない。
- 作業ごとに `navigator.locks` のロックを持ち続ける。起動時には、ロックを誰も持っていない `tmp/*` だけを消す。停止中や旧版のタブもロックを持ち続けるので、誤って消されない。
- 容量不足（`QuotaExceededError`）が起きたら、その実行を理由付きで失敗にし、過去の結果は残す。使用量は、操作を妨げない表示で示す。
- 次回の起動時に自動で復元することはしない（設計書 T01、`REQ-20`）。
- S10 の判断（`docs/evidence/losat_web_w2/README.md`）：確定はファイルを動かさず、Data worker の登録簿が持つ。所有の印は作業ごとの Web Lock だけで、ロックを持てないタブは OPFS を使わず Memory にする（他のタブがそのデータを放置と見なさないため）。OPFS が開けないとき（放置されたデータでクォータが一杯のときなど）は、先に回収してからもう一度開く。

### 5.7 結果画面・抽出・候補

- 画面の流れ：Run → Query → Subject 一覧 → HSP →（ドットプロット / 詳細）。選択の中心は HSP の ID で、行番号やヘッダーは使わない。
- Subject 一覧の列は、outfmt 6 の値と、outfmt 0 の見出し部の値から始める。Total score や Query cover のように NCBI Web にしか無い集約列は、列定義表（設計書 §11.2）で採用が決まったものだけを、NCBI の `align_format` からエンジンへ移植する（S13 で判断し、必要ならセッションを挿入する）。TS で独自に計算しない。
- ドットプロットは Canvas で描く。軸の単位（nt / aa）と frame・向きを示す。未採用の全対全計算を裏で行うことはしない。
- 抽出の座標は `domain` の 1 か所で計算する。既定は元のレコードの向き。flank は左右を別々に指定でき、端で切り詰めたときは要求した範囲と実際の範囲を両方記録する。複数の HSP は、別々の配列にするか、間を含む 1 区間にするかを選べる。ギャップ付きアラインメントは `PairwiseHit` の整列文字列から作り、原配列の抽出とは区別する。
- 候補トレイは application 層の状態で、元の run・HSP への参照とメモを持つ。

### 5.8 出力・セッション・再現性

- 互換出力（outfmt 0/6/7 のうち、その program が対応するもの）：保存したテキストを、1 バイトも変えずに書き出す。範囲は実行の結果全体だけで、フィルターや選択を反映した互換出力は作らない（S01 で書き出しを実装）。
- CSV・JSON・レポート：HSP レコードから作り、ViewState のフィルターや選択を反映できる。アプリ独自の形式と明記し、NCBI と一致するというラベルは付けない。レポートは外部スクリプトも計測も含まない静的 HTML にし、入力はすべてエスケープする。
- 設定ファイル：検索条件（program、argv のうち入力名以外、スレッドの設定）だけを書き出し、読み込むとフォームに入る（`REQ-15`）。
- セッションファイル：manifest（schema の版、app / engine の版、RunSnapshot、レコード表、入力の SHA-256）と、テキスト・レコードのブロックを長さ付きで並べた単純なコンテナにし、`CompressionStream` で gzip する（依存ライブラリ無し）。読み込んでも検索は始めない。flank や全長の抽出には元の FASTA を選び直してもらい、SHA-256 が一致したときだけつなぎ直す（`REQ-23`）。候補とメモを含めるかは、S15 の前に保守者が決める。
- CLI コマンドと NCBI との比較用コマンドは、RunSnapshot の argv から、出力形式ごとに作る。TBLASTX や TBLASTN で既定以外の subject 遺伝暗号を使う場合は、承認済みの例外（`AGENTS.md`）であることを注記する。

### 5.9 オフラインと明示更新

- 資産は `/r/<release-id>/…` に版ごとに置く。`index.html` は小さな起動用ページにし、固定した release ID の資産を読み込む。release ID は、研究データとは別のアプリ設定として `localStorage` に置き、読み書きは try/catch で囲む。
- Service Worker は、固定した release の資産（2 つの wasm、Worker、例題を含む）を取得してハッシュを確かめ、それから「Offline ready」を出す。
- 更新の流れ：新しい release を知らせる → 利用者が選ぶ → 資産を取得して確かめる → 実行中・未保存の作業への影響を案内する → 固定 ID を切り替えて再読み込みする。取得が不完全なら今の release のままにする。`skipWaiting` を無条件には使わない。旧版の資産は、それを使っているタブがある間は消さない。

### 5.10 GA4 と、研究データを送らない境界

- S16 で次の 2 案を比べて決める。(a) 設計書の案：別 origin の計測文書を iframe で埋め込む。COEP の下で埋め込めるか、gtag を読み込めるかを、対象の全ブラウザで確かめる。(b) 代替案：アプリは自分の Cloudflare Worker へ固定のイベントだけを送り、Worker がそれを GA4 の Measurement Protocol へ転送する。どの文書も Google のスクリプトを読み込まない。
- どちらの案でも、初期のイベントはアプリの訪問と release ID だけ。同意を得る前は何も読み込まず、何も送らない。拡張計測、自由形式の dataLayer、自動のエラー送信は使わない。
- CSP は `web/app/public/_headers` にある（S01）。(a) を採るときは `frame-src` に計測用の origin を足す。
- 研究データを送らないことの試験（V-PRIV）：固有の試験文字列を、配列・ヘッダー・ファイル名・メモ・不正な入力に入れる。正常・取消・失敗・書き出し・再読込・更新・計測の各経路で、Playwright がすべての要求を捕まえ、許可した origin とだけ通信していること、どの要求にも試験文字列が含まれないことを確かめる。

---

## 6. 検証戦略

### 6.1 3 段階の同一性

```text
NCBI BLAST+（oracle） ─[既存の認証]─► ネイティブ LOSAT の凍結バイト ─[新設のゲート]─► ブラウザの LOSAT Web
```

- ブラウザのゲートでは、NCBI の正解を新しく作らない。同じ argv と入力に対して、ブラウザの出力が期待値と一致することを確かめる。
- 期待値の出どころは、升目ごとに TD-5 のとおりに決め、`docs/web/verification_cells.tsv`（S02 で作る）に一覧する。認証済みの凍結バイトと一致した升目だけが、認証を引き継ぐ。
- 検証バッジの判定表は、この一覧と既存の認証記録から生成し、手では書かない。認証済みプロファイルの外の設定でも実行はできるが、その場合の表示は「Engine-supported, outside certified profile」にする。

### 6.2 試験の種類と実行場所

| ID | 内容 | 実行場所 |
|---|---|---|
| V-NAT | `run_local` の各形式 = 期待値。既存の CLI の回帰ゲートも含む | `cargo test` と既存の比較スクリプト |
| V-ABI | adapter の serial / threaded reactor × 升目 × スレッド 1/2/4 が期待値と一致する。Subject を保持した連続実行と、TBLASTX の `wasm-threads` 専用の出力箇所を通る場合を含む | Node（`LOSAT/tests/wasi_thread_host.js` と同じ方式、CI） |
| V-BR | 実アプリの EngineGateway を通した V-ABI の一部。取消 → 回復、instance の作り直しも含む | Playwright：Chromium・Firefox・WebKit（CI では一部、公開前は全部） |
| V-APP | domain と application の単体試験、層の規則の試験、BlockStore の契約試験、主な操作の E2E | Vitest と Playwright（CI） |
| V-PRIV | §5.10 の試験 | Playwright（CI） |
| V-OFF | オフラインでの起動、明示更新、新旧の資産が混ざらないこと | Playwright（CI） |
| V-MOB | iOS Safari と Android Chrome の実機 | 手動で記録する |
| V-PERF | 1 回の暖機と 3 回の計測（`AGENTS.md` の手順）。エンジンのまとめ直し（S02〜S04、SX）では、変更した program ごとに 1 つの fixture で、ネイティブと command-WASI の中央値が基準の中央値の +5% 以内であること（超えたら原因を調べて記録し、保守者に示す）。そのほかは cold / warm、serial / threaded、規模別に測る | ローカル。公開する上限値と DW-8 の判断の根拠 |

### 6.3 証拠の置き方

`docs/evidence/losat_web_<stage>/` に、ゲート記録の `README.md`、`evidence.sha256`、再現スクリプト、実行ごとに変更しない `run-<UTC>/` を置く。エンジン側の主張は、独立監査の観点（セッション README の「レビュー」）で確かめる。アプリ側は、段階ごとに E2E の記録と画面の記録を残す。回帰の基準は凍結バイトを優先し、無ければ、S02 の最初にこの worktree で全 program について 1 回だけ取った出力を使う（後の変更で基準が汚れないようにする）。そのために別の worktree を作らない。

---

## 7. 段階計画

セッションは、エンジン側とアプリ側のそれぞれの中で表の順に 1 つずつ実行し、2 本までを並行する（DW-7）。アプリ側の S10 は、S08 を入口の条件とする S09 より先に行う（S10 の本物の reactor への接続は S09 が行う）。指示書は `docs/losat_web_gui_sessions/` にある。この表の完了条件が正本である（§0.3 の 4）。

| セッション | 段階 | 内容 | 完了条件（証拠） |
|---|---|---|---|
| S01 | **W0** 契約と骨格 | PD、`web/AGENTS.md`、ルートの `AGENTS.md` への範囲の追記、要求トレース表、ABI v2 の下書き、`web/app` の骨格（層、FakeEngine、キューと状態機械、書き出し）、CI `web.yml`、`_headers`。TBLASTX の v1 で outfmt を黙って置き換える不具合の修正 | `npm run check` と `npm run e2e`（FakeEngine の一連の操作、`crossOriginIsolated`）が通る。TBLASTX の修正の単体試験が wasm32-wasip1 で通り、作り直した serial reactor で outfmt 0/7 が拒否される。ゲート記録と `evidence.sha256` がある。**完了（2026-09-29）** |
| S02 | **E1a** 核の入口：共通部と BLASTP | 最初に全 program の基準を取る。`run_local`、`ReportOutputs`（形式ごとの writer、`diagnostics`、`hits`、observer）、形式ごとのオプション解決と食い違いの検出（TD-4）を作り、BLASTP の CLI・v1 をそこに通す。`docs/web/verification_cells.tsv` を作る | 変更したコードを使う全 program の既存ゲートと Gate A の該当ハッシュが変わらない。v1 の serial / threaded reactor の検査が通る。V-NAT。V-PERF の非退行。独立監査。**完了（2026-09-29）** |
| S03 | **E1b** 核の入口：TBLASTN | 同じことを TBLASTN に行う | TLOSAN 計画の Stage G のゲートが変わらない。v1 の reactor の検査。V-NAT。V-PERF の非退行。独立監査。**完了（2026-09-29）** |
| S04 | **E1c** 核の入口：BLASTN と TBLASTX | 同じことを BLASTN と TBLASTX に行う（出力は 6/7 と 6 のまま）。TBLASTX の出力箇所を 1 つにまとめる | 既存のゲートと Gate A のハッシュが変わらない。v1 の reactor の検査。V-NAT。V-PERF の非退行。独立監査。**完了（2026-09-29）** |
| S05 | **E1d** アダプタと ABI v2 | `web/adapter` の crate、ABI v2 の確定、2 つの reactor、ビルドの同一性の検査（TD-6）、V-ABI の Node の仕組み、`scan` の性質試験（TD-8） | V-ABI（BLASTP・TBLASTN・BLASTN・TBLASTX の、その時点で対応する全形式 × スレッド 1/2/4）が期待値と一致。同一性の検査が通る。**完了（2026-09-29）** |
| S06 | **E2a-1** BLASTN outfmt 0：権威と fixture | NCBI の呼出し経路（`blast_format.cpp` → `align_format`）の記録。比較する fixture と NCBI の出力の固定 | 経路の対応表、固定した fixture と SHA-256。**完了（2026-09-29）** |
| S07 | **E2a-2** BLASTN outfmt 0：実装とゲート | S06 で見つかった BLASTN の panic の修正、座標の桁数の共有の関数（TD-9）、移植、`PairwiseHit` の作成、`run_local` と観測者への接続 | 固定した fixture（S06 の manifest の BLASTN・BLASTP・TBLASTN の 34 件、stderr を記録したものは stderr も）で NCBI とバイト一致。既存の 6/7 に退行なし。TD-9 で変わる BLASTP・TBLASTN の凍結出力は、すべて NCBI の出力と一致する。BLASTN の全升目（0/6/7 × スレッド 1/2/4）の V-ABI。独立監査。**完了（2026-09-29）** |
| S07+ | **E2c** BLASTN の得点のオプションと入力の読み方 | NCBI のオプションの検査（同じ拒否と文言）の移植。既定以外の得点で NCBI と違う原因の調査と修正。直せない組合せの明示的な拒否（TD-10）。S07 で見つかった、query の組成に依存する ungapped Karlin block と、FASTA の読み方・警告の時点の差 | S06 の `scoring_sweep.py` の全組合せが、NCBI と同じ拒否、outfmt 6 のバイト一致、明示的な拒否のどれかになる。直した組合せの fixture で NCBI とバイト一致。既存の BLASTN のゲートと S07 の fixture に退行なし。V-PERF の非退行。独立監査 |
| S07++ | **E2f** BLASTN の query の batch | NCBI の query の batch（`CBatchSizeMixer`、batch ごとの query の塊・lookup table・対角線の表）と query の分割（`CQuerySplitter`、塊ごとの予備の段階と HSP の合わせ方）の移植と、batch に依存する S07+ の拒否の置き換え（TD-14、DW-12） | 複数の query の sweep が NCBI とバイト一致。S07+ の検査、既存の BLASTN のゲート、S07 の fixture に退行なし。V-PERF の記録。BLASTN の V-ABI。独立監査 |
| S07+++ | **E2g** BLASTN の経路の棚卸しと一括の移植 | アプリが出す BLASTN のオプションの範囲で、NCBI の実行経路（`blastn_app` から検索・traceback・整形まで）に現れる関数の棚卸しと、未移植・差のある移植の一括の transpile（DW-12） | 棚卸しの表（関数ごとの NCBI のファイル・行、LOSAT の対応、状態）が経路を網羅し、未移植と差のある移植が残らない（残すものは明示的な拒否か承認済みの例外）。S07+ と S07++ の全検査、既存の BLASTN のゲート、fixture に退行なし。V-PERF の非退行。棚卸しの表を基準にした独立監査 |
| S08 | **E2b** TBLASTX outfmt 0/7 | 権威の記録、fixture の固定、移植、`PairwiseHit` の作成 | 固定した fixture で NCBI とバイト一致（承認済みの遺伝暗号の例外を除く）。既存の 6 に退行なし。TBLASTX の全升目の V-ABI。独立監査 |
| SD | **E2i** BLASTN の dc-megablast と blastn-short | NCBI の discontiguous megablast の経路（task の既定値、`-template_type`・`-template_length`、`s_DiscWordOptionsValidate`、discontiguous の template の lookup と走査、拡張）と blastn-short の task の既定値の棚卸しと一括の移植（DW-18） | 固定した fixture で NCBI とバイト一致（outfmt 0/6/7、スレッド 1/2/4、dc-megablast は template の種類・長さ・word size の全組、blastn-short は短い query と既定値の上書き）。NCBI が受け付ける両 task の組合せの sweep が、同じ拒否、バイト一致、明示的な拒否のどれか。BLASTN の既存のゲートと fixture に退行なし。V-PERF の非退行（megablast と blastn）。両 task の升目の V-ABI。独立監査 |
| S08+ | **E2e** BLASTP・TBLASTN・TBLASTX の既定以外のオプション | NCBI のオプションの検査（同じ拒否と文言）の移植。既定以外の値で NCBI と違う原因の調査と修正。直せない値の明示的な拒否（TD-13） | 各 program の sweep の全組合せが、NCBI と同じ拒否、outfmt 0/6/7 のバイト一致、明示的な拒否のどれかになる。直した組合せの fixture で NCBI とバイト一致。各 program の既存のゲートと S07・S08 の fixture に退行なし。V-PERF の非退行。変えた program の V-ABI。独立監査。**完了（2026-10-05、[E2e のゲート記録](evidence/losat_web_e2e/README.md)）** |
| S09 | **W1** ブラウザでの実行基盤 | Engine worker、WASI shim、ThreadHost、機能の確認、serial への切り替え、取消、instance の作り直し、R1。DW-8 のための前処理の割合の実測。S10 の port（`RecordScanner`、`EngineInput`、run の出力の経路）を本物の reactor と Engine worker につなぐ | V-BR（BLASTX を除く 4 program、3 ブラウザ、n=1/2/4）。取消の後の実行が成功する。メモリの推移と前処理の割合を記録する。S10 の契約試験（`record-scanner`・`engine-input`・`run-output`）が本物の reactor と Engine worker で通る。**完了（2026-10-03、[W1 のゲート記録](evidence/losat_web_w1/README.md)）** |
| S09+ | **R2**（条件付き） | DW-8 の条件を満たした program ごとに、前処理キャッシュを移植する。S09 の実測で TBLASTN が条件を満たした（subject だけの前処理が命令数の 40.9%。その 76% は codon ごとの `GeneticCode::get`）。進め方は DW-22（翻訳の表引きを先にし、20% 未満になればキャッシュは入れない） | 連続実行の出力が CLI と一致。独立監査。README の表に行を足して実施する |
| S10 | **W2** データ層 | Data worker、SourceStore、索引、DatasetRevision、BlockStore（OPFS / Memory）、登録簿、Web Locks、回収、容量不足 | 両方の実装が契約試験を通る。強制終了の後の回収、2 つのタブの保護、容量不足の試験が通る。**完了（2026-09-30、[ゲート記録](evidence/losat_web_w2/README.md)）** |
| S11 | **E2d** `-query_loc` / `-subject_loc` | BLASTN・BLASTP・TBLASTN・TBLASTX への移植（BLASTX は SX） | 固定した fixture で NCBI とバイト一致。独立監査。**完了（2026-10-06、[E2d のゲート記録](evidence/losat_web_e2d/README.md)）** |
| S12 | **W3** 検索画面 | program のタブ、入力、レコード一覧と除外、領域の指定（DW-9）、`describe` から作るフォーム、キュー、スレッド、段階と経過時間、診断、モバイル、Wake Lock など | 研究作業と境界条件の E2E。**完了（2026-10-06、[W3 のゲート記録](evidence/losat_web_w3/README.md)）** |
| S13 | **W4** 結果画面 | Subject 一覧、HSP、詳細、ドットプロット、Run の詳細、検証バッジ、フィルター、列定義表 | 対応済みの全 program の E2E。HSP と行・節との対応の試験 |
| S14 | **W5** 抽出と候補 | flank、ヒット区間、全長、複数 HSP、multi-FASTA、ギャップ付きアラインメント、候補トレイ、メモ | 領域指定、端、逆向き、翻訳の各場合で原配列と一致 |
| S15 | **W6** 出力と再現性 | CSV / JSON、レポート、設定ファイル、セッション、CLI コマンド、NCBI との比較用コマンド、使った入力 FASTA の書き出し | セッションを再読込しても再計算しない。元配列が無いとできない操作を明示する。勝手につなぎ直さない |
| S16 | **W7** 配信の仕上げ | オフライン、明示更新、GA4 の方式の決定と実装、V-PRIV、例題 5 つ、最短のチュートリアル、任意のブラウザ自己試験、Cloudflare のプレビュー（デプロイは保守者） | V-OFF と V-PRIV が通る。プレビューで `crossOriginIsolated === true` |
| SF | **E2h** FASTA の読み方（`CFastaReader` の移植） | NCBI の読み込みの経路（`CBlastFastaInputSource`・`CBlastInputReader`・`CFastaReader`・`CStreamLineReader` と、読んだ ID と題を使う報告）の棚卸しと一括の移植で、TD-12 と、E2b・E2e で TBLASTX・TBLASTN・BLASTP に入れた同じ種類の拒否をなくす（BLASTX は SX）。アダプタの `scan` に NCBI の読み込み器の種類を足し、`register` をそれで読む（TD-8、DW-13） | 棚卸しの表が経路を網羅し、未移植と差のある移植が残らない（残すものは明示的な拒否か承認済みの例外）。入力の読み方の fixture（abi_v2 の `register` の拒否の一覧と E2c の §E・§G・§M・§N の各場合、outfmt 0/6/7、スレッド 1/2/4、stdout・stderr・終了コード）で NCBI とバイト一致。入力の sweep が、一致、同じ誤り、明示的な拒否のどれか。新しい `scan` の種類がエンジンの読み込み器と性質試験で一致。各 program の既存のゲート、Gate A、TLOSAN の Stage G、fixture、sweep、v1 の検査に退行なし。変えた program の V-ABI。V-PERF の非退行。独立監査 |
| SX | **BLASTX の統合**（条件付き） | 入口の条件：LOSATX 計画の v0.2.0 の認証が `main` に入っていること。`main` を取り込み、BLASTX を `run_local` に通し（formatter の中の警告の整理、形式ごとの表示数の解決を含む）、ABI v2・V-ABI・V-BR に加え、`-query_loc` / `-subject_loc` を移植する（DW-11） | LOSATX 計画の比較ゲートが変わらない（移植した範囲指定を除く）。範囲指定の fixture で NCBI とバイト一致。BLASTX の全升目の V-ABI と V-BR。V-PERF の非退行。独立監査。条件が満たされた後の最初のセッションの区切りで実施し、S17 の前に必ず終える |
| S17 | **G** 公開判定 | 受入表、対応環境の表、実機の記録、実測した上限値の公開、既知の例外の一覧、独立レビュー。最終の commit で 5 program の認証のゲートを再実行する | 要求トレース表の initial がすべて満たされている（設計書 §2.2「初期必須を省いたものを正式公開しない」） |

1 つのセッションで段階を終えられない場合は、同じ段階の続きを次のセッションとし（例：S07b）、README の表に行を足す。

---

## 8. リスク

| リスク | 影響 | 対策 |
|---|---|---|
| 核の入口へのまとめ直しで、認証済みの出力が変わる | 既存の認証 | 出力を変えない変更だけを S02〜S04 で行い、最初に取った基準と Gate A のハッシュ、v1 の reactor の検査、性能の非退行を完了条件にする |
| 共有の formatter の変更が、BLASTX の進行中の作業とぶつかる | 衝突、LOSATX の認証の前提が崩れる | BLASTX が使う関数の本体と引数は SX まで変えない（§4.2）。`main` の取り込みのたびに、該当するゲートを再実行する |
| BLASTN の outfmt 0 の移植が大きい（NCBI の `align_format`） | 公開の時期 | 権威と fixture の固定（S06）と移植（S07）を分ける |
| BLASTX の認証が HOLD のまま | 5 program での公開 | LOSATX 計画で完了させる。SX と S17 の前提にする |
| COEP の下で GA4 の iframe が動かない | GA4 の方式 | S16 で (b) 案と比べる |
| Chromium での共有メモリの境界の問題（Node には `LOSAT/tests/wasi_shared_memory.js` という回避策がある） | threaded の経路 | S09 で判定した：Chromium 149 で起きる（60 回に 1〜2 回）。Node と同じ guard をビルドの時に threaded の module にかける（`build/reactors.ts`） |
| 304（条件付きの要求への応答）に COOP / COEP / CORP が付かない | WebKit が、作り直す Engine worker の読み込みを拒否する（取消の後の検索が失敗する） | 開発とプレビューのサーバーは、すべての応答に付ける（S12）。本番の Cloudflare は S16 で確かめる |
| reactor を繰り返し呼ぶとメモリが増える | 長時間の作業 | instance を作り直す方針。メモリの推移を RunRecord に記録する |
| メモリの上限（1 GiB、端末、モバイル） | 扱える規模 | 実測して上限を公開する。上限を超えたら理由を示して失敗にする。拡大は §2.3 |
| WebKit での OPFS の同期ハンドルと SAB の対応の差 | 機能の経路 | 機能を実際に確かめて経路を選ぶ。Memory / serial に切り替える |
| アダプタのビルドが認証済みの command-WASI とずれる | 性能と来歴 | TD-6 の同一性の検査 |
| `scan` とエンジンの解析器の食い違い | 抽出 | TD-8 の性質試験と、`register` の時点での照合 |
| 設計書の未決事項による範囲の膨張 | 期間 | 要求トレース表で initial とされた項目だけを実装する |

---

## 9. ワークフロー

### 9.1 ブランチと作業場所

規則はセッション README の規則 1〜2 にある（ブランチ・worktree・上流・`main` の取り込み方）。Rust のビルドの出力先は README の規則 4 にある。

### 9.2 セッション

共通の規則とレビューの観点はセッション README にある。毎セッションでコミットとプッシュを行い、次のセッションの指示書を、実測や結果に合わせて更新してから終える（§7 の完了条件は緩めない）。

### 9.3 レビューと判断

- 計画の改訂：大きな改訂の後に、計画レビューの観点で確認する。これまでのレビューの結果と対応は、[W0 のゲート記録](evidence/losat_web_w0/README.md)の「計画レビュー」の節にある。
- エンジン側（S02〜S08、S11、R2、SX）：独立監査。
- 画面（S12〜S16）：画面レビュー。
- 保守者が決めるもの：公開、デプロイ、GA4 の設定、§10 の項目。

---

## 10. 未決事項

設計書 §18 から持ち越すもの：表示フィルターと列の定義表、元の核酸文字列を並べて表示する時期、SVG / PNG、トラック、all-vs-all、ダークモード、NCBI との比較用コマンドのシェル形式、外部リンク、GA4 の詳細（同意・保存期間）、オフラインの範囲、モバイルの規模、旧版の提供、担当と期日。

この計画で新しく生じたもの：

| 項目 | 決める時点 |
|---|---|
| program の表示名（BLASTN 系か、LOSATN 系か、併記か） | 決定（2026-10-06、DW-21：BLASTN 系） |
| GA4 の方式（(a) iframe か、(b) Measurement Protocol への転送か） | S16 |
| §2.3 の閾値（取消の後の再準備、描画件数） | 取消の後の再準備は決定（DW-22：1 秒、S17 で対応規模の最大の入力で測り直す）。描画件数は S13 の実測の後 |
| instance を作り直す条件 | 決定（S09：検索の後の linear memory が 512 MiB 以上。§5.5） |
| セッションファイルに候補とメモを含めるか | S15 の前 |
| Cloudflare に残す旧版の数、ドメイン名 | S16 の前 |
| BLASTX の範囲の拡大（DW-11）を LOSATX 計画の範囲の記録に書くこと | SX の前（保守者） |
| v1 ABI の廃止 | gbdraw が v2 へ移る時点（TD-1） |
| BLASTN の FASTA の読み方を NCBI の `CFastaReader` に合わせる（TD-12 の拒否をなくす。アダプタの索引の走査の解析器の種類を足す）。NCBI の後の query batch の大きさの再現は S07++（TD-14）に移した | S17 の前の専用のセッション SF（2026-10-02 に保守者が決定、DW-13）。[指示書](losat_web_gui_sessions/session_sf_e2h_blastn_fasta_reader.md)は 2026-10-06 に作った（範囲は BLASTN・TBLASTX・TBLASTN・BLASTP の全入力。BLASTX は SX。DW-23） |
| `-num_threads` の NCBI の警告（CPU の数を超えると「Number of threads was reduced to N …」、`-subject` があると「'num_threads' is currently ignored when 'subject' is specified.」、`blast_args.cpp:3203-3236`）。LOSAT はどの program も `-subject` でスレッドを使い、警告を出さない（出力は同じ、stderr だけが違う。以前から。S07+ の第 4 回の監査、`docs/evidence/losat_web_e2c/AUTHORITY.md` §I）。承認済みの例外にするか、警告を出すか | S17 の前（保守者と相談） |
| 公開する版の名前 | S17 |
