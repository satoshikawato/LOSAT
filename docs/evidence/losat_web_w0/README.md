# LOSAT Web W0（Session S01）ゲート記録

- 段階：W0 契約と骨格（[総合計画書](../../losat_web_gui_plan.md) §7 の S01、[指示書](../../losat_web_gui_sessions/session_s01_w0_contract_skeleton.md)）
- ブランチ：`feature/losat-web-gui`（基点 `origin/main` の `0627c88f5`）
- 実行記録：[`run-20260928T153516Z/`](run-20260928T153516Z/)（作成後は書き換えない）
- 判定：**完了条件を満たした**。`PD-LOSAT-WEB-APP-BOUNDARY` は 2026-09-29 に保守者が承認した（計画 DW-12）

## コミット

| コミット | 内容 |
|---|---|
| `d9a114aa9` | エンジン：TBLASTX の Web ABI v1 で、未対応の outfmt を拒否する（`LOSAT/src/web_api.rs`） |
| `35e34eee0` | アプリ：`web/app` の骨格、`web/tools/wasi-test-runner.mjs`、`.github/workflows/web.yml` |
| `f64106f5f` | アプリ：2 回目の計画レビューへの対応（表示する CLI コマンドに `-outfmt N` を付ける、出力形式の型を `domain` へ移す、`ProgramDescription.formats`、層の規則の試験、CI の対象パスの拡大）。このコミットで `npm run check`（単体試験 4 ファイル 37 件）と E2E 2 件が通ることを手元で確かめた |
| このゲート記録を含むコミット | 文書：総合計画書、セッション指示書、設計書、要求トレース表、ABI v2 の下書き、PD（Accepted）、`web/AGENTS.md`、ルートの `AGENTS.md` の追記、このゲート記録 |

## 完了条件と結果

| 完了条件（計画 §7 の S01） | 結果 | 証拠 |
|---|---|---|
| `npm run check`（typecheck、lint、単体試験、本番ビルド）が通る | 通過。単体試験 3 ファイル 19 件 | `run-20260928T153516Z/npm-check.log` |
| `npm run e2e`：FakeEngine で「貼り付け → キュー → 実行 → 結果 → 書き出し」、`vite preview` で `crossOriginIsolated === true` | 通過。2 件（Chromium 149、Playwright 1.61.1） | `run-20260928T153516Z/npm-e2e.log` |
| TBLASTX の修正の単体試験が wasm32-wasip1 で通る | 通過。`web_api::tests` の 4 件（新しい試験 `tblastx_web_args_reject_unimplemented_outfmt` を含む） | `run-20260928T153516Z/wasm32-web-api-tests.log` |
| 作り直した serial reactor で outfmt 0/7 が拒否される | 通過。outfmt 6 は status 0（1,937 バイト）、`0`・`7`・`-outfmt 0`・`6 qseqid` は status −1（`unsupported TBLASTX outfmt: only 6 without custom fields is implemented`） | `run-20260928T153516Z/serial-reactor-check.{txt,json}`、スクリプト `check_tblastx_outfmt_reactor.js` |
| 層の依存の向きを lint が検出する | 確認済み。その後、規則そのものを試験する `web/app/tests/unit/layers.test.ts`（禁止 12 件、許可 6 件）を足した | 2 回目のレビューへの対応のコミットで、`npm run check` の単体試験は 4 ファイル 37 件になった |

環境は `run-20260928T153516Z/environment.txt` にある（Node 26.8.2、npm 11.19.1、Rust/Cargo 1.92.0、WSL2 の Linux x86_64）。

GitHub Actions でも、`35e34eee0` の push で Web ワークフロー（run `36445003713`）が成功した。`engine-web-api`（wasm32-wasip1 の `web_api::tests`）と `app`（`npm run check` と E2E 2 件）の両方が通った。`actions/checkout@v4` と `actions/setup-node@v4` が Node 20 を対象にしているという非推奨の注記が出たが、既存の `ci.yml` と同じ版なので、そろえて上げるかは別に判断する。

## 計画レビュー

- 1 回目：`plan_critic` の観点による読み取り専用のレビューで、10 件の指摘と「2 つの判断が必要」という結論を受けた。判断は保守者が下し（計画の DW-8、DW-9）、技術的な指摘は計画の TD-1〜TD-6 と本文の改訂で反映した。
- 2 回目：改訂版に対するレビューは、1 回目の 10 件のうち大半が解決したとし、新たに 10 件（N1〜N10）を挙げて「BLASTX / TBLASTN の担当と順序の判断が必要」と結論した。
  - 判断は保守者が下した：TBLASTN は今、BLASTX は LOSATX の認証の後（DW-10）、BLASTX も領域の指定の対象にする（DW-11）、PD を承認する（DW-12）。
  - 技術的な指摘への対応：ゲート記録と `evidence.sha256` を置いた（N1）。形式ごとに別に検索する仕組みをやめ、検索オプションが食い違えば明示的なエラーにした（N3、TD-4）。S02 の最初に全 program の基準を 1 回だけ取り、共有のコードを使う全 program のゲート、v1 の reactor の検査、性能の非退行を S02〜S04 の完了条件にした（N4）。`scan` に解析器の種類を渡し、例外として位置付けた（N5、TD-8）。共有メモリの最大値を 1 GiB のまま据え置いた（N6、TD-7）。`describe` が対応する形式を返し、`run` はそれだけを出すことにし、S05・S07・S08 の升目を埋めた（N7）。CI の対象パスに `LOSAT/` を足した（N8）。表示する CLI コマンドに `-outfmt N` を付けた（N9）。細かな指摘（型の同期の規則、PD の入力の定義、ルートの `AGENTS.md` の段階名、ui の import の規則、ブランチの規則の重複、複数レコードの `-query_loc` の fixture）を直した（N10）。
  - 引用の行番号は、レビューで確かめられたもののほか、`LOSAT/src/cli.rs` の `-query_loc` の扱いを直した（187〜193 行は TBLASTN の一覧で、BLASTN・BLASTP・TBLASTX はこの引数を定義していない）。

## 保守者の判断待ち

1. Cloudflare のプレビューへのデプロイ（S16 で設定を用意する。W0 では `vite preview` で同じヘッダーを確かめた）。
2. BLASTX の範囲の拡大（DW-11）を LOSATX 計画の範囲の記録に書くこと（SX の前）。

## 残件と注意

- Firefox と WebKit の E2E は、ブラウザを入れる S09 から行う。W0 は Chromium だけ。
- v1 の変更の影響：gbdraw が TBLASTX に 6 以外の outfmt を渡すと、これまでは outfmt 6 が返っていたが、エラーになる。outfmt 6（v1 の既定値）の出力は変わらない。
- Rust のビルドの出力は worktree の外（`/home/kawato/.cache/losat-web-gui-target/`）に置いた。
