# Session S01 — W0：契約と骨格

## INSTRUCTION PROMPT

LOSAT Web の段階 W0 を実行する。先に [セッション README](README.md) の「全セッションに共通する規則」を読み、それに従う。完了条件の正本は、総合計画書 `docs/losat_web_gui_plan.md` の §7 の S01 の行である。

目的：エンジンがまだ無い状態でアプリ側を進められるように、境界の規則・契約・骨格を先に固定する（計画 §2.1 の D）。

0. ブランチの上流を設定する：`git push -u origin feature/losat-web-gui`。以後、`git pull --ff-only` はこのブランチからだけ取り込む。
1. `docs/product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md` を書く。範囲は `web/` と、Web のために置くエンジンの配管。内容は、アプリ機能を許す条件（エンジンに渡すのは CLI と同じ argv と原 FASTA バイトだけ、互換出力のバイトを変えたり作り直したりしない、BLAST の値を `web/` で計算しない、アプリ独自の出力であることを明記する、研究データを送らない）、許可する LOSAT 側の配管の一覧（計画 §4.1〜§4.5）、Subject を保持するときのバイト同一性の条件と、過去の「検索間 cache を持たない」制約との関係、v1 ABI の扱い（fail-fast の修正を除いて凍結）、非目標。状態は「Proposed」とし、承認は保守者に委ねる。
2. `web/AGENTS.md` を書き、ルートの `AGENTS.md` の「Repository Layout」に `web/` と規約の適用範囲を追記する。
3. `docs/web/requirements_trace.tsv` を作る。設計書 §2.1 の各行と、その他の確定要求（一時ファイル、非送信、Wake Lock、セッションの読み込み）を、状態（initial / withdrawn / deferred / out-of-scope / undecided）・段階・試験 ID に対応付ける。計画の決定で取り下げた要求は withdrawn にする。
4. ABI v2 の下書き `docs/web/abi_v2.md`（S05 が確定する）と、アプリが依存する port の型 `web/app/src/ports/engine.ts` を書く。
5. `web/app` に Vite + Vue 3 + TypeScript の骨格を作る。層は計画 §3.2 と §3.3 のとおりにし、依存の向きを ESLint の `no-restricted-imports` で検査する（違反を検出することも確かめる）。アプリ層（キュー、状態機械、RunSnapshot の固定、取消と確定の順序）は Vue に依存しない TS で書き、FakeEngine とメモリ上の DataGateway で動かす。FakeEngine の出力には「検索結果ではない」と書き、画面に警告の帯を出す。
6. 応答ヘッダーは `web/app/public/_headers` の 1 か所だけで定義し（COOP / COEP / CORP / CSP）、Vite の dev と preview、試験がそれを読む。
7. `.github/workflows/web.yml` を作る。`web/**` の変更でアプリの `npm run check` と E2E を実行し、`LOSAT/src/web_api.rs` の変更で、wasm32 のときだけコンパイルされる v1 の試験を wasm32-wasip1 で実行する（runner は `web/tools/wasi-test-runner.mjs`）。
8. エンジンの不具合を 1 件直す：TBLASTX の v1（`LOSAT/src/web_api.rs` の `parse_tblastx_args`）が outfmt を検証せず、常に outfmt 6 を出す。CLI と同じ検証関数 `blastinput::value_parsers::tblastx_outfmt` を通し、NCBI の参照コメントと wasm32 の単体試験を付ける。gbdraw への影響（6 以外を渡すとエラーになる）をコミットメッセージに書く。アプリの変更とは別のコミットにする。v1 の serial reactor を既存と同じ方法で作り直し（`cargo +1.92.0 build --release --locked --lib --target wasm32-wasip1 --no-default-features`）、outfmt 0・7・カスタム列が拒否され、6 が通ることを確かめる。

完了条件は計画 §7 の S01 の行による。結果は `docs/evidence/losat_web_w0/README.md` に、実行したコマンド・結果・コミット SHA・保守者の判断待ちの項目とともに記録する。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S02 — 核の入口：共通部と BLASTP](session_s02_e1a_core_entry_blastp.md)。S02 の入口の条件は PD の承認なので、承認されていなければ、そのことを最終回答の最初に書く。
