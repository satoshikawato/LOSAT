# Session S16 — W7：配信の仕上げ

## INSTRUCTION PROMPT

LOSAT Web の段階 W7 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S16 の行である。設計は計画 §5.9〜§5.10、設計書 §13〜§14 にある。Cloudflare へのデプロイと GA4 の設定は保守者が行う（README の規則 7）。

1. release ごとの資産：資産を `/r/<release-id>/…` に置き、小さな起動用の `index.html` が、`localStorage` に固定した release ID の資産を読み込む（読み書きは try/catch で囲む）。資産の一覧とハッシュを ReleaseManifest として生成する。
2. Service Worker：固定した release の資産（2 つの wasm、Worker、例題を含む）を取得してハッシュを確かめ、それから「Offline ready」を出す。研究データは Service Worker のキャッシュに入れない。
3. 明示更新：新しい release を知らせる → 利用者が選ぶ → 資産を取得して確かめる → 実行中・未保存の作業への影響を案内する → 固定 ID を切り替えて再読み込みする。取得が不完全なら今の release のままにする。`skipWaiting` を無条件には使わない。旧版の資産は、それを使っているタブがある間は消さない。
4. GA4：計画 §5.10 の (a)（別 origin の計測文書を iframe で埋め込む）と (b)（自分の Cloudflare Worker から Measurement Protocol へ転送する）を、COEP の下で対象の全ブラウザで試し、結果とともに保守者に推奨を示す。決まった方式を実装する。初期のイベントはアプリの訪問と release ID だけで、同意の前は何も読み込まず、何も送らない。計測が塞がれていても、未同意でも、オフラインでも、検索・結果・抽出・保存が動くこと。
5. V-PRIV（計画 §5.10）と V-OFF（計画 §6.2）を Playwright で作り、CI に入れる。V-PRIV には、配信する Wasm と JavaScript にビルドした機械のパス（checkout のパス、`CARGO_HOME`、ユーザー名）が含まれないことの検査を含める（計画 TD-11）。
6. 導入支援：「例を試す」、5 program の小さな例題、最短のチュートリアル、任意のブラウザ自己試験（期待値の SHA-256 と比べる、V-BR の小さな版）。
7. Cloudflare のプレビューの設定（`web/app/public/_headers` をそのまま使う）を用意し、保守者がデプロイした後に、プレビューで `crossOriginIsolated === true` と CSP を確かめる。

完了条件は計画 §7 の S16 の行による。画面の記録を画面レビューに見せ、`docs/evidence/losat_web_w7/README.md` に記録する。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S17 — 公開判定](session_s17_g_release_decision.md)。
