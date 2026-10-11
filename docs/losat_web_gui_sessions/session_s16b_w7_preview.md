# Session S16b — W7 の続き：プレビューでの確認

## INSTRUCTION PROMPT

LOSAT Web の段階 W7 の残りを実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S16 の行（V-OFF と V-PRIV が通る、プレビューで `crossOriginIsolated === true`）である。S16 は 3 つめの条件を確かめずに終わった。Cloudflare へのデプロイは保守者が行い（README の規則 7）、S16 では行われなかったからである。V-OFF と V-PRIV は S16 のゲート記録（[`docs/evidence/losat_web_w7/README.md`](../evidence/losat_web_w7/README.md)）にある。このセッションは、プレビューでの確認を最初に行い、W7 を閉じる。

始める前に、ゲート記録の「保守者の判断待ち」を読む。保守者がプレビューをデプロイしたか、プレビューの URL（`https://<branch>.<project>.pages.dev` か、Wrangler が出した 1 回ごとの URL）を渡したかを確かめる。まだなら、4. と 5.（アプリの変更を含む）を先に行い、1.〜3. が要る点で止まって、最終回答に「保守者のデプロイ待ち」と書く。デプロイ、GA4 の設定、Cloudflare へのログイン、外部サービスへの送信はこのセッションの中で行わない。確認の script をプレビューへ向けるのは、保守者が URL を渡し、実行を許したときだけにする（読み取りだけ。GA4 へのイベントは送らない：`LOSAT_WEB_PREVIEW_SEND_VISIT` を付けない）。許されなければ、保守者が同じ command を実行し、結果の JSON 5 つを渡す。

1. **保守者がプレビューをデプロイする**：[`docs/web/deploy_cloudflare.md`](../web/deploy_cloudflare.md) の手順（§2 の初回の設定、§3 の build、§4 のプレビューへのデプロイ）。release は、未コミットの変更の無い clean な木から作り、release ID が `-dirty` でも `-fake` でもないことを確かめる。サイトは `LOSAT_WEB_ANALYTICS_ENDPOINT=/api/visit` で build する（プレビューには GA4 の secret が無く、Function は何も転送しない）。保守者の判断（ゲート記録の「保守者の判断待ち」）：Pages のプロジェクト名と production branch、旧版をいくつ残すか（`npm run assemble` は渡された release をすべて残す）、GA4 の方式 (b) の確認、同意の文言。決まらないものは、推奨の案で進めて記録する（Owner-delegated、判断の記録に日付を付ける）。保守者の判断を要するもの（デプロイそのもの、GA4 の設定、同意の法的な確認、本番）は進めない。

2. **プレビューの確認を 3 ブラウザで行う**：`web/app` で次を実行する（`LOSAT_WEB_PREVIEW_HEADERS` はデプロイした `dist/_headers`。`LOSAT_WEB_WEBKIT_EXECUTABLE` は、この機械の WebKit の launcher `/home/kawato/.cache/losat-web-gui-target/app-s09-webkit-deps/pw_run.sh`）。
   ```bash
   cd web/app
   LOSAT_WEB_PREVIEW_URL=https://preview.<project>.pages.dev \
   LOSAT_WEB_PREVIEW_HEADERS=dist/_headers \
   LOSAT_WEB_PREVIEW_OUT=$HOME/losat-preview-check/<release-id> \
   LOSAT_WEB_WEBKIT_EXECUTABLE=/home/kawato/.cache/losat-web-gui-target/app-s09-webkit-deps/pw_run.sh \
   npx playwright test -c ../../docs/evidence/losat_web_w7/check_preview/playwright.config.ts
   ```
   手順の詳細は [`check_preview/README.md`](../evidence/losat_web_w7/check_preview/README.md)。確認の script は見出しの `offline-status` の `data-state` が「ready」になるのを待つ（`35970fa1`。ゲート 2 で、なくなった `offline-ready` の test ID を待って落ちた）。アプリが変わったら、本物のサイトへ向ける前に、ローカルの server への dry run（`run_gate.sh` の `all` の一段階と同じ）を通す。エンジン入りのサイトでないと、取消の後の検索は skip される（`-fake` のサイトに向けない）。記録の場所は S16 のゲートのように、実行ごとに `docs/evidence/losat_web_w7/preview-<UTC>/` へ `summary.json`・`http.json`・`chromium.json`・`firefox.json`・`webkit.json` を写す（URL・状態・応答の header だけで、研究データは無い。作成後は書き換えない）。ゲート記録に次を書く：
   - 3 ブラウザで `crossOriginIsolated === true`、`SharedArrayBuffer` があること、起動ページが最新の release を始めたこと、見出しの release ID が `/releases.json` の `latest` と同じこと。
   - 画面に出る CSP と、`/`・`/boot.js`・`/boot.css`・`/releases.json`・`/sw.js`・manifest・entry の script と style・Engine worker・2 つの engine module の **200、304、404** のそれぞれで、COOP / COEP / CORP / CSP が `public/_headers` の規則と等しいこと。**304 に COEP と CORP が無いのは失敗**（WebKit は、304 に COEP が無いと、作り直す Engine worker の読み込みを拒否する。[W3 のゲート記録](../evidence/losat_web_w3/README.md)の判断 12）。304 と 404 の他の header（`Cache-Control`、`X-Content-Type-Options`、`Referrer-Policy`）は `notes` に残る。
   - `/r/*` の 200 の `Cache-Control: public, max-age=31536000, immutable` と、サイトの file の `no-cache`。**`/r/` の下の存在しない file の 404 に immutable が付くか**（付くなら、取り上げる release のデプロイを、それを尋ねる前にする運用か、それを防ぐ手がいる。`deploy_cloudflare.md` §6）。
   - 配った bytes が manifest と shell の一覧と byte で一致すること（Cloudflare が minify や script の注入をしていないこと）。起動ページに他の origin の参照、Cloudflare の beacon（Web Analytics）、Rocket Loader、inline の script が無いこと。
   - WebKit（と他の 2 つ）で、検索を始め、実行中に取り消し、新しい検索が完了すること。Playwright の WebKit は WPE の MiniBrowser で、Safari ではない。実機の Safari は S17 の V-MOB（macOS・iOS）で見る。
   - 「Offline ready」が出ること。Chromium では、offline にした後の再読み込みで、`crossOriginIsolated === true` で始まること。Firefox と WebKit は Playwright の offline の emulation が確かでない（ローカルの V-OFF は server を止めて確かめた。プレビューの確認では確かめない）。
   - `/api/visit`：不正な body が 400、GET が 405、Function 自身の header（`_headers` は Function に適用されない）。

3. **Pages が `_headers` と違うことをしたら直す**：プレビューの確認が失敗したら、`http.json` の `failures` と `notes` がどの header がどの応答で違うかを言う。直す前に、違いが Pages の挙動か（Cloudflare の文書は `_headers` が 304 と 404 に届くかを言わない）、デプロイの手違い（`_headers` が上がっていない、`functions/` の `_routes.json`）かを分ける。Pages が 304 に COEP / CORP を付けない、または他の点でローカルの server と違うなら、選択肢と費用を記録して、保守者に示し、推奨を進めて確認を再実行する：
   - **Pages Functions の middleware で全応答に付ける**：`_headers` の規則から作った同じ表を使い、`functions/_middleware.ts` で付ける。費用：すべての要求が Function を通る（`_routes.json` の include を `/*` にする）ので、静的 file の要求が無料枠（1 日 100,000 要求）に数えられる。headers の定義が 2 か所になる（`web/AGENTS.md` の規則 4）ので、`build/headers.ts` を共有して 1 つの source にする。
   - **304 を作らせない**：影響する file に `Cache-Control: no-store` を付ける。費用：サイトの file（起動ページ、`boot.js`、`boot.css`、`releases.json`、`sw.js`）が毎回全部届く（どれも小さい）。`/r/*` の immutable は、ブラウザが条件付きの要求を出さないことに頼る設計で、出す経路（強制再読み込みなど）では 304 が残り得る（確かめていない）。
   - **Workers の静的 assets（Pages の外）に移す**：headers をコードで完全に制御できる。費用：配信の構成の変更（`deploy_cloudflare.md` を書き直す）、保守者の判断。
   直したら、ローカルの server の側（`build/site-server.ts`、`tests/unit/site-server.test.ts`）も実際の Pages に合わせる。確認を再実行し、前の記録は書き換えずに新しい `preview-<UTC>/` を足す。

4. **wrangler が `functions/api/visit.ts` を bundle できるか**：Function は `../../src/domain/analytics`（`.ts` の指定子）を import する。**保守者を待たずに最初に行う**：`web/app` で `npx wrangler@4.149.0 pages functions build --outdir <tmp>`（アカウントも送信も要らない）。失敗したら、Function を自己完結にする（共有の `parseVisitBody` を Function の側にも置き、同じ試験で両者が一致することを確かめる）。デプロイの出力（Pages がどの file を上げ、Function を compile したか）と生成された `_routes.json` を確かめる：`/api/*` だけを Function に通すべきで、`/*` なら 3. の費用の話になる。`analytics.spec.ts`・`visit-function.test.ts` が通ること。

5. **アプリの変更（最初に行う）と、S17 の前にしておく残件**：画面レビュー 2 回目（ゲート 3 の記録、合格）の Low 4 件が、この順で**最初のアプリの変更**になる（判断 19。見た目と文言の変更なので、プレビューの確認が済んだ後、ゲートの再実行の前に入れる）。画面レビュー 2 回目の表はゲート記録の「画面レビュー」にある。
   - **L1**：見出しに収まらない release ID の全文を、ポインタ無しで読めるようにする。Getting started に 1 行「This is LOSAT Web release `<id>`」。任意で、見出しの箱を約 32em にして既定の ID が全部入るようにする（状態 46・46c）。
   - **L2**：見出しの短い状態の文「Offline: 4.3 MB of 8.2 MB…」「Offline: waiting」を「Offline copy: …」などにする（同じ短さで。アプリがオフラインだと読めないように）。状態 46・59。
   - **L3**：他の tab の知らせ（`update-other-tab-loss`、`UpdateNotice.vue`）に、`UpdateImpact.vue` と同じ数（Run・候補・note）と、フォームが空で始まることを言う。状態 51。
   - **L4**：起動ページ 55 の `first-visit-offline`（`src/boot/plan.ts`）で、保存された release が出るときは「the first time it is opened」と言わない文にする（例「LOSAT Web cannot reach this site, and this browser has not chosen a release to start offline.」）。状態 55。
   - 参考 I2（Firefox の失敗の理由の末尾「.).」。括弧に入れる前に最後の句点を除く）と I3（電話の 43 の質問の枠に 6 px の上余白）は、同じ変更に入れてよい。

   その後、残りをゲート記録の「残件と注意」から選ぶ。推奨：
   - コードレビュー 1 回目 I1：shell の保存を atomic にする（新しい cache に書いて入れ替える）。または、`/`・`/boot.js`・`/boot.css` が release をまたいで互換であることを契約として `release_manifest.md` に書く。
   - コードレビュー 1 回目 I4：自己試験の「Stop」（case の Run ID を取り消す）。1 つの case が戻らないと利用者のキューが始まらない。
   - コードレビュー 2 回目 I1（V-OFF 3b を決定的にする：`cleaned` の返事を待つ）、I3（fixture の teardown でもう一度 `findings()`）、I2（準備が検索と重なる場面を 1 つ測る）。
   - コードレビュー 3 回目 I2（V-PRIV が request header の符号化を試験文字列の 3 つの形でしか探さない）。header の許可表を作るか、記録する。
   - コードレビュー 4 回目 I1〜I3（click の時点の再確認は別の tab の掃除との競合を閉じない、cookie を塞いだ WebKit の probe、保管場所の確認の余裕）。記録のままでよい。
   次は、S17 の前に直さなくてよい：画面レビューの参考のうち受け入れたもの（電話の見出しの増加）、WebKit の scroll anchoring（W6 の判断 17）、WebKit の書き出しの開始の間隔、Chromium の JSON の書き出しの最長のフレームの間隔 250 ms（W6 の A/B で退行ではない。判断 20）、隠れた tab の timer、`query.fa` / `subject.fa`、W5 コードレビュー L3、グラフの行 8 px、BLASTX（SX）、`abi_v2.md` の kind 0、E2j の申し送り、メモリ不足・保存の上限・`runInputs`・iOS Safari（S17）。

   アプリを変えたら、`run_gate.sh` を最後の木で `all` で走らせる（アプリの本体が変わるので、`after-review` では足りない。V-OFF・V-PRIV・TD-11・画面の記録・計測が新しい木の値になる）。実行は `docs/evidence/losat_web_w7/run-<UTC>/` に残し、ゲート記録の「ゲートの実行」に足す。変えた状態（46・46c・51・55・59、Getting started）の画面レビューを受ける。

最後に、ゲート記録 `docs/evidence/losat_web_w7/README.md` を更新し（プレビューの確認の結果、3 つめの条件、直したもの、判定）、保守者の判断待ちを更新し、[S17 の指示書](session_s17_g_release_decision.md)の「S16（W7）から引き継ぐこと」を実測と結果に合わせて直す。完了条件のうち 1 つでも満たさなければ、W7 は終わらず、S17 に進まない。

## S16（W7）から引き継ぐこと

[W7 のゲート記録](../evidence/losat_web_w7/README.md)の要点。上の 1.〜5. の条件はこの節で緩めない。

- **最後の木と件数**：最後のアプリの木は `35970fa1`（ゲート 3 の木。アプリの本体は修正 4 回目の `72ae82c9` から変わっていない。S16 の起点は `e45fbd09`、83 コミット、merge 15 を含む）。release は FakeEngine のビルドが `35970fa18885-fake`、エンジン入りのビルドが `35970fa18885-e2fb7128`。ゲートは 3 回走らせ（ゲート 1 `ba7fa464`：Firefox の V-OFF のケース 1 の保管場所の下限が偽で止まった、ゲート 2 `9cb4e61e`：プレビューの確認の dry run が、なくなった test ID を待って止まった、どちらも試験の側の誤りを直した）、ゲート 3（[`run-20261011T080950Z/`](../evidence/losat_web_w7/run-20261011T080950Z/)）はすべての段階が通った（約 48 分）：V-ABI quick（60 実行がネイティブの CLI と一致、凍結ハッシュ 16/16）、`npm run check` 905 件（65 件 skip）、reactor 付きの単体試験 969 件（1 件 skip）、NCBI との比較 42 件（33 一致、承認済みの例外 6、NCBI が拒む 3）、FakeEngine のビルドの E2E 297 件（Chromium 104・Firefox 97・WebKit 96）、エンジン入りのビルドの E2E 324 件（113・106・105。V-BR を含む）、結果画面の E2E の繰り返し 90 件、プレビューの確認の dry run 25 件（1 件 skip）、画面の記録 396 枚（90 件）。コードレビュー 4 回（High 無し。3 回目の Medium 1 件はゲート 1 の落ちた試験、Low は直した。4 回目の Low 1 件は文書で `2b553a4b` が直した）、画面レビュー 2 回とも合格（2 回目は Low 4 件で、上の 5. の最初の変更）。最後のコミットは文書だけの `2b553a4b`（アプリの木は `35970fa1`）。
- **ゲートの script**：[`run_gate.sh`](../evidence/losat_web_w7/run_gate.sh)（W6 の script に W7 の試験・記録を足したもの）と W6 の `check_commands.py`。段階：release ID の確認（未コミットの変更があれば最初に止まる）、V-ABI quick、`npm ci`、`npm run check`、reactor 付きの単体試験、`check_commands.py`（oracle の lock の下）、FakeEngine の E2E、エンジン入りの E2E（サイトの索引と manifest を実行に残す）、プレビューの確認の dry run（ローカルの server）、結果画面の E2E の繰り返し、計測（結果と配信）、画面の記録。`LOSAT_WEB_GATE_STEPS` は `all`・`after-review`・`measure`。`all` には `LOSAT_WEB_NCBI_BIN`（`~/micromamba/bin`）、`LOSAT_WEB_REACTORS`（`$BUILD_ROOT/s14-reactors`）、`LOSAT_WEB_NATIVE`、`LOSAT_WEB_SCREENS_OUT` が要る。Node は 22.18 以上。`NODE_ENV` と `VITE_USER_NODE_ENV` は外される。2 つ目の checkout からのバイト同一の試験が skip されるとゲートは落ちる。各段階の前に V-PERF の lock を待つ。全段階で約 48 分（ゲート 3：2026-10-11 17:09〜17:58 JST）。1 つの失敗で止まるので、直したら最初からやり直すことになる（`after-review`・`measure` で一部だけ走らせられる）。
- **判断**：20 件、すべて Owner-delegated（2026-10-11、2026-09-29 の常設の指示、2026-10-07 に再掲）。要点：release ID は `<commit 12>-<エンジンの digest 8>`、サイトの索引 `/releases.json` はサイト単位、固定は `localStorage['losat-web:release']`、Direct Upload でデプロイ、GA4 は方式 (b)（`/api/visit`、保守者の確認待ち）、V-OFF はサーバーを止めてオフラインにする、オフラインと更新は `/r/<id>/` のときだけ、修正 1 回目から更新は切り替えの隣で影響を言う、修正 4 回目の `vite build` は Vite が production を解決したときしか進まない、Chromium が 1 つのページの 11 個目の download を捨てるので試験は download を 150 ms あける（判断 17、18。判断 19：画面レビュー 2 回目の Low 4 件は S16b に回した。判断 20：Chromium の JSON の 250 ms は W6 との A/B で退行ではない）。全文はゲート記録の「判断」。
- **残件**：ゲート記録の「残件と注意」。主なもの：コードレビュー 1 回目 I1・I4・I6、2 回目 I1〜I4、3 回目 I2、4 回目 I1〜I3、Chromium の JSON の書き出しの最長のフレームの間隔 250 ms（W6 の A/B で退行ではない）、WebKit の scroll anchoring と遅れて読まれる Alignments（W6 の判断 17）、WebKit の書き出しの開始の 128〜155 ms、隠れた tab、`query.fa` / `subject.fa`、W5 コードレビュー L3・W5 画面レビューの I3・I4、BLASTX（SX）、`abi_v2.md` の kind 0、E2j の申し送り、メモリ不足・WebKit の保存の上限・`runInputs`、WebKit の `navigator.onLine`、Playwright の offline の emulation の限界、Chromium が `/r/` の 404 を immutable として覚えること、V-PERF の lock の規則を WP-U の V-OFF の 1 回の実行で破ったこと（記録済み）。
- **画面の記録と画面レビュー**：`$BUILD_ROOT/s16-gate3-screens/`（`/home/kawato/.cache/losat-work/s16-gate3-screens/`、ゲート 3 の木、396 枚）。SHA-256 は実行記録の `screens.sha256`。状態（66 枚 / ブラウザ / サイズ）：01〜28 は W5 と同じ並び、29〜40 と 30b が W6、41〜60（と 46b・46c・55b）が W7（導入支援、オフラインの状態、更新、同意、起動ページ）。S16b で見た目を変えたら、この記録と比べる。
- **E2j（S13+）**：W7 でも merge されておらず、E2j の指示書の 6.（アプリへの申し送り）は行っていない。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S17 — 公開判定](session_s17_g_release_decision.md)。
