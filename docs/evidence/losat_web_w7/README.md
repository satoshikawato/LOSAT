# LOSAT Web W7（Session S16）ゲート記録

- 段階：W7 配信の仕上げ（[総合計画書](../../losat_web_gui_plan.md) §7 の S16、[指示書](../../losat_web_gui_sessions/session_s16_w7_delivery.md)、[release の配置と ReleaseManifest](../../web/release_manifest.md)、[計測](../../web/analytics.md)、[例題](../../web/examples.md)、[Cloudflare へのデプロイ](../../web/deploy_cloudflare.md)、[対応表](../../web/ncbi_ui_mapping.md)の §5「案内と配信」）。完了条件は計画 §7 の S16 の行：V-OFF が通る、V-PRIV が通る、プレビューで `crossOriginIsolated === true`
- ブランチ：`feature/losat-web-gui-app`（アプリ側。worktree `$WORK_ROOT/.worktrees/web-gui-app`、Linux の clone）。起点は `e45fbd09`（W6 のゲート記録と S16 の指示書の引き継ぎ）。エンジン側のブランチ `feature/losat-web-gui` は `141615ec` で、S16 の間に動いたもの（SFe は未コミットの変更を残して止まっている）はこのブランチに取り込んでいない。reactor は `$BUILD_ROOT/s14-reactors`、ネイティブの CLI は `$BUILD_ROOT/s14-gate-native/LOSAT`（sha256 `0236291b…`）で、どちらも `141615ec` の木から作った（W6 のものをそのまま使った。実行記録の `environment.txt` に wasm のハッシュ）。エンジン（`LOSAT/`、`web/adapter/`）、計画、セッションの README、README の表、依存ライブラリは変えていない。`web/app/package.json` は `engines.node` を `>=22.18` にし、script に `assemble` を足し、`preview` と `typecheck` を替えただけ。`docs/web/` は `release_manifest.md`・`analytics.md`・`examples.md`・`deploy_cloudflare.md`（新規）と `ncbi_ui_mapping.md` の §5 を足し、`web/AGENTS.md` は E2E が `build/site-server.ts` で配ることを、`session_file.md` は 1 行を直した。`.github/workflows/web.yml` は browser-engine の job に自己試験の期待値の試験を足した
- 実行記録：[`run-20261011T061357Z/`](run-20261011T061357Z/)（ゲート 1、`ba7fa464`。エンジン入りのビルドの E2E の 1 件で止まった）、[`run-20261011T073637Z/`](run-20261011T073637Z/)（ゲート 2、`9cb4e61e`。プレビューの確認の dry run の 3 件で止まった）、[`run-20261011T080950Z/`](run-20261011T080950Z/)（ゲート 3、`35970fa1`。すべての段階が通った）。どれも作成後は書き換えない。再現は [`run_gate.sh`](run_gate.sh)、プレビューの確認の script は [`check_preview/`](check_preview/README.md)、GA4 の方式の比較は [`ga4_probe/`](ga4_probe/RESULTS.md)、NCBI との比較は W6 の [`check_commands.py`](../losat_web_w6/check_commands.py)（変えていない）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)。画面の記録（PNG）はリポジトリに入れない（下の「残件と注意」）
- 判定：**完了条件を満たしていない（W7 は未完）**。計画 §7 の S16 の行の 3 条件のうち、V-OFF と V-PRIV は**最後の木（`35970fa1`）のゲート 3 で通った**。**3 つめの「プレビューで `crossOriginIsolated === true`」は確かめていない**。Cloudflare Pages のプレビューは保守者のデプロイが要り（README の規則 7）、このセッションでは行われなかった。指示書 7. のプレビューでの確認（CSP、304 に付く COOP / COEP / CORP / CSP、WebKit の取消の後の検索、immutable のキャッシュ）も同じ理由で残る。プレビューの確認の script（`check_preview/`）を本物のサイトに向けた記録は無い。ゲートが走らせた `check_preview/` の dry run はローカルの server（`build/site-server.ts`）に向けたもので、確認の script 自体が動くことの証拠であり、プレビューの証拠ではない。S16b（[指示書](../../losat_web_gui_sessions/session_s16b_w7_preview.md)）が、保守者のデプロイの後に、これを最初に行う。コードレビューは 4 回（1 回目：High 無し、Medium 3 件（M1〜M3）と Low 5 件、情報 6 件。2 回目：High・Medium 無し、Low 5 件、情報 5 件。3 回目：High 無し、Medium 1 件（M1。ゲート 1 の落ちた試験と同じ）、Low 3 件、情報 3 件。4 回目：High・Medium 無し、Low 1 件（文書。直した）、情報 3 件。Medium と Low はすべて直したか、直さない理由を記した）。画面レビュー 1 回目は合格（High・Medium 無し、Low 6 件を直した）、2 回目（ゲート 3 の記録）は合格（High・Medium 無し、Low 4 件（L1〜L4）、参考 7 件。Low 4 件は見た目と文言の変更なので、新しいゲートの実行が要り、S16b が最初の変更として直す。判断 19）
- 保守者の判断待ち（判断 1〜20 はすべて推奨案で進めた Owner-delegated で、保守者が予約した判断は無い。下の「判断」）：
  1. **プレビューのデプロイと確認**：Pages のプロジェクト名と production branch（Direct Upload のプロジェクトは後から dashboard で変えられない）を決め、[`deploy_cloudflare.md`](../../web/deploy_cloudflare.md) の手順でプレビューへデプロイし、`check_preview` の結果を S16b に渡す（URL を渡して実行を許すか、同じ command を実行して JSON 5 つを返す）。
  2. **GA4 の方式の確認**：(b)（アプリが自分の `/api/visit` へ固定のイベントだけを送り、Pages Function が Measurement Protocol へ転送する）を、probe（50 件、3 ブラウザ）の後に Owner-delegated で採った（判断 11）。指示書は推奨を保守者に示すよう求めているので、確認を求める。比較は [`ga4_probe/RESULTS.md`](ga4_probe/RESULTS.md)、実装は [`analytics.md`](../../web/analytics.md)。
  3. **GA4 の property と Cloudflare の secret**：web data stream、`GA4_MEASUREMENT_ID`、Measurement Protocol の API secret（`GA4_API_SECRET`）、custom dimension（event scope の `release_id`）、データの保持期間（短いもの。2 か月を推奨）、Google signals と広告のリンクを切ること。2 つの変数は Cloudflare の Production 環境にだけ置く（プレビューの Function は転送しない）。手順は `analytics.md` の「What the maintainer sets」。
  4. **Cloudflare に残す旧版の数とドメイン名**（計画 §10。「S16 の前」とされていた）：`npm run assemble` は渡された release をすべて残す。減らすなら `--previous` に渡さない（そのときそれを固定したタブは「もう配っていない」と言われ、最新を勧められる）。上限は 1 デプロイ 20,000 ファイル（無料）、1 ファイル 25 MiB。1 つの release は約 15 ファイル・8 MB で、どちらにも遠い。
  5. **同意の文言**：法的な確認が要る。設計は、この記録が法的な判断をしないとする。文言は `ui/AnalyticsConsent.vue` と `analytics.md` の「Consent and withdrawal」。
  6. **本番へのデプロイ**：プレビューで確認した同じ `dist/` を production branch に出す（再 build しない）。S17 の公開判定の後。
- `feature/losat-web-gui` への merge：このセッションでは行わない。エンジン側が行う（下の「合流」）
- このセッションは取りまとめ役（Opus）と、作業単位 WP-R（配置・ReleaseManifest・起動ページ・headers・server・組み立て）、WP-U（Service Worker・Offline ready・明示更新・V-OFF）、WP-G（GA4 の probe と実装）、WP-O（例題・導入支援）、WP-P（V-PRIV・TD-11）、WP-C（プレビューの手順と確認、Sonnet）、WP-GP（ゲートの準備）の agent、修正 3 回（fix 1〜3）で行った。契約（release ID、サイトの索引）は取りまとめ役が先に入れ、各 agent は一時的な worktree と作業用の branch で作業し、取りまとめ役が `feature/losat-web-gui-app` に merge した。コードレビューと画面レビューは別の agent（`losat-reviewer`）。ゲートの script の実行は取りまとめ役が背景で行った

## コミット

`e45fbd09..2b553a4b` の 84 コミット（merge 15 を含む）と、この記録のコミット 2 つ（古い順）。`9467e56f`・`bf2c1b7c`・`a49cedc1`・`0ed740cf`・`8fea746a`・`52a7d6de` は agent の作業用 branch にアプリ側の先頭を取り込んだ merge。修正 3 回目（`c458e142`〜`9cc355e6`）、修正 4 回目（`3824aa9c`〜`b2a6e066`）、修正 5 回目（`9cb4e61e`）と、ゲートで見つけた 2 件の試験の修正（`4c5bf901`、`35970fa1`）はアプリ側に fast-forward で入った。

| コミット | 内容 |
|---|---|
| `55c924d8` | 契約（アプリ）：ビルドの release ID（`__LOSAT_RELEASE_ID__`）。release の directory、固定、計測が使う |
| `f2f6af8c` | 契約（アプリ）：サイトの索引 `/releases.json`（domain、port、ブラウザの reader）。更新の確認と計測の切り替えが共有する |
| `31e2df5a` | アプリ（WP-R）：release の directory `/r/<release-id>/`、ReleaseManifest、起動ページ（指示書 1） |
| `f67d1b74` | アプリ（試験だけ、WP-P 第 1 段階）：V-PRIV と TD-11 |
| `4aa222a7` | アプリ（WP-U 第 1 段階）：オフラインと明示更新：Service Worker の論理、ReleaseManager、見出しの状態、更新の知らせ（指示書 2・3） |
| `21dffad4` | アプリ（WP-R）：パスごとの `_headers`（`/r/*` は immutable、サイトのファイルは no-cache）と、Pages と同じに配るローカルの server |
| `585388e5` | merge：WP-R（配置の部分） |
| `9467e56f` | merge：配置の部分を WP-U の branch に取り込み |
| `94fd8d50` | アプリ（WP-R）：古い release を持つサイトの組み立て（`npm run assemble`） |
| `b04f0813` | アプリ（WP-O）：例題 5 つ、「Try an example」、Getting started、ブラウザ自己試験（指示書 6） |
| `973b8b2c` | 文書：LOSAT Web の例題（file、元のレコード、領域、SHA-256、自己試験の期待値） |
| `2258c184` | アプリ（WP-O）：Getting started の枠は閉じている間も mount したまま（見出しの `aria-controls` が指す先を残す） |
| `f8c5720b` | アプリ（WP-U 第 2 段階）：サイトの Service Worker `/sw.js`：サイトのビルドが作り、起動ページが release の開始の後に登録する（指示書 2、第 2 段階） |
| `b7331f61` | 文書（WP-R）：`release_manifest.md` にサイトのビルドが消さない directory |
| `05ea4a00` | merge：WP-R（残り：組み立て、文書） |
| `64ebe0c4` | 文書：`web/AGENTS.md` に `npm run e2e` が `build/site-server.ts` で配ること |
| `115d8fa7` | 文書（WP-G）：GA4 の probe：iframe の計測文書 (a) と同一 origin の endpoint (b) を COEP の下で 3 ブラウザで比べる（指示書 4） |
| `bf2c1b7c` | merge：アプリ側の先頭を WP-U の branch に取り込み |
| `fec62539` | アプリ（WP-G）：同意の後の利用統計：サイト自身の `/api/visit` へ（計画 §5.10 (b)、probe の結果） |
| `8d150b56` | 文書（WP-G）：`analytics.md`：何を送るか、同意、Pages Function、保守者の設定 |
| `423ce55d` | merge：WP-O |
| `2bdbc377` | アプリ（CI）：browser-engine の job で自己試験の期待値をシリアルの reactor に対して検査する |
| `d24a0a78` | アプリ（試験だけ、WP-U）：V-OFF：オフラインの開始、混ざらない明示更新、不完全な更新、旧版のタブのデータ、サイトがもう配らない固定（計画 §6.2、指示書 2・3・7） |
| `6bd5b131` | アプリ（試験だけ、WP-C）：デプロイしたサイトの確認：200・304・404 の headers、bytes、起動ページ、`/api/visit`、3 ブラウザの開始・取消・オフライン |
| `10ba843e` | 文書（WP-C）：Cloudflare Pages へのデプロイ（build、assemble、Direct Upload、プレビューの確認、Pages の挙動） |
| `def52250` | merge：WP-C |
| `739508dd` | merge：WP-U |
| `68891222` | 文書：サイトの配置と headers の表に `/sw.js` |
| `a49cedc1` | merge：アプリ側の先頭（WP-R・U・O・C）を WP-P の branch に取り込み（V-PRIV の第 2 段階のため） |
| `fb383193` | アプリ（試験だけ、WP-G）：計測のオフラインの場面は入力を先に読み、WebKit は offline の emulation ではなく HTTP の要求を中断する |
| `acd0f00b` | merge：WP-G。**コードレビュー 1 回目の木** |
| `1ec45632` | アプリ（試験だけ、WP-P 第 2 段階）：配置に載せた V-PRIV：`build/site-server.ts` でサイトを配って捕まえる |
| `0ed740cf` | merge：WP-G を WP-P の branch に取り込み（計測の経路のため） |
| `9c6ecea1` | アプリ（試験だけ、WP-GP）：W7 の画面の記録（状態 41〜56）。見出しが「Offline ready」に落ち着いてから撮る |
| `5e860b74` | アプリ：`register()` が何も返さないときは「service workers are blocked in this browser」（状態 56 で見つけた） |
| `7ad7b28a` | アプリ（試験だけ、WP-P）：V-PRIV の配信の経路：更新、「Try an example」、自己試験、計測 |
| `a0183a56` | アプリ（試験だけ）：結果の計測は「Offline ready」の後に、Getting started を閉じて始める |
| `d607a2f5` | アプリ（試験だけ）：配信の計測：初回・2 回目・オフラインの訪問、準備、自己試験、例題 |
| `1a089486` | アプリ：`src/sw`・`src/boot` の ESLint の層の規則（コードレビュー 1 回目 L5） |
| `232800d6` | アプリ：release ID の commit と `app.build` を固定の長さにし、ビルドの入力の全体で `-dirty` を見る（コードレビュー 1 回目 M3） |
| `5375c139` | merge：WP-P |
| `a19a7e11` | アプリ：更新は切り替えの隣で、その時点で数えた「再読み込みが何を失うか」を言う（コードレビュー 1 回目 M1） |
| `8fea746a` | merge：WP-P を修正 1 回目の branch に取り込み |
| `96f9bc84` | 文書：W7 のゲートの script（`run_gate.sh`） |
| `3f74058b` | アプリ：取得して確かめた release は切り替えまで保持し、切り替えで保存されていることを確かめる（コードレビュー 1 回目 M2） |
| `405ad0ca` | アプリ：検査なしの切り替えは、Service Worker の無いブラウザだけに出す（コードレビュー 1 回目 L4） |
| `788c23ac` | アプリ：「Offline ready」は読み戻した固定がこの tab の release を指すときだけ（コードレビュー 1 回目 I5） |
| `895d2949` | アプリ：見出しのオフラインの状態はどの状態でも同じ大きさ（コードレビュー 1 回目 M4、WP-P の全体の実行で見つけた） |
| `56e8ee4a` | アプリ（試験だけ）：利用統計の E2E は Service Worker を塞ぐ |
| `1ca293df` | アプリ（試験だけ）：`analytics.spec.ts` は Service Worker 無しで走る（WebKit は制御下のページの要求を `page.route` に渡さない） |
| `daca999d` | アプリ：HTTP cache が覚えた `/r/` の 404 は、それを越えて読み直す（コードレビュー 1 回目 L2） |
| `70442212` | merge：WP-GP。**画面レビュー 1 回目の木** |
| `d4736ac1` | アプリ（試験だけ）：V-OFF は release の file を no-store で配り、オフラインの開始が Service Worker の cache を証明するようにする（コードレビュー 1 回目 L1） |
| `52a7d6de` | merge：WP-GP を修正 1 回目の branch に取り込み |
| `3db6c216` | アプリ：塞がれた登録は Service Worker 無しと数える（コードレビュー 1 回目 L4 の続き） |
| `71150eef` | アプリ：「Try an example」は設定と Job Title を置き換えるときも尋ねる（コードレビュー 1 回目 L3） |
| `c283324a` | 文書：`analytics.md` と `deploy_cloudflare.md` がプレビューの計測の endpoint で一致する（コードレビュー 1 回目 I3） |
| `08cd5359` | アプリ（試験だけ）：利用統計の E2E は Function が要る `Origin` を確かめる（コードレビュー 1 回目 I2 の前半） |
| `a9507eac` | アプリ：見出しのオフラインの状態は 1 行、詳細は知らせへ、他の tab の知らせは再読み込みで失うものを言う（画面レビュー 1 回目 L1・L3、コードレビュー M4・M1）。**修正 1 回目の終わり、コードレビュー 2 回目の木** |
| `0877440d` | アプリ：release を保存した後に保管場所の見積りを読み直す（画面レビュー 1 回目 L2） |
| `af821677` | アプリ：失敗した更新は、どの file がなぜかを言い、エラーの枠に出す（画面レビュー 1 回目 L5） |
| `3b457430` | アプリ：リンクのボタンは 24 px の標的、無効な「Try an example」は無効に見える、電話の小さい幅の折り返し（画面レビュー 1 回目 L6・I2・I3） |
| `9e9ee11e` | アプリ：読み込んだ Run の知らせは build を commit として言う。見出しの release ID の頭と同じ桁（画面レビュー 1 回目 I5） |
| `cf017d9d` | アプリ：起動ページは「Try again」とこのブラウザに保存された release を出し、自分では切り替えない（画面レビュー 1 回目 I6） |
| `163a22fb` | アプリ（試験だけ）：画面の状態 57〜60：更新の取得中、壊れた file でのオフライン失敗、待ち、測るサイトの初回（画面レビュー 1 回目 L4） |
| `79f724ad` | 文書：`ncbi_ui_mapping.md` に W7 の部品（画面レビュー 1 回目 I4） |
| `19ede7f0` | アプリ：サイトに Service Worker を許さないブラウザ（`register()` が `SecurityError`・`NotSupportedError`）は、毎回の失敗ではなくオフライン無しとする（コードレビュー 2 回目 L1） |
| `c458e142` | アプリ（試験だけ）：V-PRIV は要求の全体を見る：サイトへの query string、query や fragment のあるアドレス、base64・パーセント符号化した試験文字列は指摘とする（コードレビュー 2 回目 L2） |
| `ac6ef3dc` | 文書：ゲートは最初の段階の前にツリーの release ID を確かめ、測ったサイトの記録を残す（コードレビュー 2 回目 L4・I5） |
| `0c5956dc` | アプリ（試験だけ）：バイト同一の試験は、別のパスの 2 つ目の checkout からも build して比べる（コードレビュー 2 回目 L5） |
| `aeca5d93` | 文書：ゲートは `NODE_ENV` を外して build する（コードレビュー 2 回目 L5 で見つけた） |
| `9cc355e6` | アプリ（試験だけ）：V-PRIV の捕獲の自己試験が、入れ子の worker と Service Worker が他の origin に届かないことを 3 ブラウザで示す（コードレビュー 2 回目 L3）。**修正 3 回目の終わり** |
| `cc5ca34f` | アプリ：`vite build` は `NODE_ENV` が production 以外なら拒む（修正 3 回目で見つけた：開発 build は絶対パスを bundle に入れる） |
| `9f0f3bac` | アプリ（試験だけ）：状態 51 は見出しの「stored」を待つ（修正 1 回目が、他の tab が切り替えた tab にこの状態を与えるため）。**修正 2 回目の終わり** |
| `ba7fa464` | merge：修正 2 回目（画面レビュー 1 回目の L2・L4〜L6・I2〜I6、コードレビュー 2 回目 L1、`NODE_ENV` の規則）。**ゲート 1 の木、コードレビュー 3 回目の木** |
| `4c5bf901` | アプリ（試験だけ）：V-OFF のケース 1 は、見出しの保管場所の見積りを、file の大きさではなくブラウザ自身の見積りと比べる（ゲート 1 の落ちた試験、コードレビュー 3 回目 M1：Firefox は Cache API の body を圧縮して保存する） |
| `3824aa9c` | アプリ：`vite build` は、`NODE_ENV`・`.env` の file・`VITE_USER_NODE_ENV` が何を言おうと、Vite が production を解決したときしか進まない（コードレビュー 3 回目 L1） |
| `d3812c2e` | アプリ：`register()` の `SecurityError` が Service Worker 無しを意味するのは、ブラウザが `getRegistration()` も拒むときだけ。JavaScript でない `sw.js` は「Try again」付きの失敗のまま（コードレビュー 3 回目 L2） |
| `0d7f439e` | アプリ：起動ページの「not served」で、出せる最新の release が無いときも「Try again」と保存された release を出し、行き止まりにしない（コードレビュー 3 回目 L3） |
| `72ae82c9` | アプリ：起動ページの「Use release <id> (stored)」は click の時点でもう一度確かめ、まだ全部保存されている release だけを固定する。そうでなければそう言って残りを出す（コードレビュー 3 回目 I1） |
| `b2a6e066` | 文書：決定性の文は `release-bytes.test.ts` が示すこと（1 台の機械、1 つの Node、どの checkout のパスでも）に絞る。2 つ目の checkout の case が skip されるときゲートは落ちる（コードレビュー 3 回目 I3）。**修正 4 回目の終わり** |
| `9cb4e61e` | アプリ（試験だけ）：1 つのページの download は 150 ms あけて行う（Chromium は、バーストの最初から 1 s 以内の 11 個目の download を黙って捨てる。`privacy.spec.ts` の flaky、修正 5 回目）。**ゲート 2 の木** |
| `35970fa1` | アプリ（試験だけ）：プレビューの確認は、見出しの offline-status の `data-state` が「ready」になるのを待つ（`offline-ready` の test ID は修正 1 回目でなくなっていた。ゲート 2 の dry run の落ちた 3 件）。**最後のアプリの木、ゲート 3 の木** |
| `2b553a4b` | 文書：起動ページの選択肢（最新が無い「not served」、click の時点の再確認、`boot-stored-gone`）を `release_manifest.md` と `ncbi_ui_mapping.md` に（コードレビュー 4 回目 L1。アプリのバイトは変わらない。アプリの木は `35970fa1` のまま） |
| ゲートの実行記録のコミット | 文書：W7 のゲートの 3 回の実行記録（`run-20261011T061357Z/`・`run-20261011T073637Z/`・`run-20261011T080950Z/`） |
| このゲート記録を含むコミット | 文書：W7 のゲート記録と `evidence.sha256` |
| S16b・S17 への引き継ぎのコミット | 文書：S16b の指示書（新規）と、S17 の指示書に「S16（W7）から引き継ぐこと」（セッションの README の表はエンジン側が直す。下の「合流」） |

最後のアプリの木は `35970fa1`（アプリの本体は修正 4 回目の `72ae82c9` から変わっていない。後は試験・文書だけ。`2b553a4b` は文書だけで、ゲート 3 の木のアプリのバイトと同じ）。エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表、セッションの README は、アプリ側のブランチでは変えていない（README 規則 1）。

## 完了条件と結果

計画 §7 の S16 の行と指示書 1.〜7.。試験の名前は `web/app/tests/` のもの。「通過」は、作業単位と修正の branch での記録（`logs/` の各 progress）に基づく。最後の木のゲート（ゲート 3）の全体の結果は「ゲートの実行」に書く。

| 条件 | 結果 | 証拠 |
|---|---|---|
| 計画 §7 S16：V-OFF が通る | 通過。ゲート 3（`35970fa1`）で、`offline.spec.ts` が FakeEngine のビルドでもエンジン入りのビルドでも、3 ブラウザとも 12 件ずつ通った | 下の 2.・3.、`offline.spec.ts`、ゲート 3 の `npm-e2e-fake-engine.log`・`npm-e2e-engine.log` |
| 計画 §7 S16：V-PRIV が通る | 通過。ゲート 3 で、`privacy.spec.ts` が両方のビルドで Chromium 7 件・Firefox 6 件・WebKit 6 件通った（TD-11 は Chromium だけ）。`privacy-scan.test.ts` を含む単体試験も通った | 下の 5.、`privacy.spec.ts`、`privacy-scan.test.ts` |
| 計画 §7 S16：プレビューで `crossOriginIsolated === true` | **未確認**。保守者のデプロイが要り、行われなかった。ローカルの server（`build/site-server.ts`）に向けた dry run では、`check_preview` の開始の確認が 3 ブラウザで `crossOriginIsolated === true` を見て通った（確認の script の動作の証拠。プレビューの証拠ではない） | S16b。`check_preview/`、ゲート 3 の `preview-check/`（dry run：25 件が通り 1 件 skip（`/api/visit` の有効な event は `LOSAT_WEB_PREVIEW_SEND_VISIT=1` のときだけ）。http 16 件、各ブラウザ 3 件） |
| 1. release ごとの資産：`/r/<release-id>/…`、小さな起動ページ、`localStorage` の固定（try/catch）、ReleaseManifest | 通過。release ID は `LOSAT_WEB_RELEASE_ID` があればそれ、無ければ `<commit 12>-<エンジンの digest 8 または fake>[-dirty]`（`^[0-9a-z][0-9a-z.-]{0,63}$`）。同じ機械・同じ Node なら、checkout のパスによらず同じ bytes（`release-bytes.test.ts`。別の機械や別の Node での同一は示していない。判断 17）。`app.build` は commit の先頭 8 桁。ビルドの入力は `build/inputs.ts` に 1 度だけ列挙し、そのどれかに未コミットの変更があれば `-dirty`。git が無ければ `LOSAT_WEB_RELEASE_ID` が要る。manifest は format 1（`release`・`app`・`engine`・`entry`・`assets` の path・size・SHA-256）。サイトの索引 `/releases.json` は release の外（サイト単位）。起動ページは固定を読み、サイトが配らないと言われても黙って切り替えない。`npm run assemble` が古い release を新しいサイトに運ぶ（manifest に照らして検査する） | `release.spec.ts`（固定、配らない固定、索引なし、`localStorage` が使えない、全応答の headers、404、見出しの大きさ）、`release-id.test.ts`、`release-manifest.test.ts`、`site-index.test.ts`、`site-build.test.ts`、`boot.test.ts`、`assemble-site.test.ts`、`release-bytes.test.ts`（同じ commit の 2 回の build と、別のパスの 2 つ目の checkout からの build が同じ bytes）、`site-server.test.ts`、`headers.test.ts` |
| 2. Service Worker：固定した release の資産（2 つの wasm、Worker、例題）を取得して SHA-256 を確かめ、それから「Offline ready」。研究データはキャッシュに入れない | 通過。全 file のサイズと SHA-256 が合ったときだけ、manifest を最後に保存し、それから shell（`/`・`/boot.js`・`/boot.css`）。壊れた取得は何も保存せず、「Try again」で保存できる。Service Worker は release の file とその shell だけを答え、書くのは準備だけ。`/releases.json`・`/api/*`・`/sw.js`・他の origin・GET 以外には答えない。`install` の handler は無く、`skipWaiting` は利用者が切り替えるときだけ。実測：初回の訪問で「Offline ready」まで約 2.2〜2.5 s（下の「実測」） | `offline.spec.ts` の 1・1b・2、`sw-prepare.test.ts`、`sw-serve.test.ts`、`offline-client.test.ts`、`privacy.spec.ts`「an update while the tab holds research data … the caches hold only release files and the shell」 |
| 3. 明示更新：知らせる → 選ぶ → 取得して確かめる → 影響を案内 → 固定を切り替えて再読み込み。不完全なら今の release のまま。旧版の資産はそれを使うタブがある間は消さない | 通過。切り替えの隣に、その時点で数えた「再読み込みが失うもの」（実行中・待ちの Run、完了した Run、候補とメモ。Save session を勧める）を出す（コードレビュー M1）。取得した release は shared lock で保持し、切り替えで保存されていることを確かめる。消えていれば「removed」で「Download and check again」（M2）。他の tab が切り替えたことは再読み込みで失うものと一緒に言う。旧版の cache は、その release の lock を持つ tab が無くなってから消す | `offline.spec.ts` の 3・3b・4・5・6・6b・6c、`release.test.ts`、`sw-prepare.test.ts` |
| 4. GA4：(a) と (b) を COEP の下で全ブラウザで試し、結果とともに推奨を示す。決まった方式を実装する。同意の前は何も読み込まず何も送らない。塞がれていても未同意でもオフラインでも検索・結果・抽出・保存が動く | 通過（方式の確認は保守者の判断待ち）。probe は 50 件が通り、(b) は 3 ブラウザで同じに動き（204、body は固定、Referer・Cookie 無し、CSP の違反無し）、(a) は M の COEP と `CORP: cross-origin`、`frame-src`、`allow-same-origin` が要り、3 ブラウザで同じには動かない。(b) を採った（判断 11）：`{"v":1,"release","client"}` を同意の後に `/api/visit` へ。Pages Function は正確にこの文だけを受け、Google のスクリプトはどの文書にも無い。未同意・「Don't allow」・storage の失敗では何も送らず、塞がれた・500・オフラインの endpoint でも検索から保存まで動く | `ga4_probe/RESULTS.md`、`analytics.spec.ts`（WP-G の記録で 10 件 × 3 ブラウザ）、`analytics.test.ts`、`visit-function.test.ts`、`privacy.spec.ts` の計測の経路 |
| 5a. V-PRIV（計画 §5.10）と TD-11：Playwright の捕獲で、どの経路でも許した origin とだけ通信し、試験文字列を含む要求が無い。配信する Wasm と JavaScript にビルドした機械のパスが無い | 通過（ゲート 3）。経路：通常の検索・結果・全書き出し・再読込・セッション・つなぎ直し、取消、拒否した入力と option、失敗した Run、更新、「Try an example」と自己試験、計測。試験文字列は defline・file 名・配列・Job Title・メモ・設定・セッションの名前・不正な入力に入れ、平文・大文字小文字・パーセント符号化・base64（3 つのずらし、2 つの字母）・空白や JSON のエスケープを除いた形で探す。サイトへの query string と `?`・`#` のあるアドレスは指摘。body は固定の 1 件の POST 以外すべて指摘。TD-11 はサイトが配る全 file（shell、`releases.json`、`sw.js`、404、各 manifest の file と wasm）を読み、checkout・`HOME`・`CARGO_HOME`・`RUSTUP_HOME`・ユーザー名のパスを探す。捕獲の自己試験が、試験文字列の URL・他の origin への fetch・入れ子の worker と Service Worker の他の origin への fetch を実際に見つける（または塞がれる）ことを 3 ブラウザで示す | `privacy.spec.ts`（6 つの試験と TD-11）、`privacy-scan.test.ts`、CI の `web` と `browser-engine` の両方で `npm run e2e` |
| 5b. V-OFF（計画 §6.2）：オフラインでの起動、明示更新、新旧の資産が混ざらない | 通過（ゲート 3）。オフラインはサーバーを止めて作る（Playwright の `setOffline` は Chromium でしか確かでない。Firefox は Service Worker の要求に効かず、WebKit はページの遷移そのものを壊す。判断 6、10）。Chromium と Firefox は `setOffline` も併用する | `offline.spec.ts`（1・1b・2・3・3b・4・5・6・6b・6c）、`site-server.test.ts` |
| 6. 導入支援：「Try an example」、5 program の小さな例題、最短のチュートリアル、任意のブラウザ自己試験（期待値の SHA-256、V-BR の小さな版） | 通過。例題は DnaA とその遺伝子（*Candidatus* Hepatoplasma / Tyloplasma の 4+1 の genome の先頭 2,000 nt とタンパク質）。5 program × outfmt 0/6/7 の期待値 15 件は、ネイティブの LOSAT（`141615ec`）が書くバイトで、NCBI BLAST+ 2.17.0 と 15 件とも一致。BLASTX の例題は release に入っているが、SX まで「Try an example」と自己試験は無効（理由を言う）。自己試験は 1 スレッドと Auto のスレッドの 2 回、キュー・結果・候補・セッションに出ない。Getting started は初回に開き、閉じれば次から閉じたまま | `onboarding.spec.ts`（例題ごとの入力と検索、BLASTX の無効、置き換える前に尋ねる、リンクの標的の大きさ、Getting started の開閉、自己試験の全 case）、`onboarding.test.ts`、`self-test.test.ts`、`self-test-engine.test.ts`（reactor 付き 4 件、CI の browser-engine job でも）、`examples.md` |
| 7a. Cloudflare のプレビューの設定（`_headers` をそのまま使う）と手順 | 通過（準備）。`public/_headers` を唯一の source にし、`/r/*` は `Cache-Control: public, max-age=31536000, immutable`、サイトの file は `no-cache`。トップの `404.html` を作る（無いと Pages は単一ページとして未知のパスに `/` を 200 で返す）。デプロイは Direct Upload（Pages の Git build には Rust が無い）、`web/app` で `wrangler pages deploy dist`（`functions/` を一緒に上げる）。手順は `deploy_cloudflare.md` | `deploy_cloudflare.md`、`headers.test.ts`、`release.spec.ts` |
| 7b. プレビューで `crossOriginIsolated === true` と CSP、304 に COOP / COEP / CORP / CSP、WebKit の取消の後の検索 | **未確認（S16b）**。確認の script は書いて、ローカルのサイトで動かした（FakeEngine のサイト：http 13 件・各ブラウザ 1 件が通り 10 件が理由付きで skip、エンジン入りのサイト：26 件が通り、COEP を 304 から外したサイトでは 11 件が意図どおり落ちた） | `check_preview/`、WP-C の dry run（`logs/wp-c/`）、ゲート 3 の `preview-check/`（25 件通過、1 件 skip） |
| 7c. ハッシュ付きの資産に immutable のキャッシュ | 通過（検討して採用）。`/r/*` に付ける。Chromium は `/r/` の 404 も immutable として覚える（WP-U の probe で確認）ので、起動ページと Service Worker は 404・410・JSON でない答えを `cache: 'reload'` で 1 度読み直してから「配っていない」と言う。Pages が 404 にも immutable を付けるかは未確認（プレビューの確認が `notes` に残す） | `release.spec.ts`、`offline.spec.ts` の 6c、コードレビュー 1 回目 L2 |
| 7d. 旧版のタブのデータを新版が消さない（`tmp/` と Web Lock `losat-web:tmp:` の名前を変えない） | 通過。名前は W6 のまま（`release-pin.test.ts` が `SESSION_LOCK_PREFIX === 'losat-web:tmp:'`・`TMP_DIRECTORY === 'tmp'` を固定）。旧版のタブは、別の tab が B に切り替えた後も A の結果と release を保つ | `offline.spec.ts` の 5、`release-pin.test.ts` |
| S15 から引き継いだ条件：研究データは Service Worker のキャッシュと計測に入れない | 通過。cache の中身は release の file と shell だけ（`expectedCaches` と byte で比べる）。書き出し（互換出力、CSV、JSON、レポート、SVG、設定、セッション、抽出した FASTA、Run の入力）は Blob の download で、`fetch` の応答にならない。計測のイベントに file の名前・内容・大きさ・件数は無い | `privacy.spec.ts`、`session-save-load.test.ts`（W6 から変えていない） |
| 3 ブラウザ・両方のビルドで既存と新しい E2E が通る | 通過（ゲート 3）。W6 の最後の木の FakeEngine 146 件、エンジン入り 170 件から、`release.spec.ts`・`offline.spec.ts`・`analytics.spec.ts`・`onboarding.spec.ts`・`privacy.spec.ts` などが増えた。FakeEngine のビルドは 297 件（Chromium 104・Firefox 97・WebKit 96）、エンジン入りのビルドは 324 件（113・106・105。V-BR を含む）が通り、落ちたもの・flaky は無い。結果画面の E2E の繰り返し 90 件（30 件ずつ）も通った | ゲート 3 の `npm-e2e-fake-engine.log`・`npm-e2e-engine.log`・`e2e-results-repeat.log` |
| HSP の対応・V-ABI・NCBI との比較コマンドが変わらない | 通過（ゲート 3）。V-ABI quick：15 の検索 × 経路とスレッド = 60 の実行がネイティブの CLI と一致、NCBI の凍結ハッシュ 16/16 一致。reactor 付きの単体試験は 58 ファイル、969 件通過・1 件 skip（`reproduce.test.ts` の check_commands.py 用の書き出しの 1 件）で、`hsp-correspondence.test.ts` は W4 から変えていない。NCBI との比較コマンドは 42 件（一致 33、承認済みの例外 6、NCBI が拒む 3）。エンジンも W6 の比較の script も変えていない | ゲート 3 の `v-abi.log`・`unit-cases.log`・`check-commands.tsv` |

### ゲートの実行

3 回走らせた。どの実行も、release は最初の段階の前に `build/release.ts` の規則で確かめ、`environment.txt` に記録した（ビルドの入力に未コミットの変更が 0 件）。段階は W6 の段階に W7 の試験・記録を足したもの。時刻は JST、終わりは実行の folder の最後の file の時刻。

**ゲート 1（`ba7fa464`）**：2026-10-11 15:13〜15:39、約 26 分で止まった。release は `ba7fa464ad4e-fake` / `ba7fa464ad4e-e2fb7128`。V-ABI quick（60 実行が一致、凍結ハッシュ 16/16）、`npm run check`（895 件通過・65 件 skip）、reactor 付きの単体試験（58 ファイル、959 件通過・1 件 skip）、`check_commands.py`（42 の比較）、FakeEngine のビルドの E2E（291 件）が通った。エンジン入りのビルドの E2E は 317 件が通り 1 件が落ちた：Firefox の `offline.spec.ts:152`（V-OFF のケース 1）で、見出しの保管場所の見積り 5,400,000 byte が「保存した release の大きさ（約 8.2 MB）の 95 %」を下回った。Firefox は Cache API の body を圧縮して保存し、`navigator.storage.estimate()` は圧縮後の大きさを数える（修正 2 回目の画面の記録にも Firefox で「about 5.4 MB」とある）。アプリの動きは正しく、試験の下限が誤りだった。コードレビュー 3 回目の M1 が同じことを指摘していた。`4c5bf901` で、見出しの値をブラウザ自身の見積りと比べる試験にし、2 つのビルド × 3 ブラウザで通した。ここで script が止まった（`set -e`）ので、サイトの記録、プレビューの確認、計測、画面の記録はこの実行に無い。

**ゲート 2（`9cb4e61e`）**：16:36〜17:09、約 32 分で止まった。release は `9cb4e61ee59e-fake` / `9cb4e61ee59e-e2fb7128`（修正 4・5 回目を含む木）。V-ABI quick、`npm run check`（905 件通過・65 件 skip）、単体試験（969 件通過・1 件 skip）、`check_commands.py`（42 の比較：一致 33、承認済みの例外 6、NCBI が拒む 3）、FakeEngine のビルドの E2E（297 件）、エンジン入りのビルドの E2E（324 件）、サイトの索引と manifest の記録が通った。プレビューの確認の dry run は 25 件のうち 3 件が落ちた（3 ブラウザの「Offline ready」の確認が 3 分待って時間切れ）：確認の script が `offline-ready` の test ID を待っていたが、この ID は修正 1 回目で見出しの `offline-status` の `data-state` に替わっていた。確認の script の側の誤りで、アプリの誤りではない。`35970fa1` で `data-state="ready"` を待つようにし、dry run が 25 件通り 1 件 skip になった。ここで止まったので、結果画面の繰り返し、計測、画面の記録はこの実行に無い。dry run の失敗の trace（`test-results/`、3 件分で 62 MB）は、リポジトリに入れず、実行の folder から出してタスクフォルダの `logs/gate2-preview-test-results/` に置いた（ゲート 3 の実行の folder は変えていない）。

**ゲート 3（`35970fa1`、最後のアプリの木）**：17:09〜17:58、約 48 分。**すべての段階が通った**（[`run-20261011T080950Z/`](run-20261011T080950Z/)、`gate run passed`）。release は FakeEngine のビルドが `35970fa18885-fake`、エンジン入りのビルドが `35970fa18885-e2fb7128`。
- V-ABI quick：15 の検索 × 経路とスレッド = 60 の実行がネイティブの CLI と一致、NCBI の凍結ハッシュ 16/16 一致
- `npm run check`：55 ファイルが通り 3 が skip、905 件通過・65 件 skip（reactor の要る分は skip）
- reactor 付きの単体試験：58 ファイル、969 件通過・1 件 skip（`release-bytes.test.ts` の 2 つ目の checkout の case も走った）
- `check_commands.py`（NCBI BLAST+ 2.17.0 との比較、oracle の lock の下）：42 の比較のうち一致 33、承認済みの例外 6、NCBI が拒む 3
- FakeEngine のビルドの E2E：297 件（Chromium 104・Firefox 97・WebKit 96）、216 skip、約 8.0 分
- エンジン入りのビルドの E2E：324 件（113・106・105。V-BR を含む）、189 skip、約 9.4 分
- プレビューの確認の dry run（エンジン入りの `dist` をローカルの server で配る）：25 件通過、1 件 skip（`/api/visit` の有効な event は `LOSAT_WEB_PREVIEW_SEND_VISIT=1` のときだけ）。http 16 件、各ブラウザ 3 件（開始、取消の後の検索、Offline ready）
- 結果画面の E2E の繰り返し：90 件（3 ブラウザ 30 件ずつ）、約 3.5 分
- 計測（結果と配信、3 ブラウザ）：通った。記録された失敗は W6 のゲート 3 と同じ見込みどおりの組だけ（4000 写しの自己検索は 3 ブラウザでエンジンのメモリ不足、3000 写しは WebKit で保存の上限 512 MB）。下の「実測」
- 画面の記録（3 ブラウザ × 2 サイズ × 66 枚 = 396 枚、状態 01〜60 と 04・30b・46b・46c・55b の枝）：90 件が通った。SHA-256 は `screens.sha256`（396 行）。置き場は `$BUILD_ROOT/s16-gate3-screens/`
- 落ちた試験、flaky、やり直し：無い

修正 4 回目の後の FakeEngine の E2E（`logs/fix4/`）で、`privacy.spec.ts:150` が Chromium で 1 回だけ「extract-download の後に download の event が無い」で落ちた（2 回落ち、その後 5 回通った）。原因は Chromium が、1 つのページで最初の download から 1 s 以内の 11 個目の download を黙って捨てること（10 個を続けて出し、11 個目を +816 ms で出すと捨てられ、+1316 ms で出すと残った。ユーザーの操作でも解けない。Firefox に制限は無い）。アプリは 1 回の click に 1 つの download しか作らないので、人は届かない。試験が 1 つのページの download を 150 ms あけて行うようにした（`tests/e2e/support/downloads.ts`、修正 5 回目、`9cb4e61e`）。繰り返し 20/20 で通った。

どの実行の E2E の log にも W4b 以来の `FAIL [opfs] storage full …` と `FAIL [run output] storage full …` が出る（ブラウザの中の契約の確認が自分の記録として印字するもので、落ちた試験ではない）。

**ゲートの script の規則（W7 で足したもの）**：`NODE_ENV` を外して build する（修正 4 回目から `VITE_USER_NODE_ENV` も外す）。最初の段階の前に、release ID を `build/inputs.ts` の入力に対して確かめ、未コミットの変更（`-dirty`）があれば止める（コードレビュー 2 回目 L4）。E2E の後と各計測の後に、配ったサイトの `releases.json` と各 release の manifest を実行に残し、別の release が配られていれば止める。2 つ目の checkout からのバイト同一の試験が skip されるときは落ちる（`LOSAT_WEB_SECOND_CHECKOUT=required`、コードレビュー 3 回目 I3）。`all` には、プレビューの確認の dry run（`check_preview/`）が入る。

**規則を破った 1 件**：WP-U の V-OFF の実行（Firefox と WebKit）が 1 回、V-PERF の `vperf.lock` が存在する間に始まった（README の規則 5）。気づいたのは後で、lock は数分で消えていた。ゲートの script は各段階の前に lock を待つ。

## 作業との対応

| 作業 | 実装 | 試験 |
|---|---|---|
| release の配置と起動ページ（指示書 1、判断 1〜3、5、14） | `build/release.ts`・`build/inputs.ts`・`build/site.ts`・`build/site-files.ts`・`build/headers.ts`・`build/site-server.ts`・`build/assemble-site.ts`、`src/domain/release-manifest.ts`・`site-index.ts`・`site-file.ts`、`src/boot/`（`plan.ts`・`boot.ts`・`index.html`・`boot.css`・`not-found.html`）、`public/_headers` | `release.spec.ts`、`release-id.test.ts`、`release-manifest.test.ts`、`site-index.test.ts`、`site-build.test.ts`、`boot.test.ts`、`assemble-site.test.ts`、`release-bytes.test.ts`、`site-server.test.ts`、`headers.test.ts` |
| オフラインと明示更新（指示書 2・3、判断 6、10、14） | `src/sw/`（`protocol.ts`・`prepare.ts`・`serve.ts`・`cleanup.ts`・`worker.ts`）、`build/service-worker.ts`、`src/ports/offline.ts`、`src/infra/browser/offline.ts`、`src/application/release.ts`（ReleaseManager）、`ui/OfflineStatus.vue`・`OfflineNotice.vue`・`UpdateNotice.vue`・`UpdateImpact.vue` | `offline.spec.ts`、`sw-prepare.test.ts`、`sw-serve.test.ts`、`release.test.ts`、`release-pin.test.ts`、`offline-client.test.ts`、`layers.test.ts` |
| 利用統計（指示書 4、判断 11、12） | `domain/analytics.ts`、`ports/analytics.ts`、`infra/browser/analytics.ts`、`application/analytics.ts`、`ui/AnalyticsConsent.vue`・`AnalyticsSetting.vue`、`functions/api/visit.ts`、`ga4_probe/` | `analytics.spec.ts`、`analytics.test.ts`、`visit-function.test.ts` |
| 導入支援（指示書 6、判断 8） | `src/examples/`（例題 10 file、`self-test.json`）、`tools/examples.py`、`domain/examples.ts`、`application/examples.ts`・`onboarding.ts`・`self-test.ts`、`Coordinator.withEngine`、`ui/TryExample.vue`・`GettingStarted.vue`・`SelfTestPanel.vue`、`ports/local-settings.ts` | `onboarding.spec.ts`、`onboarding.test.ts`、`self-test.test.ts`、`self-test-engine.test.ts` |
| V-PRIV と TD-11（指示書 5、判断 7、13、15） | `tests/e2e/support/privacy-scan.ts`・`privacy.ts` | `privacy.spec.ts`、`privacy-scan.test.ts` |
| プレビューの手順と確認（指示書 7、判断 9） | `docs/web/deploy_cloudflare.md`、`docs/evidence/losat_web_w7/check_preview/` | WP-C の dry run、ゲートの dry run |
| 計測と画面の記録（判断 16） | `delivery-measure.spec.ts`（新規）、`results-measure.spec.ts`（「Offline ready」の後に始める）、`screens.spec.ts`（状態 41〜60）、[`run_gate.sh`](run_gate.sh) | ゲートの実行 |

## 作ったもの

- `src/domain/`：`release-manifest.ts`、`site-index.ts`、`site-file.ts`、`analytics.ts`、`examples.ts`。
- `src/application/`：`release.ts`（ReleaseManager）、`analytics.ts`、`examples.ts`、`onboarding.ts`、`self-test.ts`。`coordinator.ts`（`withEngine`）、`draft.ts`（`loadExample`）、`composition.ts` は変えた。
- `src/ports/`・`src/infra/`：`ports/offline.ts`・`site.ts`・`analytics.ts`・`examples.ts`・`local-settings.ts`、`infra/browser/offline.ts`・`site.ts`・`analytics.ts`・`examples.ts`・`local-settings.ts`。
- `src/sw/`・`src/boot/`：Service Worker とその準備・答え・掃除、起動ページ。ESLint の層の規則に入れた。
- `src/ui/`：`OfflineStatus.vue`、`OfflineNotice.vue`、`UpdateNotice.vue`、`UpdateImpact.vue`、`AnalyticsConsent.vue`、`AnalyticsSetting.vue`、`TryExample.vue`、`GettingStarted.vue`、`SelfTestPanel.vue`。`AppView.vue`、`SearchPanel.vue`、`styles.css` を変えた。
- `build/`・`functions/`・`public/`：release の build、`_headers`、ローカルの server、組み立て、Service Worker の build、Pages Function `functions/api/visit.ts`。
- 文書：`docs/web/release_manifest.md`・`analytics.md`・`examples.md`・`deploy_cloudflare.md`（新規）、`ncbi_ui_mapping.md` の追記、[`run_gate.sh`](run_gate.sh)、[`check_preview/`](check_preview/README.md)、[`ga4_probe/`](ga4_probe/RESULTS.md)、この記録、S16b の指示書、S17 の指示書への引き継ぎ。

## 判断

判断 1〜20 はすべて推奨案で進めたもの（Owner-delegated、2026-09-29 の常設の指示、2026-10-07 に再掲。2026-10-11 JST に取った）。保守者が予約した判断は含まない。出典は取りまとめ役の `DECISIONS.md`（タスクフォルダ `/home/kawato/losat-baselines/s16-w7-20261011/`）と各担当の記録（`logs/wp-*/progress.md`、`logs/fix*/progress.md`）。

1. 契約（取りまとめ役）：release ID は `LOSAT_WEB_RELEASE_ID` があればそれ、無ければ `<commit 12>-<エンジンの digest 8 または "fake">[-dirty]`（`build/release.ts`、`__LOSAT_RELEASE_ID__`）。形は `^[0-9a-z][0-9a-z.-]{0,63}$`。ビルドのバイトを決めるものから決まり、path の 1 区間にも GA4 のパラメータにも使え、ビルドの前に分かる（Vite の `base` に入れられる）。
2. 契約：サイトの索引 `/releases.json`（format 1：`latest`、`releases[]`、任意の `analytics.endpoint`（`/api/` の下））はサイト単位で、ビルドか組み立ての時に書き、release には入れない。未知の field は無視する。同じ release の bytes が、計測のあるサイトにもないサイトにも使え、1 度の読みが更新の確認と計測の切り替えに足りる。
3. 契約：固定した release は `localStorage['losat-web:release']`。アプリのすべての key は `losat-web:` で始める（計画 §5.9：アプリの設定は研究データと別）。
4. 作業の分け方：WP-R・U・G・O・P・C・GP。アプリ側は `feature/losat-web-gui` に merge しない（S15 のとおりエンジン側が行う）。
5. WP-R：FakeEngine の manifest は `engine: null`。起動ページは 2 つ目の Vite の build（domain を import）。shell は `/`・`/boot.js`・`/boot.css`（`404.html` は Pages が `/404` へ redirect するので入れない）。起動ページが manifest と一緒に `/releases.json` を読み、固定が開始できないときだけその答えを使う。404・410・JSON でない・不正な manifest は「配っていない」で、「Use release <latest>」を出す。自前のローカル server（`build/site-server.ts`、依存無し）が E2E と `npm run preview` に使われ、Pages の `.html` の redirect と既定の `Cache-Control` を真似る。`engines.node >=22.18`（組み込みの TypeScript）。組み立ては、環境の計測 endpoint が優先、無ければ新しいサイトのもの。最新は `--latest` が無ければ新しい build。
6. WP-U 第 1 段階：検査済みの file はメモリに持ち、全部合うまで何も cache しない（manifest が最後）。失敗したら shell を戻す。release の file は既定の cache mode で 1 度だけ `reload` で読み直す。shell は `no-cache`、索引は `no-store`。経路は URL だけで決める。shell は network 優先（3 s、5xx は保存した写し）。`clients.claim()` は activate で、`skipWaiting` は利用者が切り替えるときだけ。オフラインと更新は `/r/<id>/` で動くときだけ。準備は load の後 2 s と idle（5 s の timeout）、Run が準備中でないとき。Service Worker 無しでは検査なしの切り替えを注付きで出す。「Not now」は release ごと・tab ごと。固定は各確認で読み直す（取り損ねた broadcast を補う）。leave guard が再読み込みを止めたときは「Reload now」。V-OFF はサーバーを止めてオフラインにする（Playwright の `setOffline` は Chromium でしか確かでない）。
7. WP-P 第 1 段階：捕獲の server は ephemeral な port。CSP に `report-uri /__privacy/csp-report`（試験だけ。worker の違反を記録する唯一の方法）を足す。サイトの file への要求は GET・HEAD だけ（`/api/visit` は許可の hook、`/favicon.ico` は Firefox が自分で尋ねる）。試験文字列は大文字小文字を問わず、パーセント符号化、行をまたぐ残基、逆相補の形を探す。TD-11 は配る全 file を読み、短い・よくあるユーザー名は `/`・`\`・`~` の隣だけ。拒否・失敗の経路はエンジン入りのビルドだけ。
8. WP-O：5 つの例題は `LOSAT/tests/fasta` の DnaA とその遺伝子（2,000 nt の領域、タンパク質は全長、subject 3 件で identity が段階的）。「Try an example」は program の一文の後。読み込みはその program の欄を消し、スレッドを Auto に、Job Title を空にし、入力を置き換える前に尋ねる（修正 1 回目から設定・スレッド・record の扱い・Job Title も尋ね、置き換えるものを言う）。自己試験は 1 スレッドと Auto が 2 以上になる大きさの検索の 2 回で、2 つの module が動く。`Coordinator.withEngine` が自己試験の間キューを止め、動作中・待ちの Run があれば拒む。BLASTX の期待値は同梱し、SX まで skip。自己試験は Getting started の閉じた節「Test this browser (optional)」。`tools/examples.py`（cut / check / expected）。`ports/local-settings.ts` は `losat-web:` 以外の key を拒む。画面 01〜40 は Getting started を閉じて撮る。15 の期待値はネイティブの LOSAT（`141615ec`）と NCBI BLAST+ 2.17.0 の bytes が同じ。
9. WP-C：Cloudflare Pages の Direct Upload（`npx wrangler@4.149.0 pages deploy dist` を `web/app` で。`functions/` を上げるため）。Pages の Git build は使わない（build image に Rust が無く、20 分の上限）。サイトは `LOSAT_WEB_ANALYTICS_ENDPOINT=/api/visit` で build し、同じ `dist/` をそのまま本番に出す。GA4 の secret は Production だけ（プレビューの Function は転送しない）。組み立ての「前のサイト」は、保守者がリポジトリの外に残した最後の本番の `dist/`。確認は `LOSAT_WEB_PREVIEW_HEADERS`（デプロイした `_headers`）を取り、COOP・COEP・CORP・CSP が 200・304・404 で規則と等しいことを要求し、他の header（`Cache-Control` など）は 304・404 では `notes`。配った bytes を manifest と shell に照らし、注入された script を探し、索引にある全 release の manifest を要求する。確認の folder に `node_modules` の symlink を置いた（ESM の解決のため）。
10. WP-U 第 2 段階：アプリの `connect()` も `register('/sw.js')` を呼ぶ（同じ URL・option、追加の取得なし）。`sw.js` はサイト単位で shell には入れない。V-OFF は worker ごとに voff-a / voff-b を build し（B は `sw.js` を変える）、port は E2E + 100 × (slot + 1)、オフラインは server を閉じて接続を切る（Chromium・Firefox は `setOffline` も）。`release.spec.ts` の reload-304 の試験は Service Worker を塞ぐ。起動ページは `src/sw/protocol.ts` の `SW_URL`・`SW_SCOPE` を import する。6b は固定の無い、保存した起動ページを使う（本当のオフラインの初回訪問にはページが無い）。WebKit は `navigator.onLine` が true のまま（起動ページの message の種類は `failed`・`no-index`）。merge の順：見出しは「Getting started」の次にオフラインの状態、UpdateNotice は GettingStarted と ResumeNotice の上。
11. WP-G：**GA4 の方式 (b)**。アプリは同意の後に固定のイベント（`{"v":1,"release","client"}`）を同一 origin の Pages Function `/api/visit` へ送り、Function が本番の secret があるときだけ `app_visit`（`release_id` 付き）を Measurement Protocol へ転送する。どの文書にも Google のスクリプトは無く、CSP と COEP は変わらない。probe（50 件、3 ブラウザ）：(b) は Chromium・Firefox・WebKit で動く。(a) は `frame-src` が要り、a1 はどこでも拒まれ、a3 は WebKit で動かず、a2 は WebKit で第三者 cookie を置かず、a4 は Chromium だけで ID を失う。指示書は推奨を保守者に示すよう求めるので、常設の指示で採ったうえで確認を求める。
12. WP-G：endpoint は page が送る正確な文だけを受ける（`parseVisitBody` を共有）。Function は自分の header を付ける（Pages は Function に `_headers` を適用しない）。同意の知らせはページの流れの中でタブの上、フッターの「Change」で開き直す。storage が失敗したら何も送らず、フッターがそう言い、次の訪問でまた尋ねる。key は `losat-web:analytics`。merge：起動時の索引の読みは計測の切り替えと更新の確認で共有（`sharedSiteReads`、60 s）。release ID の無い build は計測しない。
13. WP-P 第 2 段階：捕獲は `build/site-server.ts` でサイトの view を配る（`onRequest` の log。別の「デプロイ」は view を別の directory に向ける、origin は同じ）。小さな front が `/api/` に Function のように 204 を返す。計測するサイトは view の `releases.json` に `analytics.endpoint` を足したもの。各経路は開始のあと「Offline ready」を待つ。`allHeaders()` は最大 2 s 待つ（server が全 header を記録する）。Firefox 自身の `/sw.js` の更新の取得は許す。TD-11 は shell、`releases.json`、`sw.js`、404、各 manifest の file を並べ、これが `_headers` を除くビルドの全 file であることを確かめる。
14. 修正 1 回目（コードレビュー 1 回目）：release ID は commit の先頭 12 桁 + エンジンの tag、`app.build` / `__LOSAT_APP_BUILD__` は commit の先頭 8 桁（clone の深さによらない。セッションファイルは 8 桁を記録する）。`-dirty` は `build/inputs.ts` に 1 度だけ列挙した入力に対して見る。commit の無い build は `LOSAT_WEB_RELEASE_ID` が要る。2 つの build を byte で比べる（`release-bytes.test.ts`）。影響の一覧は切り替えの隣でもう一度、数え直して出し、「Switch to release …」とする。取得から切り替えか取消まで、tab は `losat-web:release:<latest>` を持ち、その間に消された release は「removed」で「Download and check again」（切り替えの中では自動でもう一度取得しない）。失敗・時間切れの登録はオフラインの失敗で「Try again」（何も返さない登録は Service Worker 無し）。「Offline ready」は読み戻した固定がこの release を指すときだけ（違えば「stored」）。見出しのオフラインの状態は固定の 1 行（30em、短い文、release ID は「…」で切り title に全文）。詳細と「Try again」は見出しの下の `OfflineNotice.vue`（まれな失敗でページが 1 度動く）。404・410・JSON でない manifest は `cache: 'reload'` で読み直す。V-OFF のサイトは `/r/*` を no-store で配る（`releasesNoStore`）ので、ケース 2 が Service Worker が答えることを証明する。`vite.config.ts` は関数の設定になった。
15. 修正 3 回目（コードレビュー 2 回目）：V-PRIV は、サイトへの query string と、アドレスの `?`・`#` を指摘とし、試験文字列を全 URL・header・body でパーセント復号、base64（3 つのずらし、2 つの字母）、空白と JSON のエスケープを除いた形で探す。query の例外表は持たない。捕獲の自己試験は、入れ子の worker と Service Worker が他の origin を試すのを仕込む（WebKit は worker から CSP の report を送らない。証拠は失敗した fetch と空の他の origin の log）。ゲートは両方の release ID を `BUILD_INPUTS` に対して最初に確かめ、`-dirty` で止まり、E2E と各計測の後に `releases.json` と manifest を写し、`NODE_ENV` を外す。`release-bytes.test.ts` は、`BUILD_INPUTS` だけを持つ一時的な 2 つ目の checkout（と `node_modules` の写し）からも build する。
16. 修正 2 回目（画面レビュー 1 回目、コードレビュー 2 回目 L1）：保管場所の見積りは release を保存した後、掃除の後、更新の取得の後に読み直す。失敗の理由は、どの file がなぜかを言う（「a file did not match its checksum in the release's file list (<path>)」）、エラーの枠（`--error-bg`）に出す。リンクのボタンは 2 px の padding と同じ負の margin で当たりを 24 px 以上にする（動かさない）。無効な BLASTX の「Try an example」は muted で下線なし。「SHA-256」は 1 行に保つ。読み込んだ Run の知らせは「commit <hex>」。起動ページは「Use release <id> (stored)」と「Try again」を出し、自分では切り替えない。`register()` が拒まれる（Chromium は `NotSupportedError`、Firefox は `SecurityError`）のは Service Worker 無しで、知らせは出さない。`vite build` は `NODE_ENV` が production 以外なら拒む（`buildEnvironmentProblem`）。画面 57〜60（更新の取得中、チェックサムの失敗、待ち、測るサイトの初回）。`ncbi_ui_mapping.md` の §5。
17. 修正 4 回目（コードレビュー 3 回目）：`vite build` は、Vite が production の build を production の mode で解決したときしか進まない（`configResolved` の `productionBuildOnly`。`--mode` が production 以外でも止まる：同じ release ID に別の `.env` が入るため）。ゲートは `VITE_USER_NODE_ENV` を外し、2 つ目の checkout からのバイト同一の試験が走ることを要求する（`LOSAT_WEB_SECOND_CHECKOUT=required`）。`register()` の `SecurityError`・`NotSupportedError` は、制御する worker が無く `getRegistration('/')` も拒まれるときだけ「Service Worker 無し」（初回訪問で HTML の `/sw.js` が返るのは「Try again」付きの失敗のまま）。出せる最新の release が無い「not served」の起動ページは「Try again」と保存された release を出す。「Use release <id> (stored)」は click の時点で保存された release を読み直す（Web Lock は取らない。固定はすぐ書く）。決定性の文は、1 台の機械・1 つの Node・どの checkout のパスでも、に絞る。
18. 修正 5 回目：`privacy.spec.ts` の Chromium の flaky は試験で直した。Chromium は 1 つのページで、最初の download から 1 s 以内の 11 個目の download を黙って捨てる（測定：10 個を続けて出し、11 個目を +816 ms で出すと捨てられ、+1316 ms で残った。ユーザーの操作でも解けない。Firefox に制限は無い）。アプリは 1 回の click に 1 つの download しか作らないので、人は届かない。試験は 1 つのページの download を 150 ms あけて行う（`tests/e2e/support/downloads.ts`、`performance.now()`）。`privacy.spec.ts` と `session.spec.ts`（1 つのページから 16 回保存）に使う。
19. 取りまとめ役：画面レビュー 2 回目（合格）の Low 4 件（release ID の全文が hover だけ、「Offline:」の読み違い、他の tab の知らせの数、起動ページ 55 の文）は、アプリの変更なので新しいゲートの実行が要る。コードレビュー 4 回目の Low 1 件（文書）は `2b553a4b` で直した。四つの Low は S16b に回す（S16b はプレビューの確認の後にゲートを再実行する）。S16 は保守者のデプロイの前には終われず、1 回のゲートの実行が両方を覆う。
20. 取りまとめ役：Chromium の 10 万 query の JSON の書き出しの最長のフレームの間隔（ゲート 3 は 250 ms、W6 のゲート 3 は 117 ms）は W7 の退行ではない：W6（`e45fbd09`）と W7（`35970fa1`）を交互に 3 回ずつ測った対の A/B で、183 / 133 / 250 ms（中央値 183）対 200 / 167 / 133 ms（中央値 167）。100 ms を超える間隔はこの機械で W6 からあり、残件のまま。ゲート 2 のプレビューの確認の trace（3 件の失敗から `trace.zip` が 62 MB）は、リポジトリの外のタスクフォルダ（`logs/gate2-preview-test-results/`）に置き、実行の folder に入れない。（出典：`logs/ab-json.log`）

## 実測

実測は `tests/e2e/delivery-measure.spec.ts`（配信。新規）と `results-measure.spec.ts`（`--workers=1`、機械は他の処理と共有）。1 回の暖機と 3 回の中央値と範囲。結果の計測は W7 から「Offline ready」の後に、Getting started を閉じて始める。値はすべてゲート 3（`35970fa1`、最後の木）の [`run-20261011T080950Z/measure/`](run-20261011T080950Z/measure/)（release `35970fa18885-e2fb7128`、保存は 17 file・8,188,494 byte）。

**配信**。ms は訪問の始まり（ナビゲーションの時刻）から。`app` は検索フォームが画面にある、`ready` は見出しが「Offline ready」と言う。中央値（範囲）：

| ブラウザ | 訪問 | app | ready | 送られた byte（準備の間の分） | 自己試験 |
|---|---|---|---|---|---|
| Chromium | 初回 | 65（63〜65） | 2,189（2,185〜2,192） | 8,201,096（3,978,574） | 226 ms（225〜230） |
| Chromium | 2 回目 | 31（31〜34） | 2,049（2,045〜2,050） | 1,836（612） | |
| Chromium | オフライン | 15（15〜15） | 2,020（2,020〜2,023） | 0 | |
| Firefox | 初回 | 116（111〜118） | 2,262（2,257〜2,262） | 8,201,520（3,978,574） | 1,690 ms（1,675〜1,754） |
| Firefox | 2 回目 | 76（73〜77） | 2,103（2,102〜2,104） | 1,836（612） | |
| Firefox | オフライン | 67（65〜69） | 2,076（2,075〜2,078） | 0 | |
| WebKit | 初回 | 83（82〜87） | 2,166（2,161〜2,167） | 12,411,628（8,189,106） | 242 ms（213〜253） |
| WebKit | 2 回目 | 45（44〜45） | 2,060（2,058〜2,062） | 33,994（11,308） | |
| WebKit | オフライン | 34（33〜39） | 2,044（2,043〜2,045） | 0 | |

- 「Offline ready」までは、どの訪問でも約 2.0〜2.3 s。準備は load の後 2 s と idle を待つ設計（判断 6）で、初回の準備（ダウンロード・検査・保存）そのものは Chromium 119 ms、Firefox 142 ms、WebKit 75 ms。WebKit は準備の間にほぼ全 file（8.19 MB）をもう一度取る（Chromium・Firefox は 3.98 MB。ページが読んだ file は HTTP cache から読む）。2 回目は 1.8〜34 kB（索引と shell の確認）、オフラインは 0 byte。
- 検索フォームが出るまで（`app`）：初回 65〜116 ms、2 回目 31〜76 ms、オフライン 15〜67 ms。engine の確認（`checked`）は Chromium 445、Firefox 573、WebKit 502 ms（初回）。
- 例題（Try an example → Run LOSAT → Open results、`firstResults`、中央値）：Chromium 35〜51 ms、Firefox 168〜247 ms、WebKit 40〜140 ms（BLASTN から TBLASTX まで）。自己試験の 1 case は、Chromium 16〜33 ms、Firefox 182〜255 ms、WebKit 9〜36 ms。
- ゲートの前の WP-GP の値（`d607a2f5`、`logs/wp-gp/measure/`）では初回の「Offline ready」が 2.2〜2.5 s、自己試験が 0.6 / 3.9 / 1.2 s で、ゲート 3 のほうが短い（機械の負荷の違い。追っていない）。

**結果と書き出し、10 万 query（149,828 HSP）**（中央値。ms は click から画面が要約か message を出す（outfmt 6 は download）まで。括弧は最長のフレームの間隔 ms。W6 のゲート 3 の値は「／」の後）：

| ブラウザ | Run を開く / 最初の HSP の詳細 | CSV | JSON | レポート | outfmt 6 | セッション保存 | セッションを開く |
|---|---|---|---|---|---|---|---|
| Chromium | 266 / 266（258 / 259） | 524（50）／409 | 1,728（250）／1,794 | 1,335（17）／1,374 | 18（17）／22 | 2,334（17）／2,302 | 1,570（43）／1,680 |
| Firefox | 473 / 472（554 / 473） | 383（83）／372 | 1,643（17）／1,751 | 1,292（17）／1,458 | 30（17）／33 | 1,771（33）／1,867 | 1,872（33）／1,944 |
| WebKit | 357 / 409（410 / 461） | 580（134）／710 | 1,641（117）／1,759 | 1,352（118）／1,419 | 18（16）／20 | 1,948（32）／2,064 | 1,345（50）／1,443 |

- バイト数は W6 のゲート 3 と同じ（CSV 16,242,367、JSON 127,249,939、レポート 176,469,035、outfmt 6 8,966,787。セッションは保存の時刻の分だけ違う）。W7 の追加（Service Worker、見出しの状態、起動ページ）で、Run を開く・書き出す・セッションの時間は W6 のゲート 3 と同じか速い。W6 より遅いのは Chromium の Run を開く（+8 ms）、CSV（+115 ms）、セッション保存（+32 ms）だけ。
- 100 ms を超えたフレームの間隔：Chromium の JSON の **250 ms**（W6 のゲート 3 は 117 ms）、WebKit の書き出しの開始の 117〜134 ms（既知）。他は 83 ms 以下。Chromium の JSON の 250 ms は W7 の退行ではない：W6（`e45fbd09`）と W7（`35970fa1`）を交互に 3 回ずつ測った対の A/B で、最長のフレームの間隔は W6 が 183 / 133 / 250 ms（中央値 183）、W7 が 200 / 167 / 133 ms（中央値 167）だった（判断 20、`logs/ab-json.log`）。100 ms を超える間隔はこの機械で W6 からある。
- 1 組の場面：1,500・3,000 写しは Chromium と Firefox で通り、WebKit は 1,500 写しだけ（3,000 は保存の上限、4,000 は 3 ブラウザでエンジンのメモリ不足 `memory allocation of 49861 bytes failed`。W6 のゲート 3 と同じ見込みどおり）。3,000 写し（5,993 HSP）のドットプロットを出す・拡大・縮小は Chromium 18 / 4 / 4 ms、Firefox 41 / 28 / 28 ms。3,000 写しのセッション（約 87 MB）は Chromium で保存 8,368 ms（間隔 17 ms）・開く 7,183 ms（33 ms）、Firefox で保存 8,233 ms（17 ms）・開く 10,476 ms（50 ms）（W6 のゲート 3 の Firefox は保存 14,863 ms・間隔 633 ms、開く 16,265 ms）。WebKit の 1,500 写し（2,995 HSP）は保存 18,448,645 byte・1,291 ms、開く 1,248 ms。
- トレイ（10 万 query の 200 候補）：Chromium の印・追加・タブを開く・scroll・並べ替え・抽出 100・全部を外すは 11 / 11 / 5 / 11 / 12 / 43 / 8 ms、WebKit の抽出は 88 ms（W6 は 91 ms）。
- 準備が利用者の作業と重なる場面の数は無い（コードレビュー 2 回目 I2）。計測の各場面は「Offline ready」の後に始める。初回の訪問ですぐ検索すると、worker が threaded module（約 4 MB）を取り、約 8 MB の SHA-256 を計算し、Cache Storage に書くのと検索が重なる。`noRunPreparing` が延ばすのは準備中の Run だけ。

## 独立レビュー（コード）

agent `losat-reviewer`（code review の役）。記録はタスクフォルダの `reviews/code-review.md`、`code-review-2.md`、`code-review-3.md`、`code-review-4.md`。

**1 回目**（`e45fbd09..acd0f00b` の `web/app`・`docs/web`・`docs/evidence/losat_web_w7`・`.github/workflows/web.yml`。V-PRIV は含まない）：High 無し、Medium 3 件、Low 5 件、情報 6 件。Medium と Low は修正 1 回目（`1a089486..a9507eac`）と修正 2 回目で直した。

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 更新の影響は取得の前にしか出ず、切り替えの時点では出ない。取得の間に完了した Run や候補が、警告なしに失われる | `a19a7e11`：切り替えの隣で、数え直した一覧を出す（判断 14）。V-OFF 3 が、確認の後に完了した Run を数えることを確かめる |
| M2 | 取得して確かめた release を、別の tab の掃除が消し得る。切り替えは保存されているものとして進む | `3f74058b`：取得から切り替えまで shared lock を持ち、worker は待ちの lock も数える。切り替えは保存を確かめ、消えていれば「removed」。V-OFF 3b |
| M3 | release ID がビルドのバイトを決めない（git の短縮の桁が clone で違う、`-dirty` が `web/app` だけ、git 無し） | `232800d6`：commit の固定の桁、入力の列挙（`build/inputs.ts`）、commit 無しは拒否、`release-bytes.test.ts`（判断 14） |
| L1 | V-OFF のケース 2 が HTTP cache で通り得る | `d4736ac1`：`/r/*` を no-store で配る |
| L2 | `/r/` の 404 が immutable として覚えられ、読み直されない | `daca999d`：`cache: 'reload'` で 1 度読み直す |
| L3 | 「Try an example」が設定と Job Title を黙って消す | `71150eef`：尋ねる |
| L4 | 登録の失敗で、Service Worker のあるブラウザに検査なしの切り替えを出す | `405ad0ca`・`3db6c216` |
| L5 | `src/sw`・`src/boot` が層の規則の外 | `1a089486` |
| M4（WP-P の全体の実行と画面レビュー 1 回目 L1） | 見出しの状態の大きさが変わり、WebKit が click を落とす | `895d2949`・`a9507eac` |

情報：I3（文書の不一致）は `c283324a`、I5（「Offline ready」の前に固定を読み戻す）は `788c23ac`、I2 の前半（`Origin` の確認）は `08cd5359` で直した。**直していない（残件）**：I1（準備が shell を保存している間に打ち切られると、古い shell と新しい shell が混ざり得る。3 つの `put` の間の窓。`/`・`/boot.js`・`/boot.css` が release をまたいで互いに互換であるかぎり安全。新しい cache に書いて入れ替える手もある）、I2 の後半（wrangler が `../../src/domain/analytics` を bundle できるか。`npx wrangler@4.149.0 pages functions build --outdir <tmp>` で、アカウント無しに確かめられる。S16b）、I4（自己試験は取り消せず、1 つの case が戻らないと利用者のキューが始まらない。再読み込みで直る）、I6（V-OFF 3 は document の読みとサーバーに届く要求でしか「B だけ」を確かめない。構成が混ざりを不可能にしている）。

**2 回目**（`acd0f00b..a9507eac`：V-PRIV、ゲートの準備、修正 1 回目）：High・Medium 無し、Low 5 件、情報 5 件。V-PRIV が §5.10 の全経路を覆うこと、修正 1 回目が M1〜M4・L1〜L3・L5・I3・I5 に正しく効いていること、ゲートの script と計測が言うとおりであることを確かめた。

| 指摘 | 内容 | 対応 |
|---|---|---|
| L1 | `register()` の `SecurityError` が毎回の「Offline use failed」と知らせになる。「Try again」は成功しない | `19ede7f0`：`SecurityError`・`NotSupportedError` は Service Worker 無し（判断 16） |
| L2 | V-PRIV は query string と、base64 などの符号化を見ない | `c458e142`：query string と `?`・`#` は指摘、符号化も探す（判断 15） |
| L3 | Firefox・WebKit で、入れ子の worker と Service Worker の他の origin への試みが見えない | `9cc355e6`：仕込んで、3 ブラウザで届かないことを示す |
| L4 | ゲートが `-dirty` を E2E の後でしか見ない。`measure` の run は release を記録しない | `ac6ef3dc`：最初に確かめ、記録する |
| L5 | バイト同一の試験が同じ clone・同じパスで 2 回 build するだけ | `0c5956dc`：別のパスの 2 つ目の checkout からも build。ここで `NODE_ENV` の問題を見つけ、`aeca5d93`・`cc5ca34f` |

**直していない（残件）**：I1（V-OFF 3b が、2 つ目の tab の掃除が走る前に通り得る。`cleaned` の返事を待つか、別の release の lock を取って除かれるのを待てば決定的になる）、I2（準備が利用者の作業と重なる場面の数）、I3（最後の `findings()` の後の要求は log されるが assert されない。close のときの beacon）、I4（更新の知らせと同意の知らせが、読み込みごとにページを 1 度動かす。設計が受け入れている）。I5（ゲートの単体の一覧の注釈）は `ac6ef3dc` で直した。

**3 回目**（`a9507eac..ba7fa464`：修正 2 回目と 3 回目。ゲート 1 の実行中に、その run folder とエラーの文脈 1 つだけを読んだ）：High 無し、Medium 1 件、Low 3 件、情報 3 件。修正 2 回目の L2・L5・L6・I2〜I6、コードレビュー 2 回目の L1、修正 3 回目の V-PRIV とゲートの script と `release-bytes.test.ts` は、指摘どおりに効き、退行が無いことを確かめた。

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | V-OFF のケース 1 の保管場所の下限（保存した byte の 95 % 以上）が、エンジン入りのビルドの Firefox で偽（Firefox は Cache API の body を圧縮して保存する）。ゲート 1 を落とした | `4c5bf901`：見出しの値をブラウザ自身の `navigator.storage.estimate()` と比べる。2 つのビルド × 3 ブラウザで通した |
| L1 | `NODE_ENV` の規則を、Vite 自身の `.env`（と `VITE_USER_NODE_ENV`）が config の関数の後に回避できる。別の値の `NODE_ENV` で開発 build が clean な ID で出来得る | `3824aa9c`：`configResolved` で、Vite が production を解決したときしか進まない（`--mode` が production 以外でも止まる）。ゲートは `VITE_USER_NODE_ENV` を外す（判断 17） |
| L2 | `register()` の `SecurityError` を全部「Service Worker 無し」とすると、JavaScript でない `/sw.js`（`200 text/html`）も同じに扱い、検査なしの切り替えを出す | `d3812c2e`：制御する worker が無く `getRegistration()` も拒まれるときだけ「無し」。そうでなければ「Try again」付きの失敗 |
| L3 | 起動ページの「not served」で最新の release を出せないとき（オフライン、索引の時間切れ）、「Try again」も保存された release も無い行き止まり | `0d7f439e`：この場合も「Try again」と保存された release を出す |

情報：I1（「Use release <id> (stored)」は、message を出した時点の一覧で固定する）は `72ae82c9`（click の時点で読み直し、消えていれば言って残りを出す）、I3（決定性の文が、1 台の機械・1 つの Node より広く言う。ゲートが 2 つ目の checkout の skip を見逃す）は `b2a6e066` で直した。**直していない（残件）**：I2（V-PRIV が、request header は試験文字列の形（平文・パーセント符号化・base64）しか探さない。hex・base32・圧縮した custom header は見えない。アプリは custom header を付けない）。

**4 回目**（`4c5bf901`、`3824aa9c..b2a6e066`、`9cb4e61e`、`35970fa1`：修正 4・5 回目とゲートの 2 件の試験の修正。worktree は変更せず、npm/npx も probe も走らせず、commit と Vite の `resolveConfig`、修正 4 回目の log、ゲート 3 の run folder を読んだ）：High・Medium 無し、Low 1 件、情報 3 件。起動ページが自分では切り替えないこと、選択肢が manifest の保存された release だけを出し click で読み直すこと、本物の `register()` の失敗が失敗のままであること、production の build の規則が Vite の最終の判断を読み開発 server と Vitest に触れないこと、試験の変更が言うとおりのことを試験していること（download の間隔は失敗を隠さず、保管場所の確認は今の大きさで広い余裕を持つ）を確かめた。

| 指摘 | 内容 | 対応 |
|---|---|---|
| L1 | `release_manifest.md` の起動ページの記述が、`0d7f439e` と `72ae82c9` の前のまま（最新が無い「not served」の選択肢、click の時点の再確認、`boot-stored-gone`、対応表の起動ページの行） | `2b553a4b`：文書を直した（アプリのバイトは変わらない） |

情報は上の「残件と注意」に記した（I1〜I3）。

## 画面レビュー

agent `losat-reviewer`（screen review の役）が、画面の記録を W6 の最後の記録、対応表、W4 の規則と並べて見た（記録はタスクフォルダの `reviews/screen-review.md`。行ごとの差の script と切り出しは `tmp/review-screens/`）。

**1 回目：合格**（High・Medium 無し）。`70442212` の記録 366 枚（3 ブラウザ × デスクトップ 1280 px・電話 390 px、状態 01〜56。SHA-256 は 366 枚とも `logs/wp-gp/screens-3browsers.sha256` と一致）。基準は W6 の最後の記録（`s15-gate-screens`、`08e4557b`、状態 01〜40 と 30b）。01〜40 は、見出しの「Getting started」と「Release … Offline ready」、検索の状態の program の一文の後の「Try an example」、保管場所の見積り、実行ごとの値の他は変わらなかった。見出しは、デスクトップで「Offline ready」のとき 0 px、検索の状態で +12 px、電話で +48〜50 px（検索の状態で +81〜82 px）。新しい状態 41〜56 は仕様と W4 の規則に従い、電話の幅 390 px に収まり、3 ブラウザで同じに描かれ（WebKit の 55 だけ `navigator.onLine` が true のまま）、同意（52・53）は `analytics.md` のとおり、更新（47〜51）は影響と選択を言う。NCBI のページに寄せた置き場所で、ロゴやブランドは無い。

| 指摘 | 内容 | 対応 |
|---|---|---|
| L1 | オフラインの状態の大きさ・位置が状態で変わる（電話で 1〜3 行、括弧の中で折れる、デスクトップで右端から左の別の行へ） | `a9507eac`：見出しの 1 行（30em、短い文、ID は「…」）、詳細と「Try again」は見出しの下の知らせ（判断 14）。46・46b・56 を撮り直し |
| L2 | 保管場所の見積りが、保存した release を含まない（「Offline ready」の後も 10.9 kB など） | `0877440d` |
| L3 | 49・51 で、再読み込みが何を失うかを言い直さない | `a19a7e11`・`a9507eac`：切り替えの隣と他の tab の知らせに |
| L4 | 取得中・失敗・待ち・測るサイトの初回の記録が無い | `163a22fb`：状態 57〜60 |
| L5 | 50：失敗の文が不明瞭、青の枠の中の赤の文字 | `af821677`：「does not match the release's file list」、エラーの枠 |
| L6 | リンクのボタンが 21〜22 px（W4 の規則は 24 px） | `3b457430` |
| 参考 I1〜I6 | I1 電話の見出しの増加（+49 px）は受け入れる。I2 無効な BLASTX の例（`3b457430`）、I3 「SHA-/256」の折れと「Example loaded」の余白（`3b457430`）、I4 対応表の W7 の行（`79f724ad`）、I5 build の名前の二重（`9e9ee11e`）、I6 起動ページの 55・55b（`cf017d9d`） | |

**2 回目：合格**（High・Medium 無し）。ゲート 3（`35970fa1`）の記録 396 枚（3 ブラウザ × 2 サイズ × 66 枚、`$BUILD_ROOT/s16-gate3-screens/`。SHA-256 は 396 枚とも実行記録の `screens.sha256` と一致）。基準は W6 のゲート 3 の記録（01〜40 と 30b）と画面レビュー 1 回目の記録（01〜56）。1 回目の指摘は、L1〜L5 が直り、L6 はコードで直した（24 px の当たり。記録には出ない）ことを確かめた：オフラインの状態はどの状態でも 1 行・同じ位置・同じ大きさで、内容が動くのは失敗の知らせのときだけ（デスクトップで 101 px、電話で 176 px）、保管場所の見積りは「Offline ready」の後に release を数える（Chromium 8.2 MB、Firefox 5.4 MB）。01〜40 と 30b は、見出し、「Try an example」、保管場所、実行ごとの値、読み込んだ Run の「commit 35970fa1」の他は W6 と変わらない。見出しの増加は 1 回目と同じ（デスクトップ 0 px、検索の状態で +12 px。電話で +48〜49 px、検索の状態で +81 px）。41〜60 は 3 ブラウザで同じに描かれ、電話の幅 390 px に収まり、同意は正直で、更新は影響と選択を言う。

| 指摘 | 内容 | 対応 |
|---|---|---|
| L1 | 見出しに収まらない release ID は hover の title でしか全文が読めない（タッチの画面は title を出さない。デスクトップでも 30em の箱が 46 で既定の ID を切る） | S16b（判断 19）：Getting started に release ID の 1 行（「This is LOSAT Web release <id>」）。任意で箱を約 32em に |
| L2 | 「Offline: 4.3 MB of 8.2 MB…」「Offline: waiting」が、アプリがオフラインだと読める | S16b：「Offline copy: …」など、同じ短さで |
| L3 | 他の tab の知らせ（51）が、48 の数えた「note」を言わず、フォームが空で始まることも言わない | S16b：`UpdateImpact.vue` と同じ数（Run・候補・note）を `update-other-tab-loss` に |
| L4 | 起動ページ 55 が「the first time it is opened」と言いながら、保存された release を出す | S16b：保存された release があるときは別の文（例「LOSAT Web cannot reach this site, and this browser has not chosen a release to start offline.」） |

参考：I1 電話の見出しの増加（受け入れ）、I2 Firefox の失敗の理由の末尾「.).」（括弧に入れる前に最後の句点を除く）、I3 電話の 43 で質問の枠がリンクの 3 px 下（`.example-confirm` に 6 px の上余白）、I4 撮っていない場面（設定と Job Title を言う「Try an example」の質問、完了した作業のある 49）、I5 48 の最初の箇条書きの印が黄色の外、I6 電話の 60 で同意が最初の画面の下、I7 ブラウザごとのダウンロード失敗の文の違い。Low 4 件と I2・I3 は S16b の最初の変更に入れる（見た目と文言の変更なので、直したら新しいゲートの実行と、変えた状態の画面レビューが要る）。記録はタスクフォルダの `reviews/screen-review-2.md`、比較の script と切り出しは `tmp/review-screens-2/`。

## 見つけて直したもの

- コードレビュー 1 回目の M1〜M3・L1〜L5（と M4）、2 回目の L1〜L5、3 回目の M1・L1〜L3・I1・I3、4 回目の L1（上）。
- 画面レビュー 1 回目の L1〜L6・I2〜I6（上）。
- WP-P の全体の実行で見つかった：`analytics.spec.ts` が WebKit で 6 件落ちる（Service Worker が転送した要求を `page.route` が見ない）。試験が Service Worker を塞ぐようにして直した（`56e8ee4a`、`1ca293df`、`fb383193`）。もう 2 件は見出しの位置の flaky で、M4 で直した。
- 状態 56 で見つかった：`register()` が何も返さない（塞がれた）とき、何も言わなかった。「service workers are blocked in this browser」（`5e860b74`）。
- 修正 3 回目で見つかった：`NODE_ENV` が `development` や `test`（Vitest の子）だと、開発 build が各 Vue の部品の絶対パスを bundle に入れる（7 か所）。ゲートが外し、`vite build` が拒む（`aeca5d93`、`cc5ca34f`）。
- 修正 2 回目で見つかった：状態 51 が「ready」を待っていた（修正 1 回目が「stored」を与えるため、6 件落ちた）。`9f0f3bac`。
- WP-U の probe：Playwright の `setOffline` は Firefox では Service Worker の要求に効かず、WebKit ではページの遷移そのものが壊れる。V-OFF は server を止める（判断 6、10）。
- WP-U・WP-R：Chromium は `/r/` の 404 を immutable として覚える。起動ページと Service Worker が読み直す。Pages がどうするかはプレビューで見る。
- WP-P：Firefox は `/favicon.ico` を自分で尋ね、`/sw.js` の更新を自分で取る（許した）。WebKit は worker から CSP の report を送らない（証拠は失敗した fetch と空の他の origin の log）。
- ゲート 1 で見つかった：V-OFF のケース 1 の保管場所の下限が Firefox で偽（Cache API の body の圧縮。アプリは正しく、試験の誤り）。`4c5bf901`（コードレビュー 3 回目 M1 と同じ指摘）。
- ゲート 2 で見つかった：プレビューの確認の script が、修正 1 回目でなくなった `offline-ready` の test ID を待っていた（dry run の 3 件、3 分の時間切れ）。`35970fa1`。ゲート 1 から dry run までの間、script がアプリの変更に追随していなかった。
- 修正 4 回目の後の E2E で見つかった：Chromium が 1 つのページの 11 個目の download を、最初から 1 s 以内なら黙って捨てる（`privacy.spec.ts:150` の flaky）。アプリは 1 回の click に 1 つしか download しないので、試験の側で 150 ms あけて直した（`9cb4e61e`、判断 18）。
- コードレビュー 3 回目で見つかった：`NODE_ENV` の規則の回避（`.env`）、`SecurityError` の取り違え、起動ページの行き止まり（上）。
- 画面レビュー 2 回目で見つかった：Low 4 件（release ID の全文、「Offline:」の読み違い、他の tab の知らせの数、起動ページ 55 の文）。S16b で直す（判断 19）。
- コードレビュー 4 回目で見つかった：文書が修正 4 回目の起動ページの動きに追随していなかった（`2b553a4b`）。

## 合流

`feature/losat-web-gui` への merge は、このセッションでは行わない。エンジン側が行う。W7 は未完なので、README の表の S16 の行は「完了」にしない。merge のときにエンジン側が行うこと：

1. `origin/feature/losat-web-gui-app` を `feature/losat-web-gui` に merge する（アプリ側の起点は `e45fbd09`。エンジン側が `141615ec` から進んでいれば通常の merge）。エンジン側のパス（`LOSAT/`、`web/adapter/`）、計画、README の表、セッションの README は、このブランチで変えていない。`web/app/package.json` の `engines.node` は `>=22.18` になった（CI は Node 24）。
2. README の表の S16 の行は「S16 完了、S16b（プレビューでの確認）待ち」とし、この記録へのリンクを置く。S16b の行を足す（指示書は `session_s16b_w7_preview.md`）。計画 §7 の S16 の行の条件は S16b が閉じる。
3. 計画に、GA4 の方式 (b) の判断（Owner-delegated、保守者の確認待ち）を書く（計画 §5.10 の「S16 で決める」と §10 の「GA4 の方式」の行）。計画 §10 の「Cloudflare に残す旧版の数、ドメイン名」は保守者の判断待ちのまま。
4. `.github/workflows/web.yml` の browser-engine の job が `tests/unit/self-test-engine.test.ts` を走らせる。V-OFF・V-PRIV・TD-11 は `web` と `browser-engine` の両方の `npm run e2e` に入っている。
5. `docs/web/abi_v2.md` の §4・§9 の kind 0 の文はエンジン側が書き直す（S14 の残件。S15 の記録のまま）。
6. 判断 1〜20 は推奨案で進めた（Owner-delegated）と記す（W4・W4b・W5・W6 の判断と同じ扱い）。

## エンジン側への申し送り

- **reactor とネイティブの CLI**：`141615ec` の木から作ったもの（W6 と同じ）。release ID の後ろの 8 桁はエンジンの digest（`e2fb7128`）で、reactor を作り直すと変わる。エンジンを替えたら、`web/app/src/examples/self-test.json` の期待値を作り直す（`tools/examples.py expected`、oracle の lock の下。`self-test-engine.test.ts` が合うまで落ちる）。
- **release ID の入力**：`web/app/build/inputs.ts` が、`web/app` に加え、ビルドが読むリポジトリの file（`docs/web/verification_cells.tsv`、認証の記録の表、`LOSAT/tests/wasi_shared_memory.js`）を列挙する。エンジン側がこれらを変えると、アプリ側の release ID の `-dirty` と bytes に関わる。
- **セッションファイルの `engineBuild`**：読み込んだ Run のバッジは、記録された build 名とこのサイトの `engineBuildName` を比べる（W6 のとおり）。release ID と `app.build`（commit の先頭 8 桁）の付け方は S16 で決めたので、変えるときはこの比較が壊れないようにする。
- **E2j（S13+）の後のアプリ**：E2j は S16 でも merge されておらず、E2j の指示書の 6.（アプリへの申し送り）は行っていない（W4b・W5・W6 の記録のとおり）。
- **BLASTX**：例題と自己試験の期待値は release に入っている（skip）。SX が merge されたら、BLASTX を有効にし（`domain/programs.ts` の `unavailable`）、自己試験の skip と「Try an example」の無効を外し、実検索の試験を足す。
- **エンジンのメモリ不足**（4000 写しの自己検索）は W6 のまま。協調取消（DW-22）の閾値の測り直しは S17。

## S16b・S17 への申し送り

[S16b の指示書](../../losat_web_gui_sessions/session_s16b_w7_preview.md)に「S16（W7）から引き継ぐこと」を書いた（最後の木と件数、ゲートの script、判断、残件）。S16b は、保守者のデプロイの後に、プレビューの確認（`crossOriginIsolated`、CSP、200・304・404 の headers、immutable、WebKit の取消の後の検索、Offline ready）を最初に行い、wrangler の bundle と `_routes.json` を見て、画面レビュー 2 回目の Low 4 件（判断 19）と残件から選んで直し、ゲートを再実行し、この記録を更新する。[S17 の指示書](../../losat_web_gui_sessions/session_s17_g_release_decision.md)に「S16（W7）から引き継ぐこと」（release ID と ReleaseManifest、オフラインと更新、利用統計の状態、V-PRIV・V-OFF・TD-11 の置き場所、自己試験と例題、配信の計測、残件）を足した。S17 は S16b が終わってから始める。

## 残件と注意

- **プレビューの確認が未実施**（上の判定）。ローカルの server は Pages の挙動を真似た（headers を 200・304・308・404 に付ける、強い ETag、`.html` の redirect）が、Cloudflare Pages が `_headers` を 304 と 404 に付けるかは Cloudflare の文書に無く、確かめていない。付かなければ WebKit は取消の後に作り直す Engine worker を拒む（W3 の判断 12）。
- **wrangler が `functions/api/visit.ts` を bundle できるか**（`../../src/domain/analytics` を import する）は試していない。生成される `_routes.json` が全要求を Function に通すと、静的 file が無料枠の 1 日 100,000 要求に数えられ得る。S16b が `wrangler pages functions build` とデプロイの出力で見る。
- **コードレビュー 1 回目の情報**：I1（shell の保存が atomic でない）、I4（自己試験を取り消せない）、I6（V-OFF 3 の限界）。**2 回目の情報**：I1（V-OFF 3b が決定的でない）、I2（準備が作業と重なる場面の数）、I3（close のときの beacon）、I4（知らせが読み込みごとに 1 度ページを動かす）。**3 回目の情報**：I2（V-PRIV は request header の符号化を試験文字列の 3 つの形でしか探さない）。**4 回目の情報**：I1（click の時点の再確認は、別の tab の掃除との競合を狭めるが閉じない。見える message になり、黙って切り替えることは無い）、I2（cookie やサイトデータを塞いだ WebKit では probe していない。「failed」と「Try again」になるはずで、安全な側）、I3（保管場所の確認の余裕は、保存した release が使用量の大部分であることに頼る）。
- **Chromium の JSON の書き出しの最長のフレームの間隔**：ゲート 3 は 250 ms（W6 のゲート 3 は 117 ms）。W6 と W7 の対の A/B（3 回ずつ）で 183 / 133 / 250 対 200 / 167 / 133 ms、中央値 183 対 167 ms で、W7 の退行ではない。100 ms を超える間隔は W6 からあり、残件（100 ms の目標を満たさない。原因は追っていない）。WebKit の書き出しの開始の 117〜134 ms は既知。
- **画面レビュー 1 回目の記録のまま置いたもの**：電話の見出しは +49 px（検索の状態で +82 px）、受け入れた（I1）。
- **WebKit の scroll anchoring**：WebKit は scroll anchoring を持たず、Alignments の節が遅れて読まれると下の内容が押し下げられる（W6 の判断 17）。Alignments の節の高さを読む間確保するのは未実施。V-OFF・V-PRIV・画面の記録の spec は、見出しが「Offline ready」に落ち着くのを待って避けている。
- **WebKit の `navigator.onLine`**：Playwright の WebKit は `true` のまま（判断 10）。起動ページの message の種類は `failed`・`no-index`。Playwright の offline の emulation は Chromium でしか確かでないので、V-OFF は server を止める。
- **Chromium は `/r/` の 404 を immutable として覚える**（WP-U、判断 14）。起動ページと Service Worker は読み直す。Pages が 404 にも immutable を付けるかは未確認。
- **V-PERF の lock の規則を 1 回破った**：WP-U の V-OFF の 1 回の実行（Firefox と WebKit）が `vperf.lock` の存在する間に始まった（上）。WP-U の記録に残っている。
- **W6 から引き継いだもの（S16 が触れていない）**：WebKit は書き出しごとの最初に 128〜155 ms フレームを描かない。隠れた tab では timer が約 1 s に絞られ、書き出しとセッションの保存が遅くなる。貼り付けた入力は session に `query.fa` / `subject.fa` として記録される。W5 コードレビュー L3（行の対応が先読みの状態を持たない）。W5 画面レビュー 2 回目の I3・I4（グラフの行が 8 px、電話のタブの行）。BLASTX の検索は SX まで動かせない。`abi_v2.md` の kind 0 の文。エンジンのメモリ不足（4000 写し）、WebKit の保存の上限（512 MB）、`runInputs` の解放、iOS Safari（V-MOB、S17）。
- **E2j の申し送りは未実施**（上）。E2j は merge されていない。
- **画面の記録（PNG）はリポジトリに入れない**：ゲート 3 の画面の記録は `$BUILD_ROOT/s16-gate3-screens/`（`/home/kawato/.cache/losat-work/s16-gate3-screens/`。396 枚、SHA-256 は実行記録の `screens.sha256`）。途中の記録は `s16-wpgp-screens`（`70442212`、366 枚、画面レビュー 1 回目）、`s16-fix2-screens`（`9f0f3bac`、186 枚）。状態の一覧（66 枚 / ブラウザ / サイズ）：01〜28 は W5 と同じ並び（04 は program ごと）、29〜40 と 30b が W6、41〜45 が導入支援（Getting started、自己試験、例題）、46〜56（と 46b・46c・55b）が配信（オフラインの状態、更新の知らせと影響、同意、起動ページ）、57〜60 が修正 2 回目の追加（更新の取得中、チェックサムの失敗、待ち、測るサイトの初回）。
- **ゲートの script の順序**：`run_gate.sh` の `all` は、結果画面の E2E の繰り返しで落ちると、計測と画面の記録に進まない（W6 のゲート 1 と同じ）。W7 では 1 回の失敗ごとに最初からやり直したので 3 回かかった（ゲート 1 はエンジン入りの E2E、ゲート 2 は dry run で止まった）。`after-review`・`measure` で一部だけ走らせられる。
- **背景の作業の記録**（担当の報告、レビュー、判断の一覧）はタスクフォルダ `/home/kawato/losat-baselines/s16-w7-20261011/`（リポジトリには入れていない）。
