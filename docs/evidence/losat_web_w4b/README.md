# LOSAT Web W4b（Session S13b）ゲート記録

- 段階：W4b 画面を NCBI BLAST Web に寄せる（[総合計画書](../../losat_web_gui_plan.md) §7 の S13 と S14 の間に足す S13b、[指示書](../../losat_web_gui_sessions/session_s13b_w4b_ncbi_style_ui.md)、[対応表](../../web/ncbi_ui_mapping.md)）
- ブランチ：`feature/losat-web-gui-app`（アプリ側。worktree `$WORK_ROOT/.worktrees/web-gui-app`、Linux の clone）。起点は `ff9ebc49`（S13 の最後）。エンジン側は W4 を `feature/losat-web-gui` に merge していない（`origin/feature/losat-web-gui` は `4596b82f`）ので、起点に merge するものは無く、reactor とネイティブの CLI は S13 のもの（`app-s12-reactors`、ネイティブ sha256 `f0b8916b…`）を使った。エンジン（`LOSAT/`、`web/adapter/`）は変えていない
- 実行記録：[`run-20261009T184556Z/`](run-20261009T184556Z/)（ゲート、commit `f097595d` の木、`21cc0ff5`）、[`run-20261009T210300Z-after-review/`](run-20261009T210300Z-after-review/)（`087eaa13`、`81e012b8`）、計測だけの実行 [`run-20261009T211154Z-measure/`](run-20261009T211154Z-measure/)（`087eaa13`、`81e012b8`）、[`run-20261009T215513Z-after-review/`](run-20261009T215513Z-after-review/)（`492fad48`、最後のアプリの木、`10443d06`）。どれも作成後は書き換えない。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)、再現は [`run_gate.sh`](run_gate.sh)
- 判定：**完了条件を満たした**。ゲートの実行（`f097595d`）のすべての段階、コードレビュー 2 回（1 回目：妨げになる指摘無し、M1・L2〜L8 を直した。2 回目（修正の後）：妨げになる指摘無し、低 2 件は S14 へ）、画面レビュー 3 回（どれも High 無しの合格。1 回目と 2 回目の中程度の指摘は直し、3 回目の低 1 と参考 1 は S14 へ）、レビューの後の実行 2 回と計測だけの実行が通った。下の「完了条件と結果」。計測は `087eaa13` の木で行い、最後の木（`492fad48`）では測り直していない（下の「実測」）
- 保守者の判断待ち：判断 1〜88（下の「判断」）は推奨案で進めた（Owner-delegated、2026-09-29 の常設の指示、2026-10-07 に再掲）。判断 7（ボタンの名前「Run LOSAT」）は保守者の指示（2026-10-10）であり、委任ではない。目に付くのは 5（Job Title を RunSnapshot に持つ）、6（BLASTN の既定は NCBI の 2 配列の画面の blastn でなく、CLI と同じ megablast）、18（Edit Search は S15）、20（Percent Identity のフィルターを採らない）、24（ドットプロットの SVG の書き出しは S15）
- `feature/losat-web-gui` への merge：このセッションでは行わない。エンジン側が W4 と W4b をまとめて merge する（下の「合流」）
- このセッションは `/home/kawato/losat-baselines` で起動した（clone の外なので skill と agent が読み込まれなかった）。skill はファイルから読み、`losat-reviewer` の定義は general-purpose の agent に渡した（S13 と同じ）

## 保守者の指示

- **2026-10-09（S13 の最中）**：「次のセッションで NCBI BLAST ウェブサイトの GUI にデザインを寄せてほしい」。範囲は保守者が確かめた：NCBI BLAST の配置・節の順・言葉に寄せ、ロゴ・ヘッダー・配色などのブランドは写さない。値はエンジンが書いたまま示す規則は変えない。同日、ドットプロットについて「もしできたら [`blast2dotplot.py`](https://github.com/satoshikawato/bio_small_scripts/blob/main/blast2dotplot.py) を参考にしてほしい」「もっと洗練されて軽快に動いてクリックしてポップアップするようなのができるならそれはそれでいい」。W4 のゲート記録の「保守者の指示」に同じ記録がある。
- **2026-10-10（このセッションの最中）**：「Run BLASTじゃなくてRun LOSATにしてね」。検索のボタンは「Run LOSAT」（キューに積むときは「Run LOSAT (add to queue)」、Run が複数なら「Run LOSAT (2 runs)」）。「BLAST」を操作の名前にしない。これは保守者の決定で、判断としては取らない（対応表に反映、`d8de88b7`）。

## コミット

`ff9ebc49..10443d06` の 36 コミットと、この記録のコミット 2 つ（古い順）。

| コミット | 内容 |
|---|---|
| `9b0f15fd` | 文書：NCBI BLAST の画面と LOSAT Web の対応表 [`docs/web/ncbi_ui_mapping.md`](../../web/ncbi_ui_mapping.md)（作業 2） |
| `d8de88b7` | 文書：検索のボタンは「Run LOSAT」（保守者の指示、2026-10-10） |
| `fa5239de` | アプリ：検索画面を NCBI BLAST の順と言葉に（作業 3） |
| `e8e31a4f` | 文書：描画した NCBI の結果から対応表への補足（Clusters タブ、当たり無しの帯、ドットプロットの画像） |
| `778989b2` | アプリ（試験だけ）：Job Title のキュー上の位置は Run が終わった後に一度だけ読む |
| `2cd0f56e` | アプリ：図の縮尺と幾何を domain の関数に（`plot-scale.ts`、`plot-geometry.ts`） |
| `95687199` | アプリ：`blast2dotplot.py` に倣ったドットプロット（hover と HSP のポップアップ） |
| `743a01de` | アプリ：NCBI に倣った Graphic Summary（コンポーネント） |
| `96dc54a9` | アプリ：W4 のドットプロットの CSS を外す（図の CSS は `ui/plots.css`） |
| `22398be6` | 文書：対応表の検索の一文は「optimized for」 |
| `9ed75293` | 文書：ドットプロットの目盛りの単位は軸の題に置く（判断 B2-1） |
| `2fc21fe1` | 文書：W4b の `run_gate.sh`（W4 のコピー。記録の場所だけ W4b） |
| `7f58c537` | アプリ：結果画面を NCBI BLAST の順と言葉に。Graphic Summary と Alignments をつなぐ（作業 4） |
| `c5019bc8` | アプリ：5 px より狭い小目盛りを省く（狭い TBLASTN の図が灰色に塗れていた） |
| `f097595d` | アプリ：ドットプロットは HSP ごとに線を描き、同じ画素の不透明な線を省く（5,993 HSP のズームが Firefox で 30〜43 ms）。**ゲートの木** |
| `21cc0ff5` | 文書：`f097595d` のゲートの実行記録（すべての段階が通った。落ちたのは見込みどおりの場面の記録） |
| `9d95c5aa` | アプリ：ドットプロットと Graphic Summary の計測の時計を、描いたフレームで止める（コードレビュー M1） |
| `81f25039` | アプリ：Alignments は同じフレームに届いた節をまとめて適用し、失敗した節を読み直し、subject の見出しの読みを共有する（コードレビュー L2〜L5） |
| `57bc7f91` | アプリ：ドットプロットは 3 nt per aa で縦横比を決め、見えているもので描いて選び、スクリーンリーダー向けに application にする（判断 26、コードレビュー L6〜L8） |
| `d949ee7a` | 文書：対応表の縮尺は 3 nt per aa |
| `da0c24ca` | アプリ（試験だけ）：計測は一覧の選択と並べ替えの前に画面が落ち着くのを待つ |
| `d4c367a2` | アプリ：選んだ HSP の outfmt 6 の行は 1280 px で全欄を Range の塊に示す（画面レビュー 1 回目 M1） |
| `e81f791f` | アプリ：ドットプロットの小格子を両軸に同じに、軸の題に単位を残して 3 nt per aa を示し、ポップアップを線の脇に（M2・M3・4・L10・L11） |
| `0fe326b8` | アプリ：TBLASTN の HSP は subject の frame だけを示す（「Subject frame」、+2）（L9） |
| `9f99f956` | アプリ：キューのカードはどの状態も同じ形（L7） |
| `596c5aa6` | アプリ：Descriptions は残りの幅を Description に、数は言葉で、間隔を直す（L8・L12・L13） |
| `0659584c` | アプリ（試験だけ）：画面の記録は状態 07・08 の HSP の行を見える所で押す（L5） |
| `087eaa13` | 文書：対応表を画面レビューに合わせる。**1 回目の修正の後の木** |
| `81e012b8` | 文書：`087eaa13` の after-review と計測の実行記録（すべて通った。開く時間は W4 の水準に戻った） |
| `cc3aa24f` | アプリ：「Open results」は結果の見出しの上に 8 px 空ける（F3-5） |
| `06c1a290` | アプリ：電話の「Results for」は query ID を離して示し、長さと数を ID の下に（画面レビュー 2 回目 M1） |
| `42742404` | アプリ：label を切る select は矢印の場所を空ける（Chromium の Run の select、L1） |
| `63098d40` | アプリ：キューのカードの Cancel と「Open results」は段階の行の右端に（L2） |
| `a7b4e0bf` | アプリ：入力のカードと「複数の入力」の選択はブロックの欄の端に揃える（L3） |
| `492fad48` | 文書：対応表を 2 回目の画面レビューに合わせる（Results for、キューの操作）。**最後のアプリの木** |
| `10443d06` | 文書：`492fad48` の最後の after-review の実行記録（すべて通った） |
| このゲート記録を含むコミット | 文書：W4b のゲート記録、NCBI の参照画面の SHA-256 と索引、`evidence.sha256` |
| S14 への引き継ぎのコミット | 文書：S14 の指示書に「S13b（W4b）から引き継ぐこと」 |

エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表は、アプリ側のブランチでは変えていない（README 規則 1）。`docs/web/` は、このセッションが作った `ncbi_ui_mapping.md` と、`results_columns.md` の 1 行（`7f58c537`、当たりの無い query の言葉）だけを変えた。依存ライブラリは増やしていない。

## 完了条件と結果

| 完了条件（指示書の末尾） | 結果 | 証拠 |
|---|---|---|
| 対応表が参照画面の要素を漏れなく挙げる | 通過。[`ncbi_ui_mapping.md`](../../web/ncbi_ui_mapping.md) は撮った画面の要素（検索画面の節・欄・ボタン、結果の見出しの塊・Results for・Filter Results・各タブ・Dot Plot の画像）を 1 つずつ挙げ、採否（採る・寄せる・採らない・エンジン待ち・後の段階・LOSAT だけ）と理由、LOSAT だけのものの置き場所を書く。撮らなかった画面（下の「参照画面」。Taxonomy タブ、Download と Select columns のメニュー、Edit Search の画面）の要素は、指示書と下調べの文書に従って「採らない」または「後の段階」とした | `docs/web/ncbi_ui_mapping.md`（`9b0f15fd` から `492fad48` まで 9 つのコミットで更新）、`ncbi_reference_index.md` |
| 検索画面と結果画面が対応表に従う（画面レビューの合格） | 通過。画面レビューは 3 回行い、どれも High 無しの合格（下の「画面レビュー」）。1 回目（ゲートの記録 150 枚）は M1〜M3・L5〜L13 を、2 回目（`087eaa13` の 150 枚）は M1・L1〜L3 を出し、すべて直した。3 回目（`492fad48` の 150 枚）は直ったことを確かめ、低 1 と参考 1 を S14 に残した。「Run LOSAT」の言葉、NCBI のブランドが無いこと、390 px で横に流れないこと、触れる対象が 24 px 以上であることも確かめた | `reviews/screen-review.md`、`screen-review-2.md`、`screen-review-3.md` の要約を下の表に載せた（記録はタスクフォルダ）。記録の PNG は `$BUILD_ROOT/s13b-gate-screens/`、`s13b-review-screens/`、`s13b-screens-f2/`、`s13b-screens-f3/`、`s13b-final-screens/` |
| 既存と新しい E2E が 3 ブラウザ・両方のビルドで通る | 通過。最後の木（`492fad48`）で FakeEngine のビルド 92 件（Chromium 35・Firefox 29・WebKit 28）、エンジン入りのビルド 131 件（48・42・41）。ゲートの木（`f097595d`）も同じ件数で、エンジン入りでは V-BR（ブラウザごとに 27 の検索の outfmt 0・6・7 と診断がネイティブの CLI と一致し、NCBI の凍結バイトとの比較 21（outfmt 0 が 18、outfmt 7 が 3）がすべて一致）を含む。結果画面の E2E の 2 回の繰り返しも通った（66 件）。新しい E2E は Job Title・Run LOSAT・Algorithm parameters（開閉、印、Restore）・Graphic Summary・Alignments の Range・ドットプロットのポップアップと縮尺・各幅の収まり（1280 px と 390 px）・キューのカードの形・「Results for」の ID | `run-20261009T184556Z/`（`npm-e2e-fake-engine.log`、`npm-e2e-engine.log`、`e2e-results-repeat.log`）、`run-20261009T215513Z-after-review/` |
| HSP の対応の試験が通る（変えずに） | 通過。`hsp-correspondence.test.ts` は `ff9ebc49` から変わっていない（`git diff ff9ebc49 HEAD` が空）。BLASTN・BLASTP・TBLASTN・TBLASTX の全 fixture（165 の検索、30,467 HSP、outfmt 0 の節 30,027、outfmt 0 に現れない HSP 440）を serial の reactor で実行し、出力は NCBI の凍結バイトと 196 が byte 一致、承認済みの `-db_gencode` の例外 2 がデータベースの検索として一致、field list の 3 は比べない。単体試験は reactor 付きで 407 件 | `run-20261009T184556Z/unit-cases.log`（393 件、ゲートの木）、`run-20261009T215513Z-after-review/unit-cases.log`（407 件、最後の木） |
| 計測が S13 から大きく悪くならない | 通過（`087eaa13` の木）。開く時間は W4 の水準に戻り、操作の時計は W4 と同じ 1〜2 フレームの水準。W4 より遅く見えるのは、Firefox の 3000 写しの「一覧で選ぶ」（26 ms、W4 は 14 ms）、Chromium の 1 万 query の「終わりまで scroll」（15 ms、W4 は 1 ms）、WebKit の 200 subject の query のクリック（10 万 query で 43 ms、W4 は 15 ms。W4 自身が 1 万 query で 47 ms と揺れ、ゲートの木は 18 ms）など、1〜2 フレームの差。ゲートの木（`f097595d`）では開く時間が 16〜59% 増えたが、交互に測り直すと W4 の木と同じ水準で（Firefox の 1500 写しだけ +15%、標本が重なる）、原因は実行ごとの揺れだった。下の「実測」 | `run-20261009T184556Z/measure/`、`run-20261009T211154Z-measure/measure/` |

### ゲートの実行（`f097595d`）

`run_gate.sh` のすべての段階が通った：V-ABI quick（15 の検索 × 経路とスレッド = 60 の実行がネイティブの CLI と一致、NCBI の凍結ハッシュ 16/16 一致）、`npm run check`（393 件のうち 364 件、reactor の要る 29 件は skip）、reactor 付きの単体試験（393 件。HSP の対応の試験と、生成した検証の表 `verification-table.json`）、FakeEngine のビルドの E2E（Chromium 35・Firefox 29・WebKit 28、計 92 件）、エンジン入りのビルドの E2E（Chromium 48・Firefox 42・WebKit 41、計 131 件）、結果画面の E2E の 2 回の繰り返し（計 66 件）、計測（3 ブラウザ）、画面の記録（3 ブラウザ × 2 サイズ、150 枚、状態 01〜23）。計測の見込みどおりの場面の失敗は、W4 と同じ：20 nt の単位を 4000 写しの自己検索は 3 ブラウザでエンジンのメモリ不足（`memory allocation of 49881 bytes failed`）、3000 写しは WebKit で保存の上限（512 MB）。

### レビューの後の実行

- **`087eaa13`**（コードレビューの指摘と画面レビュー 1 回目の指摘を直した木）：`LOSAT_WEB_GATE_STEPS=after-review` が通った（[`run-20261009T210300Z-after-review/`](run-20261009T210300Z-after-review/)）。`npm run check` 378 件（29 件 skip）、reactor 付きの単体試験 407 件、FakeEngine のビルドの E2E 92 件、エンジン入りのビルドの E2E 131 件、画面の記録 150 枚。続けて `LOSAT_WEB_GATE_STEPS=measure`（[`run-20261009T211154Z-measure/`](run-20261009T211154Z-measure/)）。コードレビューの M1 で時計の意味が変わった（下の「実測」）ので、測り直した。
- **`492fad48`**（画面レビュー 2 回目の指摘を直した最後の木）：`after-review` がもう一度通った（[`run-20261009T215513Z-after-review/`](run-20261009T215513Z-after-review/)）。件数は `087eaa13` と同じ（check 378 件・29 件 skip、単体 407 件、FakeEngine 92 件、エンジン入り 131 件、画面の記録 150 枚 `$BUILD_ROOT/s13b-final-screens/`）。計測は測り直していない（F3 の変更は検索画面の入力カード、「Results for」の行、キューのカード、select の余白で、開く経路とプロットの描画には触れない）。
- どの実行の E2E の log にも `FAIL [opfs] storage full …` と `FAIL [run output] storage full …`（FakeEngine 3 行、エンジン入り 6 行）が出る。ブラウザの中の契約の確認が自分の記録として印字するもので、それを印字する Playwright の試験は通っている（落ちた試験は無い）。ゲートの実行にも同じ行がある。

## 参照画面

- **取得**：2026-10-10 01:28〜02:00（JST）、UTC では 2026-10-09 16:28:31〜17:00:26。Playwright の headless Chromium で `https://blast.ncbi.nlm.nih.gov/` の画面を撮った（1280×900。blastn の検索と結果だけ 390×844 も）。検索は 6 件で、NCBI の上限（6）に達した。
- **送った配列**：公開の RefSeq の accession だけ（NM_000518.5、NM_000519.4、NM_000558.5、NP_000509.1、NP_000510.1）。利用者のデータは送っていない。

| 検索 | RID | query | subject |
|---|---|---|---|
| blastn | `CJ5K9KV3114` | NM_000518.5 | NM_000519.4 |
| blastp | `CJ5W8ZCC114` | NP_000509.1 | NP_000510.1 |
| tblastn | `CJ6FDRHH114` | NP_000509.1 | NM_000519.4 |
| tblastx | `CJ6G5X6T114` | NM_000518.5 | NM_000519.4 |
| multi（複数 query） | `CJ6GXH5H114` | NM_000518.5 と NM_000519.4 | NM_000558.5（query 1 は当たり無し、query 2 に弱い当たり 1 つ） |
| multihit（追加） | `CJ7C8J2Y114` | NM_000518.5 と NM_000519.4 | NM_000519.4（両方に当たり） |

- **ファイル**：174 個（PNG 65・txt 61・html 48、合計約 19 MB。PNG のうち 5 つは NCBI が返したドットプロットの画像）。`$BUILD_ROOT/s13b-ncbi-reference/`（`/home/kawato/.cache/losat-work/s13b-ncbi-reference/`）に置き、リポジトリには入れていない。検索画面（blastn・blastp・tblastn・tblastx、Algorithm parameters を開いた状態、blastn の discontiguous megablast と既定と違う値）、2 配列の結果（results-top、Search Summary、Descriptions、Graphic Summary と hover と click、Alignments、Dot Plot とその画像）、複数 query の結果（「Results for」）。txt は `document.body.innerText`、html は描画後の DOM の `page.content()`（実行していない）。
- **SHA-256 と索引**：174 個のハッシュは [`ncbi_reference.sha256`](ncbi_reference.sha256)（`sha256sum` の形式。`INDEX.md` と自分自身を除く。置き場所で `sha256sum -c` が 174 個とも一致することを確かめた）。URL、時刻（UTC）、RID、program、画面の大きさ、写したもの、ハッシュを 1 ファイルずつ書いた索引は [`ncbi_reference_index.md`](ncbi_reference_index.md)（`INDEX.md` のコピー）。どちらも公開の accession、URL、RID、時刻、ハッシュ、撮影の script の置き場所だけで、利用者のデータは含まない。
- **撮れなかったもの**：Taxonomy タブ、Download と Select columns のメニュー、Edit Search の画面。headless Chromium が描かないもの（Graphic Summary の hover の title の吹き出し、「Results for」の select の選択肢の一覧）は、click の popover（`*-graphic-click`）と選択肢の txt（`multi*-queryList-options-1280.txt`）で代えた。ドットプロットで鎖の色を見分けられる例は無かった（NCBI の線はすべて濃い灰色）。390 px は blastn だけ（NCBI の検索画面は幅 2330 px、結果は約 817 px のまま折り返さず、PNG が 390 px より広い）。
- **索引の食い違い**：`ncbi_reference_index.md` の「Searches submitted」の表は blastn の RID を `None`、query と subject を空欄としている（撮影 script が最初の検索の RID を拾えなかった）。blastn の RID `CJ5K9KV3114` は各ファイルの URL と撮影担当の報告にあり、query と subject は対応表と報告の記録のとおり。索引は撮ったままのコピーなので直していない。

## 作業との対応

| 作業 | 実装 | 試験 |
|---|---|---|
| 1. 参照画面を固定する | 上の「参照画面」。撮影の script と log はタスクフォルダ（リポジトリには入れていない） | ハッシュの照合（上） |
| 2. 対応表 | [`docs/web/ncbi_ui_mapping.md`](../../web/ncbi_ui_mapping.md)：検索画面、結果の見出しの塊、タブ、Descriptions、Graphic Summary、Alignments、Dot Plot（blast2dotplot.py との対応）、LOSAT だけのものの置き場所、見た目。採否は判断 1〜34 | 画面レビュー（対応表と並べて見る） |
| 3. 検索画面 | NCBI の順（program のタブと一文、Enter Query Sequence、Enter Subject Sequence、Program Selection、ボタンと一文、Algorithm parameters、2 つ目のボタン）。Job Title（`RunSnapshot.title`、argv に入れない）、Optimize for / Algorithm の radio、「Run LOSAT」、既定と違う値の黄色と ♦、Restore default search parameters、General / Scoring / Filters and Masking / Discontiguous Word / Other Parameters の節、Query・Subject の遺伝暗号の位置。`src/ui/` の `SearchPanel`・`InputPanel`・`SourceCard`・`RegionPicker`・`RecordList`・`ParameterForm`・`ParameterField` | `search.spec.ts`、単体 `search-form.test.ts`・`draft.test.ts`・`coordinator.test.ts`（title）、画面の記録 01〜06 |
| 4. 結果画面 | 見出しの塊（Job Title、Run、Program と検証バッジ、Options、Query / Subject ID と長さ、Download All）と右の Filter Results、「Results for」、タブ Descriptions / Graphic Summary / Alignments / Dot Plot / Run details / Outputs。`SubjectTable`（Descriptions）、`GraphicSummary`、`AlignmentsView`（Range、Previous / Next subject、選んだ HSP の outfmt 6 の行と outfmt 0 の節）、`DotPlot`（`blast2dotplot.py` の約束、ポップアップ、3 nt per aa）、`ResultFilters`。図の縮尺と幾何は `src/domain/plot-scale.ts`・`plot-geometry.ts` | `results.spec.ts`（4 program、FakeEngine とエンジン入り、3 ブラウザ）、単体 `results.test.ts`・`plot-scale.test.ts`・`plot-geometry.test.ts`、計測 |
| 5. 見た目 | 節の枠と見出し、Algorithm parameters の帯、タブの帯と道具の行、表の字（0.875em）、キューのカードの形。色は LOSAT の配色（Graphic Summary の凡例とドットプロットの線の色だけ約束）。W4 の基準（1280 px で表が収まる、390 px で横に流れない、24 px 以上の対象）を保つ | E2E が 4 program × 2 幅で測る。画面レビュー |
| 6. 試験 | 既存の spec の期待値と補助を新しい画面に合わせ、test ID を保った（動いた ID は S14 の指示書の節）。HSP の対応の試験は変えない | 上の「完了条件と結果」 |
| 7. ゲートとレビュー | [`run_gate.sh`](run_gate.sh)（W4 のコピー）、コードレビューと画面レビュー 3 回 | 下の「実測」「独立レビュー」「画面レビュー」 |

## 作ったもの

- `src/ui/`：`SearchPanel`・`InputPanel`・`SourceCard`・`RegionPicker`・`RecordList`・`ParameterForm`・`ParameterField`（検索画面）、`ResultsPanel`（見出しの塊とタブ）、`QueryPicker`（Results for）、`SubjectTable`、`GraphicSummary`、`AlignmentsView`、`DotPlot`、`HspTable`、`ResultFilters`、`QueuePanel`（カードの形）、`useNarrow.ts`・`useSideScroll.ts`、`plots.css`（図の CSS）。`HspDetail.vue` は無くなった（詳細は選んだ Range の塊に）。
- `src/domain/`：`plot-scale.ts`（目盛りと縮尺。`tick_size` の表、5 px の規則）、`plot-geometry.ts`（線の描画範囲、選ぶ線、3 nt per aa の縦横比）。`run.ts` の `RunSnapshot.title`（任意）。
- `src/application/results.ts`：Range の順、節の読み（1 つの Run に 1 つの promise、失敗した節は読み直せる）、前後の subject。
- 文書：`docs/web/ncbi_ui_mapping.md`、[`run_gate.sh`](run_gate.sh)、この記録、[`ncbi_reference.sha256`](ncbi_reference.sha256)、[`ncbi_reference_index.md`](ncbi_reference_index.md)、S14 の指示書の引き継ぎ。

## 判断（推奨案で進めたもの）

Owner-delegated（2026-09-29 の常設の指示、2026-10-07 に再掲）。すべて 2026-10-10 に取った（判断 7 を除く）。由来ごとに分け、番号は通しで付けた。出典は取りまとめ役の `DECISIONS.md` と各担当の報告。

### 取りまとめ役（1〜25 は対応表の設計、26 は B1 の疑問を受けた判断、27〜32 は画面レビュー 1 回目、33・34 は 2 回目、35 は 3 回目の後）

1. 参照画面は headless Chromium、1280 px（blastn の検索と結果だけ 390 px）。検索は公開の RefSeq だけ。
2. 検索画面の枠を NCBI の順（query、subject、Program Selection、ボタン、Algorithm parameters）に積む。W3 の query と subject の左右の並びはやめる。
3. 「Enter FASTA sequence(s)」にする（NCBI の accession・gi の言葉は採らない。ネットワーク参照をしない、`web/AGENTS.md` 規則 4）。
4. 「Align two or more sequences」の印は置かない（LOSAT は常に subject を与える）。
5. Job Title を採る：`RunSnapshot.title`（任意）、argv に入れない、キューと結果の見出しに示す。自動では埋めない。
6. Program Selection は radio：BLASTN の「Optimize for」に 4 つ目（blastn-short）を足す。既定は CLI と同じ megablast（NCBI の 2 配列の画面は blastn）。BLASTP・TBLASTN の「Algorithm」は `describe` から。TBLASTX に節は無い。
7. **保守者の指示（2026-10-10）**：ボタンは「Run LOSAT」／「Run LOSAT (add to queue)」（Run の数つき）。以前は「BLAST: run locally」だった。ボタンの横に「… Runs in this browser.」で終わる一文（選んだ program と task）を置く。Algorithm parameters を開いたとき 2 つ目のボタンが下に出る。「Show results in a new window」は採らない。委任ではない。
8. Algorithm parameters は閉じて始め、argv に値が書かれるときは開く。見出しに数。
9. 既定と違う印は argv に書く値（`formParameters`）。黄色（LOSAT の注意の色）と ♦。
10. Restore default search parameters は Algorithm parameters の節だけを消す（task と遺伝暗号は残す）。
11. NCBI の select の選択肢（Word size、Match/Mismatch、Gap Costs、Max target sequences）は写さない。NCBI の名前と並びで自由な欄にし、エンジンが検査する。
12. Compositional adjustments の選択肢に NCBI の名前を添える。argv は 0〜3 のまま。
13. Query の遺伝暗号は Enter Query Sequence の中（NCBI のとおり）。Subject の遺伝暗号は Enter Subject Sequence の中（LOSAT だけ。承認済みの例外の注を保つ）。
14. TBLASTX の `-culling_limit` は「Max matches in a query range」（NCBI の `HSP_RANGE_MAX`、その説明が culling を指す）。
15. NCBI に欄の無い LOSAT の option は最後の「Other Parameters」に。
16. Short queries の自動調整、Reset page、Bookmark、Species-specific repeats、PSI/PHI/DELTA は採らない。
17. Discontiguous Word Options は task が dc-megablast のとき、または template の値があるとき出す。
18. 結果：Edit Search は S15。Search Summary は Run details タブへ、Download All は Outputs タブへ。Query / Subject Descr は見出しに置かない（outfmt 0 にある）。
19. 「Results for」は LOSAT の仮想化した query 一覧で、複数 query の Run だけ（10 万 query は select に入らない）。
20. Filter Results は E value ≤、Bit score ≥、Subject ID。Percent Identity は採らない（原値が無い）。Query Coverage は E2j 待ち。
21. Descriptions の列は #、Description、Score (bits)（Max Score の位置）、E value、HSPs、Length、Subject ID（最後。NCBI の Accession の位置）。既定の順はエンジンの順。
22. Graphic Summary は subject ごとに 1 行、NCBI の得点の階級と色を `bit_score` に当てる。最初の 100 subject と「Show all」。
23. Alignments は選んだ subject の塊と Previous / Next subject。Range n は subject の座標を小さい方から（NCBI の参照：Range 2: 588 to 608 は Sbjct 588–608）。HSP の表は塊の見出しの下。選んだ Range の前後 25 個ずつだけ示し、節は見える所だけ読む。
24. ドットプロットは `blast2dotplot.py` の約束（原点は左上、同じ縮尺、`tick_size` の表、`#D3D3D3` の格子、`#1f77b4` / `#ff7f0e`、`int(pident)` で不透明度）。短い辺は最小 120 px（「Axes not to scale」）。aa の軸は aa / kaa / Maa。SVG の書き出しは S15。選ぶとポップアップ。
25. 2 つの主タブで列の幅を同じにする（キューが動かない。W4 の画面レビュー低 4）。
26. TBLASTN・BLASTX のドットプロット：縦横比を決めるときだけ 1 aa を 3 nt と数える（目盛りとラベルは各軸の単位のまま）。コードレビュー L8 が実装の抜けを見つけた。
27. 電話の表が列の途中で終わるのは受け入れる（W4 の判断 17：横に流れる旨の一文つき）。
28. Descriptions の HSPs は E value の後に置く（E2j が Score と E value の間に Total Score・Query Cover を入れても動かない）。対応表を直した。
29. 片方の配列だけが翻訳されるときは、その frame だけを示す（TBLASTN は「Subject frame」）。
30. ドットプロットの軸の単位はスクリプトのとおり bp / kbp / Mbp（aa / kaa / Maa）。対応表の題の行の「(nt)」を直した。説明の行とポップアップは Run の単位（nt / aa）のまま。
31. キューのカードは 1 つの形（1 行目 Run の番号・program・バッジ、2 行目 段階と時間、3 行目 操作を右）。
32. ポップアップは HSP の脇（電話では図の下）。protein の軸の題に「drawn at 3 nt per aa」と示す。
33. 電話の「Results for」：ID は 12ch 以上を保ち、長さと数は 2 行目に。
34. キューのカードの操作は 2 行目の右端（終わったカードに余分な行を作らない）。
35. 3 回目の画面レビューは合格。低 1（電話の入力のカードで長いファイル名が 2 行増える）と参考 1（電話の radio が 14 px）は S14 に残し、もう 1 回ゲートを回すほどではないとした。

### A：検索画面（`logs/report-a.md`）

36. A1：argv の順は task、遺伝暗号（query、subject）、画面の順の節（TBLASTX の `-query_gencode` `-db_gencode` は `-evalue` の前）。期待値を直した。
37. A2：argv に書く欄はすべて黄色と ♦（task の radio と遺伝暗号も）。「(N changed)」の数は Algorithm parameters の欄だけ。
38. A3：Algorithm parameters は閉じて始まり、値のある program に切り替えると開き、自分では閉じない。
39. A4：レコード一覧と Combined / Separate は各ブロックの Genetic code と Job Title の行の後。
40. A5：ブロック全体がファイルの drop 先（`*-dropzone` は「Or, upload file」の箱のまま）。
41. A6：RegionPicker の凡例は「Query / Subject subrange」と「Record <id>」の行。レコードが 0 のときの一文。
42. A7：`search-message` は押したボタンの下に出す。
43. A8：一文は「optimized for …」（対応表は「optimize for」だった）。
44. A9：エンジン入りの E2E は、blastn + template でなく megablast + template の拒否（「word size must be either 11 or 12」）を確かめる（画面からは前者はもう作れない）。
45. A10：WebKit の 390 px で select の label を、折り返す span で切る（S13 の H2 と同じ）。
46. A11：Job Title は Run を積んだ後も下書きに残す。空白だけは「名前なし」。

### B2：プロット（`logs/report-b2.md`）

47. B2-1：目盛りのラベルは数字だけで、単位は軸の題に（390 px で場所が足りない）。
48. B2-2：重なるラベルは間引く。
49. B2-3：pident が数でない（FakeEngine の `FAKE`）ときは不透明で描く。
50. B2-4：長い辺は subject の方が長くても min(1000, 使える幅)（図は縦に最大 1000 px）。
51. B2-5：長さつきの説明の行は「Plot of …」の下に残す（`results.spec` の正規表現）。
52. B2-6：線でない所のクリックでポップアップを閉じる。
53. B2-7：Graphic Summary は行 12 px・棒 6 px。40 行を超えると query の棒の下で行だけ scroll。
54. B2-8：Graphic Summary の canvas は `role="application"`（矢印キー）。ドットプロットは `role="img"` のまま（コードレビュー L7 を受け、F1 で `application` にした）。
55. B2-9：query の棒は LOSAT の強調色（NCBI の teal ではない）。
56. B2-10：色と不透明度ごとに 1 つの path にし、重なりを濃くしない。**B1-10（判断 66）が覆した**（スクリプトの SVG の線は重なると濃くなる）。

### B1：結果画面（`logs/report-b1.md`）

57. B1-1：Program の欄は「BLASTN (task megablast)」。長さは単位つき。Query / Subject ID は snapshot のレコードから。
58. B1-2：「Results for」と notices はタブの上に置き、どのタブでも見える。
59. B1-3：選んだタブは Run を替えても、Search / Results を切り替えても保つ（`AppView`）。
60. B1-4：Range n はその subject の全 HSP を数える。フィルターは塊を隠すだけで番号は変わらない。a、b は素の整数。
61. B1-5：Range の窓は、選択が窓の外へ出たとき、または subject とフィルターが変わったときだけ中心を移す。`overflow-anchor` はコードで扱う（WebKit）。
62. B1-6：HSP の表のクリックは focus を動かさず Range へ scroll。Next / Previous / First Match と Show alignment は Range の label に focus。
63. B1-7：Firefox の電話の横 scroll（W4 から回された 6a）は本物のマウスの押下では再現せず、Playwright の `click()` が原因。focus を scroll 無しで当てる guard（`focusPressed`）と、マウスの押下で `scrollLeft` が 0 のままの E2E を足した。
64. B1-8：HSP の表は 1280 px に 5 px の間隔で収める。電話は Frames を Query / Subject の前に。Query の範囲は 390 px で切れてよい。
65. B1-9：キューのカードはどこでも同じ形（container query をやめた）。E2E が両タブで同じことを確かめる。
66. B1-10：ドットプロットは HSP ごとに 1 本の線、同じ画素の不透明な線は省く、端は square。B2-10 を覆す（重ねると濃くなるのはスクリプトと同じ。1 本の path にまとめると 2〜3 倍遅かった）。
67. B1-11：仮想化した一覧は描き直したとき選んだ行を見せる。
68. B1-12：計測の写しは 1500 と 3000 だけ。選択の時計は Alignments タブで。`selectWide` は詳細と Alignments の subject を待つ。

### F1：コードレビューの修正（`logs/report-f1.md`）

69. F1-1：Graphic Summary の階級の境界は、印字した精度の半分以内の HSP を E2E が見ない（アプリに hook を足さない）。
70. F1-2：節を読んだ回数は、試験が worker のメッセージを包んで数える。
71. F1-3：失敗した節は「Try again」（`range-retry`）と、見る画面に戻ったときの読み直し。自動の繰り返しはしない。
72. F1-4：M1 の印は要素に直に書く（`data-drawn`、`data-drawn-selected`）。どの時計も見始めの変化から測る。
73. F1-5：画面が落ち着くのを待つのは「一覧で選ぶ」と HSP の並べ替えだけ（Firefox の 3000 写しの並べ替えは、前の layout を測って 38〜42 ms、落ち着かせると 8〜14 ms だった）。
74. F1-6：タンパク質の軸が核酸の軸と並ぶときは 3 nt per aa の重み（BLASTX も同じ規則）。

### F2：画面レビュー 1 回目の修正（`logs/report-f2.md`）

75. F2-1：outfmt 6 の行は scroll でなく折り返す。
76. F2-2：5 px の規則は `plot-scale` の `axisTicks` に置き、単体試験をつけた。
77. F2-3：「drawn at 3 nt per aa」は縮尺どおりのときだけ書く。
78. F2-4：細い図の canvas は舞台の中で広げ、x 軸の題が切れないようにする。
79. F2-5：ポップアップは 8 か所を試し、覆いが少なく動きが小さい所を選ぶ。600 px 未満、または中点を避けられないときは図の下。タップで見える所へ scroll。
80. F2-6：キューの 1 行目は「Run n PROGRAM バッジ」、題・入力・option は操作の下。終わったカードは「Took m:ss」／「Not started」、「Cancel the group」は 3 行目。
81. F2-7：Subject ID の幅は最長の ID から決める（電話は変えない）。
82. F2-8：コミットは 1 つの検証済みの木から hunk で切り、E2E は最後の木で。

### F3：画面レビュー 2 回目の修正（`logs/report-f3.md`）

83. F3-1：電話の query の行は 2 行（46 px）、ID は幅に依らず 12ch 以上、狭い desktop は数に省略記号。`useNarrow`（`RecordList`、`QueryPicker`）。
84. F3-2：Run の select とパラメータの select の右に 24 px の余白（WebKit はさらに早く切れる）。
85. F3-3：「Cancel the group」は 3 行目のまま。
86. F3-4：Combined / Separate の箱を 6 px 内に寄せる。
87. F3-5：結果の見出しの scroll margin は 8 px。
88. F3-6：コミットは hunk で切り、どの木も `npm run check` が通る。

## 実測

計測は `tests/e2e/results-measure.spec.ts`。1 回の暖機と 3 回の計測の中央値で、ページの中の `performance.now()` による、操作から結果を示す最初のフレームまでの ms（分解能は 1 フレーム、約 17 ms）。W4 と同じ作り方で、括弧の中は W4 のゲート記録の「実測」の数値。

### ゲートの実行（`f097595d`）と W4

`run-20261009T184556Z/measure/`。**「ドットプロットを出す」の時計は描画を含まない**（コードレビュー M1。コンポーネントが frame の poll より後に mount するので、時計が描画の前に止まる）ので、W4 の値（描画を含む）と比べられない。Graphic Summary も同じ。他の時計は W4 と同じ意味。

**多数の query**：

| query | ブラウザ | 検索（s） | HSP | 開く | 最初の HSP まで | 終わりまで scroll | ID で探す（選択と詳細まで） | 絞り込みを外す | 別の query をクリック | With hits only | 200 subject の query をクリック | subject の並べ替え（最大） |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10,000 | Chromium | 0.9（0.8） | 14,981 | 41（39） | 53（43） | 1（1） | 11（11） | 12（12） | 7（8） | 11（12） | 16（9） | 12（12） |
| 10,000 | Firefox | 3.3（2.8） | 14,981 | 83（67） | 90（74） | 1（1） | 10（10） | 10（10） | 18（19） | 7（10） | 7（13） | 10（13） |
| 10,000 | WebKit | 1.4（1.4） | 14,981 | 65（61） | 101（85） | 8（10） | 40（7） | 7（15） | 39（15） | 8（7） | 24（47） | 9（11） |
| 100,000 | Chromium | 5.3（6.7） | 149,828 | 419（300） | 371（313） | 1（1） | 11（12） | 15（16） | 12（18） | 13（15） | 14（9） | 12（12） |
| 100,000 | Firefox | 25.5（22.4） | 149,828 | 676（590） | 601（660） | 1（1） | 10（10） | 18（17） | 20（21） | 17（16） | 9（8） | 10（11） |
| 100,000 | WebKit | 8.9（7.1） | 149,828 | 489（521） | 537（502） | 6（11） | 48（34） | 18（18） | 16（15） | 18（17） | 18（15） | 10（13） |

**1 組の多数の HSP**：

| 写し | ブラウザ | 検索（s） | HSP | 開く | ドットプロットを出す | 拡大 | 縮小 | Zoom to HSP | 全体 | n で選ぶ | クリックで選ぶ | 一覧で選ぶ | HSP の並べ替え |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,500 | Chromium | 4.3（3.8） | 2,995 | 294（185） | 10（42） | 4（28） | 4（38） | 3（40） | 3（38） | 13（92） | 29（99） | 23（16） | 7（11） |
| 1,500 | Firefox | 24.5（18.4） | 2,995 | 436（326） | 9（96） | 36（34） | 37（46） | 41（45） | 32（47） | 20（109） | 31（130） | 24（15） | 8（10） |
| 1,500 | WebKit | 4.3（4.3） | 2,995 | 146（139） | 9（11） | 8（10） | 7（12） | 26（4） | 9（13） | 45（151） | 50（229） | 43（14） | 11（12） |
| 3,000 | Chromium | 16.9（15.7） | 5,993 | 871（754） | 10（86） | 3（55） | 9（79） | 4（72） | 4（71） | 20（169） | 28（163） | 36（15） | 10（9） |
| 3,000 | Firefox | 92.1（76.0） | 5,993 | 1688（1308） | 11（179） | 36（67） | 34（89） | 43（89） | 31（87） | 15（186） | 30（199） | 43（14） | 10（9） |
| 3,000 | WebKit | 失敗：結果が 512 MB の保存の上限を超えた | | | | | | | | | | | |

- 開く時間が W4 より 16〜59% 遅い（1 組 1500 写しの Chromium で 294 対 185、10 万 query の Chromium で 419 対 300）。一覧で選ぶ時間も 23〜43 ms と W4（14〜16 ms）より遅い。次の節で原因を調べた。
- 3000 写しの WebKit と、4000 写しの 3 ブラウザは見込みどおりの失敗：WebKit は結果が 512 MB の保存の上限を超えた（OPFS の無い WebKit は結果をメモリに置く）、4000 写しはエンジンのメモリ不足。

### 開く時間と一覧で選ぶ時間の調べ（F1）

同じ機械で W4 の木（`50d57bba`）を作り直し、HEAD（`d949ee7a`）と 10 回、交互に測った（`logs/measure-f1/ab/`、タスクフォルダ）。括弧の中はゲートの実行（`f097595d`）の値。

| 場面 | ブラウザ | HEAD | W4 の木 | ゲートの実行 |
|---|---|---:|---:|---:|
| 開く、1 組 1500 写し | Chromium | 219 | 214 | 294 |
| 開く、1 組 1500 写し | Firefox | 371 | 322 | 436 |
| 開く、1 組 3000 写し | Chromium | 792 | 817 | 871 |
| 開く、1 組 3000 写し | Firefox | 1274 | 1434 | 1688 |
| 開く、1 万 query | Chromium | 40 | 41 | 41 |
| 開く、1 万 query | Firefox | 71 | 88 | 83 |
| 開く、10 万 query | Chromium | 347 | 395 | 419 |
| 開く、10 万 query | Firefox | 560 | 659 | 676 |
| 一覧で選ぶ、1500 写し | Chromium | 17 | 14 | 23 |
| 一覧で選ぶ、1500 写し | Firefox | 20 | 14 | 24 |
| 一覧で選ぶ、3000 写し | Chromium | 17 | 26 | 36 |
| 一覧で選ぶ、3000 写し | Firefox | 28 | 14 | 43 |

- **開く時間**：Chrome の trace（1 組 1500 写し）で、main thread は約 358 ms のうち約 300 ms が idle。Data worker が HSP のレコードと outfmt 6 を読む 250〜330 ms が支配的で、W4 から変わっていない。W4 の木を作り直して測ると 270〜309 ms で、W4 のゲートの 185 ms は速い側の揺れだった。W4b の増分は main thread で 1 組 +6 ms、10 万 query +10 ms。query ごとの表を、開くときに全 HSP の map を作らず必要になったとき作るようにし（10 万 query の main thread の task が 105 ms から 47 ms）、subject の見出しの読みを 2 回から 1 回にした。ゲートの増加は大部分が実行ごとの揺れで、残る差は Firefox の 1 組 1500 写しの +15%（標本は 325〜391 ms 対 292〜389 ms で重なる。範囲内）だけ。
- **一覧で選ぶ時間**：原因は、長い outfmt 0 の節が届くたびに 10〜25 ms の layout を強いることだった（trace で 134 ms のうち 98 ms）。1 フレームに届いた節をまとめて適用するようにした（`81f25039`）。Firefox の 3000 写しは 1 フレーム分遅いまま（28 対 14）。
- 注意：Chromium の 3000 写しで長い節がまとめて届くと、時計が止まった後、同じ layout の仕事がまとめて 1 フレームに入り、最大約 100 ms かかる。仕事の量は同じで、1 フレームに集まった。

### 時計の意味（M1 の後）

`9d95c5aa` の後、ドットプロットと Graphic Summary の描画の後に `data-drawn`（選択は `data-drawn-selected`）を要素に書き、「出す」「拡大」「縮小」「Zoom to HSP」「Whole sequences」はその印を待つ。「n で選ぶ」「クリックで選ぶ」は選択の印、ポップアップはそのフレームを待つ。つまり show と zoom の時計は**描画を含み**、W4 の時計（描画を含む）と比べられる。見始めは、どの時計も操作の前の状態からの変化で決める。Chromium の拡大と縮小の時計（2〜9 ms）は GPU の描画を含まない（`paintMs` が 29〜51 ms でそれを覆う）。WebKit は最初のフレームが早く、canvas の描画が次のフレームに出るので、`paintMs` の列も記録する（括弧の中が描画のフレームまで）。「n で選ぶ」「クリックで選ぶ」は Playwright の往復を含む上限（W4 の判断 14）。

### レビューの後の計測（`087eaa13`）と W4

`run-20261009T211154Z-measure/measure/`。ゲートの実行とは時計の意味が違う（上）ので、ドットプロットの時計は W4 と比べられる。

**多数の query**：

| query | ブラウザ | 検索（s） | HSP | 開く | 最初の HSP まで | 終わりまで scroll | ID で探す（選択と詳細まで） | 絞り込みを外す | 別の query をクリック | With hits only | 200 subject の query をクリック | subject の並べ替え（最大） |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10,000 | Chromium | 0.8（0.8） | 14,981 | 38（39） | 41（43） | 15（1） | 12（11） | 12（12） | 8（8） | 13（12） | 9（9） | 12（12） |
| 10,000 | Firefox | 2.8（2.8） | 14,981 | 64（67） | 76（74） | 1（1） | 10（10） | 10（10） | 20（19） | 11（10） | 7（13） | 11（13） |
| 10,000 | WebKit | 0.9（1.4） | 14,981 | 45（61） | 77（85） | 10（10） | 8（7） | 6（15） | 15（15） | 6（7） | 40（47） | 16（11） |
| 100,000 | Chromium | 6.7（6.7） | 149,828 | 292（300） | 266（313） | 1（1） | 12（12） | 13（16） | 14（18） | 12（15） | 10（9） | 12（12） |
| 100,000 | Firefox | 21.4（22.4） | 149,828 | 500（590） | 524（660） | 1（1） | 10（10） | 16（17） | 21（21） | 16（16） | 9（8） | 11（11） |
| 100,000 | WebKit | 4.8（7.1） | 149,828 | 415（521） | 438（502） | 6（11） | 36（34） | 14（18） | 15（15） | 14（17） | 43（15） | 15（13） |

**1 組の多数の HSP**：

| 写し | ブラウザ | 検索（s） | HSP | 開く | ドットプロットを出す | 拡大 | 縮小 | Zoom to HSP | 全体 | n で選ぶ | クリックで選ぶ | 一覧で選ぶ | HSP の並べ替え |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,500 | Chromium | 3.8（3.8） | 2,995 | 181（185） | 20（42） | 3（28） | 3（38） | 6（40） | 7（38） | 7（92） | 14（99） | 9（16） | 13（11） |
| 1,500 | Firefox | 17.9（18.4） | 2,995 | 308（326） | 37（96） | 28（34） | 30（46） | 40（45） | 28（47） | 13（109） | 30（130） | 20（15） | 11（10） |
| 1,500 | WebKit | 3.8（4.3） | 2,995 | 131（139） | 10（11） | 16（10） | 6（12） | 23（4） | 7（13） | 15（151） | 47（229） | 8（14） | 6（12） |
| 3,000 | Chromium | 15.3（15.7） | 5,993 | 672（754） | 20（86） | 4（55） | 9（79） | 7（72） | 6（71） | 16（169） | 14（163） | 16（15） | 13（9） |
| 3,000 | Firefox | 75.1（76.0） | 5,993 | 1109（1308） | 38（179） | 27（67） | 34（89） | 25（89） | 27（87） | 15（186） | 30（199） | 26（14） | 10（9） |
| 3,000 | WebKit | 失敗：結果が 512 MB の保存の上限を超えた | | | | | | | | | | | |

- WebKit の 1500 写しの描画のフレームまでの時間（`paintMs`）：ドットプロットを出す 61、拡大 60、縮小 35、Zoom to HSP 57、全体 39 ms（W4 は 224・59・120・103・115）。
- 10 万 query の Chromium の JS heap（開く前、最初に開いた後、繰り返しの後）：14.4・43.0・43.3 MB（ゲートの木は 14.4・49.6・49.9、W4 は 14.4・49.6・49.8）。Outputs で outfmt 6 の全文を出すのは Chromium 2,450 ms、Firefox 26 ms、WebKit 4 ms（ゲートの木は 2,979・51・8 ms、W4 は Chromium 2.7 秒、Firefox と WebKit は 1 フレーム）。
- 5,993 HSP のドットプロットは、出すのに Chromium 20 ms・Firefox 38 ms、拡大・縮小・移動は 4〜34 ms。W4（出すのに 86〜179 ms、拡大・縮小・移動 55〜89 ms）より速い。
- 1 万 query の Chromium の「終わりまで scroll」は 15 ms（W4 は 1 ms）。1 フレーム分で、10 万 query の run は 1 ms。
- 検索の時間と HSP の数は W4 と同じ水準（検索の揺れは同じ機械で並行した処理による）。
- 最後の木（`492fad48`）では測り直していない（上の「レビューの後の実行」）。

## 独立レビュー（コード）

agent `losat-reviewer`（code review の役。このセッションは clone の外で起動したので、定義を general-purpose の agent に渡した）が `ff9ebc49..f097595d` の `web/app`（44 ファイル）を読んだ。結論は「妨げになる指摘は無い」。指摘と対応：

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 計測の時計「ドットプロットを出す」と Graphic Summary が描画を含まない（poll が mount より先に動く）。S13 の 86〜179 ms（描画を含む）と比べられない | `9d95c5aa`：描画の後に `data-drawn` を書き、時計がそれを待つ（判断 72）。Firefox の「出す」が 42〜44 ms（描画を含む）になった |
| L2 | E2E が Graphic Summary の階級を丸めた outfmt 6 の文字列で決めるが、コンポーネントは原値で分ける。階級の境目で食い違い得る（39.96 は「40.0」と印字され、棒は < 40） | `81f25039`：印字の精度の半分以内の HSP は見ない（判断 69） |
| L3 | `expectRanges` が空の節でも通る（`toContain('')`） | `81f25039`：節が `/^ Score = /` に合うことを求める |
| L4 | 「節は見える所だけ読む」の試験が、読み終えた数を数え、要求した数を数えない | `81f25039`：`readOutputRange` の要求を worker のメッセージを包んで数える（判断 70） |
| L5 | 失敗した節を Alignments が読み直さない（`null` を保存する） | `81f25039`：「The alignment could not be read.」と「Try again」、見る画面に戻ったときの読み直し。エンジン入りの E2E が両方の経路を確かめる（判断 71） |
| L6 | ドットプロットの端：描画範囲の判定が線の幅を数えず、選ぶ処理が画面の外の線を取り得る | `57bc7f91`：`lineNearView`（1.5 px の余白）、`clipToBox`・`pickLine`。単体試験 |
| L7 | ドットプロットの canvas が `role="img"` なのに n・p・Enter を受ける。ブラウズモードのスクリーンリーダーが届かない | `57bc7f91`：`role="application"`、キーを名前に書く。E2E が `getByRole('application')` と `getByRole('dialog')` を使う（判断 54） |
| L8 | 判断 26（3 nt per aa）が HEAD で未実装 | `57bc7f91`：`codonWeights` と `plotSize`、`data-plot` 属性。E2E が縦横比を測る（判断 26、74） |

確かめて指摘の無かったもの：検索画面の印と argv が同じ関数を使うこと、Restore が Algorithm parameters の節だけを消すこと、Discontiguous Word Options で隠れた欄が argv に入らないこと、task の radio が既定のとき何も書かないこと、Job Title が argv に入らないこと、Range の番号と窓、Graphic Summary の階級が半開区間であること、ドットプロットの目盛り・縮尺・向き・色・不透明度・ポップアップの中身、BLAST の値を TS で計算・整形していないこと、画面の文字がすべて英語であること、FakeEngine の印・ネットワークと URL の規則・層の規則。

修正（F1〜F3）の後にコードレビューはもう一度行っていない。修正には試験をつけ、画面レビュー 2 回と 3 回が見える変化を確かめた。

### 2 回目（修正の後）

agent `losat-reviewer`（code review の役）が修正の差分 `f097595d..492fad48`（`web/app`、23 ファイル）を読んだ（記録は作業フォルダの `reviews/code-review-2.md`）。結論は「妨げになる指摘は無い」。1 回目の 8 件（M1 の時計、L2〜L8、判断 26 の 3 nt per aa）は報告どおり直っている。確かめて指摘の無かったもの：query ごとの表を必要なときに作る（`rowOf`、`selectHsp` は索引の前の query の HSP にも効き、古い表は残らない、Run を開く主スレッドに HSP ごとの仕事が無い）、見出しの読み取りの共有（`sIdx` の鍵、失敗は読み直せる）、節を 1 フレームでまとめて当てる（重複の読み取り無し、Run の切り替え、読み直しは繰り返さない）、ドットプロットの比と間引きと選び方、ポップアップの置き場所、TBLASTN の frame、計測の時計（それぞれ名前どおり）、規則（`web/` で BLAST の値を計算しない、層、英語、FakeEngine、ネットワーク）。低 2 件は下の「残件と注意」。

## 画面レビュー

agent `losat-reviewer`（screen review の役）が、画面の記録を対応表、NCBI の参照画面、前の記録と並べて見た。3 回とも「合格（High 無し）」だが、レビュー自身が「完全な合格ではない」と書いた回があり、中程度の指摘は次の回の前に直した。

**1 回目**：ゲートの記録 150 枚（`f097595d`、`$BUILD_ROOT/s13b-gate-screens/`）を S13 の記録（状態 01〜18）と NCBI の参照画面と並べた。検索画面の順と言葉は NCBI のとおり、結果画面は見出しの塊・タブ・道具の行が対応表のとおり。NCBI のブランドは無く、ボタンは「Run LOSAT」、390 px の 75 枚はどれも 390 px 幅。S13 の 2 回目の持ち越し 6 件のうち、3 件は直り（電話のドットプロットの説明、x 軸の端のラベル、TBLASTN の記録）、3 件は一部（Firefox の電話の記録はアプリが直ったが記録に残る、電話の表が列の途中で終わる、キューのカードの変化）。

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 1280 px で Range の塊の outfmt 6 の行が bit score を失う（S13 からの退行。desktop-09・21） | `d4c367a2`：行を折り返して全欄を示す（判断 75）。E2E が 1280 px で収まりを測る |
| M2 | 1 つのドットプロットで小格子が片方の軸だけに出る（個数が 50 を超えると省く規則） | `e81f791f`：数の上限をやめ、5 px の規則を両軸に（判断 76） |
| M3 | 電話の Y 軸の題が、長い ID で単位を失う | `e81f791f`：ID を「…」で切り、単位を残す。E2E が `data-titles` を確かめる |
| 4 | TBLASTN のドットプロット（既知。判断 26）。電話では「Axes not to scale」が残る。題に縮尺を書く | `57bc7f91`、`e81f791f`：「(aa; drawn at 3 nt per aa)」、E2E が desktop に注記が無く電話にあることを確かめる（判断 77） |
| L5 | Firefox の電話の記録 07・08 で HSP の表が横に流れている（Playwright の scroll-into-view） | `0659584c`：見える所をマウスで押す（アプリ側は B1-7 で直っていた） |
| L6 | 電話の表が列の途中で終わる | 受け入れ（判断 27） |
| L7 | キューのカードの折り返しが状態で違う（「Open results」が 1 行を取る、Cancel の位置） | `9f99f956`：どのカードも同じ形（判断 31、80） |
| L8 | Descriptions：Description が切れるのに Subject ID に余りがある、HSPs の位置 | `596c5aa6`：幅を Description に、HSPs は E value の後（判断 28） |
| L9 | TBLASTN の frame の「–/+2」が TBLASTX の「-2/+2」と並ぶと鎖に読める | `0fe326b8`：「Subject frame」と「+2」だけ（判断 29） |
| L10 | ポップアップが自分の HSP の線を覆う | `e81f791f`：線の脇、電話では図の下（判断 32、79） |
| L11 | Firefox でポップアップの focus が見えない | `e81f791f`：`:focus` にも輪郭。E2E が 3 ブラウザで確かめる |
| L12 | 単位の言葉（軸は bp、説明は nt）、「1 subj., 1 HSPs」 | `e81f791f`、`596c5aa6`、`087eaa13`：対応表を直し（判断 30）、数を言葉に（「1 subject, 1 HSP」） |
| L13 | Run details の見出しがタブの下線に触れる、電話の入力欄の右端 | `596c5aa6` |

**2 回目**：`087eaa13` の記録 150 枚（`$BUILD_ROOT/s13b-screens-f2/`）をゲートの記録と並べた。1 回目の 13 件はすべて直っていた（L6 は受け入れ）。

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 電話の「Results for」で ID が切れて 2 つが同じに見える（L12 の言葉の長さによる退行。`BDT6256…` が 2 つ） | `06c1a290`：長さと数を ID の下の 2 行目に、ID は 12ch 以上（判断 33、83）。E2E が 1280 px と 390 px で ID が切れていないことを測る |
| L1 | Chromium の Run の select の字が矢印の下に入る | `42742404`：右に 24 px（判断 84） |
| L2 | 終わったキューのカードが 1 行増え、ページが 134〜183 px 伸びた | `63098d40`：操作を段階の行の右端に（判断 34） |
| L3 | 入力のカードが欄より 6 px 外に出る | `a7b4e0bf` |
| L4 | TBLASTN の小格子が密（情報） | 対応表の規則のまま。もっと粗くしたければ最小の間隔を 8 px などに |

**3 回目**：`492fad48` の記録 150 枚（`$BUILD_ROOT/s13b-screens-f3/`。最後の after-review 実行の記録は同じ木の `$BUILD_ROOT/s13b-final-screens/`）。2 回目の 4 件はすべて直り、状態 18 は「Results」の見出しの上に 8 px ある。01・02・07・08・14・17〜23 に退行は無く、150 枚はどれも 1280 px または 390 px 幅。新しい指摘は低 1 つで、残りは参考。

| 指摘 | 内容 | 対応 |
|---|---|---|
| 低 1 | 電話で入力のカードの長いファイル名が 2 行に伸びる（`NZ_CP006932.faa` の「Remove」が 2 行目に。カードが約 48 px 伸びる。L3 の副作用） | **S14 へ**（判断 35）。任意の直し：800 px 以下でサイズを名前の下に折り返す |
| 参考 O1 | 電話の Chromium と Firefox で「More dissimilar sequences」の radio が 14 px に縮む（ゲートの記録から。`.radio input` に `flex: none` が無い）。ラベルの行は 24 px 以上 | **S14 へ**（判断 35） |
| 参考 I1〜I3 | 全パラメータの select が 18 px 広い、WebKit の Run の select は label が 20 px 早く切れる、電話の Records の行の折り返しが変わった | 直さない |

3 回とも撮らなかった状態：Alignments の「Try again」（E2E は確かめた）、TBLASTN のポップアップの subject frame、「Cancel the group」、801〜1000 px の幅。任意の状態 24 で「Try again」を記録できる。

## 見つけて直したもの

- `c5019bc8`：狭い TBLASTN の電話の図（120 px）が、5 px より狭い小目盛りの格子で灰色に塗れていた。
- `f097595d`：重なる多数の線を 1 本の path にまとめると 2〜3 倍遅かったので、HSP ごとに 1 本の線にした（5,993 HSP のズームが Firefox で 190〜274 ms から 30〜43 ms）。判断 66。
- `7f58c537` の中：Firefox の電話の横 scroll は Playwright の `click()` の scroll で、本物のマウスの押下では起きない（判断 63）。WebKit の `overflow-anchor` はコードで扱った（判断 61）。
- 検索画面の担当（A）が見つけた、両方のタブの主な格子を `minmax(0,1fr) minmax(260px,300px)` にすると 1280 px で結果の列が約 40 px 狭くなり、TBLASTN と TBLASTX の HSP の表が 2 px はみ出すこと：結果画面の担当（B1）が表を 5 px の間隔で収めて直した（判断 64）。
- `778989b2`：Job Title のキュー上の位置の試験の競合（Run が終わる前に読んでいた）。
- 3 回の画面レビューとコードレビューの指摘（上の表）。
- F3 の最初のエンジン入りの実行で 4 件落ちた：待ちのカードの Cancel が折り返した、Firefox の電話の状態 18 で見出しが 0.992 しか見えなかった（`cc3aa24f`）。直して 87 件が通った。

## 合流

`feature/losat-web-gui` への merge は、このセッションでは行わない。エンジン側が W4 と W4b を合わせて行う。merge のときにエンジン側が行うこと：

1. `origin/feature/losat-web-gui-app` を `feature/losat-web-gui` に merge する（起点 `ff9ebc49` を含み、W4 と W4b の両方が入る。W4 のゲート記録の「合流」の 1.〜3. がまだなら、同時に行う）。
2. README の表に S13（完了。[W4](../losat_web_w4/README.md)）と S13b（完了。この記録）の行を置く。
3. 計画 §7 に S13b「画面を NCBI BLAST Web に寄せる」（段階 W4b、S13 と S14 の間、アプリ側）を足す。計画 §0.4 に、保守者の指示 2 件を DW として記す：2026-10-09「NCBI BLAST ウェブサイトの GUI にデザインを寄せてほしい」（範囲は配置・節の順・言葉、ブランドは写さない。ドットプロットは `blast2dotplot.py` を参考に）、2026-10-10「Run LOSAT」（ボタンの名前）。判断 1〜88（保守者の確認待ち）は、W4 の判断と同じ扱いで、推奨案で進めたと記す。

## エンジン側への申し送り

- **E2j（S13+）の後のアプリ**：S13b は E2j が merge される前に終わったので、E2j の指示書の 6.（アプリへの申し送り）を行っていない。S14 が最初に行う（S14 の指示書の開始の条件）。列の置き場所は対応表に決めてある：Descriptions の Max Score は Score (bits) の位置、Total Score・Query Cover は Score と E value の間、Per. Ident と Acc. Len は E value の後。HSPs は E value の後のまま動かない（判断 28）。Filter Results の Query Coverage、BLASTN の 1 文字の HSP の鎖（ドットプロットの灰色の線と凡例、`detail-strand-note`）もそのときに足す。
- **検証バッジ**：既定の BLASTN は、E2j が既定の outfmt 7 の fixture を足すまで「outside」のまま（W4 の判断 3）。
- **エンジンのメモリ不足**（W4 の申し送りの続き）：20 nt の単位を 4000 写しの自己検索（約 8,000 HSP）が 3 ブラウザで `memory allocation of 49881 bytes failed`。W4b でも再現した（ゲートの実行、`087eaa13` の計測）。
- アプリ側は `describe.json` とエンジンの検査の文言を、画面に出す以外には解釈していない（検索画面の欄は `describe` から、task の radio の選択肢も `describe` の choices）。

## S14 への申し送り

[S14 の指示書](../../losat_web_gui_sessions/session_s14_w5_extraction_candidates.md)の「S13b（W4b）から引き継ぐこと」に書いた（結果画面の新しい構成、動いた test ID、検索画面の順と Job Title、ドットプロット、対応表の使い方、画面の記録の場所 `$BUILD_ROOT/s13b-final-screens/`、ゲートの script と計測の数値、残った低い指摘 2 件）。

## 残件と注意

- **TBLASTN・BLASTX の未決**：BLASTX は SX まで画面で動かせないので、3 nt per aa の規則（判断 74）は単体試験で確かめただけで、実際の画面の記録は無い（SX で BLASTX を入れるとき、ドットプロットの subject の軸と frame の「Query frame」の見え方を確かめる）。TBLASTN の電話では、query の辺が 120 px に満たないので「Axes not to scale」が残る（画面レビューは妥当と判断）。TBLASTN の desktop は小格子が query の軸で 5.2 px 間隔と密（もっと粗くしたければ最小の間隔を 8 px などにする。画面レビュー 2 回目の L4）。
- **Firefox の開く時間 +15%**（1 組 1500 写し。371 ms 対 W4 の木の 322 ms、交互に 10 回。標本が重なり、揺れの範囲内）。Firefox の 3000 写しの「一覧で選ぶ」は 1 フレーム分遅い（26〜28 ms、W4 は 14 ms）。追わない。
- **計測の対象の木**：計測は `087eaa13`。`492fad48` では測り直していない。S14 以降の計測は、今の木から取り直す。
- **S13 から残る残件**：エンジンのメモリ不足（上）、WebKit の保存の上限（OPFS の無い WebKit は結果をメモリに置き、512 MB を超える Run は「Not enough temporary storage」で失敗する。W1 の設計どおりで、V-MOB の実機（S17）で確かめる）、`runInputs` の解放（Run の削除が入るまで解放しない。W4 の判断 5）、iOS Safari での確認（V-MOB、S17）。
- **画面レビューの残り**：低 1（電話の入力のカード）と参考 1（電話の radio の 14 px）。直す場合は電話の 02・03・06 を 3 ブラウザで撮り直す。撮らなかった状態（Alignments の「Try again」、TBLASTN のポップアップ、「Cancel the group」）は状態 24 以降で足せる。
- **対応表の外**：撮らなかった NCBI の画面（Taxonomy タブ、Download と Select columns のメニュー、Edit Search の画面）の採否は、指示書と下調べの文書による。Edit Search は S15。
- **列定義表**：`results_columns.md` の列の表は W4 の並びのまま（結果画面の列の並びの正本は対応表）。E2j の後に両方を合わせる。
- **コードレビュー 2 回目の低 2 件**（`f097595d..492fad48`、妨げになる指摘無し）：(1) 「Try again」は選んでいない Range の塊だけで、選んだ HSP の詳細の読み取りが失敗したときは同じ HSP を選び直しても読み直さない（W4 からの振る舞い）。また、失敗した塊は、同じ節が詳細として読めた後も、見えたままなら誤りを示し続ける（`AlignmentsView.vue:115,188-191,282,316-326`、`results.ts:308`）。(2) 横の広い組（例：1280 px で 100 kbp 対 5 kbp の BLASTN、図の高さ約 180 px）で、最初のクリックでポップアップが図の下に回ったとき、見える所まで scroll しない（`DotPlot.vue:531,575,625,657`。電話は scroll する）。どちらも S14 で直す（S14 の指示書の残件）。
- 画面の記録（PNG）はリポジトリに入れない：ゲートの 150 枚は `$BUILD_ROOT/s13b-gate-screens/`、レビューの後の 150 枚は `s13b-review-screens/`（`087eaa13`）、`s13b-screens-f2/`、`s13b-screens-f3/`、最後の after-review の 150 枚は `$BUILD_ROOT/s13b-final-screens/`（`/home/kawato/.cache/losat-work/s13b-final-screens/`。S14 はこれと比べる）。SHA-256 は実行記録の `screens.sha256`（`s13b-screens-f2` と `s13b-screens-f3` は実行記録を残していない）。
- NCBI の参照画面（174 個）もリポジトリに入れない。SHA-256 は [`ncbi_reference.sha256`](ncbi_reference.sha256)。
