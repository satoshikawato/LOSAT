# LOSAT Web W4（Session S13）ゲート記録

- 段階：W4 結果画面（[総合計画書](../../losat_web_gui_plan.md) §7 の S13、[指示書](../../losat_web_gui_sessions/session_s13_w4_results_ui.md)）
- ブランチ：`feature/losat-web-gui-app`（アプリ側。worktree `$WORK_ROOT/.worktrees/web-gui-app`、Linux の clone）。起点は `f3048ffd`（`feature/losat-web-gui` の DW-23 の記録）。エンジン（`LOSAT/`、`web/adapter/`）は W3 のゲートの木 `4b3418992` から変わっていないので、reactor とネイティブの CLI は W3 のもの（`app-s12-reactors`、ネイティブ sha256 `f0b8916b…`）を使った
- 実行記録：[`run-20261009T122850Z/`](run-20261009T122850Z/)（commit `b2d8fb24` の木。作成後は書き換えない）。レビューの後の実行 [`run-20261009T135520Z-after-review/`](run-20261009T135520Z-after-review/)（`7ee77f8d`）と [`run-20261009T142446Z-after-review/`](run-20261009T142446Z-after-review/)（`b565eab8`、最後の木）、計測だけの実行 [`run-20261009T140235Z-measure/`](run-20261009T140235Z-measure/)（`7ee77f8d`、一部の場面を落とした。判断 22）と [`run-20261009T141300Z-measure/`](run-20261009T141300Z-measure/)（`50d57bba`）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)、再現は [`run_gate.sh`](run_gate.sh)
- 判定：**完了条件を満たした**。ゲートの実行（`b2d8fb24`）のすべての段階、コードレビュー（妨げになる指摘無し。M1・M2・L1〜L6 を直した）、画面レビュー（1 回目不合格の H1・H2 などを直し、2 回目合格）、レビューの後の実行、計測だけの実行が通った。下の「完了条件と結果」
- 保守者の判断待ち：判断 1〜23（下の「判断」）は推奨案で進めた（Owner-delegated、2026-09-29 の常設の指示、2026-10-07 に再掲）。とくに 1（集約列の採用とエンジン側の S13+ の挿入）と 2（BLASTN の 1 文字の HSP の鎖を ABI に足す）は、エンジン側の計画と README の表に行を足す（下の「エンジン側への申し送り」）
- `feature/losat-web-gui` への merge：このセッションでは行わない。エンジン側の SF が `.worktrees/web-gui` を使っているので、merge はエンジン側に任せる（下の「合流」）

## コミット

| コミット | 内容 |
|---|---|
| `e757434b` | 文書：結果画面の列定義表 [`docs/web/results_columns.md`](../../web/results_columns.md)（作業 1） |
| `f27132fc` | アプリ：Run の結果を HSP の表とバイトの範囲で読む（`readHitTable`・`readOutputRange`）。FakeEngine が検索の形の出力を書く |
| `691ce4d3` | アプリ：検証バッジの判定表をビルドのときに生成する（`build/verification.ts`、計画 §6.1） |
| `1cc9cfb5` | アプリ：`ResultsBrowser` と Run の結果の索引（`domain/result-index.ts`） |
| `be55e1eb` | アプリ：結果画面（`ResultsPanel` を置き換え、キューから結果を開く） |
| `e649391f` | アプリ（試験だけ）：結果画面の単体試験と、全 outfmt 0 fixture の HSP レコードと出力の対応の試験 |
| `1ddb0ced` | 文書：W4 の `run_gate.sh` と S13+（E2j）の指示書の下書き |
| `19154366` | アプリ：W4 のファイルの実行権を外す（内容は変えない） |
| `8f826809` | アプリ：frame は翻訳する配列にだけ示す（BLASTP が +1/+1 を示していた。エンジン入りの E2E が見つけた） |
| `028babc2` | アプリ（試験だけ）：結果画面の E2E と、既存の spec を新しい結果タブに合わせる |
| `e200e60b` | 文書：S14 の指示書の下書き（開始の条件と S13 からの引き継ぎ） |
| `bd7e4900` | アプリ：結果画面の表の見出しを行の列に揃える（画面の記録で見つけた） |
| `740a7225` | アプリ（試験だけ）：結果画面の計測と画面の記録（状態 07〜16） |
| `b2d8fb24` | アプリ：並べ替えは一覧の先頭を見せる。query のフィルターと別の Run は、隠れた選択や古いフィルターを残さない（判断 9〜11）。**ゲートの木** |
| `d7ade9e9` | アプリ：コードレビューの指摘（M1・M2・L1〜L4。判断 12・13・15） |
| `d78f1448` | アプリ（試験だけ）：計測の時計を名前どおりにする（L5・L6、判断 14）。`run_gate.sh` に `LOSAT_WEB_GATE_STEPS=measure` |
| `3c601a46` | 文書：`results_columns.md` に判断 12 の補足 |
| `057b18f1` | アプリ：画面レビューの指摘（H1・H2・M1〜M3・L1〜L8。判断 16〜21） |
| `7ee77f8d` | アプリ（試験だけ）：画面の記録の 01 を option の検査の後に、状態 17（翻訳する program のドットプロット）と 18（「Open results」の直後の画面）を足す |
| `50d57bba` | アプリ（試験だけ）：ドットプロットの計測のクリックを、選ばれた線から離れた HSP に向ける（判断 22。図の描き方と選び方は変えない） |
| `b565eab8` | アプリ：電話の Subject 一覧で Subject の列を 12ch 以上に保つ（画面レビューの 2 回目の中 1） |
| このゲート記録を含むコミット | 文書：このゲート記録、実行記録、S13b の指示書（保守者の指示、下の「保守者の指示」）、NCBI の画面の下調べ、S14 の指示書の数値と画面の記録の場所 |

エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表は、アプリ側のブランチでは変えていない（README 規則 1）。`docs/web/` は、このセッションが作った `results_columns.md` だけを変えた。依存ライブラリは増やしていない。

## 完了条件と結果

| 完了条件（計画 §7 の S13） | 結果 | 証拠 |
|---|---|---|
| 対応済みの全 program の E2E | 通過。BLASTN・BLASTP・TBLASTN・TBLASTX（BLASTX は SX まで対象外。W3 の判断 2 のとおり画面で「not available」）。3 ブラウザ、FakeEngine とエンジン入りの両方のビルド、エンジン入りは 2 回の繰り返しも。program ごとに、キューから Run を開き、Subject と HSP の一覧のセルがその HSP の outfmt 6 の欄と等しく、詳細の行が outfmt 6 の行で、見出しと節が outfmt 0 にそのまま現れ（`Strand=`・`Frame=` を一覧と照らす）、単位と frame、Run の詳細の CLI コマンドが Outputs の表示と等しく、検証バッジが生成した表どおり（エンジン：既定の option で BLASTN は「outside」（E2j まで）、BLASTP・TBLASTN・TBLASTX は「certified」。FakeEngine：「development」） | `npm-e2e-engine.log`、`e2e-results-repeat.log`、`npm-e2e-fake-engine.log`（`tests/e2e/results.spec.ts`） |
| HSP と行・節との対応の試験 | 通過。`LOSAT/tests/outfmt0_manifest.tsv` の BLASTN・BLASTP・TBLASTN・TBLASTX の全 fixture（165 の検索、30,467 HSP、outfmt 0 の節 30,027、outfmt 0 に現れない HSP 440）を serial の reactor で Node 上で実行し、出力が NCBI の凍結バイト（196 が byte 一致、承認済みの `-db_gencode` の例外 2 が NCBI のデータベース検索とデータベースの行の外で一致、field list の 3 は比べない。判断 6）、`out6` の範囲が outfmt 6 を記録の順に隙間なく覆い座標が一致、`out0` の範囲がその HSP の節（最初と最後の Query・Sbjct の座標が一致）、見出しは節の前、節が無いのは見出しが無いときだけ、outfmt 0 の節と見出しの数が記録の指す数と等しい（`d7ade9e9` で追加、L4）、結果の索引がすべての HSP をまとめる | `unit-cases.log`（`tests/unit/hsp-correspondence.test.ts`） |

### ゲートの実行（`b2d8fb24`）

`run_gate.sh` のすべての段階が通った：V-ABI quick（15 の検索 × 経路とスレッド = 60 の実行がネイティブの CLI と一致、NCBI の凍結ハッシュ 16/16 一致）、`npm run check`（309 件、reactor の要る 29 件は skip）、reactor 付きの単体試験（338 件。HSP の対応の試験と、生成した検証の表 `verification-table.json`）、FakeEngine のビルドの E2E（Chromium 31・Firefox 25・WebKit 24、計 80 件）、エンジン入りのビルドの E2E（Chromium 44・Firefox 38・WebKit 37、計 119 件。V-BR はブラウザごとに 27 の検索の outfmt 0・6・7 と診断がネイティブの CLI と一致し、NCBI の凍結バイトとの比較 21（outfmt 0 が 18、outfmt 7 が 3）がすべて一致）、結果画面の E2E の 2 回の繰り返し（3 ブラウザ × 18、54 件）、計測（3 ブラウザ。下の「実測」）、画面の記録（3 ブラウザ × 2 サイズ × 18 枚 = 108 枚）。

ゲートの準備で 2 つの失敗があり、どちらも環境の問題で記録は残していない：`run_gate.sh` が mode 644 で実行できなかった（`bash` で呼ぶ）、ネイティブの CLI を `LOSAT-native` の名前で渡すと V-ABI quick の CLI の文言（clap の `Usage: LOSAT-native …`）が reactor の `LOSAT` と食い違った（同じバイナリを `LOSAT` の名前に写して渡した。W3 は `release/LOSAT` を渡していた）。

### レビューの後の実行

コードレビューの指摘（`d7ade9e9`、`d78f1448`）と画面レビューの 1 回目の指摘（`057b18f1`、`7ee77f8d`）を直した `7ee77f8d` の木で、`LOSAT_WEB_GATE_STEPS=after-review` の `run_gate.sh`（npm ci、check、単体、両方のビルドの E2E、画面の記録）が通った（[`run-20261009T135520Z-after-review/`](run-20261009T135520Z-after-review/)）：`npm run check`（314 件、reactor の要る 29 件は skip）、reactor 付きの単体試験（343 件。HSP の対応の試験は outfmt 0 の節と見出しの数も比べる）、FakeEngine のビルドの E2E（Chromium 32・Firefox 26・WebKit 25、計 83 件）、エンジン入りのビルドの E2E（Chromium 45・Firefox 39・WebKit 38、計 122 件。1280 px で表が収まること、390 px で横に流れないこと、「Open results」の後に結果の見出しが見えることを含む）、画面の記録（3 ブラウザ × 40 枚 = 120 枚。状態 17・18 を足した）。push の前の standard（FakeEngine の E2E 83 件）も通った。

### 計測だけの実行

コードレビューの L5・L6 で計測の時計を直した（判断 14）ので、`LOSAT_WEB_GATE_STEPS=measure`（npm ci と計測だけ）で測り直した。1 回目（[`run-20261009T140235Z-measure/`](run-20261009T140235Z-measure/)、`7ee77f8d`）は、1500 写しのドットプロットの場面を Firefox と WebKit で落とした。線が約 0.6 px 間隔で並ぶので、選ばれた線のすぐ横を狙ったクリックが選択を変えなかった（計測の的の選び方の問題で、画面の振る舞いは変えていない）。`50d57bba` で的を選ばれた線から最も遠い HSP にして、2 回目（[`run-20261009T141300Z-measure/`](run-20261009T141300Z-measure/)、`50d57bba`、check 314 件も通った）がすべての場面を記録した（見込みどおりの失敗を除く。判断 22）。下の「実測」はこの 2 回目の数値。ゲートの実行（`b2d8fb24`）の `measure/` は古い時計で、`selectFar`・`selectWide` はすでに選ばれた query のクリックを測っている。

## 作業との対応

| 作業 | 実装 | 試験 |
|---|---|---|
| 1. 列定義表 | [`docs/web/results_columns.md`](../../web/results_columns.md)：Subject 一覧・HSP 一覧・Query の一覧・表示用のフィルター・状態の区別の、名称・単位・program・値の出どころ・集約の単位・並べ替え・採否。NCBI Web の説明表の集約列（Max score・Total score・Query cover・Per. ident・表形式の E value）は採用し、エンジン側の S13+（E2j）で NCBI の `align_format` から移植するまで出さない（判断 1）。TS では計算しない | 単体 `results.test.ts`（表の値は outfmt 6 の欄の文字列のまま） |
| 2. 選択と一覧 | Run の選択（番号・program・入力名・`Options:`）、仮想化した Query の選択（10 万 query）、Subject 一覧と HSP 一覧（outfmt 6 の行を `out6` の範囲で読む。並べ替えはエンジンの原値、同値はエンジンの順）。選択の中心は HSP の ID（`runId`・`qIdx`・`rank`）で、並べ替えやフィルターを通して保つ | E2E「selection by HSP identity…」、単体 `results.test.ts` |
| 3. 詳細 | 選んだ HSP の outfmt 6 の行、subject の見出し（`out0_subject`）、outfmt 0 の節（`out0`）を原文のまま。`out0` が null の HSP は、outfmt 0 が query ごとに先頭の `num_alignments` 個の subject だけを示すことを理由として示す。選択が移った後に届いた読み取りは捨てる | E2E「an HSP that outfmt 0 does not show…」、単体（遅れた詳細、M2 で試験を直した） |
| 4. ドットプロット | Canvas に選んだ query と subject の組の HSP を描く。ズーム（ホイール、ボタン、+ / −）、パン（ドラッグ、矢印キー）、HSP の選択（線分のクリック、n / p、一覧）、「Zoom to HSP」「Whole sequences」。向きは色と凡例で区別し、軸に単位（nt / aa）と frame を示す。BLASTN の 1 文字の HSP は向きを「Not in the record」とし、outfmt 0 の節の `Strand=` を指す（判断 2。鎖の欄は E2j で ABI に足す） | E2E「the dot plot…」「BLASTN: an HSP of one letter…」、計測（下の「実測」） |
| 5. Run の詳細と検証バッジ | RunSnapshot と RunRecord（経路、スレッド、build、入力）、CLI コマンド。検証バッジの判定表は `build/verification.ts` がビルドのときに `docs/web/verification_cells.tsv` と認証の記録（outfmt 0 の manifest、fixture の manifest、BLASTP・TBLASTX の Gate A、TLOSAN の行列、E2e の option sweep）から生成し、手で書かない（判断 3） | 単体 `verification-table.test.ts`、`results.test.ts`、E2E（バッジの段階） |
| 6. 表示用のフィルター | E value ≤、Bit score ≥、Subject ID、Query ID、「Queries with hits only」。再検索とは別の操作（キューに Run が増えない）。表示 0 件・ヒット無し・失敗・未完了・取消・結果数の上限に達した可能性を区別して示し、上限に達しただけでは取れなかったヒットがあると言わない | E2E「many queries; view filters…」、単体 `results.test.ts` |
| 7. 試験 | 上の「完了条件と結果」 | `hsp-correspondence.test.ts`、`results.spec.ts`（9 件 × 3 ブラウザ × 2 ビルド） |

W3 からの引き継ぎ：キューから完了した Run の結果を開く「Open results」（W3 の画面レビュー 1 回目の L7）、キューの並び（実行中、待ち、終わったものを新しい順）と、待ちと取消の印の区別（2 回目の L-e）。

## 作ったもの

- `src/application/results.ts`（`ResultsBrowser`）：選択（run、query、subject、HSP）、表示用のフィルター（ViewState）、並べ替え、詳細と見出しの読み取り、検証バッジ。Run が完了・取消・失敗したときの表示の切り替え。
- `src/domain/`：`hsp-table.ts`（HSP レコードの列の形。typed array を Data worker から transfer）、`result-index.ts`（query と subject ごとのまとめ、並べ替え、フィルター、上限、向き）、`outfmt6.ts`（12 欄の分割）、`outfmt0.ts`（見出しの title）、`verification.ts`（argv の option の鍵、バッジの判定）。
- `RunStore.readHitTable`・`readOutputRange`、Data worker の対応、worker の RPC の typed array の transfer。FakeEngine の検索の形の出力（200 query まで、3 subject まで、3 つ目は outfmt 0 に無い、逆向きの HSP、BLASTN の 1 文字の HSP。値は `FAKE`）。
- `build/verification.ts`（判定表の生成。Vite の plugin `virtual:losat-verification`）。
- `src/ui/`：`ResultsPanel`（置き換え）、`QueryPicker`、`SubjectTable`、`HspTable`、`HspDetail`、`DotPlot`、`ResultFilters`、`ResultNotices`、`RunDetails`、`VerificationBadge`、`OutputsView`（S01 の原文と書き出し）、`VirtualRows`、`SortButton`、`QueuePanel` の変更。
- 試験：単体 `results.test.ts`、`hsp-correspondence.test.ts`、`verification-table.test.ts`。E2E `results.spec.ts`、計測 `results-measure.spec.ts`（`support/synthetic.ts`）、画面の記録 `screens.spec.ts` の状態 07〜16。既存の spec の補助は `support/search.ts` に移した。
- 文書：`docs/web/results_columns.md`、[`run_gate.sh`](run_gate.sh)（`LOSAT_WEB_GATE_STEPS=after-review`・`measure`）、S13+ の指示書の下書き [`session_s13p_e2j_subject_aggregates.md`](../../losat_web_gui_sessions/session_s13p_e2j_subject_aggregates.md)、S14 の指示書の引き継ぎ。

## 判断（推奨案で進めたもの）

Owner-delegated（2026-09-29 の常設の指示、2026-10-07 に再掲）。1〜5 は 2026-10-06 の中断したセッション、6〜8 と 9〜11 は 2026-10-08、12〜21 は 2026-10-09、22・23 は 2026-10-10 に取った。

1. **集約列**（Max score・Total score・Query cover・subject の Per. ident・表形式の E value）を採用し、SF の後にエンジン側の新しいセッション S13+（E2j）で NCBI の `align_format`（`showdefline.cpp`）からエンジンへ移植する。それまで S13 は outfmt 6 / outfmt 0 の値だけを示す。TS で計算しない。
2. **BLASTN の 1 文字の HSP の鎖**：E2j で ABI に足す（エンジンの HSP の query の frame から、HSP レコードの新しい欄）。それまでアプリは「not in the record」と示し、outfmt 0 の節の `Strand=` の行を指す。
3. **検証バッジ**：option の組の完全一致（`optionKey`）で、NCBI と比べた記録と照らす。ブラウザの経路とスレッド数は V-BR の `checked` の升目（1/2/4）に入ること。範囲は `<range>` に正規化。既定の BLASTN は、E2j が既定の outfmt 7 の fixture を足すまで「outside」。
4. **データの読み方**：`readHitTable`（列の形、typed array を transfer）と `readOutputRange`（節。見出しは見えている行の分だけ後から読む）。
5. **S13 では Run を消さない**：`runInputs` の解放は残件のまま（W3 のゲート記録）。
6. **HSP の対応の試験と凍結バイト**：承認済みの `-db_gencode` の例外の行は、NCBI のデータベース検索（`<fixture_id>.db.out`、`without_posted_date` の後のハッシュ）とデータベースの行の外で比べる（`docs/evidence/losat_web_e2a/check_losat.py` と同じ）。outfmt 6/7 の field list の fixture は比べない（reactor は既定の形式を書く）が、その検索の対応は確かめる。
7. **エンジン入りの E2E の「outfmt 0 に無い」**：検索画面に `-num_alignments` が無いので、many.blastn の fixture（`-task blastn`、260 subject、outfmt 0 に 250）を使う。
8. **失敗した Run の表示**：エンジン入りのビルドだけで確かめる（FakeEngine は Run を失敗させられない）。
9. **並べ替えは一覧の先頭を見せる**：選んだ行を見せ続けると、260 subject を「#」の降順にしたとき一覧の終わりに移った。並べ替え（`orderKey`）では先頭へ、選択の変更（`revealKey`）と行数の変更（フィルター）では選んだ行へ。選択は ID で保つ。
10. **選んだ query を隠す query のフィルターは、選択を最初に並ぶ query に移す**（subject のフィルターと同じ規則）。並ぶ query が無ければ選ばない。
11. **別の Run はフィルター無しで開く**（選択と並べ替えと同じ）。同じ Run は保つ。Run ごとのフィルターの記憶はしない（計画 §5.7 は求めていない）。
12. **Subject 一覧は示す値で並べる**（L2）：行はフィルターの前の最初の HSP と HSP 数を示し、並べ替えも同じ値を使う（前はフィルターを通った最初の HSP で並べた）。`results_columns.md` に補足を書いた。
13. **保存した結果を読めなかった完了 Run は一度だけ読む**（L3）：誤りを示し、キューが変わるたびに読み直さない。Run を開き直せば読む。
14. **計測の時計は名前どおり**（L5・L6）：判断 10 の後は query のフィルターがその query を選ぶので、`findFar` はフィルターから見つけた query の最初の HSP を読み終えるまで、`selectFar` は別の query（その前の query）のクリック、`selectWide` は 2 つ目の 200 subject の query のクリック、`firstHsp` は Run を開き直して 1 つの時計で測る、`indexMs` は query のファイルだけ。`selectByKey`・`selectByClick` は Playwright の往復を含む上限（spec に書いた）。ゲートの実行の `measure/`（古い時計）は残し、下の「実測」の数値は直した spec の計測だけの実行（`LOSAT_WEB_GATE_STEPS=measure`）から取る。
15. **フィルターの欄は効いている値を示す**（M1）：フォームを描き直すと（結果タブに戻る、同じ Run を開き直す）、E value と bit score の欄は効いている数を JavaScript の書き方（`1e-10`、`5000`）で示す。打った文字は、その数が効いている間はそのまま残す。

16. **表は 1280 px の枠に収める**（画面レビュー H1・M1）：列は値に合わせ、見出しは 2 行まで折り返してよく、数の見出しは右寄せ。収まらない値は省略記号で切り、全体を `title` に出す。横の scroll は枠がもっと狭いときだけ。E2E が 3 ブラウザで測る。
17. **電話では大事な値を先に**（M2）：800 px 未満で、Subject 一覧は Score・E value・HSPs を Length と Description の前に、HSP 一覧は Bit score・E value・Query・Subject・Frames・Orientation を他の列の前に並べ、表が横に流れるときだけそう言う。2 段の行は採らなかった（仮想化した一覧は一覧ごとに 1 つの行の高さ）。
18. **ドットプロットはページの scroll を奪わない**（L4）：ホイールは Ctrl か ⌘ を押したときだけ拡大し、電話では図の上でもページを縦に動かせる（`touch-action: pan-y`）。ボタン、+ / −、ドラッグ、矢印キーは残す。
19. **Run の日付は ISO 8601**（L5）：`YYYY-MM-DD HH:MM:SS`（その機械の時刻）。en-GB の `09/10/2026` は日と月のどちらが先か読めない。
20. **結果タブの配置と表の字**：結果タブを見ている間（801 px 以上）はページの列を 4fr / 260 px 以上にする（検索タブは 3fr / 260 px 以上のまま）。表は 0.875em、列は `ch`（E value は `1.02e-129` のため 8.25ch）、見出しは 0.9em の太字で折り返してよい。見出しは行の scroll bar の幅を空ける（VirtualRows が測る）。切れて示す値が 2 種類残る（7 桁の HSP のアラインメントの長さ、Orientation の「Not in the record」）。
21. **画面の記録**：状態 01 は option の検査が終わってから撮る。状態 17（翻訳する program のドットプロット）と 18（電話で「Open results」の直後の画面）を足した。WebKit の Run の select は label を select の内側で切る（WebKit は select の `overflow` を無視する）。
22. **1 回目の計測だけの実行は残し、測り直した**：`run-20261009T140235Z-measure`（`7ee77f8d`）は、線が 0.6 px 間隔の 1500 写しの場面を Firefox と WebKit で落とした（選ばれた線の横を狙ったクリック）。`50d57bba` で的を選ばれた線から最も遠い HSP にし、2 回目の数値を使う。実行記録は書き換えない。
23. **画面レビューの 2 回目の残りは S13b へ**：中 1（電話で Subject の ID が 4 文字に切れる。M2 の修正の退行）はここで直した（`b565eab8`、試験つき）。中 2（Firefox の電話で行のクリックの後に表が横に流れる。1 回目の記録も同じで、原因は決まっていない）と低 3〜8 は、表・キュー・ドットプロットを NCBI の形に作り直す S13b に回し、その指示書に書いた。

## 実測

計測は `tests/e2e/results-measure.spec.ts`（計測だけの実行 `run-20261009T141300Z-measure/measure/`、commit `50d57bba`）。1 回の暖機と 3 回の計測の中央値で、ページの中の `performance.now()` による、操作から結果を示す最初のフレームまでの ms（分解能は 1 フレーム、約 17 ms）。

**多数の query**（BLASTN、100 文字の query × N、203 subject。70% は 1 subject、20% は 3 subject、1000 個に 1 個は 200 subject、10% はヒット無し）：

| query | ブラウザ | 検索（s） | HSP | 開く | 最初の HSP まで | 終わりまで scroll | ID で探す（選択と詳細まで） | 絞り込みを外す | 別の query をクリック | With hits only | 200 subject の query をクリック | subject の並べ替え（最大） |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10,000 | Chromium | 0.8 | 14,981 | 39 | 43 | 1 | 11 | 12 | 8 | 12 | 9 | 12 |
| 10,000 | Firefox | 2.8 | 14,981 | 67 | 74 | 1 | 10 | 10 | 19 | 10 | 13 | 13 |
| 10,000 | WebKit | 1.4 | 14,981 | 61 | 85 | 10 | 7 | 15 | 15 | 7 | 47 | 11 |
| 100,000 | Chromium | 6.7 | 149,828 | 300 | 313 | 1 | 12 | 16 | 18 | 15 | 9 | 12 |
| 100,000 | Firefox | 22.4 | 149,828 | 590 | 660 | 1 | 10 | 17 | 21 | 16 | 8 | 11 |
| 100,000 | WebKit | 7.1 | 149,828 | 521 | 502 | 11 | 34 | 18 | 15 | 17 | 15 | 13 |

- 10 万 query（13.9 MB、149,828 HSP）の Run を開くのは 0.30〜0.59 秒（HSP レコードの読み取りと transfer、索引）、最初の HSP の詳細まで 0.31〜0.66 秒。その後の操作（scroll、ID で探す、絞り込み、クリック、並べ替え）はどれも 1〜2 フレーム（34 ms 以内）。Query の選択は 20 行、Subject 一覧は 22 行だけを描く（仮想化）。
- Chromium の JS heap：開く前 14.4 MB、最初に開いた後 49.6 MB、繰り返しの後 49.8 MB（増え続けない）。
- Outputs の表示で outfmt 6 の全文（9.0 MB）を出すのは、Chromium で 2.7 秒、Firefox と WebKit は 1 フレーム（ブラウザの `<pre>` の描き方の差）。

**1 組の多数の HSP**（20 文字の単位を 1% の違いで繰り返す配列の自己検索。写しごとに約 2 HSP）：

| 写し | ブラウザ | 検索（s） | HSP | 開く | ドットプロットを出す | 拡大 | 縮小 | Zoom to HSP | 全体 | n で選ぶ | クリックで選ぶ | 一覧で選ぶ | HSP の並べ替え |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,500 | Chromium | 3.8 | 2,995 | 185 | 42 | 28 | 38 | 40 | 38 | 92 | 99 | 16 | 11 |
| 1,500 | Firefox | 18.4 | 2,995 | 326 | 96 | 34 | 46 | 45 | 47 | 109 | 130 | 15 | 10 |
| 1,500 | WebKit | 4.3 | 2,995 | 139 | 11（描画 224） | 10（59） | 12（120） | 4（103） | 13（115） | 151 | 229 | 14 | 12 |
| 3,000 | Chromium | 15.7 | 5,993 | 754 | 86 | 55 | 79 | 72 | 71 | 169 | 163 | 15 | 9 |
| 3,000 | Firefox | 76.0 | 5,993 | 1,308 | 179 | 67 | 89 | 89 | 87 | 186 | 199 | 14 | 9 |
| 3,000 | WebKit | 失敗：結果が 512 MB の保存の上限を超えた（OPFS の無い WebKit は結果をメモリに置く） | | | | | | | | | | | |
| 4,000 | 3 ブラウザ | 失敗：エンジンのメモリ不足（`memory allocation of 49881 bytes failed`） | | | | | | | | | | | |

- 5,993 HSP のドットプロットは、出すのに 86〜179 ms、拡大・縮小・移動は 55〜89 ms（Chromium と Firefox）。HSP の一覧は 20 行だけを描き、選択・並べ替え・scroll は 1 フレーム。
- WebKit は最初のフレームが早く、canvas の描画が次のフレームに出る（括弧の中が描画のフレームまで）。
- 「n で選ぶ」「クリックで選ぶ」は、Playwright の往復を含む上限（判断 14）。
- 計画 §2.3 の「図と一覧の集約表示」の閾値は、保守者がこの数値から決める（S14 以降）。
- 同じ機械でエンジン側のセッションが並行していた時間があり（32 論理プロセッサ）、検索の時間の揺れはそのため（例：10 万 query の検索は Chromium で 4.8〜6.7 秒）。

## 独立レビュー（コード）

agent `losat-reviewer`（code review の役。このセッションは clone の外で起動したので、定義を general-purpose の agent に渡した）が `f3048ffd..b2d8fb24` の `web/app` を読んだ。結論は「妨げになる指摘は無い」。指摘と対応：

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 結果タブに戻る、または同じ Run を開き直すとフィルターのフォームが描き直され、E value と bit score の欄が空なのにフィルターは効いたまま。その後 Apply すると、利用者が外していないフィルターが外れる | `d7ade9e9`：欄に効いている値を示す（判断 15）。E2E（古いコードで落ちる） |
| M2 | 遅れた詳細を捨てる単体試験が、読み取りを順に答えるので、確かめる行を消しても通る | `d7ade9e9`：後の選択の読み取りを先に答える（行を消すと落ちることを確かめた） |
| L1 | 前の Run の見出しの読み手が、次に開いた Run の待ち行列を消す（Description が「…」のまま） | `d7ade9e9`：読み手を Run の load token に結ぶ。単体試験（古いコードで落ちる） |
| L2 | Subject の並べ替えがフィルターを通った最初の HSP を使い、行はフィルターの前の最初の HSP を示す | `d7ade9e9`：示す値で並べる（判断 12）。単体試験（古いコードで落ちる） |
| L3 | 保存した結果を読めなかった完了 Run を、キューが変わるたびに読み直す | `d7ade9e9`：一度だけ読む（判断 13）。単体試験（古いコードで落ちる） |
| L4 | 対応の試験が outfmt 0 の節と見出しを数えない。判断 6 の数は記録するだけ | `d7ade9e9`：節（`^ Score =`）と見出し（`^>`）の数を記録と比べ、byte 一致とデータベース検索の比較がそれぞれ 1 つ以上あることを求める |
| L5 | 判断 10 の後、`selectFar`・`selectWide` はすでに選ばれた query をクリックしている | `d78f1448`：判断 14 |
| L6 | `indexMs` が subject を含む、`selectByKey`・`selectByClick` が Playwright の往復を含む、`firstHsp` が 2 つの段階の間を含まない | `d78f1448`：判断 14 |

確かめて指摘の無かったもの：バイトの範囲と null の欄、同値のエンジンの順と 10 万 query の費用、向きと frame（BLASTN の 1 文字、翻訳する配列だけの frame）、ドットプロットの単位とズーム・パン・選択の写像、判断 9〜11、状態の区別（上限に達しただけで「取れなかった」と言わない）、FakeEngine の印、英語の UI、層の規則、検証の表（バイトの比較だけを取り、証拠より多くを認証しない。V-BR の升目の「serial は 1、threaded は 2 と 4」を経路とスレッドの集合に分けている点は、今は `runtime.ts` が serial を 1 スレッドにするので害が無い。将来の升目が組を変えるなら組のまま持つ）。

## 画面レビュー

agent `losat-reviewer`（screen review の役）が、ゲートの画面の記録（108 枚、`$BUILD_ROOT/s13-gate-screens/`）を W3 の記録（`/home/kawato/.cache/losat-web-gui-target/app-s12-gate-screens/`）と並べて見た。

**1 回目：不合格**（High 2 件）。検索画面（01〜06）に S13 による退行は無い（キューの「Open results」と並び、待ちと取消の印の区別は予定どおり。狭い画面のレコードの行が 2 段なのは W3 の `524ea6a3f` で、W3 のゲートの記録がその前のものだから。`bd7e4900` は検索画面に触れていない）。

| 指摘 | 内容 | 対応（`057b18f1`、`7ee77f8d`） |
|---|---|---|
| H1 | 1280 px で HSP 一覧が枠より広く（約 1,200 px に約 918 px）、Subject の値が半分の字で切れ、見出しが左寄せなので Query の範囲を Subject の値と読み違える。状態 14 は frame の値が見えない | 表が枠に収まる（判断 16、20）。E2E が 4 program で `.table-scroll` の横のはみ出しが無いことを測る（3 ブラウザ、両方のビルド） |
| H2 | WebKit の 390 px で、結果画面のすべての状態がページを 522〜675 px に広げる（Run の select が選んだ label の長さに伸びる） | select の label を切る（判断 21）。E2E が 390 px で横のスクロールが無いことを結果画面の状態で測る |
| M1 | 数の列の見出しが左寄せで、右寄せの値の上に無い | 見出しを右寄せ（判断 16） |
| M2 | 390 px の表は ID しか見えない | 大事な列を先に並べ、横に流れるときだけ「Scroll the table sideways for more columns →」（判断 17） |
| M3 | キューの「Open results」が結果を見える所に出さない（電話ではキューが約 2,000 px 下） | 結果の見出しへ移ってフォーカスする。状態 18 が直後の画面を撮る |
| L1〜L8 | 「With hits only」の checkbox と label の間、option の名前が `-` の後で折れる、「Open results」が常に 1 行を取る、ドットプロット（11 px の字、選択の印の凡例、ホイール、`touch-action`、翻訳する program の図が無い）、日付の書き方と列の位置、01 の撮る時、取消の打ち消し線、電話の並べ替えのボタンの高さ | すべて直した（判断 18、19、21）。L2・L3・L7 は E2E でなく目で確かめた |

**2 回目：合格**（High 無し）。レビューの後の実行の 120 枚（`7ee77f8d`、`$BUILD_ROOT/s13-review-screens/`）を 1 回目と並べた。H1（1280 px で 3 ブラウザとも全部の列が収まり、見出しがその値の上にある）、H2（WebKit の電話の 30 枚がすべて 780 px 幅）、M1、M3（状態 18 で「Results」の見出しが画面の上にある）、L1、L2、L5〜L8 は直った。L3・L4 は副作用つきで直った。検索画面（01〜06）は予定の変化（01 の option の行、キュー）だけ。残りの指摘と対応：

| 指摘 | 内容 | 対応 |
|---|---|---|
| 中 1 | 電話の Subject 一覧で ID が 4 文字ほどに切れる（「LC73…」。M2 の修正の退行。列の最小幅の合計が 390 px を超え、Subject が 6.5ch に縮む）。共通の接頭辞の ID が同じ行に見える | `b565eab8`：電話でも Subject を 12ch 以上にする。E2E が 390 px で測る（前の CSS では 49 px 足りずに落ちる） |
| 中 2 | Firefox の電話の記録で、行のクリックの後に表が横に流れている（1 回目の記録も同じ。Playwright か Firefox のフォーカスの scroll かは決まらない） | S13b へ（判断 23） |
| 低 3〜8 | 電話の HSP 一覧が subject の範囲の途中で切れる、キューのカードの折り返しが 3 通りでタブでキューが動く、ドットプロットの電話の説明、x 軸の端のラベル、TBLASTN・BLASTX のドットプロットの記録が無い | S13b へ（判断 23）。S13b は表、キューの置き場所、ドットプロットを NCBI の形に作り直す |

`b565eab8` の後の実行：[`run-20261009T142446Z-after-review/`](run-20261009T142446Z-after-review/)（`b565eab8`）で `LOSAT_WEB_GATE_STEPS=after-review` がもう一度通った：`npm run check`（314 件、29 件 skip）、reactor 付きの単体試験（343 件）、FakeEngine のビルドの E2E（32・26・25、計 83 件）、エンジン入りのビルドの E2E（45・39・38、計 122 件）、画面の記録 120 枚（`$BUILD_ROOT/s13-review2-screens/`）。push の前の standard（check と FakeEngine の E2E 83 件）も通った。電話の状態 14 で Subject が 12 文字（「LC738875_2…」）見えることを目で確かめた（3 回目の画面レビューは行っていない）。

## 見つけて直したもの

- `8f826809`：BLASTP の HSP 一覧に Frames（+1/+1）、ドットプロットに「frames +1 / +1」が出ていた。エンジンの BLASTP の HSP レコードが両方の frame に 1 を持つため。frame は翻訳する配列（TBLASTN の subject、TBLASTX の両方、SX からの BLASTX の query）にだけ示す（エンジン入りの E2E が見つけた）。
- `bd7e4900`：表の見出しの小さい字を grid（列は em）に付けていたので、見出しの列が行の列より狭かった（1280 px で「Orientation」が x=871、セルは x=1009）。字の大きさを見出しのセルに付けた（画面の記録で見つけた）。
- `b2d8fb24`：E2E の第 1 部が見た 3 つの振る舞い（判断 9〜11）。
- 復元のとき：W4 の変更が `coordinator.test.ts` の保存の上限の試験を壊していた（FakeEngine の outfmt 6 が FAKE の印を失った。印の行を足した、`f27132fc`）。
- コードレビューの指摘（上の表）。
- 画面レビューの修正の試験（エンジン入り）で、9 文字の E value（`1.02e-129`）が列に収まらず切れていた。列を 8.25ch にした（`057b18f1`）。

## 合流

`feature/losat-web-gui` への merge は、このセッションでは行わない。エンジン側の SF（E2h）が `.worktrees/web-gui` を使っているので、エンジン側がそのセッションの区切りで行う（`.worktrees/web-gui` には触れていない）。merge のときにエンジン側が行うこと：

1. `origin/feature/losat-web-gui-app` を `feature/losat-web-gui` に merge する（このブランチの起点 `f3048ffd` はそのブランチの上にある）。
2. README の表の S13 の行を「完了」にし、このゲート記録を指す。
3. 下の「エンジン側への申し送り」の行を足す。

## 保守者の指示（2026-10-09、このセッションの最中）

- 「次のセッションで NCBI BLAST ウェブサイトの GUI にデザインを寄せてほしい」。範囲は保守者が確かめた：NCBI BLAST の配置・節の順・言葉に寄せ、NCBI のロゴ・ヘッダー・配色などのブランドは写さない。値はエンジンが書いたまま示す規則は変えない。
- 「Dotplot は、もしできたら [`blast2dotplot.py`](https://github.com/satoshikawato/bio_small_scripts/blob/main/blast2dotplot.py) を参考にしてほしい」「もっと洗練されて軽快に動いてクリックしてポップアップするようなのができるならそれはそれでいい」。
- 対応：S14 の前にアプリ側の新しいセッション S13b（段階 W4b）を入れ、指示書 [`session_s13b_w4b_ncbi_style_ui.md`](../../losat_web_gui_sessions/session_s13b_w4b_ncbi_style_ui.md) を書いた。下調べは [`docs/web/ncbi_blast_gui_survey_20261009.md`](../../web/ncbi_blast_gui_survey_20261009.md)（その agent は公開の配列 NM_000518.5 と NM_000519.4 の blastn の 2 配列の検索を NCBI で 1 件実行した。利用者のデータは送っていない）。S14 の開始の条件に S13b の完了を足した。

## エンジン側への申し送り

- **S13b（W4b）の行**：README の表と計画 §7 に、S13 と S14 の間のアプリ側のセッション S13b「画面を NCBI BLAST Web に寄せる」を足す（上の「保守者の指示」）。完了条件は指示書のとおり（NCBI の参照画面との対応表、検索画面と結果画面がそれに従うこと（画面レビュー）、E2E と HSP の対応の試験、計測が S13 から大きく悪くならないこと）。保守者の指示を計画 §0.4 に DW として記す。

- **S13+（E2j）の行**：README の表と計画 §7 に、SF の後のエンジン側のセッション S13+「subject ごとの集約の値と HSP の鎖」（段階 E2j、指示書 [`session_s13p_e2j_subject_aggregates.md`](../../losat_web_gui_sessions/session_s13p_e2j_subject_aggregates.md)）を足す。判断 1（NCBI の説明表の Max score・Total score・Query cover・Per. ident・表形式の E value を `align_format` から移植し、ABI v2 で subject ごとの値を NCBI の文字列のまま渡す）と判断 2（HSP レコードに鎖の欄）による。既定の BLASTN の outfmt 7 の fixture を足せば、既定の BLASTN の検証バッジが「certified」になる（判断 3）。
- **DW の行**：計画 §0.4 に、S13 の判断 1〜3（集約列の採用と移植、鎖の ABI への追加、検証バッジの規則）を保守者の確認待ちとして記す。保守者が確かめたら DW にする。
- **エンジンのメモリ不足**：20 nt の単位を 4000 回繰り返す配列（81 KB）の自己検索（約 8,000 HSP）で、3 ブラウザともエンジンが「memory allocation of 49881 bytes failed」で止まった（3000 回、5,993 HSP は Chromium と Firefox で完了）。エンジン側の残件として、メモリの上限と HSP の数の関係を調べる。

## S14 への申し送り

[S14 の指示書](../../losat_web_gui_sessions/session_s14_w5_extraction_candidates.md)の「S13（W4）から引き継ぐこと」に書いた（結果画面の構成、データの読み方、座標と向き、E2E の補助、FakeEngine、検証バッジ、Run の入力、画面の記録の場所 `$BUILD_ROOT/s13-review2-screens/`、ゲートの script と計測の数値）。

## 残件と注意

- **エンジンのメモリ不足**（上の「エンジン側への申し送り」）。
- **WebKit の保存の上限**：OPFS の無い WebKit（この機械の Playwright の WebKit）は結果をメモリに置き、512 MB を超える Run を「Not enough temporary storage」として失敗させる（3000 回の繰り返しの自己検索）。W1 の設計どおりで、V-MOB の実機（S17）で確かめる。
- **`runInputs` の解放**：Run の削除が入るまで解放しない（判断 5、W3 の残件のまま）。
- **集約の列**：S13+（E2j）の後に、列定義表の「採用・エンジン待ち」の列を NCBI の表形式の並びで出す（S14 の指示書に書いた）。
- **切れて示す値**：表の値のうち 2 種類は、収まらないとき省略記号で切り、全体を `title` に出す（7 桁の HSP のアラインメントの長さ、Orientation の「Not in the record」。判断 20）。S13b で表を NCBI の形に寄せるときに見直す。
- 画面の記録（PNG）はリポジトリに入れない：ゲートの 108 枚は `$BUILD_ROOT/s13-gate-screens/`（`/home/kawato/.cache/losat-work/s13-gate-screens/`）、レビューの後の 120 枚は `$BUILD_ROOT/s13-review-screens/`（`7ee77f8d`）と `$BUILD_ROOT/s13-review2-screens/`（`b565eab8`。S14 と S13b はこれと比べる）。SHA-256 はそれぞれの実行記録の `screens.sha256`。
- 同じ機械でエンジン側のセッションが並行していた時間があり（32 論理プロセッサ）、計測の絶対値には揺れがある。
