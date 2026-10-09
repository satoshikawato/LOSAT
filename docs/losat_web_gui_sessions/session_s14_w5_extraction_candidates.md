# Session S14 — W5：抽出と候補

## INSTRUCTION PROMPT

LOSAT Web の段階 W5 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S14 の行である。設計は計画 §5.7、設計書 §11.3〜§11.4 にある。

1. 座標の変換を `web/app/src/domain/` の 1 か所に置く（1 始まりで両端を含む座標、strand、単位 nt / aa）。outfmt 6 の座標（start > end は minus 鎖）から、元のレコードの区間を得る。
2. 配列の回収：ヒット区間、左右を別々に指定する flank、Subject の全長、選んだ候補の multi-FASTA。既定は元のレコードの向き。端で切り詰めたときは、要求した範囲と実際の範囲を両方記録する。核酸の長さとタンパク質の残基数は単位を示す。原配列は Data worker が元の File から読み、大文字小文字を保つ。
3. 同じ Subject の複数の HSP は、別々の配列にするか、間を含む 1 区間にするかを選べる。重なりを勝手にまとめない。CDS の復元やスプライシングの推定はしない。
4. ギャップ付きアラインメントの書き出しは、HSP レコードの整列文字列から作り、原配列の抽出とは区別する。逆向きのヒットは、検索で得た向きのまま示す。
5. 候補トレイ：複数の実行から候補を集め、元の結果（run、`q_idx`、`rank`）へ戻れる。メモ、並べ替え、削除、由来の一覧。候補の登録・抽出・書き出しは、完了した実行だけで行える（`REQ-10`。アプリ層で拒否する）。
6. 試験：領域を指定した実行、端、逆向き、翻訳（BLASTX・TBLASTN・TBLASTX）の各場合で、抽出した配列が元のファイルの配列と一致すること。重複 ID のレコードを取り違えないこと。

完了条件は計画 §7 の S14 の行による。記録は `docs/evidence/losat_web_w5/README.md`。

## 開始の条件と最初の作業

- 場所と手順：Linux の clone の worktree `/home/kawato/losat-work/.worktrees/web-gui-app`（ブランチ `feature/losat-web-gui-app`）。場所は clone の `CLAUDE.local.md`、手順は skill `losat-campaign`・`losat-worktree`・`losat-gates`・`losat-ship-pr` に従う（README の規則 1・5 の `/mnt/c` のパスは `CLAUDE.local.md` の表で読み替える）。
- 開始の条件：S13b（W4b、画面を NCBI BLAST Web に寄せる。保守者の指示、2026-10-09）が完了し、エンジン側が W4 と W4b（このブランチ）を `feature/losat-web-gui` に merge し、SF（E2h、`CFastaReader` の移植）がそのブランチに入っていること（DW-23 (6)：抽出は SF の読み方に乗る）。満たさないときは始めず、保守者に伝える。
- 最初に `git merge origin/feature/losat-web-gui` を行い（衝突の解消は独立したコミット）、reactor とネイティブの CLI をその木から作り直す（`$BUILD_ROOT/s14-reactors`・`$BUILD_ROOT/native`）。FakeEngine の `describe.json` は `LOSAT_WEB_REACTORS=<dir> npx vitest run tests/unit/engine-runtime.test.ts -u` で書き直す。
- S13+（E2j、subject ごとの集約の値と HSP の鎖）が merge されていて、S13b がまだその指示書の 6.（アプリへの申し送り）を行っていなければ、最初に行う：列定義表の「採用・エンジン待ち」の列を NCBI の表形式の並びで出し、鎖の欄で BLASTN の 1 文字の HSP の向きを示す。

## S13（W4）から引き継ぐこと

[W4 のゲート記録](../evidence/losat_web_w4/README.md)の要点：

- **結果画面の構成**：`src/application/results.ts` の `ResultsBrowser` が、選択（run、query、subject、HSP）・表示用のフィルター（ViewState）・並べ替えを持つ。選択の中心は HSP の ID（`runId`・`qIdx`・`rank`、`HspId`）で、表の行番号ではない。候補トレイ（5.）の「元の結果へ戻る」は、この ID で `ResultsBrowser.open(runId)` と `selectQuery`・`selectSubject`・`selectHsp` を呼べばよい。 S13b が見出し・タブ・配置を NCBI BLAST Web に寄せた（Descriptions・Graphic Summary・Alignments・Dot Plot。下の「S13b（W4b）から引き継ぐこと」）。`ResultsBrowser` の規則と HSP の ID は保たれ、画面の test ID の一部は動いた。候補の操作の置き場所は S13b の対応表（`docs/web/ncbi_ui_mapping.md`）に合わせる。
- **データの読み方**：`RunStore.readHitTable(runId)` は HSP レコードを列の形（`src/domain/hsp-table.ts` の `HspTable`、typed array を Data worker から transfer）で返し、整列文字列（`query_aligned`・`subject_aligned`）を含まない。ギャップ付きアラインメントの書き出し（4.）は `readHits(runId)` のレコードの整列文字列から作る。outfmt 6 の行と outfmt 0 の節は `readOutputRange(runId, format, start, end)` で、HSP の `out6`・`out0`・`out0_subject` の範囲を読む。
- **座標と向き**：表示用の向きは `src/domain/result-index.ts` の `orientation`（`forward`・`reverse`・`unknown`）、frame は `hsp-table.ts` の `frame`（0 は無し）。単位は `src/domain/programs.ts` の `residueUnit`。1.（座標の変換を 1 か所に置く）は、これらと重ねずに domain にまとめ、結果画面もそれを使うように寄せる。BLASTN の 1 文字の HSP は start = end で向きが座標から分からない（S13 は「not in the record」と outfmt 0 の `Strand=` を示す。S13+ が鎖の欄を足すまで、抽出では向きを決めずに理由を示す）。
- **E2E の補助**：検索の入力を作って Run を完了させる補助は `tests/e2e/support/search.ts`（`paste`・`openFiles`・`settled`・`program`・`submit`・`waitStatus`・`result`）。結果画面の操作は `tests/e2e/results.spec.ts` を見る。キューの「Open results」（`run-N-open`）で完了した Run の結果を開ける。
- **FakeEngine**：検索の形をした出力（200 query まで、3 subject まで、3 つ目の subject は outfmt 0 に無い、逆向きの HSP、BLASTN の 1 文字の HSP）を書く。値は `FAKE`。FakeEngine の outfmt 6 は先頭に印の行がある（行は範囲で読む）。
- **検証バッジ**：`build/verification.ts` が `docs/web/verification_cells.tsv` と認証の記録から表を生成する（手で書かない）。認証の記録を足せば、表は次のビルドで変わる。
- **Run の入力**：`Coordinator.runInputs` は解放しない（Run の削除は S13 で入れなかった。W4 の判断 5）。候補トレイの削除は Run の削除ではない。
- **画面の記録と画面レビュー**：`tests/e2e/screens.spec.ts` が検索画面と結果画面の状態 01〜23 を撮る（W4 は 01〜18：17 は翻訳する program のドットプロット、18 は「Open results」の直後の電話の画面。S13b が 19〜23 を足し、07〜18 を新しい画面で撮り直した。状態の一覧は下の「S13b（W4b）から引き継ぐこと」）。W4 の最後の記録は `$BUILD_ROOT/s13-review2-screens/`（`/home/kawato/.cache/losat-work/s13-review2-screens/`）で、画面が NCBI の形に変わったので、W5 の画面レビューの基準は S13b の最後の記録 `$BUILD_ROOT/s13b-final-screens/`（`/home/kawato/.cache/losat-work/s13b-final-screens/`、SHA-256 は [`run-20261009T215513Z-after-review/screens.sha256`](../evidence/losat_web_w4b/run-20261009T215513Z-after-review/screens.sha256)）である。W5 の画面レビューはそれと比べ、検索画面と結果画面が変わらないことを確かめる。
- **ゲートの script と計測**：W4b の [`run_gate.sh`](../evidence/losat_web_w4b/run_gate.sh)（W4 のコピー。記録の場所だけ違う）を元にできる（`LOSAT_WEB_GATE_STEPS` は `all`・`after-review`・`measure`）。ネイティブの CLI は `LOSAT` の名前で渡す（`LOSAT-native` の名前では V-ABI quick の CLI の文言が食い違う）。W4 の計測（`results-measure.spec.ts`、計測だけの実行の数値は W4 のゲート記録の「実測」。S13b の数値は下の節）：10 万 query（149,828 HSP）の Run を開くのは 0.30〜0.59 秒、最初の HSP の詳細まで 0.31〜0.66 秒、その後の選択・絞り込み・並べ替えは 1〜2 フレーム。5,993 HSP の 1 組のドットプロットは出すのに 86〜179 ms、拡大と移動は 55〜89 ms（Chromium と Firefox）。20 nt の単位の 4000 写しの自己検索はエンジンのメモリ不足、WebKit は 512 MB を超える結果を保てない。抽出と候補トレイは、この規模で同じ操作の速さを保つ。

## S13b（W4b）から引き継ぐこと

[W4b のゲート記録](../evidence/losat_web_w4b/README.md)の要点。画面の対応表は [`docs/web/ncbi_ui_mapping.md`](../web/ncbi_ui_mapping.md)。

- **結果画面の新しい構成**（`src/ui/ResultsPanel.vue`）：上から、見出しの塊 `results-summary`（Job Title `results-job-title`、Run の select `results-run` と「Download All」`results-download-all`、Program `results-program` と検証バッジ、Options `results-run-options`、Query / Subject の ID と長さ `results-query-id` ほか。右に Filter Results）、複数 query の Run だけの「Results for」（`QueryPicker`。`query-list`・`query-filter`・`filter-hits-only`。どのタブでも見える）、notices（`results-notice`）、タブ。タブは Descriptions（`results-view-hits`、`SubjectTable`・`subject-list`）、Graphic Summary（`results-view-graphic`、`GraphicSummary`・`graphic-*`）、Alignments（`pane-alignment`、`AlignmentsView`）、Dot Plot（`pane-dotplot`、`DotPlot` の下に `HspTable`）、Run details（`results-view-details`）、Outputs（`results-view-outputs`）。選んだタブは Run を替えても、検索画面と結果画面を行き来しても保たれる（`AppView`）。どのタブも `ResultsBrowser` の同じ選択（HSP の ID）に従う。
- **HSP の表と選んだ HSP の詳細の場所**：W4 は Hits タブに Subject 一覧・HSP 一覧・「Alignment / Dot plot」の切り替え（`pane-alignment`・`pane-dotplot`）・`HspDetail` があった。今は HSP の表（`hsp-table`・`hsp-list`）が、Alignments では選んだ subject の塊の見出しの下に、Dot Plot では図の下にあり、Descriptions には無い。`pane-alignment`・`pane-dotplot` はトップレベルのタブのボタン（`aria-pressed`）になった。`HspDetail.vue` は無く、選んだ HSP の詳細（`hsp-detail`・`detail-row`・`detail-section`・`detail-not-in-outfmt0`・`detail-strand-note`）は Alignments の選んだ Range の塊（`range-<q>-<r>`）の中。Alignments の他の ID：`alignments`・`alignments-subject`（`data-subject`）・`alignments-summary`・`alignments-prev-subject` / `alignments-next-subject`・`alignments-descriptions`・`range-label`・`range-next` / `range-previous` / `range-first`・`range-section`（`data-state`）・`range-retry`・`alignments-show-earlier` / `alignments-show-later`。Range n は subject の座標（小さい方から）で、その subject の全 HSP を数える。候補の操作（行の選択、候補への追加）の置き場所は、対応表が「行の選択は S14 で検討」とした Descriptions の行と、Alignments の塊が候補になる。決めたら対応表に書き足す。
- **検索画面の新しい順と Job Title**：program のタブと一文、Enter Query Sequence（`query-panel`：FASTA の欄 `query-input`、Clear、Query subrange `query-region`、Or, upload file `query-open`、Genetic code、Job Title `job-title`、レコード一覧、Combined / Separate）、Enter Subject Sequence（`subject-panel`、同じ形、subject の遺伝暗号と例外の注）、Program Selection（`program-selection`、`param-task` の radio）、ボタン「Run LOSAT」（`add-to-queue`）と一文（`search-summary-line`）、開閉する Algorithm parameters（`algorithm-parameters`、既定と違う値の黄色と ♦、Restore default search parameters `restore-defaults`、2 つ目のボタン `add-to-queue-bottom`）。Job Title は `RunSnapshot.title`（任意、前後の空白を除く、argv に入れない）で、キューのカード（`run-<n>-title`）と結果の見出しに出る。候補トレイの「由来」の一覧は、Run の名前としてこれを使える。
- **ドットプロット**：`blast2dotplot.py` の約束（原点は左上、query が X・subject が Y、両軸同じ縮尺、`tick_size` の表と `#D3D3D3` の格子、`#1f77b4` と `#ff7f0e`、`int(pident)` で不透明度）。aa の軸と nt の軸の組は 1 aa を 3 nt と数えて縦横比を決める（目盛りは各軸の単位）。短い辺は最小 120 px で、足りなければ「Axes not to scale」。HSP を選ぶ（クリック、Enter、n / p、表）とポップアップ（`dotplot-popup`。outfmt 6 の値、frame、向き、outfmt 0 に有るか、「Show alignment」）。座標と向きの変換は 1. で domain にまとめるとき、`plot-geometry.ts` と重ねない。SVG の書き出しは S15。
- **対応表は新しい NCBI 風の要素を書き留める場所**：検索画面や結果画面に要素を足す（候補トレイ、抽出の操作）ときは、先に `docs/web/ncbi_ui_mapping.md` に NCBI の対応する場所・採否・理由を書く（LOSAT だけのものは「LOSAT だけ」とし、置き場所を 3. の表に）。
- **画面の記録**：`$BUILD_ROOT/s13b-final-screens/`（`/home/kawato/.cache/losat-work/s13b-final-screens/`、`492fad48`、3 ブラウザ × desktop 1280 px / phone 390 px で 150 枚、状態 01〜23。01〜06 が検索画面、07〜23 が結果画面で、12 が Run details、13 が Outputs、14 と 17 が TBLASTX（17 は翻訳する program のドットプロット）、15 と 23 が TBLASTN、18 が電話の「Open results」の直後、19 が Descriptions、20 が Graphic Summary、21 が 2 つの Range、22 がドットプロットのポップアップ）。SHA-256 は [`run-20261009T215513Z-after-review/screens.sha256`](../evidence/losat_web_w4b/run-20261009T215513Z-after-review/screens.sha256)。S14 の画面レビューはこれと比べる。S13b の画面レビューが合格として残した低い指摘が 2 つ：電話の入力のカードで長いファイル名が 2 行に伸びる（低 1）、電話の Chromium と Firefox の radio が 14 px（参考 1。`.radio input` に `flex: none`）。直すなら電話の状態 02・03・06 を 3 ブラウザで撮り直す。
- **ゲートの script と数値**：[`run_gate.sh`](../evidence/losat_web_w4b/run_gate.sh)（`LOSAT_WEB_GATE_STEPS` は `all`・`after-review`・`measure`）。最後の木 `492fad48` の件数：`npm run check` 378 件（29 件 skip）、reactor 付きの単体 407 件、FakeEngine のビルドの E2E 92 件、エンジン入りのビルドの E2E 131 件、画面の記録 150 枚。計測は `087eaa13`（[`run-20261009T211154Z-measure/`](../evidence/losat_web_w4b/run-20261009T211154Z-measure/)）：10 万 query（149,828 HSP）の Run を開くのは Chromium 292・Firefox 500・WebKit 415 ms、最初の HSP の詳細まで 266・524・438 ms、その後の選択・絞り込み・並べ替えは 1〜2 フレーム。5,993 HSP の 1 組のドットプロットは出すのに Chromium 20・Firefox 38 ms（描画を含む）、拡大・縮小・移動は 4〜34 ms。W4 より遅いのは Firefox の 3000 写しの「一覧で選ぶ」（26 ms、W4 は 14 ms）などの 1 フレーム分。時計の意味は W4 の判断 14 に加えて、ドットプロットと Graphic Summary の「出す」「拡大」などが描画のフレームまで測る（W4b のコードレビュー M1）。抽出と候補トレイは、この規模で同じ操作の速さを保つ。
- **残件**（W4b のゲート記録の「残件と注意」）：エンジンのメモリ不足（4000 写し）、WebKit の 512 MB の上限、`runInputs` の解放、iOS Safari（V-MOB、S17）。BLASTX は SX まで画面で動かせないので、3 nt per aa の規則は単体試験だけで確かめてある。E2j が先に merge されていれば、この指示書の最初の作業の前に、列定義表の「採用・エンジン待ち」の列と BLASTN の鎖を対応表の置き場所（Descriptions の Max Score は Score (bits) の位置、Total Score・Query Cover は Score と E value の間）に出す。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S15 — 出力と再現性](session_s15_w6_export_session.md)。
