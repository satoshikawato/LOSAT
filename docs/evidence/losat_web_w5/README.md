# LOSAT Web W5（Session S14）ゲート記録

- 段階：W5 抽出と候補（[総合計画書](../../losat_web_gui_plan.md) §7 の S14、[指示書](../../losat_web_gui_sessions/session_s14_w5_extraction_candidates.md)、[対応表](../../web/ncbi_ui_mapping.md)の「候補（S14）」）。完了条件は計画 §7 の S14 の行：領域指定、端、逆向き、翻訳の各場合で原配列と一致
- ブランチ：`feature/losat-web-gui-app`（アプリ側。worktree `$WORK_ROOT/.worktrees/web-gui-app`、Linux の clone）。起点は `141615ec`。これは `origin/feature/losat-web-gui` の先頭（エンジン側が W4b を `dcbb149d` で merge し、SF から SFd までが入っている）に fast-forward したもので、衝突の解消のコミットは無い。reactor は `$BUILD_ROOT/s14-reactors`、ネイティブの CLI は `$BUILD_ROOT/s14-gate-native/LOSAT`（sha256 `0236291b…`）で、どちらも `141615ec` の木から作った。FakeEngine の `describe.json` は `141615ec` の reactor で書き直しても変わらなかった（44 件通過、判断 2）。エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表、依存ライブラリ（`package.json`）は変えていない。`docs/web/` は `ncbi_ui_mapping.md`（`7b1d3b3b`）だけを変えた
- 実行記録：[`run-20261010T090350Z/`](run-20261010T090350Z/)（ゲート、commit `b8641668` の木、全段階と計測、`9a96c6f7`）、[`run-20261010T100446Z-after-review/`](run-20261010T100446Z-after-review/)（`608f6f75`、最後のアプリの木、`563bff7d`）。どちらも作成後は書き換えない。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)、再現は [`run_gate.sh`](run_gate.sh)。画面の記録（PNG）はリポジトリに入れない（下の「残件と注意」）
- 判定：**完了条件を満たした**。ゲートの実行（`b8641668`）のすべての段階、コードレビュー 1 回（妨げになる指摘無し。M1・L2 は直し、L1・L3 は既知の限界として残す）、画面レビュー 2 回（1 回目は不合格：M1・M2・L1〜L4 を `608f6f75` で直した。2 回目は合格：低 1 と参考 4 件を S15 へ）、レビューの後の実行（`608f6f75`）が通った。下の「完了条件と結果」。計測は `b8641668` の木で行い、最後の木（`608f6f75`）では測り直していない。コードレビュー M1 の修正（`a39e19e7`）は計測の木に入っていない（下の「実測」）
- 保守者の判断待ち：判断 1〜24（下の「判断」）は推奨案で進めた（Owner-delegated、2026-09-29 の常設の指示、2026-10-07 に再掲）。すべて 2026-10-10 に取った。保守者が予約した判断は含まない。目に付くのは 1（索引の読み方を kind 0 から NCBI の読み方の kind 1・2 に替える）、4（NCBI の読み方が拒む入力は索引の時点で失敗し、除外の按は出ない）、6（抽出の header は LOSAT Web の形式）、10（コードレビュー L1・L3 は直さない）、12（プログラムを替えて kind が変わると除外が消える。保守者への問い）
- `feature/losat-web-gui` への merge：このセッションでは行わない。エンジン側が行う（下の「合流」）
- このセッションは、1 回目（S14 (1)）を保守者が中止し、Opus で最初からやり直したもの（S14 (2)）。1 回目の下書き（座標の module、`sequence-layout.ts`、型の差分。未検証でコミットしていない）は出発点の材料にした。作業は取りまとめ役と、作業単位 A・B1・B2・D・F・G の agent で行った。B1（`s14-extract`）と F は一時的な worktree で行い、`2c09616f` と `fc3a04d7` で戻した。コードレビューの修正 1 回目も別の作業用 branch（`s14-fix`）で行って載せた

## コミット

`141615ec..563bff7d` の 23 コミットと、この記録のコミット 2 つ（古い順）。

| コミット | 内容 |
|---|---|
| `2e677cea` | アプリ：HSP の区間・鎖・単位・frame の座標を 1 つの module に（`domain/coordinates.ts`、指示書 1） |
| `85f905a5` | アプリ：Data worker が原配列と HSP のレコードを索引で読む（`readResidues`、`readHspRecords`、指示書 2） |
| `0d86b9a7` | アプリ：抽出の計画と FASTA の形式（指示書 2〜4） |
| `284425c0` | アプリ：実検索（serial の reactor）に対する抽出の試験。kind 1 は ASCII だけを保つ（指示書 6） |
| `8fafb351` | アプリ：索引の読み方をエンジンの reader（scan の kind 1・2）に替え、行を名指しする拒否の記録を探す（判断 1・3〜5） |
| `2c09616f` | merge：抽出の核（座標、原配列の読み、抽出の計画）をアプリのブランチへ（B1、一時的な worktree から） |
| `969c3cb9` | アプリ：選んだ HSP の詳細の「Try again」、読めた節は塊の失敗に代わる、デスクトップのポップアップは最初のクリックで見える所へ（W4b の残件、判断 8） |
| `f3d0b891` | アプリ：電話で入力のカードの「Remove」を 1 行目に、radio は 20 px（W4b の残件） |
| `541823f0` | アプリ（試験だけ）：デスクトップのポップアップの試験は、エンジン入りのビルドで 1 px 未満の端を許す |
| `62dad0fd` | アプリ：候補トレイとその抽出をアプリ層に（指示書 2〜5） |
| `88236d65` | アプリ（試験だけ）：defline の無い最初のレコードの抽出と、実検索に対する候補トレイ |
| `fc3a04d7` | merge：W4b のレビューの残件の修正（F）をアプリのブランチへ |
| `7b1d3b3b` | 文書：対応表に候補トレイと、それに加える操作の置き場所（S14） |
| `cb875707` | アプリ：「Candidates」のタブ、行の印、結果の「Add to candidates」、トレイの 2 つのダウンロード（D） |
| `19f4b117` | アプリ（試験だけ）：候補トレイの画面の記録 24〜28 |
| `81b4b6ed` | 文書：W5 のゲートの script（W4b のコピー。記録の場所だけ W5） |
| `b8641668` | アプリ（試験だけ）：結果の計測に、候補トレイの操作を足す（G）。**ゲートの木** |
| `a39e19e7` | アプリ：抽出は区間の byte だけを読む（チャンクを読まない。コードレビュー M1） |
| `e928b64c` | アプリ：貼り付けの defline の按は、最初の行の拒否のときだけ（コードレビュー L2） |
| `010ba28c` | アプリ（試験だけ）：エンジンの抽出の試験は、flank を付けた区間を自分で計算する |
| `9a96c6f7` | 文書：`b8641668` のゲートの実行記録（すべての段階が通った） |
| `608f6f75` | アプリ：画面レビューの修正（端で切った抽出の単位、ドットプロットのポップアップの幅、トレイの列、状態 15 と 02 の記録）。**最後のアプリの木** |
| `563bff7d` | 文書：`608f6f75` の after-review の実行記録（すべての段階が通った） |
| このゲート記録を含むコミット | 文書：W5 のゲート記録と `evidence.sha256` |
| S15 への引き継ぎのコミット | 文書：S15 の指示書に「S14（W5）から引き継ぐこと」 |

エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表、セッションの README は、アプリ側のブランチでは変えていない（README 規則 1）。`docs/web/` は `ncbi_ui_mapping.md` の追記（S14 の行と置き場所）だけ。

## 完了条件と結果

指示書 6.「領域指定、端、逆向き、翻訳の各場合で、抽出した配列が元のファイルの配列と一致すること。重複 ID のレコードを取り違えないこと」の各場合。実検索の試験は reactor を載せた単体試験で、生成した入力の文字（試験が自分で持つ元の配列）と、エンジンが書いた整列文字列（`query_aligned`・`subject_aligned`）の両方に突き合わせる。

| 場合 | 結果 | 証拠 |
|---|---|---|
| 領域指定（`-query_loc` / `-subject_loc`） | 通過。座標は元のレコード全体の座標のまま、抽出は領域の外にも届く | `extraction-engine.test.ts`「BLASTN with -subject_loc and with -query_loc: coordinates stay those of the whole record」 |
| 端（レコードの両端で flank を切る） | 通過。要求した範囲と実際の範囲の両方を記録する。8 検索・45 HSP・240 の断片のうち、位置 1 で切れたもの 19、末尾で切れたもの 31。トレイ経由では 49 配列のうち 15 が端で切れた | `extraction-engine.test.ts`（6 件）。flank の計算そのものは文字どおりの期待値（`coordinates.test.ts`、`extraction.test.ts`）と、エンジンの試験が自分で計算する区間（`010ba28c`）で確かめる。FakeEngine と実エンジンの E2E `candidates.spec.ts`（flank が端で切れる。状態 28） |
| 逆向き（BLASTN の minus 鎖、TBLASTN の minus の frame） | 通過。既定は元のレコードの向きで、ヒットの鎖は header の `hit_strand=` に書く。整列文字列は検索で得た向きのまま。BLASTN の minus は逆相補、TBLASTN は frame で翻訳して、エンジンの整列行と一致する。8287 文字の整列行を比べた | `extraction-engine.test.ts`（BLASTN の plus と minus、TBLASTN の minus の frame）、`candidates.spec.ts` のエンジン入りの 2 件（BLASTN の minus、TBLASTN の subject を nt、query を aa） |
| 翻訳（TBLASTN、TBLASTX、BLASTX） | TBLASTN と TBLASTX は通過（実検索）。BLASTX は SX まで検索を動かせないので、合成のレコードで座標・単位・frame を単体試験だけで確かめた（実検索の試験は無い） | `extraction-engine.test.ts`（BLASTP・TBLASTN・TBLASTX）、`coordinates.test.ts`・`extraction.test.ts`の BLASTX（合成） |
| 重複 ID のレコード | 通過。位置（`{role}_record=K`）で区別する。除外したレコードの後は位置が詰まる | `extraction-engine.test.ts`「two subject records with the same ID: the hit on the second, then the first left out」、`data-extraction.test.ts`（除外と重複 ID）、`candidates.spec.ts`（エンジン入りの BLASTN） |
| defline の無い最初のレコード | 通過（offset 0/0。式からも checkpoint からも読む） | `extraction-engine.test.ts`「a subject file that starts with residues」、`data-extraction.test.ts` |
| 大きなレコード（checkpoint）、小文字、`U`、コメント行、CR LF、`eol` 3 | 通過 | `extraction-engine.test.ts`（BLASTN の 1 件目）、`sequence-layout.test.ts`（16 件）、`data-extraction.test.ts`（9 件。`a39e19e7` の後に 1 件増えた） |
| トレイを通した抽出（完了した Run の HSP を集め、各 region と join で抽出して読み戻す） | 通過。7 候補（2 つの Run）、49 配列、整列 7 つ | `extraction-engine.test.ts`（アプリ層の 1 件）、`candidates.test.ts`（23 件） |
| W4b の残件（コードレビュー 2 回目の低 2 件、画面レビュー 3 回目の低 1 と参考 1） | 通過。詳細の「Try again」（`detail-retry`）、読めた節は塊の失敗に代わる、ポップアップを置いてから scroll、電話の入力のカードは「Remove」が 1 行目、radio は 20 px。画面レビュー 1 回目の「W4b の残件」の確認で、電話 02・03・06 の入力のカードは 3 ブラウザとも直り、radio は Chromium と Firefox の電話で 20 px になった（WebKit は元から変わらない） | `results.test.ts`、`results.spec.ts`（`detail-retry`、デスクトップのポップアップ 1 件）。変えた試験が古い DotPlot / 節の監視で落ちることを確かめた（F の報告）。画面は `reviews/screen-review.md` の「Check 2」 |
| 速さ（指示書「抽出と候補トレイは、この規模で同じ操作の速さを保つ」） | 通過。トレイの操作は 200〜5,993 候補で、Chromium と Firefox では 1 フレーム（15 ms 以下）。抽出 100 件は 27〜137 ms（WebKit の 2,995 候補だけ 289 ms）。下の「実測」。アプリ層の試験も 20,000 候補の追加・並べ替え・移動・メモ・選択・削除が候補の数に比例することを確かめる | `run-20261010T090350Z/measure/`、`candidates.test.ts`（20,000 候補） |
| 画面（W4b の最後の記録との比較、画面レビュー） | 通過。2 回目の画面レビューが合格（High・Medium 無し）。状態 01〜23 は、新しい要素（「Candidates」のタブ、行の印、「Add to candidates」）と、実行ごとに変わる値の他は変わらない | `reviews/screen-review.md`、`screen-review-2.md`（下の「画面レビュー」）。記録は `$BUILD_ROOT/s14-gate-screens/`、`s14-review-screens/` |
| 既存と新しい E2E が 3 ブラウザ・両方のビルドで通る | 通過。FakeEngine のビルド 110 件（Chromium 41・Firefox 35・WebKit 34）、エンジン入りのビルド 140 件（51・45・44）。W4b の最後の木の 92 件、131 件から、`candidates.spec.ts`（FakeEngine 5 件 × 3、エンジン入り 2 件 × 3）とデスクトップのポップアップの試験（× 3）が増えた。結果画面の E2E の 2 回の繰り返しも通った（72 件）。エンジン入りには V-BR（ブラウザごとに 27 の検索の outfmt 0・6・7 と診断がネイティブの CLI と一致）を含む。最後の木（`608f6f75`）でも同じ件数 | `run-20261010T090350Z/`（`npm-e2e-fake-engine.log`、`npm-e2e-engine.log`、`e2e-results-repeat.log`、`records/v-br-*.json`）、`run-20261010T100446Z-after-review/` |
| HSP の対応の試験が通る（変えずに） | 通過。`hsp-correspondence.test.ts` は `141615ec` から変わっていない（`git diff 141615ec..HEAD` が空）。BLASTN・BLASTP・TBLASTN・TBLASTX の全 fixture（165 の検索、30,467 HSP、outfmt 0 の節 30,027、outfmt 0 に現れない HSP 440）の結果は W4b と同じで、NCBI の凍結バイトと 196 が byte 一致、承認済みの `-db_gencode` の例外 2 がデータベースの検索として一致、field list の 3 は比べない。単体試験は reactor 付きで 548 件（最後の木で 550 件） | `run-20261010T090350Z/unit-cases.log`、`run-20261010T100446Z-after-review/unit-cases.log` |

### ゲートの実行（`b8641668`）

2026-10-10 18:03〜18:35（JST）、約 32 分。`run_gate.sh` のすべての段階が通った：V-ABI quick（15 の検索 × 経路とスレッド = 60 の実行がネイティブの CLI と一致、NCBI の凍結ハッシュ 16/16 一致）、`npm run check`（27 ファイルのうち 25 が通り 2 が skip、488 件通過・60 件 skip。reactor の要る分は skip）、reactor 付きの単体試験（27 ファイル、548 件通過、skip 無し。HSP の対応の試験と、生成した検証の表 `verification-table.json`）、FakeEngine のビルドの E2E（110 件）、エンジン入りのビルドの E2E（140 件）、結果画面の E2E の 2 回の繰り返し（72 件）、計測（3 ブラウザ、各 3 件通過）、画面の記録（3 ブラウザ × 2 サイズ、180 枚、状態 01〜28、12 件通過）。落ちた試験、flaky、やり直しは無い。計測の見込みどおりの場面の失敗は W4b と同じ：20 nt の単位を 4000 写しの自己検索は 3 ブラウザでエンジンのメモリ不足（`memory allocation of 49861 bytes failed`）、3000 写しは WebKit で保存の上限（512 MB）。この場面ではトレイを測っていない。

### レビューの後の実行（`608f6f75`）

`LOSAT_WEB_GATE_STEPS=after-review` が通った（[`run-20261010T100446Z-after-review/`](run-20261010T100446Z-after-review/)）。`npm run check` 490 件（60 件 skip）、reactor 付きの単体試験 550 件、FakeEngine のビルドの E2E 110 件、エンジン入りのビルドの E2E 140 件、画面の記録 180 枚（`$BUILD_ROOT/s14-review-screens/`）。V-ABI と計測はこの mode では行わない。

どの実行の E2E の log にも `FAIL [opfs] storage full …` と `FAIL [run output] storage full …` が出る（FakeEngine に 3 行、エンジン入りに 6 行）。ブラウザの中の契約の確認が自分の記録として印字するもので、それを印字する Playwright の試験は通っている（落ちた試験は無い）。W4b と同じ。

## 作業との対応

| 作業 | 実装 | 試験 |
|---|---|---|
| 読み方の切り替え（指示書に無い前提。判断 1・3〜5、DW-23 (6)） | 索引の scan を kind 0 から、エンジンの reader の kind 1（核酸）・2（タンパク質）に替えた。kind は program と役割で決まる（`indexParser(program, role)`、BLASTX は query が 1・subject が 2）。`FastaParserKind = 1 \| 2`。最初のレコードの defline が無い場合（offset 0/0）、正の一様な `eol`、レコードの無いコメントだけの入力を受ける。プログラムを替えて kind が変わる役割の source は索引し直す。FakeScanner と FakeEngine も kind 1・2。`message-line.ts` は、拒否が名指す行をレコードに戻す。アダプタには kind 0 が残る | `draft.test.ts`、`coordinator.test.ts`、`tests/contract/record-scanner.contract.ts`、`data-service.test.ts`、`engine-runtime.test.ts`（拒否の行 → レコードの対応を実エンジンで）、`dataset.test.ts`、`search.spec.ts` |
| 1. 座標の変換を 1 か所に | `src/domain/coordinates.ts`：区間（1 始まり・両端を含む）、鎖、単位（nt / aa）、frame、flank、端の切り詰め（要求した範囲と実際の範囲）、複数の区間のまたぎ。`orientation`、`subjectSpan`、`residueUnit`（`Unit` を返す）、Graphic Summary の query の区間は、これを使うように寄せた。結果画面の見た目は変わらない（コードレビューが同値を確認） | `coordinates.test.ts`（10 件）、`results.test.ts`、既存の `plot-geometry.test.ts` |
| 2. 配列の回収（ヒット区間、左右別の flank、全長、multi-FASTA） | `src/domain/extraction.ts`（計画 `planExtraction`、header、60 文字の行）、`src/domain/sequence-layout.ts`（一様な `eol` と checkpoint から開始位置を求め、kind の規則で前向きに読む `ForwardReader`）、Data worker の `readResidues(revisionIds, position, intervals)`（元の File から区間の byte だけを読む。大文字小文字、`U` はそのまま。レコードの ID・長さ・残基の数を表と照らす）。端で切ったときの要求した範囲と実際の範囲、核酸の長さと残基数の単位は header と画面の要約に出る | `extraction.test.ts`（16）、`sequence-layout.test.ts`（16）、`data-extraction.test.ts`（9）、`extraction-engine.test.ts`（6）、`candidates.spec.ts` |
| 3. 同じ Subject の複数の HSP | 「Several HSPs on one record」：別々の配列（Separate sequences、重なりは残す）か、間を含む 1 区間（One region spanning them）。Run・役割・レコードの位置が同じものの間だけまたぐ。重なりを勝手にまとめない。CDS の復元もスプライシングの推定もしない。鎖の混ざった区間は「mixed」、BLASTN の 1 文字の HSP の鎖は決めず理由を示す | `extraction.test.ts`（separate、spanning、whole、別の Run / 役割をまたがない、鎖を決めない、HSP の食い違いを拒む） |
| 4. ギャップ付きアラインメントの書き出し | `alignmentHeaders`・`alignmentFasta`：HSP レコードの整列文字列（`query_aligned`・`subject_aligned`）を検索で得た向きのまま（start > end のまま）。原配列の抽出とは別の file（`losat-candidates-aligned.fa`）。整列行の無い HSP は「missing」に載せ、全部が無ければ失敗して何も保存しない。Data worker の `readHspRecords(runId, indices)` | `extraction.test.ts`（整列行、frame、行が無い HSP）、`candidates.test.ts`、`extraction-engine.test.ts`（整列 7 つ）、`candidates.spec.ts`（FakeEngine は byte 単位） |
| 5. 候補トレイ | `src/application/candidates.ts` の `CandidateTray`：複数の Run から集める、元の結果（run、`q_idx`、`rank`）へ戻る（`results.reveal`）、メモ、並べ替え、移動、削除、由来の一覧（Origins）。REQ-10：登録・抽出・書き出しは完了した Run だけ（アプリ層で拒否）。`ResultsBrowser` に `candidateSources`・`hspIdsOfSubject`・`markSubjects`・`reveal`。UI は主のタブ「Candidates」（`CandidatesPanel.vue`）、Descriptions の行の印と「Add to candidates」、Alignments の「Add all matches to candidates」と Range の「Add to candidates」、ドットプロットのポップアップの「Add to candidates」。トレイは app の状態だけで、保存も送信もしない | `candidates.test.ts`（状態ごとの REQ-10 の拒否 7 件、20,000 候補ほか 23 件）、`candidates.spec.ts`（FakeEngine 5 件、エンジン入り 2 件）、状態 24〜28 |
| 6. 試験 | 上の「完了条件と結果」 | 同上 |
| W4b の残件 | 判断 8。WP-F | 上の表 |
| 画面の対応表 | 候補の置き場所を先に書いた（`7b1d3b3b`）。NCBI の「Download」は参照画面に撮っていないので、言葉だけ寄せた | 画面レビュー |

## 作ったもの

- `src/domain/`：`coordinates.ts`、`sequence-layout.ts`、`extraction.ts`。`dataset.ts`（`FastaParserKind = 1 | 2`、`isFirstLineRefusal`）、`programs.ts`（`indexParser`）、`result-index.ts`・`hsp-table.ts` は座標の module を使うように寄せた。
- `src/ports/data.ts`：`readResidues`・`readHspRecords`。`src/infra/data/`：`data-service.ts`（読みと区間に限った byte の範囲）、`message-line.ts`（拒否の行 → レコード）。`infra/data-worker/rpc.ts` は配列の中の typed array も transfer する。
- `src/application/`：`candidates.ts`（`CandidateTray`）、`results.ts`（印、`reveal`、`candidateSources`）、`draft.ts`・`coordinator.ts`（kind を役割で選ぶ）。
- `src/ui/`：`CandidatesPanel.vue`、`SubjectTable.vue`（印の列）、`AlignmentsView.vue`、`DotPlot.vue`、`ResultsPanel.vue`（追加の確認の toast）、`SourceCard.vue`、`styles.css`、`plots.css`。
- `src/infra/fake/`：FakeEngine は整列行（`FAKE`）を、各組の 1 つ目の HSP に書く。
- 文書：`docs/web/ncbi_ui_mapping.md` の追記、[`run_gate.sh`](run_gate.sh)、この記録、S15 の指示書への引き継ぎ。

## 判断（推奨案で進めたもの）

Owner-delegated（2026-09-29 の常設の指示、2026-10-07 に再掲）。すべて 2026-10-10 に取った。1〜12 は取りまとめ役の `DECISIONS.md` の行、13 以降は各担当の報告にある判断のうち、その行に入っていないもの。出典は `DECISIONS.md` と各担当の報告（`logs/wp-*/report.md`、`logs/fix1/report.md`、`logs/fix2/report.md`。タスクフォルダ）。

### 取りまとめ役（`DECISIONS.md`）

1. 索引の scan は kind 0 から、エンジンの reader の kind 1（核酸）・kind 2（タンパク質）に替える。program と役割で選ぶ（`indexParser(program, role)`）。表の読み方が SF 以降の検索と同じになり（DW-23 (6)）、defline の無い最初のレコード、正の一様な `eol`、レコードの無いコメントだけの入力を受ける。
2. Step 0：`describe.json` を `$BUILD_ROOT/s14-reactors`（`141615ec`）で書き直して変化は無かった（44 件通過、ファイルは同一）。
3. アプリは scan の kind 0 を使わない。FakeScanner は kind 1・2（常に checkpoint の layout、失敗しない）。アダプタは kind 0 を持ち続ける。`abi_v2.md` §4・§9 の「アプリが切り替えるまで残す」の文はエンジン側の文書なので、直すのをこの記録の「合流」で頼む。
4. NCBI の reader が拒む入力（`Near line N`、`>?` の gap の行、Seq-id かもしれない最初の行）は、索引の時点で失敗する（「This input cannot be read: … line N」）。このときレコード表が無く、除外の按は出ない。利用者が file を直す。`checkInput` は拒否の行を引き続きレコードに戻す（アダプタの `Lines` の規則、NCBI の `line_reader.cpp:219-281`）。`screens.spec` の記録 02 は gap の行の file の拒否を示して外す。
5. WP-A：BLASTX は専用の kind（query 1、subject 2）。レコードの無い source はエンジンが検査する（コメントだけの query → `Empty CBlastQueryVector`）。`buildRunInput` は、他のレコードの後に続く defline の無い最初のレコードを拒む（エンジンは前のレコードに続けてしまう）。
6. 抽出の形式（WP-B1）：配列の header は `>{id}:{from}-{to} run={N} {role}_record={K} length={L} unit={nt|aa} hsps={Q.R,...} hit_strand={...}[ requested={a}-{b}]`、整列の header は `>{id}:{start}-{end} run={N} {role}_record={K} hsp={Q.R} aligned[ frame={f}]`、1 行 60 文字、空の ID は `Query_K` / `Subject_K`（outfmt 6 の名前）、`length=` はレコードの長さ、`whole` はレコードごとに 1 配列でトレイの HSP を並べる、frame は翻訳する配列だけ。NCBI の形式ではなく LOSAT Web の形式。`K` で重複 ID を区別する。
7. トレイ（WP-B2）：加えた候補は選ばれた状態、トレイにある HSP は位置とメモを保つ、file 名は `losat-candidates.fa` / `losat-candidates-aligned.fa`（`text/plain`、入力の名前を使わない）、`reveal` は `ResultsState.revealed` に報告、Descriptions の印は表示中の subject だけで query / Run の変更と filter で消える、`sortBy('subject')`＝レコードの label・座標・Run・`s_idx`・index、`sortBy('run')`＝Run の番号・エンジンの index、整列行が 1 つも無い書き出しは失敗して何も保存しない（一部の欠けは一覧に出す）、`extract` は各レコードの ID と長さを Run の snapshot と照らす。
8. W4b の残件（WP-F）：詳細の節は Range の塊の鍵に保存して読み直さない、失敗した詳細に「Try again」（`detail-retry`）で自動の繰り返しはしない、ドットプロットのポップアップは置いてから scroll、電話では入力のカードの名前と大きさを 1 つの折り返す箱に、`.radio input { flex: none }`。
9. 画面（WP-D）：FakeEngine は各組の 1 つ目の HSP に整列行を書く（書き出しと「missing」を試すため）、「In candidates」と端の Up / Down は `aria-disabled`、追加の確認は 4 秒の toast（読み上げ、操作を妨げない）、「Show in results」は Alignments のその Range に focus して `revealed.message` を示す、トレイの Strand は subject の鎖（`hit_strand`）、Run の欄は「Run N」と題、トレイの行はデスクトップで 2 行・電話で 4 行のカード、パネルは mount したまま（選択が残る）、flank の既定は 100/100、並べ替えの select は動かした後「Your order」、ポップアップは最大 24rem（電話は図の中）、Descriptions の選択の強調は印の列まで。
10. コードレビュー（High 無し）：M1 は直した（読みを区間に限る。`sequence-layout.ts` の `readLimit`：一様な layout は区間の終わりの byte、checkpoint は次の checkpoint まで）。L2 は直した（「Add the line >pasted_…」の按は最初の行が Seq-id の拒否のときだけ、`domain/dataset.ts` の `isFirstLineRefusal`）。flank の試験は独立にした。L1（「Complete sequence」の抽出がレコードの写しを 3〜4 つ持つ。250 Mbp で約 1 GB）と L3（行の対応が行 reader の先読みの状態を持たない。混ざった行末で除外の按が隣のレコードを名指し得る）は既知の限界として S15（streaming の Writer、設計書 §12.1）とこの記録に残す。
11. 画面レビュー 1 回目（不合格：M1・M2、L1〜L4）は `608f6f75` で直した：端で切った要約に単位とレコードの長さ（値は折り返さない）、図の下のポップアップは図の列の幅（値は途中で折り返さない、ボタンは 1 行、6 px の余白と 4 px の間隔）、状態 15 は Run 4 のドットプロットを開き直す、電話のトレイは 3.5 行、トレイの Subject / Run の列は値に合わせる（10〜24ch / 6〜24ch）、電話の 24 の印の間隔 4 px。L2：索引が成功した後にレコード単位で拒まれる入力はもう無い（scan と `register` が同じ reader）ので、状態 02 は利用者が除外したレコードを示し、拒否の印は FakeEngine の E2E だけが確かめる。
12. 画面レビュー 2 回目は合格（`608f6f75`、`s14-review-screens`）。低 L-a（電話：トレイの Subject ID が任意の文字で折れる、`overflow-wrap: anywhere`、`styles.css` の約 1800 行）と参考 I1〜I4 は S15 に残す（S13b が最後の低い指摘を残したのと同じく、ゲートをもう一巡しない）。I2：プログラムを替えて役割の kind が変わると索引し直して、その役割の除外が消える（WP-A、判断 5・13）。`screens.spec.ts:84` のコメント（「left out of this and the following searches」）は不正確で、S15 が直す。除外を kind の変更の後も保つべきかは保守者への問い。

### A：読み方の切り替え（`logs/wp-a/report.md`）

13. A-1：プログラムの替えで `indexParser(前, role) != indexParser(後, role)` になるとき、その役割のすべての source を索引し直す（検査と除外は消える）。

### B1：抽出の核（`logs/wp-b1/report.md`）

14. B1-1：一様な layout のレコードも checkpoint のレコードも同じ読み方（開始 byte は式または checkpoint、kind の規則で前向き、チャンクは `readChunkBytes` 以下）。全長の読みは長さと `residue_counts` に照らし、食い違えば「no longer matches its record table … Add the file again.」。
15. B1-2：`readHspRecords` は stream 1 の行を Run ごとに 1 度だけ索引し、index を検査し、index の順でなければ全体の map に落ちる。

### B2：トレイ（`logs/wp-b2/report.md`）

16. B2-1：アプリ層の実エンジンの試験は `extraction-engine.test.ts`（検索の補助を共有）。`CandidateSource` は `candidates.ts` で宣言（`results.ts` は型だけ import）。busy はトレイに 1 つ（同時に出す出力は 1 つ）。

### F：W4b の残件（`logs/wp-f/report.md`）

17. F-1：デスクトップのポップアップの E2E は 1280×520 の平たい図で `toBeInViewport` 0.95 とし、エンジン入りのビルドでは 1 px 未満の端を許す（`541823f0`）。

### G：ゲートの準備（`logs/wp-g/report.md`）

18. G-1：トレイの計測は各測定の自分の記録の後に行う（query：outfmt 6 の数の後。組：写しの loop の後、測定が終わった 1,500 と 3,000）。W4 と W4b の記録の順と名前を保つ。`record.tray = {warmup, samples, summary}` または `{error}`。トレイの失敗は測定の記録を落とさない。
19. G-2：繰り返しのたびに空のトレイから始め、全部を外して終える（外すのも測る）。トレイにある HSP をもう一度加えても何も起きないので、追加は空のトレイでしか測れない。
20. G-3：抽出する 100 候補は、Subject で並べ替えた最初の 100（測らない操作で選ぶ）。抽出は既定（Subject、Hit region、Separate）。
21. G-4：抽出の時計は新しい `extract-summary` までのページの時計（前の要約は置き換わるまで残るので、試験が前のものに `data-measured=before` を付ける）と、行為の前から Playwright の download event までの `downloadMs`（上限）。
22. G-5：200 subject の query は最初の広い query（`query-filter` が選ぶもの）、Descriptions は filter を外した状態。
23. G-6：時計の補助に行為 `select` と条件 `visible`（`getClientRects`。トレイは `v-show` で mount したまま）を足す。
24. fix 1：最初の行の拒否を見分ける matcher は domain に置いて共有する（UI は infra を import できない）。試験の補助 `source()` に任意の `every` を足す。

## 実測

計測は `tests/e2e/results-measure.spec.ts`。1 回の暖機と 3 回の計測の中央値で、ページの中の `performance.now()` による、操作から結果を示す最初のフレームまでの ms（分解能は 1 フレーム、約 17 ms）。W4b と同じ作り方で、括弧の中は W4b の「レビューの後の計測（`087eaa13`）」の数値。ドットプロットと Graphic Summary の時計はどちらも描画のフレームまで測るので比べられる（W4b のゲートの木 `f097595d` の時計は描画を含まず比べられなかった）。**計測の木は `b8641668`**。単一の実行で、機械は他の処理と共有している。`run-20261010T090350Z/measure/`。

**多数の query**：

| query | ブラウザ | 検索（s） | HSP | 開く | 最初の HSP まで | 終わりまで scroll | ID で探す | 絞り込みを外す | 別の query をクリック | With hits only | 200 subject の query をクリック | subject の並べ替え（最大） |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10,000 | Chromium | 0.8（0.8） | 14,981 | 39（38） | 38（41） | 12（15） | 12（12） | 12（12） | 7（8） | 12（13） | 9（9） | 12（12） |
| 10,000 | Firefox | 2.3（2.8） | 14,981 | 56（64） | 63（76） | 2（1） | 11（10） | 11（10） | 18（20） | 11（11） | 5（7） | 11（11） |
| 10,000 | WebKit | 0.9（0.9） | 14,981 | 49（45） | 85（77） | 8（10） | 8（8） | 11（6） | 12（15） | 7（6） | 20（40） | 12（16） |
| 100,000 | Chromium | 4.8（6.7） | 149,828 | 281（292） | 261（266） | 1（1） | 13（12） | 16（13） | 15（14） | 13（12） | 9（10） | 12（12） |
| 100,000 | Firefox | 24.6（21.4） | 149,828 | 487（500） | 494（524） | 1（1） | 9（10） | 14（16） | 22（21） | 13（16） | 15（9） | 12（11） |
| 100,000 | WebKit | 5.3（4.8） | 149,828 | 444（415） | 450（438） | 15（6） | 39（36） | 17（14） | 20（15） | 16（14） | 46（43） | 12（15） |

**1 組の多数の HSP**：

| 写し | ブラウザ | 検索（s） | HSP | 開く | ドットプロットを出す | 拡大 | 縮小 | Zoom to HSP | 全体 | n で選ぶ | クリックで選ぶ | 一覧で選ぶ | HSP の並べ替え |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,500 | Chromium | 3.3（3.8） | 2,995 | 165（181） | 20（20） | 3（3） | 4（3） | 7（6） | 6（7） | 18（7） | 14（14） | 10（9） | 13（13） |
| 1,500 | Firefox | 20.1（17.9） | 2,995 | 310（308） | 41（37） | 27（28） | 26（30） | 25（40） | 26（28） | 16（13） | 31（30） | 18（20） | 12（11） |
| 1,500 | WebKit | 5.3（3.8） | 2,995 | 196（131） | 13（10） | 15（16） | 13（6） | 32（23） | 16（7） | 45（15） | 48（47） | 11（8） | 11（6） |
| 3,000 | Chromium | 12.8（15.3） | 5,993 | 634（672） | 22（20） | 14（4） | 4（9） | 7（7） | 6（6） | 18（16） | 27（14） | 13（16） | 12（13） |
| 3,000 | Firefox | 74.6（75.1） | 5,993 | 1342（1109） | 37（38） | 31（27） | 34（34） | 26（25） | 28（27） | 11（15） | 30（30） | 24（26） | 11（10） |
| 3,000 | WebKit | 失敗：結果が 512 MB の保存の上限を超えた | | | | | | | | | | | |

**候補トレイ**（W4b に対応する数値は無い）：「印」は Descriptions の全部を印す、「追加」は印した行（または組の全 HSP）をトレイに加える、「抽出」は Subject で並べ替えた最初の 100 候補を既定の設定で FASTA にした時間で、括弧は download の event まで（上限）。200 subject の query の 200 候補：

| query | ブラウザ | 印 | 追加 | タブを開く | 終わりまで scroll | 並べ替え（Subject） | 抽出 100（download） | 全部を外す |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 10,000 | Chromium | 12 | 12 | 10 | 13 | 13 | 45（44） | 12 |
| 10,000 | Firefox | 12 | 12 | 14 | 13 | 15 | 28（26） | 12 |
| 10,000 | WebKit | 11 | 11 | 10 | 13 | 12 | 91（89） | 11 |
| 100,000 | Chromium | 12 | 12 | 11 | 13 | 14 | 45（42） | 12 |
| 100,000 | Firefox | 11 | 11 | 15 | 13 | 14 | 27（27） | 12 |
| 100,000 | WebKit | 10 | 9 | 10 | 11 | 11 | 137（141） | 6 |

1 組の全 HSP（写しが 1,500 で 2,995 候補、3,000 で 5,993 候補）。100 配列は 1 つの subject から取る（114,675 B、117,321 B）：

| 写し | ブラウザ | 候補 | 全部を追加 | タブを開く | 終わりまで scroll | 並べ替え（Subject） | 抽出 100（download） | 全部を外す |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1,500 | Chromium | 2,995 | 11 | 13 | 12 | 13 | 44（43） | 11 |
| 1,500 | Firefox | 2,995 | 10 | 14 | 13 | 15 | 27（20） | 11 |
| 1,500 | WebKit | 2,995 | 12 | 11 | 13 | 28 | 289（306） | 11 |
| 3,000 | Chromium | 5,993 | 12 | 13 | 13 | 14 | 45（44） | 11 |
| 3,000 | Firefox | 5,993 | 11 | 14 | 14 | 15 | 27（20） | 12 |
| 3,000 | WebKit | （3,000 写しは保存の上限で失敗） | | | | | | |

- トレイの操作は、200〜5,993 候補のどの場面でも Chromium と Firefox で 1 フレーム（15 ms 以下）。WebKit は 1 つを除いて 1 フレームで、2,995 候補の並べ替えが 28 ms（2 フレーム）。抽出は Data worker の読みが入るので 27〜137 ms で、WebKit の 2,995 候補だけ 289 ms と遅い（原因は調べていない）。各場面で download された file は 100 配列だった。
- W4b より遅く見える値は 1〜3 フレームの差（Chromium の 1 万 query の「終わりまで scroll」12 は W4b の計測 15、W4b のゲートの木 1 の間の揺れ）：WebKit の 10 万 query の「終わりまで scroll」15（6）と「200 subject の query をクリック」46（43。W4b のゲートの木は 18）、Chromium の 1500 写しの「n で選ぶ」18（7）、3000 写しの「拡大」14（4）と「クリックで選ぶ」27（14）、WebKit の 1500 写しの「n で選ぶ」45（15）と「縮小」13（6）、「全体」16（7）、Firefox の 3000 写しの「開く」1342（1109）。WebKit の 1500 写しの「開く」196（131。W4b のゲートの木は 146）は +34〜50%。W4b も同じ種類の揺れを経験している（W4b の「開く時間と一覧で選ぶ時間の調べ」）。単一の実行なので、これ以上追わない。
- W4b より速い値：10 万 query の「開く」（Chromium 281・Firefox 487 ms）と「最初の HSP まで」、Chromium の 3000 写しの「開く」634（672）。検索の時間は、機械の共有による揺れの範囲。
- 10 万 query の Chromium の「Outputs で outfmt 6 の全文を出す」は 2,557 ms（W4b の計測 2,450）、Firefox 13 ms（26）、WebKit 2 ms（4）。Chromium の JS heap（開く前、最初に開いた後、繰り返しの後）は 14.6・43.2・43.5 MB（W4b 14.4・43.0・43.3）。
- **M1 の修正（`a39e19e7`）は計測の木に入っていない**。計測の subject は小さく（257 KB 以下）、チャンク（8 MiB）を読んでも読む量は file の全体に等しいので、M1 が見た費用（大きな 1 レコードから多数の区間を取る）は、この計測では現れない。M1 の効果は `data-extraction.test.ts` が byte 数で確かめる（最後の木で通過）。トレイの計測は `608f6f75` では測り直していない。`608f6f75` の変更（要約の文言、ポップアップの幅、列の幅）は開く経路と抽出の読みに触れない。
- 3000 写しの WebKit と、4000 写しの 3 ブラウザは見込みどおりの失敗（上）。この場面ではトレイを測れない。

## 独立レビュー（コード）

agent `losat-reviewer`（code review の役）が `141615ec..19f4b117` の `web/app` と対応表（54 ファイル、+7037/−482）を読んだ（記録はタスクフォルダの `reviews/code-review.md`）。結論は「妨げになる指摘は無い」。`81b4b6ed`・`b8641668`（ゲートの script と計測の試験）は読んでいない。指摘と対応：

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 抽出が、区間の大きさに関わらず、区間ごとに File から最大 8 MiB を読む（`readInterval`）。100 Mbp の 1 レコードから 1,000 候補を取ると 8 GiB 分の slice になり、指示書の速さの条件（この規模で操作の速さを保つ）に触れる。試験は `readChunkBytes` を小さくしていたので見えなかった | `a39e19e7`：最初の読みを区間で限る（一様な layout は区間の終わりの byte、checkpoint は次の checkpoint）。既定のチャンクで、大きなレコードの短い区間の slice の大きさを数える試験（一様と checkpoint のレコード）。旧コードでは 305000 対 304 で落ちる（判断 10） |
| L1 | 「Complete sequence」の抽出は、UI thread にレコードの写しを 3〜4 つ持つ（worker 側の読み、`fastaRecord`、`concat`、Blob）。250 Mbp で約 0.75〜1 GB。2,147,483,647 文字に近いレコードは抽出できない（「could not be extracted」） | **直さない**（判断 10）。Writer の契約（設計書 §12.1）は S15。この記録の「残件と注意」と S15 の指示書に書いた |
| L2 | 「Add the line >pasted_query」の按が、直らない貼り付けの拒否（gap の行、`Near line N`）にも出る | `e928b64c`：最初の行の拒否のときだけ（`isFirstLineRefusal`）。`search.spec.ts` に gap の行の貼り付けで按が出ない確認（判断 10・24） |
| L3 | Data worker の行番号の対応が、アダプタの `Lines` の先読みの状態（`first_push`、`held`）を持たない。混ざった行末で、除外の按が隣のレコードを名指し得る | **直さない**（判断 10）。影響は除外の按が名指すレコードだけ（メッセージはエンジンのまま、scan が多くを索引の時点で拒む）。一般の場合は `engine-runtime.test.ts` が実エンジンで確かめる |

確かめて指摘の無かったもの：`q_idx` / `s_idx` が指すレコードと `buildRunInput` が同じ順であること、重複 ID の位置での区別、defline の無い最初のレコード、`-query_loc` / `-subject_loc` の座標、abi §9 に従った残基の読み（行の開始の空白・コメント、`;`、kind 1 は ASCII だけ、`U` はそのまま）、座標の module と、それに寄せた `orientation`・`subjectSpan`・Graphic Summary が従来と同値であること、header の形式が判断 6 と一致すること、BLAST の値を `web/` で計算も整形もしないこと、索引の切り替え（役割ごとの kind、`setProgram` の再索引、`checkRecordTable`）、REQ-10 の拒否（追加・抽出・書き出しの全状態）、`reveal` が隠す filter だけを外すこと、印の規則、トレイの線形性、file 名に入力の名前が無いこと、保存・送信を加えていないこと、層の規則、UI が英語であること、W4b の残件の直し、試験が元の配列と突き合わせること（コードの結果と比べない）。

修正（`a39e19e7`、`e928b64c`、`010ba28c`）の後にコードレビューはもう一度行っていない。修正には試験をつけた（M1：slice の大きさ、L2：`dataset.test.ts` と `search.spec.ts`、flank：独立に計算した区間で、`flanked()` を 1 ずらすと落ちることを確かめた）。

## 画面レビュー

agent `losat-reviewer`（screen review の役）が、画面の記録を W4b の最後の記録（`$BUILD_ROOT/s13b-final-screens/`、状態 01〜23）、対応表、前のレビューと並べて見た（記録はタスクフォルダの `reviews/screen-review.md`、`screen-review-2.md`。行ごとの差の script と切り出しは `tmp/review-screens*/`）。

**1 回目：不合格**（High 無し）。ゲートの記録 180 枚（`b8641668`、`$BUILD_ROOT/s14-gate-screens/`。SHA-256 は 180 枚とも `screens.sha256` と一致）。状態 01〜23 は新しい要素と、W4b の残件の修正、実行ごとの値の他は変わらず、W4b の残件（低 1、参考 1）は 3 ブラウザの電話 02・03・06 で直っていた。

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 状態 28：端で切った抽出の要約に単位とレコードの長さが無く、電話で数が「7.4 / kB」「written 1 to / 700」と折れる | `608f6f75`：「requested -2682 to 8241 nt, written 1 to 5000 nt of 5000 nt」。値と単位、「a to b」を一緒にする（判断 11） |
| M2 | 電話の状態 26：TBLASTN のドットプロットのポップアップが図の幅（約 190 px）で、値が「8.83e- / 101」などと単語の途中で折れる。Firefox では E value が別の数に読める | `608f6f75`：図の下では図の列の幅（358 px）、値は折り返さない、ボタンは 1 行。22 も同じ（判断 11） |
| L1 | 状態 15 が W4b の記録と違う内容になった（24〜28 が 15 より先に走る） | `608f6f75`：Run 4 のドットプロットを開き直してから撮る |
| L2 | 状態 02〜05：レコード単位で拒まれる入力・除外されたレコードがどの記録にも無い（gene_x が 19 nt と読まれる） | `608f6f75`：状態 02・03 に利用者が除外したレコード（取り消し線）を撮る。拒否の印は FakeEngine の E2E だけ（判断 11） |
| L3 | 電話の 27・28：トレイの箱がちょうど 3 行で、4 行目が隠れ、scroll の手がかりが無い | `608f6f75`：3.5 行 |
| L4 | 状態 27・28：Run の列に余りがあるのに Subject ID が切れる | `608f6f75`：Subject / Run の列を値に合わせる（デスクトップは完全に直った。電話は L-a が残る） |
| 参考 I1〜I6 | I1：見える変化の一覧（02 の文言、06 のエンジンの文言、12 の build の値、18 の「All matches in candidates」）。I2：gene_x が 19 nt で、5 文字を読み飛ばしたと言わない。I3：電話の Descriptions の印の列で HSPs が見えなくなる（横 scroll の按あり）。I4：電話 24 の選択の棒が 1 つ目の印に触れる。I5：「HSP 1.3」の意味が画面に無い。I6：撮っていない状態（空のトレイ、toast、flank の誤り、単位の混ざった flank の label、鎖が決まらない、整列の「Not written」の一覧、W4b の「Try again」） | I4 は `608f6f75` で直した（印との間隔 3 px）。残りは直さない（I2・I5 は S15 の候補、I6 は下の「残件と注意」） |

**2 回目：合格**（High・Medium 無し）。`608f6f75` の記録 180 枚（`$BUILD_ROOT/s14-review-screens/`、SHA-256 は `run-20261010T100446Z-after-review/screens.sha256` と一致）をゲートの記録と並べた。1 回目の M1・M2・L1〜L4 と I4 は 3 ブラウザ・両サイズで直っていた。状態 01・18 はゲートの記録と同一、04〜14・16・17・20・21・23・25 は実行ごとの値（キューの時間、build の時刻、OPFS の大きさ）の他は変わらない。幅はどの PNG も 1280 または 390 CSS px。

| 指摘 | 内容 | 対応 |
|---|---|---|
| 低 L-a | 電話の 27・28：トレイの Subject ID が任意の文字で折れ、行 1 が「LvMJNV_93001_9800」、行 2 が「0」だけになる（別の ID に見える）。値全体は画面にあり、title がある | **S15 へ**（判断 12）。直すなら `<wbr>` を `_ . | :` の後に置き `overflow-wrap: break-word`、または ID が入らないとき範囲を次の行に。電話の 27・28 を撮り直す |
| 参考 I1 | 要約の長さは「of 5000 nt」と素の数、トレイの title とドットプロットの見出しは「5,000 nt」 | 直さない（`requested=` と同じ書き方。NCBI の `Length=5000` も素の数） |
| 参考 I2 | 状態 04：プログラムを替えると索引し直して除外が消え、通知が無い。`screens.spec.ts:84` のコメントは program を替えるまでしか正しくない | **S15 へ**（コメントを直す。除外を保つかは保守者への問い。判断 12） |
| 参考 I3 | 除外した後の状態 02 に BLASTN の「1 record looks like protein sequences…」の通知が無い（同じ部品は 04 に残る） | 直さない |
| 参考 I4 | デスクトップのトレイで Subject と Range の間が約 8 px | 直さない（列の間隔のまま） |

2 回目の確認：1 回目の指摘が直っていること、検索画面の状態 01〜06 の差が WP-A の読み方の変更と F の修正だけであること、新しい要素が対応表と W4 の基準（触れる対象 24 px 以上、コントラスト、focus の輪、英語、「Run LOSAT」、NCBI のブランド無し）に従うこと。

## 見つけて直したもの

- コードレビュー M1、L2（上）。
- 画面レビュー 1 回目の M1・M2・L1〜L4・I4（上）。
- WP-A：新しい reader は gap の行、`Near line N`、Seq-id の最初の行を索引の時点で拒むので、これらの入力には除外の按が出なくなる（判断 4）。状態 02 の記録を、file 全体の拒否を示す形に替えた。
- WP-D：ポップアップが電話の図の中で 394 px はみ出した（`results.spec.ts` の「narrow screens」が検出）ので、図の下では border-box と `max-width: 100%` にした（その後、画面レビュー M2 で図の列の幅に替えた）。

## 合流

`feature/losat-web-gui` への merge は、このセッションでは行わない。エンジン側が行う。merge のときにエンジン側が行うこと：

1. `origin/feature/losat-web-gui-app` を `feature/losat-web-gui` に merge する（起点 `141615ec` を含む。エンジン側が `141615ec` より先へ進めていれば通常の merge）。エンジン側のパス（`LOSAT/`、`web/adapter/`）、計画、README の表はこのブランチで変えていない。
2. README の表の S14 の行を「完了」にして、この記録へのリンクを置く（計画 §7 の S14 の行と、セッションの README の S14 の行も、W4b のときと同じ扱いで）。
3. `docs/web/abi_v2.md` の §4（`losat_web2_scan_begin` の行）と §9（*scan* の kind 0 の項）は、kind 0 を「アプリが kind 1・2 に切り替えるまで残す」と書いている。アプリは S14 で切り替えたので、エンジン側がその文を書き直す。アダプタは kind 0 を持ったままで、アプリは使わない（判断 3）。残すか外すかはエンジン側が決める。
4. 判断 1〜24（保守者の確認待ち）は、W4・W4b の判断と同じ扱いで、推奨案で進めたと記す。

## エンジン側への申し送り

- **reactor とネイティブの CLI**：`141615ec` の木から作った。アプリの試験（`extraction-engine.test.ts`、`hsp-correspondence.test.ts`、`engine-runtime.test.ts`）は、エンジンの reader（SF）と整列文字列に依存する。SFe など、reader・scan・`register`・HSP のレコードを変える変更の後は、アプリ側で reactor を作り直してこれらを実行する。FakeEngine の `describe.json` は変わらなかった。
- **E2j（S13+）の後のアプリ**：S14 は E2j が merge される前に終わったので、E2j の指示書の 6.（アプリへの申し送り）を行っていない。E2j が merge された後のアプリ側の最初の作業が行う（W4b の記録の「エンジン側への申し送り」の内容のとおり：列の置き場所は対応表に決めてある。Filter Results の Query Coverage、BLASTN の 1 文字の HSP の鎖もそのときに足す）。鎖が分かれば、トレイの「Strand not decided」と抽出の「unknown strand」の扱いが変わる。
- **検証バッジ**：既定の BLASTN は、E2j が既定の outfmt 7 の fixture を足すまで「outside」のまま（W4 の判断 3。S14 で変えていない）。
- **エンジンのメモリ不足**（W4・W4b の申し送りの続き）：20 nt の単位を 4000 写しの自己検索（約 8,000 HSP）が 3 ブラウザで `memory allocation of 49861 bytes failed`。W5 の計測でも再現した。
- **BLASTX**：アプリの座標・単位・frame の規則は合成のレコードの単体試験だけで確かめてある。SX が merge されたら、BLASTX の検索画面（BLASTX 専用の kind は用意済み、判断 5）と、実検索の抽出の試験、ドットプロットの subject の軸と「Query frame」の見え方の確認を足す。

## S15 への申し送り

[S15 の指示書](../../losat_web_gui_sessions/session_s15_w6_export_session.md)の「S14（W5）から引き継ぐこと」に書いた（トレイと研究データの扱い、抽出の形式、Data worker の読み、reader の kind と原 FASTA のつなぎ直し、Writer の契約、画面の記録の場所、ゲートの script と数値、トレイの test ID、残件）。

## 残件と注意

- **コードレビュー L1**：「Complete sequence」の抽出はレコードの写しを UI thread に 3〜4 つ持つ（250 Mbp で約 1 GB）。2,147,483,647 文字に近いレコードは抽出できない。S15 の Writer の契約（設計書 §12.1）で、ブロックを順に書く形にする。それまでは「読みは区間に限る（8 MiB 以下の slice）」であって、メモリに全体を組み立てないとは言えない。
- **コードレビュー L3**：Data worker の行の対応（`message-line.ts`）はアダプタの先読みの状態を持たない。混ざった行末で、除外の按が隣のレコードを名指し得る。メッセージ自体はエンジンのまま。
- **索引の時点の拒否**（判断 4）：NCBI の reader が拒む入力（gap の行、`Near line N`、Seq-id の最初の行）は索引の時点で失敗し、除外の按は出ない。利用者が file を直す。
- **画面レビューの残り**：低 L-a（電話のトレイの Subject ID の折れ）、参考 I2（除外が kind の変更で消える、`screens.spec.ts:84` のコメント）、I1・I3・I4。撮っていない状態（空のトレイ、確認の toast、flank の誤り、単位の混ざった flank の label、「Strand not decided」、整列の「Not written」の一覧、W4b の「Try again」、TBLASTN のポップアップの subject frame、「Cancel the group」、801〜1000 px の幅）は状態 29 以降で足せる。
- **BLASTX の画面**：SX まで動かせない。3 nt per aa の規則と抽出の座標は単体試験だけ（上）。
- **W4b から残る残件**：エンジンのメモリ不足（4000 写し）、WebKit の保存の上限（OPFS の無い WebKit は結果をメモリに置き、512 MB を超える Run は「Not enough temporary storage」で失敗する。W1 の設計どおり。V-MOB の実機で確かめる）、`runInputs` の解放（Run の削除が入るまで解放しない。W4 の判断 5。トレイの削除は Run の削除ではない）、iOS Safari での確認（V-MOB、S17）、Firefox の「開く」+15%（W4b）と 3000 写しの「一覧で選ぶ」の 1 フレーム分、TBLASTN の電話の「Axes not to scale」と小格子の密さ（W4b の画面レビュー 2 回目 L4）、`results_columns.md` の列の表が W4 の並びのまま（E2j の後に対応表と合わせる）、撮らなかった NCBI の画面（Taxonomy、Download と Select columns のメニュー、Edit Search。Edit Search は S15）。
- **計測の対象の木**：計測は `b8641668`。`608f6f75` と、コードレビュー M1 の修正（`a39e19e7`）の後では測り直していない。S15 以降の計測は、今の木から取り直す。
- 画面の記録（PNG）はリポジトリに入れない：ゲートの 180 枚は `$BUILD_ROOT/s14-gate-screens/`、レビューの後の 180 枚は `$BUILD_ROOT/s14-review-screens/`（`/home/kawato/.cache/losat-work/s14-review-screens/`。`608f6f75`、最後のアプリの木。**S15 の画面レビューはこれと比べる**）。SHA-256 は実行記録の `screens.sha256`（180 枚。画面レビューが一致を確かめた）。途中の記録（`s14-wpa-screens`、`s14-wpd-screens`、`s14-wpf-*`、`s14-fix2-screens`）は実行記録を残していない。状態の一覧：01〜06 が検索画面、07〜23 が結果画面（W4b と同じ）、24 が Descriptions の印と「Add to candidates」、25 が Alignments の「In candidates」、26 が TBLASTN のドットプロットのポップアップ、27 が 2 つの Run の候補とメモと Origins、28 が端で flank を切った抽出と要約。
- 背景の作業の記録（担当の報告、レビュー、判断の一覧）はタスクフォルダ `/home/kawato/losat-baselines/s14-w5-20261010/`（リポジトリには入れていない）。
