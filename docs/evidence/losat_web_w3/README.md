# LOSAT Web W3（Session S12）ゲート記録

- 段階：W3 検索画面（[総合計画書](../../losat_web_gui_plan.md) §7 の S12、[指示書](../../losat_web_gui_sessions/session_s12_w3_search_ui.md)）
- ブランチ：`feature/losat-web-gui-app`（アプリ側。worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui-app`）。最初に `origin/feature/losat-web-gui` の `89613b0ec`（S11 の完了と DW-20）を merge した（`c33927577`。衝突なし）。reactor とネイティブの CLI は、その木から作り直した
- 実行記録：[`run-20261006T022027Z/`](run-20261006T022027Z/)（commit `4b3418992` の木。作成後は書き換えない）。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)、再現は [`run_gate.sh`](run_gate.sh)
- 追加の記録：[`run-20261006T024356Z-vbr/`](run-20261006T024356Z-vbr/)（commit `bd63f29ad` の木。TBLASTX の outfmt 0/7 を NCBI の凍結バイトと比べるようにした後の V-BR）と [`run-20261006T025212Z-after-review/`](run-20261006T025212Z-after-review/)（commit `524ea6a3f` の木。2 回目の画面レビューの指摘を直した後の check、E2E、画面の記録）。下の「V-BR の追加の記録」と「レビューの後の実行」
- 判定：**完了条件を満たした**。ゲートの実行（`run_gate.sh`）のすべての段階が通った：V-ABI quick（15 の検索 × 経路とスレッド = 60 の実行がネイティブの CLI と一致、NCBI の凍結ハッシュ 16/16 一致）、`npm run check`（277 件、reactor の要る 28 件は skip）、reactor 付きの単体試験（305 件）、FakeEngine のビルドの E2E（3 ブラウザ、53 件）、エンジン入りのビルドの E2E（3 ブラウザ、92 件。うち検索画面 36 件、V-BR など 30 件）、検索画面の E2E の 2 回の繰り返し（72 件）、多数のレコードの計測（3 ブラウザ）、画面の記録（3 ブラウザ × 2 サイズ × 8 状態）
- 保守者の判断（2026-10-06、このセッションの最初に確認）：program の表示名は **BLASTN 系**（BLASTN・BLASTP・BLASTX・TBLASTN・TBLASTX）。計画 §10 の項目を決定済みにする（下の「合流」で計画 §0.4 に DW-21 として足した）
- 保守者の判断待ち：このセッションで新しく生じたものは無い（W1 の 2 件は変わらず：TBLASTN の R2 の進め方、協調取消の閾値）

## コミット

| コミット | 内容 |
|---|---|
| `c33927577` | `origin/feature/losat-web-gui`（S11、E2d）の merge |
| `5f47cf5aa` | アプリ：検索画面（W3）。`SearchDraft`、入力とレコード一覧、エンジンによる入力の検査、領域、フォーム、Combined / Separate、キュー、Attention、E2E と計測と画面の記録の spec |
| `4b3418992` | アプリ：コードレビューと画面レビューの指摘の対応、計測で見つけた不具合の修正（ゲートの記録の木） |
| `bd63f29ad` | アプリ（試験だけ）：V-BR の `native.ts` の `UNCERTIFIED` から TBLASTX の outfmt 0/7 を外す（S08 が `checked` にした。W1 の合流の項目 5） |
| `524ea6a3f` | アプリ：狭い画面のレコードの行を 2 段にする（2 回目の画面レビューの M-A） |
| `5df008939` | アプリ（試験だけ）：FakeEngine の `describe.json` を単体試験の file snapshot にする（`vitest -u` で書き直す） |
| このゲート記録を含むコミット | 文書：このゲート記録、`run_gate.sh`、実行記録、S13 の指示書への申し送り |

エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表、`docs/web/` は、アプリ側のブランチでは変えていない（下の「合流」で `feature/losat-web-gui` に反映した）。依存ライブラリは増やしていない。

## 完了条件と結果

| 完了条件（計画 §7 の S12） | 結果 | 証拠 |
|---|---|---|
| 研究作業の E2E | 通過（3 ブラウザ、エンジン入りと FakeEngine の両方のビルド、エンジン入りは 2 回の繰り返しも）。Subject を保持したまま Query を変えて繰り返す（2 回目の Run の詳細が「Subject: kept from the previous search」= `subjectRetained`、R1）。キューに複数積む（長い BLASTP の実行中に、別の program・query・option と 2 つの subject の Separate を積み、グループの番号、`Options:`、グループの取消、もう一度積む、取消の後に続く Run の完了）。実行中に次のジョブを編集する（Run 1 は自分の program・入力・option を保つ） | `npm-e2e-engine.log`、`e2e-search-repeat.log`、`npm-e2e-fake-engine.log`（`tests/e2e/search.spec.ts`） |
| 境界条件の E2E | 通過（同上）。空の入力（query、subject の順に理由）、定義行の無い貼り付け（scan の誤りと「Add the line ">pasted_query" before the sequence」）、すべてのレコードの除外、エンジンが拒否するレコード（TD-12 の文言、レコードの番号、キューに入らない）、除外の後の再実行（除外したレコードが出力に無い）、取消（次の検索は新しい runtime で成功し、Subject を登録し直す）、option の拒否（エンジンの文言のまま、キューに入らない）、BLASTX（積めない）、領域の範囲（1〜長さ、2 つ以上のレコードでは指定できない）、狭い画面（横のスクロールが無い）、離れる前の注意と Wake Lock と復帰時の確認 | 同上 |

### V-BR の追加の記録

S08 が TBLASTX の outfmt 0/7 の升目を `checked` にしたので（`docs/web/verification_cells.tsv`）、V-BR の `native.ts` から未認証の升目の一覧を除き、凍結バイトのある升目はすべて一致を求めるようにした（`bd63f29ad`。W1 の合流の項目 5）。`tests/e2e/engine.spec.ts` を 3 ブラウザで実行し、30 件が通った。ブラウザごとに 27 の検索（9 の検索 × 1/2/4 スレッド）の outfmt 0・6・7 と診断がネイティブの CLI と一致し、NCBI の凍結バイトとの比較 21（outfmt 0 が 18、outfmt 7 が 3。TBLASTX の 0/7 の 6 を含む）がすべて一致した（`records/v-br-*.json` の `frozenMatch`）。取消、R1、作り直しの記録も同じ実行にある。

### レビューの後の実行

2 回目の画面レビューの指摘 M-A を直した `524ea6a3f` の木で、`LOSAT_WEB_GATE_STEPS=after-review` の `run_gate.sh`（npm ci、check、単体、両方のビルドの E2E、画面の記録）を実行した。`npm run check`（277 件、reactor の要る 28 件は skip）、reactor 付きの単体試験（305 件）、FakeEngine のビルドの E2E（3 ブラウザ、53 件）、エンジン入りのビルドの E2E（3 ブラウザ、92 件。狭い画面の E2E は、ID が切れず「refused」の印が一覧の中にあることも確かめる）、画面の記録（48 枚、`/home/kawato/.cache/losat-web-gui-target/app-s12-review-screens/`、SHA-256 は `screens.sha256`）がすべて通った。その後の `5df008939` は単体試験と注釈だけを変え、`npm run check` と reactor 付きの単体試験（305 件）を通した。

## 作業との対応

| 作業 | 実装 | 試験 |
|---|---|---|
| 1. タブと入力 | program のタブ（radio。BLASTX は「not available」）。Query / Subject ごとに、貼り付け（Data worker の source になる。300 ms の間の後に索引）、複数ファイル（`<input type=file multiple>`）、ドラッグ＆ドロップ（入力の枠の外に落としたファイルはアプリを置き換えない）。ファイルは textarea に入れず、カードにファイル名・大きさ・レコード数・総長・先頭 6 行（2 KiB まで）・レコード一覧・警告を示す | E2E「a file is summarized…」（setInputFiles と DataTransfer のドロップ） |
| 2. レコード一覧と除外 | 仮想化した一覧（行 30 px、800 px 以下では 2 段の 52 px。ID で絞り込み、行ごとの包含、絞り込んだ行の一括の除外）。重複 ID は `#番号` で区別し、警告に数を出す。不正なレコードはエンジンの `register` の結果（下の判断 3）で、文言の「query record N」を包含したレコードの番号に読み替え、「Exclude record #N」のボタンを出す。配列種別はアプリの推定（文字だけを数え、ACGTUN が 90% 以上なら核酸）で警告だけを出し、program は変えない | E2E「records: duplicate IDs…」「the sequence kind warning…」、単体 `draft.test.ts`、`search-form.test.ts` |
| 3. 領域 | その役割で包含したレコードが全体で 1 つのときだけ、プレビューの帯（ドラッグ）と 2 つの欄（From / To、1 始まり・両端を含む）。アプリは整数と 1〜レコードの長さだけを確かめ（指示書の推奨）、`start-stop` の文法（`5-5`・逆向き）は `validate` の文言を見せる。値は空白と先頭の 0 を除いて `-query_loc` / `-subject_loc` に書く。領域は選んだレコード（source、revision、番号）に属し、レコードが変われば外れる | E2E「the region of a role with one record…」（ドラッグ、欄、範囲の外、argv `-query_loc 3-20 -subject_loc 5-60`）、単体 `search-form.test.ts`・`draft.test.ts` |
| 4. フォーム | `ProgramDescriptor` に節・欄・ラベル（`src/domain/programs.ts`）。既定値（placeholder「default: N」）・help（最初の段落）・遺伝暗号の許可リストは `describe` から。既定値と同じ値と空の欄は argv に書かない。真偽の flag は印を付けたときだけ書く。`describe` の載せない欄は出さない | E2E「the options form…」、単体 `search-form.test.ts`、`engine-runtime.test.ts`（フォームの欄がすべて `describe` にある） |
| 5. Combined / Separate | 役割ごとの切り替え。Separate は source ごとの Run（query と subject の両方が Separate なら組の直積）を、同じグループ ID と「group run i of n」で積み、「Cancel the group」で待っているものから取消す。すべての検索の `validate` が通ったときだけ積む（一部だけは積まない） | E2E「several runs in the queue…」、単体 `coordinator.test.ts`・`draft.test.ts` |
| 6. スレッド・段階・診断 | Auto と 1〜min(16, 論理プロセッサ数)。Run の一覧に段階（Preparing / Searching / Organizing results）と経過時間、詳細に経路とスレッド数（要求の値）、serial の理由、Subject の保持、build、runtime の世代、エンジンのメモリ、段階の時間。query ごとの途中経過は出さない（DW-5） | E2E「several runs…」「the same subject…」「cancel…」 |
| 7. モバイルと離れること | 800 px 以下で縦に並べ（入力、フォーム、ボタン、キュー、注意、保存）、触れる対象を大きくする。Wake Lock は利用者の選択（既定は切）、離れる前の注意と復帰時の確認（下の判断 10） | E2E「narrow screens…」「leaving and returning…」、単体 `attention.test.ts` |
| 8. BLASTN | task 4 つ、template の種類と長さ（task を blastn・blastn-short にすると空にする）、得点（`-reward`・`-penalty`・`-gapopen`・`-gapextend`）、`-word_size`、`-evalue`、`-dust`、`-lcase_masking`、`-perc_identity`、`-subject_besthit` など。既定値は書かず、拒否は `validate` の文言のまま。`register` のハンドルは program ごと（エンジンの R1 の鍵が program を含み、program を変えると登録し直す）。TD-12 の拒否はレコードの警告として見せる | E2E「the options form…」（template の空、`-word_size 3`、blastn の template の拒否）、「records…」（`not supported by LOSAT's BLASTN`） |
| 9. BLASTP・TBLASTN・TBLASTX | program ごとの欄（行列、gap、`-comp_based_stats`、`-seg`、`-soft_masking`、`-db_gencode`、`-query_gencode`、`-xdrop_gap*`、`-sum_stats`、`-culling_limit` など）。NCBI と LOSAT の拒否（toolkit の語を含む）は `validate` の文言のまま、キューには入らない。入力は S10 の `DatasetStore`（`addSource`・`indexSource`・`reviseDataset`・`buildRunInput`）。BLASTX は SX まで「not available」（判断 2）。多数のレコードの索引とハッシュの時間を測った（下の「実測」） | E2E「the engine's refusals…」、計測 `search-measure.spec.ts` |
| 10. E2E | 上の「完了条件と結果」 | `tests/e2e/search.spec.ts`（12 件 × 3 ブラウザ） |

## 作ったもの

- `src/application/draft.ts`（`SearchDraft`）：検索画面の下書き。program、役割ごとの入力（貼り付けと source の一覧、Combined / Separate、領域）、program ごとの欄の値、スレッド、`describe` と `validate` の結果。索引、エンジンの検査、`validate` は、入力が止まって 300 ms の後に行い、古い結果は捨てる。`submit()` は待っている仕事を終えてから、一つの状態から要求を作り、`Coordinator.enqueueAll` に渡す。「Add to queue」を止める理由を短い文で示す（`readiness()`）。
- `src/application/coordinator.ts`：`enqueueAll`（全部の `validate` の後に積む。2 つ以上ならグループ）、`cancelGroup`、入力が `DatasetRevision` の ID の列のとき `buildRunInput` で作り、同じ ID の列の入力は同じバイトを共有する。
- `src/application/attention.ts`（`Attention`）と `src/ports/page.ts`・`src/infra/browser/page.ts`：Wake Lock、`beforeunload`、`visibilitychange`。
- `src/domain/`：`programs.ts`（`ProgramDescriptor` の節と欄、説明、BLASTX の `unavailable`）、`parameters.ts`（`describe` と欄から選択肢と argv を作る規則）、`region.ts`、`sequence-kind.ts`、`genetic-codes.ts`（NCBI の `gc.prt` の名前）、`argv.ts` の入力名。
- `src/ports/input-check.ts`、`src/infra/reactor/checker.ts`（`ReactorInputChecker`）、`DataService.checkInput`・`previewSource`、`src/infra/fake/describe.json`（FakeEngine の `describe`。エンジンの写し）。
- `src/ui/`：`SearchPanel`、`InputPanel`、`SourceCard`、`RecordList`、`RegionPicker`、`ParameterForm`、`QueuePanel`（作り直し）、`AttentionPanel`、`ResumeNotice`、`styles.css`（作り直し）。
- `vite.config.ts`：開発とプレビューのサーバーが、304 を含むすべての応答に COOP / COEP / CORP を付ける（判断 12）。
- 試験：単体 `draft.test.ts`（21）、`search-form.test.ts`（64）、`attention.test.ts`（29）、`coordinator.test.ts`・`data-service.test.ts`・`engine-runtime.test.ts` への追加。E2E `search.spec.ts`（12）、計測 `search-measure.spec.ts`、画面の記録 `screens.spec.ts`。
- `describe.json` は `tests/unit/engine-runtime.test.ts` の file snapshot（Vitest の `toMatchFileSnapshot`）：reactor があるときに `describe` と比べ、エンジンの option を変えたら `LOSAT_WEB_REACTORS=<dir> npx vitest run tests/unit/engine-runtime.test.ts -u` で書き直す（TypeScript の reactor の結合をそのまま使う）。
- 文書：[`run_gate.sh`](run_gate.sh)（`LOSAT_WEB_GATE_STEPS=after-review` で、レビューの後の実行の段階だけ）。

## 判断（推奨案で進めたもの）

1. **表示名**：保守者の判断どおり BLASTN 系（DW-21）。
2. **BLASTX は SX まで「not available」**：タブは残し（5 program の画面の形を保つ）、選ぶと「BLASTX joins LOSAT Web after its certification. Until then, use BLASTX of the LOSAT command line.」を示し、ボタンは理由とともに使えない。入力と欄は残る（他の program に戻ると使える）。BLASTX の検査と索引は行わない（レコード一覧は BLASTN の解析器の種類で作る）。
3. **不正なレコードの判定はエンジンの `register`**（計画 §5.4）：Data worker の serial の reactor で、program と役割と revision ごとに `register` してすぐ `release` する。拒否（export の −1）はエンジンの文言のままレコード一覧に出し、キューに入れない。結果は (program, 役割, revision) ごとに覚え、同じ選び方は一度だけ検査する。入力が reactor のメモリに入らない（`losat_web2_alloc` が 0）や instance が止まったときは「検査できなかった」として示し、キューは止めない（実行が自分の誤りを返す。コードレビュー M2）。検査は SHA-256 を計算しない。TS で検証の規則を作らない。
4. **貼り付けは Data worker の source**：S10 では `enqueue` が索引を作り、貼り付けは File と snapshot の 2 つの写しだった。貼り付けは入力が止まって 300 ms の後に source にし、索引と検査は file と同じ経路を通る。定義行の無い貼り付けは scan の誤りを見せ、「Add the line ">pasted_query" before the sequence」のボタンを出す（アプリが勝手に足さない）。
5. **入力名**（計画 §5.3）：貼り付けは `query.fa` / `subject.fa`、ファイル 1 つはその名前、Combined で 2 つ以上は `combined_query.fa` / `combined_subject.fa`（計画は `combined_subject.fa` だけだった。query も同じ規則にした）。Run の一覧の貼り付けの名前は「Pasted sequences (query.fa)」。
6. **Separate の組**：query と subject の source の組の直積（両方が Combined なら 1 つ）。2 つ以上ならグループにし、取消は待っているものを先に取消してから実行中のものを取消す（実行中の終わりが次のメンバーを始めないように）。一つでも `validate` が拒否すれば何も積まない（番号も進めない）。
7. **R1**：除外の組が同じなら同じ revision を使う（`reviseDataset` を (base, 除外) ごとに覚える）。Subject を変えずに Query を変えると、エンジンは Subject を保持する（`subjectRetained` = true。W1 の判断 9 の条件を満たした）。
8. **領域**：指示書の推奨どおり、欄も 1〜レコードの長さに制限する。`start == stop` と逆向きは `validate` の文言（NCBI の誤り）を見せる（アプリで BLAST の規則を作らない）。帯のクリックだけでは領域を作らない（1 文字の範囲になるため。コードレビュー L6）。領域は選んだレコードに属し、レコードが変わると適用しない（試験で見つけた。下の「見つけて直したもの」）。
9. **フォームの規則**：task の選択肢（4 つ）と template（種類 3、長さ 16/18/21）は指示書の一覧から（`describe` は help の文にだけ書く）。task を blastn・blastn-short にすると template の 2 欄を空にする。値が `describe` の既定値と同じ（前後の空白を除いて）なら書かない。task に依存する既定値（blastn-short の e-value 1000 など）は `describe` に無いので、明示した 10 はそのまま書く（NCBI と同じ結果）。TBLASTX の `-max_target_seqs` は `describe` に既定値が無いので placeholder は「default」。program を切り替えたときは、その program の `describe` を同期して適用してから argv を作る（既定値を書かないため）。
10. **離れることへの備え**（設計書 §9.4、`REQ-22`）：Wake Lock は利用者が選び（既定は切）、Run が実行中か待っている間、ページが見えている間だけ持つ（キューの Run の間で放さない。コードレビュー L4）。ブラウザが放したら、見えるようになったときに取り直す。Run が実行中か待っている間は、離れる前の注意（`beforeunload`）。実行中か待っている Run があるまま 1 秒以上隠れたページが戻ると、隠れていた時間、Run ごとの前と今の状態、Data worker の応答（10 秒で答えなければ「not responding」、遅れた答えは直す）を示す。LOSAT Web は、タブを閉じても検索が続くとは言わない。
11. **`validate` はその場で**：欄を変えて 300 ms の後に Data worker の reactor が argv を確かめ、文言をそのまま見せる（live region）。「Add to queue」でも `Coordinator` がもう一度確かめる。
12. **304 にも COOP / COEP / CORP**：WebKit で、取消の後に Engine worker を作り直すと、`vite preview` の 304 に COEP が無く、worker の読み込みが拒否された（「Worker load was blocked by Cross-Origin-Embedder-Policy」）。開発とプレビューのサーバーで、すべての応答にヘッダーを付ける plugin を足した。本番の Cloudflare の `_headers` が 304 に付くかは S16 で確かめる（下の「合流」で S16 の指示書に書いた）。
13. **レコード表はメモリのまま**（計画 §2.3）：下の「実測」のとおり、10 万レコード（34 MB）の索引は 3〜10 秒、表は 1 つの写しで約 51 MiB（Node の V8 での見積り、レコードあたり約 536 バイト。大部分は残基の内訳）。表は Data worker と画面の 2 つにある。保存（`index.blocks`）に移す条件は満たさないと判断した。S17 で公開する対応規模が 10 万レコードを大きく超えるとき、または V-MOB の実機で負担が見えたときに見直す。
14. **画面の記録はリポジトリに入れない**：ゲートの PNG 48 枚（合計 31 MB）は `/home/kawato/.cache/losat-web-gui-target/app-s12-gate-screens/`、撮り直しの 48 枚は `app-s12-review-screens/` に置き、SHA-256 をそれぞれの実行記録の `screens.sha256` に残した（大きなファイルを commit しない規則）。

## 実測

計測は `tests/e2e/search-measure.spec.ts`（ゲートの実行の `measure/`）。BLASTP の query に、300 残基・60 文字の行の蛋白の FASTA を N レコード選び、アプリを通して、1 回の暖機と 3 回の計測の中央値（ms）。「待ち時間」は試験から見た時間、「アプリ内」は下書きが Data worker の呼出しの前後で測った時間。

- 索引：ファイルを選んでからレコード表が出るまで（Data worker の serial の reactor の `scan` と、レコードごとの SHA-256）
- 検査：エンジンの `register`（Data worker の reactor）
- 除外：1 レコードを除いてから、新しい revision の検査が終わるまで（入力が止まって検査を始めるまでの 300 ms を含む）
- 積む：「Add to queue」から Run が一覧に出るまで（実行の入力とその SHA-256、`validate`）

| レコード数（大きさ） | ブラウザ | 索引（待ち時間） | 索引（アプリ内） | 検査（アプリ内） | 除外（待ち時間） | 除外の後の検査（アプリ内） | 積む（待ち時間） |
|---|---|---:|---:|---:|---:|---:|---:|
| 1,000（0.3 MB） | Chromium | 62 | 17 | 6 | 831 | 8 | 63 |
| 1,000（0.3 MB） | Firefox | 107 | 83 | 30 | 859 | 33 | 55 |
| 1,000（0.3 MB） | WebKit | 149 | 104 | 23 | 927 | 20 | 105 |
| 10,000（3.4 MB） | Chromium | 380 | 109 | 43 | 827 | 59 | 58 |
| 10,000（3.4 MB） | Firefox | 1,382 | 875 | 312 | 869 | 357 | 121 |
| 10,000（3.4 MB） | WebKit | 1,006 | 673 | 100 | 923 | 136 | 141 |
| 100,000（34.1 MB） | Chromium | 3,064 | 1,140 | 425 | 1,039 | 604 | 191 |
| 100,000（34.1 MB） | Firefox | 9,605 | 8,078 | 3,078 | 3,982 | 3,267 | 224 |
| 100,000（34.1 MB） | WebKit | 7,651 | 6,223 | 753 | 1,958 | 1,154 | 633 |

- 1 万レコード（3.4 MB）までは、どのブラウザでも索引は 1.4 秒以内、積むのは 0.15 秒以内。10 万レコード（34 MB）では索引が 3〜10 秒（ファイルを選んだときに 1 回）、除外の後の検査が 0.6〜3.3 秒、積むのは 0.2〜0.6 秒。Firefox は、W1 と同じくこの機械で Wasm が遅い（W1 の「実測」）。
- 試験から見た時間には、Playwright の確かめの間隔（100・250・500 ms…）が入る。1,000 レコードの除外が約 830 ms なのは、300 ms の間とこの間隔による。アプリの時間はアプリ内の列で読む。
- 10 万レコードの入力の検査は、このゲートの前には失敗していた（下の「見つけて直したもの」）。この計測が見つけた。
- 同じ機械で他の作業が並行していた時間がある（32 論理プロセッサ）ので、絶対値には揺れがある。

## 独立レビュー（コード）

`5f47cf5aa` を、別のエージェント（Sonnet）が読み取り専用でレビューした（判定「ready after fixes」、高い指摘は無し）。対応はすべて `4b3418992`。

| 重さ | 指摘 | 対応 |
|---|---|---|
| 中 M1 | `submit()` が例外を捨てる（worker が止まっていると、何も示さずに終わる） | 失敗を文で示す。`Coordinator` は投げた `validate` を「The options could not be checked: …」にする。単体試験 |
| 中 M2 | 入力の検査が重い（全体の SHA-256、除外のたびの検査）。確保の失敗を「拒否」と示し、検査できないだけでキューを止める | 判断 3：SHA-256 を計算しない、(program, 役割, revision) ごとに覚える、`EngineAllocationError` は検査の失敗、検査の失敗はキューを止めない。単体試験 |
| 低 L1 | `prepareAndEnqueue` が二つの時点の状態を読む | 最初の待ちの前に、一つの状態から要求を作る |
| 低 L2 | `idle()` の待ちに表示が無い | 待つ間「Preparing the inputs…」を示す |
| 低 L3 | 間（300 ms）の間に古い `validate` の結果を示すことがある | 予定した時点で世代を進める |
| 低 L4 | Wake Lock が Run の間で放され、取り直される | 判断 10：実行中か待っている間は持つ |
| 低 L5 | 遅れた Data worker の答えが無視される。短い隠れでも報告する | 遅れた答えで直す。1 秒未満の隠れは報告しない |
| 低 L6 | 帯のクリックで 1 文字の領域ができる | 判断 8 |
| 低 L7 | `-matrix` の候補と `-comp_based_stats` の選択肢が固定の一覧 | 変えない（候補と選択肢で、計算した値ではない。LOSAT の拒否は `validate` が文言で示す）。残件 |
| 低 L8 | `web/AGENTS.md` の規則 5（ui は domain の型と定数だけを import）と、ui が domain の関数を使うこと（S01 の `ResultsPanel` から） | 変えない。規則の文を実際に合わせるかを残件にした（ESLint は ui → infra だけを禁じている） |
| 細 N1〜N6 | `runInputs` の解放、入力名の重複、使っていない再索引、`dragleave`、live region、`datalist` の `v-if` の連なり、復帰の文 | 入力名を domain に移した、live region、`<template>` で分けた、復帰の文を直した、枠の外のドロップ。`runInputs` は Run の削除ができたとき（S13 以降）に見直す |

レビューで正しいとされた点：層の規則、web/ で BLAST の値を計算しないこと（領域は数字と 1〜長さだけ、種別は警告だけ）、既定値を書かないこと、task による template の消去、予約した引数、DW-9 の条件、レコード番号の読み替え、古い結果の破棄、revision の再利用、グループの番号と取消の順、Attention の競合、仮想化した一覧、信頼しない文字列を `v-html` にしないこと、FakeEngine の文言と帯、304 の plugin。

## 画面レビュー

画面の記録（`tests/e2e/screens.spec.ts`。1280×900 と 390×844、Chromium・Firefox・WebKit、全ページ）を、別のエージェント（Sonnet）が見た。前の記録は無い（W3 が最初の検索画面）。

1 回目（`4b3418992` の前の記録、8 状態）：判定「pass with minor findings」（高い指摘は無し、はみ出し・重なり・横のスクロールは無し）。対応は `4b3418992`。

| 重さ | 指摘 | 対応 |
|---|---|---|
| 中 M1 | 拒否されたレコードが含まれたままでも「The engine accepts these options.」が緑で出て、ボタンが使える | ボタンの上に、積めない理由（拒否されたレコード、入力が無い、すべて除外）を 1 行で示す。押すと同じ理由を示し、積まない |
| 中 M2 | 狭い画面でキューがずっと下にあり、押した結果が見えない | ボタンの横にキューの要約（「Run 4 is running, 2 waiting」など）と「Show the queue」のリンク |
| 中 M3 | 狭い画面でレコードの行の ID が詰まり、長さの列が揃わない | ID と長さの列を固定した（印の折り返しが無く、2 回目の M-A になった） |
| 中 M4 | WebKit の「メモリに保存」の注意が細かい灰色の字 | 警告の色と平易な文（再読込で消える） |
| 中 M5 | 選んだ BLASTX のタブが使えない状態に見える。ボタンの理由が離れている | 選んだ状態の見た目、ボタンの近くに理由 |
| 中 M6 | 狭い画面の触れる対象が小さい | 800 px 以下で checkbox・radio・帯・タブ・ボタンを大きくした |
| 低 L1〜L3、L5、L7、L8 | 復帰の文、領域の説明、遺伝暗号の注、placeholder、Run の状態の色、貼り付けの名前 | 直した |
| 低 L4、L6、L9、L10 | 詳細の列の幅、空の状態、狭い画面の見出し、ブラウザごとの部品の見た目 | L6 は M1 の行で対応。他は変えない（診断の表示、ブラウザの部品） |
| 未解決 | 待っている・失敗・取消の Run の記録が無い | 状態 06 を足した（実行中、待ち、取消、失敗の Run と復帰の報告） |
| — | WebKit の「Keep the screen on」が理由なく使えない | 誤読だった（2 回目のレビューが確かめた）。WebKit 26.5 には Screen Wake Lock があり、選択肢は使える（灰色は WebKit の部品の見た目）。E2E「leaving and returning」は 3 ブラウザで Wake Lock を持つ。API の無いブラウザでは、選択肢の代わりに「This browser cannot keep the screen on.」を示す |

2 回目（ゲートの記録 `app-s12-gate-screens/`、8 状態 × 2 サイズ × 3 ブラウザ = 48 枚、1 回目の記録と比べた）：判定「pass with minor findings」（高い指摘は無し、24 枚の狭い画面のすべてで横のはみ出しが無い）。1 回目の指摘は M1・M2・M4・M5・L1・L2・L7・L8 と未解決の状態（待ち・失敗・取消の Run）が解決、M6 はおおむね解決、L3・L5 は一部。新しい指摘と対応：

| 重さ | 指摘 | 対応 |
|---|---|---|
| 中 M-A | 狭い画面のレコードの行で、ID が「contig…」と切れ、「refused」の印が一覧の端の外に出る（M3 が満たされず、1 回目より悪い。固定の高さの行で印が折り返せない） | `524ea6a3f`：800 px 以下では行を 2 段（レコード、その下に印）にし、仮想化した一覧も高い行を使う。狭い画面の E2E が、ID が切れず印が一覧の中にあることを確かめる（前の配置では失敗することを確かめた）。撮り直して 3 回目のレビューに見せた |
| 低 L-a | Subject の遺伝暗号の select が、閉じた状態で名前の終わりを切る | 変えない（開くと全文。残件） |
| 低 L-b | 拒否されたレコードを含むままでもボタンが押せる。押した結果の記録が無い | 変えない。押すと理由を示して積まないことを E2E「records…」が確かめる（`search-message` と、Run が増えないこと） |
| 低 L-c | WebKit の「not supported」の表示が記録に無い | 上の表の最後の行（WebKit は API を持つ） |
| 低 L-d | 行の高さ 30 px、文字だけのリンク | 行は M-A で 52 px（狭い画面）。リンクは変えない |
| 低 L-e | 実行中の Run が一覧の下にある。待ちと取消の印が同じ灰色 | 変えない（ボタンの横の要約とリンクがある）。S13 の Run の一覧の作り直しで扱う（S13 の指示書） |
| 低 L-f | 「default」だけの placeholder、help の `(default: …)` と `[default: …]` の違い、欄の高さの不揃い | 変えない（`describe` に既定値の無い欄は「default」。help はエンジンの文のまま） |

3 回目（`524ea6a3f` の撮り直し `app-s12-review-screens/`、上の「レビューの後の実行」、2 回目の記録と比べた）：判定「pass」。M-A は 3 ブラウザで解決（行は 2 段、ID は全部、長さの列は揃い、「refused」の印は一覧の中）。デスクトップの大きさのレコードの行は変わらない（印が約 2 px 左に動いただけ）。退行は無い：狭い画面でレコード一覧のある状態は、行の高さの分（3 行 × 22 px = 66 px）だけ縦に伸び、一覧の上は画素まで同じ、下はその分だけ動く。24 枚の狭い画面のすべてで横のはみ出しが無い。低い指摘：印の無い行も 52 px で 2 段目が空く（仮想化した一覧の固定の高さのため。変えない）、2 回目の L-a・L-b・文字だけのリンク（上の表のとおり変えない）。

## 見つけて直したもの

- **WebKit の 304**（判断 12）：エンジン入りの E2E の最初の実行で、WebKit の取消の後の検索が 3 件失敗した。
- **10 万レコードの入力の検査が「Maximum call stack size exceeded」**：`buildRunInput` がレコード表を `push(...array)` に広げていた（約 12 万レコードから）。`Coordinator` のキューも同じ。ループにし、15 万レコードの単体試験を足した（修正の前に失敗することを確かめた）。
- **狭い画面のレコードの行**：2 回目の画面レビューの M-A（上の「画面レビュー」）。
- **領域が別のレコードに残る**：subject を変えても前の領域が残り、キューを止めていた。領域をレコードに属させた。
- 試験の側：計測の擬似乱数の精度（ほぼ `A` になっていた）、Run の状態の locator が状態の span にも当たっていた、WebKit の既定の時間切れ（30 s）が短い、thread host の単体試験の時間の揺れ（負荷の高いとき）。

## 保守者の判断待ち

このセッションで新しく生じたものは無い。W1 の 2 件（TBLASTN の R2 の進め方、協調取消の閾値）は変わらない（[W1 のゲート記録](../losat_web_w1/README.md)の「保守者の判断待ち」）。

## 合流（エンジン側の作業がこの時点で無いので、このセッションがコーディネーターとして行った）

`feature/losat-web-gui-app`（`d4503b355`）を `feature/losat-web-gui` に merge し（`57b4b4330`）、W1 と W3 のゲート記録の「合流のときにエンジン側が行うこと」を、その次の文書のコミットで反映した（このファイルのこの行もそのコミットで書いた）：

1. README の表：S09 を完了（2026-10-03）、S09+ を R2 TBLASTN（条件付き、保守者の判断 1 の後）、S12 を完了（2026-10-06）にした。
2. 計画：状態の行、§0.4 に DW-21（表示名）、§7 の S09・S12 を完了、§10 の表示名と instance の作り直しを決定済み・取消の閾値を判断待ちに、§2.3（Memory の上限 512 MB、レコード表の実測、取消の後の再準備の実測）、§3.2 のコード配置、§5.3 の `combined_query.fa`、§5.5 の Auto と作り直しの値と S12 の Attention、§8 のリスク（Chromium の共有メモリの境界、WebKit の 304）。
3. `docs/web/verification_cells.tsv`：W1 の V-BR のブラウザの升目。
4. `docs/web/abi_v2.md`：§5 に stream 2 も 1 MiB ごとに届くこと、§9 に空白だけの入力はアプリが `scan` の誤りで拒否すること（W1 の判断 8。`scan` は変えない）。
5. W1 の項目 5（`native.ts` の `UNCERTIFIED`）は、このブランチで行った（`bd63f29ad`）。
6. S16 の指示書：本番の 304 に COOP / COEP / CORP / CSP が付くこと、WebKit / Safari で取消の後の検索が動くことを確かめる。

## S13 への申し送り

[S13 の指示書](../../losat_web_gui_sessions/session_s13_w4_results_ui.md)の「S12（W3）から引き継ぐこと」に書いた（画面の構成と `ResultsPanel` の置き換え、Run の識別とグループ、E2E の補助、`describe.json` の書き直し、キューから結果を開くことと Run の一覧の残りの指摘、Run の入力の解放、多数のレコード、画面の記録の場所、ゲートの script）。

## 残件と注意

- `web/AGENTS.md` の規則 5 の文（ui は domain の型と定数だけ）と実際（ui が domain の純粋な関数を使う。S01 から）の食い違い（コードレビュー L8）。文を実際に合わせるか、関数を application から渡すかを、次に規則を見直すときに決める。
- `-matrix` の候補と `-comp_based_stats` の選択肢は `programs.ts` の固定の一覧（L7）。`describe` が選択肢を返すようになれば、それに替える。
- `describe` の help には NCBI の「Reference:」などの文が含まれ、フォームは最初の段落だけを示す。
- 除外のたびに、Data worker は新しい revision をレコード表ごと画面に返す（10 万レコードで数十 MiB の受け渡し）。画面は ID だけを使うので、負担になれば `reviseDataset` の応答を ID だけにする。
- `Coordinator.runInputs`（Run の入力のバイト）は解放しない。Run の削除を入れるとき（S13 以降）に見直す。
- 本番の 304 のヘッダー（S16）。WebKit は Playwright の WPE の MiniBrowser で試した（Safari と iOS は V-MOB と S17）。Playwright の WebKit には OPFS が無い。
- Subject の遺伝暗号の select は、閉じた状態で長い名前の終わりを切る（画面レビュー 2 の L-a）。待ちと取消の Run の印が同じ灰色（L-e、S13 の指示書）。狭い画面では印の無い行も 2 段の高さ（画面レビュー 3 の L-1）。
- 離れる前の注意（`beforeunload` のダイアログ）は Chromium だけで E2E を行った（Firefox と WebKit の Playwright は `runBeforeUnload` のダイアログを同じように出さない）。
- V-PERF の lock（`vperf.lock`）は、このセッションの間に一度も現れなかった。
