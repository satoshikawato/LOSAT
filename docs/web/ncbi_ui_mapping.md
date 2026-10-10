# NCBI BLAST の画面と LOSAT Web の対応表

段階 W4b（Session S13b、[指示書](../losat_web_gui_sessions/session_s13b_w4b_ncbi_style_ui.md)）で、保守者の指示（2026-10-09）「NCBI BLAST ウェブサイトの GUI にデザインを寄せてほしい」を具体にする。設計書 §2.1 の「NCBI Web BLAST に馴染みのある操作」と §11.2 の「NCBI Web に準じる入力項目・列を固定した参照画面と対応付ける」の参照画面と対応表である。

- 作成：2026-10-10（S13b）。採否は保守者の確認待ち（推奨案で進めた。Owner-delegated、[W4b のゲート記録](../evidence/losat_web_w4b/README.md)の判断）。
- 追記：2026-10-10（S14）。Descriptions の行の印、候補への追加（Descriptions・Alignments・Dot Plot）、候補トレイ（§2「候補（S14）」、§3）。採否は同じく保守者の確認待ち（Owner-delegated、W5 のゲート記録 `docs/evidence/losat_web_w5/README.md` の判断）。
- 追記：2026-10-10（S15、W6 の WP-H）。結果画面の最初のタブを NCBI の classic（Traditional）結果ページの 1 ページにした（Graphic Summary、その直下に Descriptions、その下に Alignments。§2「タブ」「1 ページの結果（classic）」）。保守者の指示（2026-10-10）「結果画面はClassic版というか、ヒットの座標とスコアの分布がダーッと出る絵が一発で見えるようにしてほしい」「Graphic Summaryの下にDescriptionが直接出てくるやつ。いちいちタブを開かないといけないのはめんどくさい」「全部一斉に開く必要はないんだよ。NCBIも適宜カットしてるでしょ？」による（S15 の DECISIONS 12）。S13b のタブ式の対応（下の「タブ」の履歴）はこれで置き換えた。
- 寄せるのは配置・節の順・言葉（見出し、欄の名前、ボタン、タブ、列、凡例）。NCBI のロゴ、NIH / NLM のヘッダーとフッター、「BLAST®」の文字、配色は写さない。LOSAT Web が NCBI のサービスに見えないようにする。
- 値はエンジンが書いたまま示す。TS で BLAST の値を計算・整形しない（`web/AGENTS.md` 規則 1、[`results_columns.md`](results_columns.md)）。表示の分類（Graphic Summary の得点の階級、ドットプロットの identity の階級）は、エンジンの原値を区切るだけで、値を作らない。
- 採否の語：**採る**（NCBI と同じ場所・同じ言葉）、**寄せる**（同じ場所に LOSAT の要素を置き、言葉を近づける）、**採らない**（理由を書く）、**エンジン待ち**（S13+（E2j）が値を出すまで出さない）、**後の段階**（その段階の名）。

## 参照画面

2026-10-10 に Chromium（Playwright、headless）で NCBI BLAST の画面を撮った。検索は公開の RefSeq の配列だけで行った（NM_000518.5、NM_000519.4、NM_000558.5、NP_000509.1、NP_000510.1）。検索は 6 件（RID CJ5K9KV3114 blastn、CJ5W8ZCC114 blastp、CJ6FDRHH114 tblastn、CJ6G5X6T114 tblastx、CJ6GXH5H114 と CJ7C8J2Y114 は 2 query）。画像はリポジトリに入れず `$BUILD_ROOT/s13b-ncbi-reference/`（`/home/kawato/.cache/losat-work/s13b-ncbi-reference/`）に置き、URL・時刻・RID・SHA-256 をその `INDEX.md` と [W4b のゲート記録](../evidence/losat_web_w4b/README.md)に残した。

| 画面 | ファイル（1280 px。†は 390 px も） |
|---|---|
| 検索画面（「Align two or more sequences」に印、BLAST_SPEC=blast2seq） | `search-blastn-1280`†、`search-blastn-params-1280`†（Algorithm parameters を開いた）、`search-blastn-dc-params-1280`（discontiguous megablast）、`search-blastn-nondefault-1280`（既定と違う値）、`search-blastp-*`、`search-tblastn-*`、`search-tblastx-*` |
| 2 配列の結果（blastn・blastp・tblastn・tblastx） | `<program>-results-top`、`-search-summary`、`-descriptions`、`-graphic`、`-graphic-hover`、`-alignments`、`-dotplot`、`-dotplot-image` |
| 複数 query の結果 | `multi-*`（「Results for」） |

## 1. 検索画面

NCBI の 2 配列の画面は、上から program のタブ、program の一文、「Enter Query Sequence」、（「Align two or more sequences」）、「Enter Subject Sequence」、「Program Selection」、BLAST のボタンと一文、開閉する「Algorithm parameters」、もう一度 BLAST のボタン、の順である。LOSAT はこの順に並べ直す（W3 は query と subject を左右に並べていた）。

| NCBI の要素 | LOSAT Web | 採否 | 理由・置き場所 |
|---|---|---|---|
| program のタブ「blastn \| blastp \| blastx \| tblastn \| tblastx」 | 同じ並びのタブ。表示名は BLASTN・BLASTP・BLASTX・TBLASTN・TBLASTX | 寄せる | 表示名は保守者の決定（DW-21）。BLASTX は SX まで「not available」（W3 の判断 2） |
| 画面の題「Align Sequences Nucleotide BLAST」など | なし | 採らない | NCBI の製品名の見出しで、ブランドに近い。program のタブと一文で足りる |
| program の一文（「BLASTN programs search nucleotide subjects using a nucleotide query.」） | タブの下の一文（「BLASTN searches nucleotide subjects using a nucleotide query.」など） | 寄せる | NCBI の 2 配列の画面の一文の形 |
| 「Reset page」「Bookmark」 | なし | 採らない | Reset page は開いた大きなファイルの選択も消す（入力は Remove で 1 つずつ外す）。Bookmark は URL に検索の状態を持たない（研究データを URL に入れない、`web/AGENTS.md` 規則 4） |
| 枠「Enter Query Sequence」 | 枠「Enter Query Sequence」（`query-panel`） | 採る | |
| 「Enter accession number(s), gi(s), or FASTA sequence(s)」と textarea、「Clear」 | 「Enter FASTA sequence(s)」と textarea（`query-input`）、「Clear」（貼った文字だけを消す） | 寄せる | accession の取得は NCBI への送信になる（`web/AGENTS.md` 規則 4）。FASTA だけを受ける |
| 「Query subrange」From / To | 「Query subrange」From / To と範囲の棒（`query-region`、`-query_loc`） | 寄せる | LOSAT は含めたレコードが 1 つのときだけ範囲を指定できる（DW-9）。それ以外は理由の一文 |
| 「Or, upload file」 | 「Or, upload file」：「Open FASTA files…」と落とす場所（`query-open`、`query-files`） | 寄せる | LOSAT は複数のファイルを受け、中身を貼る欄に写さない |
| 「Genetic code」（blastx・tblastx の query） | 「Genetic code」（`param-query_gencode`）を Enter Query Sequence の中に | 採る | W3 は Algorithm parameters の「Genetic codes」に置いていた |
| 「Job Title」と「Enter a descriptive title for your BLAST search」 | 「Job Title」（`job-title`）。Run の名前として RunSnapshot の `title` に持ち、argv には入れない。キューの Run と結果の見出しに示す | 採る | NCBI は最初の query の defline で自動で埋めるが、LOSAT は埋めない（空なら名前なし） |
| 「Align two or more sequences」の印 | なし | 採らない | LOSAT は常に subject を与える検索。NCBI で印を付けた形（Enter Subject Sequence を出し Choose Search Set を隠す）を常に取る |
| 枠「Enter Subject Sequence」（textarea、Clear、Subject subrange、Or, upload file） | 枠「Enter Subject Sequence」（`subject-panel`）。query と同じ形 | 寄せる | |
| （無い）subject の遺伝暗号 | 「Genetic code」（`param-db_gencode`、TBLASTN・TBLASTX）を Enter Subject Sequence の中に。承認済みの例外の注（`subject-gencode-note`）を保つ | LOSAT だけ | NCBI の 2 配列の画面に subject の遺伝暗号は無い |
| 枠「Choose Search Set」（Database、Organism、Exclude、Entrez Query） | なし | 採らない | データベースと Taxonomy は LOSAT に無い |
| LOSAT のレコード一覧と除外、検査の結果、「First lines」 | 各枠の中、ファイルの選択の下（`*-source-*`） | LOSAT だけ | |
| LOSAT の Combined / Separate | 各枠の中、レコード一覧の下（`*-mode`） | LOSAT だけ | |
| 枠「Program Selection」：BLASTN の「Optimize for」の radio | 「Optimize for」の radio（`param-task`、各 radio は `param-task-<task>`）：Highly similar sequences (megablast)、More dissimilar sequences (discontiguous megablast)、Somewhat similar sequences (blastn)、Short sequences (blastn-short) | 寄せる | blastn-short は LOSAT が足す 4 つ目（NCBI は「Short queries」の自動調整で代える）。既定は LOSAT の CLI と同じ megablast（NCBI の 2 配列の画面の既定は blastn）。出力を CLI と比べるため |
| BLASTP の「Algorithm」の radio（2 配列では blastp だけ） | 「Algorithm」の radio：エンジンの task（blastp (protein-protein BLAST)、Quick BLASTP (blastp-fast)、blastp-short） | 寄せる | task の一覧は `describe` の choices。名前の補いは表示だけ |
| tblastn・tblastx（選ぶものが無く、枠を出さない） | TBLASTN は「Algorithm」（tblastn、tblastn-fast）。TBLASTX は枠を出さない | 寄せる | LOSAT の TBLASTN には task が 2 つある |
| BLAST のボタンと「Search nucleotide sequence using Megablast (Optimize for highly similar sequences)」 | ボタン「Run LOSAT」（キューが空）／「Run LOSAT (add to queue)」（実行中か待ちがある）、複数の Run は「Run LOSAT (2 runs)」「Run LOSAT (2 runs, add to queue)」（`add-to-queue`）。横に「Search nucleotide subjects with BLASTN, task megablast (optimized for highly similar sequences). Runs in this browser.」の形の一文（`search-summary-line`） | 寄せる | ボタンの名前は保守者の指示（2026-10-10）で「Run LOSAT」（「BLAST」を操作の名前にしない）。キューに積むことはボタンの括弧、この機械で動くことは横の一文で示す（指示書 3.） |
| 「Show results in a new window」 | なし | 採らない | 結果は同じタブの「Results」に出る。新しい窓は保存と worker を分ける |
| （無い） | Threads の select、エンジンの option の検査の一文、キューの一文（「Show the queue」）、入力の不足の一文 | LOSAT だけ | ボタンの近く（ボタンの行の下） |
| （無い） | 「Search settings」の行（S15）：「Save settings」（`settings-save`、フォームの program・option・スレッドを `losat-settings-{program}.json` に。入力・名前・Job Title は書かない）、「Load settings…」（`settings-load`）、何を読んだか・入れなかったか・なぜ拒んだかの一文（`settings-message`、「Edit Search」もここに書く） | LOSAT だけ | Threads の行の下。LOSAT Web の形式であり、NCBI の「Save Search」（アカウントへの保存）ではない |
| 「Algorithm parameters」（開閉、閉じて始まる）と「Note: Parameter values that differ from the default are highlighted in yellow and marked with ♦ sign」 | 開閉する「Algorithm parameters」（`algorithm-parameters`、閉じて始める。既定と違う値があれば開いて始め、見出しに数）。同じ注 | 採る | |
| 「Restore default search parameters」 | 同じ言葉のボタン（`restore-defaults`）。Algorithm parameters の節の値を消す（task と遺伝暗号は残す） | 採る | |
| 既定と違う値の黄色と ♦ | argv に書く値（`formParameters` が書く値）の欄を黄色（LOSAT の注意の色）と ♦ で示し、読み上げに「changed from the default」 | 採る | LOSAT は既定と同じ値を argv に書かないので「書く値」がそのまま「既定と違う値」。task で変わる既定（`describe` に既定が無い）も、書く値なら印を付ける |
| 「General Parameters」：Max target sequences（select）、Short queries、Expect threshold、Word size（select）、Max matches in a query range | 「General Parameters」：Max target sequences、Expect threshold、Word size。TBLASTX の Max matches in a query range（`-culling_limit`） | 寄せる | select の選択肢は写さない（下の「select の選択肢」）。Short queries は採らない（LOSAT は CLI のとおりで自動の調整をしない。blastn-short を手で選ぶ）。「Max matches in a query range」は NCBI の `HSP_RANGE_MAX` で、その説明が culling の論文を指す（TBLASTX の `-culling_limit`） |
| 「Scoring Parameters」：Match/Mismatch Scores（select）、Gap Costs（select）、Matrix、Compositional adjustments | 「Scoring Parameters」：Match/Mismatch Scores（Match・Mismatch の 2 欄、`-reward`・`-penalty`）、Gap Costs（Existence・Extension の 2 欄、`-gapopen`・`-gapextend`）、Matrix、Compositional adjustments（0〜3 に NCBI の名前を添える） | 寄せる | 組の値を 1 行に。どの組が許されるかはエンジンの検査が決める |
| 「Filters and Masking」：Filter（Low complexity regions、Species-specific repeats）、Mask（Mask for lookup table only、Mask lower case letters） | 「Filters and Masking」：Low complexity regions filter（`-dust` / `-seg`）、Mask for lookup table only（TBLASTN の `-soft_masking`）、Mask lower case letters（`-lcase_masking`） | 寄せる | Species-specific repeats は採らない（NCBI のリピートのデータベース）。DUST・SEG は値（no、yes、引数）を書ける欄のまま |
| 「Discontiguous Word Options」：Template length、Template type | 同じ節（`-template_length`、`-template_type`）。task が dc-megablast のとき、または値があるときに出す | 採る | NCBI は discontiguous megablast のときだけ出す。値が隠れたまま argv に入らないよう、値があれば出す |
| 「PSI/PHI/DELTA BLAST」の節 | なし | 採らない | LOSAT に無い |
| （無い） | 「Other Parameters」：NCBI の画面に無い LOSAT の option（BLASTN の Max HSPs per subject・Percent identity・Subject best hit、BLASTP・TBLASTN・TBLASTX の Neighboring words threshold・Two-hit window size、TBLASTN の X-dropoff・Sum statistics） | LOSAT だけ | Algorithm parameters の最後の節 |
| 2 つ目の BLAST のボタン | Algorithm parameters を開いたとき、その下に 2 つ目の「Run LOSAT」（`add-to-queue-bottom`） | 採る | |
| フッター、ヘッダー、ロゴ | LOSAT の見出し（「LOSAT Web」と「Local processing…」） | 採らない | ブランド |

**select の選択肢**：NCBI の Word size・Match/Mismatch・Gap Costs・Max target sequences の select は、program と task ごとに許す値を画面に埋めている。LOSAT はこれを写さない。欄の名前と並びを NCBI に合わせ、値は自由な欄のまま（Matrix は候補つき）、エンジンの `validate` が判断する（W3 の判断 3・9、TS で検証の規則を作らない）。既定の値は `describe` の placeholder で示す。

## 2. 結果画面

NCBI の結果は、上から操作の行（Edit Search、Save Search、Search Summary）、見出しの塊（Job Title、RID、Program、Query ID…）と右の「Filter Results」、タブ（Descriptions、Graphic Summary、Alignments、（Taxonomy）、Dot Plot）、タブの中身、の順である。NCBI の classic（Traditional）結果ページ（新しいページの「Back to Traditional Results Page」で戻る形）は、タブの代わりに Graphic Summary、Descriptions、Alignments を 1 ページに上から並べる。LOSAT は操作の行・見出しの塊・Filter Results・「Results for」・知らせをこの順に置き、その下のタブの最初を classic の 1 ページ（2026-10-10 から。下の「タブ」）にし、Dot Plot、Run details、Outputs を続ける。

### 見出しの塊

| NCBI の要素 | LOSAT Web | 採否 | 理由・置き場所 |
|---|---|---|---|
| 「Edit Search」 | 「Edit Search」（`edit-search`、見出しの塊の上の行）：Run の設定（program、argv の `slice(5)`（範囲を含む）、スレッド）と Job Title を検索フォームに入れ、Search タブを設定の行で開き、「The search form has the settings of Run N. The inputs are the form's own.」と言う（`settings-message`）。入力はフォームのまま。検索は始めない | 寄せる | S15：設定ファイルの読み込みと同じ規則（`SearchDraft.applySettings`）。フォームに欄の無い option、入力に当たらない範囲（範囲は 1 レコードの役割だけ）、この browser に無いスレッド数は入れずに「Not applied: …」と並べる。Run の入力を戻すのは、セッションファイルのつなぎ直し（REQ-23）の役目で、ここではしない |
| 「Save Search」「How to read this report?」「BLAST Help Videos」「Back to Traditional Results Page」 | なし | 採らない | NCBI のアカウントと文書。classic のページへの切り替えは無い：2026-10-10 から最初のタブが classic の 1 ページ（「タブ」） |
| 「Search Summary」（Search Parameters・Karlin-Altschul statistics・Results Statistics） | タブ「Run details」（`results-view-details`）：RunSnapshot、RunRecord、検証バッジ、「Reproduce this run」（S15、`run-reproduce`：LOSAT のコマンド `run-command-<format>`、NCBI BLAST+ 2.17.0 の比較用のコマンド `run-ncbi-command-<format>` か比べられない理由 `run-ncbi-unavailable`、承認済みの例外、入力 FASTA の保存 `run-input-save-<query|subject>`、設定ファイル `run-settings-save`） | 寄せる | 統計の値（Lambda、K、H、有効探索空間）は outfmt 0 の末尾の原文にある（Outputs）。TS で抜き出さない |
| 「Job Title」 | 「Job Title」（Run の `title`。無ければ行を出さない） | 採る | |
| 「RID」と「Search expires on」 | 「Run」：Run の select（`results-run`、Run の番号・program・入力名） | 寄せる | RID に当たるのは Run の番号。期限は無い（作業セッションの中） |
| 「Download All」 | 「Download All」：タブ「Outputs」へ移る（`results-download-all`） | 寄せる | 原文の outfmt 0/6/7 と診断を書き出す（Outputs） |
| 「Program」「Blast 2 sequences」と「Citation」 | 「Program」：BLASTN（task megablast）など。横に検証バッジ（`verification-badge`） | 寄せる | Citation は採らない。検証バッジは LOSAT だけ（NCBI と比べた記録） |
| （無い） | 「Options」：argv の option（`results-run-options`） | LOSAT だけ | |
| 「Query ID」「Query Descr」「Query Length」 | 「Query ID」と「Query Length」（1 query の Run） | 寄せる | Query Descr は採らない（Run の記録が持つのは ID と長さ。outfmt 0 の `Query=` の行は Alignments と Outputs で読める） |
| 「Subject ID」「Subject Descr」「Subject Length」 | 「Subject ID」と「Subject Length」（subject が 1 レコードの Run）。複数なら「Subjects」：入力名とレコード数 | 寄せる | Subject Descr は Descriptions と Alignments に（outfmt 0 の見出し） |
| 「Other reports」（MSA viewer、Distance tree） | なし | 採らない | LOSAT に無い |
| 「Results for」（複数 query の select） | 「Results for」：Query の一覧（`query-list`、仮想化、ID の絞り込み `query-filter`、「With hits only」`filter-hits-only`）。複数 query の Run だけ。行は ID・長さ・subject と HSP の数。ID は 12ch 以上を保ち（共通の接頭辞を持つ ID が見分けられる）、800 px 以下では長さと数を ID の下の 2 行目に置く | 寄せる | 10 万 query を 1 つの select には入れられない。仮想化した一覧を同じ場所と名前で置く |
| 当たりの無い query：黄色の帯「No significant similarity found.」、タブが消える | 「No significant similarity found for this query.」（`results-notice` の `no-hits`）。タブは残す（Run details と Outputs は query に依らない） | 寄せる | outfmt 0 の「No hits found」と同じ事実。値は作らない |
| 「Filter Results」：Percent Identity、E value、Query Coverage の from / to、「Filter」「Reset」 | 「Filter Results」：E value ≤、Bit score ≥、Subject ID contains、「Filter」（`filter-apply`）「Reset」（`filter-clear`） | 寄せる | Percent Identity は採らない（レコードに原値が無い、`results_columns.md`）。Query Coverage はエンジン待ち（E2j）。表示だけを変え、検索をやり直さない（ViewState の規則は変えない） |

### タブ

2026-10-10 から（保守者の指示、S15 の DECISIONS 12）：「Descriptions」（`results-view-hits`。中身は classic の 1 ページ：Graphic Summary、Descriptions、Alignments。下の「1 ページの結果（classic）」）、「Dot Plot」（`pane-dotplot`）、その後に LOSAT の「Run details」（`results-view-details`）と「Outputs」（`results-view-outputs`）。最初のタブの名前は NCBI の最初のタブの名前「Descriptions」のままにした（1 ページの中心は Descriptions で、Graphic Summary はその概観、Alignments はその行を開いたもの。「Results」は主のタブと結果の見出しに重なる）。結果を開くとき、Run を替えるとき、候補の「Show in results」、Graphic Summary と Dot Plot の「Show alignment」は、どれもこの 1 ページに来る。「Taxonomy」は採らない（ローカルの FASTA に無い）。NCBI のタブの形（選んだタブを塗り、帯の下に道具の行）を LOSAT の色で写す。どのタブも同じ選択（中心は HSP の ID、`ResultsBrowser`）に従う。

履歴：S13b（2026-10-10）では NCBI の新しいタブ式の結果ページに寄せ、NCBI と同じ順で「Descriptions」（`results-view-hits`）、「Graphic Summary」（`results-view-graphic`）、「Alignments」（`pane-alignment`）、「Dot Plot」、「Run details」、「Outputs」の 6 つのタブにしていた。タブを開き直さないと概観と一覧と整列を見比べられないので、保守者の指示で上の形に替えた（`results-view-graphic` と `pane-alignment` は無くなった）。

### 1 ページの結果（classic）

NCBI の classic の結果ページのように、選んだ query について上から次の 3 つの節を並べる（`results-classic`）。節ごとに NCBI の見出しを置く。タブを開かずに、概観（ヒットの query 上の座標と得点の階級）と一覧が続けて見える。

| 節 | LOSAT Web | 切り方（NCBI のように、全部を一度に開かない） |
|---|---|---|
| 「Graphic Summary」（`results-graphic`） | 下の「Graphic Summary」の図。1 行 8 px（棒は 4 px）の詰めた帯で、30 行を一度に見せ、それより多いと図の中で行が流れる（S13b は 1 行 12 px、40 行） | 最初の 100 subject を描く（Descriptions の順と表示用のフィルターに従う）。「Show all N」で全部（W4b のまま）。描くのは図の中で見えている行だけ |
| 「Descriptions」（`results-descriptions`） | 下の「Descriptions」の表（見出し「Sequences producing significant alignments」、行の印、select all、Add to candidates） | 最初の 100 subject（今の並び順で。既定はエンジンの順）を載せ、「The first 100 of N are listed.」（`descriptions-listed`）と「Show all N」（`descriptions-show-all`）。NCBI の Descriptions は 1 ページに 100 行。「Show all」は Run か query を替えるまで続く（結果のタブを離れて戻っても続く。Graphic Summary の「Show all」も同じ）。載せた一覧は仮想化（見えている行だけ描く）。「select all」は載せた subject に印を付ける。選んだ subject が載せた範囲の外（並べ替えた後、Alignments の「Next」、候補の「Show in results」）でも、一覧は広げない（Alignments と Graphic Summary は従う） |
| 「Alignments」（`results-alignments`） | 下の「Alignments」：選んだ subject の塊 | 選んだ subject の塊だけ（全部の subject の整列を並べない）。その中の Range は選んだ HSP の前後 25 個ずつ、節の原文は見える所に来た Range だけ読む（W4b のまま） |

選び方とスクロール：Descriptions の行を選ぶと（クリック、Enter・Space）、その subject が選ばれ、Alignments の見出し（`results-alignments-heading`）を画面の上に出して focus を移す（NCBI の説明の行のリンクが整列へ飛ぶのと同じ。focus が見えている所にある）。Graphic Summary の棒のクリックと Enter は、その HSP の Range を画面に出して focus を移す。Alignments の「Descriptions」は Descriptions の見出しに戻り、選んだ行（描かれていれば）に focus を移す。ページが動くのは選んだときだけ：Graphic Summary の矢印キーは図の中で注目の HSP を動かすだけで選ばず（行は図の中で流れる）、Descriptions の行の間は Tab で動き選ばない。並べ替え、表示用のフィルター、query の切り替え、印ではページは動かない。矢印キーごとにページが飛ぶと一覧を読めないので、こうした。

速さ：開いたときに描くのは 3 つの節の見えている部分だけ（図は見えている行、一覧は見えている行、Alignments は選んだ subject の塊と、見える所に来た Range の原文）。100,000 query の Run を開く速さを W5 と比べた（S15 の WP-H、`results-measure.spec.ts`）。

### Descriptions

| NCBI の要素 | LOSAT Web | 採否 | 理由 |
|---|---|---|---|
| 見出し「Sequences producing significant alignments」 | 同じ見出し | 採る | |
| BLASTP の最初のタブ「Clusters」（「Clusters producing significant alignments」、Cluster Representative Sequence） | 「Descriptions」のまま | 採らない | NCBI のデータベースの配列のまとまり。ローカルの subject はまとめない |
| 「select all」と行の印（「Select for downloading or viewing reports」）、「N sequences selected」 | 「select all」（`descriptions-select-all`）、行の印（`subject-mark-<sIdx>`）、「N sequences selected」（`descriptions-selected`）。印はその query の subject に付き、query や Run を替えると外れる（表示用のフィルターが隠した subject の印も外れる） | 寄せる（S14） | 候補への追加に使う（NCBI は Download・GenBank・Graphics に使う）。行は選択の 1 つのボタンのまま（W4b）で、印はその左の別の check box |
| 選んだ行への操作「Download」「GenBank」「Graphics」「MSA Viewer」 | 「Add to candidates」（`descriptions-add-candidates`）：印の subject の、その query の全 HSP（表示用のフィルターに関わらず） | LOSAT だけ（S14） | 書き出しは候補トレイから（「候補（S14）」）。GenBank・Graphics・MSA Viewer は NCBI のデータベースと viewer のもの。outfmt 0/6/7 の原文の Download は Outputs |
| 「Select columns」「Show」（1 ページの行数） | 「Show」は「The first 100 of N are listed.」と「Show all N」（`descriptions-show-all`）に寄せる（2026-10-10）。「Select columns」はなし | 寄せる・採らない | 列は固定（`results_columns.md`）。NCBI の 1 ページの既定の 100 行で切り、残りは「Show all」で同じ仮想化した一覧に載せる（ページを分けない） |
| 列 Description（リンク） | 「Description」：outfmt 0 の見出しの題（`results_columns.md`）。行を選ぶとその subject が選ばれ、Alignments・Graphic Summary・Dot Plot が従い、同じページの Alignments が画面に来る（2026-10-10 から。S13b では「Alignments」タブで） | 寄せる | 行は 1 つのボタン（中にリンクを入れない） |
| 列 Scientific Name、Common Name、Taxid | なし | 採らない | ローカルの FASTA に無い |
| 列 Max Score | E2j の前：「Score (bits)」（その subject の最初の HSP の outfmt 6 の bitscore）をこの位置に | エンジン待ち | 列名は今の意味のまま（Max Score と名乗らない） |
| 列 Total Score、Query Cover、Per. Ident | なし | エンジン待ち（E2j） | |
| 列 E value | 「E value」（最初の HSP） | 寄せる | E2j の後は表形式の E value |
| （無い） | 「HSPs」：HSP の数（E value の後、Length の前） | LOSAT だけ | アプリの数え上げ。E2j で Score と E value の間に Total Score・Query Cover が入っても、この列は動かない |
| 列 Acc. Len | 「Length (nt / aa)」：Run のレコードの長さ | 寄せる | |
| 列 Accession | 「Subject ID」：outfmt 6 の `sseqid`（最後の列）。列の幅は値に合わせ（10ch から 24ch、長い ID は「…」で切り title に全文）、Description に残りの幅 | 寄せる | データベースの accession ではないので「Subject ID」と呼ぶ。幅は W4b の画面レビュー L8 |
| （無い） | 「#」：エンジンの順（最初の列、NCBI の印の位置） | LOSAT だけ | |
| 並べ替え（Accession 以外、既定は E value） | 「#」・Score・E value・HSPs・Length で並べ替え。既定はエンジンの順（`#`） | 寄せる | エンジンの順は outfmt 0 の説明表の順で、NCBI の既定の E value の順と同じ考え（同値の扱いがエンジンのまま）。Description の並べ替えは採らない（見出しは見えている行の分だけ読む） |
| 電話の幅 | 大事な列（Subject ID、Score、E value、HSPs）を先に、横に流れるときだけ知らせる（W4 の判断 17） | LOSAT の基準 | |

### Graphic Summary

| NCBI の要素 | LOSAT Web | 採否 | 理由 |
|---|---|---|---|
| 道具の行「hover to see the title」「click to show alignments」 | 同じ言葉 | 採る | |
| 「Alignment Scores」の凡例：< 40 黒、40 - 50 青、50 - 80 緑、80 - 200 マゼンタ、>= 200 赤 | 同じ階級と色（#000000、#0020e9、#75ea4c、#db3de9、#db3324） | 採る | 読み慣れた機能上の約束（ブランドではない）。HSP レコードの `bit_score`（原値）を区切るだけ |
| 「Show Conserved Domains」 | なし | 採らない | NCBI の CDD |
| 「Distribution of the top N Blast Hits on M subject sequences」 | 「Distribution of N HSPs on M subject sequences」（見せている subject の数） | 寄せる | |
| query の棒と目盛り | 同じ：棒「Query」と目盛り（1 と query の長さ、その間の切りのよい位置）、単位は query の nt / aa | 採る | 座標は Run の記録の長さ |
| subject ごとの行、HSP ごとの細い棒（同じ subject の HSP は 1 行、灰色の細い線でつなぐ） | 同じ。Descriptions の順と表示用のフィルターに従い、最初の 100 subject（「Show all M」で全部）。2026-10-10 から 1 行 8 px・棒 4 px・30 行を一度に（classic の詰めた帯。S13b は 12 px・6 px・40 行） | 採る | 棒の位置は HSP レコードの query の座標 |
| hover：題、Score、Evalue | hover と focus：Subject ID、Description（読めていれば）、outfmt 6 の bitscore と evalue | 寄せる | |
| click：その alignment へ | click（と Enter）：その HSP を選び、同じページの Alignments のその「Range」へ（focus も） | 採る | キーボード：上下で subject、左右で HSP（選ばず、ページも動かさない） |

### Alignments

| NCBI の要素 | LOSAT Web | 採否 | 理由 |
|---|---|---|---|
| 「Alignment view」（Pairwise、Query-anchored など）、「CDS feature」、「Restore defaults」、「Line length」 | なし | 採らない | エンジンが書くのは outfmt 0 の pairwise だけ。原文を書き換えない（`web/AGENTS.md` 規則 2） |
| subject ごとの塊（全部の subject を続けて） | 選んだ subject の塊（`alignments-subject`）。「Previous」「Next」で前後の subject、「Descriptions」で同じページの一覧へ戻る（2026-10-10 から。S13b ではタブを替えた） | 寄せる | 数千の subject を 1 ページに並べず、選択に従う（NCBI も整列を既定の数で切る） |
| 「Download」（subject の塊の道具の行） | 「Add all matches to candidates」（`alignments-add-subject`）：その subject の、その query の全 HSP | LOSAT だけ（S14） | 書き出しは候補トレイから。原文の Download は Outputs |
| 「Range n: a to b」の横の「GenBank」「Graphics」 | 「Add to candidates」（`range-add-<q>-<r>`）。候補にあれば「In candidates」 | LOSAT だけ（S14） | NCBI のリンクの場所に置く |
| 「Sort by」 | なし | 採らない | Range はエンジンの順。並べ替えは HSP の表で |
| 題（subject の説明）、「Sequence ID」「Length」「Number of Matches」 | 題（outfmt 0 の見出しの題）、「Sequence ID」（sseqid）、「Length」（レコードの長さ）、「Number of Matches」（HSP の数） | 採る | |
| （無い） | その subject の HSP の表（`hsp-table`、outfmt 6 の欄、並べ替え）。frame は翻訳した配列のものだけ：TBLASTN は「Subject frame」（`+2`）、BLASTX は「Query frame」、TBLASTX は「Frames (q/s)」（`-2/+2`） | LOSAT だけ | 塊の見出しの下。行を選ぶとその Range へ。TBLASTN の「–/+2」は TBLASTX の「-2/+2」と並ぶとマイナス鎖に読めた（W4b の画面レビュー L9）。NCBI の節も「Frame = +2」 |
| 「Range n: a to b」と「Next Match」「Previous Match」「First Match」 | 「Range n: a to b」（n はその subject の中のエンジンの順、a と b は HSP レコードの subject の座標を小さい方から）、「Next Match」「Previous Match」「First Match」 | 採る | NCBI の Range は subject の座標（参照画面の blastn：Range 2: 588 to 608 は Sbjct 588–608） |
| Score・Expect・Identities・Gaps・Strand（Frame、Method、Positives）の表 | outfmt 0 の節を原文のまま（その中に `Score =`・`Expect =`・`Identities =`・`Gaps =`・`Strand=` / `Frame =` の行）（`detail-section`） | 寄せる | 原文を表に組み直さない |
| Query / Sbjct の等幅の段 | 同じ（outfmt 0 の原文） | 採る | |
| 「Related Information」 | なし | 採らない | NCBI のデータベースへのリンク |
| （無い） | 選んだ HSP の塊に outfmt 6 の行（`detail-row`）、outfmt 0 に無い HSP の理由（`detail-not-in-outfmt0`）、1 文字の HSP の鎖の注（`detail-strand-note`） | LOSAT だけ | W4 の詳細（`hsp-detail`）を選んだ Range の塊に移す |

多数の HSP（例：5,993 HSP の 1 組）では、塊の Range を選んだ HSP の前後 25 個ずつだけ示し、「Show earlier matches」「Show later matches」で広げる。節の原文は見える所に来た Range だけを読む。

### Dot Plot

NCBI の Dot Plot は、題「Plot of <query> vs <subject>」とサーバーが描いた 600×300 の画像（薄い緑の地、灰色の格子、HSP はすべて濃い灰色の線で鎖を色で分けない、題は画像の外）である。タブの名前と位置を NCBI に寄せ、図の描き方は保守者の [`blast2dotplot.py`](https://github.com/satoshikawato/bio_small_scripts/blob/main/blast2dotplot.py)（commit `3b55d116`、2023-12-25、sha256 `4f0731f769eaa755e940f82f373a9bb86c8fa311df030c8810d272ed66f23f84`）を出発点にする。NCBI と違うところはスクリプトを先にする（保守者の指示、2026-10-09）。

| 要素 | スクリプト | LOSAT Web | 理由 |
|---|---|---|---|
| 題 | query の題を上、subject の題を左に 90° | 上に「Plot of <query ID> vs <subject ID>」（NCBI の言葉）、軸の名前は上に「Query <ID> (<単位>)」、左に 90° 回して「Subject <ID> (<単位>)」。単位は見えている長さの bp・kbp・Mbp（aa の軸は aa・kaa・Maa）。aa の軸と nt の軸の組で縮尺どおりのとき、aa の軸は「(aa; drawn at 3 nt per aa)」。場所が足りないときは ID を「…」で切り、単位を残す | 単位はスクリプトの bp・kbp・Mbp（保守者のスクリプトを先にする）。説明の行・ポップアップ・「Selected:」の行は Run の単位（nt / aa）のまま |
| 軸 | X が query、Y が subject。原点は左上、subject は下へ増える。目盛りとラベルは上と左 | 同じ | スクリプトどおり（W4 は原点が左下だった） |
| 縮尺 | 両軸同じ縮尺。長い方を 1000 px | 同じ縮尺。長い方を描く枠の幅（最大 1000 px）に合わせる。aa の軸と nt の軸の組（TBLASTN の query、BLASTX の subject が aa）では、縦横の比を決めるときだけ 1 aa を 3 nt と数える（目盛りとラベルは各軸の単位のまま。S13b の判断 26）。短い辺は最小 120 px（それより細いときは縮尺を変え、図の下に「Axes not to scale」と示す） | 長さの比が大きいと短い辺が 1 px 以下になる。1 aa を 1 文字と数えると HSP の線の傾きが 3 になり、aa の辺が 3 分の 1 に縮む |
| 目盛り | `tick_size` の大小の刻み（1 kbp 以下 100 / 10 … 5 Mbp より上 1 M / 500 k）。ラベルは大きい刻みだけ、5000 未満 bp・1 Mbp 未満 kbp・それ以上 Mbp | 同じ表。ズームしたときは見えている長さで刻みと単位を選び直す。単位はラベルごとに書かず軸の名前に置く（「Query q0 (kbp)」。390 px で場所が足りない）。重なるラベルは間引く。aa の軸は単位を aa・kaa・Maa にした同じ刻み。小さい刻みは間隔が 5 px 未満のときだけ省く（両軸に同じ規則。数が 50 を超えると省く規則では、縮尺が同じでも片方の軸だけ消えた。W4b の画面レビュー M2） | スクリプトは BLASTN・TBLASTX だけ |
| 格子と枠 | 大小の刻みに薄い灰色 `#D3D3D3`、周りに黒い枠 | 同じ | |
| HSP の線 | 幅 2。両方の frame の符号が同じなら青 `#1f77b4`、違えばオレンジ `#ff7f0e` | 同じ。符号は HSP レコードの座標と frame から（W4 の向き）。BLASTN の 1 文字の HSP（鎖がレコードに無い）は灰色 `#7f7f7f` と凡例（E2j まで） | |
| 不透明度 | identity（`int()` で切り捨て）60 以下 0.4、70 以下 0.6、80 以下 0.8、それより上 1 | outfmt 6 の `pident` の文字列を数として読み、整数に切り捨てて同じ階級に分ける（TS で identities / 長さを計算しない） | 丸めた `pident`（小数 3 桁）の切り捨てで、スクリプトの切り捨てと階級の境界が同じになる |
| 出力 | SVG のファイル | 画面（Canvas）と、図の脇の「Download SVG」（`dotplot-svg`）で書き出す SVG のファイル（S15）。今の拡大・移動の view と、view の絞り込みが残す HSP を、画面と同じ軸・目盛り・ラベル・格子・色・不透明度の階級・線の規則（`domain/plot-layout.ts` と `plot-scale.ts` を Canvas と共有）で、図の CSS px の大きさに書く。軸の題（query と subject の ID と単位）、`<title>` と `<desc>`（LOSAT Web の run N の dot plot で、LOSAT Web の形式であり NCBI の図ではないこと）を付ける。選択した HSP の halo と hover は画面の補助なので書かない。file 名は `losat-run{N}-dotplot-q{Q}-s{S}.svg`（Q と S は 1 始まりのレコードの位置。ID や入力の名前は使わない）、`image/svg+xml`。HSP は 1 本の `<line>` ずつ Writer のブロックに書く（1 つの文字列に組まない）。script、event 属性、外部参照（`href`、`url(`）は無く、枠の外は入れ子の `<svg>` で切る。文字はすべて escape し、XML に書けない文字は U+FFFD にする | 出力の形式と名前は S15 が決めた。書く値は HSP レコードの座標と outfmt 6 の `pident` の階級だけで、配列から新しく計算しない |
| HSP のポップアップの操作 | （無い。NCBI の Dot Plot にも無い） | 「Show alignment」の横に「Add to candidates」（`dotplot-popup-add`）。候補にあれば「In candidates」 | LOSAT だけ（S14） |

LOSAT が足すもの（W4 の判断 18 を保つ）：ズーム（ボタン、+ / −、Ctrl / ⌘ とホイール）、パン（ドラッグ、矢印キー）、HSP の選択（線のクリック、n / p、表）、「Zoom to HSP」「Whole sequences」、`touch-action: pan-y`。新しく、hover で線を太くし、HSP を選ぶ（クリック、Enter）とその HSP のポップアップ（`dotplot-popup`）を出す：outfmt 6 の行の値をそのまま（bit score、E value、identity、query と subject の範囲、翻訳した配列の frame、向き、outfmt 0 に有るか）と「Show alignment」（Alignments のその Range へ）。ポップアップは線の端の脇に置き、線の中点を覆わない（そこに場所が無いときと、描く枠が 600 px 未満の電話の幅では図の下に）。ポップアップは Escape と「Close」で閉じ、キーボードとスクリーンリーダーで届く（focus を移し、閉じたら図に戻す）。描き方は層を分ける：格子と HSP の層（色と不透明度ごとにまとめ、線は 1 本ずつ描き、見えない線（表示の外の線、両端が前の線と同じ画素に来る不透明な線）を省く。重なる多数の線を 1 本の path にすると 2〜3 倍遅かった（W4b の B1）。薄い線の重なりはスクリプトの SVG の線と同じく濃くなる）と、hover・選択の層。ズームとパンは 1 フレームに 1 回だけ描き直す。

### 候補（S14）

候補トレイは LOSAT だけの主のタブ「Candidates」（§3）にあり、結果画面の Descriptions・Alignments・Dot Plot の「Add to candidates」で、完了した Run の HSP を集める（`REQ-10`）。抽出と書き出しは、トレイで選んだ候補から作る（設計書 §11.4）。

| NCBI の要素 | LOSAT Web | 採否 | 理由 |
|---|---|---|---|
| Descriptions の「Download」の「FASTA (complete sequence)」 | Region「Complete sequence」の「Download FASTA」（`extract-download`、`losat-candidates.fa`） | 寄せる（言葉） | 原配列は利用者の File から読む（NCBI はデータベースから）。NCBI の Download のメニューは参照画面に撮っていない：メニューの項目は NCBI の言葉として知られているもので、撮った画面と照らしていない |
| Descriptions の「Download」の「FASTA (aligned sequences)」 | 「Download aligned sequences (FASTA)」（`extract-aligned`、`losat-candidates-aligned.fa`）：HSP レコードの整列文字列を検索で得た向きのまま。原配列の抽出と別の file | 寄せる（言葉） | 同上（参照画面に無い） |
| （無い） | Region「Hit region」「Hit region with flanks」（Left・Right、レコードの単位）・「Complete sequence」、「Several HSPs on one record」（Separate sequences・One region spanning them）、Sequence（Subject・Query）（`extract-form`）。書き出しの後の要約（`extract-summary`：配列の数、端で切った配列の要求した範囲と実際の範囲、鎖が分からない HSP とその理由） | LOSAT だけ | 設計書 §11.4 |
| （無い） | 候補の表（`candidate-list`、行 `candidate-<n>`）、メモ（`candidate-note-<n>`）、並べ替え（Order added・Run・Subject）、上下の移動、削除（Remove・Remove selected）、「Show in results」（`candidate-reveal-<n>`、元の Run・query・HSP へ）、由来の一覧（Origins、`candidate-origins`） | LOSAT だけ | REQ-14 |

## 3. LOSAT だけのものの置き場所

| LOSAT の要素 | 置き場所 |
|---|---|
| キュー（Run の状態、取消、「Open results」） | 右の補助の領域（どのタブでも同じ幅と位置）。カードはどの状態も同じ形：1 行目 Run の番号・program・状態、2 行目 段階と時間（終わった Run は「Took」）と、同じ行の右端に操作（Cancel・Open results。入らなければその下の右）、グループの Run には 3 行目の右に「Cancel the group」、その下に題・入力・option |
| 離れることへの備え（wake lock、戻ったときの確認）、保存の状態 | 右の補助の領域、キューの下 |
| Combined / Separate、レコード一覧と除外 | 検索画面の各枠の中 |
| Subject の保持（R1） | キューの Run の詳細 |
| 検証バッジ | 結果の見出しの塊の Program の横、Run details |
| Outputs（outfmt 0/6/7 の原文、診断、書き出し） | タブ「Outputs」と見出しの「Download All」 |
| FakeEngine の帯 | 画面の上（W0） |
| 候補トレイ（一覧、メモ、並べ替え、削除、元の結果へ、由来、抽出と書き出し） | 主のタブ「Candidates」（`tab-candidates`、Search・Results の後、候補の数を示す）。どの Run を見ていても同じ |

## 4. 見た目

- 節：検索画面の枠は NCBI と同じく、薄い灰色の面に見出しのつまみ（legend）を載せる。見出しは LOSAT の青（`--accent`）の太字。
- 「Algorithm parameters」は幅いっぱいの帯のボタン（`+` / `−` の印、`aria-expanded`）。
- 結果のタブは帯のボタンで、選んだタブを LOSAT の青で塗り白の字。タブの下に薄い色の道具の行。1 ページの節の見出し（Graphic Summary、Descriptions、Alignments）は LOSAT の青の太字に細い下線（2026-10-10）。
- 表：0.875em、見出しは薄い青の面に太字（2 行まで折り返し、数の見出しは右寄せ）、行の境に細い線。リンクは青の下線。
- 色は LOSAT の配色のまま（ブランドを写さない）。Graphic Summary の凡例とドットプロットの線の色だけは上の約束。
- S13 の基準を保つ：1280 px で表が切れない、390 px で横に流れない、触れる対象は 24 px 以上、コントラスト、キーボード操作（W4 の判断 16〜23）。
