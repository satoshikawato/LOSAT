# NCBI BLAST の画面と LOSAT Web の対応表

段階 W4b（Session S13b、[指示書](../losat_web_gui_sessions/session_s13b_w4b_ncbi_style_ui.md)）で、保守者の指示（2026-10-09）「NCBI BLAST ウェブサイトの GUI にデザインを寄せてほしい」を具体にする。設計書 §2.1 の「NCBI Web BLAST に馴染みのある操作」と §11.2 の「NCBI Web に準じる入力項目・列を固定した参照画面と対応付ける」の参照画面と対応表である。

- 作成：2026-10-10（S13b）。採否は保守者の確認待ち（推奨案で進めた。Owner-delegated、[W4b のゲート記録](../evidence/losat_web_w4b/README.md)の判断）。
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
| program のタブ「blastn \| blastp \| blastx \| tblastn \| tblastx」 | 同じ並びのタブ。表示名は BLASTN・BLASTP・BLASTX・TBLASTN・TBLASTX。LOSAT ではこれを mode と呼ぶ | 寄せる | 表示名は保守者の決定（DW-21）。呼び方の mode も保守者の決定（DW-32）。BLASTX は SX まで「not available」（W3 の判断 2） |
| 画面の題「Align Sequences Nucleotide BLAST」など | なし | 採らない | NCBI の製品名の見出しで、ブランドに近い。program のタブと一文で足りる |
| program の一文（「BLASTN programs search nucleotide subjects using a nucleotide query.」） | タブの下の一文（「BLASTN mode searches nucleotide subjects using a nucleotide query.」など） | 寄せる | NCBI の 2 配列の画面の一文の形。「programs」は「mode」に（DW-32） |
| （無い） | 「Example」の行（`example-line`、線 WF-2）：mode の組ごとのボタン（BLASTN「16S rRNA (bacteria)」「ITS (fungi)」、BLASTP・BLASTX「E. coli K-12 proteins (UniProt reviewed)」「PDB proteins (pdbaa)」、TBLASTN・TBLASTX「E. coli K-12 genome」）と、その場の一文「Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov」（保存後は「Saved in this browser (2026-10-27, 41.8 MB)」と「Save to folder…」（Chromium / Edge、他は「Download」、線 WF-3）と「Delete」）、「?」 | LOSAT だけ | 押すと query の欄に表の query、subject に組の写しだけを置く。取得は押したときだけ、ファイルの URL だけを送る（`fetch_design.md` §9）。写しを端末のファイルにするのは要求を出さない（§10） |
| 「Reset page」「Bookmark」 | なし | 採らない | Reset page は開いた大きなファイルの選択も消す（入力は Remove で 1 つずつ外す）。Bookmark は URL に検索の状態を持たない（研究データを URL に入れない、`web/AGENTS.md` 規則 4） |
| 枠「Enter Query Sequence」 | 枠「Enter Query Sequence」（`query-panel`） | 採る | |
| 「Enter accession number(s), gi(s), or FASTA sequence(s)」と textarea、「Clear」 | 「Enter accession number(s) or FASTA sequence(s)」と textarea（`query-input`）、「Clear」（貼った文字だけを消す）。欄が accession だけ（1 行 1 つ）なら下に「Fetch」「from」と選択（Fastest / NCBI / EBI / CNCB、`query-fetch-source`。既定の Fastest は線 WF-3）、「Sends the accession to NCBI」、「?」の行（`query-fetch-row`、線 WF-1）。押すと取得した FASTA を欄に写し、普通の貼り付けとして読む。gi は受けない。Run は取得せず「Fetch the query accession first.」 | 寄せる | 2026-11-03 まで「accession の取得は NCBI への送信になる。FASTA だけを受ける」（`web/AGENTS.md` 規則 4）。保守者の決定（2026-11-03、`fetch_design.md` §7）で、明示に頼んだ取得は accession と URL だけを送るものとして許す |
| 「Query subrange」From / To | 「Query subrange」From / To と範囲の棒（`query-region`、`-query_loc`） | 寄せる | LOSAT は含めたレコードが 1 つのときだけ範囲を指定できる（DW-9）。それ以外は理由の一文。BLASTX も同じ（query は nt、subject は aa。DW-11・DW-39） |
| 「Or, upload file」 | 「Or, upload file」：「Open FASTA files…」と落とす場所（`query-open`、`query-files`） | 寄せる | LOSAT は複数のファイルを受け、中身を貼る欄に写さない |
| 「Genetic code」（blastx・tblastx の query） | 「Genetic code」（`param-query_gencode`）を Enter Query Sequence の中に | 採る | W3 は Algorithm parameters の「Genetic codes」に置いていた |
| 「Job Title」と「Enter a descriptive title for your BLAST search」 | 「Job Title」（`job-title`）。Run の名前として RunSnapshot の `title` に持ち、argv には入れない。キューの Run と結果の見出しに示す | 採る | NCBI は最初の query の defline で自動で埋めるが、LOSAT は埋めない（空なら名前なし） |
| 「Align two or more sequences」の印 | なし | 採らない | LOSAT は常に subject を与える検索。NCBI で印を付けた形（Enter Subject Sequence を出し Choose Search Set を隠す）を常に取る |
| 枠「Enter Subject Sequence」（textarea、Clear、Subject subrange、Or, upload file） | 枠「Enter Subject Sequence」（`subject-panel`）。query と同じ形 | 寄せる | |
| （無い）subject の遺伝暗号 | 「Genetic code」（`param-db_gencode`、TBLASTN・TBLASTX）を Enter Subject Sequence の中に。承認済みの例外の注（`subject-gencode-note`）を保つ | LOSAT だけ | NCBI の 2 配列の画面に subject の遺伝暗号は無い |
| 枠「Choose Search Set」（Database、Organism、Exclude、Entrez Query） | なし | 採らない | データベースと Taxonomy は LOSAT に無い |
| LOSAT のレコード一覧と除外、検査の結果、「First lines」 | 各枠の中、ファイルの選択の下（`*-source-*`）。「Records」と「First lines」は閉じて出す（NCBI はファイル名だけ。保守者 2026-11-02「レコードがずらっと出てくるのいらない」、DW-34） | LOSAT だけ | 以前は 2〜20 レコードのとき「Records」が開いていた |
| LOSAT の Combined / Separate | 各枠の中、レコード一覧の下（`*-mode`） | LOSAT だけ | |
| 枠「Program Selection」：BLASTN の「Optimize for」の radio | 枠「Task Selection」（DW-32）の「Optimize for」の radio（`param-task`、各 radio は `param-task-<task>`）：Highly similar sequences (megablast)、More dissimilar sequences (discontiguous megablast)、Somewhat similar sequences (blastn)、Short sequences (blastn-short) | 寄せる | blastn-short は LOSAT が足す 4 つ目（NCBI は「Short queries」の自動調整で代える）。既定は LOSAT の CLI と同じ megablast（NCBI の 2 配列の画面の既定は blastn）。出力を CLI と比べるため |
| BLASTP の「Algorithm」の radio（2 配列では blastp だけ） | 「Algorithm」の radio：エンジンの task（blastp (protein-protein BLAST)、Quick BLASTP (blastp-fast)、blastp-short） | 寄せる | task の一覧は `describe` の choices。名前の補いは表示だけ |
| tblastn・tblastx（選ぶものが無く、枠を出さない）。blastx の画面は未取得 | BLASTX は「Algorithm」（blastx、blastx-fast）、TBLASTN は「Algorithm」（tblastn、tblastn-fast）。TBLASTX は枠を出さない | 寄せる | LOSAT の BLASTX と TBLASTN には task が 2 つある |
| BLAST のボタンと「Search nucleotide sequence using Megablast (Optimize for highly similar sequences)」 | ボタン「Run LOSAT」（キューが空）／「Run LOSAT (add to queue)」（実行中か待ちがある）、複数の Run は「Run LOSAT (2 runs)」「Run LOSAT (2 runs, add to queue)」（`add-to-queue`）。横に「Search nucleotide subjects with BLASTN, task megablast (optimized for highly similar sequences). Runs in this browser.」の形の一文（`search-summary-line`） | 寄せる | ボタンの名前は保守者の指示（2026-10-10）で「Run LOSAT」（「BLAST」を操作の名前にしない）。キューに積むことはボタンの括弧、この機械で動くことは横の一文で示す（指示書 3.）。一文のうち何を検索するか（「nucleotide subjects」）と program・task（「BLASTN, task megablast (…)」）は、NCBI の「database core_nt」「Megablast (…)」と同じく色付きの太字（line WR、2026-11-02） |
| 「Show results in a new window」 | なし | 採らない | 結果は同じタブの「Results」に出る。新しい窓は保存と worker を分ける |
| （無い） | Threads の select、エンジンの option の検査の一文、キューの一文（「Show the queue」）、入力の不足の一文 | LOSAT だけ | ボタンの近く（ボタンの行の下） |
| 「Algorithm parameters」（開閉、閉じて始まる）と「Note: Parameter values that differ from the default are highlighted in yellow and marked with ♦ sign」 | 開閉する「Algorithm parameters」（`algorithm-parameters`、閉じて始める。既定と違う値があれば開いて始め、見出しに数）。同じ注 | 採る | |
| 「Restore default search parameters」 | 同じ言葉のボタン（`restore-defaults`）。Algorithm parameters の節の値を消す（task と遺伝暗号は残す） | 採る | |
| 既定と違う値の黄色と ♦ | argv に書く値（`formParameters` が書く値）の欄を黄色（LOSAT の注意の色）と ♦ で示し、読み上げに「changed from the default」 | 採る | LOSAT は既定と同じ値を argv に書かないので「書く値」がそのまま「既定と違う値」。既定は task と他の option で決まる値（下の「既定値」）。hit list の 3 つ（`-max_target_seqs`・`-num_descriptions`・`-num_alignments`）は BLAST+ が「与えたか」を区別するので、既定と同じ数でも打てば書き、印を付ける |
| 各欄の値（NCBI はどの欄にも効く値を入れて見せる） | 欄が空のあいだ、検索が使う既定値を欄の値として示す（下の「既定値」）。select は既定の項目を選んだ形、true / false の option は印の付いた checkbox。既定値を選び直すと欄は空に戻り、task を変えるとその task の既定値に変わる | 寄せる | 2026-11-02 までは clap の既定のある option だけ placeholder に出し、ほかは空の欄だった（保守者「NCBI と比べるとやっぱ見づらい」） |
| 「General Parameters」：Max target sequences（select）と横の一文「Select the maximum number of aligned sequences to display」、Short queries、Expect threshold、Word size（select）、Max matches in a query range | 「General Parameters」：Max target sequences（NCBI の select 10・50・100・250・500・1000・5000 と「Other」、横に NCBI の同じ一文）、Expect threshold、Word size（program と task ごとの NCBI の select と「Other」）、Max matches in a query range（`-culling_limit`、5 つの mode すべて、DW-42） | 寄せる | 一覧は下の「select の選択肢」。Short queries は採らない（LOSAT は CLI のとおりで自動の調整をしない。blastn-short を手で選ぶ）。「Max matches in a query range」は NCBI の `HSP_RANGE_MAX` で、その説明が culling の論文を指す（BLAST+ の `-culling_limit`）。2026-11-03 までは BLASTX と TBLASTX だけだった |
| 「Scoring Parameters」：Match/Mismatch Scores（select）、Gap Costs（select）、Matrix、Compositional adjustments | 「Scoring Parameters」：Match/Mismatch Scores（NCBI の組の select `param-matchscores`、「1,-2」が `-reward 1 -penalty -2`）、Gap Costs（組の select `param-gapcosts`、「Existence: 11 Extension: 1」が `-gapopen 11 -gapextend 1`、megablast の「Linear」が 0 0）、Matrix（NCBI の select の 8 つと「Other」）、Compositional adjustments（0〜3 に NCBI の名前を添え、既定の 2 を選んだ形） | 寄せる | 組の select の「Other」と、一覧に無い組は 2 つの欄（Match・Mismatch、Existence・Extension）を出す。Gap Costs の一覧は BLASTN では組ごと、タンパク質の mode では matrix ごと（下の「組の select」）。どの組が許されるかは今までどおりエンジンが決める（DW-40） |
| 「Filters and Masking」：Filter（Low complexity regions、Species-specific repeats）、Mask（Mask for lookup table only、Mask lower case letters） | 「Filters and Masking」：行「Filter」に NCBI と同じ checkbox の Low complexity regions（`-dust` / `-seg`。印がある間は横に引数の欄）、行「Mask」に NCBI と同じ 2 つの checkbox：Mask for lookup table only（`-soft_masking`、BLASTN は印あり）、Mask lower case letters（`-lcase_masking`） | 寄せる | Species-specific repeats は採らない（NCBI のリピートのデータベース）。Low complexity regions の印は検索が使う値が `no` でないこと（既定の印は task の既定に従う：BLASTN は印あり・blastn-short は無し、BLASTP は無し、ほかは有り）。印を外すと `no`、付けると既定の引数（既定が `no` の task では `yes`）。横の欄は `-dust` / `-seg` の値で、CLI の引数（DUST の LEVEL WINDOW LINKER、SEG の WINDOW LOCUT HICUT）を書ける（DW-40） |
| 「Discontiguous Word Options」：Template length、Template type | 同じ節（`-template_length`、`-template_type`）。task が dc-megablast のとき、または値があるときに出す | 採る | NCBI は discontiguous megablast のときだけ出す。値が隠れたまま argv に入らないよう、値があれば出す |
| 「PSI/PHI/DELTA BLAST」の節 | なし | 採らない | LOSAT に無い |
| （無い） | 「Other Parameters」：NCBI の画面に無い LOSAT の option（4 program の Number of descriptions・Number of alignments（`-num_descriptions` / `-num_alignments`。hit list の大きさが outfmt 0 と outfmt 6 / 7 で違うと検索を 2 回行う、`docs/web/abi_v2.md` §7）、BLASTN の Max HSPs per subject・Percent identity・Subject best hit、BLASTP・TBLASTN・TBLASTX の Neighboring words threshold・Two-hit window size、TBLASTN の X-dropoff・Sum statistics）。最後の行に「Additional options」（`param-extra`）：フォームに無い CLI の option を CLI のとおりに書く 1 行（DW-42） | LOSAT だけ | Algorithm parameters の最後の節。Additional options の語は shell の規則（`'…'`・`"…"`・`\`）で分け、フォームの option の後に argv へそのまま足す。読むのはエンジンの `validate`（CLI と同じ解析と同じ文）。アプリが拒むのは、LOSAT が持つ option（`-query`・`-subject`・`-out`・`-outfmt`・`-num_threads`）、フォームにある option と範囲（`-query_loc`・`-subject_loc`。「-evalue is set in the form: Expect threshold」：フォームは既定と同じ値を書かず、その既定をフォームの task や matrix から決めるので、ここの `-task` は省いた値を黙って変える）、検証と実行で読み方の違う語（`''`、`--`、NUL）と閉じていない引用符。語があれば黄色と ♦、「(N changed)」に 1 つと数え、「Restore default search parameters」で空にする。Run の option の行と Run details の argv は語を CLI のとおりに引用して示す（そのまま打ち直せる） |
| 2 つ目の BLAST のボタン | Algorithm parameters を開いたとき、その下に 2 つ目の「Run LOSAT」（`run-losat-bottom`） | 採る | |
| BLAST のボタンを押すと結果のページへ移り、終わるまで状態（Searching…）を示す | 「Run LOSAT」（`run-losat`）で Results のタブへ移り、始めた Run の結果のページを開く。終わるまで状態の行（`run-status`：Waiting / Preparing / Searching / Organizing results、経過時間、Cancel）を示し、終わると結果が出る。保守者の指示（2026-11-02）「検索開始したら結果に直接飛ぶようにして。NCBI BLAST みたいに」 | 採る | 2026-11-02 まではキューに加えるだけで、「Open results」で開いていた |
| フッター、ヘッダー、ロゴ | LOSAT の見出し（「LOSAT」とリポジトリの一文「A miniaturized reimplementation of NCBI BLAST algorithm」。2026-11-03 まで「LOSAT Web」と鍵の「Nothing leaves this browser」、DW-41） | 採らない | ブランド |

**select の選択肢**（line WR、2026-11-02）：NCBI の画面（2026-11-02 に取得した blastn・blastp・tblastn・tblastx の 2 配列の画面）の select の一覧を写し、最後に「Other」を置く。「Other」は横に文字の欄を開き、CLI が受ける他の値を書ける（どの値が通るかは今までどおりエンジンの `validate` が決める。TS で検証の規則を作らない、W3 の判断 3・9）。

| 欄 | 一覧（`domain/programs.ts` の `list`） |
|---|---|
| Max target sequences | 10、50、100、250、500、1000、5000（NCBI の既定は 100、LOSAT は CLI の 500） |
| Word size | megablast 16、20、24、28、32、48、64／dc-megablast 11、12／blastn と blastn-short 7、11、15／BLASTP・TBLASTN 2、3、5、6／TBLASTX 2、3／BLASTX 3、5（LOSAT の BLASTX が読む値だけ、DW-36） |
| Matrix | PAM30、PAM70、BLOSUM90、BLOSUM80、BLOSUM62、BLOSUM50、BLOSUM45、PAM250（NCBI の順）／BLASTX は BLOSUM62 だけ（LOSAT の BLASTX が読む行列。Gap Costs も 11 / 1 だけ、Compositional adjustments は 0 と 2 だけ、DW-36） |

**組の select**（line WR の SURVEY #13、DW-40）：Match/Mismatch Scores と Gap Costs は NCBI の画面の組の select（`domain/programs.ts` の `PairSelect`、節の `pairs`）。値は 2 つの欄（`-reward`・`-penalty`、`-gapopen`・`-gapextend`）に入り、argv の作り方は変わらない（既定と同じ値は書かない）。

| 欄 | 一覧（NCBI の画面、2026-11-02 に取得） | 既定（NCBI の「Deflt」） |
|---|---|---|
| Match/Mismatch Scores（BLASTN） | 1,-2、1,-3、1,-4、2,-3、4,-5、1,-1 | task の既定（megablast 1,-2、ほか 2,-3、blastn-short 1,-3） |
| Gap Costs（BLASTN、組ごと） | 1,-2：Linear、5 2、2 2、1 2、0 2、3 1、2 1、1 1／1,-3：Linear、5 2、2 2、1 2、0 2、2 1、1 1／1,-4：Linear、5 2、1 2、0 2、2 1、1 1／2,-3：Linear、4 4、2 4、0 4、3 3、6 2、5 2、4 2、2 2／4,-5：Linear、12 8、6 5、5 5、4 5、3 5／1,-1：5 2、3 2、2 2、1 2、0 2、4 1、3 1、2 1 | Linear（0 0）は megablast の時だけ出し、その時の既定。ほかは 5 2、4,-5 は 12 8 |
| Gap Costs（BLASTP・TBLASTN、matrix ごと） | BLOSUM62：11 2、10 2、9 2、8 2、7 2、6 2、13 1、12 1、11 1、10 1、9 1 ほか（PAM30、PAM70、BLOSUM80、BLOSUM45、BLOSUM50、PAM250、BLOSUM90 の一覧は `programs.ts`） | BLOSUM62 11 1、PAM30 9 1、PAM70 10 1、BLOSUM80 10 1、BLOSUM45 14 2、BLOSUM50 13 2、PAM250 15 2、BLOSUM90 10 1（エンジンが matrix ごとに解決する値と同じ。unit test が `describe.json` と照らす） |
| Gap Costs（BLASTX） | 11 1 だけ（LOSAT の BLASTX が読む値、DW-36） | 11 1 |

task・組・matrix を選んで Gap Costs の一覧が変わると、Gap Costs はその一覧の既定になる（NCBI の画面の `changeOptionList` と同じ。BLAST+ は一覧に無い組を止めることがある：task blastn の 4,-5 に既定の 5 2 は「Gap existence and extension values 5 and 2 are not supported for substitution scores 4 and -5」、megablast の 1,-1 に 0 0 も止まる。CLI 2.16.0 と LOSAT で確かめた）。一覧が同じ変更（blastp と blastp-fast の切り替え）や、一覧の無い値（「Other」で打った matrix や組）では書いた値を残す。ただし一覧が同じでも、使う組が一覧から外れるとき（task blastn で 1,-1 を選び megablast に替えると既定の 0 0）は一覧の既定にする。組を選ぶと 2 つの値を書く（既定の組を選んだときは書かない）ので、task を替えても選んだ組のまま（megablast で選んだ 1,-4 は blastn でも 1,-4。argv には既定と違う値だけが出る）。

**横の一文**：NCBI の画面で欄の横に見えている文は Max target sequences の「Select the maximum number of aligned sequences to display」だけ（ほかの文は「?」の helpbox の中、checkbox の名前は「Low complexity regions filter」「Mask for lookup table only」などの label）。この一文を同じ位置に出す（`ParameterField.hint`）。ほかの行に LOSAT の文は足さない（§4）。「?」は今までどおり CLI の help を開く。

**既定値**（line WR、2026-11-02。同日の DW-37 で `describe` が返す形に）：欄に示す既定値はエンジンが解決した値で、`describe` の `defaults`（`web/adapter/src/defaults.rs`、`abi_v2.md` §9）と clap の既定から来る。アダプタはエンジン自身の解決（BLASTP・TBLASTN・BLASTX の `resolve()`、BLASTN の `coordination` の関数、TBLASTX の getter）を task ごとに、また `-matrix` の 9 つ（`ScoringMatrix::ALL`）と `-word_size` 1〜32 の値ごとに呼び、値の変わる option だけを表にする。アプリ（`parameters.ts` の `knownDefault`）は表を引くだけで、program ごとの値も規則も持たない：task の値、書いた `-matrix` か `-word_size` のその値での値（両方を書くと、両方で変わる option は不明）、`unless` の option を書くと不明。各 task の値が CLI の help（または clap の既定）の言う値と同じことは、アダプタの試験（`defaults.rs` の `every_task_default_is_a_value_that_the_help_names`、全 option）が照らす。与えたかを区別する option（`given`）は task ごと：hit list の 3 つ（全 task）と、fast の task の `-word_size`（blastp-fast・tblastn-fast・blastx-fast は 5 を与えると threshold 19.3、与えなければ 20）。select の既定の項目を選び直すと欄は空に戻る（与えない）ので、fast の task で 5 を与えるには「Other」に 5 を打つ。

| program | task で変わるもの | 他の option で変わるもの |
|---|---|---|
| BLASTN | `-evalue` 10（blastn-short 1000）、`-word_size` 28 / 11 / 11 / 7、`-reward` 1 / 2 / 2 / 1、`-penalty` -2 / -3、`-gapopen` 0 / 5、`-gapextend` 0 / 2、`-dust` 20 64 1（blastn-short no）、`-soft_masking` true | — |
| BLASTP | `-word_size` 3（blastp-fast 5）、`-evalue` 10、`-matrix` BLOSUM62、`-seg` no、`-soft_masking` false | `-gapopen`・`-gapextend`・`-window_size`・`-threshold` は `-matrix` ごとのエンジンの値（PAM30 で 9・1・15・16）。`-threshold` は与えた `-word_size` での値（5 で 19.3、6 で 21）。`-word_size` は与えたかを区別する（blastp-fast で 5 を与えると threshold は 20 でなく 19.3） |
| TBLASTN | `-word_size` 3（tblastn-fast 5）、`-evalue` 10、`-matrix` BLOSUM62、`-seg` 12 2.2 2.5、`-soft_masking` false、`-xdrop_gap` 15、`-xdrop_gap_final` 25、`-sum_stats` true | BLASTP と同じ（`-threshold` 13、tblastn-fast 20） |
| TBLASTX | `-word_size` 3、`-matrix` BLOSUM62、`-soft_masking` false、`-culling_limit` 0 | `-threshold`・`-window_size` は `-matrix` ごと（BLOSUM45 で 16・60） |
| BLASTX | `-word_size` 3（blastx-fast 5）、`-matrix` BLOSUM62、`-gapopen` 11、`-gapextend` 1、`-window_size` 40、`-culling_limit` 0、`-sum_stats` true（`-evalue` 10、`-seg` 12 2.2 2.5、`-soft_masking` false、`-comp_based_stats` 2、`-strand` both は clap の既定。help の文は線 SX が CLI に足した、DW-36） | `-threshold` 12（blastx-fast 20）は与えた `-word_size` での値（5 で 19.3）。`-matrix` は BLOSUM62 だけで何も変えない。`-word_size` は blastx-fast で与えたかを区別する |
| 5 つとも | `-max_target_seqs` 500、`-num_descriptions` 500、`-num_alignments` 250 | `-max_target_seqs` は他の 2 つが、他の 2 つは `-max_target_seqs` が空の間だけ（BLAST+ は与えたかを区別し、互いに除く） |

Template length・Template type は既定を示さない（2 つを共に与える option で、片方だけ書くと CLI が拒む）。Max HSPs・Percent identity は既定が無い（空 = 制限なし）。2026-11-02 の版はアプリの `domain/programs.ts` に program ごとの `defaults` の表（CLI の help の写し）と「その option が既定の間だけ示す」の規則を持ち、`-matrix` を変えると関わる欄は空になった（Owner-delegated 2026-11-02 の仮の出どころ）。

## 2. 結果画面

NCBI の今の結果は、上から操作の行（Edit Search、Save Search、Search Summary）、見出しの塊（Job Title、RID、Program、Query ID…）と右の「Filter Results」、タブ（Descriptions、Graphic Summary、Alignments、（Taxonomy）、Dot Plot）、タブの中身、の順である。LOSAT は見出しの塊までこの順に置き、その下は NCBI の以前の（classic の）結果のページのとおり、Graphic Summary、Descriptions、Alignments を 1 ページに上から下へ並べる（保守者の指示 2026-10-10「Graphic Summary の下に Description が直接出てくるやつ」、2026-11-02「タブで分けるんじゃなくて上下に配置して（NCBI の以前の方式）」）。Dot Plot と LOSAT の Run details・Outputs はタブのまま。

### 見出しの塊

| NCBI の要素 | LOSAT Web | 採否 | 理由・置き場所 |
|---|---|---|---|
| 「Edit Search」 | なし | 後の段階（S15） | Run の入力と条件を検索画面に戻すのは、再現性（再読込とつなぎ直し）の段階で扱う |
| 「Save Search」「How to read this report?」「BLAST Help Videos」「Back to Traditional Results Page」 | なし | 採らない | NCBI のアカウントと文書 |
| 「Search Summary」（Search Parameters・Karlin-Altschul statistics・Results Statistics） | タブ「Run details」（`results-view-details`）：CLI のコマンド、RunSnapshot、RunRecord、検証バッジ | 寄せる | 統計の値（Lambda、K、H、有効探索空間）は outfmt 0 の末尾の原文にある（Outputs）。TS で抜き出さない |
| 「Job Title」 | 「Job Title」（Run の `title`。無ければ行を出さない） | 採る | |
| 「RID」と「Search expires on」 | 「Run」：Run の select（`results-run`、Run の番号・program・入力名） | 寄せる | RID に当たるのは Run の番号。期限は無い（作業セッションの中） |
| 「Download All」 | 「Download All」：タブ「Outputs」へ移る（`results-download-all`） | 寄せる | 原文の outfmt 0/6/7 と診断を書き出す（Outputs） |
| 「Program」「Blast 2 sequences」と「Citation」 | 「Mode」（DW-32）：BLASTN（task megablast）など。横に検証バッジ（`verification-badge`） | 寄せる | Citation は採らない。検証バッジは LOSAT だけ（NCBI と比べた記録） |
| （無い） | 「Options」：argv の option（`results-run-options`） | LOSAT だけ | |
| 「Query ID」「Query Descr」「Query Length」 | 「Query ID」と「Query Length」（1 query の Run） | 寄せる | Query Descr は採らない（Run の記録が持つのは ID と長さ。outfmt 0 の `Query=` の行は Alignments と Outputs で読める） |
| 「Subject ID」「Subject Descr」「Subject Length」 | 「Subject ID」と「Subject Length」（subject が 1 レコードの Run）。複数なら「Subjects」：入力名とレコード数 | 寄せる | Subject Descr は Descriptions と Alignments に（outfmt 0 の見出し） |
| 「Other reports」（MSA viewer、Distance tree） | なし | 採らない | LOSAT に無い |
| 「Results for」（複数 query の select） | 「Results for」：Query の一覧（`query-list`、仮想化、ID の絞り込み `query-filter`、「With hits only」`filter-hits-only`）。複数 query の Run だけ | 寄せる | 10 万 query を 1 つの select には入れられない。仮想化した一覧を同じ場所と名前で置く |
| 当たりの無い query：黄色の帯「No significant similarity found.」、タブが消える | 「No significant similarity found for this query.」（`results-notice` の `no-hits`）。タブは残す（Run details と Outputs は query に依らない） | 寄せる | outfmt 0 の「No hits found」と同じ事実。値は作らない |
| 「Filter Results」：Percent Identity、E value、Query Coverage の from / to、「Filter」「Reset」 | 「Filter Results」：E value ≤、Bit score ≥、Subject ID contains、「Filter」（`filter-apply`）「Reset」（`filter-clear`） | 寄せる | Percent Identity は採らない（レコードに原値が無い、`results_columns.md`）。Query Coverage はエンジン待ち（E2j）。表示だけを変え、検索をやり直さない（ViewState の規則は変えない） |

### タブと 1 ページの並び

タブは「Hits」（`results-view-hits`：Graphic Summary・Descriptions・Alignments を 1 ページに）、「Dot Plot」（`pane-dotplot`）、その後に LOSAT の「Run details」（`results-view-details`）と「Outputs」（`results-view-outputs`）。2026-11-02 までは NCBI の今のページと同じ 3 つのタブ（Descriptions、Graphic Summary、Alignments）だった。「Taxonomy」は採らない（ローカルの FASTA に無い）。NCBI のタブの形（選んだタブを塗り、帯の下に道具の行）を LOSAT の色で写す。どのタブも同じ選択（中心は HSP の ID、`ResultsBrowser`）に従う。

「Hits」の中は上から「Graphic Summary」（`section-graphic`）、「Descriptions」（`section-descriptions`）、「Alignments」（`section-alignments`、選んだ subject の塊）の節で、Descriptions の行か Graphic Summary の棒を選ぶと Alignments がその subject に変わる。Graphic Summary の棒のクリックはその Range を見える所へ送る。

### Descriptions

| NCBI の要素 | LOSAT Web | 採否 | 理由 |
|---|---|---|---|
| 見出し「Sequences producing significant alignments」 | 同じ見出し | 採る | |
| BLASTP の最初のタブ「Clusters」（「Clusters producing significant alignments」、Cluster Representative Sequence） | 「Descriptions」のまま | 採らない | NCBI のデータベースの配列のまとまり。ローカルの subject はまとめない |
| 「Download」「Select columns」「Show」、「select all」と行の印、「GenBank」「Graphics」「MSA Viewer」 | なし | 採らない（Download は Outputs、行の選択は S14 で検討） | 抽出（S14）が行の選択を求めれば、その段階で足す |
| 列 Description（リンク） | 「Description」：outfmt 0 の見出しの題（`results_columns.md`）。行を選ぶとその subject が選ばれ、Alignments・Graphic Summary・Dot Plot が従う | 寄せる | 行は 1 つのボタン（中にリンクを入れない）。選んだ subject の Alignments は同じページの下に出る |
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
| 「Alignment Scores」の凡例：< 40 黒、40 - 50 青、50 - 80 緑、80 - 200 マゼンタ、>= 200 赤 | 同じ階級。色は Okabe–Ito のカラーユニバーサルデザインの 5 色（黒 #000000、青 #0072b2、青緑 #009e73、赤紫 #cc79a7、朱 #d55e00）：保守者の指示（2026-11-02）「案 4（色覚多様性対応）でいきましょう」。赤と緑の対比に頼らない | 寄せる | 読み慣れた機能上の約束（ブランドではない）。HSP レコードの `bit_score`（原値）を区切るだけ |
| 「Show Conserved Domains」 | なし | 採らない | NCBI の CDD |
| 「Distribution of the top N Blast Hits on M subject sequences」 | 「Distribution of N HSPs on M subject sequences」（見せている subject の数） | 寄せる | |
| query の棒と目盛り | 同じ：棒「Query」と目盛り（1 と query の長さ、その間の切りのよい位置）、単位は query の nt / aa | 採る | 座標は Run の記録の長さ |
| subject ごとの行、HSP ごとの細い棒（同じ subject の HSP は 1 行、灰色の細い線でつなぐ） | 同じ。Descriptions の順と表示用のフィルターに従い、最初の 100 subject（「Show all M」で全部） | 採る | 棒の位置は HSP レコードの query の座標 |
| hover：題、Score、Evalue | hover と focus：Subject ID、Description（読めていれば）、outfmt 6 の bitscore と evalue | 寄せる | |
| click：その alignment へ | click（と Enter）：その HSP を選び、Alignments のその「Range」へ | 採る | キーボード：上下で subject、左右で HSP |

### Alignments

| NCBI の要素 | LOSAT Web | 採否 | 理由 |
|---|---|---|---|
| 「Alignment view」（Pairwise、Query-anchored など）、「CDS feature」、「Restore defaults」、「Line length」 | なし | 採らない | エンジンが書くのは outfmt 0 の pairwise だけ。原文を書き換えない（`web/AGENTS.md` 規則 2） |
| subject ごとの塊（全部の subject を続けて） | 選んだ subject の塊（`alignments-subject`）。「Previous」「Next」で前後の subject（一覧はすぐ上） | 寄せる | 数千の subject を 1 ページに並べず、選択に従う |
| 「Download」「GenBank」「Graphics」「Sort by」 | なし | 採らない（Download は Outputs） | |
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
| 出力 | SVG のファイル | 画面（Canvas）。SVG の書き出しは採らない（S15 の出力の段階で検討） | 出力の形式と名前は S15 が決める |

**描く subject の選び方**（line WR、2026-11-02、DW-35。保守者「Dotplotの対象配列選択の挙動が謎なんだけどこれわかりづらくない？」）：Dot Plot の道具の帯の先頭に「Subject」の select（`dotplot-subject`）。選んだ query の subject を Descriptions と同じ順に「<sseqid> · <長さ> <単位> · <n> HSPs」で並べる。Descriptions の行と同じ選択で、どちらで選んでももう一方が従う。query は今までどおり結果の上の「Results for」で選び、題「Plot of <query> vs <subject>」に両方が出る。以前は Dot Plot のタブに選ぶ所が無く、Hits のタブで最後に選んだ行の subject を描いていた。

LOSAT が足すもの（W4 の判断 18 を保つ）：ズーム（ボタン、+ / −、Ctrl / ⌘ とホイール）、パン（ドラッグ、矢印キー）、HSP の選択（線のクリック、n / p、表）、「Zoom to HSP」「Whole sequences」、`touch-action: pan-y`。新しく、hover で線を太くし、HSP を選ぶ（クリック、Enter）とその HSP のポップアップ（`dotplot-popup`）を出す：outfmt 6 の行の値をそのまま（bit score、E value、identity、query と subject の範囲、翻訳した配列の frame、向き、outfmt 0 に有るか）と「Show alignment」（Alignments のその Range へ）。ポップアップは線の端の脇に置き、線の中点を覆わない（そこに場所が無いときと、描く枠が 600 px 未満の電話の幅では図の下に）。ポップアップは Escape と「Close」で閉じ、キーボードとスクリーンリーダーで届く（focus を移し、閉じたら図に戻す）。描き方は層を分ける：格子と HSP の層（色と不透明度ごとにまとめ、線は 1 本ずつ描き、見えない線（表示の外の線、両端が前の線と同じ画素に来る不透明な線）を省く。重なる多数の線を 1 本の path にすると 2〜3 倍遅かった（W4b の B1）。薄い線の重なりはスクリプトの SVG の線と同じく濃くなる）と、hover・選択の層。ズームとパンは 1 フレームに 1 回だけ描き直す。

## 3. LOSAT だけのものの置き場所

| LOSAT の要素 | 置き場所 |
|---|---|
| Run の状態、取消 | 結果のページの中：終わっていない Run は「Run」の select の下に状態の行（`run-status`：段階、経過時間、待っている理由、Cancel・Cancel the group）。Run の一覧は「Run」の select（状態つき）。2026-11-02 まで右にあったキューの領域（カードと「Open results」）は保守者の指示「検索の queue 要らない」で無くした。検索は今も 1 つずつ順に走る（coordinator は変えない） |
| 離れることへの備え（wake lock、戻ったときの確認）、保存の状態 | ページの下の注の行（`notes`、広い画面で 2 列）。wake lock の選択「Keep searching while you're away」は既定 ON で、手持ちの端末（`(pointer: coarse)`）だけに出す。理由は「?」の後ろ（DW-41、§4） |
| Combined / Separate、レコード一覧と除外 | 検索画面の各枠の中 |
| Subject の保持（R1） | Run details |
| 検証バッジ | 結果の見出しの塊の Program の横、Run details |
| Outputs（outfmt 0/6/7 の原文、診断、書き出し） | タブ「Outputs」と見出しの「Download All」 |
| FakeEngine の帯 | 画面の上（W0） |

## 4. 文言の方針（stop-describing-ui、2026-11-02）

保守者の指示（2026-11-02）「基本これでいこう。NCBI 由来の説明文は残す」で、LOSAT Web の文言は [artouc/stop-describing-ui](https://github.com/artouc/stop-describing-ui) の 9 の規則に従う：1 機能でなく操作を名乗る、2 具体的な動詞、3 最も重い帰結だけ残す、4 内部の仕組みは書かない、5 同じ操作で戻せることは説明しない、6 画面が示していることを繰り返さない、7 無いものを言わない（試して失敗する場合を除く）、8 有効にする方法でなく操作そのものを見せる、9 短い行、記号の区切りや長い散文を避ける。規則に通らず何も残らない文は何も出さない。

- **残す（NCBI の文）**：mode のタブの下の一文（「BLASTN mode searches…」）、Job Title の「Enter a descriptive title for your search」、Algorithm parameters の「Parameter values that differ from the default are highlighted in yellow and marked with ♦ sign」と「Restore default search parameters」、Run LOSAT の横の「Search nucleotide subjects with BLASTN, task megablast (…)」（LOSAT が足していた「Runs in this browser.」は消す）、Graphic Summary の「hover to see the title」「click to show alignments」と「Distribution of N HSPs on M subject sequences」、Descriptions・Alignments の見出しと Range の語、「No significant similarity found for this query.」。
- **消す（LOSAT の説明）**：見出しの「Local processing: your sequences stay in this browser.」（鍵のアイコンと「Nothing leaves this browser」に）、Threads の説明、「The engine accepts these options.」と「Checking the options…」（拒否だけ出す）、下書きの readiness の行（押したときの誤りだけ出す）、「Preparing the inputs…」「Search started.」（ボタンが「Starting…」になり、ページが Run へ移る）、subrange の「The whole record is searched…」と「A subject region can be set when…」（当てはまらないときは何も出さない）、「Files are read in place, not copied into the box.」、「The BLASTN engine reads N records.」（チェックのアイコン。拒否の文は残す）、Program Selection の task の列挙（`describe` の help は各欄の「?」の後ろ）、各 option の help の段落（「?」ボタンで開く。欄は既定値を値として示す、§1「既定値」）、Filter Results の注、Graphic Summary の説明文、dot plot の Mouse / Touch の説明と「Each line is an HSP…」、表の「Scroll the table sideways…」、Outputs と Run details の説明の行、「While searching」の段落（検索中の Run のページに「Closing or reloading this tab ends the search」の 1 行だけ）、「Temporary storage」の見出しと使用量（メモリへの退避と掃除の不能の警告だけ）、Run の状態ページの説明（「Waiting for run N」と段階の時刻だけ）、結果タブの空の文（「Run a search」のボタン）。
- **短くする**：検証バッジ「Engine-supported, outside certified profile」→「Outside certified profile」（2026-11-02 に「Supported」／「Not supported」へ。保守者「Certify って別に公式じゃないじゃん。使えます、だけでいいじゃん」「使えれば Support、使えなければ Not supported でいいと思う」、DW-36 (7)）、「Keep the screen on while a search runs」→ Run LOSAT の行の「Keep screen on while searching」（「The screen stays on until the search ends.」は消す）、subject の遺伝暗号の注、Combined / Separate の説明、重複 ID と塩基・アミノ酸の見た目の注（2 文目を消す）、BLASTX の「not available」の文。
- **2026-11-03（DW-41、保守者の 3 点）**：(1) 見出しは「LOSAT」と一文「A miniaturized reimplementation of NCBI BLAST algorithm」（リポジトリの説明）。鍵と「Nothing leaves this browser」は外した（accession や例データの取得の線が入れば無条件には言えない）。(2) 「?」の help は NCBI の helpbox と同じく、マウスやペンで乗せる・キーボードで focus すると出て、クリックで留まる（もう一度のクリック、Escape、外のクリックで閉じる。タッチはタップで開閉）。2026-11-03 まではクリックだけだった（保守者「NCBI はマウスオーバーするだけで何かしら出てくる」）。欄の下に浮かぶ箱（`HelpTip.vue`、`.help-tip`）で、ページは動かない。(3) wake lock の選択は「Keep searching while you're away」（保守者「"スクリーンが ON" が検索の駆動とどうつながるのか直感的に理解できる人は少ない」。「タブを見ていないときも検索を続ける」とは書かない：wake lock は画面を消さないだけで、別のタブやアプリへ移ると検索は止まり得る、設計 §9.4）、既定 ON（保守者「規定ONにしとこうぜ」）、手持ちの端末だけに出す（机上の PC では意味が無い。PC でも既定どおり lock は持つ）。「?」：「Keeps the screen on until the search ends. A phone or tablet ends the search when its screen turns off; switching to another app can end it too.」
- **2026-11-03（DW-42、SURVEY #4）**：Other Parameters の最後の行の名前は「Additional options」（欄だけ、placeholder なし：既定は「何も足さない」）。「?」：「Other options of the losat command line, written as there: -qcov_hsp_perc 50 -ungapped」。拒否はエンジンの文のまま（「error: unknown option or argument '-foo'; use -help for CLI v2 syntax」）。アプリの文は、フォームの欄の名前を言う「-evalue is set in the form: Expect threshold」（操作の場所を示す、規則 8）、「-num_threads is managed by LOSAT」、「Additional options: a quote (') is not closed」、「Additional options: '' is not accepted」だけ。

- **2026-11-03（線 WF-2、Example、`fetch_design.md` §9）**：行の名前は「Example」、ボタンは組の名前（操作の対象を名乗る、規則 1）。その場の一文は保守者の決定の「Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov」（実測の大きさと host）。保存後は「Saved in this browser (2026-10-27, 41.8 MB)」（版：FTP は Last-Modified の日付、UniProt は「UniProt release 2026_03」）と「Delete」。取得中は「Downloading 3.1 of 8.0 MB from ftp.ncbi.nlm.nih.gov…」（最初の bytes の前は「Downloading from ftp.ncbi.nlm.nih.gov…」）と「Cancel」、写しを開く間は「Opening the saved copy…」、保存中は「Saving in this browser…」。失敗は host と理由の一文（「rest.uniprot.org answered HTTP 500: try again in a moment.」「ftp.ncbi.nlm.nih.gov could not be reached (…).」「The file from ftp.ncbi.nlm.nih.gov is not FASTA.」）、保存を断られたら「Not kept for the next time: …」。「?」：「Example puts a query in the query box and a public set in the subject list. It asks only for the file named, from the site named, and keeps the file in this browser: the next time, it is read from there without a request.」と組ごとの出典と利用条件。設計の Temporary storage の表示（使用量と Delete）は 11-02 に外したので、写しの有無と Delete は Example の行に置く。

- **2026-11-03（線 WF-3、競走と保存先、`fetch_design.md` §10）**：取得先の選択の先頭は「Fastest」（既定。競走の後は「Fastest (NCBI)」）。その場の一文は競走の前は「Sends the accession to the fastest of NCBI, EBI and CNCB」（その accession を扱える所だけを並べる：RefSeq の核酸なら「of NCBI and CNCB」）、後は WF-1 と同じ「Sends the accession to NCBI」。競走の間は「Finding the fastest site…」。どこも全部は扱えない欄は「No one site has all of these accessions」（押すと「No one site has all of these accessions: C_AA004835.1 only from CNCB, P12345 only from NCBI or EBI.」）。「?」に「Fastest: the first fetch asks each site for one public record (MH011443.1) and sends the accessions to the first that answers; this tab keeps that order.」を足した（タンパク質は「NP_414542.1, at UniProt P0AD86」）。Example の保存した写しの横に「Save to folder…」（フォルダを選ぶ API のあるブラウザ）か「Download」（他）。書いた後は「Saved bacteria.16SrRNA.fna in Downloads.」、失敗は「Not saved: …」。「?」に「Save to folder… writes that file outside this browser's storage, which clearing the site data empties; "Open FASTA files…" opens it again.」。

- **2026-11-03（線 WF-1、accession の取得、`fetch_design.md` §8）**：欄の見出しは NCBI の「Enter accession number(s), gi(s), or FASTA sequence(s)」から gi を除いた「Enter accession number(s) or FASTA sequence(s)」。欄が accession だけのときの行は「Fetch」（操作を名乗る、規則 1・2）「from [NCBI]」と、保守者の決定のその場の一文「Sends the accession to NCBI」（複数なら accessions、選んだ先の名前）。約束は「?」の後ろ：「Fetch sends only the accessions, to the site chosen here, and writes the FASTA it gets into the box. Sequences, files and results stay in this browser. EBI answers from ENA for nucleotides and from UniProt for proteins, CNCB from GenBase.」。取得中は「Fetching 1 of 2…」、失敗は取得先と accession を名乗る一文（「NCBI has no nucleotide record AB999999.1.」「EBI (ENA) has no RefSeq records: fetch NR_074540.1 from NCBI or CNCB.」など）。

## 5. 見た目

保守者の指示（2026-11-02）「色とかボタンを gbdraw の GUI っぽく」「フォントも Inter で」。配置・節の順・言葉は NCBI のまま（§1〜§4）、色・字・部品の形は gbdraw の Web GUI（`gbdraw/web/index.html` の `--ui-*` と `css/tailwind.input.css` の `.btn` `.card` `.form-input`）に合わせる（`web/app/src/ui/styles.css` の `:root`）。

| 要素 | gbdraw | LOSAT Web |
|---|---|---|
| 字 | Inter（同梱の woff2、Noto Sans JP が後ろ） | 同じ。`web/app/public/fonts/inter/`（latin・latin-ext の 400〜700、SIL OFL 1.1、`LICENSE.txt`）。CSP は `self` だけなので外部の font は読まない |
| 地と文字 | ページ `#eef2f6`、カード `#fff`、文字 `#1e293b` / `#334155`、薄い文字 `#64748b`、線 `#e2e8f0` / `#cbd5e1` | 同じ値（`--page` `--bg` `--fg` `--fg-2` `--muted` `--border` `--border-strong`） |
| アクセント | 青 `#2563eb`、hover `#1d4ed8`、薄い青 `#eff6ff` | 同じ（`--accent` `--accent-hover` `--accent-bg`） |
| カード | 角 12 px、1 px の線、薄い影 | NCBI の各ブロック（Enter Query Sequence …）、Run LOSAT の行、Algorithm parameters、Results がカード。Algorithm parameters の中は下の「Algorithm parameters の形」 |
| ボタン | `.btn-secondary`：白、線、角 8 px、太字、影。`.btn-primary`：青地に白 | 既定のボタンが secondary、`Run LOSAT` と `Filter` が primary |
| 切り替え | mode switch：灰色の丸い枠の中で、選んだものが白くなり青い字 | Search / Results と結果の Hits / Dot Plot / Run details / Outputs は灰青の枠の中で、選んだものをアクセントで塗り白い字（案 B、DW-33）。mode のタブは下の「案 B」 |
| 入力欄 | `.form-input`：角 8 px、薄い灰色の地、focus で青い枠 | 同じ |
| 通知 | info は青、warning は琥珀 | `.notice` が青、`.warning` `.banner` と ♦ の注が琥珀 |

アイコンは gbdraw の Phosphor を取り込まず、LOSAT の少数のインラインの SVG のまま。

**画面全体の層**（line WR、2026-11-02。保守者「Enter Query Sequence とか General Parameters とかの階層の box、背景白いとちょっと明るすぎる」「背景とかも含めてデザインしなおそう」）：NCBI の画面の 3 層（青みの地、灰色の組、白い入力）に揃える。地 `--page` `#d6e3f1`、組 `--group` `#eef1f4`（Enter Query Sequence・Enter Subject Sequence・Program Selection・Run LOSAT の行・Algorithm parameters の帯とその組・結果の見出しの塊・Filter Results・Graphic Summary・Descriptions・Alignments・Dot Plot・Run details・Outputs）、組の線 `--border-strong`、組の中の帯（表の見出し、道具の行）`--band` `#e2e7ec`、入力・表・出力・レコードの一覧は白。組の名前は全部アクセントの太字のつまみ（`legend`、`.section-legend`、`.section-heading`）。カードの影はやめる。検索のページは 14 px、「Or, upload file」の点線の枠は 801 px 以上で中身の幅（落とす先はブロック全体）。地と組の色は 4 案（`web-refresh-20261102/screens/tones/`：t1 明るい、t2 濃い、t3 NCBI の青、t4 無彩色）を実画面で比べ、t3 を既定にした（Owner-delegated、3 つの token で替えられる）。

**案 B**（line WR、2026-11-02、DW-33。保守者「灰色で見ずらい。ていうか全体的に見ずらい」→ 7 案のキャンバス https://claude.ai/artifact/3mMN84B4YYgL1NG66Pmr4t から「Bにします」）：上の「画面全体の層」の色を替える。見出しの帯は紺 `--header` `#0b3a5e` に白い字（鍵と「Nothing leaves this browser」は `--header-muted` `#d4e6f4`）。地 `--page` `#e9eef3`、組 `--group` は白、組の中の帯 `--band` `#f1f4f7`、線 `--border` `#d3dbe4` / `--border-strong` `#a9b8c7`、入力欄の線 `--control-line` `#7b8a9a`。文字は `--fg` `#0b1a2a`、`--fg-2` `#1c2b3a`、薄い文字 `--muted` `#46576a`（白の上で約 7:1。使える文字はどれも 4.5:1 以上）。mode のタブは NCBI の以前のタブの列の形：灰青のつまみ（`--tab` `#d6dfe8`、字 `--tab-fg` `#14263a`）がアクセントの 3 px の線に載り、選んだタブはアクセントで塗り白い字（使えないタブの `--tab-off` `#e6ebf0` / `--tab-off-fg` `#7d8896` は、BLASTX が使えるようになって外した、DW-36）。800 px 以下では 5 つのタブが 1 行を等分する（折り返すとつまみが線の下に垂れる）。組の名前のつまみは灰青 `--legend` `#e3eaf1` に濃い字。タブの下の一文は `--fg-2`。

**Algorithm parameters の形**（line WR、2026-11-02。保守者「NCBI と比べるとやっぱ見づらい」）：NCBI の画面の読み方（太い名前、節の名前のつまみ、名前のすぐ横の内容の幅の欄、残りの幅に一文）を gbdraw の見た目の中で取る。開いた中身は地の上に組として並び（上の「画面全体の層」）、節の名前はつまみ（NCBI の legend の位置）にアクセントの色・太字・1.1 em（文字は NCBI のまま「General Parameters」、大文字にしない）。字は 14 px、行は 40 px 以下（1 行の名前で 33〜38 px）。名前は太字（600、`--fg`）、flag はその後ろに小さい薄い字（等幅の 2 行目にしない。W4b：NCBI の名前と LOSAT の flag）。全部の節が 1 つの名前の列（最も長い名前の幅、10〜22 em、CSS subgrid。fieldset には subgrid が効かないので節は `role="group"` の div）を共有し、欄はその横から内容の幅（数は 8 em、文字の値は 12 em、select は選択肢の幅）。欄の線は `--control-line` `#8b9bab`、地は白。「?」は NCBI と同じくアクセントで塗った円に白い「?」（20 px、800 px 以下は 24 px）。Run LOSAT の横の一文は、何を検索するかと program・task を NCBI と同じく色付きの太字にする。800 px 以下は W4b のまま 1 列（名前の下に欄）。1920 px と 1280 px の検査は `web-refresh-20261102/screens/params-*.png`。

色は gbdraw の青ではなく Okabe–Ito のカラーユニバーサルデザインの色：アクセント `#0072b2`（青）、OK `#007a5a`（青緑）、エラー `#b84f00`（朱）、注意 `#faf6d2` / `#5b4b00`（黄）。Graphic Summary のリボン（§3）も同じ 5 色。ページの背景は無彩色の薄い灰色 `#eef1f4`、面（`--panel`）と線（`--border`）は少し青みのある灰色。保守者の指示（2026-11-02）：青緑・黄・臙脂・案内標識の緑（`#007356`）を背景やアクセントにした案と比べて「案 4（色覚多様性対応）でいきましょう」。Dot Plot の線の青 `#1f77b4` とオレンジ `#ff7f0e`（blast2dotplot.py の色、§3）は青と橙の対比で色覚に優しいのでそのまま。

2026-11-02 までの見た目（W4b）：薄い灰色の面に青い見出しのつまみ、幅いっぱいの青い帯の「Algorithm parameters」、青く塗った結果のタブ、薄い青の表の見出し。gbdraw の見た目に置き換えた。

**ロゴとローダー**（DW-43、2026-11-03。保守者「Web UIの左上のロゴを…使ってほしい」「ロード中/ラン中に、gbdrawみたくAの三角形を左から一筆書きするアニメーションを付けてほしい」）：見出しの帯の「LOSAT」の字を保守者の描いたロゴ（Inter Bold の輪郭。A はまっすぐな 2 本の脚が平らな底辺に載り、底辺は右脚の外まで伸びる。2026-11-03 20:12 の版。その前の版の A は左下の開いた三角で、底辺の左端の斜めの切り方が Adobe の A に似るとして替えた）に替える。高さ 24 px、帯の白（`web/app/src/ui/logo.ts`、`LosatLogo.vue`）。帯の中の並びは縦の中央揃えで、ロゴの中ほどが説明文の大文字の中ほどに来る（保守者「説明文より浮いてるから中央揃えか下揃えに」。説明文のベースラインにロゴの下端を合わせると、下に出る字のぶんロゴが浮いて見えた）。見出しの名前は隠した字の「LOSAT」。ロゴは動かさない（検索は数分続くことがある）。ローダーは gbdraw の `.gbdraw-loader` と同じ作りで、マスクの中の白いペンが A の線の中ほどを左足の下 → 頂点 → 右脚と進み、底辺に着いたところで左右に分かれて底辺の両端まで伸び、塗った A を見せる（保守者「最後は左右に分岐するかたちで」。`LosatLoader.vue`、`.losat-loader`：ペン 2 本の区間を gbdraw と同じく `@property` の 1 つの位置で進める）。1 周 2.22 s（描く 0.9 s、止める 0.26 s、書き始めから消す 0.9 s、休む 0.16 s）、「動きを減らす」の設定では描き終えた形で止まる。色はアクセント、塗ったボタンの中はボタンの字の色。出すのは進行中の状態の言葉の前だけ：結果のページの Status（Waiting・Preparing・Searching・Organizing results）、Run LOSAT の「Starting…」、「Loading the engine…」（mode のエンジンの説明を待つ間、Run LOSAT の行。初回はエンジンの取得で数秒かかる。すでにあるエンジンではすぐ終わるので 0.4 s 続いたときだけ出す）、「Reading the results…」、開いた HSP の「Reading the alignment…」、「Reading the records…」、入力の「Checking…」、「Fetching … from …」、「Checking the stored results…」、Example の取得・読み出し・保存の間。Alignments の HSP ごとの読み待ちの欄（1 ページに数百ありうる）は字だけ。タブの icon は帯の紺の角丸に白い A（`public/favicon.svg`）。起動の画面は付けない：JS は gzip で 106 kB、CSS は 7 kB で、マウントまでの白紙は普通の回線で 0.2 s に満たない。

変わらない基準（S13、W4 の判断 16〜23）：1280 px で表が切れない、390 px で横に流れない、触れる対象は 24 px 以上、コントラスト、キーボード操作。Graphic Summary の凡例とドットプロットの線の色は §2 の約束のまま。

