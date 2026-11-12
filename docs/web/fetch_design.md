# LOSAT Web の取得（accession と Example data）の設計（線 WF）

作成 2026-11-03（線 WF、保守者の GO は同日 15:3x、§7）。保守者の 2026-11-03 の 2 点（5）モードごとの Example data を固定アドレスから取って
キャッシュし次回も使う、（6）accession を query / subject にして一番近い DB から取る（座標指定も）に対する設計。事実の記録（CORS の実測）は作業フォルダ web-fetch-20261103/NOTES.md。

## 1. 今の設計との関係

- 設計書 §13（`docs/web/losat_web_design_v0.1.md` 447-449 行）：「研究データは端末内で…」「『研究データを送信しない』と『ネット
  ワーク要求が一切ない』は異なる」。つまり、**配列・題・ファイル名・結果を送らない**ことが約束で、通信ゼロは約束ではない。
- `web/AGENTS.md` 規則 4 と `ncbi_ui_mapping.md` §1 の行「accession の取得は NCBI への送信になる。FASTA だけを受ける」は、この約束を
  広く取った決定。accession は研究データではない（公開 DB の識別子）が、「何を調べているか」は漏れる。
- CSP の `connect-src` は自 origin だけ（§13.3）。取得先の host を足す必要がある。
- SURVEY #11：Web は `-db` を持たない（2 配列だけ）。Example の「DB」は **subject の FASTA** として持つ（tier の makeblastdb は不要）。
- 見出しの「Nothing leaves this browser」は 11-03 に外した（DW-41）。取得が入るなら「送るのは accession と取得先の URL だけ」を取得の
  その場で言う。

**方針案**：約束を「研究データは送らない。利用者が明示に頼んだ取得（accession、Example）のときだけ、その識別子と公開ファイルの URL
を取得先（NCBI、EBI、CNCB）へ送る」に書き直す。ページを開いただけ・Run を押しただけでは通信しない。オフラインは今までどおり動く。

## 2. 取得先（ブラウザから直接。2026-11-03 実測、FETCH_NOTES.md）

| 取得先 | CORS | 用途 | 備考 |
|---|---|---|---|
| NCBI E-utilities `efetch` | `*` | accession → FASTA、`seq_start`/`seq_stop` で部分配列 | 鍵無し 3 req/s。`db=nuccore` / `protein` |
| ENA browser API `/ena/browser/api/fasta/ACC?range=a-b` | `*` | INSDC accession → FASTA（範囲も。題に `Location:a..b`） | RefSeq（NR_, NM_…）は 400 |
| UniProt REST `/uniprotkb/ACC.fasta` | `*` | タンパク質 accession | `stream?query=` は実測で 500（構文を要確認） |
| NCBI FTP `refseq/TargetedLoci/*` | `*` | Example の subject セット | 下の表 |
| NCBI FTP `blast/db/FASTA/pdbaa.gz` | `*` | タンパク質の Example 候補 | 36.8 MB gz（境界） |
| DDBJ getentry | 無し | — | CORS ヘッダ無し（curl には本文も空）。DDBJ Search の API は CORS を返すが、索引は BioProject / BioSample / DRA / GEA などで塩基配列のレコードを持たない（`/search/api/entries/nucleotide/…` は「Unknown database type」）。公開ファイル（`public/ddbj_database/16S/16S.fasta.gz`、222 MB）も CORS 無し。DDBJ が getentry に `Access-Control-Allow-Origin: *` を付ければ取得先に足せる |
| CNCB GenBase REST `https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=ACC` | 返す（origin をそのまま。`restapihelp` に記載） | accession → FASTA（全長。範囲指定の引数は無い） | GenBase の `C_AA…` に加え INSDC（MH011443.1、AB000001.1）と RefSeq（NR_074540.1）も返った。日本から 1.3〜2.6 s。`api/file/gbf` で GBF も |

DDBJ は DDBJ 側のヘッダ追加（依頼）があるまで取得先に入れない。それまで DDBJ 利用者は NCBI / EBI / CNCB から取るか、DDBJ から手で落として「Open FASTA files…」で開く。
CNCB の GenBase は取れる（2026-11-03 17:3x 実測）ので取得先に入れる（保守者 17:4x「取得先に足す」、§7 の 3）。
「一番アクセスしやすい所」= 候補（NCBI、EBI、CNCB）へ同時に Range 要求を出し、先に応えた方を取る。選んだ先はそのセッション中は固定し、
取得元の名前を画面に出す（手で NCBI / EBI / CNCB を選べる）。WF-3 では Range 要求ではなく固定の公開レコードで測り、accession の取得だけに使う（§10）。

## 3. Example セット（mode ごと、各 50 MB 以下、固定 URL、SHA-256 を表に持つ）

WF-2 で確定した中身と持ち方は §9（表に SHA-256 は持たない。TBLASTN / TBLASTX は E. coli K-12 のゲノム）。下は設計の時点の案。

| mode | subject | 大きさ（gz） | query の例 |
|---|---|---|---|
| BLASTN | RefSeq Targeted Loci bacteria 16S（`bacteria.16SrRNA.fna.gz`） | 8.0 MB | `NR_074540.1`（B. cereus 16S）の部分配列 |
| BLASTN（2 つ目） | fungi ITS（`fungi.ITS.fna.gz`） | 3.0 MB | ITS の accession |
| TBLASTN / TBLASTX | ファージのゲノム数本（efetch: λ `NC_001416.1` 49 KB、T4、T7） | < 1 MB | λ のタンパク質 1 本（TBLASTN）、T7 の断片（TBLASTX） |
| BLASTP / BLASTX | `pdbaa.gz`（36.8 MB）か、UniProt の reviewed の分類群セット（要再確認） | 37 MB / 数 MB | PDB の 1 鎖 / 核酸の CDS 断片 |

小さいもの（archaea 16S 0.3 MB、bacteria 23S 1.0 MB、fungi 18S 0.7 MB・28S 2.3 MB）は追加の候補。表はアプリ内の静的 JSON
（名前、URL、大きさ、SHA-256、出典と利用条件、mode）。

## 4. 画面（NCBI 寄せ）

- Query / Subject の欄の見出しを NCBI の「Enter accession number(s), gi(s), or FASTA sequence(s)」に戻す（gi は受けない：
  「Enter accession number(s) or FASTA sequence(s)」）。欄の中身が accession の形（`^[A-Z]{1,2}_?\d+(\.\d+)?$`、1 行 1 つ）なら
  「Fetch from NCBI」の操作が現れ、押すと取得して **FASTA を欄に写す**（写した後は普通の貼り付けと同じ。何が検索されるか見える）。
  範囲は既存の「Query subrange」From / To（NCBI と同じ）。取得時だけ efetch の `seq_start`/`seq_stop` を付ける案もあるが、
  欄に写す方式なら範囲は今までどおり `-query_loc` で足りる。
- 「Example」：mode のタブの下に「Example: 16S rRNA (RefSeq)」のような 1 行。押すと query を欄に、subject を「Or, upload file」の
  列に「16S_ribosomal_RNA (RefSeq, 8 MB)」として登録（キャッシュから、無ければ取得して進捗を出す）。
- 通信の表示：取得の行に「Sends the accession to NCBI」／「Downloads 8 MB from ftp.ncbi.nlm.nih.gov」。約束の文は "?" に。

## 5. データの持ち方

- 取得した FASTA は既存の source の経路に乗せる：小さい accession は貼り付けと同じ bytes、Example セットは Data worker の
  OPFS block store（`infra/data/opfs-block-store.ts`）に **キャッシュ**（URL + SHA-256 を鍵に。WF-2 では URL が鍵、§9）。次回は取得せずに読む。
  OPFS は origin 専用で、サイトデータの消去で消える。使用量と「Delete」は既存の Temporary storage の表示に載せる。
- 「任意のフォルダに保存」= File System Access API（`showDirectoryPicker`）。Chromium / Edge だけ。Firefox / Safari は
  「Download」で落として次回「Open FASTA files…」から開く（ファイル参照なので大きくても速い）。WF-3 では Example が読むのは
  OPFS の写しのまま、フォルダの写しは利用者のファイルで LOSAT は読み戻さない（§10）。
- 展開：`.gz` は `DecompressionStream('gzip')`（3 ブラウザとも可）で OPFS に書く。

## 6. 進め方（PR の切り方）

1. **WF-0 設計**：このメモを `docs/web/` に入れ、設計書 §13 と `web/AGENTS.md` 規則 4、`ncbi_ui_mapping.md` §1 の accession の行を
   書き換える。CSP `connect-src` に取得先を足す。
2. **WF-1 accession**：欄の accession 判定、「Fetch from NCBI / EBI / CNCB」、efetch と GenBase と ENA の port（`ports/fetch.ts`）と fake、
   e2e はローカルの stub サーバ。
3. **WF-2 Example**：静的表、OPFS キャッシュ、gz 展開、進捗と失敗の表示、mode ごとの 1 行。
4. **WF-3 保存先と mirror**：「Save to folder…」（Chromium）、NCBI / EBI / CNCB の競走、取得元の選択。

## 7. 決定（保守者、2026-11-03 15:3x）

1. **方針 = GO**：§1 の約束の書き換えと、WF-0（この文書と設計書 §13、`web/AGENTS.md` 規則 4、対応表 §1 の accession の行、CSP の
   `connect-src`）→ WF-1 accession → WF-2 Example とキャッシュ → WF-3 保存先と mirror の順。
2. **タンパク質の Example = 両方**：既定は UniProt の reviewed の分類群セット（数 MB。REST の stream 構文を再確認してから表に入れる）、
   pdbaa.gz（36.8 MB）も表に入れて選べる。
3. **取得先 = NCBI、EBI、CNCB**：15:3x は「NCBI と EBI だけ、CNCB は API を確認してから」。17:3x に GenBase の REST が CORS 付きで INSDC / RefSeq の accession も返すと分かり、17:4x「取得先に足す」。DDBJ は CORS が無いので入れない。
4. **文言 = その場の一文だけ**：取得の行に「Sends the accession to NCBI」「Downloads 8 MB from ftp.ncbi.nlm.nih.gov」。約束の文は
   "?" の後ろ。見出しは「LOSAT」と一文のまま。

## 8. WF-1 の実装（2026-11-03、Owner-delegated）

保守者の常設の指示（2026-09-29、2026-10-07 に再掲）により、推奨の案を採って記録する。作業記録は web-fetch-20261103/DECISIONS.md の 7〜13。

1. **取得先は手で選ぶ**（NCBI が既定、選んだ先はセッション中保ち、query と subject で共通。WF-3 で既定は「Fastest」になった、§10 の 3）。§2 の競走は §6 のとおり WF-3。利用者の
   accession を 3 か所へ同時に送る競走はしない（その場の一文「Sends the accession to NCBI」と食い違う）。WF-3 の競走は固定の公開
   accession か Range 要求で測る。
2. **accession の形**は §4 の `^[A-Z]{1,2}_?\d+(\.\d+)?$` を広げる：INSDC（1〜2 文字 + 5〜8 桁、タンパク質の 3 文字 + 5〜7 桁、WGS の
   4〜6 文字）、RefSeq（`NR_`、`NZ_CP…`）、GenBase（`C_AA…`）、UniProt（`P12345`、`Q9Y261`、`A0A023GPI8`）、末尾の版は任意。大文字小文字は
   問わず、大文字にして送る。数字の無い行（配列）と数字だけの行（gi）は当たらない（`web/app/src/domain/accession.ts`）。
3. **accession だけの欄は source を作らない**：読めない入力の誤りと「Add the line ">pasted_query"」の案内を出さない。Run は取得せず
   「Fetch the query accession first.」で止まる（ページを開いた・Run を押しただけでは通信しない）。取得中の Run も待たずに
   「The query accession is still being fetched.」で止まる。欄を書き換えると走っている取得は捨てる。取得中は取得先の選択を止め、
   「Fetching 1 of 2 from NCBI…」と送り先を出す。
   accession だけでない欄（線 WA の WA2-OP-04、2026-11-07）：1 行目が Seq-id として読まれる欄（accession や数字だけの行の後に
   別の文字）の読めない理由は adapter の文のまま（「…which is not supported in the browser (start the input with a '>' defline)」。
   engine が Seq-id と見ても Fetch が取らない GI、PDB の id、`sp|…|` があるので、adapter には Fetch の案内を持たせない）。
   card は、貼った文の先頭の空白と空行を除いて 1 行目の頭に `>` を足すだけのボタン「Put '>' before the first line」を出し（Seq-id の場面では「Add the line」を出さない。
   1 行目が最初の data 行になって読めない）、1 行目が `ACCESSION` に合う時だけ一文「To fetch accessions, put only the
   accessions in the box, one per line.」を足す（`pointsToFetch`）。空白・コンマで区切った accession だけの欄は 1 行 1 つでない
   ので Fetch を出さず（§4）、reader が読めるか（版の無い `AB000001 MH011443` は 4 nt の配列）に依らず card に
   「Put each accession on its own line」（押すと 1 行 1 つに書き換え、Fetch が出る。`accessionsOnALine`）。
   この時は一文と「Put '>'」は出さない（同じ事を分けるボタンがする。`>` を足すと 2 行目以降が data として読まれる）。欄を書き換えるボタン（「Add the line」「Put '>'」「Remove the byte order mark」）は、押すと card の読みが変わってボタンが消えるので、focus を card（`li`、tabindex -1）に移す（WA2-AX-04）。分けるボタンは Fetch へ移す。
4. **1 accession に 1 要求、順に**。NCBI は 340 ms 以上空ける（鍵無し 3 req/s）。1 回に 100 個まで、合計 10 MB まで（超えたら
   「Open FASTA files…」で開く案内）。全部取れたときだけ欄に写し、失敗は欄をそのままにして取得先と accession を名乗る。
   取得先が扱えない accession（5）が 1 つでもあれば、どれも送らずに断る。時間切れ（60 s）は応答の見出しまで（本文の読み込みは
   欄の書き換えで止まる）。
5. **EBI は核酸 = ENA、タンパク質 = UniProt**（ENA の browser API はタンパク質の id に CDS の核酸を返す、2026-11-03 実測）。ENA に RefSeq、
   UniProt に UniProt 以外の accession は送らずに断る。CNCB は両方（`NP_414542.1` が返る）だが、種類を指定する引数が無く、記録の
   種類のまま返す（タンパク質の欄に核酸が入り得る）。そのときは欄の入力の既存の注「1 record looks like nucleotide sequences, but
   BLASTP reads the subject as protein」が出る。
6. **要求は accession だけ**：cookie と referrer を付けない（`credentials: 'omit'`、`Referrer-Policy: no-referrer`）。NCBI の `tool`・`email`
   も付けない。
7. **e2e は Playwright の route で実アドレスを stub**（`web/app/tests/e2e/fetch.spec.ts`）。localhost の stub サーバでは配備の CSP の
   `connect-src` に例外が要る。route なら配備の CSP のまま試せ、ページから出た要求を全部数えられる。

## 9. WF-2 の実装（2026-11-03、Owner-delegated）

保守者の常設の指示（2026-09-29、2026-10-07 に再掲）により、推奨の案を採って記録する。作業記録は web-fetch-20261103/DECISIONS.md の 15〜19、
実測は同 NOTES.md「WF-2 facts」。

1. **表**（`web/app/src/domain/examples.ts`）。mode のタブの下の「Example」の行に、その mode の組をボタンで並べる。

   | mode | 組（ボタン） | 取得先（gz の大きさ、2026-11-03） | 展開後 | query（表に持つ） |
   |---|---|---|---|---|
   | BLASTN | 16S rRNA (bacteria) | NCBI FTP `refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz`（8.0 MB） | 41.8 MB、26,888 本 | NR_102783.2（B. subtilis 168 の 16S） |
   | BLASTN | ITS (fungi) | NCBI FTP `refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz`（3.0 MB） | 13.9 MB、20,390 本 | NR_111007.1（S. cerevisiae の ITS） |
   | BLASTP / BLASTX | E. coli K-12 proteins (UniProt reviewed) | UniProt REST `uniprotkb/stream?query=reviewed:true+AND+organism_id:83333&format=fasta&compressed=true`（0.98 MB） | 1.9 MB、4,531 本 | P16971（B. subtilis RecA）/ その CDS AL009126.3:1764645-1765691 |
   | BLASTP / BLASTX | PDB proteins (pdbaa) | NCBI FTP `blast/db/FASTA/pdbaa.gz`（36.8 MB） | 101 MB、156,894 本 | 同上 |
   | TBLASTN / TBLASTX | E. coli K-12 genome | NCBI FTP `genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/…_genomic.fna.gz`（1.4 MB） | 4.7 MB、1 本 | 同上（TBLASTN は RecA、TBLASTX は CDS） |

   query は表に埋め込む（数百 B〜1.6 kB）。1 回の操作で送る要求は 1 つで、その場の一文が名乗る host だけになる。§3 のファージ（efetch）は
   accession を 3 つ eutils へ送り、Content-Length も無いので、固定の版の E. coli K-12 のゲノムにした。話は 1 本：B. subtilis 168 の
   16S・RecA を RefSeq と E. coli K-12 に当てる（RecA の CDS の翻訳は P16971 と一致、単体試験で確かめる）。
2. **SHA-256 は表に持たない**。RefSeq Targeted Loci は release ごとに変わり（Last-Modified 2026-10-27）、UniProt も release ごとに
   変わる（2026_03）ので、固定した hash は次の release で Example を壊す。表は URL、その場の一文の大きさ（「Downloads 8.0 MB from
   ftp.ncbi.nlm.nih.gov」、実測の値）、ファイル名、出典と利用条件だけ。保存した写しが来たものを記録する：URL、版（FTP は
   Last-Modified の日付、UniProt は X-UniProt-Release）、取得と保存の大きさ、日時。鍵は表の id と URL（URL が変われば別物）。
   検索の入力は今までどおり SHA-256 で識別される（RunInput.sha256 とレコードごとの hash）。写しはそのまま使う（要求しない、
   オフラインで動く）。新しい release は「Delete」の後に取り直す。
3. **UniProt の stream の構文**は 11-03 の記録どおりで正しい。500 は UniProt 側の間欠的なもの（同じ URL が 500 を 3 回の後 12/12 回
   200、他の分類群も初回 500 の後 200）。失敗は「rest.uniprot.org answered HTTP 500: try again in a moment.」と出し、ボタンを
   残す。自動の再試行はしない（1 回の操作に 1 要求）。
4. **取得と展開**（`infra/browser/example-downloader.ts`）。ページ（main thread）が配備の CSP の下で取得する（`credentials: 'omit'`、
   referrer 無し。WF-1 と同じ）。NCBI の .gz は `application/x-gzip`、UniProt の `compressed=true` は Content-Encoding 無しで来るので
   ブラウザは展開しない。本文が gzip の magic で始まるときだけ展開する（途中で展開済みの本文はそのまま）。展開は 2026-11-11 から
   エンジンの decoder で行う（swebok-improve SW8 D2、Owner-delegated）：取得した本文（約 8 MiB ごとの Blob）を `<fileName>.gz`
   （サイト上のファイル名）の File として Data worker に渡し、source の `.gz` と同じ道（`infra/data/inflate.ts`、ABI v2 の
   `inflate_*`）で展開する。Web の展開器はエンジンの 1 つだけになる。それまでの `DecompressionStream('gzip')` は連結した
   member（bgzip、`cat a.gz b.gz`）を 3 ブラウザとも拒み、壊れたときの文もブラウザごとに違った。展開した部分は約 8 MiB ごとの
   Blob にまとめる。取得した本文が 256 MiB を超えるか（gzip の本文が展開後より大きいのは高々 0.03 % ほど）、展開後が
   256 MiB を超えれば止める。展開の間は圧縮した本文（最大の pdbaa で 36.8 MB）も Blob として持つ。エンジンの無い開発用の
   build は gzip の本文を展開できない（「could not be decompressed (this build has no engine)」）。
   時間切れは 60 s（応答の見出しまで）。進捗は「Downloading 3.1 of 8.0 MB from ftp.ncbi.nlm.nih.gov…」（Content-Length の無い
   UniProt は全体を言わない）と「Cancel」。失敗は host と理由の一文（HTTP の状態、届かない、gzip が壊れている（括弧の中はコマンドラインと同じエンジンの文）、大きすぎる、
   FASTA でない：`>` で始まらない 200 の答え（空や HTML）は写しにしない）で、欄はそのまま。
5. **保存**（`infra/data/example-store.ts`、Data worker）。OPFS の `examples/<id>/` に `data`、最後に `meta.json`（meta の無い
   ディレクトリは書きかけで、次の保存が消す）。書き込みは block store と同じ同期 access handle、保存と削除は Web Lock
   `losat-web:example:<id>` の下（2 つのタブが 1 つの写しを同時に書かない。保存中に他のタブが同じ URL を保存していれば、それを使う）。
   `tmp/`（作業セッション）の外なので、起動時の掃除は触れない。subject の source になるのは OPFS のファイルの File で、メモリの
   写しは捨てる。OPFS か Web Locks の無いブラウザ（結果をメモリに置くタブ）では何も残さず、そのタブでだけ使う。保存を断られた
   （容量が一杯）ときもそのタブでは使い、「Not kept for the next time: …」と出す。
6. **画面**（`ui/ExampleLine.vue`）。ボタンの横のその場の一文が、押すと何を送るか（「Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov」）、
   保存した写し（「Saved in this browser (2026-10-27, 41.8 MB)」と「Delete」）、取得中の進捗と「Cancel」を言う。設計の「Temporary
   storage の表示に載せる」（§5）は、その表示が stop-describing-ui で外れた（`ncbi_ui_mapping.md` §4）ので、Example の行に置く。
   「Delete」は写しと、それから作った subject の source を外す（File が読めなくなる）。他のタブにも BroadcastChannel
   `losat-web:examples` で知らせ（中身は表の id だけ）、そのタブも同じ source を外す。削除の途中に同じ組を押すと、読み込みは削除の
   後に始まる（消える写しを draft に渡さない）。保存中に Cancel しても写しは残るので、行は「Saved in this browser」になる。
   他の説明は「?」の後ろ。
7. **Example は両方の欄を置き換える**（EBI の「Use example sequence」と同じ）：query の欄に表の query、subject は組の写しだけ。
   他の source と範囲（subrange）は外れる。読み込みは 1 つずつで、別のボタン、Cancel、その組の無い mode への切り替え、失敗で
   終わる。query は組が届いたときの mode のもの（BLASTP から BLASTX へ移れば CDS）。
8. **e2e**（`tests/e2e/examples.spec.ts`）は WF-1 と同じく Playwright の route で実アドレスを stub し、配備の CSP のまま、ページから
   出た要求を数える。読み込み直しの後は要求が増えないこと、Delete で写しと source が消えることを OPFS のブラウザで見る。

## 10. WF-3 の実装（2026-11-03、Owner-delegated）

保守者の常設の指示（2026-09-29、2026-10-07 に再掲）により、推奨の案を採って記録する。作業記録は web-fetch-20261103/DECISIONS.md の 20〜22、
実測は同 NOTES.md「WF-3 facts」。

1. **競走は accession の取得だけ**。Example の組はどれも CORS のある host が 1 つ（NCBI の FTP、rest.uniprot.org）で、ENA と GenBase は
   レコードを返すが RefSeq Targeted Loci・pdbaa・assembly のファイルを持たない。Example はその host とその場の一文のまま。
2. **測るもの**：取得先と kind ごとに、取得が使う accession の入口（efetch、ENA browser API / UniProt REST、GenBase）へ固定の公開
   レコードを 1 つ求め、答えの全体が届くまでの時間。§2 の Range 要求は使わない（入口は API でファイルではなく、FTP の Range の時間は
   efetch の速さを言わない）。レコードは核酸が MH011443.1（TP53 の exon 5、192 B、3 か所とも）、タンパク質が大腸菌の thr operon leader
   peptide（NCBI と CNCB は NP_414542.1、UniProt は P0AD86。CNCB は UniProt の accession に空の本文を返す）。利用者の accession は
   送らない。要求は取得と同じ（`credentials: 'omit'`、referrer 無し、NCBI は 340 ms 空ける）に `cache: 'no-store'`（UniProt は
   `max-age=43200` を返し、ブラウザの写しの速さは何も測らない）。時間は要求が出た時から測る（NCBI の 340 ms の待ちは数えない、
   `FetchRequest.onStart`）。10 s 以内に答えない所は失敗と数える（ENA は 20〜25 s かかったことがある）。EBI の host は kind で
   違う（ENA / UniProt）ので kind ごとに測る。
   2026-11-03 の実測（日本から）：ENA は同じレコードに 20.1 s、25.6 s、1.3 s。アプリからの競走では 3 ブラウザとも NCBI が先（0.3〜0.6 s）。
3. **いつ走るか**：取得先の選択の既定を「Fastest」にする（保守者の 6「一番アクセスしやすい DB」、§2）。kind の最初の「Fetch」のときだけ
   走り、ページを開いた・mode を変えた・Run を押した・取得先を手で選んでいるときは走らない。1 つのタブで kind ごとに 1 回で、走っている
   間の別の Fetch は同じ競走を待つ。accession は、その全部を扱える取得先のうち最初に答えた所へ、その答えが来た時点で送る。扱えない所は
   問い合わせずに外す（§8 の 5 の EBI に加え、2026-11-03 の実測で：CNCB は UniProt の accession に空の本文を返す、GenBase の `C_…` は
   NCBI も ENA も持たない）。扱える所が 1 つなら競走しない（`C_…` は CNCB へ直に）。どこも全部は扱えない欄は何も送らず
   「No one site has all of these accessions: C_AA004835.1 only from CNCB, P12345 only from NCBI or EBI.」で止まる（WF-5 で §11 の
   組ごとのボタンに変わった）。他の答えは
   後ろで届いて順位を埋める。失敗した所は外し、扱える所のどれも答えなければその最初（ふつう NCBI）へ送る（失敗はその取得が名乗る）。どこも答えなかった競走は残さない（次の Fetch でもう一度、
   たとえばオフラインから戻った後）。順位はタブの間（メモリ）だけ保ち、読み込み直すとまた測る。手で選んだ取得先は §8 の 1 のとおり保つ
   （`web/app/src/application/source-race.ts`）。
4. **画面**：選択は「Fastest」「NCBI」「EBI」「CNCB」。競走の後は「Fastest (NCBI)」（その欄の accession を扱える最初の所）。その場の一文は
   競走の前は「Sends the accession to the fastest of NCBI, EBI and CNCB」（扱える所だけを並べる）、後は「Sends the accession to NCBI」。
   競走の間は「Finding the fastest site…」で、Run は「The query accession is still being fetched.」で止まる。選択は query と subject で
   共通なので、どちらかが取得している間は両方の欄で止める（取得中の欄の一文が送り先と食い違わないように）。「?」に競走が各取得先へ
   送るレコードの名前を足した。
5. **保存先**：Example が読むのは OPFS の写し（§9 の 5）のまま（確認の窓が出ない、オフラインで動く、3 ブラウザ）。保存した写しの横に
   「Save to folder…」（`showDirectoryPicker` のある Chromium / Edge。`mode: 'readwrite'`、`id: 'losat-examples'` で前回のフォルダから
   開く）を置き、選んだフォルダへ写しを `fileName`（展開した FASTA、subject の source と同じ名前なので、後で「Open FASTA files…」で
   開けば入力の名前も同じ）として OPFS のファイルから流して書く。同じ名前のファイルは置き換える（`createWritable` は閉じたときに
   入れ替わる。新しく作ったファイルは、写しに失敗したら消す：空のファイルを残さない）。API の無いブラウザ（Firefox、Safari、Android の Chromium）は「Download」（OPFS の File の object URL。メモリに
   写さない）。LOSAT はフォルダの handle を持たず、フォルダを自分では読まない：handle は site のデータに置くしかなく、それは OPFS の
   写しと一緒に消える（サイトデータの消去、容量不足の追い出し）。フォルダの写しが要るのはまさにそのときなので、「Open FASTA files…」で
   開く（3 ブラウザ）。メモリだけの経路（OPFS の無いタブ）には出さない。要求は出さない（ファイルはこのブラウザのもの）
   （`web/app/src/ports/download.ts` の `FileSaver`、`infra/browser/platform.ts`）。
6. **e2e**：競走は route の遅れで EBI を先にし（NCBI 1.5 s、CNCB 1.0 s）、最初の Fetch で 3 か所へ MH011443.1 だけ、accession は EBI へ、
   次の Fetch は競走しないことを数える（`tests/e2e/fetch.spec.ts`）。フォルダの picker は動かせないので、Chromium では OPFS のフォルダを
   返す picker に置き換えて書いた中身を読み、Firefox では download を受けて中身を比べる（`tests/e2e/examples.spec.ts`）。

## 11. WF-5 の実装：どこも全部は扱えない欄の組ごとのボタン（2026-11-03）

保守者（2026-11-03 21:5x）：どこも全部は扱えない欄は止めずに、取得先ごとのボタンを出す。押すとその組の accession だけを、その取得先
だけへ送る（1 回の操作で 1 か所。1 回の押しで欄を分けて送る #771 は退けた、web-fetch-20261103/DECISIONS.md の 26）。欄は全部の組が
届いてから一度に書く。細部は保守者の常設の指示（2026-09-29、2026-10-07 に再掲）により推奨の案を採った（同 DECISIONS.md の 28〜32）。

1. **組**：どこか 1 か所が全部を扱える欄はこれまでどおり 1 つの「Fetch」（§8、§10）。どこも全部は扱えない欄だけを、各 accession を
   扱える取得先の組（`servingSources`）で分ける：GenBase の `C_…` は CNCB だけ、UniProt のタンパク質は NCBI と EBI、RefSeq の核酸は
   EBI を除く。組の並びは欄の中の最初の accession の順（`domain/accession.ts` の `fetchGroups`）。CNCB はどの核酸の accession も
   扱うので、組ができるのはタンパク質の欄（`C_…` と UniProt の accession が混ざる欄）だけ。
2. **送り先**：取得先が 1 つの組はそこ。いくつかある組はその取得先だけの選択（「Fastest」は無い：ボタンは送る先を 1 つ名乗る）を持ち、
   初めは手で選んだ取得先（その組を扱えるなら）、無ければこのタブの競走の順の最初（走っている競走にも従う）、それも無ければ NCBI、EBI、
   CNCB の順で競走に落ちていない最初（`application/source-race.ts` の `groupSite`）。組は競走を始めない（読み込み、Run、組のボタンのどれでも）。
3. **手で選んだ取得先**：ボタンは選択に関わらず、どこも全部は扱えない欄に出る（手で選んだ所で止めると、その文が名指す所も別の
   accession で止まる）。1 か所が全部を扱えるのに選んだ所が扱えない欄は §8 のとおり要求の前に断る（利用者が取得先を名指した）。
4. **届いた組**：全部が届くまで draft のメモリ（`RoleDraft.fetched`、欄の行ごと）に持ち、最後の組が届いたら欄の行の順（組の順ではない）
   で一度に書く。欄の編集、mode を別の kind に変えること、欄全体の Fetch で捨てる。1 つの欄で取得は 1 組ずつ（他のボタンは待つ）。
   失敗した組は取得先と accession を名乗り（文に accession が無ければ前に付ける：「P12345: NCBI could not be reached (offline).」）、
   届いた組は残すので、そのボタンでその組だけをもう一度取れる。組の選択を変えると失敗の文は消える。§8 の大きさと数の上限は欄全体で
   数える。欄が書かれるまで Run は「Fetch the subject accessions first.」で、組を取っている間は「The subject accessions are still being
   fetched.」で止まる。
5. **画面**：行は「No one site has all of these accessions」（dev のまま）と組ごとの行：「[Fetch C_AA004835.1] from CNCB」、
   「[Fetch P12345 and 1 more] from [NCBI ▾]」（1 つの Fetch の行「Fetch from [Fastest ▾]」と同じ読み方。いくつかの accession の組は
   最初の accession を名乗る：「2 accessions」では 2 つの組が同じ名前になりうる。全部は title に並べる）。ボタンの名前（支援技術）は
   「Fetch C_AA004835.1 from CNCB」まで。届いた組は押せない「Fetched C_AA004835.1」で、取得先は届いた所のまま（後の競走や手の選択で
   変わらない。押した所は取得の間も行に残る）。取得中は §8 の「Fetching 1 of 2 from CNCB…」。その場の一文は足さず、「?」は §10 のまま
   （送るのは accession だけ、ここで選んだ所へ）。
6. **e2e**：blastp の subject に「C_AA004835.1 / P12345 / P99999」を入れ、CNCB の組は CNCB へ `C_` だけ、NCBI の組は NCBI へ 2 つ
   （P99999 で失敗、CNCB の組は残る）、EBI を選んで UniProt へ、最後に欄が欄の順で書かれること、競走の要求が出ないことを数える。
   別の欄の Fetch の競走で EBI が先になった後も、NCBI から届いた組は「Fetched P12345 from NCBI」のままで、新しい欄の組は EBI から
   始まる（`tests/e2e/fetch.spec.ts`）。
