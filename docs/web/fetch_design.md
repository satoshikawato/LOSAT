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
を取得先へ送る」に書き直す。ページを開いただけ・Run を押しただけでは通信しない。オフラインは今までどおり動く。

## 2. 取得先（ブラウザから直接。2026-11-03 実測、FETCH_NOTES.md）

| 取得先 | CORS | 用途 | 備考 |
|---|---|---|---|
| NCBI E-utilities `efetch` | `*` | accession → FASTA、`seq_start`/`seq_stop` で部分配列 | 鍵無し 3 req/s。`db=nuccore` / `protein` |
| ENA browser API `/ena/browser/api/fasta/ACC?range=a-b` | `*` | INSDC accession → FASTA（範囲も。題に `Location:a..b`） | RefSeq（NR_, NM_…）は 400 |
| UniProt REST `/uniprotkb/ACC.fasta` | `*` | タンパク質 accession | `stream?query=` は実測で 500（構文を要確認） |
| NCBI FTP `refseq/TargetedLoci/*` | `*` | Example の subject セット | 下の表 |
| NCBI FTP `blast/db/FASTA/pdbaa.gz` | `*` | タンパク質の Example 候補 | 36.8 MB gz（境界） |
| DDBJ getentry | 無し | — | CORS ヘッダ無し（curl には本文も空）。DDBJ Search の API は CORS を返すが、索引は BioProject / BioSample / DRA / GEA などで塩基配列のレコードを持たない（`/search/api/entries/nucleotide/…` は「Unknown database type」）。公開ファイル（`public/ddbj_database/16S/16S.fasta.gz`、222 MB）も CORS 無し。DDBJ が getentry に `Access-Control-Allow-Origin: *` を付ければ取得先に足せる |
| CNCB GenBase | header は返す | — | API の URL 未確認 |

DDBJ は DDBJ 側のヘッダ追加（依頼）があるまで取得先に入れない。それまで DDBJ 利用者は NCBI / EBI から取るか、DDBJ から手で落として「Open FASTA files…」で開く。
「一番アクセスしやすい所」= 候補（NCBI、EBI）へ同時に Range 要求を出し、先に応えた方を取る。選んだ先はそのセッション中は固定し、
取得元の名前を画面に出す（手で NCBI / EBI を選べる）。

## 3. Example セット（mode ごと、各 50 MB 以下、固定 URL、SHA-256 を表に持つ）

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
  OPFS block store（`infra/data/opfs-block-store.ts`）に **キャッシュ**（URL + SHA-256 を鍵に）。次回は取得せずに読む。
  OPFS は origin 専用で、サイトデータの消去で消える。使用量と「Delete」は既存の Temporary storage の表示に載せる。
- 「任意のフォルダに保存」= File System Access API（`showDirectoryPicker`）。Chromium / Edge だけ。Firefox / Safari は
  「Download」で落として次回「Open FASTA files…」から開く（ファイル参照なので大きくても速い）。
- 展開：`.gz` は `DecompressionStream('gzip')`（3 ブラウザとも可）で OPFS に書く。

## 6. 進め方（PR の切り方）

1. **WF-0 設計**：このメモを `docs/web/` に入れ、設計書 §13 と `web/AGENTS.md` 規則 4、`ncbi_ui_mapping.md` §1 の accession の行を
   書き換える。CSP `connect-src` に取得先を足す。
2. **WF-1 accession**：欄の accession 判定、「Fetch from NCBI / EBI」、efetch と ENA の port（`ports/fetch.ts`）と fake、
   e2e はローカルの stub サーバ。
3. **WF-2 Example**：静的表、OPFS キャッシュ、gz 展開、進捗と失敗の表示、mode ごとの 1 行。
4. **WF-3 保存先と mirror**：「Save to folder…」（Chromium）、NCBI / EBI の競走、取得元の選択。

## 7. 決定（保守者、2026-11-03 15:3x）

1. **方針 = GO**：§1 の約束の書き換えと、WF-0（この文書と設計書 §13、`web/AGENTS.md` 規則 4、対応表 §1 の accession の行、CSP の
   `connect-src`）→ WF-1 accession → WF-2 Example とキャッシュ → WF-3 保存先と mirror の順。
2. **タンパク質の Example = 両方**：既定は UniProt の reviewed の分類群セット（数 MB。REST の stream 構文を再確認してから表に入れる）、
   pdbaa.gz（36.8 MB）も表に入れて選べる。
3. **取得先 = NCBI と EBI だけ**：CNCB は API を確認してから取得先の表に足せる形にしておく。DDBJ は CORS が無いので入れない。
4. **文言 = その場の一文だけ**：取得の行に「Sends the accession to NCBI」「Downloads 8 MB from ftp.ncbi.nlm.nih.gov」。約束の文は
   "?" の後ろ。見出しは「LOSAT」と一文のまま。
