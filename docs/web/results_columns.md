# LOSAT Web 結果画面の列定義表

設計書 §11.2 の対応表のうち、結果画面（W4、Session S13）の列とフィルターを定める。計画 §4.4・§5.7 と `PD-LOSAT-WEB-APP-BOUNDARY` の不変条件 2.3 に従い、表に出す BLAST の値は outfmt 6 の該当行（HSP レコードの `out6` の範囲）か outfmt 0 の原文から取り、並べ替えとフィルターにはエンジンの原値（HSP レコードの数値）を使う。`web/` で BLAST の値を計算・整形しない。

- 作成：2026-10-06（S13）。状態の「採用」は S13 で画面に出すもの、「採用・エンジン待ち」は採用を決めたがエンジンが値を出すまで画面に出さないもの（設計書 §11.2「未対応項目は通常画面に出さず」）。
- 採否は保守者の確認待ち（S13 は推奨の案で進めた。[W4 のゲート記録](../evidence/losat_web_w4/README.md)の判断 1）。
- 単位の nt / aa は、その役割のレコードの種類（`src/domain/programs.ts` の `query` / `subject`）。outfmt 6 の座標は、TBLASTN の subject、TBLASTX の両方、BLASTX の query で核酸の座標である（レコードの長さと同じ単位）。

## Subject 一覧（選んだ query の subject ごとに 1 行）

| 列 | 単位 | program | 値の出どころ | 集約の単位 | 並べ替え | 状態 |
|---|---|---|---|---|---|---|
| # | — | 全部 | エンジンの順（その subject の最初の HSP の `index`） | subject | `index` | 採用 |
| Subject | — | 全部 | outfmt 6 の `sseqid`（2 列目）：その subject の最初の HSP の行 | subject の最初の HSP | — | 採用 |
| Description | — | 全部 | outfmt 0 の subject の見出し（`out0_subject` の範囲。`> ` から `Length=` の前まで、NCBI の折り返しのまま）。outfmt 0 にアラインメントの無い subject（`out0_subject` が null）は「not in outfmt 0」 | subject | — | 採用 |
| Length | nt / aa | 全部 | Run の subject のレコードの長さ（エンジンの `register` が読んだもの、RunSnapshot） | subject | 原値 | 採用 |
| Score (bits) | bits | 全部 | outfmt 6 の `bitscore`（12 列目）：その subject の最初の HSP。outfmt 0 の説明表の「Score (Bits)」と同じ HSP の同じ文字列（NCBI `showdefline.cpp:698-700` の `x_GetScoreInfo`（`:1483-1515`）が subject の最初のアラインメントから取り、outfmt 6 と同じ `CAlignFormatUtil::GetScoreString` で書く） | subject の最初の HSP | `bit_score` | 採用 |
| E value | — | 全部 | outfmt 6 の `evalue`（11 列目）：同上（説明表の「E Value」） | subject の最初の HSP | `e_value` | 採用 |
| HSPs | 個 | 全部 | この query のこの subject の HSP レコードの数（アプリの数え上げ。BLAST の値ではない） | subject | 個数 | 採用 |
| Max score | bits | 全部 | NCBI `CAlignFormatUtil::GetSeqAlignSetCalcParams`（`align_format_util.cpp:4247-4316`）の `highest_bits` を `GetScoreString` で書いたもの | subject の全 HSP | 原値 | 採用・エンジン待ち（S13+） |
| Total score | bits | 全部 | 同 `total_bits`（全 HSP の bit score の和）の `total_bit_score_buf`（`showdefline.cpp:1518-1560` の `x_GetScoreInfoForTable`） | subject の全 HSP | 原値 | 採用・エンジン待ち（S13+） |
| Query cover | % | 全部 | 同 `percent_coverage`（`GetSeqAlignCoverageParams` の query 側の被覆長 ÷ query 長、四捨五入、100 を超えれば 99） | subject の全 HSP | 原値 | 採用・エンジン待ち（S13+） |
| Per. ident | % | 全部 | 同 `percent_identity`（bit score が最も高い HSP の identity、WB-1175） | bit score が最も高い HSP | 原値 | 採用・エンジン待ち（S13+） |
| E value（表形式） | — | 全部 | 同 `lowest_evalue`（bit score が最も高い HSP の E 値） | bit score が最も高い HSP | 原値 | 採用・エンジン待ち（S13+）。エンジンが出したら、上の「Score (bits)」「E value」を NCBI の表形式（BLAST+ の `-sorthits` の説明表、NCBI Web BLAST の「Max Score / Total Score / Query Cover / E value / Per. Ident / Acc. Len」）の列と並びに置き換える |
| Scientific name、Common name、Taxid、Accession | — | — | NCBI Web の BLAST DB の taxonomy と accession。ローカルの FASTA には無い | — | — | 不採用（情報が無い。設計書 §11.2「存在しない情報を推測して列を埋めない」） |

BLAST+ は、`-sorthits` を与えたときだけ outfmt 0 の説明表に Total Score・Query Coverage・Percent Ident を出す（`blast_format.cpp:518-522` が `eShowTotalScore`・`eShowQueryCoverage`・`eShowPercentIdent` を立てる）。NCBI Web BLAST の説明表も同じ `x_GetScoreInfoForTable` の値である。したがってこれらは NCBI のソースに定義があり、BLAST+ の CLI を比較の oracle にできる（outfmt 6 の `qcovs` も同じ被覆）。移植はエンジン側のセッション S13+（段階 E2j）で行い、ABI v2 で subject ごとの値を NCBI の文字列のまま渡す。TS では計算しない。

## HSP 一覧（選んだ subject の HSP ごとに 1 行）

| 列 | 単位 | program | 値の出どころ | 並べ替え | 状態 |
|---|---|---|---|---|---|
| # | — | 全部 | その query の中の順（`rank` + 1） | `rank` | 採用 |
| Bit score | bits | 全部 | outfmt 6 の `bitscore` | `bit_score` | 採用 |
| E value | — | 全部 | outfmt 6 の `evalue` | `e_value` | 採用 |
| Identity | % | 全部 | outfmt 6 の `pident` | — | 採用 |
| Length | 列 | 全部 | outfmt 6 の `length`（アラインメントの長さ） | — | 採用 |
| Mismatches | 個 | 全部 | outfmt 6 の `mismatch` | — | 採用 |
| Gap opens | 個 | 全部 | outfmt 6 の `gapopen` | — | 採用 |
| Query start / end | nt / aa | 全部 | outfmt 6 の `qstart` / `qend` | `q_start` | 採用 |
| Subject start / end | nt / aa | 全部 | outfmt 6 の `sstart` / `send` | `s_start` | 採用 |
| Frames | — | TBLASTN・TBLASTX（BLASTX は SX） | HSP レコードの `query_frame` / `subject_frame`（整数、符号付きで示す） | — | 採用 |
| Orientation | — | 全部 | HSP レコードの座標（start > end が逆向き、ABI v2 §8）と frame の符号。BLASTN の 1 文字の HSP（start = end）は座標から分からないので「not in the HSP record」とし、詳細の outfmt 0 の節（`Strand=` の行）を見るよう示す | — | 採用。BLASTN の 1 文字の HSP の鎖は S13+ が HSP レコードに足す（判断 2） |
| In outfmt 0 | — | 全部 | `out0` が null かどうか | — | 採用 |

## Query の一覧

| 列 | 単位 | 値の出どころ | 状態 |
|---|---|---|---|
| # | — | レコードの番号（`q_idx` + 1） | 採用 |
| Query | — | Run の query のレコードの ID（エンジンの `register`、RunSnapshot） | 採用 |
| Length | nt / aa | 同じレコードの長さ | 採用 |
| Subjects / HSPs | 個 | その query の subject と HSP レコードの数（アプリの数え上げ） | 採用 |

## 表示用のフィルター（ViewState）

ViewState は表示だけを変え、検索をやり直さない（設計書 §11.2）。互換出力（outfmt 0/6/7）の書き出しには反映しない（計画 §5.8）。

| フィルター | 比べる値 | 状態 |
|---|---|---|
| E value ≤ x | HSP レコードの `e_value`（原値） | 採用 |
| Bit score ≥ x | HSP レコードの `bit_score`（原値） | 採用 |
| Subject ID を含む | Run の subject のレコードの ID と、outfmt 6 の `sseqid`（原文） | 採用 |
| Query ID を含む、Queries with hits only | Run の query のレコードの ID、その query の HSP レコードの有無 | 採用 |
| Identity ≥ x | outfmt 6 の `pident` は丸めた文字列で、レコードに原値が無い | 不採用（原値が無い） |

## 状態の区別（設計書 §11.2）

| 状態 | 条件 | 表示 |
|---|---|---|
| 表示 0 件 | HSP はあるが、フィルターがすべてを隠した | 「No HSPs of this query match the view filters」と、隠した数、フィルターを外す操作 |
| ヒット無し | Run は完了し、その query の HSP レコードが無い | 「No hits for this query.」（outfmt 0 の「No hits found」と同じ事実。値は作らない） |
| 未完了・取消・失敗 | Run の状態 | 結果を出さず、状態と失敗の文言（RunRecord の `error`） |
| 結果数の上限に達した可能性 | その query の subject の数が hit list の上限（argv の `-max_target_seqs`、無ければ `describe` の help の「default: N」）と同じ。または subject の HSP の数が `-max_hsps` と同じ | 上限に達したことと、より大きい `-max_target_seqs` で検索すれば示されることだけを示す。取れなかったヒットがあるとは言わない |
| outfmt 0 にアラインメントが無い | `out0` が null（outfmt 0 は query ごとに先頭の `num_alignments` 個の subject だけアラインメントを示す。既定は 250（`format_flags.cpp:221`）、`-max_target_seqs` を与えるとその値（`blast_args.cpp:2909-2928`）） | 詳細に理由と、outfmt 6 の行を示す |
