# E2f（S07++）：BLASTN の query の batch — NCBI の経路

NCBI のソースは固定 commit `598d8ae6`（`/mnt/c/Users/genom/GitHub/ncbi-blast/c++`）、振る舞いの確認は NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin/blastn`）で行った。LOSAT の移植箇所には、それぞれの直上に NCBI のファイル・行と断片を書いてある。

## A. NCBI の batch の経路

| NCBI | 振る舞い |
|---|---|
| `blastn_app.cpp:261-300` | `CBatchSizeMixer mixer(SplitQuery_GetChunkSize(program) - 1000)`。`BATCH_SIZE` が無ければ、subject の総文字数 `total_len` が正のとき `SetTargetHits(total_len / 3000)`、最初の batch の大きさは `mixer.GetBatchSize()`。batch ごとに `CLocalBlast` を作って検索し、`input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()))` で次の batch の大きさを決める。結果は batch ごとに query の順に書く |
| `local_blast.cpp:54-73` | chunk の大きさは blastn 1000000、megablast 5000000（環境変数 `CHUNK_SIZE` が無いとき） |
| `blast_app_util.hpp:57-83`、`blast_app_util.cpp:65-81` | `CBatchSizeMixer`：目標の hit 数は既定 1000000（最小 1000000）、最初の batch は目標の 1/200（普通は 5000）。`GetBatchSize(hits)` は `ratio = (hits + 1) / m_BatchSize`、前の比があれば 0.3 と 0.7 で混ぜ、`m_BatchSize = (Int4)(m_TargetHits / m_Ratio)`。上限（chunk − 1000）を超えれば上限、100 未満なら 100 にし、どちらのときも比を −1 に戻す |
| `blast_input.cpp:138-170`（`GetNextSeqBatch`） | 読んだ文字数が batch の大きさより少ない間、query を足す（batch の大きさに達した query は、その batch に入る） |
| `local_blast.cpp:301-310`（`GetNumExtensions`） | batch の診断の `ungapped_stat->good_init_extends`。検索しなかった batch（有効な context が無い、`local_blast.cpp:177-180`）では 0 |
| `blast_engine.c:491-500`、`na_ungapped.c:1688-1689`、`blast_diagnostics.c:102-115` | subject の chunk ごとに初期の hit の一覧を空にし、word を探した後に、その chunk で保存した初期の hit の数（`init_hitlist->total`）を `good_init_extends` に足す（lookup の hit が 0 の chunk は足さない。足す数も 0） |

`(Int4)` の変換は、オラクルの x86-64 の実行ファイルでは `int` の外の値を `INT_MIN` にする（E2c の §R）。そのため、batch の初期の hit が 0 か 1 のとき（`1000000 / ((hits + 1) / 5000)` が 2^31 を超える）、次の batch は最小の 100 文字になる。

## B. batch ごとに決まるもの

NCBI は batch ごとに query の塊（2 つの鎖と区切りの文字）を作り、その塊から次を決める。

- lookup table の種類と幅（`BlastChooseNaLookupTable` の見積もり）、word の延長の方法（E2c の §Q）、scan の歩幅
- 1 つの query の対角線の表の大きさ（`query->length`、E2c の §J）と、hash か配列か
- 近似の ungapped 伸長が止まる query の塊の端（`na_ungapped.c:164,286,317`。E2c の §K）
- gapped の X-drop（batch の有効な context の gapped Lambda の最小値、`blast_parameters.c:92-116,455-463`）
- 予備の段階の `-subject_besthit`（chunk ごとに合わせた、batch の全 query の HSP の一覧、E2c の §K）
- 検索するかどうか（有効な context が無い batch は検索しない）

query ごとの値（ungapped と gapped の block、長さの補正、探索空間、cutoff、e-value、予備の hit list の大きさ）は batch によらない。

## C. LOSAT の対応

| LOSAT | 内容 |
|---|---|
| `blastinput/query_batch.rs` の `BatchSizeMixer`、`next_query_batch_end` | `CBatchSizeMixer` と `GetNextSeqBatch` の移植。`(Int4)` は `core/blast_util.rs` の `ncbi_int4_from_double` |
| `run.rs` の `run_in_pool` | subject を 1 度だけ読み符号化し（`PreparedSubjects`）、batch を NCBI と同じく決めて `search_query_batch` を呼ぶ。batch の HSP の query の番号を全体の番号に直し、全 query の結果を 1 度に書く |
| `run.rs` の `search_query_batch` | batch の query の塊、lookup table、対角線の表、X-drop、`-subject_besthit` で検索する（S07+ までの全 query の検索と同じ処理を batch に対して行う）。subject の chunk ごとの初期の hit の数を足す |

S07+ で batch に依存するために置いた次のものを、本物の batch で置き換えた。

- 最初の batch だけを推定した `first_query_batch` と、検索されなかった batch を推定した `unsearched_queries`（最初の batch の後の無効な query の outfmt 0・7 の拒否、E2c の §F・§G）
- 表を超える gap で、context の gapped X-drop が違い、query が最初の batch に収まらない場合の拒否（`gap_x_dropoffs`、E2c の §D・§G）
- Karlin-Altschul の表の失敗の文言の繰り返しの回数（batch の query の数）

残る明示的な拒否は 1 つ：得点の表が無い得点で、最初の batch が無効な query だけのとき（NCBI は、その batch の結果と警告を書いた後に、次の batch で誤りを出す）。

## D. query の分割

NCBI は、batch の query の総文字数が 2 ×（chunk − 100）以上のとき、batch を重なり 100 の塊に分けて予備の段階を行い、HSP を合わせる（`split_query_cxx.cpp:50-62`、`split_query_aux_priv.cpp:51-145`、`prelim_stage.cpp:230-260`）。LOSAT は分割しない。EDL933（5,528,445 文字）を query にした `-task blastn`（NCBI は 5 つの塊に分ける。塊の大きさ 1,105,770、重なり 100）で、塊の境の前後に置いた 63 の subject の窓（EDL933 と Sakai、300〜20000 文字）を比べると、62 は NCBI とバイト一致し、1 つ（境でちょうど終わる 3000 文字の窓）で LOSAT に 17 文字の HSP が 1 本多かった。NCBI に `CHUNK_SIZE=10000000` を与えて分割させないと、LOSAT と同じになる。分割の移植（塊ごとの query の塊と lookup table、塊ごとの予備の段階、`BlastHSPStreamMerge` での重なりの HSP の合わせ方）は大きいので、NCBI が分割する batch（query の片方の鎖の合計が 2 ×（chunk − 100）以上、`-task blastn` で 1,999,800、megablast で 9,999,800 文字以上）は明示的に拒否する（`run.rs` の `run_in_pool`）。1,999,799 文字の query は NCBI とバイト一致し、1,999,800 文字は拒否する（`check_inputs.py` の `s07pp.query_unsplit`・`s07pp.query_split`）。Gate A の EDL933 × Sakai の megablast（5.5 Mb）は分割の対象ではない。分割の移植は計画の未決事項にした。
