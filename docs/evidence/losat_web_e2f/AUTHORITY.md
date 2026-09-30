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
| `blastn/query_split.rs` | `CQuerySplitter` の塊の数・範囲・部分・context・補正と、mask の制限（§D） |
| `run.rs` の `search_query_chunks`・`merge_query_chunk`・`merge_prelim_hit_lists` | 分割する batch の塊ごとの予備の段階（`search_query_batch` の `BatchStage::ChunkPrelim`）と、`BlastHSPStreamMerge`・`Blast_HitListMerge`・`Blast_HSPListsMerge` の移植（§D） |

S07+ で batch に依存するために置いた次のものを、本物の batch で置き換えた。

- 最初の batch だけを推定した `first_query_batch` と、検索されなかった batch を推定した `unsearched_queries`（最初の batch の後の無効な query の outfmt 0・7 の拒否、E2c の §F・§G）
- 表を超える gap で、context の gapped X-drop が違い、query が最初の batch に収まらない場合の拒否（`gap_x_dropoffs`、E2c の §D・§G）
- Karlin-Altschul の表の失敗の文言の繰り返しの回数（batch の query の数）

残る明示的な拒否は 1 つ：得点の表が無い得点で、最初の batch が無効な query だけのとき（NCBI は、その batch の結果と警告を書いた後に、次の batch で誤りを出す）。

## D. query の分割

NCBI は、batch の query の総文字数を「chunk の大きさ − 重なり 100」で割った商が 2 以上のとき（`-task blastn` で 1,999,800、megablast で 9,999,800 文字以上）、batch を query の塊に分けて予備の段階を塊ごとに行い、HSP を batch の context に移して合わせ、traceback は batch 全体で行う。LOSAT はこれを移植した（`query_split.rs`、`run.rs` の `search_query_chunks`・`merge_query_chunk`・`merge_prelim_hit_lists`）。

| NCBI | 振る舞い |
|---|---|
| `split_query_aux_priv.cpp:105-145`（`SplitQuery_CalculateNumChunks`） | 塊の数は総文字数 /（chunk − 100）。2 以上なら、塊の大きさを（総文字数 +（塊の数 − 1）× 100）/ 塊の数とし、塊の数 < 大きさ − 100 なら 1 を足す |
| `split_query_cxx.cpp:135-180` | 塊の範囲：始まりを「大きさ − 100」ずつ進め、最後の塊は総文字数まで |
| `split_query_cxx.cpp:214-293`、`196-210` | query の範囲は前の query の終わりから続く（`COpenRange`）。塊と交わる query ごとに、塊の中の部分を作る |
| `split_query_cxx.cpp:392-410` | 塊の context は、部分ごとに batch の 2 本の鎖の context |
| `split_query_cxx.cpp:580-665`、`split_query_aux_priv.cpp:251-284` | 塊の context の補正：部分が batch の context で始まる位置（minus 鎖は query の終わりから）。式は `size_t` で、部分の長さと重なりの小さい方を引く。部分が続く限り正しい位置になる（`query_split.rs` の試験） |
| `split_query_cxx.cpp:278-281`、`seqlocinfo.cpp:98-124`、`blast_setup.c:1030-1062` | 塊の mask：batch の設定のときに query の vector に書かれた mask（DUST と小文字）を部分に制限する。`RestrictToSeqInt` は部分の最後の文字を範囲から外し、交わりの終わりを 1 文字伸ばす。そのため部分の最後の文字から始まる mask は落ち、ほかの mask は 1 文字長くなる（部分の終わりまで） |
| `dust_filter.cpp:166-186`、`blast_objmgr_tools.cpp:176-228` | 塊の部分にも DUST を行い、制限した mask と合わせる |
| `split_query_aux_priv.cpp:150-182`、`blast_setup.c:670-697,774-850`、`query_data.cpp:52-58` | 塊の探索空間：options に batch の context の探索空間を入れ、塊の context は、batch の context でなく塊の中の番号の値を使う。値が 0 なら塊の長さから計算する。batch の query data は cache されるので、計算は batch の設定の無効な context を見て 0 にする |
| `prelim_stage.cpp:225-300`、`split_query_aux_priv.cpp:185-210` | 塊ごとに query の塊・lookup table・Karlin-Altschul の block（塊の組成）・X-drop・対角線の表・診断を作って予備の段階を行う（`-subject_besthit` と予備の e-value の刈り込みは塊の値で行う）。塊の警告は出さず、有効な context の無い塊の誤りは無視する |
| `blast_hspstream.c:399-536`、`hspfilter_collector.c:97-160` | 塊の HSP を query ごとの一覧に分け、context と補正で batch に移し、補正を分割点にして `Blast_HitListMerge` で合わせ、最後に得点で並べる |
| `blast_hits.c:2119-2219`、`2757-2804`、`2809-3034` | 分割点のどれかが正なら `Blast_HSPListsMerge`（query の分割の枝。重なりの帯の判定は鎖で向きが変わり、帯の HSP を前に入れ替えてから、対角線の差が 10 未満のものを `s_BlastMergeTwoHSPs` で合わせる）、そうでなければ `Blast_HSPListAppend` |
| `local_blast.cpp:301-310` | 分割した batch の `GetNumExtensions` は batch の診断で 0（塊の診断は別）。次の batch の大きさは初期の hit を 0 として決まる |

分割した batch には lookup table が無い（`blast_aux_priv.cpp:206-207`。塊ごとに作る）。`Blast_HSPListsMerge` の帯の HSP の入れ替えは、subject の分割の枝にも入れた（S07+ までは入れ替えず、得点が同じ HSP の順が違いえた）。

確かめたこと（NCBI BLAST+ 2.17.0 とのバイト比較。`split_check/`）：EDL933（5,528,445 文字）を query にした `-task blastn`（5 つの塊）の、塊の境の前後の 63 の窓（S07++ の最初の比較で 1 つが違った窓を含む）。EDL933 と Sakai をつないだ 1 つの query（11,027,023 文字）の megablast（2 つの塊）の、塊の境と 2 つのゲノムの境の窓。2 つの query にした EDL933 と Sakai（NCBI は 2 つの batch にする）の `-task blastn` の、2 つ目の batch の塊の境の窓。無効な query を分割する batch の前と後に置いた outfmt 0・6・7。短い query の後の長い query（2 つ目の塊の context が短い query の探索空間を使う）。小文字の mask の 1 文字の伸びと、塊の最後の文字から始まる mask の脱落（`-word_size 8`）。複数の subject と `-max_target_seqs`・`-subject_besthit`・`-num_threads 4`。分割した batch の後の batch。塊の境の周りの乱数の窓（分岐、挿入と欠失、逆向き）。NCBI に `CHUNK_SIZE=20000000` を与えて分割させない出力とも比べ、分割で NCBI の出力が変わる case（探索空間、mask の 2 つの細部、境の HSP）を LOSAT が再現することを確かめた。

## E. 予備の段階の hit list（S07++ の独立監査の第 1 回）

| NCBI | 振る舞い |
|---|---|
| `blast_engine.c:1409-1554`、`hspfilter_collector.c:83-161` | 予備の段階は subject を順に検索し、subject の HSP の一覧（予備の e-value の刈り込みの後）を collector に書く。collector は一覧を query ごとに分け（順を保つ）、query の hit list に入れる |
| `blast_hits.c:44-68`、`3243-3300` | hit list の大きさは `prelim_hitlist_size`（`MIN(MAX(2 × hitlist_size, 10), hitlist_size + 50)`。既定で 550、`-max_target_seqs` 1〜5 で 10）。あふれると heap にし、予備の e-value（一覧の最良）・最初の HSP の得点・subject の番号で最も悪い一覧を捨てる（`Blast_HitListUpdate`、`s_EvalueCompareHSPLists`）。heap にするとき、一覧の HSP を e-value で並べる |
| `blast_hspstream.c:133-206`、`blast_traceback.c:1500-1707` | traceback は残った一覧だけを、subject ごとに読み、query ごとの hit list（大きさ `hitlist_size`、`Blast_HSPResultsInsertHSPList`）に入れる |
| `blast_hits.c:2119-2217` | 分割した batch では、塊ごとの collector と、合わせるときの新しい hit list（`Blast_HitListMerge` の `Blast_HitListNew(hitlist1->hsplist_max)`）が、それぞれ大きさを守る。合わせた HSP の e-value は前の HSP のまま（`s_BlastMergeTwoHSPs`） |

S07+ までの LOSAT は、subject ごとに予備の段階と traceback を続けて行い、hit list の大きさ（`prelim_hitlist_size`）を traceback の後の e-value に適用していた。そのため、予備の hit list からあふれる subject があると結果が違った（分割しない batch でも。独立監査は、既定の設定で 560 の subject、`-max_target_seqs` 1 と 3、分割しない 300,000 文字の query で再現した）。LOSAT は、全 subject の予備の段階の後に collector を移植し（`run.rs` の `collect_prelim_hit_lists`。hit list は `hsp.rs` の `HitList` を予備の HSP の一覧にも使う）、残った一覧だけを traceback する（`search_subjects`）。予備の HSP は e-value を持ち（`PrelimHit::prelim_evalue`）、分割した batch は `merge_prelim_hit_list`（`Blast_HitListMerge`）で合わせる。traceback の後の hit list の大きさは NCBI と同じく `hitlist_size` にした。
