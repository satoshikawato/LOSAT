# E2c（S07+）：BLASTN の得点のオプションと入力の読み方 — NCBI の経路

NCBI のソースは固定 commit `598d8ae6`（`/mnt/c/Users/genom/GitHub/ncbi-blast/c++`）、振る舞いの確認は NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin/blastn`）で行った。LOSAT の移植箇所には、それぞれの直上に NCBI のファイル・行と断片を書いてある。この文書は、どの NCBI の経路が LOSAT のどこに当たるかと、S07+ の前に何が違っていたかをまとめる。

## A. 得点のオプションの値

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_nucl_options.cpp:137-150,198-221`（`SetLookupTableDefaults`・`SetMBLookupTableDefaults`・`SetScoringOptionsDefaults`・`SetMBScoringOptionsDefaults`）、`blast_options.h:67-68,86-95` | task の既定値：blastn は word size 11、reward 2、penalty −3、gap 5/2。megablast は 28、1、−2、0/0 | `coordination.rs` の `task_defaults` |
| `blast_args.cpp:166-170,175-182,262-271,288-300,648-660,674-679` | `-word_size`・`-gapopen`・`-gapextend`・`-reward`・`-penalty` は省略できるキーで、与えたときだけ task の既定値を置き換える。reward は 0 以上、penalty は 0 以下、gap は制約なし | `args.rs` のこれらを `Option` にした。reward は 0 以上、penalty は 0 以下（NCBI の引数の制約）、word size は 4 以上で、100 以下は §B の検査 |
| `blast_options_local_priv.hpp:1628-1643`、`blast_options.h:465-466` | reward と penalty は 16 ビット（`Int2`）に入れるので、範囲の外の値は検査の前に折り返す（例：`-penalty -40000` は 25536 になり「penalty must be negative」、65537/−65538 は 1/−2） | `coordination.rs` の `determine_scoring_params` が同じく 16 ビットに折り返す |

**S07+ の前の差（根本原因 1）：** LOSAT は clap の既定値を megablast の値にし、task が blastn のときに「既定値と同じ値」を blastn の既定値に置き換えていた。そのため、`-task blastn -reward 1` のように megablast の既定値と同じ値を明示すると、その値が無視された（S06 の sweep の reward 1 の 16 件の DIFF の原因）。

## B. オプションの検査（検索の前）

| NCBI | 振る舞い |
|---|---|
| `blast_options.c:881-906`（`BlastScoringOptionsValidate`） | penalty が 0 以上なら「BLASTN penalty must be negative」、gap open が正で gap extend が 0 なら「BLASTN gap extension penalty cannot be 0」 |
| `blast_options.c:1326-1333`（`LookupTableOptionsValidate`） | word size が 100 を超えると「Word-size must be less than or equal to 100」 |
| `blast_options.c:1518-1523`（`BlastHitSavingOptionsValidate`） | e-value が 0 以下なら「expect value or cutoff score must be greater than zero」 |
| `blast_options.c:1699-1711`（`s_BlastExtensionScoringOptionsValidate`） | gap が 0/0 で貪欲な伸長でない（blastn の task）なら「Greedy extension must be used if gap existence and extension options are zero」 |
| `blast_options.c:1759-1776`（`BLAST_ValidateOptions`） | 上の検査を、得点、lookup table、hit saving、伸長と得点の順に調べる |
| `blast_args.cpp:3631-3633`、`blastn_args.cpp:63-70`、`blast_args.cpp:2553-2557` | 引数の処理は、引数の組ごとに順に行う。`CBlastDatabaseArgs`（`-subject`）は `CStdCmdLineArgs`（`-query`・`-out`）と形式の引数より前で、ここで subject を読む（§E。S07+ の第 3 回の監査） |
| `blast_args.cpp:3474-3480` | `-out` のファイルは検査より前に開く（検査で止まっても空のファイルが残る） |
| `blast_args.cpp:3609-3612`、`blast_app_util.hpp:172-175` | 検査の失敗は `CInputException` になり、`BLAST query/options error: <msg>` と `Please refer to the BLAST+ user manual.` を出して終了コード 1 |
| `blast_args.cpp:2800-2803,2975-2977` | 出力形式の解析と「Examining 5 or more matches is recommended」の警告は、検査より前（引数の処理の時点） |
| `ncbiargs.cpp:464`、`ncbistr.cpp:1312-1318,1328-1339` | `-evalue` は `NStr::StringToDouble` で読む。最初の文字が数字、小数点、符号でなければ引数の誤り（`inf`、`nan`、` 1`）。それ以外は `strtod` でも読むので、`+inf`・`-inf`・`+nan`・`0x10`（16）も受け付ける（2.17.0 で、`+inf` と `+nan` は検索し、`0x10` は `-evalue 16` と同じ出力） |

LOSAT：`scoring.rs` の `check_scoring_options`。CLI の `run` は、NCBI と同じく、入力のファイルを 1 度ずつ開き、subject を読んでレコードが無ければ engine error で止まり、`-out` のファイルを作り、形式の解析、警告、検査（`process_options`）をこの順に、query を読む前に行う（`run_local` も、subject のレコードの検査を `process_options` より前に行う）。CLI の終了コードは、BLASTX と共有する `cli.rs` の `NativeError`（S07+ で `blastx/native.rs` から移した）で NCBI と同じにした。LOSAT の上限（§G）は NCBI の検査の後に調べる。ただし、得点の範囲の上限は、NCBI の Karlin-Altschul の表の失敗（§C。その前に行う計算の範囲を LOSAT が確かめていない）より前で、greedy の gap の上限は表の検査の後（エンジン）で調べる。環境変数 `BATCH_SIZE`・`CHUNK_SIZE` は、NCBI の検査の後で拒否する。

## C. Karlin-Altschul の表

| NCBI | 振る舞い |
|---|---|
| `ncbi_math.c:405-419`（`BLAST_Gcd`） | reward と penalty の最大公約数 |
| `blast_stat.c:3237-3371`（`s_GetNuclValuesArray`） | reward と penalty を最大公約数で割ってから、12 の表のどれかを選ぶ。表の最初の行が 0/0 なら、それは線形（megablast）の行として分ける（`s_SplitArrayOf8`）。表の外の組は「Substitution scores %d and %d are not supported」（割った後の値）。`round_down` は割った後の 2/−7、2/−5、2/−3、3/−4 |
| `blast_stat.c:3186-3219`（`s_AdjustGapParametersByGcd`） | 割った場合、表の gap の値と上限に約数を掛け、Lambda と alpha を約数で割る |
| `blast_stat.c:3845-3947`（`Blast_KarlinBlkNuclGappedCalc`） | gap 0/0 で線形の行があればその行、なければ表の行を探す。無ければ、gap が上限以上なら ungapped block を写す。それ以外は「Gap existence and extension values … are not supported for substitution scores …」に、表の組、上限、「Any values more stringent than …」を続けた文言（与えた値の reward と penalty） |
| `blast_stat.c:3965-4031`（`Blast_GetNuclAlphaBeta`）、`blast_stat.c:3955-3963` | 長さの補正の alpha と beta：表の行の値。写した場合は alpha = Lambda/H（ungapped）、beta は与えた値が 1/−1 か 2/−3 なら −2、それ以外は 0 |
| `blast_setup.c:74-107`（`Blast_ScoreBlkKbpGappedCalc`） | 有効な context ごとに計算し、最初の失敗で止まる。有効な context が無ければ計算も失敗もしない |
| `blast_aux_priv.cpp:94-102`、`blast_aux.cpp:1013-1024`、`blast_types.hpp:180-183`、`setup_factory.cpp:170-184`、`blast_app_util.hpp:225-227` | context を持たないメッセージは、その batch の全 query に付く。例外の文言は query ごとに `Error: <msg> ` をつないだもので、`BLAST engine error: ` を付けて終了コード 3 |
| `setup_factory.cpp:170-172`、`blast_stat.c:2784-2791` | 例外は、最初のメッセージが error のときだけ投げる。無効な context があると、その warning が先に来るので、NCBI は gapped block の無いまま続ける（S07+ の独立監査で、NCBI が segfault することを確かめた） |
| `blastn_app.cpp:258-261` | outfmt 0 の prolog（program、参照文献、Database の行）は検索の前に出るので、この失敗でも stdout に残る |

LOSAT：`tables.rs` の `NuclValues`（`new`・`check_gaps`・`gapped`）が 3 つの関数の移植、`scoring.rs` の `context_blocks` と `karlin_error`。失敗は、最初の batch の query がすべて有効なら NCBI と同じ文言（繰り返しの回数は最初の batch の query の数）と prolog で終わる。最初の batch に無効な query があれば、NCBI は失敗を報告しない（上の行）ので、LOSAT の文言で明示的に失敗する。

**S07+ の前の差（根本原因 2）：** LOSAT は最大公約数で割らず、表に無い組は 1/−2 の値を、表に無い gap は表の最初の行を黙って使い、拒否もしなかった。

## D. query の context ごとの block と HSP の順序

| NCBI | 振る舞い |
|---|---|
| `blast_stat.c:2762-2783`（`Blast_ScoreBlkKbpUngappedCalc`） | context（query の鎖）ごとに、その配列の組成から ungapped block を計算する。失敗した context は無効 |
| `blast_parameters.c:218-221,342-374` | ungapped の X-drop と gap trigger は context の ungapped block から |
| `blast_setup.c:786-846` | 長さの補正と探索空間は context の gapped block から |
| `blast_parameters.c:92-116,455-463` | gapped の X-drop は、batch の有効な context の gapped Lambda の最小値から |
| `blast_seqalign.cpp:1571-1577,1801-1804` | Seq-align を作るとき、subject ごとの HSP を e-value（同じなら `ScoreCompareHSPs`）で並べ直す |

2 つの鎖の組成は同じ残基の並べ替えなので、block は最後の桁で違いうる。表の gap では gapped block が表の値で全 context 同じなので差は出ないが、表を超える gap（ungapped block を写す）では、同じ得点の HSP の e-value が鎖で違い、HSP の順序が変わる（NCBI では minus の鎖の HSP が先に来る例を確かめた）。

LOSAT：S07+ の前は、エンジンは ACGT だけを仮定した 1 つの ungapped block と 1 つの gapped block をすべての query に使い、報告だけが query の組成の block を使っていた（S07 の §G.2）。S07+ で、context ごとの ungapped と gapped の block（`ContextKarlin`）を、cutoff、X-drop、gap trigger、探索空間、e-value、bit score に使うようにした。無効な context は、NCBI と同じく無限の cutoff（`blast_parameters.c:324-332`）、0 の X-drop と探索空間（`blast_setup.c:778-786`）を持ち、ヒットを持たない。有効な context が無ければ、NCBI と同じく検索しない（`local_blast.cpp:177-180`。以前は、例えば `-reward 4 -penalty -1` の 100 kb の query で 300 秒以上かかった）。gapped の X-drop は、全 context で同じ値になるとき、または最初の batch が全 query を含むときに NCBI と同じ。そうでない場合（表を超える gap で、組成の違う query が最初の batch に収まらない）は明示的に失敗する（`scoring.rs` の `gap_x_dropoffs`）。HSP の e-value の並べ直しは `post_process_hits_and_write` に足した。表の gap では並びは変わらない。

## E. 入力の読み方

| NCBI | 振る舞い |
|---|---|
| `blast_args.cpp:2553-2557`、`objmgr_query_data.cpp:375-380` | 引数の処理の中で（§B）subject を読む。レコードが無い（空、空白だけ、またはディレクトリの）subject は「BLAST engine error: Empty CBlastQueryVector」で終了コード 3。query が空でも、オプションの検査の誤り（`-penalty 0` など）があっても、`-out` を作る前で、この誤りが先に出る |
| `blast_app_util.cpp:856-866` | 次に、query のファイルに空白以外の文字が無ければ「Warning: [blastn] Query is Empty!」で成功する（ディレクトリも空）。位置の無いストリーム（パイプ）は空と見なさない |
| `fasta_reader_utils.cpp:168-225`（`ParseDefLine`） | `>` の後の空白を飛ばし、ID は `' '` 以下のバイトで終わり、title は最初の制御文字（`' '` 未満）で終わる |
| `fasta.cpp:856-935`（`ParseDataLine`） | 行の中の空白は無視、`-` は無視して警告、`;` から行末は注釈、IUPAC 以外は取り除いて「FASTA-Reader: Ignoring invalid residues」。`U` は残基として残り、BLAST+ は `T` として検索し表示する（大文字・小文字、`T` との混在のどれでも。oracle で確認） |
| `showalign.cpp` などの折り返し | title の折り返しはバイト単位で、多バイト文字を分けることがある |

LOSAT は FASTA を `bio` で読む（Web のアダプタの索引の走査も `bio` を再現する）。空白だけのファイルとディレクトリは、NCBI と同じくレコードの無いファイルとして読む（`input.rs` の `is_blank`、`run.rs` の `read_blastn_fasta_bytes`。CLI とアダプタの `register`）。標準入力（`/dev/stdin`）の query と subject は、NCBI と同じく読む（`check_inputs.py` の `audit3.piped_*`）。`run_local` は、レコードの無い subject に NCBI の engine error を、レコードの無い query に「Query is Empty!」を出す。S07+ では、`U` を `T` として読み（`input.rs` の `with_u_as_t`）、`bio` と NCBI で読み方が違う入力を明示的に拒否する：`>` の直後の空白、制御文字（tab など）や非 ASCII のバイトを含む定義行（`check_deflines`。`bio` は ID の後の区切りの文字を捨てるので、ファイルのバイトで調べる）、IUPAC の文字以外の残基（`check_residues`）。拒否の文言には「not supported by LOSAT's BLASTN」を含む。NCBI の `CFastaReader` の移植は計画の後回しの項目にした。

## F. query の batch と表形式

| NCBI | 振る舞い |
|---|---|
| `blast_app_util.hpp:60-76`、`blast_app_util.cpp:65-81` | `CBatchSizeMixer` の batch は 100 残基以上（`k_MinBatchSize`） |
| `tabular.cpp:1264-1284`、`blast_format.cpp:762` | 検索されなかった batch の query には alignment の集合が無いので、outfmt 7 の見出しに「# N hits found」が無い |
| `blast_format.cpp:129-138,336-344` | 報告と batch の大きさは、subject の総文字数を `int` で受け取る |

LOSAT：最初の batch の後の、100 残基未満の無効な query の連なりで後ろに有効な query があるものは、有効な query と同じ batch で検索される（`unsearched_queries`）。outfmt 7 は outfmt 0 と同じ規則で検索されなかった query を決め、決まらない場合は明示的に失敗する。subject の総文字数が 2^31 以上なら明示的に失敗する。

## G. LOSAT が明示的に拒否するもの

| 入力 | 文言の要点 |
|---|---|
| 16 ビットに折り返した reward が 0 以下（`-reward 0` を含む。NCBI は rmblastn の行列の得点か、すべての query が無効） | `a reward of N (NCBI BLAST+'s 16-bit value; …) is not supported by LOSAT's BLASTN`（NCBI の検査の後） |
| 無限の e-value（`+inf`、`1e400` など）、NaN（`+nan`） | `an infinite or NaN e-value is not supported by LOSAT's BLASTN`（NCBI の検査の後）、`a NaN e-value …`（引数） |
| 10 進数でない書き方（NCBI の `0x10` など） | `expected a decimal number (other forms that NCBI BLAST+ reads are not supported by LOSAT's BLASTN)`。NCBI が引数の誤りにする書き方（`inf`、` 1`）は、LOSAT も引数の誤り（`expected a number`）にする |
| 空の定義行、先頭の空白、制御文字、非 ASCII | `has a defline that …; … not supported by LOSAT's BLASTN` |
| 残基の無いレコード（NCBI は「Sequence contains no data」） | `has no residues; … not supported by LOSAT's BLASTN` |
| `bio` が読めない FASTA：最初の定義行の前の文字（空行、`;` の注釈、BOM）、UTF-8 でないバイト（NCBI は読む） | `failed to read … FASTA … (FASTA that bio cannot read (such as text before the first defline or bytes that are not UTF-8), which NCBI BLAST+ may read, is not supported by LOSAT's BLASTN)` |
| 位置の無いストリーム（パイプ）からの空の query（NCBI は空と見なさない） | `an empty query from a stream without a position … is not supported by LOSAT's BLASTN` |
| 環境変数 `BATCH_SIZE`・`CHUNK_SIZE`（NCBI の query の batch を変える。NCBI の検査の後） | `the environment variable … is not supported by LOSAT's BLASTN` |
| 16 ビットの reward − penalty が 3000 を超える得点（LOSAT の Karlin-Altschul の計算は 1000/−2000 まで NCBI と比べた。NCBI の表の失敗より前に調べる） | `… span more than 3000 score units, which is not supported by LOSAT's BLASTN` |
| greedy な伸長（`-task megablast`）の 32767 を超える gap で、表が受け付けるもの（NCBI の 32 ビットの距離が溢れうる。表の検査の後、有効な query があるとき） | `gap costs above 32767 with greedy extension … are not supported by LOSAT's BLASTN` |
| IUPAC の文字以外の残基（`X`、`-`、数字、空白など） | `has 'X' at residue N, which is not an IUPAC nucleotide letter; … not supported by LOSAT's BLASTN` |
| 得点の表の失敗で、最初の batch に無効な query がある（NCBI は報告せず、続けて落ちる） | `… the first query batch … has an invalid query; NCBI BLAST+ does not report the error then, which is not supported by LOSAT's BLASTN` |
| 表を超える gap で、context の gapped X-drop が違い、query が最初の batch に収まらない | `the gapped X-drop of these gap costs depends on … which LOSAT does not reproduce` |
| 最初の batch の後の無効な query で、その batch が決まらないもの（outfmt 0 と 7） | `the outfmt 0 and 7 reports of the invalid query N depend on NCBI BLAST+'s adaptive query batches, which LOSAT does not reproduce` |
| subject の総文字数が 2^31 以上 | `… 32-bit int, which is not supported by LOSAT's BLASTN` |

## H. greedy な伸長のメモリ（S07+ の独立監査）

| NCBI | 振る舞い |
|---|---|
| `blast_gapalign.c:240-251` | affine の greedy な伸長の配列は検索ごとに 1 度 `calloc` で取る。`last_seq2_off` の行は距離 `max_cost` までは全対角線の幅 |
| `greedy_align.c:920-942` | 毎回初期化するのは、先に読まれる最初の `xdrop_offset` の最大得点と `max_penalty` の対角線の範囲だけ |
| `greedy_align.c:1181-1197` | それより先の距離の行は、traceback が無ければ `d - max_penalty - 1` の行を使い回し、あれば memory pool から、その距離で調べる対角線の範囲だけを取る |

LOSAT は、すべての距離に全幅の行（`scaled_max_dist + max_penalty + 2` 行）を取り、呼出しごとに全部を初期化していた。そのため、gap のコストに比例して時間とメモリがかかった（multi の入力で 1000/500 に 4 秒、4000/2000 に 17 秒、10^6 でメモリ不足。NCBI は 0.02 秒）。また、1 億要素を超えると黙って収束しなかったことにしていた。S07+ で、行を「その距離で調べる対角線の範囲」だけにし（traceback が無ければ `max_penalty + 1` 行の輪、あれば距離ごと）、NCBI が初期化する要素だけを初期化するようにした（`greedy.rs` の `AffineRows`）。行は、その距離の対角線の範囲の中でだけ読まれるので、結果は変わらない。32767/32767 で 0.04 秒になり、出力は NCBI と同じ（10^6 では 0.7 秒・213 MB だったが、32767 を超える greedy の gap は、NCBI の 32 ビットの計算が溢れうるので LOSAT が拒否する。§G）。影響を受けない古い greedy の関数（`greedy_align_one_direction_ex` と `affine_greedy_align_one_direction_with_max_dist`。後者は 1 億要素を超えると黙って収束しなかったことにする）は残っているが、検索からは呼ばれない（`extend_gapped_heuristic_with_scratch` は DP の経路だけで呼ばれる）。
