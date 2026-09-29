# E2c（S07+）：BLASTN の得点のオプションと入力の読み方 — NCBI の経路

NCBI のソースは固定 commit `598d8ae6`（`/mnt/c/Users/genom/GitHub/ncbi-blast/c++`）、振る舞いの確認は NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin/blastn`）で行った。LOSAT の移植箇所には、それぞれの直上に NCBI のファイル・行と断片を書いてある。この文書は、どの NCBI の経路が LOSAT のどこに当たるかと、S07+ の前に何が違っていたかをまとめる。

## A. 得点のオプションの値

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_nucl_options.cpp:137-150,198-221`（`SetLookupTableDefaults`・`SetMBLookupTableDefaults`・`SetScoringOptionsDefaults`・`SetMBScoringOptionsDefaults`）、`blast_options.h:67-68,86-95` | task の既定値：blastn は word size 11、reward 2、penalty −3、gap 5/2。megablast は 28、1、−2、0/0 | `coordination.rs` の `task_defaults` |
| `blast_args.cpp:166-170,175-182,262-271,288-300,648-660,674-679` | `-word_size`・`-gapopen`・`-gapextend`・`-reward`・`-penalty` は省略できるキーで、与えたときだけ task の既定値を置き換える。reward は 0 以上、penalty は 0 以下、gap は制約なし | `args.rs` のこれらを `Option` にした。reward は 0 以上、penalty は 0 以下（NCBI の引数の制約）、word size は 4 以上で、100 以下は §B の検査 |
| `ncbiargs.cpp:118-129,383-390`、`ncbistr.cpp:788-791` | 整数の引数（`-reward`・`-penalty`・`-gapopen`・`-gapextend`・`-word_size`・`-max_target_seqs`・`-max_hsps`・`-num_threads`）は、符号付きの 10 進数として読めなければ、`0x`・`0X` で始まるものを 16 進数として読む。`int` の範囲の外は引数の誤り。`0x` だけは 0 だが、制約のある引数（`-reward`・`-penalty`・`-word_size`・`-max_target_seqs`・`-max_hsps`・`-num_threads`）は、制約が値を `NStr::StringToDouble` で読み直し（`blast_input_aux.hpp:110-113`、`ncbiargs.cpp:1231-1243`）、`0x` だけは読めないので引数の誤り（第 5 回の監査） | `value_parsers.rs` の `ncbi_integer`（S07+ の第 4 回の監査。BLASTN の引数だけが使う。ほかの program は S08+） |
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
| `blast_args.cpp:3624-3627,2745-2748,2801-2851`、`blast_args.hpp:1067-1072`、`blast_app_util.hpp:260-263` | `-outfmt` は、引数の組の処理より前に解析する。前後の空白（C の `isspace`）を除き、最初の空白までを符号付きの 10 進数として読む。読めなければ `BLAST query/options error: '<値>' is not a valid output format`（終了コード 1）、0 から 21 の外は `Error: Formatting choice is out of range`（255）。残り（独自の指定）は表形式（6・7 など）でだけ使い、0 では無視する（`0 qaccver` は 0）。ただし、指定が `delim` で始まれば、形式の数より前に、最初の語に `=` が無いと `Delimiter format is invalid. Valid format is delim=<delimiter value>`（終了コード 1。`0 delim` も。`blast_args.cpp:2810-2825`、第 5 回の監査）。`6 delim=,` は区切りの文字の指定、`delim=` は既定（tab） |
| `blast_args.cpp:3631-3633`、`blastn_args.cpp:63-70`、`blast_args.cpp:2553-2557` | 引数の組の処理は順に行う。`CBlastDatabaseArgs`（`-subject`）が最初で、ここで subject のファイルを開いて読む（§E。S07+ の第 3 回の監査） |
| `blast_args.cpp:3456-3481`、`ncbiargs.cpp:717-735,777-781`、`ncbiargs.cpp:615-619` | 次の `CStdCmdLineArgs` が query、`-out` の順にファイルを開く（ファイルは handler が求めたときに開く）。開けなければ `Command line argument error: Argument "<名前>". File is not accessible:  `<値>'`（終了コード 1）。`-` は標準入力・標準出力。`-query` の既定値は `-`（`cmdline_flags.cpp:47`）。`-out` は検査より前に作るので、検査で止まっても空のファイルが残る |
| `blast_args.cpp:375-384,410-426`、`blast_objmgr_tools.cpp:195-198`、`symdust.cpp:170-175` | その後の `CFilteringArgs` が `-dust` を読む。`no`・`yes` 以外は 1 つの空白で区切り（空白が続けば空の欄）、3 つでなければ `Invalid number of arguments to filtering option`、10 進数でなければ `Invalid input for filtering parameters`（どちらも query/options error、終了コード 1）。負の値や 0 も受け付け、`Uint4` にして、DUST の範囲（level 2〜64、window 8〜64、linker 1〜32）の外は既定値（20、64、1）になる |
| `blast_args.cpp:2975-2977` | その後の形式の引数の処理で「Examining 5 or more matches is recommended」を出す |
| `blast_args.cpp:3636-3639`、`blast_app_util.hpp:172-175` | 最後に検査する。失敗は `CInputException` になり、`BLAST query/options error: <msg>` と `Please refer to the BLAST+ user manual.` を出して終了コード 1 |
| `ncbiargs.cpp:464`、`ncbistr.cpp:1312-1318,1328-1339`、`ncbiargs.cpp:4821-4830` | `-evalue`・`-perc_identity` は `NStr::StringToDouble` で読む。最初の文字が数字、小数点、符号でなければ引数の誤り（`inf`、`nan`、` 1`）。それ以外は `strtod` でも読むので、`+inf`・`-inf`・`+nan`・`0x10`（16）・`1e` も受け付ける（2.17.0 で、`+inf` と `+nan` は検索し、`0x10` は `-evalue 16` と同じ出力）。`-perc_identity` は 0 から 100（NaN は範囲の外） |

LOSAT：CLI の `run` は、NCBI と同じ順に、形式の解析（`hsp.rs` の `parse_blastn_output_format`、NCBI の文言と終了コード）、subject のファイルを開いて読む（レコードが無ければ engine error）、query のファイルを開く、`-out` を開く（開けないファイルの文言は BLASTX と共有する `cli.rs` の `inaccessible`。空の名前には値の部分が無い）、`-dust` を読む（`args.rs` の `resolve_dust`、`value_parsers.rs` の `parse_dust_filtering`）、警告と検査（`process_options`、`scoring.rs` の `check_scoring_options`）を行い、その後で query を読む。各ファイルは 1 度だけ開く（名前付きパイプは 1 つの読み手にだけバイトを渡す。`-out` も、第 5 回の監査で 1 度にした）。`-out` の名前（NCBI の `CDirEntry::GetName` と同じく、末尾の区切りを除いた最後の部分。`…/.` は `.`）が 256 バイト以上なら引数の誤り（`value_parsers.rs` の `blastn_output_path`、`ncbifile.cpp:298-312,358-363,465-472`）。`run_local` も、形式の解析、subject のレコードの検査、`-dust`、`process_options` の順。CLI の終了コードは `cli.rs` の `NativeError`（S07+ で `blastx/native.rs` から移した）で NCBI と同じにした。LOSAT の上限（§G。`scoring.rs` の `check_losat_limits` と、環境変数 `BATCH_SIZE`・`CHUNK_SIZE`）は、NCBI が検索を始める所（「Query is Empty!」の成功の後、Karlin-Altschul の表の失敗の前）で調べる（第 4 回の監査：以前は検査の直後にあり、空の query で NCBI が成功する `-reward 0` などを拒否していた）。greedy の gap の上限は表の検査の後（エンジン）で調べる。subject の定義行の拒否（`check_deflines`）と残基の無いレコードの拒否（`input.rs` の `check_records_have_residues`）も、NCBI が何も出さずに読むので、「Query is Empty!」の後に行う。環境変数 `CHUNK_SIZE` は、NCBI と同じく空白だけなら無いものとする（`IsBlank`）。IUPAC 以外の残基の拒否と `bio` が読めない FASTA の拒否は、NCBI が読み込みの時点で警告を出すか、レコードが無いと判断しうるので（`bio` はレコードを区切れない）、subject を読む時点のまま（定義行の UTF-8 でないバイトを含む）。LOSAT の thread の検査（`validate_threads`）は形式の解析の後、入力を開く前。

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
| `blast_app_util.cpp:856-866` | 次に、query のファイルに空白以外の文字が無ければ「Warning: [blastn] Query is Empty!」で成功する（ディレクトリも空）。位置の無いストリーム（パイプ）は空と見なさない。`-query -` と `-subject -` では、subject を読んだ `cin` が終わりに達して位置を返さないので、query も位置の無いストリームになる（2.17.0 は警告なしに空の報告を出す） |
| `fasta_reader_utils.cpp:168-225`（`ParseDefLine`） | `>` の後の空白を飛ばし、ID は `' '` 以下のバイトで終わり、title は最初の制御文字（`' '` 未満）で終わる |
| `fasta.cpp:856-935`（`ParseDataLine`） | 行の中の空白は無視、`-` は無視して警告、`;` から行末は注釈、IUPAC 以外は取り除いて「FASTA-Reader: Ignoring invalid residues」。`U` は残基として残り、BLAST+ は `T` として検索し表示する（大文字・小文字、`T` との混在のどれでも。oracle で確認） |
| `showalign.cpp` などの折り返し | title の折り返しはバイト単位で、多バイト文字を分けることがある |

LOSAT は FASTA を `bio` で読む（Web のアダプタの索引の走査も `bio` を再現する）。空白だけのファイルとディレクトリは、NCBI と同じくレコードの無いファイルとして読む（`input.rs` の `is_blank`、`run.rs` の `read_blastn_fasta_bytes`。CLI とアダプタの `register`）。標準入力（`/dev/stdin` と `-`。`-` は標準入力を複製したファイルで、位置を共有する）と名前付きパイプの query と subject は、NCBI と同じく読む（`check_inputs.py` の `audit3.piped_*`・`audit4.stdin_*`・`audit4.fifo.*`）。`run_local` は、レコードの無い subject に NCBI の engine error を、レコードの無い query に「Query is Empty!」を出す。S07+ では、`U` を `T` として読み（`input.rs` の `with_u_as_t`）、`bio` と NCBI で読み方が違う入力を明示的に拒否する：`>` の直後の空白、制御文字（tab など）や非 ASCII のバイトを含む定義行（`check_deflines`。`bio` は ID の後の区切りの文字を捨てるので、ファイルのバイトで調べる）、IUPAC の文字以外の残基（`check_residues`）。拒否の文言には「not supported by LOSAT's BLASTN」を含む。NCBI の `CFastaReader` の移植は計画の後回しの項目にした。

## F. query の batch と表形式

| NCBI | 振る舞い |
|---|---|
| `blast_app_util.hpp:60-76`、`blast_app_util.cpp:65-81` | `CBatchSizeMixer` の batch は 100 残基以上（`k_MinBatchSize`） |
| `tabular.cpp:1264-1284`、`blast_format.cpp:762` | 検索されなかった batch の query には alignment の集合が無いので、outfmt 7 の見出しに「# N hits found」が無い |
| `blast_format.cpp:129-138,336-344` | 報告と batch の大きさは、subject の総文字数を `int` で受け取る |
| `blast_hits.c:44-68` | 予備の hit list の大きさ `MIN(MAX(2 * size, 10), size + 50)` は `Int4` で、2^30 から 2^31 − 51 では 2 倍が溢れて 10 になり、それより大きいと負になって NCBI は落ちる（第 5 回の監査） |

LOSAT：最初の batch の後の、100 残基未満の無効な query の連なりで後ろに有効な query があるものは、有効な query と同じ batch で検索される（`unsearched_queries`）。outfmt 7 は outfmt 0 と同じ規則で検索されなかった query を決め、決まらない場合は明示的に失敗する。subject の総文字数が 2^31 以上なら明示的に失敗する。予備の hit list の大きさは NCBI と同じく 32 ビットで折り返し（`hsp.rs` の `get_prelim_hitlist_size`）、負になる `-max_target_seqs`（2147483598 以上）は明示的に失敗する（`check_losat_limits`）。

## G. LOSAT が明示的に拒否するもの

| 入力 | 文言の要点 |
|---|---|
| 16 ビットに折り返した reward が 0 以下（`-reward 0` を含む。NCBI は rmblastn の行列の得点か、すべての query が無効） | `a reward of N (NCBI BLAST+'s 16-bit value; …) is not supported by LOSAT's BLASTN`（NCBI の検査の後） |
| 無限の e-value（`+inf`、`1e400` など）、NaN（`+nan`・`-nan`） | `an infinite or NaN e-value is not supported by LOSAT's BLASTN`（「Query is Empty!」の後） |
| `-evalue`・`-perc_identity` の 10 進数でない書き方（NCBI の `0x10`、`0x1p-3`、`1e` など） | `expected a decimal number (other forms, which NCBI BLAST+ may read, are not supported by LOSAT's BLASTN)`（引数）。NCBI が引数の誤りにする書き方（`inf`、` 1`）は、LOSAT も引数の誤り（`expected a number`）にする |
| 出力形式 0・6・7 以外（NCBI が受け付ける 1〜5・8〜21）と、表形式の独自の指定（欄や `delim=`） | `output format N is not supported by LOSAT's BLASTN`、`the custom output format specification … is not supported by LOSAT's BLASTN` |
| `-query -` と `-subject -` の組、またはパイプからの空の query（NCBI は位置の無いストリームを空と見なさない） | 上の「位置の無いストリーム」の文言 |
| 予備の hit list の大きさが負になる `-max_target_seqs`（NCBI は落ちる。「Query is Empty!」の後） | `… whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTN` |
| Rayon の上限（65535）を超える `-num_threads`（NCBI は CPU の数に減らす） | `… exceeds Rayon maximum 65535, which is not supported by LOSAT` |
| システムが起動できない数の `-num_threads`（NCBI は CPU の数に減らす） | `failed to build blastn pool with N threads; a thread count that the system cannot start is not supported by LOSAT` |
| NCBI の blastn の他の task（`blastn-short`・`dc-megablast`・`rmblastn`。`blast_options_handle.cpp:211-222`、`blastn_args.cpp:57-59`）と、LOSAT に無い NCBI の blastn のオプション（2.17.0 の `-help` の 53 個：`-db`、`-strand`、`-ungapped`、`-xdrop_gap`、`-num_alignments`、`-h`、`-version` など。`cli.rs` の `is_unported_blastn_arg`） | `the task … is not supported by LOSAT's BLASTN`、`the NCBI BLAST+ option -… is not supported by LOSAT's BLASTN`（引数。第 6 回の監査。以前は理由の無い引数の誤りだった） |
| 空の定義行、先頭の空白、制御文字、非 ASCII | `has a defline that …; … not supported by LOSAT's BLASTN` |
| 残基の無いレコード（NCBI は「Sequence contains no data」） | `has no residues; … not supported by LOSAT's BLASTN` |
| `bio` が読めない FASTA：最初の定義行の前の文字（空行、`;` の注釈、BOM）、UTF-8 でないバイト（NCBI は読む。subject では、読む時点で、NCBI のオプションの検査より前） | `failed to read … FASTA … (FASTA that bio cannot read (such as text before the first defline or bytes that are not UTF-8), which NCBI BLAST+ may read, is not supported by LOSAT's BLASTN)` |
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

## I. NCBI blastn に無いオプション（S07+ の第 4 回の監査）

`-verbose`・`-limit_lookup`・`-max_db_word_count`・`-min_hit_length` は、LOSAT の BLASTN の CLI にだけあった（`-limit_lookup` と `-max_db_word_count` は magicblast のオプション、`blast_args.cpp:1492-1501`。`-verbose` と `-min_hit_length` は NCBI に無い）。NCBI blastn は「Unknown argument」で終了コード 1 にするが、LOSAT は受け付けて実行していた。`AGENTS.md` の規則 5 により CLI から除いた（`args.rs` の値は NCBI の既定値のまま、`#[arg(skip)]`）。`-limit_lookup` の検査（`check_blastn_lookup_options`）は到達しなくなったので除いた。

NCBI は `-num_threads` が CPU の数を超えると「Number of threads was reduced to N to match the number of available CPUs」、`-subject` があると「'num_threads' is currently ignored when 'subject' is specified.」と警告し、1 スレッドで検索する（`blast_args.cpp:3203-3236`）。LOSAT はどの program も `-subject` でスレッドを使い、警告を出さない（出力は同じ）。この stderr の差は以前からで全 program に及ぶので、計画の未決事項に記録した（S07+ では変えない）。

LOSAT が明示的に拒否し、NCBI の文言を再現していないもの（どれも LOSAT の文言で失敗する）：NCBI が受け付ける出力形式のうち 0・6・7 以外（NCBI 自身の誤りになる 21「FASTA output format is only applicable to magicblast」と、`-out` の無い 13・14「Please provide a file name…」を含む）。全 program に共通の CLI の誤りの経路のうち、出力の書き込みの失敗（`-out /dev/full`。NCBI は outfmt 0 で「BLAST failed to write output」と終了コード 6、表形式は abort。LOSAT の BLASTN は終了コード 1）と、UTF-8 でないファイル名の表示（`inaccessible` は U+FFFD にする。NCBI は元のバイト）と、UTF-8 でない引数の値（`-dust`・`-outfmt` などに。LOSAT は clap の誤り、NCBI は自身の誤り）は、BLASTN 以外の program と合わせて S08+ で扱う（第 5・6 回の監査）。

## J. 1 つの query の two-hit の対角線の表（S07+ の第 6 回の監査、以前からの不具合）

| NCBI | 振る舞い |
|---|---|
| `blast_engine.c:1002-1003`、`blast_extend.c:52-61` | two-hit の対角線の表（`BLAST_DiagTable`）の大きさは、`query->length`（batch の query の塊。1 つの query では 2 つの鎖と区切りの 2L+1）に window を足した長さ以上の 2 の累乗 |
| `na_ungapped.c:631-700`（`s_BlastnDiagTableExtendInitialHit`） | 対角線 `s_off + 表の長さ − q_off` を表の長さで丸める。表が query の塊より window 以上長いので、同じ要素に入る 2 つの対角線は、subject で window 以上離れ、前の記録は期限切れとして扱われる |

LOSAT は、1 つの query の表（`use_array_indexing`。query の塊が 8000 以下のとき。それより長いか query が複数なら NCBI の hash と同じ扱いの hash）の大きさを、一方の鎖の長さ L から求めていた。そのため、plus 鎖と minus 鎖の対角線が、subject の近い位置で同じ要素を使い、一方の鎖の伸長の `last_hit` が他方の鎖の seed を落とした（例：`strand_query` × `strand_subject`、`-task blastn -word_size 5 -evalue 100` で NCBI の 103 行に 102 行。既定のオプションでも、逆位の反復が表の長さだけ離れた query で HSP が 1 つ欠けた）。S07+ で表の大きさを query の塊の長さから求めるようにした（`run.rs`）。query が複数の場合の hash は、期限切れの記録の扱いが NCBI の表とほぼ同じなので変えていない（one-hit の hash はセルを再利用し、初期の hit が重複しうる：`na_ungapped.c:396-449,939-941`。第 7 回の監査の、複数の query の塊が 8000 以下の 6000 件で出力の差は無かった）。`word_size_sweep.py`（8 つの fixture の組 × 2 つの task × word size 4〜28 × e-value 10・1e5 の 256 件）は、変更前に 16 件が NCBI と違い、変更後は 0 件。

## K. 第 7 回の独立監査で見つかった、以前からの差

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_gapalign.c:2679-2681,2881-2882`（`s_ReduceGaps`） | greedy な traceback の後の gap の縮小で、12 以上の一致の前を後ろへ調べる。query の側は整列の始まりで止め、subject の側は整列の始まりより前も読み、subject の先頭の sentinel（どの query の文字とも一致しない）で止まる | subject の側の範囲を確かめず、subject の先頭を越えて panic した（既定の megablast で、subject が query より少ない回数の縦列反復の中で始まるとき）。subject の先頭で止めるようにした（`greedy.rs` の `reduce_gaps`） |
| `blast_engine.c:586-591`、`blast_hits.c:2537-2606` | `-subject_besthit` の `Blast_HSPListSubjectBestHit` は、traceback の後に加えて、予備の段階でも、subject の chunk ごとに合わせた HSP の一覧（全 query の HSP、context ごとの長さ）に使う | 予備の段階で使っていなかった（traceback の前に消える HSP が残り、最終の HSP が増えた）。`subject_best_hit.rs` を HSP の欄を与える形（`subject_best_hit_by`）にし、`run.rs` の `prelim_subject_best_hit` で chunk を合わせるたびに使う |
| `align_format_util.cpp:986-990` | 99.9 を超え 99999 以下の bit score は `%3.0ld` の `(long)bit`。(99.9, 100) は ` 99`（幅 3） | 幅を付けていなかった（`outfmt6.rs` の `format_bitscore_ncbi`・`write_bitscore_ncbi`。全 program が共有する。既定の得点では出ない） |

`slice_sweep.py`（リポジトリのゲノムから固定の seed で切り出した query 2〜20 kb と subject 5〜30 kb の組を、与えたオプションで NCBI と比べる）は、`-task blastn -subject_besthit` の 120 組で変更前に 4 組が違い、変更後は 0 組。既定、`-task blastn`、`-subject_besthit`、`-word_size 8`、`-reward 3 -penalty -4 -gapopen 10 -gapextend 3` でも 0 組。

**S07++ に移したもの（TD-14）：** NCBI は複数の query を batch に分けて検索し、近似の ungapped 伸長は query の塊（batch）の端で止まる（`na_ungapped.c:164,286,317`）。最初の batch は約 5000 残基で、その後の大きさは前の batch の成功した ungapped 伸長の数で決まる（`blastn_app.cpp:261-300`、`blast_app_util.cpp:65-81`、`local_blast.cpp:301-307`）。LOSAT は全 query を 1 つの塊で検索するので、batch の境の query の端の seed で、既定のオプションでも HSP が増えうる（監査の例：5339 bp と 10 bp の 2 つの query、`-task blastn` で NCBI の 62 行に 63 行。NCBI に `BATCH_SIZE=100000000` を与えると LOSAT と同じ）。これは以前からの差で、S07+ の範囲（得点のオプションと入力の読み方）の外のエンジンの構造の変更なので、S07++ で NCBI の batch を移植する（S12 の前）。それまで、最初の batch を超える複数の query の結果は、この点で NCBI と違いうる。

## L. 第 8 回の独立監査で見つかった、以前からの差

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_setup_cxx.cpp:842-847,945-950`、`blast_objmgr_tools.cpp:424-478,515-520`、`random_gen.hpp:224-241` | 予備の段階（lookup の走査、ungapped 伸長、greedy・DP の得点だけの伸長）が読む圧縮した subject（ncbi2na）では、IUPAC の曖昧な文字を、subject の長さを種にした `CRandom` から引いた、両立する塩基にする（`N` と gap は 4 つから、2〜3 塩基の記号はその中から）。traceback は blastna で、曖昧な文字をそのまま使う | 固定の対応（blastna の下位 2 ビット：Y→C、K→T、R→A、N→G…）で詰めていた。そのため、曖昧な文字の近くの seed と予備の得点が NCBI と違い、既定のオプションでも 1 つの query で HSP の有無が違った。TBLASTN が既に移植していた `CRandom`（`NcbiRandom`）を `core/blast_encoding.rs` に移して共有し（`resolve_ncbi4na_to_ncbi2na`、`encode_subject_ncbi2na_packed`）、BLASTN の subject の圧縮に使う。TBLASTN の出力は変わらない |
| `lookup_util.c:190-203` | lookup table の大きさの見積もり（`EstimateNumTableEntries`）は、区間ごとに `right - left`（`right` は最後の文字を含む） | 区間の長さ（`end - start`、`end` は含まない）を足し、1 つ多かった。境の長さ（例：megablast の 4250 残基の query）で、NCBI と違う幅の表を選んでいた（出力の差は監査の 1380 件で見つかっていない）。`end - start - 1` にした |
| ABI v1（plan TD-1） | S07+ の前の v1 は、thread の数を最初に確かめ、空の query に空の報告を返した | 空の query でも定義行と LOSAT の上限を調べていた。空の query では、NCBI と同じ誤り（空の subject、オプションの検査）だけを残し、空の報告を返す（`run_web_pair`） |

`slice_sweep.py` に `ambiguity` の pool（EDL933 の曖昧な文字の周りの窓を subject に、その中の曖昧な文字を両立する塩基にした配列を query に）を足した。各 100 組で、変更前は既定で 1 組、`-task blastn` で 4 組、`-task blastn -word_size 7` で 5 組、`-subject_besthit` で 1 組が NCBI と違い、変更後はどれも 0 組。以前の `slice_sweep.py` の 15 のゲノムには曖昧な文字が無く、5 組は同じ配列だった（区別できる 12 にした）。

## M. 第 9 回の独立監査で見つかった、以前からの差

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `fasta.cpp:350-362,1094-1140` | `>?` で始まる行は定義行ではなく、配列の gap の行（`>?100` は 100 文字の gap、`>?unk100` は長さの分からない gap、`>?`・`>?abc` は警告して 1 文字）。`>?_` は接頭辞を除いた定義行 | `bio` は `?100` などの名前の新しいレコードにし、HSP の subject と座標が NCBI と違った。`check_deflines` が `?` で始まる定義行を拒否する（TD-12 と同じく、NCBI の読み込み器を移植せずに明示的に拒否する） |
| `fasta.cpp:670-686,1624-1643`、`blast_fasta_input.cpp:329-330` | `>` の後の文字列が 20 バイトより長く、最後の 20 バイトがすべて A/C/G/T（大文字・小文字）なら、そのレコードを読む時点で `FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?` を出す（CRLF の CR は除き、行末の空白は残す。2.17.0 で確認） | 出していなかった（stderr だけが違った）。`input.rs` の `write_title_warnings` が、subject を読む時点（CLI の `run`、`run_local` の始め）と query を読む時点（`search` の「Query is Empty!」の後）に出す。`bio` は行末の空白を捨てるので、行末に空白があり、それを除くと警告の条件を満たす定義行は拒否する |
| `blast_filter.c:1081-1115`、`blast_setup.c:633-653` | 最初の文字から最後の文字まで覆われた context は、lookup の区間 `(end + 1, end)` を持ち、見積もりは 1 つ減り、`end` は `max_off` に入る。区間は ungapped block が context を無効にする前に作るので、無効な context も数える | 覆われた context を数えていなかった。`compute_lookup_query_stats` を NCBI に合わせた（出力の差は監査で見つかっていない） |

ABI v1 は `run_local` で CLI と同じエンジンを使うので、エンジンを NCBI に合わせる修正（例：§L の曖昧な文字）は v1 の検索結果にも及ぶ。計画の TD-1 に、凍結するのは v1 の引数・形式・誤りの扱いであることを書き足した。

## N. 第 10 回の独立監査で見つかった、以前からの差

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `fasta.cpp:376,710-758,1003-1012` | 配列の行の端で除くのは ASCII の空白だけで、残る UTF-8 のバイト（U+00A0 など）は不正な残基（警告して除く。レコードの最初の配列の行で文字が少なければ誤り） | `bio` は行末の Unicode の空白（U+0085、U+00A0、U+2000〜200A、U+3000 など）も除き、黙って受け付けていた。`input.rs` の `check_sequence_lines` が、配列の行の非 ASCII のバイトを、読む時点で拒否する（CLI、`register`、ABI v1 の fail-fast） |
| `showalign.cpp:2273`、`showdefline.cpp:498`、`create_defline.cpp:219-312,3431-3446,3952-3960,4050-4095` | outfmt 0 の subject の見出しと説明の一覧の文字列は、`-parse_deflines` の無い subject では定義行全体（ID を含む）を title とした `CDeflineGenerator::GenerateDefline`：末尾の `.,;~ ` を除き、見出しでは `TPA:`・`MAG ` などの接頭辞を除き（説明は `fLeavePrefixSuffix` で残す）、HTML の文字参照を戻し、末尾の `,;~ ` を除き、`x_CleanAndCompress`（空白の連なり、` ,`、`,,`、`( `、` )` などの整理）。`sseqid`（outfmt 6/7）と query の title は定義行のまま | 定義行をそのまま出していた。`report/defline.rs` の `ncbi_nucleotide_title` を、BLASTN の見出し（`write_blastn_pairwise_report`）と説明の一覧（`write_blastn_description_table`）に使う。HTML の文字参照になりうるもの（`NStr::HtmlDecode` が変える `&名前;`・`&#数;`・`&#x16進;` を含む、それより広い形。NCBI の表に無い名前や、NCBI の trim で消える `;` も含む）を持つ subject の定義行は、outfmt 0 を求めたとき明示的に拒否する。定義行が `,;~` と空白だけで、最後の `, `（または `; `）の後が空白と同じ記号だけのとき、NCBI の `x_CleanAndCompress` は文字列の終わりを越えて読み（`size_t left` が折り返す、create_defline.cpp:221,278-281,288-291,301）、2.17.0 は outfmt 0 で落ちる（第 11 回の監査）。LOSAT は NCBI の `left` の減り方をそのまま数え、折り返す定義行を見分ける（`title_sweep.py`：先頭が空白でない長さ 1〜5 の 1023 通りで、NCBI が落ちる 66 通りをすべて拒否し、残りの 957 通りは NCBI とバイト一致）。この定義行も、outfmt 0 を求めたとき「which LOSAT does not reproduce」で拒否する。BLASTP・TBLASTN・BLASTX の outfmt 0 は変えていない（S08+。蛋白の subject は `x_AdjustProteinTitleSuffix` などの規則も持つ） |

第 9 回の `check_deflines` の、行末の空白の拒否の文言を直した（行末が制御文字でも NCBI は警告するので、「空白が無いときだけ警告する」は誤りだった）。

引数の誤り（NCBI は USAGE と誤りを出して終了コード 1、LOSAT は clap の誤りで終了コード 2。例：`-word_size 3`、`-evalue abc`、未知のオプション）は、どちらも引数の誤りで、`check_inputs.py` の `arg-error` の区分にしている。全 program の CLI に共通なので、S08+ で扱いを決める（承認済みの例外にするか、NCBI の文言と終了コードにするか）。
