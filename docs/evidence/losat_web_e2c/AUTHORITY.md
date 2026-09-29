# E2c（S07+）：BLASTN の得点のオプションと入力の読み方 — NCBI の経路

NCBI のソースは固定 commit `598d8ae6`（`/mnt/c/Users/genom/GitHub/ncbi-blast/c++`）、振る舞いの確認は NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin/blastn`）で行った。LOSAT の移植箇所には、それぞれの直上に NCBI のファイル・行と断片を書いてある。この文書は、どの NCBI の経路が LOSAT のどこに当たるかと、S07+ の前に何が違っていたかをまとめる。

## A. 得点のオプションの値

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_nucl_options.cpp:137-150,198-221`（`SetLookupTableDefaults`・`SetMBLookupTableDefaults`・`SetScoringOptionsDefaults`・`SetMBScoringOptionsDefaults`）、`blast_options.h:67-68,86-95` | task の既定値：blastn は word size 11、reward 2、penalty −3、gap 5/2。megablast は 28、1、−2、0/0 | `coordination.rs` の `task_defaults` |
| `blast_args.cpp:166-170,175-182,262-271,288-300,648-660,674-679` | `-word_size`・`-gapopen`・`-gapextend`・`-reward`・`-penalty` は省略できるキーで、与えたときだけ task の既定値を置き換える。reward は 0 以上、penalty は 0 以下、gap は制約なし | `args.rs` のこれらを `Option` にした。reward 0 は LOSAT が拒否する（`positive_i32`。NCBI の rmblastn の 0/0 と、すべての query が無効になる場合を実装しない） |

**S07+ の前の差（根本原因 1）：** LOSAT は clap の既定値を megablast の値にし、task が blastn のときに「既定値と同じ値」を blastn の既定値に置き換えていた。そのため、`-task blastn -reward 1` のように megablast の既定値と同じ値を明示すると、その値が無視された（S06 の sweep の reward 1 の 16 件の DIFF の原因）。

## B. オプションの検査（検索の前）

| NCBI | 振る舞い |
|---|---|
| `blast_options.c:881-906`（`BlastScoringOptionsValidate`） | penalty が 0 以上なら「BLASTN penalty must be negative」、gap open が正で gap extend が 0 なら「BLASTN gap extension penalty cannot be 0」 |
| `blast_options.c:1699-1711`（`s_BlastExtensionScoringOptionsValidate`） | gap が 0/0 で貪欲な伸長でない（blastn の task）なら「Greedy extension must be used if gap existence and extension options are zero」 |
| `blast_options.c:1759-1776`（`BLAST_ValidateOptions`） | 上の 2 つはこの順に調べる |
| `blast_args.cpp:3609-3612`、`blast_app_util.hpp:172-175` | 検査の失敗は `CInputException` になり、`BLAST query/options error: <msg>` と `Please refer to the BLAST+ user manual.` を出して終了コード 1 |
| `blast_args.cpp:2800-2803,2975-2977` | 出力形式の解析と「Examining 5 or more matches is recommended」の警告は、検査より前（引数の処理の時点） |

LOSAT：`scoring.rs` の `check_scoring_options`。`run.rs` の `process_options` が、形式の解析、警告、検査をこの順に、入力を読む前に行う。CLI の終了コードは、BLASTX と共有する `cli.rs` の `NativeError`（S07+ で `blastx/native.rs` から移した）で NCBI と同じにした。

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
| `blastn_app.cpp:258-261` | outfmt 0 の prolog（program、参照文献、Database の行）は検索の前に出るので、この失敗でも stdout に残る |

LOSAT：`tables.rs` の `NuclValues`（`new`・`check_gaps`・`gapped`）が 3 つの関数の移植、`scoring.rs` の `context_blocks` と `karlin_error`。失敗は、最初の batch に有効な query があれば NCBI と同じ文言（繰り返しの回数は最初の batch の query の数）と prolog で終わる。最初の batch がすべて無効で後に有効な query がある場合は、2 番目の batch の大きさが NCBI の batch の計算（§F）に依存するので、LOSAT の文言で明示的に失敗する。

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

LOSAT：S07+ の前は、エンジンは ACGT だけを仮定した 1 つの ungapped block と 1 つの gapped block をすべての query に使い、報告だけが query の組成の block を使っていた（S07 の §G.2）。S07+ で、context ごとの ungapped と gapped の block（`ContextKarlin`）を、cutoff、X-drop、gap trigger、探索空間、e-value、bit score に使うようにした。invalid な query は NCBI が検索しないので、その結果は報告の前に捨てる。gapped の X-drop は、全 context で同じ値になるとき、または最初の batch が全 query を含むときに NCBI と同じ。そうでない場合（表を超える gap で、組成の違う query が最初の batch に収まらない）は明示的に失敗する（`scoring.rs` の `gap_x_dropoffs`）。HSP の e-value の並べ直しは `post_process_hits_and_write` に足した。表の gap では並びは変わらない。

## E. 入力の読み方

| NCBI | 振る舞い |
|---|---|
| `blastn_app.cpp:192-212`、`blast_app_util.cpp:856-873` | subject を読んでから、query のファイルに空白以外の文字が無ければ「Warning: [blastn] Query is Empty!」で成功する |
| `fasta_reader_utils.cpp:168-225`（`ParseDefLine`） | `>` の後の空白を飛ばし、ID は `' '` 以下のバイトで終わり、title は最初の制御文字（`' '` 未満）で終わる |
| `fasta.cpp:856-935`（`ParseDataLine`） | 行の中の空白は無視、`-` は無視して警告、`;` から行末は注釈、IUPAC 以外は取り除いて「FASTA-Reader: Ignoring invalid residues」。`U` は残基として残り、BLAST+ は `T` として検索し表示する（大文字・小文字、`T` との混在のどれでも。oracle で確認） |
| `showalign.cpp` などの折り返し | title の折り返しはバイト単位で、多バイト文字を分けることがある |

LOSAT は FASTA を `bio` で読む（Web のアダプタの索引の走査も `bio` を再現する）。S07+ では、`U` を `T` として読み（`input.rs` の `with_u_as_t`）、`bio` と NCBI で読み方が違う入力を明示的に拒否する：`>` の直後の空白、制御文字（tab など）や非 ASCII のバイトを含む定義行（`check_deflines`。`bio` は ID の後の区切りの文字を捨てるので、ファイルのバイトで調べる）、IUPAC の文字以外の残基（`check_residues`）。拒否の文言には「not supported by LOSAT's BLASTN」を含む。NCBI の `CFastaReader` の移植は計画の後回しの項目にした。

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
| `-reward 0` | clap の値の検査（1 以上） |
| 定義行の先頭の空白、制御文字、非 ASCII | `has a defline that …; … not supported by LOSAT's BLASTN` |
| IUPAC の文字以外の残基（`X`、`-`、数字、空白など） | `has 'X' at residue N, which is not an IUPAC nucleotide letter; … not supported by LOSAT's BLASTN` |
| 最初の batch がすべて無効で、後の batch の得点の表の失敗 | `the error of these scoring options depends on NCBI BLAST+'s adaptive query batches, which LOSAT does not reproduce` |
| 表を超える gap で、context の gapped X-drop が違い、query が最初の batch に収まらない | `the gapped X-drop of these gap costs depends on … which LOSAT does not reproduce` |
| 最初の batch の後の無効な query で、その batch が決まらないもの（outfmt 0 と 7） | `the outfmt 0 and 7 reports of the invalid query N depend on NCBI BLAST+'s adaptive query batches, which LOSAT does not reproduce` |
| subject の総文字数が 2^31 以上 | `… 32-bit int, which is not supported by LOSAT's BLASTN` |
