# E2h（SF）：FASTA の読み方 — NCBI の経路

NCBI のソースは固定 commit `598d8ae6`（引用の file:line は `c++/` からの相対。作業では `$NCBI_SRC/c++` を読んだ）、振る舞いの確認は NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin`）の実行で行った。経路の関数の棚卸しは [`INVENTORY.tsv`](INVENTORY.tsv)（5 つの範囲、222 行。範囲ごとの指示・記録・oracle の観察は [`inventory/`](inventory/)）で、表の最後の列の `LR-12` は範囲と行の番号を指す。oracle の実行の名前（`n_blastn_two_mixed_o6`、`net2`、`e01` など）は棚卸しの `evidence` 列にあり、入力と出力は `$BUILD_ROOT/sf-e2h/inventory/scratch_<範囲>/` にある。「LOSAT 今」は変更前の commit `f3048ffde`（`status_before`）。この記録を書くために新しく NCBI を実行してはいない。棚卸しが未決のまま残した事実（空の `DATA_LOADERS=`）はソースで決めた（§G、§L）。

範囲の記号：LR（行の読み込みと入力の流れ）、RD（レコード：`ReadOneSeq` から `AssembleSeq` まで）、BI（blastinput の層：旗、data loader、batch、中身の無いレコード）、RP（ID と題の報告）、AD（アダプタと ABI）。program の記号：N（BLASTN）、X（TBLASTX）、T（TBLASTN）、P（BLASTP）。

## 要約：4 つの program は 1 つの経路を通り、違うのは分子の旗だけ

- `CBlastFastaInputSource::x_InitInputReader`（`blast_fasta_input.cpp:316-368`）が、役割（query・subject）と program から `fAssumeNuc` か `fAssumeProt` を決める。ほかの旗（`fNoParseID|fDLOptional`、`fNoSplit`、`fHyphensIgnoreAndWarn`、`fDisableNoResidues`、`fQuickIDCheck`）と無視する問題 3 つ、局所 ID の生成器は全員同じ（§A）。`-task`・`-num_threads`・`BLASTINPUT_GEN_DELTA_SEQ` は読み方を変えない（LR-25、LR-27、RD-42、BI-10、BI-45）。
- **subject** は引数の処理の中で（`-query` を開く前、`Query is Empty!` の前、prolog の前）全部読む。読み込みの警告は読んだ時点で stderr に出て、誤りは `BLAST query error: …`（終了コード 1）で止まる。**query** は batch ごとに、prolog と前の batch の報告の後で読む。後の batch の誤りは前の batch の報告が出た後に出る（§B）。
- 読み込みの警告は `LOG_POST` の素の行で、`Warning: [prog]` も ID も付かない。読み込みの誤りは例外で、`BLAST query error: ` が付く。中身の無いレコード・空の batch・全部が空の batch は読み込み器の外（`blast_setup_cxx.cpp`、`objmgr_query_data.cpp`）で警告か終了コード 3 になる（§E、§H）。
- 行の終わりは最初の行で決まり、途中で切り替わる。LF の中の単独の CR は行を切り、その行の終わりを**失う**（次の行に繋がる）。まれに最後の行の尾が消える（pushback の eofbit）。LOSAT は BLASTX の `FastaStream` の移植と同じ算術でこれを再現する（§C）。
- レコードは、題の byte（制御文字 0x20 未満で切れるが、それ以外の byte はそのまま）、局所 ID（`Query_N`・`Subject_N`、レコードごとに 1 つ進む）、残基（無効な byte は取り除いて警告）、小文字の mask の区間を持つ。`>?` の行は gap（N か X の連なり）、`;` 以降と `!`・`#`・`;` で始まる行は注釈、最初の定義行の前の注釈・空行は読み飛ばす（§D）。
- 報告は題の byte をそのまま使う：tabular の ID は題の最初の `' '` までの語、空なら `Query_N`・`Subject_N`、`# Query:` と `Query=` は生の byte、説明の一覧と見出しは `CDeflineGenerator` で作り直す（HTML の復号、非 UTF-8 の変換）（§F）。
- data loader が使える設定（既定）では、最初の行（または Seq-id の行の次の行）が Seq-id として解析できると NCBI は GenBank か BLAST DB から配列を取り寄せる。LOSAT はこれを**明示的に拒否**する。どの行が Seq-id として試されるかは `CSeq_id` の解析で決める（§G）。
- LOSAT が残す明示的な拒否は §J、NCBI の不具合に見える挙動は §K。

## A. 旗：program × 役割 × 分子

### A1. 分子の旗と読み手

| program | 役割 | 旗の元 | 分子 | 局所 ID | batch の大きさ | 読む時点 | 行 |
|---|---|---|---|---|---|---|---|
| N | query | `QueryIsProtein()` が偽 | `fAssumeNuc` | `Query_` | 混合（目標 `SplitQuery_GetChunkSize - 1000`、最初は 5 x 1000 nt）。`BATCH_SIZE` が優先 | batch ごと、prolog の後 | BI-9、BI-32、BI-33 |
| N | subject | `CBlastDatabaseArgs::IsProtein()` が偽（`Blast_SubjectIsNucleotide`、`blast_args.cpp:2467-2470`） | `fAssumeNuc` | `Subject_` | なし（全部） | 引数の処理の中 | BI-2、BI-3 |
| X | query・subject | 両方核酸 | `fAssumeNuc` | `Query_`・`Subject_` | 10002 | query は batch、subject は引数 | BI-9、BI-36 |
| T | query | `QueryIsProtein()` が真 | `fAssumeProt` | `Query_` | 20000 | batch（`-in_pssm` なしのときだけ読む） | BI-37 |
| T | subject | 核酸 | `fAssumeNuc` | `Subject_` | なし | 引数 | BI-9 |
| P | query・subject | 両方蛋白 | `fAssumeProt` | `Query_`・`Subject_` | 10000 | query は batch、subject は引数 | BI-38 |

`BATCH_SIZE`（環境変数）は 4 つとも守る（整数でないと `Error: NCBI C++ Exception:` と build の path、終了コード 255。LOSAT は TD-15 で拒否済み）。blastn は batch の大きさをあとから `SetBatchSize` し、tblastx・tblastn・blastp は `CBlastInput` の生成時に与えるので、不正な `BATCH_SIZE` の誤りは blastn では prolog の後、ほかの 3 つでは前に出る（BI-1、BI-32、BI-35）。batch の残基数は**レコード全体**の長さの和（残基の無いレコードは 0）（BI-29、E2d §B）。

### A2. 全員に共通の旗

| 旗・設定 | NCBI | 効果 | 行 |
|---|---|---|---|
| `fNoParseID`、`fDLOptional` | `blast_fasta_input.cpp:320-322`。`-parse_deflines` のときは `fParseRawID` だけ | 定義行から ID を作らない（題だけ）。最初の定義行が無くても最初のレコードを読む（`ParseDefLine(">")`）。`-parse_deflines` は LOSAT が拒否 | BI-9、RD-5、RD-8、RD-48 |
| `fAssumeNuc` / `fAssumeProt` | `:328-330` | 取る残基の表（§D3）、題の警告の種類、gap の文字（N か X）。`AssignMolType` は旗のとおりに分子を決める（`seqlen_thresh2guess` が `UINT_MAX` なので推測しない）。蛋白の配列を BLASTN に渡すと、ほとんどが取り除かれて長い警告になる（誤りにならない） | RD-43、RD-44 |
| `fNoSplit` | `:331-334`。`BLASTINPUT_GEN_DELTA_SEQ` が空か未設定のとき | 出力は変わらない（`>?` の gap の行でも。BI-10 の 5 組が byte 一致）。gap の行を移植した後に確かめ直す | RD-42、BI-10 |
| `fHyphensIgnoreAndWarn` | `:338` | `-` を取り除き、行ごとに 1 回警告（§E1） | RD-27 |
| `fDisableNoResidues` | `:340` | 残基の無いレコードも誤りにせず、そのまま残す | RD-40 |
| `fQuickIDCheck` | `:343` | ID を作らないので効果なし | RD-8 |
| `fSkipCheck` | 付けない（`GetSkipSeqCheck()` は常に偽） | `CheckDataLine` が毎レコード走る | RD-18 |
| 無視する問題 | `:358-360`：`ModifierFoundButNoneExpected`、`TooLong`、`TooManyAmbiguousResidues` | 題の `[key=value]`・1000 文字を超える題・最初の行の曖昧な文字が 40% を超える警告は出ない。`InvalidResidue` と `IgnoredResidue` は無視しない | RD-15、RD-14、RD-35、RD-34 |
| 局所 ID の生成器 | `CSeqIdGenerator(1, "Query_" または "Subject_")`（`:364-367`、`fasta_reader_utils.cpp:439-486`） | ファイルごとに 1 つ、すべての batch で共有。レコードごとに 1 つ進む（題が空でも、残基が無くても、`>?_x` の疑似定義行でも、後で `-query_loc` が飛ばしても） | RD-9、BI-11、RP-1 |
| 読み手の class | `UseDataLoaders()` ? `CBlastInputReader`（Seq-id を試す）: `CCustomizedFastaReader`（`:345-356`） | 後者は `x_CloseGap` を空にし、`AssignMolType` を旗に従わせる。前者は §G | BI-11、RD-43 |
| mask | `-lcase_masking` のとき `SaveMask()` を `ReadOneSeq()` の前に呼ぶ（`:383-386`）。query と subject の両方 | 小文字の区間を mask にする（§D5）。4 つの program とも NCBI は受け付ける。LOSAT の BLASTP・TBLASTX は `cli.rs:375,637` で拒否のまま（読み込み器は mask を返し、使うかは program の判断） | BI-24、RD-1、RD-31、RD-46 |

## B. 時点と順序

bl2seq（`-query Q -subject S`）の 4 つの program で同じ。blastn は `blastn_app.cpp:190-327`、tblastx は `tblastx_app.cpp:100-217`、tblastn は `tblastn_app.cpp:186-351`、blastp は `blastp_app.cpp:194-304`（BI-1、BI-35〜38）。

| 順 | NCBI | 出力・終了コード | 行 |
|---|---|---|---|
| 1 | `-subject` の file を開く（`CArg_InputFile::x_Open`、`ncbiargs.cpp:695-737`）。`-` は `std::cin`。失敗は `CArgException` | `Command line argument error: Argument "subject". File is not accessible:  \`NAME'`（空白 2 つ、backquote、quote）、終了コード 1。directory は開けて 0 byte に読める | LR-1、LR-2 |
| 2 | `SetOptions`（`blast_args.cpp:3629-3641`）が `m_Args` の順に handler を走らせる。`CBlastDatabaseArgs`（3 番目）が **subject を全部読む**（`blast_args.cpp:2525-2557` → `ReadSequencesToBlast`、`blast_input_aux.cpp:221-247` → `GetAllSeqs`、`blast_input.cpp:198-220`）。`-subject_loc` はここで文法を調べ、範囲の外は E2d の誤り | 読み込みの警告は読んだ時点で stderr。誤りは `BLAST query error: <msg>`、終了コード 1。subject の誤りは、存在しない `-query`・`-evalue -5` の誤りより先。subject が 0 個なら `BLAST engine error: Empty CBlastQueryVector`、終了コード 3（`objmgr_query_data.cpp:376-380`、prolog の前なので outfmt 0 でも stdout は空）。`GetAllSeqs` は `catch (const exception&) { continue; }` を持たない：Seq-id の失敗も致命的 | LR-10、LR-26、BI-1、BI-2、BI-30 |
| 3 | `-query` を開く（`blast_args.cpp:3454-3470`。省略すると既定 `-` で標準入力）。gzip は magicblast だけ | 開けないと `Argument "query"` の誤り、終了コード 1 | LR-3、LR-4 |
| 4 | ほかの option の handler（`-evalue` など。`-query_loc` は query options の handler） | option の誤り（`BLAST query/options error: …` + `Please refer to the BLAST+ user manual.`）、終了コード 1 | BI-1 |
| 5 | `InitializeSubject`、`InitializeQueryDataLoaderConfiguration`（`blast_app_util.cpp:127-215`）、`IsIStreamEmpty(query)`（`:845-875`） | 空なら `Warning: [prog] Query is Empty!`（stderr）、終了コード 0、stdout は空（outfmt 0・7 でも） | LR-9、BI-7、BI-34 |
| 6 | `CBlastFastaInputSource`、`CBlastInput`（tblastx・tblastn・blastp はここで `BATCH_SIZE` を解く） | 不正な `BATCH_SIZE` は 255 | BI-1、BI-32 |
| 7 | `CBlastFormat` の生成：subject の set-up（`SetupSubjects_OMF`）。残基の無い subject の警告は**ここ**（`Query is Empty!` の後、prolog の前） | `Warning: [prog] Subject_N <title>: Subject sequence contains no data`（§E3） | BI-43、BI-49 |
| 8 | `PrintProlog`（outfmt 0 の header と database の統計） | stdout | BI-35 |
| 9 | batch ごと：`GetNextSeqBatch`（読み込みの警告・誤り・Seq-id の取り寄せはここ）→ `CObjMgr_QueryFactory`（レコードの無い batch は `Empty CBlastQueryVector`）→ `CLocalBlast::Run`（query の set-up の警告・誤り）→ `PrintOneResultSet`（query ごとの警告 → 報告）。`for (; !input.End(); …)`（`blastn_app.cpp:277-282` ほか） | 読み込みの誤りはその時点で終了コード 1：前の batch の報告は書き終わって flush 済み、outfmt 0 の prolog は出ており、epilog は出ない | BI-29、BI-47、BI-51、LR-24 |
| 10 | `PrintEpilog` | 正常終了 | BI-35 |

- `CATCH_ALL`（`blast_app_util.hpp:149-262`）：`CObjReaderParseException` は `BLAST query error: <msg>` 1。`CInputException` は `BLAST query/options error: <msg>` + `Please refer to the BLAST+ user manual.` 1。`CArgException` は `Command line argument error: <msg>` 1。`CBlastException` は eInvalidOptions が `BLAST options error: <msg>` 1、メモリ不足が 4、それ以外が `BLAST engine error: <msg>` 3。`CException`・`std::exception` は `Error: <what>` 255（BI-39）。
- 診断の順（stdout と stderr を混ぜて観察）：`Query is Empty!` → 空の subject の警告 → prolog → [batch ごと：読み込みの警告・誤り → query ごとの警告 → 報告] → epilog。全部の query が空の誤りは prolog の後（BI-35、BI-41）。
- 標準入力が `-query -` と `-subject -` の両方なら 1 つの `cin`：subject が全部読み、query は EOF の stream を見る。`IsIStreamEmpty` は `tellg() < 0` で「空でない」と答え、batch は 0 個で、`Query is Empty!` は出ず、終了コード 0（outfmt 0 は header と epilog、7 は `# BLAST processed 0 queries`、6 は何も出さない）。regular file を `<` で渡しても同じ（LR-1、BI-46）。
- LOSAT 今：subject を先に読む順は同じ。query は file 全体を読んで検査してから検索するので、後の batch の誤りが前の batch の出力の後に出ない。移植後は query を batch ごとに読む（BLASTX の `emit_fasta_batch` と同じ形）（LR-24、RD-49、BI-29）。

## C. 入力の stream と行の読み込み（範囲 LR）

### C1. 空の入力と標準入力

`IsIStreamEmpty`（`blast_app_util.cpp:845-875`、非 Windows の枝 855-874）：`tellg() < 0`（pipe・FIFO・tty、または fail/eof が立った stream）なら「空でない」。そうでなければ skipws の `in >> c`（C locale の空白：space、`\t`、`\n`、`\v`、`\f`、`\r`）が読めなければ空で、位置と状態を戻す。NUL・0x1A・0xA0 は空白ではない。**subject には使わない**。

| `-query` の入力 | NCBI | 行 |
|---|---|---|
| 0 byte・`/dev/null`・空白だけ・空行だけの file・directory・`< empty`・`<&-` | `Warning: [prog] Query is Empty!`、終了コード 0、stdout は空 | LR-6、LR-9 |
| 空の pipe・空の FIFO | 何も出ず終了コード 0。outfmt 6・10 は無出力、7 は `# BLAST processed 0 queries`、0 は header と epilog（`Database: …`、`Matrix: …`、`Gap Penalties: …`）。batch が 0 個 | LR-7、BI-34 |
| 空白・空行だけの pipe、注釈だけの file | `BLAST engine error: Empty CBlastQueryVector`、終了コード 3（outfmt 0 は header が出た後）。同じ byte でも seekable な file なら `Query is Empty!` の 0 | LR-8、BI-50 |
| `-query - -subject -` | 上のとおり、batch が 0 個の成功（subject が `cin` を使い切る） | LR-1、BI-46 |
| subject が空・空白だけ・directory・空の pipe | `BLAST engine error: Empty CBlastQueryVector`、終了コード 3、stdout は空 | LR-10 |
| 空の query と空の subject | subject が先：終了コード 3 | LR-9 |

LOSAT 今：seekable な file は同じ。pipe の空・空白、両方 `-`、注釈だけの file は `an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's <PROGRAM>`（終了コード 1）などで拒否。移植後は零 batch の成功（「query の無い報告」）と終了コード 3 の経路を作り、`stream_position` が失敗する stream を空でないとする（`stream_is_empty`）。Windows の `IsIStreamEmpty`（`:848-854`）は空の pipe を空とする。oracle は Linux なので Linux の意味を全 platform で採る（LR-30、§L）。

### C2. 行の読み方（`CStreamLineReader`、`line_reader.cpp:77-297`）

| 項目 | NCBI | 行 |
|---|---|---|
| 状態 | stream ごとに 1 つ：行番号・押し戻した行の印・EOL の型・pushback の streambuf・eof/fail。record や batch で戻さない。subject と query の reader は別（ただし両方 `-` なら同じ `cin`） | LR-11 |
| 型の決定 | 最初の行だけで決める（`x_AdvanceEOLUnknown`、`:219-243`）：`NcbiGetline("\r\n")` の後、読み戻した区切りが CR（LF が続かない）なら CR 型、LF（CRLF の LF も）なら「CRLF」型（LF か CRLF かは遅延）、区切りの無い 1 行の file は未決 | LR-16 |
| CR 型 | CR で切る。直後の LF は食う（CRLF は 1 つの区切り） | LR-17 |
| LF・CRLF 型（`x_AdvanceEOLCRLF`） | LF で切り、直前の CR は落とす。行の途中の単独の CR はそこで切り、型を CR 型・mixed に切り替える。**その物理行の終わりは失われ**、尾は押し戻されて次の行に繋がる | LR-17 |
| mixed 型 | CR、LF、CRLF がそれぞれ 1 つの区切り。LFCR は 2 つ（行番号が倍になるが record は正しい） | LR-18 |
| 失う尾 | 押し戻しが既存の pushback buffer（256 byte 以下）に書き戻され、最後の改行の無い行が EOF で eofbit を立てていると、`AtEOF()` が真になり尾が返らない（`n_lost_g55`）。埋め込みの EOL が 1 つだけなら失わない（`n_lost_a..d`）。fill の大きさは file（残りの byte）と pipe・cin（FIONREAD か 1 byte）で違う | LR-20、LR-21 |
| `AtEOF` / `End()` | 押し戻した行が無く、`eof()` か `peek() == EOF`。`peek` が EOF を見つけると eofbit が残る。最後の record が後続の空行・注釈を読み飛ばすので、末尾の空 batch はできない（`BATCH_SIZE` が小さくても） | LR-12 |
| 行番号 | 物理行ごとに 1 つ進む（空行・注釈・定義行も）。`UngetLine` は戻す。EOF では進まない。警告の「line N」はこの番号 | LR-14、RD-45 |
| 行の長さ・byte | 制限なし（3 MB の行を確認）。BOM・NUL・0x1A・0xA0・0xFF は行の普通の byte。区切りは CR と LF だけ | LR-19、LR-22 |
| 最後の行 | 改行の無い最後の行も返る。LF で終わる file に余分な空行は無い。連続する空行はそれぞれ数える。LF の file の末尾の単独の CR は落ちる | LR-23 |

観察：LF・CRLF・CR だけ・それぞれ最後の改行なしは byte 一致。LF/CRLF/CR の混在 file では q2 の配列が題に繋がる（`Query_2 q2 secondTGGC..`、データなし、20 nt の題の警告）。`\r\r\n` は区切り 1 つと空行 1 つ。hyphen の警告の行番号は LF・CRLF・CR で 2,3,5,7,8、LFCR で 3,4,7,11,13（LR notes）。

LOSAT 今：行の読み込み器は無い（`bio` は LF で切り、末尾の CR を落とす）。CR だけの file は拒否。新しい `blastinput/fasta_reader/stream.rs` は BLASTX の `FastaStream` と同じ算術（`push_back`・`raw_peek`・`failed`）で、BLASTX の C++ oracle（`LOSAT/tests/unit/blastx_stage_e_io_stream_expected.tsv`、12 180 行）と Python の模型（`scratch_LR/lrmodel.py`、`file_` の 1413 行すべてと fuzz で一致）が基準（LR-11〜21、S0）。

## D. レコード（範囲 RD、`fasta.cpp`）

### D1. 1 行ごとの分岐（`ReadOneSeq`、`:312-440`）

| 先頭の形 | NCBI | 行 |
|---|---|---|
| 行頭が `>?_` | `>` に書き換えて普通の定義行として扱う（`>?_x` は題 `x` の新しい record。record の途中なら今の record を終え、元の行を `UngetLine`） | RD-3 |
| 行頭が `>?`（その他） | `UngetLine` して data 行として `ParseGapLine`（§D6） | RD-3、RD-36 |
| 行頭が `>` | 定義行（`need_defline` なら `ParseDefLine`、そうでなければ今の record を終えて `UngetLine`） | RD-3 |
| それ以外 | `TruncateSpaces_Unsafe`（C locale の `isspace`。U+00A0 は含まない）で両端を落とす。空の行は黙って飛ばし、`!`・`#`・`;` で始まる行は注釈として飛ばす（どこでも：最初の定義行の前、配列の途中、gap の後）。飛ばした行も行番号に数える。先頭に空白のある `>`（` >q1`）は data 行 | RD-4、RD-5 |
| 最初の定義行の前の data 行 | `need_defline` と `fDLOptional`：`ParseDefLine(">")`（題なし、ID を 1 つ使う）の後、同じ行を data 行として処理。配列だけの file は record `Query_1` 1 つ。BOM が配列の前なら 3 byte の無効な残基、`>q1` の前なら `CheckDataLine` の誤り | RD-5 |
| EOF（`need_defline` が真のまま） | 注釈・空行だけが残った：`CFastaReader: Expected defline around line N`（eEOF）。`CBlastInput` が握りつぶし、batch は空（`Empty CBlastQueryVector`、終了コード 3）。後続の注釈・空行は最後の record に属する | RD-7、BI-50 |

### D2. 定義行と題（`CFastaDeflineReader::ParseDefline`、`fasta_reader_utils.cpp:146-226`）

| 項目 | NCBI | 行 |
|---|---|---|
| 空の定義行 | `>` だけ、または `>` と空白だけ：題なし・ID なしで record は残る（`Query_N`。outfmt 0 は `Query= ` の後に `Length=`）。`>\x01` は空ではない | RD-10 |
| 題の始まり | `>` の後の `isspace` を飛ばす。`fNoParseID` なので残り全部が題（最初の語も ID にならない）。題の最初の byte は制御文字でも残る | RD-11 |
| 題の終わり | 2 byte 目から、最初の 0x20 未満の byte で切る（tab・CR・NUL・その他の制御文字。以降は黙って捨てる）。0x20〜0xFF（DEL と 0x80 以上を含む）は残る。末尾の空白は題の警告の判定までは残り、`x_ApplyMods` で落とす | RD-12 |
| 非 ASCII | UTF-8 の検査も正規化もしない。非 UTF-8 の byte も題に残り、報告にそのまま書く | RD-13 |
| `[key=value]` | `fAddMods` が無いので解析しない。題に括弧のまま残る（警告は無視される） | RD-15 |
| 題の警告 | `AssembleSeq` の `ParseTitle`（`:670-686`、`:1617-1679`）で **record の最後**、その record の data 行の警告の後に出る。核酸：題の長さが 20 byte を超え、最後の 20 byte がすべて `ACGTacgt`。蛋白：50 byte を超え、最後の 50 byte がすべて ASCII の英字。末尾の空白・数字・非 ASCII・`U`・`N` は判定を崩す。題の長さの判定は末尾の空白を落とす前 | RD-14、RD-16、RD-17 |
| 1000 文字を超える題 | 警告なし（`TooLong` は無視）、題はそのまま | RD-14 |

### D3. 残基（`ParseDataLine`、`:774-1016`）

トリムした行を前から 1 byte ずつ見る。

| byte | 核酸（`fAssumeNuc`） | 蛋白（`fAssumeProt`） | 行 |
|---|---|---|---|
| `A B C D G H K M N R S T U V W Y`（大小） | 残す（大文字で保存。`U` は `U` のまま Seq-data に入り、検索と表示は `T`。LOSAT は `U`→`T`、`u`→`t`） | 英字は全部残す（`U O J B Z X` も） | RD-22、RD-29 |
| `E F I J L O P Q Z X`（大小）と `*` | 無効な残基（取り除き、行ごとに 1 回の警告）。mask は開閉しない | `*` は残基（mask を閉じる）。英字は残す | RD-23、RD-29 |
| 数字・`.`・`@`・`#`・`!`・`>`（列 0 以外）・DEL・NUL・0x80 以上の byte・その他 | 無効な残基（1 byte ごとに 1 つ。U+00A0 は 2 つ、BOM は 3 つ） | 同じ | RD-24、RD-25 |
| `-` | 取り除く。mask は開閉しない。**警告は行ごとに 1 回**（`-` が何個あっても。無効な残基の警告の前） | 同じ | RD-27 |
| space・`\t`・`\v`・`\f`・`\r`・`\n` | 黙って飛ばす（警告なし、mask の run は続く）。後の無効な残基の「位置」には数える | 同じ | RD-26 |
| `;` | その行の残りを黙って無視（mask の run は次の行に続く） | 同じ | RD-28 |

警告の位置は、トリムした行の 1 始まりの byte の位置で、連続する位置は `a-b` に併合し、`, ` で繋ぐ。行ごとに範囲は最大 1000 で、超えた分は黙って切る（RD-33）。

### D4. `CheckDataLine`（`:710-772`）

record の最初の data 行（残基がまだ 1 つも無い間の、gap の行以外の各 data 行）で、トリムした行の最初の 70 byte を数える：英字と `*` が good、`-` は中立（`fHyphensIgnoreAndWarn`）、空白と数字は中立、`;` で打ち切り、それ以外（0x80 以上を含む）が bad。`bad >= good / 3 && (len_to_check > 3 || good == 0 || bad > good)`（整数の除算）なら誤り：

`CFastaReader: Near line N, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.`（`CObjReaderParseException` eFormat。stderr は `BLAST query error: ` + 本文、終了コード 1。subject でも同じ）。

例：`@@A` は誤り、`@A` は警告だけ、`ACGTACGT@@` は誤り（2 >= 2）、`AC;@@@@@@@` は誤り、`---`・`12345` は誤り（good == 0）、BOM + `>q1` は誤り。2 行目以降は検査しない（無効な残基の警告だけ）。核酸で 40% を超える曖昧な文字の警告は、`bIsNuc` が偽のため数が 0 のままで、無視される警告でもあり、出ない（RD-18〜20、RD-35）。

### D5. 小文字と mask（`x_OpenMask`・`x_CloseMask`、`:1062-1092`）

`-lcase_masking` のときだけ、残した小文字が `GetCurrentPos(ePosWithGapsAndSegs)`（残基の数 + それまでの gap の長さ）で区間を開き、次の**大文字**の残基（蛋白は `*` も）と record の終わりで `[開始, 位置-1]`（0 始まり）を閉じる。閉じないもの：空白、`-`、無効な byte（小文字の `e` も）、数字、`;` の後の文字、行の終わり、注釈・空行、gap の行。したがって「小文字、無効な byte、小文字」は 1 つの区間で、「小文字、無効な byte、大文字」は大文字で閉じる。gap の塩基は run に入る（`>?20` の前後が小文字なら 1 つの区間で `n` と表示。gap が小文字の run の前にあれば区間は gap の後から）。mask の無い実行でも小文字は大文字になる。subject も同じ reader なので同じ mask を持つ（表示だけに使う）（RD-30〜32、RD-39、RD-46）。

LOSAT 今：完成した `bio` の配列の小文字の最大の run から作る（受け付ける入力では同じ）。新しい record は区間と同じ場所を小文字、ほかを大文字にして持ち、`collect_lowercase_masks` はそのまま使える（port_plan §1）。

### D6. gap の行（`ParseGapLine`、`:1094-1349`、`AssembleSeq` の N/X への置き換え、`:1387-1432`）

行頭が `>?`（`>?_` 以外。先頭の空白があっても data 行として届く）。`remaining = TruncateSpaces(line[2:])`、大文字小文字を区別する接頭辞 `unk`（以降の空白は落とす）は未知の長さの印だが、`fParseGaps` が無いので常に数字の長さを使う。先頭の数字（`isdigit`）を `StringToUInt` で読み、空・0・4294967295 超は `CFastaReader: Bad gap size at line N` で長さ 1 とする。残りから空白を落とし、何か残れば `[key=value]` として解析し、`[` で始まらない・`=` か `]` が無いと `CFastaReader: Problem parsing gap mods at line N`。長さ N の gap は record の**今の位置**に N 個の `N`（核酸）か `X`（蛋白）を入れる（連続する行は足し算、record の先頭・末尾も可、長さの上限なし）。gap の行は `CheckDataLine` を受けず、残基と数えない。修飾子（`gap-type`・`linkage-evidence`）は残基を変えず、stderr の警告だけ（§E1）。

### D7. record の終わり、中身の無い record、ID

`AssembleSeq`（`:1351-1516`）は `CloseMask`、`AssignMolType`（旗のとおり）、`ParseTitle` の順。残基も gap も無ければ長さ 0 の record のまま返る（誤りなし、局所 ID と題は残る）。record に残る情報：題の byte、局所 ID（レコードの通し番号）、残基、mask、警告（行の番号つき、読んだ順）。`CheckDataLine` の誤りは例外で、それまでの record の警告は stderr に出済みで、その batch の結果は何も出ない。

## E. 文言（byte）と終了コード

### E1. 読み込みの警告（stderr、改行で終わる、`Warning:` も ID も無い。`LOG_POST_X(1, Warning << Message())`）

| 文言 | 出る条件・順 | 行 |
|---|---|---|
| `CFastaReader: Hyphens are invalid and will be ignored around line N` | `-` を含む data 行ごとに 1 回。その行の無効な残基の警告の前 | RD-27 |
| `FASTA-Reader: Ignoring invalid residues at position(s): On line N: 51-54, 62` | 無効な byte を含む data 行ごとに 1 回（「nucleotide」「protein」の語は入らない） | RD-33 |
| `FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?`（ピリオドの後に空白 2 つ）／`… at least 50 valid amino acid characters.  Was …` | record の最後（その record の data 行の警告の後）。行番号なし | RD-16、RD-17 |
| `CFastaReader: Bad gap size at line N`、`CFastaReader: Problem parsing gap mods at line N`、`Unknown gap-type: V`、`Unknown linkage-evidence: V`（値の全体）、`Unknown gap modifier name(s): K`、`There were conflicting gap-types around line N`、`FASTA-Reader: Unknown gap-type can have linkage-evidence of type 'unspecified' only.`、`CFastaReader: This gap-type should have at least one specified linkage-evidence.` | `>?` の行。修飾子は key の順に並ぶ。`This gap-type cannot have any linkage-evidence specified, so any will be ignored.` と `Gap mods are ignored because …` は無視される問題なので出ない。`gap-type` の値は `CSeq_gap::NameToGapTypeInfo` が `CanonicalizeString(sName).c_str()` で `const char*` の表から引く（`Seq_gap.cpp:157-172`、`Seq_gap.hpp:93`）ので、最初の NUL byte までで照合する（`telomere\0zz` は telomere、警告の `V` は値の全体）。key と `linkage-evidence` は byte 列の全体で比べる（SFc の監査 A-1 で移植） | RD-36、RD-37、RD-38 |

読み込みの警告は、record を読む時点で出る。subject の警告は引数の処理の中（最初）、query の警告は batch を読む時点：blastn の最初の batch（5 query）の後の record 6 の警告は q1〜q5 の報告の後、q6 の報告の前に出る。tblastx（10002 で 1 batch）は全報告の前（BI-47）。`cerr` は `cout` に結ばれているので、ある batch の報告は次の batch の警告と誤りの前に標準出力に届く。止まる batch（警告の後に読み込みの誤り）でも同じで、BLASTP の outfmt 6・7 は報告を書いた時点で flush する（`blast_input.cpp:134-146`。SFc の監査 A-3・C-1 で直した。outfmt 0 と他の program は前から同じ）。

### E2. 読み込みの誤り

| 文言（stderr） | 終了コード | 条件 | 行 |
|---|---|---|---|
| `BLAST query error: CFastaReader: Near line N, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.` | 1 | §D4。subject なら引数の処理の中、query ならその batch を読む時点 | RD-20、BI-51 |
| `BLAST engine error: Empty CBlastQueryVector` | 3 | batch に query が 1 つも無い（注釈・空白だけ、または全部が `-query_loc` などで飛ばされた）。subject が 0 個。最初の batch なら prolog の後、後の batch なら前の batch の報告（epilog なし）の後 | BI-29、BI-50、LR-8、LR-10 |

### E3. アプリの層の文言

| 文言 | 条件・時点・終了コード | 行 |
|---|---|---|
| `Warning: [prog] Query is Empty!` | `IsIStreamEmpty`、終了コード 0（§C1） | LR-9 |
| `Warning: [prog] Query_N[ <title>]: Sequence contains no data `（末尾の空白 1 つ） | 残基の無い query。その query の報告の直前。ID は `Query_N`（題があれば ` ` + 題）で、35 byte を超えると最初の 25 byte + `.. `。query は残り、報告は「ヒットなし」（outfmt 0：`Query= …`・`Length=0`・`***** No hits found *****`・`Effective search space used: 0`。outfmt 7：`# 0 hits found`、`Fields` 行なし。outfmt 6：何も出ない）。`# BLAST processed N queries` に数える | BI-40、BI-42 |
| `BLAST engine error: Warning: Sequence contains no data Warning: Sequence contains no data ` | batch の全 query が空。空の query ごとに 1 つ（ID なし）。終了コード 3、prolog の後。**batch ごと** | BI-41 |
| `Warning: [prog] Subject_N <title>: Subject sequence contains no data`（末尾の空白なし。題が空なら `Subject_N : …`） | 残基の無い subject。`CBlastFormat` の生成時。subject は database の統計に残り、ヒットしない | BI-43、RP-44 |
| `BLAST engine error: The average subject length is too short` | N・X だけ。全 subject が空。prolog の後、終了コード 3。P・T は終了コード 0 で全 query が「ヒットなし」 | BI-44 |
| `BLAST query/options error: Invalid from coordinate (greater than sequence length)` | `-subject_loc` が record の長さを越える（致命的）。query では黙って飛ばす | BI-26、E2d §B |

### E4. 終了コード

0：警告だけ、正常、`Query is Empty!`、零 batch の成功。1：CArgs・option・読み込み・範囲の誤り。3：engine の誤り（空の batch、全部が空、subject が短すぎる）。4：メモリ。5：network（data loader）。6：出力の失敗。255：その他の例外（`BATCH_SIZE` の文字列、NCBI の内部の例外）（BI-39）。

## F. 局所 ID・題・報告（範囲 RP、byte のまま）

| 出力 | NCBI | 行 |
|---|---|---|
| tabular の ID の既定の列（`qaccver saccver`）と `qseqid`・`sseqid`・`qacc`・`sacc`・`sallseqid`・`sallacc` | 題の最初の `' '`（0x20 だけ。NBSP・U+3000・0x0B・0x0C・BOM は語に残る）までの語（`s_ReplaceLocalId`、`tabular.cpp:474-504`）。題が空なら `Query_N`・`Subject_N`。byte をそのまま書く | RP-2〜6 |
| `sseqid`・`sacc`・`saccver` | `GetSeqIdList`（`showdefline.cpp:219-247`）：置き換えた ID が `lcl|Subject_` を含むとき（題が空、または最初の語が `Subject_` で始まる）、`GenerateDefline` の最初の語に置き換える。例：題 `Subject_&amp;x desc` → `Subject_&x`、`Subject_4.` → `Subject_4`。蛋白の subject で題が空なら `unnamed`（`sallseqid`・`sallacc` は `Subject_N` のまま）。query には掛からない。LOSAT 今は受け付ける入力でここが違う（`v01`・`v09`） | RP-17〜19 |
| 題が `Subject_` の語で非 UTF-8（0x81・0x8D・0x8F・0x90・0x9D を含む） | outfmt 6 で `Subject_\x81\t` を書いた後に `CCoreException` eNullPtr（build の path つき）、終了コード 255 | RP-20（§K） |
| `stitle`・`salltitles` / `qgi`・`sgi`・`sallgi` | FASTA の subject は `N/A` / `0` | RP-13、RP-15 |
| outfmt 7 `# Query: ` | ID なしで題の byte（`TruncateSpaces`、長さ制限・折り返しなし）。題が空なら `# Query: `。`# Database: User specified sequence set (Input: <-subject の引数のまま>)` | RP-21、RP-25 |
| outfmt 0 `Query= ` | 題を `NStr::Wrap` で 68 byte に折り返す（byte 単位で UTF-8 の途中も切れる）。題が空なら `Query= ` の直後に `\nLength=N`（空行なし）。HTML の復号も変換もしない | RP-26〜28 |
| outfmt 0 の説明の一覧 | `GenerateDefline(fLeavePrefixSuffix)`。ID は出さない。68 byte を超えると最初の 65 byte + `...`、足りなければ 68 byte まで空白で埋める | RP-29、RP-30 |
| outfmt 0 の見出し | `> ` + `GenerateDefline(0)`（TPA/MAG などの接頭辞を取る）を `s_WrapOutputLine` で 60 byte ごとに折り返し（次の空白の後で改行）、`\nLength=N`。蛋白で題が空なら `> unnamed protein product` | RP-32、RP-39、RP-40 |
| `GenerateDefline` | 末尾の `.,;~ ` を落とす（全部がそれなら落とさない）→ 接頭辞を取る（見出しだけ）→ `HtmlDecode`（実体、`&#N;`・`&#xH;` は 2^32 で巻き、検査なしで UTF-8 に書く）→ 末尾の `,;~ ` → `x_CleanAndCompress`（題の最後の byte が 0x80 以上なら**落ちる**：signed char） | RP-31〜38 |
| 非 UTF-8 の題 | `GuessEncoding`：構造が正しい UTF-8 ならそのまま、そうでなければ ISO-8859-1（0x80〜0x9F が無いとき）か Windows-1252（0x81・0x8D・0x8F・0x90・0x9D が無いとき）で 1 byte ずつ UTF-8 に。それ以外は例外：説明が `Unknown`、見出しなし、HSP ごとに stdout へ `Sequence with id Subject_N no longer exists in database...alignment skipped`、終了コード 0 | RP-33〜36、RP-41 |
| query の警告の ID・題 | `Query_N`（題があれば ` ` + 題）、35 byte を超えると 25 byte + `.. `。生の byte（制御文字・非 UTF-8 も） | RP-43 |

LOSAT 今：受け付ける入力（ASCII で制御文字なし、先頭の空白なし）では一致するが、IDは bio の最初の Unicode 空白までの語、題は `String`。拒否している入力：空・先頭が空白・制御文字・非 ASCII の定義行、実体のある subject の題、非 UTF-8。移植後は record が題の byte を持ち、`PairwiseHit.subject_title` は使わず 4 つの program が題を `&[Arc<[u8]>]`（`q_idx`・`s_idx` 順）で渡す。BLASTX の関数と `PairwiseHit` の型は変えない（RP-49〜52、port_plan §1、§2 S2）。

## G. data loader と Seq-id の経路（範囲 BI）

### G1. 読み手の選択（`DATA_LOADERS`）

`SDataLoaderConfig(is_prot)` は blastdb と genbank を有効で始め、`app->GetConfig()`（環境変数 `NCBI_CONFIG__BLAST__DATA_LOADERS` > `<prog>.ini` > `.ncbirc`。`ncbireg.cpp:1567-1581`）の `[BLAST] DATA_LOADERS` があれば `blast_scope_src.cpp:76-94` に従う：値に `blastdb`（大小無視の部分文字列）が無ければ blastdb を、`genbank` が無ければ genbank を切り、`none` を含めば両方切る。`UseDataLoaders() = blastdb || genbank`。したがって**空でない値で `blastdb` も `genbank` も含まないもの（例：`x`）も両方切る**。既定（項目なし）は query・subject とも有効（BI-4、BI-54）。**空の値 `DATA_LOADERS=`（空白だけ・`""` も同じ）は、その項目を持つ層では項目がある（両方の loader を切る）**：`CCompoundRegistry::FindByContents` は各層に `fCountCleared` を付けて尋ね（`ncbireg.cpp:1235-1246`）、最初に「ある」と答えた層の値を使う。`fCountCleared` があると `CMemoryRegistry::x_HasEntry`（`<prog>.ini` の層）も `CCompoundRWRegistry::x_HasEntry`（`.ncbirc` の層。空の値は cleared として記録される、`:2041-2045`）も空の値を数える（`:984-991`、`:1939-1947`）。**ただし `.ncbirc` の空の値は、既定ではアプリの registry に届かない**：BLAST の各 program は `CBlastUsageReport` を member に持ち（`blastn_app.cpp:81`）、その構築子は環境変数 `BLAST_USAGE_REPORT`（環境だけ、`getenv`）が偽の Boolean でなければ `fWithNcbirc` の registry を作り、`.ncbirc` を `CMetaRegistry` の cache に読む（`blast_usage_report.cpp:195-211`、`ncbireg.cpp:1622-1625`・`:1650-1653`）。アプリは `<prog>.ini` が無ければ同じ鍵で `.ncbirc` を読み、cache の registry を `Write` と `Read` で自分の registry に写す（`ncbiapp.cpp:1278-1279`、`metareg.cpp:152-171`・`:200-208`）。`Write` は `fCountCleared` なしで項目を数えるので、空の値は写らない（`ncbireg.cpp:226-234`・`:1997-2004`。他の値は `Printable` と `ParseEscapes` が逆なので変わらない）。`<prog>.ini` があると、`.ncbirc` は `CNcbiRegistry::x_Read` から `fJustCore` を付けて（`SEntry::Reload`、`metareg.cpp:89-90`、`ncbireg.cpp:1688-1691`）cache に無い鍵で読まれ、アプリの registry に直接入る。したがって `.ncbirc` の空の値が項目になる（loader を切る）のは、`BLAST_USAGE_REPORT` が偽（`0`・`false`・`f`・`no`・`n`・`off`、大小無視。`NStr::StringToBool`）のときか `<prog>.ini` があるときで、それ以外では項目なし（次の層が決め、無ければ有効。SFb の oracle の観察はこの場合）。`NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT` や `.ncbirc` 自身の `BLAST_USAGE_REPORT` は usage report の読み込みの後に読まれるので効かない。以前のこの節は `ncbireg.cpp:984-991` を `fCountCleared` なしで読み、file の空の値を一律に項目なしとしていた（SFd の再監査 A-1・A-2・B-3 が oracle で見つけた）。**探索 path**（`metareg.cpp:331-396`、`ncbi_environment.rs` の `registry_search_path`）：`NCBI_CONFIG_PATH` が空の値で設定されていると、`NStr::Split` は token を作らず（`ncbistr_util.hpp:268-273`）、探索 path は空になり、`<prog>.ini` も `.ncbirc` も読まれない（再監査 A-3・B-4）。`NCBI_CONFIG_PATH`・`NCBI`・`HOME` の値は byte のまま使う（UTF-8 でない directory の名前も探す。再監査の第 2 回 R-3）。home は `CDir::GetHome`（`ncbifile.cpp:3586-3657`）：`HOME` があればその値（空なら home なし）、無ければ passwd の利用者の項目の directory（`getpwuid(getuid())`。空なら home なし。R-2）。passwd に利用者の項目が無いと NCBI は login 名（`USER`・`LOGNAME`・`getlogin()`）の項目を `getpwnam` で引く。LOSAT は passwd の項目を Rust の `std::env::home_dir`（`getpwuid_r` を `sysconf(_SC_GETPW_R_SIZE_MAX)` の大きさ、glibc で 1024 byte の buffer で 1 回だけ呼ぶ）で読むので、項目が無いときも、項目がその buffer に収まらないとき（NCBI の `getpwuid` は読める）も home が分からず、探索がその home の位置に届くときだけ拒否する（§J-8。再監査の第 3 回 S-3。この機械では nss-systemd が root の項目を足すので、`passwd: files` だけの `nsswitch.conf` を mount 名前空間に置いて確かめた）。**program の directory**（R-1）：探索 path は `CMetaRegistry` を最初に使うときに一度だけ作られる（`metareg.hpp:257-261`）。既定ではそれは usage report が `.ncbirc` を読むとき（アプリの構築中、`AppMain` が引数を渡す前）で、program の名前はまだ無く（`ncbi`）、Linux では `/proc/<pid>/exe` の link を解いた directory だけを探す（`ncbienv.cpp:373-393`。macOS・Windows では program の directory を探さない）。usage report が `.ncbirc` を読まないとき（環境変数 `BLAST_USAGE_REPORT` が偽、`NCBI_DONT_USE_NCBIRC` か `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` が在る）は `LoadConfig` で作られ、起動したときの名前の directory が先、link を解いた directory がそれと違えば次に入る（`metareg.cpp:375-388`）。起動したときの名前は `FindProgramExecutablePath`（`ncbiapp.cpp:1426-1609`）：`argv[0]` が絶対ならそのまま、作業 directory から見て file なら（末尾の `/` は落として見る。`CDirEntry::Reset`、`ncbifile.cpp:298-313`。再監査の第 3 回 S-5）その前に作業 directory を付け、そうでなければアプリの環境（下の `CNcbiEnvironment`）の `PATH` の各 directory で base 名（最後の拡張子を除く）を探し、link を辿らず文字の上で正規化する（`NormalizePath`、`ncbifile.cpp:820-971`）。**`<prog>.ini` の名前**：NCBI ではこの名前の base 名（無ければ link を解いた名前の base 名、`ncbiapp.cpp:1243-1258`）で、`argv[0]` が空なら `ncbi`（`ncbi.ini` を先に読む、`ncbiapp.cpp:886-909`）。LOSAT は NCBI の program が自分の名前（`blastn`・`blastp`・`tblastn`・`tblastx`）で起動された場合を写す：`<prog>.ini` は LOSAT の program 名で、LOSAT の `argv[0]` は directory だけを与える。別の名前で起動した NCBI（`foo` なら `foo.ini`）はこの範囲の外。LOSAT の `argv[0]` が空なら写す NCBI の起動が無いので拒否する（§J-8。S-2）。directory は名前の最後の `/`・`\`・`:` までの部分（`ncbienv.cpp:408-415`）。したがって symlink から起動した program では、symlink の側の directory の registry の file を読むかどうかがこの 3 つの環境変数で決まる（LOSAT も同じ。以前は解決済みの directory だけを探していた）。`NCBI_USAGE_REPORT_ENABLED`・`DO_NOT_TRACK` の有効な値は usage report の `.ncbirc` の読み込みを変えない（`x_CheckBlastUsageEnv` はどちらも読まない。R-8 の確認、oracle で一致）。registry の file の読み方（C locale の白空白、名前の文字、引用符、escape、継続行、BOM、構文の誤り）は `ncbi_environment.rs` の `registry_entries` が `IRWRegistry::x_Read` を移す（§J-8）。BOM の検査（`GetTextEncodingForm`）は `.ncbirc` で 1 回、`<prog>.ini` で 2 回（`IRWRegistry::Read` を 2 度通る、`ncbireg.cpp:605-627`・`:1681-1692`。R-4）。探索で選ばれた file が開けないと registry は無く、探索も先へ進まない（`SEntry::Reload`、`metareg.cpp:62-84`・`:228-231`。開けない `<prog>.ini` は無いのと同じで、`.ncbirc` は cache の写しになる。R-6）。**同じ名前の環境変数が複数あるとき**（再監査の第 3 回 S-1）：NCBI の `getenv` の読み手（`BLAST_USAGE_REPORT`、`NCBI_CONFIG_OVERRIDES`、`NCBI_CONFIG_PATH`、`NCBI_DONT_USE_LOCAL_CONFIG`、`NCBI_DONT_USE_NCBIRC`、`NCBI`、`HOME`）は最初の項目を、アプリの `CNcbiEnvironment`（registry の環境変数の層 `NCBI_CONFIG__…`、`FindProgramExecutablePath` の `PATH`）は最後の項目を使う（`Reset` が `environ` の順に map に入れる、`ncbienv.cpp:87-105`）。LOSAT も同じに読み、値で決まる拒否（`BLAST_USAGE_REPORT`・`NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT` の Boolean、空の `NCBI_CONFIG_OVERRIDES`）は NCBI が使う方の値だけを見る。**環境変数は空の値でも項目がある**：`CEnvironmentRegistry::x_HasEntry`（`env_reg.cpp:157-167`）は変数が在れば `found` を返し（`CNcbiEnvironment` は `NAME=` も値の在る変数として持つ、`ncbienv.cpp:97-104,117-124`）、空の `NCBI_CONFIG__BLAST__DATA_LOADERS=` は最上位の層として答えて blastdb も genbank も切る（`.ncbirc` の `DATA_LOADERS=genbank` も上書きする）。以前のこの節は file の規則を環境変数の層にも当てはめていた。SFc の監査 B-1 が oracle で見つけ、`ncbi_environment.rs` の `data_loaders_of` を直した。`BLASTDB_NUCL_DATA_LOADER`・`BLASTDB_PROT_DATA_LOADER` は blastdb の loader が開く DB の名前だけを決め（`:98-125`）、Seq-id を拒否する限り出力に効かない（BI-5）。`DATA_LOADERS` だけは、Seq-id の行を配列として読むか拒否するかを決めるので、LOSAT の `ncbi_environment.rs` の「出力を変えない」扱いを改める（BI-54）。

### G2. 試される行（`CBlastInputReader::ReadOneSeq`、`blast_fasta_input.cpp:130-161`）

読み手が `CBlastInputReader` のとき、**`ReadOneSeq` のたびに**、次の 1 行を `TruncateSpaces_Unsafe` して読み、空でなく先頭の byte が `isalnum`（`[A-Za-z0-9]`。BOM・`>`・`;`・`#`・`!`・非 ASCII は試さない）なら `CSeq_id(line, fParse_AnyRaw | fParse_ValidLocal)` を作り、局所 ID になった場合は `lcl|` で始まらなければ既定の旗で作り直す（ここで `Malformatted ID` が投げられる）。`CSeqIdException` で文言に `Malformatted ID` を含むものは FASTA に戻り（`UngetLine`）、それ以外の `CSeqIdException`・`std::exception` は投げ直し、`catch (...)` も FASTA に戻る。試される行は、stream の最初の行か、Seq-id の行の直後の行だけ（FASTA の record の後は次が `>` か EOF）。Seq-id の行はちょうど 1 行を消費し、`Query_N` の番号は進めない。空行・注釈・BOM が先にあると試されない（BI-13、BI-14、BI-20）。

### G3. `CSeq_id` の解析で Seq-id になる形（`Seq_id.cpp:2457-2551`）

| 順 | 形 | 行 |
|---|---|---|
| 1 | tag：3 byte 目が `\|` で長さ > 3 なら先頭 2 byte、4 byte 目が `\|` で長さ > 4 なら先頭 3 byte が、`gb gi sp tr`／`bbm bbs dbj emb gim gnl gpp lcl nat pat pdb pgp pir prf ref tpd tpe tpg`（大小無視）のどれか。tag が合えば `Malformatted ID` は決して出ず、成功すれば取り寄せ、他の例外なら投げ直し | BI-15 |
| 2 | 生の accession（`IdentifyAccession`）。`.` の前を大文字にして、数字だけ（先頭が `0` でなく版なし）なら GI。数字で始まる PDB・PRF の形、`[A-Z]\d[A-Z0-9]{3}\d`（6 byte、先頭 O/P/Q または 3 byte 目が英字）か先頭が O/P/Q でない 10 byte の Swiss-Prot の形。英字（`_` 可）+ 数字の形は accession guide（`accguide2.inc`）の 26 の（英字数, 数字数）— 1+5、2+6、2+8、2+10、3+5、3+6、3+7、3+8、3+9、3+11、4+5、4+6、4+8、4+9、4+10、5+6、5+7、6+9、6+10、6+11、7+8、7+9、7+10、9+9、9+10、9+11 — のとき、接頭辞が約 1450 の規則に合えば accession | BI-16、BI-17 |
| 3 | `DB:` の接頭辞（`:` の前を大文字にして）が `ATGC BCMHGSC BERKELEY CELERA GSDB HOOD LANLCHGS LRG MIPS NCBI_EXT_ACC NCBI_GENOMES NCBI_MITO PGEC PID SGD SHGC SRA TIGR UOKNOR UWGC WASHU WIBR WUGSC`（`dbGSS`・`dbSTS` は大文字化した後なので一致しない） | BI-18 |
| 4 | それ以外は局所 ID か `Malformatted ID`：**FASTA の data として読む**。`ACGTACGT`、`ACGT ACGT`、40 nt の行、`MKVLAAGIVGLLLAQ`、`0123`、`12.5`、`XP_123`、`contig1`、`ACGT1234`（4+4 は guide に無い）、`gb\|`、`xx\|abc`。配列だけの file は record `Query_1`（題は空） | BI-19 |

取り寄せの結果（oracle `net1`〜`net6`、blastn と blastp。ネットワークに出た 6 回だけ）：存在する accession（`AB123456`）は `AB123456.1` として GenBank の題で検索される（`# Query: AB123456.1 Verrucosispora sp. 431C09 gene for 16S rRNA, partial sequence`）。存在しない accession（`AB000000`）の query は**黙って飛ばされる**（`CInputException` eSeqIdNotFound を `GetNextSeqBatch` が `continue` で握る。数にも入らない）が、subject なら致命的：`BLAST query/options error: Sequence ID not found: 'dbj|AB000000|'` + `Please refer to the BLAST+ user manual.`、終了コード 1。分子の型が違う accession（blastp に核酸、blastn に蛋白）も飛ばされる。`lcl|foo` は tag の経路で、見つからず飛ばされる（BI notes §0、BI-21、BI-52）。

### G4. LOSAT の判定の規則（決定済み：`$TASK_DIR/DECISIONS.md` の 2026-10-08、保守者の判断 2）

loader が有効（G1）のとき、stream の最初の行と各 Seq-id の行の次の行について、トリムした行が `[A-Za-z0-9]` で始まり、(a) G3-1 の tag、(b) G3-2 の表の要らない形（GI・PDB・PRF・Swiss-Prot）、(c) 英字と数字の形で（英字数, 数字数）が 26 の guide の形のどれか（`.数字` や S/P の旗の形を含む）、(d) G3-3 の `DB:` のどれかなら、明示的に拒否する（`… may be a sequence identifier that NCBI BLAST+ fetches through a data loader, which is not supported by LOSAT's <PROGRAM>`。`DATA_LOADERS=none` を案内する）。それ以外（G3-4）は FASTA の data として読む。(c) だけは表なしでは忠実に移せず、guide が未認識の型に写す接頭辞（例：`ZZ123456`）と、guide の形の核酸の配列行（例：最初の行 `AB123456` は NCBI では本物の accession）を広めに拒否する。研究の FASTA の最初の行はこの形ではないので費用は無い。判定は読む時点（遅延）で行い、subject にも同じ規則と文言を使う（NCBI も致命的）。`scratch_BI/seqid_probe.py` が規則を実装し、6 回の network の観察と `DATA_LOADERS=none` の実行を再現する（BI-55）。loader が無効なら（`none`・`blastdb`/`genbank` を含まない値）NCBI と同じく全部 FASTA として読む（`ncbi_environment.rs` が `DATA_LOADERS` の実効値を `ApplicationSettings` で返す。環境変数 > `<prog>.ini` > `.ncbirc`）。現在の `seq_id.rs` の暫定の規則（英字だけの行以外は拒否）は、NCBI が data として読む `0123`・`ACGT ACGT`・`contig1` も拒否するので、この規則で置き換える。

## H. 中身の無いレコードと batch の端の場合（BI、RD、AD）

| 場合 | NCBI | LOSAT 今 → 移植後 | 行 |
|---|---|---|---|
| 残基の無い query | record は batch に残り、警告（§E3）→ 「ヒットなし」の報告。`size_read` には 0 を足す | 全部拒否 → 4 つの program で NCBI と同じ | BI-40、BI-42 |
| batch の全 query が空 | `BLAST engine error: Warning: …`（§E3）、終了コード 3、batch ごと。前の batch の報告は完全 | 拒否 → 移植 | BI-41 |
| 残基の無い subject | 警告を prolog の前に。database の統計と有効な検索空間には数える。局所 ID の番号を使う。ヒットしない | BLASTP だけ警告の後で拒否 → 移植 | BI-43 |
| 全 subject が空 | N・X：`The average subject length is too short`（3）。P・T：終了コード 0 | 拒否 → program ごと | BI-44 |
| `-query_loc` が範囲外の query | 黙って飛ばす。局所 ID の番号は使う。`# BLAST processed` に数えない。最後の batch が全部飛ばされると `Empty CBlastQueryVector`（3） | 忠実（`cut_queries`、E2d） | BI-26、BI-52 |
| 範囲の始まりが record の長さちょうど | 文字の無い区間。query は警告つきで残り、subject は警告つきで飛ばされる | 拒否のまま（E2d 判断 R2） | BI-27 |
| 注釈だけの file（pipe でも seekable でも） | 空の batch：`Empty CBlastQueryVector`（3） | 拒否 → 移植 | BI-50 |
| 2 つ目以降の batch の読み込みの誤り | 前の batch の報告は出た後。epilog なし、終了コード 1 | 何も出さず拒否 → 遅延読みで一致 | BI-51、RD-49 |

## I. アダプタ：`scan` の種類 1・2、`register`、`run`（範囲 AD、port_plan §1、S9〜S10）

NCBI の読み方は program と役割だけで決まる（AD-9、AD-21）。`scan` の解析器の種類を次のとおり足す（`docs/web/abi_v2.md` §4 の「1 は BLASTX、SX」を書き直す）：

| 種類 | 旗 | 使う program・役割 |
|---|---|---|
| 0 | `bio::io::fasta` 1.6.0 の再現（残す。アプリが切り替えるまで使う） | 今のアプリの 4 program |
| 1 | `fAssumeNuc` | BLASTN の query・subject、TBLASTX の query・subject、TBLASTN の subject |
| 2 | `fAssumeProt` | BLASTP の query・subject、TBLASTN の query |

（BLASTX の query は 1、subject は 2 になるが SX で行う。AD-9）

- **scan**：chunk の大きさによらずエンジンの読み込み器と同じ record・ID・長さを出す。byte 単位（UTF-8 の誤りなし）。行の分け方は §C2 と同じ（EOL の型、CR の先読み、落とす LF の印、行番号）で、chunk の終わりの CR は次の byte か EOF で決める（AD-1、AD-12）。record の始まりは行頭の `>`（先頭の空白は data 行）。`>?_x` は定義行、`>?` で始まる他の行は**拒否**（§J-3）。最初の定義行の前の data は定義行の無い最初の record で、`header_offset == sequence_offset`（アプリが連結のとき `>` を足す。S09 の判断 Q5）。空行・注釈だけの file は record `[]`（旧 `Expected > at record start.` から変わる）。残基の無い record も record として出す（`length` 0、`end_offset` は次の定義行）。
- **`line_layout` の前へ読む規則**（`abi_v2.md` §9 に書く）：checkpoint から前へ読むとき、空白（space・`\t`・`\v`・`\f`・`\r`・`\n`）、`-`、`;` 以降（その行の残り）、無効な残基、注釈の行（トリムして最初の byte が `!`・`#`・`;`）、行の終わりを飛ばす。`residue_counts` は残した残基を**大文字にした**形（`U`→`T` の後）の byte を数え、エンジンの record も同じ数を返す（`check_scan` の網）。`CheckDataLine` は残基のできるまでの各行に効くので、行の先頭 70 byte の窓を持つ（AD-6）。
- **register**：program と役割の読み込み器（loader は CLI の既定と同じく有効：Seq-id の行を拒否）で読み、同じ種類の `scan` と照合する。JSON の `id` は題の最初の `' '` までの語（題が無ければ空。`Query_N` は入れない）で、非 UTF-8 の byte は JSON でだけ `U+FFFD` にする（出力の stream 0/6/7/3 は byte のまま。`abi_v2.md` の stream 3「UTF-8」を直す）。`q_idx`・`s_idx` は残基の無い record も数え、`-query_loc` で飛ばされた record があっても登録順の index に写す（`ordinals`）。読み込みの誤りと注釈だけの query は register で早く失敗し、NCBI の文言を使う。失敗した run は stream 3 を送らないので、後の batch の誤りの前の警告は見えない（成功する run だけが CLI と一致）。reader の警告は record ごとに持ち、run のたびに NCBI の時点で再生する（AD-14〜19）。
- **register が今拒否するもの**（AD-15、AD-16）：空の定義行、先頭の空白、制御文字・非 ASCII・非 UTF-8 の定義行、非 IUPAC の残基、残基の無い record、最初の定義行の前の text は、すべて NCBI の読み方で record になる。残る拒否は `CheckDataLine` の誤り、`>?` の行（gap）、Seq-id の行、record の長さの上限、I/O だけ。
- **ABI v1**：変えない（TD-1）。v1 は今の検査を先に行い、通った入力を `from_bio` で新しい record にして `run_local` に渡す（§J-5）。v1 の BLASTN は `bio` と NCBI が同じ。v1 の TBLASTX・BLASTP は定義行・配列行の検査を持たず（AD-24、AD-25）、新しい読み込み器にすると受け付ける入力の出力が変わる（空・tab・非 ASCII・先頭の空白の定義行、`>?` の gap）。それを避けるための `from_bio` である。ただし新しい報告の層で出力が S11 とも NCBI とも違う byte になる定義行は、SFd の再監査 D-1 とその第 2 回から検索の後で拒否する（§J-5）。
- **性質試験**：`web/adapter/tests/scan_ncbi_properties.rs`（新しい種類。oracle はエンジンの読み込み器、record の始まりと残基の byte 位置は独立した基準の走査）。生成する入力と端の場合は指示書の手順 4 のとおり（LF・CRLF・CR だけ・混在、最後の改行なし、注釈、定義行の無い最初の record、`>?` の行、`-`・`;`・空白・無効な残基・非 ASCII・非 UTF-8、空の record、BOM、checkpoint を越える長さ、任意の chunk）。

## J. LOSAT が残す明示的な拒否

文言は `not supported by LOSAT's <PROGRAM>` を含む。

| # | 拒否 | 範囲 | 理由 | 行 |
|---|---|---|---|---|
| 1 | loader が有効なとき Seq-id として試される最初の行（§G4） | 4 program、query・subject、CLI・Web | NCBI の結果が GenBank か BLAST DB と network に依る。同じ入力が機械ごとに違う出力になる（BI §9-6）。LOSAT は network も BLAST DB も使わない。保守者の判断 2。読む時点（遅延）に出す。`DATA_LOADERS=none`（実効値）なら NCBI と同じく FASTA として読む | BI-13〜23、BI-55 |
| 2 | `-parse_deflines`（`cli.rs:322,386,593,649,713`） | 全 program | 定義行の最初の語を Seq-id として解析する（局所 ID・`Could not construct seq-id`・data loader・題の文言が変わり、subject の誤りの文言に NUL が入る）。範囲の外 | RD-48、RP-7、BI-9 |
| 3 | `>?` の gap の行（`>?_x` は除く） | ABI v2 の `register` と `scan` だけ | gap の残基には元の byte が無く、索引（`line_layout`）から切り出せない。CLI と `run_local` は gap の行を移植する。保守者の判断 3 | AD-7 |
| 4 | 2^31 - 1（2147483647）残基を超える record（gap の長さを含む） | 4 program | NCBI の位置は `TSeqPos`（32 bit 符号なし、`ncbimisc.hpp:879`）、BLAST の配列の長さは `Int4`（`blast_def.h:246`）で、`size_read` も 32 bit で巻く（BI-29）。この長さの結果は桁あふれに依り、検証できる値が定義できない。`reader.rs` の `check_length`（`a <role> record longer than 2147483647 letters is not supported by LOSAT's <PROGRAM>`） | BI-29 |
| 5 | ABI v1 は `from_bio` を通る | v1 の BLASTN・TBLASTX・BLASTP | 計画 TD-1 の凍結（gbdraw が使う）。今の検査と文言と順を変えない。通った入力は新しい読み込み器でなく `bio` の record を `from_bio` で包む。保守者の判断 4 の前提（「受け付ける入力は `bio` と NCBI で同じ読み方」）が TBLASTX・BLASTP の v1 では成り立たない（AD-24、AD-25）ので、文字どおり凍結を採る（`DECISIONS.md` 2026-10-08）。v1 の BLASTP・TBLASTX は、新しい報告の層で出力が S11 とも NCBI とも違う byte になる定義行を検索の後で拒否する（`check_bio_deflines_of`、README 判断 12・26・28）：`bio` が ID を空にし `from_bio` の record が NCBI の読み込み器と違うもの（非 ASCII の白空白、行の途中の単独の CR、行頭の白空白の後の制御文字。題の無い record は名前の出るところだけ）、`>?` の行の後の題の無い record で番号の出るところ、BLASTP の outfmt 0 の subject で `bio` の読み方が NCBI と違い題が非 ASCII で終わるもの | AD-23〜27、SFd 再監査 D-1 |
| 6 | 題が `Subject_` の語で非 UTF-8 の byte を含む subject の tabular（RP-20） | N・X・T・P | NCBI は部分的な行と build の path つきの例外で終了コード 255。byte 一致で再現できない。§K の一括の問いの回答待ちの間は拒否 | RP-20 |
| 7 | そのほか既存の拒否で読み込みと別の事柄 | — | `-lcase_masking`（BLASTP・TBLASTX。`cli.rs:375,637`）、N・X・T の tabular の custom の field（RP-9〜11）、`-html`・`-show_gis`・`-line_length`・`-num_descriptions`（RP-42）、`BATCH_SIZE` の非整数・`BL2SEQ_LEGACY`（BI-7、BI-32）、UTF-8 でない file 名 | 各行 |
| 8 | NCBI が構文の誤りとする registry の file（`.ncbirc`・`<prog>.ini`）、UTF-16 の registry の file、開けたが読めない file、home の分からない探索 | 4 program、CLI | NCBI の読み手（`IRWRegistry::x_Read`、`ncbireg.cpp:652-781`）は次の行で `CRegistryException` を投げる：`]` の無い section の行、名前の無い section、section・項目の名前に `[A-Za-z0-9_./-]` 以外の byte（行頭・名前の中の NBSP・U+3000・U+0085、名前の中の空白・NUL）、`=` の無い行、値の中ほどの引用符、`NStr::ParseEscapes` が拒む `\`（末尾の `\`、桁の無い `\x`、`kMax_UInt` を超える `\x`）。`.ncbirc` では `Critical: (110.6) CNcbiRegistry: Syntax error in system-wide configuration file: NCBI C++ Exception:` と、NCBI の build の source の path と行を含む例外の文を、file を読む者ごとに出す（usage report とアプリで 2 回、片方だけなら 1 回）。誤りの行より前の項目だけを持って続ける。`<prog>.ini` では `Error: (CRegistryException::…) …` と `Error: (106.15) Application's initialization failed …` を 2 回ずつ出して終了コード 2。文言が NCBI の build に依り byte で再現できないので、LOSAT は file を拒否する（`the registry file <path> has <理由> on line <n>, which NCBI BLAST+ reports as a syntax error; this is not supported by LOSAT`、終了コード 1。program 名を含まないのは他の registry・環境変数の拒否と同じ）。範囲は NCBI がその file を読む場合だけ（探索 path は §G1。symlink から起動したときの directory と passwd の home を含む。`<prog>.ini` は NCBI の program が自分の名前で起動された場合の名前、§G1）：`NCBI_CONFIG_PATH` の空の値・`NCBI_DONT_USE_NCBIRC`・`NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` で読まれない file は調べない。`<prog>.ini` の `[NCBI] DONT_USE_NCBIRC` でアプリが読まない `.ncbirc` も、`BLAST_USAGE_REPORT` が偽でなければ usage report が読むので、構文の誤りと `[BLAST] BLAST_USAGE_REPORT` の Boolean でない値（空を含む。NCBI は `StringToBool` の例外で abort）を拒否する（他の項目は調べない）。BOM の検査は NCBI の `GetTextEncodingForm`（`ncbistre.cpp:782-826`）のとおりで、`.ncbirc` は 1 回、`<prog>.ini` は 2 回検査する（`ncbireg.cpp:605-627`・`:1681-1692`）：先頭が `EF`・`FE`・`FF` のどれかで次が `BB BF` なら UTF-8 の BOM として読み飛ばす（`<prog>.ini` は 2 つまで）。検査が UTF-16 の byte order mark（`FF FE`・`FE FF`）を見つける file（`<prog>.ini` では UTF-8 の BOM の後の UTF-16 の BOM も）は NCBI が UTF-8 に変換して読む（`ReadIntoUtf8`、`ncbireg.cpp:617-623`）が、LOSAT は変換を移さず拒否する（`… is UTF-16, which NCBI BLAST+ converts to UTF-8 before reading it; this is not supported by LOSAT`）。`EF`・`FE`・`FF` の 1 byte だけの file（`<prog>.ini` は UTF-8 の BOM の後の 1 byte も）は stream が失敗した状態で空として読まれ、NCBI は `Error: (110.4) Error reading the registry after line 1: ` の 1 行を file を読む者ごとに出して続ける（`ncbireg.cpp:813-816`。1 byte だけの `<prog>.ini` と、ちょうど 2 byte の `EF`・`FE`・`FF` + `BB` の `<prog>.ini` は何も出さずに空）。文言が固定で回数が読む者（§G1 の usage report とアプリ）で決まるので、LOSAT は同じ行を同じ回数出す（拒否ではない。再監査の第 2 回 R-5）。開けない file は NCBI と同じく registry なし（`SEntry::Reload`、`metareg.cpp:62-84`。R-6）。開けたが読めない file は拒否（`NCBI BLAST+ would read the registry file <path>, which LOSAT could open but not read (…); this is not supported by LOSAT`）。section より前の項目は NCBI と同じく読んで捨てる（`IRWRegistry::Set`、`:833-838`）。**home の分からない探索**：`HOME` が無く、`std::env::home_dir` が passwd の項目を返さないとき（利用者の項目が無い：NCBI は login 名 `USER`・`LOGNAME`・`getlogin()` の home を `getpwnam` で引く、`ncbifile.cpp:3599-3619`。または項目が `getpwuid_r` の buffer（glibc で 1024 byte）に収まらない：NCBI の `getpwuid` は読む。再監査の第 3 回 S-3）、LOSAT は home を決められず、registry の file の探索がその home の位置に届くときだけ拒否する（`HOME is not set and LOSAT could not determine the home directory where NCBI BLAST+ would look for its registry file <name> (std::env::home_dir found no passwd entry for the user: there is none, or it does not fit its getpwuid_r buffer; NCBI uses getpwuid, then the login name's entry); this is not supported by LOSAT`）。`NCBI_DONT_USE_LOCAL_CONFIG` があるとき、`NCBI_CONFIG_PATH` に空の要素が無いとき、home より前（`NCBI_CONFIG_PATH` の前の部分、`.`）で file が見つかるときは拒否しない。**空の `argv[0]`**：NCBI は自分を `ncbi` と名付け（`ncbi.ini` を `<prog>.ini` より先に読み、作業 directory と `PATH` で `ncbi` を探す）、LOSAT が写す起動（§G1）に当たらないので、LOSAT は常に拒否する（`the program was started with an empty argv[0]; NCBI BLAST+ then names itself ncbi (…), which is not supported by LOSAT`。`execve` で `argv[0]` を空にしたときだけ起きる） | SFd 再監査 A-4・B-2、第 2 回 R-2・R-4〜R-6、第 3 回 S-2・S-3 |

## K. NCBI の不具合に見える挙動（`PD-LOSAT-NCBI-DEFECTS` の一括の問いの候補）

規則（同文書の Rule）：1 NCBI が失敗して LOSAT が有効な結果を出せる＝承認済みの例外、2 NCBI の決まった結果は byte で再現、3 NCBI が失敗して有効な結果を定義できない＝明示的な拒否、4 NCBI が決まって受け付けるなら移植。以下は推奨を付けた候補で、保守者の答えを待つ。

| # | 挙動 | NCBI の結果 | 推奨 | 行 |
|---|---|---|---|---|
| 1 | tabular の `Subject_` の最初の語 + 非 UTF-8（0x81・0x8D・0x8F・0x90・0x9D）の題：`GenerateDefline` が投げ、`SetFields` が捕まえない | `Subject_\x81\t` を書いた後、`CCoreException` eNullPtr と build の path つきの stack trace、終了コード 255（outfmt 6 を観察。7・custom は同じ code と読んだ） | 規則 3（拒否）。堅い結果ではなく部分的な行 | RP-20 |
| 2 | 同じ題の outfmt 0：説明が `Unknown`、見出しなし、HSP ごとに `Sequence with id Subject_N no longer exists in database...alignment skipped`、終了コード 0 | 決まった結果。誤解を招く文言 | 規則 2（再現） | RP-36、RP-41 |
| 3 | LF/CRLF の file の途中の単独の CR（または CR の file の途中の LF）が行を切り、**その行の終わりを失う**：次の行が繋がり、題に配列が入って record が空になる | 決まった結果（混在の file のすべて） | 規則 2（再現）。BLASTX の移植と同じ | LR-17 |
| 4 | pushback の eofbit が残り、最後の改行の無い行を読み終えた後に、まれに最後の行の尾が失われる。fill の大きさが file（残りの byte）と pipe・cin で違うので、**同じ byte が届き方で失われたり残ったりする** | 決まった結果だが届き方に依る（regular file が基準。`from_bytes` は file と同じ） | 規則 2（file の挙動を基準に再現）。pipe との差は保守者に知らせる。Web は常に file の挙動 | LR-20、LR-21、AD-32 |
| 5 | `x_CleanAndCompress` で題の最後の byte が 0x80 以上なら落ちる（signed char） | 決まった結果（説明・見出し、4 program） | 規則 2（再現）。承認済みの例外 2 の punctuation と別 | RP-37 |
| 6 | `HtmlDecode` が `&#xD800;`・`&#0;`・`&#1114112;` などを検査なしで UTF-8 に書く（NUL や不正な列）。`GuessEncoding` が E0/ED/F0/F4 の 2 byte 目を検査しない | 決まった結果 | 規則 2（再現） | RP-33、RP-34 |
| 7 | `CheckDataLine` の 40% の曖昧な文字の警告は到達しない（`bIsNuc` が偽）。`AssignMolType` は推測せず、蛋白の配列を BLASTN に渡すと静かに取り除く | 決まった結果。到達しないコード | 規則 2（再現。移植は不要） | RD-35、RD-43 |
| 8 | 最初の行が `AB123456` のような accession の file は GenBank に接続する。取り寄せの失敗は query では黙って飛ばされ、subject では致命的。同じ配列の 2 つ目の綴りも黙って落ちる（原因は未確認：BI-23） | 機械と network に依る | 規則 3（拒否、§J-1） | BI-21、BI-23 |
| 9 | 同じ byte でも seekable な file なら `Query is Empty!`（0）、pipe なら `Empty CBlastQueryVector`（3）。空の pipe は成功で空の報告、`-query - -subject -` は黙って零 batch | 決まった結果 | 規則 2（再現） | LR-7、LR-8、BI-46 |
| 10 | Windows の `IsIStreamEmpty` は tellg を使わず、空の pipe を空とする | platform の差 | Linux の意味を全 platform で採る | LR-30 |
| 11 | 無効な残基の位置の一覧を 1 行あたり最大 1000 の範囲で黙って切る。`>?` の長さの上限が無く、非常に大きな gap は巨大な確保になる | 決まった結果／確保の失敗 | 切る方は規則 2。大きな gap は §J-4 の拒否 | RD-33、RD-36 |
| 12 | `-parse_deflines` と定義行の無い subject の誤りに NUL byte が入る | 範囲の外（§J-2） | 拒否のまま | BI notes §9-4 |
| 13 | NCBI の parameter の環境変数（`DIAG_`・`NCBI_CONFIG_` と別の名前）に NCBI が解析できない値：`DEBUG_CATCH_UNHANDLED_EXCEPTIONS`・`NCBI_STATIC_ARRAY_COPY_WARNING`・`NCBI_USAGE_REPORT_ENABLED`・`NCBI_USAGE_REPORT_MAXQUEUESIZE`・`OBJMGR_SCOPE_AUTORELEASE`・`OBJMGR_SCOPE_POSTPONE_DELETE`・`OBJMGR_BLOB_CACHE`（例 `maybe`・`abc`） | `Error reading CParam value …` か `terminate called after throwing …` の後に SIGABRT、`OBJMGR_BLOB_CACHE` は `CParamException` で終了コード 255（parameter の定義：`objmgr/scope_info.cpp:67`・`:89-90`、`objmgr/data_source.cpp:94-95`、`util/static_set.cpp:40-41`、`connect/ncbi_usage_report.cpp:75`・`:84`、`corelib/ncbiapp.cpp:396-397`）。どれも NCBI 自身の parameter の解析で、FASTA の読み方と関係しない。有効な値は出力を変えない（`NCBI_USAGE_REPORT_ENABLED`・`DO_NOT_TRACK` は usage report の `.ncbirc` の読み込みも変えない、§G1）。`NCBI_USAGE_REPORT_CONN_TIMEOUT`・`_CONN_MAX_TRY`・`_WAIT_TIMEOUT`・`DO_NOT_TRACK`・`THREAD_CATCH_UNHANDLED_EXCEPTIONS`・`NCBI_STATIC_ARRAY_UNSAFE_TYPE_WARNING` の解析できない値は network の無い実行では影響が無かった | LOSAT は異常終了を再現せず普通に走る（規則 1 に近い。保守者への問い：記録に留めるか、変数を §J の拒否に加えるか） | SFd 再監査の第 2 回 R-8 |
| 14 | 環境に `=` を含まない項目（`FOOBAR`、空の項目。`execve` で envp を作ったときだけ起きる） | NCBI 2.17.0 は SIGSEGV で出力なしに落ちる（blastn・blastp・tblastx、`-version` も。`CNcbiEnvironment::Reset` はその項目に `bad string` を記録して飛ばす、`ncbienv.cpp:99-103` が、どこで落ちるかは未確認） | LOSAT は異常終了を再現せず普通に走る（Rust の `std::env::vars_os` もその項目を飛ばす。保守者への問い：記録に留めるか拒否にするか。libc なしでは環境の生の項目を見られない） | SFd 再監査の第 3 回 S-4 |

## L. 判断と、棚卸しの食い違い・未解決

### L1. 作業の中で採った判断（Owner-delegated、2026-10-08、`DECISIONS.md`。保守者が 2026-10-06 に決めた 5 つ（DW-23）に加えて）

| 判断 | 内容 |
|---|---|
| 読み込み器の状態（Seq-id の拒否の後） | Seq-id の行 1 行を消費し、局所 ID の番号は進めない。次の呼び出しは次の行から読む。program は拒否で止まる（`read_all`、`read_queries`）。NCBI 自身の状態（BI の net2、net5） |
| 計画 | `port_plan.md` の S0〜S10。`AUTHORITY.md` は S1 の前 |
| ABI v1（Q1） | `bio` と今の検査を残し `from_bio` で入る（§J-5） |
| Seq-id の判定 | BI-15〜19、BI-55 の分類（§G4）。表なしで広めに拒否 |
| Web の `register`（Q2、Q3） | loader 有効で読み、読み込みの誤りと注釈だけの query を早く NCBI の文言で拒否。`run_local` は `QueryEnd::Input` |
| そのほか（Q4〜Q6、Q8、Q10、Q11） | byte の `HtmlDecode`・`GuessEncoding` を S2 で。定義行の無い record は `header_offset == sequence_offset`。`residue_counts` は大文字にした残基。`q_idx` は `ordinals` 経由。gap の行の修飾子は移植。`NCBI_CONFIG__BLAST__DATA_LOADERS` は registry の項目と同様に読む |
| BI §10 | query は遅延で batch ごとに読む。実効の `DATA_LOADERS` が読み手を決める。subject にも同じ規則。空の record・空の batch・零 batch の pipe は 4 program で移植。`-parse_deflines` は拒否のまま。空の `DATA_LOADERS=` は環境変数・`<prog>.ini` では項目があり loader 無効、`.ncbirc` では `BLAST_USAGE_REPORT` が偽か `<prog>.ini` があるときだけ項目（§G1。SFc B-1、SFd A-1・A-2・B-3） |
| registry の層（SFd の再監査の第 2 回、R-1〜R-9、2026-10-10） | 移植：探索 path を作る時点で決まる program の directory（usage report が `.ncbirc` を読まないときは起動したときの名前の directory が先。Linux で oracle と照合、macOS・Windows は source から）、`HOME` が無いときの passwd の home（Windows は `APPDATA`、無ければ `USERPROFILE`。source から）、byte のままの `NCBI_CONFIG_PATH`、`<prog>.ini` の 2 回の BOM の検査（UTF-16 は拒否のまま）、1 byte の file の `Error reading the registry after line 1: `（`.ncbirc` も。文言が固定で、回数は読む者で決まる）、開けない file は registry なし、空の `NCBI_CONFIG_OVERRIDES=` は無いのと同じ。拒否：passwd に利用者の項目が無く探索が home に届くとき（login 名の `getpwnam` を移さない）、開けたが読めない file（§J-8）。R-8 は §K-13（LOSAT は普通に走る） |
| registry の層（SFd の再監査の第 3 回、S-1〜S-5、2026-10-10） | 移植：重複した環境変数は NCBI の読み手ごとに最初（`getenv`）か最後（`CNcbiEnvironment`）、相対の `argv[0]` の末尾の `/` は `CFile` の検査で落とす。範囲の記録：`<prog>.ini` の名前は NCBI の program が自分の名前で起動された場合のもの（LOSAT の `argv[0]` は directory だけ）。拒否：空の `argv[0]`、`std::env::home_dir` が passwd の項目を返さないとき（項目なし、または buffer に収まらない長さ）に探索が home に届くとき（libc への依存は足さない）。S-4 は §K-14（LOSAT は普通に走る） |

### L2. 棚卸しの食い違いと未解決（この記録を書く際に見つけたもの）

1. **空の `DATA_LOADERS=`**：BI の notes（§10-3）は未確認で保守側（有効）としていた。SFb はソース（`ncbireg.cpp:984-991`、`:1300-1304`）から registry の file の空の値を項目なしと決め、`.ncbirc` について oracle で確かめた。SFd の再監査（A-1・A-2・B-3）で、この読みは `CCompoundRegistry::FindByContents` が `fCountCleared` を付けること（`:1235-1246`）を見落としていたと分かった：空の値は `<prog>.ini` でも `.ncbirc` でも項目で、`.ncbirc` の空の値が既定で項目なしになるのは usage report の cache の写し（`Write` が空の値を落とす）のためだけ（§G1）。環境変数の空の値（`NCBI_CONFIG__BLAST__DATA_LOADERS=`）は項目があり、両方の loader を切る（`env_reg.cpp:157-167`。SFc の監査 B-1 で訂正、§G1）。BI-54 の式「`loaders_on = 項目なし || blastdb/genbank を含み none を含まない`」は、空でない値の場合の正確な書き方。開けない `<prog>.ini` は無いのと同じ（`.ncbirc` は cache の写し）で、空や 1 byte の `<prog>.ini` は在る（`.ncbirc` は直接読まれ、空の値が項目になる。再監査の第 2 回 R-5・R-6）。
2. **AD-22 と判断 Q3**：AD-22 は「Web には data loader が無い（`UseDataLoaders` が偽）」とするが、`register` は loader 有効で読み Seq-id の行を拒否すると決めた。AD-22 の前提は古い。
3. **RD-40 と BI-42**：RD-40 は残基の無い query を「skipped」とするが、BI-42 は record が batch と報告に残り「ヒットなし」で数えられるとする。BI-42（`e01`・`e02`・`e11`・`e12` で観察）に従う。
4. **RP-20 の対象**：「outfmt 6/7/custom」と書くが、観察（`v06`）は outfmt 6 だけ。7 と custom はソースの読み（`SetFields` の経路）。一括の問いに出すとき明記する。
5. **loader 有効の実行**：Seq-id の観察は blastn（`net1`〜`net3`、`net5`、`net6`）と blastp（`net4`）だけ。tblastn・tblastx は同じコード。subject の Seq-id は blastn の `net6` 1 回。tblastn の query（蛋白）の `SDataLoaderConfig(true)` は読んで確認しただけ。この記録のために新しい network の実行はしていない。
6. **BI-23（未確認）**：既に取り寄せた accession の別の綴り（`gb|AB123456.1|`、`ab123456`）が黙って落ちる原因は未確認。LOSAT は全部拒否するので影響なし。
7. **AD-20（未確認）**：BLASTN の `q_idx` が batch をまたいで通しかは未確認。
8. **LR-1 の「推測」**：両方 `-` のとき `tellg()` が失敗するのは LR-1 が「推測」とし、BI-46 が eofbit による sentry の failbit で説明している。観察（零 batch、`Query is Empty!` なし）は 4 program とも一致。
9. **RP の重複行**：RP-2〜5、RP-9〜12、RP-21〜24、RP-49〜52 は program ごとに同じ内容の行（実際は約 38 の別の事柄）。件数 222 は重複を含む。
10. **`status_before` の古い箇所**：LR・BI の LOSAT 今の欄は新しい `fasta_reader/` の完成前（`new:` の注）。移植の後で再確認する（S0 のテストが基準）。
11. **`LOSAT/tests/unit/blastx_stage_e_io_stream_expected.tsv`**：行番号・`PeekChar`・`UngetLine`・`End()`・batch・標準入力の共有は覆わない（LR notes）。S0 の `from_bytes` と `from_file` の性質試験が補う。
12. **BLASTX**：この記録は BLASTX の関数を変えない（SX まで。計画 §4.2、DW-10）。BLASTX の読み方との差（題を `String` で持つ、警告が定義行で出る、`BATCH_SIZE` を無視、標準入力を拒否、`GenerateDefline` なし）は SX の引き継ぎ（RD、BI-32、RP の BLASTX の欄）。
