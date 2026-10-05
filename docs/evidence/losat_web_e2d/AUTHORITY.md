# E2d（S11）：`-query_loc` / `-subject_loc` — NCBI の経路

NCBI のソースは固定 commit `598d8ae6`（`/mnt/c/Users/genom/GitHub/ncbi-blast/c++`。作業では LF に揃えた複製 `~/.cache/losat-web-gui-target/s08p/ncbi/c++` を読んだ）、振る舞いの確認は NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin`）で行った。範囲の経路の関数の棚卸しは [`INVENTORY.tsv`](INVENTORY.tsv)（5 つの範囲、248 行。範囲ごとの指示・記録・oracle の観察は [`inventory/`](inventory/)。列 `status_before` は S11 の前の LOSAT `4fab73fdb`）。LOSAT の移植箇所には、それぞれの直上に NCBI のファイル・行と断片を書いてある。

範囲の記号：LA（引数）、QN（核酸の query：blastn・tblastx）、QP（蛋白の query：blastp・tblastn）、SB（subject）、RF（結果と報告）。

## 要約：NCBI の範囲は「区間の文字だけを検索し、報告はレコードに戻す」

NCBI は、範囲を、その役割（query か subject）の入力元が読む**すべてのレコード**に同じに当てる（1 つの範囲。レコードごとの範囲は無い）。各レコードの区間 `[from, min(to, len-1)]`（0 始まり・閉区間）の文字だけを検索し（統計・DUST・SEG・lowercase の mask・lookup・分割・翻訳の frame はすべて区間のもの）、結果の座標を区間の始まりだけ動かしてレコードの座標で報告する。報告の数のうち、`Length=`・`slen`・frame はレコード全体から、`qlen`・`qcovs` は**入力どおりの範囲**（レコードの端で切らない）から、database の合計は区間の長さの和から取る。LOSAT は、各レコードを区間に切って検索し（`seq_range.rs` の `cut_subjects`・`cut_queries`）、レコードの中の位置（`RecordPlacement`：区間の始まりとレコードの長さ）を報告に渡す。区間に切ったレコードを NCBI の範囲付きの検索と比べると、すべての program で検索と統計がバイト一致する（[`inventory/QN_notes.md`](inventory/QN_notes.md) §0、[`inventory/QP_notes.md`](inventory/QP_notes.md) の `evidence_slice_cmp`）。違いは下の §B〜§E の点だけで、どれも移植した。

## A. 引数（範囲 LA）

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_args.cpp:1945-1949`、`2372-2386` | `-query_loc`・`-subject_loc` は文字列の key。`-subject_loc` は `-db` の一族と `-remote` を除外（CArgs の誤りは USAGE、終了コード 1） | 4 つの program の `*Args` の `query_loc`・`subject_loc`。`cli.rs` の `is_unported_{blastn,blastp,tblastn,tblastx}_arg` から 2 つを除いた（BLASTX は SX まで拒否のまま）。CArgs の誤り（key の重複・値の欠落）は clap（承認済みの例外 1）。`-db`・`-remote` は LOSAT が拒否する |
| `blast_input_aux.cpp:145-179` `ParseSequenceRange` | `-` で `NStr::Split`（空の token を残す）。2 つの空でない token でなければ `(Format: start-stop)`。始め、終わりの順に `NStr::StringToInt`。0 以下は `(range elements cannot be less than or equal to 0)`、同じ値は `(range cannot be empty)`、始め > 終わりは `(start cannot be larger than stop)`。どれも `CBlastException` で `BLAST engine error: Invalid specification of query location (…)`（subject は `subject location`）、終了コード 3 | `seq_range.rs` の `parse_sequence_range`（`ncbi_string_to_int` は `i32::from_str` と同じ文字列を読む：先頭の `+`、先頭の 0 は可、空白・`0x`・`,`・小数点・2^31 以上は不可） |
| `ncbistr.cpp:635-643,798-880` | `StringToInt` が読めない部分は `CStringException`。`Error: NCBI C++ Exception:` と build のソースのパス・行（`/opt/conda/conda-bld/…/ncbistr.cpp`）、終了コード 255 | **明示的な拒否**（判断 R1、計画 TD-15 と同じ）：「the -query_loc value … has a part … that NCBI BLAST+ cannot convert to an int …; this is not supported by LOSAT's <PROGRAM>」、終了コード 1 |
| `blast_args.cpp:2537-2545` | subject の範囲は、subject のファイルを開いた後、読む前に読む（`CBlastDatabaseArgs` は 3 番目の handler）。`-subject` が無ければ範囲を読まずに `Either a BLAST database or subject sequence(s) must be specified` | BLASTN `run`、BLASTP `run`、TBLASTN・TBLASTX は `blastn/input.rs` の `read_nucleotide_subjects`（`subject_loc` を受ける）。`run_local` は最初に読む |
| `blast_args.cpp:1995-1999`、`blastn_args.cpp:44-120` ほか | query の範囲は query options の handler：filtering（`-dust`・`-seg`）・window size の後、formatting（`-max_target_seqs` の警告）と `Validate` の前 | BLASTN `search_cli`・`run_local` の `parse_query_range`（`resolve_dust` と `process_options` の間）、BLASTP・TBLASTN `check_options`（window size と `parse_formatting_string` の間）、TBLASTX `check_ncbi_options`（`seg_spec` と `formatting_handler_check` の間） |
| ABI v2 | — | `web/adapter/src/run.rs` の `validate`：subject の範囲、program の検査（query の範囲を含む）の順。レコードとの関係（範囲の始まりがレコードを越える）は `run`（`run_local`）が調べる |

誤りの順（oracle `inventory/scratch_LA` の `ord*`、固定 fixture `range_regression` の `parse_*`）：CArgs の誤り → `-outfmt` → subject のファイルを開く → **subject の範囲の文法** → subject のレコードを読む（**範囲の始まりがレコードを越える誤り**） → `-query`・`-out` を開く → filtering → **query の範囲の文法** → formatting → `Validate` → `Query is Empty!`。

## B. 入力のレコード（範囲 QN・QP・SB）

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_fasta_input.cpp:433-460` `x_FastaToSeqLoc` | レコード全体を読んでから区間を作る。0 始まりの始め `from` が長さを**越える**（`from > seqlen`、1 始まりの始めが長さ + 2 以上）と `CInputException eInvalidRange` `Invalid from coordinate (greater than sequence length)`。`from == seqlen`（1 始まりで長さ + 1）は通り、文字の無い区間 `[seqlen, seqlen-1]` になる。終わりが長さ以上なら黙って最後の文字に | `seq_range.rs` の `record_interval`（`Letters`・`Empty`・`PastEnd`） |
| `blast_input.cpp:198-219` `GetAllSeqs` | subject はすべて先に読み、例外を捕まえない：範囲を越えるレコードで `BLAST query/options error: Invalid from coordinate (greater than sequence length)` と `Please refer to the BLAST+ user manual.`、終了コード 1、出力の前。それまでに読んだレコードの題の警告は出る | `cut_subjects`（`app::options_error`）、題の警告は `subjects_read` までのレコード |
| `blast_input.cpp:144-155` `GetNextSeqBatch` | query は batch ごとに読み、`catch (const exception&) { continue; }`（SB-2307）：範囲を越える query のレコードは**黙って飛ばす**（警告も報告も無い。outfmt 7 の `# BLAST processed N queries` にも数えない） | `cut_queries`（`QueryInput` の `skipped`）。読み手の通し番号 `Query_<n>` は飛ばしたレコードも数える（`QueryInput::ordinal`、警告の `query_warning`・`invalid_query_warning`） |
| `blast_input.cpp:157-168` | batch の残基数は**レコード全体**の長さ（`sequence::GetLength(id)`）。飛ばしたレコードは数えない | `next_ranged_batch_end`・`QueryInput::batches`（BLASTN の `run_in_pool`、BLASTP の分割と skip の印、TBLASTN の batch、TBLASTX の `run_in_pool`）。範囲が無ければ `next_query_batch_end` と同じ |
| `objmgr_query_data.cpp:375-380` | 検索するレコードの無い batch は `BLAST engine error: Empty CBlastQueryVector`、終了コード 3。最初の batch なら outfmt 0 の prolog の後、後の batch なら前の batch の報告（epilog なし）の後。その batch の題の警告は前に出る | BLASTN `run_in_pool`、BLASTP `check_blastp_environment`（最初）と `run_resolved_in_pool` の最後、TBLASTN `search`、TBLASTX `run_in_pool`。後の batch のときは各 writer の epilog を書かない（BLASTP・TBLASTN の `write_*_pairwise_report(…, epilog)`、tabular の `# BLAST processed` も） |
| `blast_setup_cxx.cpp:632-639,653-657,733-850` | 文字の無い区間：query は `Sequence contains no data` の警告（batch の唯一の query なら終了コード 3 の誤り）、subject は `Subject sequence contains no data` の警告（BLASTN・TBLASTX で subject がそれだけなら `The average subject length is too short`、終了コード 3） | **明示的な拒否**（判断 R2）：query は 4 つの program、subject は BLASTN・TBLASTN・TBLASTX（`check_no_empty_interval`、文字の無いレコードの既存の拒否と同じ類）。BLASTP の subject は、LOSAT が文字の無い subject を NCBI と同じに扱う（警告・合計）ので、そのまま NCBI と一致する（fixture `e2d.blastp.sloc_empty.0`） |
| CFastaReader | 読み手の警告と誤り（無効な文字、ハイフン、題）はレコード全体に対して | LOSAT はそれらの入力を範囲の前に（レコード全体で）拒否する（以前から） |

## C. 検索（区間の文字だけ）

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_setup_cxx.cpp:153-281,485-658`、`blast_objmgr_tools.cpp:324` | context の長さ・バッファ・lookup・統計（長さの補正、有効な検索空間、Karlin のブロック、組成による統計）は区間の長さと文字から | 区間に切ったレコードを、そのまま既存の検索に渡す |
| `dust_filter.cpp:166-186` | DUST は区間の文字に（端の効果も区間の端） | 同上（検証：`inventory/scratch_QN/dust`、fixture `e2d.blastn.dust_edge.0`） |
| `blast_filter.c:1262` | SEG は区間の文字（tblastx は区間の翻訳）に | 同上 |
| `blast_setup.c:1029-1062` `BlastSeqLoc_RestrictToInterval`、`blast_setup_cxx.cpp:660-731` | lowercase の mask はレコード全体で読み、区間に切って移す | 区間に切ったレコードの lowercase をそのまま使う（同じになる） |
| `seqsrc_multiseq.cpp`、`blast_setup.c:896-980` | subject の長さ・合計・数は区間の長さ（文字の無い区間のレコードも数に入る）。tblastn・tblastx の合計は塩基数 | 区間に切った subject から（`Number of letters in database`、検索空間） |
| `split_query_cxx.cpp`、`split_query_aux_priv.cpp:104-145` | query の分割は区間の長さの和で決め、chunk の座標は区間の中 | 区間の長さ（BLASTN `split_query_batch`、BLASTP・TBLASTN `split_protein_batch`）。TBLASTX は以前から分けない（E2e R2c-1） |
| `blast_engine.c` | 長い subject の分割（`MAX_DBSEQ_LEN`）は区間の文字に | 同上 |

## D. 報告（範囲 RF）

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_seqalign.cpp:1520-1534` `RemapToQueryLoc`、`1541-1547` `s_RemapToSubjectLoc`、`seq_align_util.cpp:52-78` | Seq-align の query の行を区間の始まりだけ動かす。subject の区間は strand both なので mapper も始まりだけ動かす（minus の strand・frame も同じ） | BLASTN `shift_to_records`（`run.rs`）、BLASTP `shift_to_records`・`shift_hsp_to_records`、TBLASTN `record_coordinates`（`stage_e_report.rs`）、TBLASTX `shift_hit`・`shift_to_records`（`report.rs`）。表示の行は区間の文字から先に作る |
| `align_format_util.cpp:742-744`、`showalign.cpp:2479` | outfmt 0 の `Length=` は query・subject のレコードの長さ | `Placements::length` |
| `tabular.cpp:837-842` | `slen` は subject のレコードの長さ | BLASTP の `write_blastp_hsp_tabular_row`・PairwiseHit の `subject_length` |
| `tabular.cpp:884-885`（`*_app.cpp` の `formatter.SetQueryRange`） | `qlen` は入力どおりの範囲の長さ（`to - from + 1`。レコードの端で切らない。すべての query で同じ） | `QueryInput::range_length`（BLASTP の `qlen`。BLASTN・TBLASTN・TBLASTX は custom の field を書かない） |
| `showalign.cpp:411-440`、`align_format_util.cpp:1233-1245` | frame（outfmt 0 の `Frame =`、tabular の qframe・sframe）は、レコードの座標とレコードの長さから計算し直す | `seq_range.rs` の `record_frame`（TBLASTN の subject、TBLASTX の query と subject） |
| `blast_format.cpp:129-138` | database の合計は区間の長さの和 | 区間に切った subject から |
| `blast_aux.cpp:825-842` `Map`、`859-958` | 表示する query の mask：区間の座標の mask を区間の始まりだけ動かし、区間の終わりで切り、区間全体を覆う mask は出さない | BLASTN `pairwise.rs` の `shown_query_masks(…, offset)`、BLASTP `masked_query_rows(…, offsets)`、TBLASTX `report.rs` の `shown_query_masks` |
| `blast_seqalign.cpp:1588-1600`、`seqinfosrc_seqvec.cpp:110-160` | 表示する subject の lowercase：**レコードの座標**の mask を、HSP の**区間の座標**の `[offset, end]`（end を含む）と比べて残す（NCBI の振る舞いとして再現） | BLASTN `pairwise.rs` の `shown_subject_masks`（レコードの lowercase、HSP list ごと）。BLASTP・TBLASTX は `-lcase_masking` を拒否する（以前から）。tblastn の翻訳した subject は lowercase を出さない |
| `blast_aux.cpp:904-958`、`showalign.cpp:2500-2521` | tblastx の表示の mask：mask の frame の label は**区間**の frame（`BLAST_ContextToFrame`）、行の frame は**レコード**の frame。`locFrame == frame` で選ぶので、範囲の始まりが 3 の倍数でないと別の frame の mask を表示する（NCBI の振る舞いとして再現） | TBLASTX `report.rs` の `lowercase_query_row`（レコードの frame・座標）と `shown_query_masks`（区間の frame の label のまま） |

## E. NCBI の振る舞いとして再現するもの（決まった結果。`PD-LOSAT-NCBI-DEFECTS` の方針「決まった結果は再現」）

1. **分割した範囲付きの query の chunk の mask**（`split_query_cxx.cpp:256-281`）：chunk の Seq-interval はレコードの座標（`q_sl_offset`）だが、mask を切る区間は区間の座標（offset 0）で作る。mask はレコードの座標なので、chunk `[pf, pt)` の mask は `[pf + from, pt - 1]` の中だけ残る。BLASTN（`-lcase_masking`、DUST の mask も同じ経路）と TBLASTN（`-lcase_masking`）。LOSAT：`query_split.rs` の `restrict_masks(…, offset)`、`protein_chunk_part(…, offset)`。fixture：`range_regression` の `blastn.split_masks`（`CHUNK_SIZE=20000`）、`e2d.tblastn.split_mask.6`。
2. **subject の lowercase の表示**（§D）：fixture `e2d.blastn.sloc_lcase.0`（表示されない）と `e2d.blastn.sloc_lcase_start.0`（表示される）。
3. **tblastx の mask の表示の frame**（§D）：fixture `e2d.tblastx.seg5.0`・`seg31.0`・`seg61.0`。

## F. 判断（推奨の案で進め、記録した。保守者の常の指示）

- **R1**：`StringToInt` が読めない範囲（NCBI は build のパスを含む文言で終了コード 255）は明示的な拒否。計画 TD-15（`BATCH_SIZE` などの環境変数）と同じ扱い。
- **R2**：文字の無い区間（範囲の始まりがレコードの長さ + 1）は明示的な拒否（BLASTP の subject を除く）。NCBI の結果は決まっているが、実用が無く、LOSAT は文字の無いレコードも拒否している（同じ類）。保守者の DW-16 の基準（実用の無い設定で再現の費用が高いものは拒否のまま）に当たる。BLASTP の subject は、既存の空の subject の経路でそのまま NCBI と一致する。
- **R3**：subject のレコードの読み込みの順：NCBI はレコードを 1 つずつ読み、範囲の誤りで止まる。LOSAT はファイル全体を読んで LOSAT の拒否（NCBI が警告して読む入力）を先に出すことがある。NCBI と LOSAT がどちらも NCBI の誤りを出す組では順が同じ（題の警告も範囲の誤りのレコードまで）。違いは LOSAT の明示的な拒否の側だけ。
- **R4**：`-strand`（BLASTN・TBLASTX）は範囲と別の option で、拒否のまま（S11 の範囲外）。

## G. 範囲の外

- BLASTX（SX、計画 DW-11）。`is_unported_blastx_arg` は変えていない。
- 保存した検索の範囲（`blast_app_util.cpp:595-600`）：`-import_search_strategy` は LOSAT が拒否する。
- `-db` との組（`-subject_loc` は `-db` を除外する）：LOSAT は `-db` を拒否する。
- `qcovs`・`qcovhsp`・`qcovus`・`sstrand`（範囲に依る field）：LOSAT はどの program もこれらの field を書かない（BLASTP は明示的な拒否）。tblastx の HSP 1 つの subject の `qcovhsp` の癖（`inventory/RF_notes.md` §6.3）も、field を移植するときの記録。
