# LA notes: the arguments `-query_loc` / `-subject_loc` (blastn, blastp, tblastn, tblastx; bl2seq only)

Oracle: NCBI BLAST+ 2.17.0 (`/home/kawato/micromamba/bin/*`, the bioconda build). Inputs, scripts and every run are in
`scratch_LA/` (`run.sh` wrapper, `spell.sh`, `q3.sh`, `q3s.sh`, `mkinputs.py`; outputs `out/N_<label>.{cmd,out,err,rc}` for NCBI,
`out/L_<label>.*` for LOSAT-before). All NCBI line numbers are relative to `c++/` (`src/...`), pinned commit 598d8ae6.
LOSAT line numbers are relative to `LOSAT/src` of the frozen base copy.

## 1. Call path (NCBI), in execution order

1. `CArgDescriptions` parse (corelib, before `Run`): `-query_loc` is `AddOptionalKey(kArgQueryLocation "query_loc", "range", eString)`
   (`blast_args.cpp:1944-1948`, in `CQueryOptionsArgs::SetArgumentDescriptions` 1935-1965); `-subject_loc` is
   `AddOptionalKey(kArgSubjectLocation "subject_loc", "range", eString)` (`blast_args.cpp:2373-2386`, in
   `CBlastDatabaseArgs::SetArgumentDescriptions`, only when not RPS/KBlast/IgBlast). No constraint on the string. Dependencies:
   `-subject_loc` excludes every `database_args` entry (`-db`, `-gilist`, `-seqidlist`, `-negative_*`, `-taxids*`, `-no_taxid_expansion`,
   `-ipglist`, `-db_soft_mask`, `-db_hard_mask`; list built at 2211-2228, `SetDependency` at 2377-2383) and `-remote` (2384-2385). Errors here print the USAGE block,
   `Error: <message>` and `Error:  (CArgException::eXxx) <message>`, exit 1 (LOSAT: exception 1 of
   `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, clap message, exit 2).
2. `CBlastnApp::Run` etc. (`blastn_app.cpp:116-170`; the same in the other three apps): `RecoverSearchStrategy` (only
   `-import_search_strategy`), then `CBlastAppArgs::SetOptions(args)` (`blast_args.cpp:3585-3650`):
   a. `m_FormattingArgs->ArchiveFormatRequested(args)` (3625) calls `ParseFormattingString` (2795-2870) FIRST: a bad `-outfmt` is
      raised before any range is read (`Error: Formatting choice is out of range`, exit 255, or `'x' is not a valid output format`).
   b. `x_CreateOptionsHandle`, then the loop `ExtractAlgorithmOptions` over `m_Args` in the order the app constructor pushed them
      (`blastn_args.cpp:44-120`, `blastp_args.cpp`, `tblastn_args.cpp:44-135`, `tblastx_args.cpp:44-125`):
      - blastn: Program, Task, **BlastDb (subject read, `-subject_loc`)**, StdCmdLine (opens `-query` then `-out`), GenericSearch, Nucl, DMB,
        Filtering (`-dust`, `-soft_masking`, `-window_masker*`), Gapped, HspFiltering, WindowSize, OffDiagonal, MbIndex,
        **QueryOpts (`-strand`, `-query_loc`)**, Formatting, MT, Remote, Debug.
      - blastp: Program, Task, BlastDb, StdCmdLine, GenericSearch, Filtering (`-seg`), Matrix, WordThreshold, HspFiltering, WindowSize,
        **QueryOpts**, Formatting, MT, Gapped, Remote, CompBasedStats, Debug.
      - tblastn: Program, Task, BlastDb, StdCmdLine, GenericSearch, GeneticCode(db), Gapped, LargestIntron, Filtering, Matrix,
        WordThreshold, HspFiltering, WindowSize, **QueryOpts**, Formatting, MT, Remote, CompBasedStats, PsiBlast, Debug.
      - tblastx: Program, BlastDb, StdCmdLine, GenericSearch, LargestIntron, Filtering (`-seg`), Matrix, WordThreshold, HspFiltering,
        WindowSize, **QueryOpts**, GeneticCode(query), GeneticCode(db), Formatting, MT, Remote, Debug.
   c. `CBlastDatabaseArgs::ExtractAlgorithmOptions` (`blast_args.cpp:2425-2565`), `-db` absent: if `-subject` given (2525-2556):
      `args[kArgSubject].AsInputFile()` (open error first) -> `ParseSequenceRange(args["subject_loc"], "Invalid specification of subject
      location")` (2540-2545, only if `-subject_loc` was given) -> `ReadSequencesToBlast(stream, IsProtein(), subj_range, ...)`
      (`blast_input_aux.cpp:222-246`; `CBlastInputSourceConfig::SetRange`, `SetSubjectLocalIdMode`; `CBlastInput::GetAllSeqs`) ->
      `CObjMgr_QueryFactory(*subjects)` (empty file -> `Empty CBlastQueryVector`). If `-subject` is absent (and no `-db`):
      `CInputException eInvalidInput "Either a BLAST database or subject sequence(s) must be specified"` (2558-2562): `-subject_loc` is
      never parsed in that case.
   d. `CStdCmdLineArgs::ExtractAlgorithmOptions` (3454-3480): `-query` stream (`AsInputFile`), then `-out` (`AsOutputFile`). A missing
      query file and an unwritable `-out` are raised here (rc 1, `Command line argument error: Argument "query". File is not
      accessible:  `x'`).
   e. `CQueryOptionsArgs::ExtractAlgorithmOptions` (1967-2004): strand, then `m_Range = ParseSequenceRange(args["query_loc"],
      "Invalid specification of query location")` (1995-1999; only if given), then `-lcase_masking`, `-parse_deflines`.
   f. After the loop `m_OptsHandle->Validate()` (3640-3648): its errors (`-evalue 0`: `expect value or cutoff score must be greater than
      zero`, bad gap costs ...) come AFTER the range errors (oracle `sc3` vs `sc4`, `val4` vs `val5`).
3. `x_RunMTBySplitDB` (`blastn_app.cpp:179-305`; blastp 183-; tblastn 183-; tblastx 112-): `InitializeSubject`; `iconfig(dlconfig, strand,
   lcase, parse_deflines, query_opts->GetRange())` (blastn_app:205-209, blastp 207-210, tblastn 209-212, tblastx 128-131);
   `IsIStreamEmpty` -> `Query is Empty!` rc 0 (blastn 209, blastp 211, tblastn 213, tblastx 132); formatter ctor and
   `formatter.SetQueryRange(query_opts->GetRange())` (blastn 250, blastp 245, tblastn 262, tblastx 165; the TYPED range, see RF);
   `PrintProlog`; then per batch `input.GetNextSeqBatch(*scope)` (`blast_input.cpp:135-176`) -> `CBlastFastaInputSource::GetNextSequence` ->
   `x_FastaToSeqLoc` (`blast_fasta_input.cpp:377-466`).
4. `ParseSequenceRange` (`blast_input_aux.cpp:146-180`): `NStr::Split(str, "-", tokens)` with flags 0 (no merge, no truncate: empty
   tokens are kept); if `tokens.size() != 2 || front empty || back empty` -> `<prefix> (Format: start-stop)`; `from = StringToInt(front)`,
   `to = StringToInt(back)` (front first, BOTH before any other check); `from<=0||to<=0` -> `(range elements cannot be less than or equal
   to 0)`; `from==to` -> `(range cannot be empty)`; `from>to` -> `(start cannot be larger than stop)`; then `from--, to--` (0-based
   closed `[from,to]`). All three are `CBlastException eInvalidArgument` -> `CATCH_ALL` (`blast_app_util.hpp:221-227`): stderr
   `BLAST engine error: <msg>`, exit `BLAST_ENGINE_ERROR` = 3. A `StringToInt` failure is a `CStringException` -> the generic
   `CException` branch (`blast_app_util.hpp:240-243`): `Error: ` + `e.what()`, exit `BLAST_UNKNOWN_ERROR` = 255.
5. `x_FastaToSeqLoc` range checks (`blast_fasta_input.cpp:433-460`): no `-query_loc`/`-subject_loc` -> `from=to=0` -> whole record
   `[0, seqlen-1]`. Else `from`,`to` as parsed (0-based); `to>0 && to<from` -> `Invalid sequence range` (unreachable: the parser guarantees
   from<to); `from > seqlen` -> `CInputException eInvalidRange "Invalid from coordinate (greater than sequence length)"` (450-453);
   `SetFrom(from)`; `SetTo((to > 0 && to < seqlen) ? to : seqlen-1)` (459-460): an end at or past the record end is silently set to the last
   letter. Hence with the 1-based start `s`: `s > seqlen+1` error; `s == seqlen+1` (from == seqlen) gives the EMPTY interval `[seqlen, seqlen-1]`;
   `s == seqlen` gives a one-letter interval.

## 2. Coordinates

- User string: 1-based, closed, nucleotide letters for blastn/tblastx queries, amino acids for blastp/tblastn queries; subject: nucleotide
  letters (blastn, tblastn, tblastx) or amino acids (blastp). One range for every record of the role (all query records share `-query_loc`; all
  subject records share `-subject_loc`).
- After parse: 0-based closed `[from, to]` (`from--, to--`); after `x_FastaToSeqLoc`: `[from, min(to, seqlen-1)]` on strand both (nucleotide) or
  unknown (protein) (unless `-strand` plus/minus, which sets the strand of the same interval).
- Output coordinates stay in record coordinates (oracle `stdin_s`: `-subject_loc 400-1200` hit at `sstart` 501; `k3`: `qstart 1`). `qlen`, `Length=`,
  `# Query:` etc. belong to RF.

## 3. The spellings (oracle: `spell.sh`, `spellings.tsv`; `out/N_sp<i>_<prog>_{q,s}.{out,err,rc}`; 35 spellings x 4 programs x 2 roles, outfmt 6)

For all four programs and both roles the result is identical (the only difference is the program-specific result: `-query_loc 1-2` on tblastx
gives `Warning: [tblastx] Query_1 nq1 len1300: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its
translation. Please verify the query sequence(s) and/or filtering options ` rc 0, no hits). Prefix `P` = `Invalid specification of query
location` (query role) or `Invalid specification of subject location` (subject role). stdout is empty on every error; every error is raised
before any output (also with `-outfmt 0`/`7`: stdout 0 bytes; `o0_sp*`).

| spelling | result |
|---|---|
| `10-20`, `1-2`, `+1-5`, `01-05`, `1-+5`, `1-2147483647` | accepted (0-based `[9,19]`, `[0,1]`, `[0,4]`, `[0,4]`, `[0,4]`, `[0,2147483646]`) |
| `1-1`, `5-5` | rc 3 `BLAST engine error: P (range cannot be empty)` |
| `20-10`, `2-1` | rc 3 `BLAST engine error: P (start cannot be larger than stop)` |
| `0-10`, `10-0` | rc 3 `BLAST engine error: P (range elements cannot be less than or equal to 0)` |
| `-10`, `10-`, `10`, `1--5`, `1-2-3`, ``, `-`, `--`, `-5-10` | rc 3 `BLAST engine error: P (Format: start-stop)` |
| ` 1-5`, `1- 5`, `1 -5`, `1-5 `, `10 -20`, `a-5`, `1-5a`, `1.0-5`, `0x10-0x20`, `1,0-50`, `5-5x`, `0-a`, `1-0x`, non-ASCII digit | rc 255 `Error: NCBI C++ Exception:` (see below) |
| `1-2147483648`, `1-4294967296`, `2147483648-1`, `0-2147483648` | rc 255, `... StringToInt() - Cannot convert string '2147483648' to int, overflow (m_Pos = 0)` (ncbistr.cpp line 640) |
| `1-9999999999999999999`, `1-9223372036854775808` | rc 255, `... StringToInt8() - Cannot convert string '9999999999999999999' to Int8, overflow (m_Pos = 18)` (line 852) |

The rc-255 stderr is exactly (bytes, `cat -A`): `Error: NCBI C++ Exception:\n    T0 "/opt/conda/conda-bld/blast_1754905306566/work/c++/src/corelib/ncbistr.cpp", line 861: Error:
(CStringException::eConvert) ncbi::NStr::StringToInt8() - Cannot convert string ' 1' to Int8 (m_Pos = 0)\n\n`. Line 861 = the first character is not a digit
(token ` 1`, `a`, `0x...` is NOT this: see 870), 870 = trailing garbage (`5 `, `5a`, `1.0`, `0x10`, `1,0`; `m_Pos` is the position of the first bad byte), 852 = Int8
overflow, 640 = int overflow. Non-printable/non-ASCII bytes are shown as octal escapes (`'\331\243'`). The text contains the build path of
the bioconda binary: not reproducible from NCBI's source, same situation as TD-15 (`docs/losat_web_gui_plan.md:96`, the `BATCH_SIZE` integer
parse, where LOSAT rejects explicitly and exits 1). **Decision D1** for the session (see section 9).

`NStr::StringToInt` with flags 0 (`ncbistr.cpp:635-643,798-880`): optional ONE leading `+` (`-` cannot occur: it is the split delimiter), then decimal digits
only; leading zeros fine; no spaces (leading or trailing), no `0x`, no `,` (fAllowCommas is off), no `.`; at least one digit; value must fit in `int` (`<= 2147483647`).
`i32::from_str` on the token accepts and rejects exactly the same strings (existing helper `blastinput/query_batch.rs:98 ncbi_string_to_int`).
The two ints are converted front first, then back, BEFORE the `<=0` / `==` / `>` checks (`0-a` reports `a`, `a-0` reports `a`). `Split` with flags 0 keeps empty tokens:
`1--5` is 3 tokens, `-10` is `["", "10"]`, the empty string yields 0 tokens (all `(Format: start-stop)`).

CArgs level (all four programs; parse time, USAGE block + `Error:` lines, rc 1; covered by approved exception 1 -> LOSAT clap text, rc 2):
- key twice (`-query_loc 1-500 -query_loc 2-600`, same for `-subject_loc`): `Error: Argument with this name is defined already: query_loc` (`k1`,`k2`).
- value missing (`-query_loc` last, or followed by end): `Error: Argument "-query_loc". Value is missing` (`k4`,`k5`).
- `--query_loc`: `Unknown argument: "-query_loc"` (`k6`); `-query_l`, `-Query_loc`: unknown argument (`abbr`,`upper`); `-query_loc=1-500` accepted (`k3`).
- the value may start with `-` and may be empty: `-query_loc -5-10`, `-query_loc=-5-10`, `-query_loc ""`, `-query_loc=` all reach `ParseSequenceRange`
  (rc 3 `(Format: start-stop)`; `k16-k19`). LOSAT's `cli.rs try_parse_from` (cli.rs:121-160) already takes the next token as the value whatever it
  looks like and passes `--name=value`, and clap accepts an empty value (`clapempty`, `clapdash` show the same for `-dust ""`/`-seg -10`), so the new `Option<String>`
  field needs only `allow_hyphen_values`-free handling (the translation already neutralises it) and no value parser.
- `-query_loc -- 1-500` (the `--` is taken as the value of nothing, `1-500` is positional): `Too many positional arguments (1), the offending value: 1-500` (`k20`).

## 4. Order of errors (oracle `ord*`, blastn unless noted)

1. CArgs parse errors (bad `-evalue`: `ord3`/`ord20`; `-num_threads 0`; `-max_target_seqs 0`; `-word_size 3`; `-perc_identity 200`; `-task bad`; `-soft_masking foo`; `-strand foo`;
   `-query_gencode 99`; `-db` with `-subject_loc`; duplicate keys): USAGE, rc 1, always before any range error (also before the subject range error `ord20`).
2. `-outfmt` errors (`ord4`: `Error: Formatting choice is out of range`, rc 255; before both range errors, `ord15` too).
3. Subject file open error (`ord21`: `Command line argument error: Argument "subject". File is not accessible:  `nonexist.fa'`, rc 1) precedes the SUBJECT RANGE parse.
4. SUBJECT range parse and subject read (all records; per record: FASTA read, then range check): `ord1`, `ord16` (before `-dust foo`), `ord17` (before `-out` error),
   `ord18` (before a missing `-query` file), `ord2` (subject error beats a bad `-query_loc`), `ord40` (but `-query_gencode 99` is a CArgs error and comes first).
5. `-query` open (`ord19`: `Command line argument error: Argument "query". File is not accessible`) and `-out` open (`ord5`, `ord22`: before `-dust foo`).
6. GenericSearch ... Filtering (`-dust foo`: `BLAST query/options error: Invalid number of arguments to filtering option`, rc 1 `ord6`; `-seg foo` same text in blastp/tblastn/tblastx
   `ord36-38`), WindowSize ...
7. QUERY range parse (`ord1` ... `ord44`). tblastx: before `-query_gencode`/`-db_gencode` handlers (CArgs constraints in any case: `ord33`).
8. Formatting (custom fields; `ord23` shows an unknown field is NOT checked before the query range error), MT, Validate (`-evalue 0` `sc3`/`sc4`; `val4`/`val5`).
9. `x_RunMTBySplitDB`: `InitializeSubject`; the empty-query check (`Query is Empty!` rc 0 even with a VALID `-query_loc`, `emp_empty_loc`; a BAD `-query_loc` on an empty file still fails: `emp_empty_badloc` rc 3);
   an empty subject FILE: `Empty CBlastQueryVector` rc 3 whatever `-subject_loc` (valid or not: a bad one fails first, `emp_subj_bad_empty`).
10. per query batch: record-dependent range effects (section 5).

`-subject_loc` without `-subject` (and without `-db`): `BLAST query/options error: Invalid ... ` is NOT raised; stderr `BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified\nPlease
refer to the BLAST+ user manual.`, rc 1, also when `-subject_loc` is malformed (`ord11`, `ord12`). With `-db`: USAGE, `Error: Argument "subject_loc". Incompatible with argument:  `db'` (`ord13`); `-subject` with `-db`:
`Argument "subject". Incompatible with argument: `db'` (`ord14`); with `-remote`: `Argument "subject_loc". Incompatible with argument:  `remote'` (`dep1`). tblastn `-in_pssm` with `-query_loc`:
`Argument "query". Incompatible with argument:  `in_pssm'` (`dep2`; LOSAT rejects `-in_pssm`).
`-num_threads N` with `-subject`: stdout identical to 1 thread, NCBI adds `Warning: [<prog>] 'num_threads' is currently ignored when 'subject' is specified.` (`nt4_*`; approved exception 2). `-mt_mode 1` with `-query_loc` is the same (`dep9`).
No other CArgs-level interaction exists between the range options and `-lcase_masking`, `-dust`, `-seg`, `-soft_masking`, `-query_gencode`, `-db_gencode`, `-comp_based_stats`, `-max_target_seqs`, `-evalue`,
`-outfmt`, `-task` (`cmb1-cmb11` all rc 0).

## 5. Start past (or at) the end of a record (oracle `q3m_*` queries, `q3s_*` subjects, `bp_*`/`q3_bp_*` batches, `bd_*`, `idx*`; `q3_matrix.txt`, `q3s_matrix.txt`)

Notation: record length `L`, 1-based start `s` (the end is chosen so that it is past the record end).

QUERIES (read in batches by `GetNextSeqBatch`; blastn, blastp, tblastn, tblastx identical):
- `s > L+1` (from > seqlen): `x_FastaToSeqLoc` throws `CInputException eInvalidRange`, but `GetNextSeqBatch` has `catch (const exception&) { continue; }`
  (`blast_input.cpp:153-155`, "SB-2307") which swallows it: the record is SILENTLY SKIPPED, no message, rc 0, whatever its position (`q3m_*_past_{ABA,BAA,AAB}`: 2 of 3 queries searched,
  outfmt 6 and 0). If every record of a batch is skipped: `CObjMgr_QueryFactory` throws `Empty CBlastQueryVector` (`objmgr_query_data.cpp:376-380`), stderr `BLAST engine error: Empty CBlastQueryVector`, rc 3
  (single bad record `*_past_B`; outfmt 0 has already printed the 12-line prolog (stdout 12 lines for blastn), outfmt 6/7 nothing printed for that batch; with outfmt 7 the earlier batches' header+rows are printed,
  `q3_bp_longpast_7`: 6 lines then rc 3). The skipped record still consumes its ordinal: the reader's local id `Query_<n>` counts it (`idx1`: file [B(900), A, C901] with `-query_loc 902-951`: B skipped, C901 is
  `Query_3`).
- `s == L+1` (from == seqlen): the interval is `[L, L-1]`, empty. `SetupQueries_OMF` (`api/blast_setup_cxx.cpp:632-639`, via `IBlastSeqVector::size()` `blast_setup.hpp:192-197`) catches the `CBlastException`, records a warning and
  invalidates the query's contexts: stderr `Warning: [<prog>] Query_<n> <title>: Sequence contains no data ` (trailing space), rc 0, no hits for it; outfmt 0 still prints its `Query=` block with `***** No hits found *****` and `Effective search space used: 0`
  (`q3_bp_long_0`, lines 79-90). If no query of the batch is valid: `BlastSetup_Validate` fails (`blast_setup_cxx.cpp:653-657`) -> `CBlastException eSetup`, stderr `BLAST engine error: Warning: Sequence contains no data`, rc 3 (outfmt 0: prolog printed first;
  `*_eqp1_B`). Whether the bad query shares a batch with good ones decides warning (rc 0) or error (rc 3): see the batch accounting below.
- `s == L` (from == seqlen-1): one-letter interval: searched, no hits, rc 0 (`*_eq_*`); tblastx adds the Karlin-Altschul warning of section 3 for a 1-nt query.
- the end past the end (`1-<L+1000>`, `1-2147483647`): silently truncated to the record end, no message (`*_endpast_*`).
- BATCH ACCOUNTING USES THE FULL RECORD LENGTH. `CBlastInput::GetNextSeqBatch` adds `sequence::GetLength(loc->GetInt().GetId(), scope)` (`blast_input.cpp:157-168`), the whole bioseq, not the interval.
  Oracle: blastp batch size 10000: file [P1 (10050 aa), P2 (400 aa)] with `-query_loc 401-800`: P1 (range of 400) alone fills batch 1 and is printed, batch 2 = [P2] -> `BLAST engine error: Warning: Sequence contains no data` rc 3 after P1's rows
  (`q3_bp_longP1P2_6`: 1 row on stdout, rc 3; outfmt 0 stdout 76 lines); the SAME file with a 1000-aa first record (`bp_S1P2`) puts both in one batch: `Warning: [blastp] Query_2 P2 len400: Sequence contains no data` rc 0
  (`q3_bp_shortS1P2_6`). With `-query_loc 402-800` P2 is swallowed: long file -> batch 2 empty -> `Empty CBlastQueryVector` rc 3 after P1's output; short file -> rc 0 (`q3_bp_longpast_*`, `q3_bp_shortpast_6`).
  So with a range, earlier queries' results ARE printed before a later batch's error (outfmt 0, 6 and 7), and batch boundaries depend on the full lengths.

SUBJECTS (all read up front by `GetAllSeqs`, `blast_input.cpp:199-220`, no catch):
- `s > L+1` for ANY record: the exception propagates out of `SetOptions`: stderr `BLAST query/options error: Invalid from coordinate (greater than sequence length)\nPlease refer to the BLAST+ user manual.`, rc 1, stdout 0 bytes (also outfmt 0/7: nothing,
  not even the prolog), position of the record in the file irrelevant (`q3s_*_past_{SLL,LSL,LLS,S}`), the first failing record in file order. FASTA read errors of later records are not reached; those of earlier records come first.
- `s == L+1`: empty interval: per-subject warning `Warning: [<prog>] Subject_<n> <title>: Subject sequence contains no data` (ordinal = file position; same text as for an empty subject record, `emp_subj_emptyrec`), the subject is
  dropped, rc 0. If the only subject is dropped: blastn and tblastx stderr `BLAST engine error: The average subject length is too short`, rc 3 (`core/blast_message.c:216`; outfmt 0: 12 prolog lines printed), blastp and tblastn rc 0 with an
  empty result (`q3s_*_eqp1_S`).
- `s == L`: one-letter subject, searched, rc 0; `s` anywhere with the end past the end: silently truncated.

EMPTY RECORDS: a query record without letters behaves as without a range (`from = 0` always passes `from > 0`): `Warning: [blastn] Query_2 emptyrec: Sequence contains no data`, rc 0, or rc 3 `BLAST engine error: Warning: Sequence contains no data` when it is the only
record (`emp_qe_*`). LOSAT rejects records without residues explicitly today (adapter doc, `input.rs:475`), unchanged by a range.

STDIN: `-query -` and `-subject -` take the range the same way (`stdin_q`, `stdin_s`); an empty stdin query with a valid `-query_loc`: `Query is Empty!` rc 0 (`stdin_q_empty`).

## 6. LOSAT side (frozen copy of base commit 4fab73fdb)

- Parsing: `cli.rs try_parse_from` (cli.rs:63-160) translates `-name value`/`-name=value` to clap `--name=value`; a name that is not an argument of the subcommand goes to `unknown_option_error` (cli.rs:334-400).
  `query_loc` and `subject_loc` are not fields of `BlastnArgs`/`BlastpArgs`/`TblastnArgs`/`TblastxArgs`, and are named in `is_unported_blastn_arg` (cli.rs:247, 256), `is_unported_blastp_arg` (313, 322),
  `is_unported_tblastn_arg` (522, 530), `is_unported_tblastx_arg` (580, 590); the message is `the NCBI BLAST+ option -query_loc is not supported by LOSAT's BLASTN` (BLASTP, TBLASTN, TBLASTX alike), clap exit 2, stdout empty
  (`out/L_lb_*`: identical for every spelling, even `5-5`, `a-5`, `-query_loc` without value, key twice; `-strand` is rejected the same way in blastn and tblastx, and is `unknown option or argument '-strand'` in blastp/tblastn).
  A rejected option is raised at translation time, before any other check. `algorithm/tblastn/args.rs:1174` is only the unit test (`["-subject_loc", "1-30"]` expects the "is not supported by LOSAT's TBLASTN" text); there is no
  separate rejection in `args.rs` (the line 364 of the brief is the `window_size` computation in this copy; the rejection is only in `cli.rs`).
- ABI v2 adapter (`web/adapter/src/run.rs:62-97`): `parse` calls the SAME `LOSAT::cli::try_parse_from`, so the rejection is shared; `validate` = `parse` + per-program checks (`blastn::scoring::check_scoring` after `resolve_dust`,
  `tblastx::check_options`, `blastp::blast_engine::check_options`, `tblastn::check_options`). Record-dependent range errors cannot be found in `validate` (no records) nor in `register` (the argv is not available; each input is checked alone);
  they belong to `run` (`run_local` with the registered records, `run_web_pair*`). The syntax errors (rc-3 messages) can be raised in `validate` in NCBI's order: subject range first, then query range.
- Run order today (CLI): blastn `run` (`algorithm/blastn/blast_engine/run.rs:5134`): `parse_output_formats` (outfmt) -> `validate_threads` -> missing `-subject` (5209) -> `open_input(subject)` (5214) -> `read_fasta_bytes` -> `read_records` (5217) ->
  `check_subjects_not_empty` (5224) -> `open_input(query)` -> create `-out` -> `search_cli` (5315): `resolve_dust` (5323) -> `process_options` (5324: few-matches warning, `check_scoring_options`) -> read query -> blank query warning -> `parse_fasta` (5365).
  blastp `run` (`algorithm/blastp/blast_engine.rs:4581`): `parse_formatting_string` -> `report_format` -> `validate_threads` -> missing subject (4627) -> `open_input(subject)` (4630) -> `read_fasta_bytes` -> `bio_records_of` (4635) -> query open -> `-out` -> `search_cli` (4729):
  `args.check_options` (4738; `blastp/args.rs:467`: seg parse, ..., `window_size` 599, `parse_formatting_string` 611, `validate_protein_options`) -> read query -> `bio_records_of` (4791).
  tblastn `run` (`algorithm/tblastn/args.rs:583`) and tblastx `run` (`algorithm/tblastx/blast_engine/run_impl.rs:621`): same shape, the subject read is `blastn/input.rs read_nucleotide_subjects` (input.rs:808: `open_input` 824, `read_fasta_bytes` 825);
  `check_options` at tblastn args.rs:704 (`TblastnArgs::check_options` 205: seg 331, `window_size` 370, `parse_formatting_string` 381, comp_based_stats, `validate_protein_options`); tblastx `check_ncbi_options` (run_impl.rs:810: `seg_spec` 814, `formatting_handler_check` 815, Validate
  messages 825-837). The structure already follows NCBI's order, so the two new steps slot in as follows.
- Where each range step must go (NCBI order, section 4): SUBJECT range parse = immediately after the subject file is opened and before it is read (blastn run.rs:5214-5215; blastp blast_engine.rs:4630-4631; tblastn/tblastx inside `read_nucleotide_subjects`
  between input.rs:824 and 825, as a parameter), i.e. after outfmt parse, `-subject` check and the subject open error, before the query open, `-out` open, option extraction. QUERY range parse = after the Filtering/WindowSize extraction and before the formatting handler and the
  Validate checks: blastn `search_cli` between `resolve_dust()` (5323) and `process_options` (5324); blastp `blastp/args.rs check_options` between `options.window_size` (599) and `parse_formatting_string` (611); tblastn `tblastn/args.rs check_options` between
  `window_size` (370) and `parse_formatting_string` (381); tblastx `check_ncbi_options` between `args.seg_spec()?` (814) and `formatting_handler_check` (815). For `validate` of the adapter the same two parses run (subject first) from `validate` itself.
- Record-dependent steps: query records are parsed by `parse_fasta` (blastn/tblastx) / `bio_records_of` + checks (blastp/tblastn); subjects by `read_records` / `bio_records_of` / `read_nucleotide_subjects`. The range cut, the `from > seqlen` skip (queries) or error
  (subjects), the empty-interval warnings and the clamp of the end belong right after a record is read (NCBI does them per record in file order: a FASTA error in record k+1 comes after a range error in record k; LOSAT parses the whole file first, a corner case to note).
- Batches: `blastinput/query_batch.rs:19-62` (`query_batches`, `next_query_batch_end`) and `report/query_warnings.rs:73` take `lengths`; they are fed `record.seq().len()` at blastn run.rs:5882, 6421, 12695, blastp blast_engine.rs:5274 (and 5116, 5437), tblastn args.rs:910, tblastx run_impl.rs:1232, 1504. With a cut record these would be
  range lengths; NCBI uses the full record length (section 5). The `Query_<n>` numbering (`report/query_warnings.rs:187`, `index + 1`) must keep counting skipped records.
- Errors to add: `engine_error` (`blastinput/app.rs:40`, rc 3) and `options_error` (app.rs:21, rc 1) already produce NCBI's texts and exit codes for `(Format: start-stop)` etc. and for `Invalid from coordinate ...`; there is no builder for rc 255 / `Error: NCBI C++ Exception`
  (only `Error: Formatting choice is out of range`, app.rs:243). `Empty CBlastQueryVector` exists (`app.rs:74 empty_subjects_error`); `Warning: [prog] Subject_n title: Subject sequence contains no data` exists (`blastn/input.rs:360-373`); the query-side `Sequence contains no data`
  does not exist (LOSAT rejects records without residues).

## 7. Surprises

1. A query record whose range start is past the end is SILENTLY DROPPED (not an error), because of `catch (const exception&) { continue; }` in `GetNextSeqBatch`; a subject record in the same situation is a hard error rc 1. Do not "fix" this.
2. The batch size accounting uses the FULL record length even with a range (so `Sequence contains no data` is a warning rc 0 or an error rc 3 depending on the batch the empty-interval query lands in; results of earlier batches are printed first).
3. `from == seqlen` (start = L+1) is accepted by the check and produces an EMPTY interval (warning/error "no data"); start = L is a one-letter query. The end past the record end is clamped silently.
4. Subject range is parsed and applied before ANY other extraction, even before `-out`/`-query` are opened; the query range after the Filtering handlers but before Formatting/Validate. `-outfmt` is parsed before both.
5. Number-parse failures (`a-5`, ` 1-5`, `1-5 `, `0x10-0x20`, overflow >= 2^31) are NOT the `BLAST engine error` frame: rc 255, `Error: NCBI C++ Exception:` with the oracle build's source path and line (build specific).
6. `+1-5`, `01-05`, `1-+5` are accepted; `5-5` and `1-1` are errors; `1-2` is valid (a 2-letter query: tblastx warns, no hits).
7. The subject range applies to every subject record, the query range to every query record; `m_QueryRange` handed to the formatter is the typed range (RF).
8. The ids of skipped/empty-interval records keep their file ordinal (`Query_3`, `Subject_2`).

## 8. Out of scope here (and why)

- Range from an imported search strategy (`blast_app_util.cpp:595-600`, `SetRange(strategy.GetQueryRange())`): needs `-import_search_strategy` (ASN.1 Blast4 request), which LOSAT rejects as unported (is_unported_*_arg); the saved-strategy path also never reads
  `-subject_loc` (the strategy carries the database). `x_IssueWarningsForIgnoredOptions` (`blast_args.cpp:3697-3710`, lists `query_loc`/`subject_loc` as overridable) is the same branch.
- `-db`, `-remote` combinations (`-subject_loc` excludes them): LOSAT rejects `-db` and `-remote`; the exclusion can never fire in LOSAT.
- IgBlast/RPS/Mapper code in `blast_args.cpp` that reads subjects with a range (1887-1899) and `CMapperQueryOptionsArgs`: no LOSAT program.
- `-strand`: out of scope (rejected today); recorded for the session: `-query_loc 1-500 -strand plus` gives the same hit as without `-strand` (`k7`), `-strand minus` finds nothing for the plus-strand hit (`k8`), `-strand both` = default (`k9`); tblastx
  accepts `-strand` in NCBI (`k12`, `k13`: plus 23 rows, minus 25, none 48 on the 1-500 range), blastp and tblastn do not (`Unknown argument: "strand"` USAGE rc 1, `ord42`, `ord43`).

## 9. Decisions the session must take

- D1: how to treat the rc-255 `CStringException` text (number-parse failures). NCBI bytes include the bioconda build path; precedent TD-15 rejects explicitly (exit 1, not byte-identical). Options: (a) reproduce the oracle build's string verbatim (path
  `/opt/conda/conda-bld/blast_1754905306566/work/c++/src/corelib/ncbistr.cpp`, lines 861/870/852/640, `m_Pos`, `PrintableString` escaping) with exit 255; (b) reject with `... is not supported by LOSAT's <PROGRAM>`, exit 1 (TD-15). The table proposes (b) as `reject`
  and (a) as the alternative.
- D2: the `-query_loc` on tblastx with `-strand` and `-lcase_masking`/`-soft_masking` stay rejected (not part of this range).
- D3: reading order of subject FASTA vs range errors (per record interleaving) needs the subject reader to stop at the first failing record; today LOSAT reads and checks the whole file and defers `subject_checks`.
