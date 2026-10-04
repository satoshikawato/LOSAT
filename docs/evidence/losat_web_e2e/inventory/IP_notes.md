# IP notes: protein FASTA input (BLASTP query and subject, TBLASTN query)

Table: `IP.tsv` (55 rows: 27 divergent, 14 unported, 11 faithful, 3 n/a, 0 rejected, 0 exception).
Evidence: `scratch_IP/`. A tag `X.Y` means the four files `X.Y.ncbi.{out,err,rc}` and `X.Y.losat.{out,err,rc}`
(`X` = the input file `X.faa`; `Y` = run kind: `P6` blastp, input is the query, subject `s1.faa`, outfmt 6;
`S6` blastp, input is the subject, query `q1.faa`; `N6` tblastn, input is the query, subject `nt1.fna` (back-translation of `s1.faa`);
`Q7`/`Q0` and `T7`/`T0` the same with outfmt 7/0 for the query role (Q) and the subject role (T); `NQ7`/`NQ0` tblastn query role outfmt 7/0;
`F0`/`F7` blastp query role outfmt 0/7 for the error cases). `scratch_IP/cmp.sh`, `runall.sh`, `runsubj.sh`, `io.sh`, `show.sh` run both binaries;
`mk.py`/`mk2.py` and the small generators in the shell history created the inputs (`cases*.txt` list them). Oracle: `/home/kawato/micromamba/bin/{blastp,tblastn}` 2.17.0, LOSAT: `bin/LOSAT-before`.

## Call path followed (NCBI, pinned commit 598d8ae6)

`blastp_app.cpp` `CBlastpApp::Run` (195-262) and `tblastn_app.cpp` (196-262), in this order:

1. Argument processing: `CBlastDatabaseArgs::ExtractAlgorithmOptions` reads the **subject** first (`blast_args.cpp:2540-2560` -> `ReadSequencesToBlast`, `blast_input_aux.cpp:222-244`:
   `CBlastInputSourceConfig` with `SetSubjectLocalIdMode()` = prefix `Subject_`, `-parse_deflines`, `-lcase_masking`, then `CBlastInput::GetAllSeqs` = `CBlastFastaInputSource::GetNextSequence` in a loop).
   Reader warnings of the subject appear here (stderr, plain text). A subject without any record: `CObjMgr_QueryFactory` throws `Empty CBlastQueryVector` (`objmgr_query_data.cpp:375-380`) -> `BLAST engine error:` rc=3.
2. `InitializeSubject` -> `SetupSubjects_OMF` (`blast_setup_cxx.cpp:733-860`): empty subject records are skipped with `ERR_POST(Warning)`; the O->X replacement of `GetSequenceProtein` also happens here (no message for subjects).
   (Observed: this warning appears only if the query is not empty, i.e. after the `Query is Empty!` check, and before any query FASTA warning.)
3. The `-query` stream is opened (`File is not accessible`, `ncbiargs.cpp:606-620`, rc=1), `IsIStreamEmpty` (`blast_app_util.cpp:846-870`): not a single non-white-space byte -> `Warning: [blastp] Query is Empty!`, exit 0. A pipe (`tellg() < 0`) is never empty.
4. `CBlastInputSourceConfig(dlconfig, strand, lcase, believe_defline=-parse_deflines, range)` (`blast_input.cpp:54-76`: prefix `Query_`, `seqlen_thresh2guess = UINT_MAX`),
   `CBlastFastaInputSource` -> `x_InitInputReader` (`blast_fasta_input.cpp:317-368`): `CCustomizedFastaReader` (class at 64-103: `x_CloseGap` is a no-op, `AssignMolType` forced to the flags because the threshold is UINT_MAX) with flags
   `fNoParseID | fDLOptional | fAssumeProt | fNoSplit | fHyphensIgnoreAndWarn | fDisableNoResidues | fQuickIDCheck` (and `fParseRawID` instead of `fNoParseID|fDLOptional` with `-parse_deflines`); ignored problems: `ModifierFoundButNoneExpected`, `TooLong`, `TooManyAmbiguousResidues`; id generator `Query_<n>` starting at 1.
5. `CBlastFormat` is built, `PrintProlog()` (outfmt 0 writes the 21-line `BLASTP 2.17.0+ / Reference` text now, before any query is read).
6. Loop: `CBlastInput::GetNextSeqBatch` (`blast_input.cpp:135-170`) reads records until the sum of lengths reaches the batch size (`GetQueryBatchSize`, `blast_input_aux.cpp:70-141`: blastp 10000, tblastn 20000 residues), `catch (const exception&) { continue; }` drops a record that throws a non-parse exception (e.g. bad_alloc);
   `CObjMgr_QueryFactory` (`Empty CBlastQueryVector` if the batch is empty), `CLocalBlast::Run` -> `SetupQueries_OMF` (`blast_setup_cxx.cpp:485-660`): per query `GetSequenceProtein` (O->X, `Sequence contains no data` for length 0), messages with id `Query_<n> <title>` (cut to 25 chars + `.. ` above 35), `BlastSetup_Validate` (all queries invalid -> `BLAST engine error: Warning: Sequence contains no data ...` rc=3), then the search and `PrintOneResultSet`.
   A later batch can therefore fail after the earlier batches were printed.

`CFastaReader::ReadOneSeq` (`fasta.cpp:312-440`), as used here, per record:
- loop over lines of `CStreamLineReader` (`line_reader.cpp:100-296`: end of line style detected from the first line: LF, CRLF, bare CR, and mixed; CR inside an LF file switches style);
- if `PeekChar()=='>'` (column 0): `>?_x` is rewritten to `>x`; `>?...` is a data line (gap line); otherwise a defline (ends the previous record);
- otherwise the line is `TruncateSpaces_Unsafe`d (ASCII `isspace`); an empty line is skipped; first char `;`, `#`, `!` -> comment, skipped; if no defline was seen yet `ParseDefLine(">")` (fDLOptional: id `Query_n`, empty title);
- `ParseDefLine` -> `CFastaDeflineReader::ParseDefline` (`fasta_reader_utils.cpp:146-226`): `len<=1` or blank after `>` -> no title; with `fNoParseID` the id is not parsed at all (so every record gets the generated local id `Query_<n>`), the title starts after `isspace` bytes and ends before the first byte `< 0x20` (the first byte is exempt); the title is parsed again in `AssembleSeq` (`ParseTitle` 670-686: >1000 and mods ignored, `CreateWarningsForSeqDataInTitle` 1617-1681, `x_ApplyMods` `TruncateSpacesInPlace`);
- `ParseDataLine` (774-1016): `>?` -> `ParseGapLine` (1094-1349); otherwise `CheckDataLine` (710-772; only while no residue has been read) and the residue switch (856-931): `A-Z a-z *` kept (lowercase upper-cased and marked in the mask), ASCII white space skipped, `-` skipped with a per-line warning, `;` ends the line, everything else removed with one `Ignoring invalid residues` warning per line (positions 1-based in the trimmed line);
- `AssembleSeq` (1351-1516): gaps become `X` runs (`!fParseGaps`), empty sequences are allowed (`fDisableNoResidues`), title warnings are emitted now.
- Output side: `tabular.cpp:474-541` `s_ReplaceLocalId`: id = first `' '`-delimited token of the title, local id (`Query_n`, `Subject_n`) if the title is empty; for subjects `showdefline.cpp:213-245` replaces ids containing `lcl|Subject_` by the first word of `CDeflineGenerator::GenerateDefline` (protein, no title: `unnamed protein product` -> `unnamed`).

## Tables the port needs

- NCBIstdaa <-> letter: 28 codes, `blast_encoding.c:104-121` (`AMINOACID_TO_NCBISTDAA` 128 bytes), also `seqport_util.cpp:6308-6341`. LOSAT already has the identical table (`utils/matrix.rs:133`). Only valid reader letters reach it (A-Z, `*`; `O` is replaced by `X` first).
- Only for gap-line modifiers (optional, row 40): 10 gap types (`src/objects/seq/Seq_gap.cpp:176-195`), 12 linkage-evidence names (`src/objects/seq/seq.asn:383-398`), `CanonicalizeString` (`fasta.cpp:2129-2146`).
- Protein branch of `CDeflineGenerator` for the outfmt 0 subject title (row 44, owned by RP/RN): `create_defline.cpp` (see AUTHORITY §N of e2c for the nucleotide twin).

## What LOSAT has today

- BLASTP and TBLASTN read both FASTA files with `bio::io::fasta` (`blast_engine.rs:4052,4058`; `tblastn/args.rs:466`), pass `fasta::Record`s to the engines and use `record.id()` (first Unicode-white-space token, `unknown` when empty in BLASTP, empty in TBLASTN) and `id + " " + desc` as the title.
  They encode `record.seq()` bytes directly (`encoding.rs:55`, `utils/matrix.rs:133`), so every byte that the NCBI reader removes (digits, punctuation, blanks, hyphens, `;` text, comment lines) becomes code 0 or X inside the sequence: silent wrong alignments (e.g. a numbered GenBank block gives a 85.0% / 140 hit instead of 99.167% / 120, g_genbank).
- `algorithm/blastx/input.rs` + `native.rs` (BLASTX, base commit) already contain an NCBI-style reader and app flow with `protein: bool`: a transpile of `CStreamLineReader` (`FastaStream`, 1321-2030), the record loop (`parse_fasta_lines`, 373-780), `check_data_line` (891), residue/hyphen/gap/title handling, batches (`emit_fasta_batch` 1091), `stream_is_empty` (1229), and in `native.rs` the NCBI order (210-470: subject first, `Empty CBlastQueryVector`, `inaccessible`, `Query is Empty!`, subject-empty warnings, prolog, lazy batches), `reader_warning` (565) and the `FastaParseError` -> `BLAST query error:` rc=1 mapping.
  I read its protein branch against the oracle outputs of this range (it was not executed: no builds are allowed in this task); by reading it agrees with the NCBI behaviour recorded in the table except for the points below, which must be confirmed by a run after the generalisation. **The port for BLASTP/TBLASTN should generalise this code (move it to a shared module, give it a role Query/Subject) rather than transpile again.**
- `algorithm/blastn/input.rs` is the older "reject what bio reads differently" model (BLASTN, TBLASTX); it does not apply to proteins (`IUPAC_NUCLEOTIDE` residue check) but supplies `open_input`/stdin (387-420), `read_fasta_bytes` (423), `check_utf8_file_name` (559).
- `report/query_warnings.rs` has `invalid_query_warning` (used by TBLASTN) and the 35-char rule.

### Points where `algorithm/blastx/input.rs` must change before reuse (found by reading it against the oracle runs)

1. `new_record` names ids `Subject_n` when `protein` is true and `Query_n` otherwise; BLASTP/TBLASTN need `Query_n` for protein queries (a role parameter, row 27 and 42).
2. Title warning order: `parse_fasta_lines` calls `warn()` for the title at the defline (input.rs:688); NCBI emits it in `AssembleSeq` after the residue/hyphen warnings of that record (row 31, `i_order.P6`).
3. A title that is not UTF-8 fails (`invalid UTF-8 FASTA title`, input.rs:605) and `Record.title: String`; NCBI works on bytes (row 15). The engines print titles via `String`; decide whether to keep bytes (`Vec<u8>`) for titles.
4. Gap lines: only `>?N` / `>?unkN` with N a positive integer that fits are accepted; everything else (size 0, no digits, modifiers, text after the number) bails with `unsupported FASTA assembly-gap size/modifiers` (rows 39-41). The explicit rejection is a legitimate choice; row 40 gives the full NCBI behaviour if the session wants to port it.
5. The residue warning (`reader_warning`) caps at 1000 ranges like `ConvertBadIndexesToString`; checked, same.

## Surprises

- LOSAT BLASTP drops ALL hits when any subject record is empty (rc=0, empty output; `t1.S6`-`t3.S6`, `a_multi_emptymid.S6`). NCBI warns and searches the others (row 25).
- Every byte that the reader deletes silently changes the LOSAT alignment (rows 35-37), and the mid-line `;`, comment lines (`;`, `#`, `!`) and blank first lines are read as sequence or fail with bio's `Expected > at record start.` (rows 17-20).
- NCBI reads a line starting with `!` inside a record as a comment even when it looks like residues (`g_shortline.P6`: hit 5-64 instead of 1-120).
- The `# Query:`/`Query=` text, outfmt 0 subject titles and the query warnings use the title (control character cut, leading blanks trimmed), not bio's id+description (row 28).
- A subject with an empty defline is called `unnamed` / `unnamed protein product`, never `Subject_n` (row 42).
- `Query_n` ordinal counts every record (empty, gap-first, titled); an empty defline gives `Query_<ordinal>` in outfmt 6/7 (row 27).
- NCBI reader messages are plain stderr lines (no `Warning:` prefix); search messages (`Warning: [blastp] Query_2 q2 empty: ...`) have a trailing space and several messages of one query are joined on one line with single spaces (`h_onlyO.P6`).
- Reader errors say `BLAST query error:` even for the subject (`a_bom.S6`).
- `>?N` (assembly gap) is common in real nucleotide FASTA but rare in protein files; `>?4000000000` depends on memory (`k_gap_*`): bad_alloc is swallowed by `GetNextSeqBatch` and surfaces as `Empty CBlastQueryVector` rc=3 (row 41).
- Protein-type guessing never happens (UINT_MAX threshold): `ACGT...` queries are accepted silently (row 49).
- The web adapter's index scan (`web/adapter/src/scan.rs`, parser kind 0 = bio) must follow the new reader (kind 1 is announced for BLASTX in `docs/web/abi_v2.md:85`); `register` compares ids and lengths against the search parser.
- The engines take `bio::fasta::Record` (`run_local(.., query_records: &[fasta::Record], ..)` and ABI v1/v2); an NCBI record has three names (local id `Query_n`, title, printed id token). Bridging with `Record::with_attrs(id, desc, seq)` loses the empty-title case (`# Query: ` empty but id `Query_1`), so the record type of `blastx/input.rs` (`internal_id`, `title`, `sequence`, `lowercase_masks`, `warnings`) or an equivalent is needed.

## Decisions for the session

1. Reuse/generalise the BLASTX reader (recommended, rows with `port:S`) and the BLASTX app flow for BLASTP and TBLASTN; keep BLASTN/TBLASTX as they are.
2. Gap-line modifiers (row 40): port the two small tables (M) or keep BLASTX's explicit rejection. Gap sizes beyond Int4 (row 41): reject with an explicit LOSAT limit.
3. Non-UTF-8 file names (row 9): reject like BLASTN (`check_utf8_file_name`) or carry `OsString` into the report.
4. Title bytes (rows 15, 30): `String` with an explicit rejection of non-UTF-8 titles, or bytes.
5. Subject titles in outfmt 0 (row 44) belong to RP/RN; the id token/`unnamed` rule (row 42) is input-side and should be done with the reader.
6. `-lcase_masking` for BLASTP (row 45) needs the engine part from PA/NS (mask joins the SEG mask for the lookup table, shown in lowercase in the `Query` lines); TBLASTN's mask already works for plain lowercase letters and should be re-derived from the reader's mask intervals after the port (row 46).
7. All-ambiguous queries (row 47): BLASTP lacks the per-query invalid-statistics handling that TBLASTN has (warning, no `# N hits found`, `-1.00` Lambda block); owner PS/RP.
8. Environment variable `BLASTINPUT_GEN_DELTA_SEQ` and `.ncbirc` are not modelled (outside the scope rules).

## Not in this range (referred)

`-parse_deflines` (fParseRawID: ids are parsed as Seq-ids; `*A` row of PA), `-query_loc`/`-subject_loc` (S11), `-in_pssm`, outfmt 0 body (RP/RN), nucleotide TBLASTN subject (IN), per-query statistics (PS/NS), `-seg` interplay with lowercase masks (NS), `-num_threads` warning text.
