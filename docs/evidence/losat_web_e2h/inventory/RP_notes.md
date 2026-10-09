# RP notes: reports that print IDs and titles (outfmt 0/6/7, warnings, ABI)

Oracle: NCBI BLAST+ 2.17.0 (`/home/kawato/micromamba/bin`), LOSAT-native of f3048ffde. Every run had a `.ncbirc` with `[BLAST]` / `DATA_LOADERS=none` in its directory (all first lines were `>` lines). Inputs and outputs: `scratch_RP/` (`gen.py` builds `cases/<case>/{q,s}_{nuc,prot}.fa` with the title variants; `run_oracle.sh`, `run_losat.sh`, `run_mode.sh`; results `out/`, `out_losat/`, `out_qvar`, `out_svar`, `out_losat_qvar`, `out_losat_svar`; `cmp.py` compares; `warn/` warning cases; `multi/` multi-record file; `pd/` -parse_deflines; `bx/` BLASTX).

## Call path (what prints an ID or title)
1. Reader (range RD): each record gets a local ID from CSeqIdGenerator (`Query_N` / `Subject_N`, N = record position in that input, every record counts) and a Title descriptor = defline after `>` and skipped white space, cut at the first byte < 0x20, trimmed; empty title -> no descriptor.
2. Warnings at setup (stderr): subjects `blast_setup_cxx.cpp:773-788` (`Subject_N <title>: Subject sequence contains no data`), queries `blast_setup_cxx.cpp:533-543` (query id string, printed by `PrintOneResultSet` at `blast_format.cpp:1450-1452` before the query's report).
3. outfmt 6/7: `x_PrintTabularReport` (blast_format.cpp:759) -> `CBlastTabularInfo::SetFields` -> `SetQueryId` (s_ReplaceLocalId) and `SetSubjectId` (s_ReplaceLocalId, then `CShowBlastDefline::GetSeqIdList`). outfmt 7 header: `PrintHeader` -> `x_PrintQueryAndDbNames` -> `AcknowledgeBlastQuery(tabular)`.
4. outfmt 0: `PrintOneResultSet` (blast_format.cpp:1411): `AcknowledgeBlastQuery` (Query=), description list `CShowBlastDefline::x_DisplayDefline` (title from `CDeflineGenerator::GenerateDefline(fLeavePrefixSuffix)`), per HSP `CDisplaySeqalign::x_ShowAlnvecInfo -> x_PrintDefLine` (heading from GenerateDefline flags 0). bl2seq apps always run with m_IsDbScan, so the `Subject=` branch (blast_format.cpp:1497) is unreachable.

## State that persists
- ID counters per input (Query_N, Subject_N), never reset by empty or skipped records (multi/: 2nd query empty -> `Query_2`, 5th subject empty -> `Subject_5`).
- No state between query batches in the report code; the invalid-query warning uses the global query index.

## Key results (all from oracle runs)
- Tabular IDs (default fields `qaccver saccver`, and custom `qseqid sseqid qacc sacc qaccver saccver sallseqid sallacc`): first token of the title split only at `' '`; empty title -> `Query_N`/`Subject_N`; bytes verbatim. `stitle`/`salltitles` = `N/A`; `qgi sgi sallgi` = `0`.
- sseqid/sacc/saccver go through `GetSeqIdList`: when the replaced ID contains `lcl|Subject_` (empty title, or a first token that starts with `Subject_`), the ID is the first token of `GenerateDefline`. Consequences: BLASTP subject with empty title -> `unnamed` (from `unnamed protein product`) while sallseqid says `Subject_N`; `Subject_&amp;x` -> `Subject_&x`; `Subject_4.` -> `Subject_4`. LOSAT is wrong today on accepted input for the last two (v01, v09; all four programs). If the title is not decodable (non-UTF-8 with 0x81/0x8D/0x8F/0x90/0x9D) NCBI dies in outfmt 6 with exit 255, `CCoreException::eNullPtr`, stdout `Subject_\x81\t` (v06): proposed reject (stack trace has NCBI build paths); needs the batch question.
- `# Query: ` + title bytes (no ID, empty title leaves `# Query: `). `# Database: User specified sequence set (Input: <-subject as given>)`.
- `Query= ` + title wrapped by `NStr::Wrap` (bytes, width 68; may split a UTF-8 character), `\nLength=`. Empty title: `Query= \nLength=330` (no blank line). The Query= line is raw: no HTML decoding, no encoding change.
- Description list and heading use `GenerateDefline`: strip trailing `.,;~ ` (not when only those bytes), TPA/MAG prefix removal (heading only), `HtmlDecode` (entities; invalid UTF-8 treated as ISO-8859-1/Windows-1252 and re-encoded, valid UTF-8 untouched; undecodable -> exception), strip trailing `,;~ `, `x_CleanAndCompress` (a final byte >= 0x80 is dropped because of signed `char`). Description: 68-byte cut `...`, no ID. Heading: `> ` + title wrapped by `s_WrapOutputLine` (60 bytes, newline after the next white space), `Length=`.
- Exception path: `Unknown` description, no heading, `Sequence with id Subject_N no longer exists in database...alignment skipped` per HSP on stdout, exit 0 (u06).
- Warnings: `Warning: [prog] Query_N<space title>: <messages each + space>` with query_id cut to 25 bytes + `.. ` when longer than 35 bytes; subject empty-data warning `Subject_N <title>: Subject sequence contains no data`; raw bytes on stderr (no escaping of 0x01, 0x7f, 0xff). The empty-query error has no ID (`BLAST engine error: Warning: Sequence contains no data `, exit 3). Range errors print no ID/title.
- Same output in all four programs, except: protein subjects (BLASTP) get `unnamed protein product`; nucleotide subjects an empty title. `-num_threads` has no effect (ignored with -subject).

## LOSAT today
Accepted inputs equal NCBI (no DIFF except `Subject_` first tokens). Everything else is rejected by `algorithm/blastn/input.rs:100-132` (empty, leading white space, control byte, non-ASCII), bio's UTF-8 reading, `report/defline.rs` check functions (entities, subjects with hits), and the custom-field rejections. All places print/hold `String`/`&str`; the port needs bytes (`Vec<u8>`/`Arc<[u8]>`) for IDs, titles, Query=, headers, warnings; PairwiseHit.subject_title must carry the full title.

## BLASTX
`local_id` (report.rs:799) and the tabular/`# Query:`/Query=/warning paths are faithful for UTF-8 titles. Divergent: no GenerateDefline at all in description/heading (raw title: trailing punctuation, TPA prefix, entities, `unnamed protein product`, `Subject_` token not applied), `subject_title` split at the first space, titles are `String` (non-UTF-8 rejected). Verified in `bx/`.

## Decisions for the session
1. Port GenerateDefline + HtmlDecode + GuessEncoding + CleanAndCompress as one byte module (rows 30-39); reuse for BLASTX.
2. Reproduce the `Unknown` / `Sequence with id ... skipped` path (deterministic). Ask in a batch about rejecting the outfmt 6/7 exit-255 crash with a `Subject_` first token and an undecodable title.
3. `-parse_deflines`, `-html`, `-show_gis`, `-line_length`, `-num_descriptions` stay rejected.
4. Custom tabular fields of BLASTN/TBLASTX/TBLASTN stay rejected (not a reader matter); if ported, values above.
5. ABI: JSON `id` lossy UTF-8 only; streams 0/6/7/3 raw bytes; adjust `docs/web/abi_v2.md` (stream 3 'UTF-8') and the scan/store ID equality check.
6. Approved punctuation exception (PD-LOSAT-NCBI-DEFECTS 2) stays valid; the crash condition is evaluated on the decoded title bytes (entities cannot produce it; non-ASCII bytes exclude it).
