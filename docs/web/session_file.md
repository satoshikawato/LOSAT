# LOSAT Web session file (container 1, schema 1)

A session file saves the completed runs of a LOSAT Web working session so that they can be opened
again later, in another tab or browser, **without searching again** (plan §5.8, design §12.2,
S15 item 4). This document is the format's contract; `web/app/src/domain/session-file.ts` writes,
checks and reads it, and `web/app/src/application/session.ts` saves, opens and re-attaches.

## What it holds, and what it never holds

It holds, for each completed run: what was searched (the RunSnapshot without the bytes of its
inputs: program, argv, title, requested threads, group, queue time), what happened (the RunRecord
fields that Run details shows), the outputs 0, 6 and 7, the HSP records (ABI stream 1,
`docs/web/abi_v2.md` §8) and the diagnostics (stream 3), **byte for byte as stored**, and the
**identity of each input**: the SHA-256 and length of the engine input, the reader kind, the record
table (ID, length, SHA-256 of each record), and the sources with their exclusions. When "Include
candidates and notes" is checked (the default; S15 decision 1, Owner), it also holds the candidate
tray's candidates of those runs with their notes, in tray order.

It never holds the input FASTA (design §12.2: the dataset is not moved), storage (OPFS) paths, the
working session's token, or the IDs of runs, groups, revisions or sources. References inside the
file are positions in it. Queued, running, cancelled and failed runs are never saved.

The aligned rows of the HSP records are research data: a session file is not "without sequences".

- File name: `losat-session-YYYYMMDD-HHMMSS.losat-session.gz` (local time of saving); it never
  names an input. MIME type `application/gzip`.
- Compression: gzip, written with the browser's `CompressionStream('gzip')` and read with
  `DecompressionStream('gzip')`; no library. One gzip member; anything after it is refused.

## The container (the bytes inside the gzip)

```
LOSAT-WEB-SESSION 1\n
manifest <n>\n <n bytes: the manifest, JSON, UTF-8> \n
run1.out0 <n>\n <n bytes> \n
run1.out6 <n>\n <n bytes> \n
run1.out7 <n>\n <n bytes> \n
run1.hits <n>\n <n bytes> \n
run1.diagnostics <n>\n <n bytes> \n
run2.out0 <n>\n ...                      (five blocks per run, runs in the manifest's order)
candidates <n>\n <n bytes: JSON> \n      (only when the manifest's "candidates" is true)
LOSAT-WEB-SESSION-END\n
```

- The first line names the container and its version. Another first line is "not a LOSAT Web
  session file"; a higher version is "saved by a newer LOSAT Web".
- Each block is a header line `<name> <length>` (a decimal length without leading zeros), then
  exactly `<length>` bytes, then one line feed. Every block is present, in this order, even when
  empty (length 0). A run block's length must equal the manifest's `blocks` value for it.
- Lines (the first line, block headers and the end line) are ASCII and at most 64 bytes.
- `runK` is the run's 1-based position in the manifest's `runs`.
- Nothing may follow the end line.

## The manifest

A JSON object. Every field below is required unless marked optional; **a field that is not listed
is refused** (a later change adds a new schema). Numbers are whole numbers from 0 to 2^53 - 1 unless
said otherwise; times are milliseconds since the epoch.

| Field | Value |
|---|---|
| `format` | `"losat-web-session"` |
| `schema` | `1` (a higher one is "saved by a newer LOSAT Web") |
| `app` | `{ "version", "build" }`: the app's `package.json` version and the git short SHA of its build (`unknown` when the build had no git) |
| `savedAt` | when the file was saved |
| `candidates` | `true` when the file has a candidates block |
| `runs` | 1 to 1,000 runs, in run-number order of the saving session |

A run:

| Field | Value |
|---|---|
| `number` | the run's number in the session that saved it (1 or more); shown as "run N there" |
| `title` | optional; the Job Title (not empty) |
| `program` | `blastn`, `blastp`, `blastx`, `tblastn` or `tblastx` |
| `argv` | the run's argv: `[program, "-query", <query name>, "-subject", <subject name>, options...]`; the names must be the inputs' `name` |
| `requestedThreads` | `"auto"` or 1 or more |
| `group` | optional; `{ "index", "position", "size" }`: runs with the same `index` were queued together as separate searches; `position` ≤ `size` |
| `queuedAt` | when it was queued |
| `record` | the RunRecord, each field optional: `runtimePath` (`threaded`, `serial`, `fake`), `threads` (1 or more), `fallbackReason`, `engineBuild`, `runtimeGeneration`, `memory` (`linearBytesBefore`, `linearBytesAfter`, `instanceRuns`), `subjectRetained` (boolean), `startedAt`, `phaseTimes` (`preparing`, `running`, `finalizing`, each optional), `endedAt` |
| `hitCount` | the number of HSP records |
| `blocks` | `{ "out0", "out6", "out7", "hits", "diagnostics" }`: the length of each block of the run |
| `query`, `subject` | the input identity (below) |

An input identity:

| Field | Value |
|---|---|
| `name` | the `-query` / `-subject` name (a file's name, `query.fa`/`subject.fa` for pasted text, `combined_query.fa`/`combined_subject.fa` for joined files) |
| `sha256` | lower-case hex SHA-256 of the engine input (`InputSnapshot.sha256`) |
| `length` | bytes of the engine input |
| `reader` | the index scan's reader kind that made the record table: 1 (nucleotide) or 2 (protein); it must be `indexParser(program, role)` |
| `records` | columns `{ "id": [...], "length": [...], "sha256": [...] }` of equal length; element `k` is the record at position `k` of the engine input (the HSP records' `q_idx` / `s_idx`): its ID (may be empty), its length in residues, and the SHA-256 of its original bytes (`DatasetRecord.sha256`) |
| `sources` | 1 or more files, in the order the input joined them: `{ "name", "size", "records", "excluded" }`, `excluded` being the 0-based indices of the records left out, ascending, each below `records`. The included records of all sources must be exactly as many as `records` lists |

## The candidates block

`{ "candidates": [ ... ] }`, in tray order. Each candidate is `{ "run", "index", "qIdx", "rank",
"note", "addedAt" }`: `run` is the run's 1-based position in the manifest's `runs`; `index` is the
HSP's index (its HSP record), below the run's `hitCount`; `qIdx` is below the run's query record
count; `qIdx` and `rank` must be those of the HSP record at `index`; the same HSP appears once;
`note` is the user's note (text). Only candidates of saved runs are written.

## Limits

| Item | Limit |
|---|---|
| manifest block, candidates block | 64 MiB each (refused from the header, before reading) |
| runs | 1,000 |
| argv | 1,000 words of at most 10,000 characters |
| title, record ID, engine build, fallback reason, app version and build | 10,000 characters |
| file and input names | 1,000 characters |
| sources of one input | 10,000 |
| candidates | 1,000,000; a note at most 100,000 characters |

Saving checks the session with the same checks as opening, so a saved file always opens; a
session over a limit is not saved, and the message says which value.

## Opening

The file is read incrementally: the decompressed chunks go through a state machine that reads the
header lines, collects only the manifest and the candidates (within their limits), and sends each
run's blocks straight to the Data worker over the run output channel (`ports/run-output.ts`), as
the engine's output would arrive, with backpressure; nothing else is held. Once the Data worker
cannot store a run (the storage ran out), the backpressure wait fails and loading stops there: the
rest of the file is not decompressed, and the file is refused with the reason. After a run's last block
the run is committed and its HSP records are checked against the manifest. A loaded run is in the
RunStore like a searched one: the results screen, the Alignments, the outputs and the alignment
export work unchanged.

**Nothing is searched.** The runs join the working session as completed runs, with new run IDs and
the next run numbers of this session; they are never queued, never validated, and the engine is
never asked. The queue, the results header and Run details say "from `<file>`, run N there".
The run's verification badge (Run details, the JSON export's `run.verification`, the report) is
this site's only when one of this site's engine builds wrote its outputs (the run record's
`engineBuild`), and then adds that the run was loaded from a session file; otherwise it is
"Written by another engine build" and names that build, the LOSAT Web that saved the file (`app`)
and this site's builds. A FakeEngine run stays a development run. A
loaded run's InputSnapshot has no engine bytes and no revisions (`bytes` is absent); Run details
shows the input's size from the manifest. Candidates are added to the tray with their notes and
times, rebuilt from the loaded runs' HSP records (read 1,000 at a time; the tray keeps no aligned
rows) and outfmt 6 rows: the file gives only the reference.

A file is **refused** with a message that says what is wrong and where, and then nothing is loaded:
every run staged or committed for it is deleted, and the runs and candidates join the working
session only after every check has passed and the tray's entries are made:

- not gzip; damaged or cut gzip data (a changed byte, a cut file, data after the gzip member);
- a wrong first line; a newer container version or schema;
- a manifest or candidates block that is not UTF-8 JSON, misses a field, has an unknown field or a
  value of the wrong type, or is over a limit (the message names the field, e.g.
  `runs[1].query.sha256`);
- blocks missing, extra, out of order, of another length than the manifest gives, longer than the
  header states, or a file that ends anywhere before its end line, or has bytes after it;
- HSP records that the run cannot have, checked in the Data worker as the JSON of their lines,
  before anything coerces them (the results' typed arrays, the exports): another count than
  `hitCount`; a line that is not a JSON object; a field missing or of the wrong type (`index`,
  `q_idx`, `s_idx` and `rank` whole numbers of 0 or more; coordinates whole numbers of 1 or more;
  frames null or -3..3 other than 0; scores numbers; `subject_length` null or a whole number;
  aligned rows null or text; `out6`/`out0`/`out0_subject` null or `[start, end]`); an index
  outside 0..count-1 or repeated; a `q_idx`/`s_idx` beyond the record tables; a rank repeated
  within a query; a range outside its output. Fields that a later ABI adds are left alone;
- candidates that name a run or an HSP that the file does not have, or whose `qIdx`/`rank`
  disagree with the HSP record at `index`.

Text from the file (names, IDs, titles, argv, notes, the file's own name) is shown as text only,
never as HTML, a script or a URL.

## Re-attaching the original FASTA (REQ-23)

Only the extraction of original residues (hit regions, flanks, complete sequences) and the download
of a run's input FASTA need the original files; the outputs, the HSP records, the results screen and
the export of aligned rows do not. For a loaded run, Run details shows per role that the original is
not attached and what that prevents, with "Choose the original query/subject FASTA…" (several files,
chosen together in any order, for a joined input). The extract form names the runs whose original is missing as soon as
such a candidate is selected, and the extraction is refused before anything is read.

Re-attachment is never automatic (a file of the same name in the search form, or one attached to
another run, attaches nothing). The chosen files must be as many as the recorded sources; each is
indexed with the recorded reader kind and matched to a recorded source by its records, whatever
order the files were chosen in (file dialogs rarely let one order a selection): a file is a
source's when it has the recorded number of records and, with the source's exclusions applied, the
recorded IDs, lengths and SHA-256s of the records that the input has from it. When a source has no
file, the message names a file's record count or the first record that differs. The run input is
then made in the recorded order with the recorded exclusions, and must have the recorded SHA-256
(which also covers lines before or between records). Only then is the input attached to that run
and role (the file names listed in the recorded order); otherwise nothing changes. The sources and
record tables of a refused attempt, and of an original that a new attachment replaces, are released
from the Data worker. The attachment is kept apart from the RunSnapshot (`RunView.attached`), which
stays as the file recorded it. The download of the run's input FASTA rebuilds the input from the
attached original and refuses it unless its SHA-256 is still the recorded one.
