# LOSAT Web ABI v2 (engine adapter ↔ application)

Status: **version 13** (§10), implemented in `web/adapter/` (version 2 in session S05). Drafted in S01
([`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md) §4, §7). The application
depends on the TypeScript ports under `web/app/src/ports/` (the searches and the HSP record in
`engine.ts`, the index scan in `scan.ts`, the inflation in `inflate.ts`), not on the exported
functions. The HSP record (§8) and the JSON responses (§9) cross into the
application, so any change to them updates the port that types them in the same
commit, and an incompatible change increments the version (§10).

Programs: BLASTP, TBLASTN, BLASTN, TBLASTX and BLASTX (since session SX); the tool `windowmasker` (`losat_web2_tool`, since version 13).

## 1. Purpose

The application (`web/app`) runs LOSAT searches through a Wasm module built from the
adapter crate `web/adapter`, which links the engine crate `LOSAT/` as a library.

- Search conditions are one CLI argv, read by the same argument reader as the CLI (NCBI's argument
  reader since line RG item K4i-2b, 2026-11-23).
- One run writes every output format the program supports (listed by `describe`),
  the HSP records and the diagnostics from the same search result (§7), and the output of
  the argv's own `-outfmt` when it is another format that the CLI writes, or 6 or 7 with fields
  (§5, §7).
- All BLAST-defined values and all compatibility text come from the engine.

ABI v1 (`losat_web_*` in `LOSAT/src/web_api.rs`, used by gbdraw) stays in the same
module unchanged, apart from fail-fast defect fixes (plan TD-1). The two ABIs share no
state.

## 2. Artifacts and module kind

| Artifact | Target | Cargo features of `LOSAT` | Notes |
|---|---|---|---|
| `losat-web-serial.wasm` | `wasm32-wasip1` | none | Serial. The engine rejects `-num_threads` above 1 |
| `losat-web-threads.wasm` | `wasm32-wasip1-threads` | `parallel`, `wasm-threads` | Needs a cross-origin isolated page and shared memory |

Both are WASI preview1 **reactors**: the host calls `_initialize` once after
instantiation, then calls the exports below any number of times. The build must match
the certified command-WASI builds in profile, rustflags, link arguments and shared
dependency versions (plan TD-6). `web/adapter/tools/check_build_identity.py` checks
this, and `web/adapter/tools/build_reactors.py` builds both artifacts and records their
identity. The reactor start-up (`crt1-reactor.o`, `--entry=_initialize`) comes from
`LOSAT/build.rs` through the `LOSAT` dependency: Cargo passes that build script's
`rustc-cdylib-link-arg` to the adapter's cdylib as well, so the adapter has no build
script of its own (a second copy would define `_initialize` twice).

The reactors' bytes depend on where they are built. `LOSAT` is a path dependency
outside the adapter's workspace, so rustc gets its sources by absolute path, and every
crate from the registry is compiled from `CARGO_HOME`. `build_reactors.py` replaces both
prefixes in the source paths that the modules embed (panic locations) with `/losat` and
`/cargo` (`--remap-path-prefix`, plan TD-11), so a module carries no path of the build
machine. The checkout path still changes the bytes: Cargo's metadata hash of the path
dependency, which appears in the symbol names of the name section, includes it. Two
builds of the same commit from the same checkout path are identical; the identity record
(`losat-web-*.json`) records the checkout path, `CARGO_HOME` and the remapping.

## 3. Imports

| Module | Name | Used by | Purpose |
|---|---|---|---|
| `wasi_snapshot_preview1` | `fd_write`, `fd_close`, `fd_fdstat_get`, `clock_time_get`, `random_get`, `sched_yield`, `environ_get`, `environ_sizes_get`, `proc_exit`, and the filesystem functions `path_open`, `path_create_directory`, `path_filestat_get`, `fd_prestat_get`, `fd_prestat_dir_name` | both | stderr, clocks, random, environment (the host gives an empty environment: variables such as `BATCH_SIZE` change the CLI's results, whose parity is defined without them). The filesystem functions read the files that the argv names, in the run's working directory, which the host preopens as `/work` (the first preopen, descriptor 3; `RUN_DIRECTORY` in `web/adapter/src/run.rs`), read-only: the files of the database of `-db`, its `taxdb.*` and `taxonomy4blast.sqlite3`, the lists of its filters (`-seqidlist`, `-taxidlist`, `-gilist` and their negatives), and the files of the database of repeats of BLASTN's `-filtering_db`, with `-db` or with `-subject` (line RG item D4; a run of `-subject` has its subject registered as always, and the repeat database is looked up by the name as `-db` is). The module makes `/work` its current directory at its first call (a host that preopens no `/work` leaves it `/`), so a relative name (`db/nt`) is a file of the run's directory, as on the CLI run from that directory, and an absolute name (`/db/nt`) names a file outside it, where the module has none: the run fails as the CLI run from that directory fails for a file that its root does not hold (`No alias or index file found ...`). The CLI's messages that name its current directory name `/work` (for example the search path of `No alias or index file found for protein database [x] in search path [/work::]`, and when BLASTN's `-filtering_db` repeat database files are missing from `/work`, the batch warning `Warning: No alias or index file found for nucleotide database [x] in search path [/work::]`). A host looks a path up one name after another, as the OS does: a name after a file is `ENOTDIR`, a name after one that is not there `ENOENT`, `..` included (`nofolder/../db/nt` fails as on the CLI), and `..` past the directory is `ENOTCAPABLE` (review of D12a, low 4, 2026-11-21). The engine opens a file when a search needs it (a database once per run, BLASTN's `-filtering_db` database at each query batch) and may read a descriptor on a thread other than the one that opened it (the threaded module's search threads), so a host of the threaded module gives every thread the same descriptors (`fd_read`, `fd_pread`, `fd_seek`, `fd_tell`, `fd_filestat_get`, `fd_close`, and `path_open`, `path_filestat_get` from `.`). A host that gives a run no file preopens an empty directory, or none: the run then fails where the CLI fails for a file it cannot open. The browser host (`web/packages/losat/src/infra/reactor/run-directory.ts`) reads each file in place, the bytes that each read asks for, from the user's `File` with `FileReaderSync`; the Node hosts of the tests preopen a directory of the OS (line WQA D12a, 2026-11-20). In the threaded module, the shared-memory guard also calls `sched_yield` between the slices of a wait once the host ended the memory's waits (or while its table of waits is full) (`LOSAT/tests/wasi_shared_memory.js`, `endWaits` in `web/packages/losat/src/infra/reactor/shared-memory.ts`): WebKit ends a worker that the page ends only when its thread is in a host function |
| `wasi` | `thread-spawn` | threads | starts a rayon worker thread |
| `env` | `memory` | threads | shared linear memory, created by the host. The module declares a maximum of 65,535 pages (4 GiB less one 64 KiB page; plan TD-7, DW-28); the host uses that maximum or, where it cannot be reserved, a smaller one (below) |
| `losat_host` | `emit(stream: u32, ptr: u32, len: u32)` | both | receives output bytes (§5). Worker instances of the threaded module import it too, but never call it |

The shared memory's maximum. The threaded module is linked with `--max-memory=4294901760`
(`LOSAT/.cargo/config.toml` and its copy `web/adapter/.cargo/config.toml`). One page less
than wasm32's 4 GiB keeps the byte size of the memory within 32 bits, in the module's C
library and in the host. A host must not provide a memory with a larger maximum than the
declared one; it may provide a smaller one. An engine reserves the maximum as address
space when the memory is created. If it cannot, it either throws `RangeError`
(JavaScriptCore without a fast memory) or keeps a smaller reservation without an error,
down to the initial pages alone, so that `memory.grow` fails (V8 without guard regions,
as on Android or under an address-space limit). The hosts of this repository therefore
try 65,535, 49,152, 32,768, 24,576 and 16,384 pages: they create each candidate one page
short of the module's initial size, grow that page, and use the first candidate for
which both steps succeed. Below 16,384 pages (1 GiB, the maximum declared until plan
DW-28) a host does what it did before: LOSAT Web runs the search on the serial module,
and the Node host fails. The serial module defines its own memory without a maximum; that
memory is not reserved in advance and can grow to wasm32's 4 GiB as well.

## 4. Exports

All functions return `i32`: `0` (or a non-negative handle) on success, `-1` on failure.
On failure, `losat_web2_last_error_ptr/len` hold a UTF-8 message. Engine errors use the
CLI's wording, which is NCBI's in LOSAT's words (`LOSAT/src/wording.rs`, AGENTS.md rule 8): every
message that this document quotes as NCBI's (`BLAST query/options error: …`, `BLAST engine error: …`)
reaches the host with LOSAT's prefix (`LOSAT query/options error: …`, `LOSAT engine error: …`), and a
printed version states LOSAT's (`blastn: LOSAT 0.2.0`) where NCBI's states `2.17.0+`. Argv errors are the command line's standard error for the argv, ending with a line end
(`losat::web_adapter::argv_message`): NCBI's short usage, a line of 72 `=`, an empty line and NCBI's
two posts (`Error: <message>` and `Error:  <message>`); with `-logfile NAME` a run's error is NCBI's
short usage, its line `Error in command-line arguments. See error logs for more details.` and the line
of `=`, and the two posts are in the run's file of `-logfile` (stream 8), as on the CLI, while
`validate`, whose start keeps no file, keeps them in its message, as the CLI writes it without
`-logfile` (line RG item K4i-2b round 2); or, where LOSAT refuses an NCBI option that it does
not run, the refusal text (`the NCBI BLAST+ option -NAME is not supported by LOSAT's PROGRAM`, without the
command line's `Error: `). The application marks the fields that the posts name. Search errors are the error and its causes joined by
`: ` (as ABI v1 reports them). The adapter reads the argv exactly as the command line
(every program accepts the CLI's default `-outfmt 0`), so an argv error of `validate`, an
unknown program's included, is the CLI's message (session S08; line RG item K4i-2b). Only the adapter's own
rules give other messages: `-out` is not accepted (an argv whose parse stops at a help word prints the help, as the CLI
does, `-out` or not; line WQA's D14 review, 2026-11-20), an `-outfmt` of a format other than 0, 6
and 7 is refused (§7). `run` finds the
handles of the argv's program before it parses the argv, so it reports an unknown program
as `unknown program '<name>'`; a host validates the argv first.

| Export | Arguments | Result |
|---|---|---|
| `losat_web2_abi_version()` | — | `13` (§10) |
| `losat_web2_alloc(len)` / `losat_web2_dealloc(ptr, len)` | — | memory for inputs (§6) |
| `losat_web2_describe(program_ptr, program_len)` | program name | emits a *describe* JSON on stream 2 |
| `losat_web2_validate(argv_ptr, argv_len)` | argv (§7) | `0` if the argv is valid, with a *validation* JSON on stream 2 (§9: whether the run reads its inputs with NCBI's data loaders on, which `register` and `scan_begin` take for the run's inputs), else `-1` with the CLI error. Since version 13 it also takes the argv of a tool run (`losat_web2_tool`: argv[0] `windowmasker`): it makes the checks of the command line before the tool reads its input or opens a file: the application start, which reads the registry files of the run's directory as the command line reads them in its directory (`.ncbirc`, the file of `-conffile`; §3), NCBI's argument reader and its checks, LOSAT's settings, windowmasker's mode (`one of '-mk_counts', '-convert' or '-ustat <stat_file>' must be specified`) and the checks of the mode that come before its files: `-mk_counts` with both `-ids` and `-exclude_ids` (NCBI's exception, exit status 3) and, without an id list, with `-infmt seqids` (LOSAT's refusal, exit status 1), and `-convert -in -` (NCBI's exception). An error that the command line finds after it opened a file is the run's, not `validate`'s: the masking opens `-in` before it reads `-ids` and `-exclude_ids`, and `-mk_counts` reads an id list before it refuses `-infmt seqids` (line NP WW3a's review, M3 and L4, D-WW3-10, 2026-11-23). It emits nothing, and its error is what the command line writes for that argv: its standard error, then the file of `-logfile` when the argv names one, without the last line end (for example the usage, the line of `=` and `Error: …` of an unknown argument); a help or version word is valid (line NP WW3a, D-WW3-7, 2026-11-22). The *validation* is emitted as soon as the start has read the argv's registry and LOSAT has taken its settings, before the checks of the arguments and the options, so a call that returns `-1` for one of those checks emits it too; one whose start or settings fail emits none (line WQA's D13 round 3, L1, 2026-11-20). For BLASTN it also checks the scoring options as NCBI does before a search (the checks that do not depend on the query), and returns NCBI's message (`BLAST query/options error: …`), or LOSAT's message with `not supported by LOSAT's BLASTN` for options that NCBI runs and LOSAT does not (`docs/evidence/losat_web_e2c/AUTHORITY.md` §G). For TBLASTX it checks the options alone, with the messages of a run: NCBI's `-evalue` check (`BLAST query/options error: expect value or cutoff score must be greater than zero`). For BLASTX it checks the options of each search's formats (§7) as the run checks them before it searches (`blastx::check_options`), with the messages of a run: LOSAT Web's refusals, LOSAT's BLASTX limits (for example `-matrix IDENTITY` with the compressed lookup table of `-ungapped -word_size 5`, with which NCBI crashes) and NCBI's errors of the options (`BLAST query/options error: …`); as for BLASTN, the argv is valid when one search is. For every program it reads `-subject_loc` and `-query_loc` as NCBI's `ParseSequenceRange` does, in the command line's order (for BLASTN, BLASTP, TBLASTN and TBLASTX the subject range before the program's option checks, the query range among them; for BLASTX both after LOSAT's limits of its options, the subject range first, then the query range, then NCBI's errors of the options, as its run reads them), with NCBI's messages (`BLAST engine error: Invalid specification of query location (…)`, `subject location` for the subject), and for a part that NCBI cannot convert to an int (`a-10`, `1-1.5`) NCBI's report of its exception (`Error: NCBI C++ Exception:` and the `CStringException` of `NStr::StringToInt`, which the CLI ends with exit status 255; env parity H). Whether a range fits the records is checked by `run`: a start past a subject record's end (NCBI's `Invalid from coordinate (greater than sequence length)`), a start just past a record's end (LOSAT's rejection; BLASTX searches that record as a sequence without data, as NCBI does, with NCBI's `Sequence contains no data` or `Subject sequence contains no data` warning), every query record of a batch skipped (NCBI's `Empty CBlastQueryVector`). The CLI reads the subjects before it checks the query options, so for an argv and inputs with both kinds of error the CLI reports the record's error and `validate` the argv's (`docs/evidence/losat_web_e2d/AUTHORITY.md` §A, §B). The checks of the inputs come with `register` and `run` |
| `losat_web2_start(argv_ptr, argv_len)` | argv as for `validate` | emits a *start* JSON on stream 2 (§9): what NCBI's application reads of the argv before its arguments (line RG item E1w-3a): the file of the last `-conffile`, the file of `-logfile` (`log_file`: NCBI's scan of it, the last `-logfile FILE` or `-logfile=FILE` before `--`, whatever it names; a run gives a stream 8 file only for a name that makes a file, not for `-`, `""`, `/dev/null`, a directory or the standard error), whether the argv has `-dryrun`, and whether a run of it reads its inputs together (`register_inputs`: an argv of `-parse_deflines` that is not a dry run; line WQA's sweep, D7, 2026-11-19). A host calls it before `validate` to learn which registry file to register (kind 2 of `file_begin`), under which name a run gives the file of `-logfile` (stream 8), whether a run is a dry run, and whether it registers the run's inputs with `losat_web2_register_inputs`. Since version 10 it also says what a run of the argv prints in place of a search (`prints`; line WQA's D14, 2026-11-20): `"version"` for a version word (`-version`, `-version-full`, `-version-full-xml`, `-version-full-json`) anywhere in the argv, which the CLI looks for before it reads a registry file (`conffile` is then null), and `"help"` for an argv whose parse stops at a help word (`-h`, `-help`, `-help-full`, `-xmlhelp`; `--help` is NCBI's `Unknown argument: "-help"` since line RG item K4i-2b): the words before it are read, so `-evalue abc -h` is the parser's error and `-h -evalue abc` the help; `dry_run` is then false, as the CLI prints before a dry run ends. It fails for a program that ABI v2 does not run (`unknown program` and `not available in the browser`, as `run`). Every `validate` and `run` makes NCBI's application start itself, as the CLI does before it reads the arguments: in an environment of its own (no variable but `NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH=windowmasker`, the directory of `-window_masker_taxid`, set where NCBI's environment layer comes before its registry files; no home directory, no program directory), with the registered registry files (kind 2) as the files of the CLI's working directory (the registry search path is `.` then `/etc`; a registered name is compared as a path of the CLI's working directory, a leading `./` being the same file): `-conffile NAME` reads the registered file `NAME`, also when the argv writes `./NAME`, in place of `<program>.ini`, and a registered `<program>.ini` or `.ncbirc` is read as the CLI reads one in its working directory. Its errors are the CLI's first ones, with the text that the CLI writes to its standard error without `-logfile`, without its last line end (§5, stream 3): a `-conffile` that names no registered file (`Critical: (106.23) Registry file "NAME" not found, aborting.`, exit 255), a registry file that NCBI reports (`Error: (CRegistryException::eSection) ...`, exit 2), a setting that LOSAT refuses. Then the argv is parsed (the three words are taken), the settings that LOSAT refuses are refused, and a dry run ends: `validate` returns 0 and `run` emits no output (stream 3 and 8 aside) and returns 0. A version word ends the start after the file of `-logfile` is made and before a registry file is read, and a help word ends the parse (the registry file, read before it, may still fail it): `validate` returns 0, and `run` emits the CLI's standard output on stream 11 (§5; streams 3, 8 and 10 aside), no output format, and returns 0 (line WQA's D14, 2026-11-20). The order is the CLI's: start, argv parse, refusal of settings, end of a dry run, and only then `-outfmt` is read (as NCBI reads it in `Run`), so `-outfmt x -dryrun` passes, as on the command line. After `-outfmt`'s number, an argv without `-subject` and `-db` fails with NCBI's error (`BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified`), before the other checks, as on the command line, whatever inputs the host registered: a help word that is an option's value is that value, as the CLI's parser reads it, so `blastp -query -h` names the query `-h` and no subject (line WQA's D14, 2026-11-21). A dry run's `run` takes the registered query and subject handles as every `run` does, but reads no record. A `-logfile` name is read as the CLI on Linux reads it in a working directory where the name's directories exist, as a `-conffile` name is a path of the CLI's working directory (line WQA's sweep, D11, 2026-11-20; a name with a directory was refused before), and the CLI reads it as NCBI 2.17.0 does (its round 2): the spaces at the ends of a relative name are cut, the name is put after the working directory and read lexically (`//` and `/./` as `/`, `..` taking the component before it away, links not followed), cut at 2047 bytes, and the `/`s at its end name no directory (`logs/` is the file `logs`); NCBI checks only that the directory of the name as written exists. `-logfile NAME` makes the file `NAME` of the run's files (stream 8), under the name as written (`run.log`, `logs/run.log`, `./run.log`, `/tmp/run.log`, ` run.log`, `logs/`, `/proc/self/cwd/run.log`): it holds what the CLI writes to that file (the start's messages, and in TBLASTX the run's messages and NCBI's error that ends the run, which then are not on stream 3), as the CLI's run of the argv writes it: what the run's other formats add to its messages goes to stream 3 (§5, stream 3; line WQA's sweep, D3 review r2, 2026-11-19). `-logfile -` is stream 3, `-logfile ""` and `-logfile /dev/null` send the start's messages (and TBLASTX's run messages) nowhere. The system then opens the name as Linux does for a user's command line with no terminal, whose standard output and error are pipes and whose standard input is the null device, and the messages go: to stream 3 for a name that NCBI cannot open, as NCBI sends a log it cannot open to `cerr` (a directory: `.`, `..`, `logs/.`, `/tmp/`, `/proc/self/cwd`; a file descriptor that the run has not opened: `/dev/fd/3`; a file of procfs or sysfs that a user may not write or that does not exist: `/proc/self/status`, `/sys/...`; a terminal: `/dev/pts/0`, `/dev/tty`; a name in a file: `/dev/null/x`, `/dev/stdout/../stderr`; a directory of `/`s only, which NCBI reads as none: `//run.log`; a name in a directory of root's, where a user makes no file, or a device that a user may not open: `/run.log`, `/etc/run.log`, `/var/log/run.log`, `/home/run.log`, `/dev/run.log`, `/dev/kmsg`, `/dev/console`; round 3 of D11, 2026-11-20; a directory that only root makes and Linux does not have, which the system does not find, also behind a link that `..` follows: `/data/run.log`, `/bin/../tmp/run.log`, as `/bin` is `/usr/bin` and `/usr/tmp` is none; a component of more than 255 bytes or a directory of 4096 bytes or more, which the system refuses; round 4 of D11, 2026-11-20); nowhere for a file in a file system that reports no free space, where NCBI's log gets no message (it checks for 20 KiB free before it writes: `/dev/fd/2`, `/dev/fd/1`, `/proc/self/fd/2`, `/proc/thread-self/fd/2`, `/dev/fd/../stderr`, whose directory is `/proc/self`, `/proc/self/comm` and the other files of procfs that a process may write, a message queue `/dev/mqueue/q`), for a device whose writes are lost (`//dev/null`, `/dev/full`, `/dev/zero`, `/dev/random`, `/dev/urandom`, `/dev/ptmx`, `/dev/fuse`, `/dev/net/tun`) and for `/dev/stdin`; to stream 3 as a file for the standard error (`/dev/stderr`, also `//dev/stderr` or `/dev/stderr/`): the messages go there as with `-logfile -`, and the file also gets the start's messages that the CLI writes to `cerr` at once (a registered `.ncbirc` that NCBI cannot read: `Error: (110.4) Error reading the registry ...`), which stream 3 then holds twice, as the CLI's standard error does. The standard output (`/dev/stdout`) is a file of that name on stream 8: the CLI writes its bytes into its standard output among the report's, where their place depends on the CLI's buffering of the report, and ABI v2 gives each output as the CLI's bytes of that format alone (§5), so the host gets the messages apart from the outputs. Other names are files of that name where a user makes files, with the directories that the name gives there: in the working directory, in a user's home (`/home/<user>/...`), and in the directories that the user writes on the engine's reference host (`/tmp`, `/var/tmp`, `/dev/shm`, `/run/shm`, `/run/lock`, also as `/var/lock`, `/var/crash`, the user's `/run/user/1000`, and WSL's `/mnt/wsl`, `/mnt/wslg`, `/mnt/c` and `/mnt/d`). The reference host is Ubuntu 24.04.3, kernel 5.15 (WSL2), as its first user (uid 1000, in that user's groups), and every list of directories and procfs entries is measured on it as that user (round 7 of D11, 2026-11-21). Every other directory is root's: only root makes entries in `/` (round 4 of D11, 2026-11-20). A name may pass a directory of root's that the engine does not list on its way back out (`/etc/default/../../tmp/run.log` is `/tmp/run.log`, as on the CLI), but must end in one that the engine knows (`/bin/../tmp/run.log` ends in `/usr/tmp`: stream 3), and looks nothing up in a directory that the user cannot search, as the system refuses it (`/root/../tmp/run.log`, `/var/cache/ldconfig/..`, another process's `/proc/1/fd/..`: stream 3; rounds 6 and 7 of D11). The engine carries no list of the system's files, so a name that climbs with `..` out of a file or out of a missing name in a directory of root's is a file of the run in the Web (stream 8), where the CLI on Linux refuses the path (`ENOTDIR`, `ENOENT`) and writes the messages to its standard error: `-logfile /etc/passwd/../../tmp/x` and `-logfile /etc/zzz/../../tmp/x` (D11, 2026-11-21). `/proc/thread-self` is a link to a task of the run's process, as on Linux, and the system follows at most 40 links in a lookup (`/proc/self/root` 21 times: stream 3; round 7). The files of procfs are those of the reference host's Linux 5.15 |
| `losat_web2_register(program_ptr, program_len, role, data_loaders, bytes_ptr, bytes_len)` | role `0` query, `1` subject; `data_loaders` `1` on, `0` off, as the *validation* of the run's argv gives it (§9); original FASTA bytes | handle ≥ 1; emits a *register* JSON on stream 2. A handle belongs to the program that registered it and to its setting of the data loaders: a run whose argv reads with them on refuses a handle registered with them off (`query handle N was read with the data loaders off (DATA_LOADERS=none), which this run has on`), and takes one registered with them on in a run with them off (the readers differ only in a first line that may be a Seq-id, which the reader with the loaders on refuses). The records are those of the program's own reader, the engine's port of NCBI BLAST+'s `CBlastFastaInputSource` (session SF): nucleotide flags for BLASTN, TBLASTX, the TBLASTN subject and the BLASTX query, protein flags for BLASTP, the TBLASTN query and the BLASTX subject, the bytes read as a file, the data loaders on or off as given (on is the CLI without `.ncbirc` or `DATA_LOADERS`; `DATA_LOADERS=none` in the `[BLAST]` section of a registry file that the run's start reads, such as the file of `-conffile`, turns them off, as in the CLI; line WQA's D16, 2026-11-20). Every record that the reader returns is registered, records without residues and a first record without a defline included; the reader's warnings stay with the records and come on stream 3 of a run, where the CLI writes them. `register` first scans the bytes with the reader of the program's input (`scan_begin`, §9) and fails with the scan's error, which is the reader's and the CLI's: LOSAT's refusals, which the CLI makes too and writes after `Error: `: with the data loaders on, a first line that NCBI's data loaders would fetch as a Seq-id (`the first line of the query ("2HHB") is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT's BLASTP (start the query with a '>' defline, or set DATA_LOADERS=none in the [BLAST] section of .ncbirc)`), and a record over 2147483647 letters (`a query record longer than 2147483647 letters is not supported by LOSAT's BLASTN`); and for a subject the reader's error, with the CLI's text, which the CLI writes as it is (`BLAST query error: CFastaReader: Near line N, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.`), as the CLI reads every subject before its first query. A query is read as the CLI's query batches read it (line WQA's D10, 2026-11-19): the records before the reader's error are registered with the error, which *register* gives as `reader_error` (§9), and a run searches the query batches before the batch that reads the record, then ends with the error, as the CLI does (the host keeps what the run's own output held at the error, §5). A `>?` gap line and the lines that NCBI's line reader joins (a CR line end in a file of LF line ends) are read as the CLI reads them (line WQA's D13, 2026-11-20; version 4 refused them as `not supported in the browser`). A query with blank and comment lines only (`!`, `#`, `;`) has no record and is registered: a run of it ends with NCBI's `BLAST engine error: Empty CBlastQueryVector`, as the CLI's first batch without a query. A query of white space only has no record: a run gives NCBI's `Query is Empty!` warning. A subject without records is registered: a run gives NCBI's `BLAST engine error: Empty CBlastQueryVector`. `register` then checks that the scan found the reader's records (§9). `register` judges one input alone, in the order the host calls it, and never looks at the other role's input. The CLI reads the subjects while it processes its arguments, before it reads any query: a subject that fails to read, and a subject without records (`BLAST engine error: Empty CBlastQueryVector`, exit 3), end the CLI before the query is read, in BLASTN, BLASTP, TBLASTN and TBLASTX (checked with NCBI BLAST+ 2.17.0; audit finding D-1 of session SFc), and in LOSAT's BLASTX (its command line and `run_local`). So for a query and a subject that both fail, the CLI reports the subject's error, while a host that registers the query first sees the query's. A host that wants the CLI's first error registers the subject first, stops at its error, and takes a subject registered without records as the run's `BLAST engine error: Empty CBlastQueryVector` before it registers the query. A registration gives bytes, never a name: the names of `-query` and `-subject` in a run's argv are labels, which the run writes where the command line writes a file name (the `Input:` line of a subject set in outfmt 0 to 4 and 7, SAM's `@PG` command line), and no engine call reads one as a file. A host that inflated a compressed file itself (`inflate_begin`) registers the inflated bytes and puts the file's own name in the argv (`q.fa.gz`, not `q.fa`), so that the outputs equal the command line's on that argv (line WQA's D21, 2026-11-23) |
| `losat_web2_register_inputs(argv_ptr, argv_len, query_ptr, query_len, subject_ptr, subject_len)` | argv as for `run` (with its `-num_threads`); the original FASTA bytes of the query and of the subject | `0`; emits an *inputs* JSON on stream 2 (§9) with a query handle and a subject handle; for an argv of `-db` (*start* `database`) the subject bytes are empty, the JSON's `subject` is null and `run` takes subject handle `0` (the subjects are the database's). It reads both inputs as the CLI reads its `-query` and `-subject` for that argv (line WQA's sweep, D7, 2026-11-19): the inputs of a run whose *start* says `register_inputs` (`-parse_deflines`), whose ids the reader parses from the deflines in one scope for the queries and the subjects, the subjects read first and the query with the argv's `-strand`, and whose ranged ids (`x:10-40`) are refused when the run writes the BLAST archive (outfmt 11). One library function does this for the CLI's reading of its inputs and for this export (`api::input::read_search_inputs`, `losat::web_adapter::read_search_inputs`). Its errors are the CLI's, in the CLI's order: a subject's error (and a subject input without records, `BLAST engine error: Empty CBlastQueryVector`) before the query is read, with the CLI's text (the reader's errors, LOSAT's refusals of ids: two records of one id, `… is not supported by LOSAT's BLASTN`). The query is read as in `register` (line WQA's D10, 2026-11-19): its records before the reader's error, or before LOSAT's refusal of a record's id, are registered with that error, which the run meets in the record's query batch, after the batches before and their errors (two query records of one id in a batch before it, for example); the refusals of the scan of `register` (§9) are the reader's own, in the reader's order. The scan's records are checked as in `register`. Only `run` with the same argv takes the two handles, together: `run` refuses them with any other argv, with another handle, and refuses handles of `register` for an argv whose *start* says `register_inputs` |
| `losat_web2_release(handle)` | handle | `0` |
| `losat_web2_file_begin(kind, name_ptr, name_len, size)` / `losat_web2_file_chunk(file, ptr, len)` / `losat_web2_file_end(file)` / `losat_web2_file_release(file)` | *named input files* (line NP WW1, WW2; kind 2 line RG E1w-3a): kind `2` = a registry file of NCBI's application start, the file of `-conffile` (the start of every later `validate` and `run` reads it, `losat_web2_start`), whose bytes are kept as given (no check at `file_end`: the start reads them as the CLI reads the file, and reports them where the CLI does); kind `1` = a WindowMasker unit counts file, the file that BLASTN's `-window_masker_db` names, or the file of `-window_masker_taxid <taxid>`; the name (UTF-8) is that option's value, or `windowmasker/<taxid>/wmasker.obinary` (the taxid in decimal: LOSAT Web's directory of `-window_masker_taxid` is `windowmasker`, where the CLI reads NCBI's registry entry `[WINDOW_MASKER] WINDOW_MASKER_PATH`, and the engine joins the names as NCBI does), compared byte for byte with the argv's or the joined one (no path is normalized: `./x` and `x` are two names); `size` is the file's length in bytes; then the file's bytes in order, in chunks of any size | `file_begin` returns a file handle ≥ 1, or fails for another kind (`unknown input file kind N`), an empty name, a name that is registered or loading (`an input file named N is registered; release it before registering another`: the host releases the old file first, nothing is replaced silently), or `size` bytes that cannot be had (`not enough memory for the WindowMasker counts N (S bytes)`, `not enough memory for the registry file N (S bytes)`; the bytes are reserved here, before any chunk). A name is one file of one kind; an option reads the files of its kind. `file_chunk` fails, and drops the file, when the chunks pass `size`. `file_end` fails, and drops the file, when fewer than `size` bytes came; otherwise it registers the file and emits a *file* JSON on stream 2 (§9). For kind 1 it checks the bytes as NCBI's `DiscoverStatType` does and, for the two binary formats that BLASTN reads, reads the counts with BLASTN's settings (`s_BuildSeqMasker`), keeps both results and drops the bytes; for kind 2 it checks nothing and keeps the bytes. An error of the contents (a format other than the two binary ones, a file too short for its tables) does not fail `file_end`: it is in the JSON's `error` and `stat_type`, and the run reports it as the CLI reports that file, with the CLI's message, exit status and place among the other errors. A registered file is seen by every later `run` of the instance, as a file of that name in the working directory of the CLI, until `file_release`, which frees its memory. A run whose `-window_masker_db` names no registered file fails as the CLI fails for a file it cannot open (`Error: … could not open N`, exit 255); a run of a taxid without a registered file is the CLI's run of a taxid without a file (each query batch reports `could not open ` with the empty name, exit 3). `validate` reads the two options where the CLI reads them (after `-dust`, before the query range and the other options: the check of the file of `-window_masker_db`, "at most one of" the options, the file of the taxid), so a host that registers the files before it validates gets the CLI's first error; for a `-window_masker_db` name that is not registered `validate` leaves the file to `run` (a host may validate without the file). A run without the options is unchanged |
| `losat_web2_scan_begin(program_ptr, program_len, role, data_loaders)` / `losat_web2_scan_chunk(scanner, ptr, len)` / `losat_web2_scan_end(scanner)` | the program, the role (`0` query, `1` subject) and the data loaders (`1` on, `0` off) of the input, as for `register`: the scan reads with the reader of that input, NCBI BLAST+'s reader with nucleotide flags (BLASTN, TBLASTX, the TBLASTN subject, the BLASTX query) or with protein flags (BLASTP, the TBLASTN query, the BLASTX subject) (version 4 took a parser kind, `1` or `2`, and read with the data loaders on; kind `0`, the `bio::io::fasta` 1.6.0 reader of the four programs before session SF, was removed after the application switched to kinds 1 and 2, plan DW-38); FASTA bytes in chunks of any size | `scan_begin` returns a scanner handle (`unknown program 'N'`, `unknown input role N`, `data_loaders must be 0 or 1, not N`); `scan_end` emits a *scan* JSON on stream 2 (§9), or fails with the reader's refusal, as `register`; the scan of an input whose reader fails ends at the reader's error, with the records before it (`reader_error`, §9) |
| `losat_web2_inflate_begin(name_ptr, name_len)` / `losat_web2_inflate_chunk(inflation, ptr, len)` / `losat_web2_inflate_end(inflation)` | the file's name (UTF-8); the file's bytes in chunks of any size | `inflate_begin` reads the name with the engine's rule, the one by which the command line inflates `-query`, `-subject` and makeblastdb's `-in` (`losat::web_adapter::compression_of_name`: the name ends in `.gz`, `.bz2` or `.zst`, in any case). A name that ends in none of them returns `0`, no handle and no response: the file is read as it is. Otherwise it returns an inflation handle and emits an *inflate* JSON on stream 2 (§9). `inflate_chunk` inflates the next bytes, and `inflate_end` ends the file; both emit the inflated bytes on stream 5 (§5) with the command line's decoder (`PushInflater`, the push form of the reader `Compression::decoder`; swebok-improve SW8 DECISIONS 44), so a file gives the bytes that the command line reads, however it is cut into chunks. A file that is not what its name says, is cut short or is empty fails `inflate_chunk` or `inflate_end` with the command line's sentence without its `-query: ` (`Compression::not_data`: `'x.fa.zst' is named .zst but is not zstd data (a zstd frame's checksum does not match its data); LOSAT reads a file named .zst through zstd`), where the decoder meets the fault, and the host discards what stream 5 gave for the file. The handle is dropped after an error and by `inflate_end`; a host that stops before the end of the file calls `inflate_end` and ignores its result. The application refuses a `.bz2` name itself (line GZ decision 11), after `inflate_begin`. A Wasm linear memory never shrinks, and the decoders grow it (a zstd frame's window is reserved at its first frame, up to 128 MiB), so the application runs each file on an instance of its own, opened for that file from the compiled module and dropped when the inflation ends, fails or is cancelled; the session's shared instance only reads the name rule (`inflate_begin`, then `inflate_end` with nothing pushed) |
| `losat_web2_run(argv_ptr, argv_len, query_handle, subject_handle)` | argv, handles registered for the argv's program; for an argv of `-db` (*start* `database`), subject handle `0`; for an argv whose *start* says `prints`, none: the host passes 0, which `run` does not read | emits the program's supported format streams (0, 6, 7), stream 4 when the argv's `-outfmt` is another format or 6 or 7 with fields (§7), stream 1 (§8) and stream 3; returns after the run ends. It searches once for each distinct hit list size of the formats (§7): one search, or two when `-num_descriptions` / `-num_alignments` give format 0 another size than formats 6 and 7. The registered records enter the engine's `run_local` as the reader's records, so a run reads them as the CLI reads the same file: the subjects' warnings first, then the query batches with their warnings; a first batch that fails (for example a batch of queries without residues only, NCBI's `BLAST engine error: Warning: Sequence contains no data ...`) fails the run after the CLI would have written the outfmt 0 prolog to its file; where outfmt 0 is the argv's own format the run sends it, and the host keeps it as the output of the failed run (§5; version 7). Subject titles are written by the report layer ported in E2h (`NStr::HtmlDecode`, `GuessEncoding`, `x_CleanAndCompress`), so a title that NCBI decodes no longer fails a run (BLASTN, TBLASTX and TBLASTN outfmt 0/6 byte-identical to NCBI 2.17.0, SFc 2026-10-10); the CLI's explicit rejections that happen after the search starts are LOSAT's refusals of one format's report (`docs/evidence/losat_web_e2h/AUTHORITY.md` §J: BLASTP's punctuation-only subject titles in outfmt 0; `sacc`, `saccver` or a first row's `staxid` of a subject whose title NCBI cannot decode, where NCBI stops, PD-LOSAT-NCBI-DEFECTS 1.37; BLASTX's ungapped `-subject` subjects of empty rows that NCBI's mapper rewrites, which the CLI refuses without NCBI's `Error pre-fetching sequence data` warning, with `-subject_loc` or without; with `-sorthits`, a description table whose last total is above 99999 bits in outfmt 0 to 4, which NCBI pads without end), and NCBI's errors of one format's report alone (for example BLASTX's `CSeq_align::GetSeqRange(): Row is empty` in outfmt 0, where its outfmt 6 completes): an error of the argv's own format fails the run, as the CLI fails, and one of a format that the argv does not ask for fails that stream alone while the run writes the others, at a check before the report's rows or in its middle (§7; a stream that failed after the host received a chunk of it is discarded, §5). A run also fails for NCBI's own errors of the search (for example `The average subject length is too short` when every subject is without residues) |
| `losat_web2_tool(argv_ptr, argv_len)` | argv (§7): argv[0] the tool's name, `windowmasker`, then the command line's words after `losat windowmasker` | runs the tool as the command line runs it in its working directory, which is the run's (§3: it opens the files that the words name there, `-in`, `-ustat`, `-fa_list` and the files it lists, `-ids`, `-exclude_ids`, a BLAST database of `-infmt blastdb`, `-conffile`), with the null device as its standard input (`-in -` reads no record): its standard output goes to stream 11 and its standard error to stream 3 as they are written, each file that it writes (`-out FILE` of every mode, `-logfile FILE`) to stream 12 as it is written, in frames of up to 1 MiB (§5), and when the run ends its exit status to stream 2 (*status*, §9); it returns `0` whatever the exit status. The bytes are those of the command line's run of the words in a folder of the run's files (exit status, standard output, standard error, every file it writes; line NP WW3a, 2026-11-22). An output file name is opened as the command line of a user on Linux opens it there: a name in a folder that the run's directory does not have, a directory and a directory of root's do not open (the command line's error or lost counts, as there); `/dev/stdout` is the standard output (stream 11), `/dev/null`, `/dev/zero`, `/dev/random`, `/dev/urandom` and `/dev/stdin` lose the bytes, `/dev/full` fails its writes; LOSAT's refusal (`… is not supported by LOSAT's WINDOWMASKER`, exit status 1 as the command line's refusals) for a name that is the standard error, another device or a file of procfs, and for a file of the run's directory, which the command line would empty. The name of `-logfile` is read as for a search (§4 `losat_web2_start`). Every name that the system opens at one path is one file of the run (`x`, `./x`, `sub/../x`), as the system shares a file that a process opens twice: with `-out x -logfile x` each opening writes at its own place, the output's opening empties the file, and the log's messages go to its end (stream 12's offsets and truncations, §5; line NP WW3a's review, M2, D-WW3-9, 2026-11-23). `-h`, `-help`, `-xmlhelp` (NCBI's XML) and `-version` print on stream 11 what the command line prints. It fails (`-1`) for an argv that is not UTF-8 (`argv is not UTF-8: …`) and for a name that is not `windowmasker`: `<name> is not available in the browser (ABI v2)` for the command line's other programs, `unknown program '<name>'` otherwise, before any stream |
| `losat_web2_last_error_ptr()` / `losat_web2_last_error_len()` | — | last error message |

`validate` has no query, so it does not check the Karlin-Altschul tables of BLASTN (scores and gap costs without a row,
LOSAT's limit of greedy gap costs). A CLI run meets them only when some query context has an ungapped
Karlin-Altschul block, which depends on the query's letters; `run` makes them and stops with the CLI's message and
exit code (`BLAST engine error: Error: …`), or writes NCBI's warning `Could not calculate ungapped Karlin-Altschul
parameters …` and no hit. What `validate` rejects, the CLI rejects.

## After a trap

The release profile sets `panic = "abort"` (`LOSAT/Cargo.toml:141-145`). The adapter has
a copy of it (`web/adapter/Cargo.toml:27-32`, checked by
`web/adapter/tools/check_build_identity.py`), and the adapter builds the shipped
reactors. So a Rust panic inside any export of this ABI or of ABI v1
([`abi_v1.md`](abi_v1.md)) does not unwind and does not return `-1`. It ends as a Wasm
trap, and the host sees an exception (in JavaScript, a `WebAssembly.RuntimeError`). The
same holds for an allocation failure and for a pointer outside the memory (§6).

After a trap:

- Do not call the instance again, not even `losat_web2_dealloc` or
  `losat_web2_last_error_ptr`. Its linear memory is in an unknown state.
- Every handle (`register`, `scan_begin`, `inflate_begin`), the stored records and the last error are
  lost. There is no error message to read for the call that trapped. The chunks of
  streams that the host received before the trap belong to a run that did not finish:
  the host discards them (§5).
- The host discards the instance and instantiates the module again, calls
  `_initialize` (§2), and registers its inputs again. The new handles differ from the
  old ones.
- In the threaded module the main instance and the worker instances share one memory
  that the host created (§3). The host discards the instance, the shared memory and the
  workers together, and creates new ones. The old shared memory is not reused: its
  contents are in an unknown state.

A failure that is reported (`-1` and a message) is not a trap. The instance stays usable
after it.

## 5. Output streams

`emit` may be called many times per stream. Chunks of one stream arrive in order. The
bytes are valid only during the call, so the host must copy them. `emit` is always
called on the thread that called the export: the engine formats on the calling thread,
which runs slot zero of the search's thread pool (`LOSAT/src/utils/threading.rs`). `emit` must
not call an export of the instance, and an exception thrown by `emit` leaves the instance as
after a trap (§4, "After a trap").

| Stream | Content |
|---|---|
| `0` | outfmt 0 text (only for programs that support outfmt 0) |
| `6` | outfmt 6 text |
| `7` | outfmt 7 text (only for programs that support outfmt 7) |
| `4` | the output of the argv's `-outfmt` when it is a format besides 0, 6 and 7 that the CLI writes (5, 8 to 16, 18, 20, BLASTN's 17; with or without words after the number), or 6 or 7 with words after the number (fields, `delim=`, §7): the CLI's output for the same argv, text, or binary ASN.1 for outfmt 9. For 13 and 14 it is a ZIP archive of the files that the CLI writes with `-out results`: the side files `results_<n>.json` (13) or `results_<n>.xml` (14), one per query, each added as its report is written, then the master file that lists them, as `results.json` or `results.xml`, then the central directory. The entries are stored (no compression) with their CRC-32 and the date 1980-01-01 00:00, so that the archive of a run is always the same bytes; ZIP64 records come only where the classic fields cannot hold a value (65,535 entries or more, or an offset of 4 GiB - 1 or more). Only such a run writes it. An output without bytes (outfmt 6 of a search without hits) comes as no chunk at all, so a host tells from the argv, as the adapter does, whether the run has this output |
| `1` | HSP records as JSON Lines, one object per HSP (§8); every program |
| `2` | JSON response of `describe`, `validate` (*validation*), `register`, `losat_web2_register_inputs` (*inputs*), `scan_end`, `inflate_begin`, `losat_web2_start` (*start*), `file_end` (*file*) or `losat_web2_tool` (*status*) |
| `3` | diagnostics that the CLI writes to stderr (warnings), UTF-8, each written once. It starts with the CLI's standard error for the argv (for its `-outfmt`, or outfmt 0 without one; line WQA's sweep, D3 review r2, 2026-11-19), as the CLI writes it: the formatting handler's lines for the argv's format (`-num_descriptions`, `-line_length` and `-sorthits` for formats above 4, `-sorthsps` for formats other than 0, `Examining 5 or more matches` for a hit list size below 5), the reader's and the queries' warnings (which the writers of the structured formats 5, 8, 9 and 11 to 17 do not write) and the warnings of the fields' writer (NCBI's `Taxonomy name lookup from taxid requires installation of taxdb database`). A run that completes adds what the CLI's standard error for the run's other formats adds to it, so that each text comes once, the argv's search first and each search's formats in its order (the tabular ones, then those up to 4). For each other format of the argv's search: the formatting handler's lines (outfmt 6's `-line_length` line where the argv asks for outfmt 0), then, where the argv's format is a structured one, what the format's run writes with its queries' reports: the queries' warnings (outfmt 0, 6 and 7 write a query's `Sequence contains no data`, outfmt 5 does not; D3 review r3 M2, 2026-11-19; BLASTX: every line of its reports' standard error, such as NCBI's `Error pre-fetching sequence data`), those of the queries whose reports it began where it failed alone, then its refusal where it failed alone (a format refused before its reports adds its refusal alone; one whose report stops after a query's warning adds that warning, which the formats after it do not add again; review of 2026-11-20). Then, in a run of two searches (§7), the lines of the other search's warnings that stream 3 does not hold yet (a line as many times as that search writes it more often than stream 3 holds it: the reader's and the queries' warnings, which both searches write, and a handler's line that both write, `Examining 5 or more matches`, come with the argv's search), followed by its error when it failed (its formats fail alone, §7) or, format by format, by the same lines of its formats; an error or a refusal that stream 3 holds (streams 6 and 7 stop at the same row) is not written again. A run that fails also sends it: what the CLI writes to its standard error for the run's argv before the error that ends the run (its warnings, as NCBI's `Subject sequence contains no data` before `The average subject length is too short`), without the error, which is the run's (`losat_web2_last_error`), so that stream 3 and then the error as the CLI writes it are the CLI's standard error for the argv: NCBI's messages as they are (with their line end), but for those of the application start (§4), which come without their last line end and begin with NCBI's severity and error code (`Critical: (106.23) Registry file "x.ini" not found, aborting.`, `Error: (CRegistryException::eSection) ...`; D3 review r2 L2), and LOSAT's own errors after `Error: ` and with a line end. In a run of two searches, whose outcome is that of the search of the argv's own format (§7), these are the warnings of that search alone, not those of the other search, which the CLI writes for formats that the argv does not ask for (line WQA's sweep, D3 review, 2026-11-19). Where the CLI writes its run's messages to the file of `-logfile` (TBLASTX), the file holds what the CLI's file holds for the argv, its warnings and the error that ends the run, and stream 3 does not (line WQA's sweep, D3, 2026-11-19); a run that completes adds to stream 3 what the other formats add, which no file of the argv's CLI run holds (D3 review r2, 2026-11-19: the file had the tabular formats' handler lines and the other search's warnings and error). Where TBLASTX's run messages go nowhere (`-logfile ""`, `-logfile /dev/null`, §4), so does what the other formats add: stream 3 has none of it, as the CLI's standard error for no format has (D3 review r3, 2026-11-19), but for LOSAT's own refusals of the outputs that failed, which the CLI writes to its standard error whatever `-logfile` is (stream 9) |
| `8` | the files of a run besides its outputs (line RG item E1w-3a): the file of `-logfile NAME` (§4, `losat_web2_start`), named `NAME`, as a ZIP archive with stored entries as stream 4's. It is emitted after the run, only when the run has such a file, and also when the run fails (the host keeps it: TBLASTX's file then holds NCBI's error, and a start that ended holds what the CLI writes to the file) |
| `5` | the inflated bytes of `inflate_chunk` and `inflate_end` (§4): each call sends what the decoder released during it, in chunks of at most 1 MiB, in order; the host copies them during the call |
| `9` | a run's signal that one output stream (0, 6, 7 or 4) of this run failed, that the host discards what it received of it, and where its error is: each `emit` on stream 9 is one signal of 8 bytes, that stream's number and then the place of its error, each an unsigned 32-bit little-endian integer. The place is where the run wrote the error, with the values of stream 10's places: `3`, stream 3; `8`, the file of `-logfile` (stream 8); `0xFFFFFFFF`, nowhere; and since version 7 `0xFFFFFFFE`, the error that ended the run, which is the run's own (§4) and on no stream, for an output that a failed run does not keep (§5, end). The error of a format besides the argv's `-outfmt` (outfmt 6 and 7 for an argv of outfmt 0) is the one that the CLI's run of that format writes; the CLI's run of the argv writes none, and its file of `-logfile` holds none (§5, stream 3). So NCBI's message goes where what the run's other formats add goes (stream 3): nowhere in TBLASTX wherever the run's messages go nowhere (stream 10: `-logfile ""`, `-logfile /dev/null` and every other name whose messages go nowhere, §4 *start*, such as `/dev/fd/2` and `/dev/full`), else stream 3, also in TBLASTX with `-logfile NAME`, whose file holds what the CLI's file holds for the argv (line WQA's sweep, D3 review r2, 2026-11-19), so a run gives no place `8`. LOSAT's own refusal (`Error: ` and the reason, which NCBI has no message for) is on stream 3 whatever `-logfile` is, as the CLI writes LOSAT's own errors to its standard error. The adapter takes each place from the writer of the error (LOSAT's `CallerSink::error`), so the place is where the error is. It comes for every stream that fails alone (§7), whether or not a chunk of it was sent, and for each stream of a search that fails in a run of two searches (§7); no chunk of that stream follows, the stream of the run is empty, and the run goes on with the other streams, except for the place `0xFFFFFFFE`, where the run has failed and goes on with nothing. A host tells from these signals alone which outputs a run could not write and where their errors are (an empty stream is not such a sign: outfmt 6 of a search without hits is empty). Since ABI version 3 after a chunk was sent (4 bytes, the stream), since version 4 for every such stream, with the place of its error (§10) |
| `10` | a run's signal of where its messages go: the warnings that the CLI writes to its standard error without `-logfile`. One `emit` of 4 bytes, the place as an unsigned 32-bit little-endian integer, with these values: `3`, stream 3; `8`, the file of `-logfile` (stream 8), where TBLASTX with `-logfile NAME` writes them, as the CLI does; `0xFFFFFFFF`, nowhere, wherever TBLASTX's `-logfile` sends them nowhere, as on the CLI: `""`, `/dev/null` and every other name whose messages go nowhere (§4 *start*: a file in a file system without free space such as `/dev/fd/2`, a device whose writes are lost such as `/dev/full`). The adapter takes it from LOSAT's start of the run (`CallerStart::messages_place`). It comes once in every run whose start does not end it, a dry run's included, before the run's outputs. The warnings that the CLI's run of the argv writes go there (§5, stream 3): in TBLASTX with `-logfile NAME` the file holds what the CLI's file holds for the argv, those of its own `-outfmt`; what the run's other formats add (the formatting handler's lines of the other formats, in a run of two searches, §7, the other search's warnings, and NCBI's errors of the outputs that failed, stream 9) goes to stream 3 there, and nowhere where the place is nowhere (line WQA's sweep, D3 review r2 and r3, 2026-11-19). Where the place is not `3`, stream 3 holds only the start's messages that NCBI writes to `cerr` at once (the report of a registry file that NCBI cannot read, `Error: (110.4) Error reading the registry ...`, which the file of `-logfile NAME` holds too), what the run's other formats add where the place is the file, and LOSAT's own refusals of the outputs that failed (stream 9). An error that ends the run is the error of `run`; in TBLASTX with `-logfile NAME` the file also holds NCBI's error that ends the search (`-evalue 0`), and in every program the file holds the posts of NCBI's argument reader for an argument error (`-evalue abc`), whose run error is then the CLI's standard error with NCBI's line that points to the error logs (§4), as on the CLI. A host tells from it whether an empty stream 3 is a run without warnings (place `3`) or one whose warnings are elsewhere, where they may be none (an empty file). Since ABI version 4 (line WQA's sweep, D7, 2026-11-20; §10) |
| `11` | the text that the CLI prints on its standard output in place of a run (line WQA's D14, 2026-11-20): for a help word the CLI's help (`-xmlhelp`: its XML), for a version word its version text, the bytes that the CLI of the same argv writes, with exit status 0. It comes in a run whose *start* says `prints`, after stream 10; such a run reads no input and emits no output format, HSP record or signal of stream 9. One producer gives the text to the CLI and to the adapter: NCBI's argument reader, whose help ends `losat::web_adapter::try_parse_cli` with `ArgvError::Printed` (the XML of `-xmlhelp` from `cli::xml_help_text`), and the start's version (`CallerStart::version`). Since ABI version 10 (§10) |
| `12` | the files that a tool run writes (`losat_web2_tool`; version 13; line NP WW3a, D-WW3-6 and D-WW3-9, 2026-11-22 and 23), each as it is written, in frames of up to 1 MiB: writes that follow each other in a file are joined into one frame until it is full, a write elsewhere in the file or a truncation comes, or the run ends, so a file of small writes (the messages of `-logfile`) may come only when the run ends. Every `emit` on stream 12 is one frame, 8 bytes and a payload: the file's number in the run (from 0, in the order the run opens its files; the file of `-logfile` opens first) and the frame's kind, each an unsigned 32-bit little-endian integer. Kind `0` is the file's beginning: the payload is its name as the argv writes it where the run opens it first, UTF-8, and the file is empty. Kind `1` is bytes of the file: the payload is the offset where they go in the file, an unsigned 64-bit little-endian integer, then the bytes (at most 1 MiB − 16 a frame); the host writes them there, over what the file has there, and past its end with zero bytes between. Kind `3` is a truncation: the payload is the file's new length, an unsigned 64-bit little-endian integer (shorter: the rest is cut; longer: zero bytes are added). Kind `2` is the file's end (no payload). A file's bytes and truncations come in the order the run makes them, which the host applies in that order: most bytes come at the end of what came before, and a run that opens one file twice (`-out x -logfile x`, §4) writes through each opening at its own offset and empties the file when its output opens it (a truncation to 0). The ends come when the run ends, after its last write, every file's in the order of their numbers, before *status*; the frames of two files may alternate. A run that fails ends its files too: the host keeps them as the command line's run leaves its files. A run that stops the instance (a trap, memory exhausted) sends neither the ends nor *status* (§4). A tool run emits no stream 8; a search emits no stream 12 |

The adapter sends each output stream in chunks of 1 MiB as they fill, and the rest when
the run ends, so an output never has to fit in linear memory as a whole. Stream 2 is
chunked the same way: a large *scan* or *register* response (thousands of records)
arrives in several chunks, which the host joins before it parses the JSON (S09). Stream 3 and
the HSP records are sent after the search. If `run` fails, the host keeps streams 3 and 8 with the failed
run, and the output of the argv's own `-outfmt` (stream 4 for an `-outfmt` of its own, else the stream of
the argv's format, 0 when it gives none) as the CLI's run of the argv leaves it before its error: the
bytes that the output's writer held are sent before the run fails, so the host holds exactly the bytes
of the CLI's file (`-out` or standard output), also those of its last, incomplete chunk (line WQA's
sweep, D15, 2026-11-20; version 7). The adapter drops every other output of the run with the signal of
stream 9 (the place `0xFFFFFFFE`, below), which it gives for each of them, the outputs of a search that
did not run too: they are not what the CLI's run of the argv writes, and a format of the search that
failed stopped where that search did. A run that fails sends no HSP records (stream 1). The archive of
outfmt 13 or 14 of a failed run holds the side files written before the error and the master file as the CLI
leaves it (empty, where the run failed before it was written). A failure of the adapter itself (an error in the output ranges) keeps no output. The
host discards what it received if the run traps (§4, "After a trap") or is stopped. A stream that fails alone
(§7), and each stream of a search that fails in a run of two searches, is followed by its signal
on stream 9: the host discards what it received of that stream (a chunk of it may have been sent
before it failed) and keeps the others, and the signal says where the stream's error is (stream
3, or nowhere). The signal of stream 10, where the run's messages go, comes before the outputs,
once in every run whose start did not end it.

## 6. Memory ownership

The host allocates input buffers with `losat_web2_alloc`, writes the bytes, passes the
pointer, and frees the buffer with `losat_web2_dealloc` after the call returns. The
adapter never keeps a pointer to a host buffer after a call: `register` copies what it
keeps, and so does `file_chunk` (one buffer can serve every chunk of a file).

A named input file holds its size in linear memory from `file_begin` to `file_end` (the bytes),
then about the size again while its counts are read (the tables), and keeps the tables until
`file_release` (the bytes are freed at `file_end`). The memory of a Wasm instance does not
shrink: the freed bytes are reused by the searches. A host that keeps counts files across runs
keeps at most what its memory allows; the application keeps one per instance (plan §4.6 R1).

Pointers and lengths are unsigned 32-bit values. JavaScript receives a Wasm `i32` (a
result of an export, an argument of an import) as a signed number, so a pointer of 2 GiB
or more arrives negative. A JavaScript host converts every pointer and length to its
unsigned value (`ptr >>> 0`) before it uses it as an offset into the memory: those of
this ABI (`losat_web2_alloc`, `losat_web2_last_error_ptr`, `emit`) and the arguments of
the WASI imports.

## 7. argv

- UTF-8 words separated by NUL (`\0`). The first word is the program name
  (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`).
- `-query <name>` and `-subject <name>` are required. The names are used only in
  output text (for example the `# Database:` line of outfmt 7); the sequences come from
  the handles.
- `-out` is rejected: the adapter writes every format the program supports (`describe` lists
  them).
- `-outfmt` (once, as on the command line; `-outfmt=value` too) is read as the CLI reads it, with
  the CLI's messages: a value that is not a format (`'abc' is not a valid output format`,
  `Formatting choice is out of range`), a format that LOSAT's program does not write, NCBI's
  formatting handler's errors (another program's format such as SAM, a custom delimiter that a
  field uses to separate its entries). Outfmt 13 and 14 do not need the CLI's `-out`: the run gives
  the engine a receiver of their side files (`OutputSink::Files`), and stream 4 is the archive of
  their files (§5). The value 0, 6 or 7 without
  words after the number adds nothing: the run writes those formats anyway; nor does 0 with words
  (`0 qseqid`), which the CLI reads as 0. Every other format that the CLI writes, with or without
  words, and 6 and 7 with words after the number (after NCBI's white space around the value is
  removed, a space follows the number: fields, `-field`, `std`, `delim=`) add stream 4, written by
  the search of its hit list size as its first format, so that it is the CLI's output for the same
  argv and adds no search: the tabular search (the first search, below) for every format above 4,
  whose hit list size is that of 6 and 7, and the format 0 search for 1 to 4, whose size is that of
  0. What the CLI's run of that format has besides, the run has too (session W11): SAM's `@PG`
  line shows the command line (`CL:losat ` and the run's argv with the host's `-num_threads`, each
  word followed by a space, as NCBI's `GetCmdlineArgs`; the CLI writes its own argv there); the BLAST archive's Bioseqs hold
  the whole records (with `-query_loc` and `-subject_loc` too) and the deflines as read (`register`
  keeps a defline that `>` and the title do not give). A custom delimiter that a field uses to separate its entries
  (`6 delim=; qseqid staxids`) is NCBI's formatting handler's error, which the CLI makes after the
  ranges and the query options: `validate` reports the error of the CLI's own checks of the argv
  with that `-outfmt` (a `-subject_loc` of bad syntax first), and `run` fails with it: the tabular
  search writes the argv's own format, whose search decides the run (below).
  TBLASTN's and TBLASTX's engine reads every format of a run before the ranges and the option
  checks that come before the handler (`-seg`, `-query_loc`), so their `run` of an argv with both
  errors reports the delimiter's, where the CLI and `validate` report the other.
  The formats that the handler rejects (SAM of another program, 19, 21) are rejected with the
  handler's message as the argv is parsed, before the ranges. The rule reads the words, not the fields that they give, so that a
  host reads it from the argv without the table of fields: `6 std` adds a stream 4 with the bytes of
  stream 6. The engine is given the run's formats, and the checks of the options read them as
  without `-outfmt` but for the argv's own format, which comes first in its search: streams 0, 6, 7
  and 1 are those of the argv without `-outfmt` (sessions W10, W11).
- A format other than the argv's own (stream 0, 6 or 7 where the argv's `-outfmt` is another) whose
  report fails, as the CLI's run of that format alone fails (the refusals and NCBI's errors of one
  format of the `run` row, §4), fails alone in the search that writes the others
  (`FormatRefusal::FailsAlone`), whether the failure comes at a check before its rows or in the
  middle of its report: its stream is empty, its error follows the search's warnings on stream 3
  (NCBI's message goes nowhere where the run's messages go nowhere, LOSAT's own refusal comes on
  stream 3 whatever `-logfile` is; §5, stream 9; each text once: streams 6 and 7 stop at the same
  row), and the other formats are written;
  without stream 6 there are no HSP records (stream 1 points into it, §8), and without stream 0
  the records have no `out0` and `out0_subject`. The bytes that it wrote stay in the adapter,
  and the run's signal on stream 9 tells the host that the stream failed, where its error is, and to
  discard what it received of it (§5), where a chunk of them was sent (a report that failed after 1 MiB of it, or,
  in BLASTX, which writes its formats query batch by query batch, a batch of 10,002 query letters,
  NCBI's, a refusal at a later batch) and where none was. A subject whose title NCBI's `GenerateDefline` cannot decode
  stops the default tabular formats (streams 6 and 7 fail alone where the argv asks for outfmt 0
  or another format; outfmt 0 writes NCBI's line in place of the alignment), and a BLASTP subject
  titled `, ,` fails stream 0 alone where the argv asks for `-outfmt 6` (whose row NCBI writes).
  So does a description table of `-sorthits` whose last total is above 99999 bits where the argv
  asks for `-outfmt 6`, 5 or another format (stream 0 fails alone; the CLI writes NCBI's outfmt 6
  and 5), an NCBI error of outfmt 0 alone, such as BLASTX's `Row is empty`, and LOSAT's refusal of
  outfmt 0 of a record that is both a query and a subject under `-parse_deflines` and has masks
  (BLASTN, BLASTP and TBLASTX: stream 0 fails alone where the argv asks for `-outfmt 6` or another
  format, whose output the CLI writes; line WQA's sweep, D7, 2026-11-19). An error of the
  argv's own format fails the run, as the CLI fails (`-outfmt 6` or `7` with such a subject: the
  CLI's exit status 1; outfmt 0 with `Row is empty`: NCBI's exception; outfmt 0, or no `-outfmt`,
  with such a record: the CLI's exit status 1).
- `-num_descriptions` and `-num_alignments` are accepted (every program; `-max_target_seqs` excludes both, as on
  the command line) and described (§9). NCBI derives the hit list size of a format from them
  (`blast_args.cpp:2894-2973`): format 0 keeps the larger of `-num_descriptions` (500) and `-num_alignments` (250) where
  either is given, else `-max_target_seqs`; formats 6 and 7 keep `-max_target_seqs`, else `-num_alignments`, else 500.
  The sizes are equal without the options, with `-max_target_seqs`, and for example with `-num_descriptions 2
  -num_alignments 4`; they differ for `-num_alignments 10` and `-num_descriptions 10`.
- `-num_threads` is appended by the host. When the host falls back to the serial
  module, it passes `-num_threads 1` and records the reason (plan §4.7). The host
  validates the argv that it will run, with its `-num_threads`, so a `-num_threads` in
  the user's words is reported by `validate` (the parser rejects the repeated option).
- Every other word is read by the engine's CLI argument reader, so options and errors behave
  exactly as on the command line. LOSAT's own debug output (environment variables such
  as `LOSAT_TIMING`) goes to the module's WASI standard error, not to stream 3; stream 3
  carries the warnings that the CLI writes.
- Options are resolved separately for each output format, as the CLI resolves its single `-outfmt`. `run` makes
  one search for each distinct hit list size of the formats 0, 6 and 7 (the engine's `run_local` writes the formats of
  one hit list size from one search, and refuses formats whose search options differ otherwise): the tabular
  formats first, then format 0 if its size differs. Without `-num_descriptions` and `-num_alignments` this is the one
  search of today. Each stream is the CLI's output for that format alone with the same options, because the
  hit list size takes part in the search (a smaller list is not a prefix of a larger one: the preliminary list evicts
  subjects before their alignments are extended). With two searches the run costs about twice the time of one. Its
  memory peak is that of the larger search plus the HSP list of the tabular formats, which stays in memory (about the size of
  stream 1) while the second search runs (measured: +21.7 MiB on 142.6 MiB for BLASTN of Sakai against EDL933 with
  `-num_alignments 3`).
- The outcome of a two-search run is that of the search of the argv's own format (the tabular search for `-outfmt 6`
  and `7` and for the formats above 4 of stream 4; the format 0 search without `-outfmt`, for `-outfmt 0` and for 1
  to 4), as the CLI runs that format alone (line WQA's sweep, D3 review, 2026-11-19): the run completes when that
  search completes, and fails when it fails, with its error and, on stream 3, the warnings that it wrote before the
  error (§5). The other search writes only formats that the argv does not ask for (the run writes formats 0, 6 and 7 whatever the
  argv): when it fails, its formats fail alone, as the CLI's run of each of them alone fails. Their streams stay
  empty (a chunk that the host received of them is discarded, stream 9); if it is the tabular search there are no
  HSP records (stream 1), else the records have no `out0`. This is `-num_alignments 0` for the tabular formats
  (NCBI's error `BLAST query/options error: No hits are being saved` for `-outfmt 6` and `-outfmt 7`, while
  `-outfmt 0` keeps its size, 500 or `-num_descriptions`): an argv without `-outfmt` completes with format 0, and
  one with `-outfmt 6` fails with that error after the formatting handler's warning, as the CLI. (A
  `-num_descriptions` of 2,147,483,598 or more, whose preliminary hit list size overflows NCBI's 32-bit `int` in
  BLASTN's format 0, runs: LOSAT computes the size without overflow, PD-LOSAT-NCBI-DEFECTS.) A fault of the adapter
  itself (the output ranges of the records) fails the run in either search. Stream 3 of a run that completes is the
  CLI's standard error for the argv (the argv's search's warnings) followed by the lines of the CLI's standard error
  for the other formats that it does not hold yet (the formatting handler's lines of the other formats, the queries'
  warnings that a structured format does not write, the other search's warnings, and its error when it failed; §5):
  the lines that the CLI user sees running the argv and then the other formats, each text once (line WQA's sweep, D3
  review r2, 2026-11-19: the tabular search's lines came first). When both hit lists are empty (`-num_descriptions 0
  -num_alignments 0`, one search), the run fails as the CLI fails every format. `validate` makes the checks of the
  argv's own search and reports their error: BLASTN's and BLASTX's with the hit list of its formats, BLASTP's,
  TBLASTN's and TBLASTX's with the argv's own format, whose hit list the CLI checks before the e-value (`No hits are
  being saved` for `-num_alignments 0 -evalue 0 -outfmt 6`; D3 review r2, 2026-11-19: they checked format 0's).

## 8. HSP record

One JSON object per HSP of the final, sorted result. Field names follow the Rust
`common::Hit` and `report::PairwiseHit` fields. Every program emits them (TBLASTX since
session S08, BLASTX since SX). BLASTN rows are nucleotide rows: the query on its
plus strand, residues of masked regions in lowercase, no frames. TBLASTX rows are the
translated rows of outfmt 0 (both sequences translated again from the nucleotides with
`-query_gencode` and `-db_gencode`, SEG-masked query residues in lowercase), with both
frames. BLASTX rows are the rows of outfmt 0: the query translated in its frame with
`-query_gencode` (`query_frame` 1 to 3 or -1 to -3; residues of the query's masked
regions in lowercase) and the protein subject (`subject_frame` null); with `-query_loc` and
`-subject_loc` the coordinates and `query_frame` are those of the whole records, as outfmt 6
prints them.

| Field | Type | Meaning |
|---|---|---|
| `index` | integer | 0-based position of the HSP in the final, sorted hit list of the formats 6 and 7 (`HspIndex` in `LOSAT/src/api/local_blast.rs`); the same HSP has the same index in every output format of the same search |
| `q_idx`, `s_idx` | integer | 0-based record index of the query and the subject in the registered inputs (the `index` of the *register* response: records without residues count; with `-query_loc`, a query record that NCBI skips because its interval starts past its end keeps its place, so `q_idx` is not the position among the searched queries). In a run of `-db`, `s_idx` is the subject's OID in the database (its 0-based place in the database's volumes), and its IDs are the database's: no registered record, and no entry of the *inputs* response, stands for it |
| `rank` | integer | 0-based position of the HSP among the HSPs of its query (derived from `index`) |
| `raw_score`, `bit_score`, `e_value` | number | engine values, not rounded |
| `q_start`, `q_end`, `s_start`, `s_end` | integer | 1-based coordinates as printed in outfmt 6 (start > end means minus strand; a BLASTN HSP of one letter has start = end on either strand, and only its `Strand=` line in the stream 0 section shows the strand, S07+ AUTHORITY §P) |
| `query_frame`, `subject_frame` | integer or null | translation frames where applicable |
| `subject_length` | integer or null | subject length |
| `query_aligned`, `subject_aligned` | string or null | aligned sequences with `-` for gaps |
| `out6` | [integer, integer] or null | byte range [start, end) of this HSP's row in the stream 6 text; null if outfmt 6 does not show it |
| `out0` | [integer, integer] or null | byte range [start, end) of this HSP's section in the stream 0 text: its score lines and alignment. The subject heading that precedes the first HSP of each subject is not part of any section; null if outfmt 0 does not show the HSP (for example BLASTX shows alignments for the first 250 subjects by default). In a run of two searches (§7) outfmt 0 has its own hit list: the range is set only where outfmt 0 shows exactly the HSPs that the records hold for that query and subject (the same count and every field of each HSP equal, in order; the subject has a heading in stream 0), else null (for example when outfmt 0 shows only some of the pair's HSPs, or other HSPs besides), and stream 0 may hold sections and headings that no record points to |
| `out0_subject` | [integer, integer] or null | byte range [start, end) of the heading of this HSP's subject in the stream 0 text: the defline and the `Length=` line that NCBI's `x_ShowAlnvecInfo` writes before the first HSP of each subject (`showalign.cpp:3613-3632`). Every HSP of the subject in that query has the same range |

The ranges come from formatter observer events: a formatter reports when it starts and
ends the row or section of an HSP, identified by its `index`, and the heading of a
subject, identified by the index of its first HSP; the adapter records the byte
positions at those moments. The formatter's output bytes do not
change (plan TD-3, §4.5). Displayed numbers are taken from the outfmt 6 row, not
formatted from the raw values (plan §4.4).

## 9. JSON responses (stream 2)

- *describe*: `{ "program", "formats": [0, 6, 7], "parameters": [{ "flag", "help", "takes_value", "default"?, "choices"? }], "query_gencodes"?, "subject_gencodes"? }`, generated from NCBI's argument tables (made from each program's NCBI `-xmlhelp` XML), which the CLI's reader reads: the options of the program's XML but those that LOSAT refuses by name, in the order of the usage line that the program's `-help` prints first, the order in which NCBI's program describes its options (blastn: `-task`, `-db`, `-dbsize`, …), not that of the list of options under it, which `-help` prints by group (blastn: `-query`, `-query_loc`, `-strand`, `-task`, …) (since line RG item K4i-3a, 2026-11-24; the order of the clap definitions before), each with `help` (NCBI's comment of the option, the text that `-help` prints for it, in LOSAT's words; since K4i-3a, the clap definitions' help before), `takes_value` (false for a flag), `choices` (since K4i-2b; a flag has none; `-template_type` and `-template_length` have theirs) and `default` (NCBI's `<default>`, since K4i-2b: none for `-evalue` of BLASTX and TBLASTX, whose help says `Default = 10`), the rows of the library API's `describe` (`losat::web_adapter::described_parameters`), without `-out`, `-outfmt` and `-num_threads`, which the adapter and the host set or refuse, or read from the argv's words as the CLI does (`-outfmt`, §7; `-num_descriptions` and `-num_alignments` are listed, with NCBI's help text, which states the defaults 500 and 250, and no `default`: NCBI tells them given from not given, and `defaults` gives the values that the engine resolves), the program's supported output formats and the genetic codes that the engine's parser accepts for `-query_gencode` and `-db_gencode`, and `"defaults"`.
  - `"defaults": { "tasks": { "<task>": { "<flag>": "<value>" } }, "by": { "<flag>": { "<value>": { "<task>": { "<flag>": "<value>" } } } }, "unless": { "<flag>": ["<flag>"] }, "given": ["<flag>"] }` are the values that the engine resolves for the options a search leaves out (`web/adapter/src/defaults.rs`), from the engine's own resolution of an argv (BLASTP's, TBLASTN's and BLASTX's `resolve`, BLASTN's `coordination` functions, TBLASTX's getters), printed as the CLI's help prints them (`"10"`, `"19.3"`, `"20 64 1"`, `"no"`, `"true"`). `tasks`: each task's values with no other option given (the key `""` for TBLASTX, which has no `-task`). `by`: for `-matrix` and `-word_size`, for each value that the engine resolves (the nine matrices of `ScoringMatrix::ALL`, word sizes 1-32), the values by task of the options that some value changes (a matrix's gap costs, threshold and window; the threshold of a given word size); a program lists only the options that change something (BLASTN neither, BLASTX the word size). `unless`: an option's default holds only while none of these options is given (`-max_target_seqs` excludes `-num_descriptions` and `-num_alignments`, §7). `given`: the options that BLAST+ tells given from not given, with the tasks where it does, so that a host writes one that is set even at its default: the hit list options in every task (§7), and each option that, given at the value shown for it (the task's, or the table's value for a given `-matrix` or `-word_size`), resolves another search (the word size of blastp-fast, tblastn-fast and blastx-fast, whose given 5 sets NCBI's threshold 19.3 where the task's own 5 keeps 20). The response of each program is computed once (some hundred probe argvs) and kept. An option without an entry has the `default` of `parameters` (NCBI's argument XML), or no default.
- *validation*: `{ "data_loaders" }` of `losat_web2_validate` (§4): `true` when the run reads its inputs with NCBI's data loaders on (the CLI's default), `false` when the `DATA_LOADERS` entry of the `[BLAST]` section of a registry file that the run's start reads (the file of `-conffile`) turns them off, by the CLI's rule (`ApplicationSettings::data_loaders`, NCBI's `SDataLoaderConfig`): a value that names neither `blastdb` nor `genbank`, or that names `none` (case folded; `DATA_LOADERS=none`, `DATA_LOADERS=foo` and `DATA_LOADERS=blastdb none` turn them off, `DATA_LOADERS=BlastDB` and `DATA_LOADERS=genbank,foo` keep them on). An empty value in a registry file is no entry, and keeps them on (NCBI's registry; the CLI's probe of 2026-11-20). The CLI reads the setting from `.ncbirc` too, which the browser has not; line WQA's D16, 2026-11-20. A host gives it to `register` and `scan_begin` for the run's inputs. It comes before the option checks (§4), so a host has it also while another option refuses the argv.
- *register*: `{ "handle", "records": [{ "index", "id", "length" }], "reader_error"? }`, from the records of the program's own reader (§4). `reader_error`, a query's only: the reader's error on the record after the records, with the CLI's text, which a run ends with in that record's query batch (§4). `id` is the record's title up to its first space (0x20), empty when the record has no title (never NCBI's `Query_N` or `Subject_N`); bytes that are not UTF-8 become U+FFFD in the JSON only (the output streams keep the title's bytes). `length` is the number of residues the reader stores. `register` also scans the same bytes with the reader of the input (`scan_begin`) and fails if the records, IDs, lengths or residue counts differ.
- *scan*: `{ "records": [{ "index", "id", "header_offset", "sequence_offset", "end_offset", "length", "residue_counts" }], "reader_error"?: { "offset", "defline", "message" } }`. Offsets are byte offsets in the scanned input and may exceed 2³²; they are encoded as JSON numbers and stay below 2⁵³. `residue_counts` maps each byte value of the stored residues to its count: the bytes 0x21-0x7E as the character itself, every other byte (including the space) as `"0xNN"`.
  - The scan (session SF; line WQA's D13 and D16, 2026-11-20) gives exactly the records that the engine's reader of the input (`scan_begin`: the program, the role and the data loaders) returns from the bytes read as a file: records without residues and a first record without a defline included. The scan writes no messages. An input of white space, blank and comment lines only has no record and scans to `{"records":[]}`.
    - `id`: the title up to its first space (0x20), as in *register*. The title is the defline after `>` (after `>?_`, which NCBI reads as `>`) and the white space that follows `>`, up to the first byte below 0x20 after its first byte, without trailing white space; empty without a title; bytes that are not UTF-8 become U+FFFD.
    - `header_offset`: the defline's `>`. A first record without a defline has `header_offset` = `sequence_offset` = 0, so the extracted bytes of the record start the input as before.
    - `sequence_offset`: the first byte that NCBI's line reader reads after the defline's end of line (CR, LF or CR LF; a CR or LF that the line reader drops right after it is skipped), or the end of the input. `end_offset`: the next record's `>`, or the end of the input (a tail that NCBI's line reader loses is inside the last record's range). A `>?` gap line is inside its record's range.
    - `length` and `residue_counts`: the stored residues, upper-cased. With nucleotide flags the reader stores `A B C D G H K M N R S T U V W Y` in either case and counts `U` as `T` (keys `A B C D G H K M N R S T V W Y`); with protein flags it stores every ASCII letter and `*` (keys `A`-`Z`, `*`). A `>?` gap line adds the letters of its gap, counted by the reader's own parser of gap lines (`N` with nucleotide flags, `X` with protein flags; `>?` alone is one letter, `>?unk100` 100). Nothing else is counted: no `-`, digits, `*` with nucleotide flags, invalid letters, white space or other bytes. Lines that NCBI's line reader joins (a CR line end in a file of LF line ends, or an LF in a file of CR line ends) are counted as it reads them.
    - The reader's error `CFastaReader: Near line N, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.` (NCBI's text without the CLI's `BLAST query error: `; N is NCBI's line number: every line counts, CR LF is one line end) ends the scan, as it ends the reader (line WQA's D10, 2026-11-19): the records are those before it, and `reader_error` gives its `message`, its `offset`, the `header_offset` of the record that the reader fails on (0 for a first record without a defline), and `defline`, whether that record starts with a defline (`false` for the input's first record read without one: a byte-order mark, a first line of other text); the bytes from `offset` on belong to no record, and a host that extracts records passes them on after the records so that the run's reader meets the error. A host that joins inputs passes such bytes without a defline after another input's records only to have them read as that input's last record: it refuses that join (review M2, 2026-11-20). `scan_end` fails, in the reader's order, with the refusals of the engine's reader, in its words, which the CLI makes too: with the data loaders on, a first line that NCBI's data loaders would fetch as a Seq-id (the CLI's `the first line of the query ("…") is not a defline and may be a sequence identifier … which is not supported by LOSAT's BLASTP (…)`, the line consumed as NCBI consumes it), and a record over 2147483647 letters (`a query record longer than 2147483647 letters is not supported by LOSAT's BLASTN`; the scan refuses it as it passes that length, since such a record does not fit in the reader's memory in the browser). With the data loaders off, the first line is data, as the CLI reads it.
    - Version 4 also gave each record's `line_layout` (where each residue is in the bytes), which no application read, and refused a `>?` gap line, the lines that NCBI's line reader joins, and the reader's two refusals with LOSAT Web's own `not supported in the browser` (line WQA's D13, 2026-11-20).
  - The scan exists only to index original records for extraction; the search always reads the input with the program's own reader (plan TD-8). `web/adapter/tests/scan_ncbi_properties.rs` checks the scan of nucleotide and protein input, with the data loaders on and off, against the engine's reader and an independent reference walk, in chunks of any size.
- *status*: `{ "status" }` of `losat_web2_tool` (§4): the tool's exit status, the command line's (an integer; since version 13). It comes last, after the run's streams 3, 11 and 12.
- *start*: `{ "conffile", "log_file", "dry_run", "register_inputs", "prints", "database" }` of `losat_web2_start` (§4): `conffile` and `log_file` are the names (strings) or null, `dry_run`, `register_inputs` and `database` Booleans, `prints` `"help"`, `"version"` or null (since version 10). `database` is true for an argv of `-db`: its run searches the database that the engine opens in the run's working directory (§3), takes no subject (subject handle `0` in `run`; `register_inputs` reads the query alone) and refuses a subject handle (line WQA D12a, 2026-11-20; since version 12).
- *inputs*: `{ "query": { "handle", "records": [{ "index", "id", "length", "shown_id" }] }, "subject": { … } }` of `losat_web2_register_inputs` (§4), whose `subject` is null for an argv of `-db`, which reads no subject input (the subject bytes must then be empty; the IDs that the run's outputs show for its subjects are the database's, by the OID that a record's `s_idx` gives, §8): the two handles and the records as in *register* (`id` is the scan's title word, which the host compares with its record table), and `shown_id`, the ID that the run's outputs show for the record: `qaccver` of a query and `saccver` of a subject, the first two columns of outfmt 6 and 7 (`AB000001.1` of `>gb|AB000001.1| …`), lossy UTF-8. A query whose reading ends with an error has `reader_error` as in *register*; the records that the scan found from the record of that error on (LOSAT's refusal of an id, which the scan does not read) are listed with the scan's length and a null `shown_id`.
- *file* of kind 2: `{ "handle", "name", "kind": "registry", "size", "sha256" }`.
- *file* of kind 1: `{ "handle", "name", "kind": "windowmasker_counts", "size", "sha256", "stat_type", "unit_size", "error" }` (`file_end`, §4). `sha256` is the lowercase hex SHA-256 of the bytes given, counted by the adapter over the chunks. `stat_type` is what NCBI's format check found: `"obinary"`, `"binary"`, `"ascii"`, `"oascii"` or `"unknown"`, or null when the check failed. `unit_size` is the counts' unit size when they were read, else null. `error` is null, or the message of the check's or the reader's error without NCBI's framing (for example `not enough data to fill the hash table`); a `stat_type` other than `"obinary"` and `"binary"` with a null `error` is a format that BLASTN refuses (`Only optimized binary windowmasker stat files are supported`). The run gives these errors with the CLI's framing (§4).
- *inflate*: `{ "compression", "plain_name" }` of `inflate_begin`: `compression` is `"gzip"`, `"bzip2"` or `"zstd"`, the format that the name's suffix announces; `plain_name` is the name without the suffix (`a.fa` for `a.fa.GZ`), or `input` when nothing is left (`.gz`), the name under which the application keeps the inflated file. `web/adapter/tests/inflate.rs` checks the three formats whole, joined, damaged, cut and empty, in chunks of one byte and whole, against the command line's reader and sentences.

## 10. Versioning

`losat_web2_abi_version` returns `13`. Any incompatible change to a signature, a stream,
or a JSON field increments the version, and the application refuses a module whose
version it does not know. The name of this ABI (v2, beside ABI v1 of gbdraw in the same module)
and the prefix `losat_web2_` of its exports stay.

Stream 4 and the `-outfmt` of §7 (session W10, 2026-11-05) kept version 2: a host that does not
know them never sends `-outfmt` (an earlier adapter refused it) and drops a stream it does not
know, and a host that sends `-outfmt` to an earlier adapter gets its refusal
(`-outfmt is not accepted`), not other outputs.

The named input files of §4 (line NP WW1, 2026-11-12) kept version 2: they are four exports more,
with no change to a signature, a stream or a JSON field of the others. The application ships with
its module and checks that the exports exist when it instantiates the module; a host that does
not know them never registers a file, and a run of `-window_masker_db` then fails as the CLI's
does for a file it cannot open.

`-window_masker_taxid` (line NP WW2, 2026-11-12) kept version 2: it is a name of the kind 1 files of
WW1 and adds no export, stream or field. A host that does not know it never registers a file of that
name, and a run of a taxid then reports the CLI's taxid without a file.

The inflation exports, the *inflate* response and stream 5 (swebok-improve SW8 D2, 2026-11-11) kept
version 2 as stream 4 did: they add calls and a stream that a host which does not know them never
calls and never receives, and change no other export or stream. The application requires them in the
modules it loads (`ABI_EXPORTS`, `web/packages/losat/src/infra/reactor/artifact.ts`).

The application start of §4 (line RG item E1w-3a, 2026-11-14) kept version 2: `losat_web2_start`, the
kind 2 files, the *start* response and stream 8 are calls, a name kind and a stream that a host which
does not know them never makes or receives. An argv without `-logfile`, `-conffile` and `-dryrun`
gives the same streams as before (the start of such a run in LOSAT Web's environment writes nothing
and changes no setting), and an argv with them, which an earlier adapter refused with the parser's
message, now runs.

Version 3 (line RG item T7-W, 2026-11-16) adds the signal of stream 9 (§5): a stream that fails
alone after the host received a chunk of it, which made the run fail in version 2, is discarded
while the run goes on. A host of version 2 does not know stream 9 and would keep the chunks of
such a stream as if they were its output, so the change is incompatible: the application refuses
a module of version 2, and a host of version 2 refuses a module of version 3. Every other run gives
the streams of version 2, except that a format which fails in the middle of its report fails alone
as one refused at a check before its rows did (§7), where the run failed, and that stream 3 of a run
of two searches holds each text once (§5), where version 2 wrote the format 0 search's warnings
again after a failure.

Version 4 (line WQA's sweep, D7, 2026-11-19) adds `losat_web2_register_inputs`, the *inputs*
response and the *start* field `register_inputs`: a run of `-parse_deflines`, which version 3 refused
(`-parse_deflines is not supported in the browser yet`), runs with the inputs that the new export reads
for its argv, and refuses the handles of `register`; in such a run LOSAT's refusal of outfmt 0 of a
record that is both a query and a subject and has masks fails stream 0 alone (§7). The signal of
stream 9 comes for every output stream that fails, alone or with its search in a run of two
searches, where version 3 sent it only after a chunk of the stream, and has 8 bytes, the stream and
the place of its error (stream 3, the file of `-logfile`, or nowhere; 2026-11-20), so that a host
knows which outputs a run could not write and where their errors are (§5). The signal of stream 10
says where each run's messages go (stream 3, the file of `-logfile`, or nowhere; 2026-11-20), so that a host does not take an
empty stream 3 for a run without warnings when TBLASTX's `-logfile` took them. A host of version 3 registers each input alone,
so its runs of the option would now fail rather than be refused at `validate`: the application
refuses a module of version 3, and a host of version 3 refuses a module of version 4. Runs without
the option give the streams of version 3. A `-logfile` name with a directory, which version 3 refused,
makes a file of that name on stream 8, or sends the messages to stream 3 (§4, *start*; line WQA's
sweep, D11, 2026-11-20): a stream 8 name may now contain `/`, which a host shows and saves as a name,
never as a path.

Stream 3 of a run that fails (line WQA's sweep D3, 2026-11-19) kept version 4: the adapter sends what
it held of the stream when the run fails, and a host that does not know it discards it with the
other streams of the failed run, as before. So did the outcome of a run of two searches, which is
that of the argv's own search (§7; D3 review, 2026-11-19): a run whose own format's search fails now
fails with that search's error, as the CLI, where it completed or failed with the tabular search's
error. This changes which runs complete, not a signature, a stream or a JSON field. So did stream 3
starting with the CLI's standard error for the argv (§5; D3 review r2, 2026-11-19): a run whose argv
asks for format 0 to 4 gives the formatting handler's lines of the tabular formats after its own
lines, where they came first, and a failed one does not give them; TBLASTX's file of `-logfile`
holds only the argv's lines, the rest coming on stream 3, NCBI's errors of the outputs that failed
included (stream 9 gives their place as `3`). A run that completes gives
the same lines, in another order where the argv asks for format 0 to 4. So did a run whose argv asks for a
structured format (§5; D3 review r3 M2, 2026-11-19): stream 3 adds the queries' warnings that the other
formats write, where it had none.

Version 7 (line WQA's sweep, D15, 2026-11-20; versions 5 and 6 are taken by other open branches, and
whichever merges later renumbers) keeps the output of a run that fails (§5): the adapter sends what the
writer of the argv's own output held when the run failed, so that the host holds the bytes that the
CLI leaves in that output before its error, and drops every other output of the run with a signal of
stream 9 whose place is `0xFFFFFFFE` (the error that ended the run). A host of version 6 or earlier
keeps nothing of a failed run, and takes the place `0xFFFFFFFE` for an error that no stream holds, so
that the run ends as an engine stop; the application therefore refuses a module of version 6 or
earlier, and a host of version 6 refuses a module of version 7. A run that completes gives the streams of version 6.

Version 8 (line WQA's D13 and D16, 2026-11-20; numbered 5 on its branch until D15's version 7 merged
first, 2026-11-21) reads the inputs of `register` and the index scan as
the CLI reads them: `scan_begin` takes the program, the role and the data loaders in place of a parser
kind, `register` takes the data loaders, and `validate` emits the *validation* response, which gives
them for the run's argv (`DATA_LOADERS=none` of a registry file, which version 7 and earlier ignored: it refused a
first line that may be a Seq-id, which the CLI then reads as data). The scan reads a `>?` gap line and
the lines that NCBI's line reader joins, which version 7 and earlier refused (`not supported in the browser`) while
the CLI reads them, gives the reader's refusals of a Seq-id first line and of a record over 2147483647
letters in the CLI's words, and no longer gives `line_layout`. A run refuses a handle registered with
the data loaders off when its argv has them on. A host of version 7 or earlier calls `scan_begin` and `register`
with other arguments, so the change is incompatible: the application refuses a module of version 7 or
earlier, and a host of version 7 refuses a module of version 8. A run of inputs that version 7 read gives
the same streams in version 8: only the responses on stream 2 changed, as above.

Version 10 (line WQA's D14, 2026-11-20; version 9 is taken by another open branch) adds stream 11 and the *start* field `prints`: an argv with a help or version
word, which version 8 and earlier refused (the parser's refusal of `-version` and `-xmlhelp`, line RG E3-2; the
application refused the help words), prints the CLI's text on stream 11, searches nothing, and its `run`
takes no input handles. A host of version 8 or earlier does not know stream 11 and would show such a run
as a search without output, so the change is incompatible: the application refuses a module of version 8
or earlier, and a host of version 8 refuses a module of version 10. Runs without these words give the
streams of version 8.

Version 11 (line WQA's D10, 2026-11-19; it kept version 4 on its branch until D13's version 8 merged
first, and version 9 until D14's version 10 merged first, 2026-11-21) reads a query as the command
line's query batches read it: `register` and `register_inputs` register the records before the reader's error of a query record with that error
(`reader_error` of *register* and *inputs*), where they failed, and the run ends with the error in
the record's query batch; `scan_end` gives the records before the reader's error and the error
(`reader_error` of *scan*, with whether the record it fails on starts with a defline), where it
failed with it; a query of blank and comment lines only is registered, and its run ends with
`Empty CBlastQueryVector`. A host of version 10 or earlier does not know `reader_error`: it would
take the records of a scan that ended at the reader's error for the whole input, and run them
without the error that ends the CLI's run, so the change is incompatible: the application refuses a
module of version 10 or earlier, and a host of version 10 refuses a module of version 11. A run of
inputs that version 10 read gives the same streams in version 11.

Version 12 (line WQA D12a, 2026-11-20; numbered 9 on its branch until D13's version 8 and D14's version 10 merged
first, 2026-11-21, and D10's version 11, 2026-11-22) runs `-db` and its filters (`-seqidlist`, `-negative_seqidlist`, `-taxids`,
`-negative_taxids`, `-taxidlist`, `-negative_taxidlist`, `-gilist`, `-negative_gilist`, BLASTN's database masks),
which versions 11 and earlier refused (`-db (a database search) is not supported in the browser yet`): the engine runs the
library's database search (`run_local_database_*`, the CLI's), opening the files that the argv names in the run's
working directory, which the host preopens as `/work`, the module's current directory (§3). The *start* field `database` says that a run takes no
subject: `run` takes subject handle `0` and refuses another, and the *inputs* `subject` is null. Its query is read as
version 8 reads every input, with the data loaders of the argv (`register` with those of the *validation*,
`register_inputs` with those of the argv's registry), and a run with them on refuses a query handle read with them off.
As in version 11, its records end at the reader's error of a query record, which the run meets in that record's query batch
(`run_local_database_*` take the query's end as `run_local_*` do).
A host of version 11 or
earlier gives no preopened directory and a subject handle to every run, so its runs of `-db` would fail at the database or
at the handle rather than be refused at `validate`: the application refuses a module of an earlier version, and a host of
an earlier version refuses a module of version 12. Runs without `-db` give the streams of version 11, except BLASTN's `-filtering_db` (a database of repeats), which version 11
refused: it runs with `-db` and with `-subject`, the engine opening the repeat database by the option's name in `/work`
(line RG item D4, 2026-11-22; no new version, since every host of version 12 gives the run directory). `-entrez_query` and
`-remote` stay refused, as on the CLI.

Version 13 (line NP WW3a, 2026-11-22) runs the tool `windowmasker` (`losat_web2_tool`, §4) with any command line that `losat windowmasker` accepts, on the run's files (§3): `-mk_counts` in each `-sformat`, `-convert`, the masking of `-ustat` in each `-outfmt`, with their options, help and version. Its standard output comes on stream 11, its standard error on stream 3 and its files on the new stream 12 as they are written (D-WW3-6: a tool's output may be larger than the browser's memory, human counts 724 MB, a genome masked as FASTA 3.1 GB), each write at its offset in its file, so that a file that the run opens twice is one file as on the command line (D-WW3-9), its exit status on stream 2 (*status*), and `validate` takes its argv (D-WW3-7, with the checks before the tool's files, D-WW3-10). A host of version 12 does not know the export, stream 12 or the *status*: the application refuses a module of an earlier version, and a host of version 12 refuses a module of version 13. Searches give the streams of version 12.

The *describe* rows of line RG item K4i-3a (2026-11-24) kept version 13: they have the same fields, and a host finds an
option by its `flag`. Their `help` is NCBI's comment of the option in place of the clap definitions' text, and they come
in the order of the usage line of the program's `-help`; psiblast's `-out_pssm` and `-out_ascii_pssm` are paths in the
library API's `describe` (`value_kind`, which this response does not carry). The application reads no value from `help`:
the hit list defaults that it read from the clap texts (`default: 500`) come from `defaults`.
