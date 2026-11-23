# Troubleshooting

Questions and answers for the problems people meet most. Each answer is based on what LOSAT v0.2.0 does. If your problem is not here, see [Where to report a problem](#where-to-report-a-problem).

## LOSAT says an option is not supported

**What does it look like?**

```text
Error: the NCBI BLAST+ option -html is not supported by LOSAT's BLASTN
```

The command exits with status 1 and runs nothing. For an output format that NCBI writes and LOSAT does not (in `psiblast` or `tblastn -in_pssm`), see [Output formats](output-formats.md) and [`psiblast` refuses some `-outfmt` values](#psiblast-refuses-some--outfmt-values).

**Why?** LOSAT stops with an error for each option and value it does not implement. It never ignores an option and never runs something else in its place, so a result you get is the result you asked for.

**What can I do?** Look the option up in [What LOSAT refuses](refusals.md) and on the page of your program. Options that need NCBI's servers or its GI and IPG databases (`-remote`, `-entrez_query`, `-gilist`, `-ipglist`) are not available. For an option that is refused, use NCBI BLAST+ for that run. The [Introduction](introduction.md) lists the programs and options planned for a later release.

A different message, NCBI's short usage followed by `Error: Unknown argument: "xyz"` (exit status 1), means the option does not exist. Check the spelling with `losat <program> -help`. LOSAT takes NCBI's single-dash names (`-evalue`, not `--evalue`).

## `losat --install-shims` stops with "nothing was installed"

```text
losat --install-shims: these files exist and are not this executable; nothing was installed: /path/to/bin/blastn
```

A file with one of the 16 program names is already in the folder, often NCBI's own `blastn`. LOSAT never overwrites a file that is not a LOSAT link or copy, and it makes none of the names when one is in the way. Exit status is 1.

You can:

- give another folder, `losat --install-shims ~/losat-bin`, and put that folder first on your `PATH`; or
- rename or remove the file that is in the way, then run the command again.

Links made by an earlier LOSAT are replaced without a message, so you can run the command again after an upgrade. See [Coming from NCBI BLAST+](from-ncbi-blast.md).

## `blastn: command not found`, or the wrong `blastn` answers

The names exist only after `losat --install-shims`. Make sure the folder you gave is on your `PATH`, and before any folder with NCBI's programs if you want LOSAT to answer. `blastn -version` shows who answers: LOSAT's first line is `blastn: LOSAT 0.2.0` and its second line starts with `Package: LOSAT`.

If you moved or replaced `losat` after making the links, run `losat --install-shims DIR` again.

## The database is not found

```text
LOSAT database error: No alias or index file found for protein database [nodb] in search path [...]
```

The message is NCBI's; where NCBI starts it with `BLAST Database error:`, LOSAT starts it with `LOSAT database error:`. Check the database name, that it is a protein or nucleotide database as the program needs (`blastp` reads protein, `blastn` reads nucleotide), and where LOSAT looks for it: the current directory, then the directories in `BLASTDB`, then the configuration in `.ncbirc`, as NCBI BLAST+ does. See [BLAST databases](databases.md).

## A `-seqidlist` warning names `blastdb_aliastool`

With a text file in `-seqidlist` or `-negative_seqidlist`, LOSAT writes NCBI's warning on the standard error (`Warning: [blastp] To obtain better run time performance, please run blastdb_aliastool -seqid_file_in <INPUT_FILE_NAME> -seqid_file_out <OUT_FILE_NAME> and use <OUT_FILE_NAME> as the argument to -seqidlist`). The search runs normally. `blastdb_aliastool` is NCBI's tool; LOSAT does not provide it. To get a binary list, run NCBI BLAST+'s `blastdb_aliastool`. LOSAT reads both the text and the binary list.

## A compressed input fails

A `-query`, `-subject` or `makeblastdb -in` file whose name ends in `.gz`, `.bz2` or `.zst` is decompressed by LOSAT, and the name alone decides. If the file is not really compressed, LOSAT stops before the search with a message such as `-query: 'x.fa.gz' is named .gz but is not gzip data` (exit status 1). Rename the file, or decompress it yourself.

A search program writes the decompressed data into a temporary file next to the input first. If that folder refuses it or has no room, LOSAT uses the folder in `TMPDIR`, so point `TMPDIR` to one with room. `makeblastdb` streams and needs no copy. A WASI command has no `TMPDIR` to fall back on: if the input's folder refuses the copy, it stops; see [WASI host problems](#wasi-host-problems).

## macOS blocks the downloaded program

Files you download with a browser can carry a quarantine mark. macOS (Gatekeeper) can then refuse to start `losat`, with a message that it cannot be opened because it comes from an unidentified developer. The release binaries are not code-signed or notarized, so this can happen.

If you trust the file you downloaded, remove the mark. In a terminal, go to the folder with `losat` and run:

```bash
xattr -d com.apple.quarantine ./losat
```

Use `xattr -dr com.apple.quarantine LOSAT-0.2.0-aarch64-apple-darwin` for the whole unpacked folder. Then run `./losat --version`.

The maintainers cannot test on Apple hardware, so this is the standard macOS step and not something LOSAT checks. To be sure the file is the one the project published, [check it](installation.md#check-a-download) against `SHA256SUMS` first.

## Windows

- **A warning when you run `losat.exe`.** The Windows binary is not code-signed, so Windows can ask you to confirm before it runs a downloaded file. This is Windows behaviour, and the same advice applies as on macOS: check the file against `SHA256SUMS` first.
- **`losat` is not found.** Unzip the file, then put the folder with `losat.exe` on your `PATH`, or run it by its full path.
- **The NCBI names.** On Windows `losat --install-shims` makes copies of `losat.exe`, not links. After you replace `losat.exe`, run the command again to refresh the copies.
- **`BLASTDB` with several folders.** Separate them with `;` on Windows (`:` elsewhere).
- **Times.** `makeblastdb` writes its times and the database date in local time, as NCBI BLAST+ does on Windows.

## WASI host problems

These apply when you run `losat.wasm` yourself. See [Platforms and WebAssembly](platforms.md) for the full list of what a host must pass.

- **"File not found" for a file that exists.** The host must preopen the directory that holds it. The release tests give the module `/`.
- **`makeblastdb` writes UTC times, or `New DB name:` shows `/...`.** The host did not pass `TZ` (or the zone files) and `PWD`. Pass both.
- **A compressed input stops the search.** The WASI command writes the decompressed copy beside the input and has no temporary directory to fall back on. Make the input's folder writable for the module, or pipe the data to `-query -`.
- **A path with `..` works when it should not.** With Node.js, `..` is read as text without looking at the directory before it, so `-query gone/../q.fa` reads `q.fa` even if `gone` does not exist. NCBI BLAST+ and native LOSAT stop with an error.
- **Calls fail with `EINVAL` on large runs.** A host that passes pointers as signed numbers breaks once the module uses 2 GiB or more of memory. Pass WASI parameters as unsigned 32-bit values; the repository's Node.js runner does.
- **The threaded module will not start.** It needs a host with shared memory and WASI threads. It declares a shared memory of up to 4 GiB less one page, and a host that cannot reserve that much address space fails when it creates the memory. The repository's Node.js runner (`LOSAT/tests/run_losat_wasi_threads.js`) retries with smaller maximums.

## Memory and threads

- **`-num_threads` sets the number of threads.** The default is 1. The output is the same at every thread count, and it is what NCBI BLAST+ writes with `-num_threads 1` (apart from the version line, which names LOSAT). `-mt_mode` is accepted and changes nothing in the output (`blastn`, `blastp`, `blastx` and `tblastn` take 0, 1 or 2; `rpsblast` and `rpstblastn` take 0 or 1; `tblastx` and `psiblast` have no `-mt_mode`, as in NCBI BLAST+ 2.17.0). `-num_threads` takes an integer of 1 or more (`rpsblast` also reads 0 as one thread for each volume of the database).
- **Threads and memory.** The page [Threads, speed and memory](performance.md) explains how threads and memory use relate.
- **A search of a compressed file needs room for the copy.** See "A compressed input fails" above.
- **Development variables.** Variables that start with `LOSAT_X_`, and `LOSAT_TIMING`, are for development; leave them unset.

## The output differs from NCBI BLAST+

First check [Differences from NCBI BLAST+](differences.md). A few differences are on purpose, for example a local `tblastx` or `tblastn` subject translated with `-db_gencode`, the last digit of some full-precision values on macOS, and the `CL:` tag of SAM output. If your case is not covered there, it is a bug: please report it.

## `psiblast` refuses some `-outfmt` values

`psiblast` writes `-outfmt 0`, `5` to `11` and `20` from a FASTA query, `-in_pssm` or `-in_msa`, at any `-num_iterations`. `tblastn -in_pssm` writes `0` to `10` and `20`.

- A `psiblast` format `1` to `4`, `12` to `16` or `18` stops before the search: `Error: output format N is not supported by LOSAT's PSIBLAST yet` with a FASTA query, and `Error: output format N with -in_pssm is not supported by LOSAT's PSIBLAST yet` (`with -in_msa` for an alignment) with `-in_pssm` or `-in_msa`. The `-in_pssm` and `-in_msa` refusals come after the file is read, so a file that NCBI BLAST+ cannot read gets NCBI's error for it first.
- A `tblastn -in_pssm` format `12` to `16` or `18` stops after the PSSM is read: `Error: output format N with -in_pssm is not supported by LOSAT's TBLASTN yet`. `-outfmt 11` stops there too, because NCBI's own run of it stops with a C++ exception whose report changes from run to run.
- Each exits with status 1.

Use NCBI BLAST+ if your workflow requires those formats.

## `rpsblast` or `rpstblastn` rejects `-subject`

Both `rpsblast` and `rpstblastn` search profile databases (`-db`, such as CDD) only and have no `-subject` option, matching NCBI BLAST+ 2.17.0. Passing `-subject` prints NCBI's usage and `Unknown argument: "subject"`, and exits with status 1, as NCBI BLAST+ does.

## `blastn -use_index true` stops with an error

`losat blastn -use_index true` searches with MegaBLAST index files built by `makembindex`. It refuses the search when an index volume is damaged or does not match the database, and a later query batch of an index of only some volumes of `-db` where NCBI BLAST+ writes to memory it has freed (see [What LOSAT refuses](refusals.md#blastn)). Omit `-use_index true` to search without an index.

## `rpsblast`, `rpstblastn`, or `rpsbproc` is missing in Python (`pylosat`), Rust (`losat::api`), or LOSAT Web

In v0.2.0, `losat::api` and `pylosat` run the six core search programs (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`) and wrap the seven tools (`makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask` and `makembindex`). The `losat` npm package and LOSAT Web run `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx` (`psiblast` is not available in the browser engine). `rpsblast`, `rpstblastn` and `rpsbproc` run on the command line only (through `losat` or its installed shims).

## A full-precision E-value on macOS differs in the 15th digit

On macOS (`-outfmt 8` or `-outfmt 11`), the last (15th) significant digit of a full-precision E-value can differ by 1 because of Apple's maths library. Tabular, JSON, XML and pairwise formats are unaffected in practice, and Linux, Windows and WASI match Linux x86_64. See [Differences from NCBI BLAST+](differences.md#last-digit-of-a-full-precision-value-on-macos) and [Platforms and WebAssembly](platforms.md#macos).

## Where to report a problem

Open an issue at [github.com/satoshikawato/LOSAT/issues](https://github.com/satoshikawato/LOSAT/issues) when LOSAT fails or its output differs from NCBI BLAST+ 2.17.0 for the same command. The issue form asks for:

1. the exact LOSAT command line and the output of `losat --version`;
2. the matching NCBI BLAST+ command and the output of `<program> -version` (write "not run" if you did not run it);
3. a small input that shows the problem (attach files if they are long);
4. the output you expected, and the output you got (for a difference, the lines that differ);
5. your platform: operating system and architecture, or the WebAssembly host (browser, Node.js or another host).

Use small public or synthetic sequences. Do not post confidential, restricted or proprietary data. Report a security vulnerability privately, through the repository's Security tab ("Report a vulnerability"), and not in a public issue.
