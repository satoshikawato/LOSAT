# Platforms and WebAssembly

LOSAT runs as a native program on five platforms, as a WebAssembly (WASI) command in two forms, and in the browser as a Web GUI. This page lists what each form is, how well it is covered, and what to know before you use it.

## Status at a glance

| Platform | File | Status in v0.2.0 |
|:---|:---|:---|
| Linux x86_64 | `LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz` | Supported: matches NCBI BLAST+ 2.17.0 on the release test set (apart from the version line, which names LOSAT) |
| Linux aarch64 | `LOSAT-0.2.0-aarch64-unknown-linux-gnu.tar.gz` | Supported |
| Windows x86_64 | `LOSAT-0.2.0-x86_64-pc-windows-msvc.zip` | Supported |
| macOS arm64 | `LOSAT-0.2.0-aarch64-apple-darwin.tar.gz` | Supported |
| macOS x86_64 | `LOSAT-0.2.0-x86_64-apple-darwin.tar.gz` | Supported |
| WASI serial (`wasm32-wasip1`) | `LOSAT-0.2.0-wasm32-wasip1.tar.gz` | Supported: matches NCBI BLAST+ 2.17.0 on the release test set (apart from the version line, which names LOSAT) |
| WASI threads (`wasm32-wasip1-threads`) | `LOSAT-0.2.0-wasm32-wasip1-threads.tar.gz` | Supported: tested in CI, outside the release's NCBI parity claim |
| Web GUI (browser) | not a release file | Preview: outside the release's NCBI parity claim |

On Linux aarch64, 43 cases were compared with NCBI BLAST+ on that platform, and every other case gives the same bytes as Linux x86_64. Windows and macOS are Supported on the same terms, apart from the last-digit case under macOS below.

## Native platforms

Download the file for your platform from the release page and follow [Installation](installation.md).

### Linux

Both Linux files are built for glibc systems and need no other software. [Threads, speed and memory](performance.md) explains `-num_threads`.

### macOS

On macOS, an E value that is printed with full precision (`-outfmt 11` and `-outfmt 8`) can differ from other platforms by 1 in its last (15th) digit. This comes from Apple's maths library, and NCBI BLAST+ on macOS is expected to show the same (not tested). Every other value is the same. The decision is recorded in [the decision record on the last digit of full-precision values](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-PLATFORM-LIBM.md).

The macOS binaries are not code-signed or notarized. If macOS stops a downloaded file, see [Troubleshooting](troubleshooting.md).

### Windows

The Windows file is a `.zip` with `losat.exe`. `makeblastdb` writes its times and the database date in local time, as NCBI BLAST+ prints them on Windows. `losat --install-shims` makes copies of `losat.exe` under the NCBI names, because Windows links need a privilege. In `BLASTDB`, separate directories with `;`.

### Output differences between platforms

NCBI's own builds for Windows, macOS and Linux aarch64 can differ slightly from its Linux x86_64 build. LOSAT's output is checked against the Linux x86_64 bytes on every platform. See [the decision record on differences between NCBI builds](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-NCBI-PLATFORM-VARIANCE.md).

## WASI: `losat.wasm`

The two WASI files are command modules for runtimes that implement WASI preview 1. They run the same programs as the native `losat`, with the same arguments. Because a module has no program name to link, the first argument is the program: `blastn -query ...`, not `losat blastn -query ...`.

| Module | Threads | Needs from the host |
|:---|:---|:---|
| `wasm32-wasip1` (serial) | one | WASI preview 1 |
| `wasm32-wasip1-threads` | several (`-num_threads`) | WASI preview 1, shared memory and the `wasi` thread-spawn call |

The serial module gives the same bytes as the native Linux x86_64 binary.

### What a WASI host must pass

A WASI module sees only what its host gives it.

- **Files.** The host must open the directories the module reads and writes (a preopened directory). The release tests give the module `/`, so absolute and relative file names work as they do natively.
- **Working directory.** Pass the directory you start in as the environment variable `PWD`. Relative file names, and the `New DB name:` that `makeblastdb` prints, start from `PWD`. Without `PWD` they start from `/`.
- **Time zone.** Pass `TZ`, or give the module the zone files (`/etc/localtime` and `/usr/share/zoneinfo`). `makeblastdb` then writes local times. Without them it writes UTC.
- **Temporary files.** A search program decompresses a `.gz`, `.bz2` or `.zst` input into a temporary file beside the input. A WASI command has no temporary directory to fall back on, so if the input's folder refuses the file, it stops with an error. Give the module a writable folder for the input, or pipe the data to `-query -`.
- **Memory.** A module can grow past 2 GiB. A host must then pass pointers and lengths to WASI calls as unsigned 32-bit numbers. Node's `node:wasi` refuses the signed values that the plain WebAssembly interface gives for 2 GiB and above, so the repository's Node.js runner (below) converts them.

### Run with Node.js

The repository's Node.js runners handle the points above. They are in `LOSAT/tests`, so clone the repository to use them. They need Node.js 24 or later.

```bash
git clone https://github.com/satoshikawato/LOSAT.git

# serial module
node LOSAT/tests/run_losat_wasi.js losat.wasm blastn -query query.fna -subject subject.fna -outfmt 6

# threaded module
node LOSAT/tests/run_losat_wasi_threads.js losat-threads.wasm blastn -query query.fna -subject subject.fna -num_threads 4 -outfmt 6
```

Both archives contain a file called `losat.wasm`. Rename them, or keep them in separate folders, as the two examples above suggest. The runner takes the module first, then the arguments of `losat`. It passes the environment, including `TZ`, sets `PWD` to the current directory, and gives the module `/`. Node.js prints a warning that WASI is experimental; this is Node's message and does not change the output.

### Other WASI hosts

Any host that implements WASI preview 1 can run `losat.wasm`, if it gives the module the working directory and the files (see above). LOSAT is tested with Node.js, through the runners in `LOSAT/tests` described above; other hosts are not part of the tests. The threaded module needs a host that supports shared memory and WASI threads; use the Node.js runner unless your host documents that support.

### A Node.js difference with `..`

When the host is Node.js (`node:wasi`), the host program opens the files, and Node reads `..` in a path as text without looking at the directory before it. So `makeblastdb -out gone/../db` builds `db`, and `-query gone/../q.fa` reads `q.fa`, even if `gone` does not exist. NCBI BLAST+ and the native `losat` on Linux stop with an error there. This only matters for paths that contain `..` after a directory that is not there.

## Web GUI (Preview)

LOSAT also runs in a web browser; see the dedicated [LOSAT Web (Browser GUI)](web.md) page for capabilities, build instructions, and browser isolation details. The Web GUI is a preview outside the release files (`web/app`), following the NCBI Web BLAST search layout (Graphic Summary and Descriptions on one page).

**Status: Preview.** It is outside the release's NCBI parity claim.

- **What each run says.** A run is marked "Supported" only when its program, task and options were compared with NCBI BLAST+ on LOSAT's test fixtures, for every output format the run wrote, and its browser runtime was compared with the native engine. Other runs are marked "Not yet supported". The mark is about the build and the options. It is not a comparison of your own inputs with NCBI.
- **Same engine.** The browser runs the same engine as the command line. It receives the same options and whole FASTA records, and the pairwise text and tabular outputs (`-outfmt 0`, `6` and `7`) are the engine's bytes. You can export the exact input of a run and repeat it with the command line.
- **Application features.** The job queue, sequence extraction, CSV, JSON and report exports, and session files belong to the Web GUI. They are labelled as LOSAT Web formats and are not NCBI output.
- **No databases.** The Web GUI rejects `-db`; database searches are native only.
- **Your data stays in your browser.** Sequences, headers, file names and results are not sent over the network.

The boundary of the Web GUI is set by [the decision record on the Web GUI boundary](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md). The same WebAssembly engine powers the browser version of [gbdraw](https://github.com/satoshikawato/gbdraw), a genome diagram tool.
