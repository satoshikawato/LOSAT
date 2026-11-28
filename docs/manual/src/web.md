# LOSAT Web (Browser GUI)

LOSAT Web is an in-browser graphical interface for local sequence searches built on LOSAT's WebAssembly engine (`losat`). It follows the familiar workflow of NCBI Web BLAST while running computations locally inside your browser.

> **Status: Preview.** LOSAT Web is outside the v0.2.0 CLI NCBI parity claim. See [the Web GUI boundary decision](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md).

## Key properties

- **Data stays on your machine.** Queries, subject FASTA files, local BLAST databases, headers and alignments never leave your browser over the network. Local database folders are read in place from your device without uploading or copying whole volumes into WebAssembly memory.
- **Same search engine.** LOSAT Web calls the WebAssembly build (`losat` npm package) of the same Rust engine used by the `losat` command line. Pairwise (`-outfmt 0`) and tabular (`-outfmt 6`, `7`) reports come directly from the engine's output bytes.
- **GUI features vs. NCBI output.** Graphic Summary, Descriptions, job queueing, sequence extraction, and CSV/JSON/session exports are LOSAT Web features labelled as LOSAT Web formats. Each run shows an exportable reproduction recipe so you can repeat the exact search with `losat` on the command line.
- **Per-run status badge.** Each search reports **Supported** when its program, task and options were compared against NCBI BLAST+ on LOSAT's test fixtures in the browser runtime; other combinations are labelled **Not yet supported**.

## What runs in LOSAT Web

- **Programs:** `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, and `psiblast` (selected via **PSI-BLAST** under BLASTP's Algorithm radios, matching NCBI BLAST), plus the `windowmasker` (`-mk_counts`, `-convert`, `-ustat`) and **Format Archive** (`blast_formatter`, `-archive`) tool tabs.
- **Format Archive (`blast_formatter`):** writes the reports of a search again from its `-outfmt 11` BLAST archive: a file, pasted ASN.1 text, or the `-outfmt 11` output of a search of this session (**Format archive…** beside **Export** on that output), in `-outfmt` `8`, `9` or `12` as the command line writes them. For the archive of a `-db` search, open the database's files under **Run directory**. `-rid` is refused.
- **Queries:** FASTA sequences (`-query`) for all six programs; in `psiblast`, also multiple sequence alignments (`-in_msa`, with `-msa_master_idx` and `-ignore_msa_master`) and PSSM checkpoints (`-in_pssm`); in `tblastn`, also PSSM checkpoints (`-in_pssm`). Generated PSSM files (`-out_pssm`, `-out_ascii_pssm`, with `-save_pssm_after_last_round` and `-save_each_pssm`) are captured as downloadable run files in Results (Outputs and Download All), matching native CLI bytes.
- **Subjects:** FASTA subjects (`-subject`) and local BLAST databases (`-db`) opened from a folder (with multi-file fallback). Database searches auto-discover valid `.nal`/`.pal`, `.nin`/`.pin` and multi-volume databases matching the selected program's molecule type, mount optional `taxdb.*` and `taxonomy4blast.sqlite3` files in the folder, and expose grouped database masking and filter controls: mutually exclusive database masking (`-db_soft_mask` or `-db_hard_mask`, by algorithm ID or name, for programs that support database masks), taxonomy filtering (`-taxids`, `-negative_taxids`, `-taxidlist`, `-negative_taxidlist`, `-no_taxid_expansion`), and ID list files (`-seqidlist`, `-negative_seqidlist`, `-gilist`, `-negative_gilist`, plus `-ipglist`/`-negative_ipglist` for protein databases). After a page reload, browsers do not keep folder handles, so LOSAT Web asks you to choose the folder again.
- **Repeats filtering in BLASTN:** local nucleotide repeats databases (`-filtering_db`), WindowMasker counts files (`-window_masker_db`), and downloaded or browser-made taxonomy ID counts (`-window_masker_taxid`), with both `-subject` and `-db` searches. When `-no_taxid_expansion` is enabled on a database search, LOSAT Web guides WindowMasker to a counts file (`-window_masker_db`) and highlights the NCBI incompatibility if `-window_masker_taxid` is also set.
- **Not available in LOSAT Web:**
  - Remote NCBI options (`-remote`, `-entrez_query`, `blast_formatter -rid`);
  - `psiblast -phi_pattern`;
  - DELTA-BLAST (`deltablast`), RPS programs (`rpsblast`, `rpstblastn`), `rpsbproc`, and database tools other than `windowmasker` and `blast_formatter` (`makeblastdb` runs in the `losat` npm package's engine, with every option and the command line's bytes, but has no tab in LOSAT Web yet; `blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask`, `makembindex`, `makeprofiledb`, `blastdb_aliastool`, `blastdbcheck`).

## Threads in the browser

LOSAT Web automatically uses multi-threaded WebAssembly (`losat-web-threads.wasm`) when the page is served with cross-origin isolation (`SharedArrayBuffer` enabled):

```text
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

Without those headers, or if the browser cannot reserve shared address space, multi-threaded requests automatically fall back to the single-threaded worker (`losat-web-serial.wasm`) and record the reason in `RuntimeInfo.fallbackReason`. See [WebAssembly and Browser](libraries/webassembly.md#what-a-page-must-do).

## Running LOSAT Web locally from the repository

Pre-built Web GUI bundles are not published separately in v0.2.0 releases. To build and run LOSAT Web from a clone of the repository (requires Rust `1.92.0` with `wasm32-wasip1` and `wasm32-wasip1-threads`, Python `3.11+`, and Node.js `22.18+`):

```bash
python3 web/adapter/tools/build_reactors.py --target-dir <build directory> --output-dir <modules directory>
cd web
npm ci
LOSAT_WEB_REACTORS=<modules directory> npm run build -w losat
npm run dev -w losat-web
```

For production bundling, `npm run build -w losat-web` outputs static assets under `web/app/dist/` (`vite preview --port 4173` serves a local preview).

## See also

- [WebAssembly and Browser](libraries/webassembly.md) — `losat` npm package and `WasmEngine` host API
- [Platforms and WebAssembly](platforms.md) — platform support matrix and WASI CLI runners
- [gbdraw](https://github.com/satoshikawato/gbdraw) — genome diagram tool whose browser edition uses the same LOSAT WebAssembly engine
