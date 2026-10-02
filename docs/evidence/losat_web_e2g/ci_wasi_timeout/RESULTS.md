# E2g WASI threaded TBLASTX p11 timeout investigation

## A1/A2. CI run facts
Note: run 36882087513 ran at headSha eb262bac2 (the "c6c3386cf" in the task is 1 commit later, E2g docs-only; eb262bac2 is its parent). Checkout ref in log = eb262bac2.

| run | headSha | branch | date | conclusion | runner image | node | frozen step |
|---|---|---|---|---|---|---|---|
| 35201795582 | 6bfb1b09b | implement-v010-distribution | 2026-09-17 | success | ubuntu-24.04 20260907.300.1 (runner 2.337.0) | node 24 (action default; matrix-job log) | 08:54:36 -> 10:23:52 (89m16s) |
| 36745543144 | 7693a9c73 | feature/losat-web-gui | 2026-09-30 | failure (blastn Sakai.MG1655 mismatch, see below) | ubuntu-24.04 20260927.320.1 | 24.x | 16:44:16 -> 18:37:17 (113m) |
| 36882087513 | eb262bac2 | feature/losat-web-gui | 2026-10-01 | failure (p11 threaded TIMEOUT) | ubuntu-24.04 20260927.320.1 | 24.21.0 | 15:17:30 -> 17:07:43 (110m) |

Command (all runs): `python3 LOSAT/tests/check_wasm_threading_regressions.py --native LOSAT/target/release/LOSAT --threaded .tmp/threading/artifacts/losat-threaded-command.wasm --oracle-dir .tmp/threading/oracle/ncbi-blast-2.17.0+/bin --jobs 3 --output-dir .tmp/threading/frozen-regressions`  (no --timeout-seconds on the command line; default 3600 s per search, seen in the runs.json record 'timeout_seconds': 3600). CPU model is not printed in the logs (ubuntu-24.04 standard hosted runner = 4 vCPU; not verified from logs).
Job step timing JSON: run-<id>.json; full logs: log-<id>.txt in this dir.

Completion timestamps of tblastx labels (UTC; lines printed when a label finishes):

### 35201795582 (green)
native p11 PASS 09:36:01; threaded p11 PASS 10:01:43 ; threaded d01 PASS 09:45:32; threaded d02 10:00:16; native d01 09:31:22
### 36745543144
native p11 PASS 17:20:07; threaded p11 NEVER printed; threaded d01 PASS 17:42:00 ; threaded d02 18:03:02
### 36882087513
native p11 PASS 16:05:15; threaded p11 FAIL TIMEOUT printed 16:32:20 (so it started ~15:32:20, right when threaded p10/p06 finished); threaded d01 PASS 16:22:55
Full label list: a2_labels.txt

## A3. Artifacts and per-case timings (frozen-regressions/runs.json)
Artifacts downloaded to artifacts-<id>/ (old threaded wasm + its host JS kept in artifacts-35201795582/artifacts/).
threaded losat-threaded-command.wasm sha256:
- green 35201795582: d49206b531da26977849857fdfc81f74ae5bd5efacb9f0e881e41486cfcbf51f (native sha fb01e0f6...)
- red 36745543144:   85483126c0ad02db871481b7678b99b9b61984304c0d023afbdf00e10c39553c (native 4721d606...)
- red 36882087513:   85483126c0ad02db871481b7678b99b9b61984304c0d023afbdf00e10c39553c (same wasm and native as 36745543144)

KEY FACTS
- threaded p11 in the GREEN run took 3358 s wall (6617 CPU-s, peak RSS 362 MB): only 242 s (7%) under the 3600 s limit.
- threaded p11 TIMED OUT (3600.0 s) in BOTH red runs (36745543144 too; its log never printed the label because the script raised on the blastn Sakai.MG1655 mismatch first; runs.json records TIMEOUT).
- Timeout records have cpu=None (process killed) so CPU of the red p11 is unknown.
- The same wasm/native binary was used in both red runs, so red-vs-red differences (native p11: 1413 s vs 2060 s; ratio 1.46) are pure runner/concurrency noise.
- Almost every other label (native AND threaded) is 1.1-1.9x slower in the red runs than in the green run: e.g. native d01 1497 -> 2327/1898 s, threaded d01 850 -> 1313/1125, threaded d05 534 -> 792/692, native p14 38 -> 73/74, but native p11 1909 -> 1413/2060.
- Concurrency identical in all three: --jobs 3; threaded p11 always overlapped (2 other jobs for its whole window, avg 2.0) native p11 + native d01 + d02 chain (see a3_overlap.md).

Full table in a3_table.md (copied below):

Wall seconds / CPU seconds (user+sys) per tblastx label, from frozen-regressions/runs.json

| label | green 35201795582 | red 36745543144 | red 36882087513 | wall ratio red1/green | red2/green |
|---|---|---|---|---|---|
| native/tblastx/d01_nz_self_code4 | 1497s / cpu 1055 | 2327s / cpu 1503 | 1898s / cpu 1347 | 1.55 | 1.27 |
| threaded/tblastx/d01_nz_self_code4 | 850s / cpu 1286 | 1313s / cpu 1809 | 1125s / cpu 1672 | 1.55 | 1.32 |
| native/tblastx/d02_ap027132_nz_code4 | 1737s / cpu 1412 | 2389s / cpu 2017 | 2356s / cpu 1797 | 1.38 | 1.36 |
| threaded/tblastx/d02_ap027132_nz_code4 | 884s / cpu 1678 | 1262s / cpu 2376 | 1149s / cpu 2173 | 1.43 | 1.30 |
| native/tblastx/d03_ap027078_ap027131_code4 | 736s / cpu 580 | 1013s / cpu 898 | 835s / cpu 768 | 1.38 | 1.13 |
| threaded/tblastx/d03_ap027078_ap027131_code4 | 321s / cpu 676 | 513s / cpu 1020 | 506s / cpu 898 | 1.60 | 1.57 |
| native/tblastx/d04_ap027131_ap027133_code4 | 174s / cpu 161 | 232s / cpu 208 | 279s / cpu 207 | 1.33 | 1.60 |
| threaded/tblastx/d04_ap027131_ap027133_code4 | 104s / cpu 217 | 128s / cpu 284 | 148s / cpu 279 | 1.23 | 1.42 |
| native/tblastx/d05_ap027133_ap027132_code4 | 869s / cpu 773 | 1309s / cpu 1186 | 1140s / cpu 1045 | 1.51 | 1.31 |
| threaded/tblastx/d05_ap027133_ap027132_code4 | 534s / cpu 1076 | 792s / cpu 1603 | 692s / cpu 1400 | 1.48 | 1.30 |
| native/tblastx/d06_ap027131_ap027133_db4 | 36s / cpu 36 | 51s / cpu 45 | 48s / cpu 48 | 1.43 | 1.34 |
| threaded/tblastx/d06_ap027131_ap027133_db4 | 37s / cpu 56 | 45s / cpu 71 | 45s / cpu 71 | 1.21 | 1.23 |
| native/tblastx/p01_ap027280_self | 48s / cpu 43 | 63s / cpu 56 | 59s / cpu 52 | 1.31 | 1.23 |
| threaded/tblastx/p01_ap027280_self | 39s / cpu 61 | 47s / cpu 76 | 54s / cpu 86 | 1.19 | 1.38 |
| native/tblastx/p02_mje_mela | 21s / cpu 21 | 31s / cpu 29 | 29s / cpu 27 | 1.46 | 1.36 |
| threaded/tblastx/p02_mje_mela | 22s / cpu 35 | 28s / cpu 44 | 31s / cpu 47 | 1.27 | 1.39 |
| native/tblastx/p03_mela_pemojnva | 3s / cpu 3 | 6s / cpu 6 | 6s / cpu 6 | 1.86 | 1.81 |
| threaded/tblastx/p03_mela_pemojnva | 6s / cpu 9 | 7s / cpu 11 | 8s / cpu 12 | 1.24 | 1.34 |
| native/tblastx/p04_pemojnva_pesemjnv | 23s / cpu 22 | 27s / cpu 26 | 33s / cpu 31 | 1.18 | 1.43 |
| threaded/tblastx/p04_pemojnva_pesemjnv | 31s / cpu 47 | 37s / cpu 63 | 41s / cpu 67 | 1.18 | 1.32 |
| native/tblastx/p05_pesemjnv_pemojnva | 61s / cpu 59 | 69s / cpu 63 | 78s / cpu 72 | 1.12 | 1.27 |
| threaded/tblastx/p05_pesemjnv_pemojnva | 68s / cpu 113 | 99s / cpu 157 | 110s / cpu 175 | 1.46 | 1.61 |
| native/tblastx/p06_pemojnva_lvmjnv | 452s / cpu 299 | 438s / cpu 299 | 567s / cpu 373 | 0.97 | 1.25 |
| threaded/tblastx/p06_pemojnva_lvmjnv | 515s / cpu 1179 | 596s / cpu 1329 | 690s / cpu 1560 | 1.16 | 1.34 |
| native/tblastx/p07_lvmjnv_trcumjnv | 6s / cpu 3 | 5s / cpu 3 | 6s / cpu 4 | 0.95 | 1.12 |
| threaded/tblastx/p07_lvmjnv_trcumjnv | 7s / cpu 6 | 8s / cpu 6 | 10s / cpu 7 | 1.07 | 1.40 |
| native/tblastx/p08_trcumjnv_mellatmjnv | 22s / cpu 14 | 27s / cpu 18 | 28s / cpu 18 | 1.22 | 1.27 |
| threaded/tblastx/p08_trcumjnv_mellatmjnv | 24s / cpu 28 | 26s / cpu 34 | 31s / cpu 38 | 1.09 | 1.27 |
| native/tblastx/p09_mellatmjnv_meenmjnv | 96s / cpu 62 | 121s / cpu 84 | 126s / cpu 84 | 1.26 | 1.31 |
| threaded/tblastx/p09_mellatmjnv_meenmjnv | 84s / cpu 106 | 108s / cpu 129 | 116s / cpu 138 | 1.28 | 1.38 |
| native/tblastx/p10_meenmjnv_mejomjnv | 157s / cpu 107 | 238s / cpu 145 | 229s / cpu 143 | 1.52 | 1.46 |
| threaded/tblastx/p10_meenmjnv_mejomjnv | 134s / cpu 242 | 156s / cpu 246 | 177s / cpu 308 | 1.16 | 1.32 |
| native/tblastx/p11_avclpv_psclpv | 1909s / cpu 1248 | 1413s / cpu 997 | 2060s / cpu 1472 | 0.74 | 1.08 |
| threaded/tblastx/p11_avclpv_psclpv | 3358s / cpu 6617 | 3600s **TIMEOUT** | 3600s **TIMEOUT** | 1.07 | 1.07 |
| native/tblastx/p12_lc738874_lc738875_default | 6s / cpu 6 | 7s / cpu 7 | 12s / cpu 8 | 1.21 | 1.92 |
| threaded/tblastx/p12_lc738874_lc738875_default | 11s / cpu 17 | 14s / cpu 17 | 17s / cpu 20 | 1.36 | 1.63 |
| native/tblastx/p13_mela_mje_reverse | 18s / cpu 18 | 35s / cpu 23 | 32s / cpu 23 | 1.98 | 1.79 |
| threaded/tblastx/p13_mela_mje_reverse | 18s / cpu 27 | 31s / cpu 34 | 34s / cpu 35 | 1.72 | 1.86 |
| native/tblastx/p14_ap027131_ap027133_query4 | 38s / cpu 38 | 73s / cpu 49 | 74s / cpu 50 | 1.93 | 1.97 |
| threaded/tblastx/p14_ap027131_ap027133_query4 | 37s / cpu 57 | 68s / cpu 73 | 52s / cpu 75 | 1.83 | 1.40 |

threaded/native wall ratio per case (green / red1 / red2):
  d01_nz_self_code4: 0.57 / 0.56 / 0.59
  d02_ap027132_nz_code4: 0.51 / 0.53 / 0.49
  d03_ap027078_ap027131_code4: 0.44 / 0.51 / 0.61
  d04_ap027131_ap027133_code4: 0.60 / 0.55 / 0.53
  d05_ap027133_ap027132_code4: 0.61 / 0.60 / 0.61
  d06_ap027131_ap027133_db4: 1.02 / 0.87 / 0.93
  p01_ap027280_self: 0.82 / 0.74 / 0.92
  p02_mje_mela: 1.04 / 0.90 / 1.06
  p03_mela_pemojnva: 1.82 / 1.21 / 1.35
  p04_pemojnva_pesemjnv: 1.36 / 1.36 / 1.25
  p05_pesemjnv_pemojnva: 1.12 / 1.45 / 1.41
  p06_pemojnva_lvmjnv: 1.14 / 1.36 / 1.22
  p07_lvmjnv_trcumjnv: 1.31 / 1.47 / 1.63
  p08_trcumjnv_mellatmjnv: 1.08 / 0.97 / 1.09
  p09_mellatmjnv_meenmjnv: 0.88 / 0.89 / 0.92
  p10_meenmjnv_mejomjnv: 0.85 / 0.65 / 0.77
  p11_avclpv_psclpv: 1.76 / 2.55 / 1.75
  p12_lc738874_lc738875_default: 1.72 / 1.93 / 1.46
  p13_mela_mje_reverse: 1.03 / 0.90 / 1.07
  p14_ap027131_ap027133_query4: 0.99 / 0.94 / 0.71

Overlap during threaded p11 (a3_overlap.md):
```
### 35201795582 threaded p11 2026-09-17T09:05:45Z -> 2026-09-17T10:01:43Z (3358s)
   09:04:12Z 09:06:26Z threaded/tblastx/p10_meenmjnv_mejomjnv 41
   09:04:12Z 09:36:01Z native/tblastx/p11_avclpv_psclpv 1816
   09:06:26Z 09:31:22Z native/tblastx/d01_nz_self_code4 1496
   09:31:22Z 09:45:32Z threaded/tblastx/d01_nz_self_code4 850
   09:36:01Z 10:04:58Z native/tblastx/d02_ap027132_nz_code4 1542
   09:45:32Z 10:00:16Z threaded/tblastx/d02_ap027132_nz_code4 884
   10:00:16Z 10:12:32Z native/tblastx/d03_ap027078_ap027131_code4 87
   other-job overlap seconds total 6716 avg concurrent others 2.0
### 36745543144 threaded p11 2026-09-30T16:57:03Z -> 2026-09-30T17:57:03Z (3600s)
   16:47:14Z 16:57:10Z threaded/tblastx/p06_pemojnva_lvmjnv 7
   16:56:34Z 17:20:07Z native/tblastx/p11_avclpv_psclpv 1384
   16:57:10Z 17:35:57Z native/tblastx/d01_nz_self_code4 2327
   17:20:07Z 17:42:00Z threaded/tblastx/d01_nz_self_code4 1313
   17:35:57Z 18:15:47Z native/tblastx/d02_ap027132_nz_code4 1266
   17:42:00Z 18:03:02Z threaded/tblastx/d02_ap027132_nz_code4 903
   other-job overlap seconds total 7200 avg concurrent others 2.0
### 36882087513 threaded p11 2026-10-01T15:32:20Z -> 2026-10-01T16:32:20Z (3600s)
   15:29:35Z 15:32:32Z threaded/tblastx/p10_meenmjnv_mejomjnv 12
   15:30:55Z 16:05:15Z native/tblastx/p11_avclpv_psclpv 1975
   15:32:32Z 16:04:10Z native/tblastx/d01_nz_self_code4 1898
   16:04:10Z 16:22:55Z threaded/tblastx/d01_nz_self_code4 1125
   16:05:15Z 16:44:31Z native/tblastx/d02_ap027132_nz_code4 1625
   16:22:55Z 16:42:04Z threaded/tblastx/d02_ap027132_nz_code4 565
   other-job overlap seconds total 7200 avg concurrent others 2.0
```

## A4. What differs between green 6bfb1b09b and run commit eb262bac2
Workflow/test infrastructure (diffstat): ci.yml, nightly.yml, wasm-threading.yml (+11), web.yml, check_wasm_threading_regressions.py (62 lines), wasi_thread_host.js (+42/-8). run_losat_wasi_threads.js, wasi_shared_memory.js, wasi_artifact.js, wasm_performance.py unchanged.
wasi_thread_host.js change: env now `{...process.env, PWD: cwd}`, and host-import plumbing for the Web adapter (WebAssembly.Module.imports scan at instantiate) - startup-only, no hot-path change.
Engine source: LOSAT/src/algorithm/tblastx changed in 4 files (blast_engine/mod.rs, blast_engine/run_impl.rs +329 lines, lookup/backbone.rs +247, lookup/mod.rs); commits touching tblastx: edf825d36 (Report outfmt 0 subject headings..., S05), 458cbe110 (Route BLASTN and TBLASTX through shared run_local entries, S04), 1c131e602 (Implement BLASTX ...), f9832b6ce (TBLASTN Stage C...). Also src/utils/threading.rs, common.rs, stats/*, core/* changed (88 files, +59169/-2506 in LOSAT/src overall). So an engine-side slowdown of threaded TBLASTX cannot be excluded from CI evidence alone -> task C.

## B. Existing local runs (wasm-threading-*)
These are the small-case `check_wasm_threading.py` gate matrices (409-423 records), NOT the frozen regressions: none of them contains a p11_avclpv_psclpv record (and no other local runs.json does either; grep over the whole target cache finds p11 only in the CI artifacts downloaded here). So there is no earlier local p11 timing.

| dir | metadata.json mtime | threaded wasm sha256 |
|---|---|---|
| wasm-threading-after | 2026-09-29 03:09 | 268f1329216d796c... |
| wasm-threading-final | 2026-09-29 03:51 | fae798b8f4d33e5608e... |
| wasm-threading-s03 | 2026-09-29 04:32 | 055ac71be4fa6... |
| wasm-threading-s04 | 2026-09-29 06:03 | 0a14aaabe3359... |
| wasm-threading-s05 | 2026-09-29 07:15 | 61f94c56e40c... |
| wasm-threading-s07 | 2026-09-29 14:31 | 03c5d5a2ce2b1... |
| wasm-threading-s07p | 2026-10-01 02:50 | 478d6fab8c700f13... (s07pg artifact = HEAD-equivalent) |

## Additional A observation: native p11 is serial, threaded p11 uses 4 threads
command.json for the CI records: native p11 runs with `-num_threads 1`, threaded p11 with `-num_threads 4` (RAYON_NUM_THREADS=1, LOSAT_WASI_THREADS_DEBUG=1, NODE_NO_WARNINGS=1, cwd = repo root). Threaded p11 in the green run consumed 6617 CPU-s vs 1248 CPU-s for the serial native run (5.3x), on ~2.0 cores average -> p11 is the one case where the threaded Wasm path burns far more CPU than native (other big cases: d01 1286 vs 1055, p06 1179 vs 299). Its wall time therefore depends strongly on how many vCPUs it actually gets while 2 other searches run (--jobs 3 on a 4-vCPU runner).
The timed-out stderr only shows the pool/stage banner (stage=linking work_items=4 parallel_selected=true), no progress.

## C. Local measurement
Machine: i9-14900HX, WSL2, 32 logical CPUs, node v26.8.2 (CI used node 24.21.0). Cores 28-31 are (probably) E-cores of the hybrid CPU - slower than CI vCPUs per core, but identical for both runs. Driver: drive_p11.py (imports authority.load_catalog/build_steps/stage_required_fixtures and wasm_performance.execute from the LOSAT-web-gui worktree, HEAD = 9c810a3d7; same env filtering, prefix [node, run_losat_wasi_threads.js, wasm], -num_threads 4, LOSAT_WASI_THREADS_DEBUG=1, LC_ALL=C, NODE_NO_WARNINGS=1, cwd = repo root, timeout 7200), run as `taskset -c 28-31 python3 drive_p11.py ...`.
Command check: local argv == CI command.json argv (node, run_losat_wasi_threads.js, wasm, tblastx -query /tmp/losat-pr5-runtime-cert-5845d22/LOSAT/tests/fasta/AvCLPV.fasta -subject .../PsCLPV.fasta -outfmt 6 -num_threads 4 -query_gencode 1 -db_gencode 1 -out ...) apart from the wasm path and output path. Env overrides: LOSAT_WASI_THREADS_DEBUG=1, NODE_NO_WARNINGS=1, RAYON_NUM_THREADS=1 (from step.environment) - same as CI.
Note: LOSAT-web-gui worktree has an uncommitted modification LOSAT/src/algorithm/blastn/blast_engine/run.rs (not mine; irrelevant, tests use prebuilt wasm).
Run 1 (HEAD-equivalent s07pg artifact 478d6fab...) started 2026-10-01T21:41:49Z (=06:41 JST), load avg before 16.65 (decaying; top showed ~97% idle at start).
Run 1 result (HEAD-equivalent 478d6fab...): PASS, wall 2163.6 s (36.1 min), cpu user 6327.4 s + sys 1.5 s, peak RSS 375.7 MB, output sha256 1eb11f5c... == expected (raw_equal true). Load average at start 16.65 (15 min 10.66), at end (22:20:59Z) 25.74 / 23.73 / 22.07 -> machine was heavily loaded by other work during the whole run (load up to ~30 on 32 logical CPUs; cores 28-31 shared with it). CPU/wall = 2.9 cores effective.
Host JS: old (green) artifact's run_losat_wasi_threads.js / wasi_shared_memory.js / wasi_artifact.js vs HEAD's: see next line.
- run_losat_wasi_threads.js identical
- wasi_thread_host.js DIFFERS (old host from green run vs HEAD host)
- wasi_shared_memory.js identical
- wasi_artifact.js identical
