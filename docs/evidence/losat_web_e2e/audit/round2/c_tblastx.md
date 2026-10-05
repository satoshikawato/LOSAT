# Round 2, angle (c) TBLASTX: auditor's report

> S08+b. The final reply of the Sonnet auditor (read-only), as returned (the harness refused its report file). Binary `d8d18ec0…cb116` (the gate native of `491292327`). Work dir `~/.cache/losat-web-gui-target/s08pb-audit/c/` (`out/hA.txt`, `hB.txt`, `hC.txt`, `misc.txt`, `tk.txt`, `thr.tsv`, `big.tsv`, `scal.tsv`, `gc.tsv`, `gc2.tsv`, `cul.tsv`, `cul2.tsv`, `cul3.tsv`, `bnd.txt`; raw outputs in `r/`, `t/`, `g/`). Brief: [`brief/COMMON.md`](brief/COMMON.md), [`brief/ANGLE_C.md`](brief/ANGLE_C.md). How the findings were handled: [`../ROUND2.md`](../ROUND2.md). The auditor stopped only its own `thr.sh` run (NCBI `-word_size 4 -threshold +inf` used 7 GB per run).

## Overall verdict: SUPPORTED

## Round 1 findings, re-run

| ID | Re-run verdict | Evidence |
|---|---|---|
| TX-1 | LOSAT-REJECTS as decided (D11) | 8 values where NCBI ends with no hits; LOSAT exit 1 with the D11 text |
| TX-2 | LOSAT-REJECTS as decided (D8) | 14 values on q1/s1 and q2/s3; NCBI not run (it does not end) |
| TX-3 | LOSAT-REJECTS as decided (D10) | regular-file stdin and pipes; `cat qm.fna | … -query - -subject sm.fna` SAME (470,951 bytes) |
| TX-4 | SAME (26) and LOSAT-REJECTS (14) | malformed locut/hicut byte-identical "Invalid input for filtering parameters"; values NCBI runs (`0x10 0X1p3 +inf -infinity -nan 1e400 1e99999`) explicitly rejected |
| TX-5 | ACCEPTED | 12 runs, both exit 1 |
| TX-6 | LOSAT-REJECTS with the phrase / ACCEPTED | thread limits; 2, 4, 200, 1000: stdout SAME, only NCBI's thread warnings differ |
| TX-7 | FIXED | `protein_options.rs:492` cites `blast_options.c:1518-1521` |
| TX-8 | ACCEPTED, not run | ABI v1 frozen |
| TX-9 | ACCEPTED; CLI side SAME | `Query is Empty!` byte-identical |
| TX-10 | FIXED | AUTHORITY.md lines 71, 115, 123 agree |

Round 1 harness (730 runs): SAME 319, ACCEPTED (thread warnings) 28, ACCEPTED (parser, exit 2 vs 1) 254, LOSAT-REJECTS 127, DIFF 2 (`-help`, `-help 1`: approved), TIMEOUT 0. The classes per id are identical to S08+a's rerun. The S08+a CLI changes (`ncbi_preparsed_toolkit_word`, `unknown_option_error`) changed no accepted TBLASTX argv (706 round 1 argv, 50 further `=`-syntax probes).

## New checks

- Threshold extremes (71 runs, `/usr/bin/time -v`): `+inf 1e300 2147483647 2147483648 4294967296 0.5 1`, word size 3 on four input pairs, outfmt 0 and 6: 53 SAME, 18 LOSAT-REJECTS (word size 2 and 4), 0 DIFF, no signal 9. Values that wrap to INT_MIN (D4) make every word a neighbour; time is within 1.0-1.3× NCBI for every run of 10 s or more; memory below (R2c-1).

| Query (`-threshold +inf`) | NCBI | LOSAT | RSS ratio |
|---|---|---|---|
| 5 kb × s1 | 961 MB, 8.2 s | 1,205 MB, 9.1 s | 1.25 |
| 10 kb | 1,811 MB, 15.7 s | 2,398 MB, 16.3 s | 1.32 |
| 20 kb | 2,509 MB, 23.1 s | 4,862 MB, 27.7 s | 1.94 |
| 30 kb | 3,817 MB, 37.7 s | 8,726 MB, 45.2 s | 2.29 |
| 20 kb × 30 kb subject | 2,510 MB, 106.6 s | 4,862 MB, 104.2 s | 1.94 |
| qbig × sbig | 2,051 MB, 58-72 s | 2,638 MB, 75-78 s | 1.29 |

- `-db_gencode`/`-query_gencode` against NCBI with the subject as a database (`makeblastdb -dbtype nucl`, `-db`): 288 + 63 combinations SAME (outfmt 6 byte-identical; outfmt 0 identical outside the database lines). Invalid codes are parse errors on both sides.
- `-culling_limit`: 354 runs SAME (1, 2, 3, 10, 0x2 × `-max_target_seqs` 1-5 × `-evalue` 10/1000; with genetic codes against the `-db` oracle; with `-seg`, `-threshold`, `-window_size`, `-evalue` crosses).
- Toolkit words in a value position (104 runs): 79 explicit rejections, 13 SAME (NCBI also reads a plain value), 10 parser class, 2 DIFF (`-help`; `-evalue 1 --`: R2c-2).
- D8/D11 boundary: at `window = 2^30 - L` NCBI finishes (about 5.3 GB) and LOSAT gives the same stdout; at +1 NCBI does not finish and LOSAT rejects (q1/s1, qm/s1, qbig/sbig, gq/gsub; the sum is exact for multi-record batches).

## New findings

**R2c-1 (medium; resource use only, no output difference).** LOSAT TBLASTX uses up to 2.3× NCBI's memory when the neighbourhood is complete (`-threshold` values that wrap to INT_MIN). NCBI `blast_aalookup.c:413,472,546-603`; LOSAT `tblastx/lookup/backbone.rs:161-162,1018,1108,1205` (per-word chains, then the overflow array). Repro: `/usr/bin/time -v … -query c/bq30000.fna -subject c/s1.fna -word_size 3 -threshold +inf -outfmt 6`: same output, NCBI 3,817 MB / 37.7 s, LOSAT 8,726 MB / 45.2 s. Explains the signal 9 of the earlier aborted gate run (a 60 kb query would need about 17 GB in LOSAT).

**R2c-2 (low; argument syntax).** A last `--` is accepted by NCBI and rejected by LOSAT's parser: `tblastx -query q2.fna -subject s3.fna -evalue 1 --` (NCBI exit 0, 44,081 bytes; LOSAT exit 2 `error: unknown option or argument '--'`). NCBI `ncbiargs.cpp:89,2867-2868` (`--` ends the options), `ncbiapp.cpp:938-941`. With words after `--`, NCBI fails with USAGE exit 1 (accepted parser class).

**R2c-3 (low; wording only, justified rejection).** The TX-3 rejection says "stream without a position (such as a pipe)" also when standard input is a regular file (`$FINAL tblastx -subject - < c/qm.fna`; NCBI's `tellg()` fails at EOF, `blast_app_util.cpp:856-860`). D10.

## Not run

Web ABI v1 (TX-8) and the adapter's `validate` (TX-9).
