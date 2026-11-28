# A database sequence of length 0 at more than one thread (PD-LOSAT-NCBI-DEFECTS exception 36)

Line RG-D2, 2026-11-27, found in review R4 of PR #1660. NCBI BLAST+ 2.17.0 (the official Linux x86-64
build) and LOSAT (commit da0ec9442 of PR #1660), on ws (2 x Xeon, Linux).

## What the runs show

`probe.py` builds each database with NCBI's `makeblastdb`, rewrites its index and sequence file so that
one sequence has length 0, and runs each search at `-num_threads` 1, 2, 4 and 8 (3 runs above 1 thread),
`-outfmt "6 std slen"`. `results.tsv` has one row per run: tool, program, database, threads, run, exit
status (Python's: `-6` is SIGABRT, the shell's 134), rows, SHA-256 of the output (8 digits), first line
of the standard error.

| Databases | Searches | NCBI at 1 thread | NCBI at 2, 4, 8 threads | LOSAT at 1 to 8 threads |
|---|---|---|---|---|
| protein: `s_first`, `s_mid`, `s_last` (one volume of 2 or 3), `m_first`, `m_mid` (one volume of 300), `a_0s`, `a_s0`, `a_s0s` (alias of 2 or 3 volumes, the sequence of length 0 in the first, last or middle one) | `blastp`, `blastp -ungapped -comp_based_stats 0`, `blastx`, `psiblast -num_iterations 2` | exit 0 | SIGABRT in all 288 runs: `BLAST engine error: Database memory map file error` (`psiblast`: `BLAST Database error: Sequence not returned.`), no output | exit 0, NCBI's 1-thread bytes |
| nucleotide: `na_0s` (alias of 2 volumes), `ns_mid` (one volume of 3) | `blastn`, `blastn -task blastn`, `tblastn`, `tblastx` | exit 0 | SIGABRT in all 72 runs: `BLAST engine error: Database memory map file error`, no output | exit 0, NCBI's 1-thread bytes |

Exception 36 (Owner 2026-11-27 09:25 JST) covers the protein rows. The nucleotide rows are the same
defect of NCBI's per-thread buffer; they wait for the Owner's answer and are not part of exception 36.

`a_0s` is the test database `dba0` of `LOSAT/tests/protein_db_32bit_lengths.rs` (the same sequences and
queries); `LOSAT/tests/fixtures/zero_length_mt/` holds NCBI's 1-thread output of its five searches.

```sh
python3 -I docs/evidence/zero_length_protein_mt/probe.py <NCBI BLAST+ 2.17.0 bin> <work dir> [losat]
```

With `SKIP_NCBI=1` and a `losat` path it runs LOSAT only. `freeze.sh` wrote NCBI's 1-thread outputs on the same
databases (2026-11-27), among them the fixture files (the `blastx -ungapped` one has no row in `results.tsv`). Every run is capped
(`prlimit --as=8589934592 --core=1 timeout 120`).
