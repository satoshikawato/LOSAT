# v0.2.0 release test set J: memory and speed

Set J of [the release test plan](../../../release/v0.2.0_certification_plan.md) (section 2) on the
candidate commit `C` = `93ca2f652617acac09d416c28c293cc8a1a501d1`, measured on the Linux x86_64
workstation (dual Intel Xeon Gold 6338, 2 NUMA nodes) on 2026-11-19 and 2026-11-20. NCBI BLAST+ 2.17.0
(official Linux x64 archive) ran on the same machine in the same windows. LOSAT was built with
`cargo build --release --locked` and Rust 1.92.0 (binary SHA-256 prefix `7343715b5e551b28`).

## Memory: passed

1,954 runs, none failed or hit a cap. Peak is the peak `RssAnon` read by a sampler (50 ms, 2 ms or
0.1 ms as named in each file); NCBI is the same command line of NCBI BLAST+ 2.17.0.

| Part | Rule | Result | File |
| --- | --- | --- | --- |
| Standard table, 1 to 16 threads | max RSS at or below NCBI on every cell | 0 of 99 cells above NCBI (highest: tblastx `-subject` 16 threads, 0.64 of NCBI) | `memory_standard_n1-16.md` |
| 24 to 64 threads | at or below NCBI | 0 of 95 cells above NCBI (highest: blastx 800 queries x 200 subjects `-subject` 64 threads, 0.986 of NCBI over 8 rounds) | `memory_n24-64.md`, `memory_bx800s200_n64_8rounds.md` |
| Memory floor rows | at most 1.05 x the recorded floor | 0 of 25 rows above (tblastx Sakai x MG1655 `-subject` 16 threads: 1,985 MiB against a floor of 1,990 MiB; the 10 ms trivial run is read against its floor scaled by the reference build's peak re-measured at 0.1 ms, 1.016) | `memory_floor_50ms.md`, `memory_floor_2ms_8rounds.md`, `memory_triv_0.1ms_16rounds.md` |

The tblastx 16-thread row failed on the earlier candidate `703caed17` (2,226 MiB); #1338 fixed it, which is
why `C` replaced that candidate. Output: every cell gives the same bytes as the earlier candidate; it differs
from NCBI only in the tblastn and tblastx `-subject` cells covered by PD-TLOSAN-LOCAL-GENCODE-32 and
PD-TBLASTX-LOCAL-SUBJECT-GENCODE.

## Speed: passed

Nine search modes (`speed_inputs.tsv` names every input and its SHA-256). Arms: `cp` = `C`; `d48` = the dev
build `6a0a8719a` (SHA-256 prefix `e05f579903a000f3`) as the previous reference; `v010` = the v0.1.0 Linux
x64 release binary (it cannot run blastp-fast, dc-megablast, blastn-short, tblastn and blastx, so those
cells are empty); `ncbi` = NCBI BLAST+ 2.17.0 with `-db` against a database it built from the same FASTA.

- 16 threads: wall time, 16 physical cores of NUMA node 1 booked for timing, 8 rounds after one discarded
  round, arm order rotated each round.
- 1 thread: `cycles:u` and `instructions:u`, median of 5 alternated runs for `cp` and `d48`, one run for
  `v010` and `ncbi`, node 0.
- LOSAT runs each mode with `-subject` (`cp`, `d48`, `v010`); `ncbi` runs the same search with `-db` (the
  plan's "NCBI `-db` same conditions"), and for blastp and blastp-fast `ncbis` runs NCBI with LOSAT's own
  `-subject` command line once.

| Rule | Result |
| --- | --- |
| `cp` at most 1.05 x `d48` | met in all 9 modes at 16 threads (0.903-1.014) and at 1 thread (cycles 0.951-1.011) |
| `cp` faster than NCBI | met in all 9 modes at 16 threads (0.0016-0.420 of NCBI's wall time) and at 1 thread (0.010-0.95 of NCBI's cycles) |
| `cp` against v0.1.0 | 0.022-0.182 of v0.1.0's wall time at 16 threads in the 4 modes v0.1.0 can run |

Every arm wrote the same output in every run (`speed_grid.md`, the sha16 columns). `cp` and `d48` give the
same bytes in every mode; v0.1.0 differs in megablast and tblastx (fixed since). NCBI's `-db` output equals
LOSAT's `-subject` output in 7 modes; in blastp and blastp-fast the two forms report different statistics
(`-subject` searches each subject on its own, `-db` the whole database as one search space), and NCBI's own `-subject` run (`ncbis`) gives LOSAT's bytes.
