# NCBI's 1-thread output on a protein database with a sequence of length 0

`<search>.dba0.out` and `.err`: the standard output and standard error of NCBI BLAST+ 2.17.0 (the official Linux x86-64 build),
`-num_threads 1 -outfmt "6 std slen"`, on `dba0` of `LOSAT/tests/protein_db_32bit_lengths.rs` (volume `db0`, one
sequence of length 0, then volume `dbs`, one sequence), frozen on 2026-11-27 (line RG-D2) from the probe
database `a_0s` of `docs/evidence/zero_length_protein_mt/` (the same sequences and queries):

| File | Search |
|---|---|
| `blastp.dba0.out` | `blastp -query q.fa` |
| `blastp_ug.dba0.out` | `blastp -query q.fa -ungapped -comp_based_stats 0` |
| `blastx.dba0.out` | `blastx -query qn.fa` |
| `blastx_ug.dba0.out` | `blastx -query qn.fa -ungapped -comp_based_stats 0` |
| `psiblast.dba0.out` | `psiblast -query q.fa -num_iterations 2` |

At 2 or more threads NCBI aborts on this database (PD-LOSAT-NCBI-DEFECTS exception 36); the test
`a_protein_of_length_0_keeps_the_one_thread_output_at_any_thread_count` checks that LOSAT prints these
bytes (both streams) at 1, 2 and 4 threads. `docs/evidence/zero_length_protein_mt/freeze.sh` froze them.
