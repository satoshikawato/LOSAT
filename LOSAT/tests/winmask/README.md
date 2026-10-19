# WindowMasker fixtures

Written by NCBI BLAST+ 2.17.0 `windowmasker` (build Aug 11 2025, conda `blast 2.17.0`) on the
benchmark workstation, 2026-10-19. The tests in `LOSAT/src/utils/winmask/tests.rs` mask each
FASTA record of `manifest.tsv` with the counts file and the flags of its row and compare the
`-outfmt interval` text with `expected/`.

## Counts

| File | Source | Command |
| --- | --- | --- |
| `counts/mg1655_u9.obinary` | `LOSAT/tests/fasta/MG1655.fna` | `windowmasker -mk_counts -in MG1655.fna -unit 9 -sformat obinary -out mg1655_u9.obinary` |
| `counts/dm2L_u9.binary` | *D. melanogaster* BDGP6.54 2L:1-3,000,000 (Ensembl `dna.toplevel`) | `windowmasker -mk_counts -in dm2L_3M.fa -unit 9 -sformat binary -out dm2L_u9.binary` |

`counts/dm_u2.obinary` (`-mk_counts -in dm2L_60k.fa -unit 2`) has a cache bit array of 0 words;
NCBI keeps it in use, and the boost would read past it (its thresholds of 500 keep the boost
from running on `dm2L_3k.fa`).

Units of 9 letters keep the files small: the default for these inputs is 10, whose optimized
binary file holds a hash table of 2 MiB. Units of 7 and 8 give thresholds stuck at the count
limit of 500. Both files have thresholds (t_low, t_extend, t_threshold, t_high) of
`[74, 147, 173, 210]` and `[43, 106, 136, 184]`.

## FASTA

- `fasta/dm2L_60k.fa`: 2L:1-60,000; `fasta/dm2L_3k.fa`: its first 3,000 bytes.
- `fasta/dm_edge.fa`: `python3 make_edge.py dm2L_1_200000.fa dm_edge.fa` (2L:1-200,000): 20 kb
  pieces with a lower-case stretch, an N run and single ambiguity letters, an N every 13
  letters, records of 1, 9, 10, 13, 14, 15 and 30 letters, an all-N record and a tail of N.
- The other rows use `LOSAT/tests/fasta/{NZ_CP006932,AP027131,PsCLPV}.fasta`.

## Expected output

`windowmasker -ustat <counts> -in <fasta> -outfmt interval [flags] -out <expected>`.

- The `p_*` and `dmp_*` rows set the masker's arguments: `-set_t_low 200` (above T_extend, so
  the cache bit array skips windows that score above it), `-t_extend 50`, `-window 9` (one
  unit), `-window 255`, and `-window 256` (`Uint1`, so the default).
- `expected/mg1655_u9.p_low50.stderr`: the warning of `-t_low 50` (below the file's 74).
- `expected/mg1655_u9.p_window5.stderr`: the error of `-window 5` (below the unit of 9).

## Checked against blastn

`blastn -window_masker_db <counts>` builds the masker with the same settings as the
application (`winmask_filter.cpp:99-154`, `win_mask_config.cpp:318-361`). On 2026-10-19 the
query masks of blastn's archive (`-outfmt 11`) equaled the application's intervals for the four
tasks, `-dust yes|no`, the obinary and binary formats and `-lcase_masking`, with counts of
units of 10 from 3 Mb of 2L and of human chromosome 21.

## Larger files (not in the repository)

The ignored test `full_size_counts_match_ncbi` reads a manifest of the same columns from
`LOSAT_WINMASK_FULL_CHECK`. On the workstation it lists counts of whole 2L (unit 12), the whole
*D. melanogaster* genome (unit 13) and human chromosome 21 (units 12, 15 and 16; obinary and
binary), with whole 2L and chromosome 21 as queries.
