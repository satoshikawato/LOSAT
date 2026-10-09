# SFb step 6 sweeps: scripts, frozen NCBI side, before counts (2026-10-08)

Scripts `fasta_sweep.py` and `check_inputs.py` in this directory (both
scripts import each other from their own directory, and `check_inputs.py` finds E2g's script in
`docs/evidence/losat_web_e2g/` next to it, or under `$WT/docs/evidence` when run from here):

- `fasta_sweep.py`: seeded FASTA mutation sweep (also the shared library: runner, classes, summary).
- `check_inputs.py`: E2g's `check_inputs.py` with the E2h expectations plus TBLASTX/TBLASTN/BLASTP input cases.
- `before-counts.md`: class counts of the S11 build (before any port) per program x role x outfmt.

Outputs: `$BUILD_ROOT/sfb-e2h/sweeps/` (`fasta_sweep/` inputs and `cases.tsv`, `work/check_inputs/`,
`ncbi/{fasta_sweep,check_inputs}/` frozen NCBI, `ncbi-run2/` the determinism run, `before/*.tsv`).
Logs: `$TASK_DIR/logs/sweeps/`. NCBI BLAST+ 2.17.0 (`$NCBI_BIN`), run inside `unshare -rn`, from the
sweep directory (relative paths), without `~/.ncbirc`, with `BATCH_SIZE` and the other report variables unset.

## Reproduce

```bash
source /home/kawato/losat-baselines/sfb-e2h-20261008/env.sh
S=$WT/docs/evidence/losat_web_e2h; R=$BUILD_ROOT/sfb-e2h/sweeps; export PYTHONDONTWRITEBYTECODE=1; cd $R
# 1. FASTA mutation sweep (2280 cases)
python3 $S/fasta_sweep.py generate --dir $R/fasta_sweep                       # deterministic, seconds
flock $BUILD_ROOT/oracle.lock python3 $S/fasta_sweep.py freeze-ncbi --dir $R/fasta_sweep --ncbi-bin $NCBI_BIN --jobs 3 --out $R/ncbi/fasta_sweep
flock $BUILD_ROOT/oracle.lock python3 $S/fasta_sweep.py check --dir $R/fasta_sweep --ncbi $R/ncbi/fasta_sweep --losat <LOSAT> --jobs 3 --out <TSV>
# 2. check_inputs (1032 cases); --work must be the same path at freeze and check (it is inside the paths of the cases)
flock $BUILD_ROOT/oracle.lock python3 $S/check_inputs.py freeze-ncbi --ncbi-bin $NCBI_BIN --work $R/work/check_inputs --out $R/ncbi/check_inputs --jobs 3
flock $BUILD_ROOT/oracle.lock python3 $S/check_inputs.py check --losat <LOSAT> --work $R/work/check_inputs --ncbi $R/ncbi/check_inputs --out <TSV> --jobs 3
# 3. tables
python3 $S/fasta_sweep.py summary <TSV> [<TSV> ...]
# 4. determinism of the frozen side
python3 $S/fasta_sweep.py compare-frozen $R/ncbi/fasta_sweep $R/ncbi-run2/fasta_sweep
```

`<LOSAT>` is the binary that has the port (the before run used `~/.cache/losat-web-gui-target/sf/bin/LOSAT-native`).
Both `check` commands exit 1 when a case is `differs`, `timeout`, an unlisted rejection, or (check_inputs)
not the expected class. Take the load average under 12 first (`freeze` and `check` wait by themselves).

## Cases and run times

| Script | Cases | NCBI freeze (jobs 3) | LOSAT check, before binary |
| --- | --- | --- | --- |
| `fasta_sweep.py` | 2280 = 456 files (114 per molecule x role) x programs x outfmt 0, 6 | 97.5 s, 80.8 s (second run) | 3.1 s |
| `check_inputs.py` | 1032 = E2g's 300 BLASTN rows + 244 each for TBLASTX, TBLASTN, BLASTP | 42.8 s, 41.7 s | 2.6 s |

Determinism: the two NCBI runs of each script give identical `manifest.tsv` (rc, stdout, stderr and `2>&1`
hashes): fasta_sweep `152104448d65aca4...`, check_inputs `0678ecaf172083cab...`; `inputs.sha256` (hash of every
generated input file and every engine file the rows read) is identical too. NCBI exit codes of the fasta_sweep:
0 x 2166, 1 x 58, 3 x 10 (all records empty), 255 x 36 (Seq-id subjects, unshare), -11 x 10 (punctuation titles).

fasta_sweep matrix: program (blastn default task = megablast, blastn `-task blastn`, tblastx, tblastn
protein query / nucleotide subject, blastp) x role (query, subject: the other role is a clean file) x
outfmt 0, 6. Families per file: `eol` (LF, CRLF, CR only, mixed, no final newline, blank lines at the end, lone CR/LF
inside a line), `structure` (empty records, headerless first record, blank/comment lines and BOM before the first
defline, `>?` gap lines, `>?_x`, very long titles 999-40000 bytes, punctuation titles, HTML entities, 20-letter
titles), `defline-insert`, `defline-replace`, `seq-insert`, `seq-replace`, `seq-run` (control bytes, tab, CR, NUL,
`>`, `?`, `;`, `-`, `*`, digits, spaces, NBSP, 0x80/0xE9/0xFF, UTF-8, BOM, residues of the molecule), `combo` (a
defline and a sequence mutation and a line-end kind), `seqid` (6 files per molecule x role whose first line NCBI
may try as a Seq-id; NCBI is run in `unshare -rn` and cannot reach the network: 24 files, 120 cases).

## Classes

| Class | Meaning |
| --- | --- |
| `same` | NCBI succeeds; LOSAT has the same exit code, stdout, stderr and `2>&1` stream |
| `same-error` | NCBI fails; LOSAT has the same exit code and NCBI's stdout, stderr and `2>&1` |
| `explicit-rejection` | LOSAT exits non-zero, stderr contains `not supported by LOSAT's <PROGRAM>`, and the message is one AUTHORITY.md section J lists (`LISTED_REJECTIONS` in `fasta_sweep.py`: Seq-id line, `-parse_deflines`, record over 2147483647 letters, non-UTF-8 `Subject_` title). In check_inputs also the option rejections of E2g (the case expects it) |
| `explicit-rejection/unlisted` | the same, with a message that J does not list: a failure of the port (turn it into `same` or add it to AUTHORITY.md) |
| `approved-exception` | NCBI dies of a signal and LOSAT exits 0 (approved exception 2 of PD-LOSAT-NCBI-DEFECTS, punctuation title; column `approval` names it), or both argument parsers reject (approved exception 1 of PD-LOSAT-CLI-NONSEARCH-DIFFERENCES; check_inputs only) |
| `differs` | anything else, with the first differing line in `detail` |
| `timeout` | a run took over 120 s |

check_inputs `expect` values: `same`, `same-error`, `explicit-rejection` and `approved-exception` as E2g (its
`losat-rejects` = `explicit-rejection`, `arg-error` and `exception-2` = `approved-exception`); `fixed` = the 31
E2g `losat-rejects` rows about reading the FASTA file (set `FIXED` in the script), and `auto` = the new rows: the
expected class is decided by the frozen NCBI exit code (0: `same`, signal: `approved-exception`, otherwise
`same-error`), so a port that reads a file differently from NCBI cannot pass by expectation.

## Not covered / limits

- `approved-exception` for punctuation titles is decided by the NCBI signal and LOSAT exit 0 only; the content of
  LOSAT's report is checked by E2g's `title_sweep.py` (stand-in subject), not here.
- The inputs of check_inputs are read from the engine worktree (`LOSAT/tests/fasta/`) and E2c's `inputs/`;
  `check` compares their hashes with `inputs.sha256` and warns when they changed.
- Not covered: `-lcase_masking`, `-num_threads`, `BATCH_SIZE` batches of the new programs, queries over one batch.
  They are in `LOSAT/tests/fasta_input_fixtures.py`.
