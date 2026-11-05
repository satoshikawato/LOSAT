# pylosat on Python gbdraw (stage 3 of the library plan)

The gate of stage 3 (`docs/losat_library_surface_plan.md` §10.1, §12): gbdraw runs without a
subprocess. gbdraw's commit is `033a3aa`; its repository is not changed.

- `gbdraw_pylosat.patch`: the local patch. `run_losat_search` searches through `pylosat` when it
  can be imported, with the words of `build_losat_command` and `outfmt="6"`. gbdraw's order of
  runtimes, its records and its tests of that order are gbdraw's product decision, so this patch is
  a check and not a proposal.
- `compare_runtimes.py`: runs gbdraw's entry points twice per scenario, through `pylosat` with
  `subprocess.Popen` refused, and through the `losat` executable of the same LOSAT commit
  (`--losat_bin`), and compares the SVG files and the raw LOSAT entries of the sessions.

## Commands

```bash
git clone https://github.com/satoshikawato/gbdraw && cd gbdraw && git checkout 033a3aa
git apply gbdraw_pylosat.patch
pip install -e . && (cd <LOSAT>/python && maturin develop)
LOSAT_BIN=<losat> PYTHONPATH=. python compare_runtimes.py <empty directory>
```

## Results

Dev builds on the machine of the stage 3 session, Rust 1.97.0, Python 3.13.18, LOSAT at the head of
`feature/pylosat`:

| Scenario | SVG | Raw searches |
|---|---|---|
| `linear-losatn` (lambda, DE3) | identical | 1 identical |
| `linear-losatn-blastn-task` | identical | 1 identical |
| `linear-tlosatx` | identical | 1 identical |
| `linear-tlosatx-code-11` | identical | 1 identical |
| `linear-losatp` | identical | 4 identical |
| `circular-tlosatx` (human mtDNA, codes 2 and 5) | identical | 3 identical |
| `circular-losatn` | identical | 3 identical |

gbdraw's own tests with a real search, through `pylosat` with process starts refused:

- Pass: `test_cli_losatn_reproduces_the_t_cli_07_figure` (the frozen SVG, byte for byte),
  `test_losatn_ring_is_identical_for_fasta_genbank_and_ddbj`,
  `test_native_losatn_reproduces_lambda_de3_bytes`,
  `test_native_tlosatx_reproduces_lambda_de3_rows_for_any_thread_count`.
- `test_native_losat_records_its_runtime_identity` fails: it expects 0.1.0, the release gbdraw pins;
  the engine here is 0.2.0.
- `test_cli_tlosatx_reproduces_the_t_cli_09_figure` and the three
  `test_native_tlosatx_reproduces_metazoan_mtdna_bytes` fail, and so does the executable of the same
  commit. Those tables were made by LOSAT v0.1.0: `losat` built at the tag reproduces them byte for
  byte. Its tables of genetic codes 2, 3, 5, 13 and 21 have `ATA` as I where NCBI's `gc.prt` has M,
  and code 31's has `TGA` as a stop where `gc.prt` has W. The 27 tables of the current engine
  (`LOSAT/src/utils/genetic_code.rs`) equal `gc.prt`. Rows across the `N` at position 3107 of the
  human sequence differ as well, also with code 1: the current engine draws ambiguous subject letters
  as NCBI's preliminary stage does (`tblastx/subjects.rs`); which commit changed them was not traced.
  gbdraw's tutorial rows for these inputs change when it moves to a newer LOSAT.
