# TBLASTN: the share of subject-only preprocessing (native callgrind, S09 DW-8)

Binary: the native LOSAT of commit `dfa65cfd8` (release, symbols, no line tables), one
thread, run from the repository root by `run_case.sh <query> <subject>` (outfmt 6). The
callgrind outputs are `cg_<query>_<subject>.out`; `ann_incl_*.txt` and `ann_excl_*.txt` are
`callgrind_annotate --inclusive=yes|no` of them.

Inputs:

- full query: `docs/evidence/tlosan_stage_g/benchmark/first_AvCLPV_protein.faa` (1 record, 1046 aa; the V-PERF TBLASTN fixture)
- qmin: `qmin.faa`, the same header and the first 30 residues
- full subject: `LOSAT/tests/fasta/AvCLPV.fasta` (416,069 nt); half subject: its header and first 208,034 nt (70 per line)

## Instruction counts (Ir)

| case | Ir | wall, median of 25 interleaved runs (`timings_ms_v2.txt`) |
|---|---:|---:|
| full query, full subject | 425,472,929 | 82.5 ms |
| qmin, full subject | 219,527,350 | 47.7 ms |
| full query, half subject | 317,978,243 | 56.6 ms |
| qmin, half subject | 121,090,295 | 27.8 ms |

## Roles (full query, full subject; inclusive Ir)

| role | Ir | share | functions |
|---|---:|---:|---|
| (a) subject only | 173,852,222 | 40.9% | `tblastx::translation::generate_frames` from the per-subject closure (153.5 M, 1 call: 6 frames of the whole subject; 132.7 M of it in `GeneticCode::get`, 832,134 codons, about 159 Ir per codon); `resolve_local_subject_ncbi2na` (15.4 M) and `resolve_ncbi4na_to_ncbi2na` (5.0 M) |
| (b) query only | 45,458,675 | 10.7% | `build_ncbi_lookup_for_profile` (34.0 M, built twice), `SegMasker::mask_sequence` / query encoding (11.5 M) |
| (c) scan | ≤ 39,306,403 | ≤ 9.2% | inlined into the per-subject closure (self cost 33.6 M) and seed vector growth |
| (d) extension, composition-based statistics, traceback, linking | 157,670,147 | 37.1% | `kappa::redo_preliminary_match` (116.2 M), `find_protein_init_hsps_by_chunk_with_mask_mode` (36.4 M), hit-window re-translation (4.4 M) |
| (e)+(f) output, FASTA reading, copies, start-up, other | about 9.2 M | 2.2% | |

- Lower bound of (a): the single `generate_frames` call, 36.1%. Upper bound, counting the
  whole closure self cost and the FASTA read as subject-only: about 50%. The wall-time
  scaling with the subject length (qmin: 47.7 ms full, 27.8 ms half) agrees (about 42%).
- For qmin, (a) is 79.2% of the Ir. (a) is exactly proportional to the subject length
  (86.9 M Ir for the half subject, 418 Ir per nucleotide).
- Where: `LOSAT/src/algorithm/tblastn/search_seed.rs:428` (`resolve_local_subject_ncbi2na`,
  defined at `:138-155`) and `:448` (`generate_frames`, `LOSAT/src/algorithm/tblastx/translation.rs:98-126`;
  `translate_sequence` `:68-88`; `GeneticCode::get`, `LOSAT/src/utils/genetic_code.rs:188-231`).
  Once per subject per search; nothing is cached across searches.
- Sensitivity (estimate, not measured): a 64-entry codon table per genetic code with the
  current function for ambiguous codons (about 10 Ir per codon) would leave (a) at about
  50 M Ir, about 16% of the full case.
