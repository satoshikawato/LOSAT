# Pure-Rust runtime authority and production delegation inventory

## Record

```text
Audit date: 2026-08-30
Branch: pure-rust/runtime-boundary
BASE_SHA: 2844c19e3f3b5f2ed29dc88be1a69ff7af2277a0
NCBI behavioral source: official NCBI BLAST+ 2.17.0
```

This record implements
[`PD-PURE-RUST-RUNTIME-AUTHORITY`](../product_decisions/PD-PURE-RUST-RUNTIME-AUTHORITY.md).
It inventories project-authored production implementation and build boundaries
without changing BLAST behavior. Line numbers below describe `BASE_SHA`; the
checker identity is the exact rule/path/symbol tuple and therefore does not
depend on line numbers.

## Audited surfaces

The live audit covered:

```text
LOSAT/src/**/*.rs
LOSAT/Cargo.toml
LOSAT/Cargo.lock
LOSAT/.cargo/**
all project build.rs files outside tests, target, and local instruction data
.github/workflows/** for build/runtime invocation classification
resolved Cargo normal/build dependencies from cargo metadata --locked
```

The source scan classified imported and implemented `extern "C"`/`extern
"system"` functions, `#[link]`/`link_name`, direct `libc` and dynamic-loader
tokens, production `Command` use for NCBI executable names, `c_void`/raw C
types, native build tools, `qsort` declarations/calls, and production callers.
Rust comments are removed before matching so embedded NCBI snippets are not
mistaken for LOSAT implementation.

## Classification summary

| Class | Live result |
|---|---|
| Production algorithm delegation | Four native `qsort` import groups containing four declarations and five direct calls; nine exact temporary findings. |
| Production non-algorithm ABI boundary | Thirteen Rust-implemented web API exports in `web_api.rs` and five Rust-implemented qsort callbacks. These export/callback functions are not imported implementation owners. |
| Build-time native-code dependency | None in project build logic. The project has no `build.rs`, direct `cc`/`cmake`/`bindgen` dependency, native algorithm library, or project linker attribute. |
| Cargo `links` review signal | Two resolved normal-dependency findings, `rayon-core 1.13.0` and `wasm-bindgen-shared 0.2.121`; both reviewed as non-algorithm integration below. |
| Production subprocess/dynamic load | None. `LOSAT/src` contains no `Command` launch of `blastn`, `blastp`, `tblastx`, or `makeblastdb`, and no `libloading`, `dlopen`, `LoadLibrary`, or direct `libc::` call. |
| Test/diagnostic/certification oracle | NCBI executable use is outside production. The release-readiness workflow's package-manager NCBI jobs are explicitly named diagnostics and are not frozen release certification. |
| Documentation/comment/false positive | Program-name strings, NCBI C snippets, internal `NonNull<c_void>` storage, and Rust `extern "C" fn` exports/callbacks do not import native behavior. |

The live audit matches the revised authoritative inventory: four production
delegation groups, four imported C `qsort` declarations, five direct C calls,
and nine exact temporary project-owned findings. It found no additional
project-authored production algorithm delegation.

## Exact adapter and caller inventory

### Group B — BLASTN preliminary-HSP adapter

Owner and mechanism:

```text
LOSAT/src/algorithm/blastn/blast_engine/run.rs:3133-3375
native_prelim_qsort
qsort_prelim_hits_by
target cfg: native adapter under cfg(not(target_arch = "wasm32"));
            Rust fallback under cfg(target_arch = "wasm32")
native: clone PrelimHit values, qsort usize indices through a thread-local
        raw-pointer comparator context, then clone/replay complete values
serial Wasm: sort_unstable_by with the same comparator
```

Concrete comparator routes and timing:

| Comparator | Direct route | Timing and downstream consumer |
|---|---|---|
| `score_compare_prelim_hits` | `sort_prelim_hits_by_score` at 3312 | After combining preliminary lists (3533), after per-chunk common-endpoint purge (8905-8906), and before traceback interval-tree traversal (9214). Ordering controls merge/containment/traceback input order. |
| `query_offset_compare_prelim_hits` | `purge_prelim_hits_with_common_endpoints` at 3344 | First preliminary common-start purge; the first record in each context/query-start/subject-start group survives. |
| `query_end_compare_prelim_hits` | `purge_prelim_hits_with_common_endpoints` at 3361 | Second preliminary common-end purge; the first record in each context/query-end/subject-end group survives. |

There is no focused component test for this adapter. BLASTN program tests and
the certified manifest exercise the route, but the stable-equality seam belongs
to PR 2. PR 2 owns this adapter together with Group A because both feed the same
BLASTN survivor, traceback, and output pipeline.

### Group A — BLASTN final-HSP adapter

Owner and mechanism:

```text
LOSAT/src/algorithm/blastn/hsp.rs:473-658
qsort_blastn_hsps_by / native_hsp_qsort
target cfg: native adapter under cfg(not(target_arch = "wasm32"));
            Rust fallback under cfg(target_arch = "wasm32")
native: clone BlastnHsp values, qsort usize indices through a thread-local
        raw-pointer comparator context, then clone/replay complete values
serial Wasm: stable Rust slice sort_by with the same comparator
```

Concrete comparator routes and timing:

| Comparator | Callers | Downstream consumer |
|---|---|---|
| `score_compare_hsps` | `sort_hsps_by_score` (616), called from BLASTN engine seams including `blast_engine/mod.rs` (119, 254) and `blast_engine/run.rs` (4412, 10149) | Score-ordered HSP lists feed hit-list update, trimming, pruning, and final reporting. Existing already-sorted checks remain before the adapter. |
| `evalue_compare_hsps` | `sort_hsplist_by_evalue` (655) | E-value-first HSP-list ordering used by hit-list/output processing after its existing already-sorted check. |
| `compare_query_offset_hsps_for_common_endpoint` and `compare_query_end_hsps_for_common_endpoint` | `filtering/purge_endpoints.rs:qsort_hsp_array_prefix_by` (753-782), invoked by both common-endpoint purge passes | First-survivor common-start/common-end purge, gap-script preservation, re-evaluation, containment pruning, then score resort/output. |

Existing focused coverage includes
`common_endpoint_comparator_ignores_edit_script`,
`source_compatible_first_survivor_accepts_either_equal_order`, and
`source_compatible_purge_never_synthesizes_edit_script`. PR 2 adds only the
missing target-neutral stable-sort seam coverage. PR 2 owns this adapter
together with Group B.

### Group D — TBLASTX `UngappedHit` search/linking adapter

Owner and mechanism:

```text
LOSAT/src/algorithm/tblastx/ncbi_qsort.rs:1-163
qsort_ungapped_hits_by / native
target cfg: native adapter under cfg(not(target_arch = "wasm32"));
            Rust fallback under cfg(target_arch = "wasm32")
native: clone UngappedHit values, qsort usize indices through a thread-local
        raw-pointer comparator context, then clone/replay complete values
serial Wasm: stable Rust index sort followed by complete-record replay
```

Concrete comparator routes and timing:

| Comparator | Callers | Downstream consumer |
|---|---|---|
| `score_compare_ungapped_hits_ncbi` | `blast_gapalign.rs:ncbi_qsort_ungapped_hits_by_score` (362-385), reached from `blast_gapalign.rs`, `blast_engine/mod.rs`, and `blast_engine/run_impl.rs` | Score-ordered HSP lists before gap-align/list processing and combined-search reduction. |
| `rev_compare_hsps_tbx` | `sum_stats_linking/linking.rs:sort_hsps_by_ncbi_link_order` at 483 and 2278 | Reverse TBLASTX linking pass and test seam. |
| `rev_compare_hsps_transl` | Same owner at 656 | Reverse translated-query linking pass before chain construction/replay. |
| `fwd_compare_hsps_transl` | Same owner at 657 | Forward translated-query linking pass before chain construction/replay. |

The linking route must preserve translated context grouping and the complete
`UngappedHit`, including `link_id`, `chain_next_link_id`, `linked_set`,
`start_of_chain`, and `ordering_method`, through replay. Existing focused tests
include `index_replay_sort_uses_ncbi_comparator_without_extra_ties`,
`test_link_qsort_tie_preserves_current_hsp_list_order_for_equal_hits`,
`test_ncbi_score_sort_preserves_current_hsp_list_order_for_equal_hits`, and
chain-replay tests. PR 4 owns the remaining target-neutral wiring/replay gate.

### Group C — Shared `Hit` / `SubjectGroup` output adapter

Owner and mechanism:

```text
LOSAT/src/common.rs:459-722
qsort_hits_by / qsort_subject_groups_by / output_qsort_native
target cfg: native adapter under cfg(not(target_arch = "wasm32"));
            Rust fallback under cfg(target_arch = "wasm32")
native: move complete Hit or SubjectGroup records into Option slots, qsort
        usize indices through thread-local raw-pointer contexts, then replay
        each complete record once
serial Wasm: stable Rust index sort followed by complete-record replay
```

Concrete comparator routes and timing:

| Comparator | Caller | Downstream consumer |
|---|---|---|
| `score_compare_hsps` | `write_output_ncbi_order_with_format_to_writer_impl` at 1133 | Score-first per-subject HSP output. |
| `evalue_compare_hsps` | Same owner at 1143 | E-value-first per-subject HSP output; current concrete caller is TBLASTX. |
| `compare_subject_groups` | Same owner at 1179 | Subject-group order before pairwise/tabular formatting. |

The public APIs expose both score-order and E-value-order output paths. Current
in-repository TBLASTX callers in `blast_engine/run_impl.rs` use only the
E-value-order entry points; the score-order entry points currently have no
concrete in-repository engine caller. BLASTN uses its program-local HSP-list
formatter. That concrete reach does not make `common.rs` TBLASTX-owned: the
public APIs and record types are shared/generic output infrastructure. PR 3
owns this adapter as an independently reviewable shared-output migration.
There is no focused component test for the native owner; existing TBLASTX
output contracts exercise the E-value route at program level.

## Direct C `qsort` call accounting

| Group | Imported declaration | Direct call site(s) |
|---|---:|---|
| A — BLASTN final HSP | 1 | `hsp.rs:qsort@native_qsort_blastn_hsps_by` |
| B — BLASTN preliminary HSP | 1 | `blast_engine/run.rs:qsort@native_qsort_prelim_hits_by` |
| C — shared output | 1 | `common.rs:qsort@qsort_hit_indices`; `common.rs:qsort@qsort_subject_group_indices` |
| D — TBLASTX search/linking | 1 | `ncbi_qsort.rs:qsort@native_qsort_ungapped_hits_by` |
| **Total** | **4** | **5 direct calls** |

Together, the four declarations and five calls are the nine exact temporary
project-owned findings frozen in
`LOSAT/tests/pure_rust_runtime_allowlist.tsv`.

## Route-family map

```text
BLASTN preliminary score/common-endpoint ordering
  -> blast_engine/run.rs native_prelim_qsort
  -> merge, first-survivor purge, traceback interval tree

BLASTN final HSP score/evalue/common-endpoint ordering
  -> hsp.rs native_hsp_qsort
  -> hit-list pruning, endpoint purge, reporting

TBLASTX UngappedHit search HSP and reverse/forward translated-link ordering
  -> tblastx/ncbi_qsort.rs native
  -> gap-align list order, sum-statistics linking, chain replay/filtering

Shared Hit and SubjectGroup output ordering (current engine reach: TBLASTX)
  -> common.rs output_qsort_native
  -> pairwise/tabular formatting
```

## Cargo `links` review

`cargo metadata --format-version 1 --locked` reports two reachable
normal-dependency packages with non-null `links` metadata:

| Package | Dependency path | Review |
|---|---|---|
| `rayon-core 1.13.0`, `links = "rayon-core"` | `LOSAT -> rayon -> rayon-core` under the optional/default `parallel` feature | Its `build.rs` states that it is not linking anything; `links` enforces that only one `rayon-core` version is used. The script emits only `rerun-if-changed`. Reviewed non-algorithm Rust concurrency integration. |
| `wasm-bindgen-shared 0.2.121`, `links = "wasm_bindgen"` | Target-dependent Wasm routes through `getrandom`/`wasm-bindgen` and `indicatif`/`web-time` | Its build script hashes the Rust schema source and may query Git for a revision; it does not compile or link native algorithm code. Reviewed Wasm schema/version coordination. |

The exact package version, `links` value, and resolved dependency kinds are in
the separate machine-reviewed dependency inventory with classification
`reviewed_non_algorithm_integration`. They are not counted among the nine
temporary production findings. A change to any identity field becomes an
unexpected dependency signal and blocks automatic acceptance pending review.

## Mechanical ratchet

Files:

```text
LOSAT/tests/check_pure_rust_runtime_boundary.py
LOSAT/tests/pure_rust_runtime_allowlist.tsv
LOSAT/tests/pure_rust_runtime_dependency_review.tsv
LOSAT/tests/test_check_pure_rust_runtime_boundary.py
```

The exact finding key is:

```text
rule + repository-relative path + imported symbol/calling function
```

The production allowlist contains only the nine temporary project-owned
delegation findings. The separate dependency-review file contains only the two
reviewed Cargo metadata signals. The ratchet rejects a new rule or symbol in an
already allowlisted path, a new path, duplicate ambiguous keys, an
unreviewed/new/changed Cargo `links` identity, and every stale entry in either
inventory. It separately reports Rust-implemented ABI exports/callbacks as
non-importing observations. Test, diagnostic, documentation, generated target,
and local instruction data are outside the production scan.

Focused tests prove imported ABI/call detection, Rust export classification,
comment and explicit test-oracle exclusion, unknown finding rejection, stale
entry rejection, and Cargo `links` review semantics.

## Migration states

| Merged state | Temporary production findings | Reviewed metadata baselines |
|---|---:|---:|
| PR 1 | 9: both BLASTN groups, shared output, and TBLASTX search/linking | 2 |
| After PR 2 | 5: shared output (3) plus TBLASTX search/linking (2) | Re-evaluate exact identities; normally 2 |
| After PR 3 | 2: TBLASTX search/linking | Retain only exact reviewed non-algorithm identities |
| After PR 4 | 0, enforced in ordinary CI | New/changed identities require review |
| PR 5 onward | 0 | New/changed identities require review |

PR 2 changes only both BLASTN adapter groups, PR 3 changes only shared output,
and PR 4 changes only TBLASTX search/linking. This PR does not change any
comparator, sort invocation, grouping, purge, linking, replay, filtering,
formatting, manifest, fixture, expected output, or classifier.
