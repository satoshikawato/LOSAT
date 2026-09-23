# TLOSAN v0.2.0 Stage C continuation — still incomplete

Branch: `feature/tlosan-tblastn-v0.2.0`; starting commit:
`ad3e8381fbdda721a41e483e11f6e6f237a18014` (clean working tree).
Authority: NCBI C/C++ commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`;
comparison executable SHA-256
`e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.
The initial [Stage C record](README.md), [Stage A](../tlosan_stage_a/README.md),
[Stage B](../tlosan_stage_b/README.md), and
[PD-TLOSAN-LOCAL-GENCODE-32](../../product_decisions/PD-TLOSAN-LOCAL-GENCODE-32.md) set the
source, option, and product boundaries. This continuation supersedes the
initial Stage C record's statements that ambiguity, lowercase masking, full
candidate order, and WordFinder initial HSP parity are still untested.

## NCBI → Rust boundary reached

The pinned `blast_objmgr_tools.cpp:427-474` constructs local-subject ncbi2na
from ncbi4na, resolves ambiguous positions using `CRandom(subject length)`,
and retains ncbi4na for later reevaluation. `random_gen.cpp:98,227-230,287-308`
defines the generator. `blast_engine.c:747-841` translates the preliminary
ncbi2na subject in frame order `+1,+2,+3,-1,-2,-3` and converts subject
lowercase mask ranges. `aa_ungapped.c:478-614` scans ordered word pairs with
one diagonal state surviving frame boundaries and saves sorted initial HSPs.
Rust now ports and compares these steps for the internal BLOSUM62, word-size-3,
one-query, one-subject-at-a-time path. NCBI snippets and line references appear
above the corresponding Rust code in `search_seed.rs` and `search_init.rs`.
The public CLI still exits with the explicit Stage C–E unimplemented error.

| Fixture / option state | NCBI trace | Rust comparison | Result |
| --- | --- | --- | --- |
| Fixed 13-subject local `-subject`, code 1, `-comp_based_stats 0 -seg no -sum_stats false` | 78 full preliminary frame byte strings; 2,030 ordered candidate pairs; 11 saved initial HSPs | Full ordered vectors; initial raw scores, seed positions, frame-relative `q_start/s_start/length`, and HSP order | Exact on these fields |
| Lowercase `-lcase_masking`, code 1 | Ordered 277 candidates; two initial HSPs | Same ordered candidates and HSP fields; uppercase, lowercase, mixed subject records | Exact on these fields |
| Eleven IUPAC ambiguity masks, code 1 | Six preliminary frame byte strings; 181 ordered candidates; one initial HSP | Same bytes, candidates, and initial HSP fields | Exact on these fields |
| Code 32, Stage A comparison-only C++ API `-db` oracle | 199 ordered candidates; two initial HSPs: `(+1, q/s start 0/0, length 120, raw 656)` and `(-3, 21/102, length 6, raw 16)` | Same 199 pairs and both initial HSPs, including order and internal coordinates | Exact for this preliminary path; DB statistics are **not** local-subject oracles |

The code-32 C++ API probe's format-6 bytes match the unprobed Stage A output
(SHA-256 `126fd98f7bea5a8bcf101d2328c9243ad90f6c5f35c3c7c5a8de18adf7c80026`).
The code-1 subject-32 control also matches Stage A. Probe `.so` files were
loaded only into comparison runs of NCBI; LOSAT build and runtime do not load
or call NCBI code. Code 32 remains an approved local-subject product exception,
not a license to change scoring, ordering, linking, or formatting.

## First unresolved stage and residual classification

| Boundary | Fixed NCBI observation | Rust state / required resolution |
| --- | --- | --- |
| WordFinder saved initial HSP | The fixed local fixture has ten raw-656 initial HSPs and one raw-647 stop HSP. The ambiguous initial HSP is raw 656. | Exact initial HSP fields and order on the fixed fixture. This is the **last proven boundary**. |
| Gapped HSP and frame aggregation | `blast_engine.c:478-599,835-899` invokes gapped construction, purges, adjusts chunk offsets, appends per-frame lists, then links/evaluates. | No TBLASTN Rust gapped construction, frame merge, chunk merge, or stage trace. First unresolved boundary is the NCBI gapped-HSP call after WordFinder. Candidate and initial-HSP parity does not establish gapped-HSP parity. |
| Traceback reevaluation and deletion | The fixed ambiguous prelim raw 656 becomes final raw 646 when the retained ncbi4na is reevaluated; the stop HSP remains raw 647 under raw-isolation options. `blast_traceback.c:640-668` is the relevant reevaluation path. | No Rust TBLASTN reevaluation/deletion trace. Cannot compare deletion order or assert that the 656→646 change is implemented. |
| Profile coverage | NCBI supports the Stage B matrix-specific defaults, word sizes, multiple queries, and long-subject chunks. | Current internal Rust path covers BLOSUM62/word 3 and one query/subject at a time. Other matrix/word-size paths, query context offsets, and long-subject chunk boundaries remain unported and unverified. |
| Stage D statistics/linking | NCBI's default and raw-isolation final values are retained below. | **Not started or accepted.** Stage C gapped/post-reevaluation HSP parity is a prerequisite. No Rust raw/bit/E-value/final rank/deletion-order comparison is possible yet. |
| Stage E display | NCBI outfmt 0/6/7 remain the output authority. | Display-only residuals cannot be isolated while HSP/statistical residuals remain. The public CLI must continue to fail explicitly. |

This is a bounded Stage C advance, **not Stage C completion or TBLASTN
certification**. In particular, the internal test's equal candidate count and
initial HSP count are backed by ordered field-by-field comparisons; they do not
substitute for the missing downstream comparisons.

## NCBI final reference values; Rust final comparison pending

[score_reference_20260923](score_reference_20260923/manifest.txt) was generated
from the same fixed *local-subject* fixture and pinned executable. It records
all ranks, raw/bit scores, and E-values for the raw-isolation and default
profiles. It verifies that the `(subject, raw score)` order matches the
original Stage C outfmt-6 files. Selected rows:

| Profile | NCBI ranks 1–9 | NCBI rank 10 | NCBI rank 11 | Rust final |
| --- | --- | --- | --- | --- |
| Raw isolation (`comp_based_stats=0`, SEG off, sum statistics off) | nine raw 656 / bit 257 / `1.66e-94`, ordered `tie_b`, `tie_a`, `partial_codon`, `minus3`, `minus2`, `minus1`, `plus3`, `plus2`, `plus1` | `internal_stop`: raw 647 / bit 253 / `3.92e-93` | `ambiguous`: raw 646 / bit 253 / `5.57e-93` | unavailable |
| Default (`comp_based_stats=2`, SEG on, sum statistics on) | same nine IDs in order, raw 640 / bit 251 / `4.42e-92` | `ambiguous`: raw 631 / bit 247 / `1.19e-90` | `internal_stop`: raw 630 / bit 247 / `1.79e-90` | unavailable |

`stop_only` yields no final row. These final values are NCBI reference data,
not passing Rust comparisons. The Stage A `-db` header/statistical values and
code-32 API output must not be substituted for the local `-subject` contract.

## Rerun and evidence integrity

From repository root, with fresh output directories:

```bash
python3 docs/evidence/tlosan_stage_c/run_six_frame_oracle.py /tmp/tlosan-c-six-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_trace.py /tmp/tlosan-c-six-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_frame_trace.py docs/evidence/tlosan_stage_c/run_20260923 /tmp/tlosan-c-frames-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_candidate_trace.py docs/evidence/tlosan_stage_c/run_20260923 /tmp/tlosan-c-candidates-recheck
python3 docs/evidence/tlosan_stage_c/run_lowercase_oracle.py /tmp/tlosan-c-lowercase-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_trace.py --lcase-masking /tmp/tlosan-c-lowercase-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_candidate_trace.py --lcase-masking /tmp/tlosan-c-lowercase-recheck /tmp/tlosan-c-lowercase-pairs-recheck
python3 docs/evidence/tlosan_stage_c/run_ambiguity_oracle.py /tmp/tlosan-c-ambiguity-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_trace.py /tmp/tlosan-c-ambiguity-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_frame_trace.py /tmp/tlosan-c-ambiguity-recheck /tmp/tlosan-c-ambiguity-frames-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_candidate_trace.py /tmp/tlosan-c-ambiguity-recheck /tmp/tlosan-c-ambiguity-pairs-recheck
python3 docs/evidence/tlosan_stage_c/run_ncbi_score_reference.py docs/evidence/tlosan_stage_c/run_20260923 /tmp/tlosan-c-score-recheck
gcc -shared -fPIC -std=c11 -O2 -o /tmp/tlosan-candidate-probe.so docs/evidence/tlosan_stage_c/ncbi_candidate_trace.c -ldl
LD_PRELOAD=/tmp/tlosan-candidate-probe.so bash docs/evidence/tlosan_stage_a/run_api_oracle.sh /tmp/tlosan-c-code32-candidate-recheck
gcc -shared -fPIC -std=c11 -O2 -o /tmp/tlosan-word-probe.so docs/evidence/tlosan_stage_c/ncbi_wordfinder_trace.c -ldl
LD_PRELOAD=/tmp/tlosan-word-probe.so bash docs/evidence/tlosan_stage_a/run_api_oracle.sh /tmp/tlosan-c-code32-word-recheck
python3 docs/evidence/tlosan_stage_c/parse_code32_api_trace.py /tmp/tlosan-c-code32-candidate-recheck /tmp/tlosan-c-code32-word-recheck /tmp/tlosan-c-code32-trace-recheck
(cd LOSAT && cargo test --lib algorithm::tblastn::)
```

The probe scripts compare final NCBI bytes against their unprobed fixture.
The retained `*.sha256` files verify artifact bytes with `sha256sum -c` from
each artifact directory. New `code32_20260923/outputs.sha256` and
`score_reference_20260923/outputs.sha256` cover their extracted tables and
manifests. A fresh output path is required by each script.
