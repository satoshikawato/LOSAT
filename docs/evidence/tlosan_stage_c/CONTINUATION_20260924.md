# TLOSAN v0.2.0 Stage C continuation — gapped boundary advanced, Stage C incomplete

Branch: `feature/tlosan-tblastn-v0.2.0`; starting LOSAT commit
`afb719c413464be2903c8bb46058260a42a7ebae` (clean tree).
NCBI authority: C/C++ commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`;
comparison CLI SHA-256
`e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.
This follows [the 2026-09-23 checkpoint](CONTINUATION_20260923.md);
[Stage A](../tlosan_stage_a/README.md), [Stage B](../tlosan_stage_b/README.md),
and [PD-TLOSAN-LOCAL-GENCODE-32](../../product_decisions/PD-TLOSAN-LOCAL-GENCODE-32.md)
remain the authority and product boundaries. This record is **not** Stage C,
D, or E acceptance and is not TBLASTN certification.

## First unresolved boundary and result

The preceding checkpoint had exact ordered candidate and WordFinder initial-HSP
comparisons but no `GetGappedScore` trace. The new comparison-only
[native probe](ncbi_gapped_trace.c) records input parameters, every initial
HSP, every returned gapped HSP, and traceback list transitions. It is loaded
only into pinned NCBI oracle runs; LOSAT never links or calls it.
The [runner](run_ncbi_gapped_trace.py) checks that probing leaves the pinned
NCBI final output bytes unchanged.

| Fixture and profile | NCBI `GetGappedScore` | Rust internal comparison | Current boundary |
| --- | --- | --- | --- |
| Fixed 13-subject local `-subject`, BLOSUM62/word 3, code 1, SEG off, composition 0, sum statistics off | 11 initial → 11 gapped HSPs; ten raw 656, one raw 647, in the recorded subject/frame order | All 11 raw scores, frame-local query/subject start/end, gapped starts, frame, and order exact | [gapped output](gapped_20260924/gapped_output.tsv) |
| Lowercase masking | Two initial → two gapped; raw 656 each | Both full gapped rows exact | [gapped output](gapped_lowercase_20260924/gapped_output.tsv) |
| Eleven IUPAC ambiguity masks | One initial → one gapped; raw 625 | Full gapped row exact | [gapped output](gapped_ambiguity_20260924/gapped_output.tsv) |
| Code 32 via Stage A comparison-only NCBI C++ API `-db` oracle | Two gapped HSPs in order: `+1` raw 656, query/subject 0–120; `−3` raw 16, query 21–27, subject 102–108 | Both raw scores, internal coordinates, gapped starts, frames, and order exact | [API trace](gapped_code32_20260924/gapped_code32_api.tsv) |

For the fixed local fixture, `Blast_TracebackFromHSPList` is called **twice**
per retained subject. The first call has `fence_hit=1`, so
`blast_traceback.c:1644-1684` refetches the whole nucleotide subject and
retries with the same HSP list. The second call has `fence_hit=0`.
The ambiguous HSP is raw **656** at WordFinder, `GetGappedScore`, and the
first traceback return; it becomes raw **646** on the second traceback return,
with frame `+1` and internal query/subject coordinates 0–120 unchanged.
The probe recorded no `Blast_HSPReevaluateWithAmbiguitiesGapped` calls in
this fixture: the score change occurs during the full-subject gapped
traceback, so it must not be attributed to the later reevaluation loop.
The [all-HSP traceback table](gapped_20260924/traceback_hsps.tsv) records all
22 calls, score/coordinate transitions, fence flags, and list sizes. Every
fixed-fixture list has one HSP before and after both calls; this fixture has
**no observed traceback deletion**. Do not infer deletion-order parity from
that absence. The NCBI source path is
`blast_hits.c:1147-1222` (target translation and fence),
`blast_traceback.c:259-721,1644-1684` (traceback, rollback/retry, reevaluation,
sorting and purge), and `blast_gapalign.c:3664-4091,4549-4650`.

The Rust internal diagnostic reuses the existing NCBI-derived BLASTP protein
gapped alignment primitive. On the fixed small HSPs, it compares all 11 first-pass fence flags and
all 11 post-fence raw scores/internal coordinates against NCBI, including
the 656→646 ambiguity change. The reused primitive now
exposes its internal fence flag without changing BLASTP's callers. Its current
TBLASTN target-range model is bounded to these full-frame small fixtures. It
does not implement general partial target windows, long-sequence chunks,
frame aggregation, or public output. The public TBLASTN CLI retains the
explicit Stage C–E unimplemented error.

The code-32 API run is a **database** oracle. The probe's code-32 format-6
bytes equal the unprobed Stage A API output; code 1 is a matching control.
Its statistics and headers are not local `-subject` oracles. The approved
subject-code exception remains limited to selected non-default translation.

## Difference classification and remaining gates

- **Resolved on bounded fixtures:** The first missing `GetGappedScore`
  boundary; full ordered gapped HSP fields for the fixed, lowercase, IUPAC,
  and code-32 API fixtures. The fixed raw 656→646 transition and first-pass
  fence retry are reproduced as an internal diagnostic.
- **Stage C remains open:** General NCBI target-range translation/fence
  lifecycle; multi-HSP endpoint purge and deletion order; six-frame list
  append and cross-frame ordering; long-subject chunk offsets/merge; multiple
  protein queries and context offsets; matrix/word-size-specific paths;
  traceback reevaluation and deletion on fixtures that actually exercise
  them. No all-candidate/all-internal-HSP comparison exists for these profiles.
- **Stage D has not begun:** No TBLASTN effective-length, linking, composition
  mode-2/Kappa redo, E-value/bit-score, filter, hitlist, rank, tie, or
  deletion-order Rust comparison. `-max_intron_length 0` must retain NCBI
  default linking; BLASTP `do_link_hsps=false` is not a TBLASTN rule.
- **Stage E is gated:** No TBLASTN local `-subject` output in format 0/6/7;
  therefore no display-only residual can be classified. No Stage E pass,
  TBLASTN completion, release certification, or performance claim is made.

## Reproduce and verification

From repository root, each trace command requires a new output path:

```bash
python3 docs/evidence/tlosan_stage_c/run_ncbi_gapped_trace.py docs/evidence/tlosan_stage_c/run_20260923 /tmp/tlosan-c-gapped-fixed-new
python3 docs/evidence/tlosan_stage_c/run_ncbi_gapped_trace.py docs/evidence/tlosan_stage_c/lowercase_20260923 /tmp/tlosan-c-gapped-lower-new --lcase-masking
python3 docs/evidence/tlosan_stage_c/run_ncbi_gapped_trace.py docs/evidence/tlosan_stage_c/ambiguity_20260923 /tmp/tlosan-c-gapped-ambiguity-new
gcc -shared -fPIC -std=c11 -O2 -o /tmp/tlosan-gapped-probe.so docs/evidence/tlosan_stage_c/ncbi_gapped_trace.c -ldl
LD_PRELOAD=/tmp/tlosan-gapped-probe.so bash docs/evidence/tlosan_stage_a/run_api_oracle.sh /tmp/tlosan-c-code32-gapped-api-new
python3 docs/evidence/tlosan_stage_c/parse_code32_gapped_trace.py /tmp/tlosan-c-code32-gapped-api-new /tmp/tlosan-c-code32-gapped-extract-new
(cd LOSAT && cargo test --lib algorithm::tblastn::search_gapped::tests)
```

All four retained evidence directories include `outputs.sha256`; every
checksum passed with `sha256sum -c`. The fresh Stage A API run's 11
format-0/6/7 output checksums also passed. The five new Rust oracle tests
passed. The comparison-only NCBI trace changed no final bytes on the three
local fixtures; the API code-1 and code-32 format-6 outputs equal Stage A's
unprobed bytes.

Program and task: TBLASTN, local `-subject` code 1 plus comparison-only
code-32 API `-db`. Target: native Linux, one thread. LOSAT owner:
`src/algorithm/tblastn/search_gapped.rs`; shared fence reporting:
`src/algorithm/blastp/gapalign.rs`. Focused first difference: gapped HSP
construction absent before this checkpoint; remaining first general boundary:
partial target translation/fence and list pruning on multi-HSP fixtures.
Accepted exception: only the approved non-default local subject genetic code.
Native/Wasm comparison, broad output-byte comparison, benchmarks, and release
audit: not applicable until the Stage C/D/E gates are met.
