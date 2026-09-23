# TLOSAN v0.2.0 Stage C work record — incomplete

**Latest checkpoint:** [CONTINUATION_20260924.md](CONTINUATION_20260924.md)
records exact gapped HSP comparisons on bounded fixtures and the still-open
Stage C gate.

**Prior checkpoint after `ad3e8381`:** [CONTINUATION_20260923.md](CONTINUATION_20260923.md) records the exact preliminary comparisons. The status tables below describe the initial checkpoint.

Branch: feature/tlosan-tblastn-v0.2.0. Starting tree: clean at b478d701.
NCBI source: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4.
NCBI executable: tblastn 2.17.0+, SHA-256
e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0.
Target: WSL2 Linux native, local subject, one thread. The NCBI source worktree has
line-ending modifications; the Stage A API rerun archives and verifies the frozen
commit. None of the following is TBLASTN search parity or release certification.

## Stage B gate

Before search work, the Stage B checker ran in a fresh /tmp directory and passed:
27 gc.prt tables, 1,728 codons, nine matrix gap tables, the code-32 API
control, and CLI accept/reject cases. Focused Rust tests for TBLASTN arguments,
genetic codes, and cli_v2 passed. No new Stage B discrepancy was found. The
Stage A C++ API oracle was rerun in
/tmp/tlosan-stage-c-code32-api-20260923; all eleven output checksums passed.
Its code-32 format-6 file matches the retained Stage A file byte for byte:
126fd98f7bea5a8bcf101d2328c9243ad90f6c5f35c3c7c5a8de18adf7c80026.
This is comparison-only code outside LOSAT's runtime and build.

## Pinned NCBI → LOSAT call path and state

All source paths below are relative to c++/src/algo/blast. The NCBI references
are from the frozen commit, verified via the Stage A archive run.

| Order | NCBI owner | Input, coordinates and effect | LOSAT owner / status |
| --- | --- | --- | --- |
| 1 | api/blast_setup_cxx.cpp:486-651; core/blast_util.c:112-127; core/blast_filter.c:337-370 | Protein query becomes NCBISTDAA with NULLB=0 sentinels. SEG masks the working query, retaining unmasked residues for later identity; protein context frame is 0 and each query context advances length+1. | blastp/encoding.rs is shared for the same protein query encoding and SEG state. Current internal scanner handles one query. |
| 2 | core/lookup_wrap.c:83-110; core/blast_engine.c:1023-1050; core/blast_options.c:1163-1169 | For word size 3, the AA lookup indexes protein query once before subjects. The protein word finder is selected for both BLASTP and TBLASTN. Word sizes above 5 select a compressed lookup. | tblastx/lookup and tblastx/blast_aascan are shared only for BLOSUM62 word size 3; search_seed.rs uses the same threshold 13 control. Other matrices/word sizes are not exposed by this internal entry. |
| 3 | api/blast_setup_cxx.cpp:736-853; api/blast_objmgr_tools.cpp:428-473; core/blast_engine.c:1458-1466 | Local subject keeps uncompressed ambiguous 4na for reevaluation, plus compressed ncbi2na without sentinels for prelim. Ambiguities in local ncbi2na are resolved with CRandom(subject base length). The local source normally supplies code 1, but PD-TLOSAN-LOCAL-GENCODE-32 requires LOSAT to apply the selected db_gencode throughout. | search_seed.rs validates the selected 27-ID table, including 32. It currently rejects ambiguous or lowercase subjects explicitly; exact ncbi2na ambiguity resolution, lowercase mask handling and preserved 4na reevaluation remain open. |
| 4 | core/blast_util.c:1045-1101; core/blast_engine.c:747-785,804-813 | Translate one subject in context order +1,+2,+3,-1,-2,-3. Each frame begins with NULLB; frame_offsets advance translated length+1 so boundaries share a sentinel. Subject AA offsets are frame-relative, while the original nucleotide length is retained. | tblastx/translation.rs is reused for unambiguous frame amino acids and sentinels; search_seed.rs visits six frames in that order. It has no NCBI packed-buffer implementation and does not claim equivalence for ambiguous DNA. |
| 5 | core/blast_engine.c:815-841; core/aa_ungapped.c:492-505,763-769 | For each frame, translate subject mask ranges at positive context 0 and negative context 3; run BlastAaWordFinder. Query offsets are in concatenated query contexts; subject offsets are frame-relative AA. One diagonal structure survives between frames. | search_seed.rs scans words for one unmasked/unambiguous subject, preserving returned pair order. Lowercase subject mask conversion, diagonal state, ungapped extension and candidate parity tracing are open. |
| 6 | core/blast_engine.c:478-599; core/blast_gapalign.c:3697-4008; core/blast_hits.c:2455-2603 | Reset init hits per chunk, get gapped HSPs, purge common endpoints, sort by score, adjust chunk offsets, merge overlapping chunks. The HSP query offset becomes context-relative when saved; subject offset stays frame-relative. | No TBLASTN Rust HSP construction or chunk merge yet. BLASTP's existing prelim function is coupled to BLASTP parameters and cannot simply be called on six translated frames. |
| 7 | core/blast_engine.c:835-899; core/blast_hits.c:2809 onward | Append each frame HSP list with kHspNumMax. After all frames restore original subject, then link or get E-values and reap by preliminary E-value. Do not prune members before linking. | No TBLASTN HSP append, linking, or reap yet. Linking/statistical formulas belong to Stage D after Stage C HSP parity. |
| 8 | core/blast_engine.c:1481-1538; core/blast_traceback.c:199-259,613-712,808-853; core/blast_kappa.c:2942-2981 | Ungapped mode reevaluates with the ambiguous 4na subject after prelim and can relink/reap. Gapped mode proceeds to traceback, composition redo where enabled, common-endpoint/containment purge, score and E-value processing. | TBLASTN reevaluation and deletion are open. Stage D owns linking, composition adjustment, statistics and their post-traceback effect once preliminary HSPs match. |

The existing BLASTP extension and gapped alignment code shares NCBI source
functions with TBLASTN but currently computes BLASTP-specific cutoffs,
subject length, hitlist, and postprocessing within one large engine. Calling it
per translated frame would change parameter-update timing and frame HSP
merging, so no such shortcut was added. TBLASTX also shares amino-acid
translation/scan primitives but has six query contexts and different
statistics, so its engine is not reused as a whole.

## Reproducible fixture and comparison

From repository root, with fresh output paths:

    python3 docs/evidence/tlosan_stage_c/run_six_frame_oracle.py /tmp/tlosan-c-six-new
    python3 docs/evidence/tlosan_stage_c/run_ncbi_trace.py /tmp/tlosan-c-six-new
    (cd /tmp/tlosan-c-six-new && sha256sum -c inputs_outputs.sha256 && sha256sum -c trace_outputs.sha256)
    bash docs/evidence/tlosan_stage_a/run_api_oracle.sh /tmp/tlosan-c-api-new
    (cd /tmp/tlosan-c-api-new && sha256sum -c outputs.sha256)
    (cd LOSAT && cargo test --lib algorithm::tblastn::search_seed)

The committed run_20260923 directory contains all inputs, full NCBI commands,
binary hash, raw score/frame tabular output, pairwise output, and hashes.
run_ncbi_trace.py compiles a temporary LD_PRELOAD probe against the NCBI
WordFinder ABI. It records the sorted ungapped init HSPs immediately after
BlastAaWordFinder; the probe is comparison-only and its final NCBI output
matches the unprobed output byte for byte. Its manifest records the probe
source hash. It observed 78 calls (six per 13 subjects) and 11 saved init HSPs.
The initial HSP seed pair (query 3, subject 3) of each of ten unambiguous
records exists in Rust's corresponding word scan. The ambiguous record is
explicitly rejected by the incomplete Rust seed path. This is a necessary
seed check, not a complete candidate-set comparison.
The main 120-aa query is placed in six subject reading frames; additional
records cover trailing partial codons, one ambiguous base, one internal stop,
two equal-score subjects, a stop-only no-hit subject, and a low-complexity
subject. A separate K-rich query/AAA subject control gives one NCBI hit with
SEG off and none with SEG on.

| Stage / field | New NCBI observation | LOSAT observation | Difference / disposition |
| --- | --- | --- | --- |
| Frame translation | Six full-length 120-aa matches: +1/+2/+3/-1/-2/-3 | All six translated AA strings contain the NCBI reported subject AA; unit test passes. | No mismatch for these unambiguous frame strings. This does not prove all NCBI2na input encodings. |
| Protein lookup and seed | NCBI WordFinder made 78 frame calls and saved 11 init HSPs; ten unambiguous saved seed pairs are (query 3, subject 3). | Each of those ten seed pairs exists in the Rust frame scan. Code 32 changes the seed list versus code 1 (199 versus 160 offsets for the Stage A subject). | Full candidate set and candidate order are unverified. Saved init HSPs do not expose every scanned word. |
| SEG | K-rich query gives 0 rows with SEG on, 1 row with SEG off. | Internal seed test checks masked count 0 and unmasked count 41,536. | Observable control agrees; full candidate trace remains missing. |
| Raw score | WordFinder init HSPs: ten score 656, one internal-stop score 647; ambiguous init score is 656, then final reported score is 646. | No TBLASTN Rust HSP/raw score. | First unverified Rust stage after seed scan: ungapped and gapped HSP construction. NCBI's ambiguous score change corresponds to blast_traceback.c:640-668 reevaluation. |
| Internal HSP coordinates | NCBI printed 1..120 query and the six frame-specific nucleotide spans. | Only seed query/subject AA offsets exist. | No internal HSP coordinate comparison yet. |
| HSP order | The eleven NCBI rows start tie_b, tie_a, partial_codon, minus3, minus2, minus1, plus3, plus2, plus1, internal_stop, ambiguous. | No TBLASTN HSP list. | Tie and frame order parity unverified. |
| No-hit | stop-only subject emits no row. | No final search result exists. | Unverified. |
| Code 32 | Fresh C++ API oracle matches Stage A code-32 format-6 bytes. | Selected code changes internal seed list; no Rust HSP. | Search/HSP parity unverified. |

## Final local verification

- cargo fmt --all -- --check: pass.
- cargo test --lib algorithm::tblastn::search_seed -- --nocapture:
  four tests passed; captured in rust_seed_test.log. This checks the six
  translated frame strings against new NCBI output, ten NCBI saved seed pairs,
  SEG on/off, and selected code 32.
- cargo clippy --all-targets -- -D warnings: pass.
- cargo build --release: pass.
- Final Stage B checker: pass on the new release binary (27 source tables,
  1,728 codons, nine matrix tables, code-32 API control, CLI cases).
- Public release CLI with -db_gencode 32 exits 1 and prints:
  "TBLASTN local search is unimplemented (Stages C-E)".
- Fixture and NCBI trace SHA-256 manifests: all entries pass. The temporary
  NCBI API output checksums pass, and code32_fmt6 matches Stage A raw bytes.
- git diff --check on Rust, scripts and documentation: pass. The unedited
  NCBI pairwise outputs and warning text retain upstream trailing spaces, so
  a whole-tree whitespace check reports those raw oracle bytes.
- No benchmark or Wasm gate was run: the TBLASTN search engine does not yet
  produce results and the partial seed function is internal only.

The public TBLASTN CLI still fails explicitly after validation. No statistical,
formatting, linking, output-byte, or performance parity is claimed. Stage C
cannot be accepted until all seed, ungapped/gapped HSP, reevaluation, deletion,
internal coordinate, and ordering rows have NCBI stage trace comparisons and
zero unexplained differences. Stage D can then address translated subject
statistics, unequal-gap/intron linking, composition redo, bit score, E-value,
and pruning timing; Stage E owns format 0/6/7 bytes.
