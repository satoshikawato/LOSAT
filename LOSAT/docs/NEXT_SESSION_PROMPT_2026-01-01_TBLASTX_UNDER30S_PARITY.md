## NEXT SESSION PROMPT — 2025-12-31

## Instruction: Strict Algorithmic Fidelity Protocol (Output is Ground Truth)

**Role:** You are an elite Systems Engineer tasked with correcting LOSAT to match NCBI BLAST behavior, with **NCBI output as ground truth**.

**Hard constraints:**

- **Output parity is non-negotiable.** Match NCBI BLAST+ results (hit set) for tblastx.
- “Simplification” that changes output is prohibited.

**Relaxed constraints (explicit user agreement):**

- **Side effects and computational cost may change if they become smaller**, as long as output parity is preserved.

**Operational rule (session discipline):**

- You must keep an always-up-to-date devlog: append to `LOSAT/docs/TBLASTX_NCBI_PARITY_DEVLOG.md` roughly **every ~10 minutes**, with a timestamp (use `date -Is`).

## Current state (what has been implemented already)

### Key code changes already landed

- NCBI lookup parity:
  - `LOSAT/src/algorithm/tblastx/lookup.rs` `build_ncbi_lookup()` was rewritten to match:
    - `BlastLookupIndexQueryExactMatches()` (group exact words first)
    - `s_AddWordHits()` semantics (**index exact hits for low self-score words** even if `self_score < threshold`)
  - New lookup invariants were unit-tested:
    - low self-score exact word indexing (AAA self-score 12 < 13 still indexed)
    - `longest_chain` correctness

- “HF suppression”:
  - `--max-hits-per-kmer` is now **optional** (`Option<usize>`). Default is **disabled** (NCBI parity).

- Stop handling:
  - Defaults are NCBI-parity:
    - stop seeds included by default (`--include-stop-seeds` default true)
    - stop-stop score default is `*-* = +1` (NCBI BLOSUM62), with opt-out available

- ScanSubject buffering parity:
  - `offset_pairs` sized to NCBI: `4096 + lookup.longest_chain`
  - scan loop “hits==0 break” was replaced with a guard that only breaks if `scan_range` fails to advance (infinite-loop safety).

### Files touched / important references

- Implementation:
  - `LOSAT/src/algorithm/tblastx/lookup.rs`
  - `LOSAT/src/algorithm/tblastx/utils.rs`
  - `LOSAT/src/algorithm/tblastx/args.rs`
  - `LOSAT/src/utils/matrix.rs`
  - `LOSAT/src/algorithm/tblastx/extension.rs`
- Tests:
  - `LOSAT/tests/unit/tblastx/lookup.rs`
  - `LOSAT/tests/unit/tblastx/extension.rs`
- Devlog:
  - `LOSAT/docs/TBLASTX_NCBI_PARITY_DEVLOG.md`

## Known performance data (single-thread)

Diagnostics summary (single-thread, qframes=all(6), sframe=1):

- `kmer_matches` ≈ **1.589e8**
- `ungapped_extensions` ≈ **3.767e6**
- `cutoff_score` failures ≈ **3.723e6**
- `hsps_before_chain` ≈ **4.40e4**
- final hits ≈ **9.55e3**
- wall time ≈ **33.27 sec**

This shows the runtime is dominated by **scan → two-hit gating → ungapped extension volume**, not linking.

## Goal for next session

- **Primary:** NZ_CP006932 self, single-thread, **< 30 sec** (LOSAT), while maintaining **output parity**.
- **Secondary:** keep overlong-tail safety checks passing.

## Highest-suspicion parity/perf gap to investigate next

### Subject masking is likely missing (and can massively affect kmer_matches)

LOSAT currently applies:

- nucleotide DUST mask: query-side only
- amino-acid SEG mask: query-side only

But NCBI’s scan code (`blast_aascan.c` + `masksubj.inl`) is built around `subject->seq_ranges[]` which represent **allowed scanning ranges** (skipping masked regions). If NCBI tblastx masks the subject by default (very likely), LOSAT is currently scanning **too much** and doing extra work (and may also be output-divergent).

## Concrete next tasks (in order)

1) **Confirm NCBI tblastx default subject filtering/masking behavior**
   - Use NCBI source references already in workspace:
     - `ncbi-blast/c++/src/algo/blast/core/blast_filter.c`
     - `ncbi-blast/c++/src/algo/blast/core/blast_setup.c`
     - `ncbi-blast/c++/src/algo/blast/core/blast_aascan.c`
     - `ncbi-blast/c++/src/algo/blast/core/masksubj.inl`
   - Determine whether tblastx applies DUST/SEG to **subject** by default and how it translates into `seq_ranges`.

2) **Implement subject-side masking in LOSAT (parity + performance)**
   - Add subject nucleotide DUST masking (mirror query path) and propagate into translated subject frames.
   - Add subject AA SEG masking (mirror query path) for translated subject frames.
   - Create `seq_ranges` equivalent for each subject frame (unmasked AA intervals) and drive scanning by those ranges (NCBI `masksubj.inl` behavior).
   - Ensure semantics match NCBI: “skip masked regions for seeding/scan” while preserving correct extension/identity computation rules.

3) **Re-measure diagnostics**
   - Run single-thread with `LOSAT_DIAGNOSTICS=1` on NZ self and confirm that:
     - `kmer_matches` and `ungapped_extensions` drop substantially
     - wall time drops under 30 sec

4) **Re-validate output parity**
   - Compare to BLAST+ output using:
     - `LOSAT/tests/analyze_tblastx_hit_gap.py`
     - `LOSAT/tests/check_tblastx_overlong_tail.py`
   - Confirm “worst bin” improvement does not regress and no overlong tail.

## Repro commands (copy/paste)

Build:

```bash
cd LOSAT/LOSAT
cargo build --release
```

Diagnostics (single-thread):

```bash
cd LOSAT/LOSAT/tests
LOSAT_DIAGNOSTICS=1 ../target/release/LOSAT tblastx \
  -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta \
  -o ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.diag.out \
  --query-gencode 4 --db-gencode 4 -n 1
```

Gap analysis + tail safety:

```bash
cd LOSAT/LOSAT/tests
python3 analyze_tblastx_hit_gap.py \
  --blast ./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out \
  --losat ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.diag.out
python3 check_tblastx_overlong_tail.py ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.diag.out
```

## Devlog rule (must follow)

- Every ~10 minutes while working, append a new section to:
  - `LOSAT/docs/TBLASTX_NCBI_PARITY_DEVLOG.md`
- Include:
  - timestamp (`date -Is`)
  - what was tried
  - what was measured (numbers)
  - why you changed direction
  - file/function pointers


