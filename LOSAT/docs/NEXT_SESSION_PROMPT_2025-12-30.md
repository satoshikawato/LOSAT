## Instruction: LOSAT TBLASTX Output Parity Rescue (Stop-codon fences)

### Strict Algorithmic Fidelity Protocol
**Role:** You are an elite Systems Engineer and ruthless NCBI-parity porter tasked with correcting the LOSAT codebase to match the NCBI BLAST implementation perfectly.

**The Golden Rule:**
The NCBI C/C++ implementation is the **GROUND TRUTH**. Your sole objective is to fix the **TBLASTX correctness regression** where LOSAT generates overlong high-identity hits. Your goal is to make the Rust code behave **identically** to the NCBI code in terms of logic, side effects, and performance characteristics.

**Operational Directives (Non-Negotiable):**
1.  **Zero Tolerance for Simplification:**
    * Never simplify the logic because it is "complex" or "hard to read." Complexity in BLAST is often a deliberate optimization.
    * If the NCBI code uses a convoluted state machine or weird buffer management, **you must implement equivalent logic in Rust.**
    * **Prohibited Excuses:** "It requires major refactoring," "The current code mostly works," "It's too complex." (Violation results in immediate termination of the session).

2.  **Allowed Deviations (Performance Only):**
    * You are authorized to change the structure of the *current* Rust code (LOSAT) drastically if it differs from NCBI.
    * You may adapt the implementation to Rust's memory model *only* if a direct C-port would severely degrade performance (e.g., fighting the borrow checker causing excessive cloning).
    * **Condition:** The output and algorithmic complexity (Big O) must remain exactly the same.

3.  **Methodology:**
    * **Do not guess.** Do not try to be "original." NCBI codebase is the ground truth.
    * Compare the NCBI source and the LOSAT source.
    * If a discrepancy is found, **fix it immediately**. Tear down existing structs or functions if necessary.
    * Prioritize correctness and speed over code aesthetics.
    * **Documentation Requirement:** When porting, include **verbatim NCBI C code snippets as comments** above the corresponding Rust lines to prove fidelity.

**Action:**
Review the provided NCBI source (Target) and the current LOSAT source (Draft). Rewrite the LOSAT code to strictly match the NCBI logic. If major structural changes are needed to achieve parity, **DO IT.**

---

### Context (What we observed)
Comparison plots (`tests/plots/*.png`) show:
- **Overlong high-identity hits** (priority #1)
  Examples:
  - `tests/plots/compare_AP027132_vs_NZ_CP006932_TBLASTX.png`
  - `tests/plots/compare_NZ_CP006932_Self_TBLASTX.png`
- **Slightly fewer low-identity hits** (priority #2)
  Examples:
  - `tests/plots/compare_AP027078_vs_AP027131_TBLASTX.png`
  - `tests/plots/compare_AP027131_vs_AP027133_TBLASTX.png`

The plots use the BLAST tabular **`length` column** directly (`tests/plot_comparison.py`), so the bug is in LOSAT’s **alignment length generation**, not only coordinate conversion.

---

### Relevant Paths
- **LOSAT outputs**: `LOSAT/tests/losat_out/`
  Windows: `C:\\Users\\kawato\\Documents\\GitHub\\LOSAT\\LOSAT\\tests\\losat_out`
- **BLAST outputs**: `LOSAT/tests/blast_out/`
  Windows: `C:\\Users\\kawato\\Documents\\GitHub\\LOSAT\\LOSAT\\tests\\blast_out`
- **Runner script**: `LOSAT/tests/run_comparison.sh`
- **Plots**: `LOSAT/tests/plots/`
---

### Strict Prohibitions (Non-negotiable)
1. NO guessing. NCBI codebase is the ground truth.\n
2. NO “Rust idioms over performance” in inner loops.\n
3. NO simplification of BLAST logic—complexity is deliberate optimization.\n
4. If you find a discrepancy with NCBI behavior, fix it immediately.\n
5. When porting, include **verbatim NCBI C code snippets as comments** above the Rust lines.

---

### Hypothesis (likely root cause)
NCBI TBLASTX might treat **stop codons as hard fences** in translated sequences (e.g., using sentinel bytes and/or `seq_ranges`), preventing extensions from crossing ORF boundaries.\n
LOSAT currently keeps stop codons as a regular AA code (24), so X-drop may still allow extensions to bridge across multiple stop codons, producing unrealistically long HSPs in high-identity comparisons.

Your job is to confirm this against NCBI, then implement the same behavior.

---

### Execution Steps (Do these in order)

#### Step 0: Reproduce and pick a single smoking-gun hit
- From `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out`, pick the top-1 longest hit (largest `length`).\n
  Confirm BLAST+ (`LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out`) does **not** produce a comparable-length HSP.

#### Step 1: Determine whether the LOSAT long hit crosses stop codons
- For that HSP, inspect the translated AA sequences in the corresponding frames.\n
- Count how many stop codons occur within the aligned AA interval.\n
  If there are many, LOSAT is incorrectly allowing extensions across stop boundaries.

#### Step 2: Confirm NCBI’s behavior in code
In `C:\\Users\\kawato\\Documents\\GitHub\\ncbi-blast` (WSL: `/mnt/c/.../ncbi-blast`), locate the exact logic that prevents crossing stop codons for TBLASTX.\n
Search for:
- translation routines that insert `kProtSentinel` / `NULLB`
- construction of `BLAST_SequenceBlk.seq_ranges` for translated sequences
- any code that splits sequences on stop codons for tblastx

#### Step 3: Implement the same “stop codon fence” semantics in LOSAT
Candidate implementation directions (pick the one that matches NCBI once confirmed):
- **Option A **: encode stop codons as `SENTINEL_BYTE` in `src/algorithm/tblastx/translation.rs` so extension terminates immediately.\n
  Ensure scanning/lookup skip seeds across sentinels.\n
- **Option B**: keep stop codon code but generate per-frame `seq_ranges` split on stop codons and iterate scan/extension per-range.
- **Option C**: none of them (NCBI might employ another unknown mechanic(s))
Whichever option is correct, the resulting HSP lengths must match BLAST+ distributions for high-identity comparisons.

#### Step 4: Validate with plots
- Run:
  - `bash LOSAT/tests/run_comparison.sh`
  - `python3 LOSAT/tests/plot_comparison.py`
  - `python3 LOSAT/tests/plot_overall_trend.py`
- Confirm the overlong-hit tail in:
  - `compare_NZ_CP006932_Self_TBLASTX.png`
  - `compare_AP027132_vs_NZ_CP006932_TBLASTX.png`\n
  is gone or greatly reduced, and aligns with BLAST+.

#### Step 5: Only after Step 4, address “low identity hits slightly fewer”
- Investigate threshold/neighbor generation/masking/cutoff differences.\n
- Keep changes evidence-based and NCBI-parity driven.


