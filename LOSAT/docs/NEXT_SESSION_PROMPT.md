# Instruction: LOSAT Performance Rescue & Strict Porting Protocol

**Role:** You are an elite Systems Engineer and Ruthless Refactorer. Your sole objective is to fix the critical performance regression in LOSAT by enforcing strict parity with NCBI BLAST.

**Context:**
We are facing a CRITICAL PERFORMANCE ISSUE in `tblastx` self-comparison. The Rust implementation is significantly slower than NCBI (timeout vs seconds).
Below is the technical status report and the constraints you must operate under.

---

## 1. TECHNICAL STATUS (From NEXT_SESSION_PROMPT.md)

**Current Status:**
Full test (657kb self-comparison) times out. Small test works (222 hits in 0.043s).

**Root Cause Identified:**
**8 million lookup table entries** for full genome self-comparison:
- 657kb genome → ~219k AA per frame × 6 frames = 1.3M query positions
- Each position generates ~50-100 neighbors (threshold 13)
- Self-comparison: every neighbor matches somewhere
- Result: 8M+ total lookup entries → O(n²) behavior

**The Fundamental Problem:**
NCBI BLAST handles this efficiently. Our Rust implementation has too much overhead per hit or is missing a critical filter mechanism present in the C code.
- **NCBI:** `blast_aascan.c` / `aa_ungapped.c` processes hits in few cycles.
- **LOSAT:** Binary search overhead and likely structural inefficiency in the scan loop.

**Relevant NCBI Files:**
- `blast_aascan.c` (s_BlastAaScanSubject)
- `aa_ungapped.c` (s_BlastAaWordFinder_TwoHit)
- `blast_aalookup.c` (BlastAaLookupFinalize)

---

## 2. OPERATIONAL PROTOCOL (Strict Algorithmic Fidelity)

**The Golden Rule:**
The NCBI C/C++ implementation is the **GROUND TRUTH**.
You must port NCBI's logic **line-by-line** to Rust.

**Strict Prohibitions (Non-Negotiable - Violation = Termination):**
1.  **NO Simplification:** Do not simplify logic because it is "complex" or "messy." Complexity in BLAST is a deliberate optimization.
2.  **NO "Rust Idioms" over Performance:** Do not use `match`, iterators, or safe wrappers if they introduce *any* overhead compared to C pointer arithmetic.
3.  **NO Excuses:** Phrases like "It requires major refactoring," "The current code mostly works," "It's too complex," or "I need to modernize this" are **forbidden**. If the structure is wrong, tear it down and rewrite it.
4.  **NO Guessing:** Do not implement "Option A/B/C" from your imagination. Look at `blast_aascan.c` and implement exactly what they do.
5.  **Immediate Action:** If you find a discrepancy, fix it immediately. Do not waste time on diagnostics if the code clearly differs from NCBI.

**Execution Steps:**

**Step 1: Architectural Analysis (Pre-computation)**
Analyze `blast_aascan.c` and `aa_ungapped.c` specifically regarding the **Lookup Table explosion**.
* How does NCBI avoid the O(n^2) penalty in self-comparison?
* Do they use a specific "Subject Masking" or "Diagonal Skipping" logic inside the inner loop?
* Identify the exact lines of C code that prevent the slowdown we are seeing.

**Step 2: Evidence-Based Implementation**
Rewrite the critical scan loop in `utils.rs` (or wherever appropriate) to match NCBI exactly.
* Use `unsafe` pointers freely.
* **Must include original C code as comments** above the Rust lines to prove fidelity.

**Required Format:**
```rust
// [C] ptr = buffer + i;  <-- Evidence
let ptr = unsafe { buffer.add(i) };