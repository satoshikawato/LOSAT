# NEXT SESSION PROMPT (2026-01-02): Fix neighbor-map two-hit/diag parity + reduce HSP explosion (NCBI ground truth)

## （日本語要約・冒頭で必ず読む）

次セッションは **NCBI実装だけ**を参照して、LOSATのTBLASTX（特に`--neighbor-map`経路）の **two-hit/diag判定のズレ**を直し、`sum_stats_linking`に流れ込むHSP数の爆発を止める。

- **NCBI C/C++が仕様（ground truth）**
- **出力に影響する“簡略化”は禁止**（greedy overlap filter / ad-hoc suppression等はNG）
- 作業中は `docs/TBLASTX_NCBI_PARITY_DEVLOG.md` を **10分ごとに追記**

下の英語セクションは「そのまま実装に落とすためのチェックリスト」として残してある。

## Goal

Make LOSAT TBLASTX **match NCBI BLAST+ output** while eliminating the current “HSP explosion” that makes LOSAT slow (especially in `--neighbor-map` mode).

The immediate target is the `NZ_CP006932 self` benchmark (full 6x6, single-thread / multi-thread as needed) with:

- Output parity with NCBI (ideally identical; at minimum: **hit distribution within ±10%** for length/identity/bitscore/evalue and no pathological differences)
- No “simplification” that changes output (Strict Algorithmic Fidelity Protocol)

## Non‑negotiables (Strict Algorithmic Fidelity Protocol)

- **NCBI C/C++ implementation is ground truth**.
- **No “heuristic simplifications”** (e.g., greedy overlap filtering, ad-hoc hit suppression).
- Big refactors are allowed **only** if output/complexity matches NCBI behavior.
- Keep sentinel / encoding behavior NCBI-compatible.
- Update `docs/TBLASTX_NCBI_PARITY_DEVLOG.md` **at least every 10 minutes** while working.

## Reference (NCBI sources in this repo)

NCBI repo location: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast`

Key ground-truth files:

- Two-hit / diagonal logic: `c++/src/algo/blast/core/aa_ungapped.c`
- Subject scan loop: `c++/src/algo/blast/core/blast_aascan.c`
- Masked scanning offsets: `c++/src/algo/blast/core/masksubj.inl` (`s_DetermineScanningOffsets`)
- Sum-statistics HSP linking: `c++/src/algo/blast/core/link_hsps.c` (`s_BlastEvenGapLinkHSPs`)
- Lookup build + neighbor generation: `c++/src/algo/blast/core/blast_aalookup.c`

## Current LOSAT status (what is wrong)

### 1) neighbor-map mode is not NCBI-equivalent upstream

In LOSAT `--neighbor-map` scan loop (`src/algorithm/tblastx/utils.rs`), the two-hit/diag behavior diverges from NCBI:

- Missing exact NCBI `diff` rules (`diff >= window`, `diff < wordsize`, query range guard).
- Incorrect `last_hit` update semantics (NCBI updates `last_hit` only at specific points).
- `flag==1` “reset” should **not** attempt extension on the same hit (NCBI resets then moves on).
- A non-NCBI “one-hit ungapped extension fallback” exists when two-hit is absent (gated by an arbitrary `seed_score` threshold). This inflates ungapped HSP counts and makes linking dominate.

Result: `sum_stats_linking` sees too many HSPs → slow, and output deviates.

### 2) Subject masked range skipping is missing

NCBI does not scan masked intervals (uses `subject->seq_ranges` via `s_DetermineScanningOffsets`). LOSAT tends to scan full translated frames, increasing seed/HSP workload.

### 3) sum_stats_linking fast-path is incomplete

LOSAT has `next_larger` and some of `changed/linked_to` behavior, but lacks NCBI’s:

- `path_changed` / `use_current_max` reuse path (avoid recomputing when removal didn’t affect best paths)

This causes more recomputation than NCBI in practice.

## Work plan (do in this order)

### Task A — Make neighbor-map two-hit/diag match NCBI exactly

File: `LOSAT/src/algorithm/tblastx/utils.rs` (neighbor-map scan loop)

Implement NCBI logic from `aa_ungapped.c` **verbatim** (translated into Rust), including:

1. **`flag` block**:
   - If `flag==1` and `(subject_offset + diag_offset) < last_hit` → `continue`
   - Else set `last_hit = subject_offset + diag_offset`, set `flag = 0`, and **continue to next hit** (do not attempt extension on this hit)
2. **`flag==0` block**:
   - `last_hit_unoffset = last_hit - diag_offset`
   - `diff = subject_offset - last_hit_unoffset`
   - If `diff >= window`: set `last_hit = subject_offset + diag_offset`; `continue`
   - If `diff < wordsize`: `continue`
   - If `query_offset - diff < frame_base`: set `last_hit = subject_offset + diag_offset`; `continue`
   - Only then call two-hit ungapped extension (`s_BlastAaExtendTwoHit` equivalent) with NCBI’s arguments.
3. **After extension**:
   - If `right_extend`: set `flag=1`; set `last_hit = s_last_off - (wordsize - 1) + diag_offset` (NCBI)
   - Else update `last_hit` exactly as NCBI does (confirm in `aa_ungapped.c`)
4. **Remove the non-NCBI one-hit extension fallback** from neighbor-map mode.
   - If supporting `--window-size 0`, implement NCBI one-hit path explicitly (NCBI has separate behavior; don’t invent thresholds like `seed_score < 30`).

Acceptance check:

- `Total raw ungapped hits` in neighbor-map mode should drop substantially toward NCBI-like behavior (order-of-magnitude sanity).
- Output should move toward NCBI (hit count and length distribution).

### Task B — Implement subject masked scanning offsets (NCBI `seq_ranges`)

File: `LOSAT/src/algorithm/tblastx/utils.rs` (both normal scan and neighbor-map scan)

- Mirror NCBI’s `s_DetermineScanningOffsets` idea:
  - Derive unmasked intervals for each **subject translated frame** (SEG masked regions, possibly also DNA masking if applicable).
  - Scan only those intervals, not the full frame.
- This is a performance parity feature: fewer seeds/HSPs → linking becomes tractable.

### Task C — Finish NCBI `sum_stats_linking` fast-path (no heuristics)

File: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`

- Add **`path_changed` / `use_current_max`** logic exactly as in `link_hsps.c`.
- Ensure tie behavior matches NCBI (pay attention to NCBI’s “-1 to preserve ordering” trick in the large-gap loop).
- Keep output identical; do not add heuristic filters.

## Commands (repro)

LOSAT neighbor-map:

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests
timeout 300 ../target/release/LOSAT tblastx \
  -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta \
  -o ./losat_out/next_session_nm.out \
  --query-gencode 4 --db-gencode 4 -n 8 --seg --neighbor-map
```

LOSAT normal mode:

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests
timeout 300 ../target/release/LOSAT tblastx \
  -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta \
  -o ./losat_out/next_session_normal.out \
  --query-gencode 4 --db-gencode 4 -n 8 --seg
```

NCBI BLAST+ tblastx (example; adjust paths/options as used locally):

```bash
tblastx -query ./fasta/NZ_CP006932.fasta -subject ./fasta/NZ_CP006932.fasta \
  -seg yes -query_gencode 4 -db_gencode 4 -outfmt 6 \
  > ./losat_out/ncbi_tblastx.out
```

## Deliverables

- neighbor-map scan/two-hit parity fixes merged (no non-NCBI one-hit fallback)
- subject scanning offsets implemented (masked region skipping)
- `sum_stats_linking` has NCBI fast-path (`path_changed/use_current_max`)
- DEVLOG updated every 10 minutes during work


