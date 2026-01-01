# NEXT SESSION PROMPT (2026-01-02): Fix hit count disparity between LOSAT and NCBI

## （日本語要約・冒頭で必ず読む）

前セッションで **two-hit/diag ロジックを NCBI 準拠に修正**し、HSP 爆発問題を解決した（14.7M → 232K）。
しかし **最終 hit 数が NCBI の半分**（32K vs 62K）という差分が残っている。

- **NCBI C/C++が仕様（ground truth）**
- **出力に影響する"簡略化"は禁止**
- 作業中は `docs/TBLASTX_NCBI_PARITY_DEVLOG.md` を **10分ごとに追記**
- NCBI BLAST+ (C/C++) が ground truth。出力を一致させる。
- 出力に影響を与える"簡略化"は禁止。
- 副作用・計算量を減らせる場合は変更を許容。
- 大規模な変更も躊躇なく行え。NCBI BLAST+ (C/C++) の ground truthに合わせろ。

## Current Status

### What was fixed (previous session)

1. **neighbor-map two-hit/diag logic** - now matches NCBI `aa_ungapped.c` verbatim
2. **cutoff_score computation** - using proper E-value and GAP_TRIGGER_BIT_SCORE
3. **sum_stats_linking fast-path** - `path_changed/use_current_max` optimization

### Remaining disparity

| Metric | LOSAT | NCBI | Issue |
|--------|-------|------|-------|
| Final hits | 32,680 | 62,053 | LOSAT ~50% of NCBI |
| Max alignment | 1,967 AA | 2,260 AA | LOSAT 300 AA shorter |
| E<=1e-10 hits | 15,396 | 56,711 | NCBI 3.7x more |
| Min bit score | 22.1 | 22.1 | Same (good) |

## Root Cause Hypotheses

### Hypothesis 1: Linked set output difference

NCBI may output **all HSPs in a linked chain** separately, while LOSAT might be filtering some out.

- Check: `src/algorithm/tblastx/sum_stats_linking.rs` - how are linked HSPs output?
- NCBI reference: `link_hsps.c` lines 961-1000 - chain marking and output

### Hypothesis 2: Extension length difference

LOSAT's max alignment is 300 AA shorter than NCBI's. This could be:

- X-drop value difference
- Extension boundary handling
- Coordinate conversion bug

- Check: `src/algorithm/tblastx/extension.rs` - `extend_hit_two_hit` implementation
- NCBI reference: `aa_ungapped.c` - `s_BlastAaExtendTwoHit`

### Hypothesis 3: E-value calculation difference

NCBI has 3.7x more hits with E <= 1e-10. This suggests:

- Sum-statistics formula difference
- Search space calculation difference

- Check: `src/stats/sum_statistics.rs` 
- NCBI reference: `link_hsps.c` - sum-statistics E-value calculation

## Investigation Plan

### Step 1: Compare normal mode vs neighbor-map mode

First, verify if the disparity exists in both modes or only neighbor-map:

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests

# Normal mode
timeout 300 ../target/release/LOSAT tblastx \
  -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta \
  -o ./losat_out/normal_mode_check.out \
  --query-gencode 4 --db-gencode 4 -n 8 --seg

# Compare hit counts
wc -l ./losat_out/normal_mode_check.out
wc -l ./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out
```

If normal mode also shows disparity → issue is in shared code (extension/linking)
If normal mode is closer to NCBI → issue is neighbor-map specific

### Step 2: Debug extension length

Compare top hit coordinates:

- NCBI top hit: 537769-530990 (2260 AA)
- LOSAT top hit: 537382-531482 (1967 AA)

The positions are close but lengths differ. Check:
1. `extend_hit_two_hit` X-drop handling
2. `convert_coords` coordinate transformation

### Step 3: Debug linked set output

Add diagnostic output to see how many HSPs are in each linked chain:

```rust
// In sum_stats_linking.rs, after marking chain as removed:
eprintln!("Chain: {} HSPs, e-value: {}", chain_length, evalue);
```

## Key Files

### LOSAT
- `src/algorithm/tblastx/utils.rs` - scan/two-hit logic
- `src/algorithm/tblastx/extension.rs` - ungapped extension
- `src/algorithm/tblastx/sum_stats_linking.rs` - HSP linking
- `src/stats/sum_statistics.rs` - E-value calculation

### NCBI (ground truth)
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c` - extension functions
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c` - HSP linking and output

## Non-negotiables (Strict Algorithmic Fidelity Protocol)

- **NCBI C/C++ implementation is ground truth**
- **No heuristic simplifications** that change output
- Update `docs/TBLASTX_NCBI_PARITY_DEVLOG.md` **every 10 minutes**

## Success Criteria

- LOSAT hit count within ±10% of NCBI (55K-68K hits)
- Top alignment lengths match (±5%)
- E-value distribution similar to NCBI

## Commands

```bash
# LOSAT neighbor-map (current default for testing)
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests
../target/release/LOSAT tblastx \
  -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta \
  -o ./losat_out/debug.out \
  --query-gencode 4 --db-gencode 4 -n 8 --seg --neighbor-map

# NCBI reference (already generated)
# Output: ./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out
```

