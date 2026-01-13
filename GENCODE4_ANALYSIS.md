# Gencode 4 Excess HSP Analysis

## Summary

LOSAT produces significantly more HSPs than NCBI BLAST when using genetic code 4 (Mold mitochondrial). The issue is **specific to gencode 4** and does NOT occur with gencode 1.

## Key Finding

| Metric | Gencode 1 | Gencode 4 |
|--------|-----------|-----------|
| LOSAT low-score HSPs (21-26 bits) | 175 | 1162 |
| NCBI low-score HSPs (21-26 bits) | 175 | 203 |
| Ratio | **1.00x** | **5.7x** |

**The excess is concentrated in low-score HSPs (21-26 bits).**

## Pipeline Stage Analysis (100k × 100k test)

Detailed comparison of each pipeline stage between gencode 1 and gencode 4:

| Pipeline Stage | Gencode 1 | Gencode 4 | Ratio |
|----------------|-----------|-----------|-------|
| **Lookup Table** | | | |
| Total entries | 1,108,192 | 1,734,018 | **1.56x** |
| Expansion factor | 8.51x | 12.84x | 1.51x |
| Neighbor entries | 1,083,016 | 1,714,019 | **1.58x** |
| Longest chain | 715 | 1,183 | 1.65x |
| **Seed Stage** | | | |
| K-mer matches | 23,331,988 | 31,823,985 | **1.36x** |
| Seeds to extension | 485,742 | 864,719 | **1.78x** |
| **Extension Stage** | | | |
| Total extensions | 485,742 | 864,719 | 1.78x |
| Ungapped HSPs | 3,577 | 8,860 | **2.48x** |
| **Final Output** | | | |
| E-value passed | 545 | 1,874 | **3.44x** |

**Observation:** The amplification grows at each stage:
- Lookup → Seeds: 1.56x → 1.78x
- Seeds → HSPs: 1.78x → 2.48x
- HSPs → Final: 2.48x → 3.44x

## Root Cause: Neighbor Word Generation

The fundamental difference is in the **lookup table construction**:

### Gencode 1 (TGA → Stop codon `*`)
- Stop codon has BLOSUM62 self-score = **+1**
- Very few neighbors pass threshold (only `*` can match `*` well)
- Words containing `*` contribute minimally to neighbor expansion

### Gencode 4 (TGA → Tryptophan `W`)
- Tryptophan has BLOSUM62 self-score = **+11**
- Many neighbors pass threshold (W can pair with F, Y, etc.)
- Words containing `W` generate significantly more neighbors

### Top K-mers Comparison

**Gencode 1 top k-mers:**
```
#1: count=715, residues 22,22,22 (SSS)
#2: count=675, residues 6,22,6  (ESE)
#3: count=648, residues 22,6,6  (SEE)
```

**Gencode 4 top k-mers:**
```
#1: count=1183, residues 20,6,20  (WEW) ← W-containing!
#2: count=1149, residues 6,20,20  (EWW) ← W-containing!
#3: count=1126, residues 20,20,6  (WWE) ← W-containing!
```

**All top 10 k-mers in gc4 contain Tryptophan (W)!**

## Truncation Test Results

Tests run on truncated AP027131 vs AP027133 sequences:

| Sequence Length | LOSAT gc1 | NCBI gc1 | Ratio gc1 | LOSAT gc4 | NCBI gc4 | Ratio gc4 |
|-----------------|-----------|----------|-----------|-----------|----------|-----------|
| 50k × 50k | 407 | 407 | **1.00x** | 850 | 464 | **1.83x** |
| 100k × 100k | 545 | 545 | **1.00x** | 1874 | 613 | **3.06x** |
| 200k × 200k | 1479 | 1479 | **1.00x** | 4911 | 1856 | **2.65x** |
| 300k × 300k | 4583 | 4533 | **1.01x** | 10874 | 5369 | **2.02x** |

## Critical Observation: NCBI vs LOSAT Amplification

| Genetic Code | LOSAT Hits | NCBI Hits | LOSAT Ratio | NCBI Ratio |
|--------------|------------|-----------|-------------|------------|
| gc1 | 545 | 545 | baseline | baseline |
| gc4 | 1874 | 613 | **3.44x** | **1.12x** |

NCBI only produces **12% more hits** with gc4 while LOSAT produces **244% more**. This is the smoking gun: NCBI has additional filtering that prevents the amplification seen in LOSAT.

## Hypothesis

Since NCBI BLAST with gc4 does NOT produce 3x more hits (only ~1.12x more than gc1), NCBI must have additional mechanisms:

1. **NCBI filtering that LOSAT is missing:**
   - NCBI may filter low-score HSPs more aggressively
   - May have composition-based filtering for W-rich sequences
   - May have hidden optimizations not in public source

2. **Lookup table construction difference:**
   - NCBI might cap neighbor expansion for high-frequency words
   - May have additional PV (position vector) filtering

3. **Two-hit detection difference:**
   - NCBI might suppress seeds differently for W-containing words

## Next Steps

1. **Compare NCBI lookup table size** - Does NCBI also have 1.56x more entries with gc4?
2. **Trace W-containing seeds** - Compare seed-to-extension ratio for W-rich regions
3. **Check NCBI composition filter** - Look for filtering based on amino acid frequency
4. **Examine NCBI's PV handling** - May filter high-frequency words differently

## Files Involved

- `LOSAT/src/algorithm/tblastx/lookup/backbone.rs` - Lookup table construction
- `LOSAT/src/algorithm/tblastx/blast_aascan.rs` - Seed scanning with PV
- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs` - Two-hit detection
- `LOSAT/src/algorithm/tblastx/sum_stats_linking/` - HSP linking

## Test Commands

```bash
# Build
cd LOSAT && cargo build --release

# Test with diagnostics (gc1)
LOSAT_DIAGNOSTICS=1 ./target/release/LOSAT tblastx \
    -q /tmp/AP027131_100k.fasta \
    -s /tmp/AP027133_100k.fasta \
    --query-gencode 1 --db-gencode 1 \
    -o /dev/null 2>&1 | grep -A50 "Phase 1"

# Test with diagnostics (gc4)
LOSAT_DIAGNOSTICS=1 ./target/release/LOSAT tblastx \
    -q /tmp/AP027131_100k.fasta \
    -s /tmp/AP027133_100k.fasta \
    --query-gencode 4 --db-gencode 4 \
    -o /dev/null 2>&1 | grep -A50 "Phase 1"
```

---
Generated: 2026-01-12
