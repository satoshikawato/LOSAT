# Phase 1: 2-Bit Packed Nucleotide Format

## Implementation Overview

This phase implements a 2-bit packed nucleotide format for LOSAT, following NCBI BLAST's approach to compress DNA sequences and enable efficient k-mer extraction. The optimization reduces memory usage by 4x and enables O(1) sliding window k-mer extraction.

## Technical Details

### 2-Bit Encoding Specification

Following NCBI BLAST's encoding scheme (from `blast_encoding.c`):

| Base | 2-Bit Code | Binary |
|------|------------|--------|
| A    | 0          | 00     |
| C    | 1          | 01     |
| G    | 2          | 10     |
| T/U  | 3          | 11     |

Ambiguous bases (N, R, Y, etc.) are tracked separately in a bit vector and are not encoded in the packed data.

### Packing Order

Following NCBI BLAST's `BLAST_PackDNA` function, bases are packed with the first base in the most significant bits:

```
Byte layout: [base0: bits 6-7] [base1: bits 4-5] [base2: bits 2-3] [base3: bits 0-1]

Example: "ACGT" -> (0 << 6) | (1 << 4) | (2 << 2) | 3 = 0b00010111 = 23
```

### Compression Ratio

`COMPRESSION_RATIO = 4` (from NCBI BLAST's `blast_def.h`)

This means 4 nucleotides are stored per byte, achieving 4x memory reduction compared to ASCII representation.

### Bit Operations

#### Base Extraction
```rust
let byte_idx = pos / 4;
let bit_offset = 6 - (pos % 4) * 2;
let base = (data[byte_idx] >> bit_offset) & 0x03;
```

#### K-mer Extraction
K-mers are extracted as `u64` values with the first base in the most significant bits:
```rust
// For k-mer at position pos with length k:
let mut kmer: u64 = 0;
for i in 0..k {
    let base = self.get_base(pos + i);
    kmer = (kmer << 2) | (base as u64);
}
```

#### Sliding Window Optimization
New k-mers can be computed in O(1) from the previous k-mer:
```rust
let mask = (1u64 << (2 * k)) - 1;
let new_kmer = ((prev_kmer << 2) | new_base as u64) & mask;
```

This is the key optimization - instead of O(k) per position, we achieve O(1) per position.

## Code Changes

### New Files Created

1. **`src/sequence/mod.rs`**
   - Module declaration file exposing the packed_nucleotide module

2. **`src/sequence/packed_nucleotide.rs`**
   - Complete implementation of 2-bit packed nucleotide format
   - Key structures:
     - `PackedSequence`: Main struct holding packed data and ambiguous base tracking
     - `KmerIterator`: Efficient iterator over k-mers with sliding window optimization
   - Key functions:
     - `PackedSequence::new(seq: &[u8]) -> Option<Self>`: Create packed sequence from ASCII
     - `PackedSequence::get_base(pos: usize) -> u8`: Extract single base (2-bit code)
     - `PackedSequence::extract_kmer(pos: usize, k: usize) -> Option<u64>`: Extract k-mer
     - `PackedSequence::sliding_kmer(prev: u64, new_base: u8, k: usize) -> u64`: O(1) sliding window
     - `PackedSequence::iter_kmers(k: usize) -> KmerIterator`: Iterator over all k-mers
     - `encode_kmer_from_ascii(seq: &[u8], start: usize, k: usize) -> Option<u64>`: Compatibility function

### Modified Files

1. **`src/lib.rs`**
   - Added `pub mod sequence;` to expose the new sequence module

2. **`src/algorithm/blastn.rs`**
   - Added import: `use crate::sequence::packed_nucleotide::PackedSequence;`
   - Modified subject scanning loop (lines 3081-3285) to:
     - Pack subject sequence once at the start of processing
     - Use `KmerIterator` with sliding window optimization for k-mer extraction
     - Fall back to original `encode_kmer` for sequences that can't be packed

## NCBI BLAST Reference Code

The implementation references the following NCBI BLAST source files:

1. **`blast_encoding.c`** (lines 40-60): Base encoding tables
2. **`blast_util.c`** (`BLAST_PackDNA` function): Packing algorithm
3. **`blast_def.h`** (line ~100): `COMPRESSION_RATIO = 4` constant
4. **`blast_nascan.c`** (lines 95-134): K-mer extraction from packed format

## Performance Considerations

### Memory Usage
- Original: 1 byte per nucleotide
- Packed: 0.25 bytes per nucleotide (4x reduction)
- Additional overhead: ~1 bit per nucleotide for ambiguous base tracking

### K-mer Extraction Performance
- Original `encode_kmer`: O(k) per position (reads k bytes, performs k operations)
- Packed `iter_kmers`: O(1) per position after initialization (single shift and OR)

### Cache Efficiency
- Packed format improves cache utilization due to 4x smaller memory footprint
- Sequential access pattern in k-mer iterator is cache-friendly

## Known Limitations

1. **Ambiguous Base Handling**: K-mers containing ambiguous bases (N, R, Y, etc.) return `None` and are skipped. This matches NCBI BLAST behavior.

2. **Maximum K-mer Size**: K-mers are stored as `u64`, limiting k to 31 (62 bits for 2-bit encoding).

3. **Fallback Path**: Sequences that fail to pack (e.g., containing only ambiguous bases) fall back to the original `encode_kmer` approach.

## Testing

The module includes 14 comprehensive test cases covering:
- Base encoding/decoding
- Sequence packing/unpacking
- K-mer extraction
- Ambiguous base handling
- Sliding window optimization
- NCBI BLAST format compatibility
- Partial byte handling

Run tests with:
```bash
cd LOSAT && cargo test sequence::packed_nucleotide
```

## Next Phase Handoff

### Phase 2: Optimized Seed Generation and Lookup

The packed format foundation enables the following Phase 2 optimizations:

1. **Presence-Vector (PV) Checks**: Implement NCBI BLAST's PV array for fast lookup filtering (reference: `blast_nascan.c` lines 41-82)

2. **Query Packing**: Extend packed format to query sequences for lookup table construction

3. **Direct Lookup Optimization**: Use packed k-mers directly as array indices for O(1) lookup

### Files to Modify in Phase 2
- `src/algorithm/blastn.rs`: `build_lookup` and `build_direct_lookup` functions
- `src/seed/na_word_finder.rs`: K-mer encoding functions

### Key NCBI BLAST References for Phase 2
- `blast_nascan.c` lines 41-82: Presence-vector implementation
- `blast_lookup.c`: Lookup table construction with packed format
