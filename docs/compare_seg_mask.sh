#!/bin/bash
# Compare SEG mask rates between LOSAT and NCBI segmasker
# Usage: ./compare_seg_mask.sh <fasta_file>

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fasta_file>"
    exit 1
fi

FASTA_FILE="$1"
OUTPUT_DIR="seg_mask_comparison"
mkdir -p "$OUTPUT_DIR"

echo "=== SEG Mask Rate Comparison ==="
echo "Input file: $FASTA_FILE"
echo ""

# Check if segmasker is available
if ! command -v segmasker &> /dev/null; then
    echo "ERROR: segmasker not found. Please install NCBI BLAST+ tools."
    echo "  On Ubuntu/Debian: sudo apt-get install ncbi-blast+"
    echo "  On macOS: brew install blast"
    exit 1
fi

# Run NCBI segmasker
echo "Running NCBI segmasker..."
segmasker -in "$FASTA_FILE" -outfmt interval -out "$OUTPUT_DIR/ncbi_segmasker.intervals" 2>&1 | tee "$OUTPUT_DIR/ncbi_segmasker.log"

# Calculate mask rate from intervals
if [ -f "$OUTPUT_DIR/ncbi_segmasker.intervals" ]; then
    echo ""
    echo "=== NCBI segmasker Results ==="
    TOTAL_BASES=$(grep -v "^>" "$FASTA_FILE" | tr -d '\n' | wc -c)
    MASKED_BASES=0
    while IFS=$'\t' read -r start end; do
        if [[ "$start" =~ ^[0-9]+$ ]] && [[ "$end" =~ ^[0-9]+$ ]]; then
            MASKED_BASES=$((MASKED_BASES + (end - start)))
        fi
    done < <(grep -v "^#" "$OUTPUT_DIR/ncbi_segmasker.intervals" | grep -v "^$")
    
    if [ "$TOTAL_BASES" -gt 0 ]; then
        MASK_RATE=$(echo "scale=2; $MASKED_BASES * 100 / $TOTAL_BASES" | bc)
        echo "Total bases: $TOTAL_BASES"
        echo "Masked bases: $MASKED_BASES"
        echo "Mask rate: ${MASK_RATE}%"
    fi
fi

echo ""
echo "=== LOSAT SEG Mask Rate ==="
echo "Run LOSAT with --seg flag and check the output for:"
echo "  'SEG masked X aa (Y%) of Z total aa across all frames'"
echo "  'Combined DNA masks: X bases (Y%)'"
echo ""
echo "Example command:"
echo "  cd LOSAT"
echo "  time ./target/release/losat tblastx \\"
echo "    -q $FASTA_FILE \\"
echo "    -s $FASTA_FILE \\"
echo "    --seg \\"
echo "    -o test.out 2>&1 | tee test.log"


