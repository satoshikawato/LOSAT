#!/bin/bash
# Compare LOSAT vs NCBI BLAST+ for long sequences (600kb+) to investigate excessive hits
#
# This script runs both LOSAT and NCBI BLAST+ on long sequences and collects
# detailed debug output to identify the root cause of excessive hits.
#
# Usage: ./compare_long_sequences_debug.sh <query.fasta> <subject.fasta> [genetic_code]
#
# Example:
#   ./compare_long_sequences_debug.sh AP027131.fna AP027133.fna 4
#
# Requirements:
#   - LOSAT must be built and in PATH or LOSAT/target/release/losat
#   - NCBI BLAST+ tblastx must be installed and in PATH
#   - Python3 for analysis script

set -e

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <query.fasta> <subject.fasta> [genetic_code]"
    echo "Example: $0 AP027131.fna AP027133.fna 4"
    exit 1
fi

QUERY_FILE="$1"
SUBJECT_FILE="$2"
GENCODE="${3:-1}"

if [ ! -f "$QUERY_FILE" ]; then
    echo "Error: Query file not found: $QUERY_FILE"
    exit 1
fi

if [ ! -f "$SUBJECT_FILE" ]; then
    echo "Error: Subject file not found: $SUBJECT_FILE"
    exit 1
fi

# Get base name for output files
QUERY_BASE=$(basename "$QUERY_FILE" | sed 's/\.[^.]*$//')
SUBJECT_BASE=$(basename "$SUBJECT_FILE" | sed 's/\.[^.]*$//')
OUTDIR="debug_results_${QUERY_BASE}_vs_${SUBJECT_BASE}"
mkdir -p "$OUTDIR"

echo "============================================="
echo "Long sequence comparison: $QUERY_BASE vs $SUBJECT_BASE"
echo "Genetic code: $GENCODE"
echo "Output directory: $OUTDIR"
echo "============================================="

# Find LOSAT executable
if command -v losat &> /dev/null; then
    LOSAT_CMD="losat"
elif [ -f "LOSAT/target/release/losat" ]; then
    LOSAT_CMD="./LOSAT/target/release/losat"
elif [ -f "../LOSAT/target/release/losat" ]; then
    LOSAT_CMD="../LOSAT/target/release/losat"
else
    echo "Error: LOSAT executable not found. Please build with 'cargo build --release'"
    exit 1
fi

echo ""
echo "Using LOSAT: $LOSAT_CMD"
echo ""

# Get sequence lengths for reference
QUERY_LEN=$(grep -v "^>" "$QUERY_FILE" | tr -d '\n' | wc -c)
SUBJECT_LEN=$(grep -v "^>" "$SUBJECT_FILE" | tr -d '\n' | wc -c)
echo "Query length: $QUERY_LEN bp"
echo "Subject length: $SUBJECT_LEN bp"
echo ""

# Run LOSAT TBLASTX with debug output
echo "Running LOSAT TBLASTX (with debug output)..."
LOSAT_DIAGNOSTICS=1 $LOSAT_CMD tblastx \
    -query "$QUERY_FILE" \
    -subject "$SUBJECT_FILE" \
    --query_gencode "$GENCODE" \
    --db_gencode "$GENCODE" \
    -out "$OUTDIR/losat_output.tsv" \
    2>"$OUTDIR/losat_debug.log"

# Extract debug information from LOSAT log
echo "Extracting debug information from LOSAT output..."
grep -E "\[DEBUG (CUTOFF|CUTOFF_CALC|CUTOFF_UPDATE|LINKING_CUTOFF|LINKING_FILTER|HSP_SAVING)" "$OUTDIR/losat_debug.log" > "$OUTDIR/losat_debug_extracted.txt" || true

# Run NCBI BLAST+ TBLASTX if available
if command -v tblastx &> /dev/null; then
    echo "Running NCBI BLAST+ TBLASTX..."
    tblastx \
        -query "$QUERY_FILE" \
        -subject "$SUBJECT_FILE" \
        -query_gencode "$GENCODE" \
        -db_gencode "$GENCODE" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out "$OUTDIR/ncbi_output.tsv" \
        2>"$OUTDIR/ncbi_debug.log"
    
    # Extract effective search space and length adjustment from NCBI log
    echo "Extracting debug information from NCBI output..."
    grep -E "(Effective search space|Length adjustment|cutoff)" "$OUTDIR/ncbi_debug.log" > "$OUTDIR/ncbi_debug_extracted.txt" || true
else
    echo "NCBI BLAST+ not found. Skipping BLAST+ comparison."
    echo "Please install NCBI BLAST+ to compare results."
fi

echo ""
echo "============================================="
echo "Results Summary"
echo "============================================="

# Count hits and report statistics
echo ""
echo "Hit counts:"
if [ -f "$OUTDIR/losat_output.tsv" ]; then
    LOSAT_COUNT=$(wc -l < "$OUTDIR/losat_output.tsv")
    echo "  LOSAT: $LOSAT_COUNT hits"
fi
if [ -f "$OUTDIR/ncbi_output.tsv" ]; then
    NCBI_COUNT=$(wc -l < "$OUTDIR/ncbi_output.tsv")
    echo "  NCBI: $NCBI_COUNT hits"
    if [ -f "$OUTDIR/losat_output.tsv" ]; then
        DIFF=$((LOSAT_COUNT - NCBI_COUNT))
        RATIO=$(awk "BEGIN {printf \"%.2f\", $LOSAT_COUNT / $NCBI_COUNT}")
        echo "  Difference: $DIFF hits (LOSAT is ${RATIO}x NCBI)"
    fi
fi

# Report score distribution
echo ""
echo "Score distribution (bit score from column 12):"
if [ -f "$OUTDIR/losat_output.tsv" ] && [ -s "$OUTDIR/losat_output.tsv" ]; then
    echo "  LOSAT:"
    awk -F'\t' 'NR>0 {
        score = $12
        if (score < 30) low++
        else if (score < 50) mid++
        else high++
        total++
        sum += score
    }
    END {
        if (total > 0) {
            printf "    Total: %d, Low (<30): %d (%.1f%%), Mid (30-50): %d (%.1f%%), High (>50): %d (%.1f%%)\n",
                total, low, (low/total)*100, mid, (mid/total)*100, high, (high/total)*100
            printf "    Average bit score: %.1f\n", sum/total
        }
    }' "$OUTDIR/losat_output.tsv"
fi
if [ -f "$OUTDIR/ncbi_output.tsv" ] && [ -s "$OUTDIR/ncbi_output.tsv" ]; then
    echo "  NCBI:"
    awk -F'\t' 'NR>0 {
        score = $12
        if (score < 30) low++
        else if (score < 50) mid++
        else high++
        total++
        sum += score
    }
    END {
        if (total > 0) {
            printf "    Total: %d, Low (<30): %d (%.1f%%), Mid (30-50): %d (%.1f%%), High (>50): %d (%.1f%%)\n",
                total, low, (low/total)*100, mid, (mid/total)*100, high, (high/total)*100
            printf "    Average bit score: %.1f\n", sum/total
        }
    }' "$OUTDIR/ncbi_output.tsv"
fi

# Display debug information
echo ""
echo "============================================="
echo "Debug Information (LOSAT)"
echo "============================================="
if [ -f "$OUTDIR/losat_debug_extracted.txt" ] && [ -s "$OUTDIR/losat_debug_extracted.txt" ]; then
    cat "$OUTDIR/losat_debug_extracted.txt"
else
    echo "No debug information found in LOSAT output."
    echo "Check $OUTDIR/losat_debug.log for full output."
fi

# Display NCBI debug information if available
if [ -f "$OUTDIR/ncbi_debug_extracted.txt" ] && [ -s "$OUTDIR/ncbi_debug_extracted.txt" ]; then
    echo ""
    echo "============================================="
    echo "Debug Information (NCBI)"
    echo "============================================="
    cat "$OUTDIR/ncbi_debug_extracted.txt"
fi

echo ""
echo "============================================="
echo "Analysis complete. Output files in: $OUTDIR"
echo "============================================="
echo ""
echo "Next steps:"
echo "1. Review debug output in: $OUTDIR/losat_debug_extracted.txt"
echo "2. Compare cutoff values between LOSAT and NCBI"
echo "3. Check HSP saving statistics and filtering rates"
echo "4. Analyze score distributions to identify low-score hits"
echo ""

