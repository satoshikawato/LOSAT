#!/bin/bash
# Self-comparison test for TBLASTX: Compare LOSAT vs NCBI BLAST+ output
#
# This script runs both LOSAT and NCBI BLAST+ on the same self-comparison
# (query == subject) and compares the results.
#
# Usage: ./compare_self_tblastx.sh <fasta_file> [genetic_code]
#
# Example:
#   ./compare_self_tblastx.sh NZ_CP006932.fna 4
#
# Requirements:
#   - LOSAT must be built and in PATH or LOSAT/target/release/losat
#   - NCBI BLAST+ tblastx must be installed and in PATH
#   - Python3 for analysis script

set -e

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <fasta_file> [genetic_code]"
    echo "Example: $0 NZ_CP006932.fna 4"
    exit 1
fi

FASTA_FILE="$1"
GENCODE="${2:-1}"

if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: File not found: $FASTA_FILE"
    exit 1
fi

# Get base name for output files
BASENAME=$(basename "$FASTA_FILE" | sed 's/\.[^.]*$//')
OUTDIR="test_results_${BASENAME}"
mkdir -p "$OUTDIR"

echo "============================================="
echo "Self-comparison test: $FASTA_FILE"
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

# Run LOSAT TBLASTX (standard mode)
echo "Running LOSAT TBLASTX (standard mode)..."
LOSAT_DIAGNOSTICS=1 $LOSAT_CMD tblastx \
    -query "$FASTA_FILE" \
    -subject "$FASTA_FILE" \
    --query_gencode "$GENCODE" \
    --db_gencode "$GENCODE" \
    -out "$OUTDIR/losat_standard.tsv" \
    2>"$OUTDIR/losat_standard.log"

# Run LOSAT TBLASTX (NCBI-compatible mode)
echo "Running LOSAT TBLASTX (NCBI-compat mode)..."
LOSAT_DIAGNOSTICS=1 $LOSAT_CMD tblastx \
    -query "$FASTA_FILE" \
    -subject "$FASTA_FILE" \
    --query_gencode "$GENCODE" \
    --db_gencode "$GENCODE" \
    --ncbi_compat \
    -out "$OUTDIR/losat_ncbi_compat.tsv" \
    2>"$OUTDIR/losat_ncbi_compat.log"

# Run NCBI BLAST+ TBLASTX if available
if command -v tblastx &> /dev/null; then
    echo "Running NCBI BLAST+ TBLASTX..."
    tblastx \
        -query "$FASTA_FILE" \
        -subject "$FASTA_FILE" \
        -query_gencode "$GENCODE" \
        -db_gencode "$GENCODE" \
        -outfmt 6 \
        -out "$OUTDIR/blast_plus.tsv" \
        2>"$OUTDIR/blast_plus.log"
else
    echo "NCBI BLAST+ not found. Skipping BLAST+ comparison."
fi

echo ""
echo "============================================="
echo "Results Summary"
echo "============================================="

# Count hits and report statistics
echo ""
echo "Hit counts:"
for f in "$OUTDIR"/*.tsv; do
    if [ -f "$f" ]; then
        COUNT=$(wc -l < "$f")
        echo "  $(basename "$f"): $COUNT hits"
    fi
done

# Report alignment length statistics
echo ""
echo "Alignment length statistics (from column 4 = alignment length):"
for f in "$OUTDIR"/*.tsv; do
    if [ -f "$f" ] && [ -s "$f" ]; then
        echo "  $(basename "$f"):"
        # Min, max, average alignment length
        awk -F'\t' 'NR==1 {min=$4; max=$4; sum=$4; count=1; next}
                    {if($4<min) min=$4; if($4>max) max=$4; sum+=$4; count++}
                    END {if(count>0) printf "    Min: %d, Max: %d, Avg: %.1f\n", min, max, sum/count}' "$f"
    fi
done

# Check for very long alignments (potential self-comparison issue)
echo ""
echo "Long alignments (>10000 bp):"
for f in "$OUTDIR"/*.tsv; do
    if [ -f "$f" ] && [ -s "$f" ]; then
        LONG_COUNT=$(awk -F'\t' '$4 > 10000 {count++} END {print count+0}' "$f")
        if [ "$LONG_COUNT" -gt 0 ]; then
            echo "  $(basename "$f"): $LONG_COUNT alignments > 10kb"
        fi
    fi
done

# Display diagnostic logs if available
echo ""
echo "============================================="
echo "Diagnostic Logs"
echo "============================================="
for f in "$OUTDIR"/*.log; do
    if [ -f "$f" ] && [ -s "$f" ]; then
        echo ""
        echo "--- $(basename "$f") ---"
        # Show diagnostics section
        grep -A 100 "Pipeline Diagnostics" "$f" 2>/dev/null || echo "(No diagnostics found)"
    fi
done

echo ""
echo "============================================="
echo "Test complete. Output files in: $OUTDIR"
echo "============================================="

