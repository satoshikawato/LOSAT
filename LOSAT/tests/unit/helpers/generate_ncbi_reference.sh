#!/bin/bash
# Generate NCBI BLAST reference data for unit tests
#
# This script runs NCBI BLAST with specific parameters to generate
# reference data for comparing LOSAT's implementation with NCBI BLAST.
#
# Usage:
#   ./generate_ncbi_reference.sh <query.fasta> <subject.fasta> <output_prefix>
#
# Requirements:
#   - blastn or tblastx must be in PATH
#   - Input FASTA files must exist

set -e

if [ $# -lt 3 ]; then
    echo "Usage: $0 <query.fasta> <subject.fasta> <output_prefix> [blastn|tblastx]"
    echo ""
    echo "Example:"
    echo "  $0 test_query.fasta test_subject.fasta test_case blastn"
    exit 1
fi

QUERY_FASTA="$1"
SUBJECT_FASTA="$2"
OUTPUT_PREFIX="$3"
BLAST_MODE="${4:-blastn}"

if [ ! -f "$QUERY_FASTA" ]; then
    echo "Error: Query file not found: $QUERY_FASTA"
    exit 1
fi

if [ ! -f "$SUBJECT_FASTA" ]; then
    echo "Error: Subject file not found: $SUBJECT_FASTA"
    exit 1
fi

# Check if blastn or tblastx is available
if ! command -v "$BLAST_MODE" &> /dev/null; then
    echo "Error: $BLAST_MODE not found in PATH"
    echo "Please install NCBI BLAST+ or add it to your PATH"
    exit 1
fi

OUTPUT_DIR="ncbi_reference_data"
mkdir -p "$OUTPUT_DIR"

OUTPUT_FILE="$OUTPUT_DIR/${OUTPUT_PREFIX}.out"
LOG_FILE="$OUTPUT_DIR/${OUTPUT_PREFIX}.log"

echo "Running NCBI $BLAST_MODE..."
echo "  Query: $QUERY_FASTA"
echo "  Subject: $SUBJECT_FASTA"
echo "  Output: $OUTPUT_FILE"

# Run NCBI BLAST with specific parameters for reference data generation
# Important flags:
#   -dust no: Disable DUST filtering (BLASTN)
#   -seg no: Disable SEG filtering (TBLASTX)
#   -comp_based_stats 0: Disable composition-based statistics
#   -reward 1 -penalty -2: Explicit scoring parameters
#   -gapopen 0 -gapextend 0: Explicit gap costs
#   -outfmt 7: Detailed output format with comments

if [ "$BLAST_MODE" = "blastn" ]; then
    "$BLAST_MODE" \
        -query "$QUERY_FASTA" \
        -subject "$SUBJECT_FASTA" \
        -dust no \
        -reward 1 \
        -penalty -2 \
        -gapopen 0 \
        -gapextend 0 \
        -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out "$OUTPUT_FILE" \
        2>&1 | tee "$LOG_FILE"
elif [ "$BLAST_MODE" = "tblastx" ]; then
    "$BLAST_MODE" \
        -query "$QUERY_FASTA" \
        -subject "$SUBJECT_FASTA" \
        -seg no \
        -comp_based_stats 0 \
        -gapopen 11 \
        -gapextend 1 \
        -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out "$OUTPUT_FILE" \
        2>&1 | tee "$LOG_FILE"
else
    echo "Error: Unsupported BLAST mode: $BLAST_MODE"
    echo "Supported modes: blastn, tblastx"
    exit 1
fi

echo ""
echo "NCBI BLAST output saved to: $OUTPUT_FILE"
echo "Log saved to: $LOG_FILE"
echo ""
echo "Next steps:"
echo "  1. Extract effective search space from the output:"
echo "     grep 'Effective search space' $OUTPUT_FILE"
echo ""
echo "  2. Extract query/database lengths from the output:"
echo "     grep -E 'Query|Database' $OUTPUT_FILE"
echo ""
echo "  3. Generate test cases using:"
echo "     python extract_ncbi_cases_enhanced.py $OUTPUT_FILE <q_len> <db_len> <db_num_seqs> $BLAST_MODE"


