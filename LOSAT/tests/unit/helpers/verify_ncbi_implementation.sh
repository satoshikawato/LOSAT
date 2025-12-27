#!/bin/bash
# Verify LOSAT's implementation against NCBI BLAST
#
# This script:
# 1. Runs NCBI BLAST to generate reference data
# 2. Extracts effective search space and length adjustment
# 3. Runs LOSAT's calculation
# 4. Compares the results
#
# Usage:
#   ./verify_ncbi_implementation.sh <query.fasta> <subject.fasta> <test_name>

set -e

if [ $# -lt 3 ]; then
    echo "Usage: $0 <query.fasta> <subject.fasta> <test_name> [blastn|tblastx]"
    exit 1
fi

QUERY_FASTA="$1"
SUBJECT_FASTA="$2"
TEST_NAME="$3"
BLAST_MODE="${4:-blastn}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Step 1: Generate NCBI BLAST reference data
echo "Step 1: Generating NCBI BLAST reference data..."
./generate_ncbi_reference.sh "$QUERY_FASTA" "$SUBJECT_FASTA" "$TEST_NAME" "$BLAST_MODE"

OUTPUT_FILE="ncbi_reference_data/${TEST_NAME}.out"

# Step 2: Extract information from NCBI BLAST output
echo ""
echo "Step 2: Extracting information from NCBI BLAST output..."

# Extract effective search space
EFFECTIVE_SPACE=$(grep "Effective search space" "$OUTPUT_FILE" | head -1 | sed -n 's/.*Effective search space: \([0-9.]*\).*/\1/p' || echo "")

# Extract length adjustment
LENGTH_ADJ=$(grep "Length adjustment" "$OUTPUT_FILE" | head -1 | sed -n 's/.*Length adjustment: \([0-9.]*\).*/\1/p' || echo "")

# Extract query length - try multiple patterns
QUERY_LEN=$(grep "^# Query:" "$OUTPUT_FILE" | head -1 | sed -n 's/.*Query: [0-9]* letters (\([0-9]*\) letters).*/\1/p')
if [ -z "$QUERY_LEN" ]; then
    QUERY_LEN=$(grep "^# Query:" "$OUTPUT_FILE" | head -1 | sed -n 's/.*Query: \([0-9]*\) letters.*/\1/p')
fi
if [ -z "$QUERY_LEN" ]; then
    # Try to get from sequence file directly
    QUERY_LEN=$(grep -v "^>" "$QUERY_FASTA" | tr -d '\n' | wc -c)
fi

# Extract database length - try multiple patterns
DB_LEN=$(grep "^# Database:" "$OUTPUT_FILE" | head -1 | sed -n 's/.*Database: [0-9]* sequences; \([0-9]*\) total letters.*/\1/p')
if [ -z "$DB_LEN" ]; then
    DB_LEN=$(grep "^# Database:" "$OUTPUT_FILE" | head -1 | sed -n 's/.*Database:.*\([0-9]*\) total letters.*/\1/p')
fi
if [ -z "$DB_LEN" ]; then
    # Try to get from sequence file directly
    DB_LEN=$(grep -v "^>" "$SUBJECT_FASTA" | tr -d '\n' | wc -c)
fi

# Extract number of sequences
DB_NUM_SEQS=$(grep "^# Database:" "$OUTPUT_FILE" | head -1 | sed -n 's/.*Database: \([0-9]*\) sequences.*/\1/p')
if [ -z "$DB_NUM_SEQS" ]; then
    DB_NUM_SEQS=1  # Default to 1 for single-subject searches
fi

echo "  Query length: $QUERY_LEN"
echo "  Database length: $DB_LEN"
echo "  Number of sequences: $DB_NUM_SEQS"
echo "  Effective search space: $EFFECTIVE_SPACE"
echo "  Length adjustment: $LENGTH_ADJ"

if [ -z "$EFFECTIVE_SPACE" ] || [ -z "$QUERY_LEN" ] || [ -z "$DB_LEN" ]; then
    echo ""
    echo "Warning: Could not extract all required information from NCBI BLAST output"
    echo "Please check the output file: $OUTPUT_FILE"
    echo ""
    echo "Available information:"
    head -50 "$OUTPUT_FILE" | grep -E "^#"
    exit 1
fi

# Step 3: Generate test case
echo ""
echo "Step 3: Generating test case..."

# Create a temporary Rust test file
TEST_CASE_FILE="ncbi_reference_data/${TEST_NAME}_test_case.rs"
cat > "$TEST_CASE_FILE" << EOF
// Test case generated from NCBI BLAST output
// Query: $QUERY_FASTA
// Subject: $SUBJECT_FASTA
// NCBI BLAST mode: $BLAST_MODE

use LOSAT::stats::search_space::SearchSpace;
use LOSAT::stats::tables::KarlinParams;

#[test]
#[ignore] // Generated test case - enable after verification
fn test_${TEST_NAME}_effective_search_space() {
    // Parameters from NCBI BLAST output
    let query_len: usize = $QUERY_LEN;
    let db_len: usize = $DB_LEN;
    let db_num_seqs: usize = $DB_NUM_SEQS;
    
    // TODO: Set correct Karlin-Altschul parameters based on scoring scheme
    // For BLASTN with reward=1, penalty=-2, gapopen=0, gapextend=0:
    let params = KarlinParams {
        lambda: 1.28,  // TODO: Verify from NCBI BLAST output
        k: 0.46,       // TODO: Verify from NCBI BLAST output
        h: 0.85,       // TODO: Verify from NCBI BLAST output
        alpha: 1.5,    // TODO: Verify from NCBI BLAST output
        beta: -2.0,    // TODO: Verify from NCBI BLAST output
    };
    
    // Calculate using LOSAT
    let search_space = SearchSpace::for_database_search(
        query_len,
        db_len,
        db_num_seqs,
        &params,
        true, // use_length_adjustment
    );
    
    // Expected values from NCBI BLAST
    let expected_effective_space: f64 = $EFFECTIVE_SPACE;
    let expected_length_adj: i64 = ${LENGTH_ADJ%.*}; // Convert to integer
    
    println!("LOSAT effective space: {}", search_space.effective_space);
    println!("NCBI effective space: {}", expected_effective_space);
    println!("Difference: {:.2}%", 
        ((search_space.effective_space - expected_effective_space).abs() / expected_effective_space) * 100.0);
    
    println!("LOSAT length adjustment: {}", search_space.length_adjustment);
    println!("NCBI length adjustment: {}", expected_length_adj);
    println!("Difference: {}", (search_space.length_adjustment - expected_length_adj).abs());
    
    // Verify with 0.1% tolerance for effective search space
    let relative_diff = (search_space.effective_space - expected_effective_space).abs() / expected_effective_space;
    assert!(
        relative_diff <= 0.001,
        "Effective search space mismatch: LOSAT={}, NCBI={}, diff={:.4}%",
        search_space.effective_space,
        expected_effective_space,
        relative_diff * 100.0
    );
    
    // Verify length adjustment (should match exactly or ±1)
    let adj_diff = (search_space.length_adjustment - expected_length_adj).abs();
    assert!(
        adj_diff <= 1,
        "Length adjustment mismatch: LOSAT={}, NCBI={}, diff={}",
        search_space.length_adjustment,
        expected_length_adj,
        adj_diff
    );
}
EOF

echo "  Test case saved to: $TEST_CASE_FILE"

# Step 4: Instructions
echo ""
echo "Step 4: Next steps"
echo "  1. Update Karlin-Altschul parameters in the test case"
echo "     (Extract lambda, K, H, alpha, beta from NCBI BLAST output)"
echo ""
echo "  2. Run the test:"
echo "     cargo test --test unit_tests test_${TEST_NAME}_effective_search_space -- --nocapture"
echo ""
echo "  3. If the test passes, integrate it into the main test suite"

