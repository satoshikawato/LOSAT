#!/bin/bash
# Compare LOSAT TBLASTX results with NCBI TBLASTX results

compare_file() {
    local losat_file=$1
    local ncbi_file=$2
    local test_name=$3
    
    echo "=========================================="
    echo "Test: $test_name"
    echo "=========================================="
    
    # Count hits (skip header lines for NCBI BLAST)
    local losat_hits=$(grep -v "^#" "$losat_file" 2>/dev/null | wc -l | awk '{print $1}')
    local ncbi_hits=$(grep -v "^#" "$ncbi_file" 2>/dev/null | wc -l | awk '{print $1}')
    
    echo "Hit counts:"
    echo "  LOSAT: $losat_hits"
    echo "  NCBI:  $ncbi_hits"
    if [ "$ncbi_hits" -gt 0 ]; then
        local ratio=$(echo "scale=2; $losat_hits * 100 / $ncbi_hits" | bc)
        echo "  Ratio:  ${ratio}%"
    fi
    echo ""
    
    # Parse identity, length, and bit score
    if [ -f "$losat_file" ] && [ -f "$ncbi_file" ]; then
        echo "Identity distribution:"
        echo "  LOSAT (<50%): $(awk -F'\t' 'NR>0 && $3<50 {count++} END {print count+0}' "$losat_file")"
        echo "  LOSAT (50-80%): $(awk -F'\t' 'NR>0 && $3>=50 && $3<80 {count++} END {print count+0}' "$losat_file")"
        echo "  LOSAT (>=80%): $(awk -F'\t' 'NR>0 && $3>=80 {count++} END {print count+0}' "$losat_file")"
        
        echo "  NCBI (<50%): $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 && $3<50 {count++} END {print count+0}')"
        echo "  NCBI (50-80%): $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 && $3>=50 && $3<80 {count++} END {print count+0}')"
        echo "  NCBI (>=80%): $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 && $3>=80 {count++} END {print count+0}')"
        echo ""
        
        echo "Length statistics:"
        echo "  LOSAT min: $(awk -F'\t' 'NR>0 {if(min=="" || $4<min) min=$4} END {print min+0}' "$losat_file")"
        echo "  LOSAT max: $(awk -F'\t' 'NR>0 {if(max=="" || $4>max) max=$4} END {print max+0}' "$losat_file")"
        echo "  LOSAT avg: $(awk -F'\t' 'NR>0 {sum+=$4; count++} END {printf "%.1f\n", (count>0 ? sum/count : 0)}' "$losat_file")"
        
        echo "  NCBI min: $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 {if(min=="" || $4<min) min=$4} END {print min+0}')"
        echo "  NCBI max: $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 {if(max=="" || $4>max) max=$4} END {print max+0}')"
        echo "  NCBI avg: $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 {sum+=$4; count++} END {printf "%.1f\n", (count>0 ? sum/count : 0)}')"
        echo ""
        
        echo "Bit score statistics:"
        echo "  LOSAT min: $(awk -F'\t' 'NR>0 {if(min=="" || $12<min) min=$12} END {printf "%.1f\n", min+0}' "$losat_file")"
        echo "  LOSAT max: $(awk -F'\t' 'NR>0 {if(max=="" || $12>max) max=$12} END {printf "%.1f\n", max+0}' "$losat_file")"
        echo "  LOSAT avg: $(awk -F'\t' 'NR>0 {sum+=$12; count++} END {printf "%.1f\n", (count>0 ? sum/count : 0)}' "$losat_file")"
        
        echo "  NCBI min: $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 {if(min=="" || $12<min) min=$12} END {printf "%.1f\n", min+0}')"
        echo "  NCBI max: $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 {if(max=="" || $12>max) max=$12} END {printf "%.1f\n", max+0}')"
        echo "  NCBI avg: $(grep -v "^#" "$ncbi_file" | awk -F'\t' 'NR>0 {sum+=$12; count++} END {printf "%.1f\n", (count>0 ? sum/count : 0)}')"
        echo ""
    fi
}

# Compare each test case
compare_file "losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out" \
             "LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out" \
             "NZ_CP006932 self"

compare_file "losat_out/AP027132.NZ_CP006932.tlosatx.n1.out" \
             "LOSAT/tests/blast_out/AP027132.NZ_CP006932.tblastx.n1.out" \
             "AP027132 vs NZ_CP006932"

compare_file "losat_out/AP027078.AP027131.tlosatx.n1.out" \
             "LOSAT/tests/blast_out/AP027078.AP027131.tblastx.n1.out" \
             "AP027078 vs AP027131"

compare_file "losat_out/AP027131.AP027133.tlosatx.n1.out" \
             "LOSAT/tests/blast_out/AP027131.AP027133.tblastx.n1.out" \
             "AP027131 vs AP027133"

compare_file "losat_out/AP027133.AP027132.tlosatx.n1.out" \
             "LOSAT/tests/blast_out/AP027133.AP027132.tblastx.n1.out" \
             "AP027133 vs AP027132"

echo "=========================================="
echo "Summary"
echo "=========================================="

