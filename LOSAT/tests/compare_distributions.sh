#!/bin/bash
echo "=== Hit Distribution Comparison ==="
echo ""

compare_dist() {
    local test=$1
    local losat_file="./losat_out/${test}.out"
    local blast_file="./blast_out/${test}.out"
    
    if [ ! -f "$losat_file" ] || [ ! -f "$blast_file" ]; then
        return
    fi
    
    echo "--- $test ---"
    
    # Extract columns: length, evalue, bitscore, identity (assuming outfmt 6)
    # Format: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
    
    echo "LOSAT stats:"
    awk -F'\t' '!/^#/ && NF>=12 {
        len=$4; eval=$11; bit=$12; ident=$3
        sum_len+=len; sum_eval+=eval; sum_bit+=bit; sum_ident+=ident; n++
        if(len<min_len||n==1) min_len=len
        if(len>max_len||n==1) max_len=len
        if(eval<min_eval||n==1) min_eval=eval
        if(eval>max_eval||n==1) max_eval=eval
        if(bit<min_bit||n==1) min_bit=bit
        if(bit>max_bit||n==1) max_bit=bit
        if(ident<min_ident||n==1) min_ident=ident
        if(ident>max_ident||n==1) max_ident=ident
    } END {
        if(n>0) {
            printf "  Count: %d\n", n
            printf "  Length: min=%d, max=%d, avg=%.1f\n", min_len, max_len, sum_len/n
            printf "  E-value: min=%.2e, max=%.2e, avg=%.2e\n", min_eval, max_eval, sum_eval/n
            printf "  Bit-score: min=%.1f, max=%.1f, avg=%.1f\n", min_bit, max_bit, sum_bit/n
            printf "  Identity: min=%.1f%%, max=%.1f%%, avg=%.1f%%\n", min_ident, max_ident, sum_ident/n
        }
    }' "$losat_file"
    
    echo "NCBI stats:"
    awk -F'\t' '!/^#/ && NF>=12 {
        len=$4; eval=$11; bit=$12; ident=$3
        sum_len+=len; sum_eval+=eval; sum_bit+=bit; sum_ident+=ident; n++
        if(len<min_len||n==1) min_len=len
        if(len>max_len||n==1) max_len=len
        if(eval<min_eval||n==1) min_eval=eval
        if(eval>max_eval||n==1) max_eval=eval
        if(bit<min_bit||n==1) min_bit=bit
        if(bit>max_bit||n==1) max_bit=bit
        if(ident<min_ident||n==1) min_ident=ident
        if(ident>max_ident||n==1) max_ident=ident
    } END {
        if(n>0) {
            printf "  Count: %d\n", n
            printf "  Length: min=%d, max=%d, avg=%.1f\n", min_len, max_len, sum_len/n
            printf "  E-value: min=%.2e, max=%.2e, avg=%.2e\n", min_eval, max_eval, sum_eval/n
            printf "  Bit-score: min=%.1f, max=%.1f, avg=%.1f\n", min_bit, max_bit, sum_bit/n
            printf "  Identity: min=%.1f%%, max=%.1f%%, avg=%.1f%%\n", min_ident, max_ident, sum_ident/n
        }
    }' "$blast_file"
    echo ""
}

for test in "AP027280.AP027280.tlosatx.n1" "MjeNMV.MelaMJNV.tlosatx.n1" "MelaMJNV.PemoMJNVA.tlosatx.n1"; do
    compare_dist "$test"
done
