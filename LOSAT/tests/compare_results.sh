#!/bin/bash
echo "=== Hit Count Comparison ==="
echo ""
for test in "AP027280.AP027280.tlosatx.n1" "AP027280.AP027280.tlosatx.n8" "MjeNMV.MelaMJNV.tlosatx.n1" "MelaMJNV.PemoMJNVA.tlosatx.n1" "PemoMJNVA.PeseMJNV.tlosatx.n1" "PeseMJNV.PemoMJNVB.tlosatx.n1" "PemoMJNVB.LvMJNV.tlosatx.n1" "LvMJNV.TrcuMJNV.tlosatx.n1" "TrcuMJNV.MellatMJNV.tlosatx.n1" "MellatMJNV.MeenMJNV.tlosatx.n1" "MeenMJNV.MejoMJNV.tlosatx.n1"; do
    losat_file="./losat_out/${test}.out"
    blast_file="./blast_out/${test}.out"
    
    if [ -f "$losat_file" ] && [ -f "$blast_file" ]; then
        losat_count=$(grep -v "^#" "$losat_file" | wc -l)
        blast_count=$(grep -v "^#" "$blast_file" | wc -l)
        diff=$((losat_count - blast_count))
        pct_diff=$(awk "BEGIN {printf \"%.2f\", ($diff / $blast_count) * 100}")
        echo "$test: LOSAT=$losat_count, NCBI=$blast_count, diff=$diff (${pct_diff}%)"
    fi
done
