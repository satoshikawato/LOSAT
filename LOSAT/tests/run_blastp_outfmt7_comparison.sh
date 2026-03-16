#!/bin/bash

set -euo pipefail

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgQuery("query");
# const string kArgOutput("out");
# const string kArgSubject("subject");
# const string kArgEvalue("evalue");
# const string kArgMaxTargetSequences("max_target_seqs");
# ```
#
# NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:37
# ```c
# const int kDfltArgOutputFormat = 0;
# ```

LOSAT_BIN="../target/debug/LOSAT"
REF_DIR="./blast_out/blastp_outfmt7"
FASTA_DIR="./fasta"

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgOutput("out");
# ```
#
# Comparison-harness note: keep LOSAT outputs in a temporary directory so
# tracked `losat_out/` artifacts cannot mask current-source divergence.
OUT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/losat-blastp-outfmt7-XXXXXX")"
trap 'rm -rf "${OUT_DIR}"' EXIT

run_case() {
    local query="$1"
    local subject="$2"
    local stem="$3"
    local ref_file="${REF_DIR}/${stem}.out"

    if [ ! -f "${ref_file}" ]; then
        echo "Missing blastp outfmt 7 reference: ${ref_file}" >&2
        exit 1
    fi

    "${LOSAT_BIN}" blastp \
        -query "${FASTA_DIR}/${query}" \
        -subject "${FASTA_DIR}/${subject}" \
        -outfmt 7 \
        -evalue 1e-5 \
        -max_target_seqs 10 \
        -out "${OUT_DIR}/${stem}.out" \
        > "${OUT_DIR}/${stem}.log" 2>&1

    diff -u "${ref_file}" "${OUT_DIR}/${stem}.out"
}

run_case "WSSV.faa" "PajaWSV.faa" "WSSV.PajaWSV.blastp"
run_case "PajaWSV.faa" "SicyWSV.faa" "PajaWSV.SicyWSV.blastp"
run_case "WSSV.faa" "SicyWSV.faa" "WSSV.SicyWSV.blastp"
run_case "WSSV.faa" "CoBV.faa" "WSSV.CoBV.blastp"
run_case "CoBV.faa" "PajaWSV.faa" "CoBV.PajaWSV.blastp"
run_case "CoBV.faa" "SicyWSV.faa" "CoBV.SicyWSV.blastp"
run_case "AP027078.faa" "AP027131.faa" "AP027078.AP027131.blastp"
run_case "AP027078.faa" "AP027132.faa" "AP027078.AP027132.blastp"
run_case "AP027078.faa" "AP027133.faa" "AP027078.AP027133.blastp"
run_case "AP027078.faa" "NZ_CP006932.faa" "AP027078.NZ_CP006932.blastp"
run_case "AP027131.faa" "AP027132.faa" "AP027131.AP027132.blastp"
run_case "AP027131.faa" "AP027133.faa" "AP027131.AP027133.blastp"
run_case "AP027131.faa" "NZ_CP006932.faa" "AP027131.NZ_CP006932.blastp"
run_case "AP027132.faa" "AP027133.faa" "AP027132.AP027133.blastp"
run_case "AP027132.faa" "NZ_CP006932.faa" "AP027132.NZ_CP006932.blastp"
run_case "AP027133.faa" "NZ_CP006932.faa" "AP027133.NZ_CP006932.blastp"
