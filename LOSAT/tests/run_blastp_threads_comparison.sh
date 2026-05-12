#!/bin/bash

set -uo pipefail

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgQuery("query");
# const string kArgOutput("out");
# const string kArgSubject("subject");
# const string kArgNumThreads("num_threads");
# const string kArgEvalue("evalue");
# const string kArgMaxTargetSequences("max_target_seqs");
# ```
#
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
# ```c
# const int kMaxValue = static_cast<int>(CSystemInfo::GetCpuCount());
# int num_threads = args[kArgNumThreads].AsInteger();
# if (num_threads > kMaxValue) {
#     m_NumThreads = kMaxValue;
# } else {
#     m_NumThreads = num_threads;
# }
# ```

LOSAT_BIN="${LOSAT_BIN:-../target/debug/LOSAT}"
case "$(uname -s 2>/dev/null || true)" in
    MINGW*|MSYS*|CYGWIN*)
        if [ -x "${LOSAT_BIN}.exe" ]; then
            LOSAT_BIN="${LOSAT_BIN}.exe"
        fi
        ;;
esac
if [ ! -x "${LOSAT_BIN}" ] && [ -x "${LOSAT_BIN}.exe" ]; then
    LOSAT_BIN="${LOSAT_BIN}.exe"
fi

THREADS="${LOSAT_BLASTP_THREADS:-8}"
REF_DIR="./blast_out/blastp"
FASTA_DIR="./fasta"
source ./blastp_parity_options.sh

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgOutput("out");
# ```
#
# Comparison-harness note: keep LOSAT outputs in a temporary directory so
# tracked `losat_out/` artifacts cannot mask current-source divergence.
OUT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/losat-blastp-threads-XXXXXX")"
trap 'rm -rf "${OUT_DIR}"' EXIT

serial_thread_failures=()
ncbi_failures=()

run_case() {
    local query="$1"
    local subject="$2"
    local stem="$3"
    local serial_out="${OUT_DIR}/${stem}.n1.out"
    local serial_log="${OUT_DIR}/${stem}.n1.log"
    local threaded_out="${OUT_DIR}/${stem}.n${THREADS}.out"
    local threaded_log="${OUT_DIR}/${stem}.n${THREADS}.log"
    local serial_thread_diff="${OUT_DIR}/${stem}.serial-vs-threaded.diff"
    local ncbi_diff="${OUT_DIR}/${stem}.threaded-vs-ncbi.diff"

    "${LOSAT_BIN}" blastp \
        -query "${FASTA_DIR}/${query}" \
        -subject "${FASTA_DIR}/${subject}" \
        -outfmt 6 \
        "${BLASTP_PARITY_STRICT_OPTIONS[@]}" \
        -num_threads 1 \
        -out "${serial_out}" \
        > "${serial_log}" 2>&1
    local status=$?
    if [ "${status}" -ne 0 ]; then
        echo "LOSATP serial failed for ${stem}; see ${serial_log}" >&2
        serial_thread_failures+=("${stem}:serial-exit-${status}")
        return
    fi

    "${LOSAT_BIN}" blastp \
        -query "${FASTA_DIR}/${query}" \
        -subject "${FASTA_DIR}/${subject}" \
        -outfmt 6 \
        "${BLASTP_PARITY_STRICT_OPTIONS[@]}" \
        -num_threads "${THREADS}" \
        -out "${threaded_out}" \
        > "${threaded_log}" 2>&1
    status=$?
    if [ "${status}" -ne 0 ]; then
        echo "LOSATP threaded failed for ${stem}; see ${threaded_log}" >&2
        serial_thread_failures+=("${stem}:threaded-exit-${status}")
        return
    fi

    if ! diff -u "${serial_out}" "${threaded_out}" > "${serial_thread_diff}"; then
        echo "serial/threaded mismatch for ${stem}; see ${serial_thread_diff}" >&2
        serial_thread_failures+=("${stem}")
    fi

    if ! diff -u "${REF_DIR}/${stem}.out" "${threaded_out}" > "${ncbi_diff}"; then
        echo "threaded/NCBI mismatch for ${stem}; see ${ncbi_diff}" >&2
        ncbi_failures+=("${stem}")
    fi
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

echo "serial/threaded mismatches: ${#serial_thread_failures[@]}"
if [ "${#serial_thread_failures[@]}" -gt 0 ]; then
    printf '  %s\n' "${serial_thread_failures[@]}"
fi

echo "threaded/NCBI mismatches: ${#ncbi_failures[@]}"
if [ "${#ncbi_failures[@]}" -gt 0 ]; then
    printf '  %s\n' "${ncbi_failures[@]}"
fi

if [ "${#serial_thread_failures[@]}" -gt 0 ] || [ "${#ncbi_failures[@]}" -gt 0 ]; then
    exit 1
fi
