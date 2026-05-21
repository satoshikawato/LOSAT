#!/bin/bash

# --- Setup ---

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Create output directories
mkdir -p losat_out blast_out

# Path to LOSAT binary
LOSAT_BIN="../target/release/LOSAT"
LOSAT_THREADS="${LOSAT_THREADS:-${LOSATP_THREADS:-${LOSAT_BLASTP_THREADS:-8}}}"
LOSATP_THREADS="${LOSATP_THREADS:-${LOSAT_THREADS}}"
LOSAT_WASM_BIN="${LOSAT_WASM_BIN:-../target/wasm32-wasip1/release/LOSAT.wasm}"
LOSAT_WASM_RUNNER="./losat_out/run_losat_wasi.js"
RUN_LOSAT_WASM="${RUN_LOSAT_WASM:-1}"
BUILD_LOSAT_WASM="${BUILD_LOSAT_WASM:-0}"
LOSAT_WASM_TARGET="${LOSAT_WASM_TARGET:-wasm32-wasip1}"
LOSAT_WASM_FEATURES="${LOSAT_WASM_FEATURES:---no-default-features}"
LOSAT_WASM_EXECUTION_MODE="${LOSAT_WASM_EXECUTION_MODE:-command-wasi-serial}"
LOSAT_WASM_THREADS_VERIFIED="${LOSAT_WASM_THREADS_VERIFIED:-0}"
TBLASTX_BIN="${TBLASTX_BIN:-tblastx}"
FASTA_DIR="${SCRIPT_DIR}/fasta"
LOSAT_OUT_DIR="${SCRIPT_DIR}/losat_out"

case "${LOSAT_WASM_BIN}" in
    /*) ;;
    *) LOSAT_WASM_BIN="${SCRIPT_DIR}/${LOSAT_WASM_BIN}" ;;
esac

wasm_has_start_export() {
    local wasm_path="$1"

    [ -f "${wasm_path}" ] || return 1
    node -e 'const fs=require("fs"); const p=process.argv[1]; const m=new WebAssembly.Module(fs.readFileSync(p)); process.exit(WebAssembly.Module.exports(m).some((e)=>e.name==="_start") ? 0 : 1);' "${wasm_path}" >/dev/null 2>&1
}

resolve_losat_wasm_bin() {
    local requested="$1"
    local release_dir
    local deps_dir
    local newest=""
    local candidate

    if wasm_has_start_export "${requested}"; then
        printf '%s\n' "${requested}"
        return 0
    fi

    release_dir="$(dirname "${requested}")"
    deps_dir="${release_dir}/deps"
    if [ -d "${deps_dir}" ]; then
        for candidate in "${deps_dir}"/LOSAT-*.wasm; do
            [ -f "${candidate}" ] || continue
            if wasm_has_start_export "${candidate}"; then
                if [ -z "${newest}" ] || [ "${candidate}" -nt "${newest}" ]; then
                    newest="${candidate}"
                fi
            fi
        done
    fi

    if [ -n "${newest}" ]; then
        printf '%s\n' "${newest}"
        return 0
    fi

    return 1
}

wasm_requested_num_threads() {
    local expect_value=0
    local arg

    for arg in "$@"; do
        if [ "${expect_value}" = "1" ]; then
            printf '%s\n' "${arg}"
            return 0
        fi
        case "${arg}" in
            -n|-num_threads|--num_threads|--num-threads)
                expect_value=1
                ;;
            --num_threads=*|--num-threads=*)
                printf '%s\n' "${arg#*=}"
                return 0
                ;;
        esac
    done

    printf 'unspecified\n'
}

wasm_rayon_compiled() {
    case " ${LOSAT_WASM_FEATURES} " in
        *" parallel "*|*"--features parallel"*|*"--features=parallel"*)
            printf 'true\n'
            ;;
        *)
            printf 'false\n'
            ;;
    esac
}

# NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
# ```c
# db_length = BlastSeqSrcGetTotLen(seq_src);
# itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
# while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
#     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
#         continue;
#     }
# }
# ```
wasm_effective_engine_threads() {
    case "${LOSAT_WASM_EXECUTION_MODE}" in
        command-wasi-serial|browser-in-memory-serial)
            printf '1\n'
            ;;
        browser-worker-parallel|future-wasi-threaded)
            if [ "${LOSAT_WASM_THREADS_VERIFIED}" = "1" ]; then
                printf '%s\n' "$1"
            else
                printf '1\n'
            fi
            ;;
        *)
            printf 'unknown\n'
            ;;
    esac
}

# NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
# ```c
# while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
#     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
#         continue;
#     }
# }
# ```
wasm_threading_status() {
    case "${LOSAT_WASM_EXECUTION_MODE}" in
        command-wasi-serial)
            printf 'serial-command-wasi\n'
            ;;
        browser-in-memory-serial)
            printf 'serial-browser-in-memory\n'
            ;;
        browser-worker-parallel)
            if [ "${LOSAT_WASM_THREADS_VERIFIED}" = "1" ]; then
                printf 'verified-browser-worker-parallel\n'
            else
                printf 'unverified-browser-worker-parallel-disabled\n'
            fi
            ;;
        future-wasi-threaded)
            if [ "${LOSAT_WASM_THREADS_VERIFIED}" = "1" ]; then
                printf 'verified-wasi-threaded\n'
            else
                printf 'unverified-wasi-threaded-disabled\n'
            fi
            ;;
        *)
            printf 'unknown\n'
            ;;
    esac
}

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:75
# ```c
# const string kArgNumThreads("num_threads");
# ```
wasm_benchmark_label() {
    local requested_threads="$1"

    case "${LOSAT_WASM_EXECUTION_MODE}" in
        command-wasi-serial)
            if [ "${requested_threads}" = "1" ]; then
                printf 'LOSAT command-WASI serial\n'
            else
                printf 'LOSAT command-WASI serial (requested n%s)\n' "${requested_threads}"
            fi
            ;;
        browser-in-memory-serial)
            if [ "${requested_threads}" = "1" ]; then
                printf 'LOSAT browser in-memory serial\n'
            else
                printf 'LOSAT browser in-memory serial (requested n%s)\n' "${requested_threads}"
            fi
            ;;
        browser-worker-parallel)
            if [ "${LOSAT_WASM_THREADS_VERIFIED}" = "1" ]; then
                printf 'LOSAT browser worker n%s\n' "${requested_threads}"
            else
                printf 'LOSAT browser worker unverified (requested n%s)\n' "${requested_threads}"
            fi
            ;;
        future-wasi-threaded)
            if [ "${LOSAT_WASM_THREADS_VERIFIED}" = "1" ]; then
                printf 'LOSAT WASI threaded n%s\n' "${requested_threads}"
            else
                printf 'LOSAT WASI threaded unverified (requested n%s)\n' "${requested_threads}"
            fi
            ;;
        *)
            printf 'LOSAT wasm unknown-mode\n'
            ;;
    esac
}

# ==========================================
# Part 1: LOSAT Execution
# ==========================================
echo "Starting LOSAT commands..."
<<COMMENTOUT
run_losatn_native_threaded_case() {
    local query="$1"
    local subject="$2"
    local stem="$3"
    local task="$4"
    local threaded_suffix="n${LOSAT_THREADS}"
    local task_args=()

    if [ "${task}" = "blastn" ]; then
        task_args=(-task blastn)
    fi

    echo "Starting ${stem} (LOSAT native ${threaded_suffix})..."
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    # ```c
    # const string kArgQuery("query");
    # const string kArgSubject("subject");
    # const string kArgOutput("out");
    # const string kArgNumThreads("num_threads");
    # const string kTask("task");
    # ```
    (time "${LOSAT_BIN}" blastn -query "./fasta/${query}" -subject "./fasta/${subject}" -out "./losat_out/${stem}.${threaded_suffix}.out" "${task_args[@]}" -num_threads "${LOSAT_THREADS}" )&>"./losat_out/${stem}.${threaded_suffix}.log"
}

# --- LOSATN Commands (Default / Megablast behavior) ---

# NZ_CP006932 self (Default/Megablast)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.losatn.megablast.out -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.losatn.megablast.log

# EDL933 vs Sakai
(time $LOSAT_BIN blastn -q ./fasta/EDL933.fna -s ./fasta/Sakai.fna -o ./losat_out/EDL933.Sakai.losatn.megablast.out -n 1 )&>./losat_out/EDL933.Sakai.losatn.megablast.log

# Sakai vs MG1655
(time $LOSAT_BIN blastn -q ./fasta/Sakai.fna -s ./fasta/MG1655.fna -o ./losat_out/Sakai.MG1655.losatn.megablast.out -n 1 )&>./losat_out/Sakai.MG1655.losatn.megablast.log


# --- LOSATN Commands (Task: blastn) ---
# NZ_CP006932 self (Task: blastn)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out --task blastn -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.log
# PesePMNV vs MjPMNV
(time $LOSAT_BIN blastn -q ./fasta/AP027152.fasta -s ./fasta/AP027202.fasta -o ./losat_out/PesePMNV.MjPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PesePMNV.MjPMNV.losatn.blastn.log

# MelaMJNV vs PemoMJNVA
(time $LOSAT_BIN blastn -q ./fasta/LC738874.fasta -s ./fasta/LC738870.fasta -o ./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.log

# SiNMV vs ChdeNMV
(time $LOSAT_BIN blastn -q ./fasta/LC738884.fasta -s ./fasta/AP027155.fasta -o ./losat_out/SiNMV.ChdeNMV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/SiNMV.ChdeNMV.losatn.blastn.log

# PmeNMV vs MjPMNV
(time $LOSAT_BIN blastn -q ./fasta/LC738869.fasta -s ./fasta/AP027202.fasta -o ./losat_out/PmeNMV.MjPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PmeNMV.MjPMNV.losatn.blastn.log

# PmeNMV vs PesePMNV
(time $LOSAT_BIN blastn -q ./fasta/LC738869.fasta -s ./fasta/AP027152.fasta -o ./losat_out/PmeNMV.PesePMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PmeNMV.PesePMNV.losatn.blastn.log

# PeseMJNV vs PemoMJNVB
(time $LOSAT_BIN blastn -q ./fasta/LC738873.fasta -s ./fasta/LC738871.fasta -o ./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.log

# PemoMJNVA vs PeseMJNV
(time $LOSAT_BIN blastn -q ./fasta/LC738870.fasta -s ./fasta/LC738873.fasta -o ./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.log

# MjeNMV vs MelaMJNV
(time $LOSAT_BIN blastn -q ./fasta/LC738868.fasta -s ./fasta/LC738874.fasta -o ./losat_out/MjeNMV.MelaMJNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MjeNMV.MelaMJNV.losatn.blastn.log

# MjPMNV vs MlPMNV
(time $LOSAT_BIN blastn -q ./fasta/AP027202.fasta -s ./fasta/LC738875.fasta -o ./losat_out/MjPMNV.MlPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MjPMNV.MlPMNV.losatn.blastn.log

# --- LOSATN multithread timing commands ---
run_losatn_native_threaded_case "NZ_CP006932.fasta" "NZ_CP006932.fasta" "NZ_CP006932.NZ_CP006932.losatn.megablast" "megablast"
run_losatn_native_threaded_case "EDL933.fna" "Sakai.fna" "EDL933.Sakai.losatn.megablast" "megablast"
run_losatn_native_threaded_case "Sakai.fna" "MG1655.fna" "Sakai.MG1655.losatn.megablast" "megablast"

run_losatn_native_threaded_case "NZ_CP006932.fasta" "NZ_CP006932.fasta" "NZ_CP006932.NZ_CP006932.losatn.blastn" "blastn"
run_losatn_native_threaded_case "AP027152.fasta" "AP027202.fasta" "PesePMNV.MjPMNV.losatn.blastn" "blastn"
run_losatn_native_threaded_case "LC738874.fasta" "LC738870.fasta" "MelaMJNV.PemoMJNVA.losatn.blastn" "blastn"
run_losatn_native_threaded_case "LC738884.fasta" "AP027155.fasta" "SiNMV.ChdeNMV.losatn.blastn" "blastn"
run_losatn_native_threaded_case "LC738869.fasta" "AP027202.fasta" "PmeNMV.MjPMNV.losatn.blastn" "blastn"
run_losatn_native_threaded_case "LC738869.fasta" "AP027152.fasta" "PmeNMV.PesePMNV.losatn.blastn" "blastn"
run_losatn_native_threaded_case "LC738873.fasta" "LC738871.fasta" "PeseMJNV.PemoMJNVB.losatn.blastn" "blastn"
run_losatn_native_threaded_case "LC738870.fasta" "LC738873.fasta" "PemoMJNVA.PeseMJNV.losatn.blastn" "blastn"
run_losatn_native_threaded_case "LC738868.fasta" "LC738874.fasta" "MjeNMV.MelaMJNV.losatn.blastn" "blastn"
run_losatn_native_threaded_case "AP027202.fasta" "LC738875.fasta" "MjPMNV.MlPMNV.losatn.blastn" "blastn"
COMMENTOUT
# --- LOSATP Commands ---
echo "Starting LOSATP commands..."
<<COMMENTOUT
run_losatp_case() {
    local query="$1"
    local subject="$2"
    local stem="$3"
    local threaded_suffix="n${LOSAT_THREADS}"

    echo "Starting ${stem} (LOSATP)..."
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    # ```c
    # const string kArgQuery("query");
    # const string kArgSubject("subject");
    # const string kArgOutput("out");
    # const string kArgNumThreads("num_threads");
    # ```
    (time "${LOSAT_BIN}" blastp -query "./fasta/${query}" -subject "./fasta/${subject}" -out "./losat_out/${stem}.losatp.out" -num_threads 1 -outfmt 6 )&>"./losat_out/${stem}.losatp.log"
    (time "${LOSAT_BIN}" blastp -query "./fasta/${query}" -subject "./fasta/${subject}" -out "./losat_out/${stem}.losatp.${threaded_suffix}.out" -num_threads "${LOSAT_THREADS}" -outfmt 6 )&>"./losat_out/${stem}.losatp.${threaded_suffix}.log"
}
COMMENTOUT
run_tlosatx_case() {
    local query="$1"
    local subject="$2"
    local stem="$3"
    local query_gencode="$4"
    local db_gencode="$5"
    local threaded_suffix="n${LOSAT_THREADS}"

    echo "Starting ${stem} (TLOSATX)..."
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
    # ```c
    # const string kArgQuery("query");
    # const string kArgOutput("out");
    # const string kArgSubject("subject");
    # const string kArgQueryGeneticCode("query_gencode");
    # const string kArgDbGeneticCode("db_gencode");
    # const string kArgNumThreads("num_threads");
    # ```
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2638-2662
    # ```c
    # arg_desc.AddDefaultKey(kArgOutputFormat, "format",
    #                        kOutputFormatDescription,
    #                        CArgDescriptions::eString,
    #                        NStr::IntToString(dft_outfmt));
    # ```
    (time "${LOSAT_BIN}" tblastx -query "./fasta/${query}" -subject "./fasta/${subject}" -out "./losat_out/${stem}.tlosatx.n1.out" --query-gencode "${query_gencode}" --db-gencode "${db_gencode}" -num_threads 1 -outfmt 6 )&>"./losat_out/${stem}.tlosatx.n1.log"
    (time "${LOSAT_BIN}" tblastx -query "./fasta/${query}" -subject "./fasta/${subject}" -out "./losat_out/${stem}.tlosatx.${threaded_suffix}.out" --query-gencode "${query_gencode}" --db-gencode "${db_gencode}" -num_threads "${LOSAT_THREADS}" -outfmt 6 )&>"./losat_out/${stem}.tlosatx.${threaded_suffix}.log"
}

LOSATP_CASES=(
    "WSSV.faa:PajaWSV.faa:WSSV.PajaWSV"
    "WSSV.faa:SicyWSV.faa:WSSV.SicyWSV"
    "WSSV.faa:CoBV.faa:WSSV.CoBV"
    "PajaWSV.faa:SicyWSV.faa:PajaWSV.SicyWSV"
    "PajaWSV.faa:CoBV.faa:PajaWSV.CoBV"
    "SicyWSV.faa:CoBV.faa:SicyWSV.CoBV"
    "AP027078.faa:AP027131.faa:AP027078.AP027131"
    "AP027078.faa:AP027132.faa:AP027078.AP027132"
    "AP027078.faa:AP027133.faa:AP027078.AP027133"
    "AP027078.faa:NZ_CP006932.faa:AP027078.NZ_CP006932"
    "AP027131.faa:AP027132.faa:AP027131.AP027132"
    "AP027131.faa:AP027133.faa:AP027131.AP027133"
    "AP027131.faa:NZ_CP006932.faa:AP027131.NZ_CP006932"
    "AP027132.faa:AP027133.faa:AP027132.AP027133"
    "AP027132.faa:NZ_CP006932.faa:AP027132.NZ_CP006932"
    "AP027133.faa:NZ_CP006932.faa:AP027133.NZ_CP006932"
)

TBLASTX_CASES=(
    "AP027280.fasta:AP027280.fasta:AP027280:AP027280.AP027280:1:1"
    "MjeNMV.fasta:MelaMJNV.fasta:MelaMJNV:MjeNMV.MelaMJNV:1:1"
    "MelaMJNV.fasta:PemoMJNVA.fasta:PemoMJNVA:MelaMJNV.PemoMJNVA:1:1"
    "PemoMJNVA.fasta:PeseMJNV.fasta:PeseMJNV:PemoMJNVA.PeseMJNV:1:1"
    "PeseMJNV.fasta:PemoMJNVB.fasta:PemoMJNVB:PeseMJNV.PemoMJNVB:1:1"
    "PemoMJNVB.fasta:LvMJNV.fasta:LvMJNV:PemoMJNVB.LvMJNV:1:1"
    "LvMJNV.fasta:TrcuMJNV.fasta:TrcuMJNV:LvMJNV.TrcuMJNV:1:1"
    "TrcuMJNV.fasta:MellatMJNV.fasta:MellatMJNV:TrcuMJNV.MellatMJNV:1:1"
    "MellatMJNV.fasta:MeenMJNV.fasta:MeenMJNV:MellatMJNV.MeenMJNV:1:1"
    "MeenMJNV.fasta:MejoMJNV.fasta:MejoMJNV:MeenMJNV.MejoMJNV:1:1"

)
#     "AvCLPV.fasta:PsCLPV.fasta:PsCLPV:AvCLPV.PsCLPV:1:1"
#     "NZ_CP006932.fasta:NZ_CP006932.fasta:NZ_CP006932:NZ_CP006932.NZ_CP006932:4:4"
#     "AP027132.fasta:NZ_CP006932.fasta:NZ_CP006932:AP027132.NZ_CP006932:4:4"
#     "AP027078.fasta:AP027131.fasta:AP027131:AP027078.AP027131:4:4"
#     "AP027131.fasta:AP027133.fasta:AP027133:AP027131.AP027133:4:4"
#     "AP027133.fasta:AP027132.fasta:AP027132:AP027133.AP027132:4:4"
for losatp_case in "${LOSATP_CASES[@]}"; do
    IFS=":" read -r query subject stem <<<"$losatp_case"
    run_losatp_case "$query" "$subject" "$stem"
done


echo "Finished LOSATP commands!"

echo "Starting TLOSATX commands..."
for tblastx_case in "${TBLASTX_CASES[@]}"; do
    IFS=":" read -r query subject _db stem query_gencode db_gencode <<<"$tblastx_case"
    run_tlosatx_case "$query" "$subject" "$stem" "$query_gencode" "$db_gencode"
done
echo "Finished TLOSATX commands!"

# --- LOSAT Wasm Commands ---
echo "Starting LOSAT Wasm commands..."
write_losat_wasm_runner() {
    cat >"${LOSAT_WASM_RUNNER}" <<'JS'
const { WASI } = require('wasi');
const fs = require('fs');

const wasmPath = process.argv[2];
const args = process.argv.slice(3);
const wasi = new WASI({
  version: 'preview1',
  args: [wasmPath, ...args],
  env: process.env,
  preopens: { '/': '/' },
});

(async () => {
  const bytes = fs.readFileSync(wasmPath);
  const module = await WebAssembly.compile(bytes);
  const instance = await WebAssembly.instantiate(module, {
    wasi_snapshot_preview1: wasi.wasiImport,
  });
  if (typeof instance.exports._start !== 'function') {
    console.error(`${wasmPath} does not export _start; use the LOSAT command Wasm artifact.`);
    process.exit(1);
  }
  wasi.start(instance);
})().catch((err) => {
  console.error(err);
  process.exit(1);
});
JS
}

run_losat_wasm() {
    local log="$1"
    local requested_threads
    local effective_threads
    shift

    requested_threads="$(wasm_requested_num_threads "$@")"
    effective_threads="$(wasm_effective_engine_threads "${requested_threads}")"
    {
        printf '[WASM-METADATA] target_triple=%s\n' "${LOSAT_WASM_TARGET}"
        printf '[WASM-METADATA] feature_set=%s\n' "${LOSAT_WASM_FEATURES}"
        printf '[WASM-METADATA] execution_mode=%s\n' "${LOSAT_WASM_EXECUTION_MODE}"
        printf '[WASM-METADATA] threading_status=%s\n' "$(wasm_threading_status)"
        printf '[WASM-METADATA] rayon_compiled=%s\n' "$(wasm_rayon_compiled)"
        printf '[WASM-METADATA] requested_num_threads=%s\n' "${requested_threads}"
        printf '[WASM-METADATA] effective_engine_threads=%s\n' "${effective_threads}"
        printf '[WASM-METADATA] benchmark_label=%s\n' "$(wasm_benchmark_label "${requested_threads}")"
        printf '[WASM-METADATA] browser_worker_parallel_verified=false\n'
        printf '[WASM-METADATA] wasi_runtime=node-wasi-preview1\n'
        printf '[WASM-METADATA] node_version=%s\n' "$(node --version 2>/dev/null || printf 'unavailable')"
        time env NODE_NO_WARNINGS=1 node "${LOSAT_WASM_RUNNER}" "${LOSAT_WASM_BIN}" "$@"
    } &>"${log}"
}

run_losatn_wasm_case() {
    local query="$1"
    local subject="$2"
    local stem="$3"
    local task="$4"
    local threaded_suffix="n${LOSAT_THREADS}"
    local out="${LOSAT_OUT_DIR}/${stem}.wasm.out"
    local log="${LOSAT_OUT_DIR}/${stem}.wasm.log"
    local threaded_out="${LOSAT_OUT_DIR}/${stem}.wasm.${threaded_suffix}.out"
    local threaded_log="${LOSAT_OUT_DIR}/${stem}.wasm.${threaded_suffix}.log"
    local task_args=()

    if [ "${task}" = "blastn" ]; then
        task_args=(-task blastn)
    fi

    echo "Starting ${stem} (LOSAT Wasm)..."
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    # ```c
    # const string kArgQuery("query");
    # const string kArgSubject("subject");
    # const string kArgOutput("out");
    # const string kArgNumThreads("num_threads");
    # const string kTask("task");
    # ```
    run_losat_wasm "${log}" blastn -query "${FASTA_DIR}/${query}" -subject "${FASTA_DIR}/${subject}" -out "${out}" "${task_args[@]}" -num_threads 1
    run_losat_wasm "${threaded_log}" blastn -query "${FASTA_DIR}/${query}" -subject "${FASTA_DIR}/${subject}" -out "${threaded_out}" "${task_args[@]}" -num_threads "${LOSAT_THREADS}"
}

run_losatp_wasm_case() {
    local query="$1"
    local subject="$2"
    local stem="$3"
    local threaded_suffix="n${LOSAT_THREADS}"

    echo "Starting ${stem} (LOSATP Wasm)..."
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    # ```c
    # const string kArgQuery("query");
    # const string kArgSubject("subject");
    # const string kArgOutput("out");
    # const string kArgNumThreads("num_threads");
    # ```
    run_losat_wasm \
        "${LOSAT_OUT_DIR}/${stem}.losatp.wasm.log" \
        blastp -query "${FASTA_DIR}/${query}" -subject "${FASTA_DIR}/${subject}" \
        -out "${LOSAT_OUT_DIR}/${stem}.losatp.wasm.out" -num_threads 1 -outfmt 6
    run_losat_wasm \
        "${LOSAT_OUT_DIR}/${stem}.losatp.wasm.${threaded_suffix}.log" \
        blastp -query "${FASTA_DIR}/${query}" -subject "${FASTA_DIR}/${subject}" \
        -out "${LOSAT_OUT_DIR}/${stem}.losatp.wasm.${threaded_suffix}.out" -num_threads "${LOSAT_THREADS}" -outfmt 6
}

run_tlosatx_wasm_case() {
    local query="$1"
    local subject="$2"
    local stem="$3"
    local query_gencode="$4"
    local db_gencode="$5"
    local threaded_suffix="n${LOSAT_THREADS}"

    echo "Starting ${stem} (TLOSATX Wasm)..."
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
    # ```c
    # const string kArgQuery("query");
    # const string kArgOutput("out");
    # const string kArgSubject("subject");
    # const string kArgQueryGeneticCode("query_gencode");
    # const string kArgDbGeneticCode("db_gencode");
    # const string kArgNumThreads("num_threads");
    # ```
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2638-2662
    # ```c
    # arg_desc.AddDefaultKey(kArgOutputFormat, "format",
    #                        kOutputFormatDescription,
    #                        CArgDescriptions::eString,
    #                        NStr::IntToString(dft_outfmt));
    # ```
    run_losat_wasm \
        "${LOSAT_OUT_DIR}/${stem}.tlosatx.wasm.log" \
        tblastx -query "${FASTA_DIR}/${query}" -subject "${FASTA_DIR}/${subject}" \
        -out "${LOSAT_OUT_DIR}/${stem}.tlosatx.wasm.out" --query-gencode "${query_gencode}" \
        --db-gencode "${db_gencode}" -num_threads 1 -outfmt 6
    run_losat_wasm \
        "${LOSAT_OUT_DIR}/${stem}.tlosatx.wasm.${threaded_suffix}.log" \
        tblastx -query "${FASTA_DIR}/${query}" -subject "${FASTA_DIR}/${subject}" \
        -out "${LOSAT_OUT_DIR}/${stem}.tlosatx.wasm.${threaded_suffix}.out" --query-gencode "${query_gencode}" \
        --db-gencode "${db_gencode}" -num_threads "${LOSAT_THREADS}" -outfmt 6
}

if [ "${RUN_LOSAT_WASM}" = "0" ]; then
    echo "Skipping LOSAT Wasm commands because RUN_LOSAT_WASM=0."
elif [ "${LOSAT_WASM_EXECUTION_MODE}" != "command-wasi-serial" ]; then
    echo "Skipping LOSAT command-Wasm commands because run_comparison.sh only executes command-wasi-serial mode."
    echo "Requested LOSAT_WASM_EXECUTION_MODE=${LOSAT_WASM_EXECUTION_MODE}; use a browser/web_api harness for browser modes or a dedicated WASI-threaded runner once verified."
else
    if [ "${BUILD_LOSAT_WASM}" = "1" ]; then
        (cd .. && cargo build --release --target "${LOSAT_WASM_TARGET}" ${LOSAT_WASM_FEATURES})
    fi

    if ! command -v node >/dev/null 2>&1; then
        echo "Skipping LOSAT Wasm commands because node is not available."
    elif [ ! -f "${LOSAT_WASM_BIN}" ]; then
        echo "Skipping LOSAT Wasm commands because ${LOSAT_WASM_BIN} does not exist."
    elif ! LOSAT_WASM_BIN_RESOLVED="$(resolve_losat_wasm_bin "${LOSAT_WASM_BIN}")"; then
        echo "Skipping LOSAT Wasm commands because no _start-exporting LOSAT command Wasm was found near ${LOSAT_WASM_BIN}."
    else
        if [ "${LOSAT_WASM_BIN_RESOLVED}" != "${LOSAT_WASM_BIN}" ]; then
            echo "Using LOSAT Wasm command artifact: ${LOSAT_WASM_BIN_RESOLVED}"
            LOSAT_WASM_BIN="${LOSAT_WASM_BIN_RESOLVED}"
        fi
        write_losat_wasm_runner
        <<COMMENTOUT
        run_losatn_wasm_case "NZ_CP006932.fasta" "NZ_CP006932.fasta" "NZ_CP006932.NZ_CP006932.losatn.megablast" "megablast"
        run_losatn_wasm_case "EDL933.fna" "Sakai.fna" "EDL933.Sakai.losatn.megablast" "megablast"
        run_losatn_wasm_case "Sakai.fna" "MG1655.fna" "Sakai.MG1655.losatn.megablast" "megablast"

        run_losatn_wasm_case "NZ_CP006932.fasta" "NZ_CP006932.fasta" "NZ_CP006932.NZ_CP006932.losatn.blastn" "blastn"
        run_losatn_wasm_case "AP027152.fasta" "AP027202.fasta" "PesePMNV.MjPMNV.losatn.blastn" "blastn"
        run_losatn_wasm_case "LC738874.fasta" "LC738870.fasta" "MelaMJNV.PemoMJNVA.losatn.blastn" "blastn"
        run_losatn_wasm_case "LC738884.fasta" "AP027155.fasta" "SiNMV.ChdeNMV.losatn.blastn" "blastn"
        run_losatn_wasm_case "LC738869.fasta" "AP027202.fasta" "PmeNMV.MjPMNV.losatn.blastn" "blastn"
        run_losatn_wasm_case "LC738869.fasta" "AP027152.fasta" "PmeNMV.PesePMNV.losatn.blastn" "blastn"
        run_losatn_wasm_case "LC738873.fasta" "LC738871.fasta" "PeseMJNV.PemoMJNVB.losatn.blastn" "blastn"
        run_losatn_wasm_case "LC738870.fasta" "LC738873.fasta" "PemoMJNVA.PeseMJNV.losatn.blastn" "blastn"
        run_losatn_wasm_case "LC738868.fasta" "LC738874.fasta" "MjeNMV.MelaMJNV.losatn.blastn" "blastn"
        run_losatn_wasm_case "AP027202.fasta" "LC738875.fasta" "MjPMNV.MlPMNV.losatn.blastn" "blastn"

        for losatp_case in "${LOSATP_CASES[@]}"; do
            IFS=":" read -r query subject stem <<<"$losatp_case"
            run_losatp_wasm_case "$query" "$subject" "$stem"
        done
        COMMENTOUT
        for tblastx_case in "${TBLASTX_CASES[@]}"; do
            IFS=":" read -r query subject _db stem query_gencode db_gencode <<<"$tblastx_case"
            run_tlosatx_wasm_case "$query" "$subject" "$stem" "$query_gencode" "$db_gencode"
        done
    fi
fi
echo "Finished LOSAT Wasm commands!"

<<COMMENTOUT
# --- TLOSATX Commands (Genetic Code: 4) ---
echo "Starting TLOSATX commands..."

# --- TLOSATX Commands (Genetic Code: 1) ---
echo "Starting AP027280 self..."
(time $LOSAT_BIN tblastx -q ./fasta/AP027280.fasta -s ./fasta/AP027280.fasta -o ./losat_out/AP027280.AP027280.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1)&>./losat_out/AP027280.AP027280.tlosatx.n1.log 
(time $LOSAT_BIN tblastx -q ./fasta/AP027280.fasta -s ./fasta/AP027280.fasta -o ./losat_out/AP027280.AP027280.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 )&>./losat_out/AP027280.AP027280.tlosatx.n8.log 

# MjeNMV vs MelaMJNV
echo "Starting MjeNMV vs MelaMJNV (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/MjeNMV.fasta -s ./fasta/MelaMJNV.fasta -o ./losat_out/MjeNMV.MelaMJNV.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/MjeNMV.MelaMJNV.tlosatx.n8.log 

# MelaMJNV vs PemoMJNVA
echo "Starting MelaMJNV vs PemoMJNVA (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/MelaMJNV.fasta -s ./fasta/PemoMJNVA.fasta -o ./losat_out/MelaMJNV.PemoMJNVA.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/MelaMJNV.PemoMJNVA.tlosatx.n8.log 

# PemoMJNVA vs PeseMJNV
echo "Starting PemoMJNVA vs PeseMJNV (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/PemoMJNVA.fasta -s ./fasta/PeseMJNV.fasta -o ./losat_out/PemoMJNVA.PeseMJNV.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/PemoMJNVA.PeseMJNV.tlosatx.n8.log 

# PeseMJNV vs PemoMJNVB
echo "Starting PeseMJNV vs PemoMJNVB (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/PeseMJNV.fasta -s ./fasta/PemoMJNVB.fasta -o ./losat_out/PeseMJNV.PemoMJNVB.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/PeseMJNV.PemoMJNVB.tlosatx.n8.log 

# PemoMJNVB vs LvMJNV
echo "Starting PemoMJNVB vs LvMJNV (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/PemoMJNVB.fasta -s ./fasta/LvMJNV.fasta -o ./losat_out/PemoMJNVB.LvMJNV.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/PemoMJNVB.LvMJNV.tlosatx.n8.log 

# LvMJNV vs TrcuMJNV
echo "Starting LvMJNV vs TrcuMJNV (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/LvMJNV.fasta -s ./fasta/TrcuMJNV.fasta -o ./losat_out/LvMJNV.TrcuMJNV.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/LvMJNV.TrcuMJNV.tlosatx.n8.log 

# TrcuMJNV vs MellatMJNV
echo "Starting TrcuMJNV vs MellatMJNV (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/TrcuMJNV.fasta -s ./fasta/MellatMJNV.fasta -o ./losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.log 

# MellatMJNV vs MeenMJNV
echo "Starting MellatMJNV vs MeenMJNV (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/MellatMJNV.fasta -s ./fasta/MeenMJNV.fasta -o ./losat_out/MellatMJNV.MeenMJNV.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/MellatMJNV.MeenMJNV.tlosatx.n8.log 

# MeenMJNV vs MejoMJNV
echo "Starting MeenMJNV vs MejoMJNV (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/MeenMJNV.fasta -s ./fasta/MejoMJNV.fasta -o ./losat_out/MeenMJNV.MejoMJNV.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 ) &> ./losat_out/MeenMJNV.MejoMJNV.tlosatx.n8.log 

# AvCLPV vs PsCLPV
echo "Starting AvCLPV vs PsCLPV (TLOSATX)..."
(time $LOSAT_BIN tblastx -q ./fasta/AvCLPV.fasta -s ./fasta/PsCLPV.fasta -o ./losat_out/AvCLPV.PsCLPV.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 )&>./losat_out/AvCLPV.PsCLPV.tlosatx.n8.log 

echo "Starting NZ_CP006932 self (TLOSATX)..."
# NZ_CP006932 self
(time $LOSAT_BIN tblastx -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n8.out --query-gencode 4 --db-gencode 4 -n 8 )&>./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n8.log 
echo "Starting AP027132 vs NZ_CP006932 (TLOSATX)..."
# AP027132 vs NZ_CP006932
(time $LOSAT_BIN tblastx -q ./fasta/AP027132.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/AP027132.NZ_CP006932.tlosatx.n8.out --query-gencode 4 --db-gencode 4 -n 8 )&>./losat_out/AP027132.NZ_CP006932.tlosatx.n8.log 
echo "Starting AP027078 vs AP027131 (TLOSATX)..."
# AP027078 vs AP027131
(time $LOSAT_BIN tblastx -q ./fasta/AP027078.fasta -s ./fasta/AP027131.fasta -o ./losat_out/AP027078.AP027131.tlosatx.n8.out --query-gencode 4 --db-gencode 4 -n 8 )&>./losat_out/AP027078.AP027131.tlosatx.n8.log 
echo "Starting AP027131 vs AP027133 (TLOSATX)..."
# AP027131 vs AP027133
(time $LOSAT_BIN tblastx -q ./fasta/AP027131.fasta -s ./fasta/AP027133.fasta -o ./losat_out/AP027131.AP027133.tlosatx.n8.out --query-gencode 4 --db-gencode 4 -n 8 )&>./losat_out/AP027131.AP027133.tlosatx.n8.log 
echo "Starting AP027133 vs AP027132 (TLOSATX)..."
# AP027133 vs AP027132
(time $LOSAT_BIN tblastx -q ./fasta/AP027133.fasta -s ./fasta/AP027132.fasta -o ./losat_out/AP027133.AP027132.tlosatx.n8.out --query-gencode 4 --db-gencode 4 -n 8 )&>./losat_out/AP027133.AP027132.tlosatx.n8.log 

echo "Finished TLOSATX commands!"
COMMENTOUT


echo "Finished LOSAT commands!"

<<COMMENTOUT
# ==========================================
# Part 2: BLAST+ Execution
# ==========================================
echo "Starting BLAST+ commands..."


# --- TBLASTX Commands (gencode 1) ---

# makeblastdb -in AP027280.fasta -dbtype nucl -title AP027280 -parse_seqids -hash_index -out AP027280
# makeblastdb -in MjeNMV.fasta -dbtype nucl -title MjeNMV -parse_seqids -hash_index -out MjeNMV
# makeblastdb -in MelaMJNV.fasta -dbtype nucl -title MelaMJNV -parse_seqids -hash_index -out MelaMJNV
# makeblastdb -in PemoMJNVA.fasta -dbtype nucl -title PemoMJNVA -parse_seqids -hash_index -out PemoMJNVA
# makeblastdb -in PeseMJNV.fasta -dbtype nucl -title PeseMJNV -parse_seqids -hash_index -out PeseMJNV
# makeblastdb -in PemoMJNVB.fasta -dbtype nucl -title PemoMJNVB -parse_seqids -hash_index -out PemoMJNVB
# makeblastdb -in LvMJNV.fasta -dbtype nucl -title LvMJNV -parse_seqids -hash_index -out LvMJNV
# makeblastdb -in TrcuMJNV.fasta -dbtype nucl -title TrcuMJNV -parse_seqids -hash_index -out TrcuMJNV
# makeblastdb -in MellatMJNV.fasta -dbtype nucl -title MellatMJNV -parse_seqids -hash_index -out MellatMJNV
# makeblastdb -in MeenMJNV.fasta -dbtype nucl -title MeenMJNV -parse_seqids -hash_index -out MeenMJNV
# makeblastdb -in MejoMJNV.fasta -dbtype nucl -title MejoMJNV -parse_seqids -hash_index -out MejoMJNV
# makeblastdb -in AvCLPV.fasta -dbtype nucl -title AvCLPV -parse_seqids -hash_index -out AvCLPV
# makeblastdb -in PsCLPV.fasta -dbtype nucl -title PsCLPV -parse_seqids -hash_index -out PsCLPV

# AP027280 self
echo "Starting AP027280_Self (BLAST)..."
(time tblastx -query ./fasta/AP027280.fasta -db ./fasta/AP027280 -out ./blast_out/AP027280.AP027280.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6 )&>./blast_out/AP027280.AP027280.tblastx.n8.log
(time tblastx -query ./fasta/AP027280.fasta -db ./fasta/AP027280 -out ./blast_out/AP027280.AP027280.tblastx.windowsize0.n1.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6 -window_size 0)&>./blast_out/AP027280.AP027280.tblastx.windowsize0.n1.log

# MjeNMV vs MelaMJNV
echo "Starting MjeNMV vs MelaMJNV (BLAST)..."
(time tblastx -query ./fasta/MjeNMV.fasta -db ./fasta/MelaMJNV -out ./blast_out/MjeNMV.MelaMJNV.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/MjeNMV.MelaMJNV.tblastx.n8.log

# MelaMJNV vs PemoMJNVA
echo "Starting MelaMJNV vs PemoMJNVA (BLAST)..."
(time tblastx -query ./fasta/MelaMJNV.fasta -db ./fasta/PemoMJNVA -out ./blast_out/MelaMJNV.PemoMJNVA.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/MelaMJNV.PemoMJNVA.tblastx.n8.log

# PemoMJNVA vs PeseMJNV
echo "Starting PemoMJNVA vs PeseMJNV (BLAST)..."
(time tblastx -query ./fasta/PemoMJNVA.fasta -db ./fasta/PeseMJNV -out ./blast_out/PemoMJNVA.PeseMJNV.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/PemoMJNVA.PeseMJNV.tblastx.n8.log

# PeseMJNV vs PemoMJNVB
echo "Starting PeseMJNV vs PemoMJNVB (BLAST)..."
(time tblastx -query ./fasta/PeseMJNV.fasta -db ./fasta/PemoMJNVB -out ./blast_out/PeseMJNV.PemoMJNVB.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/PeseMJNV.PemoMJNVB.tblastx.n8.log

# PemoMJNVB vs LvMJNV
echo "Starting PemoMJNVB vs LvMJNV (BLAST)..."
(time tblastx -query ./fasta/PemoMJNVB.fasta -db ./fasta/LvMJNV -out ./blast_out/PemoMJNVB.LvMJNV.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/PemoMJNVB.LvMJNV.tblastx.n8.log

# LvMJNV vs TrcuMJNV
echo "Starting LvMJNV vs TrcuMJNV (BLAST)..."
(time tblastx -query ./fasta/LvMJNV.fasta -db ./fasta/TrcuMJNV -out ./blast_out/LvMJNV.TrcuMJNV.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/LvMJNV.TrcuMJNV.tblastx.n8.log

# TrcuMJNV vs MellatMJNV
echo "Starting TrcuMJNV vs MellatMJNV (BLAST)..."
(time tblastx -query ./fasta/TrcuMJNV.fasta -db ./fasta/MellatMJNV -out ./blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.log

# MellatMJNV vs MeenMJNV
echo "Starting MellatMJNV vs MeenMJNV (BLAST)..."
(time tblastx -query ./fasta/MellatMJNV.fasta -db ./fasta/MeenMJNV -out ./blast_out/MellatMJNV.MeenMJNV.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/MellatMJNV.MeenMJNV.tblastx.n8.log

# MeenMJNV vs MejoMJNV
echo "Starting MeenMJNV vs MejoMJNV (BLAST)..."
(time tblastx -query ./fasta/MeenMJNV.fasta -db ./fasta/MejoMJNV -out ./blast_out/MeenMJNV.MejoMJNV.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/MeenMJNV.MejoMJNV.tblastx.n8.log

# AvCLPV vs PsCLPV
echo "Starting AvCLPV vs PsCLPV (BLAST)..."
(time tblastx -query ./fasta/AvCLPV.fasta -db ./fasta/PsCLPV -out ./blast_out/AvCLPV.PsCLPV.tblastx.n8.out -query_gencode 1 -db_gencode 1 -num_threads 8 -outfmt 6) &> ./blast_out/AvCLPV.PsCLPV.tblastx.n8.log

COMMENTOUT


<<COMMENTOUT
# --- TBLASTX Commands (gencode 4) ---

# makeblastdb -in NZ_CP006932.fasta -dbtype nucl -title NZ_CP006932 -parse_seqids -hash_index -out NZ_CP006932
# makeblastdb -in AP027078.fasta -dbtype nucl -title AP027078 -parse_seqids -hash_index -out AP027078
# makeblastdb -in AP027131.fasta -dbtype nucl -title AP027131 -parse_seqids -hash_index -out AP027131
# makeblastdb -in AP027132.fasta -dbtype nucl -title AP027132 -parse_seqids -hash_index -out AP027132
# makeblastdb -in AP027133.fasta -dbtype nucl -title AP027133 -parse_seqids -hash_index -out AP027133


# NZ_CP006932 self
(time tblastx -query ./fasta/NZ_CP006932.fasta -db ./fasta/NZ_CP006932 -out ./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n8.out -query_gencode 4 -db_gencode 4 -num_threads 8 -outfmt 6 )&>./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n8.log &

# AP027132 vs NZ_CP006932
(time tblastx -query ./fasta/AP027132.fasta -db ./fasta/NZ_CP006932 -out ./blast_out/AP027132.NZ_CP006932.tblastx.n8.out -query_gencode 4 -db_gencode 4 -num_threads 8 -outfmt 6 )&>./blast_out/AP027132.NZ_CP006932.tblastx.n8.log &

# AP027078 vs AP027131
(time tblastx -query ./fasta/AP027078.fasta -db ./fasta/AP027131 -out ./blast_out/AP027078.AP027131.tblastx.n8.out -query_gencode 4 -db_gencode 4 -num_threads 8 -outfmt 6 )&>./blast_out/AP027078.AP027131.tblastx.n8.log &

# AP027131 vs AP027133
(time tblastx -query ./fasta/AP027131.fasta -db ./fasta/AP027133 -out ./blast_out/AP027131.AP027133.tblastx.n8.out -query_gencode 4 -db_gencode 4 -num_threads 8 -outfmt 6 )&>./blast_out/AP027131.AP027133.tblastx.n8.log &

# AP027133 vs AP027132
(time tblastx -query ./fasta/AP027133.fasta -db ./fasta/AP027132 -out ./blast_out/AP027133.AP027132.tblastx.n8.out -query_gencode 4 -db_gencode 4 -num_threads 8 -outfmt 6 )&>./blast_out/AP027133.AP027132.tblastx.n8.log &
COMMENTOUT
<<COMMENTOUT
run_tblastx_blast_case() {
    local query="$1"
    local db="$2"
    local stem="$3"
    local query_gencode="$4"
    local db_gencode="$5"
    local threaded_suffix="n${LOSAT_THREADS}"

    echo "Starting ${stem} (TBLASTX BLAST+)..."
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
    # ```c
    # const string kArgQuery("query");
    # const string kArgOutput("out");
    # const string kArgDb("db");
    # const string kArgQueryGeneticCode("query_gencode");
    # const string kArgDbGeneticCode("db_gencode");
    # const string kArgNumThreads("num_threads");
    # ```
    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2638-2662
    # ```c
    # arg_desc.AddDefaultKey(kArgOutputFormat, "format",
    #                        kOutputFormatDescription,
    #                        CArgDescriptions::eString,
    #                        NStr::IntToString(dft_outfmt));
    # ```
    (time "${TBLASTX_BIN}" -query "${FASTA_DIR}/${query}" -db "${FASTA_DIR}/${db}" -out "./blast_out/${stem}.tblastx.${threaded_suffix}.out" -query_gencode "${query_gencode}" -db_gencode "${db_gencode}" -num_threads "${LOSAT_THREADS}" -outfmt 6 )&>"./blast_out/${stem}.tblastx.${threaded_suffix}.log"
}

echo "Starting TBLASTX BLAST+ commands..."
for tblastx_case in "${TBLASTX_CASES[@]}"; do
    IFS=":" read -r query _subject db stem query_gencode db_gencode <<<"$tblastx_case"
    run_tblastx_blast_case "$query" "$db" "$stem" "$query_gencode" "$db_gencode"
done
echo "Finished TBLASTX BLAST+ commands!"
COMMENTOUT
# --- BLASTN Commands (Default / Megablast) ---

# NZ_CP006932 self (Default/Megablast)
(time blastn -query ./fasta/NZ_CP006932.fasta -subject ./fasta/NZ_CP006932.fasta -out ./blast_out/NZ_CP006932.NZ_CP006932.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/NZ_CP006932.NZ_CP006932.blastn.log

# EDL933 vs Sakai
(time blastn -query ./fasta/EDL933.fna -subject ./fasta/Sakai.fna -outfmt 7 -out ./blast_out/EDL933.Sakai.blastn.megablast.out )&>./blast_out/EDL933.Sakai.blastn.megablast.log

# Sakai vs MG1655
(time blastn -query ./fasta/Sakai.fna -subject ./fasta/MG1655.fna -outfmt 7 -out ./blast_out/Sakai.MG1655.blastn.megablast.out )&>./blast_out/Sakai.MG1655.blastn.megablast.log


# --- BLASTN Commands (Task: blastn) ---
<<COMMENTOUT
# NZ_CP006932 self (Task: blastn)
(time blastn -task blastn -query ./fasta/NZ_CP006932.fasta -subject ./fasta/NZ_CP006932.fasta -out ./blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/NZ_CP006932.NZ_CP006932.task_blastn.log

# PesePMNV vs MjPMNV
(time blastn -task blastn -query ./fasta/AP027152.fasta -subject ./fasta/AP027202.fasta -out ./blast_out/PesePMNV.MjPMNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PesePMNV.MjPMNV.blastn.log

# MelaMJNV vs PemoMJNVA
(time blastn -task blastn -query ./fasta/LC738874.fasta -subject ./fasta/LC738870.fasta -out ./blast_out/MelaMJNV.PemoMJNVA.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/MelaMJNV.PemoMJNVA.blastn.log

# SiNMV vs ChdeNMV
(time blastn -task blastn -query ./fasta/LC738884.fasta -subject ./fasta/AP027155.fasta -out ./blast_out/SiNMV.ChdeNMV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/SiNMV.ChdeNMV.blastn.log

# PmeNMV vs MjPMNV
(time blastn -task blastn -query ./fasta/LC738869.fasta -subject ./fasta/AP027202.fasta -out ./blast_out/PmeNMV.MjPMNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PmeNMV.MjPMNV.blastn.log

# PmeNMV vs PesePMNV
(time blastn -task blastn -query ./fasta/LC738869.fasta -subject ./fasta/AP027152.fasta -out ./blast_out/PmeNMV.PesePMNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PmeNMV.PesePMNV.blastn.log

# PeseMJNV vs PemoMJNVB
(time blastn -task blastn -query ./fasta/LC738873.fasta -subject ./fasta/LC738871.fasta -out ./blast_out/PeseMJNV.PemoMJNVB.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PeseMJNV.PemoMJNVB.blastn.log

# PemoMJNVA vs PeseMJNV
(time blastn -task blastn -query ./fasta/LC738870.fasta -subject ./fasta/LC738873.fasta -out ./blast_out/PemoMJNVA.PeseMJNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PemoMJNVA.PeseMJNV.blastn.log

# MjeNMV vs MelaMJNV
(time blastn -task blastn -query ./fasta/LC738868.fasta -subject ./fasta/LC738874.fasta -out ./blast_out/MjeNMV.MelaMJNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/MjeNMV.MelaMJNV.blastn.log

# MjPMNV vs MlPMNV
(time blastn -task blastn -query ./fasta/AP027202.fasta -subject ./fasta/LC738875.fasta -out ./blast_out/MjPMNV.MlPMNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/MjPMNV.MlPMNV.blastn.log


wait
COMMENTOUT
<<COMMENTOUT
# --- BLASTP Commands ---
echo "Starting WSSV vs PajaWSV (BLASTP)..."
(time  blastp -query ./fasta/WSSV.faa -subject ./fasta/PajaWSV.faa -out ./blast_out/WSSV.PajaWSV.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/WSSV.PajaWSV.BLASTP.log
echo "Starting WSSV vs SicyWSV (BLASTP)..."
(time  blastp -query ./fasta/WSSV.faa -subject ./fasta/SicyWSV.faa -out ./blast_out/WSSV.SicyWSV.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/WSSV.SicyWSV.BLASTP.log
echo "Starting WSSV vs CoBV (BLASTP)..."
(time  blastp -query ./fasta/WSSV.faa -subject ./fasta/CoBV.faa -out ./blast_out/WSSV.CoBV.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/WSSV.CoBV.BLASTP.log
echo "Starting PajaWSV vs SicyWSV (BLASTP)..."
(time  blastp -query ./fasta/PajaWSV.faa -subject ./fasta/SicyWSV.faa -out ./blast_out/PajaWSV.SicyWSV.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/PajaWSV.SicyWSV.BLASTP.log
echo "Starting PajaWSV vs CoBV (BLASTP)..."
(time  blastp -query ./fasta/PajaWSV.faa -subject ./fasta/CoBV.faa -out ./blast_out/PajaWSV.CoBV.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/PajaWSV.CoBV.BLASTP.log
echo "Starting SicyWSV vs CoBV (BLASTP)..."
(time  blastp -query ./fasta/SicyWSV.faa -subject ./fasta/CoBV.faa -out ./blast_out/SicyWSV.CoBV.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/SicyWSV.CoBV.BLASTP.log
echo "Starting AP027078 vs AP027131 (BLASTP)..."
(time  blastp -query ./fasta/AP027078.faa -subject ./fasta/AP027131.faa -out ./blast_out/AP027078.AP027131.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027078.AP027131.BLASTP.log
echo "Starting AP027078 vs AP027132 (BLASTP)..."
(time  blastp -query ./fasta/AP027078.faa -subject ./fasta/AP027132.faa -out ./blast_out/AP027078.AP027132.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027078.AP027132.BLASTP.log
echo "Starting AP027078 vs AP027133 (BLASTP)..."
(time  blastp -query ./fasta/AP027078.faa -subject ./fasta/AP027133.faa -out ./blast_out/AP027078.AP027133.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027078.AP027133.BLASTP.log
echo "Starting AP027078 vs NZ_CP006932 (BLASTP)..."
(time  blastp -query ./fasta/AP027078.faa -subject ./fasta/NZ_CP006932.faa -out ./blast_out/AP027078.NZ_CP006932.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027078.NZ_CP006932.BLASTP.log
echo "Starting AP027131 vs AP027132 (BLASTP)..."
(time  blastp -query ./fasta/AP027131.faa -subject ./fasta/AP027132.faa -out ./blast_out/AP027131.AP027132.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027131.AP027132.BLASTP.log
echo "Starting AP027131 vs AP027133 (BLASTP)..."
(time  blastp -query ./fasta/AP027131.faa -subject ./fasta/AP027133.faa -out ./blast_out/AP027131.AP027133.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027131.AP027133.BLASTP.log
echo "Starting AP027131 vs NZ_CP006932 (BLASTP)..."
(time  blastp -query ./fasta/AP027131.faa -subject ./fasta/NZ_CP006932.faa -out ./blast_out/AP027131.NZ_CP006932.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027131.NZ_CP006932.BLASTP.log
echo "Starting AP027132 vs AP027133 (BLASTP)..."
(time  blastp -query ./fasta/AP027132.faa -subject ./fasta/AP027133.faa -out ./blast_out/AP027132.AP027133.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027132.AP027133.BLASTP.log
echo "Starting AP027132 vs NZ_CP006932 (BLASTP)..."
(time  blastp -query ./fasta/AP027132.faa -subject ./fasta/NZ_CP006932.faa -out ./blast_out/AP027132.NZ_CP006932.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027132.NZ_CP006932.BLASTP.log
echo "Starting AP027133 vs NZ_CP006932 (BLASTP)..."
(time  blastp -query ./fasta/AP027133.faa -subject ./fasta/NZ_CP006932.faa -out ./blast_out/AP027133.NZ_CP006932.BLASTP.out -num_threads 1 -outfmt 6 )&>./blast_out/AP027133.NZ_CP006932.BLASTP.log
COMMENTOUT
echo "All done."
