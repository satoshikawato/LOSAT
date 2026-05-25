#!/bin/bash

# --- Setup ---

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

mkdir -p losat_out

LOSAT_THREADS="${LOSAT_THREADS:-${LOSATP_THREADS:-${LOSAT_BLASTP_THREADS:-8}}}"
LOSATP_THREADS="${LOSATP_THREADS:-${LOSAT_THREADS}}"
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
# ```c
# int num_threads = args[kArgNumThreads].AsInteger();
# if (num_threads > kMaxValue) {
#     m_NumThreads = kMaxValue;
# } else {
#     m_NumThreads = num_threads;
# }
# ```
LOSAT_WASM_BIN="${LOSAT_WASM_BIN:-../target/wasm32-wasip1/release/LOSAT.wasm}"
LOSAT_WASM_THREADED_BIN="${LOSAT_WASM_THREADED_BIN:-../target/wasm32-wasip1-threads/release/LOSAT.wasm}"
LOSAT_WASM_RUNNER="./losat_out/run_losat_wasi.js"
LOSAT_WASM_THREADED_RUNNER="${SCRIPT_DIR}/run_losat_wasi_threads.js"
RUN_LOSAT_WASM="${RUN_LOSAT_WASM:-1}"
RUN_LOSAT_WASM_THREADED="${RUN_LOSAT_WASM_THREADED:-1}"
BUILD_LOSAT_WASM="${BUILD_LOSAT_WASM:-0}"
BUILD_LOSAT_WASM_THREADED="${BUILD_LOSAT_WASM_THREADED:-${BUILD_LOSAT_WASM}}"
FASTA_DIR="${SCRIPT_DIR}/fasta"
LOSAT_OUT_DIR="${SCRIPT_DIR}/losat_out"

case "${LOSAT_WASM_BIN}" in
    /*) ;;
    *) LOSAT_WASM_BIN="${SCRIPT_DIR}/${LOSAT_WASM_BIN}" ;;
esac
case "${LOSAT_WASM_THREADED_BIN}" in
    /*) ;;
    *) LOSAT_WASM_THREADED_BIN="${SCRIPT_DIR}/${LOSAT_WASM_THREADED_BIN}" ;;
esac

wasm_has_start_export() {
    local wasm_path="$1"

    [ -f "${wasm_path}" ] || return 1
    node -e 'const fs=require("fs"); const p=process.argv[1]; const m=new WebAssembly.Module(fs.readFileSync(p)); process.exit(WebAssembly.Module.exports(m).some((e)=>e.name==="_start") ? 0 : 1);' "${wasm_path}" >/dev/null 2>&1
}

wasm_imports_wasi_thread_spawn() {
    local wasm_path="$1"

    [ -f "${wasm_path}" ] || return 1
    node -e 'const fs=require("fs"); const p=process.argv[1]; const m=new WebAssembly.Module(fs.readFileSync(p)); process.exit(WebAssembly.Module.imports(m).some((e)=>e.module==="wasi" && e.name==="thread-spawn") ? 0 : 1);' "${wasm_path}" >/dev/null 2>&1
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
    "AvCLPV.fasta:PsCLPV.fasta:PsCLPV:AvCLPV.PsCLPV:1:1"
    "NZ_CP006932.fasta:NZ_CP006932.fasta:NZ_CP006932:NZ_CP006932.NZ_CP006932:4:4"
    "AP027132.fasta:NZ_CP006932.fasta:NZ_CP006932:AP027132.NZ_CP006932:4:4"
    "AP027078.fasta:AP027131.fasta:AP027131:AP027078.AP027131:4:4"
    "AP027131.fasta:AP027133.fasta:AP027133:AP027131.AP027133:4:4"
    "AP027133.fasta:AP027132.fasta:AP027132:AP027133.AP027132:4:4"
)

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
    shift

    if [ "${LOSAT_WASM_AVAILABLE}" != "1" ]; then
        echo "Skipping ${log} because serial LOSAT Wasm is unavailable."
        return 0
    fi

    (time env NODE_NO_WARNINGS=1 node "${LOSAT_WASM_RUNNER}" "${LOSAT_WASM_BIN}" "$@" )&>"${log}"
}

run_losat_wasm_threaded() {
    local log="$1"
    shift

    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
    # ```c
    # int num_threads = args[kArgNumThreads].AsInteger();
    # if (num_threads > kMaxValue) {
    #     m_NumThreads = kMaxValue;
    # } else {
    #     m_NumThreads = num_threads;
    # }
    # ```
    if [ "${LOSAT_WASM_THREADED_AVAILABLE}" != "1" ]; then
        echo "Skipping ${log} because threaded LOSAT Wasm is unavailable."
        return 0
    fi

    (time env NODE_NO_WARNINGS=1 node "${LOSAT_WASM_THREADED_RUNNER}" "${LOSAT_WASM_THREADED_BIN}" "$@" )&>"${log}"
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
    run_losat_wasm_threaded "${threaded_log}" blastn -query "${FASTA_DIR}/${query}" -subject "${FASTA_DIR}/${subject}" -out "${threaded_out}" "${task_args[@]}" -num_threads "${LOSAT_THREADS}"
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
    run_losat_wasm_threaded \
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
    run_losat_wasm_threaded \
        "${LOSAT_OUT_DIR}/${stem}.tlosatx.wasm.${threaded_suffix}.log" \
        tblastx -query "${FASTA_DIR}/${query}" -subject "${FASTA_DIR}/${subject}" \
        -out "${LOSAT_OUT_DIR}/${stem}.tlosatx.wasm.${threaded_suffix}.out" --query-gencode "${query_gencode}" \
        --db-gencode "${db_gencode}" -num_threads "${LOSAT_THREADS}" -outfmt 6
}

echo "Starting LOSAT Wasm commands..."
if [ "${RUN_LOSAT_WASM}" = "0" ]; then
    echo "Skipping LOSAT Wasm commands because RUN_LOSAT_WASM=0."
else
    if [ "${BUILD_LOSAT_WASM}" = "1" ]; then
        (cd .. && cargo build --release --target wasm32-wasip1 --no-default-features)
    fi
    if [ "${BUILD_LOSAT_WASM_THREADED}" = "1" ]; then
        (cd .. && cargo build --release --target wasm32-wasip1-threads --features wasm-threads)
    fi

    if ! command -v node >/dev/null 2>&1; then
        echo "Skipping LOSAT Wasm commands because node is not available."
    else
        LOSAT_WASM_AVAILABLE=0
        if [ ! -f "${LOSAT_WASM_BIN}" ]; then
            echo "Skipping LOSAT Wasm n1 commands because ${LOSAT_WASM_BIN} does not exist."
        elif ! LOSAT_WASM_BIN_RESOLVED="$(resolve_losat_wasm_bin "${LOSAT_WASM_BIN}")"; then
            echo "Skipping LOSAT Wasm n1 commands because no _start-exporting LOSAT command Wasm was found near ${LOSAT_WASM_BIN}."
        else
            if [ "${LOSAT_WASM_BIN_RESOLVED}" != "${LOSAT_WASM_BIN}" ]; then
                echo "Using LOSAT Wasm command artifact: ${LOSAT_WASM_BIN_RESOLVED}"
                LOSAT_WASM_BIN="${LOSAT_WASM_BIN_RESOLVED}"
            fi
            # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
            # ```c
            # int num_threads = args[kArgNumThreads].AsInteger();
            # if (num_threads > kMaxValue) {
            #     m_NumThreads = kMaxValue;
            # } else {
            #     m_NumThreads = num_threads;
            # }
            # ```
            LOSAT_WASM_AVAILABLE=1
        fi

        LOSAT_WASM_THREADED_AVAILABLE=0
        if [ "${RUN_LOSAT_WASM_THREADED}" != "0" ]; then
            if [ ! -f "${LOSAT_WASM_THREADED_BIN}" ]; then
                echo "Skipping LOSAT Wasm n${LOSAT_THREADS} commands because ${LOSAT_WASM_THREADED_BIN} does not exist."
            elif ! LOSAT_WASM_THREADED_BIN_RESOLVED="$(resolve_losat_wasm_bin "${LOSAT_WASM_THREADED_BIN}")"; then
                echo "Skipping LOSAT Wasm n${LOSAT_THREADS} commands because no _start-exporting threaded LOSAT command Wasm was found near ${LOSAT_WASM_THREADED_BIN}."
            elif ! wasm_imports_wasi_thread_spawn "${LOSAT_WASM_THREADED_BIN_RESOLVED}"; then
                echo "Skipping LOSAT Wasm n${LOSAT_THREADS} commands because ${LOSAT_WASM_THREADED_BIN_RESOLVED} does not import wasi/thread-spawn."
            elif [ ! -f "${LOSAT_WASM_THREADED_RUNNER}" ]; then
                echo "Skipping LOSAT Wasm n${LOSAT_THREADS} commands because ${LOSAT_WASM_THREADED_RUNNER} does not exist."
            else
                if [ "${LOSAT_WASM_THREADED_BIN_RESOLVED}" != "${LOSAT_WASM_THREADED_BIN}" ]; then
                    echo "Using LOSAT threaded Wasm command artifact: ${LOSAT_WASM_THREADED_BIN_RESOLVED}"
                    LOSAT_WASM_THREADED_BIN="${LOSAT_WASM_THREADED_BIN_RESOLVED}"
                fi
                # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
                # ```c
                # int num_threads = args[kArgNumThreads].AsInteger();
                # if (num_threads > kMaxValue) {
                #     m_NumThreads = kMaxValue;
                # } else {
                #     m_NumThreads = num_threads;
                # }
                # ```
                LOSAT_WASM_THREADED_AVAILABLE=1
            fi
        fi
        write_losat_wasm_runner

        # run_losatn_wasm_case "NZ_CP006932.fasta" "NZ_CP006932.fasta" "NZ_CP006932.NZ_CP006932.losatn.megablast" "megablast"
        # run_losatn_wasm_case "EDL933.fna" "Sakai.fna" "EDL933.Sakai.losatn.megablast" "megablast"
        # run_losatn_wasm_case "Sakai.fna" "MG1655.fna" "Sakai.MG1655.losatn.megablast" "megablast"

        # run_losatn_wasm_case "NZ_CP006932.fasta" "NZ_CP006932.fasta" "NZ_CP006932.NZ_CP006932.losatn.blastn" "blastn"
        # run_losatn_wasm_case "AP027152.fasta" "AP027202.fasta" "PesePMNV.MjPMNV.losatn.blastn" "blastn"
        # run_losatn_wasm_case "LC738874.fasta" "LC738870.fasta" "MelaMJNV.PemoMJNVA.losatn.blastn" "blastn"
        # run_losatn_wasm_case "LC738884.fasta" "AP027155.fasta" "SiNMV.ChdeNMV.losatn.blastn" "blastn"
        # run_losatn_wasm_case "LC738869.fasta" "AP027202.fasta" "PmeNMV.MjPMNV.losatn.blastn" "blastn"
        # run_losatn_wasm_case "LC738869.fasta" "AP027152.fasta" "PmeNMV.PesePMNV.losatn.blastn" "blastn"
        # run_losatn_wasm_case "LC738873.fasta" "LC738871.fasta" "PeseMJNV.PemoMJNVB.losatn.blastn" "blastn"
        # run_losatn_wasm_case "LC738870.fasta" "LC738873.fasta" "PemoMJNVA.PeseMJNV.losatn.blastn" "blastn"
        # run_losatn_wasm_case "LC738868.fasta" "LC738874.fasta" "MjeNMV.MelaMJNV.losatn.blastn" "blastn"
        # run_losatn_wasm_case "AP027202.fasta" "LC738875.fasta" "MjPMNV.MlPMNV.losatn.blastn" "blastn"

        for losatp_case in "${LOSATP_CASES[@]}"; do
            IFS=":" read -r query subject stem <<<"$losatp_case"
            run_losatp_wasm_case "$query" "$subject" "$stem"
        done

        # for tblastx_case in "${TBLASTX_CASES[@]}"; do
        #     IFS=":" read -r query subject _db stem query_gencode db_gencode <<<"$tblastx_case"
        #     run_tlosatx_wasm_case "$query" "$subject" "$stem" "$query_gencode" "$db_gencode"
        # done
    fi
fi
echo "Finished LOSAT Wasm commands!"
