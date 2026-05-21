#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUT_DIR="${OUT_DIR:-"${SCRIPT_DIR}/losat_out"}"

QUERY="${QUERY:-"${SCRIPT_DIR}/fasta/small_test.fasta"}"
SUBJECT="${SUBJECT:-"${SCRIPT_DIR}/fasta/small_test.fasta"}"
STEM="${STEM:-small_test.small_test.tlosatx}"
QUERY_GENCODE="${QUERY_GENCODE:-1}"
DB_GENCODE="${DB_GENCODE:-1}"
OUTFMT="${OUTFMT:-6}"
THREADS="${THREADS:-1}"

WASM_TARGET="${WASM_TARGET:-wasm32-wasip1}"
WASM_FEATURES="${WASM_FEATURES:---no-default-features}"
BUILD_NATIVE="${BUILD_NATIVE:-0}"
BUILD_WASM_SIMD="${BUILD_WASM_SIMD:-0}"
BUILD_WASM_SCALAR="${BUILD_WASM_SCALAR:-0}"
ALLOW_LOSAT_INTERNAL_DIFFS="${ALLOW_LOSAT_INTERNAL_DIFFS:-0}"
REQUIRE_NCBI_TBLASTX_PARITY="${REQUIRE_NCBI_TBLASTX_PARITY:-0}"

NATIVE_BIN="${NATIVE_BIN:-"${CRATE_DIR}/target/release/LOSAT"}"
WASM_SIMD_BIN="${WASM_SIMD_BIN:-"${OUT_DIR}/LOSAT.wasm.simd.wasm"}"
WASM_SCALAR_BIN="${WASM_SCALAR_BIN:-"${OUT_DIR}/LOSAT.wasm.scalar.wasm"}"
NODE_RUNNER="${OUT_DIR}/run_losat_wasi.js"

mkdir -p "${OUT_DIR}"

case "${NATIVE_BIN}" in
    /*) ;;
    *) NATIVE_BIN="${CRATE_DIR}/${NATIVE_BIN}" ;;
esac

case "${WASM_SIMD_BIN}" in
    /*) ;;
    *) WASM_SIMD_BIN="${CRATE_DIR}/${WASM_SIMD_BIN}" ;;
esac

case "${WASM_SCALAR_BIN}" in
    /*) ;;
    *) WASM_SCALAR_BIN="${CRATE_DIR}/${WASM_SCALAR_BIN}" ;;
esac

if ! command -v node >/dev/null 2>&1; then
    echo "node is required for TBLASTX Wasm SIMD/scalar benchmarks." >&2
    exit 2
fi

cat >"${NODE_RUNNER}" <<'JS'
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

build_wasm_copy() {
    local output="$1"
    shift

    (
        cd "${CRATE_DIR}"
        cargo build --release --target "${WASM_TARGET}" "$@"
    )

    local built="${CRATE_DIR}/target/${WASM_TARGET}/release/LOSAT.wasm"
    local resolved
    if ! resolved="$(resolve_losat_wasm_bin "${built}")"; then
        echo "No _start-exporting command Wasm artifact found after build." >&2
        exit 2
    fi
    cp "${resolved}" "${output}"
}

if [ "${BUILD_NATIVE}" = "1" ]; then
    (cd "${CRATE_DIR}" && cargo build --release)
fi

if [ "${BUILD_WASM_SIMD}" = "1" ]; then
    build_wasm_copy "${WASM_SIMD_BIN}" ${WASM_FEATURES}
fi

if [ "${BUILD_WASM_SCALAR}" = "1" ]; then
    build_wasm_copy "${WASM_SCALAR_BIN}" ${WASM_FEATURES} --features tblastx-wasm-scalar
fi

if [ ! -x "${NATIVE_BIN}" ]; then
    echo "Missing native LOSAT binary: ${NATIVE_BIN}" >&2
    echo "Set BUILD_NATIVE=1 to build it." >&2
    exit 2
fi

if ! wasm_has_start_export "${WASM_SIMD_BIN}"; then
    echo "Missing command-WASI SIMD artifact: ${WASM_SIMD_BIN}" >&2
    echo "Set BUILD_WASM_SIMD=1 to build and copy it." >&2
    exit 2
fi

if ! wasm_has_start_export "${WASM_SCALAR_BIN}"; then
    echo "Missing command-WASI scalar artifact: ${WASM_SCALAR_BIN}" >&2
    echo "Set BUILD_WASM_SCALAR=1 to build and copy it." >&2
    exit 2
fi

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgQuery("query");
# const string kArgSubject("subject");
# const string kArgOutput("out");
# const string kArgNumThreads("num_threads");
# ```
run_native() {
    local out="$1"
    local log="$2"

    {
        printf '[WASM-SIMD-BENCH] execution_mode=native-release\n'
        printf '[WASM-SIMD-BENCH] requested_num_threads=%s\n' "${THREADS}"
        printf '[WASM-SIMD-BENCH] effective_engine_threads=%s\n' "${THREADS}"
        printf '[WASM-SIMD-BENCH] query=%s\n' "${QUERY}"
        printf '[WASM-SIMD-BENCH] subject=%s\n' "${SUBJECT}"
        time env LOSAT_TIMING="${LOSAT_TIMING:-1}" "${NATIVE_BIN}" tblastx \
            -query "${QUERY}" -subject "${SUBJECT}" -out "${out}" \
            --query-gencode "${QUERY_GENCODE}" --db-gencode "${DB_GENCODE}" \
            -num_threads "${THREADS}" -outfmt "${OUTFMT}"
    } &>"${log}"
}

# NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:846-921
# ```c
# status =
#     s_BlastAaExtendTwoHit(query, subject, word_params, ext_params,
#                           hit_params, init_hitlist, hsp_list);
# ```
run_wasm() {
    local mode="$1"
    local wasm_bin="$2"
    local out="$3"
    local log="$4"

    {
        printf '[WASM-SIMD-BENCH] execution_mode=command-wasi-serial\n'
        printf '[WASM-SIMD-BENCH] wasm_build_mode=%s\n' "${mode}"
        printf '[WASM-SIMD-BENCH] target_triple=%s\n' "${WASM_TARGET}"
        printf '[WASM-SIMD-BENCH] feature_set=%s%s\n' "${WASM_FEATURES}" "$([ "${mode}" = "scalar" ] && printf ' --features tblastx-wasm-scalar' || true)"
        printf '[WASM-SIMD-BENCH] requested_num_threads=%s\n' "${THREADS}"
        printf '[WASM-SIMD-BENCH] effective_engine_threads=1\n'
        printf '[WASM-SIMD-BENCH] benchmark_label=LOSAT Wasm %s\n' "${mode}"
        printf '[WASM-SIMD-BENCH] wasi_runtime=node-wasi-preview1\n'
        printf '[WASM-SIMD-BENCH] node_version=%s\n' "$(node --version 2>/dev/null || printf 'unavailable')"
        time env LOSAT_TIMING="${LOSAT_TIMING:-1}" NODE_NO_WARNINGS=1 node "${NODE_RUNNER}" "${wasm_bin}" \
            tblastx -query "${QUERY}" -subject "${SUBJECT}" -out "${out}" \
            --query-gencode "${QUERY_GENCODE}" --db-gencode "${DB_GENCODE}" \
            -num_threads "${THREADS}" -outfmt "${OUTFMT}"
    } &>"${log}"
}

NATIVE_OUT="${OUT_DIR}/${STEM}.native.out"
NATIVE_LOG="${OUT_DIR}/${STEM}.native.log"
SIMD_OUT="${OUT_DIR}/${STEM}.wasm.simd.out"
SIMD_LOG="${OUT_DIR}/${STEM}.wasm.simd.log"
SCALAR_OUT="${OUT_DIR}/${STEM}.wasm.scalar.out"
SCALAR_LOG="${OUT_DIR}/${STEM}.wasm.scalar.log"

run_native "${NATIVE_OUT}" "${NATIVE_LOG}"
run_wasm "SIMD" "${WASM_SIMD_BIN}" "${SIMD_OUT}" "${SIMD_LOG}"
run_wasm "scalar" "${WASM_SCALAR_BIN}" "${SCALAR_OUT}" "${SCALAR_LOG}"

run_ncbi=0
if command -v "${NCBI_TBLASTX:-tblastx}" >/dev/null 2>&1; then
    run_ncbi=1
    NCBI_OUT="${OUT_DIR}/${STEM}.ncbi.out"
    NCBI_LOG="${OUT_DIR}/${STEM}.ncbi.log"
    {
        printf '[WASM-SIMD-BENCH] execution_mode=ncbi-oracle\n'
        printf '[WASM-SIMD-BENCH] requested_num_threads=%s\n' "${THREADS}"
        time "${NCBI_TBLASTX:-tblastx}" -query "${QUERY}" -subject "${SUBJECT}" \
            -out "${NCBI_OUT}" -query_gencode "${QUERY_GENCODE}" \
            -db_gencode "${DB_GENCODE}" -num_threads "${THREADS}" -outfmt "${OUTFMT}"
    } &>"${NCBI_LOG}"
else
    printf 'NCBI tblastx not found; set NCBI_TBLASTX to capture oracle output.\n' \
        >"${OUT_DIR}/${STEM}.ncbi.SKIPPED"
fi

# NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
# ```c
# while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
#     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
#         continue;
#     }
# }
# ```
diff_pair() {
    local left_label="$1"
    local left="$2"
    local right_label="$3"
    local right="$4"
    local diff_file="${OUT_DIR}/${STEM}.${left_label}_vs_${right_label}.diff"

    if diff -u "${left}" "${right}" >"${diff_file}"; then
        printf '0\n'
    else
        wc -l <"${diff_file}"
    fi
}

native_vs_simd_lines="$(diff_pair native "${NATIVE_OUT}" wasm_simd "${SIMD_OUT}")"
native_vs_scalar_lines="$(diff_pair native "${NATIVE_OUT}" wasm_scalar "${SCALAR_OUT}")"
scalar_vs_simd_lines="$(diff_pair wasm_scalar "${SCALAR_OUT}" wasm_simd "${SIMD_OUT}")"

if [ "${run_ncbi}" = "1" ]; then
    ncbi_vs_native_lines="$(diff_pair ncbi "${NCBI_OUT}" native "${NATIVE_OUT}")"
    ncbi_vs_simd_lines="$(diff_pair ncbi "${NCBI_OUT}" wasm_simd "${SIMD_OUT}")"
    ncbi_vs_scalar_lines="$(diff_pair ncbi "${NCBI_OUT}" wasm_scalar "${SCALAR_OUT}")"
else
    ncbi_vs_native_lines="skipped"
    ncbi_vs_simd_lines="skipped"
    ncbi_vs_scalar_lines="skipped"
fi

{
    printf 'out_dir\t%s\n' "${OUT_DIR}"
    printf 'query\t%s\n' "${QUERY}"
    printf 'subject\t%s\n' "${SUBJECT}"
    printf 'query_gencode\t%s\n' "${QUERY_GENCODE}"
    printf 'db_gencode\t%s\n' "${DB_GENCODE}"
    printf 'outfmt\t%s\n' "${OUTFMT}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'native_log\t%s\n' "${NATIVE_LOG}"
    printf 'wasm_simd_log\t%s\n' "${SIMD_LOG}"
    printf 'wasm_scalar_log\t%s\n' "${SCALAR_LOG}"
    printf 'native_vs_wasm_simd_raw_diff_lines\t%s\n' "${native_vs_simd_lines}"
    printf 'native_vs_wasm_scalar_raw_diff_lines\t%s\n' "${native_vs_scalar_lines}"
    printf 'wasm_scalar_vs_wasm_simd_raw_diff_lines\t%s\n' "${scalar_vs_simd_lines}"
    printf 'ncbi_vs_native_raw_diff_lines\t%s\n' "${ncbi_vs_native_lines}"
    printf 'ncbi_vs_wasm_simd_raw_diff_lines\t%s\n' "${ncbi_vs_simd_lines}"
    printf 'ncbi_vs_wasm_scalar_raw_diff_lines\t%s\n' "${ncbi_vs_scalar_lines}"
} | tee "${OUT_DIR}/${STEM}.wasm-simd-modes.summary.tsv"

losat_internal_diff=$(( native_vs_simd_lines + native_vs_scalar_lines + scalar_vs_simd_lines ))
if [ "${losat_internal_diff}" -ne 0 ] && [ "${ALLOW_LOSAT_INTERNAL_DIFFS}" != "1" ]; then
    echo "LOSAT native/Wasm scalar/Wasm SIMD outputs differ; inspect ${OUT_DIR}/*.diff" >&2
    exit 1
fi

if [ "${REQUIRE_NCBI_TBLASTX_PARITY}" = "1" ] && [ "${run_ncbi}" = "1" ]; then
    if [ "${ncbi_vs_native_lines}" -ne 0 ] || [ "${ncbi_vs_simd_lines}" -ne 0 ] || [ "${ncbi_vs_scalar_lines}" -ne 0 ]; then
        echo "NCBI oracle output differs; inspect ${OUT_DIR}/*.diff" >&2
        exit 1
    fi
fi
