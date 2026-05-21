#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUT_DIR="${OUT_DIR:-"${SCRIPT_DIR}/losat_out/wasm_serial_modes"}"

PROGRAM="${PROGRAM:-tblastx}"
QUERY="${QUERY:-"${SCRIPT_DIR}/fasta/AvCLPV.fasta"}"
SUBJECT="${SUBJECT:-"${SCRIPT_DIR}/fasta/PsCLPV.fasta"}"
OUTFMT="${OUTFMT:-6}"

COMMAND_WASM_TARGET="${COMMAND_WASM_TARGET:-wasm32-wasip1}"
WEB_WASM_TARGET="${WEB_WASM_TARGET:-wasm32-unknown-unknown}"
COMMAND_WASM_BIN="${COMMAND_WASM_BIN:-"${CRATE_DIR}/target/${COMMAND_WASM_TARGET}/release/LOSAT.wasm"}"
WEB_WASM_BIN="${WEB_WASM_BIN:-"${CRATE_DIR}/target/${WEB_WASM_TARGET}/release/LOSAT.wasm"}"
WASM_FEATURES="${WASM_FEATURES:---no-default-features}"
BUILD_COMMAND_WASM="${BUILD_COMMAND_WASM:-0}"
BUILD_WEB_WASM="${BUILD_WEB_WASM:-0}"

mkdir -p "${OUT_DIR}"

# NCBI reference: ncbi-blast/c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
# ```c
# #ifdef NCBI_NO_THREADS
#     m_NumThreads = CThreadable::kMinNumThreads;
#     m_MTMode = eNotSupported;
# #endif
# ```
wasm_rayon_compiled() {
    case " ${WASM_FEATURES} " in
        *" parallel "*|*"--features parallel"*|*"--features=parallel"*)
            printf 'true\n'
            ;;
        *)
            printf 'false\n'
            ;;
    esac
}

case "${PROGRAM}" in
    tblastx)
        EXTRA_ARGS=(--query-gencode "${QUERY_GENCODE:-1}" --db-gencode "${DB_GENCODE:-1}")
        ;;
    blastn|blastp)
        EXTRA_ARGS=()
        ;;
    *)
        echo "Unsupported PROGRAM=${PROGRAM}; expected tblastx, blastn, or blastp." >&2
        exit 2
        ;;
esac

if [ "${BUILD_COMMAND_WASM}" = "1" ]; then
    (cd "${CRATE_DIR}" && cargo build --release --target "${COMMAND_WASM_TARGET}" ${WASM_FEATURES})
fi

if [ "${BUILD_WEB_WASM}" = "1" ]; then
    (cd "${CRATE_DIR}" && cargo build --release --lib --target "${WEB_WASM_TARGET}" ${WASM_FEATURES})
fi

if ! command -v node >/dev/null 2>&1; then
    echo "node is required for Wasm serial mode benchmarks." >&2
    exit 2
fi

COMMAND_RUNNER="${OUT_DIR}/run_losat_wasi.js"
WEB_RUNNER="${OUT_DIR}/run_losat_web_pair.js"

cat >"${COMMAND_RUNNER}" <<'JS'
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

cat >"${WEB_RUNNER}" <<'JS'
const fs = require('fs');

const [wasmPath, program, queryPath, subjectPath, outPath, outfmt, ...extraArgs] =
  process.argv.slice(2);

function encoderBytes(value) {
  return new TextEncoder().encode(value);
}

(async () => {
  const imports = { env: {} };
  const bytes = fs.readFileSync(wasmPath);
  const { instance } = await WebAssembly.instantiate(bytes, imports);
  const e = instance.exports;
  const memory = e.memory;

  for (const name of [
    'losat_web_alloc',
    'losat_web_dealloc',
    'losat_web_run_pair',
    'losat_web_result_ptr',
    'losat_web_result_len',
    'losat_web_error_ptr',
    'losat_web_error_len',
  ]) {
    if (typeof e[name] !== 'function') {
      throw new Error(`${wasmPath} does not export ${name}; use the LOSAT web_api Wasm artifact.`);
    }
  }

  function put(bytes) {
    const ptr = e.losat_web_alloc(bytes.length);
    new Uint8Array(memory.buffer, ptr, bytes.length).set(bytes);
    return { ptr, len: bytes.length };
  }

  function read(ptr, len) {
    return new Uint8Array(memory.buffer, ptr, len).slice();
  }

  const programMem = put(encoderBytes(program));
  const queryMem = put(encoderBytes(fs.readFileSync(queryPath, 'utf8')));
  const subjectMem = put(encoderBytes(fs.readFileSync(subjectPath, 'utf8')));
  const outfmtMem = put(encoderBytes(outfmt || '6'));
  const extraMem = put(encoderBytes(extraArgs.join('\0')));

  const started = process.hrtime.bigint();
  const status = e.losat_web_run_pair(
    programMem.ptr,
    programMem.len,
    queryMem.ptr,
    queryMem.len,
    subjectMem.ptr,
    subjectMem.len,
    outfmtMem.ptr,
    outfmtMem.len,
    extraMem.ptr,
    extraMem.len,
  );
  const elapsed = Number(process.hrtime.bigint() - started) / 1e9;

  if (status !== 0) {
    const error = new TextDecoder().decode(read(e.losat_web_error_ptr(), e.losat_web_error_len()));
    throw new Error(error || `losat_web_run_pair returned ${status}`);
  }

  fs.writeFileSync(outPath, read(e.losat_web_result_ptr(), e.losat_web_result_len()));
  console.error(`[TIMING] web_api_in_memory_total: ${elapsed.toFixed(3)}s`);

  for (const item of [programMem, queryMem, subjectMem, outfmtMem, extraMem]) {
    e.losat_web_dealloc(item.ptr, item.len);
  }
})().catch((err) => {
  console.error(err);
  process.exit(1);
});
JS

case "${COMMAND_WASM_BIN}" in
    /*) ;;
    *) COMMAND_WASM_BIN="${CRATE_DIR}/${COMMAND_WASM_BIN}" ;;
esac

case "${WEB_WASM_BIN}" in
    /*) ;;
    *) WEB_WASM_BIN="${CRATE_DIR}/${WEB_WASM_BIN}" ;;
esac

if [ ! -f "${COMMAND_WASM_BIN}" ]; then
    echo "Missing command-WASI artifact: ${COMMAND_WASM_BIN}" >&2
    echo "Set BUILD_COMMAND_WASM=1 to build it." >&2
    exit 2
fi

if [ ! -f "${WEB_WASM_BIN}" ]; then
    echo "Missing browser/web_api artifact: ${WEB_WASM_BIN}" >&2
    echo "Set BUILD_WEB_WASM=1 to build it." >&2
    exit 2
fi

COMMAND_OUT="${OUT_DIR}/${PROGRAM}.command-wasi.out"
COMMAND_LOG="${OUT_DIR}/${PROGRAM}.command-wasi.log"
WEB_OUT="${OUT_DIR}/${PROGRAM}.browser-in-memory.out"
WEB_LOG="${OUT_DIR}/${PROGRAM}.browser-in-memory.log"

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgQuery("query");
# const string kArgSubject("subject");
# const string kArgOutput("out");
# const string kArgNumThreads("num_threads");
# ```
{
    printf '[WASM-METADATA] execution_mode=command-wasi-serial\n'
    printf '[WASM-METADATA] target_triple=%s\n' "${COMMAND_WASM_TARGET}"
    printf '[WASM-METADATA] feature_set=%s\n' "${WASM_FEATURES}"
    printf '[WASM-METADATA] rayon_compiled=%s\n' "$(wasm_rayon_compiled)"
    printf '[WASM-METADATA] requested_num_threads=1\n'
    printf '[WASM-METADATA] effective_engine_threads=1\n'
    printf '[WASM-METADATA] benchmark_label=LOSAT command-WASI serial\n'
    printf '[WASM-METADATA] browser_worker_parallel_verified=false\n'
    printf '[WASM-METADATA] wasi_runtime=node-wasi-preview1\n'
    printf '[WASM-METADATA] node_version=%s\n' "$(node --version 2>/dev/null || printf 'unavailable')"
    time env NODE_NO_WARNINGS=1 node "${COMMAND_RUNNER}" "${COMMAND_WASM_BIN}" \
        "${PROGRAM}" -query "${QUERY}" -subject "${SUBJECT}" -out "${COMMAND_OUT}" \
        "${EXTRA_ARGS[@]}" -num_threads 1 -outfmt "${OUTFMT}"
} &>"${COMMAND_LOG}"

# NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
# ```c
# db_length = BlastSeqSrcGetTotLen(seq_src);
# while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
#     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
#         continue;
#     }
# }
# ```
{
    printf '[WASM-METADATA] execution_mode=browser-in-memory-serial\n'
    printf '[WASM-METADATA] target_triple=%s\n' "${WEB_WASM_TARGET}"
    printf '[WASM-METADATA] feature_set=%s\n' "${WASM_FEATURES}"
    printf '[WASM-METADATA] rayon_compiled=%s\n' "$(wasm_rayon_compiled)"
    printf '[WASM-METADATA] requested_num_threads=1\n'
    printf '[WASM-METADATA] effective_engine_threads=1\n'
    printf '[WASM-METADATA] benchmark_label=LOSAT browser in-memory serial\n'
    printf '[WASM-METADATA] browser_worker_parallel_verified=false\n'
    printf '[WASM-METADATA] node_version=%s\n' "$(node --version 2>/dev/null || printf 'unavailable')"
    time env NODE_NO_WARNINGS=1 node "${WEB_RUNNER}" "${WEB_WASM_BIN}" \
        "${PROGRAM}" "${QUERY}" "${SUBJECT}" "${WEB_OUT}" "${OUTFMT}" "${EXTRA_ARGS[@]}"
} &>"${WEB_LOG}"

if diff -u "${COMMAND_OUT}" "${WEB_OUT}" >"${OUT_DIR}/${PROGRAM}.serial-modes.diff"; then
    printf 'serial_mode_output_diff=0\n' | tee "${OUT_DIR}/${PROGRAM}.serial-modes.summary"
else
    printf 'serial_mode_output_diff=1\n' | tee "${OUT_DIR}/${PROGRAM}.serial-modes.summary"
    echo "Output differs; inspect ${OUT_DIR}/${PROGRAM}.serial-modes.diff" >&2
    exit 1
fi

printf 'command_wasi_log=%s\n' "${COMMAND_LOG}" | tee -a "${OUT_DIR}/${PROGRAM}.serial-modes.summary"
printf 'browser_in_memory_log=%s\n' "${WEB_LOG}" | tee -a "${OUT_DIR}/${PROGRAM}.serial-modes.summary"
