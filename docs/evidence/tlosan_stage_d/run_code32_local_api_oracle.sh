#!/usr/bin/env bash
set -euo pipefail

# Comparison-only pinned NCBI C++ local-subject oracle. Never used by LOSAT.
root=$(cd "$(dirname "$0")" && pwd)
stage_a="$root/../tlosan_stage_a"
out=${1:?usage: run_code32_local_api_oracle.sh NEW_OUTPUT_DIRECTORY}
mkdir "$out"
out=$(cd "$out" && pwd)
source_repo=${NCBI_SOURCE_REPO:-/mnt/c/Users/genom/GitHub/ncbi-blast}
source_commit=598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4
blast_lib=${NCBI_BLAST_LIB:-/home/kawato/micromamba/lib/ncbi-blast+}
datatool_bin=${NCBI_DATATOOL_BIN:-/home/kawato/micromamba/bin/datatool}
tblastn_bin=${TBLASTN_BIN:-/home/kawato/micromamba/bin/tblastn}
[[ $(git -C "$source_repo" rev-parse HEAD) == "$source_commit" ]]
[[ $(sha256sum "$tblastn_bin" | cut -d' ' -f1) == e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0 ]]

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir "$work/source"
git -C "$source_repo" archive "$source_commit" | tar -x -C "$work/source"
(cd "$work/source" && sha256sum -c "$stage_a/ncbi_source.sha256") > "$out/source_verified.log"
src="$work/source/c++"
chmod +x "$src/src/build-system/configure" "$src/src/build-system/config.sub" "$src/src/build-system/config.guess" "$src/scripts/common/impl/get_lock.sh"
find "$src/scripts" -type f \( -name '*.sh' -o -name '*.awk' \) -exec chmod +x {} +
(
  cd "$work"
  bash "$src/configure.orig" --with-build-root="$work/build" --without-debug --with-mt \
    --without-vdb --without-gnutls --without-gcrypt --with-experimental=Int8GI
) > "$out/configure.log" 2>&1
find "$work/build/build" -type f -name '*.sh' -exec chmod +x {} +
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects/seq" sources > "$out/seq_sources.log" 2>&1
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects" sources > "$out/objects_sources.log" 2>&1
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects/genomecoll" sources > "$out/genomecoll_sources.log" 2>&1

# NCBI api/blast_setup_cxx.cpp:800-812 selects the subject's code string.
# core/blast_engine.c:1460-1466 consumes a code supplied by BlastSeqSrc.
g++ -std=c++17 -O2 -I"$work/build/inc" -I"$src/include" \
  "$root/tblastn_code32_local_oracle.cpp" -L"$blast_lib" -Wl,-rpath,"$blast_lib" \
  -lxblastformat -lalign_format -lxformat -lxcleanup -lgbseq -lblastinput \
  -lxblast -lblast -lseqdb -lxobjmgr -lseqset -lseq -lgeneral -lxser -lxncbi \
  -o "$work/tblastn_code32_local_oracle"
gcc -std=c11 -O0 -shared -fPIC "$root/ncbi_kappa_traceback_trace.c" -ldl -o "$work/kappa_trace.so"
gcc -std=c11 -O0 -shared -fPIC "$root/ncbi_kappa_mode2_trace.c" -ldl -o "$work/mode2_trace.so"
gcc -std=c11 -O0 -shared -fPIC "$root/ncbi_kappa_composition_trace.c" -ldl -o "$work/composition_trace.so"
for code in 1 32; do
  "$work/tblastn_code32_local_oracle" \
    "$stage_a/fixtures/query.faa" "$stage_a/fixtures/subject_code32.fna" "$code" 6 \
    > "$out/code${code}_fmt6.out" 2> "$out/code${code}.stderr"
  LD_PRELOAD="$work/kappa_trace.so" "$work/tblastn_code32_local_oracle" \
    "$stage_a/fixtures/query.faa" "$stage_a/fixtures/subject_code32.fna" "$code" 6 \
    > "$out/code${code}_traced.out" 2> "$out/code${code}_traced.stderr"
  python3 - "$out" "$code" <<'PY'
from pathlib import Path
import sys
out, code = Path(sys.argv[1]), sys.argv[2]
prefix = b"K_TRACE_"
plain_out = (out / f"code{code}_fmt6.out").read_bytes()
traced_out = (out / f"code{code}_traced.out").read_bytes()
plain_err = (out / f"code{code}.stderr").read_bytes()
traced_err = (out / f"code{code}_traced.stderr").read_bytes()
selected = [line for line in traced_err.splitlines(keepends=True) if line.startswith(prefix)]
other = b"".join(line for line in traced_err.splitlines(keepends=True) if not line.startswith(prefix))
assert plain_out == traced_out
assert plain_err == other
assert any(line.startswith(b"K_TRACE_HEAP_HSP\t") for line in selected)
(out / f"code{code}.trace").write_bytes(b"".join(selected))
PY
  for probe in mode2 composition; do
    LD_PRELOAD="$work/${probe}_trace.so" "$work/tblastn_code32_local_oracle" \
      "$stage_a/fixtures/query.faa" "$stage_a/fixtures/subject_code32.fna" "$code" 6 \
      > "$out/code${code}_${probe}.out" 2> "$out/code${code}_${probe}.stderr"
    python3 - "$out" "$code" "$probe" <<'PY'
from pathlib import Path
import sys
out, code, probe = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
prefixes = (
    (b"K_CALL\t", b"K_ALIGN\t", b"K_TRANSLATED\t") if probe == "mode2"
    else (b"K_COMP_CALL\t", b"K_MATRIX\t", b"K_COMP\t", b"K_COMP_RESULT\t", b"K_ADJUSTED\t")
)
traced_out = (out / f"code{code}_{probe}.out").read_bytes()
traced_err = (out / f"code{code}_{probe}.stderr").read_bytes()
plain_out = (out / f"code{code}_fmt6.out").read_bytes()
plain_err = (out / f"code{code}.stderr").read_bytes()
selected = [line for line in traced_err.splitlines(keepends=True) if line.startswith(prefixes)]
other = b"".join(line for line in traced_err.splitlines(keepends=True) if not line.startswith(prefixes))
assert traced_out == plain_out
assert other == plain_err
assert selected
(out / f"code{code}_{probe}.trace").write_bytes(b"".join(selected))
PY
  done
done
"$tblastn_bin" -task tblastn -query "$stage_a/fixtures/query.faa" \
  -subject "$stage_a/fixtures/subject_code32.fna" -db_gencode 1 -num_threads 1 \
  -outfmt 6 > "$out/cli_code1_fmt6.out" 2> "$out/cli_code1.stderr"
cmp "$out/cli_code1_fmt6.out" "$out/code1_fmt6.out"
LD_PRELOAD="$work/kappa_trace.so" "$tblastn_bin" -task tblastn \
  -query "$stage_a/fixtures/query.faa" \
  -subject "$stage_a/fixtures/subject_code32.fna" -db_gencode 1 -num_threads 1 \
  -outfmt 6 > "$out/cli_code1_traced.out" 2> "$out/cli_code1_traced.stderr"
python3 - "$out" <<'PY'
from pathlib import Path
import sys
out = Path(sys.argv[1])
traced = (out / "cli_code1_traced.stderr").read_bytes()
prefix = b"K_TRACE_"
selected = b"".join(x for x in traced.splitlines(keepends=True) if x.startswith(prefix))
ordinary = b"".join(x for x in traced.splitlines(keepends=True) if not x.startswith(prefix))
assert (out / "cli_code1_fmt6.out").read_bytes() == (out / "cli_code1_traced.out").read_bytes()
assert (out / "cli_code1.stderr").read_bytes() == ordinary
assert (out / "code1.trace").read_bytes() == selected
(out / "cli_code1.trace").write_bytes(selected)
PY

printf 'NCBI source commit: %s\nAPI mode: local subject adapter with explicit FindGeneticCode on BlastSeqSrc retrieval\nCode 1 API output and Kappa call trace match local -subject CLI bytes.\n' "$source_commit" > "$out/manifest.txt"
sha256sum "$root/tblastn_code32_local_oracle.cpp" "$root/ncbi_kappa_traceback_trace.c" \
  "$root/ncbi_kappa_mode2_trace.c" "$root/ncbi_kappa_composition_trace.c" \
  "$stage_a/fixtures/query.faa" "$stage_a/fixtures/subject_code32.fna" \
  "$work/tblastn_code32_local_oracle" "$tblastn_bin" > "$out/inputs.sha256"
(cd "$out" && sha256sum ./*.out ./*.trace > outputs.sha256)
