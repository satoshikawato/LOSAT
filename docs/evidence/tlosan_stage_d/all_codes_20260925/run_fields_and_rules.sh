#!/usr/bin/env bash
set -euo pipefail
# Comparison-only NCBI API oracle. It is never invoked by LOSAT builds/runtime.
# Pinned source: c++/src/algo/blast/core/blast_engine.c:1460-1466;
# api/seqsrc_multiseq.cpp:140-164; core/gencode_singleton.c:FindGeneticCode.
root=$(cd "$(dirname "$0")" && pwd)
stage_d=$(cd "$root/.." && pwd)
stage_a=$(cd "$root/../../tlosan_stage_a" && pwd)
out=${1:?usage: run_ncbi_all_codes.sh NEW_OUTPUT_DIRECTORY}
mkdir "$out"
out=$(cd "$out" && pwd)
source_repo=${NCBI_SOURCE_REPO:-/mnt/c/Users/genom/GitHub/ncbi-blast}
source_commit=598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4
blast_lib=${NCBI_BLAST_LIB:-/home/kawato/micromamba/lib/ncbi-blast+}
datatool_bin=${NCBI_DATATOOL_BIN:-/home/kawato/micromamba/bin/datatool}
tblastn_bin=${TBLASTN_BIN:-/home/kawato/micromamba/bin/tblastn}
[[ $(git -C "$source_repo" rev-parse HEAD) == "$source_commit" ]]
[[ $(sha256sum "$tblastn_bin" | cut -d' ' -f1) == e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0 ]]
(cd "$root" && sha256sum -c fixtures.sha256) > "$out/fixtures_verified.log"
work=$(mktemp -d /tmp/tlosan-allcodes-XXXXXX)
trap 'rm -rf "$work"' EXIT
mkdir "$work/source"
git -C "$source_repo" archive "$source_commit" | tar -x -C "$work/source"
(cd "$work/source" && sha256sum -c "$stage_a/ncbi_source.sha256") > "$out/source_verified.log"
src="$work/source/c++"
chmod +x "$src/src/build-system/configure" "$src/src/build-system/config.sub" "$src/src/build-system/config.guess" "$src/scripts/common/impl/get_lock.sh"
find "$src/scripts" -type f \( -name '*.sh' -o -name '*.awk' \) -exec chmod +x {} +
(cd "$work" && bash "$src/configure.orig" --with-build-root="$work/build" --without-debug --with-mt --without-vdb --without-gnutls --without-gcrypt --with-experimental=Int8GI) > "$out/configure.log" 2>&1
find "$work/build/build" -type f -name '*.sh' -exec chmod +x {} +
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects/seq" sources > "$out/seq_sources.log" 2>&1
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects" sources > "$out/objects_sources.log" 2>&1
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects/genomecoll" sources > "$out/genomecoll_sources.log" 2>&1
g++ -std=c++17 -O2 -I"$work/build/inc" -I"$src/include" \
  "$stage_d/tblastn_code32_local_oracle.cpp" -L"$blast_lib" -Wl,-rpath,"$blast_lib" \
  -lxblastformat -lalign_format -lxformat -lxcleanup -lgbseq -lblastinput \
  -lxblast -lblast -lseqdb -lxobjmgr -lseqset -lseq -lgeneral -lxser -lxncbi \
  -o "$work/tblastn_local_api_oracle"
gcc -std=c11 -O0 -shared -fPIC "$stage_d/ncbi_kappa_traceback_trace.c" -ldl -o "$work/kappa_trace.so"
gcc -std=c11 -O0 -shared -fPIC "$stage_d/ncbi_d_call_trace.c" -ldl -o "$work/d_trace.so"
# Derive a comparison-only custom-fields formatter from the same CLI-calibrated
# API search source. Only the CBlastFormat constructor receives the custom spec.
python3 "$root/make_fields_oracle.py" "$stage_d/tblastn_code32_local_oracle.cpp" "$out/fields_oracle.cpp"
g++ -std=c++17 -O2 -I"$work/build/inc" -I"$src/include" \
  "$out/fields_oracle.cpp" -L"$blast_lib" -Wl,-rpath,"$blast_lib" \
  -lxblastformat -lalign_format -lxformat -lxcleanup -lgbseq -lblastinput \
  -lxblast -lblast -lseqdb -lxobjmgr -lseqset -lseq -lgeneral -lxser -lxncbi \
  -o "$work/tblastn_fields_oracle"
gcc -std=c11 -O0 -shared -fPIC "$stage_d/ncbi_kappa_rule_trace.c" -ldl -o "$work/rule_trace.so"
python3 "$root/run_fields_and_rules.py" "$out" "$work/tblastn_fields_oracle" "$work/rule_trace.so" "$work/d_trace.so" "$tblastn_bin"
