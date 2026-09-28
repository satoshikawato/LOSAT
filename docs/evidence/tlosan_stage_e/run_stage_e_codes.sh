#!/usr/bin/env bash
set -euo pipefail
# Comparison-only API search and formatting oracle. Never used by LOSAT
# runtime/build/fallback. NCBI source: tblastn_app.cpp:242-265,288-301;
# blast_setup_cxx.cpp:800-812; blast_format.cpp:759-835,1411-1458.
root=$(cd "$(dirname "$0")" && pwd)
stage_a=$(cd "$root/../tlosan_stage_a" && pwd)
out=${1:?usage: run_stage_e_codes.sh NEW_OUTPUT_DIRECTORY}
source_repo=${NCBI_SOURCE_REPO:-/mnt/c/Users/genom/GitHub/ncbi-blast}
source_commit=598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4
blast_lib=${NCBI_BLAST_LIB:-/home/kawato/micromamba/lib/ncbi-blast+}
datatool_bin=${NCBI_DATATOOL_BIN:-/home/kawato/micromamba/bin/datatool}
[[ $(git -C "$source_repo" rev-parse HEAD) == "$source_commit" ]]
mkdir "$out"
out=$(cd "$out" && pwd)
work=$(mktemp -d /tmp/tlosan-stagee-codes-XXXXXX)
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
  "$root/tblastn_stage_e_local_oracle.cpp" -L"$blast_lib" -Wl,-rpath,"$blast_lib" \
  -lxblastformat -lalign_format -lxformat -lxcleanup -lgbseq -lblastinput \
  -lxblast -lblast -lseqdb -lxobjmgr -lseqset -lseq -lgeneral -lxser -lxncbi \
  -o "$work/tblastn_stage_e_local_oracle"
cp "$work/tblastn_stage_e_local_oracle" "$out/tblastn_stage_e_local_oracle"
python3 "$root/run_stage_e_codes.py" "$out/tblastn_stage_e_local_oracle" "$out/comparison"
