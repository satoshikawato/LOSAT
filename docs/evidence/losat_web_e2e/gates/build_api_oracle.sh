#!/usr/bin/env bash
# Build TLOSAN Stage E's comparison-only NCBI C++ API local-subject oracle
# (docs/evidence/tlosan_stage_e/tblastn_stage_e_local_oracle.cpp, never used by LOSAT)
# into a kept directory, for the TBLASTN non-default -db_gencode rows of the S08+ sweep
# (AGENTS.md: PD-TLOSAN-LOCAL-GENCODE-32 and the approved local-subject db_gencode
# behavior). Same steps as docs/evidence/tlosan_stage_e/run_stage_e_codes.sh.
set -euo pipefail
W=$(cd "$(dirname "$0")/../../../.." && pwd)
out=${1:?usage: build_api_oracle.sh OUTPUT_DIRECTORY}
source_repo=${NCBI_SOURCE_REPO:-/mnt/c/Users/genom/GitHub/ncbi-blast}
source_commit=598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4
blast_lib=${NCBI_BLAST_LIB:-/home/kawato/micromamba/lib/ncbi-blast+}
datatool_bin=${NCBI_DATATOOL_BIN:-/home/kawato/micromamba/bin/datatool}
[[ $(git -C "$source_repo" rev-parse HEAD) == "$source_commit" ]]
mkdir -p "$out"
out=$(cd "$out" && pwd)
work="$out/work"
rm -rf "$work"; mkdir -p "$work/source"
git -C "$source_repo" archive "$source_commit" | tar -x -C "$work/source"
(cd "$work/source" && sha256sum -c "$W/docs/evidence/tlosan_stage_a/ncbi_source.sha256") > "$out/source_verified.log"
src="$work/source/c++"
chmod +x "$src/src/build-system/configure" "$src/src/build-system/config.sub" "$src/src/build-system/config.guess" "$src/scripts/common/impl/get_lock.sh"
find "$src/scripts" -type f \( -name '*.sh' -o -name '*.awk' \) -exec chmod +x {} +
(cd "$work" && bash "$src/configure.orig" --with-build-root="$work/build" --without-debug --with-mt --without-vdb --without-gnutls --without-gcrypt --with-experimental=Int8GI) > "$out/configure.log" 2>&1
find "$work/build/build" -type f -name '*.sh' -exec chmod +x {} +
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects/seq" sources > "$out/seq_sources.log" 2>&1
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects" sources > "$out/objects_sources.log" 2>&1
PREBUILT_DATATOOL_EXE="$datatool_bin" make -C "$work/build/build/objects/genomecoll" sources > "$out/genomecoll_sources.log" 2>&1
g++ -std=c++17 -O2 -I"$work/build/inc" -I"$src/include" \
  "$W/docs/evidence/tlosan_stage_e/tblastn_stage_e_local_oracle.cpp" -L"$blast_lib" -Wl,-rpath,"$blast_lib" \
  -lxblastformat -lalign_format -lxformat -lxcleanup -lgbseq -lblastinput \
  -lxblast -lblast -lseqdb -lxobjmgr -lseqset -lseq -lgeneral -lxser -lxncbi \
  -o "$out/tblastn_stage_e_local_oracle"
rm -rf "$work"
sha256sum "$out/tblastn_stage_e_local_oracle" "$W/docs/evidence/tlosan_stage_e/tblastn_stage_e_local_oracle.cpp" > "$out/oracle.sha256"
cat "$out/oracle.sha256"
