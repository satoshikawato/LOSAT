#!/usr/bin/env bash
set -euo pipefail

# Comparison-only NCBI API build. Never call from a LOSAT build/runtime path.
# NCBI c++/src/algo/blast/api/blast_aux.cpp:588-613,629-646:
#   FindGeneticCode(id) reads the NCBI genetic-code table;
#   CAutomaticGenCodeSingleton registers the selected code.
# NCBI c++/src/app/blast/tblastn_app.cpp:189-193,287-305:
#   InitializeSubject -> CLocalBlast::Run -> CBlastFormat::PrintOneResultSet.
root=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: run_api_oracle.sh NEW_OUTPUT_DIRECTORY}
mkdir "$out"
out=$(cd "$out" && pwd)
source_repo=${NCBI_SOURCE_REPO:-/mnt/c/Users/genom/GitHub/ncbi-blast}
source_commit=598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4
blast_lib=${NCBI_BLAST_LIB:-/home/kawato/micromamba/lib/ncbi-blast+}
datatool_bin=${NCBI_DATATOOL_BIN:-/home/kawato/micromamba/bin/datatool}
tblastn_bin=${TBLASTN_BIN:-/home/kawato/micromamba/bin/tblastn}
makeblastdb_bin=${MAKEBLASTDB_BIN:-/home/kawato/micromamba/bin/makeblastdb}
[[ $(sha256sum "$tblastn_bin" | cut -d' ' -f1) == e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0 ]]
[[ $(sha256sum "$makeblastdb_bin" | cut -d' ' -f1) == 407bf59ec36ee6444104a90f28eab282e374a0dfe4c7ffee89a65d82dc51c4bf ]]
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir "$work/source"
git -C "$source_repo" archive "$source_commit" | tar -x -C "$work/source"
(cd "$work/source" && sha256sum -c "$root/ncbi_source.sha256") > "$out/source_verified.log"
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

# NCBI c++/src/algo/blast/api/blast_setup_cxx.cpp:800-810:
#   the subject's genetic-code ID selects FindGeneticCode for translation.
g++ -std=c++17 -O2 -I"$work/build/inc" -I"$src/include" \
  "$root/oracle/tblastn_code32_oracle.cpp" -L"$blast_lib" -Wl,-rpath,"$blast_lib" \
  -lxblastformat -lalign_format -lxformat -lxcleanup -lgbseq -lblastinput \
  -lxblast -lblast -lseqdb -lxobjmgr -lseqset -lseq -lgeneral -lxser -lxncbi \
  -o "$work/tblastn_code32_oracle"

"$tblastn_bin" -version > "$out/tblastn.version"
"$datatool_bin" -version > "$out/datatool.version" 2>&1 || true
sha256sum "$tblastn_bin" "$makeblastdb_bin" "$datatool_bin" \
  "$root/oracle/tblastn_code32_oracle.cpp" "$root"/fixtures/* \
  "$work/tblastn_code32_oracle" > "$out/inputs.sha256"
for lib in xblastformat align_format xformat xcleanup gbseq blastinput \
  xblast blast seqdb xobjmgr seqset seq general xser xncbi; do
  sha256sum "$blast_lib/lib${lib}.so"
done > "$out/libraries.sha256"
printf '%s\n' "NCBI source commit: $source_commit" \
  'target: x86_64-unknown-linux-gnu, NCBI C++ API DB search, one thread' \
  'options: CTBlastnOptionsHandle defaults; SEG 12 2.2 2.5; soft_masking=false; SetDbGeneticCode(code); outfmt 0/6/7' \
  > "$out/manifest.txt"
for code in 1 4 32; do
  "$makeblastdb_bin" -in "$root/fixtures/subject_code${code}.fna" -dbtype nucl \
    -parse_seqids -out "$work/db_${code}" > "$out/makeblastdb_code${code}.log"
  for fmt in 0 6 7; do
    "$work/tblastn_code32_oracle" "$root/fixtures/query.faa" \
      "$root/fixtures/subject_code${code}.fna" "$work/db_${code}" "$code" "$fmt" \
      > "$out/code${code}_fmt${fmt}.out" 2> "$out/code${code}_fmt${fmt}.stderr"
  done
done
for code in 1 4; do
  "$work/tblastn_code32_oracle" "$root/fixtures/query.faa" \
    "$root/fixtures/subject_code32.fna" "$work/db_32" "$code" 6 \
    > "$out/subject32_g${code}_fmt6.out" 2> "$out/subject32_g${code}_fmt6.stderr"
done
(cd "$out" && sha256sum ./*.out > outputs.sha256)
