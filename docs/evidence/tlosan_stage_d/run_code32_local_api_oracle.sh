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
# NCBI core/blast_engine.c:870-899 links and reaps before
# core/blast_kappa.c:3577-3585 converts the retained preliminary HSPs.
gcc -std=c11 -O0 -shared -fPIC "$root/ncbi_d_call_trace.c" -ldl -o "$work/d_trace.so"
# NCBI core/blast_engine.c:1407,1434-1443 reads TotLen and conditionally
# updates per-subject search parameters.
gcc -std=c11 -O0 -shared -fPIC "$root/ncbi_seqsrc_callstate_trace.c" -ldl -o "$work/seqsrc_trace.so"
# NCBI core/blast_parameters.c:902-999 and core/blast_gapalign.c:3924-3927
# expose the actual setup and per-context search cutoffs.
gcc -std=c11 -O0 -shared -fPIC "$root/ncbi_parameter_trace.c" -ldl -o "$work/parameter_trace.so"
gcc -std=c11 -O0 -shared -fPIC "$root/../tlosan_stage_c/ncbi_context_cutoff_trace.c" -ldl -o "$work/hit_cutoff_trace.so"
gcc -std=c11 -O0 -shared -fPIC "$root/../tlosan_stage_c/ncbi_wordfinder_context_cutoff_trace.c" -ldl -o "$work/word_cutoff_trace.so"
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
  # NCBI core/blast_engine.c:870-899; core/blast_kappa.c:3577-3585:
  # retain the complete preliminary link/reap input and output call order.
  LD_PRELOAD="$work/d_trace.so" "$work/tblastn_code32_local_oracle" \
    "$stage_a/fixtures/query.faa" "$stage_a/fixtures/subject_code32.fna" "$code" 6 \
    > "$out/code${code}_d.out" 2> "$out/code${code}_d.stderr"
  python3 - "$out" "$code" <<'PY_D_TRACE'
from pathlib import Path
import sys
out, code = Path(sys.argv[1]), sys.argv[2]
traced_out = (out / f"code{code}_d.out").read_bytes()
traced_err = (out / f"code{code}_d.stderr").read_bytes()
selected = b"".join(line for line in traced_err.splitlines(keepends=True) if line.startswith(b"D_"))
ordinary = b"".join(line for line in traced_err.splitlines(keepends=True) if not line.startswith(b"D_"))
assert traced_out == (out / f"code{code}_fmt6.out").read_bytes()
assert ordinary == (out / f"code{code}.stderr").read_bytes()
assert b"\tlink_before\t" in selected and b"\treap_after\t" in selected
(out / f"code{code}_d.trace").write_bytes(selected)
PY_D_TRACE
  # NCBI core/blast_engine.c:1407,1434-1443: capture TotLen and the
  # conditional OneSubjectUpdateParameters call in the API local path.
  LD_PRELOAD="$work/seqsrc_trace.so" "$work/tblastn_code32_local_oracle" \
    "$stage_a/fixtures/query.faa" "$stage_a/fixtures/subject_code32.fna" "$code" 6 \
    > "$out/code${code}_seqsrc.out" 2> "$out/code${code}_seqsrc.stderr"
  python3 - "$out" "$code" <<'PY_SEQSRC_TRACE'
from pathlib import Path
import sys
out, code = Path(sys.argv[1]), sys.argv[2]
traced = (out / f"code{code}_seqsrc.stderr").read_bytes()
selected = b"".join(line for line in traced.splitlines(keepends=True) if line.startswith((b"D_SEQSRC_", b"D_ONE_SUBJECT_")))
ordinary = b"".join(line for line in traced.splitlines(keepends=True) if not line.startswith((b"D_SEQSRC_", b"D_ONE_SUBJECT_")))
assert (out / f"code{code}_seqsrc.out").read_bytes() == (out / f"code{code}_fmt6.out").read_bytes()
assert ordinary == (out / f"code{code}.stderr").read_bytes()
assert b"D_SEQSRC_TOTLEN\t" in selected
(out / f"code{code}_seqsrc.trace").write_bytes(selected)
PY_SEQSRC_TRACE
  # NCBI core/blast_setup.c:1011-1024, core/blast_gapalign.c:3924-3927,
  # core/aa_ungapped.c:547-583: preserve initial and call-site cutoffs.
  for cutoff_probe in parameter hit_cutoff word_cutoff; do
    LD_PRELOAD="$work/${cutoff_probe}_trace.so" "$work/tblastn_code32_local_oracle" \
      "$stage_a/fixtures/query.faa" "$stage_a/fixtures/subject_code32.fna" "$code" 6 \
      > "$out/code${code}_${cutoff_probe}.out" 2> "$out/code${code}_${cutoff_probe}.stderr"
    python3 - "$out" "$code" "$cutoff_probe" <<'PY_CUTOFF_TRACE'
from pathlib import Path
import sys
out, code, probe = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
prefix = {"parameter": b"D_PARAM_", "hit_cutoff": b"GAPPED_CONTEXT_CUTOFF\t", "word_cutoff": b"WORD_CONTEXT_CUTOFF\t"}[probe]
assert (out / f"code{code}_{probe}.out").read_bytes() == (out / f"code{code}_fmt6.out").read_bytes()
selected = [line for line in (out / f"code{code}_{probe}.stderr").read_bytes().splitlines(keepends=True) if line.startswith(prefix)]
assert selected
(out / f"code{code}_{probe}.trace").write_bytes(b"".join(selected))
PY_CUTOFF_TRACE
  done
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

# NCBI core/blast_engine.c:870-899,1407,1434-1443: code-1 API and CLI
# must have the same preliminary D calls and sequence-source update state.
for probe in d seqsrc; do
  LD_PRELOAD="$work/${probe}_trace.so" "$tblastn_bin" -task tblastn \
    -query "$stage_a/fixtures/query.faa" \
    -subject "$stage_a/fixtures/subject_code32.fna" -db_gencode 1 -num_threads 1 \
    -outfmt 6 > "$out/cli_code1_${probe}.out" 2> "$out/cli_code1_${probe}.stderr"
  python3 - "$out" "$probe" <<'PY_CLI_D_TRACE'
from pathlib import Path
import sys
out, probe = Path(sys.argv[1]), sys.argv[2]
prefixes = (b"D_",) if probe == "d" else (b"D_SEQSRC_", b"D_ONE_SUBJECT_")
traced = (out / f"cli_code1_{probe}.stderr").read_bytes()
selected = b"".join(line for line in traced.splitlines(keepends=True) if line.startswith(prefixes))
ordinary = b"".join(line for line in traced.splitlines(keepends=True) if not line.startswith(prefixes))
assert (out / f"cli_code1_{probe}.out").read_bytes() == (out / "cli_code1_fmt6.out").read_bytes()
assert ordinary == (out / "cli_code1.stderr").read_bytes()
assert selected == (out / f"code1_{probe}.trace").read_bytes()
(out / f"cli_code1_{probe}.trace").write_bytes(selected)
PY_CLI_D_TRACE
done
# NCBI core/blast_engine.c:1407,1434-1443; core/blast_kappa.c:3577-3595:
# this pinned fixture must use the CLI local-subject call state and preserve
# both code-32 preliminary HSPs through Kappa entry.
python3 - "$out" <<'PY_LOCAL_STATE'
from pathlib import Path
import sys
out = Path(sys.argv[1])
for code in (1, 32):
    lines = (out / f"code{code}_seqsrc.trace").read_text().splitlines()
    assert lines and all(line == "D_SEQSRC_TOTLEN\t360" for line in lines)
code32_d = (out / "code32_d.trace").read_text().splitlines()
prelim = [line.split("\t") for line in code32_d if line.startswith("D_HSP\t0\tlink_before\t")]
assert len(prelim) == 2
assert [int(row[10]) for row in prelim] == [656, 16]
code32_kappa = (out / "code32.trace").read_text().splitlines()
assert len([line for line in code32_kappa if line.startswith("K_TRACE_PRELIM\t")]) == 2
assert (out / "code32_hit_cutoff.trace").read_text().splitlines()[0].endswith("\t9\t9")
assert all(line.endswith("\t9") for line in (out / "code32_word_cutoff.trace").read_text().splitlines())
PY_LOCAL_STATE
printf 'NCBI source commit: %s\nAPI mode: local -subject dbscan_mode=true; explicit FindGeneticCode on BlastSeqSrc retrieval.\nCode 1 API outfmt 6, Kappa, Stage D, and TotLen call traces match the local -subject CLI bytes.\nCode 32 TotLen=360 skips OneSubjectUpdateParameters; initial word/gapped cutoff=9; two Kappa preliminary HSPs score 656 and 16.\n' "$source_commit" > "$out/manifest.txt"
sha256sum "$root/tblastn_code32_local_oracle.cpp" "$root/ncbi_kappa_traceback_trace.c" \
  "$root/ncbi_kappa_mode2_trace.c" "$root/ncbi_kappa_composition_trace.c" \
  "$root/ncbi_d_call_trace.c" "$root/ncbi_parameter_trace.c" \
  "$root/ncbi_seqsrc_callstate_trace.c" \
  "$root/../tlosan_stage_c/ncbi_context_cutoff_trace.c" \
  "$root/../tlosan_stage_c/ncbi_wordfinder_context_cutoff_trace.c" \
  "$stage_a/fixtures/query.faa" "$stage_a/fixtures/subject_code32.fna" \
  "$work/tblastn_code32_local_oracle" "$tblastn_bin" > "$out/inputs.sha256"
(cd "$out" && sha256sum ./*.out ./*.trace > outputs.sha256)
