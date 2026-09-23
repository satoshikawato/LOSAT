#!/usr/bin/env bash
set -euo pipefail

# Comparison-only NCBI oracle. Run from any directory with a new output path.
# NCBI c++/src/algo/blast/blastinput/blast_args.cpp:997-1056:
#   static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33};
#   opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
# NCBI c++/src/algo/blast/blastinput/blast_args.cpp:2538-2557:
#   m_Scope = ReadSequencesToBlast(... subjects ...);
#   m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
root=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: run_cli_oracles.sh NEW_OUTPUT_DIRECTORY}
mkdir "$out"
out=$(cd "$out" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
tblastn_bin=${TBLASTN_BIN:-$(command -v tblastn)}
makeblastdb_bin=${MAKEBLASTDB_BIN:-$(command -v makeblastdb)}
[[ $(sha256sum "$tblastn_bin" | cut -d' ' -f1) == e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0 ]]
[[ $(sha256sum "$makeblastdb_bin" | cut -d' ' -f1) == 407bf59ec36ee6444104a90f28eab282e374a0dfe4c7ffee89a65d82dc51c4bf ]]
"$tblastn_bin" -version > "$out/tblastn.version"
"$tblastn_bin" -help > "$out/tblastn.help"
"$makeblastdb_bin" -version > "$out/makeblastdb.version"
sha256sum "$tblastn_bin" "$makeblastdb_bin" "$root"/fixtures/* > "$out/inputs.sha256"
printf '%s\n' 'task=tblastn evalue=10 matrix=BLOSUM62 word_size=3 threshold=13 window_size=40 gapopen=11 gapextend=1 comp_based_stats=2 seg="12 2.2 2.5" sum_stats=true max_intron_length=0 num_threads=1; outfmt=0,6,7' > "$out/options.txt"

for subject in code1 code4 code32; do
  subject_file="$root/fixtures/subject_${subject}.fna"
  "$makeblastdb_bin" -in "$subject_file" -dbtype nucl -out "$work/$subject" > "$out/makeblastdb_${subject}.log"
  for code in 1 4; do
    for mode in subject db; do
      for fmt in 0 6 7; do
        name="${subject}_g${code}_${mode}_fmt${fmt}"
        if [[ $mode == subject ]]; then
          source_args=(-subject "$subject_file")
        else
          source_args=(-db "$work/$subject")
        fi
        "$tblastn_bin" -task tblastn -query "$root/fixtures/query.faa" "${source_args[@]}" \
          -evalue 10 -matrix BLOSUM62 -word_size 3 -threshold 13 -window_size 40 \
          -gapopen 11 -gapextend 1 -comp_based_stats 2 -seg "12 2.2 2.5" -sum_stats true \
          -max_intron_length 0 -db_gencode "$code" -num_threads 1 -outfmt "$fmt" \
          -out "$out/$name.out" 2> "$out/$name.stderr"
      done
      # NCBI c++/src/objtools/align_format/tabular.cpp:907-1095:
      #   score and sframe expose raw-score and subject-frame differences.
      "$tblastn_bin" -task tblastn -query "$root/fixtures/query.faa" "${source_args[@]}" \
        -evalue 10 -matrix BLOSUM62 -word_size 3 -threshold 13 -window_size 40 \
        -gapopen 11 -gapextend 1 -comp_based_stats 2 -seg "12 2.2 2.5" -sum_stats true \
        -max_intron_length 0 -db_gencode "$code" -num_threads 1 \
        -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score sframe qseq sseq' \
        -out "$out/${subject}_g${code}_${mode}_detail.out" \
        2> "$out/${subject}_g${code}_${mode}_detail.stderr"
    done
  done
done
"$tblastn_bin" -task tblastn -query "$root/fixtures/query.faa" \
  -subject "$root/fixtures/subject_code1.fna" -outfmt 6 -num_threads 1 \
  > "$out/code1_default_subject_fmt6.out"
cmp "$out/code1_default_subject_fmt6.out" "$out/code1_g1_subject_fmt6.out"
set +e
"$tblastn_bin" -task tblastn -query "$root/fixtures/query.faa" \
  -subject "$root/fixtures/subject_code32.fna" -db_gencode 32 -outfmt 6 \
  > "$out/code32_cli.stdout" 2> "$out/code32_cli.stderr"
status=$?
printf '%s\n' "$status" > "$out/code32_cli.status"
"$tblastn_bin" -task tblastn -query "$root/fixtures/query.faa" \
  -db "$work/code32" -db_gencode 32 -outfmt 6 \
  > "$out/code32_db_cli.stdout" 2> "$out/code32_db_cli.stderr"
db_status=$?
set -e
printf '%s\n' "$db_status" > "$out/code32_db_cli.status"
if [[ $status -eq 0 || $db_status -eq 0 ]]; then
  echo 'Unexpected acceptance of genetic code 32 by pinned NCBI CLI' >&2
  exit 1
fi
(cd "$out" && sha256sum ./*.out ./*.stderr ./*.stdout > outputs.sha256)
