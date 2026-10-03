#!/usr/bin/env bash
# Standard output closed at the start (`>&-`): LOSAT and NCBI BLAST+ 2.17.0, five programs,
# outfmt 0/6/7 (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 6, Session S08b).
# NCBI's first write to cout fails: outfmt 0 "BLAST failed to write output", exit 6;
# outfmt 6/7 an uncaught std::ios_base::failure, abort (134). The Rust runtime opens
# /dev/null on a closed standard descriptor before main, so LOSAT writes the report there
# and exits 0. "open-stdout bytes" is the report's size with stdout open (LOSAT and NCBI).
# Usage: closed_stdout.sh LOSAT NCBI_BIN_DIR WORK_DIR
set -u
L=$1; NC=$2; D=$3
F=$(cd "$(dirname "$0")/../../../../LOSAT/tests/fasta" && pwd)
mkdir -p "$D" && cd "$D" || exit 1
awk '/^>/{n++} n<=2' "$F/SicyWSV.faa" > p.faa
head -c 3000 "$F/LC738884.fasta" > n.fna
tail -n +2 "$F/LC738884.fasta" | head -c 12000 | sed '1i>s1' > s.fna
run() {
  prog=$1; shift
  for f in 0 6 7; do for who in LOSAT NCBI; do
    if [ $who = LOSAT ]; then cmd="$L $prog"; else cmd="$NC/$prog"; fi
    sh -c "exec $cmd $* -outfmt $f >&-" 2>err.txt; ec=$?
    bytes=$(sh -c "exec $cmd $* -outfmt $f" 2>/dev/null | wc -c)
    printf '%s\toutfmt %s\t%s\texit %s\topen-stdout bytes %s\tstderr: %s\n' \
      "$prog" "$f" "$who" "$ec" "$bytes" "$(head -c 70 err.txt | tr '\n' '|')"
  done; done
}
{
  run blastn "-query n.fna -subject s.fna"
  run tblastx "-query n.fna -subject s.fna"
  run blastp "-query p.faa -subject p.faa"
  run tblastn "-query p.faa -subject s.fna"
  run blastx "-query n.fna -subject p.faa"
} 2>/dev/null
