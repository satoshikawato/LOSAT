#!/bin/bash
# Builds the fixtures of the rows db.swgap.* (line RG item X4b-gap, 2026-11-24) with NCBI BLAST+ 2.17.0's
# makeblastdb (comparison-only tooling, as the fixtures' NCBI outputs; LOSAT never runs NCBI):
#   tests/fasta/swgap/f101.faa   residues 96..285 of tests/fasta/rps/nest.faa (WSSV fragment 101, which nest
#                                holds twice with its parts interleaved), 60 a line
#   tests/fasta/swgap/f101.fna   bases 286..855 of tests/fasta/rps/nest_p.fna (the same residues' codons)
#   tests/fasta/swgap/f101m.fna  its reverse complement
#   nestp, nestm                 version 4 nucleotide databases of tests/fasta/rps/nest_p.fna and nest_m.fna
#   nestaa                       a version 4 protein database of tests/fasta/rps/nest.faa
# With -use_sw_tback the Smith-Waterman alignment of f101 against nest starts with 95 gaps in the query's row
# (its edit script starts with a deletion): BLASTP, BLASTX and TBLASTN report it as NCBI does
# (`report::serial::seqalign::row_end_gaps`, `mapped_script`).
# makeblastdb writes the build date into the .nin/.pin, so a rebuild differs from the committed files: the
# committed files are the fixtures, this script records how they were made.
#
# Usage: build.sh NCBI_BIN
set -euo pipefail
N=$1
HERE=$(cd "$(dirname "$0")" && pwd)
LOSAT=$(cd "$HERE/../../.." && pwd)
RPS=$LOSAT/tests/fasta/rps
Q=$LOSAT/tests/fasta/swgap
mkdir -p "$Q"
python3 -I - "$RPS" "$Q" <<'PY'
import sys
rps, q = sys.argv[1], sys.argv[2]
def seq(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">"))
def write(path, title, letters):
    with open(path, "w") as out:
        out.write(f">{title}\n")
        for at in range(0, len(letters), 60):
            out.write(letters[at:at + 60] + "\n")
bases = seq(f"{rps}/nest_p.fna")[285:855]
write(f"{q}/f101.faa", "f101", seq(f"{rps}/nest.faa")[95:285])
write(f"{q}/f101.fna", "f101n", bases)
write(f"{q}/f101m.fna", "f101m", bases[::-1].translate(str.maketrans("ACGT", "TGCA")))
PY
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
cp "$RPS/nest_p.fna" "$RPS/nest_m.fna" "$RPS/nest.faa" "$W/"
cd "$W"
"$N/makeblastdb" -in nest_p.fna -dbtype nucl -out "$HERE/nestp" -blastdb_version 4
"$N/makeblastdb" -in nest_m.fna -dbtype nucl -out "$HERE/nestm" -blastdb_version 4
"$N/makeblastdb" -in nest.faa -dbtype prot -out "$HERE/nestaa" -blastdb_version 4
