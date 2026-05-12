#!/bin/bash

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kTask("task");
# const string kArgMatrixName("matrix");
# const string kArgEvalue("evalue");
# const string kArgMaxTargetSequences("max_target_seqs");
# const string kArgGapOpen("gapopen");
# const string kArgGapExtend("gapextend");
# ```
#
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:101-211
# ```c
# const string kArgWindowSize("window_size");
# const string kArgWordSize("word_size");
# const string kArgCompBasedStats("comp_based_stats");
# const string kArgSegFiltering("seg");
# ```
#
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:570-622
# ```c
# const string kArgWordThreshold("threshold");
# ```
#
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:190-210
# ```c
# const string kArgMaxHSPsPerSubject("max_hsps");
# ```
BLASTP_PARITY_COMMON_OPTIONS=(
    -task blastp
    -matrix BLOSUM62
    -gapopen 11
    -gapextend 1
    -word_size 3
    -threshold 11
    -window_size 40
    -comp_based_stats 2
    -seg no
    -max_hsps 0
)

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgEvalue("evalue");
# const string kArgMaxTargetSequences("max_target_seqs");
# ```
BLASTP_PARITY_STRICT_OPTIONS=(
    "${BLASTP_PARITY_COMMON_OPTIONS[@]}"
    -evalue 1e-5
    -max_target_seqs 10
)

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgEvalue("evalue");
# const string kArgMaxTargetSequences("max_target_seqs");
# ```
BLASTP_PARITY_DEFAULT_OPTIONS=(
    "${BLASTP_PARITY_COMMON_OPTIONS[@]}"
    -evalue 10
    -max_target_seqs 500
)
