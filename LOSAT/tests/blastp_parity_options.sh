#!/bin/bash

# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
#
# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:101-211
#
# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:109
#
# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:208
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

# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
BLASTP_PARITY_STRICT_OPTIONS=(
    "${BLASTP_PARITY_COMMON_OPTIONS[@]}"
    -evalue 1e-5
    -max_target_seqs 10
)

# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
BLASTP_PARITY_DEFAULT_OPTIONS=(
    "${BLASTP_PARITY_COMMON_OPTIONS[@]}"
    -evalue 10
    -max_target_seqs 500
)
