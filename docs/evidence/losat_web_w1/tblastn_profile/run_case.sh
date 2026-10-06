#!/bin/bash
# usage: run_case.sh <query> <subject> -> prints stdout path
BIN=/home/kawato/.cache/losat-web-gui-target/app-s09-native/release/LOSAT
exec $BIN tblastn -task tblastn -query "$1" -subject "$2" -outfmt 6 -num_threads 1
