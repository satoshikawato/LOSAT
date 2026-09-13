#!/bin/bash
# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
# const string kArgQuery("query");
# const string kArgSubject("subject");
# Reuse exactly the native comparison's inputs/options for command-Wasm.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export RUN_NATIVE=0 RUN_NCBI=0
exec bash "$SCRIPT_DIR/run_comparison.sh" "$@"
