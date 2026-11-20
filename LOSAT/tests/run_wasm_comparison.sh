#!/bin/bash
# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
# Reuse exactly the native comparison's inputs/options for command-Wasm.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export RUN_NATIVE=0 RUN_NCBI=0
exec bash "$SCRIPT_DIR/run_comparison.sh" "$@"
