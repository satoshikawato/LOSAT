#!/usr/bin/env bash
# Line WO: mutants of the bulk guard (memory.fill, memory.copy) against a test file. Each must fail.
# usage: mutants_bulk.sh <worktree> <test file> <node>...
set -u
WT=$(realpath "$1"); TESTS=$(realpath "$2"); shift 2; NODES="$*"
TMP=$(mktemp -d "${TMPDIR:-$HOME}/wo-mutants.XXXXXX")
cp "$TESTS" "$TMP/test_wasi_shared_memory.js"
mut() { # name, old text, new text
  cp "$WT/LOSAT/tests/wasi_shared_memory.js" "$TMP/wasi_shared_memory.js"
  python3 - "$TMP/wasi_shared_memory.js" "$2" "$3" <<'P' || { echo "$1: PATCH FAILED"; return; }
import sys
p, a, b = sys.argv[1:4]; s = open(p).read()
assert s.count(a) == 1, a
open(p, 'w').write(s.replace(a, b))
P
  for node in $NODES; do
    timeout 280 "$node" --test "$TMP/test_wasi_shared_memory.js" > "$TMP/log" 2>&1
    printf '%s %s: %s | %s\n' "$1" "$("$node" -v)" "$(grep -E 'ℹ (pass|fail) ' "$TMP/log" | tr -s ' ' | tr '\n' ' ')" "$(grep -E '^✖ ' "$TMP/log" | grep -v 'failing tests' | sed 's/ (.*//; s/^✖ //' | sort -u | cut -c1-70 | tr '\n' ';' | cut -c1-300)"
  done
}
mut none "const GUARD_VERSION = 3;" "const GUARD_VERSION = 3;"
mut no-bulk-compare "const outside = (index) => [0x20, index, 0xad, 0x20, 2, 0xad, 0x7c, 0x20, 3, 0x56];" "const outside = (index) => [0x41, 0];"
mut no-source-compare "...outside(0), ...(copy ? [...outside(1), 0x72] : [])," "...outside(0),"
mut start-only "const outside = (index) => [0x20, index, 0xad, 0x20, 2, 0xad, 0x7c, 0x20, 3, 0x56];" "const outside = (index) => [0x20, index, 0xad, 0x20, 3, 0x56];"
rm -rf "$TMP"
