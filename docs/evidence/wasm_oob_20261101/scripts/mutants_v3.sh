#!/usr/bin/env bash
# Line WO: guard v3 mutants. Each must fail the guard's tests. usage: mutants_v3.sh <worktree> <node>...
set -u
WT=$1; shift
TMP=$(mktemp -d "${TMPDIR:-$HOME}/wo-mutants.XXXXXX")
cp "$WT/LOSAT/tests/test_wasi_shared_memory.js" "$TMP/"
mut() { # name, python expression on s
  cp "$WT/LOSAT/tests/wasi_shared_memory.js" "$TMP/wasi_shared_memory.js"
  python3 - "$TMP/wasi_shared_memory.js" "$2" "$3" <<'P' || { echo "$1: PATCH FAILED"; return; }
import sys
p,a,b=sys.argv[1:4]; s=open(p).read()
assert a in s, a
open(p,'w').write(s.replace(a,b,1))
P
  for n in "$@"; do :; done
  for node in $NODES; do
    timeout 280 "$node" --test "$TMP/test_wasi_shared_memory.js" > "$TMP/log" 2>&1
    printf '%s %s: %s | %s\n' "$1" "$("$node" -v)" "$(grep -E 'ℹ (pass|fail) ' "$TMP/log" | tr -s ' ' | tr '\n' ' ')" "$(grep -E '^✖ ' "$TMP/log" | grep -v 'failing tests' | sed 's/ (.*//; s/^✖ //' | sort -u | cut -c1-60 | tr '\n' ';' | cut -c1-260)"
  done
}
NODES="$*"
mut none "const GUARD_VERSION = 3;" "const GUARD_VERSION = 3;"
mut no-refresh-after "const atomicAcquires = (sub) => sub !== 0 && !(sub >= 0x17 && sub <= 0x1d);" "const atomicAcquires = (sub) => false;"
mut grow-no-raise "0x41, 0x7f, 0x47, // (old = memory.grow(pages)) != -1" "0x41, 0x7f, 0x46, // MUTANT: only when it failed"
mut no-bound-check "0x20, 0, 0xad, 0x42, ...signedLeb(offset + size), 0x7c, 0x20, bytes, 0x56, // address + offset + size > bytes" "0x41, 0, 0xad, 0x41, 1, 0xad, 0x56, // MUTANT: never"
mut raise-any "0x3f, 0, 0x49, // seen < memory.size" "0x3f, 0, 0x47, // MUTANT: seen != memory.size"
rm -rf "$TMP"
