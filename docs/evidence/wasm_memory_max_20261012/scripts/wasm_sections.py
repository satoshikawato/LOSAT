#!/usr/bin/env python3
"""Compare two Wasm modules section by section.
Usage: wasm_sections.py A.wasm B.wasm   (prints each section's id, size, SHA-256 prefix and whether it is equal)
"""
import hashlib
import sys

NAMES = {0: "custom", 1: "type", 2: "import", 3: "function", 4: "table", 5: "memory", 6: "global", 7: "export",
         8: "start", 9: "element", 10: "code", 11: "data", 12: "datacount"}


def uleb(data, pos):
    value = shift = 0
    while True:
        byte = data[pos]
        pos += 1
        value |= (byte & 127) << shift
        shift += 7
        if not byte & 128:
            return value, pos


def sections(path):
    data = open(path, "rb").read()
    pos, out = 8, []
    while pos < len(data):
        sid = data[pos]
        size, start = uleb(data, pos + 1)
        body = data[start:start + size]
        name = NAMES.get(sid, str(sid))
        if sid == 0:
            n, p = uleb(body, 0)
            name = "custom:" + body[p:p + n].decode()
        out.append((name, body))
        pos = start + size
    return out


def memory_import(body):
    count, pos = uleb(body, 0)
    for _ in range(count):
        for _ in range(2):
            n, pos = uleb(body, pos)
            pos += n
        kind = body[pos]
        pos += 1
        if kind == 0:
            _, pos = uleb(body, pos)
        elif kind == 1:
            pos += 1
            flags, pos = uleb(body, pos)
            _, pos = uleb(body, pos)
            if flags & 1:
                _, pos = uleb(body, pos)
        elif kind == 2:
            flags, pos = uleb(body, pos)
            initial, pos = uleb(body, pos)
            maximum = None
            if flags & 1:
                maximum, pos = uleb(body, pos)
            return flags, initial, maximum
        elif kind == 3:
            pos += 2
    return None


a, b = sections(sys.argv[1]), sections(sys.argv[2])
assert [n for n, _ in a] == [n for n, _ in b], "different section lists"
differing = []
for (name, x), (_, y) in zip(a, b):
    equal = x == y
    if not equal:
        differing.append(name)
    print(f"{name:28s} {len(x):9d} {len(y):9d} {hashlib.sha256(x).hexdigest()[:12]} {hashlib.sha256(y).hexdigest()[:12]} {'equal' if equal else 'DIFFERENT'}")
    if name == "import":
        print("  memory import (flags, initial, maximum):", memory_import(x), memory_import(y))
        if not equal:
            diff = [i for i in range(min(len(x), len(y))) if x[i] != y[i]]
            print(f"  differing byte offsets in the section: {diff[:16]} (count {len(diff)}, lengths {len(x)} {len(y)})")
print("differing sections:", differing or "none")
