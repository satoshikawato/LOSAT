#!/usr/bin/env python3
"""The profile database rps_wide and its variant rps_wide_ug (line PS, series PR, item PS12m).

No NCBI counterpart: test inputs for build.sh. Two subcommands:

  build_wide.py tile F101.smp K OUT.smp
      Tiles the PSSM of f101 (build.sh: psiblast PSSM with the scores dropped) K times: a profile of
      K x 190 columns, its consensus letters and frequency ratios repeated. build.sh makes rps_wide of
      K = 700 (133,000 columns: 131,072 rows are 64 buckets of 2,048, the most that one group of
      buckets of LOSAT's scan holds, so the volume has 66 buckets and two groups) with makeprofiledb
      -threshold 20 -force: only the words that match exactly are in the table (9 words, each with
      700 hits, one in every tile, so a cell's hits lie in 65 buckets and both groups).

  build_wide.py ungroup IN.loo OUT.loo
      Writes the table with the hits of its cells not in ascending order (never written by
      makeprofiledb, which appends the hits of the profile's rows in order): the lists of the nine
      cells reversed (every hit of the later group before every hit of the earlier one), rotated by
      a third, and with the first and the last hit swapped, in turn; the cell's first hit stays in
      the cell, the others in the overflow run (the layout is unchanged).
"""

import re
import struct
import sys

CELLS = 32768
HEADER = 40


def tile(src: str, k: int, out: str) -> None:
    text = open(src).read()
    # The final data (scores) dropped, the query id CDD:901, as build.sh does for the fragments.
    start = text.index(",\n    finalData {")
    end = text.index("\n    }", start) + len("\n    }")
    text = text[:start] + text[end:]
    text = re.sub(
        r'id \{\s*local str "[^"]*"\s*\}',
        'id {\n        general {\n          db "CDD",\n          tag id 901\n        }\n      }',
        text,
        count=1,
    )
    match = re.search(r"seq-data ncbistdaa '([0-9A-F\s]+)'H", text)
    letters = re.sub(r"\s+", "", match.group(1))
    columns = len(letters) // 2
    letters *= k
    lines = "\n".join(letters[i : i + 70] for i in range(0, len(letters), 70))
    text = text[: match.start(1)] + lines + text[match.end(1) :]
    text = text.replace(f"length {columns},", f"length {columns * k},", 1)
    text = text.replace(f"numColumns {columns},", f"numColumns {columns * k},", 1)
    first = text.index("freqRatios {") + len("freqRatios {")
    last = text.index("\n      }", first)
    rows = [row.strip().rstrip(",") for row in text[first:last].strip("\n").split("\n")]
    assert len(rows) == 28 * columns, (len(rows), columns)
    text = text[:first] + "\n" + ",\n".join("        " + row for row in rows * k) + text[last:]
    open(out, "w").write(text)


def ungroup(src: str, out: str) -> None:
    data = bytearray(open(src, "rb").read())
    assert struct.unpack_from("<i", data, 32)[0] == HEADER
    cells_at = HEADER
    overflow_at = HEADER + 16 * (CELLS + 1)
    changed = 0
    for index in range(CELLS):
        at = cells_at + 16 * index
        num_used, first, offset, _ = struct.unpack_from("<4i", data, at)
        if num_used <= 3:
            continue
        run_at = overflow_at + offset
        run = list(struct.unpack_from(f"<{num_used - 1}i", data, run_at))
        hits = [first] + run
        kind = changed % 3
        if kind == 0:
            hits.reverse()
        elif kind == 1:
            cut = len(hits) // 3
            hits = hits[cut:] + hits[:cut]
        else:
            hits[0], hits[-1] = hits[-1], hits[0]
        struct.pack_into("<i", data, at + 4, hits[0])
        struct.pack_into(f"<{num_used - 1}i", data, run_at, *hits[1:])
        changed += 1
    assert changed > 0
    open(out, "wb").write(bytes(data))
    print(changed, "cells permuted")


if __name__ == "__main__":
    if sys.argv[1] == "tile" and len(sys.argv) == 5:
        tile(sys.argv[2], int(sys.argv[3]), sys.argv[4])
    elif sys.argv[1] == "ungroup" and len(sys.argv) == 4:
        ungroup(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
