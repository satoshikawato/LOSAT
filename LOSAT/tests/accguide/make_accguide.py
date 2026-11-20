#!/usr/bin/env python3
"""Writes LOSAT's copy of NCBI's built-in accession guide (`kBuiltInGuide`).

Input: NCBI BLAST+ 2.17.0 source (598d8ae6), c++/src/objects/seqloc/accguide2.inc (the C array
`kBuiltInGuide` that `SAccGuide::SAccGuide` feeds line by line to `AddRule`, Seq_id.cpp:1462-1471).
Output, under LOSAT/src/seqid/accguide/ (read by seqid::accguide at first use):

  rules.txt   every line of the guide that AddRule does not ignore, its comment (from '#') removed
              and its tokens joined by one space, except the ':' continuation lines of `special2`
              ranges: each run of them is one line ':'.
  ranges.bin  the tokens of those continuation lines, packed two per byte, high nibble first:
              0-9 digit, 10 '-', 11 next token, 12 end of line, 13 end of run, 15 padding.

AddRule strips the comment and splits at blanks with merged delimiters (Seq_id.cpp:1037-1042), so
the dropped bytes do not change the guide; only `count` (its line counter, used in diagnostics only)
differs. The script prints the FNV-1a 64 hash of the compacted lines (continuation lines included,
as ': <tokens>'), which seqid::accguide's test recomputes from the two files.

Usage: python3 -I make_accguide.py <ncbi c++ dir> <LOSAT/src/seqid/accguide>
"""
import os
import re
import sys


def fnv1a64(data):
    h = 0xCBF29CE484222325
    for b in data:
        h ^= b
        h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def main():
    src_dir, out_dir = sys.argv[1], sys.argv[2]
    src = open(os.path.join(src_dir, "src/objects/seqloc/accguide2.inc"), encoding="ascii").read()
    start = src.index("kBuiltInGuide[] = {")
    body = src[start:src.index("};", start)]
    literals = re.findall(r'"((?:[^"\\]|\\.)*)"', body)
    if any(re.search(r'\\[^"]', lit) for lit in literals):
        sys.exit("an escape other than \\\" in a string literal: not expected")
    literals = [lit.replace('\\"', '"') for lit in literals]
    compact = []
    for line in literals:
        tokens = line.split("#", 1)[0].split()
        if tokens:
            compact.append(tokens)
    rules, nibbles = [], []
    in_run = False
    for tokens in compact:
        if tokens[0] == ":":
            if not all(re.fullmatch(r"\d+(-\d+)?", t) for t in tokens[1:]):
                sys.exit("unexpected continuation token: %r" % tokens)
            if not in_run:
                rules.append(":")
                in_run = True
            for i, token in enumerate(tokens[1:]):
                if i:
                    nibbles.append(11)
                nibbles.extend(10 if c == "-" else int(c) for c in token)
            nibbles.append(12)
        else:
            if any(c in t for t in tokens for c in "\\["):
                sys.exit("unexpected mask character: %r" % tokens)
            if in_run:
                nibbles.append(13)
                in_run = False
            rules.append(" ".join(tokens))
    if in_run:
        nibbles.append(13)
    if len(nibbles) % 2:
        nibbles.append(15)
    packed = bytes((nibbles[i] << 4) | nibbles[i + 1] for i in range(0, len(nibbles), 2))
    with open(os.path.join(out_dir, "rules.txt"), "w", encoding="ascii", newline="\n") as f:
        f.write("\n".join(rules) + "\n")
    with open(os.path.join(out_dir, "ranges.bin"), "wb") as f:
        f.write(packed)
    canonical = "".join(" ".join(t) + "\n" for t in compact).encode("ascii")
    print("lines", len(compact), "rules.txt lines", len(rules), "ranges.bin bytes", len(packed))
    print("tokens", sum(len(t) - 1 for t in compact if t[0] == ":"))
    print("canonical bytes", len(canonical), "fnv1a64 0x%016x" % fnv1a64(canonical))


main()
