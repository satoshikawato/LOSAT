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

Every range of a `special2` block must be <digits>[-<digits>] with at most as many digits a side
as x_SplitSpecial reads for the bit (Seq_id.cpp:1599-1605), and each block must follow its
`special2` line: seqid::accguide reads that form only, and the script refuses any other.

The script also loads the guide as SAccGuide's constructor does (AddRule, x_AddSpecial,
x_SplitSpecial and the cleanup, Seq_id.cpp:1030-1053,1054-1131,1133-1247,1292-1323,1443-1498,
1513-1594,1596-1625), with the values of the type names in acc_info.rs of the output directory,
and prints the fingerprint of the loaded guide that seqid::accguide's test compares
(GUIDE_FINGERPRINT): FNV-1a 64 of every format, prefix, wildcard, big special and small special
with its options and ranges, each list after its length.

Usage: python3 -I make_accguide.py <ncbi c++ dir> <LOSAT/src/seqid/accguide>
"""
import bisect
import os
import re
import struct
import sys


def fnv1a64(data):
    h = 0xCBF29CE484222325
    for b in data:
        h ^= b
        h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


MAX_SMALL_SPECIAL_DIGITS = 9


def bit_digits(digits):
    """The number of last digits that x_SplitSpecial reads as the bit (Seq_id.cpp:1599-1605)."""
    if digits == MAX_SMALL_SPECIAL_DIGITS + 1:
        return digits - 2
    return min(digits, MAX_SMALL_SPECIAL_DIGITS)


def to_uint(text):
    """NStr::StringToUInt with fConvErr_NoThrow: 0 for anything but decimal digits."""
    return int(text) if text.isdigit() and text.isascii() else 0


def acc_info_values(path):
    """The value of each type name of acc_info.rs (ACC_INFO_NAMES), its constants evaluated
    with `|`, `<<` and acc(a, b, c) = a | b | c."""
    text = open(path, encoding="utf-8").read()
    consts = dict(re.findall(r"pub(?:\(super\))? const (\w+): u32 = ([^;]+);", text))
    values = {}

    def evaluate(expr):
        tokens = re.findall(r"0x[0-9a-fA-F_]+|\d[\d_]*|\w+|<<|[|(),]", expr)
        pos = 0

        def atom():
            nonlocal pos
            token = tokens[pos]
            pos += 1
            if token == "acc":
                assert tokens[pos] == "("
                pos += 1
                result = alternatives()
                while tokens[pos] == ",":
                    pos += 1
                    result |= alternatives()
                assert tokens[pos] == ")"
                pos += 1
                return result
            if token[0].isdigit():
                return int(token.replace("_", ""), 0)
            if token not in values:
                values[token] = evaluate(consts[token])
            return values[token]

        def shifted():
            nonlocal pos
            result = atom()
            if pos < len(tokens) and tokens[pos] == "<<":
                pos += 1
                result <<= atom()
            return result

        def alternatives():
            nonlocal pos
            result = shifted()
            while pos < len(tokens) and tokens[pos] == "|":
                pos += 1
                result |= shifted()
            return result

        result = alternatives()
        assert pos == len(tokens), expr
        return result

    table = text[text.index("ACC_INFO_NAMES"):]
    table = table[table.index("= [") + 3:table.index("\n];")]
    names, depth, item = {}, 0, ""
    for c in table:
        if c == "(":
            depth += 1
            if depth == 1:
                item = ""
                continue
        elif c == ")":
            depth -= 1
            if depth == 0:
                name, expr = re.fullmatch(r'\s*"([a-z0-9_]+)",\s*(.*?),?\s*', item, re.S).groups()
                names[name] = evaluate(expr)
                continue
        if depth >= 1:
            item += c
    return names, evaluate("F_SPECIALS"), evaluate("F_FALLBACK")


class Bits:
    """A bm::bvector as sorted, disjoint, non-adjacent closed ranges; set and clear take the two
    ends in either order (BitMagic 8.0.1, NCBI 598d8ae6 copy: c++/include/util/bitset/bm.h:1216-1217,
    2377-2378)."""

    def __init__(self):
        self.starts, self.ends = [], []

    def set(self, left, right):
        left, right = min(left, right), max(left, right)
        starts, ends = self.starts, self.ends
        first = bisect.bisect_left(ends, left - 1)
        last = bisect.bisect_right(starts, right + 1)
        if first < last:
            left, right = min(left, starts[first]), max(right, ends[last - 1])
        starts[first:last] = [left]
        ends[first:last] = [right]

    def clear(self, left, right):
        left, right = min(left, right), max(left, right)
        starts, ends = self.starts, self.ends
        first = bisect.bisect_left(ends, left)
        last = bisect.bisect_right(starts, right)
        if first >= last:
            return
        kept = []
        if starts[first] < left:
            kept.append((starts[first], left - 1))
        if ends[last - 1] > right:
            kept.append((right + 1, ends[last - 1]))
        starts[first:last] = [k[0] for k in kept]
        ends[first:last] = [k[1] for k in kept]


def guide_fingerprint(compact, names, f_specials, f_fallback):
    """The fingerprint of the guide that AddRule makes of the compacted lines."""
    rules, hints = {}, {"submap": None, "small": None, "special2": None}

    def known(name):
        value = names.get(name.lower())
        if value is None:
            sys.exit("unknown type name %r: the fallbacks are not modelled" % name)
        return value

    def find_sub_map(fmt):
        if hints["submap"] != fmt:
            rules.setdefault(fmt, ({}, [], {}, {}))
            hints["submap"] = fmt
            hints["small"] = None

    def split_special(acc, fmt):
        letters, raw = fmt >> 16, fmt & 0xFFFF
        digits = bit_digits(raw)
        if len(acc) == letters + raw:
            pos = len(acc) - digits
            return acc[:pos], to_uint(acc[pos:])
        pos = letters + 2
        result = to_uint(acc[pos + 1:])
        if digits == raw:
            pos -= 2
            result += to_uint(acc[pos:pos + 2]) * 10 ** (digits - 2)
        return acc[:pos], result & 0xFFFFFFFF

    def same(a, b):
        return a & ~f_fallback == b & ~f_fallback

    def add_special(fmt, first, last, value):
        _, _, big, small = rules[fmt]
        key, left = split_special(first, fmt)
        to_key, right = split_special(last, fmt)
        if key != to_key:
            big.setdefault(last, [first, value])[1] = value
            return
        options = small.setdefault(key, [])
        hint = hints["small"]
        if hint is not None and hint[0] == key:
            index = hint[1]
            options[index][0].clear(left, right)
            while not same(options[index][1], value):
                if index == 0:
                    index = hint[1] + 1
                    break
                index -= 1
        else:
            index = 0
        while index < len(options) and not same(options[index][1], value):
            index += 1
        if index == len(options):
            options.append([Bits(), value])
        options[index][0].set(left, right)
        options[index][1] = value
        hints["small"] = (key, index)

    for tokens in compact:
        first = tokens[0]
        if len(tokens) == 2 and first.lower() == "version":
            hints["special2"] = None
        elif "+" in first and (len(tokens) == 3 or (len(tokens) == 4 and tokens[3] == "*")):
            hints["special2"] = None
            letters, digits = first.split("+", 1)
            fmt = (to_uint(letters) & 0xFFFF) << 16 | (to_uint(digits) & 0xFFFF)
            value = known(tokens[2])
            find_sub_map(fmt)
            prefixes, wildcards, _, _ = rules[fmt]
            if len(tokens) == 4:
                value |= f_specials
            prefix = tokens[1]
            if "?" not in prefix and "*" not in prefix:
                prefixes[prefix] = value
            else:
                for wildcard in wildcards:
                    if wildcard[0] == prefix:
                        wildcard[1] = value
                        break
                else:
                    wildcards.append([prefix, value])
        elif len(tokens) == 3 and first.lower() == "special":
            hints["special2"] = None
            accs = tokens[1]
            pos = next((i for i, c in enumerate(accs) if c.isdigit()), len(accs))
            dash = accs.find("-", pos)
            end = len(accs) if dash < 0 else dash
            fmt = (pos & 0xFFFF) << 16 | ((end - pos) & 0xFFFF)
            value = known(tokens[2])
            find_sub_map(fmt)
            add_special(fmt, accs[:end], accs if dash < 0 else accs[dash + 1:], value)
        elif len(tokens) >= 3 and first.lower() == "special2":
            letters, _, digits = tokens[1].partition("+")
            fmt = (len(letters) & 0xFFFF) << 16 | (to_uint(digits) & 0xFFFF)
            find_sub_map(fmt)
            hints["special2"] = [fmt, list(letters + "0" * to_uint(digits)), known(tokens[2])]
        elif len(tokens) >= 2 and first == ":":
            fmt, acc, value = hints["special2"]
            for token in tokens[1:]:
                low, _, high = token.partition("-")
                acc[len(acc) - len(low):] = low
                start = "".join(acc)
                if high:
                    acc[len(acc) - len(high):] = high
                add_special(fmt, start, "".join(acc), value)
        elif len(tokens) == 3 and first.lower() in ("gnl", "fallback"):
            hints["special2"] = None
        else:
            sys.exit("a line that AddRule ignores with a warning: %r" % tokens)

    out = bytearray()

    def put(data):
        out.extend(struct.pack("<I", len(data)))
        out.extend(data)

    def u32(number):
        return struct.pack("<I", number & 0xFFFFFFFF)

    put(u32(len(rules)))
    for fmt in sorted(rules):
        prefixes, wildcards, big, small = rules[fmt]
        put(u32(fmt))
        for items in (sorted(prefixes.items(), key=lambda kv: kv[0].encode()), wildcards):
            put(u32(len(items)))
            for key, value in items:
                put(key.encode())
                put(u32(value))
        put(u32(len(big)))
        for last in sorted(big, key=str.encode):
            first, value = big[last]
            put(last.encode())
            put(first.encode())
            put(u32(value))
        kept = {}
        for key, options in small.items():
            options = [option for option in options if option[0].starts]
            if options:
                kept[key] = options
        put(u32(len(kept)))
        for key in sorted(kept, key=str.encode):
            put(key.encode())
            put(u32(len(kept[key])))
            for bits, value in kept[key]:
                put(u32(value))
                put(u32(len(bits.starts)))
                for left, right in zip(bits.starts, bits.ends):
                    put(u32(left))
                    put(u32(right))
    return fnv1a64(out)


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
    block_bits = None
    for tokens in compact:
        if tokens[0] == ":":
            if block_bits is None:
                sys.exit("a continuation line that does not follow a special2 line: %r" % tokens)
            side = r"\d{1,%d}" % block_bits if block_bits else r"(?!)"
            if not all(re.fullmatch(r"%s(-%s)?" % (side, side), t) for t in tokens[1:]):
                sys.exit("a range of more digits than the bit has, or not digits: %r" % tokens)
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
            block_bits = None
            if len(tokens) >= 3 and tokens[0].lower() == "special2":
                block_bits = bit_digits(to_uint(tokens[1].partition("+")[2]))
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
    names, f_specials, f_fallback = acc_info_values(os.path.join(out_dir, "acc_info.rs"))
    print("guide fingerprint 0x%016x" % guide_fingerprint(compact, names, f_specials, f_fallback))


if __name__ == "__main__":
    main()
