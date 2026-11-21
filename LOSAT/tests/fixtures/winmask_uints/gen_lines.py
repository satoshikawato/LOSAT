#!/usr/bin/env python3
"""Input lines of uint_probe: systematic (one field kind at each position of "M k roff bc") and
seeded random token sequences. Prints one escaped line per row."""
import random
base = [3694, 17, 0, 2]
U32, U64 = 1 << 32, 1 << 64
def kinds(t):
    z = lambda s: s
    out = ["0", "00", "-0", "+0", "-00", "1", "-1", "-2", "+1", str(U32 - 1), str(U32), str(U32 + 1),
           "-" + str(U32 - 1), "-" + str(U32), "-" + str(U32 + 1), "+" + str(U32 - 1), "+" + str(U32),
           str(U64 - 1), str(U64), str(U64 + 1), "-" + str(U64 - 1), "-" + str(U64), "-" + str(U64 + 1),
           "-" + str(U64 - U32 + 1), "-" + str(U64 - U32), "-" + str(U64 - U32 - 1), "-" + str(U64 - U32 + 2),
           "-", "+", "+-1", "-+1", "--1", "++1", "- 1", "-\t1", "+ 1", "-x", "x", "1x", "-1x", "+1x",
           "0x10", "010", "-010", "1e3", "1.5", "-1.5", "1,5", "1_0", "9" * 30, "-" + "9" * 30,
           "0" * 25 + "17", "-" + "0" * 25 + "17", "0" * 40, "-" + "0" * 40,
           "-" + str(U32 - t) if t else "-0", "-" + str(U64 - t) if t else "-0",
           "+" + str(t), str(t), "-" + str(t), "-" + str(U32 + t), "-" + str(U64 + t) ]
    return out
def freeze_kinds(t):
    # the kinds of the blastn freezes (tests/winmask/text_counts.tsv, taxids 101-160)
    return [("minus1", "-1"), ("minus0", "-0"), ("minus_u32max", f"-{U32 - 1}"), ("minus_u32", f"-{U32}"),
            ("wrap32", f"-{U32 - t}" if t else "-00"), ("wrap64", f"-{U64 - t}"), ("minus_u64max", f"-{U64 - 1}"),
            ("u32max", f"{U32 - 1}"), ("u32", f"{U32}"), ("u64", f"{U64}"), ("plus", f"+{t}"), ("zeros", f"{t:025d}"),
            ("sign_only", "-"), ("plus_only", "+"), ("junk", f"{t}x")]
rows = []
for p in range(4):
    for v in kinds(base[p]):
        f = [str(b) for b in base]
        f[p] = v
        rows.append(" ".join(f))
        rows.append(" ".join(f[:p + 1]))          # the line ends at the field
        if p == 3:
            rows.append(v)
        rows.append(" ".join(f) + " 99 x")
# several kinds in one line
for v in kinds(5)[:40]:
    for w in kinds(5)[:12]:
        rows.append(f"{v} {w}")
        rows.append(f"{v} {w} {v} {w}")
# separators
for sep in [" ", "  ", "\\t", " \\t ", "\\r", "\\v", "\\f", "\\r\\n", "\\x0b\\x0c", "\\x85", "\\xa0", "\\xff", "\\0", ",", ";"]:
    rows.append(sep.join(["3694", "17", "0", "2"]))
    rows.append(sep + "3694" + sep + "17" + sep + "0" + sep + "2" + sep)
    rows.append("-" + sep + "1" + sep + "2")
rows += ["", " ", "\\t", "\\r", "\\0", "\\0 1 2 3 4", "1\\0 2 3 4", "1 2\\r3 4", "\\xef\\xbb\\xbf3694 17 0 2", "3694 17 0 2\\r",
         "3694\\r17 0 2", "+", "-", "+ ", "- ", "-\\0", "4294967295 4294967295 4294967295 4294967295",
         "-4294967295 -4294967295 -4294967295 -4294967295"]
# the lines of the blastn freezes (tests/winmask/text_counts.tsv taxids 101-160, line 507)
for p in range(4):
    for _tag, v in freeze_kinds(base[p]):
        f = [str(b) for b in base]
        f[p] = v
        rows.append(" ".join(f))
# random
rnd = random.Random(20261120)
pool = ["0", "1", "2", "17", "3694", "9", "31", "32", "-1", "-2", "-17", "+5", "4294967295", "4294967296", "-4294967295",
        "-4294967296", "18446744073709551615", "18446744073709551616", "-18446744073709551615", "-18446744073709551616",
        "-", "+", "x", "1x", "0x5", "abc", "1.5", "", "99999999999999999999999", "-99999999999999999999999", "00017",
        "-00017", "2147483647", "2147483648", "-2147483648", "-2147483649", "65535", "65536", "-65536"]
seps = [" ", " ", " ", "  ", "\\t", " \\t", "\\r", "\\v", "\\f"]
for _ in range(1500):
    n = rnd.randint(1, 6)
    s = rnd.choice(["", "", "", " ", "\\t"])
    for j in range(n):
        s += rnd.choice(pool)
        if j < n - 1:
            s += rnd.choice(seps)
    s += rnd.choice(["", "", "", " ", "\\r"])
    rows.append(s)
seen = set()
for r in rows:
    if r not in seen:
        seen.add(r)
        print(r.replace("\t", "\\t"))
