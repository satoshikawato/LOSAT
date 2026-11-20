#!/usr/bin/env python3
"""expected.tsv for the Rust tests, from out/cmd_<taxid>_<variant>.{out,err,rc} (probe_taxdb.sh stage cmd).

Columns: variant, taxid, expect, messages. `expect` is the escaped `%T|%S|%L|%B|%K` line of
`blastdbcmd -tax_info` (bytes outside 0x20..0x7e and the backslash as \\xNN), `N/A` when NCBI found nothing
(blastdbcmd prints either N/A fields or `T||||`), `CRASH` for a process killed by a signal,
`messages` are the stderr lines without the `Error: [blastdbcmd] ` prefix, joined with U+001F.
Usage: make_expected.py OUTDIR VARIANTS_TSV > expected.tsv
"""
import os
import re
import sys

OVERRIDES = {}


def esc(b):
    return ''.join(chr(c) if 0x20 <= c < 0x7f and c != 0x5c else '\\x%02x' % c for c in b)


def main(outdir, variants):
    for v in [l.split('\t')[0] for l in open(variants)]:
        for t in (9606, 10090, 99990, 7777, 32630, 12908, 562, 777, 3, 100000, 0):
            base = os.path.join(outdir, 'cmd_%d_%s' % (t, v))
            rc = int(open(base + '.rc').read())
            out = open(base + '.out', 'rb').read()
            err = open(base + '.err', 'rb').read().decode('latin-1')
            msgs = [re.sub(r'^Error: \[blastdbcmd\] ', '', l) for l in err.split('\n') if l]
            if rc != 0:
                exp = 'CRASH'
            else:
                head, _, rest = out.partition(b'\n')
                assert head.startswith(b'# of Tax IDs in Database: '), (v, t, head)
                assert rest.endswith(b'\n'), (v, t)
                payload = rest[:-1]
                if payload in (b'%d|N/A|N/A|N/A|N/A' % t, b'%d||||' % t):
                    exp = 'N/A'
                else:
                    exp = esc(payload)
            exp = OVERRIDES.get((v, t), exp)
            print('%s\t%d\t%s\t%s' % (v, t, exp, '\x1f'.join(msgs)))


main(sys.argv[1], sys.argv[2])
