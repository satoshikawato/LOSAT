#!/usr/bin/env python3
"""Write a small NCBI taxdb (taxdb.bti + taxdb.btd) and damaged variants of it.

Format (c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp, 598d8ae6), all integers big-endian u32:
  taxdb.bti  magic 0x8739, taxid count, 4 reserved words, then `count` pairs (taxid, offset)
  taxdb.btd  the records back to back: "scientific\\tcommon\\tblast\\tkingdom"; the record of an
             entry ends where the next entry's offset starts (the last: at the end of the file).

Usage:
  make_taxdb.py write TSV OUTDIR           one good taxdb from TSV (taxid, sci, common, blast, kingdom)
  make_taxdb.py all TSV OUTROOT            every variant in OUTROOT/<name>/ and OUTROOT/variants.tsv
  make_taxdb.py list                       the variant names and what each is

Python 3 standard library only.
"""
import os
import struct
import sys

MAGIC = 0x8739


def read_tsv(path):
    rows = []
    with open(path, 'rb') as f:
        for line in f.read().split(b'\n'):
            if not line:
                continue
            cols = line.split(b'\t')
            assert len(cols) == 5, line
            rows.append((int(cols[0]), cols[1], cols[2], cols[3], cols[4]))
    return rows


def u32(v):
    return struct.pack('>I', v & 0xFFFFFFFF)


def record(row):
    return b'\t'.join(row[1:5])


def layout(rows, records=None, term=b''):
    """(entries, btd): entries = [(taxid, offset)], records back to back (each followed by `term`)."""
    if records is None:
        records = [record(r) for r in rows]
    entries, btd = [], b''
    for r, rec in zip(rows, records):
        entries.append((r[0], len(btd)))
        btd += rec + term
    return entries, btd


def bti(entries, magic=MAGIC, count=None, extra=b''):
    if count is None:
        count = len(entries)
    out = u32(magic) + u32(count) + u32(0) * 4
    for t, o in entries:
        out += u32(t) + u32(o)
    return out + extra


def good(rows):
    e, d = layout(rows)
    return bti(e), d


def with_row(rows, taxid, fn):
    return [fn(r) if r[0] == taxid else r for r in rows]


def with_record(rows, taxid, rec):
    e, d = layout(rows, [rec if r[0] == taxid else record(r) for r in rows])
    return bti(e), d


# name -> (description, builder(rows) -> (bti bytes | None, btd bytes | None))
VARIANTS = {}


def variant(name, desc):
    def deco(fn):
        VARIANTS[name] = (desc, fn)
        return fn
    return deco


@variant('good', 'sorted, 7 entries; the last entry is taxid 99990')
def _(rows): return good(rows)


@variant('bad_magic', 'magic 0x8738')
def _(rows):
    e, d = layout(rows); return bti(e, magic=0x8738), d


@variant('bad_magic_le', 'magic written little-endian')
def _(rows):
    e, d = layout(rows); b = bytearray(bti(e)); b[0:4] = struct.pack('<I', MAGIC); return bytes(b), d


@variant('count_larger', 'count = entries + 1 (room for fewer)')
def _(rows):
    e, d = layout(rows); return bti(e, count=len(e) + 1), d


@variant('count_smaller', 'count = entries - 1 (room for more)')
def _(rows):
    e, d = layout(rows); return bti(e, count=len(e) - 1), d


@variant('count_zero', 'count 0')
def _(rows):
    e, d = layout(rows); return bti(e, count=0), d


@variant('index_short', 'index file 31 bytes (header + 7)')
def _(rows):
    b, d = good(rows); return b[:31], d


@variant('index_header_only', 'index file 24 bytes (header, no entry), count 0')
def _(rows):
    e, d = layout(rows); return bti([], count=0), d


@variant('index_one_entry', 'one entry (index 32 bytes, the smallest valid)')
def _(rows):
    e, d = layout(rows[2:3]); return bti(e), d


@variant('index_trailing', 'valid index followed by 7 stray bytes (integer division ignores them)')
def _(rows):
    e, d = layout(rows); return bti(e, extra=b'\x01' * 7), d


@variant('index_trailing_8', 'valid index followed by 8 stray bytes (room for one more: count mismatch)')
def _(rows):
    e, d = layout(rows); return bti(e, extra=b'\x00' * 8), d


@variant('no_btd', 'taxdb.btd absent')
def _(rows):
    b, d = good(rows); return b, None


@variant('no_bti', 'taxdb.bti absent')
def _(rows):
    b, d = good(rows); return None, d


@variant('btd_empty', 'taxdb.btd exists with 0 bytes')
def _(rows):
    b, d = good(rows); return b, b''


@variant('btd_truncated', 'taxdb.btd cut inside the record of taxid 9606 (index unchanged)')
def _(rows):
    e, d = layout(rows)
    off = dict(e)[9606]
    return bti(e), d[:off + 5]


@variant('empty_kingdom_mid', 'taxid 9606 record ends with a tab and an empty kingdom')
def _(rows): return with_record(rows, 9606, b'Homo sapiens\thuman\tprimates\t')


@variant('empty_kingdom_last', 'last entry (99990) record ends with a tab and an empty kingdom')
def _(rows): return with_record(rows, 99990, b'Last entry sp.\tlast common\tlastblast\t')


@variant('empty_fields', 'taxid 9606: empty scientific, common and blast names, kingdom X')
def _(rows): return with_record(rows, 9606, b'\t\t\tX')


@variant('high_bytes_utf8', 'taxid 9606 scientific name with UTF-8 bytes (>= 0x80)')
def _(rows): return with_record(rows, 9606, b'Caf\xc3\xa9 sapiens\thuman\tprimates\tEukaryota')


@variant('high_bytes_latin1', 'taxid 10090 common name with a Latin-1 byte 0xE9 and kingdom with 0xFF')
def _(rows): return with_record(rows, 10090, b'Mus musculus\tmus\xe9\trodents\tEuk\xffaryota')


@variant('tabs0', 'taxid 9606 record without tab')
def _(rows): return with_record(rows, 9606, b'Homo sapiens')


@variant('tabs1', 'taxid 9606 record with 1 tab')
def _(rows): return with_record(rows, 9606, b'Homo sapiens\thuman')


@variant('tabs2', 'taxid 9606 record with 2 tabs (no kingdom, no fourth field)')
def _(rows): return with_record(rows, 9606, b'Homo sapiens\thuman\tprimates')


@variant('tabs4', 'taxid 9606 record with 4 tabs (kingdom = "Eukaryota\\textra")')
def _(rows): return with_record(rows, 9606, b'Homo sapiens\thuman\tprimates\tEukaryota\textra')


@variant('tabs0_last', 'last entry record without tab')
def _(rows): return with_record(rows, 99990, b'Last entry sp.')


@variant('term_nul', 'every record followed by a NUL byte')
def _(rows):
    e, d = layout(rows, term=b'\x00'); return bti(e), d


@variant('term_nl', 'every record followed by a newline')
def _(rows):
    e, d = layout(rows, term=b'\n'); return bti(e), d


@variant('offsets_swapped', 'offsets of 9606 and 10090 swapped (9606 would end before it begins)')
def _(rows):
    e, d = layout(rows)
    i = [t for t, _ in e].index(9606)
    e[i], e[i + 1] = (e[i][0], e[i + 1][1]), (e[i + 1][0], e[i][1])
    return bti(e), d


@variant('offsets_equal', 'taxid 10090 has the offset of 9606 (9606 gets an empty record)')
def _(rows):
    e, d = layout(rows)
    i = [t for t, _ in e].index(9606)
    e[i + 1] = (e[i + 1][0], e[i][1])
    return bti(e), d


@variant('offset_last_past_end', 'offset of the last entry is past the end of taxdb.btd')
def _(rows):
    e, d = layout(rows)
    e[-1] = (e[-1][0], len(d) + 10)
    return bti(e), d


@variant('offset_last_at_end', 'offset of the last entry equals the size of taxdb.btd (empty record)')
def _(rows):
    e, d = layout(rows)
    e[-1] = (e[-1][0], len(d))
    return bti(e), d


@variant('offset_mid_past_end', 'offset of taxid 9606 is past the end of taxdb.btd')
def _(rows):
    e, d = layout(rows)
    i = [t for t, _ in e].index(9606)
    e[i] = (e[i][0], len(d) + 10)
    return bti(e), d


@variant('offset_next_past_end', 'offset of taxid 10090 is past the end of taxdb.btd (the record of 9606 runs past the end)')
def _(rows):
    e, d = layout(rows)
    i = [t for t, _ in e].index(10090)
    e[i] = (e[i][0], len(d) + 10)
    return bti(e), d


@variant('unsorted_adjacent', 'entries 9606 and 10090 in the wrong order')
def _(rows):
    r = list(rows)
    i = [x[0] for x in r].index(9606)
    r[i], r[i + 1] = r[i + 1], r[i]
    return good(r)


@variant('dup_taxid', 'taxid 9606 twice (second copy has other names)')
def _(rows):
    r = list(rows)
    i = [x[0] for x in r].index(9606)
    r.insert(i + 1, (9606, b'Homo duplicatus', b'dup', b'dupblast', b'DupK'))
    return good(r)


@variant('taxid_high_bit', 'an entry with taxid 0x80000005 (negative as Int4) after 99990')
def _(rows):
    r = list(rows) + [(0x80000005, b'Negative sp.', b'neg', b'negblast', b'NegK')]
    e, d = layout(r)
    return bti(e), d


def write_variant(outdir, b, d):
    os.makedirs(outdir, exist_ok=True)
    for name, data in (('taxdb.bti', b), ('taxdb.btd', d)):
        p = os.path.join(outdir, name)
        if os.path.exists(p):
            os.remove(p)
        if data is not None:
            with open(p, 'wb') as f:
                f.write(data)


def main(argv):
    if len(argv) >= 2 and argv[1] == 'list':
        for n, (desc, _) in VARIANTS.items():
            print('%s\t%s' % (n, desc))
        return 0
    if len(argv) == 4 and argv[1] == 'write':
        b, d = good(read_tsv(argv[2]))
        write_variant(argv[3], b, d)
        return 0
    if len(argv) == 4 and argv[1] == 'all':
        rows = read_tsv(argv[2])
        with open(os.path.join(argv[3], 'variants.tsv'), 'w') as m:
            for n, (desc, fn) in VARIANTS.items():
                b, d = fn(rows)
                write_variant(os.path.join(argv[3], n), b, d)
                m.write('%s\t%s\n' % (n, desc))
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv))
