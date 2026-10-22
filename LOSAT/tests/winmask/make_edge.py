#!/usr/bin/env python3
"""Edge-case query records for the WindowMasker oracle, cut from a FASTA slice (deterministic)."""
import sys

src, out = sys.argv[1], sys.argv[2]
s = ''.join(l.strip() for l in open(src) if not l.startswith('>'))
recs = []
a = s[0:20000]
recs.append(('edge_plain', a))
# lowercase stretch inside a repeat-rich piece
b = s[20000:40000]
recs.append(('edge_lower', b[:5000] + b[5000:10000].lower() + b[10000:]))
# N run and single ambiguity letters
c = s[40000:60000]
c = c[:3000] + 'N' * 1000 + c[4000:8000] + 'R' + c[8001:9000] + 'Y' + c[9001:9002] + 'K' + c[9003:]
recs.append(('edge_ambig', c))
# ambiguity every 13 letters (shorter than a window of unit+4)
d = list(s[60000:62000])
for i in range(0, len(d), 13):
    d[i] = 'N'
recs.append(('edge_dense_n', ''.join(d)))
for n in (1, 9, 10, 13, 14, 15, 30):
    recs.append(('edge_len%d' % n, s[70000:70000 + n]))
recs.append(('edge_alln', 'N' * 500))
recs.append(('edge_tail_n', s[80000:81000] + 'N' * 7))
with open(out, 'w') as o:
    for name, seq in recs:
        o.write('>%s\n' % name)
        for i in range(0, len(seq), 60):
            o.write(seq[i:i + 60] + '\n')
