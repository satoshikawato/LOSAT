"""Hand-edit a version 4 makeblastdb nucleotide database so that the stored 2-bit base at every ambiguous
position differs from makeblastdb's (the ambiguity data stay as they are).

makeblastdb stores a random base (CRandom(length) draws, the ones `-subject` uses) at each ambiguous letter and the
ambiguity data beside it. The rule here: the stored base becomes the first base compatible with the letter
(R: A G, Y: C T, ..., N: A C G T) that is not the stored one (the next base when there is none). A `-db` search that translates the
stored bases then differs from one that re-resolves the letters as `-subject` does (fixture rows db.tblastx.flipped.*).

Usage (from LOSAT/): python3 tests/db_search/flip_stored_bases.py FASTA NIN NSQ OUT_NSQ
(FASTA is the file makeblastdb read; NIN and NSQ the database to read; OUT_NSQ may be NSQ itself).
"""
import struct, sys
fasta, nin, nsq, out = sys.argv[1:5]
seqs, cur = [], None
for line in open(fasta):
    if line.startswith(">"):
        cur = []; seqs.append(cur)
    else:
        cur.append(line.strip())
seqs = ["".join(s) for s in seqs]
b = open(nin, "rb").read()
p = 0
def i32():
    global p; v = struct.unpack(">i", b[p:p+4])[0]; p += 4; return v
ver = i32(); typ = i32(); assert ver == 4 and typ == 0, (ver, typ)
tl = i32(); p += tl
dl = i32(); p += dl
n = i32(); p += 8; i32()
hdr = [i32() for _ in range(n + 1)]
seq = [i32() for _ in range(n + 1)]
amb = [i32() for _ in range(n)]
data = bytearray(open(nsq, "rb").read())
COMPAT = {"R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG", "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT", "-": "ACGT"}
flipped = 0
for k, s in enumerate(seqs):
    assert len(s) <= (amb[k] - seq[k]) * 4, (k, len(s))
    for i, c in enumerate(s.upper()):
        if c in "ACGTU":
            continue
        byte = seq[k] + i // 4; shift = 6 - 2 * (i % 4)
        old = (data[byte] >> shift) & 3
        opts = [x for x in COMPAT.get(c, "ACGT") if "ACGT".index(x) != old] or ["ACGT"[(old + 1) % 4]]
        new = "ACGT".index(opts[0])
        data[byte] = (data[byte] & ~(3 << shift) & 0xFF) | (new << shift)
        flipped += 1
open(out, "wb").write(data)
print("sequences", n, "flipped", flipped)
