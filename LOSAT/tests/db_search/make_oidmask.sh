#!/bin/bash
# Builds the fixture directory LOSAT/tests/db_search/oidmask/ (series D, item D1b-2: the alias key OID_MASK_TYPE)
# with NCBI BLAST+ 2.17.0 (comparison-only tooling, as the fixtures' NCBI outputs; LOSAT never runs NCBI).
# The committed files in oidmask/ are the fixtures, this script records how they were made:
#   om.faa      (kept, the source) 8 protein records; record 4 has two definition lines joined by byte 0x01
#               (ref|XP_000004.1| and ref|NP_000004.1|), records 2, 3 and 6 are predicted RefSeq (XP_)
#   omp, omp4   makeblastdb -parse_seqids -oid_masks 1 over om.faa, version 5 and 4 (omp.pxm, omp4.pxm: the mask)
#   ompx, omp4x blastdb_aliastool -oid_masks 1 aliases over them (OIDLIST + OID_MASK_TYPE 1)
#   omn.fna     the back-translation of om.faa (one fixed codon per amino acid), titles of the same kinds
#               (XM_ in place of XP_, NM_ in place of NP_, XR_/NR_ for records 6 and 8)
#   omn         makeblastdb -parse_seqids -oid_masks 1 over omn.fna in two volumes (omn.00, omn.01, omn.nal, .nxm each)
#   omnx        blastdb_aliastool -oid_masks 1 over omn (omnx.nal, omnx.00.nal, omnx.01.nal)
#   q.faa, q.fna  queries: record 1 of om.faa / omn.fna
#   omb, ombx   a version 4 database with the byte 0x80 in two predicted lines' titles, and its alias (below)
#   *.pal, *.nal, *.msk  hand-written alias files and OID list files (below)
# makeblastdb writes the build date into the .pjs/.njs files, so a rebuild differs from the committed files.
#
# Usage: make_oidmask.sh NCBI_BIN    (writes into the oidmask directory next to this script)
set -euo pipefail
N=$1
HERE=$(cd "$(dirname "$0")" && pwd)
O=$HERE/oidmask
cd "$O"
# om.faa (committed) is the source; every other file is made again.
find . -maxdepth 1 -type f ! -name om.faa -delete

# --- protein databases
"$N/makeblastdb" -in om.faa -dbtype prot -parse_seqids -oid_masks 1 -out omp > /dev/null
"$N/makeblastdb" -in om.faa -dbtype prot -parse_seqids -oid_masks 1 -blastdb_version 4 -out omp4 > /dev/null
"$N/blastdb_aliastool" -db omp -dbtype prot -oid_masks 1 -out ompx > /dev/null
"$N/blastdb_aliastool" -db omp4 -dbtype prot -oid_masks 1 -out omp4x > /dev/null

# --- nucleotide records and queries
python3 - <<'PY'
codon = dict(A="GCC", R="CGC", N="AAC", D="GAC", C="TGC", Q="CAG", E="GAG", G="GGC", H="CAC", I="ATC",
             L="CTG", K="AAG", M="ATG", F="TTC", P="CCC", S="AGC", T="ACC", W="TGG", Y="TAC", V="GTG",
             B="GAC", Z="GAG", X="AAA", U="TGC", J="CTG", O="AAG")
codon["*"] = "TAA"
protein = []
name = None
for line in open("om.faa"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        protein.append("")
    else:
        protein[-1] += line
assert len(protein) == 8
titles = ["ref|NM_000001.1| normal one",
          "ref|XM_000002.1| predicted two",
          "ref|XM_000003.1|gb|AY000003.1| mixed three",
          "ref|XM_000004.1| predicted four\x01ref|NM_000004.1| normal four",
          "lcl|nloc5 local five",
          "ref|XR_000006.1| predicted six",
          "emb|AJ000007.1| embl seven",
          "ref|NR_000008.1| noncoding eight"]
dna = ["".join(codon[c] for c in p) for p in protein]
with open("omn.fna", "w", newline="") as out:
    for t, d in zip(titles, dna):
        out.write(">" + t + "\n")
        for i in range(0, len(d), 80):
            out.write(d[i:i + 80] + "\n")
with open("q.faa", "w") as out:
    out.write(">q1\n" + protein[0] + "\n")
with open("q.fna", "w") as out:
    out.write(">qn1\n")
    for i in range(0, len(dna[0]), 80):
        out.write(dna[0][i:i + 80] + "\n")
PY

# --- nucleotide database in two volumes
"$N/makeblastdb" -in omn.fna -dbtype nucl -parse_seqids -oid_masks 1 -max_file_sz "${OM_MAX_FILE_SZ:-500B}" -out omn > /dev/null
"$N/blastdb_aliastool" -db omn -dbtype nucl -oid_masks 1 -out omnx > /dev/null

# --- OID list files: big-endian u32 (last OID) then the bit bytes
python3 - <<'PY'
import struct
open("ones.msk", "wb").write(struct.pack(">I", 7) + b"\xff")
open("onesn.msk", "wb").write(struct.pack(">I", 3) + b"\xff")
# OIDs 0, 3, 4, 6 and 7: every sequence keeps a line under the type (ompart.pal, makeblastdb -input_type blastdb)
open("nol.msk", "wb").write(struct.pack(">I", 7) + b"\x9b")
PY

# --- hand-written alias files over the protein database
al() { local name=$1; shift; printf '%s\n' "$@" > "$name"; }
for t in 0 2 3 max plus over abc neg empty; do
  case $t in
    max) v=4294967295 ;;
    plus) v=+1 ;;
    over) v=4294967296 ;;
    neg) v=-1 ;;
    empty) v= ;;
    *) v=$t ;;
  esac
  al "t$t.pal" "DBLIST omp" "OIDLIST omp.pxm" "OID_MASK_TYPE $v" "NSEQ 6" "LENGTH 851"
done
al tnone.pal "DBLIST omp" "OIDLIST omp.pxm" "NSEQ 6" "LENGTH 851"
al tnol.pal "DBLIST omp" "OID_MASK_TYPE abc"
al tfirst.pal "DBLIST omp" "FIRST_OID abc" "OIDLIST omp.pxm" "OID_MASK_TYPE abc"
al tscan.pal "DBLIST omp" "OIDLIST omp.pxm" "OID_MASK_TYPE 1"
al j2.pal "DBLIST ompx" "OIDLIST ones.msk"
al ompy.pal "DBLIST omp" "OIDLIST ones.msk" "NSEQ 8" "LENGTH 1133"
al j5.pal "DBLIST ompx ompy"
al j1.pal "DBLIST omp" "FIRST_OID 2" "LAST_OID 6" "OIDLIST omp.pxm" "OID_MASK_TYPE 1"
al tall.pal "DBLIST omp" "OIDLIST ones.msk" "OID_MASK_TYPE 1" "NSEQ 8" "LENGTH 1133"
al ompart.pal "TITLE omp partly erased" "DBLIST omp" "OIDLIST nol.msk" "OID_MASK_TYPE 1"
# --- over the nucleotide database: a per-volume list (bits past each volume's end are ignored)
al ntall.nal "DBLIST omn" "OIDLIST onesn.msk" "OID_MASK_TYPE 1"
# --- omb (version 4, -oid_masks 1) and ombx (blastdb_aliastool -oid_masks 1): OID 0 has a normal line and a
# predicted one whose title holds the byte 0x80, OID 1 a normal line, OID 2 a predicted one with 0x80 (masked out by
# omb.pxm). makeblastdb writes '#' for such a byte, so the stored "bad#" is patched back to "bad\x80" (as a header that
# NCBI's decoder reads as '#': it writes an error for each, also in a line OID_MASK_TYPE then erases; review A of D1b-2).
python3 - <<'PY'
seqs = []
for line in open("om.faa"):
    if line.startswith(">"):
        seqs.append("")
    else:
        seqs[-1] += line.strip()
with open("omb.faa", "w") as out:
    out.write(">ref|NP_000011.1| ok eleven\x01ref|XP_000012.1| bad#twelve\n" + seqs[0] + "\n")
    out.write(">ref|NP_000013.1| plain thirteen\n" + seqs[6] + "\n")
    out.write(">ref|XP_000014.1| bad#fourteen\n" + seqs[7] + "\n")
PY
"$N/makeblastdb" -in omb.faa -dbtype prot -parse_seqids -oid_masks 1 -blastdb_version 4 -out omb > /dev/null
python3 - <<'PY'
data = open("omb.phr", "rb").read()
assert data.count(b"bad#") == 2
open("omb.phr", "wb").write(data.replace(b"bad#", b"bad\x80"))
PY
"$N/blastdb_aliastool" -db omb -dbtype prot -oid_masks 1 -out ombx > /dev/null
rm omb.faa
# --- the seqid lists of rows db.om.bp.tall.neg.6 and .pos.6
al neg.txt WP_000008.1
al pos.txt XP_000003.1 NP_000001.1
ls -l
