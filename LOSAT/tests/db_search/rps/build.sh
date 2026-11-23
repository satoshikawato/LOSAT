#!/bin/bash
# Builds the RPS-BLAST fixture databases of LOSAT/tests/db_search/rps/ with NCBI BLAST+ 2.17.0
# (comparison-only tooling, as the fixtures' NCBI outputs; LOSAT never runs NCBI):
#   rps_mini   version 5, three profiles of WSSV protein fragments (the frequency ratios of two
#              psiblast rounds against tests/db_search/paja), Seq-ids gnl|CDD|101..103, scores at
#              makeprofiledb's default scale 100 (as CDD's)
#   rps_mini4  the same profiles, version 4
#   rps_mini3v the same profiles as three volumes of one profile each (rps_mini3v.00 to .02, version 5,
#              each made by makeprofiledb from one profile, and the alias file rps_mini3v.pal that
#              lists them, written below: makeprofiledb -max_smp_vol never makes one profile per
#              volume of three profiles)
#   rps_tie    two volumes of one profile each (rps_tie.00, .01 and rps_tie.pal), both the ribonucleotide
#              reductase profile, gnl|CDD|102 and gnl|CDD|202: every hit of one ties exactly with the other's
#              (the order of the merge at equal e-values and bit scores)
#   rps_mini26 version 4 with a 26-column .rps (magic 0x1e16 in .rps and .loo, rows cut to their
#              first 26 columns), its .freq left at 28 columns: the layout of NCBI's own test
#              database rpstest_le
# makeprofiledb writes the build date into the index files, so a rebuild differs from the
# committed files: the committed .gz files are the fixtures, this script records how they were made.
#
# Not a fixture (the fixture budget is tight), but how the review of PS11c-2 (REVIEW_ps11c2 S2) filled the
# preliminary hit list past its capacity of 10 (`-max_target_seqs` 1 to 5): a database of 18 profiles, copies of 6
# of these PSSMs saved by psiblast (`-out_pssm`, as above) under 18 Seq-ids gnl|CDD|300..317 (p0.smp .. p17.smp, the
# query id rewritten as above, the scores dropped), made by makeprofiledb -in list.txt once as one volume (big18) and
# once as four volumes plus a .pal (big18v), so that one query hits up to 15 profiles. NCBI and LOSAT
# agreed on every run (`-max_target_seqs` 1, 2, 5, 10, 11, 30, `-searchsp`, `OLD_FSC`, CBS 0/1/3); the unit test
# `a_full_prelim_list_of_equal_e_values_keeps_the_best_ten_of_the_total_order` (rps_search/prelim.rs) holds the order.
#
# rps_cdd (series PR, item PS15a-1): the profiles of rps_mini under three real CDD PSSM ids, for the pipeline test (below).
#
# xp101.phr, .pin, .pdb and .pos (series D, item D1b-2): rps_mini made with profile 101 under the predicted RefSeq Seq-id
#   ref|XP_000101.1| (below), which the alias key OID_MASK_TYPE 1 erases; the other files of that build are the same as
#   rps_mini's, so the rows replace only these four (`rps-gz:rps_mini!phr=xp101.phr.gz...`, gzip -9 -n each).
# xp3v01.phr, .pin, .pdb and .pos (series D, item D1b-2): rps_mini3v.01 made from profile 102 under ref|XP_000102.1|
#   (below; the other files of that build are the same as rps_mini3v.01's), the second of three volumes searched bare
#   (`rps-gz:rps_mini3v!01.phr=xp3v01.phr.gz...@om3vrps`).
# tx3v00..02.{phr,pin,pot,ptf} (line PS, PS15b-2): rps_mini3v.00 to .02 made as below with makeprofiledb -taxid 9606,
#   10090 and 562 (`rps-gz:rps_mini3v!00.phr=tx3v00.phr.gz...`, row db.ps15b2.rps.tax3v.18).
#
# rps_wide and the .loo rps_wide_ug.loo (series PR, item PS12m): one profile of 133,000 columns (f101 tiled 700 times by
#   build_wide.py, version 5, gnl|CDD|901) with makeprofiledb -threshold 20 -force, so that its table holds the nine words that
#   match exactly (700 hits each, one per tile): 65 buckets of 2,048 rows, a cell's hits in all of them, in both groups of
#   LOSAT's scan (64 buckets and one). Its .freq is left out (the rows use -comp_based_stats 0). rps_wide_ug.loo is its table with
#   the hits of those nine cells not in ascending order (reversed, rotated, first and last swapped; build_wide.py), which
#   makeprofiledb never writes: the cell-order read of a table like it would lose hits of the later group.
#
# Usage: build.sh NCBI_BIN OUT_DIR    (then: gzip -9 -n OUT_DIR/rps_mini*.* and copy the .gz files here)
set -euo pipefail
N=$1
O=$2
HERE=$(cd "$(dirname "$0")" && pwd)
LOSAT=$(cd "$HERE/../../.." && pwd)
mkdir -p "$O"
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
cd "$W"
# The fragments: (CDD tag, WSSV record, first residue, last residue), 1-based.
cat > fragments.tsv <<'TSV'
101	YP_009220601.1	1	190
102	YP_009220529.1	301	480
103	YP_009220649.1	1101	1290
TSV
python3 - "$LOSAT/tests/fasta/WSSV.faa" <<'PY'
import sys
records, name = {}, None
for line in open(sys.argv[1]):
    line = line.rstrip("\n")
    if line.startswith(">"):
        name = line[1:].split()[0]
        records[name] = [line[1:], ""]
    else:
        records[name][1] += line
for row in open("fragments.tsv"):
    tag, acc, first, last = row.split()
    title, seq = records[acc]
    with open(f"f{tag}.faa", "w") as out:
        out.write(f">f{tag} {title.split(' ', 1)[1]} residues {first}-{last}\n{seq[int(first) - 1:int(last)]}\n")
PY
for tag in 101 102 103; do
  "$N/psiblast" -query f$tag.faa -db "$LOSAT/tests/db_search/paja" -num_iterations 2 \
    -out_pssm f$tag.smp -save_pssm_after_last_round -num_threads 1 > /dev/null
  # The PSSM's query id becomes gnl|CDD|<tag> (CDD's ids), its title stays the fragment's; its
  # scores (`finalData`) are dropped, so that makeprofiledb computes them from the frequency
  # ratios at its -scale (default 100), as for CDD's profiles (a PSSM with scores keeps its own
  # scaling factor, 1 for psiblast's: makeprofiledb.cpp x_RPSAddFirstSequence).
  python3 - f$tag.smp $tag <<'PY'
import re, sys
path, tag = sys.argv[1], sys.argv[2]
text = open(path).read()
new, count = re.subn(r'id \{\s*local str "[^"]*"\s*\}',
                     'id {\n        general {\n          db "CDD",\n          tag id ' + tag + '\n        }\n      }',
                     text, count=1)
assert count == 1, path
start = new.index(",\n    finalData {")
end = new.index("\n    }", start) + len("\n    }")
new = new[:start] + new[end:]
open(path, "w").write(new)
PY
  echo f$tag.smp >> list.txt
done
"$N/makeprofiledb" -in list.txt -out "$O/rps_mini" -title rps_mini -dbtype rps
"$N/makeprofiledb" -in list.txt -out "$O/rps_mini4" -title rps_mini4 -dbtype rps -blastdb_version 4
i=0
for tag in 101 102 103; do
  echo f$tag.smp > list3v.txt
  "$N/makeprofiledb" -in list3v.txt -out "$O/rps_mini3v.0$i" -title rps_mini3v -dbtype rps
  i=$((i + 1))
done
# NSEQ and LENGTH are the sums over the volumes, as makeprofiledb writes them into its alias files.
nseq=0; len=0
for i in 00 01 02; do
  info=$("$N/blastdbcmd" -db "$O/rps_mini3v.$i" -info | grep ' sequences;')
  nseq=$((nseq + $(echo "$info" | awk '{print $1}')))
  len=$((len + $(echo "$info" | awk '{gsub(",", "", $3); print $3}')))
done
printf '#\n# Alias file created by build.sh\n#\nTITLE rps_mini3v\nDBLIST "rps_mini3v.00" "rps_mini3v.01" "rps_mini3v.02" \nNSEQ %s\nLENGTH %s\n' "$nseq" "$len" > "$O/rps_mini3v.pal"
sed 's/tag id 102/tag id 202/' f102.smp > g102.smp
echo f102.smp > tie0.txt; echo g102.smp > tie1.txt
"$N/makeprofiledb" -in tie0.txt -out "$O/rps_tie.00" -title rps_tie -dbtype rps
"$N/makeprofiledb" -in tie1.txt -out "$O/rps_tie.01" -title rps_tie -dbtype rps
tie_len=$(for i in 00 01; do "$N/blastdbcmd" -db "$O/rps_tie.$i" -info | grep ' sequences;' | awk '{gsub(",", "", $3); print $3}'; done | paste -sd+ | bc)
printf '#\n# Alias file created by build.sh\n#\nTITLE rps_tie\nDBLIST "rps_tie.00" "rps_tie.01" \nNSEQ 2\nLENGTH %s\n' "$tie_len" > "$O/rps_tie.pal"
"$N/makeprofiledb" -in list.txt -out "$O/rps_mini26" -title rps_mini26 -dbtype rps -blastdb_version 4
# rps_mini26: the .rps rows cut to 26 columns and the 26-letter magic in .rps and .loo.
python3 - "$O/rps_mini26" <<'PY'
import struct, sys
base = sys.argv[1]
rps = open(base + ".rps", "rb").read()
magic, n = struct.unpack("<2i", rps[:8])
assert magic == 0x1e17
offsets = struct.unpack(f"<{n + 1}i", rps[8:8 + 4 * (n + 1)])
start = 8 + 4 * (n + 1)
rows = offsets[n]
assert len(rps) == start + rows * 28 * 4
out = bytearray(struct.pack("<2i", 0x1e16, n)) + rps[8:start]
for i in range(rows):
    row = start + i * 28 * 4
    out += rps[row:row + 26 * 4]
open(base + ".rps", "wb").write(bytes(out))
loo = bytearray(open(base + ".loo", "rb").read())
assert struct.unpack("<i", loo[:4])[0] == 0x1e17
loo[:4] = struct.pack("<i", 0x1e16)
open(base + ".loo", "wb").write(bytes(loo))
PY
# rps_cdd (series PR, item PS15a-1; the pipeline test `tests/rps_pipeline.rs`): the three profiles of rps_mini under the
# PSSM ids of three rows of CDD v3.22 that the data files of tests/fixtures/rpsbproc hold (177038 CHL00113, 439874 COG0104,
# 440288 COG0522), so that NCBI's and LOSAT's rpsbproc find them in cddid.tbl.
n=0
: > listcdd.txt
for tag in 177038 439874 440288; do
  n=$((n + 1))
  sed "s/tag id 10$n\$/tag id $tag/" f10$n.smp > c$tag.smp
  grep -q "tag id $tag\$" c$tag.smp
  echo c$tag.smp >> listcdd.txt
done
"$N/makeprofiledb" -in listcdd.txt -out "$O/rps_cdd" -title rps_cdd -dbtype rps
# xp101: profile 101 under `other { accession "XP_000101", version 1 }` (ref|XP_000101.1|), then 102 and 103 as in rps_mini.
python3 - f101.smp x101.smp <<'PY'
import re, sys
text = open(sys.argv[1]).read()
new, count = re.subn(r'id \{\s*general \{\s*db "CDD",\s*tag id 101\s*\}\s*\}',
                     'id {\n        other {\n          accession "XP_000101",\n          version 1\n        }\n      }', text, count=1)
assert count == 1
open(sys.argv[2], "w").write(new)
PY
printf '%s\n' x101.smp f102.smp f103.smp > listxp.txt
mkdir xp
"$N/makeprofiledb" -in listxp.txt -out xp/rps_mini -title rps_mini -dbtype rps
for ext in phr pin pdb pos; do cp "xp/rps_mini.$ext" "$O/xp101.$ext"; done
# xp3v01: profile 102 alone under `other { accession "XP_000102", version 1 }`, as volume rps_mini3v.01.
python3 - f102.smp x102.smp <<'PY'
import re, sys
text = open(sys.argv[1]).read()
new, count = re.subn(r'id \{\s*general \{\s*db "CDD",\s*tag id 102\s*\}\s*\}',
                     'id {\n        other {\n          accession "XP_000102",\n          version 1\n        }\n      }', text, count=1)
assert count == 1
open(sys.argv[2], "w").write(new)
PY
echo x102.smp > listxp3v.txt
mkdir xp3v
"$N/makeprofiledb" -in listxp3v.txt -out xp3v/rps_mini3v.01 -title rps_mini3v -dbtype rps
for ext in phr pin pdb pos; do cp "xp3v/rps_mini3v.01.$ext" "$O/xp3v01.$ext"; done
# rps_wide and its table with ungrouped cells.
python3 "$HERE/build_wide.py" tile f101.smp 700 wide.smp
echo wide.smp > wide.txt
"$N/makeprofiledb" -in wide.txt -out "$O/rps_wide" -title rps_wide -dbtype rps -threshold 20 -force
rm "$O/rps_wide.freq"
python3 "$HERE/build_wide.py" ungroup "$O/rps_wide.loo" "$O/rps_wide_ug.loo"
ls -l "$O"
