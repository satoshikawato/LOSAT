#!/bin/bash
# Does NCBI BLAST+ 2.17.0 read a single-volume `makeblastdb -gi_mask` database once the alias file's DBLIST line
# holds the database's file name without a directory instead of NCBI's absolute path? (PD-LOSAT-NCBI-DEFECTS
# approved exception 23; line NP item M4c-1, 2026-11-20.)
# Builds single-volume -gi_mask databases with NCBI (version 5 and 4, nucleotide and protein, -out in a
# subdirectory, the default GI-mask name), rewrites only the DBLIST line, then reads them with NCBI's blastdbcmd
# (-info, %m0, -mask_sequence_with 0, -metadata) and blastn/blastp -db_soft_mask 0, from the database's directory
# and from another directory, with a database of several volumes (which NCBI reads as it writes it) and the
# mask-column database as references; also NCBI's own rebuilds into the changed database.
# Usage: probe.sh NCBI_BIN WORK OUT  (WORK and OUT are made anew; the inputs are those of
# LOSAT/tests/seqdb/make_gi_mask_fixtures.py). Results of the run of 2026-11-20 on Linux x86-64: README.md.
set -u
export BLAST_USAGE_REPORT=false
NB=$(cd "${1:?NCBI_BIN}" && pwd)
export PATH=$NB:$PATH
W=${2:?WORK}
OUT=${3:?OUT}
GEN=$(cd "$(dirname "$0")/../../../LOSAT/tests/seqdb" && pwd)
rm -rf $W $OUT; mkdir -p $W/other $OUT
LOG=reads.txt

python3 -I - <<EOF
import sys
sys.path.insert(0, "$GEN")
import make_gi_mask_fixtures as m
from pathlib import Path
d = Path("$W/in"); d.mkdir(parents=True)
m.write_inputs(d); m.write_masks(Path("$NB"), d)
EOF
IN=$W/in

{ echo "### host $(hostname) date $(date '+%F %T %Z') TZ=${TZ:-unset}"; makeblastdb -version; blastdbcmd -version; blastn -version; } > $OUT/00_env.txt 2>&1

# r NAME DIR CMD...: run CMD in DIR, record command, status, stdout, stderr (the work path shown as <W>).
r() {
  local name=$1 dir=$2; shift 2
  ( cd "$dir" && "$@" ) > $W/o.txt 2> $W/e.txt; local rc=$?
  { echo "### $name"; echo "### CWD ${dir#$W/}"; echo "### CMD $*"; echo "### RC $rc"
    echo "### STDOUT"; sed "s#$W#<W>#g" $W/o.txt; echo "### STDERR"; sed "s#$W#<W>#g" $W/e.txt; echo; } >> $OUT/$LOG
  cp $W/o.txt $W/last.out
  return $rc
}

# show NAME FILE: the bytes of an alias file (od -c), work path shown as <W>.
show() { { echo "### $1 ($2)"; sed "s#$W#<W>#g" "$2" | od -c; echo; } >> $OUT/$LOG; }

# repair ALIAS: the one change of decision B: DBLIST <path> -> DBLIST <file name of path>.
repair() {
  python3 -I - "$1" <<'EOF'
import sys, os
p = sys.argv[1]
lines = open(p, "rb").read().split(b"\n")
out = []
for l in lines:
    if l.startswith(b"DBLIST "):
        l = b"DBLIST " + os.path.basename(l[len(b"DBLIST "):])
    out.append(l)
open(p, "wb").write(b"\n".join(out))
EOF
}

NARGS="-dbtype nucl -parse_seqids -title t"
PARGS="-dbtype prot -parse_seqids -title pt"
mk() { mkdir -p $W/$1; for f in "${@:2}"; do cp $IN/$f $W/$1/; done; }

# ---- build with NCBI
mk v5n nucl.fa dust_nucl.asnb
r build_v5n $W/v5n makeblastdb -in nucl.fa $NARGS -out db -mask_data dust_nucl.asnb -gi_mask -gi_mask_name gm1
mk v4n nucl.fa dust_nucl.asnb
r build_v4n $W/v4n makeblastdb -in nucl.fa $NARGS -out db -blastdb_version 4 -mask_data dust_nucl.asnb -gi_mask -gi_mask_name gm1
mk v5p prot.fa seg_prot.asnb
r build_v5p $W/v5p makeblastdb -in prot.fa $PARGS -out db -mask_data seg_prot.asnb -gi_mask -gi_mask_name pm1
mk v4p prot.fa seg_prot.asnb
r build_v4p $W/v4p makeblastdb -in prot.fa $PARGS -out db -blastdb_version 4 -mask_data seg_prot.asnb -gi_mask -gi_mask_name pm1
mk v5sub nucl.fa dust_nucl.asnb
r build_v5sub $W/v5sub makeblastdb -in nucl.fa $NARGS -out sub/db -mask_data dust_nucl.asnb -gi_mask -gi_mask_name gm1
mk v5def nucl.fa dust_nucl.asnb
r build_v5def $W/v5def makeblastdb -in nucl.fa $NARGS -out db -mask_data dust_nucl.asnb -gi_mask
# references: several volumes (NCBI's readable form) and the OID column
mk mv nucl.fa dust_nucl.asnb
r build_mv $W/mv makeblastdb -in nucl.fa $NARGS -out db -max_file_sz 400 -mask_data dust_nucl.asnb -gi_mask -gi_mask_name gm1
mk oid nucl.fa dust_nucl.asnb
r build_oid $W/oid makeblastdb -in nucl.fa $NARGS -out db -mask_data dust_nucl.asnb
mk oidp prot.fa seg_prot.asnb
r build_oidp $W/oidp makeblastdb -in prot.fa $PARGS -out db -mask_data seg_prot.asnb

# ---- before the repair: NCBI cannot read its own single-volume database
r before_info_v5n $W/v5n blastdbcmd -db db -dbtype nucl -info
r before_metadata_v5n $W/v5n blastdbcmd -db db -dbtype nucl -metadata

# ---- repair (decision B) and read
for c in v5n:db.nal:nucl v4n:db.nal:nucl v5p:db.pal:prot v4p:db.pal:prot v5sub:sub/db.nal:nucl v5def:db.nal:nucl; do
  d=${c%%:*}; rest=${c#*:}; a=${rest%%:*}; t=${rest#*:}; db=${a%.?al}
  show "alias_before_$d" $W/$d/$a
  cp $W/$d/$a $W/$d/$a.ncbi
  repair $W/$d/$a
  show "alias_after_$d" $W/$d/$a
  r info_$d $W/$d blastdbcmd -db $db -dbtype $t -info
  r fmt_m0_$d $W/$d blastdbcmd -db $db -dbtype $t -entry all -outfmt "%o %a %m0"
  cp $W/last.out $W/$d.m0
  r msw0_$d $W/$d blastdbcmd -db $db -dbtype $t -entry all -mask_sequence_with 0
  r metadata_$d $W/$d blastdbcmd -db $db -dbtype $t -metadata
  cp $W/last.out $W/$d.metadata
  if [ $t = nucl ]; then
    r bn_soft0_$d $W/$d blastn -query nucl.fa -db $db -outfmt 6 -db_soft_mask 0
    cp $W/last.out $W/$d.bn_soft0
    r bn_hard0_$d $W/$d blastn -query nucl.fa -db $db -outfmt 6 -db_hard_mask 0
    cp $W/last.out $W/$d.bn_hard0
    r bn_nomask_$d $W/$d blastn -query nucl.fa -db $db -outfmt 6
    cp $W/last.out $W/$d.bn_nomask
  else
    r bp_soft0_$d $W/$d blastp -query prot.fa -db $db -outfmt 6 -db_soft_mask 0
    cp $W/last.out $W/$d.bp_soft0
  fi
done

# ---- references from their own directories
r fmt_m0_mv $W/mv blastdbcmd -db db -dbtype nucl -entry all -outfmt "%o %a %m0"; cp $W/last.out $W/mv.m0
r metadata_mv $W/mv blastdbcmd -db db -dbtype nucl -metadata
r fmt_m11_oid $W/oid blastdbcmd -db db -dbtype nucl -entry all -outfmt "%o %a %m11"; cp $W/last.out $W/oid.m11
r bn_soft11_oid $W/oid blastn -query nucl.fa -db db -outfmt 6 -db_soft_mask 11; cp $W/last.out $W/oid.bn_soft11
r bn_hard11_oid $W/oid blastn -query nucl.fa -db db -outfmt 6 -db_hard_mask 11; cp $W/last.out $W/oid.bn_hard11
r fmt_m21_oidp $W/oidp blastdbcmd -db db -dbtype prot -entry all -outfmt "%o %a %m21"; cp $W/last.out $W/oidp.m21
r bp_soft21_oidp $W/oidp blastp -query prot.fa -db db -outfmt 6 -db_soft_mask 21; cp $W/last.out $W/oidp.bp_soft21

# ---- from another directory (D5): repaired single volume and NCBI's multi-volume database side by side
for d in v5n mv; do
  r other_info_$d $W/other blastdbcmd -db $W/$d/db -dbtype nucl -info
  r other_m0_$d $W/other blastdbcmd -db $W/$d/db -dbtype nucl -entry all -outfmt "%o %a %m0"
  r other_bn_soft0_$d $W/other blastn -query $IN/nucl.fa -db $W/$d/db -outfmt 6 -db_soft_mask 0
  r other_blastdb_m0_$d $W/other env BLASTDB=$W/$d blastdbcmd -db db -dbtype nucl -entry all -outfmt "%o %a %m0"
  r other_metadata_$d $W/other blastdbcmd -db $W/$d/db -dbtype nucl -metadata
done
r abs_m0_v5n $W/v5n blastdbcmd -db $W/v5n/db -dbtype nucl -entry all -outfmt "%o %a %m0"
r sub_parent_m0_v5sub $W/v5sub blastdbcmd -db sub/db -dbtype nucl -entry all -outfmt "%o %a %m0"

# ---- NCBI's rebuild into the repaired database (x01_2 with -gi_mask, x01_3 without)
cp -a $W/v5n $W/re1; cp -a $W/v5n $W/re2
r rebuild_gi_re1 $W/re1 makeblastdb -in nucl.fa $NARGS -out db -mask_data dust_nucl.asnb -gi_mask -gi_mask_name gm1
r rebuild_files_re1 $W/re1 ls -l --time-style=+ .
r rebuild_oid_re2 $W/re2 makeblastdb -in nucl.fa $NARGS -out db -mask_data dust_nucl.asnb
r rebuild_files_re2 $W/re2 ls -l --time-style=+ .

# ---- comparisons
{
  echo "### diff GI %m0 (repaired v5n) vs OID %m11 (oid): expected only the GI-less OIDs 4 and 8 (D4)"
  diff $W/v5n.m0 $W/oid.m11
  for d in v4n v5sub v5def mv; do echo "### diff %m0 v5n vs $d"; diff $W/v5n.m0 $W/$d.m0; done
  echo "### diff %m0 v5p vs v4p"; diff $W/v5p.m0 $W/v4p.m0
  echo "### diff %m0 v5p vs OID %m21 (oidp)"; diff $W/v5p.m0 $W/oidp.m21
  echo "### diff blastn soft0 v5n vs OID soft11"; diff $W/v5n.bn_soft0 $W/oid.bn_soft11
  echo "### diff blastn hard0 v5n vs OID hard11"; diff $W/v5n.bn_hard0 $W/oid.bn_hard11
  echo "### diff blastn soft0 v5n vs unmasked v5n (the masks act)"; diff $W/v5n.bn_soft0 $W/v5n.bn_nomask | head -20
  echo "### diff blastp soft0 v5p vs OID soft21"; diff $W/v5p.bp_soft0 $W/oidp.bp_soft21
  for d in v5n v4n v5p v4p v5sub v5def; do echo "### metadata $d"; cat $W/$d.metadata; done
} > $OUT/compare.txt 2>&1

echo "PROBE END"
