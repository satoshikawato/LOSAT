#!/bin/bash
# runmotif.sh LABEL MOTIF [AS_BYTES] [cat|dim]: NCBI rpsbproc 0.5.1 on mfirst.asn with the catalytic-site motif of
# cd00438 replaced by MOTIF (cat: row 2, 3 mapped residues; dim: row 3, 24 positions, row 2's motif cleared) (the data of tests/fixtures/rpsbproc `sites` + the motif); logs the
# allocations of 16 bytes or more with the malloc shim (mallog.so) to mal_LABEL.log; stdout/stderr/rc
# to out_LABEL.{out,err,rc}.
cd ~/libo-o1c/rps || exit 2
label=$1; motif=$2; as=${3:-0}; site=${4:-cat}
rm -rf d_$label; cp -r data d_$label; cp data_with/sites/* d_$label/
python3 - "$label" "$motif" "$site" <<'PY'
import sys
label, motif, site = sys.argv[1:4]
p = f"d_{label}/cddannot.dat"
lines = open(p, newline="").read().split("\n")
f = lines[1].split("\t"); assert f[5] == "HKY", f
g = lines[2].split("\t"); assert g[4] == "homodimer interface", g
if site == "cat":
    f[5] = motif
else:  # the 24-position site is the first motif site of the process
    f[5] = "0"; g[5] = motif
lines[1] = "\t".join(f); lines[2] = "\t".join(g)
open(p, "w", newline="").write("\n".join(lines))
PY
if [ "$as" != 0 ]; then lim="prlimit --as=$as"; else lim=""; fi
LD_PRELOAD=$HOME/libo-o1c/mallog.so MALLOG_FILE=mal_$label.log MALLOG_MIN=${MALLOG_MIN:-16} timeout 120 $lim \
  $HOME/losat-bench/tools/env_rpsbproc-0.5.1/bin/rpsbproc -i arch/mfirst.asn -d d_$label \
  > out_$label.out 2> out_$label.err < /dev/null
echo $? > out_$label.rc
