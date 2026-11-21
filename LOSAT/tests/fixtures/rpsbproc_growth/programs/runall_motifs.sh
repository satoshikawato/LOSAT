#!/bin/bash
# runall_motifs.sh: runs every row of motifs.tsv (site, motif, baseline motif of the same shape) with runmotif.sh: the
# motif as run mNN and its baseline as bNN (labels of equal length, so that paths do not change the small allocations).
cd ~/libo-o1c
i=0
while IFS=$'\t' read -r site motif basem; do
  i=$((i+1))
  MALLOG_MIN=8 ./runmotif.sh "$(printf m%02d $i)" "$motif" 0 "$site"
  MALLOG_MIN=8 ./runmotif.sh "$(printf b%02d $i)" "$basem" 0 "$site"
done < motifs.tsv
touch done_motifs
