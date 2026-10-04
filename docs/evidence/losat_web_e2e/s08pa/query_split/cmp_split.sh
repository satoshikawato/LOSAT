#!/bin/bash
# cmp_split.sh PROGRAM QUERY SUBJECT "ENV" "ARGS..."  -> prints one verdict line
P=$1; Q=$2; S=$3; E=$4; shift 4; A="$*"
L=${LOSAT:-/home/kawato/.cache/losat-web-gui-target/s08pa/bin/LOSAT-split}
env $E timeout 600 /home/kawato/micromamba/bin/$P -query $Q -subject $S $A > n.out 2> n.err; nrc=$?
env $E timeout 600 $L $P -query $Q -subject $S $A > l.out 2> l.err; lrc=$?
if cmp -s n.out l.out && cmp -s n.err l.err && [ $nrc = $lrc ]; then v=same; elif grep -q "not supported by LOSAT" l.err; then v="rejects: $(head -c 120 l.err)"; else v="DIFF rc $nrc/$lrc out $(diff n.out l.out | grep -c '^[<>]') err $(diff n.err l.err | grep -c '^[<>]')"; fi
echo "$P [$E] $(basename $Q) $(basename $S) $A :: $v"
