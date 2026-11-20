#!/bin/bash
# dbsearch-20261012 session 14: build LOSAT at SHA (origin/dev) into db3/bin/LOSAT-db6-<sha8> on 4 idle CPUs of node 1
# (untimed; SHARED marker naming them; cpu_busy.py pick). The LOSAT crate builds mostly on one core (codegen-units=1).
set -u
SHA=$1
B=$HOME/losat-bench; ST=$B/state; D3=$B/db3; R=$D3/repo
pick=$(python3 $D3/task/scripts/cpu_busy.py pick 3 4 $ST)
[ -n "$pick" ] || { echo "no 4 idle cores"; exit 1; }
read -r node CPUS SIBS <<< "$pick"
MK=$ST/NODE${node}_SHARED_dbsearch-s14-build.txt
echo "dbsearch-20261012 session 14 (local losat-baselines-4d) since $(date '+%F %T'): node $node CPUs $CPUS ONLY (UNTIMED): cargo build --release of $SHA (~6 min); driver PID $$; removed at its end." > $MK
trap "rm -f $MK" EXIT
export PATH=$HOME/.cargo/bin:/usr/bin:/bin CARGO_BUILD_JOBS=4 TMPDIR=$D3/tmp RUSTUP_TOOLCHAIN=1.92.0
N=$D3/bin/LOSAT-db6-${SHA:0:8}
git -C $R fetch -q origin dev && git -C $R checkout -q --detach "$SHA" || { echo "FAILED checkout"; exit 1; }
touch $R/LOSAT/src/main.rs
(cd $R/LOSAT && numactl --physcpubind=$CPUS --membind=$node cargo build --release --locked --bin losat --target-dir $D3/target_db6) > $D3/logs/s14build.log 2>&1 || { echo "FAILED build"; tail -20 $D3/logs/s14build.log; exit 1; }
grep -q "Compiling LOSAT" $D3/logs/s14build.log || { echo "FAILED: no Compiling line"; exit 1; }
cp $D3/target_db6/release/losat $N && echo "$(date '+%F %T') built $N $(sha256sum $N | cut -c1-16) from $(git -C $R rev-parse HEAD) on node $node CPUs $CPUS"
