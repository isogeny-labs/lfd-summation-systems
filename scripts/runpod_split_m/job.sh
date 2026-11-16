#!/usr/bin/env bash
# Runs on the pod: one cell of experiments/ic-bounds-split-m (arguments in cell.env).
set -euo pipefail
cd /work
. ./cell.env
sha256sum bundle.tar.gz
mkdir -p repo && tar xzf bundle.tar.gz -C repo
python3 -m pip install -q numpy==2.5.3
cc --version | head -1; nproc; grep MemTotal /proc/meminfo; python3 -c "import numpy; numpy.show_config()" | grep -i -m3 "name\|openblas" || true
cd repo
python3 -c "from csd import elim; print('compiled elimination:', elim.compiled())"
python3 scripts/run_ic_split_m.py cell --m "$M" --k "$K" --primes $PRIMES --seeds "$SEEDS" --cap "$CAP" &
pid=$!
while kill -0 $pid 2>/dev/null; do cp -f experiments/ic-bounds-split-m/cell-m${M}-k${K}.jsonl /work/out/ 2>/dev/null || true; sleep 30; done
wait $pid; rc=$?
cp -f experiments/ic-bounds-split-m/cell-m${M}-k${K}.jsonl /work/out/ 2>/dev/null || true
echo "job exit $rc"; exit $rc
