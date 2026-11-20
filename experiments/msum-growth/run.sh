#!/usr/bin/env bash
# One closure per process (clean peak RSS), address space capped at 14 GB.
# usage: run.sh OUT.jsonl LOG args...
set -u
out=$1; log=$2; shift 2
here=$(cd "$(dirname "$0")/../.." && pwd)
prlimit --as=14000000000 -- "$here/target/release/msum_direct" "$@" --out "$out" >/dev/null 2>>"$log"
echo "exit=$? args=$*" >>"$log"
