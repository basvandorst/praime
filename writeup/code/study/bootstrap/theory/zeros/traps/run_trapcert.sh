#!/bin/sh
# D3 for the nodes given as arguments (default 2..12), sequentially
cd "$(dirname "$0")"
for k in ${@:-2 3 4 5 6 7 8 9 10 11 12}; do
  (ulimit -v 1500000; timeout 2400 nice python trapcert.py $k > trapcert_k$k.stdout 2>&1)
  grep RESULT trapcert_k$k.log || (echo "k=$k FAILED"; tail -3 trapcert_k$k.stdout)
done
