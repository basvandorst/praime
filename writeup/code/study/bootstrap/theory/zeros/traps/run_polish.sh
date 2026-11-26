#!/bin/sh
cd "$(dirname "$0")"
for d in 2 4 6 8 10 12 14 16 18 20 22 24; do
  (ulimit -v 1500000; timeout 1200 nice python polish_hp.py $d > pol_d$d.log 2>&1)
  tail -1 pol_d$d.log
done
