#!/bin/bash
# every odd n in [2^k, 2^(k+1)) on the norm-2 wheel, k = 1..31: how many come home
for k in $(seq 1 31); do
  lo=$((1<<k)); hi=$((1<<(k+1)))
  echo -n "k $k " ; ./norm2 range $lo $hi 2>&1 | head -n 1
done
echo done
