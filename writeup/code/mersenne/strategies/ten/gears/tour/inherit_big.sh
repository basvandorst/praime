#!/bin/bash
# rounds of random large primes on machine 23; one JSON line per size
PY=/home/bas/praime/.venv/bin/python3
for spec in "8 3000" "10 3000" "12 3000" "15 3000" "20 3000" "25 2000" "30 1500" "35 800" "40 500" "50 200" "60 80"; do
  set -- $spec
  timeout 1500 nice $PY inherit_big.py $1 $2 2>/dev/null | grep digits
done
echo done
