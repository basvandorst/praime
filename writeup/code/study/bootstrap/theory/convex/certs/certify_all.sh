#!/bin/sh
# Rigorous certificates for every class-A optimum: eps = 10^(floor(log10 cost) - 3), Taylor order 16,
# U = 5 and b_j rounded to 25 digits for deg <= 32; U = 6 and 60 digits for deg >= 36 (polish_A_d*.json must exist).
PY=python
for d in ${DLIST:-2 4 6 8 10 12 14 16 18 20 22 24}; do
  [ -f polish_A_d$d.json ] || continue
  den=$($PY -c "import json,math; c=float(json.load(open('polish_A_d$d.json'))['cost']); print(10**(3 - math.floor(math.log10(c))))")
  if [ $d -le 16 ]; then U=5; SIG=25; else U=6; SIG=60; fi
  (export SIG; ulimit -v 1500000; timeout 900 nice $PY certify.py $d 1/$den $U 16 > cert_d$d.out 2>&1)
  echo "deg $((2*d)): $(tail -1 cert_d$d.out)"
done
