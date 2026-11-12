#!/bin/sh
# The numbers with 8 decision histories (M-algorithm) in RIEMANN.md section 12 (about 40 minutes).
# Run from the repository root with the virtualenv active; can run next to band_validate.sh.
set -e
T=$(mktemp -d)
D=writeup/data/riemann
R=writeup/code/riemann
BD="python3 $R/band_decoder.py"
ALL=$D/odlyzko_zeros3.txt+$D/odlyzko_zeros4.txt+$D/odlyzko_zeros5.txt
node $R/zeros.js 1000000 $T/z1m.f64
echo "== with 8 decision histories"
for z in 3 4 5; do $BD $D/odlyzko_zeros$z.txt 20000 --M 8 --out $T/m$z.txt; done
$BD $ALL 40000 --M 8 --out $T/m345.txt
$BD $ALL 70000 --M 8 --wheel 2,3,5,7 --out $T/m345w210.txt
$BD low:1342:$T/z1m.f64 12000 --M 8
$BD band:500000:1342:$T/z1m.f64 12000 --M 8
echo "== cold start with 8 decision histories"
for s in 3300 4000 4600 5200; do $BD $D/odlyzko_zeros5.txt $((s + 600)) --cold $s --M 8 --stop 1 --out $T/mc5_$s.txt; done
for s in 9500 10000 11000; do $BD $ALL $((s + 1000)) --cold $s --M 8 --stop 1 --out $T/mc345_$s.txt; done
echo "== independent checks of the lists (sympy.isprime)"
python3 $R/band_check.py list $T/m3.txt $T/m4.txt $T/m5.txt $T/m345.txt $T/m345w210.txt $T/mc5_3300.txt $T/mc5_4000.txt $T/mc345_10000.txt
rm -r $T
