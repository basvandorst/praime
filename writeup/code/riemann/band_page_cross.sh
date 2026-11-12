#!/bin/sh
# Transformation 141 of prime.html against band_decoder.py (RIEMANN.md section 12, "On the page").
# DIR holds the decided lists written by   OUT=DIR node writeup/code/riemann/band_page_check.js reach
# For k = 1 and k = 4 on every band of the band slider, band_decoder.py (numpy, mpmath phases) runs on the same zeros
# up to the same limit and must decide exactly the same primes, one by one.
# usage: sh writeup/code/riemann/band_page_cross.sh DIR   (from the repository root, virtualenv active; about 45 minutes)
set -e
T=$1
D=writeup/data/riemann
R=writeup/code/riemann
[ -f $T/z.f64 ] || node $R/zeros.js 1010000 $T/z.f64
i=0
for spec in zeros:0:10000:$T/z.f64 zeros:10000:10000:$T/z.f64 zeros:100000:10000:$T/z.f64 zeros:1000000:10000:$T/z.f64 \
    $D/odlyzko_zeros3.txt $D/odlyzko_zeros4.txt $D/odlyzko_zeros5.txt $D/odlyzko_zeros3.txt+$D/odlyzko_zeros4.txt+$D/odlyzko_zeros5.txt; do
  for k in 1 4; do
    w=2; [ $k = 4 ] && w=2,3,5,7
    N=$(head -1 $T/band${i}_k$k.txt | cut -d= -f2)
    python3 $R/band_decoder.py $spec $N --M 8 --wheel $w --out $T/py.txt > $T/py.log
    tail -n +2 $T/py.txt > $T/a.txt; tail -n +2 $T/band${i}_k$k.txt > $T/b.txt
    if cmp -s $T/a.txt $T/b.txt; then r="the same $(wc -l < $T/a.txt) primes"; else r=DIFFERENT; fi
    echo "band $i, k = $k: $(sed 's/^.*: \([0-9]* primes decided; first wrong \[[0-9, ]*\]\).*/\1/' $T/py.log); page and band_decoder.py decide $r"
  done
  i=$((i + 1))
done
