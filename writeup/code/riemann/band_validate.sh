#!/bin/sh
# Every number of RIEMANN.md sections 12 and 13, from scratch (about 40 minutes on 2 cores).
# Run from the repository root with the virtualenv active.
set -e
T=$(mktemp -d)
D=writeup/data/riemann
R=writeup/code/riemann
BD="python3 $R/band_decoder.py"
ALL=$D/odlyzko_zeros3.txt+$D/odlyzko_zeros4.txt+$D/odlyzko_zeros5.txt
node $R/zeros.js 2000000 $T/z2m.f64
echo "== bands of Odlyzko's zeros, odd candidates"
for z in 3 4 5; do $BD $D/odlyzko_zeros$z.txt 16000 --out $T/p$z.txt | tee -a $T/log.txt; done
$BD $ALL 32000 --out $T/p345.txt
echo "== wheels"
$BD $D/odlyzko_zeros5.txt 16000 --wheel 2,3 --out $T/p5w6.txt
$BD $D/odlyzko_zeros5.txt 16000 --wheel 2,3,5,7 --out $T/p5w210.txt
$BD $ALL 60000 --wheel 2,3,5,7 --out $T/p345w210.txt
echo "== same width, low and high (own zeros)"
for W in 1342 3000 10000; do
  $BD low:$W:$T/z2m.f64 $((W * 8))
  $BD band:500000:$W:$T/z2m.f64 $((W * 6)) | tee -a $T/log.txt
  $BD band:1000000:$W:$T/z2m.f64 $((W * 6)) | tee -a $T/log.txt
done
$BD low:1342:$T/z2m.f64+$D/odlyzko_zeros5.txt 16000
echo "== the phase rule at the reach of each band"
python3 $R/band_phase.py $T/log.txt
echo "== cold start (nothing known below S - 6 scales), and the oracle start (true values from S - 6 scales to S)"
for s in 1300 2000 2700 3300; do $BD $D/odlyzko_zeros5.txt $((s + 600)) --cold $s --stop 1 --out $T/c5_$s.txt; done
for s in 5000 8000 11000; do $BD $ALL $((s + 1000)) --cold $s --stop 1 --out $T/c345_$s.txt; done
for s in 3300 4000 5000; do $BD $D/odlyzko_zeros5.txt $((s + 600)) --cold $s --warm --stop 1; done
echo "== independent checks of the lists (sympy.isprime)"
python3 $R/band_check.py list $T/p3.txt $T/p4.txt $T/p5.txt $T/p345.txt $T/p5w6.txt $T/p5w210.txt $T/p345w210.txt \
  $T/c5_1300.txt $T/c5_2000.txt $T/c5_2700.txt $T/c345_5000.txt $T/c345_8000.txt
echo "== the formula, zero side against prime side"
python3 $R/band_check.py formula $D/odlyzko_zeros5.txt
echo "== the primes check the zeros: Odlyzko's tables and our own zeros"
python3 $R/zeros_vs_primes.py $D/odlyzko_zeros3.txt $D/odlyzko_zeros4.txt $D/odlyzko_zeros5.txt \
  band:2000:1000:$T/z2m.f64 band:20000:1000:$T/z2m.f64 band:200000:1000:$T/z2m.f64 band:1100000:1000:$T/z2m.f64
echo "== H4 at height 2.68e11"
python3 $R/zeros_high.py 2.0,1.6,1.3,1.0 500,2500,4500,6500,8500
rm -r $T
