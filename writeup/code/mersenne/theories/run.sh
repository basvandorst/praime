#!/bin/sh
set -e
PY=${PY:-python3}
$PY sieve1e9.py 1e7 s1e7.npz
$PY sieve1e9.py 1e9 s1e9.npz > sieve1e9.txt
$PY hlbytes.py > hlbytes.txt
$PY info.py > info.txt
$PY mlane.py > mlane.txt
$PY t2oos.py > t2oos.txt
$PY gold.py 4e6 > gold.txt
$PY gold8.py 2e6 > gold8.txt
$PY ord2.py 1e7 > ord2.txt
$PY t7form.py > t7form.txt
$PY zrace.py > zrace.txt
$PY goldrace.py 4e6 12 > goldrace12.txt
$PY goldrace.py 4e6 8 > goldrace8.txt
$PY -c "import cypari2; P = cypari2.Pari(); P.allocatemem(4*10**8); [open(f, 'w').write('\\n'.join('%.12f' % float(x) for x in P.lfunzeros(P.lfuncreate(d), 400))) for d, f in ((-4, 'zeros_m4.txt'), (12, 'zeros_12.txt'))]"
$PY goldzeros.py 4e6 > goldzeros.txt
