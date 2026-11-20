"""Compare the Landau-Ramanujan formulas with the counts of sums of two squares up to 10^n from OEIS A164775 (offset 1;
up to 10^16, Shiu 1986 and later), read from the local OEIS copy.  usage: tiles_oeis.py stripped.gz out.txt"""
import sys, gzip, math
from mpmath import mp, quad, exp, sqrt, log
mp.dps = 30
K, C = 0.7642236535892206, 0.5819486593
row = next(l for l in gzip.open(sys.argv[1], 'rt') if l.startswith('A164775 '))
terms = [int(t) for t in row.split(',')[1:] if t.strip()]
out = open(sys.argv[2], 'w')
def pr(s): print(s); out.write(s + '\n')
pr('   x        B(x) (A164775)        B / (K x/sqrt ln x)   B / (same * (1 + C/ln x))   B / (K int_2^x dt/sqrt ln t)')
for n, b in enumerate(terms, start=1):
    x = mp.mpf(10) ** n; L = log(x)
    f1 = K * x / sqrt(L); f2 = f1 * (1 + C / L)
    ri = K * quad(lambda u: exp(u) / sqrt(u), [log(2), L / 2, L])
    pr('10^%-3d %20d %16.5f %22.5f %24.5f' % (n, b, float(b / f1), float(b / f2), float(b / ri)))
