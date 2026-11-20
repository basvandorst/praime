"""Residuals of the Landau expansion against OEIS A164775: r2 = (B/(K x/sqrt L) - 1 - c1/L) L^2 should tend to c2, and the
ratio of B to the four-term formula K x/sqrt L (1 + c1/L + c2/L^2 + c3/L^3).  usage: delange_fit.py stripped.gz out.txt"""
import sys, gzip, math
K, c1, c2, c3 = 0.764223653589221, 0.5819486593, 1.2909972064, 4.87469162949
row = next(l for l in gzip.open(sys.argv[1], 'rt') if l.startswith('A164775 '))
terms = [int(t) for t in row.split(',')[1:] if t.strip()]
out = open(sys.argv[2], 'w')
for n, b in enumerate(terms, start=1):
    if n < 4: continue
    L = n * math.log(10); f = K * 10 ** n / math.sqrt(L)
    r2 = (b / f - 1 - c1 / L) * L * L; four = f * (1 + c1 / L + c2 / L ** 2 + c3 / L ** 3)
    s = '10^%-3d  (B/Landau - 1 - c1/L) L^2 = %.3f   (c2 + c3/L = %.3f)   B / four-term formula = %.6f' % (n, r2, c2 + c3 / L, b / four); print(s); out.write(s + '\n')
