"""Parameter-free prediction of the real error of the local criterion, from the random-phase model of the zeros above T.
E_rms(n)^2 = (c sqrt n C0)^2 / ((q-1) K0)^2 * sum over all zeros with |gamma| > T of Z(c|gamma|)^2,
the zero count of the family being (q-1)(1/pi) ln(q t/2 pi) dt; worst of N numbers ~ E_rms sqrt(2 ln N)."""
import json, math, glob
import numpy as np
from scipy.integrate import quad
from scipy.special import iv
from kbm import Z as Zkb
rows = []
for f in ('fc3_101.json', 'fc3_101b.json', 'fc3_1009.json', 'fc3_1009c.json', 'fc3_10007.json', 'fc3_10007b.json', 'fc3_20011.json', 'fc3test.json'):
    try: rows += json.load(open(f))
    except FileNotFoundError: pass
X = 10000; m = 4; Nn = 10000
seen = {}
for r in rows: seen[(r['q'], r['beta'])] = r          # later files override (the corrected q = 1009 rerun)
print('q       beta  cT     measured worst/(ln n/4)   predicted worst/(ln n/4)   ratio')
for (q, beta), r in sorted(seen.items()):
    c = r['c']; T = r['T']
    C0 = math.sqrt(2 * math.pi) * beta ** m; K0 = iv(m, beta)
    integ = quad(lambda s: math.log(q * s / (2 * math.pi * c)) * float(Zkb(s, beta, m)) ** 2, c * T, c * T + 400, limit=400)[0]
    n = X
    Erms = (c * math.sqrt(n) * C0 / ((q - 1) * K0)) * math.sqrt((q - 1) / (math.pi * c) * integ)
    worst = Erms * math.sqrt(2 * math.log(Nn)) / (0.25 * math.log(n))
    meas = r['max_err_over_quarter']
    print(f'{q:<7d} {beta:4.0f}  {c * T:5.1f}  {meas:22.2e}   {worst:22.2e}   {meas / worst:6.2f}')
