"""Decomposition of the family mean of Re(G - pred) at low T (sharp window (0,T]) by the explicit formula:
 G = E0 + (primes part) + (prime squares + 4th powers) ; E0 = (1/2 pi i)(e^{2i theta(T)} - 1) exactly per curve;
 prime squares/4th powers: family mean from r3 (USp(4) means); primes part = remainder (zero side - E0 - squares)."""
import math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
rows = [r for r in load() if r['cls'] == 'g2']
K = np.load('r3_K.npz')['g2']; KT = [5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
for T in (5, 10, 15, 20, 30, 60, 150):
    zs = np.array([(np.exp(2j * theta(r['g'][(r['g'] > 0) & (r['g'] <= T)], r['Q'], r['mu'])).sum() + r['eps'] * r['c'] * T / (2 * math.pi)) for r in rows])
    E0 = np.array([(np.exp(2j * theta(T, r['Q'], r['mu'])) - 1) / (2j * math.pi) for r in rows])
    sq = K[:, KT.index(T)].mean()
    pr = zs - E0 - sq
    print('T=%3d  residual %+5.2f%+5.2fi (jk %.2f)  E0 %+5.2f%+5.2fi  squares %+5.2f%+5.2fi  => primes part (minus law) %+5.2f%+5.2fi ; eps+1 %+5.2f eps-1 %+5.2f' % (
        T, zs.real.mean(), zs.imag.mean(), jk_mean(zs.real)[1], E0.real.mean(), E0.imag.mean(), sq.real, sq.imag, pr.real.mean(), pr.imag.mean(),
        pr.real[np.array([r['eps'] for r in rows]) > 0].mean(), pr.real[np.array([r['eps'] for r in rows]) < 0].mean()))
