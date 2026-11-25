# Localise a zeros-vs-primes discrepancy: G and R with smooth bumps (Hann, width 2h) at many centres.
import os, sys
import numpy as np
import cypari2
from scipy.special import loggamma
SC = os.path.dirname(os.path.abspath(__file__))+'/'
CU = np.load(SC + 'rv_curves.npz')
s = sys.argv[1]; XMAX = int(float(sys.argv[2])); h = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
k = list(CU['ids']).index(s); N = int(CU['N'][k])
pari = cypari2.Pari(); pari.allocatemem(600 * 10**6)
primes = np.array([int(x) for x in pari('primes(primepi(%d))' % XMAX)], dtype=np.int64)
a = [int(x.replace('m', '-')) for x in s.split('_')[1:]]
pari('global(EE); EE = ellinit(%s)' % str(a))
ap = np.array([int(x) for x in pari('my(v=vectorsmall(%d),i=0);forprime(p=2,%d,i++;v[i]=ellap(EE,p));v' % (len(primes), XMAX))], float)
pf = primes.astype(float); lam = ap / np.sqrt(pf); bad = (N % primes == 0); th = np.arccos(np.clip(lam / 2, -1, 1))
ns, Ls = [], []
for kk in range(1, 25):
    m = pf ** kk <= XMAX
    if not m.any():
        break
    coef = lam[m] if kk == 1 else np.where(bad[m], lam[m] ** kk, 2 * np.cos(kk * th[m]))
    ns.append(pf[m] ** kk); Ls.append(coef * np.log(pf[m]))
n = np.concatenate(ns); L = np.concatenate(Ls)
t = 2 * np.pi * np.sqrt(n / N)
for it in range(6):
    g = np.log(N * t * t / (4 * np.pi ** 2)) + 1 / (6 * t * t) + 1 / (60 * t ** 4) - np.log(n)
    t = t - g / (2 / t - 1 / (3 * t ** 3) - 1 / (15 * t ** 5))
ok = t > 3
n, L, t = n[ok], L[ok], t[ok]
f = t * np.log(N) - 2 * t * np.log(2 * np.pi) + 2 * loggamma(1 + 1j * t).imag - t * np.log(n)
term = -(1 / (2 * np.pi)) * L / np.sqrt(n) * np.sqrt(2 * np.pi / (2 / t - 1 / (3 * t ** 3))) * np.exp(1j * (f + np.pi / 4))
z = np.sort(np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r144/fresh_zeros.npz')[s])
e2 = np.exp(1j * (z * np.log(N) - 2 * z * np.log(2 * np.pi) + 2 * loggamma(1 + 1j * z).imag))
print(s, 'N', N, 'zeros', len(z), 'max t from primes %.1f' % t.max())
for c0 in np.arange(20 + h, 401 - h, h / 2):
    wz = np.where(np.abs(z - c0) < h, np.cos(np.pi * (z - c0) / (2 * h)) ** 2, 0)
    wp = np.where(np.abs(t - c0) < h, np.cos(np.pi * (t - c0) / (2 * h)) ** 2, 0)
    G = np.sum(e2 * wz); R = np.sum(term * wp)
    flag = ' <==' if abs(G - R) > 0.3 else ''
    print('  centre %6.1f: zeros %8.3f%+8.3fi primes %8.3f%+8.3fi  |diff| %.3f%s' % (c0, G.real, G.imag, R.real, R.imag, abs(G - R), flag))
gaps = np.diff(z); i = np.argsort(-gaps)[:5]
print('largest gaps (t, gap, mean spacing there):', [(round(z[j], 3), round(gaps[j], 3), round(2 * np.pi / np.log(N * z[j] ** 2 / (4 * np.pi ** 2)), 3)) for j in i])
print('smallest gaps:', [(round(z[j], 4), float('%.2e' % gaps[j])) for j in np.argsort(gaps)[:5]])
