"""Independent check of an information limit found by lpchain.py (mpmath and sympy; shares only the budget delta).

lpchain.py stops with 'ambiguous' when, for the next unknown b_j, it finds two 0/1 patterns, one with b_j = 0 and one
with b_j = 1, that both satisfy every sample of the window. It writes them to <run>_ambiguous.json. This script
rebuilds the window in mpmath (zero side from the zero strings, prime side exact, kernel from mpmath's Bessel
functions) and checks that both patterns satisfy |A b - y| <= delta in every row. If they do, the samples cannot
decide b_j: the data, not the search, has run out.
Usage: ZFILE=... EPS=... info_law_check.py run_ambiguous.json
"""
import json, math, os, sys
import mpmath as mp
import numpy as np
from sympy import factorint, isprime
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cert import sample   # the rigorous budget delta

mp.mp.dps = 30
ZS = [l.strip() for l in open(os.environ.get('ZFILE', os.path.join(HERE, '../../../data/riemann/zeros_hp2419.txt'))) if l.strip()]
Zm = [mp.mpf(s) for s in ZS]
Gld = np.array(ZS, dtype=np.longdouble)
EPS = float(os.environ.get('EPS', '2e-16'))

def Kmp(v, a, beta, m):
    t = 1 - (v / a) ** 2
    return t ** (mp.mpf(m) / 2) * mp.besseli(m, beta * mp.sqrt(t)) if t > 0 else mp.mpf(0)

def Khmp(w, a, beta, m):
    nu = m + mp.mpf(1) / 2; s = abs(a * w)
    if s < beta:
        q = mp.sqrt(beta ** 2 - s ** 2); z = mp.besseli(nu, q) / q ** nu if q > 0 else 1 / (2 ** nu * mp.gamma(nu + 1))
    else:
        y = mp.sqrt(s ** 2 - beta ** 2); z = mp.besselj(nu, y) / y ** nu if y > 0 else 1 / (2 ** nu * mp.gamma(nu + 1))
    return a * mp.sqrt(2 * mp.pi) * beta ** m * z

def Pmp(x, a, beta, m):
    nu = m + mp.mpf(1) / 2; r = a / (2 * x); q = mp.sqrt(beta ** 2 + r ** 2)
    return a * mp.sqrt(2 * mp.pi) * beta ** m * mp.besseli(nu, q) / q ** nu

def lam(n):
    f = factorint(n); return mp.log(next(iter(f))) if len(f) == 1 else mp.mpf(0)

def zero_side(x, a, beta, m, T):
    x = mp.mpf(x); L = mp.log(x); s = mp.mpf(0)
    for g in Zm:
        if g > T: break
        s += mp.cos(g * L) * Khmp(g / x, a, beta, m)
    triv = mp.quad(lambda v: Kmp(v, a, beta, m) * mp.sqrt(x * mp.e ** (v / x)) / (x * ((x * mp.e ** (v / x)) ** 3 - x * mp.e ** (v / x))), [-a, 0, a])
    return Pmp(x, a, beta, m) - 2 / mp.sqrt(x) * s - mp.sqrt(x) * triv

R = json.load(open(sys.argv[1]))
U, xs, j = R['U'], R['xs'], R['bit']
a, beta, m, cc, k, h = R['params']
fits = {int(v): p for v, p in R['fits'].items()}
k0 = mp.besseli(m, beta)
worst = {0: mp.mpf(-1), 1: mp.mpf(-1)}; ok = {0: True, 1: True}
for x in xs:
    T = cc * beta * x / a
    xm = mp.mpf(x)
    lo, hi = int(mp.floor(xm * mp.e ** (-a / xm))), int(mp.ceil(xm * mp.e ** (a / xm)))
    known = sum(lam(n) * mp.sqrt(xm / n) * Kmp(xm * mp.log(mp.mpf(n) / xm), a, beta, m) for n in range(max(2, lo), hi + 1) if n not in U)
    y = (zero_side(x, a, beta, m, T) - known) / k0
    d = mp.mpf(sample(x, Gld, T, a, beta, m, EPS)[1]) / k0
    col = [mp.log(u) * mp.sqrt(xm / u) * Kmp(xm * mp.log(mp.mpf(u) / xm), a, beta, m) / k0 for u in U]
    for v, p in fits.items():
        r = abs(sum(c for c, b in zip(col, p) if b) - y)
        worst[v] = max(worst[v], r / d)
        if r > d: ok[v] = False
truth = [1 if isprime(u) else 0 for u in U]
diff = [u for u, p0, p1 in zip(U, fits[0], fits[1]) if p0 != p1]
print(f'window at {U[0]} (a = {a}, beta = {beta}, m = {m}, height at x: {cc * beta / a:.4f} x, {len(xs)} samples, {len(U)} unknowns), bit {j} = {U[j]} (prime: {isprime(U[j])})')
print(f'pattern with b = 0: fits every sample: {ok[0]} (largest |A b - y| / delta = {mp.nstr(worst[0], 4)}); equals the truth: {fits[0] == truth}')
print(f'pattern with b = 1: fits every sample: {ok[1]} (largest |A b - y| / delta = {mp.nstr(worst[1], 4)}); equals the truth: {fits[1] == truth}')
print(f'the two patterns differ at {diff}')
print('information limit confirmed' if ok[0] and ok[1] else 'NOT confirmed')
