"""Rain of the primes: the numbers 2..x that no prime p <= y divides ('dry'), Phi(x, y) - 1, against the forecast of
independence x prod_{p <= y} (1 - 1/p) (Mertens), and against Buchstab: the ratio dry / forecast tends to
e^gamma omega(u), u = ln x / ln y; the refined count is integral_y^x omega(ln t / ln y) / ln y dt.
Predictions written before the run: at y = sqrt(x) the ratio is near e^gamma / 2 = 0.8905 only for the asymptotic
formula; with the refined count the page's ratio is matched; near u = 3 the ratio exceeds 1 (omega(3) > e^-gamma).
usage: rain.py x out.txt"""
import sys, math
import numpy as np
x = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
g = 0.5772156649015329
h = 1e-4; U = np.arange(1, 40 + h / 2, h); om = np.zeros_like(U)
m = U <= 2; om[m] = 1 / U[m]
k1 = int(round(1 / h))
for i in range(len(U)):
    if U[i] <= 2: continue
    om[i] = (U[i - 1] * om[i - 1] + h * (om[i - k1] + om[i - 1 - k1]) / 2) / U[i]
def omega(u): return float(np.interp(u, U, om))
pr('Buchstab omega: omega(2) = %.5f, omega(3) = %.5f, max on [2,4] %.5f at u = %.3f, e^-gamma = %.5f' % (omega(2), omega(3), om[(U >= 2) & (U <= 4)].max(), U[(U >= 2) & (U <= 4)][om[(U >= 2) & (U <= 4)].argmax()], math.exp(-g)))
alive = np.ones(x + 1, bool); alive[:2] = False
sv = np.ones(x + 1, bool); sv[:2] = False
for i in range(2, int(x ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
prod = 1.0; rows = []
for p in P:
    alive[p::p] = False; prod *= 1 - 1 / p
    dry = int(alive.sum()); fc = x * prod; u = math.log(x) / math.log(p)
    rows.append((int(p), dry, fc, u))
    if p * p > x and dry == int(sv[p + 1:].sum()) and p > 3 * x ** .5: break
def refined(y, u):
    if u <= 1: return 0.0
    t = np.exp(np.linspace(math.log(y), math.log(x), 4001)); f = np.array([omega(math.log(v) / math.log(y)) for v in t]) / math.log(y)
    return float(np.sum((f[1:] + f[:-1]) / 2 * np.diff(t)))
pr('x = %d: dry numbers in [2, x] after the clouds up to y' % x)
pr('     y        u      dry   forecast  ratio  e^g*omega(u)  refined/forecast')
pick = sorted(set([r for r in rows if r[0] in (2, 3, 5, 7, 11, 13, 17, 23, 31, 47, 71, 97, 101, 151, 211, 307, 401, 503, 701, 997, 1009, 1499, 1999)] + [min(rows, key=lambda r: abs(r[3] - u0)) for u0 in (1.5, 2, 2.5, 3, 3.5, 4, 5, 6)]))
for p, dry, fc, u in pick:
    pr('%6d %8.3f %8d %10.1f %6.4f %8.4f %12.4f' % (p, u, dry, fc, dry / fc, math.exp(g) * omega(u), refined(p, u) / fc))
r2 = min(rows, key=lambda r: abs(r[3] - 2)); r3 = min(rows, key=lambda r: abs(r[3] - 3))
pr('at the largest y <= sqrt(x): see the row near u = 2; the asymptotic value e^g/2 = %.4f' % (math.exp(g) / 2))
above = [(p, u, dry / fc) for p, dry, fc, u in rows if dry > fc]
pr('clouds y where dry > forecast: %d of %d; first %s; the y of the largest ratio %s' % (len(above), len(rows), above[:3], max(rows, key=lambda r: r[1] / r[2])[:1]))
