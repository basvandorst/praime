"""Step 3 (development curves only): the lines of the one-point spectrum G_E(xi) = sum e^{2 i theta(gamma) - i xi gamma} at
xi = ln r, r = a/b, against the recipe amplitude -w c_E rho(r)/sqrt(ab) T/2pi (sharp; Hann: half); off-line control
values of xi.  Output: s3_out.txt."""
import os, math, sys, json
import numpy as np
sys.path.insert(0, '.')
import elib as L
T1, T2 = 20.0, 400.0; T = T2 - T1
C = np.load('curves.npz'); ids = [str(s) for s in C['ids']]; sets = [str(s) for s in C['set']]
dev = [i for i in range(len(ids)) if sets[i] == 'dev']
DEVZ = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r148', 'dev_zeros.npz'))
pr = C['primes']

def rho(i, a, b):
    """rho(a/b) for curve i (a, b coprime)"""
    bad = json.loads(str(C['bad'][i])); out = 1.0
    for n_, sgn in ((a, 1), (b, -1)):
        m = n_; q = 2
        while m > 1:
            if m % q == 0:
                v = 0
                while m % q == 0: m //= q; v += 1
                k = int(np.searchsorted(pr, q))
                out *= float(L.rho_pv(q, np.array([float(C['ap'][i][k])]), np.array([q in bad]), sgn * v)[0])
            q += 1
    return out

R = [(2, 1), (1, 2), (3, 1), (1, 3), (5, 1), (1, 5), (7, 1), (1, 7), (3, 2), (2, 3), (4, 1), (1, 4)]
out = open('s3_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
P('development curves (20): lines of G(xi) at xi = ln(a/b); law -w c rho(a/b)/sqrt(ab) T/2pi; slope over curves (jackknife)')
Gs = {}; Ls = {}
for i in dev:
    g = np.sort(DEVZ[ids[i]]); b = g[(g > T1) & (g <= T2)]
    ph = np.exp(2j * L.theta_E(int(C['N'][i]), b)); x = (b - T1) / T; tap = np.sin(math.pi * x) ** 2
    for (a, bb) in R + [(1, 1)]:
        xi = math.log(a / bb); e = ph * np.exp(-1j * xi * b)
        Gs.setdefault((a, bb), []).append((e.sum(), (e * tap).sum()))
        Ls.setdefault((a, bb), []).append(-C['w'][i] * C['c'][i] * rho(i, a, bb) / math.sqrt(a * bb) * T / (2 * math.pi))
    for xi in (0.35, 0.55, 0.9, 1.3, -0.35, -0.55, -0.9, -1.3):     # control: no line (not ln of a small rational)
        e = ph * np.exp(-1j * xi * b); Gs.setdefault(('ctl', xi), []).append((e.sum(), (e * tap).sum()))
def slope(y, l):
    c = np.dot(l, y) / np.dot(l, l); n = len(y)
    reps = [np.dot(np.delete(l, k), np.delete(y, k)) / np.dot(np.delete(l, k), np.delete(l, k)) for k in range(n)]
    return c, math.sqrt((n - 1) / n * np.sum((np.array(reps) - np.mean(reps)) ** 2))
allG = []; allL = []
for key in [(1, 1)] + R:
    G = np.array(Gs[key]); l = np.array(Ls[key])
    s, se = slope(G[:, 0].real, l); sh, seh = slope(G[:, 1].real, l / 2)
    P('  r = %d/%d: law range [%6.1f, %6.1f], sd %.1f | slope sharp %.3f +- %.3f, Hann %.3f +- %.3f | rms(Re G - law) %.1f, mean Im G %.1f' % (
        key[0], key[1], l.min(), l.max(), l.std(), s, se, sh, seh, np.sqrt(np.mean((G[:, 0].real - l) ** 2)), G[:, 0].imag.mean()))
    if key != (1, 1): allG.append(G[:, 0].real); allL.append(l)
s, se = slope(np.concatenate(allG), np.concatenate(allL))
P('  pooled r != 1 (12 lines x 20 curves): slope %.4f (naive jackknife over cells %.4f)' % (s, se))
P('control frequencies (no line expected): rms Re G sharp %.1f, Hann %.1f; mean Re G %.1f' % (
    np.sqrt(np.mean([np.mean(np.array(v)[:, 0].real ** 2) for k, v in Gs.items() if k[0] == 'ctl'])),
    np.sqrt(np.mean([np.mean(np.array(v)[:, 1].real ** 2) for k, v in Gs.items() if k[0] == 'ctl'])),
    np.mean([np.mean(np.array(v)[:, 0].real) for k, v in Gs.items() if k[0] == 'ctl'])))
out.close()
