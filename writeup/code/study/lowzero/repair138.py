"""Stricter completeness check (round 138's registration) of saved zero sets: N(t) - theta_chi(t)/pi on 200 heights and
the largest unfolded gap; a flagged real-character set is searched again around the first height where the count drops,
with 256 subdivisions, and the missing zeros are added. usage: repair138.py FILE.npz LO HI"""
import sys, math, os
import numpy as np, cypari2
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
f = sys.argv[1]; lo, hi = float(sys.argv[2]), float(sys.argv[3]); Z = dict(np.load(f)); changed = []
def theta(t, q, a): return np.imag(loggamma(0.25 + a / 2 + 0.5j * t)) + t / 2 * math.log(q / math.pi)
def check(g, D):
    q = abs(D); a = 0 if D > 0 else 1; t = np.linspace(lo, hi, 200); g = np.sort(g)
    s = np.searchsorted(g, t) - (theta(t, q, a) - theta(0.0, q, a)) / math.pi
    b = g[(g >= lo) & (g <= hi)]; u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e))
    return s, t, float(np.diff(u).max()), float(b[np.argmax(np.diff(u))])
for k in sorted(Z, key=lambda k: abs(int(k[1:]))):
    D = int(k[1:]); s, t, gap, where = check(Z[k], D)
    if np.abs(s).max() > 1.5 or gap > 3.0:
        g = np.sort(Z[k]); i = int(np.argmax(np.abs(s) > 1.0)); t0 = max(0.5, (t[i - 1] if i else lo) - 6)
        new = [float(x) for x in pari('lfunzeros(lfuncreate(%d),[%f,%f],256)' % (D, min(t0, where - 3), max(t[i], where) + 6))]
        add = [x for x in new if np.min(np.abs(g - x)) > 1e-7]
        Z[k] = np.sort(np.concatenate([g, add])); s2, _, gap2, _ = check(Z[k], D)
        changed.append((D, round(float(np.abs(s).max()), 2), round(gap, 2), [round(x, 4) for x in add], round(float(np.abs(s2).max()), 2), round(gap2, 2)))
for c in changed: print('D = %d: max |N - theta/pi| %.2f, largest gap %.2f; added %s; now %.2f, %.2f' % c)
print('%d sets, %d flagged' % (len(Z), len(changed)))
if any(c[3] for c in changed): np.savez_compressed(f, **Z); print('saved')
